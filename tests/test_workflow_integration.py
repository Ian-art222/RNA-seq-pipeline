#!/usr/bin/env python3
"""Opt-in, real-tool integration tests using isolated synthetic paired-end BAMs.

Run with a Python containing pysam. Tool overrides are RNASEQ_TEST_SNAKEMAKE,
RNASEQ_TEST_RSEQC, RNASEQ_TEST_STRINGTIE, and RNASEQ_TEST_WORKFLOW_PYTHON.
The workflow Python must contain pandas/openpyxl; it may differ from the test
Python used for BAM validation and RSeQC. No environment installation occurs.
Set RNASEQ_TEST_KEEP=1 to retain temporary workdirs and command logs.
"""

import csv
import gzip
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import unittest

try:
    import pysam
except ImportError:
    pysam = None


ROOT = Path(__file__).resolve().parents[1]


def executable(variable, name):
    candidate = os.environ.get(variable) or shutil.which(name)
    return str(Path(candidate).resolve()) if candidate else None


class WorkflowIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.snakemake = executable("RNASEQ_TEST_SNAKEMAKE", "snakemake")
        cls.rseqc = executable("RNASEQ_TEST_RSEQC", "infer_experiment.py")
        cls.stringtie = executable("RNASEQ_TEST_STRINGTIE", "stringtie")
        cls.workflow_python = os.environ.get("RNASEQ_TEST_WORKFLOW_PYTHON", sys.executable)
        if pysam is None or not all((cls.snakemake, cls.rseqc, cls.stringtie)):
            raise unittest.SkipTest("Requires pysam and existing Snakemake, RSeQC, and StringTie executables")
        cls.root = Path(tempfile.mkdtemp(prefix="rnaseq-strand-integration-"))
        print("Synthetic integration workdir: {}".format(cls.root), flush=True)
        bindir = cls.root / "bin"
        bindir.mkdir()
        for name, target in (("infer_experiment.py", cls.rseqc), ("stringtie", cls.stringtie)):
            (bindir / name).symlink_to(target)
        # Mirror the workflow's separate quantification/strand environments
        # without asking Snakemake to create or download a conda environment.
        dispatcher = bindir / "python"
        dispatcher.write_text(
            "#!/bin/sh\ncase \"$1\" in\n"
            "  */strandness.py|*/prepare_strand_bed.py|*/validate_bam_reference.py) exec {} \"$@\" ;;\n"
            "  *) exec {} \"$@\" ;;\nesac\n".format(
                shlex.quote(sys.executable), shlex.quote(cls.workflow_python)
            )
        )
        dispatcher.chmod(0o755)
        cls.env = dict(os.environ)
        cls.env["PATH"] = str(bindir) + os.pathsep + cls.env.get("PATH", "")
        cls.env["PYTHONDONTWRITEBYTECODE"] = "1"
        cls.env["OMP_NUM_THREADS"] = "1"

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, "root") and os.environ.get("RNASEQ_TEST_KEEP") != "1":
            shutil.rmtree(cls.root)

    def run_command(self, command, cwd, timeout=180):
        process = subprocess.Popen(
            command, cwd=cwd, env=self.env, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, text=True, start_new_session=True,
        )
        try:
            output, _ = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
            output, _ = process.communicate()
            self.fail("Command timed out: {}\n{}".format(shlex.join(command), output[-10000:]))
        return process.returncode, output

    def make_case(self, name, samples, stranded):
        workdir = self.root / name
        workdir.mkdir()
        annotation = workdir / "reference.gtf"
        rows = []
        for gene, strand, start in (("gene_plus", "+", 100), ("gene_minus", "-", 2000)):
            attributes = 'gene_id "{}"; transcript_id "{}.t1";'.format(gene, gene)
            for feature in ("transcript", "exon"):
                rows.append("chrSynthetic\tfixture\t{}\t{}\t{}\t.\t{}\t.\t{}\n".format(
                    feature, start + 1, start + 900, strand, attributes
                ))
        annotation.write_text("".join(rows))
        for i in range(1, 9):
            index_file = workdir / f"unused_index.{i}.ht2"
            index_file.write_text("fixture only: these tests start from BAM\n")
            os.utime(index_file, (time.time() - 120, time.time() - 120))
        for sample, direction in samples.items():
            self.make_sample(workdir, sample, direction)
        config = {
            "stranded": stranded,
            "strandness": {"sample_size": 200000, "mapq": 30, "min_determinable": 0.5,
                           "min_dominance": 0.8, "bed12": ""},
            "paths": {"clean_dir": "cleandata", "align_dir": "alignment",
                      "quant_dir": "quantification", "logs_dir": "logs"},
            "reference": {"hisat2_index_prefix": str(workdir / "unused_index"),
                          "annotation": str(annotation)},
            "sample_scan": {"r1_suffix": "_1.clean.fastq.gz", "r2_suffix": "_2.clean.fastq.gz", "strict": True},
            "alignment": {"mapq_filter": 0, "hisat2_extra": ""},
            "stringtie": {"use_ballgown": False, "extra": ""},
            "threads": {"hisat2": 1, "samtools_view": 1, "samtools_sort": 1,
                        "samtools_index": 1, "stringtie": 1},
            "output": {"xlsx": "quantification/gene_expression_matrices.xlsx"},
        }
        # JSON is valid YAML and avoids an additional test-only PyYAML dependency.
        (workdir / "config.yaml").write_text(json.dumps(config, indent=2))
        return workdir

    def make_sample(self, workdir, sample, direction):
        clean = workdir / "cleandata" / sample
        clean.mkdir(parents=True, exist_ok=True)
        old_time = time.time() - 120
        for mate in (1, 2):
            fastq = clean / "{}_{}.clean.fastq.gz".format(sample, mate)
            with gzip.open(fastq, "wt") as handle:
                handle.write("@placeholder/{}\n{}\n+\n{}\n".format(mate, "A" * 75, "I" * 75))
            os.utime(fastq, (old_time, old_time))
        # Imported BAM fixture: measure its FASTQ fixture before publishing BAM.
        import importlib.util
        spec = importlib.util.spec_from_file_location("read_lengths", ROOT / "workflow/scripts/read_lengths.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        report = module.measure(sample, clean / (sample + "_1.clean.fastq.gz"), clean / (sample + "_2.clean.fastq.gz"))
        module.atomic_json(workdir / "alignment/read_lengths" / (sample + ".json"), report)
        directory = workdir / "alignment" / sample
        directory.mkdir(parents=True, exist_ok=True)
        unsorted = directory / "synthetic.unsorted.bam"
        bam_path = directory / (sample + ".sorted.bam")
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": "chrSynthetic", "LN": 5000}],
                  "PG": [{"ID": "synthetic_fixture", "PN": "synthetic_fixture"}]}
        # 600 fragments per isolated gene, 2 genes: 2400 informative alignments.
        with pysam.AlignmentFile(str(unsorted), "wb", header=header) as bam:
            for gene, strand, gene_start in (("gene_plus", "+", 100), ("gene_minus", "-", 2000)):
                for number in range(600):
                    forward = direction == "fr" or (direction == "balanced" and number % 2 == 0)
                    r1_plus = (strand == "+") == forward
                    left = gene_start + number % 500
                    right = left + 200
                    for mate in (1, 2):
                        read = pysam.AlignedSegment()
                        read.query_name = "{}_{}_{}".format(sample, gene, number)
                        read.query_sequence = "A" * 75
                        read.query_qualities = pysam.qualitystring_to_array("I" * 75)
                        read.flag = (99 if r1_plus else 83) if mate == 1 else (147 if r1_plus else 163)
                        read.reference_id = read.next_reference_id = 0
                        plus = r1_plus if mate == 1 else not r1_plus
                        read.reference_start = left if plus else right
                        read.next_reference_start = right if plus else left
                        read.template_length = 275 if plus else -275
                        read.mapping_quality = 60
                        read.cigarstring = "75M"
                        read.set_tags([("NH", 1), ("NM", 0)])
                        bam.write(read)
        pysam.sort("-@", "1", "-o", str(bam_path), str(unsorted))
        unsorted.unlink()
        pysam.index(str(bam_path))
        (directory / (sample + ".align_summary.txt")).write_text(
            "1200 reads; of these:\n"
            "  1200 (100.00%) were paired; of these:\n"
            "    0 (0.00%) aligned concordantly 0 times\n"
            "    1200 (100.00%) aligned concordantly exactly 1 time\n"
            "    0 (0.00%) aligned concordantly >1 times\n"
            "100.00% overall alignment rate\n"
        )

    def run_workflow(self, workdir, dry_run=False, override=None):
        command = [self.snakemake, "--snakefile", str(ROOT / "workflow" / "Snakefile"),
                   "--directory", str(workdir), "--configfile", str(workdir / "config.yaml"),
                   "--cores", "2", "--printshellcmds"]
        if dry_run:
            command.append("--dry-run")
        if override is not None:
            command.extend(["--config", "stranded=" + override])
        status, output = self.run_command(command, cwd=workdir)
        (workdir / ("snakemake.dry-run.log" if dry_run else "snakemake.run.log")).write_text(output)
        self.assertNotIn("localrule hisat2_align_sort:", output, output[-10000:])
        return status, output

    def assert_completed_matrices(self, workdir, samples, stranded):
        quant = workdir / "quantification"
        if stranded:
            quant /= "stranded"
        for sample in samples:
            gtf = quant / sample / (sample + ".gtf")
            self.assertIn('\ttranscript\t', gtf.read_text())
        for filename in ("gene_count_matrix.csv", "gene_fpkm_matrix.csv", "gene_tpm_matrix.csv"):
            with (quant / filename).open() as handle:
                table = list(csv.reader(handle))
            self.assertEqual(table[0][1:], sorted(samples))
            self.assertEqual({row[0] for row in table[1:]}, {"gene_plus", "gene_minus"})
            self.assertTrue(all(float(value) > 0 for row in table[1:] for value in row[1:]))
        with (quant / "gene_tpm_matrix.csv").open() as handle:
            table = list(csv.reader(handle))
        for index in range(len(samples)):
            self.assertAlmostEqual(sum(float(row[index + 1]) for row in table[1:]), 1e6, delta=10)
        workbook = quant / "gene_expression_matrices.xlsx"
        command = [self.workflow_python, "-c",
                   "import openpyxl,sys; w=openpyxl.load_workbook(sys.argv[1],read_only=True); "
                   "assert w.sheetnames==['alignment_report','gene_count','gene_fpkm','gene_tpm']; "
                   "assert w['gene_count'].max_row==3", str(workbook)]
        status, output = self.run_command(command, workdir, timeout=30)
        self.assertEqual(status, 0, output)

    def test_dry_run_selects_only_requested_branch(self):
        for enabled in (False, True):
            with self.subTest(stranded=enabled):
                workdir = self.make_case("dry-{}".format(enabled), {"sample": "fr"}, not enabled)
                status, output = self.run_workflow(workdir, dry_run=True, override=str(enabled).lower())
                self.assertEqual(status, 0, output[-10000:])
                self.assertIn("prepare_reference", output)
                for rule in ("infer_strandness", "strandness_gate"):
                    if enabled:
                        self.assertIn(rule, output)
                    else:
                        self.assertNotIn(rule, output)
                self.assertEqual((workdir / "quantification" / "strandness").exists(), False)

    def test_invalid_config_and_direction_overrides_are_rejected(self):
        cases = (
            ("invalid-switch", None, None, "maybe", "stranded must be true or false"),
            ("conflicting-hisat", "alignment", "hisat2_extra", "--rna-strandness RF", "remove --rna-strandness"),
            ("conflicting-stringtie", "stringtie", "extra", "--rf", "remove --fr/--rf"),
        )
        for name, group, field, value, expected in cases:
            with self.subTest(case=name):
                workdir = self.make_case(name, {"sample": "fr"}, True)
                override = value if group is None else None
                if group is not None:
                    config_path = workdir / "config.yaml"
                    config = json.loads(config_path.read_text())
                    config[group][field] = value
                    config_path.write_text(json.dumps(config))
                status, output = self.run_workflow(workdir, dry_run=True, override=override)
                self.assertNotEqual(status, 0)
                self.assertIn(expected, output)

    def test_unstranded_preserves_original_outputs(self):
        workdir = self.make_case("unstranded", {"ordinary": "balanced"}, False)
        status, output = self.run_workflow(workdir)
        self.assertEqual(status, 0, output[-10000:])
        self.assert_completed_matrices(workdir, ["ordinary"], False)
        self.assertFalse((workdir / "quantification" / "stranded").exists())
        self.assertNotIn("infer_strandness", output)

    def test_fr_and_rf_run_real_inference_and_quantification(self):
        for direction in ("fr", "rf"):
            with self.subTest(direction=direction):
                workdir = self.make_case("success-" + direction, {"sample": direction}, True)
                expected_reads = 2400
                if direction == "fr":
                    config_path = workdir / "config.yaml"
                    config = json.loads(config_path.read_text())
                    config["strandness"]["sample_size"] = 1000
                    config_path.write_text(json.dumps(config))
                    expected_reads = 1000
                status, output = self.run_workflow(workdir)
                self.assertEqual(status, 0, output[-10000:])
                self.assert_completed_matrices(workdir, ["sample"], True)
                quant = workdir / "quantification" / "stranded"
                report = json.loads((quant / "strandness" / "sample.json").read_text())
                self.assertEqual(report["status"], "PASS")
                self.assertEqual(report["stringtie_flag"], "--" + direction)
                self.assertEqual(report["usable_reads"], expected_reads)
                self.assertAlmostEqual(report["dominant_fraction"], 1.0)
                self.assertIn("STRINGTIE_STRAND=--" + direction,
                              (quant / "strandness" / "sample.env").read_text())
                gtf = (quant / "sample" / "sample.gtf").read_text()
                self.assertIn(" --" + direction + " ", gtf)
                self.assertFalse((workdir / "quantification" / "gene_count_matrix.csv").exists())

    def read_cohort_diagnostics(self, workdir):
        diagnostic = workdir / "logs" / "stranded" / "strandness" / "strandness_diagnostic.tsv"
        with diagnostic.open() as handle:
            return {row["sample"]: row for row in csv.DictReader(handle, delimiter="\t")}

    def test_opposite_protocols_use_each_samples_own_direction(self):
        samples = {"forward": "fr", "reverse": "rf"}
        workdir = self.make_case("mixed-protocols", samples, True)
        status, output = self.run_workflow(workdir)
        self.assertEqual(status, 0, output[-10000:])
        self.assert_completed_matrices(workdir, samples, True)
        quant = workdir / "quantification" / "stranded"
        diagnostics = self.read_cohort_diagnostics(workdir)
        for sample, direction in samples.items():
            with self.subTest(sample=sample):
                env = quant / "strandness" / (sample + ".env")
                self.assertEqual(env.read_text().strip(), "STRINGTIE_STRAND=--" + direction)
                gtf = (quant / sample / (sample + ".gtf")).read_text()
                self.assertIn(" --" + direction + " ", gtf)
                self.assertEqual(diagnostics[sample]["status"], "PASS")
                self.assertEqual(diagnostics[sample]["stringtie_flag"], "--" + direction)
        self.assertFalse((quant / "strandness" / "strandness.env").exists())

    def test_ambiguous_sample_does_not_block_good_sample_with_default_profile(self):
        # Do not pass --keep-going: the workflow's default profile must keep
        # independent work running when a different sample's gate fails.
        workdir = self.make_case("ambiguous", {"good": "fr", "balanced": "balanced"}, True)
        status, output = self.run_workflow(workdir)
        self.assertNotEqual(status, 0, output[-10000:])
        quant = workdir / "quantification" / "stranded"
        self.assertIn(" --fr ", (quant / "good" / "good.gtf").read_text())
        self.assertTrue((quant / "good" / "good_abundance.tab").is_file())
        self.assertIn("STRINGTIE_STRAND=--fr", (quant / "strandness" / "good.env").read_text())
        self.assertFalse((quant / "strandness" / "balanced.env").exists())
        self.assertFalse((quant / "balanced" / "balanced.gtf").exists())
        self.assertFalse((quant / "gene_count_matrix.csv").exists())
        self.assertFalse((quant / "gene_expression_matrices.xlsx").exists())
        report = json.loads((quant / "strandness" / "balanced.json").read_text())
        self.assertEqual(report["status"], "FAIL")
        self.assertAlmostEqual(report["dominant_fraction"], 0.5)
        for sample in ("good", "balanced"):
            logdir = workdir / "logs" / "stranded" / "strandness" / sample
            for filename in ("stdout.txt", "stderr.txt", "gate.log", "strandness_diagnostic.tsv"):
                self.assertTrue((logdir / filename).is_file(), str(logdir / filename))
        diagnostics = self.read_cohort_diagnostics(workdir)
        self.assertEqual(diagnostics["good"]["status"], "PASS")
        self.assertEqual(diagnostics["balanced"]["status"], "FAIL")

    def test_changed_sample_does_not_rerun_completed_other_sample(self):
        workdir = self.make_case("resume-isolated", {"good": "fr", "changed": "fr"}, True)
        status, output = self.run_workflow(workdir)
        self.assertEqual(status, 0, output[-10000:])
        self.assert_completed_matrices(workdir, ["good", "changed"], True)
        quant = workdir / "quantification" / "stranded"
        preserved = [
            quant / "good" / "good.gtf",
            quant / "good" / "good_abundance.tab",
            quant / "strandness" / "good.json",
            quant / "strandness" / "good.env",
        ]
        original = {path: (path.stat().st_mtime_ns, path.read_bytes()) for path in preserved}
        self.make_sample(workdir, "changed", "balanced")
        status, output = self.run_workflow(workdir)
        self.assertNotEqual(status, 0, output[-10000:])
        changed = json.loads((quant / "strandness" / "changed.json").read_text())
        self.assertEqual(changed["status"], "FAIL")
        self.assertFalse((quant / "strandness" / "changed.env").exists())
        for path in preserved:
            with self.subTest(file=path.name):
                self.assertEqual((path.stat().st_mtime_ns, path.read_bytes()), original[path])
        diagnostics = self.read_cohort_diagnostics(workdir)
        self.assertEqual(diagnostics["good"]["status"], "PASS")
        self.assertEqual(diagnostics["changed"]["status"], "FAIL")


if __name__ == "__main__":
    unittest.main(verbosity=2)
