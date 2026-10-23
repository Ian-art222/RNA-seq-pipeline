#!/usr/bin/env python3
"""Opt-in FASTQ-to-Excel regressions using installed tools and isolated copies.

Enable with RNASEQ_TEST_FULL_FASTQ=1. Run with a Python containing pysam.
Tools are found on PATH or through RNASEQ_TEST_{SNAKEMAKE,HISAT2,HISAT2_BUILD,
SAMTOOLS,STRINGTIE,RSEQC,GFFREAD}. RNASEQ_TEST_WORKFLOW_PYTHON selects the
pandas/openpyxl interpreter. RNASEQ_TEST_KEEP=1 retains all fixtures and logs.
No installation, download, production source mutation, or biological data is used.
"""
import csv
import gzip
import json
import math
import os
from pathlib import Path
import random
import re
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
SAMPLES = {"001": 75, "NA": 150}
GENES = (("001", "FirstGene", 1000, "+"), ("Gene|literal", "SecondGene", 5000, "-"))
CTABS = ("e2t.ctab", "e_data.ctab", "i2t.ctab", "i_data.ctab", "t_data.ctab")


def executable(variable, name):
    candidate = os.environ.get(variable) or shutil.which(name)
    return str(Path(candidate).resolve()) if candidate else None


def execute(command, cwd, env, log, timeout=180):
    """Retain exact argv, combined output and timing, also for failures."""
    started = time.monotonic()
    process = subprocess.Popen(command, cwd=cwd, env=env, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True, start_new_session=True)
    try:
        output, _ = process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        output, _ = process.communicate()
        output += "\nTEST COMMAND TIMED OUT\n"
    Path(log).write_text(output)
    Path(str(log) + ".command.json").write_text(json.dumps({
        "argv": command, "cwd": str(cwd), "returncode": process.returncode,
        "seconds": time.monotonic() - started,
    }, indent=2) + "\n")
    return process.returncode, output


@unittest.skipUnless(os.environ.get("RNASEQ_TEST_FULL_FASTQ") == "1", "Set RNASEQ_TEST_FULL_FASTQ=1 for installed-tool FASTQ regressions")
class FullFastqRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        names = {"snakemake": "SNAKEMAKE", "hisat2": "HISAT2", "hisat2-build": "HISAT2_BUILD",
                 "samtools": "SAMTOOLS", "stringtie": "STRINGTIE",
                 "infer_experiment.py": "RSEQC", "gffread": "GFFREAD"}
        cls.tools = {name: executable("RNASEQ_TEST_" + variable, name) for name, variable in names.items()}
        if pysam is None or not all(cls.tools.values()):
            raise unittest.SkipTest("Requires pysam and installed tools: " + ", ".join(k for k, v in cls.tools.items() if not v))
        cls.workflow_python = os.environ.get("RNASEQ_TEST_WORKFLOW_PYTHON", sys.executable)
        cls.root = Path(tempfile.mkdtemp(prefix="rnaseq-fastq-regressions-"))
        print("FASTQ regression artifacts: " + str(cls.root), flush=True)
        bindir = cls.root / "bin"
        bindir.mkdir()
        for name, target in cls.tools.items():
            (bindir / name).symlink_to(target)
        dispatcher = bindir / "python"
        dispatcher.write_text(
            "#!/bin/sh\ncase \"$1\" in\n"
            "  */strandness.py|*/prepare_strand_bed.py|*/validate_bam_reference.py) exec "
            + shlex.quote(sys.executable) + " \"$@\" ;;\n"
            "  *) exec " + shlex.quote(cls.workflow_python) + " \"$@\" ;;\nesac\n"
        )
        dispatcher.chmod(0o755)
        cls.env = dict(os.environ)
        cls.env.update(PATH=str(bindir) + os.pathsep + cls.env.get("PATH", ""),
                       PYTHONDONTWRITEBYTECODE="1", OMP_NUM_THREADS="1")
        rng = random.Random(20261022)
        cls.genome = "".join(rng.choices("ACGT", k=10000))
        cls.reference = cls.root / "reference"
        cls.reference.mkdir()
        fasta = cls.reference / "mini.fa"
        fasta.write_text(">chrMini\n" + cls.genome + "\n")
        code, output = execute([cls.tools["hisat2-build"], "-p", "2", str(fasta), str(cls.reference / "mini_index")],
                               cls.root, cls.env, cls.root / "hisat2-build.log")
        if code:
            raise RuntimeError("Synthetic HISAT2 index construction failed: " + output[-8000:])

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, "root") and os.environ.get("RNASEQ_TEST_KEEP") != "1":
            shutil.rmtree(cls.root)

    def make_case(self, name, stranded, samples=None):
        work = self.root / name
        work.mkdir()
        source = work / "pipeline source"
        shutil.copytree(ROOT, source, ignore=shutil.ignore_patterns("__pycache__", ".git", ".snakemake"))
        references = work / "references"
        shutil.copytree(self.reference, references)
        annotation = references / "reference.gtf"
        rows = []
        for gene, symbol, start, strand in GENES:
            attrs = f'gene_id "{gene}"; transcript_id "{gene}.t1"; gene_name "{symbol}";'
            for feature in ("transcript", "exon"):
                rows.append(f"chrMini\tfixture\t{feature}\t{start + 1}\t{start + 1400}\t.\t{strand}\t.\t{attrs}\n")
        annotation.write_text("".join(rows))
        for sample, length in (samples or SAMPLES).items():
            self.make_fastq(work, sample, length)
        config = {
            "stranded": stranded,
            "strandness": {"sample_size": 200000, "mapq": 30, "min_determinable": 0.5, "min_dominance": 0.8, "bed12": ""},
            "paths": {"clean_dir": "cleandata", "align_dir": "alignment", "quant_dir": "quantification", "logs_dir": "logs"},
            "reference": {"hisat2_index_prefix": str(references / "mini_index"), "annotation": str(annotation)},
            "sample_scan": {"r1_suffix": "_1.clean.fastq.gz", "r2_suffix": "_2.clean.fastq.gz", "strict": True},
            "alignment": {"mapq_filter": 30, "hisat2_extra": ""},
            "stringtie": {"use_ballgown": True, "extra": ""},
            "threads": {"hisat2": 2, "samtools_view": 1, "samtools_sort": 1, "samtools_index": 2, "stringtie": 2},
            "output": {"xlsx": "quantification/gene_expression_matrices.xlsx"},
        }
        (work / "config.json").write_text(json.dumps(config, indent=2))
        return work, source

    def make_fastq(self, work, sample, length):
        directory = work / "cleandata" / sample
        directory.mkdir(parents=True)
        rng = random.Random(31)
        complement = str.maketrans("ACGT", "TGCA")
        reverse = lambda sequence: sequence.translate(complement)[::-1]
        with gzip.open(directory / (sample + "_1.clean.fastq.gz"), "wt") as r1, gzip.open(directory / (sample + "_2.clean.fastq.gz"), "wt") as r2:
            for gene, _, start, strand in GENES:
                for number in range(100):
                    left = start + rng.randrange(1000)
                    right = left + 200
                    first, second = self.genome[left:left + length], reverse(self.genome[right:right + length])
                    if strand == "-":
                        first, second = second, first
                    for mate, handle, sequence in ((1, r1, first), (2, r2, second)):
                        handle.write(f"@{sample}_{gene}_{number}/{mate}\n{sequence}\n+\n" + "I" * length + "\n")

    def run_workflow(self, work, source, cores, name, dry=False):
        command = [self.tools["snakemake"], "--snakefile", str(source / "workflow" / "Snakefile"),
                   "--directory", str(work), "--configfile", str(work / "config.json"),
                   "--cores", str(cores), "--printshellcmds"]
        if dry:
            command.append("--dry-run")
        return execute(command, work, self.env, work / (name + ".log"))

    @staticmethod
    def table(path):
        with path.open(newline="") as handle:
            rows = list(csv.reader(handle))
        return rows[0], {row[0]: row[1:] for row in rows[1:]}

    def assert_outputs(self, work, stranded, cores, ballgown=True):
        quant = work / "quantification"
        if stranded:
            quant /= "stranded"
        expected_genes = {gene for gene, _, _, _ in GENES}
        raw_header, raw = self.table(quant / "gene_count_matrix.raw.csv")
        final_header, final = self.table(quant / "gene_count_matrix.csv")
        self.assertEqual(raw_header[1:], list(SAMPLES))
        self.assertEqual(final_header, ["GeneID"] + list(SAMPLES))
        self.assertEqual(set(raw), expected_genes)
        self.assertEqual(raw, final)
        self.assertEqual(json.loads((work / "alignment" / "read_lengths.json").read_text()), SAMPLES)
        for sample_index, (sample, length) in enumerate(SAMPLES.items()):
            lengths = json.loads((work / "alignment" / "read_lengths" / (sample + ".json")).read_text())
            self.assertEqual(lengths["read_pairs"], 200)
            self.assertEqual(lengths["mean_read_length"], length)
            bam = work / "alignment" / sample / (sample + ".sorted.bam")
            with pysam.AlignmentFile(str(bam), "rb") as handle:
                self.assertEqual(sum(1 for _ in handle), 400)
            resource_lines = [line for line in (work / "logs" / sample / "hisat2_align_sort.log").read_text().splitlines() if line.startswith("[RESOURCES] ")]
            self.assertEqual(len(resource_lines), 1)
            budget = json.loads(resource_lines[0].split(" ", 1)[1])
            self.assertEqual(budget["budget"], cores)
            self.assertEqual(budget["mode"], "sequential" if cores == 1 else "pipeline")
            totals = [budget[key] for key in ("hisat2_total", "view_total", "sort_total")]
            self.assertLessEqual(max(totals) if cores == 1 else sum(totals), cores)
            # Independently recompute the published coverage/exon-length formula.
            coverage, exon_lengths, gene_ids = {}, {}, {}
            for line in (quant / sample / (sample + ".gtf")).read_text().splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.split("\t")
                attrs = dict(re.findall(r'(\w+) "([^\"]*)"', fields[8]))
                tid = attrs["transcript_id"]
                gene_ids[tid] = attrs["gene_id"]
                if fields[2] == "transcript":
                    coverage[tid] = float(attrs["cov"])
                elif fields[2] == "exon":
                    exon_lengths[tid] = exon_lengths.get(tid, 0) + int(fields[4]) - int(fields[3]) + 1
            expected = {gene: 0 for gene in expected_genes}
            fixed_75 = {gene: 0 for gene in expected_genes}
            for tid, cov in coverage.items():
                expected[gene_ids[tid]] += math.ceil(cov * exon_lengths[tid] / length)
                fixed_75[gene_ids[tid]] += math.ceil(cov * exon_lengths[tid] / 75)
            observed = {gene: int(final[gene][sample_index]) for gene in expected_genes}
            self.assertEqual(observed, expected)
            self.assertTrue(all(value > 0 for value in observed.values()))
            if length == 150:
                self.assertNotEqual(observed, fixed_75)
            for name in CTABS:
                self.assertEqual((quant / sample / name).is_file(), ballgown)
            validation = json.loads((quant / "reference" / "bam_validation" / (sample + ".json")).read_text())
            self.assertEqual(validation["status"], "PASS")
            if stranded:
                report = json.loads((quant / "strandness" / (sample + ".json")).read_text())
                self.assertEqual(report["status"], "PASS")
                self.assertEqual(report["stringtie_flag"], "--fr")
        _, tpm = self.table(quant / "gene_tpm_matrix.csv")
        for column in range(len(SAMPLES)):
            self.assertAlmostEqual(sum(float(row[column]) for row in tpm.values()), 1e6, delta=10)
        code = (
            "import openpyxl,sys; w=openpyxl.load_workbook(sys.argv[1],read_only=True); "
            "s=w['gene_count']; assert [c.value for c in s[1]]==['GeneID','001','NA']; "
            "assert all(c.data_type=='s' for c in s[1]); "
            "assert {r[0].value for r in list(s.rows)[1:]}=={'001','Gene|literal'}; "
            "assert all(r[0].data_type=='s' for r in list(s.rows)[1:]); "
            "assert all(c.data_type=='n' for r in list(s.rows)[1:] for c in r[1:]); "
            "assert {r[0].value for r in list(w['alignment_report'].rows)[1:]}=={'001','NA'}; w.close()"
        )
        status, output = execute([self.workflow_python, "-c", code, str(quant / "gene_expression_matrices.xlsx")],
                                 work, self.env, work / "check-workbook.log", timeout=30)
        self.assertEqual(status, 0, output)
        return quant

    def assert_mutation_schedules(self, work, source, cores, path, name, rule, replacement):
        original, metadata = path.read_bytes(), path.stat()
        try:
            if replacement is None:
                path.unlink()
            else:
                path.write_bytes(replacement(original))
            status, output = self.run_workflow(work, source, cores, name, dry=True)
            self.assertEqual(status, 0, output[-10000:])
            self.assertRegex(output, r"(?:localrule|rule) " + re.escape(rule) + r":")
        finally:
            path.write_bytes(original)
            os.utime(path, ns=(metadata.st_atime_ns, metadata.st_mtime_ns))

    def test_fastq_to_excel_both_modes_and_default_rerun_tracking(self):
        for stranded, cores, name in ((False, 1, "unstranded work 'quoted'"), (True, 4, "stranded work")):
            with self.subTest(stranded=stranded, cores=cores):
                work, source = self.make_case(name, stranded)
                status, output = self.run_workflow(work, source, cores, "workflow")
                self.assertEqual(status, 0, output[-12000:])
                quant = self.assert_outputs(work, stranded, cores)
                status, output = self.run_workflow(work, source, cores, "unchanged", dry=True)
                self.assertEqual(status, 0, output[-10000:])
                self.assertIn("Nothing to be done", output)
                self.assert_mutation_schedules(work, source, cores, quant / "001" / "t_data.ctab", "missing-ballgown", "stringtie_quant", None)
                self.assert_mutation_schedules(work, source, cores, source / "workflow" / "scripts" / "build_expression_matrices.py",
                                               "changed-matrix-source", "finalize_expression_matrices",
                                               lambda data: data + b"\nAUDIT_SOURCE_CHANGE = True\n")
                self.assert_mutation_schedules(work, source, cores, work / "references" / "mini_index.1.ht2",
                                               "changed-index", "hisat2_align_sort",
                                               lambda data: data[:-1] + bytes([data[-1] ^ 1]))

    def test_ballgown_disabled_with_complete_reference_still_runs(self):
        work, source = self.make_case("ballgown disabled", False)
        config_path = work / "config.json"
        config = json.loads(config_path.read_text())
        config["stringtie"]["use_ballgown"] = False
        config["alignment"].pop("mapq_filter")
        config_path.write_text(json.dumps(config, indent=2))
        status, output = self.run_workflow(work, source, 3, "no-ballgown")
        self.assertEqual(status, 0, output[-12000:])
        self.assert_outputs(work, False, 3, ballgown=False)
        self.assertNotIn("--ballgown", output)
        self.assertEqual(re.findall(r"--mapq (\d+) --hisat2-extra", output), ["0", "0"])

    def test_alignment_mapq_default_fallback_and_override_preserve_strand_qc_cutoff(self):
        for mode, expected in (("shipped default", 0), ("missing key fallback", 0), ("explicit override", 30)):
            with self.subTest(mode=mode):
                work, source = self.make_case("mapq " + mode, True, {"001": 75})
                config_path = work / "config.json"
                config = json.loads(config_path.read_text())
                if mode != "explicit override":
                    config["alignment"].pop("mapq_filter")
                if mode == "missing key fallback":
                    # Snakemake merges the user's config with the shipped YAML;
                    # omit the key from both to exercise Snakefile's fallback.
                    shipped = source / "config/config.yaml"
                    text, count = re.subn(r"(?m)^  mapq_filter:.*\n", "", shipped.read_text())
                    self.assertEqual(count, 1)
                    shipped.write_text(text)
                config_path.write_text(json.dumps(config, indent=2))
                status, output = self.run_workflow(work, source, 3, "mapq-dry-run", dry=True)
                self.assertEqual(status, 0, output[-12000:])
                self.assertEqual(re.findall(r"--mapq (\d+) --hisat2-extra", output), [str(expected)])
                self.assertRegex(output, r"--sample-size 200000 --mapq 30")

    def test_bad_fastq_isolated_in_both_input_modes_and_repair_reuses_good_sample(self):
        fastp = executable("RNASEQ_TEST_FASTP", "fastp")
        if fastp is None:
            self.skipTest("Requires fastp for the raw input isolation case")
        fastp_link = self.root / "bin" / "fastp"
        if not fastp_link.exists():
            fastp_link.symlink_to(fastp)
        for raw in (False, True):
            with self.subTest(raw=raw):
                work, source = self.make_case("isolated " + ("raw" if raw else "clean"), raw)
                config_path = work / "config.json"
                config = json.loads(config_path.read_text())
                if raw:
                    for sample in SAMPLES:
                        destination = work / "rawdata" / sample
                        destination.mkdir(parents=True)
                        for mate in (1, 2):
                            (work / "cleandata" / sample / f"{sample}_{mate}.clean.fastq.gz").rename(
                                destination / f"{sample}_{mate}.fastq.gz")
                    config["fastp"] = {"enabled": True, "trim_poly_g": False}
                    config["threads"]["fastp"] = 2
                    config_path.write_text(json.dumps(config, indent=2))
                bad = work / ("rawdata" if raw else "cleandata") / "001" / (
                    "001_1.fastq.gz" if raw else "001_1.clean.fastq.gz")
                original = bad.read_bytes()
                # Both files exist and are nonempty, but this read has mismatched
                # sequence/quality lengths: a sample job failure, not a DAG error.
                with gzip.open(bad, "wt") as handle:
                    handle.write("@broken/1\nAAAA\n+\nIII\n")
                status, output = self.run_workflow(work, source, 3, "one-bad-sample")
                self.assertNotEqual(status, 0, output[-12000:])
                self.assertIn("profiles/default", output)
                quant = work / "quantification"
                if raw:
                    quant /= "stranded"
                good_paths = [
                    quant / "NA" / "NA.gtf",
                    quant / "NA" / "NA_abundance.tab",
                    work / "alignment" / "NA" / "NA.sorted.bam",
                    work / "alignment" / "read_lengths" / "NA.json",
                ]
                if raw:
                    good_paths += [quant / "strandness" / "NA.env",
                                   work / "qc" / "fastp" / "NA" / "NA.fastp.json"]
                for path in good_paths:
                    self.assertTrue(path.is_file(), str(path) + "\n" + output[-12000:])
                self.assertFalse((quant / "001" / "001.gtf").exists())
                self.assertFalse((quant / "gene_count_matrix.csv").exists())
                self.assertFalse((quant / "gene_expression_matrices.xlsx").exists())
                failure_log = work / "logs" / "001" / ("fastp.log" if raw else "read_lengths.log")
                self.assertIn("invalid/truncated FASTQ record", failure_log.read_text())
                before = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in good_paths}
                bad.write_bytes(original)
                status, output = self.run_workflow(work, source, 3, "repaired-sample")
                self.assertEqual(status, 0, output[-12000:])
                self.assert_outputs(work, raw, 3)
                for path, expected in before.items():
                    self.assertEqual((path.read_bytes(), path.stat().st_mtime_ns), expected, str(path))

    def test_missing_index_is_rejected_before_alignment(self):
        work, source = self.make_case("missing-index", False, {"001": 75})
        (work / "references" / "mini_index.8.ht2").unlink()
        status, output = self.run_workflow(work, source, 1, "missing-index", dry=True)
        self.assertNotEqual(status, 0)
        self.assertIn("complete, nonempty set of 8", output)
        self.assertFalse((work / "alignment").exists())

    def test_unstranded_wrong_contig_blocks_quantification(self):
        work, source = self.make_case("wrong-contig", False, {"001": 75})
        annotation = work / "references" / "reference.gtf"
        annotation.write_text(annotation.read_text().replace("chrMini", "chrWrong"))
        status, output = self.run_workflow(work, source, 4, "wrong-contig")
        self.assertNotEqual(status, 0)
        log = work / "logs" / "001" / "validate_bam_reference.log"
        self.assertIn("absent from BAM SQ", log.read_text())
        self.assertNotIn("localrule stringtie_quant:", output)
        self.assertFalse((work / "quantification" / "001" / "001.gtf").exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
