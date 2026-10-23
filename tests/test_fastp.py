#!/usr/bin/env python3
"""Real fastp regression tests on synthetic paired FASTQ; never install tools.

RNASEQ_TEST_FASTP selects an existing fastp executable. The raw FASTQ-to-Excel
case also requires RNASEQ_TEST_FULL_FASTQ=1 and the tool variables documented
in test_workflow_regressions.py. RNASEQ_TEST_KEEP=1 preserves artifacts.
"""
import gzip
import hashlib
import json
import os
from pathlib import Path
import random
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "workflow/scripts/run_fastp.py"


def checksum(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_fastq(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    if len(lines) % 4:
        raise AssertionError("Truncated output FASTQ: " + str(path))
    return {lines[i].split()[0][1:]: (lines[i + 1], lines[i + 3])
            for i in range(0, len(lines), 4)}


class FastpTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        candidate = os.environ.get("RNASEQ_TEST_FASTP") or shutil.which("fastp")
        if not candidate:
            raise unittest.SkipTest("Requires an existing fastp executable")
        cls.fastp = str(Path(candidate).resolve())
        cls.root = Path(tempfile.mkdtemp(prefix="rnaseq-fastp-tests-"))
        print("fastp regression artifacts: " + str(cls.root), flush=True)

    @classmethod
    def tearDownClass(cls):
        if os.environ.get("RNASEQ_TEST_KEEP") != "1":
            shutil.rmtree(cls.root)

    def setUp(self):
        self.work = self.root / (self._testMethodName + " work 'quoted'")
        self.work.mkdir()
        bindir = self.work / "bin"
        bindir.mkdir()
        (bindir / "fastp").symlink_to(self.fastp)
        self.env = dict(os.environ)
        self.env["PATH"] = str(bindir) + os.pathsep + self.env.get("PATH", "")
        self.env["PYTHONDONTWRITEBYTECODE"] = "1"

    def write_pairs(self, cases, directory="raw", bad_ids=False):
        """cases contains (name, R1 sequence, R2 sequence, R1 quality, R2 quality)."""
        folder = self.work / directory
        folder.mkdir()
        paths = [folder / "sample_1.fastq.gz", folder / "sample_2.fastq.gz"]
        for mate, path in enumerate(paths, 1):
            with gzip.open(path, "wt") as handle:
                for index, (name, first, second, quality1, quality2) in enumerate(cases, 1):
                    # A NextSeq instrument ID deliberately triggers fastp's
                    # instrument-dependent auto polyG mode unless disabled.
                    key = f"NS500123:1:HG:1:1:{index}:{name}"
                    if bad_ids and mate == 2:
                        key += "mismatch"
                    sequence, quality = (first, quality1) if mate == 1 else (second, quality2)
                    handle.write(f"@{key} {mate}:N:0:ACGT\n{sequence}\n+\n{quality}\n")
        return paths

    @staticmethod
    def case(name, length=80, quality=None, n=0, poly_g=0):
        rng = random.Random("synthetic-pair-" + name)
        sequences = ["".join(rng.choices("ACGT", k=length)) for _ in range(2)]
        if n:
            sequences[0] = "N" * n + sequences[0][n:]
        if poly_g:
            sequences = [seq[:length - poly_g - 1] + "A" + "G" * poly_g for seq in sequences]
        quality = quality if quality is not None else "I" * length
        return name, sequences[0], sequences[1], quality, "I" * length

    def run_fastp(self, paths, trim=False, name="clean", extra=()):
        folder = self.work / name
        folder.mkdir(exist_ok=True)
        outputs = {"r1": folder / "sample_1.clean.fastq.gz", "r2": folder / "sample_2.clean.fastq.gz",
                   "html": folder / "sample.fastp.html", "json": folder / "sample.fastp.json"}
        command = [sys.executable, str(SCRIPT), "--sample", "sample", "--r1", str(paths[0]),
                   "--r2", str(paths[1]), "--out-r1", str(outputs["r1"]), "--out-r2", str(outputs["r2"]),
                   "--html", str(outputs["html"]), "--json", str(outputs["json"]),
                   "--trim-poly-g", str(trim).lower(), *extra]
        result = subprocess.run(command, cwd=self.work, env=self.env, capture_output=True, text=True, timeout=90)
        (folder / "command.json").write_text(json.dumps(command, indent=2) + "\n")
        (folder / "command.log").write_text(result.stdout + result.stderr)
        return result, outputs

    def assert_no_outputs(self, outputs):
        self.assertTrue(all(not path.exists() for path in outputs.values()), outputs)

    def test_default_filters_keep_exact_boundaries_and_reports(self):
        cases = [self.case("good"), self.case("q20", quality="5" * 80),
                 self.case("q19", quality="4" * 80),
                 self.case("low30percent", quality="!" * 24 + "I" * 56),
                 self.case("low_over30percent", quality="!" * 25 + "I" * 55),
                 self.case("N5", n=5), self.case("N6", n=6),
                 self.case("length40", length=40), self.case("length39", length=39)]
        paths = self.write_pairs(cases)
        before = [checksum(path) for path in paths]
        result, outputs = self.run_fastp(paths)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([checksum(path) for path in paths], before)
        expected = {"good", "q20", "low30percent", "N5", "length40"}
        first, second = read_fastq(outputs["r1"]), read_fastq(outputs["r2"])
        self.assertEqual(set(first), set(second))
        self.assertEqual({key.rsplit(":", 1)[1] for key in first}, expected)
        # The requested -q/-u flags filter reads; they do not imply cutting
        # every low-Q base from otherwise acceptable reads.
        kept_low = next(value for key, value in first.items() if key.endswith(":low30percent"))
        self.assertEqual(kept_low[1], "!" * 24 + "I" * 56)
        report = json.loads(outputs["json"].read_text())
        self.assertEqual(report["summary"]["before_filtering"]["total_reads"], len(cases) * 2)
        self.assertEqual(report["summary"]["after_filtering"]["total_reads"], len(expected) * 2)
        self.assertIn("<html", outputs["html"].read_text().lower())
        arguments = shlex.split(report["command"])
        self.assertIn("--detect_adapter_for_pe", arguments)
        self.assertIn("--disable_trim_poly_g", arguments)
        self.assertNotIn("--trim_poly_g", arguments)
        for alternatives, expected_value in ((('-w', '--thread'), '8'), (('-q', '--qualified_quality_phred'), '20'),
                                             (('-u', '--unqualified_percent_limit'), '30'),
                                             (('-n', '--n_base_limit'), '5'), (('-l', '--length_required'), '40')):
            option = next(flag for flag in alternatives if flag in arguments)
            self.assertEqual(arguments[arguments.index(option) + 1], expected_value)

    def test_polyg_enabled_trims_and_disabled_preserves_nextseq_reads(self):
        paths = self.write_pairs([self.case("polyG", poly_g=20)])
        before = [checksum(path) for path in paths]
        for enabled, length in ((False, 80), (True, 60)):
            with self.subTest(trim_poly_g=enabled):
                result, outputs = self.run_fastp(paths, enabled, name="clean-" + str(enabled))
                self.assertEqual(result.returncode, 0, result.stderr)
                for mate in ("r1", "r2"):
                    values = list(read_fastq(outputs[mate]).values())
                    self.assertEqual(len(values), 1)
                    self.assertEqual(len(values[0][0]), length)
                    self.assertEqual(values[0][0].endswith("G" * 20), not enabled)
                arguments = shlex.split(json.loads(outputs["json"].read_text())["command"])
                self.assertIn("--trim_poly_g" if enabled else "--disable_trim_poly_g", arguments)
        self.assertEqual([checksum(path) for path in paths], before)

    def test_mismatched_ids_rejected_even_if_low_quality_would_hide_them(self):
        for quality in ("I" * 80, "!" * 80):
            with self.subTest(quality=quality[0]):
                name = "high" if quality[0] == "I" else "low"
                paths = self.write_pairs([self.case("bad_pair", quality=quality)], directory=name, bad_ids=True)
                before = [checksum(path) for path in paths]
                result, outputs = self.run_fastp(paths, name="clean-" + name)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("read IDs disagree", result.stderr)
                self.assert_no_outputs(outputs)
                self.assertEqual([checksum(path) for path in paths], before)

    def test_mismatched_record_counts_rejected(self):
        paths = self.write_pairs([self.case("one"), self.case("two")])
        with gzip.open(paths[1], "wt") as handle:
            case = self.case("one")
            handle.write("@NS500123:1:HG:1:1:1:one 2:N:0:ACGT\n" + case[2] + "\n+\n" + case[4] + "\n")
        result, outputs = self.run_fastp(paths)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("different record counts", result.stderr)
        self.assert_no_outputs(outputs)

    def test_truncated_input_rejected_before_fastp(self):
        paths = self.write_pairs([self.case("one")])
        with gzip.open(paths[0], "at") as handle:
            handle.write("@broken/1\nACTG\n+\n")
        result, outputs = self.run_fastp(paths)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("invalid/truncated FASTQ", result.stderr)
        self.assert_no_outputs(outputs)

    def test_input_output_alias_rejected_and_raw_preserved(self):
        paths = self.write_pairs([self.case("one")])
        before = [checksum(path) for path in paths]
        result, outputs = self.run_fastp(paths, extra=("--out-r1", str(paths[0])))
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual([checksum(path) for path in paths], before)
        self.assert_no_outputs(outputs)

    def test_symlink_and_hardlink_to_raw_are_rejected(self):
        paths = self.write_pairs([self.case("one")])
        before = [checksum(path) for path in paths]
        for kind in ("symlink", "hardlink"):
            with self.subTest(kind=kind):
                alias = self.work / (kind + ".fastq.gz")
                if kind == "symlink":
                    alias.symlink_to(paths[0])
                else:
                    os.link(paths[0], alias)
                result, outputs = self.run_fastp(paths, extra=("--out-r1", str(alias)), name=kind)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("aliases raw input", result.stderr)
                self.assertEqual([checksum(path) for path in paths], before)
                self.assertEqual(checksum(alias), before[0])
                self.assert_no_outputs(outputs)

    def test_input_mates_and_output_files_must_be_distinct(self):
        paths = self.write_pairs([self.case("one")])
        result, outputs = self.run_fastp([paths[0], paths[0]], name="same-mates")
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("different files", result.stderr)
        self.assert_no_outputs(outputs)
        result, outputs = self.run_fastp(paths, name="same-outputs", extra=(
            "--html", str(self.work / "same-outputs" / "sample_1.clean.fastq.gz")))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("distinct", result.stderr)
        self.assert_no_outputs(outputs)

    def test_fastp_failure_does_not_publish_partial_outputs(self):
        paths = self.write_pairs([self.case("one")])
        before = [checksum(path) for path in paths]
        fastp = self.work / "bin/fastp"
        fastp.unlink()
        fastp.write_text("#!" + sys.executable + "\n"
                         "import sys\nfrom pathlib import Path\n"
                         "if '--version' in sys.argv:\n    print('fastp failure fixture')\n    sys.exit(0)\n"
                         "for flag in ('-o', '-O', '-h', '-j'):\n"
                         "    Path(sys.argv[sys.argv.index(flag) + 1]).write_text('partial output')\n"
                         "sys.exit(7)\n")
        fastp.chmod(0o755)
        result, outputs = self.run_fastp(paths)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("7", result.stderr)
        self.assert_no_outputs(outputs)
        self.assertEqual(list(self.work.glob("**/.fastp-*")), [])
        self.assertEqual([checksum(path) for path in paths], before)

    def test_invalid_filter_values_are_rejected(self):
        paths = self.write_pairs([self.case("one")])
        for option, value in (("--threads", "0"), ("--qualified-quality-phred", "94"),
                              ("--unqualified-percent-limit", "101"), ("--n-base-limit", "-1"),
                              ("--n-base-limit", "51"),
                              ("--length-required", "0"), ("--trim-poly-g", "maybe")):
            with self.subTest(option=option):
                result, outputs = self.run_fastp(paths, extra=(option, value), name=option[2:])
                self.assertNotEqual(result.returncode, 0)
                self.assert_no_outputs(outputs)

    def test_all_filtered_reads_keep_diagnostic_reports(self):
        paths = self.write_pairs([self.case("too_short", length=39)])
        result, outputs = self.run_fastp(paths)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(read_fastq(outputs["r1"]), {})
        self.assertEqual(read_fastq(outputs["r2"]), {})
        report = json.loads(outputs["json"].read_text())
        self.assertEqual(report["summary"]["after_filtering"]["total_reads"], 0)
        self.assertIn("<html", outputs["html"].read_text().lower())
        measured = subprocess.run([sys.executable, str(ROOT / "workflow/scripts/read_lengths.py"), "measure",
                                   "--sample", "sample", "--r1", str(outputs["r1"]), "--r2", str(outputs["r2"]),
                                   "--output", str(self.work / "read-lengths.json")],
                                  capture_output=True, text=True, timeout=20)
        self.assertNotEqual(measured.returncode, 0)
        self.assertIn("contains no reads", measured.stderr)
        self.assertFalse((self.work / "read-lengths.json").exists())


@unittest.skipUnless(os.environ.get("RNASEQ_TEST_FULL_FASTQ") == "1", "Set RNASEQ_TEST_FULL_FASTQ=1 for raw-to-Excel regression")
class FastpWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import test_workflow_regressions as regression
        candidate = os.environ.get("RNASEQ_TEST_FASTP") or shutil.which("fastp")
        if not candidate:
            raise unittest.SkipTest("Requires an existing fastp executable")
        cls.regression = regression
        regression.FullFastqRegressionTests.setUpClass()
        cls.fixture = regression.FullFastqRegressionTests()
        (cls.fixture.root / "bin/fastp").symlink_to(Path(candidate).resolve())

    @classmethod
    def tearDownClass(cls):
        cls.regression.FullFastqRegressionTests.tearDownClass()

    def test_raw_path_aliases_are_rejected_before_jobs_start(self):
        for kind in ("clean-symlink", "report-hardlink", "log-hardlink"):
            with self.subTest(kind=kind):
                work, source = self.fixture.make_case("raw-protect-" + kind, False, {"001": 75})
                raw_directory = work / "rawdata/001"
                raw_directory.mkdir(parents=True)
                for mate in (1, 2):
                    (work / "cleandata/001" / f"001_{mate}.clean.fastq.gz").rename(
                        raw_directory / f"001_{mate}.fastq.gz")
                shutil.rmtree(work / "cleandata")
                raw = raw_directory / "001_1.fastq.gz"
                before = checksum(raw)
                config_path = work / "config.json"
                config = json.loads(config_path.read_text())
                config["fastp"] = {"enabled": True, "trim_poly_g": False}
                config["paths"].update(raw_dir="rawdata", qc_dir="qc/fastp")
                config_path.write_text(json.dumps(config, indent=2))
                if kind == "clean-symlink":
                    (work / "cleandata").mkdir()
                    (work / "cleandata/001").symlink_to(raw_directory, target_is_directory=True)
                    expected = "output resolves inside a raw input directory"
                else:
                    output = work / ("qc/fastp/001/001.fastp.html" if kind == "report-hardlink" else "logs/001/fastp.log")
                    output.parent.mkdir(parents=True)
                    os.link(raw, output)
                    expected = "linked to the same file"
                status, output = self.fixture.run_workflow(work, source, 3, "raw-protect", dry=True)
                self.assertNotEqual(status, 0)
                self.assertIn(expected, output)
                self.assertNotIn("localrule fastp_clean:", output)
                self.assertEqual(checksum(raw), before)

    def test_raw_fastq_to_excel_and_fastp_dependency_tracking(self):
        work, source = self.fixture.make_case("raw work 'quoted'", False)
        raw_inputs = []
        for sample in self.regression.SAMPLES:
            directory = work / "rawdata" / sample
            directory.mkdir(parents=True)
            for mate in (1, 2):
                raw = directory / f"{sample}_{mate}.fastq.gz"
                (work / "cleandata" / sample / f"{sample}_{mate}.clean.fastq.gz").rename(raw)
                raw_inputs.append(raw)
        shutil.rmtree(work / "cleandata")
        before = {str(path): checksum(path) for path in raw_inputs}
        config_path = work / "config.json"
        config = json.loads(config_path.read_text())
        config["fastp"] = {"enabled": True, "trim_poly_g": False}
        config["paths"].update(raw_dir="rawdata", qc_dir="qc/fastp")
        config["raw_sample_scan"] = {"r1_suffix": "_1.fastq.gz", "r2_suffix": "_2.fastq.gz", "strict": True}
        config["threads"]["fastp"] = 8
        config_path.write_text(json.dumps(config, indent=2))
        status, output = self.fixture.run_workflow(work, source, 3, "raw-dry-run", dry=True)
        self.assertEqual(status, 0, output[-12000:])
        self.assertIn("fastp_clean", output)
        self.assertFalse((work / "cleandata").exists())
        status, output = self.fixture.run_workflow(work, source, 3, "raw-workflow")
        self.assertEqual(status, 0, output[-12000:])
        self.fixture.assert_outputs(work, False, 3)
        baseline_counts = (work / "quantification/gene_count_matrix.csv").read_bytes()
        for sample in self.regression.SAMPLES:
            report_dir = work / "qc/fastp" / sample
            report = json.loads((report_dir / f"{sample}.fastp.json").read_text())
            self.assertEqual(report["summary"]["before_filtering"]["total_reads"], 400)
            self.assertEqual(report["summary"]["after_filtering"]["total_reads"], 400)
            self.assertTrue((report_dir / f"{sample}.fastp.html").is_file())
            log = (work / "logs" / sample / "fastp.log").read_text()
            command = next(line for line in log.splitlines() if line.startswith("[INFO] Command: "))
            arguments = shlex.split(command.split(": ", 1)[1])
            self.assertEqual(arguments[arguments.index("-w") + 1], "3")
            self.assertIn("--disable_trim_poly_g", arguments)
        self.assertEqual({str(path): checksum(path) for path in raw_inputs}, before)
        status, output = self.fixture.run_workflow(work, source, 3, "raw-unchanged", dry=True)
        self.assertEqual(status, 0, output[-12000:])
        self.assertIn("Nothing to be done", output)
        self.fixture.assert_mutation_schedules(
            work, source, 3, source / "workflow/scripts/run_fastp.py", "changed-fastp-script", "fastp_clean",
            lambda data: data + b"\nAUDIT_FASTP_SOURCE_CHANGE = True\n")
        self.fixture.assert_mutation_schedules(
            work, source, 3, work / "qc/fastp/001/001.fastp.html", "missing-fastp-html", "fastp_clean", None)
        config["fastp"]["trim_poly_g"] = True
        config_path.write_text(json.dumps(config, indent=2))
        status, output = self.fixture.run_workflow(work, source, 3, "changed-polyg", dry=True)
        self.assertEqual(status, 0, output[-12000:])
        self.assertIn("fastp_clean", output)
        self.assertIn("--trim-poly-g true", output)
        # Now execute the changed mode with synthetic polyG-contaminated raw
        # reads. Append tails only after six terminal non-G bases: fastp's
        # polyG detector tolerates mismatches, so a single preceding non-G
        # base would not guarantee preservation of the original sequence.
        # Raw data mix 75/95 bp and 150/170 bp reads; the expected clean
        # lengths remain exactly 75/150 bp and still match the reference.
        expected_clean = {}
        raw_bases = {sample: 0 for sample in self.regression.SAMPLES}
        tail_reads = {sample: 0 for sample in self.regression.SAMPLES}
        for raw in raw_inputs:
            original = read_fastq(raw)
            expected_clean[raw] = original
            with gzip.open(raw, "wt") as handle:
                for key, (sequence, quality) in original.items():
                    tail = 20 if "G" not in sequence[-6:] else 0
                    tail_reads[raw.parent.name] += bool(tail)
                    raw_bases[raw.parent.name] += len(sequence) + tail
                    handle.write(f"@{key}\n{sequence + 'G' * tail}\n+\n{quality + 'I' * tail}\n")
        contaminated_hashes = {str(path): checksum(path) for path in raw_inputs}
        status, output = self.fixture.run_workflow(work, source, 3, "polyg-actual-rerun")
        self.assertEqual(status, 0, output[-12000:])
        self.fixture.assert_outputs(work, False, 3)
        self.assertEqual((work / "quantification/gene_count_matrix.csv").read_bytes(), baseline_counts)
        for raw, expected in expected_clean.items():
            sample = raw.parent.name
            mate = 1 if raw.name.endswith("_1.fastq.gz") else 2
            clean = work / "cleandata" / sample / f"{sample}_{mate}.clean.fastq.gz"
            self.assertEqual(read_fastq(clean), expected)
        for sample, length in self.regression.SAMPLES.items():
            report = json.loads((work / "qc/fastp" / sample / f"{sample}.fastp.json").read_text())
            self.assertGreater(tail_reads[sample], 0)
            self.assertEqual(report["summary"]["before_filtering"]["total_bases"], raw_bases[sample])
            self.assertEqual(report["summary"]["after_filtering"]["total_bases"], length * 400)
            self.assertIn("--trim_poly_g", shlex.split(report["command"]))
        self.assertEqual({str(path): checksum(path) for path in raw_inputs}, contaminated_hashes)
        # The raw-data entrance also schedules the existing strand-specific
        # branch; its real RSeQC execution is covered by the original suite.
        config["stranded"] = True
        config["fastp"]["trim_poly_g"] = False
        config_path.write_text(json.dumps(config, indent=2))
        status, output = self.fixture.run_workflow(work, source, 3, "raw-stranded-dry-run", dry=True)
        self.assertEqual(status, 0, output[-12000:])
        self.assertIn("fastp_clean", output)
        self.assertIn("infer_strandness", output)


if __name__ == "__main__":
    unittest.main(verbosity=2)
