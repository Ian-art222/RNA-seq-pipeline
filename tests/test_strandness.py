"""Synthetic strand evidence tests; no biological data or external tools needed."""

import argparse
import contextlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts"))
import strandness


def output(mode1="0.04", mode2="0.91", failed="0.05", protocol="PairEnd"):
    return "\n".join([
        "This is {} Data".format(protocol),
        "{}: {}".format(strandness.FRACTION_LABELS["failed_fraction"], failed),
        "{}: {}".format(strandness.FRACTION_LABELS["mode1"], mode1),
        "{}: {}".format(strandness.FRACTION_LABELS["mode2"], mode2),
        "",
    ])


def parsed(sample="sample1", **kwargs):
    result = strandness.parse_rseqc(output(**kwargs), "Total 200000 usable reads were sampled\n")
    result.update(sample=sample, schema_version=1, tool={"returncode": 0})
    return result


class ParseTests(unittest.TestCase):
    def test_both_directions(self):
        for first, second, flag in [("0.04", "0.91", "--rf"), ("0.91", "0.04", "--fr")]:
            with self.subTest(flag=flag):
                report = parsed(mode1=first, mode2=second)
                self.assertEqual(report["status"], "PASS")
                self.assertEqual(report["stringtie_flag"], flag)
                self.assertEqual(report["usable_reads"], 200000)
                self.assertAlmostEqual(report["dominant_fraction"], 0.91 / 0.95)

    def test_thresholds_are_inclusive(self):
        self.assertEqual(parsed(mode1="0.1", mode2="0.4", failed="0.5")["status"], "PASS")
        for first, second, failed in [("0.1", "0.3999", "0.5001"), ("0.2001", "0.7999", "0")]:
            with self.subTest(values=(first, second, failed)):
                self.assertEqual(parsed(mode1=first, mode2=second, failed=failed)["status"], "FAIL")

    def test_missing_or_duplicate_each_fraction(self):
        for key, label in strandness.FRACTION_LABELS.items():
            with self.subTest(field=key):
                line = next(line for line in output().splitlines(True) if label in line)
                for text in [output().replace(line, ""), output() + line]:
                    report = strandness.parse_rseqc(text, "Total 200000 usable reads were sampled")
                    self.assertEqual(report["status"], "FAIL")
                    self.assertIsNone(report["stringtie_flag"])
                    self.assertIn("Expected exactly one", "; ".join(report["errors"]))

    def test_invalid_numbers(self):
        for value in ["nan", "NaN", "inf", "-inf", "1.01", "-0.1", "wrong", ""]:
            for key in ["mode1", "mode2", "failed"]:
                with self.subTest(value=value, field=key):
                    report = parsed(**{key: value})
                    self.assertEqual(report["status"], "FAIL")
                    json.dumps(report, allow_nan=False)

    def test_protocol_single_mixed_unknown_and_duplicate(self):
        for text in [output(protocol="SingleEnd"), "Unknown Data type\n", output() + "This is SingleEnd Data\n", output() * 2]:
            report = strandness.parse_rseqc(text, "Total 200000 usable reads were sampled")
            self.assertEqual(report["status"], "FAIL")
            self.assertIsNone(report["stringtie_flag"])

    def test_invalid_or_missing_sampling_counter(self):
        for count in ["0", "200001", "-1", "1.2", "unknown"]:
            with self.subTest(count=count):
                result = strandness.parse_rseqc(output(), "Total {} usable reads were sampled".format(count))
                self.assertEqual(result["status"], "FAIL")
        for stderr in ["", "Total 200000 usable reads were sampled\n" * 2]:
            self.assertEqual(strandness.parse_rseqc(output(), stderr)["status"], "FAIL")

    def test_short_positive_sample_is_recorded_and_warned(self):
        for actual in (1, 999):
            result = strandness.parse_rseqc(output(), "Total {} usable reads were sampled".format(actual))
            self.assertEqual(result["status"], "PASS")
            self.assertEqual(result["usable_reads"], actual)
            self.assertTrue(result["warnings"])

    def test_sampling_counter_may_be_on_stdout(self):
        report = strandness.parse_rseqc(output() + "Total 200000 usable reads were sampled", "")
        self.assertEqual(report["status"], "PASS")

    def test_sampling_counter_after_progress_message(self):
        # RSeQC does not insert a newline after this progress message when the
        # requested count is reached before EOF.
        stderr = "Reading reference gene model a.bed ... Done\nLoading SAM/BAM file ... Total 200000 usable reads were sampled\n"
        self.assertEqual(strandness.parse_rseqc(output(), stderr)["status"], "PASS")
        self.assertEqual(strandness.parse_rseqc(output(), stderr + stderr)["status"], "FAIL")

    def test_zero_and_inconsistent_fraction_sums(self):
        self.assertEqual(parsed(mode1="0", mode2="0", failed="1")["status"], "FAIL")
        self.assertEqual(parsed(mode1="0.1", mode2="0.8", failed="0.5")["status"], "FAIL")
        self.assertEqual(parsed(mode1="0.04", mode2="0.91", failed="0.0501")["status"], "PASS")


class FileTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name)

    def write_report(self, name, report):
        path = self.path / (name + ".json")
        path.write_text(json.dumps(report))
        return str(path)

    def gate(self, reports, samples=None):
        return strandness.gate(argparse.Namespace(
            samples=samples or ["sample1", "sample2"], reports=reports,
            diagnostic=str(self.path / "diagnostic.tsv"), output=str(self.path / "strandness.env"),
        ))

    def test_gate_all_pass(self):
        reports = [self.write_report(sample, parsed(sample)) for sample in ("sample1", "sample2")]
        self.assertEqual(self.gate(reports), 0)
        self.assertEqual((self.path / "strandness.env").read_text(), "STRINGTIE_STRAND=--rf\n")
        diagnostic = (self.path / "diagnostic.tsv").read_text()
        self.assertIn("200000", diagnostic)
        self.assertIn("0.91", diagnostic)

    def test_gate_mixed_direction_removes_stale_env(self):
        (self.path / "strandness.env").write_text("STRINGTIE_STRAND=--rf\n")
        reports = [self.write_report("one", parsed("sample1")), self.write_report("two", parsed("sample2", mode1="0.91", mode2="0.04"))]
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(self.gate(reports), 1)
        self.assertFalse((self.path / "strandness.env").exists())
        self.assertIn("disagree", (self.path / "diagnostic.tsv").read_text())

    def test_gate_missing_wrong_duplicate_and_failed_samples(self):
        good = self.write_report("one", parsed("sample1"))
        failed = self.write_report("bad", parsed("sample2", mode1="0.5", mode2="0.45"))
        for reports, samples in [([good], None), ([good, good], None), ([good, failed], None), ([good, good], ["sample1", "sample1"]), ([good, "missing.json"], None)]:
            with self.subTest(reports=reports, samples=samples), contextlib.redirect_stderr(io.StringIO()):
                (self.path / "strandness.env").write_text("stale")
                self.assertEqual(self.gate(reports, samples), 1)
                self.assertFalse((self.path / "strandness.env").exists())
                self.assertTrue((self.path / "diagnostic.tsv").is_file())

    def test_gate_revalidates_pass_evidence(self):
        good = self.write_report("one", parsed("sample1"))
        for mutation in [{"usable_reads": 0}, {"stringtie_flag": "--fr"}, {"mode2": 0.5}, {"protocol": "SingleEnd"}, {"tool": {"returncode": 1}}]:
            with self.subTest(mutation=mutation), contextlib.redirect_stderr(io.StringIO()):
                report = parsed("sample2")
                report.update(mutation)
                bad = self.write_report("two", report)
                self.assertEqual(self.gate([good, bad]), 1)

    def test_option_allows_only_the_exact_assignment(self):
        env = self.path / "strandness.env"
        for value in ["STRINGTIE_STRAND=--fr\n", "STRINGTIE_STRAND=--rf"]:
            env.write_text(value)
            with contextlib.redirect_stdout(io.StringIO()) as captured:
                self.assertEqual(strandness.option(argparse.Namespace(env=env)), 0)
                self.assertIn(captured.getvalue().strip(), strandness.FLAGS)
        for value in ["", "STRINGTIE_STRAND=NA\n", "STRINGTIE_STRAND=--rf\nCOMMAND=bad\n", "STRINGTIE_STRAND=--rf; bad\n"]:
            env.write_text(value)
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(strandness.option(argparse.Namespace(env=env)), 1)

    def infer_args(self):
        return argparse.Namespace(sample="sample1", bam="missing.bam", bai="missing.bam.bai", bed12="missing.bed", output=str(self.path / "report.json"), log_dir=str(self.path / "logs"), sample_size=200000, mapq=30, min_determinable=0.50, min_dominance=0.80)

    def test_input_failure_writes_report_and_logs(self):
        args = self.infer_args()
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(strandness.infer(args), 0)
        report = json.loads(Path(args.output).read_text())
        self.assertEqual(report["status"], "FAIL")
        self.assertIn("BAM", report["errors"][0])
        self.assertTrue((self.path / "logs" / "command.json").is_file())

    def test_invalid_threshold_still_writes_diagnostic_json(self):
        args = self.infer_args()
        args.min_dominance = float("nan")
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(strandness.infer(args), 0)
        report = json.loads(Path(args.output).read_text())
        self.assertEqual(report["status"], "FAIL")
        self.assertEqual(report["thresholds"]["min_dominance"], "nan")

    def test_runtime_failure_retains_actual_reads_and_raw_fractions(self):
        args = self.infer_args()
        completion = subprocess.CompletedProcess(["infer_experiment.py"], 2, output(), "Total 50 usable reads were sampled\n")
        version = subprocess.CompletedProcess(["infer_experiment.py", "--version"], 0, "infer_experiment.py 5.0.4\n", "")
        with mock.patch.object(strandness, "validate_inputs", return_value={}), mock.patch.object(strandness.shutil, "which", return_value="/tool/infer_experiment.py"), mock.patch.object(strandness.subprocess, "run", side_effect=[version, completion]), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(strandness.infer(args), 0)
        report = json.loads(Path(args.output).read_text())
        self.assertEqual(report["status"], "FAIL")
        self.assertEqual(report["usable_reads"], 50)
        self.assertEqual(report["mode2"], 0.91)
        self.assertEqual(report["tool"]["returncode"], 2)
        self.assertEqual(report["tool"]["version"], "infer_experiment.py 5.0.4")
        self.assertIsNone(report["stringtie_flag"])
        self.assertEqual((self.path / "logs" / "stdout.txt").read_text(), output())

    def test_inference_output_cannot_alias_input(self):
        args = self.infer_args()
        bam = self.path / "input.bam"
        bam.write_bytes(b"protected input")
        args.bam = str(bam)
        for path in [bam, self.path / "hardlink.bam", self.path / "symlink.bam"]:
            if path.name == "hardlink.bam":
                os.link(bam, path)
            elif path.name == "symlink.bam":
                path.symlink_to(bam)
            args.output = str(path)
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(strandness.infer(args), 1)
            self.assertEqual(bam.read_bytes(), b"protected input")

    def test_inference_log_cannot_alias_input(self):
        args = self.infer_args()
        raw = self.path / "logs" / "stdout.txt"
        raw.parent.mkdir()
        raw.write_bytes(b"protected input")
        args.bam = str(raw)
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(strandness.infer(args), 1)
        self.assertEqual(raw.read_bytes(), b"protected input")

    def test_gate_cannot_delete_a_report_aliased_by_env(self):
        report = self.write_report("one", parsed("sample1"))
        args = argparse.Namespace(samples=["sample1"], reports=[report], diagnostic=str(self.path / "diagnostic.tsv"), output=report)
        before = Path(report).read_bytes()
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(strandness.gate(args), 1)
        self.assertEqual(Path(report).read_bytes(), before)

    def test_actual_rseqc_full_requested_sample(self):
        executable = shutil.which("infer_experiment.py")
        if not executable:
            self.skipTest("RSeQC is only required in the RNA-seq runtime environment")
        try:
            import pysam
        except ImportError:
            self.skipTest("pysam is only required in the RNA-seq runtime environment")
        args = self.infer_args()
        bam, bed = self.path / "fixture.bam", self.path / "fixture.bed"
        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
        with pysam.AlignmentFile(str(bam), "wb", header=header) as handle:
            # Synthetic R1 reverse / R2 forward on a plus transcript -> --rf.
            for index in range(20):
                read = pysam.AlignedSegment()
                read.query_name = "synthetic{}".format(index // 2)
                read.query_sequence = "A" * 20
                read.flag = 81 if index % 2 == 0 else 161
                read.reference_id = 0
                read.reference_start = 10 + index * 2
                read.mapping_quality = 60
                read.cigarstring = "20M"
                read.next_reference_id = 0
                read.next_reference_start = read.reference_start + (2 if index % 2 == 0 else -2)
                handle.write(read)
        pysam.index(str(bam))
        bed.write_text("chr1\t0\t1000\tt1\t0\t+\t0\t1000\t0\t1\t1000,\t0,\n")
        args.bam, args.bai, args.bed12 = str(bam), str(bam) + ".bai", str(bed)
        args.sample_size = 10
        self.assertEqual(strandness.infer(args), 0)
        report = json.loads(Path(args.output).read_text())
        self.assertEqual(report["status"], "PASS", report["errors"])
        self.assertEqual(report["stringtie_flag"], "--rf")
        self.assertEqual(report["usable_reads"], 10)
        self.assertIn("Loading SAM/BAM file ...  Total 10", Path(report["logs"]["stderr"]).read_text())

    def test_bam_annotation_validation(self):
        try:
            import pysam
        except ImportError:
            self.skipTest("pysam is only required in the RNA-seq runtime environment")
        bam = self.path / "input.bam"
        bed = self.path / "annotation.bed"
        bed.write_text("chr1\t10\t30\tt1\t0\t+\t10\t30\t0\t1\t20,\t0,\n")
        header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 100}]}
        with pysam.AlignmentFile(str(bam), "wb", header=header):
            pass
        pysam.index(str(bam))
        self.assertEqual(strandness.validate_inputs(bam, str(bam) + ".bai", bed)["bed12_transcripts"], 1)
        bed.write_text("chr2\t10\t30\tt1\t0\t+\t10\t30\t0\t1\t20,\t0,\n")
        with self.assertRaisesRegex(ValueError, "absent from BAM"):
            strandness.validate_inputs(bam, str(bam) + ".bai", bed)
        bed.write_text("chr1\t90\t110\tt1\t0\t+\t90\t110\t0\t1\t20,\t0,\n")
        with self.assertRaisesRegex(ValueError, "exceeds BAM"):
            strandness.validate_inputs(bam, str(bam) + ".bai", bed)
        header["PG"] = [{"ID": "hisat2", "PN": "hisat2", "CL": "hisat2 --rna-strandness RF -x ref"}]
        with pysam.AlignmentFile(str(bam), "wb", header=header):
            pass
        pysam.index(str(bam))
        with self.assertRaisesRegex(ValueError, "--rna-strandness"):
            strandness.validate_inputs(bam, str(bam) + ".bai", bed)


if __name__ == "__main__":
    unittest.main()
