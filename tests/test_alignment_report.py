"""Regression tests for standard HISAT2 paired-end summary validation."""

import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import build_alignment_report


SUMMARY = """100 reads; of these:
  100 (100.00%) were paired; of these:
    20 (20.00%) aligned concordantly 0 times
    70 (70.00%) aligned concordantly exactly 1 time
    10 (10.00%) aligned concordantly >1 times
    ----
    20 pairs aligned concordantly 0 times; of these:
      5 (25.00%) aligned discordantly 1 time
    ----
    15 pairs aligned 0 times concordantly or discordantly; of these:
      30 mates make up the pairs; of these:
        10 (33.33%) aligned 0 times
        18 (60.00%) aligned exactly 1 time
        2 (6.67%) aligned >1 times
95.00% overall alignment rate
"""


class AlignmentReportTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="alignment-report-test-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def parse(self, text=SUMMARY, name="sample.align_summary.txt"):
        path = self.root / name
        path.write_text(text)
        return build_alignment_report.parse_hisat2_summary(path)

    def test_standard_mixed_mapping_report_uses_correct_denominators(self):
        result = self.parse()
        self.assertEqual(result["Total_Read_Pairs"], "100")
        self.assertEqual(result["Concordant_Unique_Pairs"], "70")
        self.assertEqual(result["Concordant_Multi_Pairs"], "10")
        self.assertEqual(result["Concordant_0_Pairs"], "20")
        self.assertEqual(result["Discordant_1_Pairs"], "5")
        # Input discordance percentage uses the 20 concordantly unmapped pairs;
        # the exported percentage consistently uses all 100 input pairs.
        self.assertEqual(result["Discordant_1(%)"], "5.00")
        self.assertEqual(result["Overall_Map(%)"], "95.00")

    def test_missing_or_duplicate_required_fields_fail(self):
        required = [line for line in SUMMARY.splitlines(True)
                    if any(label in line for label in (
                        "reads; of these:", "were paired", "aligned concordantly 0 times\n",
                        "aligned concordantly exactly", "aligned concordantly >", "overall alignment"))]
        for line in required:
            with self.subTest(line=line):
                for changed in (SUMMARY.replace(line, ""), SUMMARY + line):
                    with self.assertRaises(ValueError):
                        self.parse(changed)

    def test_empty_and_unrecognized_reports_fail(self):
        for text in ("", " \n", "HISAT2 failed: missing index\n", "NaN% overall alignment rate\n"):
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.parse(text)

    def test_invalid_nonfinite_and_out_of_range_percentages_fail(self):
        for value in ("nan", "NaN", "inf", "-inf", "-1", "100.01", "not-a-number"):
            for original in ("70.00%", "95.00%"):
                with self.subTest(value=value, field=original), self.assertRaises(ValueError):
                    self.parse(SUMMARY.replace(original, value + "%"))

    def test_negative_and_fractional_counts_fail(self):
        for value in ("-70", "70.5", "NaN"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.parse(SUMMARY.replace("70 (70.00%)", value + " (70.00%)"))

    def test_counts_and_percentages_must_agree(self):
        for old, new in (("70.00%", "71.00%"), ("25.00%", "5.00%"),
                         ("33.33%", "34.00%"), ("95.00%", "94.00%")):
            with self.subTest(field=old), self.assertRaises(ValueError):
                self.parse(SUMMARY.replace(old, new))

    def test_concordant_counts_must_partition_total(self):
        with self.assertRaisesRegex(ValueError, "do not sum"):
            self.parse(SUMMARY.replace("10 (10.00%) aligned concordantly", "11 (11.00%) aligned concordantly"))

    def test_discordant_count_cannot_exceed_concordant_zero(self):
        with self.assertRaises(ValueError):
            self.parse(SUMMARY.replace("5 (25.00%) aligned discordantly", "21 (105.00%) aligned discordantly"))

    def test_inconsistent_nested_totals_and_partial_sections_fail(self):
        variants = [
            SUMMARY.replace("20 pairs aligned concordantly", "19 pairs aligned concordantly"),
            SUMMARY.replace("15 pairs aligned 0 times", "14 pairs aligned 0 times"),
            SUMMARY.replace("30 mates make up", "29 mates make up"),
            SUMMARY.replace("      5 (25.00%) aligned discordantly 1 time\n", ""),
            SUMMARY.replace("      30 mates make up the pairs; of these:\n", ""),
            SUMMARY.replace("        18 (60.00%) aligned exactly 1 time\n", ""),
            SUMMARY.replace("2 (6.67%) aligned >1 times", "3 (10.00%) aligned >1 times"),
        ]
        for number, text in enumerate(variants):
            with self.subTest(variant=number), self.assertRaises(ValueError):
                self.parse(text)

    def test_optional_section_counts_cannot_lose_their_headers(self):
        headers = ["    20 pairs aligned concordantly 0 times; of these:\n",
                   "    15 pairs aligned 0 times concordantly or discordantly; of these:\n"]
        for header in headers:
            with self.subTest(header=header), self.assertRaisesRegex(ValueError, "incomplete"):
                self.parse(SUMMARY.replace(header, ""))

    def test_no_mixed_output_omits_mate_section(self):
        text = SUMMARY.split("    15 pairs aligned")[0] + "85.00% overall alignment rate\n"
        self.assertEqual(self.parse(text)["Overall_Map(%)"], "85.00")

    def test_no_discordant_output_omits_discordant_section(self):
        text = SUMMARY.replace("    20 pairs aligned concordantly 0 times; of these:\n", "")
        text = text.replace("      5 (25.00%) aligned discordantly 1 time\n", "")
        text = text.replace("15 pairs aligned", "20 pairs aligned").replace("30 mates", "40 mates")
        text = text.replace("10 (33.33%)", "10 (25.00%)").replace("18 (60.00%)", "28 (70.00%)")
        text = text.replace("2 (6.67%)", "2 (5.00%)")
        self.assertEqual(self.parse(text)["Discordant_1_Pairs"], "0")

    def test_zero_pairs_report_is_explicitly_validated(self):
        text = """0 reads; of these:
0 (0.00%) were paired; of these:
0 (0.00%) aligned concordantly 0 times
0 (0.00%) aligned concordantly exactly 1 time
0 (0.00%) aligned concordantly >1 times
0.00% overall alignment rate
"""
        self.assertEqual(self.parse(text)["Total_Read_Pairs"], "0")
        with self.assertRaises(ValueError):
            self.parse(text.replace("0.00% overall", "100.00% overall"))

    def test_nonpaired_input_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "all input reads"):
            self.parse(SUMMARY.replace("100 (100.00%) were paired", "90 (90.00%) were paired"))

    def test_sample_id_strips_only_terminal_suffix(self):
        name = "cohort.align_summary.txt.replicate.align_summary.txt"
        self.assertEqual(self.parse(name=name)["SampleID"], "cohort.align_summary.txt.replicate")

    def test_invalid_encoding_is_not_silently_discarded(self):
        path = self.root / "sample.align_summary.txt"
        path.write_bytes(SUMMARY.encode() + b"\xff")
        with self.assertRaises(UnicodeError):
            build_alignment_report.parse_hisat2_summary(path)

    def run_cli(self, paths):
        output = self.root / "report.csv"
        result = subprocess.run(
            [sys.executable, str(SCRIPTS / "build_alignment_report.py"), "--output", str(output)]
            + [str(path) for path in paths], text=True, capture_output=True,
        )
        return result, output

    def test_cli_rejects_bad_input_without_writing_report(self):
        path = self.root / "bad.align_summary.txt"
        path.write_text("")
        result, output = self.run_cli([path])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(str(path), result.stderr)
        self.assertFalse(output.exists())

    def test_cli_rejects_duplicate_sample_ids(self):
        paths = []
        for directory in ("first", "second"):
            parent = self.root / directory
            parent.mkdir()
            path = parent / "sample.align_summary.txt"
            path.write_text(SUMMARY)
            paths.append(path)
        result, output = self.run_cli(paths)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Duplicate sample IDs", result.stderr)
        self.assertFalse(output.exists())

    def test_cli_exports_valid_report_with_unchanged_columns(self):
        path = self.root / "sample.align_summary.txt"
        path.write_text(SUMMARY)
        result, output = self.run_cli([path])
        self.assertEqual(result.returncode, 0, result.stderr)
        with output.open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0], self.parse())


if __name__ == "__main__":
    unittest.main(verbosity=2)
