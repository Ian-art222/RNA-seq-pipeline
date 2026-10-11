"""Scientific validation tests for transcript BED12 preparation (standard library only)."""
import importlib.util
import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/prepare_strand_bed.py"
SPEC = importlib.util.spec_from_file_location("prepare_strand_bed", SCRIPT)
bed = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bed)


def gtf(feature, start, end, tid="tx1", strand="+", chrom="chr1", gene="g1"):
    return f'{chrom}\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{gene}"; transcript_id "{tid}";\n'


class StrandBedTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.annotation = self.root / "reference.gtf"
        self.output = self.root / "reference.bed12"
        self.report = self.root / "report.json"

    def run_prepare(self, text, *arguments):
        self.annotation.write_text(text)
        result = bed.main(["--annotation", str(self.annotation), "--output", str(self.output),
                           "--report", str(self.report), *arguments])
        return result, json.loads(self.report.read_text())

    def test_plus_minus_multiexon_coordinates(self):
        text = gtf("exon", 101, 150) + gtf("exon", 201, 250)
        text += gtf("exon", 401, 450, "tx2", "-") + gtf("exon", 301, 350, "tx2", "-")
        code, report = self.run_prepare(text)
        self.assertEqual(code, 0)
        records = bed.read_bed12(self.output)
        self.assertEqual(records[0]["blocks"], [(100, 150), (200, 250)])
        self.assertEqual(records[1]["blocks"], [(300, 350), (400, 450)])
        self.assertEqual(records[1]["strand"], "-")
        self.assertEqual(report["transcript_count"], 2)
        self.assertEqual(self.annotation.read_text(), text)

    def test_duplicate_and_overlapping_exons_are_merged(self):
        code, report = self.run_prepare(gtf("exon", 1, 20) * 2 + gtf("exon", 10, 30))
        self.assertEqual(code, 0)
        self.assertEqual(bed.read_bed12(self.output)[0]["blocks"], [(0, 30)])
        self.assertEqual(report["input_exon_records"], 3)
        self.assertEqual(report["merged_exon_blocks"], 1)

    def test_conflicting_transcript_identity_fails(self):
        for change in ({"strand": "-"}, {"chrom": "chr2"}, {"gene": "g2"}):
            with self.subTest(change=change):
                code, report = self.run_prepare(gtf("exon", 1, 20) + gtf("exon", 101, 120, **change))
                self.assertEqual(code, 1)
                self.assertIn("inconsistent", report["error"])

    def test_missing_exons_are_not_silently_omitted_or_inferred(self):
        for feature in ("transcript", "CDS", "five_prime_utr"):
            with self.subTest(feature=feature):
                code, report = self.run_prepare(gtf("exon", 1, 20) + gtf(feature, 101, 120, "incomplete"))
                self.assertEqual(code, 1)
                self.assertIn("lack exon", report["error"])

    def test_empty_annotation_fails(self):
        code, report = self.run_prepare("# annotation without models\n")
        self.assertEqual(code, 1)
        self.assertIn("no transcript exon", report["error"])

    def test_invalid_coordinates_and_strand_fail(self):
        for record in (gtf("exon", 0, 20), gtf("exon", 20, 10), gtf("exon", 1, 20, strand=".")):
            with self.subTest(record=record):
                code, report = self.run_prepare(record)
                self.assertEqual(code, 1)
                self.assertEqual(report["status"], "failed")

    def test_exons_outside_declared_transcript_fail(self):
        code, report = self.run_prepare(gtf("transcript", 10, 20) + gtf("exon", 1, 20))
        self.assertEqual(code, 1)
        self.assertIn("outside declared", report["error"])

    def test_missing_identifier_fails(self):
        code, report = self.run_prepare('chr1\ttest\texon\t1\t20\t.\t+\t.\tgene_id "g1";\n')
        self.assertEqual(code, 1)
        self.assertIn("requires transcript_id", report["error"])

    def test_external_bed_exact_match_passes_mismatch_fails(self):
        external = self.root / "provided.bed"
        external.write_text("chr1\t0\t20\ttx1\t0\t+\t0\t0\t0\t1\t20,\t0,\n")
        code, report = self.run_prepare(gtf("exon", 1, 20), "--bed12", str(external))
        self.assertEqual(code, 0)
        external.write_text(external.read_text().replace("\t+\t", "\t-\t"))
        code, report = self.run_prepare(gtf("exon", 1, 20), "--bed12", str(external))
        self.assertEqual(code, 1)
        self.assertIn("1 inconsistent", report["error"])

    def test_bed_rejects_bad_blocks_and_empty_file(self):
        examples = ["", "chr1\t0\t20\ttx1\t0\t+\t0\t0\t0\t2\t10,\t0,\n",
                    "chr1\t0\t20\ttx1\t0\t+\t0\t0\t0\t2\t15,10,\t0,10,\n",
                    "chr1\t0\t20\ttx1\t0\t+\t0\t0\t0\t1\t19,\t0,\n"]
        for text in examples:
            with self.subTest(text=text):
                self.output.write_text(text)
                with self.assertRaises(ValueError):
                    bed.read_bed12(self.output)

    def test_gff3_conversion_is_recorded_and_models_checked(self):
        self.annotation = self.root / "reference.gff3"
        source = "##gff-version 3\nchr1\ttest\tmRNA\t1\t20\t.\t+\t.\tID=tx1;Parent=g1\n"
        source += "chr1\ttest\tCDS\t1\t20\t.\t+\t0\tParent=tx1\n"
        def convert(command, **kwargs):
            Path(command[-1]).write_text(gtf("exon", 1, 20))
            return type("Result", (), {"returncode": 0})()
        with patch.object(bed.subprocess, "run", side_effect=convert):
            code, report = self.run_prepare(source)
        self.assertEqual(code, 0)
        self.assertIn("CDS/UTR", report["exon_model_source"])
        self.assertTrue(Path(report["conversion_log"]).is_file())
        self.assertEqual(self.annotation.read_text(), source)

    def test_gff3_normalization_cannot_silently_drop_transcripts(self):
        self.annotation = self.root / "reference.gff3"
        source = "chr1\ttest\tmRNA\t1\t20\t.\t+\t.\tID=tx1\nchr1\ttest\tmRNA\t30\t40\t.\t+\t.\tID=dropped\n"
        source += "chr1\ttest\texon\t1\t20\t.\t+\t.\tParent=tx1\nchr1\ttest\texon\t30\t40\t.\t+\t.\tParent=dropped\n"
        def convert(command, **kwargs):
            Path(command[-1]).write_text(gtf("exon", 1, 20))
            return type("Result", (), {"returncode": 0})()
        with patch.object(bed.subprocess, "run", side_effect=convert):
            code, report = self.run_prepare(source)
        self.assertEqual(code, 1)
        self.assertIn("omitted 1 source", report["error"])

    def test_missing_gffread_retains_failure_report_and_log(self):
        self.annotation = self.root / "reference.gff3"
        with patch.object(bed.subprocess, "run", side_effect=FileNotFoundError("gffread")):
            code, report = self.run_prepare("chr1\ttest\tmRNA\t1\t20\t.\t+\t.\tID=tx1\nchr1\ttest\texon\t1\t20\t.\t+\t.\tParent=tx1\n")
        self.assertEqual(code, 1)
        self.assertIn("requires gffread", report["error"])
        self.assertTrue(Path(report["conversion_log"]).is_file())


    def test_gff3_incomplete_gene_and_transcript_not_silently_omitted(self):
        self.annotation = self.root / "reference.gff3"
        complete = "chr1\ttest\tgene\t1\t20\t.\t+\t.\tID=g1\n"
        complete += "chr1\ttest\tmRNA\t1\t20\t.\t+\t.\tID=tx1;Parent=g1\nchr1\ttest\texon\t1\t20\t.\t+\t.\tParent=tx1\n"
        incomplete = ["chr1\ttest\tgene\t101\t200\t.\t+\t.\tID=incomplete\n",
                      "chr1\ttest\tmRNA\t101\t200\t.\t+\t.\tID=incomplete\n"]
        for record in incomplete:
            with self.subTest(record=record), patch.object(bed.subprocess, "run") as converter:
                code, report = self.run_prepare(complete + record)
                self.assertEqual(code, 1)
                self.assertIn("incomplete", report["error"])
                converter.assert_not_called()

    @unittest.skipUnless(shutil.which("gffread"), "gffread is unavailable on PATH")
    def test_real_gffread_plus_minus_multiexon_conversion(self):
        self.annotation = self.root / "reference.gff3"
        source = "##gff-version 3\n"
        for tid, gene, strand, offset in (("tx1", "g1", "+", 0), ("tx2", "g2", "-", 300)):
            source += f"chr1\ttest\tgene\t{offset + 1}\t{offset + 150}\t.\t{strand}\t.\tID={gene}\n"
            source += f"chr1\ttest\tmRNA\t{offset + 1}\t{offset + 150}\t.\t{strand}\t.\tID={tid};Parent={gene}\n"
            for left, right in ((1, 50), (101, 150)):
                source += f"chr1\ttest\texon\t{offset + left}\t{offset + right}\t.\t{strand}\t.\tParent={tid}\n"
        code, report = self.run_prepare(source)
        self.assertEqual(code, 0, report)
        records = bed.read_bed12(self.output)
        self.assertEqual(records[0]["blocks"], [(0, 50), (100, 150)])
        self.assertEqual(records[1]["blocks"], [(300, 350), (400, 450)])
        self.assertEqual(records[1]["strand"], "-")
        self.assertEqual(report["source_models_with_explicit_exons"], 2)
        self.assertEqual(report["source_models_without_explicit_exons"], 0)
        self.assertTrue(Path(report["normalized_annotation"]).is_file())
        self.assertEqual(self.annotation.read_text(), source)

    @unittest.skipUnless(shutil.which("gffread"), "gffread is unavailable on PATH")
    def test_real_gffread_cdsonly_requires_normalized_exons_and_records_origin(self):
        self.annotation = self.root / "reference.gff3"
        source = "chr1\ttest\tgene\t1\t150\t.\t+\t.\tID=g1\n"
        source += "chr1\ttest\tmRNA\t1\t150\t.\t+\t.\tID=tx1;Parent=g1\n"
        source += "chr1\ttest\tCDS\t1\t50\t.\t+\t0\tParent=tx1\nchr1\ttest\tCDS\t101\t150\t.\t+\t1\tParent=tx1\n"
        code, report = self.run_prepare(source)
        self.assertEqual(report["source_models_with_explicit_exons"], 0)
        self.assertEqual(report["source_models_without_explicit_exons"], 1)
        self.assertIn("gffread-normalized GFF3", report["exon_model_source"])
        self.assertEqual(report["conversion_exit_code"], 0)
        normalized = Path(report["normalized_annotation"]).read_text()
        # gffread versions differ in whether CDS-only GFF3 produces exon records.
        # Accept only an actual normalized exon model, never silently infer one here.
        if "\texon\t" in normalized:
            self.assertEqual(code, 0, report)
            self.assertEqual(bed.read_bed12(self.output)[0]["blocks"], [(0, 50), (100, 150)])
        else:
            self.assertEqual(code, 1)
            self.assertIn("lack exon records", report["error"])
            self.assertFalse(self.output.exists())
        self.assertTrue(Path(report["conversion_log"]).is_file())
        self.assertEqual(self.annotation.read_text(), source)


    def test_gff3_source_children_must_match_declared_transcript(self):
        self.annotation = self.root / "reference.gff3"
        declaration = "chr1\ttest\tmRNA\t1\t50\t.\t+\t.\tID=tx1;Parent=g1\n"
        supported = "chr1\ttest\texon\t1\t50\t.\t+\t.\tParent=tx1\n"
        for feature in ("exon", "CDS", "five_prime_UTR"):
            for chrom, strand, start, end, message in (
                ("chr1", "+", 51, 60, "outside declared"),
                ("chr2", "+", 1, 20, "chromosome/strand disagrees"),
                ("chr1", "-", 1, 20, "chromosome/strand disagrees"),
            ):
                child = f"{chrom}\ttest\t{feature}\t{start}\t{end}\t.\t{strand}\t.\tParent=tx1\n"
                for source in (declaration + supported + child, supported + child + declaration):
                    with self.subTest(feature=feature, source=source), patch.object(bed.subprocess, "run") as converter:
                        code, report = self.run_prepare(source)
                        self.assertEqual(code, 1)
                        self.assertIn(message, report["error"])
                        converter.assert_not_called()

    def test_gff3_valid_child_before_transcript_is_supported(self):
        self.annotation = self.root / "reference.gff3"
        source = "chr1\ttest\texon\t1\t20\t.\t+\t.\tParent=tx1\n"
        source += "chr1\ttest\tmRNA\t1\t20\t.\t+\t.\tID=tx1;Parent=g1\n"
        def convert(command, **kwargs):
            Path(command[-1]).write_text(gtf("exon", 1, 20))
            return type("Result", (), {"returncode": 0})()
        with patch.object(bed.subprocess, "run", side_effect=convert):
            code, report = self.run_prepare(source)
        self.assertEqual(code, 0, report)
        self.assertEqual(bed.read_bed12(self.output)[0]["blocks"], [(0, 20)])


if __name__ == "__main__":
    unittest.main()
