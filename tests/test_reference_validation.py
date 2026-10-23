"""BAM/reference validation shared by stranded and unstranded quantification."""
import contextlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts"))
import validate_bam_reference as validation

try:
    import pysam
except ImportError:
    pysam = None


@unittest.skipUnless(pysam, "pysam is required for real BAM validation")
class ReferenceValidationTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.bam = self.root / "sample.bam"
        self.bai = Path(str(self.bam) + ".bai")
        self.bed = self.root / "reference.bed12"
        self.report = self.root / "validation.json"
        self.bed.write_text("chr1\t0\t20\ttx1\t0\t+\t0\t0\t0\t1\t20,\t0,\n")
        self.make_bam()

    def make_bam(self, sort_order="coordinate", program=None):
        header = {"HD": {"VN": "1.6", "SO": sort_order}, "SQ": [{"SN": "chr1", "LN": 100}]}
        if program:
            header["PG"] = [program]
        with pysam.AlignmentFile(str(self.bam), "wb", header=header) as bam:
            read = pysam.AlignedSegment()
            read.query_name = "read1"
            read.query_sequence = "A" * 20
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 0
            read.mapping_quality = 60
            read.cigarstring = "20M"
            bam.write(read)
        pysam.index(str(self.bam))

    def run_validation(self):
        with contextlib.redirect_stderr(io.StringIO()):
            code = validation.main(["--bam", str(self.bam), "--bai", str(self.bai),
                                    "--bed12", str(self.bed), "--output", str(self.report)])
        return code, json.loads(self.report.read_text())

    def test_matching_reference_passes_without_enforcing_library_direction(self):
        self.make_bam(program={"ID": "hisat2", "CL": "hisat2 --rna-strandness RF -x ref"})
        code, report = self.run_validation()
        self.assertEqual(code, 0, report)
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["validation"]["bed12_transcripts"], 1)

    def test_missing_contig_replaces_stale_pass_with_failure(self):
        self.assertEqual(self.run_validation()[0], 0)
        self.bed.write_text(self.bed.read_text().replace("chr1", "chr2"))
        code, report = self.run_validation()
        self.assertEqual(code, 1)
        self.assertEqual(report["status"], "FAIL")
        self.assertIn("absent from BAM", report["errors"][0])

    def test_out_of_bounds_coordinates_fail(self):
        self.bed.write_text("chr1\t90\t110\ttx1\t0\t+\t90\t90\t0\t1\t20,\t0,\n")
        code, report = self.run_validation()
        self.assertEqual(code, 1)
        self.assertIn("exceeds BAM", report["errors"][0])

    def test_coordinate_sort_declaration_required(self):
        self.make_bam(sort_order="queryname")
        code, report = self.run_validation()
        self.assertEqual(code, 1)
        self.assertIn("coordinate sorting", report["errors"][0])

    def test_corrupt_index_and_truncated_bam_fail(self):
        self.bai.write_bytes(b"invalid-index")
        self.assertEqual(self.run_validation()[0], 1)
        self.make_bam()
        self.bam.write_bytes(self.bam.read_bytes()[:-28])
        self.assertEqual(self.run_validation()[0], 1)

    def test_hardlinked_output_cannot_overwrite_bam(self):
        before = self.bam.read_bytes()
        self.report.hardlink_to(self.bam)
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()):
            self.run_validation()
        self.assertEqual(self.bam.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
