"""Reference-model contract across normalization, BED12, StringTie and prepDE.

Unit cases use only synthetic GTF. Opt-in installed-tool integrations use the
same RNASEQ_TEST_* variables as test_workflow_regressions.py; no downloads.
"""
import contextlib
import io
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))
import prepare_strand_bed as reference
import prepDE


def record(feature, start, end, gene="Gene|literal", tid="tx|one", strand="+",
           chrom="chr1", name="Gene Alpha", coverage=None):
    attrs = f'gene_id "{gene}";'
    if tid is not None:
        attrs += f' transcript_id "{tid}";'
    if name is not None:
        attrs += f' gene_name "{name}";'
    if coverage is not None:
        attrs += f' cov "{coverage}";'
    phase = "0" if feature == "CDS" else "."
    return f"{chrom}\tfixture\t{feature}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs}\n"


def source_gtf(flavor):
    rows = ["# Synthetic reference with two spliced transcript models\n"]
    for gene, tid, chrom, strand, offset, name in (
            ("Gene|literal", "tx|one", "chr1", "+", 0, "Gene Alpha"),
            ("001", "NA", "chr2", "-", 900, "Gene Beta")):
        metadata = dict(gene=gene, tid=tid, chrom=chrom, strand=strand, name=name)
        if flavor == "full":
            rows.append(record("gene", offset + 101, offset + 400,
                               **dict(metadata, tid=None)))
        if flavor in {"basic", "full"}:
            rows.append(record("transcript", offset + 101, offset + 400, **metadata))
        # Reverse input order ensures all consumers derive genomic exon order.
        rows.append(record("exon", offset + 301, offset + 400, **metadata))
        rows.append(record("exon", offset + 101, offset + 200, **metadata))
        if flavor in {"full", "exon_cds"}:
            rows.append(record("CDS", offset + 121, offset + 200, **metadata))
            rows.append(record("CDS", offset + 301, offset + 380, **metadata))
        if flavor == "full":
            lower_utr, upper_utr = (("five_prime_utr", "three_prime_utr") if strand == "+"
                                    else ("three_prime_utr", "five_prime_utr"))
            lower_codon, upper_codon = (("start_codon", "stop_codon") if strand == "+"
                                       else ("stop_codon", "start_codon"))
            rows.append(record(lower_utr, offset + 101, offset + 120, **metadata))
            rows.append(record(upper_utr, offset + 381, offset + 400, **metadata))
            rows.append(record(lower_codon, offset + 121, offset + 123, **metadata))
            rows.append(record(upper_codon, offset + 378, offset + 380, **metadata))
    return "".join(rows)


class ReferenceModelContractTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="rnaseq-reference-contract-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / "source.gtf"
        self.normalized = self.root / "reference.gtf"
        self.bed = self.root / "reference.bed12"
        self.report = self.root / "reference.json"

    def prepare(self, text):
        self.source.write_text(text)
        before = self.source.read_bytes()
        stderr = io.StringIO()
        with contextlib.redirect_stderr(stderr):
            status = reference.main([
                "--annotation", str(self.source), "--output", str(self.bed),
                "--normalized-gtf", str(self.normalized), "--report", str(self.report)])
        self.prepare_stderr = stderr.getvalue()
        self.assertEqual(self.source.read_bytes(), before)
        return status, json.loads(self.report.read_text())

    def test_four_reference_flavors_share_normalized_bed_and_count_models(self):
        baseline = None
        for flavor in ("basic", "full", "exon_cds", "exon_only"):
            with self.subTest(flavor=flavor):
                status, report = self.prepare(source_gtf(flavor))
                self.assertEqual(status, 0, report)
                self.assertEqual(report["status"], "passed")
                direct = prepDE.parse_gtf(self.source, gene_id_only=True, reference=True)
                normalized = prepDE.parse_gtf(self.normalized, gene_id_only=True, reference=True)
                self.assertEqual(direct, normalized)
                if baseline is None:
                    baseline = direct
                self.assertEqual(direct, baseline)
                self.assertEqual(set(direct), {"tx|one", "NA"})
                for tid, gene, chrom, strand, offset in (
                        ("tx|one", "Gene|literal", "chr1", "+", 0),
                        ("NA", "001", "chr2", "-", 900)):
                    self.assertEqual(direct[tid], {
                        "gene": gene, "canonical_gene": gene, "chrom": chrom,
                        "strand": strand, "start": offset + 101, "end": offset + 400,
                        "coverage": 0.0, "exons": ((offset + 101, offset + 200),
                                                   (offset + 301, offset + 400)),
                        "length": 200})
                beds = {row["name"]: row for row in reference.read_bed12(self.bed)}
                self.assertEqual(beds["tx|one"]["blocks"], [(100, 200), (300, 400)])
                self.assertEqual(beds["NA"]["blocks"], [(1000, 1100), (1200, 1300)])
                for tid, model in direct.items():
                    bed = beds[tid]
                    self.assertEqual((bed["chrom"], bed["strand"], bed["start"], bed["end"]),
                                     (model["chrom"], model["strand"], model["start"] - 1, model["end"]))
                    self.assertEqual(tuple((left + 1, right) for left, right in bed["blocks"]),
                                     model["exons"])
                features = [line.split("\t")[2] for line in self.normalized.read_text().splitlines()
                            if line and not line.startswith("#")]
                self.assertEqual(features.count("transcript"), 2)
                self.assertEqual(features.count("exon"), 4)
                self.assertEqual(set(features), {"transcript", "exon"})
                named_direct = prepDE.parse_gtf(self.source, reference=True)
                named_normalized = prepDE.parse_gtf(self.normalized, reference=True)
                self.assertEqual(named_direct, named_normalized)
                self.assertEqual(named_normalized["tx|one"]["gene"], "Gene|literal|Gene Alpha")
                self.assertEqual(named_normalized["NA"]["gene"], "001|Gene Beta")

    def test_overlapping_duplicate_and_adjacent_exons_use_one_validated_union(self):
        text = (record("transcript", 90, 410) + record("exon", 101, 180)
                + record("exon", 150, 200) + record("exon", 150, 200)
                + record("exon", 201, 220) + record("exon", 301, 400))
        status, report = self.prepare(text)
        self.assertEqual(status, 0, report)
        direct = prepDE.parse_gtf(self.source, gene_id_only=True, reference=True)
        normalized = prepDE.parse_gtf(self.normalized, gene_id_only=True, reference=True)
        self.assertEqual(direct, normalized)
        self.assertEqual(normalized["tx|one"]["exons"], ((101, 220), (301, 400)))
        self.assertEqual(normalized["tx|one"]["length"], 220)
        self.assertEqual((normalized["tx|one"]["start"], normalized["tx|one"]["end"]), (101, 400))
        self.assertEqual(reference.read_bed12(self.bed)[0]["blocks"], [(100, 220), (300, 400)])

    def test_subfeatures_outside_exon_union_fail_both_reference_entry_points(self):
        valid = record("transcript", 101, 400) + record("exon", 101, 200) + record("exon", 301, 400)
        for feature in ("CDS", "UTR", "five_prime_utr", "three_prime_utr", "start_codon", "stop_codon"):
            for start, end in ((250, 252), (199, 301), (401, 403)):
                with self.subTest(feature=feature, interval=(start, end)):
                    # Publish a success first: a subsequent error must not leave
                    # a usable BED/normalized GTF from the preceding attempt.
                    self.assertEqual(self.prepare(valid)[0], 0)
                    status, report = self.prepare(valid + record(feature, start, end))
                    self.assertEqual(status, 1)
                    self.assertEqual(report["status"], "failed")
                    self.assertIn("not contained in exon blocks", report["error"])
                    self.assertFalse(self.bed.exists())
                    self.assertFalse(self.normalized.exists())
                    with self.assertRaisesRegex(ValueError, "not contained in exon blocks"):
                        prepDE.parse_gtf(self.source, gene_id_only=True, reference=True)

    def test_stringtie_short_intron_boundary_is_shared_by_all_model_consumers(self):
        for gap in (1, 4, 5):
            with self.subTest(intron_bases=gap):
                text = record("exon", 101, 150) + record("exon", 151 + gap, 200 + gap)
                status, report = self.prepare(text)
                self.assertEqual(status, 0, report)
                expected_exons = ((101, 200 + gap),) if gap < 5 else ((101, 150), (151 + gap, 200 + gap))
                direct = prepDE.parse_gtf(self.source, gene_id_only=True, reference=True)
                normalized = prepDE.parse_gtf(self.normalized, gene_id_only=True, reference=True)
                self.assertEqual(direct, normalized)
                self.assertEqual(normalized["tx|one"]["exons"], expected_exons)
                self.assertEqual(normalized["tx|one"]["length"], 100 + gap if gap < 5 else 100)
                self.assertEqual(reference.read_bed12(self.bed)[0]["blocks"],
                                 [(start - 1, end) for start, end in expected_exons])
                expected_merges = [{"transcript_id": "tx|one", "chrom": "chr1", "strand": "+",
                                    "start_1based": 151, "end_1based": 150 + gap,
                                    "added_exonic_bases": gap}] if gap < 5 else []
                self.assertEqual(report["short_intron_merges"], expected_merges)
                self.assertEqual(report["short_intron_added_exonic_bases"], gap if gap < 5 else 0)
                if gap < 5:
                    self.assertIn("[WARNING]", self.prepare_stderr)
                    self.assertIn(f"adding {gap} exonic bases", self.prepare_stderr)
                    self.assertIn("Source annotation is unchanged", self.prepare_stderr)
                else:
                    self.assertNotIn("[WARNING]", self.prepare_stderr)

    def test_short_intron_normalization_does_not_hide_invalid_source_subfeatures(self):
        text = record("exon", 101, 150) + record("exon", 155, 204)
        for feature in ("CDS", "UTR", "start_codon", "stop_codon"):
            with self.subTest(feature=feature):
                status, report = self.prepare(text + record(feature, 151, 153))
                self.assertEqual(status, 1)
                self.assertIn("not contained in exon blocks", report["error"])
                with self.assertRaisesRegex(ValueError, "not contained in exon blocks"):
                    prepDE.parse_gtf(self.source, gene_id_only=True, reference=True)

    def test_cds_only_gtf_is_not_mistaken_for_a_complete_exon_model(self):
        status, report = self.prepare(record("transcript", 101, 400) + record("CDS", 121, 380))
        self.assertEqual(status, 1)
        self.assertIn("lack exon", report["error"])
        with self.assertRaisesRegex(ValueError, "lack exon"):
            prepDE.parse_gtf(self.source, gene_id_only=True, reference=True)

    def test_unsupported_stringtie_attributes_fail_before_publishing_models(self):
        for field in ("gene", "tid", "name"):
            for value in ('A"B', 'A\\"B', "A;B", "A\x01B"):
                with self.subTest(field=field, value=value):
                    status, report = self.prepare(record("exon", 101, 200, **{field: value}))
                    self.assertEqual(status, 1, report)
                    self.assertEqual(report["status"], "failed")
                    self.assertFalse(self.normalized.exists())
                    self.assertFalse(self.bed.exists())
                    with self.assertRaises(ValueError):
                        prepDE.parse_gtf(self.source, gene_id_only=True, reference=True)
        gene, tid, name = "Gene|id\\literal", "tx|id\\literal", "Gene|name\\literal"
        status, report = self.prepare(record("exon", 101, 200, gene=gene, tid=tid, name=name))
        self.assertEqual(status, 0, report)
        for path in (self.source, self.normalized):
            canonical = prepDE.parse_gtf(path, gene_id_only=True, reference=True)
            named = prepDE.parse_gtf(path, reference=True)
            self.assertEqual(set(canonical), {tid})
            self.assertEqual(canonical[tid]["canonical_gene"], gene)
            self.assertEqual(canonical[tid]["gene"], gene)
            self.assertEqual(named[tid]["gene"], gene + "|" + name)

    def test_gene_name_fallback_and_transcript_priority_preserve_canonical_gene_id(self):
        for gene_name, exon_name, transcript_name, expected_name in (
                ("GeneFallback", None, None, "GeneFallback"),
                ("GeneFallback", "ExonName", None, "ExonName"),
                ("GeneFallback", "ExonName", "TranscriptName", "TranscriptName"),
                (None, None, None, None)):
            with self.subTest(gene=gene_name, exon=exon_name, transcript=transcript_name):
                text = (record("gene", 101, 400, tid=None, name=gene_name)
                        + record("exon", 101, 200, name=exon_name)
                        + record("transcript", 101, 400, name=transcript_name)
                        + record("exon", 301, 400, name=exon_name))
                status, report = self.prepare(text)
                self.assertEqual(status, 0, report)
                for path in (self.source, self.normalized):
                    canonical = prepDE.parse_gtf(path, gene_id_only=True, reference=True)["tx|one"]
                    named = prepDE.parse_gtf(path, reference=True)["tx|one"]
                    self.assertEqual(canonical["gene"], "Gene|literal")
                    self.assertEqual(canonical["canonical_gene"], "Gene|literal")
                    self.assertEqual(named["gene"], "Gene|literal" + ("|" + expected_name if expected_name else ""))

    def test_gene_id_cannot_span_chromosomes_or_strands_across_exon_only_transcripts(self):
        first = record("exon", 101, 200, tid="tx1")
        for second_chrom, second_strand in (("chr2", "+"), ("chr1", "-"), ("chr2", "-")):
            with self.subTest(chrom=second_chrom, strand=second_strand):
                text = first + record("exon", 501, 600, tid="tx2", chrom=second_chrom,
                                      strand=second_strand)
                status, report = self.prepare(text)
                self.assertEqual(status, 1, report)
                self.assertIn("gene_id", report["error"])
                self.assertIn("inconsistent chromosome or strand across transcripts", report["error"])
                self.assertFalse(self.normalized.exists())
                self.assertFalse(self.bed.exists())
                with self.assertRaisesRegex(ValueError, "inconsistent chromosome or strand across transcripts"):
                    prepDE.parse_gtf(self.source, gene_id_only=True, reference=True)
        # Disjoint transcript intervals alone are valid for the same gene and
        # must remain distinct models when chromosome and strand agree.
        status, report = self.prepare(first + record("exon", 501, 600, tid="tx2"))
        self.assertEqual(status, 0, report)
        direct = prepDE.parse_gtf(self.source, gene_id_only=True, reference=True)
        normalized = prepDE.parse_gtf(self.normalized, gene_id_only=True, reference=True)
        self.assertEqual(direct, normalized)
        self.assertEqual(set(normalized), {"tx1", "tx2"})
        self.assertEqual({model["canonical_gene"] for model in normalized.values()}, {"Gene|literal"})
        self.assertEqual(normalized["tx1"]["exons"], ((101, 200),))
        self.assertEqual(normalized["tx2"]["exons"], ((501, 600),))
        self.assertEqual([row["blocks"] for row in reference.read_bed12(self.bed)],
                         [[(100, 200)], [(500, 600)]])

    def test_stringtie_sample_parser_remains_strict_after_reference_split(self):
        provenance = "# stringtie input.bam -e -G reference.gtf\n"
        rows = (record("transcript", 101, 400, coverage="2.5")
                + record("exon", 101, 200, coverage="2.5")
                + record("exon", 301, 400, coverage="2.5"))
        self.source.write_text(provenance + rows)
        model = prepDE.parse_gtf(self.source, gene_id_only=True)["tx|one"]
        self.assertEqual(model["coverage"], 2.5)
        self.assertEqual(model["length"], 200)
        cases = {
            "missing provenance": rows,
            "missing coverage": provenance + rows.replace(' cov "2.5";', ''),
            "nonfinite coverage": provenance + rows.replace('cov "2.5"', 'cov "NaN"'),
            "negative coverage": provenance + rows.replace('cov "2.5"', 'cov "-1"'),
            "unexpected CDS": provenance + rows + record("CDS", 121, 200),
            "unexpected codon": provenance + rows + record("start_codon", 121, 123),
            "missing transcript": provenance + "".join(line for line in rows.splitlines(True)
                                                        if "\ttranscript\t" not in line),
            "duplicate exon": provenance + rows + record("exon", 101, 200),
        }
        for label, text in cases.items():
            with self.subTest(case=label):
                self.source.write_text(text)
                with self.assertRaises(ValueError):
                    prepDE.parse_gtf(self.source, gene_id_only=True)


@unittest.skipUnless(os.environ.get("RNASEQ_TEST_FULL_FASTQ") == "1",
                     "Set RNASEQ_TEST_FULL_FASTQ=1 for real reference-to-Excel regressions")
class ReferenceWorkflowContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import test_workflow_regressions as regression
        candidate = os.environ.get("RNASEQ_TEST_FASTP") or shutil.which("fastp")
        if not candidate:
            raise unittest.SkipTest("Requires an existing fastp executable")
        cls.regression = regression
        # Composition avoids unittest discovering inherited regression methods.
        regression.FullFastqRegressionTests.setUpClass()
        cls.fixture = regression.FullFastqRegressionTests()
        (cls.fixture.root / "bin" / "fastp").symlink_to(Path(candidate).resolve())

    @classmethod
    def tearDownClass(cls):
        cls.regression.FullFastqRegressionTests.tearDownClass()

    def test_full_and_exon_only_references_reach_excel_in_both_modes(self):
        counts = {}
        for stranded in (False, True):
            for flavor in ("full", "exon_only"):
                with self.subTest(stranded=stranded, flavor=flavor):
                    work, source = self.fixture.make_case(f"reference-{stranded}-{flavor}", stranded)
                    annotation = work / "references" / "reference.gtf"
                    rows = []
                    for gene, symbol, start, strand in self.regression.GENES:
                        metadata = dict(gene=gene, tid=gene + ".t1", chrom="chrMini",
                                        strand=strand, name=symbol)
                        if flavor == "full":
                            lower_utr, upper_utr = (("five_prime_utr", "three_prime_utr") if strand == "+"
                                                    else ("three_prime_utr", "five_prime_utr"))
                            lower_codon, upper_codon = (("start_codon", "stop_codon") if strand == "+"
                                                       else ("stop_codon", "start_codon"))
                            rows += [record("gene", start + 1, start + 1400, **dict(metadata, tid=None)),
                                     record("transcript", start + 1, start + 1400, **metadata),
                                     record(lower_utr, start + 1, start + 100, **metadata),
                                     record("CDS", start + 101, start + 1300, **metadata),
                                     record(upper_utr, start + 1301, start + 1400, **metadata),
                                     record(lower_codon, start + 101, start + 103, **metadata),
                                     record(upper_codon, start + 1298, start + 1300, **metadata)]
                        rows.append(record("exon", start + 1, start + 1400, **metadata))
                    annotation.write_text("".join(rows))
                    original = annotation.read_bytes()
                    status, output = self.fixture.run_workflow(work, source, 3, "reference-workflow")
                    self.assertEqual(status, 0, output[-16000:])
                    quant = self.fixture.assert_outputs(work, stranded, 3)
                    self.assertEqual(annotation.read_bytes(), original)
                    normalized = quant / "reference" / "reference.gtf"
                    self.assertEqual(prepDE.parse_gtf(annotation, gene_id_only=True, reference=True),
                                     prepDE.parse_gtf(normalized, gene_id_only=True, reference=True))
                    counts[(stranded, flavor)] = (quant / "gene_count_matrix.csv").read_bytes()
            self.assertEqual(counts[(stranded, "full")], counts[(stranded, "exon_only")])

    def test_invalid_shared_reference_stops_before_cleaning_or_alignment(self):
        for raw in (False, True):
            with self.subTest(raw=raw):
                work, source = self.fixture.make_case("bad-shared-reference-" + str(raw), raw, {"001": 75})
                annotation = work / "references" / "reference.gtf"
                annotation.write_text(annotation.read_text() + record(
                    "CDS", 2501, 2503, gene="001", tid="001.t1", chrom="chrMini", name="FirstGene"))
                if raw:
                    raw_dir = work / "rawdata" / "001"
                    raw_dir.mkdir(parents=True)
                    for mate in (1, 2):
                        (work / "cleandata" / "001" / f"001_{mate}.clean.fastq.gz").rename(
                            raw_dir / f"001_{mate}.fastq.gz")
                    shutil.rmtree(work / "cleandata")
                    config_path = work / "config.json"
                    config = json.loads(config_path.read_text())
                    config["fastp"] = {"enabled": True, "trim_poly_g": False}
                    config["paths"].update(raw_dir="rawdata", qc_dir="qc/fastp")
                    config_path.write_text(json.dumps(config, indent=2))
                inputs = list((work / ("rawdata" if raw else "cleandata")).rglob("*.fastq.gz"))
                original = {path: path.read_bytes() for path in inputs}
                status, output = self.fixture.run_workflow(work, source, 4, "invalid-reference")
                self.assertNotEqual(status, 0, output[-16000:])
                logs = work / "logs"
                if raw:
                    logs /= "stranded"
                report = json.loads((logs / "reference" / "prepare_bed.json").read_text())
                self.assertEqual(report["status"], "failed")
                self.assertIn("not contained in exon blocks", report["error"])
                self.assertNotIn("localrule fastp_clean:", output)
                self.assertNotIn("localrule measure_read_lengths:", output)
                self.assertNotIn("localrule hisat2_align_sort:", output)
                self.assertFalse((work / "logs" / "001" / "fastp.log").exists())
                self.assertFalse((work / "logs" / "001" / "hisat2_align_sort.log").exists())
                self.assertFalse((work / "alignment" / "001" / "001.sorted.bam").exists())
                if raw:
                    self.assertFalse(list((work / "cleandata").rglob("*.fastq.gz")))
                    self.assertFalse((work / "qc" / "fastp" / "001" / "001.fastp.json").exists())
                self.assertEqual({path: path.read_bytes() for path in inputs}, original)

    def test_reference_only_rerun_preserves_raw_clean_read_lengths_and_bams(self):
        work, source = self.fixture.make_case("reference-only-rerun", False)
        raw_inputs = []
        for sample in self.regression.SAMPLES:
            raw_dir = work / "rawdata" / sample
            raw_dir.mkdir(parents=True)
            for mate in (1, 2):
                raw = raw_dir / f"{sample}_{mate}.fastq.gz"
                (work / "cleandata" / sample / f"{sample}_{mate}.clean.fastq.gz").rename(raw)
                raw_inputs.append(raw)
        shutil.rmtree(work / "cleandata")
        config_path = work / "config.json"
        config = json.loads(config_path.read_text())
        config["fastp"] = {"enabled": True, "trim_poly_g": False}
        config["paths"].update(raw_dir="rawdata", qc_dir="qc/fastp")
        config["alignment"].pop("mapq_filter")
        config_path.write_text(json.dumps(config, indent=2))
        status, output = self.fixture.run_workflow(work, source, 3, "initial-raw-workflow")
        self.assertEqual(status, 0, output[-16000:])
        quant = self.fixture.assert_outputs(work, False, 3)
        expected_counts = (quant / "gene_count_matrix.csv").read_bytes()
        artifacts = raw_inputs + [work / "alignment" / "read_lengths.json"]
        for sample in self.regression.SAMPLES:
            artifacts.extend((work / "cleandata" / sample).glob("*.fastq.gz"))
            artifacts += [work / "alignment" / "read_lengths" / f"{sample}.json",
                          work / "alignment" / sample / f"{sample}.sorted.bam",
                          work / "alignment" / sample / f"{sample}.sorted.bam.bai"]
        self.assertEqual(len(artifacts), 15)
        before = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in artifacts}
        annotation = work / "references" / "reference.gtf"
        annotation.write_text(annotation.read_text() + "# Added reference annotation comment for rerun regression\n")
        status, output = self.fixture.run_workflow(work, source, 3, "changed-reference-workflow")
        self.assertEqual(status, 0, output[-16000:])
        self.assertIn("localrule prepare_reference:", output)
        for rule in ("fastp_clean", "measure_read_lengths", "hisat2_align_sort"):
            self.assertNotIn("localrule " + rule + ":", output)
        self.fixture.assert_outputs(work, False, 3)
        self.assertEqual((quant / "gene_count_matrix.csv").read_bytes(), expected_counts)
        self.assertEqual({path: (path.read_bytes(), path.stat().st_mtime_ns) for path in artifacts}, before)


if __name__ == "__main__":
    unittest.main(verbosity=2)
