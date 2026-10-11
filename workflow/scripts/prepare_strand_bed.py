#!/usr/bin/env python3
"""Prepare checked transcript BED12 models for RSeQC from the quantification annotation."""
import argparse
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from urllib.parse import unquote


GTF_ATTRIBUTE = re.compile(r'\s*([^\s;]+)\s+(?:"([^"\\]*(?:\\.[^"\\]*)*)"|([^;\s]+))\s*(?:;|$)')
MODEL_FEATURES = {"exon", "cds", "utr", "five_prime_utr", "three_prime_utr", "transcript", "mrna"}


def annotation_rows(path):
    """Yield checked nine-column annotation records, preserving line numbers."""
    with Path(path).open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{number}: expected nine annotation columns")
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise ValueError(f"{path}:{number}: annotation coordinates must be integers") from exc
            if start < 1 or end < start:
                raise ValueError(f"{path}:{number}: invalid one-based annotation interval {start}-{end}")
            if not fields[0] or any(char.isspace() for char in fields[0]):
                raise ValueError(f"{path}:{number}: invalid reference sequence name")
            yield number, fields, start - 1, end


def parse_gtf_attributes(text, location):
    attributes = {}
    position = 0
    while position < len(text):
        if not text[position:].strip():
            break
        match = GTF_ATTRIBUTE.match(text, position)
        if not match:
            raise ValueError(f"{location}: malformed GTF attributes near {text[position:position + 60]!r}")
        key = match.group(1)
        value = match.group(2) if match.group(2) is not None else match.group(3)
        # Repeated metadata such as tag is valid; model identifiers must be unique.
        if key in {"gene_id", "transcript_id"} and key in attributes:
            raise ValueError(f"{location}: repeated {key}")
        attributes[key] = value
        position = match.end()
    return attributes


def read_gtf_models(path):
    """Read exon models; never infer missing exons from CDS or transcript spans."""
    models = {}
    input_exons = 0
    for number, fields, start, end in annotation_rows(path):
        feature = fields[2].lower()
        location = f"{path}:{number}"
        attributes = parse_gtf_attributes(fields[8], location)
        tid = attributes.get("transcript_id")
        if not tid:
            if feature in MODEL_FEATURES or feature.endswith("_utr"):
                raise ValueError(f"{location}: {fields[2]} requires transcript_id and gene_id")
            continue
        gene = attributes.get("gene_id")
        if not gene or any(char.isspace() for char in tid) or not tid:
            raise ValueError(f"{location}: transcript requires a nonempty gene_id and a whitespace-free transcript_id")
        if fields[6] not in {"+", "-"}:
            raise ValueError(f"{location}: transcript {tid!r} must have strand '+' or '-'")
        identity = (fields[0], fields[6], gene)
        if tid not in models:
            models[tid] = {"identity": identity, "exons": [], "spans": []}
        model = models[tid]
        if model["identity"] != identity:
            raise ValueError(f"{location}: transcript_id {tid!r} has inconsistent chromosome, strand, or gene_id")
        if feature == "exon":
            model["exons"].append((start, end))
            input_exons += 1
        elif feature in {"transcript", "mrna"}:
            model["spans"].append((start, end))
    if not models:
        raise ValueError(f"{path}: no transcript exon models found; provide a GTF with exon records")
    missing = sorted(tid for tid, model in models.items() if not model["exons"])
    if missing:
        examples = ", ".join(missing[:5])
        raise ValueError(f"{path}: {len(missing)} transcript(s) lack exon records ({examples}); supply explicit exon models")
    records = []
    for tid, model in models.items():
        blocks = []
        for start, end in sorted(model["exons"]):
            if blocks and start <= blocks[-1][1]:
                blocks[-1] = (blocks[-1][0], max(end, blocks[-1][1]))
            else:
                blocks.append((start, end))
        start, end = blocks[0][0], blocks[-1][1]
        if any(start < left or end > right for left, right in model["spans"]):
            raise ValueError(f"{path}: exon blocks extend outside declared transcript {tid!r}")
        chrom, strand, _ = model["identity"]
        records.append({"chrom": chrom, "start": start, "end": end, "name": tid, "strand": strand, "blocks": blocks})
    records.sort(key=lambda row: (row["chrom"], row["start"], row["end"], row["name"]))
    return records, {"input_exon_records": input_exons,
                     "merged_exon_blocks": sum(len(row["blocks"]) for row in records)}


def read_bed12(path):
    """Return validated BED12 records with absolute, zero-based half-open blocks."""
    records, names = [], set()
    with Path(path).open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith(("#", "track ", "browser ")):
                continue
            location = f"{path}:{number}"
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"{location}: expected at least 12 tab-separated BED columns")
            try:
                start, end = int(fields[1]), int(fields[2])
                thick_start, thick_end = int(fields[6]), int(fields[7])
                count = int(fields[9])
                sizes = [int(item) for item in fields[10].rstrip(",").split(",")]
                offsets = [int(item) for item in fields[11].rstrip(",").split(",")]
            except ValueError as exc:
                raise ValueError(f"{location}: invalid integer in BED coordinates or exon blocks") from exc
            if start < 0 or end <= start or not (start <= thick_start <= thick_end <= end):
                raise ValueError(f"{location}: invalid BED transcript or thick interval")
            if not fields[0] or any(char.isspace() for char in fields[0]):
                raise ValueError(f"{location}: invalid reference sequence name")
            if fields[3] in {"", "."} or any(char.isspace() for char in fields[3]):
                raise ValueError(f"{location}: BED requires a whitespace-free transcript name")
            if fields[3] in names:
                raise ValueError(f"{location}: duplicate transcript name {fields[3]!r}")
            if fields[5] not in {"+", "-"}:
                raise ValueError(f"{location}: BED strand must be '+' or '-'")
            if count < 1 or len(sizes) != count or len(offsets) != count:
                raise ValueError(f"{location}: BED blockCount does not match positive exon block count")
            blocks = [(start + offset, start + offset + size) for offset, size in zip(offsets, sizes)]
            if any(size <= 0 or offset < 0 for size, offset in zip(sizes, offsets)):
                raise ValueError(f"{location}: exon sizes must be positive and offsets nonnegative")
            if blocks[0][0] != start or blocks[-1][1] != end:
                raise ValueError(f"{location}: first/last exon block must match transcript bounds")
            if any(left < start or right > end for left, right in blocks) or any(
                blocks[index][0] < blocks[index - 1][1] for index in range(1, count)
            ):
                raise ValueError(f"{location}: exon blocks must be sorted, nonoverlapping, and inside the transcript")
            names.add(fields[3])
            records.append({"chrom": fields[0], "start": start, "end": end,
                            "name": fields[3], "strand": fields[5], "blocks": blocks})
    if not records:
        raise ValueError(f"{path}: no BED12 transcript records found")
    return records


def write_bed12(records, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, prefix=path.name + ".", delete=False) as handle:
            temporary = Path(handle.name)
            for row in records:
                blocks = row["blocks"]
                fields = [row["chrom"], row["start"], row["end"], row["name"], 0, row["strand"],
                          row["start"], row["start"], 0, len(blocks),
                          ",".join(str(right - left) for left, right in blocks) + ",",
                          ",".join(str(left - row["start"]) for left, _ in blocks) + ","]
                handle.write("\t".join(map(str, fields)) + "\n")
        read_bed12(temporary)
        temporary.replace(path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def annotation_format(path):
    for _, fields, _, _ in annotation_rows(path):
        attributes = fields[8]
        if re.search(r'(?:^|;)\s*(?:gene_id|transcript_id)\s+', attributes):
            return "gtf"
        if "=" in attributes:
            return "gff3"
        break
    if Path(path).suffix.lower() in {".gff", ".gff3"}:
        return "gff3"
    return "gtf"


def gff3_transcripts(path, report):
    """Audit source model support before gffread can synthesize or omit records."""
    features, gene_ids, expected = [], set(), set()
    exon_models, cds_utr_models, identities = set(), set(), {}
    declared_transcripts, child_records = {}, []
    for number, fields, start, end in annotation_rows(path):
        attributes = {}
        for item in fields[8].split(";"):
            if not item or item == ".":
                continue
            if "=" not in item:
                raise ValueError(f"{path}:{number}: malformed GFF3 attributes")
            key, value = item.split("=", 1)
            if key in attributes and key in {"ID", "Parent"}:
                raise ValueError(f"{path}:{number}: repeated GFF3 {key}")
            attributes[key] = value
        feature = fields[2].lower()
        tid = unquote(attributes.get("ID", ""))
        parents = {unquote(value) for value in attributes.get("Parent", "").split(",") if value}
        is_gene = feature in {"gene", "pseudogene"} or feature.endswith("_gene")
        is_transcript = feature in {"transcript", "mrna"} or feature.endswith("rna") or feature.endswith("_transcript")
        is_exon_feature = feature in {"exon", "cds", "utr"} or feature.endswith("_utr")
        if is_gene or is_transcript:
            if not tid:
                raise ValueError(f"{path}:{number}: {fields[2]} requires a GFF3 ID")
            identity = (fields[0], fields[6], frozenset(parents))
            if tid in identities and identities[tid] != identity:
                raise ValueError(f"{path}:{number}: GFF3 ID {tid!r} has inconsistent chromosome, strand, or Parent")
            identities[tid] = identity
        if (is_transcript or is_exon_feature) and fields[6] not in {"+", "-"}:
            raise ValueError(f"{path}:{number}: transcript/exon models require strand '+' or '-'")
        if is_gene:
            gene_ids.add(tid)
        if is_transcript:
            expected.add(tid)
            declared_transcripts.setdefault(tid, []).append((fields[0], fields[6], start, end))
        if is_exon_feature:
            if not parents:
                raise ValueError(f"{path}:{number}: {fields[2]} requires a transcript Parent")
            expected.update(parents)
            child_records.append((number, fields[2], fields[0], fields[6], start, end, parents))
            if feature == "exon":
                exon_models.update(parents)
            else:
                cds_utr_models.update(parents)
        features.append((tid, parents))
    expected.update(tid for tid, parents in features if tid and parents & gene_ids)
    genes_with_models = gene_ids & expected
    for tid, parents in features:
        if tid in expected:
            genes_with_models.update(parents & gene_ids)
    missing_genes = gene_ids - genes_with_models
    unsupported = expected - exon_models - cds_utr_models
    report.update(source_gene_count=len(gene_ids), source_transcript_model_count=len(expected),
                  source_models_with_explicit_exons=len(exon_models),
                  source_models_without_explicit_exons=len(expected - exon_models))
    # Validate source coordinates before normalization can repair a malformed
    # transcript span. GFF3 permits child records before their parent record.
    for number, feature, chrom, strand, start, end, parents in child_records:
        for parent in parents:
            for declared_chrom, declared_strand, left, right in declared_transcripts.get(parent, []):
                if chrom != declared_chrom or strand != declared_strand:
                    raise ValueError(f"{path}:{number}: {feature} chromosome/strand disagrees with declared transcript {parent!r}")
                if start < left or end > right:
                    raise ValueError(f"{path}:{number}: {feature} interval extends outside declared transcript {parent!r}")
    if missing_genes:
        raise ValueError(f"{path}: {len(missing_genes)} gene(s) lack transcript/exon/CDS/UTR models "
                         f"({', '.join(sorted(missing_genes)[:5])}); provide complete exon annotation")
    if unsupported:
        raise ValueError(f"{path}: {len(unsupported)} transcript(s) lack exon/CDS/UTR records "
                         f"({', '.join(sorted(unsupported)[:5])}); transcript spans alone cannot define exons")
    return expected


def prepare_annotation(annotation, output, report, external_bed=None):
    source = Path(annotation)
    report["annotation_format"] = annotation_format(source)
    normalized = source
    expected_ids = None
    if report["annotation_format"] == "gff3":
        expected_ids = gff3_transcripts(source, report)
        normalized = Path(str(output) + ".normalized.gtf")
        conversion_log = Path(str(output) + ".gffread.log")
        normalized.parent.mkdir(parents=True, exist_ok=True)
        command = ["gffread", str(source), "-T", "-o", str(normalized)]
        report.update(normalized_annotation=str(normalized), conversion_log=str(conversion_log), conversion_command=command)
        with conversion_log.open("w", encoding="utf-8") as log:
            log.write("GFF3 to GTF normalization; gffread may construct exon models from GFF3 CDS/UTR.\n")
            log.flush()
            try:
                result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
            except FileNotFoundError as exc:
                raise ValueError("GFF3 conversion requires gffread in PATH; use the pipeline environment or provide exon-containing GTF") from exc
        report["conversion_exit_code"] = result.returncode
        if result.returncode:
            raise ValueError(f"gffread conversion failed with exit code {result.returncode}; see {conversion_log}")
        report["exon_model_source"] = "gffread-normalized GFF3 (may include CDS/UTR-derived exons)"
    else:
        report["exon_model_source"] = "explicit GTF exon records"
    records, counts = read_gtf_models(normalized)
    report.update(counts)
    if expected_ids is not None:
        missing = expected_ids - {row["name"] for row in records}
        if missing:
            raise ValueError(f"gffread omitted {len(missing)} source transcript model(s): {', '.join(sorted(missing)[:5])}; provide complete exon annotation")
    if external_bed:
        provided = read_bed12(external_bed)
        signature = lambda row: (row["chrom"], row["start"], row["end"], row["strand"], tuple(row["blocks"]))
        expected = {row["name"]: signature(row) for row in records}
        actual = {row["name"]: signature(row) for row in provided}
        missing = set(expected) - set(actual)
        extra = set(actual) - set(expected)
        changed = {name for name in set(expected) & set(actual) if expected[name] != actual[name]}
        if missing or extra or changed:
            raise ValueError(f"external BED12 does not match annotation: {len(missing)} missing, {len(extra)} extra, {len(changed)} inconsistent transcript model(s)")
        report["external_bed12_validation"] = "exact transcript ID, chromosome, strand, interval, and exon-block match"
    write_bed12(records, output)
    report.update(status="passed", transcript_count=len(records), reference_sequence_count=len({row["chrom"] for row in records}))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation", required=True, help="The same GTF/GFF3 used by StringTie")
    parser.add_argument("--output", required=True, help="Validated BED12 output")
    parser.add_argument("--report", required=True, help="JSON validation report, including failures")
    parser.add_argument("--bed12", help="Optional external BED12; must exactly match annotation models")
    args = parser.parse_args(argv)
    # Never overwrite source annotations, or confuse report and result paths.
    sources = {Path(args.annotation).resolve()}
    if args.bed12:
        sources.add(Path(args.bed12).resolve())
    outputs = [Path(path).resolve() for path in (args.output, args.report, args.output + ".normalized.gtf", args.output + ".gffread.log")]
    if len(set(outputs)) != len(outputs) or sources.intersection(outputs):
        parser.error("annotation, external BED12, output, report, and normalization paths must be distinct")
    report = {"status": "failed", "annotation": str(Path(args.annotation)), "bed12": str(Path(args.output)),
              "external_bed12": args.bed12}
    try:
        prepare_annotation(args.annotation, args.output, report, args.bed12)
    except (ValueError, OSError) as exc:
        report["error"] = str(exc)
        print(f"[ERROR] {exc}", file=sys.stderr)
    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
