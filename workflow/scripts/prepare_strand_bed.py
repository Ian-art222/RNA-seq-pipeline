#!/usr/bin/env python3
"""Prepare checked transcript BED12 models for RSeQC from the quantification annotation."""
import argparse
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path


GTF_ATTRIBUTE = re.compile(r'\s*([^\s;]+)\s+(?:"([^"\\]*(?:\\.[^"\\]*)*)"|([^;\s]+))\s*(?:;|$)')
MODEL_FEATURES = {"exon", "cds", "utr", "five_prime_utr", "three_prime_utr", "transcript", "mrna",
                  "start_codon", "stop_codon"}


def protect_paths(inputs, outputs):
    """Reject input/output aliases, including symlinks and hard links, before writing."""
    inputs, outputs = [Path(path) for path in inputs], [Path(path) for path in outputs]

    def same_path(first, second):
        return first.resolve() == second.resolve() or (
            first.exists() and second.exists() and first.samefile(second)
        )

    for index, output in enumerate(outputs):
        for source in inputs:
            if same_path(output, source):
                raise ValueError(f"Output {output} aliases protected input {source}")
        for previous in outputs[:index]:
            if same_path(output, previous):
                raise ValueError(f"Outputs {output} and {previous} alias each other")


def write_json(path, report):
    """Atomically replace a report without following an existing output hard link."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, prefix=path.name + ".", delete=False) as handle:
            temporary = Path(handle.name)
            json.dump(report, handle, ensure_ascii=False, indent=2, allow_nan=False)
            handle.write("\n")
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


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


def check_model_attribute(value, label, location):
    # StringTie's GTF loader cannot preserve quoted/semicolon-containing IDs
    # or names, including backslash-escaped quotes. Reject before conversion.
    if value is not None and any(char in '\";' or ord(char) < 32 or ord(char) == 127 for char in value):
        raise ValueError(f"{location}: {label} contains a character StringTie cannot preserve")


def read_gtf_models(path):
    """Read explicit exon models using StringTie 2.2.3's close-exon policy.

    Never infer missing exon models from CDS or transcript spans.
    """
    models, declared_genes, gene_names, gene_loci = {}, {}, {}, {}
    short_introns = []
    input_exons = 0
    for number, fields, start, end in annotation_rows(path):
        feature = fields[2].lower()
        location = f"{path}:{number}"
        attributes = parse_gtf_attributes(fields[8], location)
        gene = attributes.get("gene_id")
        if feature in {"gene", "pseudogene"} or feature.endswith("_gene"):
            if not gene:
                raise ValueError(f"{location}: gene requires a nonempty gene_id")
            declared_genes.setdefault(gene, []).append((fields[0], fields[6], start, end))
            if attributes.get("gene_name"):
                gene_names.setdefault(gene, attributes["gene_name"])
        tid = attributes.get("transcript_id")
        if not tid:
            if feature in MODEL_FEATURES or feature.endswith("_utr"):
                raise ValueError(f"{location}: {fields[2]} requires transcript_id and gene_id")
            continue
        if not gene or any(char.isspace() for char in tid) or not tid:
            raise ValueError(f"{location}: transcript requires a nonempty gene_id and a whitespace-free transcript_id")
        check_model_attribute(tid, "transcript_id", location)
        check_model_attribute(gene, "gene_id", location)
        if fields[6] not in {"+", "-"}:
            raise ValueError(f"{location}: transcript {tid!r} must have strand '+' or '-'")
        identity = (fields[0], fields[6], gene)
        locus = (fields[0], fields[6])
        if gene in gene_loci and gene_loci[gene] != locus:
            raise ValueError(f"{location}: gene_id {gene!r} has inconsistent chromosome or strand across transcripts")
        gene_loci[gene] = locus
        if tid not in models:
            models[tid] = {"identity": identity, "exons": [], "spans": [], "subfeatures": [],
                           "gene_name": None, "name_priority": -1}
        model = models[tid]
        if model["identity"] != identity:
            raise ValueError(f"{location}: transcript_id {tid!r} has inconsistent chromosome, strand, or gene_id")
        symbol = attributes.get("gene_name") or attributes.get("ref_gene_name")
        priority = 2 if feature in {"transcript", "mrna"} else 1
        if symbol and priority > model["name_priority"]:
            model.update(gene_name=symbol, name_priority=priority)
        if feature == "exon":
            model["exons"].append((start, end))
            input_exons += 1
        elif feature in {"transcript", "mrna"}:
            model["spans"].append((start, end))
        elif feature in {"cds", "utr", "start_codon", "stop_codon"} or feature.endswith("_utr"):
            model["subfeatures"].append((number, fields[2], start, end))
    if not models:
        raise ValueError(f"{path}: no transcript exon models found; provide a GTF with exon records")
    missing = sorted(tid for tid, model in models.items() if not model["exons"])
    if missing:
        examples = ", ".join(missing[:5])
        raise ValueError(f"{path}: {len(missing)} transcript(s) lack exon records ({examples}); supply explicit exon models")
    missing_genes = set(declared_genes) - {model["identity"][2] for model in models.values()}
    if missing_genes:
        raise ValueError(f"{path}: {len(missing_genes)} gene(s) lack transcript/exon models "
                         f"({', '.join(sorted(missing_genes)[:5])}); provide complete exon annotation")
    records = []
    for tid, model in models.items():
        blocks = []
        for start, end in sorted(model["exons"]):
            if blocks and start <= blocks[-1][1]:
                blocks[-1] = (blocks[-1][0], max(end, blocks[-1][1]))
            else:
                blocks.append((start, end))
        start, end = blocks[0][0], blocks[-1][1]
        chrom, strand, gene = model["identity"]
        for gene_chrom, gene_strand, left, right in declared_genes.get(gene, []):
            if chrom != gene_chrom or (gene_strand in {"+", "-"} and strand != gene_strand):
                raise ValueError(f"{path}: transcript {tid!r} chromosome/strand disagrees with declared gene {gene!r}")
            if start < left or end > right or any(a < left or b > right for a, b in model["spans"]):
                raise ValueError(f"{path}: transcript {tid!r} extends outside declared gene {gene!r}")
        if any(start < left or end > right for left, right in model["spans"]):
            raise ValueError(f"{path}: exon blocks extend outside declared transcript {tid!r}")
        for number, feature, left, right in model["subfeatures"]:
            if not any(exon_start <= left and right <= exon_end for exon_start, exon_end in blocks):
                raise ValueError(f"{path}:{number}: {feature} interval is not contained in exon blocks of transcript {tid!r}")
        # StringTie 2.2.3 joins exon blocks separated by fewer than five bases.
        # Apply this same loading behavior to BED12 and prepDE, after validating
        # metadata against the original blocks. Record every added base.
        compatible = []
        for left, right in blocks:
            gap = left - compatible[-1][1] if compatible else 0
            if compatible and 0 < gap < 5:
                short_introns.append({"transcript_id": tid, "chrom": chrom, "strand": strand,
                                      "start_1based": compatible[-1][1] + 1,
                                      "end_1based": left, "added_exonic_bases": gap})
                compatible[-1] = (compatible[-1][0], right)
            else:
                compatible.append((left, right))
        symbol = model["gene_name"] or gene_names.get(gene)
        check_model_attribute(symbol, "gene_name", f"{path}: transcript {tid!r}")
        records.append({"chrom": chrom, "start": start, "end": end, "name": tid, "strand": strand,
                        "blocks": compatible, "gene_id": gene,
                        "gene_name": symbol})
    records.sort(key=lambda row: (row["chrom"], row["start"], row["end"], row["name"]))
    return records, {"input_exon_records": input_exons,
                     "merged_exon_blocks": sum(len(row["blocks"]) for row in records),
                     "short_intron_merge_policy": "StringTie 2.2.3: merge exon gaps of 1-4 bp",
                     "short_intron_merges": short_introns,
                     "short_intron_added_exonic_bases": sum(row["added_exonic_bases"] for row in short_introns)}


def write_model_gtf(records, path):
    """Publish the same explicit exon models used by BED12 and prepDE.

    The source annotation remains untouched. CDS/UTR/codon rows have already
    been checked; quantification needs only transcript/exon rows. IDs and
    optional gene names retain the source GTF attribute escaping.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent,
                                         prefix=path.name + ".", delete=False) as handle:
            temporary = Path(handle.name)
            handle.write("# Validated reference exon models for BED12, StringTie and prepDE\n")
            for row in records:
                attrs = f'gene_id "{row["gene_id"]}"; transcript_id "{row["name"]}";'
                if row["gene_name"]:
                    attrs += f' gene_name "{row["gene_name"]}";'
                for feature, start, end in [("transcript", row["start"], row["end"])] + [
                    ("exon", start, end) for start, end in row["blocks"]
                ]:
                    handle.write(f'{row["chrom"]}\treference\t{feature}\t{start + 1}\t{end}\t.\t'
                                 f'{row["strand"]}\t.\t{attrs}\n')
        # Validate serialization with the shared parser before publishing it.
        checked, _ = read_gtf_models(temporary)
        if checked != records:
            raise ValueError("Normalized GTF changed the validated reference exon models")
        temporary.replace(path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


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
    declared_transcripts, declared_genes, child_records, exon_intervals = {}, {}, [], {}
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
        # gffread preserves percent escapes in identifiers. Keep the same spelling
        # throughout validation, normalized GTF, BED12 and downstream matrices.
        tid = attributes.get("ID", "")
        parents = {value for value in attributes.get("Parent", "").split(",") if value}
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
            declared_genes.setdefault(tid, []).append((fields[0], fields[6], start, end))
        if is_transcript:
            # GTF/StringTie assigns each transcript to one gene. gffread 0.12.9
            # silently keeps only the first gene when Parent has multiple IDs.
            if len(parents) > 1:
                raise ValueError(f"{path}:{number}: transcript {tid!r} has multiple Parent IDs; "
                                 "GTF quantification requires one gene per transcript")
            expected.add(tid)
            declared_transcripts.setdefault(tid, []).append((fields[0], fields[6], start, end))
        if is_exon_feature:
            if not parents:
                raise ValueError(f"{path}:{number}: {fields[2]} requires a transcript Parent")
            expected.update(parents)
            child_records.append((number, fields[2], fields[0], fields[6], start, end, parents))
            if feature == "exon":
                exon_models.update(parents)
                for parent in parents:
                    exon_intervals.setdefault((parent, fields[0], fields[6]), []).append((start, end))
            else:
                cds_utr_models.update(parents)
        features.append((tid, parents))
    # A gene may also parent regulatory features (e.g. TF_binding_site). Only
    # declared transcripts and parents of exon/CDS/UTR records are model IDs.
    genes_with_models = gene_ids & expected
    for tid, parents in features:
        if tid in expected:
            genes_with_models.update(parents & gene_ids)
            # Check declared genes before gffread discards their feature rows.
            for parent in parents & gene_ids:
                for chrom, strand, start, end in declared_transcripts.get(tid, []):
                    for gene_chrom, gene_strand, left, right in declared_genes[parent]:
                        if chrom != gene_chrom or (gene_strand in {"+", "-"} and strand != gene_strand):
                            raise ValueError(f"{path}: transcript {tid!r} chromosome/strand disagrees with declared gene {parent!r}")
                        if start < left or end > right:
                            raise ValueError(f"{path}: transcript {tid!r} extends outside declared gene {parent!r}")
    missing_genes = gene_ids - genes_with_models
    unsupported = expected - exon_models - cds_utr_models
    report.update(source_gene_count=len(gene_ids), source_transcript_model_count=len(expected),
                  source_models_with_explicit_exons=len(exon_models),
                  source_models_without_explicit_exons=len(expected - exon_models),
                  source_models_without_explicit_exon_ids=sorted(expected - exon_models))
    # Validate source coordinates before normalization can repair a malformed
    # transcript span. GFF3 permits child records before their parent record.
    for intervals in exon_intervals.values():
        intervals.sort()
    for number, feature, chrom, strand, start, end, parents in child_records:
        for parent in parents:
            for declared_chrom, declared_strand, left, right in declared_transcripts.get(parent, []):
                if chrom != declared_chrom or strand != declared_strand:
                    raise ValueError(f"{path}:{number}: {feature} chromosome/strand disagrees with declared transcript {parent!r}")
                if start < left or end > right:
                    raise ValueError(f"{path}:{number}: {feature} interval extends outside declared transcript {parent!r}")
            if feature.lower() != "exon" and parent in exon_models:
                exons = exon_intervals.get((parent, chrom, strand), [])
                covered = start
                for left, right in exons:
                    if left > covered:
                        break
                    covered = max(covered, right)
                    if covered >= end:
                        break
                if covered < end:
                    raise ValueError(f"{path}:{number}: {feature} interval is not contained in exon blocks of transcript {parent!r}")
    if missing_genes:
        raise ValueError(f"{path}: {len(missing_genes)} gene(s) lack transcript/exon/CDS/UTR models "
                         f"({', '.join(sorted(missing_genes)[:5])}); provide complete exon annotation")
    if unsupported:
        raise ValueError(f"{path}: {len(unsupported)} transcript(s) lack exon/CDS/UTR records "
                         f"({', '.join(sorted(unsupported)[:5])}); transcript spans alone cannot define exons")
    return expected


def prepare_annotation(annotation, output, report, external_bed=None, normalized_gtf=None):
    source = Path(annotation)
    report["annotation_format"] = annotation_format(source)
    normalized = source
    expected_ids = None
    if report["annotation_format"] == "gff3":
        expected_ids = gff3_transcripts(source, report)
        normalized = Path(normalized_gtf) if normalized_gtf else Path(str(output) + ".normalized.gtf")
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
    if counts["short_intron_merges"]:
        print(f'[WARNING] Matched StringTie 2.2.3 reference loading: merged '
              f'{len(counts["short_intron_merges"])} intron(s) of 1-4 bp, adding '
              f'{counts["short_intron_added_exonic_bases"]} exonic bases; '
              'see short_intron_merges in the reference report. Source annotation is unchanged.', file=sys.stderr)
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
    if normalized_gtf or report["annotation_format"] == "gff3":
        destination = Path(normalized_gtf) if normalized_gtf else normalized
        write_model_gtf(records, destination)
        report["normalized_annotation"] = str(destination)
        report["normalized_model_format"] = "explicit transcript/exon GTF from validated exon blocks"
    write_bed12(records, output)
    report.update(status="passed", transcript_count=len(records), reference_sequence_count=len({row["chrom"] for row in records}))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--annotation", required=True, help="The same GTF/GFF3 used by StringTie")
    parser.add_argument("--output", required=True, help="Validated BED12 output")
    parser.add_argument("--report", required=True, help="JSON validation report, including failures")
    parser.add_argument("--bed12", help="Optional external BED12; must exactly match annotation models")
    parser.add_argument("--normalized-gtf", help="Canonical transcript/exon GTF for StringTie and matrices; source IDs are preserved")
    args = parser.parse_args(argv)
    # Never overwrite source annotations, or confuse report and result paths.
    sources = [args.annotation]
    if args.bed12:
        sources.append(args.bed12)
    normalized = args.normalized_gtf or args.output + ".normalized.gtf"
    outputs = [args.output, args.report, normalized, args.output + ".gffread.log"]
    try:
        protect_paths(sources, outputs)
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    report = {"status": "failed", "annotation": str(Path(args.annotation)), "bed12": str(Path(args.output)),
              "external_bed12": args.bed12}
    try:
        # Remove old success artifacts only after proving none aliases an input.
        for path in outputs:
            Path(path).unlink(missing_ok=True)
        prepare_annotation(args.annotation, args.output, report, args.bed12, args.normalized_gtf)
    except (ValueError, OSError) as exc:
        report["status"] = "failed"
        report["error"] = str(exc)
        Path(args.output).unlink(missing_ok=True)
        Path(normalized).unlink(missing_ok=True)
        print(f"[ERROR] {exc}", file=sys.stderr)
    write_json(args.report, report)
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    sys.exit(main())
