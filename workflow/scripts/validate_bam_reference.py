#!/usr/bin/env python3
"""Validate BAM/BAI and annotation coordinates before either quantification mode."""
import argparse
from pathlib import Path
import sys

from prepare_strand_bed import protect_paths, read_bed12, write_json


def validate_reference(bam_path, bai_path, bed_path):
    for path, label in ((bam_path, "BAM"), (bai_path, "BAI"), (bed_path, "BED12")):
        path = Path(path)
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"{label} must be an existing, nonempty file: {path}")
        with path.open("rb") as handle:
            handle.read(1)
    import pysam

    records = read_bed12(bed_path)
    # quickcheck verifies the header and terminal EOF block without rescanning
    # the entire BAM; the coordinate index is produced by the upstream rule.
    pysam.quickcheck(str(bam_path))
    with pysam.AlignmentFile(str(bam_path), "rb", index_filename=str(bai_path), require_index=True) as bam:
        if not bam.is_bam:
            raise ValueError("The input is not a BAM file")
        if not bam.check_index():
            raise ValueError("BAM index is not usable")
        if bam.header.to_dict().get("HD", {}).get("SO") != "coordinate":
            raise ValueError("BAM header must declare coordinate sorting (HD:SO:coordinate)")
        lengths = dict(zip(bam.references, bam.lengths))
        if not lengths or len(lengths) != len(bam.references) or any(length <= 0 for length in lengths.values()):
            raise ValueError("BAM SQ reference names/lengths are missing, duplicated or invalid")
        index_statistics = bam.get_index_statistics()
        if [row.contig for row in index_statistics] != list(bam.references):
            raise ValueError("BAM index reference list does not match BAM SQ")
        for record in records:
            chrom = record["chrom"]
            if chrom not in lengths:
                raise ValueError(f"BED12 contig {chrom!r} is absent from BAM SQ")
            if record["start"] < 0 or record["end"] > lengths[chrom]:
                raise ValueError(f"BED12 transcript {record['name']!r} exceeds BAM SQ coordinates")
    return {
        "bam_contigs": len(lengths), "bed12_transcripts": len(records),
        "annotation_contigs": len({row["chrom"] for row in records}),
        "sort_order": "coordinate", "bam_quickcheck": "passed", "index_check": "passed",
        "assembly_check": "Sequence names and coordinate bounds only; matching assembly provenance must be established separately",
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("bam", "bai", "bed12", "output"):
        parser.add_argument("--" + name, required=True)
    args = parser.parse_args(argv)
    try:
        protect_paths([args.bam, args.bai, args.bed12], [args.output])
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    report = {
        "schema_version": 1, "status": "FAIL", "errors": [],
        "inputs": {key: str(Path(getattr(args, key)).resolve()) for key in ("bam", "bai", "bed12")},
    }
    Path(args.output).unlink(missing_ok=True)
    try:
        report["validation"] = validate_reference(args.bam, args.bai, args.bed12)
        report["status"] = "PASS"
    except Exception as exc:
        report["errors"].append(f"{type(exc).__name__}: {exc}")
        print(f"[ERROR] BAM/reference validation failed: {exc}", file=sys.stderr)
    write_json(args.output, report)
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
