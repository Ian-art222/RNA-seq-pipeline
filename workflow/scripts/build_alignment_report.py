#!/usr/bin/env python3
import argparse
import csv
import re
from decimal import Decimal, InvalidOperation
from pathlib import Path


def pct(count, total):
    return f"{count / total * 100:.2f}" if total > 0 else "0.00"


def parse_hisat2_summary(summary_file: Path):
    """Read the standard HISAT2 PE summary, rejecting incomplete or corrupt data.

    Counts in concordant rows describe pairs; the overall rate describes mates.
    HISAT2 omits the discordant/mate sections with --no-discordant/--no-mixed.
    Their absence therefore means zero alignments in that optional category.
    """
    summary_file = Path(summary_file)
    lines = summary_file.read_text(encoding="utf-8").splitlines()

    def fail(message):
        raise ValueError(f"Invalid HISAT2 paired-end summary {summary_file}: {message}")

    def field(pattern, name, optional=False):
        matches = [match for line in lines if (match := re.fullmatch(pattern, line.strip()))]
        if not matches and optional:
            return None
        if len(matches) != 1:
            fail(f"expected exactly one {name} field, found {len(matches)}")
        return matches[0].groups()

    def count(value, name):
        if not re.fullmatch(r"[0-9]+", value):
            fail(f"{name} must be a nonnegative integer: {value!r}")
        return int(value)

    def percentage(value, name, numerator, denominator):
        try:
            result = Decimal(value)
        except InvalidOperation:
            fail(f"invalid {name} percentage: {value!r}")
        if not result.is_finite() or not 0 <= result <= 100:
            fail(f"{name} percentage must be finite and between 0 and 100")
        expected = Decimal(numerator) * 100 / denominator if denominator else Decimal(0)
        # HISAT2 prints two decimal places. Allow rounding of its displayed rate.
        if abs(result - expected) > Decimal("0.0051"):
            fail(f"{name} percentage {value} disagrees with counts {numerator}/{denominator}")
        return result

    def count_percentage(label, name, denominator, optional=False):
        values = field(r"(\S+)\s+\((\S+)%\)\s+" + re.escape(label), name, optional)
        if values is None:
            return None
        number = count(values[0], name)
        if number > denominator:
            fail(f"{name} count exceeds its denominator {denominator}")
        percentage(values[1], name, number, denominator)
        return number

    total = count(field(r"(\S+)\s+reads; of these:", "total reads")[0], "total reads")
    paired = count_percentage("were paired; of these:", "paired reads", total)
    if paired != total:
        fail("this pipeline requires all input reads to be paired-end")
    concordant_0 = count_percentage("aligned concordantly 0 times", "concordant zero", total)
    unique = count_percentage("aligned concordantly exactly 1 time", "concordant unique", total)
    multi = count_percentage("aligned concordantly >1 times", "concordant multiple", total)
    if concordant_0 + unique + multi != total:
        fail("concordant zero/unique/multiple counts do not sum to total pairs")

    discordant_1 = count_percentage(
        "aligned discordantly 1 time", "discordant unique", concordant_0, optional=True
    )
    discordant_header = field(
        r"(\S+)\s+pairs aligned concordantly 0 times; of these:", "discordant pairs", optional=True
    )
    if (discordant_header is None) != (discordant_1 is None):
        fail("discordant alignment section is incomplete")
    if discordant_header is not None:
        if count(discordant_header[0], "discordant pairs") != concordant_0:
            fail("discordant section pair count differs from concordant zero")
    if discordant_1 is None:
        discordant_1 = 0

    remaining_pairs = concordant_0 - discordant_1
    mate_header = field(r"(\S+)\s+mates make up the pairs; of these:", "mate total", optional=True)
    remaining_header = field(
        r"(\S+)\s+pairs aligned 0 times concordantly or discordantly; of these:",
        "remaining pairs", optional=True,
    )
    if (remaining_header is None) != (mate_header is None):
        fail("remaining pairs/mate alignment section is incomplete")
    if remaining_header is not None:
        if count(remaining_header[0], "remaining pairs") != remaining_pairs:
            fail("remaining pair count is inconsistent with concordant/discordant counts")

    mate_total = remaining_pairs * 2
    mate_counts = [count_percentage(label, name, mate_total, optional=True) for label, name in (
        ("aligned 0 times", "unmapped mates"),
        ("aligned exactly 1 time", "unique mates"),
        ("aligned >1 times", "multiple mates"),
    )]
    if mate_header is not None or any(value is not None for value in mate_counts):
        if mate_header is None or any(value is None for value in mate_counts):
            fail("mate alignment section is incomplete")
        if count(mate_header[0], "mate total") != mate_total or sum(mate_counts) != mate_total:
            fail("mate counts do not match the remaining pairs")
        mapped_mates = mate_counts[1] + mate_counts[2]
    else:
        mapped_mates = 0
    overall_text = field(r"(\S+)%\s+overall alignment rate", "overall alignment rate")[0]
    overall = percentage(overall_text, "overall alignment", 2 * (unique + multi + discordant_1) + mapped_mates, 2 * total)
    suffix = ".align_summary.txt"
    sample_id = summary_file.name[:-len(suffix)] if summary_file.name.endswith(suffix) else summary_file.name
    if not sample_id:
        fail("empty sample ID")
    return {
        "SampleID": sample_id,
        "Total_Read_Pairs": str(total),
        "Concordant_Unique_Pairs": str(unique),
        "Concordant_Unique(%)": pct(unique, total),
        "Concordant_Multi_Pairs": str(multi),
        "Concordant_Multi(%)": pct(multi, total),
        "Concordant_0_Pairs": str(concordant_0),
        "Concordant_0(%)": pct(concordant_0, total),
        "Discordant_1_Pairs": str(discordant_1),
        "Discordant_1(%)": pct(discordant_1, total),
        "Overall_Map(%)": f"{overall:.2f}",
    }


def main():
    parser = argparse.ArgumentParser(description="Build alignment_report.csv from HISAT2 summary files.")
    parser.add_argument("--output", required=True, help="Output CSV path")
    parser.add_argument("summary_files", nargs="+", help="HISAT2 summary files")
    args = parser.parse_args()
    summary_paths = [Path(x) for x in args.summary_files]
    summary_paths = sorted(summary_paths)
    if not summary_paths:
        raise SystemExit("[ERROR] No HISAT2 summary files provided.")
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        rows = [parse_hisat2_summary(p) for p in summary_paths]
        if len({row["SampleID"] for row in rows}) != len(rows):
            raise ValueError("Duplicate sample IDs in HISAT2 summary inputs")
    except (OSError, UnicodeError, ValueError) as error:
        raise SystemExit(f"[ERROR] {error}") from error
    fieldnames = [
        "SampleID",
        "Total_Read_Pairs",
        "Concordant_Unique_Pairs",
        "Concordant_Unique(%)",
        "Concordant_Multi_Pairs",
        "Concordant_Multi(%)",
        "Concordant_0_Pairs",
        "Concordant_0(%)",
        "Discordant_1_Pairs",
        "Discordant_1(%)",
        "Overall_Map(%)",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"[INFO] Alignment report written to: {output_path}")
if __name__ == "__main__":
    main()
