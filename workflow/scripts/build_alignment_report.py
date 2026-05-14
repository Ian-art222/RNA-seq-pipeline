#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path
def pct(count, total):
    return f"{count / total * 100:.2f}" if total > 0 else "0.00"

def parse_hisat2_summary(summary_file: Path):
    text = summary_file.read_text(encoding="utf-8", errors="ignore")
    total_match = re.search(r"(\d+)\s+reads; of these:", text)
    unique_match = re.search(r"(\d+)\s+\([\d\.]+%\)\s+aligned concordantly exactly 1 time", text)
    multi_match = re.search(r"(\d+)\s+\([\d\.]+%\)\s+aligned concordantly >1 times", text)
    concordant_0_match = re.search(r"(\d+)\s+\([\d\.]+%\)\s+aligned concordantly 0 times", text)
    discordant_1_match = re.search(r"(\d+)\s+\([\d\.]+%\)\s+aligned discordantly 1 time", text)
    overall_match = re.search(r"([\d\.]+)%\s+overall alignment rate", text)
    total = int(total_match.group(1)) if total_match else 0
    unique = int(unique_match.group(1)) if unique_match else 0
    multi = int(multi_match.group(1)) if multi_match else 0
    concordant_0 = int(concordant_0_match.group(1)) if concordant_0_match else 0
    discordant_1 = int(discordant_1_match.group(1)) if discordant_1_match else 0
    overall = overall_match.group(1) if overall_match else "0.00"
    sample_id = summary_file.name.replace(".align_summary.txt", "")
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
        "Overall_Map(%)": overall,
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
    rows = [parse_hisat2_summary(p) for p in summary_paths]
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
