#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
def main():
    parser = argparse.ArgumentParser(description="Export RNA-seq result CSVs into one Excel workbook.")
    parser.add_argument("--alignment-report", required=True, help="alignment_report.csv")
    parser.add_argument("--count", required=True, help="gene_count_matrix.csv")
    parser.add_argument("--fpkm", required=True, help="gene_fpkm_matrix.csv")
    parser.add_argument("--tpm", required=True, help="gene_tpm_matrix.csv")
    parser.add_argument("--output", required=True, help="Output XLSX file")
    args = parser.parse_args()
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df_align = pd.read_csv(args.alignment_report)
    df_count = pd.read_csv(args.count)
    df_fpkm = pd.read_csv(args.fpkm)
    df_tpm = pd.read_csv(args.tpm)
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        df_align.to_excel(writer, sheet_name="alignment_report", index=False)
        df_count.to_excel(writer, sheet_name="gene_count", index=False)
        df_fpkm.to_excel(writer, sheet_name="gene_fpkm", index=False)
        df_tpm.to_excel(writer, sheet_name="gene_tpm", index=False)
    print(f"[INFO] Excel workbook written to: {output_path}")
if __name__ == "__main__":
    main()
