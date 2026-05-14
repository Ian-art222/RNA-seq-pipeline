#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path
def extract_ordered_gene_ids(annotation_path: Path):
    ordered = []
    seen = set()
    with annotation_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            attrs = parts[8]
            gene_id = None
            # GFF/GFF3: 优先 gene feature 的 ID=
            if feature_type == "gene":
                m = re.search(r"(?:^|;)ID=([^;]+)", attrs)
                if m:
                    gene_id = m.group(1)
            # GTF: gene_id "xxx"
            if gene_id is None:
                m = re.search(r'gene_id "([^"]+)"', attrs)
                if m:
                    gene_id = m.group(1)
            # 兜底：某些 GFF 可能在非 gene feature 上带 Parent/ID，但这里不强行解析 Parent
            if gene_id and gene_id not in seen:
                seen.add(gene_id)
                ordered.append(gene_id)
    if not ordered:
        raise ValueError(f"[ERROR] No gene IDs extracted from annotation file: {annotation_path}")
    return ordered
def load_raw_count_matrix(raw_count_path: Path):
    with raw_count_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        if not header:
            raise ValueError(f"[ERROR] Empty raw count matrix: {raw_count_path}")
        sample_names = header[1:]
        data = {}
        for row in reader:
            if not row:
                continue
            gene_id = row[0]
            data[gene_id] = row[1:]
    return sample_names, data
def parse_abundance_file(abundance_path: Path):
    data_fpkm = {}
    data_tpm = {}
    with abundance_path.open("r", encoding="utf-8", errors="ignore") as fh:
        header = next(fh, None)
        if header is None:
            return data_fpkm, data_tpm
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                continue
            gene_id = parts[0]
            fpkm = parts[7]
            tpm = parts[8]
            data_fpkm[gene_id] = fpkm
            data_tpm[gene_id] = tpm
    return data_fpkm, data_tpm
def write_count_matrix(out_path: Path, ordered_genes, sample_names, count_data):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["GeneID"] + sample_names)
        for gene_id in ordered_genes:
            row = count_data.get(gene_id)
            if row is None:
                row = ["0"] * len(sample_names)
            writer.writerow([gene_id] + row)
def write_float_matrix(out_path: Path, ordered_genes, sample_names, matrix_dict):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["GeneID"] + sample_names)
        for gene_id in ordered_genes:
            row = [matrix_dict.get(gene_id, {}).get(sample, "0.000000") for sample in sample_names]
            writer.writerow([gene_id] + row)
def main():
    parser = argparse.ArgumentParser(description="Finalize count/FPKM/TPM matrices from StringTie outputs.")
    parser.add_argument("--annotation", required=True, help="Reference annotation file (GTF/GFF)")
    parser.add_argument("--raw-count", required=True, help="prepDE raw gene count matrix")
    parser.add_argument("--out-count", required=True, help="Final gene count matrix CSV")
    parser.add_argument("--out-fpkm", required=True, help="Final gene FPKM matrix CSV")
    parser.add_argument("--out-tpm", required=True, help="Final gene TPM matrix CSV")
    parser.add_argument("abundance_tabs", nargs="+", help="Per-sample StringTie abundance.tab files")
    args = parser.parse_args()
    annotation_path = Path(args.annotation)
    raw_count_path = Path(args.raw_count)
    out_count = Path(args.out_count)
    out_fpkm = Path(args.out_fpkm)
    out_tpm = Path(args.out_tpm)
    abundance_paths = [Path(x) for x in args.abundance_tabs]
    ordered_genes = extract_ordered_gene_ids(annotation_path)
    sample_names, count_data = load_raw_count_matrix(raw_count_path)
    fpkm_matrix = {}
    tpm_matrix = {}
    for abundance_path in sorted(abundance_paths):
        sample = abundance_path.name.replace("_abundance.tab", "")
        data_fpkm, data_tpm = parse_abundance_file(abundance_path)
        for gene_id, value in data_fpkm.items():
            fpkm_matrix.setdefault(gene_id, {})[sample] = value
        for gene_id, value in data_tpm.items():
            tpm_matrix.setdefault(gene_id, {})[sample] = value
    write_count_matrix(out_count, ordered_genes, sample_names, count_data)
    write_float_matrix(out_fpkm, ordered_genes, sample_names, fpkm_matrix)
    write_float_matrix(out_tpm, ordered_genes, sample_names, tpm_matrix)
    print(f"[INFO] Wrote final count matrix: {out_count}")
    print(f"[INFO] Wrote final FPKM matrix: {out_fpkm}")
    print(f"[INFO] Wrote final TPM matrix: {out_tpm}")
if __name__ == "__main__":
    main()
