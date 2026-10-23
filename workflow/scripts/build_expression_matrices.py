#!/usr/bin/env python3
"""Validate and combine reference-only StringTie estimates without inventing zeros."""
import argparse
import csv
from decimal import Decimal, InvalidOperation
import math
from pathlib import Path

from prepDE import attributes, atomic_csv


def finite_nonnegative(value, label):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f'{label}: expected a number, got {value!r}')
    if not math.isfinite(number) or number < 0:
        raise ValueError(f'{label}: expected a finite nonnegative number, got {value!r}')
    return number


def nonnegative_integer(value, label):
    # Decimal checks integrality before float rounding can turn 1.000...001
    # into an integer or lose large integer counts.
    try:
        number = Decimal(value)
    except (InvalidOperation, TypeError, ValueError):
        raise ValueError(f'{label}: counts must be nonnegative integers')
    if not number.is_finite() or number < 0 or number != number.to_integral_value():
        raise ValueError(f'{label}: counts must be nonnegative integers')
    return int(number)


def checked_ids(values, label):
    if not values or any(not value or value != value.strip() for value in values):
        raise ValueError(f'{label}: IDs must be nonempty, without surrounding whitespace')
    if len(set(values)) != len(values):
        raise ValueError(f'{label}: duplicate IDs')
    return values


def require_same_ids(actual, expected, label):
    missing, extra = set(expected) - set(actual), set(actual) - set(expected)
    if missing or extra:
        raise ValueError(f'{label}: missing IDs {sorted(missing)[:10]}; unexpected IDs {sorted(extra)[:10]}')


def annotation_models(path):
    """Read GTF or GFF3 gene/transcript/exon relationships, including ncRNA_gene.

    A gene is quantifiable only when it has a validated transcript/exon model.
    Missing output for such a model is an error: StringTie -e reports explicit
    zero-valued records for unexpressed models, so absence is not zero evidence.
    """
    records = []
    with Path(path).open(encoding='utf-8') as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'{path}:{number}: expected nine annotation fields')
            start, end = int(fields[3]), int(fields[4])
            if start < 1 or end < start:
                raise ValueError(f'{path}:{number}: invalid coordinates')
            records.append((number, fields, start, end))
    if not records:
        raise ValueError(f'{path}: empty annotation')
    is_gtf = bool(attributes(records[0][1][8]))
    models, gene_order, seen_genes = {}, [], set()

    def add_gene(gid):
        if not gid:
            raise ValueError(f'{path}: missing gene ID')
        if gid not in seen_genes:
            seen_genes.add(gid)
            gene_order.append(gid)

    def add_model(tid, gid, chrom, strand, exon=None):
        if not tid or not gid or strand not in ('+', '-'):
            raise ValueError(f'{path}: invalid transcript ID, gene ID or strand')
        add_gene(gid)
        metadata = (gid, chrom, strand)
        if tid in models and models[tid]['metadata'] != metadata:
            raise ValueError(f'{path}: conflicting transcript model {tid}')
        model = models.setdefault(tid, {'metadata': metadata, 'exons': []})
        if exon:
            model['exons'].append(exon)

    if is_gtf:
        for number, fields, start, end in records:
            attrs = attributes(fields[8])
            feature = fields[2]
            if feature == 'gene':
                add_gene(attrs.get('gene_id'))
            elif feature in ('transcript', 'exon'):
                add_model(attrs.get('transcript_id'), attrs.get('gene_id'), fields[0], fields[6],
                          (start, end) if feature == 'exon' else None)
    else:
        nodes, parsed = {}, []
        for number, fields, start, end in records:
            attrs = {}
            for field in fields[8].split(';'):
                if not field:
                    continue
                key, separator, value = field.partition('=')
                if not separator or key in attrs:
                    raise ValueError(f'{path}:{number}: invalid or duplicate GFF3 attribute')
                attrs[key] = value
            identifier = attrs.get('ID', '')
            # Match StringTie/gffread: preserve encoded ID text (e.g. %3B).
            parents = [p for p in attrs.get('Parent', '').split(',') if p]
            record = dict(id=identifier, parents=parents, feature=fields[2], chrom=fields[0], strand=fields[6], start=start, end=end)
            parsed.append(record)
            if identifier:
                if identifier in nodes:
                    raise ValueError(f'{path}: duplicate GFF3 ID {identifier}')
                nodes[identifier] = record
            if fields[2] == 'gene' or fields[2].endswith('_gene') or fields[2] == 'pseudogene':
                add_gene(identifier)
        def gene_for(identifier, trail=()):
            if identifier in trail:
                raise ValueError(f'{path}: cyclic GFF3 Parent hierarchy')
            node = nodes.get(identifier)
            if node is None:
                raise ValueError(f'{path}: missing GFF3 Parent {identifier}')
            if identifier in seen_genes:
                return identifier
            genes = {gene_for(parent, trail + (identifier,)) for parent in node['parents']}
            if len(genes) != 1:
                raise ValueError(f'{path}: transcript {identifier} must resolve to one gene')
            return next(iter(genes))
        for record in parsed:
            if record['feature'] != 'exon':
                continue
            if not record['parents']:
                raise ValueError(f'{path}: exon without Parent')
            for tid in record['parents']:
                parent = nodes.get(tid)
                if parent is None:
                    raise ValueError(f'{path}: missing exon Parent {tid}')
                gid = gene_for(tid)
                if tid in seen_genes or (record['chrom'], record['strand']) != (parent['chrom'], parent['strand']):
                    raise ValueError(f'{path}: inconsistent GFF3 exon parent {tid}')
                if record['start'] < parent['start'] or record['end'] > parent['end']:
                    raise ValueError(f'{path}: exon outside transcript {tid}')
                add_model(tid, gid, record['chrom'], record['strand'], (record['start'], record['end']))
        # Refuse transcript models with no exon rather than silently losing CDS-only models.
        for record in parsed:
            if record['feature'] in ('mRNA', 'transcript', 'ncRNA', 'lnc_RNA', 'tRNA', 'rRNA') and record['id'] not in models:
                raise ValueError(f"{path}: transcript {record['id']} has no exon model")
    if not models:
        raise ValueError(f'{path}: no quantifiable transcript/exon models')
    gene_loci = {}
    for tid, model in models.items():
        if not model['exons']:
            raise ValueError(f'{path}: transcript {tid} has no exons')
        exons = sorted(model['exons'])
        if any(current[0] <= previous[1] for previous, current in zip(exons, exons[1:])):
            raise ValueError(f'{path}: overlapping or duplicate exons for {tid}')
        gid, chrom, strand = model['metadata']
        locus = gene_loci.setdefault(gid, (chrom, strand))
        if locus != (chrom, strand):
            raise ValueError(f'{path}: gene {gid} spans different contigs or strands')
    require_same_ids(gene_loci, gene_order, f'{path}: genes without a quantifiable exon model')
    return gene_order, gene_loci


def extract_ordered_gene_ids(annotation_path):
    return annotation_models(annotation_path)[0]


def load_raw_count_matrix(raw_count_path):
    with Path(raw_count_path).open(encoding='utf-8', newline='') as handle:
        reader = csv.reader(handle)
        header = next(reader, None)
        if not header or header[0] not in ('gene_id', 'GeneID'):
            raise ValueError(f'{raw_count_path}: missing gene count header')
        sample_names = checked_ids(header[1:], f'{raw_count_path}: samples')
        data = {}
        for number, row in enumerate(reader, 2):
            if len(row) != len(header) or not row[0] or row[0] in data:
                raise ValueError(f'{raw_count_path}:{number}: invalid width, missing or duplicate gene ID')
            values = []
            for sample, value in zip(sample_names, row[1:]):
                values.append(nonnegative_integer(value, f'{raw_count_path}:{number}/{sample}'))
            data[row[0]] = values
        if not data:
            raise ValueError(f'{raw_count_path}: count matrix has no genes')
    return sample_names, data


def parse_abundance_file(abundance_path, expected_loci=None):
    data_fpkm, data_tpm, intervals = {}, {}, {}
    expected_header = ['Gene ID', 'Gene Name', 'Reference', 'Strand', 'Start', 'End', 'Coverage', 'FPKM', 'TPM']
    with Path(abundance_path).open(encoding='utf-8', newline='') as handle:
        reader = csv.reader(handle, delimiter='\t')
        if next(reader, None) != expected_header:
            raise ValueError(f'{abundance_path}: missing or invalid StringTie abundance header')
        for number, parts in enumerate(reader, 2):
            if len(parts) != 9 or not parts[0]:
                raise ValueError(f'{abundance_path}:{number}: malformed abundance record')
            gid, chrom, strand = parts[0], parts[2], parts[3]
            start, end = int(parts[4]), int(parts[5])
            if start < 1 or end < start or strand not in ('+', '-'):
                raise ValueError(f'{abundance_path}:{number}: invalid gene locus')
            if expected_loci is not None and expected_loci.get(gid) != (chrom, strand):
                raise ValueError(f'{abundance_path}:{number}: gene/locus absent from validated annotation: {gid}')
            finite_nonnegative(parts[6], f'{abundance_path}:{number}: Coverage')
            fpkm = finite_nonnegative(parts[7], f'{abundance_path}:{number}: FPKM')
            tpm = finite_nonnegative(parts[8], f'{abundance_path}:{number}: TPM')
            for prior_chrom, prior_strand, prior_start, prior_end in intervals.get(gid, []):
                if (chrom, strand) != (prior_chrom, prior_strand) or start <= prior_end and prior_start <= end:
                    raise ValueError(f'{abundance_path}:{number}: overlapping duplicate or inconsistent locus for gene {gid}')
            intervals.setdefault(gid, []).append((chrom, strand, start, end))
            # StringTie can emit disjoint bundles for one gene. These are distinct
            # isoform contributions; add them instead of retaining the last row.
            data_fpkm[gid] = data_fpkm.get(gid, 0.0) + fpkm
            data_tpm[gid] = data_tpm.get(gid, 0.0) + tpm
    if not data_tpm:
        raise ValueError(f'{abundance_path}: no quantified genes')
    return data_fpkm, data_tpm


def validate_tpm(values, counts, fpkm, label):
    total = math.fsum(values)
    if total == 0:
        if any(counts) or any(fpkm):
            raise ValueError(f'{label}: TPM sum is zero despite nonzero count/FPKM')
    elif not math.isclose(total, 1e6, rel_tol=1e-5, abs_tol=10.0):
        raise ValueError(f'{label}: TPM sum {total:.9g} is inconsistent with 1,000,000')
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--annotation', required=True)
    parser.add_argument('--raw-count', required=True)
    parser.add_argument('--out-count', required=True)
    parser.add_argument('--out-fpkm', required=True)
    parser.add_argument('--out-tpm', required=True)
    parser.add_argument('--samples', nargs='+', help='Expected sample IDs, in output column order')
    parser.add_argument('abundance_tabs', nargs='+')
    args = parser.parse_args()
    ordered_genes, loci = annotation_models(args.annotation)
    sample_names, count_data = load_raw_count_matrix(args.raw_count)
    require_same_ids(count_data, ordered_genes, 'Raw counts versus annotation (use prepDE --gene-id-only)')
    if args.samples:
        checked_ids(args.samples, '--samples')
        require_same_ids(sample_names, args.samples, 'Raw count samples')
        indices = [sample_names.index(sample) for sample in args.samples]
        count_data = {gid: [values[index] for index in indices] for gid, values in count_data.items()}
        sample_names = args.samples
    abundances = {}
    for value in args.abundance_tabs:
        path = Path(value)
        suffix = '_abundance.tab'
        if not path.name.endswith(suffix):
            raise ValueError(f'Unexpected abundance filename: {path.name}')
        sample = path.name[:-len(suffix)]
        if sample in abundances:
            raise ValueError(f'Duplicate abundance sample: {sample}')
        abundances[sample] = parse_abundance_file(path, loci)
    require_same_ids(abundances, sample_names, 'Abundance versus count samples')
    for index, sample in enumerate(sample_names):
        fpkm, tpm = abundances[sample]
        require_same_ids(tpm, ordered_genes, f'{sample}: abundance versus annotation')
        validate_tpm(tpm.values(), (count_data[gid][index] for gid in ordered_genes), fpkm.values(), sample)
    atomic_csv(args.out_count, ['GeneID'] + sample_names, ([gid] + count_data[gid] for gid in ordered_genes))
    for target, metric in ((args.out_fpkm, 0), (args.out_tpm, 1)):
        atomic_csv(target, ['GeneID'] + sample_names,
                   ([gid] + [format(abundances[sample][metric][gid], '.12g') for sample in sample_names] for gid in ordered_genes))
    print(f'[INFO] Validated {len(ordered_genes)} genes and {len(sample_names)} samples; wrote count/FPKM/TPM matrices')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as error:
        raise SystemExit(f'[ERROR] {error}')
