#!/usr/bin/env python3
"""StringTie coverage-derived read estimates, with explicit sample/read-length checks.

Formula and legacy CLI are based on gpertea/stringtie prepDE.py3:
https://github.com/gpertea/stringtie/blob/master/prepDE.py3
Counts are ceil(transcript coverage * summed exon length / mean single-read length).
They are estimates of reads, not directly counted paired-end fragments.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile
from collections import defaultdict

from prepare_strand_bed import read_gtf_models

ATTR = re.compile(r'(?:^|;)\s*([^\s;]+)\s+"([^"]*)"\s*(?=;|$)')


def attributes(text):
    result = {}
    for key, value in ATTR.findall(text):
        if key in result and key in {'gene_id', 'transcript_id', 'cov', 'gene_name', 'ref_gene_name', 'FPKM', 'TPM'}:
            raise ValueError(f"Duplicate GTF attribute: {key}")
        # Repeated descriptive metadata such as GENCODE tag is valid GTF.
        result.setdefault(key, value)
    return result


def positive_length(value, label):
    if isinstance(value, bool):
        raise ValueError(f"{label} must be a positive finite read length")
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{label} must be a positive finite read length")
    return value


def sample_list(path, pattern):
    path = Path(path)
    samples = []
    if path.is_file():
        with path.open(encoding="utf-8") as handle:
            for number, line in enumerate(handle, 1):
                if not line.strip() or line.lstrip().startswith('#'):
                    continue
                fields = line.strip().split(None, 1)
                if len(fields) != 2:
                    raise ValueError(f"{path}:{number}: expected sample ID and GTF path")
                samples.append((fields[0], Path(fields[1])))
    elif path.is_dir():
        for directory in sorted(path.iterdir()):
            if directory.is_dir() and re.search(pattern, directory.name):
                samples.extend((directory.name, p) for p in sorted(directory.glob('*.gtf')))
    else:
        raise ValueError(f"Input does not exist: {path}")
    if not samples:
        raise ValueError(f"No sample GTF files found: {path}")
    seen = set()
    for sample, gtf in samples:
        if sample in seen:
            raise ValueError(f"Duplicate sample ID: {sample}")
        if not sample or any(char.isspace() for char in sample):
            raise ValueError(f"Invalid sample ID: {sample!r}")
        seen.add(sample)
        if not gtf.is_file():
            raise ValueError(f"GTF file not found: {gtf}")
    return sorted(samples)


def unique_json_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"Duplicate read-length sample: {key}")
        result[key] = value
    return result


def read_lengths(args, samples):
    default = positive_length(args.length, '-l/--length')
    if not args.read_lengths:
        return {sample: default for sample, _ in samples}
    with open(args.read_lengths, encoding='utf-8') as handle:
        values = json.load(handle, object_pairs_hook=unique_json_object)
    expected = {sample for sample, _ in samples}
    if not isinstance(values, dict) or set(values) != expected:
        raise ValueError('Read-length JSON sample IDs must exactly match the sample list')
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) for value in values.values()):
        raise ValueError('Read-length JSON values must be numbers, not strings or booleans')
    return {sample: positive_length(value, f'read length for {sample}') for sample, value in values.items()}


def parse_reference_gtf(path, gene_id_only=False):
    """Use the validated exon model parser shared with reference preparation.

    Reference GTF may be exon-only or contain CDS/UTR/codon metadata. Coverage
    and StringTie command provenance belong exclusively to quantified output.
    """
    records, _ = read_gtf_models(path)
    models = {}
    for row in records:
        gid, symbol = row['gene_id'], row['gene_name']
        exons = tuple((start + 1, end) for start, end in row['blocks'])
        models[row['name']] = dict(
            gene=gid if gene_id_only or not symbol else gid + '|' + symbol,
            canonical_gene=gid, chrom=row['chrom'], strand=row['strand'],
            start=row['start'] + 1, end=row['end'], coverage=0.0,
            exons=exons, length=sum(end - start + 1 for start, end in exons))
    return models


def parse_gtf(path, gene_id_only=False, reference=False):
    if reference:
        return parse_reference_gtf(path, gene_id_only)
    path = Path(path)
    transcripts, exons = {}, defaultdict(list)
    command_seen = False
    with path.open(encoding='utf-8') as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            if line.startswith('#'):
                if re.search(r'(?:^|\s)-e(?:\s|$)', line) and re.search(r'(?:^|\s)-G\S*(?:\s|$)', line):
                    command_seen = True
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                raise ValueError(f'{path}:{number}: expected nine GTF fields')
            if fields[2] not in ('transcript', 'exon'):
                raise ValueError(f'{path}:{number}: unexpected StringTie feature {fields[2]!r}')
            start, end = int(fields[3]), int(fields[4])
            if start < 1 or end < start or fields[6] not in ('+', '-'):
                raise ValueError(f'{path}:{number}: invalid coordinates or strand')
            attrs = attributes(fields[8])
            tid, gid = attrs.get('transcript_id'), attrs.get('gene_id')
            if not tid or not gid:
                raise ValueError(f'{path}:{number}: transcript_id and gene_id are required')
            if fields[2] == 'transcript':
                if tid in transcripts:
                    raise ValueError(f'{path}:{number}: duplicate transcript {tid}')
                if 'cov' not in attrs:
                    raise ValueError(f'{path}:{number}: missing transcript coverage')
                coverage = float(attrs['cov'])
                if not math.isfinite(coverage) or coverage < 0:
                    raise ValueError(f'{path}:{number}: coverage must be finite and nonnegative')
                symbol = attrs.get('gene_name') or attrs.get('ref_gene_name')
                output_gene = gid if gene_id_only or not symbol else gid + '|' + symbol
                transcripts[tid] = dict(gene=output_gene, canonical_gene=gid, chrom=fields[0], strand=fields[6],
                                        start=start, end=end, coverage=coverage)
            else:
                exons[tid].append((fields[0], fields[6], start, end, gid))
    if not command_seen:
        raise ValueError(f'{path}: missing StringTie -e -G command provenance')
    if not transcripts:
        raise ValueError(f'{path}: no quantified transcript models')
    if set(exons) != set(transcripts):
        raise ValueError(f'{path}: transcript/exon IDs differ; every transcript requires exons')
    for tid, transcript in transcripts.items():
        ordered = sorted(exons[tid], key=lambda exon: exon[2])
        previous_end = 0
        for chrom, strand, start, end, gid in ordered:
            if (chrom, strand, gid) != (transcript['chrom'], transcript['strand'], transcript['canonical_gene']):
                raise ValueError(f'{path}: inconsistent exon metadata for {tid}')
            if start < transcript['start'] or end > transcript['end'] or start <= previous_end:
                raise ValueError(f'{path}: overlapping or out-of-bounds exons for {tid}')
            previous_end = end
        transcript['exons'] = tuple((x[2], x[3]) for x in ordered)
        transcript['length'] = sum(end - start + 1 for start, end in transcript['exons'])
    return transcripts


def match_reference(models, reference, label):
    extra = set(models) - set(reference)
    if extra:
        raise ValueError(f'{label}: unexpected reference transcript IDs {sorted(extra)[:10]}')
    fields = ('gene', 'canonical_gene', 'chrom', 'strand', 'start', 'end', 'exons')
    for tid, model in models.items():
        if any(model[field] != reference[tid][field] for field in fields):
            raise ValueError(f'{label}: inconsistent reference model for transcript {tid}')


def complete_from_ballgown(models, reference, path, sample):
    """Recover GTF omissions only from explicit zero estimates in StringTie -B.

    StringTie 2.2.3 can omit a redundant, zero-estimated isoform from its GTF
    inside an expressed bundle, while retaining its zero in t_data.ctab.
    This is an estimator-assigned zero, not evidence of biological absence.
    """
    expected_header = ['t_id', 'chr', 'strand', 'start', 'end', 't_name',
                       'num_exons', 'length', 'gene_id', 'gene_name', 'cov', 'FPKM']
    records, row_ids = {}, set()
    with path.open(encoding='utf-8', newline='') as handle:
        reader = csv.reader(handle, delimiter='\t')
        if next(reader, None) != expected_header:
            raise ValueError(f'{path}: invalid Ballgown transcript header')
        for number, row in enumerate(reader, 2):
            if len(row) != len(expected_header):
                raise ValueError(f'{path}:{number}: invalid Ballgown transcript row width')
            rid, tid = int(row[0]), row[5]
            if rid < 1 or rid in row_ids or not tid or tid in records or tid not in reference:
                raise ValueError(f'{path}:{number}: missing, duplicate or unexpected transcript ID')
            row_ids.add(rid)
            model = reference[tid]
            observed = (row[8], row[1], row[2], int(row[3]), int(row[4]), int(row[6]), int(row[7]))
            expected = (model['canonical_gene'], model['chrom'], model['strand'], model['start'],
                        model['end'], len(model['exons']), model['length'])
            if observed != expected:
                raise ValueError(f'{path}:{number}: Ballgown/reference model mismatch for {tid}')
            coverage, fpkm = float(row[10]), float(row[11])
            if any(not math.isfinite(value) or value < 0 for value in (coverage, fpkm)):
                raise ValueError(f'{path}:{number}: non-finite or negative Ballgown estimate')
            if tid in models and coverage != models[tid]['coverage']:
                raise ValueError(f'{path}:{number}: Ballgown/GTF coverage mismatch for {tid}')
            records[tid] = (coverage, fpkm)
    if set(records) != set(reference):
        raise ValueError(f'{path}: Ballgown transcript IDs do not exactly match the reference')
    missing = sorted(set(reference) - set(models))
    for tid in missing:
        if records[tid] != (0.0, 0.0):
            raise ValueError(f'{sample}: GTF omitted nonzero Ballgown transcript {tid}; cannot restore zero')
        models[tid] = dict(reference[tid], coverage=0.0)
    if missing:
        print(f'[INFO] {sample}: restored {len(missing)} reference transcript(s) omitted from GTF '
              f'using explicit zero cov/FPKM in {path}: {", ".join(missing)}', file=sys.stderr)
    return models


def atomic_csv(path, header, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name + '.', dir=path.parent)
    try:
        with os.fdopen(fd, 'w', newline='', encoding='utf-8') as handle:
            writer = csv.writer(handle)
            if header is not None:
                writer.writerow(header)
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def clustered_genes(models, prefix, key):
    """Keep legacy -c: connected same-strand exon overlaps form clusters."""
    mapping = {tid: model['gene'] for tid, model in models.items()}
    parent = {tid: tid for tid in models if not models[tid]['gene'].startswith(prefix)}
    def find(tid):
        while parent[tid] != tid:
            parent[tid] = parent[parent[tid]]
            tid = parent[tid]
        return tid
    loci = defaultdict(list)
    for tid in parent:
        model = models[tid]
        loci[(model['chrom'], model['strand'])].append(tid)
    for tids in loci.values():
        tids.sort(key=lambda tid: models[tid]['start'])
        for index, first in enumerate(tids):
            for second in tids[index + 1:]:
                if models[second]['start'] > models[first]['end']:
                    break
                if any(a <= d and c <= b for a, b in models[first]['exons'] for c, d in models[second]['exons']):
                    parent[find(second)] = find(first)
    groups = defaultdict(list)
    for tid in parent:
        groups[find(tid)].append(tid)
    clusters = sorted(sorted(group) for group in groups.values() if len(group) > 1)
    legend = []
    for index, group in enumerate(clusters, 1):
        gid = key + str(index)
        legend.append([gid] + group)
        for tid in group:
            mapping[tid] = gid
    return mapping, legend


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', '--in', default='.')
    parser.add_argument('-g', default='gene_count_matrix.csv')
    parser.add_argument('-t', default='transcript_count_matrix.csv')
    parser.add_argument('-l', '--length', type=float, default=75.0, help='Mean single-read length (default: 75)')
    parser.add_argument('--read-lengths', help='JSON mapping every sample ID to its positive finite mean read length')
    parser.add_argument('--reference', type=Path, help='Normalized reference GTF; require all its transcript models')
    parser.add_argument('--ballgown', action='store_true', help='Validate each GTF against adjacent t_data.ctab; restore only explicitly zero reference models (requires --reference)')
    parser.add_argument('--gene-id-only', action='store_true', help='Use exact annotation gene_id; never append or split gene symbols')
    parser.add_argument('-p', '--pattern', default='.')
    parser.add_argument('-c', '--cluster', action='store_true')
    parser.add_argument('-s', '--string', default='MSTRG')
    parser.add_argument('-k', '--key', default='prepG')
    parser.add_argument('--legend', default='legend.csv')
    parser.add_argument('-v', action='store_true')
    args = parser.parse_args()
    if args.ballgown and not args.reference:
        parser.error('--ballgown requires --reference')
    samples = sample_list(args.input, args.pattern)
    lengths = read_lengths(args, samples)
    reference = parse_gtf(args.reference, args.gene_id_only, reference=True) if args.reference else None

    def load_models(sample, path):
        models = parse_gtf(path, args.gene_id_only)
        if reference is not None:
            match_reference(models, reference, sample)
            if args.ballgown:
                models = complete_from_ballgown(models, reference, path.with_name('t_data.ctab'), sample)
            missing = set(reference) - set(models)
            if missing:
                raise ValueError(f'{sample}: GTF is missing reference transcripts {sorted(missing)[:10]}; '
                                 'enable StringTie Ballgown output and --ballgown for explicit zero verification, '
                                 'or inspect the annotation and StringTie output')
            models = {tid: models[tid] for tid in reference}
        return models

    first = load_models(*samples[0])
    gene_ids = {tid: model['gene'] for tid, model in first.items()}
    if args.cluster:
        gene_ids, legend = clustered_genes(first, args.string, args.key)
        atomic_csv(args.legend, None, legend)
    # Retain one reference model set and the numeric matrices, not a complete
    # copy of every sample's exon/attribute objects in memory.
    names = [sample for sample, _ in samples]
    transcript_counts = {tid: [] for tid in first}
    gene_counts = {gid: [0] * len(samples) for gid in dict.fromkeys(gene_ids.values())}
    for sample_index, (sample, path) in enumerate(samples):
        sample_models = first if sample_index == 0 else load_models(sample, path)
        if set(sample_models) != set(first):
            raise ValueError(f'{sample}: transcript IDs do not match other samples')
        for tid, model in sample_models.items():
            fields = ('gene', 'canonical_gene', 'chrom', 'strand', 'exons')
            if any(model[field] != first[tid][field] for field in fields):
                raise ValueError(f'{sample}: inconsistent reference model for transcript {tid}')
            estimate = model['coverage'] * model['length'] / lengths[sample]
            if not math.isfinite(estimate):
                raise ValueError(f'{sample}/{tid}: non-finite count estimate')
            value = math.ceil(estimate)
            transcript_counts[tid].append(value)
            gene_counts[gene_ids[tid]][sample_index] += value
    atomic_csv(args.t, ['transcript_id'] + names, ([tid] + values for tid, values in transcript_counts.items()))
    atomic_csv(args.g, ['gene_id'] + names, ([gid] + values for gid, values in gene_counts.items()))
    if args.v:
        print(f'Wrote {len(gene_counts)} genes and {len(transcript_counts)} transcripts for {len(samples)} samples')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as error:
        raise SystemExit(f'[ERROR] {error}')
