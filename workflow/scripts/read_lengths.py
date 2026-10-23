#!/usr/bin/env python3
"""Validate paired clean FASTQ and measure the exact mean length of individual reads."""
import argparse
import gzip
import json
import math
from pathlib import Path
import tempfile


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile('w', dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        json.dump(data, handle, indent=2, allow_nan=False)
        handle.write('\n')
    try:
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def records(path, mate):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='ascii') as handle:
        number = 0
        while True:
            name = handle.readline()
            if not name:
                return
            number += 1
            seq, plus, quality = (handle.readline().rstrip('\r\n') for _ in range(3))
            if not name.startswith('@') or not plus.startswith('+') or not seq or len(seq) != len(quality):
                raise ValueError(f'{path}: invalid/truncated FASTQ record {number}')
            fields = name[1:].split()
            if not fields:
                raise ValueError(f'{path}: missing read ID at record {number}')
            key = fields[0]
            if key.endswith(('/1', '/2')):
                if key[-1] != str(mate):
                    raise ValueError(f'{path}: incorrect mate suffix at record {number}')
                key = key[:-2]
            if len(fields) > 1 and fields[1].startswith(('1:', '2:')) and fields[1][0] != str(mate):
                raise ValueError(f'{path}: incorrect mate field at record {number}')
            if any(ord(c) < 33 or ord(c) > 126 for c in quality) or any(c.isspace() for c in seq):
                raise ValueError(f'{path}: invalid sequence/quality at record {number}')
            yield key, len(seq)


def measure(sample, r1, r2):
    from itertools import zip_longest
    if Path(r1).resolve() == Path(r2).resolve() or Path(r1).samefile(r2):
        raise ValueError('R1 and R2 must be different files')
    bases, pairs, minimum, maximum = 0, 0, None, 0
    for first, second in zip_longest(records(r1, 1), records(r2, 2)):
        if first is None or second is None:
            raise ValueError('R1 and R2 have different record counts')
        if first[0] != second[0]:
            raise ValueError(f'R1/R2 read IDs disagree at pair {pairs + 1}: {first[0]!r} vs {second[0]!r}')
        lengths = (first[1], second[1])
        bases += sum(lengths)
        pairs += 1
        minimum = min(lengths) if minimum is None else min(minimum, *lengths)
        maximum = max(maximum, *lengths)
    if not pairs:
        raise ValueError('FASTQ pair contains no reads')
    return {'sample': sample, 'read_pairs': pairs, 'total_reads': pairs * 2, 'total_bases': bases,
            'mean_read_length': bases / (pairs * 2), 'min_read_length': minimum,
            'max_read_length': maximum, 'method': 'all clean FASTQ records, both mates',
            'r1': str(Path(r1).resolve()), 'r2': str(Path(r2).resolve())}


def combine(samples, reports):
    if not samples or len(set(samples)) != len(samples):
        raise ValueError('Expected sample IDs must be nonempty and unique')
    lengths = {}
    for path in reports:
        report = json.loads(Path(path).read_text())
        sample, value = report.get('sample'), report.get('mean_read_length')
        if sample in lengths or sample not in samples:
            raise ValueError(f'Duplicate or unexpected sample: {sample!r}')
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError(f'Invalid mean read length for {sample}')
        lengths[sample] = value
    if set(lengths) != set(samples):
        raise ValueError('Read length reports do not match expected sample IDs')
    return {sample: lengths[sample] for sample in samples}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('measure')
    for name in ('sample', 'r1', 'r2', 'output'):
        p.add_argument('--' + name, required=True)
    p = sub.add_parser('combine')
    p.add_argument('--samples', nargs='+', required=True)
    p.add_argument('--reports', nargs='+', required=True)
    p.add_argument('--output', required=True)
    args = parser.parse_args()
    inputs = [args.r1, args.r2] if args.command == 'measure' else args.reports
    out = Path(args.output)
    if any(out.resolve() == Path(p).resolve() or (out.exists() and out.samefile(p)) for p in inputs):
        parser.error('Output aliases input')
    out.unlink(missing_ok=True)
    result = measure(args.sample, args.r1, args.r2) if args.command == 'measure' else combine(args.samples, args.reports)
    atomic_json(out, result)


if __name__ == '__main__':
    main()
