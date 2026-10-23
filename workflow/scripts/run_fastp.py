#!/usr/bin/env python3
"""Validate raw FASTQ pairs and run fastp with explicit polyG behavior."""
import argparse
from contextlib import ExitStack
import json
from pathlib import Path
import shlex
import subprocess
import tempfile

from read_lengths import measure


def bounded_int(minimum, maximum=None):
    def parse(value):
        try:
            number = int(value)
        except ValueError:
            raise argparse.ArgumentTypeError('must be an integer')
        if number < minimum or (maximum is not None and number > maximum):
            limit = f'{minimum}..{maximum}' if maximum is not None else f'>={minimum}'
            raise argparse.ArgumentTypeError(f'must be {limit}')
        return number
    return parse


def aliases(first, second):
    return first.resolve() == second.resolve() or (
        first.exists() and second.exists() and first.samefile(second)
    )


def validate_paths(inputs, outputs):
    for path in inputs:
        if not path.is_file() or not path.stat().st_size:
            raise ValueError(f'Input FASTQ is missing or empty: {path}')
    if aliases(*inputs):
        raise ValueError('R1 and R2 must be different files')
    for index, output in enumerate(outputs):
        if any(aliases(output, source) for source in inputs):
            raise ValueError(f'Output aliases raw input: {output}')
        if any(aliases(output, other) for other in outputs[:index]):
            raise ValueError(f'Output paths must be distinct: {output}')
        if output.exists() and not output.is_file():
            raise ValueError(f'Output is not a regular file: {output}')
        if any(output.resolve() in other.resolve().parents or other.resolve() in output.resolve().parents
               for other in outputs[:index]):
            raise ValueError('One output path cannot be the parent of another output')


def validate_report(path, raw_reads):
    report = json.loads(path.read_text())
    before = report['summary']['before_filtering']['total_reads']
    after = report['summary']['after_filtering']['total_reads']
    if any(isinstance(n, bool) or not isinstance(n, int) for n in (before, after)):
        raise ValueError('fastp report contains invalid read counts')
    if before != raw_reads or not 0 <= after <= before or after % 2:
        raise ValueError(f'fastp report read counts disagree with validated paired input: {before}, {after}, {raw_reads}')
    if report['filtering_result']['passed_filter_reads'] != after:
        raise ValueError('fastp passed-filter count disagrees with after-filtering count')
    return after


def run(args):
    inputs = [Path(args.r1).absolute(), Path(args.r2).absolute()]
    outputs = [Path(p).absolute() for p in (args.out_r1, args.out_r2, args.html, args.json)]
    validate_paths(inputs, outputs)
    # fastp itself can accept mismatched read IDs. Check all pairs before trimming
    # so that filtering cannot conceal a mismatch in a discarded read pair.
    raw = measure(args.sample, *inputs)
    print(f'[INFO] Validated raw sample {args.sample}: {raw["read_pairs"]} pairs, '
          f'{raw["total_bases"]} bases', flush=True)
    version = subprocess.run(['fastp', '--version'], check=True, capture_output=True, text=True)
    print('[INFO] ' + (version.stdout + version.stderr).strip(), flush=True)
    with ExitStack() as stack:
        staged = []
        for index, output in enumerate(outputs):
            output.parent.mkdir(parents=True, exist_ok=True)
            directory = Path(stack.enter_context(tempfile.TemporaryDirectory(prefix='.fastp-', dir=output.parent)))
            # Preserve compression semantics even when the configured suffix is plain FASTQ.
            filename = f'read{index + 1}.fastq' + ('.gz' if str(output).endswith('.gz') else '')
            if index >= 2:
                filename = 'report.html' if index == 2 else 'report.json'
            staged.append(directory / filename)
        command = ['fastp', '-i', str(inputs[0]), '-I', str(inputs[1]),
                   '-o', str(staged[0]), '-O', str(staged[1]),
                   '--detect_adapter_for_pe',
                   '--trim_poly_g' if args.trim_poly_g == 'true' else '--disable_trim_poly_g',
                   '-q', str(args.qualified_quality_phred),
                   '-u', str(args.unqualified_percent_limit),
                   '-n', str(args.n_base_limit), '-l', str(args.length_required),
                   '-h', str(staged[2]), '-j', str(staged[3]), '-w', str(args.threads)]
        print('[INFO] Command: ' + shlex.join(command), flush=True)
        subprocess.run(command, check=True)
        for output in staged:
            if not output.is_file():
                raise ValueError(f'fastp did not produce required output: {output}')
        if not staged[2].stat().st_size:
            raise ValueError('fastp produced an empty HTML report')
        after = validate_report(staged[3], raw['total_reads'])
        # Publish only after fastp succeeds and its report matches the raw input.
        # Keep reports even when all reads are filtered; downstream read-length
        # validation then prevents empty clean data from reaching alignment.
        for temporary, destination in zip(staged, outputs):
            temporary.replace(destination)
        print(f'[INFO] Published clean FASTQ and HTML/JSON reports for {args.sample}; '
              f'{after // 2} pairs retained; polyG={args.trim_poly_g}', flush=True)
        if not after:
            print('[WARN] No clean read pairs remain; downstream analysis will stop.', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ('sample', 'r1', 'r2', 'out-r1', 'out-r2', 'html', 'json'):
        parser.add_argument('--' + option, required=True)
    parser.add_argument('--threads', type=bounded_int(1), default=8)
    parser.add_argument('--trim-poly-g', choices=('true', 'false'), default='false')
    parser.add_argument('--qualified-quality-phred', type=bounded_int(0, 93), default=20)
    parser.add_argument('--unqualified-percent-limit', type=bounded_int(0, 100), default=30)
    parser.add_argument('--n-base-limit', type=bounded_int(0, 50), default=5)
    parser.add_argument('--length-required', type=bounded_int(1), default=40)
    args = parser.parse_args()
    try:
        run(args)
    except (ValueError, OSError, KeyError, TypeError, subprocess.CalledProcessError) as error:
        parser.exit(1, f'[ERROR] fastp failed for {args.sample}: {error}\n')


if __name__ == '__main__':
    main()
