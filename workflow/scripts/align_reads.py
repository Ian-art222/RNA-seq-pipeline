#!/usr/bin/env python3
"""Run HISAT2/filter/sort within the actual Snakemake CPU allocation."""
import argparse
import json
import re
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile


def extra_arguments(value):
    """Keep alignment inputs, read lengths and output files owned by the DAG.

    HISAT2 accepts long aliases, abbreviated long names and bundled short
    options. Check those forms as well as the usual spelling; otherwise, for
    example, --trim3 or -q3 silently changes reads after their length is measured.
    """
    tokens = shlex.split(value)
    managed = {
        '-p', '--threads', '-x', '--index', '-1', '--1', '-2', '--2',
        '-U', '--unpaired', '--reads', '--12', '-S', '--output',
        '--summary-file', '--new-summary', '--sra-acc', '--interleaved',
        '--tab5', '--tab6', '--qseq', '-c', '-f', '-F', '-r', '-b', '-z', '-C',
        '--color', '-Q', '--quals', '--Q1', '--Q2',
        '--ref-string', '--reference-string', '--novel-splicesite-outfile',
        '--met-file', '--metrics-file', '--log-file', '--temp-directory',
        '--bam', '--passthrough', '--wrapper', '--hadoopout',
        '--al', '--un', '--al-gz', '--un-gz', '--al-bz2', '--un-bz2',
        '--al-conc', '--un-conc', '--al-conc-gz', '--un-conc-gz',
        '--al-conc-bz2', '--un-conc-bz2', '--al-conc-disc',
        '--al-conc-disc-gz', '--al-conc-disc-bz2',
        '-3', '-5', '--trim3', '--trim5', '-s', '--skip', '-u', '--upto',
        '--qupto', '--sample', '--read-lengths', '--rna-strandness',
        '--sam-no-head', '--sam-nohead', '--sam-noHD', '--sam-no-hd',
        '--sam-nosq', '--sam-no-sq', '--sam-noSQ', '--no-head', '--no-hd',
        '--no-sq', '--no-HD', '--no-SQ',
    }
    path_options = {'--known-splicesite-infile', '--novel-splicesite-infile'}
    protected_long = {option for option in managed | path_options if option.startswith('--')}
    # Remaining short options that consume a value. An attached value is not
    # another option (e.g. -k5 is valid); flags can be bundled (e.g. -at).
    short_values = set('vPkMINXiLDROgw')
    short_flags = set('qat')
    for token in tokens:
        option = token.split('=', 1)[0]
        incompatible = option in managed or token == '--'
        if option.startswith('--'):
            incompatible |= option.startswith(('--1', '--2'))
            incompatible |= any(full.startswith(option) and full != option for full in protected_long)
        elif token.startswith('-') and len(token) > 1:
            # Numeric negative values may be arguments to a preceding option.
            if not re.fullmatch(r'-[0-9]+(?:\.[0-9]+)?', token) or token[:2] in managed:
                for character in token[1:]:
                    if '-' + character in managed:
                        incompatible = True
                        break
                    if character in short_values:
                        break
                    if character not in short_flags:
                        incompatible = True
                        break
        if incompatible:
            raise ValueError(
                f'HISAT2 option {token!r} is managed by the workflow or incompatible with paired FASTQ; '
                'keep trimming in fastp, control library direction with stranded, and use full long option names'
            )
    return tokens


def extra_input_files(tokens, base_dir=None):
    paths = []
    index = 0
    base = Path.cwd() if base_dir is None else Path(base_dir)
    while index < len(tokens):
        option, separator, value = tokens[index].partition('=')
        if option in {"--known-splicesite-infile", "--novel-splicesite-infile"}:
            if not separator:
                index += 1
                if index >= len(tokens):
                    raise ValueError(f'Missing file for {option}')
                value = tokens[index]
            path = Path(value).expanduser()
            paths.append(str((base / path).resolve() if not path.is_absolute() else path.resolve()))
        index += 1
    return paths


def prepare_extra_paths(tokens, temporary):
    """Keep supported annotation inputs relative to the invocation workdir.

    HISAT2's own wrapper reconstructs a shell command; alias these files before
    changing cwd rather than handing it names requiring nested escaping.
    """
    path_options = {"--known-splicesite-infile", "--novel-splicesite-infile"}
    converted = []
    index = 0
    while index < len(tokens):
        token = tokens[index]
        option, separator, value = token.partition('=')
        if option in path_options:
            if not separator:
                index += 1
                if index >= len(tokens):
                    raise ValueError(f'Missing file for {option}')
                value = tokens[index]
            source = Path(value).expanduser().resolve()
            if not source.is_file():
                raise ValueError(f'{option}: input file does not exist: {source}')
            alias = f'extra_input_{index}.txt'
            (temporary / alias).symlink_to(source)
            converted.extend([option, alias])
            print('[EXTRA_INPUT] ' + json.dumps({'option': option, 'source': str(source), 'alias': alias}), file=sys.stderr)
        else:
            # Non-file values (e.g. a read group containing whitespace) reach
            # hisat2-align via the wrapper's shell command.
            converted.append(shlex.quote(token))
        index += 1
    return converted


def allocate_threads(total, hisat2, view, sort):
    """view/sort settings are additional threads; count every process main thread."""
    if total < 1 or hisat2 < 1 or view < 0 or sort < 0:
        raise ValueError('Invalid thread budget')
    requested = [hisat2, view + 1, sort + 1]
    if total < 3:
        return 'sequential', [min(total, cap) for cap in requested]
    allocation = [1, 1, 1]
    for _ in range(min(total, sum(requested)) - 3):
        candidates = [i for i in range(3) if allocation[i] < requested[i]]
        chosen = max(candidates, key=lambda i: (requested[i] - allocation[i]) / requested[i])
        allocation[chosen] += 1
    return 'pipeline', allocation


def run_pipeline(commands, cwd=None):
    processes = []
    try:
        previous = None
        for i, command in enumerate(commands):
            print('[COMMAND] ' + shlex.join(command), file=sys.stderr, flush=True)
            process = subprocess.Popen(command, cwd=cwd, stdin=previous, stdout=subprocess.PIPE if i < len(commands) - 1 else None)
            if previous is not None:
                previous.close()
            previous = process.stdout
            processes.append(process)
        codes = [process.wait() for process in processes]
        if any(codes):
            raise RuntimeError(f'HISAT2/view/sort pipeline failed: exit codes {codes}')
    finally:
        for process in processes:
            if process.poll() is None:
                process.terminate()
        for process in processes:
            process.wait()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('r1', 'r2', 'index', 'bam', 'summary'):
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--threads', required=True, type=int)
    parser.add_argument('--hisat2-threads', required=True, type=int)
    parser.add_argument('--view-threads', required=True, type=int)
    parser.add_argument('--sort-threads', required=True, type=int)
    parser.add_argument('--mapq', type=int, default=0)
    parser.add_argument('--hisat2-extra', default='')
    args = parser.parse_args()
    if not 0 <= args.mapq <= 255:
        parser.error('mapq must be in [0,255]')
    extra = extra_arguments(args.hisat2_extra)
    mode, allocation = allocate_threads(args.threads, args.hisat2_threads, args.view_threads, args.sort_threads)
    print('[RESOURCES] ' + json.dumps({'mode': mode, 'budget': args.threads,
          'hisat2_total': allocation[0], 'view_total': allocation[1], 'sort_total': allocation[2]}), file=sys.stderr)
    bam, summary = Path(args.bam).resolve(), Path(args.summary).resolve()
    bam.parent.mkdir(parents=True, exist_ok=True)
    summary.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='alignment-', dir=bam.parent) as directory:
        temporary = Path(directory)
        temp_summary, temp_bam = temporary / 'summary.txt', temporary / 'sorted.bam'
        # HISAT2 2.2.1's Perl wrapper joins its own shell command without escaping.
        # Use safe relative aliases for paths it handles, even when cwd contains quotes.
        aliases = {}
        for mate, source in ((1, args.r1), (2, args.r2)):
            alias = f'reads_{mate}.fastq' + ('.gz' if source.endswith('.gz') else '')
            (temporary / alias).symlink_to(Path(source).resolve())
            aliases[f'r{mate}'] = alias
        groups = [[Path(f'{args.index}.{i}.{extension}') for i in range(1, 9)] for extension in ('ht2', 'ht2l')]
        complete = [group for group in groups if all(path.is_file() and path.stat().st_size for path in group)]
        if len(complete) != 1:
            raise ValueError('HISAT2 requires one complete set of 8 index files')
        for index, source in enumerate(complete[0], 1):
            (temporary / f'reference.{index}{source.suffix}').symlink_to(source.resolve())
        print('[INPUTS] ' + json.dumps({'r1': str(Path(args.r1).resolve()), 'r2': str(Path(args.r2).resolve()),
              'index': str(Path(args.index).resolve()), 'aliases': aliases}), file=sys.stderr)
        extra = prepare_extra_paths(extra, temporary)
        hisat = ['hisat2', '--dta', '--summary-file', 'summary.txt', '-p', str(allocation[0]),
                 *extra, '-x', 'reference', '-1', aliases['r1'], '-2', aliases['r2']]
        view = ['samtools', 'view', '-u', '-@', str(allocation[1] - 1)]
        if args.mapq:
            view += ['-q', str(args.mapq)]
        sort = ['samtools', 'sort', '-@', str(allocation[2] - 1), '-o', str(temp_bam)]
        if mode == 'pipeline':
            run_pipeline([hisat, view + ['-'], sort + ['-']], cwd=temporary)
        else:
            sam, filtered = temporary / 'alignment.sam', temporary / 'filtered.bam'
            for command in [hisat + ['-S', 'alignment.sam'], view + ['-o', str(filtered), str(sam)], sort + [str(filtered)]]:
                print('[COMMAND] ' + shlex.join(command), file=sys.stderr, flush=True)
                subprocess.run(command, cwd=temporary, check=True)
        subprocess.run(['samtools', 'quickcheck', '-v', str(temp_bam)], check=True)
        if not temp_summary.is_file() or not temp_summary.stat().st_size:
            raise ValueError('HISAT2 did not produce a nonempty summary')
        temp_bam.replace(bam)
        temp_summary.replace(summary)


if __name__ == '__main__':
    main()
