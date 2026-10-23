"""Validate optional StringTie arguments and protect workflow input files."""
from pathlib import Path
import re
import shlex


def stringtie_arguments(value, workdir):
    """Keep extra arguments within single-sample, reference-only BAM estimation."""
    if not isinstance(value, str):
        raise ValueError('stringtie.extra must be a string')
    tokens = shlex.split(value)
    switches = {'-t', '-v', '-u', '--conservative'}
    values = {'-l', '-f', '-m', '-a', '-j', '-c', '-s', '-g', '-M', '-E', '-x', '--ptf'}
    forbidden = {'-G', '-o', '-A', '-p', '-B', '-b', '-e', '-L', '-R',
                 '--mix', '--merge', '--fr', '--rf', '--ref', '--cram-ref',
                 '-h', '--help', '--version', '--viral'}
    parsed, inputs = [], []
    position = 0
    while position < len(tokens):
        token = tokens[position]
        option, sep, argument = token.partition('=')
        if token.startswith('-') and not token.startswith('--') and len(token) > 2:
            option, sep, argument = token[:2], '=', token[2:]
        if option in {'--fr', '--rf'}:
            raise ValueError('Library direction is managed by stranded and RSeQC; remove --fr/--rf from stringtie.extra.')
        if option in forbidden:
            raise ValueError(f'stringtie.extra cannot override workflow input, output, threads, strand or mode: {token}')
        if token in switches:
            parsed.append(token)
            position += 1
            continue
        if option not in values:
            raise ValueError(f'Unsupported stringtie.extra argument {token!r}; extra BAM inputs and mode overrides are not allowed')
        if not sep:
            position += 1
            if position == len(tokens):
                raise ValueError(f'stringtie.extra option {option} requires a value')
            argument = tokens[position]
        if not argument or (argument.startswith('-') and not re.fullmatch(r'-\d+(?:\.\d*)?(?:[eE][+-]?\d+)?', argument)):
            raise ValueError(f'stringtie.extra option {option} requires a value, not another option')
        if option == '--ptf':
            path = Path(argument).expanduser()
            if not path.is_absolute():
                path = Path(workdir) / path
            path = path.resolve()
            if not path.is_file() or not path.stat().st_size:
                raise ValueError(f'StringTie point-feature input is missing or empty: {path}')
            argument = str(path)
            inputs.append(argument)
        parsed.extend((option, argument))
        position += 1
    return parsed, inputs


def protect_inputs(inputs, destinations):
    """Reject aliases before Snakemake removes outputs or a shell opens logs.

    destinations is an iterable of (human-readable producer, path) pairs.
    Check all samples together to catch cross-sample links as well.
    """
    sources, inodes = set(), set()
    for item in inputs:
        path = Path(item)
        sources.add(path.resolve())
        if path.exists():
            stat = path.stat()
            inodes.add((stat.st_dev, stat.st_ino))
    targets, target_inodes = {}, {}
    for label, item in destinations:
        path = Path(item)
        resolved = path.resolve()
        if resolved in sources:
            raise ValueError(f'Workflow output/log aliases an input: {label}: {path}')
        if resolved in targets:
            raise ValueError(f'Workflow output/log paths collide: {label} and {targets[resolved]}: {path}')
        if path.exists():
            if not path.is_file():
                raise ValueError(f'Workflow output/log is not a regular file: {label}: {path}')
            stat = path.stat()
            inode = (stat.st_dev, stat.st_ino)
            if inode in inodes:
                raise ValueError(f'Workflow output/log is linked to an input: {label}: {path}')
            if inode in target_inodes:
                raise ValueError(f'Workflow output/log files are linked: {label} and {target_inodes[inode]}: {path}')
            target_inodes[inode] = label
        targets[resolved] = label
    for target, label in targets.items():
        for parent in target.parents:
            if parent in sources or parent in targets:
                raise ValueError(f'Workflow output/log has a file as its parent: {label}: {target}')
