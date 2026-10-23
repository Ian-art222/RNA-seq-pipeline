#!/usr/bin/env python3
"""Collect RSeQC evidence and gate paired-end StringTie quantification.

An inference failure is a diagnostic result, not a successful strand call.
Only ``gate`` creates the file that unlocks downstream quantification.
``summarize`` reports cohort evidence without gating individual samples.
"""

import argparse
import csv
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys


FRACTION_LABELS = {
    "mode1": 'Fraction of reads explained by "1++,1--,2+-,2-+"',
    "mode2": 'Fraction of reads explained by "1+-,1-+,2++,2--"',
    "failed_fraction": "Fraction of reads failed to determine",
}
FLAGS = {"--fr": "fr-secondstrand", "--rf": "fr-firststrand"}


def protect_paths(inputs, outputs):
    """Refuse output aliases before cleanup can overwrite an input or log."""
    inputs, outputs = [Path(path) for path in inputs], [Path(path) for path in outputs]

    def same_path(first, second):
        return first.resolve() == second.resolve() or (first.exists() and second.exists() and first.samefile(second))

    for index, output in enumerate(outputs):
        for source in inputs:
            if same_path(output, source):
                raise ValueError("Output {} aliases protected input {}".format(output, source))
        for previous in outputs[:index]:
            if same_path(output, previous):
                raise ValueError("Outputs {} and {} alias each other".format(output, previous))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def finite_fraction(value):
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(value)
        and 0 <= value <= 1
    )


def evaluate_evidence(report):
    """Validate parsed values, preserving all evidence on failure."""
    errors = report.setdefault("errors", [])
    warnings = report.setdefault("warnings", [])
    report["status"] = "FAIL"
    report["stringtie_flag"] = None
    report["library_type"] = None
    report["determinable_fraction"] = None
    report["dominant_fraction"] = None
    if report.get("protocol") != "PairEnd":
        errors.append("RSeQC must identify exactly one PairEnd dataset; single-end, mixed or unknown data are unsupported")
    usable = report.get("usable_reads")
    requested = report.get("requested_reads")
    if not isinstance(usable, int) or isinstance(usable, bool) or usable <= 0:
        errors.append("Actual usable reads must be a positive integer")
    elif not isinstance(requested, int) or requested <= 0 or usable > requested:
        errors.append("Actual usable reads exceed the requested sample size, or the request is invalid")
    elif usable < requested:
        warnings.append("Only {} of {} requested usable reads were sampled; no additional minimum-read threshold was applied".format(usable, requested))
    if not all(finite_fraction(report.get(key)) for key in FRACTION_LABELS):
        errors.append("All three fractions must be finite numbers in [0, 1]")
        return report
    mode1, mode2 = report["mode1"], report["mode2"]
    determinable = mode1 + mode2
    report["determinable_fraction"] = determinable
    if abs(determinable + report["failed_fraction"] - 1.0) > 0.0002 + 1e-12:
        errors.append("The two mode fractions and failed fraction do not sum to 1 (tolerance 0.0002)")
    if determinable <= 0:
        errors.append("No reads have a determinable strand direction")
        return report
    dominance = max(mode1, mode2) / determinable
    report["dominant_fraction"] = dominance
    thresholds = report.get("thresholds", {})
    minimum_determinable = thresholds.get("min_determinable")
    minimum_dominance = thresholds.get("min_dominance")
    if not finite_fraction(minimum_determinable) or not finite_fraction(minimum_dominance):
        errors.append("Thresholds must be finite numbers in [0, 1]")
        return report
    if determinable + 1e-12 < minimum_determinable:
        errors.append("Determinable fraction {:.6f} is below {:.6f}".format(determinable, minimum_determinable))
    if dominance + 1e-12 < minimum_dominance:
        errors.append("Dominant fraction {:.6f} is below {:.6f}".format(dominance, minimum_dominance))
    if mode1 == mode2:
        errors.append("The two strand modes are tied")
    if not errors:
        flag = "--fr" if mode1 > mode2 else "--rf"
        report.update(status="PASS", stringtie_flag=flag, library_type=FLAGS[flag])
    return report


def parse_rseqc(stdout, stderr, requested_reads=200000, min_determinable=0.50, min_dominance=0.80):
    """Parse unique, complete RSeQC evidence from both output streams."""
    combined = stdout + "\n" + stderr
    report = {
        "requested_reads": requested_reads,
        "thresholds": {"min_determinable": min_determinable, "min_dominance": min_dominance},
        "errors": [], "warnings": [], "raw_values": {}, "usable_reads": None,
    }
    protocols = re.findall(r"^\s*This is (\S+) Data\s*$", combined, re.MULTILINE)
    report["protocol"] = protocols[0] if len(protocols) == 1 and "Unknown Data type" not in combined else None
    report["raw_values"]["protocol"] = protocols
    for key, label in FRACTION_LABELS.items():
        values = re.findall(r"^\s*" + re.escape(label) + r":\s*([^\r\n]*)$", combined, re.MULTILINE)
        report["raw_values"][key] = values
        report[key] = None
        if len(values) != 1:
            report["errors"].append("Expected exactly one {} fraction, found {}".format(key, len(values)))
            continue
        try:
            value = float(values[0])
        except ValueError:
            report["errors"].append("{} is not numeric: {!r}".format(key, values[0]))
            continue
        if finite_fraction(value):
            report[key] = value
        else:
            report["errors"].append("{} is nonfinite or outside [0, 1]: {!r}".format(key, values[0]))
    # RSeQC's progress message uses end=' ', so a full sample reports
    # "Loading SAM/BAM file ... Total N usable reads were sampled" on one line.
    counts = re.findall(r"(?<!\S)Total[ \t]+(\S+)[ \t]+usable reads were sampled(?=$|\s)", combined)
    report["raw_values"]["usable_reads"] = counts
    if len(counts) != 1 or not re.fullmatch(r"[0-9]+", counts[0]):
        report["errors"].append("Expected exactly one integer 'Total N usable reads were sampled' record")
    else:
        report["usable_reads"] = int(counts[0])
    return evaluate_evidence(report)


def require_readable_file(path, label):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError("{} must be an existing, nonempty file: {}".format(label, path))
    with path.open("rb") as handle:
        handle.read(1)


def validate_inputs(bam_path, bai_path, bed_path):
    """Check the index, reference coordinates and absence of strand constraints."""
    for path, label in ((bam_path, "BAM"), (bai_path, "BAI"), (bed_path, "BED12")):
        require_readable_file(path, label)
    import pysam
    from prepare_strand_bed import read_bed12

    records = read_bed12(bed_path)
    if not records:
        raise ValueError("BED12 contains no transcript records")
    with pysam.AlignmentFile(str(bam_path), "rb", index_filename=str(bai_path), require_index=True) as bam:
        if not bam.is_bam:
            raise ValueError("The input is not a BAM file")
        if not bam.check_index():
            raise ValueError("BAM index is not usable")
        header = bam.header.to_dict()
        if header.get("HD", {}).get("SO") != "coordinate":
            raise ValueError("BAM header must declare coordinate sorting (HD:SO:coordinate)")
        lengths = dict(zip(bam.references, bam.lengths))
        if not lengths or len(lengths) != len(bam.references) or any(length <= 0 for length in lengths.values()):
            raise ValueError("BAM SQ reference names/lengths are missing, duplicated or invalid")
        for program in header.get("PG", []):
            command = program.get("CL", "")
            identity = " ".join(str(program.get(key, "")) for key in ("ID", "PN", "CL"))
            if "hisat2" in identity.lower() and re.search(r"(?:^|\s)--rna-strandness(?:=|\s|$)", command):
                raise ValueError("BAM was generated with HISAT2 --rna-strandness; independent inference requires an unconstrained alignment")
        for record in records:
            chrom = record["chrom"]
            if chrom not in lengths:
                raise ValueError("BED12 contig {!r} is absent from BAM SQ".format(chrom))
            if record["start"] < 0 or record["end"] > lengths[chrom]:
                raise ValueError("BED12 transcript {!r} exceeds BAM SQ coordinates".format(record["name"]))
    return {"bam_contigs": len(lengths), "bed12_transcripts": len(records), "sort_order": "coordinate"}


def infer(args):
    report = {
        "schema_version": 1, "sample": args.sample, "status": "FAIL",
        "errors": [], "warnings": [], "requested_reads": args.sample_size,
        "mapq": args.mapq, "usable_reads": None, "protocol": None,
        "mode1": None, "mode2": None, "failed_fraction": None,
        "determinable_fraction": None, "dominant_fraction": None,
        "stringtie_flag": None, "library_type": None,
        "thresholds": {"min_determinable": args.min_determinable, "min_dominance": args.min_dominance},
        "inputs": {"bam": str(Path(args.bam).resolve()), "bai": str(Path(args.bai).resolve()), "bed12": str(Path(args.bed12).resolve())},
        "tool": {"name": "RSeQC infer_experiment.py", "version": None, "returncode": None,
                 "command": ["infer_experiment.py", "-i", args.bam, "-r", args.bed12, "-s", str(args.sample_size), "-q", str(args.mapq)]},
    }
    log_dir = Path(args.log_dir)
    stdout_path, stderr_path = log_dir / "stdout.txt", log_dir / "stderr.txt"
    output_path, command_path = Path(args.output), log_dir / "command.json"
    try:
        protect_paths(
            [args.bam, args.bai, args.bed12],
            [output_path, output_path.with_name(output_path.name + ".tmp"), stdout_path, stderr_path,
             command_path, command_path.with_name(command_path.name + ".tmp")],
        )
    except (OSError, ValueError) as exc:
        print("Unsafe inference output paths; no files were changed: {}".format(exc), file=sys.stderr)
        return 1
    log_dir.mkdir(parents=True, exist_ok=True)
    stdout_path.write_text("")
    stderr_path.write_text("")
    report["logs"] = {"stdout": str(stdout_path.resolve()), "stderr": str(stderr_path.resolve()), "command": str((log_dir / "command.json").resolve())}
    Path(args.output).unlink(missing_ok=True)
    try:
        if not args.sample.strip() or any(char in args.sample for char in "\t\r\n"):
            raise ValueError("Sample ID must be nonempty and cannot contain tabs or newlines")
        if args.sample_size <= 0 or not 0 <= args.mapq <= 255:
            raise ValueError("Sample size must be positive and MAPQ must be in [0, 255]")
        if not finite_fraction(args.min_determinable) or not finite_fraction(args.min_dominance):
            raise ValueError("Thresholds must be finite numbers in [0, 1]")
        report["input_validation"] = validate_inputs(args.bam, args.bai, args.bed12)
        executable = shutil.which("infer_experiment.py")
        if not executable:
            raise ValueError("infer_experiment.py is not available in the active environment")
        command = [executable, "-i", args.bam, "-r", args.bed12, "-s", str(args.sample_size), "-q", str(args.mapq)]
        report["tool"]["command"] = command
        version_command = [executable, "--version"]
        version = subprocess.run(version_command, capture_output=True, text=True, errors="replace", check=False)
        report["tool"].update(
            version=(version.stdout + version.stderr).strip() or "unavailable",
            version_command=version_command, version_returncode=version.returncode,
            version_stdout=version.stdout, version_stderr=version.stderr,
        )
        if re.search(r"\b5\.0\.4\b", report["tool"]["version"]):
            # Checked against qcmodule.SAM.ParseBAM.configure_experiment in the
            # pinned RSeQC version. Do not describe these as random read pairs,
            # or interpret determinable_fraction as a whole-library exon rate.
            report["sampling"] = {
                "implementation": "RSeQC 5.0.4 configure_experiment",
                "selection": "First eligible alignment records in BAM order; not random genome-wide sampling",
                "unit": "Individual alignment records; paired-end mates are counted separately",
                "overlap": "BED transcript spans, including introns; exon blocks are not consulted",
                "read_interval": "Alignment start plus query alignment length; not the full spliced genomic interval",
                "limitation": "Fractions describe the selected annotation-overlapping records, not all input reads",
            }
        completed = subprocess.run(command, capture_output=True, text=True, errors="replace", check=False)
        stdout_path.write_text(completed.stdout)
        stderr_path.write_text(completed.stderr)
        report["tool"]["returncode"] = completed.returncode
        report.update(parse_rseqc(completed.stdout, completed.stderr, args.sample_size, args.min_determinable, args.min_dominance))
        if completed.returncode != 0:
            report["errors"].append("infer_experiment.py exited with status {}".format(completed.returncode))
            report.update(status="FAIL", stringtie_flag=None, library_type=None)
    except Exception as exc:
        report["errors"].append("{}: {}".format(type(exc).__name__, exc))
        report.update(status="FAIL", stringtie_flag=None, library_type=None)
        # JSON has no NaN/Infinity values, including when invalid CLI thresholds
        # themselves caused the failure. Preserve their literal forms instead.
        report["thresholds"] = {
            key: value if finite_fraction(value) else str(value)
            for key, value in report["thresholds"].items()
        }
    write_json(log_dir / "command.json", report["tool"])
    write_json(args.output, report)
    if report["status"] != "PASS":
        print("Strand inference failed for {}: {}".format(args.sample, "; ".join(report["errors"])), file=sys.stderr)
    return 0  # Let the gate collect diagnostics for every sample.


def check_report(report):
    if not isinstance(report, dict):
        raise ValueError("Report is not a JSON object")
    if report.get("status") != "PASS":
        raise ValueError("Inference failed: " + "; ".join(str(error) for error in report.get("errors", ["no PASS result"])))
    if report.get("errors") or report.get("tool", {}).get("returncode") != 0:
        raise ValueError("PASS report has errors or lacks a successful RSeQC exit status")
    checked = dict(report, errors=[], warnings=[])
    evaluate_evidence(checked)
    if checked["status"] != "PASS":
        raise ValueError("Invalid PASS evidence: " + "; ".join(checked["errors"]))
    if report.get("stringtie_flag") != checked["stringtie_flag"]:
        raise ValueError("Stored StringTie direction disagrees with the recorded fractions")


def diagnostic_rows(samples, reports):
    """Distinguish bad evidence from incomplete or incorrectly assigned reports."""
    rows, failures, structural, flags, observed = [], [], [], set(), set()
    if len(samples) != len(set(samples)):
        structural.append("Expected sample IDs are duplicated")
    if len(samples) != len(reports):
        structural.append("Expected sample and report counts differ")
    if not samples:
        structural.append("No expected samples were provided")
    for index, sample in enumerate(samples):
        report_path = reports[index] if index < len(reports) else ""
        row = {"sample": sample, "status": "FAIL", "report": report_path}
        structure_valid = False
        try:
            if not sample.strip() or any(char in sample for char in "\t\r\n"):
                raise ValueError("Sample ID must be nonempty and cannot contain tabs or newlines")
            if not report_path:
                raise ValueError("Missing report path for expected sample")
            report = json.loads(Path(report_path).read_text())
            if not isinstance(report, dict):
                raise ValueError("Report is not a JSON object")
            for key in ("mode1", "mode2", "failed_fraction", "determinable_fraction", "dominant_fraction", "usable_reads", "requested_reads", "stringtie_flag"):
                row[key] = report.get(key)
            row["warnings"] = "; ".join(str(warning) for warning in report.get("warnings", []))
            if report.get("sample") != sample:
                raise ValueError("Report sample {!r} does not match expected {!r}".format(report.get("sample"), sample))
            if sample in observed:
                raise ValueError("Duplicate sample report")
            observed.add(sample)
            structure_valid = True
            check_report(report)
            flags.add(report["stringtie_flag"])
            row["status"] = "PASS"
        except Exception as exc:
            row["reason"] = str(exc)
            failures.append("{}: {}".format(sample, exc))
            if not structure_valid:
                structural.append("{}: {}".format(sample, exc))
        rows.append(row)
    return rows, failures, structural, flags


def write_diagnostic(path, rows):
    diagnostic = Path(path)
    diagnostic.parent.mkdir(parents=True, exist_ok=True)
    fields = ["sample", "status", "mode1", "mode2", "failed_fraction", "determinable_fraction", "dominant_fraction", "usable_reads", "requested_reads", "stringtie_flag", "warnings", "reason", "report"]
    with diagnostic.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def gate(args):
    output = Path(args.output)
    try:
        protect_paths(args.reports, [output, args.diagnostic])
    except (OSError, ValueError) as exc:
        print("Unsafe gate output paths; no files were changed: {}".format(exc), file=sys.stderr)
        return 1
    output.unlink(missing_ok=True)
    rows, failures, structural, flags = diagnostic_rows(args.samples, args.reports)
    failures = list(dict.fromkeys(structural + failures))
    if len(flags) > 1:
        failures.append("Samples disagree on strand direction; expected one shared library protocol")
        for row in rows:
            if row["status"] == "PASS":
                row.update(status="FAIL", reason="Sample directions disagree across the cohort")
    if failures:
        rows.append({"sample": "__gate__", "status": "FAIL", "reason": "; ".join(failures)})
    write_diagnostic(args.diagnostic, rows)
    if failures:
        print("Stranded quantification blocked. See {}: {}".format(args.diagnostic, "; ".join(failures)), file=sys.stderr)
        return 1
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("STRINGTIE_STRAND={}\n".format(next(iter(flags))))
    return 0


def summarize(args):
    """Report failures and mixed protocols without creating a quantification gate.

    Failed or invalid strand evidence produces a FAIL row and a warning, with
    exit status zero. Missing/corrupt reports, incorrect sample assignments,
    duplicate IDs and count mismatches write diagnostics but return one.
    """
    try:
        protect_paths(args.reports, [args.diagnostic])
    except (OSError, ValueError) as exc:
        print("Unsafe summary output paths; no files were changed: {}".format(exc), file=sys.stderr)
        return 1
    rows, failures, structural, flags = diagnostic_rows(args.samples, args.reports)
    if structural:
        rows.append({"sample": "__summary__", "status": "FAIL", "reason": "; ".join(structural)})
    write_diagnostic(args.diagnostic, rows)
    if len(flags) > 1:
        print("WARNING: Samples have mixed FR/RF strand directions. Each passing sample uses its own inferred StringTie direction; verify the library protocols and sample metadata. See {}".format(args.diagnostic), file=sys.stderr)
    if failures:
        print("WARNING: Failed strand evidence is recorded in {}. Only the corresponding per-sample gates block quantification: {}".format(args.diagnostic, "; ".join(failures)), file=sys.stderr)
    if structural:
        print("Incomplete strand summary in {}: {}".format(args.diagnostic, "; ".join(structural)), file=sys.stderr)
        return 1
    return 0


def option(args):
    try:
        value = Path(args.env).read_text()
        match = re.fullmatch(r"STRINGTIE_STRAND=(--fr|--rf)\n?", value)
        if not match:
            raise ValueError("Expected exactly one STRINGTIE_STRAND=--fr or STRINGTIE_STRAND=--rf assignment")
        print(match.group(1))
        return 0
    except (OSError, ValueError) as exc:
        print("Invalid strand gate file: {}".format(exc), file=sys.stderr)
        return 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    inference = commands.add_parser("infer")
    for name in ("sample", "bam", "bai", "bed12", "output", "log-dir"):
        inference.add_argument("--" + name, required=True)
    inference.add_argument("--sample-size", type=int, default=200000)
    inference.add_argument("--mapq", type=int, default=30)
    inference.add_argument("--min-determinable", type=float, default=0.50)
    inference.add_argument("--min-dominance", type=float, default=0.80)
    inference.set_defaults(func=infer)
    gating = commands.add_parser("gate")
    gating.add_argument("--samples", nargs="+", required=True)
    gating.add_argument("--reports", nargs="+", required=True)
    gating.add_argument("--diagnostic", required=True)
    gating.add_argument("--output", required=True)
    gating.set_defaults(func=gate)
    summary = commands.add_parser("summarize")
    summary.add_argument("--samples", nargs="+", required=True)
    summary.add_argument("--reports", nargs="+", required=True)
    summary.add_argument("--diagnostic", required=True)
    summary.set_defaults(func=summarize)
    selection = commands.add_parser("option")
    selection.add_argument("--env", required=True)
    selection.set_defaults(func=option)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
