"""Exercise alignment MAPQ selection through the wrapper and real samtools."""
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]


class AlignmentMapqTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.samtools = os.environ.get("RNASEQ_TEST_SAMTOOLS") or shutil.which("samtools")
        if not cls.samtools:
            raise unittest.SkipTest("Requires installed samtools")

    def test_default_zero_preserves_low_mapq_in_both_execution_modes(self):
        # Use a controlled SAM source because HISAT2's confidence for an ambiguous
        # fragment depends on its other candidate loci. The real filter/sort/BAM
        # reader are exercised; each cutoff includes a boundary observation.
        root = Path(tempfile.mkdtemp(prefix="rnaseq-mapq-regression-"))
        if os.environ.get("RNASEQ_TEST_KEEP") == "1":
            print("MAPQ regression artifacts: " + str(root), flush=True)
        else:
            self.addCleanup(shutil.rmtree, root)
        bindir = root / "bin"
        bindir.mkdir()
        (bindir / "samtools").symlink_to(Path(self.samtools).resolve())
        sam = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n"
        for index, mapq in enumerate((0, 29, 30, 60)):
            sam += f"mapped{mapq}\t0\tchr1\t{101 + index * 10}\t{mapq}\t4M\t*\t0\t0\tACGT\tIIII\n"
        sam += "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n"
        stub = bindir / "hisat2"
        stub.write_text(
            "#!" + sys.executable + "\n"
            "from pathlib import Path\nimport sys\n"
            "sam = " + repr(sam) + "\n"
            "args = sys.argv[1:]\n"
            "Path(args[args.index('--summary-file') + 1]).write_text('synthetic aligner summary\\n')\n"
            "if '-S' in args:\n"
            "    Path(args[args.index('-S') + 1]).write_text(sam)\n"
            "else:\n"
            "    sys.stdout.write(sam)\n"
        )
        stub.chmod(0o755)
        for mate in (1, 2):
            (root / f"reads{mate}.fastq").write_text(f"@one/{mate}\nACGT\n+\nIIII\n")
        for number in range(1, 9):
            (root / f"reference.{number}.ht2").write_bytes(b"synthetic stub index")
        env = dict(os.environ, PATH=str(bindir) + os.pathsep + os.environ.get("PATH", ""))
        evidence = []
        for threads in (1, 3):
            for cutoff in (None, 0, 30):
                with self.subTest(threads=threads, cutoff=cutoff):
                    label = f"threads{threads}-mapq{cutoff}"
                    bam = root / (label + ".bam")
                    command = [
                        sys.executable, str(ROOT / "workflow/scripts/align_reads.py"),
                        "--r1", str(root / "reads1.fastq"), "--r2", str(root / "reads2.fastq"),
                        "--index", str(root / "reference"), "--bam", str(bam),
                        "--summary", str(root / (label + ".summary.txt")),
                        "--threads", str(threads), "--hisat2-threads", "1",
                        "--view-threads", "0", "--sort-threads", "0",
                    ]
                    if cutoff is not None:
                        command += ["--mapq", str(cutoff)]
                    process = subprocess.run(command, env=env, capture_output=True, text=True, timeout=30)
                    (root / (label + ".log")).write_text(shlex.join(command) + "\n" + process.stdout + process.stderr)
                    self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
                    observed = subprocess.run([self.samtools, "view", str(bam)], check=True,
                                              capture_output=True, text=True, timeout=30)
                    rows = [line.split("\t") for line in observed.stdout.splitlines()]
                    expected = {"mapped30", "mapped60"} if cutoff == 30 else {
                        "mapped0", "mapped29", "mapped30", "mapped60", "unmapped"}
                    self.assertEqual({row[0] for row in rows}, expected)
                    self.assertEqual(len(rows), len(expected))
                    evidence.append({"threads": threads, "mapq_argument": cutoff,
                                     "retained": {row[0]: int(row[4]) for row in rows}})
        (root / "results.json").write_text(json.dumps(evidence, indent=2) + "\n")


if __name__ == "__main__":
    unittest.main(verbosity=2)
