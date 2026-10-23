#!/usr/bin/env python3
"""Independent regressions for reference completeness and Excel exact integers."""
import csv
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'workflow' / 'scripts'
sys.path.insert(0, str(SCRIPTS))
import prepDE


class QuantAuditTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='rnaseq-quant-audit-')
        self.directory = Path(self.temporary.name)
        self.reference = self.directory / 'reference.gtf'
        lines = []
        for tid, gid, exons in [('T1', 'G1', [(101, 200), (301, 400)]),
                                ('T2', 'G1', [(101, 200), (301, 500)]),
                                ('T3', 'G1', [(601, 800)]),
                                ('T4', 'G2', [(1001, 1200)])]:
            attrs = f'gene_id "{gid}"; transcript_id "{tid}";'
            lines.append(f'chr1\tref\ttranscript\t{exons[0][0]}\t{exons[-1][1]}\t.\t+\t.\t{attrs}\n')
            lines.extend(f'chr1\tref\texon\t{start}\t{end}\t.\t+\t.\t{attrs}\n' for start, end in exons)
        self.reference.write_text(''.join(lines))
        self.models = prepDE.parse_gtf(self.reference, True, reference=True)

    def tearDown(self):
        self.temporary.cleanup()

    def ballgown(self, missing='T1', cov=0, fpkm=0):
        path = self.directory / 't_data.ctab'
        with path.open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['t_id', 'chr', 'strand', 'start', 'end', 't_name', 'num_exons',
                             'length', 'gene_id', 'gene_name', 'cov', 'FPKM'])
            for index, (tid, model) in enumerate(self.models.items(), 1):
                writer.writerow([index, model['chrom'], model['strand'], model['start'], model['end'],
                                 tid, len(model['exons']), model['length'], model['canonical_gene'], '.',
                                 cov if tid == missing else model['coverage'], fpkm if tid == missing else 0])
        return path

    def test_missing_isoform_restored_only_with_explicit_zero(self):
        observed = {tid: dict(model) for tid, model in self.models.items() if tid != 'T1'}
        restored = prepDE.complete_from_ballgown(observed, self.models, self.ballgown(), 'S')
        self.assertEqual(set(restored), set(self.models))
        self.assertEqual(restored['T1']['coverage'], 0)
        self.assertEqual(restored['T1']['exons'], ((101, 200), (301, 400)))

    def test_nonzero_missing_isoform_is_never_imputed(self):
        for cov, fpkm in ((1, 0), (0, 1), (1, 1)):
            with self.subTest(cov=cov, fpkm=fpkm):
                observed = {tid: dict(model) for tid, model in self.models.items() if tid != 'T1'}
                with self.assertRaisesRegex(ValueError, 'nonzero Ballgown'):
                    prepDE.complete_from_ballgown(observed, self.models, self.ballgown(cov=cov, fpkm=fpkm), 'S')

    def test_ballgown_missing_wrong_model_and_stale_coverage_rejected(self):
        path = self.ballgown(); original = path.read_text()
        variants = [original.replace('1\tchr1\t+\t101\t400\tT1', '1\tchr1\t+\t102\t400\tT1'),
                    '\n'.join(original.splitlines()[:-1]) + '\n',
                    original.replace('T2\t2\t300\tG1', 'T2\t2\t301\tG1'),
                    original.replace('T2\t2\t300\tG1', 'T2\t2\t300\twrong')]
        for text in variants:
            path.write_text(text)
            with self.assertRaises(ValueError):
                prepDE.complete_from_ballgown(dict(self.models), self.models, path, 'S')
        path.write_text(original)
        observed = {tid: dict(model) for tid, model in self.models.items()}
        observed['T2']['coverage'] = 1
        with self.assertRaisesRegex(ValueError, 'coverage mismatch'):
            prepDE.complete_from_ballgown(observed, self.models, path, 'S')

    def test_reference_catches_models_missing_from_every_sample(self):
        lines = [line for line in self.reference.read_text().splitlines(True) if '"T1"' not in line]
        gtf = self.directory / 'sample.gtf'
        gtf.write_text('# stringtie input.bam -G reference.gtf -e\n' +
                       ''.join(line.rstrip('\n') + ' cov "0.0";\n' for line in lines))
        listing = self.directory / 'samples.txt'; listing.write_text(f'S\t{gtf}\n')
        command = [sys.executable, str(SCRIPTS / 'prepDE.py'), '-i', str(listing), '--gene-id-only',
                   '--reference', str(self.reference), '-g', str(self.directory / 'gene.csv'),
                   '-t', str(self.directory / 'tx.csv')]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('enable StringTie Ballgown output', result.stderr)
        self.assertFalse((self.directory / 'gene.csv').exists())
        result = subprocess.run(command + ['--ballgown'], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('t_data.ctab', result.stderr)
        self.ballgown()
        result = subprocess.run(command + ['--ballgown'], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('T1', result.stderr)

    def test_real_stringtie_multi_isoform_two_sample_counts_and_tpm(self):
        stringtie = os.environ.get('RNASEQ_TEST_STRINGTIE') or shutil.which('stringtie')
        samtools = os.environ.get('RNASEQ_TEST_SAMTOOLS') or shutil.which('samtools')
        if not stringtie or not samtools:
            self.skipTest('Requires existing StringTie and samtools')
        gtfs, abundance = [], []
        for sample, pos, cigar in [('A', 101, '100M100N50M'), ('B', 1501, '150M')]:
            directory = self.directory / sample; directory.mkdir()
            sam, bam, gtf = directory / 'reads.sam', directory / 'reads.bam', directory / 'sample.gtf'
            sam.write_text('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n' + ''.join(
                f'r{i}\t0\tchr1\t{pos}\t60\t{cigar}\t*\t0\t0\t' + 'A' * 150 + '\t' + 'I' * 150 + '\tXS:A:+\n'
                for i in range(20)))
            commands = [[samtools, 'view', '-bo', str(bam), str(sam)],
                        [stringtie, str(bam), '-G', str(self.reference), '-e', '-B', '-p', '1',
                         '-o', str(gtf), '-A', str(directory / (sample + '_abundance.tab'))]]
            for command in commands:
                result = subprocess.run(command, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
            gtfs.append((sample, gtf)); abundance.append(directory / (sample + '_abundance.tab'))
        self.assertNotIn('T1', prepDE.parse_gtf(gtfs[0][1], True))
        self.assertIn('T1', prepDE.parse_gtf(gtfs[1][1], True))
        listing = self.directory / 'samples.txt'; listing.write_text(''.join(f'{s}\t{g}\n' for s, g in gtfs))
        command = [sys.executable, str(SCRIPTS / 'prepDE.py'), '-i', str(listing), '--gene-id-only',
                   '--reference', str(self.reference), '--ballgown', '-l', '150',
                   '-g', str(self.directory / 'raw.csv'), '-t', str(self.directory / 'tx.csv')]
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / 'raw.csv').open() as handle:
            self.assertEqual(list(csv.reader(handle)), [['gene_id', 'A', 'B'], ['G1', '20', '0'], ['G2', '0', '0']])
        with (self.directory / 'tx.csv').open() as handle:
            self.assertEqual(list(csv.reader(handle))[1:], [['T1', '0', '0'], ['T2', '20', '0'], ['T3', '0', '0'], ['T4', '0', '0']])
        command = [sys.executable, str(SCRIPTS / 'build_expression_matrices.py'), '--annotation', str(self.reference),
                   '--raw-count', str(self.directory / 'raw.csv'), '--out-count', str(self.directory / 'count.csv'),
                   '--out-fpkm', str(self.directory / 'fpkm.csv'), '--out-tpm', str(self.directory / 'tpm.csv')] + list(map(str, abundance))
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / 'tpm.csv').open() as handle:
            self.assertEqual(list(csv.reader(handle))[1:], [['G1', '1000000', '0'], ['G2', '0', '0']])

    def test_ballgown_invalid_estimates_rejected(self):
        for cov, fpkm in ((-1, 0), (float('nan'), 0), (float('inf'), 0),
                          (0, -1), (0, float('nan')), (0, float('inf'))):
            with self.subTest(cov=cov, fpkm=fpkm):
                path = self.ballgown(cov=cov, fpkm=fpkm)
                with self.assertRaisesRegex(ValueError, 'non-finite or negative'):
                    prepDE.complete_from_ballgown(dict(self.models), self.models, path, 'S')

    def test_real_two_expressed_splice_isoforms_aggregate_read_estimates(self):
        stringtie = os.environ.get('RNASEQ_TEST_STRINGTIE') or shutil.which('stringtie')
        samtools = os.environ.get('RNASEQ_TEST_SAMTOOLS') or shutil.which('samtools')
        if not stringtie or not samtools:
            self.skipTest('Requires existing StringTie and samtools')
        reference, sam = self.directory / 'branch.gtf', self.directory / 'branch.sam'
        reference.write_text(''.join(
            f'chr1\tref\t{feature}\t{start}\t{end}\t.\t+\t.\tgene_id "G"; transcript_id "{tid}";\n'
            for tid, second in [('T1', 301), ('T2', 501)]
            for feature, start, end in [('transcript', 101, second + 99), ('exon', 101, 200), ('exon', second, second + 99)]))
        sam.write_text('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:2000\n' + ''.join(
            f'r{index}_{read}\t0\tchr1\t101\t60\t{cigar}\t*\t0\t0\t' + 'A'*150 + '\t' + 'I'*150 + '\tXS:A:+\n'
            for index, cigar in enumerate(['100M100N50M', '100M300N50M']) for read in range(20)))
        bam, gtf = self.directory / 'branch.bam', self.directory / 'branch.actual.gtf'
        listing = self.directory / 'branch.samples'; listing.write_text(f'S\t{gtf}\n')
        commands = [[samtools, 'view', '-bo', str(bam), str(sam)],
                    [stringtie, str(bam), '-G', str(reference), '-e', '-B', '-p', '1', '-o', str(gtf)],
                    [sys.executable, str(SCRIPTS / 'prepDE.py'), '-i', str(listing), '--gene-id-only',
                     '--reference', str(reference), '--ballgown', '-l', '150',
                     '-g', str(self.directory / 'counts.csv'), '-t', str(self.directory / 'tx.csv')]]
        for command in commands:
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
        with (self.directory / 'counts.csv').open() as handle:
            self.assertEqual(list(csv.reader(handle)), [['gene_id', 'S'], ['G', '40']])
        with (self.directory / 'tx.csv').open() as handle:
            self.assertEqual(list(csv.reader(handle)), [['transcript_id', 'S'], ['T1', '20'], ['T2', '20']])

    def test_excel_rejects_changed_integer_and_preserves_previous_output(self):
        python = os.environ.get('RNASEQ_TEST_WORKFLOW_PYTHON', sys.executable)
        if subprocess.run([python, '-c', 'import openpyxl'], capture_output=True).returncode:
            self.skipTest('Requires existing openpyxl')
        output = self.directory / 'prior.xlsx'; output.write_bytes(b'existing output')
        code = ("import sys; sys.path.insert(0, sys.argv[1]); from export_expression_xlsx import write_workbook; "
                "write_workbook(sys.argv[2], {'gene_count': (['GeneID', 'S'], [['G', 10000000000000001]])})")
        result = subprocess.run([python, '-c', code, str(SCRIPTS), str(output)], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('numeric value changed', result.stderr)
        self.assertEqual(output.read_bytes(), b'existing output')

    def test_excel_rejects_columns_beyond_excel_limit(self):
        python = os.environ.get('RNASEQ_TEST_WORKFLOW_PYTHON', sys.executable)
        if subprocess.run([python, '-c', 'import openpyxl'], capture_output=True).returncode:
            self.skipTest('Requires existing openpyxl')
        code = ("import sys; sys.path.insert(0, sys.argv[1]); from export_expression_xlsx import write_workbook; "
                "write_workbook(sys.argv[2], {'gene_count': (['GeneID']+[f'S{i}' for i in range(16384)], [['G']+[1]*16384])})")
        result = subprocess.run([python, '-c', code, str(SCRIPTS), str(self.directory / 'too-wide.xlsx')], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Excel worksheet', result.stderr)
        self.assertFalse((self.directory / 'too-wide.xlsx').exists())


if __name__ == '__main__':
    unittest.main()
