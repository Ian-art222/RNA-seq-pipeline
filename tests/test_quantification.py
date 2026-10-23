#!/usr/bin/env python3
"""Regression checks for count estimation, gene/sample identity, and Excel export.

Real StringTie tests use RNASEQ_TEST_STRINGTIE when set (no installation).
Excel tests use RNASEQ_TEST_WORKFLOW_PYTHON when set.
"""
import csv
import importlib.util
import json
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
import build_expression_matrices as matrices
import prepDE


class QuantificationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='rnaseq-quant-tests-')
        self.directory = Path(self.temporary.name)
        self.env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')

    def tearDown(self):
        self.temporary.cleanup()

    def run_script(self, script, arguments, python=None):
        return subprocess.run([python or sys.executable, str(SCRIPTS / script)] + list(map(str, arguments)),
                              text=True, capture_output=True, env=self.env, timeout=60)

    def assert_success(self, result):
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def write_gtf(self, name='S1.gtf', genes=None, symbol=False):
        genes = genes or [('G1', 1, 150, 10), ('G2', 201, 350, 0)]
        path = self.directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        lines = ['# stringtie input.bam -e -G reference.gtf\n']
        for gid, start, end, coverage in genes:
            attrs = f'gene_id "{gid}"; transcript_id "{gid}.t1";'
            if symbol:
                attrs += ' ref_gene_name "Symbol";'
            lines.extend(f'chr1\tStringTie\t{feature}\t{start}\t{end}\t.\t+\t.\t{attrs} cov "{coverage}";\n' for feature in ('transcript', 'exon'))
        path.write_text(''.join(lines))
        return path

    def prep(self, samples, extra=()):
        listing = self.directory / 'samples.txt'
        listing.write_text(''.join(f'{sample}\t{path}\n' for sample, path in samples))
        return self.run_script('prepDE.py', ['-i', listing, '-g', self.directory / 'raw.csv',
                                            '-t', self.directory / 'tx.csv'] + list(extra))

    def read_csv(self, filename):
        with (self.directory / filename).open() as handle:
            return list(csv.reader(handle))

    def abundance(self, sample='S1', rows=None):
        rows = rows if rows is not None else [
            ['G1', '-', 'chr1', '+', 1, 150, 10, 1, 1000000],
            ['G2', '-', 'chr1', '+', 201, 350, 0, 0, 0]]
        path = self.directory / (sample + '_abundance.tab')
        with path.open('w', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t')
            writer.writerow(['Gene ID', 'Gene Name', 'Reference', 'Strand', 'Start', 'End', 'Coverage', 'FPKM', 'TPM'])
            writer.writerows(rows)
        return path

    def finalize(self, annotation, abundances, samples=None):
        command = ['--annotation', annotation, '--raw-count', self.directory / 'raw.csv',
                   '--out-count', self.directory / 'count.csv', '--out-fpkm', self.directory / 'fpkm.csv',
                   '--out-tpm', self.directory / 'tpm.csv']
        if samples:
            command += ['--samples'] + samples + ['--']
        return self.run_script('build_expression_matrices.py', command + abundances)

    def basic(self):
        annotation = self.write_gtf()
        self.assert_success(self.prep([('S1', annotation)], ['--gene-id-only', '-l', '150']))
        return annotation, self.abundance()

    def test_per_sample_fractional_read_lengths_and_path_spaces(self):
        a = self.write_gtf('directory with spaces/A.gtf', [('G', 1, 301, 10)])
        b = self.write_gtf('directory with spaces/B.gtf', [('G', 1, 301, 10)])
        lengths = self.directory / 'lengths.json'
        lengths.write_text(json.dumps({'A': 150.5, 'B': 75.25}))
        self.assert_success(self.prep([('B', b), ('A', a)], ['--gene-id-only', '--read-lengths', lengths]))
        self.assertEqual(self.read_csv('raw.csv'), [['gene_id', 'A', 'B'], ['G', '20', '40']])

    def test_gene_id_only_preserves_literal_pipe_and_legacy_default(self):
        gtf = self.write_gtf(genes=[('literal|G', 1, 150, 10)], symbol=True)
        self.assert_success(self.prep([('S1', gtf)], ['--gene-id-only', '-l', '150']))
        self.assertEqual(self.read_csv('raw.csv')[1], ['literal|G', '10'])
        self.assert_success(self.prep([('S1', gtf)], ['-l', '150']))
        self.assertEqual(self.read_csv('raw.csv')[1], ['literal|G|Symbol', '10'])

    def test_float_global_read_length(self):
        gtf = self.write_gtf(genes=[('G', 1, 301, 10)])
        self.assert_success(self.prep([('S1', gtf)], ['-l', '150.5']))
        self.assertEqual(self.read_csv('raw.csv')[1][1], '20')

    def test_duplicate_sample_rejected(self):
        gtf = self.write_gtf()
        result = self.prep([('S1', gtf), ('S1', gtf)])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Duplicate sample', result.stderr)

    def test_invalid_read_length_maps_rejected(self):
        gtf = self.write_gtf()
        for value in ({'wrong': 150}, {'S1': 0}, {'S1': -1}, {'S1': True}, {'S1': '150'}, {'S1': float('nan')}, {'S1': float('inf')}):
            with self.subTest(value=value):
                path = self.directory / 'lengths.json'; path.write_text(json.dumps(value))
                self.assertNotEqual(self.prep([('S1', gtf)], ['--read-lengths', path]).returncode, 0)
        path.write_text('{"S1":100,"S1":150}')
        self.assertNotEqual(self.prep([('S1', gtf)], ['--read-lengths', path]).returncode, 0)

    def test_missing_transcript_in_second_sample_rejected(self):
        a = self.write_gtf('A.gtf')
        b = self.write_gtf('B.gtf', [('G1', 1, 150, 0)])
        self.assertNotEqual(self.prep([('A', a), ('B', b)]).returncode, 0)

    def test_coverage_missing_nonfinite_negative_rejected(self):
        for value in ('NaN', 'inf', '-1'):
            with self.subTest(value=value):
                self.assertNotEqual(self.prep([('S1', self.write_gtf(genes=[('G', 1, 150, value)]))]).returncode, 0)
        gtf = self.write_gtf(); gtf.write_text(gtf.read_text().replace(' cov "10";', ''))
        self.assertNotEqual(self.prep([('S1', gtf)]).returncode, 0)

    def test_empty_gtf_missing_exons_and_duplicate_transcript_rejected(self):
        gtf = self.write_gtf(); original = gtf.read_text()
        for content in ('', '# stringtie input.bam -e -G reference.gtf\n',
                        ''.join(line for line in original.splitlines(True) if '\texon\t' not in line),
                        original + original.splitlines(True)[1]):
            with self.subTest(content=content[:50]):
                gtf.write_text(content)
                self.assertNotEqual(self.prep([('S1', gtf)]).returncode, 0)

    def test_explicit_zero_models_and_exact_sample_option(self):
        annotation, abundance = self.basic()
        self.assert_success(self.finalize(annotation, [abundance], ['S1']))
        self.assertEqual(self.read_csv('count.csv'), [['GeneID', 'S1'], ['G1', '10'], ['G2', '0']])

    def test_valid_all_zero_sample(self):
        annotation = self.write_gtf(genes=[('G', 1, 150, 0)])
        self.assert_success(self.prep([('S1', annotation)], ['--gene-id-only']))
        abundance = self.abundance(rows=[['G', '-', 'chr1', '+', 1, 150, 0, 0, 0]])
        self.assert_success(self.finalize(annotation, [abundance]))

    def test_gff3_ncrna_gene_is_retained(self):
        annotation, abundance = self.basic()
        gff = self.directory / 'reference.gff3'
        gff.write_text('##gff-version 3\n' + ''.join(
            f'chr1\tref\t{feature}\t{start}\t{end}\t.\t+\t.\t{attrs}\n'
            for feature, start, end, attrs in [('gene', 1, 150, 'ID=G1'), ('mRNA', 1, 150, 'ID=T1;Parent=G1'),
            ('exon', 1, 150, 'Parent=T1'), ('ncRNA_gene', 201, 350, 'ID=G2'),
            ('ncRNA', 201, 350, 'ID=T2;Parent=G2'), ('exon', 201, 350, 'Parent=T2')]))
        self.assert_success(self.finalize(gff, [abundance]))
        self.assertEqual([row[0] for row in self.read_csv('tpm.csv')[1:]], ['G1', 'G2'])

    def test_gff_escaped_ids_remain_literal(self):
        gff = self.directory / 'escaped.gff3'
        gff.write_text('chr1\tref\tgene\t1\t150\t.\t+\t.\tID=G%3B1\n'
                       'chr1\tref\tmRNA\t1\t150\t.\t+\t.\tID=T%3B1;Parent=G%3B1\n'
                       'chr1\tref\texon\t1\t150\t.\t+\t.\tParent=T%3B1\n')
        self.assertEqual(matrices.extract_ordered_gene_ids(gff), ['G%3B1'])
        (self.directory / 'raw.csv').write_text('gene_id,S1\nG%3B1,10\n')
        abundance = self.abundance(rows=[['G%3B1', '-', 'chr1', '+', 1, 150, 10, 1, 1000000]])
        self.assert_success(self.finalize(gff, [abundance]))

    def test_repeated_descriptive_gtf_tags_allowed(self):
        annotation = self.write_gtf()
        annotation.write_text(annotation.read_text().replace('cov "10";', 'cov "10"; tag "basic"; tag "canonical";'))
        self.assert_success(self.prep([('S1', annotation)], ['--gene-id-only']))
        self.assert_success(self.finalize(annotation, [self.abundance()]))

    def test_missing_gene_not_filled_with_zero(self):
        annotation, abundance = self.basic()
        (self.directory / 'raw.csv').write_text('gene_id,S1\nG1,10\n')
        result = self.finalize(annotation, [abundance])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('missing IDs', result.stderr)
        self.assertFalse((self.directory / 'count.csv').exists())

    def test_missing_abundance_gene_or_sample_rejected(self):
        annotation, abundance = self.basic()
        abundance.write_text('\n'.join(abundance.read_text().splitlines()[:2]) + '\n')
        self.assertNotEqual(self.finalize(annotation, [abundance]).returncode, 0)
        abundance = self.abundance('wrong')
        self.assertNotEqual(self.finalize(annotation, [abundance]).returncode, 0)

    def test_duplicate_abundance_sample_rejected(self):
        annotation, abundance = self.basic()
        self.assertNotEqual(self.finalize(annotation, [abundance, abundance]).returncode, 0)

    def test_suffix_removed_once(self):
        sample = 'S_abundance.tab_extra'
        annotation = self.write_gtf()
        self.assert_success(self.prep([(sample, annotation)], ['--gene-id-only']))
        self.assert_success(self.finalize(annotation, [self.abundance(sample)]))
        self.assertEqual(self.read_csv('tpm.csv')[0][1], sample)

    def test_disjoint_gene_bundles_are_added(self):
        annotation = self.write_gtf(genes=[('G', 1, 150, 10)])
        annotation.write_text(annotation.read_text() + 'chr1\tref\ttranscript\t201\t350\t.\t+\t.\tgene_id "G"; transcript_id "other"; cov "10";\nchr1\tref\texon\t201\t350\t.\t+\t.\tgene_id "G"; transcript_id "other";\n')
        self.assert_success(self.prep([('S1', annotation)], ['--gene-id-only']))
        abundance = self.abundance(rows=[['G', '-', 'chr1', '+', 1, 150, 10, 1, 500000],
                                        ['G', '-', 'chr1', '+', 201, 350, 10, 1, 500000]])
        self.assert_success(self.finalize(annotation, [abundance]))
        self.assertEqual(self.read_csv('tpm.csv')[1], ['G', '1000000'])
        self.assertEqual(self.read_csv('fpkm.csv')[1], ['G', '2'])

    def test_overlapping_abundance_duplicates_rejected(self):
        annotation, abundance = self.basic()
        line = abundance.read_text().splitlines()[1]
        with abundance.open('a') as handle:
            handle.write(line + '\n')
        self.assertNotEqual(self.finalize(annotation, [abundance]).returncode, 0)

    def test_raw_invalid_numeric_width_and_duplicates_rejected(self):
        annotation, abundance = self.basic()
        for text in ('gene_id,S1\nG1,NaN\nG2,0\n', 'gene_id,S1\nG1,-1\nG2,0\n',
                     'gene_id,S1\nG1,1.5\nG2,0\n', 'gene_id,S1\nG1,1.00000000000000001\nG2,0\n', 'gene_id,S1\nG1,\nG2,0\n',
                     'gene_id,S1\nG1,10\nG1,20\nG2,0\n', 'gene_id,S1,S1\nG1,1,2\nG2,0,0\n'):
            with self.subTest(text=text):
                (self.directory / 'raw.csv').write_text(text)
                self.assertNotEqual(self.finalize(annotation, [abundance]).returncode, 0)

    def test_abundance_invalid_numeric_empty_and_tpm_sum_rejected(self):
        annotation, abundance = self.basic()
        original = abundance.read_text()
        for content in ('', original.splitlines()[0] + '\n', original.replace('\t1000000', '\tNaN'),
                        original.replace('\t1000000', '\t-1'), original.replace('\t1000000', '\t500000'),
                        original.replace('\t1000000', '\t0')):
            with self.subTest(content=content[-20:]):
                abundance.write_text(content)
                self.assertNotEqual(self.finalize(annotation, [abundance]).returncode, 0)

    def test_legacy_combined_gene_ids_rejected(self):
        annotation = self.write_gtf(symbol=True)
        self.assert_success(self.prep([('S1', annotation)]))
        result = self.finalize(annotation, [self.abundance()])
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('--gene-id-only', result.stderr)

    def excel_fixture(self):
        python = os.environ.get('RNASEQ_TEST_WORKFLOW_PYTHON', sys.executable)
        available = subprocess.run([python, '-c', 'import openpyxl'], capture_output=True)
        if available.returncode:
            self.skipTest('Excel tests need existing openpyxl or RNASEQ_TEST_WORKFLOW_PYTHON')
        samples = ['001', 'NA']
        genes = ['001', 'NA', 'literal|G', '=1+1']
        (self.directory / 'alignment.csv').write_text('SampleID,Total_Read_Pairs\n001,20\nNA,20\n')
        for metric, value in [('count', 10), ('fpkm', 1), ('tpm', 250000)]:
            with (self.directory / (metric + '.csv')).open('w', newline='') as handle:
                writer = csv.writer(handle); writer.writerow(['GeneID'] + samples)
                writer.writerows([gid, value, value] for gid in genes)
        command = ['--alignment-report', self.directory / 'alignment.csv', '--count', self.directory / 'count.csv',
                   '--fpkm', self.directory / 'fpkm.csv', '--tpm', self.directory / 'tpm.csv', '--output', self.directory / 'result.xlsx']
        return python, genes, command

    def test_excel_preserves_numeric_na_and_formula_identifiers_as_text(self):
        python, genes, command = self.excel_fixture()
        self.assert_success(self.run_script('export_expression_xlsx.py', command, python))
        code = ('import openpyxl,sys,json; w=openpyxl.load_workbook(sys.argv[1]); '
                'print(json.dumps({"genes":[[c.value,c.data_type] for c in list(w["gene_count"].columns)[0]][1:],'
                '"samples":[r[0] for r in list(w["alignment_report"].values)[1:]],'
                '"numeric_type":w["gene_count"]["B2"].data_type}))')
        result = subprocess.run([python, '-c', code, str(self.directory / 'result.xlsx')], capture_output=True, text=True)
        self.assert_success(result)
        observed = json.loads(result.stdout)
        self.assertEqual(observed['genes'], [[gid, 's'] for gid in genes])
        self.assertEqual(observed['samples'], ['001', 'NA'])
        self.assertEqual(observed['numeric_type'], 'n')

    def test_excel_mismatched_samples_and_genes_rejected(self):
        python, genes, command = self.excel_fixture()
        tpm = self.directory / 'tpm.csv'; original = tpm.read_text()
        for content in (original.replace('GeneID,001,NA', 'GeneID,NA,001'), original.replace('literal|G', 'wrong')):
            with self.subTest(content=content[:30]):
                tpm.write_text(content)
                self.assertNotEqual(self.run_script('export_expression_xlsx.py', command, python).returncode, 0)
                self.assertFalse((self.directory / 'result.xlsx').exists())

    def test_real_stringtie_names_and_explicit_zero_models(self):
        stringtie = os.environ.get('RNASEQ_TEST_STRINGTIE') or shutil.which('stringtie')
        samtools = os.environ.get('RNASEQ_TEST_SAMTOOLS') or shutil.which('samtools')
        if not stringtie or not samtools:
            self.skipTest('Real-tool regression needs existing StringTie and samtools')
        annotation = self.write_gtf(symbol=True)
        sam = self.directory / 'reads.sam'
        sam.write_text('@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n' + ''.join(
            f'read{i}\t0\tchr1\t1\t60\t150M\t*\t0\t0\t' + 'A' * 150 + '\t' + 'I' * 150 + '\n' for i in range(20)))
        bam = self.directory / 'reads.bam'
        result = subprocess.run([samtools, 'view', '-b', '-o', str(bam), str(sam)], capture_output=True, text=True)
        self.assert_success(result)
        gtf, abundance = self.directory / 'actual.gtf', self.directory / 'S1_abundance.tab'
        result = subprocess.run([stringtie, str(bam), '-G', str(annotation), '-e', '-p', '1', '-o', str(gtf), '-A', str(abundance)], capture_output=True, text=True)
        self.assert_success(result)
        self.assert_success(self.prep([('S1', gtf)], ['--gene-id-only', '-l', '150']))
        self.assert_success(self.finalize(annotation, [abundance]))
        observed = self.read_csv('count.csv')
        self.assertEqual(observed[1], ['G1', '20'])
        self.assertEqual(observed[2], ['G2', '0'])


if __name__ == '__main__':
    unittest.main()
