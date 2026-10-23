import importlib.util
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / 'workflow/scripts'
def module(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / (name + '.py'))
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj
lengths, align = module('read_lengths'), module('align_reads')

class ReadLengthTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='rna lengths ')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
    def fastq(self, name, records):
        path = self.root / name
        with gzip.open(path, 'wt') as handle:
            for key, length in records:
                handle.write('@' + key + '\n' + 'A' * length + '\n+\n' + 'I' * length + '\n')
        return path
    def test_exact_mean_both_mates_and_variable_lengths(self):
        r1 = self.fastq('r1.gz', [('a/1',150),('b/1',100)])
        r2 = self.fastq('r2.gz', [('a/2',149),('b/2',101)])
        result = lengths.measure('001', r1, r2)
        self.assertEqual(result['mean_read_length'],125)
        self.assertEqual(result['total_bases'],500)
        self.assertEqual(result['total_reads'],4)
        self.assertEqual(result['sample'],'001')
    def test_fractional_read_length_preserved(self):
        r1 = self.fastq('r1.gz',[('a',150)])
        r2 = self.fastq('r2.gz',[('a',149)])
        self.assertEqual(lengths.measure('s',r1,r2)['mean_read_length'],149.5)
    def test_bad_pairs_rejected(self):
        for first,second in [([('a',75)],[('b',75)]),([('a',75)],[]),([('a/2',75)],[('a/2',75)]),([],[])]:
            with self.subTest(first=first,second=second):
                r1=self.fastq('r1.gz',first);r2=self.fastq('r2.gz',second)
                with self.assertRaises(ValueError):lengths.measure('s',r1,r2)
    def test_same_file_rejected(self):
        p=self.fastq('r.gz',[('a',75)])
        with self.assertRaises(ValueError):lengths.measure('s',p,p)
    def test_truncation_rejected_and_old_result_removed(self):
        r1=self.fastq('r1.gz',[('a',75)])
        r2=self.root/'r2.fastq';r2.write_text('@a\nAAA\n+\nII\n')
        out=self.root/'result.json';out.write_text('{"old":"success"}')
        process=subprocess.run([sys.executable,str(SCRIPTS/'read_lengths.py'),'measure','--sample','s','--r1',str(r1),'--r2',str(r2),'--output',str(out)],capture_output=True)
        self.assertNotEqual(process.returncode,0)
        self.assertFalse(out.exists())
    def test_combine_requires_exact_sample_set_and_numeric_lengths(self):
        a=self.root/'a.json';b=self.root/'b.json'
        a.write_text(json.dumps({'sample':'NA','mean_read_length':75}))
        b.write_text(json.dumps({'sample':'001','mean_read_length':149.5}))
        self.assertEqual(list(lengths.combine(['001','NA'],[a,b])),['001','NA'])
        for samples,reports in [(['001','NA'],[a]),(['NA'],[a,a]),(['other'],[a])]:
            with self.assertRaises(ValueError):lengths.combine(samples,reports)
        for value in [True,0,-1,float('nan'),float('inf'),'150']:
            a.write_text(json.dumps({'sample':'NA','mean_read_length':value}))
            with self.assertRaises(ValueError):lengths.combine(['NA'],[a])

class ResourceTests(unittest.TestCase):
    def test_pipeline_budgets_respect_total_and_tool_caps(self):
        for cores in range(3,65):
            mode,allocation=align.allocate_threads(cores,16,8,8)
            self.assertEqual(mode,'pipeline')
            self.assertEqual(sum(allocation),min(cores,34))
            self.assertTrue(all(1<=n<=cap for n,cap in zip(allocation,[16,9,9])))
    def test_one_or_two_cores_runs_sequentially(self):
        for cores in [1,2]:
            mode,allocation=align.allocate_threads(cores,16,8,8)
            self.assertEqual(mode,'sequential')
            self.assertLessEqual(max(allocation),cores)
    def test_managed_options_cannot_override_budget_or_inputs(self):
        for value in ['-p 99','-p99','--threads=99','-S overwrite.sam','-xother','--new-summary','-U single.fastq']:
            with self.subTest(value=value):
                with self.assertRaises(ValueError):align.extra_arguments(value)
        self.assertEqual(align.extra_arguments('--no-mixed --no-discordant'),['--no-mixed','--no-discordant'])
    def test_hisat_aliases_and_short_bundles_cannot_bypass_owned_inputs_outputs(self):
        values = [
            '--output alternate.sam', '--out alternate.sam', '--reads other.fastq',
            '--1 other.fastq', '--1other.fastq', '--2 other.fastq', '--2other.fastq',
            '--reference-string ACTG', '--ref-str ACTG', '--thread 99', '-qp99',
            '--metrics-file metrics.txt', '--log-file other.log',
            '--al-conc-disc paired.fastq', '--al-conc-disc-gz paired.gz',
            '--known-splicesite known.ss', '--novel-splicesite-i novel.ss', '--',
        ]
        for value in values:
            with self.subTest(value=value), self.assertRaises(ValueError):
                align.extra_arguments(value)

    def test_hisat_cannot_trim_or_subsample_after_clean_read_length_measurement(self):
        values = ['-3 10', '-310', '-q310', '--trim3=10', '--trim5 10', '-5 10',
                  '-510', '-q510', '--trim 10', '--skip 10', '-s10', '--upto 10',
                  '-u10', '--qupto 10', '--sample 0.5', '--read-lengths 75',
                  '--rna-strandness RF', '--rna-strand RF']
        for value in values:
            with self.subTest(value=value), self.assertRaises(ValueError):
                align.extra_arguments(value)

    def test_supported_alignment_tuning_and_quoted_values_are_preserved(self):
        import shlex
        value = '--no-mixed --no-discordant --score-min L,0,-0.2 -k5 -N 0 --rg-id "sample x"'
        self.assertEqual(align.extra_arguments(value), shlex.split(value))

    def test_relative_splice_input_remains_bound_to_original_workdir(self):
        with tempfile.TemporaryDirectory(prefix="splice ' input ") as directory:
            root=Path(directory);source=root/'known sites.ss';source.write_text('chr1\t1\t10\t+\n')
            staging=root/'stage';staging.mkdir()
            tokens=align.prepare_extra_paths(['--known-splicesite-infile='+str(source)],staging)
            self.assertEqual(align.extra_input_files(['--known-splicesite-infile','known sites.ss'], root),[str(source)])
            self.assertEqual(tokens[0],'--known-splicesite-infile')
            self.assertEqual((staging/tokens[1]).resolve(),source)
            self.assertFalse(any(c.isspace() for c in tokens[1]))

    def test_upstream_pipeline_failure_propagates(self):
        with self.assertRaises(RuntimeError):
            align.run_pipeline([[sys.executable,'-c','import sys;sys.exit(3)'],[sys.executable,'-c','import sys;sys.stdin.read()']])

if __name__=='__main__':unittest.main()
