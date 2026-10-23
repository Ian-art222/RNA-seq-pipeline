"""Whole-DAG input protection and StringTie extra-input regressions."""
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
spec = importlib.util.spec_from_file_location('workflow_options', ROOT / 'workflow/scripts/workflow_options.py')
options = importlib.util.module_from_spec(spec)
spec.loader.exec_module(options)


class WorkflowOptionsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="workflow audit ' ")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_stringtie_extra_bams_and_modes_cannot_bypass_sample_inputs(self):
        for value in ('other.bam', '-c 1 other.bam', '--merge', '-R', '--mix', '--fr', '--rf',
                      '-Gother.gtf', '-oother.gtf', '-p12', '-B', '--ref ref.fa', '--version', '--'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                options.stringtie_arguments(value, self.root)
        args, inputs = options.stringtie_arguments('-c1.5 -g 50 -t -u --conservative -l label', self.root)
        self.assertEqual(args, ['-c', '1.5', '-g', '50', '-t', '-u', '--conservative', '-l', 'label'])
        self.assertEqual(inputs, [])

    def test_point_feature_input_is_resolved_and_declared(self):
        ptf = self.root / "point features ' test.tsv"
        ptf.write_text('chr1\t10\t+\tTSS\n')
        import shlex
        for flag in ('--ptf ' + shlex.quote(ptf.name), '--ptf=' + shlex.quote(ptf.name)):
            args, inputs = options.stringtie_arguments(flag, self.root)
            self.assertEqual(args, ['--ptf', str(ptf)])
            self.assertEqual(inputs, [str(ptf)])
        for value in ('--ptf absent.tsv', '--ptf', '-c', '-l --version'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                options.stringtie_arguments(value, self.root)

    def test_all_original_inputs_are_protected_from_outputs_and_logs(self):
        source = self.root / 'input.fastq'; source.write_text('original')
        link = self.root / 'log'; link.symlink_to(source)
        hard = self.root / 'hard'; os.link(source, hard)
        for path in (source, link, hard):
            with self.subTest(path=path), self.assertRaises(ValueError):
                options.protect_inputs([source], [('test log', path)])
            self.assertEqual(source.read_text(), 'original')

    def test_cross_sample_output_collisions_and_file_parents_are_rejected(self):
        first = self.root / 'out'; first.write_text('old')
        second = self.root / 'other'; os.link(first, second)
        cases = [(first, first), (first, second), (self.root/'new', self.root/'new/child')]
        for one, two in cases:
            with self.subTest(one=one,two=two), self.assertRaises(ValueError):
                options.protect_inputs([], [('sample1', one), ('sample2', two)])
        options.protect_inputs([self.root/'source'], [('sample1',self.root/'a/out'), ('sample2',self.root/'b/out')])


class WorkflowDagAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.snakemake = os.environ.get('RNASEQ_TEST_SNAKEMAKE') or shutil.which('snakemake')
        if not cls.snakemake:
            raise unittest.SkipTest('Requires installed Snakemake')

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='rnaseq dag audit ')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        reference = self.root/'reference'; reference.mkdir()
        for i in range(1,9): (reference/f'idx.{i}.ht2').write_text('synthetic index for dry-run only')
        (reference/'ref.gtf').write_text('chr1\ttest\texon\t1\t40\t.\t+\t.\tgene_id "g"; transcript_id "t";\n')
        sample = self.root/'cleandata/read_lengths'; sample.mkdir(parents=True)
        for name,mate in (('read_lengths.log',1),('read_lengths_2.fastq',2)):
            (sample/name).write_text('@one/'+str(mate)+'\n'+'A'*40+'\n+\n'+'I'*40+'\n')
        self.raw = sample/'read_lengths.log'
        self.config = {'reference':{'hisat2_index_prefix':str(reference/'idx'),'annotation':str(reference/'ref.gtf')},
                       'sample_scan':{'r1_suffix':'.log','r2_suffix':'_2.fastq','strict':True}}

    def run_dag(self, dry=True, target=None):
        cfg=self.root/'config.json';cfg.write_text(json.dumps(self.config))
        cmd=[self.snakemake,'--snakefile',str(ROOT/'workflow/Snakefile'),'--directory',str(self.root),
             '--configfile',str(cfg),'--cores','1','--printshellcmds']
        if dry:cmd.append('--dry-run')
        if target:cmd.append(str(target))
        return subprocess.run(cmd,capture_output=True,text=True,timeout=45)

    def test_log_same_path_as_fastq_rejected_before_snakemake_deletes_input(self):
        self.config['paths']={'logs_dir':'cleandata'}
        before=self.raw.read_bytes()
        process=self.run_dag(dry=False,target=self.root/'alignment/read_lengths/read_lengths.json')
        self.assertNotEqual(process.returncode,0)
        self.assertIn('output/log aliases an input',process.stdout+process.stderr)
        self.assertEqual(self.raw.read_bytes(),before)

    def test_point_feature_is_in_stringtie_dependency_list(self):
        ptf=self.root/'point features.tsv';ptf.write_text('chr1\t10\t+\tTSS\n')
        self.config['stringtie']={'extra':'--ptf "point features.tsv"'}
        process=self.run_dag()
        self.assertEqual(process.returncode,0,process.stdout+process.stderr)
        self.assertTrue(any(line.strip().startswith('input:') and str(ptf) in line
                            for line in (process.stdout+process.stderr).splitlines()))
        ptf.unlink()
        process=self.run_dag()
        self.assertNotEqual(process.returncode,0)
        self.assertIn('point-feature input is missing',process.stdout+process.stderr)

    def test_extra_bam_and_output_collision_fail_dag(self):
        self.config['stringtie']={'extra':'other.bam'}
        process=self.run_dag()
        self.assertNotEqual(process.returncode,0)
        self.assertIn('extra BAM inputs',process.stdout+process.stderr)
        self.config.pop('stringtie')
        self.config['output']={'xlsx':'quantification/gene_count_matrix.csv'}
        process=self.run_dag()
        self.assertNotEqual(process.returncode,0)
        self.assertIn('output/log paths collide',process.stdout+process.stderr)


if __name__ == '__main__':
    unittest.main()
