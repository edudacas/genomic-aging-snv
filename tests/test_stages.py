import tempfile
import unittest
from pathlib import Path

from aging.common import atomic, sha
from aging.stages import STAGES, execute_stage, resource_budget


class StageTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = Path(self.tmp.name)
        self.rows = [{'sample_id': 'a'}, {'sample_id': 'b'}]
        self.calls = []

    def callback(self, row, directory):
        self.calls.append((directory.parent.name, row['sample_id']))
        atomic(directory / 'result.txt', row['sample_id'])

    def run_stage(self, stage, callback=None, fp='frozen'):
        return execute_stage(self.out, stage, self.rows, fp, 2, callback or self.callback)

    def test_stage_order_and_no_automatic_successors(self):
        with self.assertRaises(FileNotFoundError):
            self.run_stage('sage')
        for stage in STAGES:
            self.run_stage(stage)
            self.assertTrue((self.out / (stage.upper() + '_COMPLETE.json')).exists())
            if stage != 'lofreq':
                self.assertFalse((self.out / 'CALLERS_COMPLETE.json').exists())
        self.assertEqual(len(self.calls), 8)
        self.run_stage('mutect')
        self.assertEqual(len(self.calls), 8)

    def test_failure_resume_reuses_completed_samples(self):
        def fails(row, directory):
            if row['sample_id'] == 'b':
                raise RuntimeError('synthetic interruption')
            self.callback(row, directory)
        with self.assertRaises(RuntimeError):
            self.run_stage('qc', fails)
        self.assertFalse((self.out / 'QC_COMPLETE.json').exists())
        digest = sha(self.out / 'qc/a/COMPLETE.json')
        receipt = self.run_stage('qc')
        self.assertEqual(digest, sha(self.out / 'qc/a/COMPLETE.json'))
        self.assertTrue(receipt['samples']['a']['reused'])
        self.assertEqual(self.calls.count(('qc', 'a')), 1)
        self.assertTrue(list((self.out / 'qc').glob('b.failed.*')))

    def test_corruption_and_config_change_fail_closed(self):
        self.run_stage('qc')
        with self.assertRaises(ValueError):
            self.run_stage('qc', fp='different')
        atomic(self.out / 'qc/a/result.txt', 'corrupted')
        with self.assertRaises(ValueError):
            self.run_stage('sage')

    def test_budget(self):
        c = dict(sample_workers=12, cpus=16, memory_gib=96, threads=2, java_heap_gb=4)
        self.assertEqual(resource_budget(c), 12)
        for key, value in [('sample_workers', 13), ('cpus', 14), ('memory_gib', 97), ('threads', 4)]:
            with self.assertRaises(ValueError):
                resource_budget(c | {key: value})


if __name__ == '__main__':
    unittest.main()
