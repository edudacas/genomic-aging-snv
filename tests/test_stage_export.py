import hashlib
import json
import tarfile
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
import threading
import time
from pathlib import Path

from aging.common import atomic, sha
from aging.stage_export import export_stage
from aging.stages import sample_task, stage_fingerprint, execute_stage
from aging.stage_contract import required_outputs


class ExportTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.run=Path(self.tmp.name)/'run';self.dest=Path(self.tmp.name)/'exports'
        self.run.mkdir()
        atomic(self.run/'QC_COMPLETE.json','frozen QC receipt fixture')

    def complete(self,sid):
        def populate(d,stage):
            for name in required_outputs(stage):atomic(d/name,'fixture')
        sample_task(self.run/'qc'/sid,stage_fingerprint('fp','qc'),lambda d: populate(d,'qc'))
        sample_task(self.run/'sage'/sid,stage_fingerprint('fp','sage',sha(self.run/'QC_COMPLETE.json')),
                    lambda d: populate(d,'sage'))

    def test_partial_incremental_read_only_and_archive_hashes(self):
        self.complete('a')
        atomic(self.run/'sage/b/in-progress.txt','partial')
        before={str(p.relative_to(self.run)):sha(p) for p in self.run.rglob('*') if p.is_file()}
        receipt=export_stage(self.run,self.dest,'sage','fp')
        self.assertEqual(set(receipt['samples']),{'a'})
        self.assertEqual(before,{str(p.relative_to(self.run)):sha(p) for p in self.run.rglob('*') if p.is_file()})
        with tarfile.open(self.dest/receipt['archive']) as archive:
            sums=archive.extractfile('SHA256SUMS').read().decode().splitlines()
            for line in sums:
                digest,name=line.split('  ',1)
                self.assertEqual(hashlib.sha256(archive.extractfile(name).read()).hexdigest(),digest)
        self.assertEqual(export_stage(self.run,self.dest,'sage','fp')['status'],'NO_NEW_SAMPLES')
        self.complete('b')
        self.assertEqual(set(export_stage(self.run,self.dest,'sage','fp')['samples']),{'b'})

    def test_reject_run_destination_and_changed_output(self):
        with self.assertRaises(ValueError):
            export_stage(self.run,self.run/'exports','sage','fp')
        self.complete('a')
        atomic(self.run/'sage/a/sage.vcf.gz','corrupt')
        with self.assertRaises(ValueError):
            export_stage(self.run,self.dest,'sage','fp')

    def test_live_export_does_not_wait_for_running_sample(self):
        # A real thread remains inside sample b while sample a is exported.
        ready=threading.Event();release=threading.Event()
        for sid in ('a','b'):
            sample_task(self.run/'qc'/sid,stage_fingerprint('fp','qc'),
                        lambda d:[atomic(d/name,'fixture') for name in required_outputs('qc')])
        def worker(sid):
            def publish(d):
                if sid=='b':
                    ready.set()
                    if not release.wait(5):raise RuntimeError('test timeout')
                for name in required_outputs('sage'):atomic(d/name,'fixture')
            return sample_task(self.run/'sage'/sid,stage_fingerprint('fp','sage',sha(self.run/'QC_COMPLETE.json')),publish)
        with ThreadPoolExecutor(max_workers=2) as pool:
            a=pool.submit(worker,'a');b=pool.submit(worker,'b')
            try:
                a.result(timeout=5);self.assertTrue(ready.wait(5));self.assertFalse(b.done())
                exported=export_stage(self.run,self.dest,'sage','fp')
                self.assertEqual(set(exported['samples']),{'a'});self.assertFalse(b.done())
            finally:release.set()
            b.result(timeout=5)
        self.assertEqual(set(export_stage(self.run,self.dest,'sage','fp')['samples']),{'b'})

    def test_required_index_missing_from_receipt_is_rejected(self):
        self.complete('a')
        marker=self.run/'sage/a/COMPLETE.json'
        receipt=json.loads(marker.read_text())
        del receipt['outputs']['sage.vcf.gz.tbi']
        atomic(marker,json.dumps(receipt))
        with self.assertRaises(ValueError):export_stage(self.run,self.dest,'sage','fp')


if __name__=='__main__':
    unittest.main()
