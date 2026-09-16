import json
import tempfile
import threading
import unittest
from pathlib import Path
from aging.callers import run_jobs

class SampleWorkersTest(unittest.TestCase):
    def test_four_jobs_overlap_and_results_keep_input_order(self):
        barrier=threading.Barrier(4,timeout=5)
        guard=threading.Lock(); active=0; maximum=0
        def job(item,stage):
            nonlocal active, maximum
            with guard:
                active+=1;maximum=max(maximum,active)
            stage('TEST')
            barrier.wait()
            with guard: active-=1
            return item['sample_id']
        with tempfile.TemporaryDirectory() as td:
            items=[{'sample_id':str(i)} for i in range(4)]
            result=run_jobs(items,4,td,'TEST',job)
            self.assertEqual(result,['0','1','2','3'])
            self.assertEqual(maximum,4)
            self.assertEqual(json.loads((Path(td)/'progress.json').read_text())['completed'],4)
    def test_failed_sample_not_reported_as_complete(self):
        with tempfile.TemporaryDirectory() as td:
            run_jobs([{'sample_id':'x'}],1,td,'TEST',lambda r,s:({'explore_status':'FAILED'},None))
            states=json.loads((Path(td)/'progress.json').read_text())['worker_states']
            self.assertEqual(next(iter(states.values()))['stage'],'FAILED')
