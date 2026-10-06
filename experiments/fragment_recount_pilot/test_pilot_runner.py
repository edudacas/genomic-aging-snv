import json,tempfile,unittest
from pathlib import Path
from run_pilot import finalize

class ReceiptTests(unittest.TestCase):
    def test_incomplete_pair_coverage_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            (out/'PILOT_VALIDATION_PENDING.PRIVATE.json').write_text(json.dumps(dict(pilot_samples=['a','b','c'],
                candidate_alleles=10,evaluated_pairs=29,evaluation_states={'EVALUATED':29})))
            with self.assertRaisesRegex(ValueError,'Incomplete pair coverage'):finalize(out)
            self.assertFalse((out/'RECOUNT_PILOT_COMPLETE.PRIVATE.json').exists())
    def test_non_evaluable_cannot_be_approved(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            (out/'PILOT_VALIDATION_PENDING.PRIVATE.json').write_text(json.dumps(dict(pilot_samples=['a','b','c'],
                candidate_alleles=1,evaluated_pairs=3,evaluation_states={'REF_MISMATCH':1,'EVALUATED':2},outputs={},
                synthetic_tests_passed=True,native_tests_passed=True,failed_evaluation_pairs=1)))
            with self.assertRaisesRegex(ValueError,'Non-evaluable'):finalize(out)
            self.assertFalse((out/'RECOUNT_PILOT_COMPLETE.PRIVATE.json').exists())

if __name__=='__main__':unittest.main()
