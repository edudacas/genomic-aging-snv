import unittest,json
from unittest.mock import patch
from aging.optimized import recount_sample
from aging.window_recount import recount_window
from aging.common import read_tsv
import test_explore

class ResumeTests(unittest.TestCase):
    def test_resume_uses_only_verified_complete_chunks(self):
        fixture=test_explore.ExploreIntegrationTests();fixture.setUp();self.addCleanup(fixture.tearDown)
        p=fixture.p;vcf=fixture.vcf(p);paths={n:vcf for n in ['mutect','lofreq','sage']}
        c=fixture.c;c['pon_by_assay']={'TEST':vcf};c['recount']={'window_bp':5,'gap_bp':5};sample=dict(fixture.row,assay_id='TEST')
        sites=[dict(chrom='chr1',pos=n,ref='A',alt='T',variant_id=f'chr1:{n}:A:T') for n in [20,30,150]];pop={r['variant_id']:{} for r in sites}
        out=p/'recount';seen=[]
        def interrupted(*args):
            seen.append(args[2][0]['pos'])
            if len(seen)==2:raise RuntimeError('synthetic interruption')
            return recount_window(*args)
        with patch('aging.optimized.recount_window',side_effect=interrupted):
            with self.assertRaises(RuntimeError):recount_sample(sample,sites,c,paths,pop,out,'fp')
        self.assertEqual(json.loads((out/'progress.json').read_text())['stage'],'FAILED')
        with patch('aging.optimized.recount_window',wraps=recount_window) as counter:
            result=recount_sample(sample,sites,c,paths,pop,out,'fp');self.assertEqual(counter.call_count,2)
        self.assertEqual(len(read_tsv(result)),3)
        timing=json.loads((out/'timing.json').read_text())
        self.assertEqual(timing['computed_windows'],2);self.assertEqual(timing['cached_windows'],1)
        for field in ('pileup_elapsed_seconds','integration_elapsed_seconds','elapsed_seconds'):
            self.assertGreaterEqual(timing[field],0)
        self.assertLessEqual(timing['pileup_elapsed_seconds']+timing['integration_elapsed_seconds'],timing['elapsed_seconds'])
        receipt=json.loads((out/'chunks/0000001.json').read_text())
        self.assertIn('pileup_elapsed_seconds',receipt);self.assertIn('integration_elapsed_seconds',receipt)
        # A tampered chunk is re-counted instead of reused.
        with (out/'chunks/0000000.tsv').open('a') as f:f.write('corrupt\n')
        with patch('aging.optimized.recount_window',wraps=recount_window) as counter:
            recount_sample(sample,sites,c,paths,pop,out,'fp');self.assertEqual(counter.call_count,1)
        with self.assertRaises(ValueError):recount_sample(sample,sites,c,paths,pop,out,'changed')
if __name__=='__main__':unittest.main()
