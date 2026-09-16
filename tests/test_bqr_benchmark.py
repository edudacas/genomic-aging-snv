import tempfile
import unittest
import argparse
from unittest.mock import patch
from pathlib import Path
from aging.common import atomic, write_json, load_config
from aging.bqr_benchmark import PILOTS, arm_config, compare_sample, raw_pass_records, vaf_bucket, run_benchmark


HEADER=('##fileformat=VCFv4.2\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR\n')


class BqrTests(unittest.TestCase):
    def test_exact_boundaries_and_missing(self):
        self.assertEqual(vaf_bucket('0.02'),'[0.02,0.03)')
        self.assertEqual(vaf_bucket('0'),'[0,0.02)')
        self.assertEqual(vaf_bucket('1'),'[0.50,1]')
        for value in ('.','NaN','Infinity',None,'-0.1','1.1'):
            self.assertEqual(vaf_bucket(value),'UNKNOWN')

    def test_raw_duplicates_pass_and_strata_migration(self):
        with tempfile.TemporaryDirectory() as temp:
            left,right=Path(temp)/'skip.vcf',Path(temp)/'enabled.vcf'
            record='chr1\t1\t.\tA\tC\t30\tPASS\tTIER=PANEL\tAF:AD\t0.02:98,2\n'
            unknown='chr1\t2\t.\tA\tG\t30\tPASS\t.\tAF:AD\t.:98,2\n'
            atomic(left,HEADER+record*2+unknown+record.replace('\tPASS\t','\t.\t'))
            atomic(right,HEADER+record.replace('0.02:','0.03:')+unknown)
            summary,rows=compare_sample('test',left,right,{})
            self.assertEqual(summary['skip_bqr']['pass_records'],3)
            self.assertEqual(summary['bqr_enabled']['pass_records'],2)
            self.assertEqual(summary['retained_pass_records'],2)
            self.assertEqual(summary['lost_pass_records'],1)
            self.assertEqual(sum(r['gained_in_bucket'] for r in rows),1)
            self.assertTrue(any(r['vaf_bin']=='UNKNOWN' and r['sage_tier']=='UNKNOWN' for r in rows))

    def test_arm_config_only_changes_bqr(self):
        base=dict(threads=2,java_heap_gb=4,sage_extra_args=['-skip_bqr','-skip_msi_jitter'],
                  profile='exploratory',callers=['sage'],evidence={'max_depth':100000})
        skip,enabled=arm_config(base,'skip_bqr'),arm_config(base,'enabled')
        self.assertEqual({k for k in skip if skip[k]!=enabled[k]},{'sage_extra_args'})
        self.assertEqual([x for x in skip['sage_extra_args'] if x!='-skip_bqr'],enabled['sage_extra_args'])

    def test_mismatched_sample_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            a,b=Path(temp)/'a.vcf',Path(temp)/'b.vcf'
            atomic(a,HEADER);atomic(b,HEADER.replace('TUMOR','OTHER'))
            with self.assertRaises(ValueError):compare_sample('x',a,b,{})

    def test_paired_runner_is_sage_only_and_resumes_all_six(self):
        import gzip
        with tempfile.TemporaryDirectory() as temp:
            root=Path(temp);config=root/'config.json'
            write_json(config,dict(threads=2,java_heap_gb=4,sage_extra_args=['-skip_bqr','-skip_msi_jitter']))
            rows=[dict(sample_id=sid,bam_path='/fixture.bam',bai_path='/fixture.bai') for sid in PILOTS]
            args=argparse.Namespace(config=str(config),workbook='fixture',linkage='fixture',out=str(root/'benchmark'))
            calls=[]
            def fake_sage(row,c,d,pon):
                self.assertEqual(c['callers'],['sage']);self.assertIsNone(pon)
                calls.append((row['sample_id'],tuple(c['sage_extra_args'])))
                with gzip.open(d/'sage.vcf.gz','wt') as stream:stream.write(HEADER)
                for name in ('sage.norm.vcf.gz','sage.norm.vcf.gz.tbi','sage.vcf.gz.tbi','sage.arguments.txt'):
                    atomic(d/name,'fixture')
                write_json(d/'sage.timing.json',{'sage_wall_elapsed_seconds':1})
            with patch('aging.bqr_benchmark.load_workbook_samples',return_value=(rows,{})), \
                 patch('aging.bqr_benchmark.freeze',return_value={}), \
                 patch('aging.bqr_benchmark.fingerprint',return_value='fixture-fingerprint'), \
                 patch('aging.bqr_benchmark.runtime.call_sample',side_effect=fake_sage), \
                 patch('aging.bqr_benchmark.runtime.run'):
                self.assertEqual(run_benchmark(args),0)
                self.assertEqual(run_benchmark(args),0)
            self.assertEqual(len(calls),6)
            self.assertIn('-skip_bqr',calls[0][1]);self.assertNotIn('-skip_bqr',calls[2][1])
            report=load_config(root/'benchmark/bqr_comparison.json')
            self.assertEqual(len(report['samples']),3)
            self.assertFalse(report['cohort_config_prepared'])
            self.assertIsNone(report['selected_bqr_variant'])


if __name__=='__main__':unittest.main()
