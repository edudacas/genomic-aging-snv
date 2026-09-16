import tempfile
import unittest
from pathlib import Path
from aging.callers import candidate_union


class CohortUnionTest(unittest.TestCase):
    def test_three_callers_resume_and_stop_before_union(self):
        import gzip,json,shutil
        import pysam
        from types import SimpleNamespace
        from unittest.mock import patch
        from test_explore import ExploreIntegrationTests
        from aging.callers import execute
        from aging.common import write_json,write_tsv
        fixture=ExploreIntegrationTests();fixture.setUp()
        self.addCleanup(fixture.tearDown)
        d=fixture.p
        pop=d/'pop.vcf'
        pop.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=200>\n'
                       '##INFO=<ID=AF,Number=A,Type=Float,Description="frequency">\n'
                       '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        gz=str(pop)+'.gz';pysam.tabix_compress(str(pop),gz);pysam.tabix_index(gz,preset='vcf')
        c=fixture.c
        c.update(callers=['mutect','lofreq','sage'],sage_extra_args=['-skip_bqr','-skip_msi_jitter'],
                 pon_by_assay={'TEST':gz}, sample_workers=4, export_parquet=True)
        c['resources'].update(gnomad={'path':gz},contamination_sites={'path':gz},
                              sage_jar={'path':'/opt/sage/sage.jar'})
        rows=[dict(fixture.row,assay_id='TEST',germline_vcf='/missing/deferred.vcf.gz'),
              dict(fixture.row,sample_id='S2',person_id='P2',analysis_id='A2',assay_id='TEST',germline_vcf='/missing/deferred.vcf.gz')]
        write_json(d/'config.json',c);write_tsv(d/'samples.tsv',rows)
        a=SimpleNamespace(command='callers-only',config=str(d/'config.json'),samples=str(d/'samples.tsv'),
                          sample_ids=None,out=str(d/'out'))
        def fake_calls(r,c,d,pon):
            source=fixture.vcf(d)
            for name in ('lofreq','sage'):
                shutil.copyfile(source,d/(name+'.norm.vcf.gz'))
        with patch('aging.callers.freeze',return_value={'synthetic_resource_validation':True}), \
             patch('aging.runtime.qc',side_effect=lambda r,c,d:write_json(d/'qc.json',{'synthetic':True})), \
             patch('aging.runtime.call_sample',side_effect=fake_calls) as calls, \
             patch('aging.callers.candidate_union',side_effect=AssertionError('Forbidden union')), \
             patch('aging.optimized.run_single_pass_evidence',side_effect=AssertionError('Forbidden recount')):
            self.assertEqual(execute(a),0)
            self.assertEqual(execute(a),0)
            self.assertEqual(calls.call_count,2)
        status=json.loads((d/'out/status.json').read_text())
        self.assertEqual(status['status'],'CALLERS_COMPLETE')
        self.assertFalse(status['union_performed'])
        self.assertFalse(status['recount_performed'])
        receipt=json.loads((d/'out/CALLERS_COMPLETE.json').read_text())
        self.assertEqual(len(receipt['samples']),2)
        self.assertFalse(receipt['pilot_cache_imported'])
        for name in ('candidate_union.tsv','evidence','joined.evidence.vcf.gz','evidence.all.tsv.gz'):
            self.assertFalse((d/'out'/name).exists())
        a.command='run'
        with self.assertRaisesRegex(ValueError,'future gated invocation'):
            execute(a)
        # A failed caller must not produce a cohort completion receipt; other
        # samples still finish and remain available for same-run resume.
        a.command='callers-only';a.out=str(d/'failed-out')
        def partial_failure(r,c,d,pon):
            if r['sample_id']=='S2':raise RuntimeError('synthetic caller failure')
            return fake_calls(r,c,d,pon)
        with patch('aging.callers.freeze',return_value={'synthetic_resource_validation':True}), \
             patch('aging.runtime.qc',side_effect=lambda r,c,d:write_json(d/'qc.json',{'synthetic':True})), \
             patch('aging.runtime.call_sample',side_effect=partial_failure), \
             patch('aging.callers.candidate_union',side_effect=AssertionError('Forbidden union')):
            self.assertEqual(execute(a),2)
        self.assertFalse((d/'failed-out/CALLERS_COMPLETE.json').exists())
        self.assertEqual(json.loads((d/'failed-out/status.json').read_text())['failed_samples'],['S2'])

    def test_foreign_candidate_is_recounted_without_origin_classification(self):
        from test_explore import ExploreIntegrationTests
        from aging.explore import observe
        fixture=ExploreIntegrationTests()
        fixture.setUp()
        self.addCleanup(fixture.tearDown)
        fixture.c['label_origins']=False
        fixture.row['germline_vcf_status']='DEFERRED'
        union={'chr1:20:A:T':dict(chrom='chr1',pos=20,ref='A',alt='T'),
               'chr1:150:A:T':dict(chrom='chr1',pos=150,ref='A',alt='T')}
        rows=observe(fixture.row,fixture.c,{},fixture.p/'evidence.tsv',{},set(),union)
        self.assertEqual(len(rows),2)
        by_key={r['variant_id']:r for r in rows}
        self.assertEqual(by_key['chr1:20:A:T']['ad_alt'],4)
        self.assertEqual(by_key['chr1:20:A:T']['caller_support'],0)
        self.assertEqual(by_key['chr1:150:A:T']['evidence_status'],'NO_COVERAGE')
        self.assertEqual(by_key['chr1:150:A:T']['vaf'],'')
        self.assertTrue(all(r['origin_label']=='NOT_CLASSIFIED_CALLER_EVIDENCE_ONLY' for r in rows))

    def test_nonpass_cross_sample_and_territory(self):
        with tempfile.TemporaryDirectory() as tmp:
            d=Path(tmp)
            bed=d/'common.bed';bed.write_text('chr1\t10\t20\n')
            header=('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n'
                    '##FILTER=<ID=artifact,Description="test">\n'
                    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
            a=d/'a.vcf';a.write_text(header+'chr1\t11\t.\tA\tT\t30\tartifact\t.\n')
            b=d/'b.vcf';b.write_text(header+'chr1\t20\t.\tG\tC\t30\tPASS\t.\n'
                                   +'chr1\t21\t.\tG\tT\t30\tPASS\t.\n')
            paths={'A':dict.fromkeys(['mutect','lofreq','sage'],str(a)),
                   'B':dict.fromkeys(['mutect','lofreq','sage'],str(b))}
            result=candidate_union(paths,str(bed))
            self.assertEqual(set(result),{'chr1:11:A:T','chr1:20:G:C'})


if __name__=='__main__':unittest.main()
