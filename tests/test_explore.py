import csv
import gzip
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
from aging.common import write_tsv,write_json,read_tsv
from aging.explore import selected,labels,load_ibd,observe,explore


class OriginTests(unittest.TestCase):
    def row(self,dp=200,alt=10):
        return dict(sample_id='S',variant_id='chr17:10:A:T',dp=dp,ad_alt=alt,alt_forward=5,alt_reverse=5)
    def test_high_vaf_is_inconclusive_even_above_half(self):
        for alt in [92,100,160,200]:
            self.assertEqual(labels(self.row(200,alt))['origin_label'],'INCONCLUSIVE')
    def test_low_vaf_and_both_reference_means(self):
        r=labels(self.row())
        self.assertEqual(r['origin_label'],'LIKELY_SOMATIC')
        self.assertLess(r['binomial_p_050'],r['binomial_p_046'])
        self.assertEqual(r['classification_scope'],'EXPLORATORY_EXTRAPOLATION')
    def test_shared_haplotype_and_vaf_both_required(self):
        ibd={('S','chr17:10:A:T'):[dict(shared_cm='1.5',source='independent_phasing')]}
        self.assertEqual(labels(self.row(200,100),ibd)['origin_label'],'LIKELY_GERMLINE')
        self.assertEqual(labels(self.row(),ibd)['origin_label'],'LIKELY_SOMATIC')
        ibd[('S','chr17:10:A:T')][0]['shared_cm']='1.49'
        self.assertEqual(labels(self.row(200,100),ibd)['origin_label'],'INCONCLUSIVE')
    def test_missing_and_zero_alt_are_not_somatic(self):
        for dp,alt in [(0,0),(200,0),(200,1),(5,3)]:
            self.assertEqual(labels(self.row(dp,alt))['origin_label'],'INSUFFICIENT_SUPPORT')
    def test_original_mutect_counts_have_explicit_provenance(self):
        r=self.row();r.update(mutect_dp=200,mutect_ad_alt=100)
        r=labels(r)
        self.assertEqual(r['origin_counts_source'],'MUTECT_FORMAT_DP_AD')
        self.assertEqual(r['origin_label'],'INCONCLUSIVE')
    def test_clinical_holds_do_not_block_but_identity_conflict_does(self):
        r=dict(sample_id='S',bam_path='/b',bai_path='/i',technical_status='HOLD',technical_reasons='STRATUM_MISSING',
               analysis_status='HOLD',hold_reasons='TREATMENT_UNKNOWN;GENETIC_PENDING;ANCESTRY_MISSING')
        self.assertEqual(selected([r])[0]['explore_status'],'SELECTED_NOT_SERVER_VALIDATED')
        r['technical_reasons']='CONFLICT_FAP'
        self.assertEqual(selected([r])[0]['explore_status'],'TECHNICAL_INPUT_PENDING')
    def test_ibd_rejects_self_sharing_and_nonoverlapping_segment(self):
        with tempfile.TemporaryDirectory() as t:
            p=Path(t)/'i.tsv'
            r=dict(sample_id='S',partner_sample_id='P',variant_id='chr17:10:A:T',build='GRCh38',chrom='chr17',
                   start_1based=1,end_1based=100,shared_cm=2,phased='YES',partner_same_variant='YES',source='independent')
            write_tsv(p,[r]);self.assertTrue(load_ibd(p))
            for change in [dict(partner_sample_id='S'),dict(start_1based=11),dict(partner_same_variant='NO')]:
                write_tsv(p,[dict(r,**change)])
                with self.assertRaises(ValueError):load_ibd(p)


class ExploreIntegrationTests(unittest.TestCase):
    def setUp(self):
        import pysam
        self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)
        self.fa=self.p/'ref.fa';self.fa.write_text('>chr1\n'+'A'*200+'\n');pysam.faidx(str(self.fa))
        self.fa.with_suffix('.dict').write_text('@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:200\n')
        self.bed=self.p/'targets.bed';self.bed.write_text('chr1\t0\t200\n')
        self.bam=self.p/'source.bam'
        header=dict(HD=dict(SO='coordinate'),SQ=[dict(SN='chr1',LN=200)],RG=[dict(ID='r',SM='synthetic')])
        with pysam.AlignmentFile(self.bam,'wb',header=header) as f:
            for i in range(20):
                r=pysam.AlignedSegment();r.query_name='r'+str(i)
                seq=list('A'*50)
                if i<4: seq[19]='T'
                if i<12: seq[29]='T'
                r.query_sequence=''.join(seq);r.flag=16 if i%2 else 0;r.reference_id=0;r.reference_start=0
                r.mapping_quality=60;r.cigar=[(0,50)];r.query_qualities=pysam.qualitystring_to_array('I'*50)
                r.set_tag('RG','r');f.write(r)
        pysam.index(str(self.bam))
        self.row=dict(sample_id='S',person_id='P',analysis_id='A',group='LFS_R337H',genetic_status='PENDING_GERMLINE_VAF',
            technical_status='LINKED_NOT_SERVER_VALIDATED',technical_reasons='',analysis_status='HOLD',
            hold_reasons='TREATMENT_UNKNOWN;GENETIC_PENDING',assay_id='',batch_id='',bam_path=str(self.bam),
            bai_path=str(self.bam)+'.bai',germline_vcf='/absent/optional.vcf',carrier_variant='')
        self.c=dict(build='GRCh38',profile='exploratory',threads=1,java_heap_gb=1,keep_derived_bam=False,
            callers=['mutect'],evidence=dict(min_baseq=20,min_mapq=30,max_depth=100000),assays={},secondary={'telseq':False},
            resources={n:dict(path=str(p),build='GRCh38') for n,p in dict(fasta=self.fa,fai=str(self.fa)+'.fai',
                dict=self.fa.with_suffix('.dict'),common_bed=self.bed).items()})
    def tearDown(self):self.tmp.cleanup()
    def vcf(self,d):
        import pysam
        p=d/'mutect.vcf'
        p.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=200>\n'
          '##FILTER=<ID=artifact,Description="synthetic rejected call">\n'
          '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="depth">\n'
          '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="allele depth">\n'
          '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n'
          'chr1\t20\t.\tA\tT\t.\tartifact\t.\tDP:AD\t20:16,4\n'
          'chr1\t30\t.\tA\tT\t.\tPASS\t.\tDP:AD\t20:8,12\n')
        out=d/'mutect.norm.vcf.gz';pysam.tabix_compress(str(p),str(out),force=True);pysam.tabix_index(str(out),preset='vcf',force=True)
        return str(out)
    def test_no_deletion_of_filtered_or_high_vaf_alleles_and_unknown_not_zero(self):
        self.row['germline_vcf_status']='NOT_LINKED'
        rows=observe(self.row,self.c,{'mutect':self.vcf(self.p)},self.p/'e.tsv',{},set())
        self.assertEqual(len(rows),2)
        self.assertEqual(rows[0]['mutect_filter'],'artifact')
        self.assertEqual(rows[1]['vaf'],.6)
        self.assertEqual(rows[1]['origin_label'],'INCONCLUSIVE')
        self.assertEqual(rows[0]['germline_overlap'],'')
        self.assertEqual(rows[0]['in_pon'],'')
    def test_complete_orchestration_and_resume_with_synthetic_caller_outputs(self):
        samples=self.p/'samples.tsv';write_tsv(samples,[self.row,dict(self.row,sample_id='M',bam_path='',bai_path='',technical_status='HOLD',technical_reasons='REGISTRY_ABSENT')])
        config=self.p/'config.json';write_json(config,self.c)
        a=SimpleNamespace(config=str(config),samples=str(samples),out=str(self.p/'run'),ibd=None,sample_ids=None,preflight_only=False)
        def fake_qc(s,c,d):write_json(d/'qc.json',{'mean_depth':1,'test':'synthetic QC; no minimum-depth exclusion'})
        def fake_calls(s,c,d,pon):self.vcf(d)
        with patch('aging.runtime.qc',side_effect=fake_qc),patch('aging.runtime.call_sample',side_effect=fake_calls) as call:
            self.assertEqual(explore(a),0);self.assertEqual(call.call_count,1)
            self.assertEqual(explore(a),0);self.assertEqual(call.call_count,1)
        out=Path(a.out)
        with gzip.open(out/'variants.all.private.tsv.gz','rt') as f: rows=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(len(rows),2)
        self.assertEqual(rows[0]['germline_vcf_status'],'PATH_UNAVAILABLE')
        self.assertEqual(len(read_tsv(out/'samples.all.private.tsv')),2)
        self.assertEqual(json.loads((out/'explore.status.json').read_text())['status'],'AVAILABLE_SUBSET_COMPLETE')
    def test_mutect_command_retains_germline_and_does_not_require_pon(self):
        from aging.runtime import call_sample
        commands=[]
        row=dict(self.row,germline_vcf='')
        with patch('aging.runtime.run',side_effect=lambda args,*rest:commands.append(list(map(str,args)))),patch('aging.runtime.normalize'):
            call_sample(row,self.c,self.p/'calls',None)
        mutect=next(x for x in commands if 'Mutect2' in x)
        self.assertIn('--genotype-germline-sites',mutect)
        self.assertNotIn('--germline-resource',mutect)
        self.assertFalse(any('GetPileupSummaries' in x for x in commands))
        self.assertTrue(any('FilterMutectCalls' in x for x in commands))
    def test_one_failed_sample_does_not_stop_other_samples(self):
        samples=self.p/'samples.tsv';write_tsv(samples,[dict(self.row,sample_id='FAIL'),self.row])
        config=self.p/'config.json';write_json(config,self.c)
        a=SimpleNamespace(config=str(config),samples=str(samples),out=str(self.p/'run'),ibd=None,sample_ids=None,preflight_only=False)
        def fake_qc(s,c,d):
            if s['sample_id']=='FAIL': raise RuntimeError('synthetic failure')
            write_json(d/'qc.json',{'synthetic':True})
        with patch('aging.runtime.qc',side_effect=fake_qc),patch('aging.runtime.call_sample',side_effect=lambda s,c,d,pon:self.vcf(d)):
            self.assertEqual(explore(a),2)
        out=Path(a.out)
        self.assertEqual([r['explore_status'] for r in read_tsv(out/'samples.all.private.tsv')],['FAILED','EVIDENCE_COMPLETE'])
        self.assertEqual(json.loads((out/'explore.status.json').read_text())['evidence_rows'],2)
