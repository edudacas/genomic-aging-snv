import json
import tempfile
import unittest
from pathlib import Path
from aging.common import *
from aging.cohort import load_truthset,prepare
from aging.contracts import bed_intervals,selected_samples,validate_gate
from aging.evidence import classify,callable_threshold,pileup_site,union
from aging.cli import task
from aging.calibration import resample,split_half,summarize

class LogicTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)
        self.tier=dict(min_callers=2,min_baseq=20,min_mapq=30,min_alt_each_strand=1,
            min_alt_each_orientation=1,max_population_af=.001,germline_alpha=.01)
        self.ev=dict(dp=200,ad_alt=10,vaf=.05,evidence_status='OBSERVED',carrier_variant_exact=0,
            caller_support=2,alt_mean_bq=35,alt_mean_mq=60,alt_forward=5,alt_reverse=5,
            alt_f1r2=5,alt_f2r1=5,germline_overlap=0,germline_binomial_p=1e-20,
            in_pon=0,driver_whitelist=0,gnomad_af='',abraom_af='')
    def tearDown(self):self.tmp.cleanup()
    def test_pass(self):self.assertEqual(classify(self.ev,self.tier,.03,5),'PASS')
    def test_pon_retained_as_flag(self):
        self.ev['in_pon']=1
        self.assertIn('PON',classify(self.ev,self.tier,.03,5))
        self.assertEqual(classify(self.ev,self.tier,.03,5,True),'PASS')
    def test_whitelist_exception(self):
        self.ev.update(in_pon=1,driver_whitelist=1)
        self.assertEqual(classify(self.ev,self.tier,.03,5),'PASS')
    def test_germline_overlap_is_not_blacklist(self):
        self.ev['germline_overlap']=1
        self.assertEqual(classify(self.ev,self.tier,.03,5),'PASS')
        self.ev['germline_binomial_p']=.4
        self.assertIn('GERMLINE_COMPATIBLE',classify(self.ev,self.tier,.03,5))
    def test_carrier_exact_excluded(self):
        self.ev['carrier_variant_exact']=1
        self.assertIn('CARRIER_ALLELE',classify(self.ev,self.tier,.03,5))
    def test_missing_depth_not_negative(self):
        self.ev.update(dp=0,ad_alt=0,vaf='',evidence_status='NO_COVERAGE')
        self.assertIn('NOT_CALLABLE',classify(self.ev,self.tier,.03,5))
    def test_saturation_rejected(self):
        self.ev['evidence_status']='DEPTH_SATURATED'
        self.assertIn('DEPTH_SATURATED',classify(self.ev,self.tier,.03,5))
    def test_population_filter(self):
        self.ev['abraom_af']=.01
        self.assertIn('COMMON_POPULATION',classify(self.ev,self.tier,.03,5))
    def test_callable_threshold(self):
        from scipy.stats import binom
        n=callable_threshold(.03,5,.95)
        self.assertGreaterEqual(binom.sf(4,n,.03),.95)
        self.assertLess(binom.sf(4,n-1,.03),.95)
    def test_zero_depth_bed(self):
        p=self.p/'x.bed';p.write_text('chr1\t0\t10\nchr1\t10\t15\n')
        self.assertEqual(sum(e-s for c,s,e in bed_intervals(p)),15)
    def test_overlapping_bed_rejected(self):
        p=self.p/'x.bed';p.write_text('chr1\t0\t10\nchr1\t9\t15\n')
        with self.assertRaises(ValueError):bed_intervals(p)
    def test_empty_bed_rejected(self):
        p=self.p/'x.bed';p.write_text('')
        with self.assertRaises(ValueError):bed_intervals(p)
    def test_resumption_and_corruption(self):
        d=self.p/'task';task(d,'fp',lambda d:atomic(d/'result','ok'))
        task(d,'fp',lambda d:self.fail('resumed task should not execute'))
        with self.assertRaises(ValueError):task(d,'changed',lambda d:None)
        (d/'result').write_text('corrupt')
        with self.assertRaises(ValueError):task(d,'fp',lambda d:None)
    def test_failed_task_archived(self):
        d=self.p/'task';d.mkdir();(d/'partial').write_text('old')
        task(d,'fp',lambda d:atomic(d/'result','new'))
        self.assertEqual(len(list(self.p.glob('task.failed.*'))),1)
        self.assertFalse((d/'partial').exists())
    def test_draft_gate_rejected(self):
        p=self.p/'gate.json';write_json(p,dict(status='DRAFT'))
        with self.assertRaises(ValueError):validate_gate(p,p,p,p)
    def test_empty_calibration_rejected(self):
        p=self.p/'m.tsv';p.write_text('component\tassay_id\tstatus\tpath\tsha256\tmethod\treviewer\n')
        with self.assertRaises(ValueError):summarize(p,self.p/'cal.json')
    def test_hold_blocks_pilot(self):
        p=self.p/'s.tsv';write_tsv(p,[dict(sample_id='GA_1',analysis_status='HOLD',person_id='p')])
        ids=self.p/'ids';ids.write_text('GA_1\n')
        with self.assertRaises(ValueError):selected_samples(p,'pilot',ids)
    def test_private_permissions(self):
        p=self.p/'private';atomic(p,'private')
        self.assertEqual(p.stat().st_mode & 0o777,0o600)

class WorkbookTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)
    def tearDown(self):self.tmp.cleanup()
    def workbook(self,rows):
        from openpyxl import Workbook
        p=self.p/'t.xlsx';w=Workbook();s=w.active;s.title='truthset_wes_v2'
        s.append(list(rows[0]));
        for r in rows:s.append(list(r.values()))
        w.save(p);return p
    def truth(self):
        return dict(cohort_person_id='person1',franklin_analysis_id='analysis1',franklin_run_id='batch1',
            assay_id='assay1',fap='123456789012',study_group='LFS_R337H',genetic_eligibility_status='PENDING_GERMLINE_VAF',
            batch_matched_design_member='YES',variant_id='chr17:7670699:C:T',age_years='30')
    def test_duplicate_person_rejected(self):
        t=self.truth();p=self.workbook([t,dict(t,franklin_analysis_id='analysis2')])
        with self.assertRaises(ValueError):load_truthset(p)
    def test_formula_rejected(self):
        t=self.truth();t['age_years']='=1+1'
        with self.assertRaises(ValueError):load_truthset(self.workbook([t]))
    def test_exact_link_and_pending_preserved(self):
        t=self.truth();reg=self.p/'r.tsv'
        write_tsv(reg,[dict(analysis_id='analysis1',batch_id='batch1',assay_id='assay1',fap='123456789012',bam_path='/not-present.bam',bai_path='/not-present.bai')])
        result=prepare(self.workbook([t]),reg,self.p/'out')
        self.assertEqual(result['linked'],1);self.assertEqual(result['metadata_ready'],0)
        self.assertEqual(result['hold_reasons']['GENETIC_PENDING'],1)
        self.assertFalse(result['server_accessed'])
    def test_no_prefix_fallback(self):
        t=self.truth();reg=self.p/'r.tsv';write_tsv(reg,[dict(analysis_id='different',fap='123456789012')])
        result=prepare(self.workbook([t]),reg,self.p/'out')
        self.assertEqual(result['linked'],0)
    def test_registry_conflict_holds(self):
        t=self.truth();reg=self.p/'r.tsv'
        write_tsv(reg,[dict(analysis_id='analysis1',batch_id='wrong',assay_id='assay1',fap=t['fap'])])
        result=prepare(self.workbook([t]),reg,self.p/'out')
        self.assertEqual(result['hold_reasons']['CONFLICT_BATCH_ID'],1)

class BAMTests(unittest.TestCase):
    def setUp(self):
        import pysam
        self.ps=pysam;self.tmp=tempfile.TemporaryDirectory();self.p=Path(self.tmp.name)
        self.fa=self.p/'ref.fa';self.fa.write_text('>chr1\n'+'A'*200+'\n');pysam.faidx(str(self.fa))
        self.bam=self.p/'input.bam'
        header=dict(HD=dict(VN='1.6',SO='coordinate'),SQ=[dict(SN='chr1',LN=200)],RG=[dict(ID='r',SM='synthetic')])
        with pysam.AlignmentFile(self.bam,'wb',header=header) as out:
            for i in range(20):
                r=pysam.AlignedSegment();r.query_name='frag'+str(i);r.query_sequence='A'*19+('T' if i<4 else 'A')+'A'*30
                r.flag=0 if i%2 else 16;r.reference_id=0;r.reference_start=0;r.mapping_quality=60;r.cigar=[(0,50)]
                r.query_qualities=pysam.qualitystring_to_array('I'*50);r.set_tag('RG','r');out.write(r)
        pysam.index(str(self.bam))
    def tearDown(self):self.tmp.cleanup()
    def evidence(self,pos=20,ref='A',alt='T'):
        with self.ps.AlignmentFile(self.bam,'rb') as b,self.ps.FastaFile(self.fa) as f:
            return pileup_site(b,f,dict(chrom='chr1',pos=pos,ref=ref,alt=alt),20,30,100000)
    def test_synthetic_snv_evidence(self):
        e=self.evidence();self.assertEqual(e['dp'],20);self.assertEqual(e['ad_alt'],4);self.assertEqual(e['vaf'],.2)
    def test_missing_site(self):
        e=self.evidence(150);self.assertEqual(e['dp'],0);self.assertEqual(e['vaf'],'')
    def test_reference_mismatch(self):
        with self.assertRaises(ValueError):self.evidence(ref='G')
    def test_resample_is_deterministic(self):
        resample(str(self.bam),self.p/'a.bam',.5,17);resample(str(self.bam),self.p/'b.bam',.5,17)
        with self.ps.AlignmentFile(self.p/'a.bam') as a,self.ps.AlignmentFile(self.p/'b.bam') as b:
            self.assertEqual([r.query_name for r in a],[r.query_name for r in b])
    def test_split_disjoint_and_complete(self):
        split_half(str(self.bam),self.p/'split',17);names=[]
        for i in [0,1]:
            with self.ps.AlignmentFile(self.p/'split'/f'half{i}.bam') as b:names.append({r.query_name for r in b})
        self.assertFalse(names[0]&names[1]);self.assertEqual(len(names[0]|names[1]),20)
    def test_exact_indel_not_any_deletion(self):
        from aging.evidence import allele_from_pileup
        from types import SimpleNamespace
        r=SimpleNamespace(query_sequence='AAA')
        p=SimpleNamespace(query_position=0,is_refskip=False,is_del=False,alignment=r,indel=-2)
        self.assertEqual(allele_from_pileup(p,'AAA','A'),'A')
        self.assertEqual(allele_from_pileup(p,'AA','A'),'OTHER')
    def test_union_forced_driver_without_calls(self):
        p=self.p/'bed';p.write_text('chr1\t0\t100\n')
        w=self.p/'drivers';write_tsv(w,[dict(chrom='chr1',pos=20,ref='A',alt='T')])
        out=self.p/'union';rows=union({},out,p,w)
        self.assertEqual(len(rows),1);self.assertEqual(rows[0]['discovery_sample_count'],0)

class ModelTests(unittest.TestCase):
    def test_primary_model_and_nonestimable_outcome(self):
        import numpy as np
        from aging.analysis import model
        rng=np.random.default_rng(42)
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);rows=[]
            for st in range(10):
                for i in range(30):
                    age=20+2*i;case=i in [8,12,16,20]
                    rows.append(dict(group='LFS_R337H' if case else 'TP53_negative',endpoint_status='EVALUABLE',
                        age_at_collection=age,PC1=rng.normal(),PC2=rng.normal(),PC3=rng.normal(),
                        callable_mb=20,somatic_count=int(rng.poisson(4*np.exp(.1*age/10+.15*case*age/10))),
                        retain_pon_count=8,treatment_status='NAIVE',sex='M' if rng.random()>.5 else 'F',
                        cancer_history='YES' if rng.random()>.5 else 'NO',stratum='s'+str(st)))
            write_tsv(p/'features.tsv',rows)
            result=model(p/'features.tsv',dict(min_controls_per_stratum=5,min_model_strata=5),p/'out')
            self.assertEqual(result[0]['status'],'ESTIMATED_EXPLORATORY_UNTIL_VALIDATED')
            self.assertEqual(result[0]['n'],300)
            self.assertGreater(result[0]['rate_ratio_per_decade'],0)
            for r in rows:r['somatic_count']=0
            write_tsv(p/'zero.tsv',rows)
            result=model(p/'zero.tsv',dict(min_controls_per_stratum=5,min_model_strata=5),p/'zero')
            self.assertEqual(result[0]['status'],'NOT_ESTIMABLE')

class ExtensionTests(unittest.TestCase):
    def test_interval_offsets(self):
        from aging.extensions import coverage_index,covered_bp
        x=coverage_index([('chr1',10,20),('chr1',30,50)])
        self.assertEqual(covered_bp(x,'chr1',15,35),10)
        self.assertEqual(covered_bp(x,'chr1',20,30),0)
        self.assertEqual(covered_bp(x,'chr1',0,100),30)
        self.assertEqual(covered_bp(x,'chr2',0,100),0)
    def test_bh_adjustment(self):
        from aging.extensions import bh
        self.assertEqual(bh([.01,.02,.9]),[.03,.03,.9])
    def test_signature_fit_with_synthetic_basis(self):
        from aging.extensions import signature_fit
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);mat=p/'matrix.tsv'
            write_tsv(mat,[dict(context='C'+str(i),S1=1/48 if i<48 else 0,S2=1/48 if i>=48 else 0) for i in range(96)])
            plan=p/'plan.json';write_json(plan,dict(status='FROZEN',reviewer='synthetic_test',
                resources=dict(signatures=dict(path=str(mat),sha256=sha(mat))),min_snvs=100,
                seed=17,bootstrap_reps=10,artifact_signatures=['S2']))
            spec=p/'spec.json';write_json(spec,dict(snvs=480,sbs96_counts={'C'+str(i):10 for i in range(48)}))
            signature_fit(spec,plan,p/'out.json');r=load_config(p/'out.json')
            self.assertAlmostEqual(r['model_artifact_fraction'],0)
            self.assertAlmostEqual(r['exposures'][0]['fraction'],1)
    def test_mosaic_missing_evidence_not_negative(self):
        from aging.extensions import mosaic_candidates
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);bed=p/'tp53.bed';bed.write_text('chr1\t10\t20\n')
            plan=p/'p.json';write_json(plan,dict(status='FROZEN',reviewer='synthetic_test',
                resources=dict(tp53_bed=dict(path=str(bed),sha256=sha(bed))),min_probes=2,min_het_sites=2,
                min_baf_depth=20,min_segment_bp=5,log2_loss=-.2,log2_gain=.2,baf_deviation=.1))
            seg=p/'seg.tsv';write_tsv(seg,[dict(chromosome='chr1',start=0,end=100,log2=-.3,probes=10)])
            baf=p/'baf.tsv';write_tsv(baf,[],['chrom','pos','dp','vaf','evidence_status'])
            mosaic_candidates(seg,baf,plan,p/'o.tsv')
            self.assertEqual(read_tsv(p/'o.tsv')[0]['screen'],'NOT_EVALUABLE')

if __name__=='__main__':unittest.main()
