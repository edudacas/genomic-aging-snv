"""Exactly-three-sample evidence pilot; no cohort launch or somatic classification."""
import argparse,csv,gzip,hashlib,json,os,subprocess,sys,time
from collections import Counter
from pathlib import Path
from fragment_recount import recount_site
from legacy_recount_snapshot import pileup_site

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()

def read(p):return json.loads(Path(p).read_text())

def finalize(out,review=None):
    report=read(out/'PILOT_VALIDATION_PENDING.PRIVATE.json')
    if len(report['pilot_samples'])!=3 or len(set(report['pilot_samples']))!=3 or report['candidate_alleles']<1:
        raise ValueError('Invalid pilot sample/site coverage')
    if report['evaluated_pairs']!=3*report['candidate_alleles'] or sum(report['evaluation_states'].values())!=report['evaluated_pairs']:
        raise ValueError('Incomplete pair coverage')
    for filename,digest in report['outputs'].items():
        if sha(out/filename)!=digest:raise ValueError('Changed pilot output')
    if report['synthetic_tests_passed'] is not True or report['native_tests_passed'] is not True:raise ValueError('Tests incomplete')
    if report['failed_evaluation_pairs']:raise ValueError('Non-evaluable candidate pairs must be resolved before final validation')
    approvals=read(review) if review else {}
    if report['discordances']:
        if approvals.get('report_sha256')!=sha(out/'PILOT_VALIDATION_PENDING.PRIVATE.json'):raise ValueError('Review report identity mismatch')
        if not approvals.get('reviewer'):raise ValueError('Named technical reviewer required')
        with gzip.open(out/'legacy_discordances.PRIVATE.jsonl.gz','rt') as f:
            for line in f:
                row=json.loads(line);entry=approvals.get('discordance_reviews',{}).get(row['discordance_id'],{})
                if entry.get('accepted') is not True or not entry.get('explanation','').strip():raise ValueError('Unexplained discordance')
    receipt=dict(report,status='RECOUNT_PILOT_COMPLETE',validation_report_sha256=sha(out/'PILOT_VALIDATION_PENDING.PRIVATE.json'),
        review_sha256=sha(review) if review else None,cohort_execution_authorized=False,somatic_origin_classified=False)
    (out/'RECOUNT_PILOT_COMPLETE.PRIVATE.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'cohort_execution_authorized':False}))

def main():
    os.umask(0o077);p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prepared',type=Path);p.add_argument('--fasta',type=Path);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--min-bq',type=int,default=30);p.add_argument('--min-mq',type=int,default=30);p.add_argument('--raw-max-depth',type=int,default=100000)
    p.add_argument('--execute-three-pilot-bams',action='store_true');p.add_argument('--finalize-with-review',type=Path)
    a=p.parse_args();out=a.out.resolve()
    if a.finalize_with_review:return finalize(out,a.finalize_with_review)
    if not a.execute_three_pilot_bams:raise ValueError('Explicit --execute-three-pilot-bams required')
    if not a.prepared or not a.fasta:raise ValueError('Prepared pilot and FASTA required')
    if out.exists():raise ValueError('Refuse overwrite')
    prepared=a.prepared.resolve();prep=read(prepared/'PREPARATION_RECEIPT.PRIVATE.json')
    if prep['status']!='PILOT_CANDIDATE_INPUTS_PREPARED_NOT_RECOUNTED' or prep['recount_performed'] is not False:raise ValueError('Invalid preparation receipt')
    for name,digest in prep['outputs'].items():
        path=prepared/name
        if Path(name).is_absolute() or '..' in Path(name).parts or prepared not in path.resolve().parents:raise ValueError('Unsafe prepared path')
        if sha(path)!=digest:raise ValueError('Preparation hash mismatch')
    with (prepared/'pilot_samples.PRIVATE.tsv').open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f,delimiter='\t'))
    ids={r['sample_id'] for r in rows}
    if len(rows)!=3 or len(ids)!=3 or ids!=set(prep['pilot_samples']):raise ValueError('Exactly three prepared samples required')
    # Only specimen/file identity columns are accessed. No age, genotype, treatment or outcome is used.
    sites=[]
    with (prepared/'candidate_union.pilot.vcf').open() as f:
        for line in f:
            if line.startswith('#'):continue
            c,pos,_,ref,alt=line.split('\t')[:5];sites.append(dict(chrom=c,pos=int(pos),ref=ref,alt=alt,variant_id=f'{c}:{pos}:{ref}:{alt}'))
    if not sites or len(sites)!=prep['candidate_snv_alleles'] or len({s['variant_id'] for s in sites})!=len(sites):raise ValueError('Empty or mismatched candidate count')
    import pysam,unittest
    suite=unittest.defaultTestLoader.discover(str(Path(__file__).parent),pattern='test_fragment_recount*.py')
    import io
    testlog=io.StringIO();tests=unittest.TextTestRunner(stream=testlog,verbosity=2).run(suite)
    if not tests.wasSuccessful() or tests.skipped:raise ValueError('All synthetic/native tests must pass without skips')
    for row in rows:
        for field in ('bam_path','bai_path'):
            if not Path(row[field]).is_file():raise ValueError('Pilot BAM/index missing')
    protected=[prepared]+[Path(row[f]).resolve() for row in rows for f in ('bam_path','bai_path')]+[a.fasta.resolve()]
    if any(out==x or x in out.parents or out in x.parents for x in protected):raise ValueError('Output overlaps immutable inputs')
    out.mkdir(parents=True);start=time.monotonic();(out/'synthetic_tests.log').write_text(testlog.getvalue())
    inputs={'prepared_receipt_sha256':sha(prepared/'PREPARATION_RECEIPT.PRIVATE.json'),'fasta_sha256':sha(a.fasta),
            'fasta_index_sha256':sha(str(a.fasta)+'.fai'),'bams':{}}
    states=Counter();discordances=0;scanned=0
    with pysam.FastaFile(str(a.fasta)) as fasta,gzip.open(out/'evidence.PRIVATE.jsonl.gz','wt') as evidence,gzip.open(out/'legacy_discordances.PRIVATE.jsonl.gz','wt') as diffs:
        for row in rows:
            sid=row['sample_id'];inputs['bams'][sid]={'bam_sha256':sha(row['bam_path']),'bai_sha256':sha(row['bai_path'])}
            with pysam.AlignmentFile(row['bam_path'],'rb',index_filename=row['bai_path']) as bam:
                if not bam.has_index():raise ValueError('Indexed BAM required')
                bam_lengths=dict(zip(bam.references,bam.lengths))
                for chrom in {s['chrom'] for s in sites}:
                    if chrom not in bam_lengths or chrom not in fasta.references or bam_lengths[chrom]!=fasta.get_reference_length(chrom):
                        raise ValueError('BAM/FASTA contig identity or length mismatch')
                for site in sites:
                    ev=recount_site(bam,fasta,site,a.min_bq,a.min_mq,a.raw_max_depth);state=ev['evaluation_status'];states[state]+=1;scanned+=1
                    evidence.write(json.dumps(dict(sample_id=sid,**site,**ev),sort_keys=True)+'\n')
                    if state in ('EVALUATED','NO_COVERAGE','NO_USABLE_FRAGMENTS'):
                        old=pileup_site(bam,fasta,site,a.min_bq,a.min_mq,a.raw_max_depth)
                        mapping={'dp':'usable_fragment_depth','ad_ref':'ref_fragments','ad_alt':'alt_fragments',
                                 'alt_forward':'alt_forward','alt_reverse':'alt_reverse','alt_f1r2':'alt_f1r2','alt_f2r1':'alt_f2r1'}
                        changed={k:{'legacy':old.get(k),'new':ev[new]} for k,new in mapping.items() if old.get(k)!=ev[new]}
                        if changed:
                            record=dict(sample_id=sid,variant_id=site['variant_id'],changed=changed,
                                diagnostics=ev['diagnostics'],explanation_status='TECHNICAL_REVIEW_REQUIRED')
                            record['discordance_id']=hashlib.sha256(json.dumps(record,sort_keys=True).encode()).hexdigest()
                            diffs.write(json.dumps(record,sort_keys=True)+'\n');discordances+=1
                    if scanned%1000==0:print(json.dumps({'stage':'FRAGMENT_RECOUNT_PILOT','completed_pairs':scanned,'total_pairs':3*len(sites)}),flush=True)
    assert scanned==3*len(sites)
    report=dict(status='PILOT_EVIDENCE_COMPLETE_VALIDATION_PENDING',pilot_samples=sorted(ids),candidate_alleles=len(sites),
        evaluated_pairs=scanned,evaluation_states=dict(states),discordances=discordances,
        failed_evaluation_pairs=sum(n for s,n in states.items() if s not in ('EVALUATED','NO_COVERAGE','NO_USABLE_FRAGMENTS')),
        synthetic_tests_passed=True,native_tests_passed=True,tests_run=tests.testsRun,
        settings={'min_bq':a.min_bq,'min_mq':a.min_mq,'raw_max_depth':a.raw_max_depth,'scientifically_frozen':False},
        inputs=inputs,pysam_version=pysam.__version__,elapsed_seconds=time.monotonic()-start,
        code_sha256={name:sha(Path(__file__).parent/name) for name in ('fragment_recount.py','run_pilot.py','legacy_recount_snapshot.py','test_fragment_recount.py','test_fragment_recount_native.py')},
        outputs={name:sha(out/name) for name in ('evidence.PRIVATE.jsonl.gz','legacy_discordances.PRIVATE.jsonl.gz','synthetic_tests.log')},
        cohort_execution_authorized=False,somatic_origin_classified=False,endpoint_applied=False,
        runtime='Pointwise raw pileup; explicit non-evaluable saturation, no fallback',
        limitation='Pointwise pilot may be slow on network storage; not a cohort execution engine')
    (out/'PILOT_VALIDATION_PENDING.PRIVATE.json').write_text(json.dumps(report,indent=2)+'\n')
    if not report['failed_evaluation_pairs'] and not discordances:finalize(out)
    else:print(json.dumps({'status':report['status'],'discordances':discordances,'failed_evaluation_pairs':report['failed_evaluation_pairs']}),flush=True)

if __name__=='__main__':main()
