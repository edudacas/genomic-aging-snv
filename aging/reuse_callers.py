"""Reintegrate locally or recount existing calls on the server, never re-run callers."""
import argparse,csv,fcntl,hashlib,json,math,os,time
from pathlib import Path
import pysam
from .common import read_tsv,write_json,sha
from .caller_records import load_records
from .optimized import (integrate,record_sidecar,atomic_tsv,FIELDS,CALLERS,
                        EV_FIELDS,run_evidence,run_single_pass_evidence)
from .joined_vcf import export_joined

def scope_artifact(samples,paths,candidates,outputs,c,out):
    from .scopes import assign_scopes
    discovery={key:set() for key in candidates}
    for sample in samples:
        for caller in CALLERS:
            for key in load_records(paths[sample['sample_id']][caller],caller):
                if key in discovery:discovery[key].add(sample['sample_id'])
    def resource_keys(name):
        resource=c.get('resources',{}).get(name)
        if not resource:return set()
        path=Path(resource['path']);keys=set()
        if path.name.endswith(('.vcf','.vcf.gz')):
            with pysam.VariantFile(str(path)) as vcf:
                for record in vcf:
                    for alt in record.alts or ():
                        if not alt.startswith('<') and alt!='*':keys.add(f'{record.chrom}:{record.pos}:{record.ref}:{alt}')
        else:
            for row in read_tsv(path):keys.add('{chrom}:{pos}:{ref}:{alt}'.format(**row))
        return keys & set(candidates)
    evidence=[]
    for path in outputs:evidence.extend(read_tsv(path))
    forced={name:resource_keys(name) for name in ('driver_whitelist','sage_hotspots','calibration_truth')}
    config=c.get('recount_scope')
    if not config:raise ValueError('v3 recount_scope config missing')
    rows=assign_scopes(candidates,evidence,discovery,forced,config,
        {sample['sample_id']:sample['assay_id'] for sample in samples},require_frozen=False)
    observed=';'.join(sorted(sample['sample_id'] for sample in samples))
    for row in rows:
        row['observed_recount_sample_ids']=observed
        row['scope_execution_mode']='PILOT_ALL_UNION_VALIDATION_OVERRIDE'
    fields=list(rows[0]) if rows else ['variant_id','recount_scope','recount_scope_reason']
    path=out/'site_scopes.tsv';atomic_tsv(path,rows,fields)
    counts={}
    for row in rows:counts[row['recount_scope']]=counts.get(row['recount_scope'],0)+1
    receipt=dict(file=path.name,sha256=sha(path),rows=len(rows),counts=counts,
        selection_uses_caller_pass=False,selection_uses_population_af=False,selection_uses_pon=False,
        cohort_gate_frozen=config.get('frozen') is True,
        scope_applied_to_reduce_pilot_matrix=False,
        scope_execution_mode='PILOT_ALL_UNION_VALIDATION_OVERRIDE')
    write_json(out/'site_scopes.manifest.json',receipt)
    return receipt

def mnv_artifact(samples,paths,candidates,c,out):
    from .mnv_groups import build_mnv_groups
    rows=[]
    for sample in samples:
        sage=load_records(paths[sample['sample_id']]['sage'],'sage')
        rows.extend(build_mnv_groups(sample['sample_id'],sage,candidates))
    fields=(list(rows[0]) if rows else
            ['sample_id','mnv_variant_id','component_variant_ids','representation_status','event_group_id'])
    path=out/'mnv_representation_groups.tsv';atomic_tsv(path,rows,fields)
    counts={}
    for row in rows:counts[row['representation_status']]=counts.get(row['representation_status'],0)+1
    receipt=dict(file=path.name,sha256=sha(path),rows=len(rows),counts=counts,
        sage_version=c.get('sage_version','UNRECORDED'),components_and_mnvs_retained_in_union=True,
        burden_deduplication_applied=False)
    write_json(out/'mnv_representation_groups.manifest.json',receipt)
    return receipt

def compare(source,outputs,out):
    results={};any_mismatch=False
    for path in outputs:
        sid=path.parent.name;old={r['variant_id']:r for r in read_tsv(source/'evidence'/sid/'evidence.tsv')};new=read_tsv(path)
        mismatch=[];changed=0
        if len(old)!=len(new) or set(old)!={r['variant_id'] for r in new}:raise ValueError('Comparison candidate sets differ')
        for r in new:
            baseline=old[r['variant_id']]
            if baseline['evidence_status']=='UNSUPPORTED_COMPLEX_ALLELE':
                if r['evidence_status']!=baseline['evidence_status'] or any(r.get(k) not in ('',None) for k in ('dp','ad_ref','ad_alt','vaf')):mismatch.append((r['variant_id'],'unsupported_semantics'))
                changed+=1;continue
            for k in EV_FIELDS:
                x,y=baseline.get(k,''),r.get(k,'')
                equal=x==y
                if not equal and x and y:
                    try:equal=math.isclose(float(x),float(y),rel_tol=1e-12,abs_tol=1e-12)
                    except ValueError:pass
                if not equal:mismatch.append((r['variant_id'],k))
        any_mismatch |= bool(mismatch)
        results[sid]=dict(rows=len(new),supported_field_mismatches=len(mismatch),first_mismatches=mismatch[:20],unsupported_zero_to_null_rows=changed)
    report=dict(status='MISMATCH' if any_mismatch else 'MATCH_SUPPORTED_EVIDENCE',samples=results)
    write_json(out/'baseline_comparison.json',report)
    return not any_mismatch

def execute(a):
    run_start=time.monotonic()
    source=Path(a.source).resolve();out=Path(a.out).resolve()
    if source==out or source in out.parents:raise ValueError('Output must be separate from baseline')
    state=json.loads((source/'status.json').read_text())
    if state['status']!='CALLERS_AND_EVIDENCE_COMPLETE' or state['failed_samples']:raise ValueError('Baseline not complete')
    samples=[r for r in read_tsv(source/'samples.private.tsv') if r.get('explore_status')=='CALLERS_COMPLETE']
    if len(samples)!=state['selected_samples']:raise ValueError('Baseline sample count mismatch')
    paths={};source_hashes={}
    for sample in samples:
        sid=sample['sample_id'];paths[sid]={};root=source/'calls'/sid
        done=json.loads((root/'COMPLETE.json').read_text())
        for caller in CALLERS:
            name=caller+'.norm.vcf.gz'
            for suffix in ('','.tbi'):
                file=root/(name+suffix);digest=sha(file)
                if done['outputs'].get(name+suffix)!=digest:raise ValueError('Baseline caller receipt mismatch: '+str(file))
                source_hashes[str(file)]=digest
            paths[sid][caller]=str(root/name)
    candidates={r['variant_id']:{**r,'pos':int(r['pos'])} for r in read_tsv(source/'candidate_union.tsv')}
    if len(candidates)!=state['evidence_rows']//len(samples):raise ValueError('Baseline candidate count mismatch')
    c=json.loads(Path(a.config).read_text()) if a.config else {}
    extra={str(source/'candidate_union.tsv'):sha(source/'candidate_union.tsv'),str(source/'samples.private.tsv'):sha(source/'samples.private.tsv')}
    for sample in samples:
        f=source/'evidence'/sample['sample_id']/'evidence.tsv';extra[str(f)]=sha(f)
    code={f.name:sha(f) for f in Path(__file__).parent.glob('*.py')}
    resource_hashes={}
    if a.mode in ('recount','single-pass-compare'):
        from .explore import freeze
        from .cli import fingerprint
        resource_hashes=freeze(c,None)
        rows=[dict(r,germline_vcf='') for r in samples]
        bam_fp=fingerprint(rows,c,resource_hashes)
    else:bam_fp=None
    fp=hashlib.sha256(json.dumps(dict(mode=a.mode,source=source_hashes,extra=extra,config=c,code=code,bam_fp=bam_fp),sort_keys=True).encode()).hexdigest()
    out.mkdir(parents=True,exist_ok=True)
    with (out/'.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        marker=out/'RUN_INPUTS.json'
        if marker.exists() and json.loads(marker.read_text())['fingerprint']!=fp:raise ValueError('Inputs/code changed; choose new output directory')
        write_json(marker,dict(fingerprint=fp,mode=a.mode,source_hashes=source_hashes,extra=extra,code=code))
        with pysam.VariantFile(paths[samples[0]['sample_id']]['mutect']) as v:
            contigs=[(k,v.header.contigs[k].length) for k in v.header.contigs]
        rank={k:i for i,(k,n) in enumerate(contigs)}
        sites=sorted(candidates.values(),key=lambda r:(rank[r['chrom']],r['pos'],r['ref'],r['alt']))
        atomic_tsv(out/'candidate_union.tsv',sites,['variant_id','chrom','pos','ref','alt','discovered_by','discovery_sample_count'])
        atomic_tsv(out/'samples.private.tsv',samples,list(samples[0]))
        if a.mode in ('recount','single-pass-compare'):
            if a.workers not in range(1,5):raise ValueError('workers must be 1..4')
            outputs=(run_single_pass_evidence(samples,c,paths,candidates,out,fp,a.workers)
                     if a.mode=='single-pass-compare' else
                     run_evidence(samples,c,paths,candidates,out,fp,a.workers))
        else:
            outputs=[]
            for sample in samples:
                sid=sample['sample_id'];directory=out/'evidence'/sid;directory.mkdir(parents=True,exist_ok=True)
                records={n:load_records(paths[sid][n],n) for n in CALLERS};record_sidecar(directory/'caller_records.tsv.gz',records,sid)
                old={r['variant_id']:r for r in read_tsv(source/'evidence'/sid/'evidence.tsv')}
                def rows():
                    for site in sites:
                        previous=old[site['variant_id']];ev={k:previous.get(k,'') for k in EV_FIELDS}
                        if ev['evidence_status']=='UNSUPPORTED_COMPLEX_ALLELE':ev.update(dp='',ad_ref='',ad_alt='',vaf='')
                        row=integrate(site,ev,sample,records,set(),{},dict(resources={}))
                        for k in ['in_pon','pon_status','gnomad_af','gnomad_resource_status','abraom_af','abraom_resource_status']:row[k]=previous.get(k,'')
                        yield row
                output=directory/'evidence.tsv';atomic_tsv(output,rows(),FIELDS);outputs.append(output)
                print(json.dumps(dict(mode=a.mode,sample_completed=sid)),flush=True)
        from .explore import export_tables
        phase_start=time.monotonic()
        counts=export_tables(outputs,out/'evidence.all.tsv.gz')
        tsv_seconds=time.monotonic()-phase_start;phase_start=time.monotonic()
        joined=export_joined(outputs,out/'joined.evidence.vcf.gz',contigs)
        join_seconds=time.monotonic()-phase_start
        parquet={}
        parquet_seconds=0
        if a.mode=='single-pass-compare':
            from .parquet import export_run
            phase_start=time.monotonic()
            parquet=export_run(out,samples,candidates,outputs)
            parquet_seconds=time.monotonic()-phase_start
            phase_start=time.monotonic()
            scopes=scope_artifact(samples,paths,candidates,outputs,c,out)
            mnv=mnv_artifact(samples,paths,candidates,c,out)
            from .validate_views import validate_views
            views=validate_views(out,outputs)
        else:scopes={};mnv={};views={}
        derived_metadata_seconds=time.monotonic()-phase_start if a.mode=='single-pass-compare' else 0
        phase_start=time.monotonic()
        matched=compare(source,outputs,out)
        compare_seconds=time.monotonic()-phase_start
        timing_path=out/'phase_timings.json'
        timing=json.loads(timing_path.read_text()) if timing_path.exists() else {}
        timing.update(tsv_export_elapsed_seconds=tsv_seconds,join_elapsed_seconds=join_seconds,
            parquet_elapsed_seconds=parquet_seconds,scope_and_mnv_elapsed_seconds=derived_metadata_seconds,
            baseline_comparison_elapsed_seconds=compare_seconds,total_elapsed_seconds=time.monotonic()-run_start)
        if a.mode=='recount':timing['recount_and_exports_elapsed_seconds']=timing['population_elapsed_seconds']+timing['sample_recount_wall_elapsed_seconds']+tsv_seconds+join_seconds
        write_json(timing_path,timing)
        status=('SINGLE_PASS_COMPARISON_PASS' if a.mode=='single-pass-compare' and matched else
                'RECOUNT_COMPARISON_PASS' if a.mode=='recount' and matched else
                'REINTEGRATION_COMPLETE' if matched else 'EVIDENCE_MISMATCH')
        write_json(out/'status.json',dict(status=status,mode=a.mode,selected_samples=len(samples),evidence_rows=counts['evidence_rows'],joined=joined,**parquet,artifact_validation=views,site_scopes=scopes,mnv_representation_groups=mnv,baseline_matched=matched,callers_rerun=False,bam_recount_performed=a.mode in ('recount','single-pass-compare'),somatic_origin_classified=False,cohort_execution_authorized=False))
        if matched and (out/'failure.json').exists():os.replace(out/'failure.json',out/'previous_failure.json')
        return 0 if matched else 2

def main():
    os.umask(0o077)
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',required=True);p.add_argument('--out',required=True);p.add_argument('--config');p.add_argument('--mode',choices=['reintegrate','recount','single-pass-compare'],required=True);p.add_argument('--workers',type=int,default=4)
    a=p.parse_args()
    if a.mode in ('recount','single-pass-compare') and not a.config:p.error('--config required for BAM recount')
    try:return execute(a)
    except Exception as e:
        out=Path(a.out).resolve();source=Path(a.source).resolve()
        if out!=source and source not in out.parents:
            out.mkdir(parents=True,exist_ok=True)
            write_json(out/'failure.json',dict(error=type(e).__name__,message=str(e)))
        raise
if __name__=='__main__':raise SystemExit(main())
