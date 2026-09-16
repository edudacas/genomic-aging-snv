"""Prepared blood BAMs: three independent callers and cohort-wide evidence, no interpretation."""
import argparse
import fcntl
import hashlib
import json
import os
import threading
import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from .common import load_config, read_tsv, write_json, write_tsv, sha
from .cli import fingerprint, task
from .explore import selected, freeze, records, observe, export_tables
from .contracts import territory_index, within_territory
from . import runtime


def candidate_union(paths_by_sample, bed):
    territory=territory_index(bed)
    result={}
    origins={};discoveries={}
    for sample_id,paths in paths_by_sample.items():
        for caller in ('mutect','lofreq','sage'):
            for key,r in records(paths[caller],caller).items():
                if within_territory(territory,r['chrom'],r['pos'],r['ref']):
                    result[key]={k:r[k] for k in ('chrom','pos','ref','alt')}
                    origins.setdefault(key,set()).add(caller)
                    discoveries.setdefault(key,set()).add(sample_id)
    for key,r in result.items():
        r.update(discovered_by=';'.join(sorted(origins[key])),discovery_sample_count=len(discoveries[key]))
    return result


def run_jobs(items, workers, out, phase, fn):
    """Bounded sample parallelism; one atomic progress writer at a time."""
    guard = threading.Lock()
    states = {}
    completed = 0
    def publish():
        snapshot = dict(phase=phase, workers=workers, completed=completed,
                        total=len(items), worker_states=dict(states),
                        updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        write_json(Path(out)/'progress.json', snapshot)
        print(json.dumps(snapshot), flush=True)
    def job(item):
        worker = threading.current_thread().name
        def stage(name):
            with guard:
                states[worker] = dict(sample_id=item['sample_id'], stage=name)
                publish()
        try:
            result = fn(item, stage)
            stage('FAILED' if isinstance(result, tuple) and result[0].get('explore_status')=='FAILED' else 'COMPLETE')
            return result
        except Exception:
            stage('FAILED')
            raise
    with guard:
        publish()
    with ThreadPoolExecutor(max_workers=workers, thread_name_prefix='worker') as pool:
        futures = {pool.submit(job, item): i for i, item in enumerate(items)}
        results = [None] * len(items)
        for future in as_completed(futures):
            results[futures[future]] = future.result()
            with guard:
                completed += 1
                publish()
    return results


def execute(a):
    if a.command not in ('preflight','callers-only'):
        raise ValueError('This release stops at CALLERS_COMPLETE; recount requires a future gated invocation')
    c=load_config(a.config)
    c.update(profile='exploratory',label_origins=False)
    if c.get('callers')!=['mutect','lofreq','sage']:
        raise ValueError('This protocol requires callers [mutect, lofreq, sage]')
    if not 1 <= c['threads'] <= int(os.environ.get('AGING_CPU_LIMIT','8')):
        raise ValueError('Thread budget exceeded')
    workers=c.get('sample_workers',1)
    if type(workers) is not int or not 1 <= workers <= 4:
        raise ValueError('sample_workers must be an integer from 1 to 4')
    for k in ('gnomad','contamination_sites'):
        if k not in c['resources']:
            raise ValueError('Protocol resource missing: '+k)
    if c.get('secondary',{}).get('cnv') or c.get('secondary',{}).get('telseq'):
        raise ValueError('This delivery runs only the SNV/indel callers and evidence')
    if c['resources']['sage_jar']['path']!='/opt/sage/sage.jar':
        raise ValueError('Use the versioned SAGE JAR embedded in this image')
    if c.get('sage_extra_args')!=['-skip_bqr','-skip_msi_jitter']:
        raise ValueError('Frozen SAGE profile requires explicit default models without Redux files')
    input_receipt=None
    if getattr(a,'workbook',None):
        from .workbook_input import load_workbook_samples
        if not a.linkage:raise ValueError('--linkage required with --workbook')
        input_rows,input_receipt=load_workbook_samples(a.workbook,a.linkage,check_files=True)
        c['study_input']=input_receipt
    else:
        input_rows=read_tsv(a.samples)
    rows=selected(input_rows,Path(a.sample_ids).read_text().split() if a.sample_ids else None)
    # Germline inputs are provenance only in this callers-only entry point.
    for r in rows:
        r['germline_vcf_deferred_source'] = r.get('germline_vcf', '')
        r['germline_vcf'] = ''
    chosen=[r for r in rows if r['explore_status']=='SELECTED_NOT_SERVER_VALIDATED']
    if not chosen: raise ValueError('No technically selectable BAMs')
    for field in ('person_id','analysis_id'):
        values=[r.get(field) for r in chosen]
        if not all(values) or len(set(values))!=len(values):
            raise ValueError('Missing or duplicate '+field+' in selected samples')
    for r in chosen:
        if not c.get('pon_by_assay',{}).get(r['assay_id']):
            raise ValueError('Assay PoN missing: '+r['assay_id'])
    out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True)
    with (out/'.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        resources=freeze(c,None)
        fp=fingerprint(chosen,c,resources)
        marker=out/'RUN_INPUTS.json'
        if marker.exists() and load_config(marker)['fingerprint']!=fp:
            raise ValueError('Inputs or code changed; use a new output directory')
        write_json(marker,dict(fingerprint=fp,profile='callers_only'))
        if input_receipt:write_json(out/'workbook_input.json',input_receipt)
        write_json(out/'resources.lock.json',resources)
        paths={}; failures=[]
        def process_sample(original, stage):
            r=dict(original); sid=r['sample_id']
            try:
                stage('VALIDATE_BAM')
                r['bam_sample']=runtime.validate_bam(r,c)
                if a.command=='preflight':
                    r['explore_status']='PREFLIGHT_PASS'
                    return r, None
                sfp=fingerprint([r],c,resources)
                stage('QC')
                task(out/'qc'/sid,sfp,lambda d:runtime.qc(r,c,d))
                sample=dict(r,germline_vcf='',carrier_variant='')
                d=out/'calls'/sid
                stage('CALLERS')
                task(d,sfp,lambda d:runtime.call_sample(sample,c,d,c['pon_by_assay'][r['assay_id']]))
                r['explore_status']='CALLERS_COMPLETE'
                return r, {n:str(d/(n+'.norm.vcf.gz')) for n in c['callers']}
            except Exception as e:
                r['explore_status']='FAILED';r['error']=str(e)
                return r, None
        processed=run_jobs(chosen,workers,out,'PREFLIGHT' if a.command=='preflight' else 'QC_AND_CALLERS',process_sample)
        by_id={r['sample_id']:r for r in rows}
        for r, sample_paths in processed:
            by_id[r['sample_id']].update(r)
            if r['explore_status']=='FAILED': failures.append(r['sample_id'])
            if sample_paths: paths[r['sample_id']]=sample_paths
        write_tsv(out/'samples.private.tsv',rows)
        counts={}
        if a.command=='callers-only' and not failures:
            receipt=dict(status='CALLERS_COMPLETE',selected_samples=len(chosen),fingerprint=fp,
                         union_performed=False,recount_performed=False,pilot_cache_imported=False,
                         config_sha256=sha(a.config),resources_lock_sha256=sha(out/'resources.lock.json'),
                         samples_manifest_sha256=sha(out/'samples.private.tsv'),workbook_input=input_receipt,
                         samples={sid:dict(qc_complete_sha256=sha(out/'qc'/sid/'COMPLETE.json'),
                                           calls_complete_sha256=sha(out/'calls'/sid/'COMPLETE.json'),
                                           normalized_vcfs={name:dict(file=str(Path(p).relative_to(out)),sha256=sha(p))
                                                            for name,p in sample_paths.items()})
                                  for sid,sample_paths in paths.items()})
            write_json(out/'CALLERS_COMPLETE.json',receipt)
            counts['receipt']='CALLERS_COMPLETE.json'
        status='FAILED' if failures else ('PREFLIGHT_PASS' if a.command=='preflight' else 'CALLERS_COMPLETE')
        write_json(out/'status.json',dict(status=status,selected_samples=len(chosen),failed_samples=failures,
                   union_performed=False,recount_performed=False,
                   functional_annotation=False,somatic_origin_classified=False,**counts))
        return 2 if failures else 0


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['preflight','callers-only'])
    for k in ('config','out'):p.add_argument('--'+k,required=True)
    inputs=p.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--samples');inputs.add_argument('--workbook')
    p.add_argument('--linkage')
    p.add_argument('--sample-ids')
    a=p.parse_args()
    try:return execute(a)
    except Exception as e:
        write_json(Path(a.out)/'failure.json',dict(error=type(e).__name__,message=str(e)))
        raise


if __name__=='__main__':
    raise SystemExit(main())
