"""Entry point. Explore generates data without clinical or statistical gates."""
import argparse
import fcntl
import json
import os
import shutil
import sys
from pathlib import Path
from .common import *
from .contracts import freeze_resources, selected_samples, validate_gate

def fingerprint(rows,c,extra=None):
    code={p.name:sha(p) for p in Path(__file__).parent.glob('*.py')}
    files={}
    for r in rows:
        for k in ['bam_path','bai_path','germline_vcf']:
            if r.get(k):
                p=Path(r[k]);st=p.stat()
                files[str(p.resolve())]=[st.st_size,st.st_mtime_ns,r.get('bam_sha256','') if k=='bam_path' else sha(p)]
    return hashlib.sha256(json.dumps(dict(rows=rows,config=c,code=code,files=files,extra=extra),sort_keys=True).encode()).hexdigest()

def task(directory,fp,fn):
    directory=Path(directory);marker=directory/'COMPLETE.json'
    if marker.exists():
        old=json.loads(marker.read_text())
        if old['fingerprint']!=fp:raise ValueError('stale task: use a new run directory')
        for name,digest in old['outputs'].items():
            if not (directory/name).is_file() or sha(directory/name)!=digest:raise ValueError('completed output corrupted')
        return
    if directory.exists():
        # Preserve failures; never mix partial callers into completed evidence.
        import time
        directory.rename(directory.with_name(directory.name+'.failed.'+str(time.time_ns())))
    directory.mkdir(parents=True)
    fn(directory)
    outputs={str(p.relative_to(directory)):sha(p) for p in directory.rglob('*') if p.is_file() and not p.is_symlink()}
    if not outputs:raise ValueError('empty task output')
    write_json(marker,dict(fingerprint=fp,outputs=outputs,status='COMPLETE'))

def pipeline(a):
    from . import runtime,evidence,analysis,secondary
    c=load_config(a.config);out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True)
    lock=(out/'.lock').open('w')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if c['threads']<1 or c['threads']>int(os.environ.get('AGING_CPU_LIMIT','8')):raise ValueError('thread budget exceeded')
    c['threads']=int(c['threads'])
    resources=freeze_resources(c,out/'resources.lock.json')
    if a.stage=='preflight':
        rows=selected_samples(a.samples,'qc');checks=[]
        for r in rows:
            try:sm=runtime.validate_bam(r,c);checks.append(dict(sample_id=r['sample_id'],status='PASS',bam_sample=sm))
            except Exception as e:checks.append(dict(sample_id=r['sample_id'],status='FAIL',reason=type(e).__name__))
        write_tsv(out/'preflight.private.tsv',checks,['sample_id','status','bam_sample','reason'])
        failures=sum(r['status']!='PASS' for r in checks)
        write_json(out/'preflight.json',dict(samples=len(checks),failures=failures,status='PASS' if checks and not failures else 'FAIL'))
        if failures or not checks:raise ValueError('BAM preflight failed')
        return
    rows=selected_samples(a.samples,'qc' if a.stage=='qc' else ('pilot' if a.stage=='pilot' else 'cohort'),a.pilot_ids)
    if any(r['assay_id'] not in c['assays'] for r in rows):raise ValueError('unconfigured assay')
    fp=fingerprint(rows,c,resources)
    if a.stage=='qc':
        for r in rows:task(out/'qc'/r['sample_id'],fp,lambda d,r=r:runtime.qc(r,c,d))
        write_json(out/'QC_COMPLETE.json',dict(samples=len(rows),fingerprint=fp));return
    if not a.pon:raise ValueError('explicit technical PoN directory required')
    ponroot=Path(a.pon)
    if not (ponroot/'pon.complete.json').is_file():raise ValueError('PoN not completed')
    pm=load_config(ponroot/'pon.complete.json')
    if pm.get('status')!='COMPLETE' or pm.get('resource_lock_sha256')!=sha(out/'resources.lock.json'):
        raise ValueError('PoN reference/resource lock changed')
    if set(pm['technical_people'])&{r['person_id'] for r in read_tsv(a.samples)}:
        raise ValueError('PoN overlaps biological people')
    for name,digest in pm['outputs'].items():
        if sha(ponroot/name)!=digest:raise ValueError('PoN output changed')
    g=None
    if a.stage in ['cohort','features','models','secondary']:
        if not a.gate:raise ValueError('frozen calibration gate required')
        g=validate_gate(a.gate,a.samples,out/'resources.lock.json',a.config)
        if set(r['assay_id'] for r in rows)-set(g['assays']):raise ValueError('assay not calibrated')
        if any(c['evidence'][k]!=g['tier'][k] for k in ['min_baseq','min_mapq']):raise ValueError('QC/evidence/tier filters differ')
    # Require successful uniform QC for every selected participant before calling.
    for r in rows:
        p=out/'qc'/r['sample_id']/'qc.json'
        if not p.is_file():raise ValueError('uniform QC missing; run qc first')
        marker=load_config(p.parent/'COMPLETE.json')
        if marker['outputs'].get('qc.json')!=sha(p):raise ValueError('QC output hash mismatch')
        q=json.loads(p.read_text())
        bp=Path(r['bam_path'])
        if q['bam_path']!=str(bp) or q['bam_size']!=bp.stat().st_size or q['bam_mtime_ns']!=bp.stat().st_mtime_ns or q['evidence_filters']!=c['evidence'] or q['common_bed_sha256']!=sha(runtime.resource(c,'common_bed')):
            raise ValueError('QC is stale for BAM or filters')
        if q['status']!='TECHNICAL_QC_COMPLETE' or q['mean_depth']<c['qc']['min_mean_depth']:
            raise ValueError('uniform QC threshold failed')
    work=out/('pilot_calls' if a.stage=='pilot' else 'cohort_calls')
    if a.stage in ['pilot','cohort']:
        calls={}
        for r in rows:
            pon=ponroot/safe_id(r['assay_id'])/'pon.vcf.gz'
            if not pon.is_file():raise ValueError('assay PoN absent')
            rfp=fingerprint([r],c,dict(resources=resources,pon=sha(pon)))
            task(work/r['sample_id'],rfp,lambda d,r=r,p=pon:runtime.call_sample(r,c,d,str(p)))
            calls[r['sample_id']]={n:str(work/r['sample_id']/(n+'.norm.vcf.gz')) for n in ['mutect','lofreq','sage','germline']}
        if a.calibration_truth and a.stage!='pilot':raise ValueError('calibration truth only permitted in pilot')
        if a.calibration_truth:
            runtime.normalize(a.calibration_truth,out/'calibration.truth.norm.vcf.gz',c,out/'calibration.normalize.log')
        candidates=evidence.union(calls,out/(a.stage+'.candidates.tsv'),runtime.resource(c,'common_bed'),runtime.resource(c,'driver_whitelist'),
            out/'calibration.truth.norm.vcf.gz' if a.calibration_truth else None)
        # Evidence has its own completion marker because the cohort union can expand.
        for r in rows:
            ed=out/(a.stage+'_evidence')/r['sample_id'];efp=fingerprint([r],c,dict(candidates=candidates,call_fingerprint=fp))
            task(ed,efp,lambda d,r=r:evidence.interrogate(r,candidates,calls[r['sample_id']],c,
                str(ponroot/safe_id(r['assay_id'])/'pon.vcf.gz'),d/'evidence.tsv'))
        write_json(out/(a.stage.upper()+'_CALLS_COMPLETE.json'),dict(samples=len(rows),candidate_sites=len(candidates),fingerprint=fp))
        if a.stage=='pilot':return
    if a.stage in ['cohort','features']:
        analysis.features(rows,c,g,out/'cohort_evidence',out/'features')
    if a.stage in ['cohort','models']:
        analysis.model(out/'features/features.private.tsv',g,out/'analysis')
    if a.stage in ['cohort','secondary']:
        for r in rows:
            task(out/'secondary'/r['sample_id'],fp,lambda d,r=r:secondary.secondary_sample(r,c,work/r['sample_id'],d))
            secondary.spectra(out/'features'/(r['sample_id']+'.classified.tsv'),runtime.resource(c,'fasta'),out/'secondary'/r['sample_id']/'spectrum.json')
    if a.stage=='cohort':
        write_json(out/'COMPUTATIONAL_RUN_COMPLETE.json',dict(samples=len(rows),fingerprint=fp,
            status='COMPUTATIONAL_RUN_COMPLETE_REQUIRES_SCIENTIFIC_REVIEW',
            not_claimed=['orthogonal validation','clonal CNV confirmation','causal aging effect']))
        from .report import report
        report(out,out/'report.aggregate.html')

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    from .explore import add_parser
    add_parser(sub)
    q=sub.add_parser('prepare');q.add_argument('--truthset',required=True);q.add_argument('--registry',required=True)
    q.add_argument('--clinical');q.add_argument('--germline-index');q.add_argument('--assay-aliases');q.add_argument('--out',required=True)
    q=sub.add_parser('run');q.add_argument('--config',required=True);q.add_argument('--samples',required=True);q.add_argument('--out',required=True)
    q.add_argument('--stage',choices=['preflight','qc','pilot','cohort','features','models','secondary'],required=True)
    q.add_argument('--pilot-ids');q.add_argument('--gate');q.add_argument('--pon');q.add_argument('--calibration-truth')
    q=sub.add_parser('pon');q.add_argument('--config',required=True);q.add_argument('--samples',required=True)
    q.add_argument('--biological-samples',required=True);q.add_argument('--out',required=True)
    q=sub.add_parser('resample');q.add_argument('--bam',required=True);q.add_argument('--out',required=True)
    q.add_argument('--fraction',type=float,required=True);q.add_argument('--seed',type=int,default=1701)
    q=sub.add_parser('split-half');q.add_argument('--bam',required=True);q.add_argument('--out',required=True);q.add_argument('--seed',type=int,default=1701)
    q=sub.add_parser('spike');q.add_argument('--bam',required=True);q.add_argument('--sites',required=True);q.add_argument('--fasta',required=True)
    q.add_argument('--out',required=True);q.add_argument('--kind',choices=['snv','indel'],required=True);q.add_argument('--seed',type=int,default=1701);q.add_argument('--threads',type=int,default=2)
    q=sub.add_parser('calibrate');q.add_argument('--manifest',required=True);q.add_argument('--out',required=True)
    q=sub.add_parser('compare-truth');
    for k in ['truth-vcf','evidence','out','vaf-bin','depth-bin','context']:q.add_argument('--'+k,required=True)
    q=sub.add_parser('power');q.add_argument('--features',required=True);q.add_argument('--out',required=True);q.add_argument('--reps',type=int,default=1000);q.add_argument('--design')
    q=sub.add_parser('pilot-score');q.add_argument('--config',required=True);q.add_argument('--samples',required=True)
    q.add_argument('--pilot-ids',required=True);q.add_argument('--rules',required=True);q.add_argument('--work',required=True);q.add_argument('--out',required=True)
    q=sub.add_parser('normalize-telseq');q.add_argument('--table',required=True);q.add_argument('--out',required=True)
    q=sub.add_parser('report');q.add_argument('--work',required=True);q.add_argument('--out',required=True)
    q=sub.add_parser('ancestry');q.add_argument('--samples',required=True);q.add_argument('--config',required=True);q.add_argument('--out',required=True)
    q=sub.add_parser('signatures');q.add_argument('--spectrum',required=True);q.add_argument('--plan',required=True);q.add_argument('--out',required=True)
    q=sub.add_parser('mosaic-screen');q.add_argument('--segments',required=True);q.add_argument('--baf',required=True);q.add_argument('--plan',required=True);q.add_argument('--out',required=True)
    q=sub.add_parser('regional-burden');q.add_argument('--features',required=True);q.add_argument('--regions',required=True);q.add_argument('--plan',required=True);q.add_argument('--out',required=True)
    q=sub.add_parser('annotate');q.add_argument('--candidates',required=True);q.add_argument('--plan',required=True);q.add_argument('--out',required=True)
    q=sub.add_parser('derive-manifest');q.add_argument('--samples',required=True);q.add_argument('--sample-id',required=True);q.add_argument('--bam',required=True);q.add_argument('--out',required=True)
    a=p.parse_args()
    try:
        from . import cohort,calibration,analysis,secondary,runtime
        if a.command=='prepare':print(json.dumps(cohort.prepare(a.truthset,a.registry,a.out,a.clinical,a.germline_index,a.assay_aliases),indent=2))
        elif a.command=='explore':
            from .explore import explore
            return explore(a)
        elif a.command=='label-origins':
            from .explore import label_table
            label_table(a)
        elif a.command=='run':pipeline(a)
        elif a.command=='pon':
            c=load_config(a.config);technical=read_tsv(a.samples);bio=read_tsv(a.biological_samples)
            for k in ['person_id','analysis_id','bam_path']:
                if not all(r.get(k) for r in technical):raise ValueError('technical PoN identity incomplete')
                if len({r[k] for r in technical})!=len(technical):raise ValueError('duplicate technical PoN identity')
                if set(r[k] for r in technical)&set(r[k] for r in bio):raise ValueError('technical PoN overlaps biological cohort')
            freeze_resources(c,Path(a.out)/'resources.lock.json');runtime.build_pon(technical,c,a.out)
        elif a.command=='resample':calibration.resample(a.bam,a.out,a.fraction,a.seed)
        elif a.command=='split-half':calibration.split_half(a.bam,a.out,a.seed)
        elif a.command=='spike':calibration.spike(a.bam,a.sites,a.fasta,a.out,a.kind,a.seed,a.threads)
        elif a.command=='calibrate':calibration.summarize(a.manifest,a.out)
        elif a.command=='compare-truth':calibration.compare_truth(a.truth_vcf,a.evidence,a.out,a.vaf_bin,a.depth_bin,a.context)
        elif a.command=='power':analysis.power(a.features,a.out,reps=a.reps,design=a.design)
        elif a.command=='pilot-score':
            rows=selected_samples(a.samples,'pilot',a.pilot_ids)
            rules=load_config(a.rules)
            if rules.get('status')!='PROVISIONAL_PILOT_ONLY':raise ValueError('explicit provisional pilot rules required')
            analysis.features(rows,load_config(a.config),rules,a.work,a.out)
            write_json(Path(a.out)/'PILOT_ONLY.json',dict(status='PROVISIONAL_PILOT_NOT_COHORT_RESULT',rules_sha256=sha(a.rules)))
        elif a.command=='normalize-telseq':secondary.normalize_telseq(a.table,a.out)
        elif a.command=='report':
            from .report import report
            report(a.work,a.out)
        elif a.command=='ancestry':
            from .ancestry import ancestry
            c=load_config(a.config);freeze_resources(c,Path(a.out)/'resources.lock.json')
            ancestry(selected_samples(a.samples,'qc'),c,a.out)
        elif a.command=='signatures':
            from .extensions import signature_fit
            signature_fit(a.spectrum,a.plan,a.out)
        elif a.command=='mosaic-screen':
            from .extensions import mosaic_candidates
            mosaic_candidates(a.segments,a.baf,a.plan,a.out)
        elif a.command=='regional-burden':
            from .extensions import regional_burdens
            regional_burdens(a.features,a.regions,a.plan,a.out)
        elif a.command=='annotate':
            from .annotation import annotate
            annotate(a.candidates,a.plan,a.out)
        elif a.command=='derive-manifest':
            rows=[r for r in read_tsv(a.samples) if r['sample_id']==a.sample_id]
            if len(rows)!=1 or rows[0]['analysis_status']!='METADATA_READY':raise ValueError('derived manifest requires one ready parent')
            r=dict(rows[0]);r.update(sample_id=safe_id(a.sample_id+'_DERIVED'),origin='COMPUTATIONAL_CALIBRATION',
                bam_path=a.bam,bai_path=a.bam+'.bai',bam_sha256=sha(a.bam),parent_sample_id=a.sample_id)
            if not Path(r['bai_path']).is_file():raise ValueError('derived BAM index missing')
            write_tsv(a.out,[r]);atomic(str(a.out)+'.pilot.ids',r['sample_id']+'\n')
    except Exception as e:
        # Full tracebacks can leak clinical paths; explicitly opt in during private debugging.
        if os.environ.get('AGING_DEBUG')=='1':raise
        print('FAILED: '+str(e),file=sys.stderr);return 2
    return 0

if __name__=='__main__':sys.exit(main())
