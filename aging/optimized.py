"""Resumable window-based cohort evidence; no clinical selection or origin classification."""
import contextlib,csv,gzip,hashlib,json,os,time
from collections import defaultdict
from pathlib import Path
import pysam
from .common import read_tsv,write_json,sha
from .caller_records import load_records
from .window_recount import ordered_sites,windows,recount_window

CALLERS=('mutect','lofreq','sage')
EV_FIELDS=['evidence_status','dp','ad_ref','ad_alt','vaf','alt_forward','alt_reverse','alt_f1r2','alt_f2r1','alt_mean_bq','alt_mean_mq']
EXTENDED_FIELDS=['dp_raw_reads','alt_raw_reads','mq0_raw_reads','duplicate_raw_reads','raw_depth_saturated','filtered_depth_saturated','flag_filtered_raw_reads',
    'low_mq_raw_reads','low_bq_eligible_reads','eligible_reads_before_bq','eligible_reads_after_bq',
    'other_fragments','discordant_raw_fragments','discordant_eligible_fragments','alt_concordant_both_mates',
    'ref_forward','ref_reverse','ref_f1r2','ref_f2r1','coverage_diagnostic']
for prefix in ('ref_','alt_ext_'):
    EXTENDED_FIELDS += [prefix+name for name in ('min_bq','median_bq','sum_bq','min_mq','median_mq','sum_mq',
        'median_distance_to_end','near_end_fraction','unique_fragment_starts','mean_nm','mean_abs_tlen')]
EXTENDED_FIELDS += ['saturation_fallback']
FIELDS=['group_curated','genetic_status_curated','clinical_hold_flags','carrier_variant_exact','sample_id','variant_id','chrom','pos','ref','alt','discovered_by','discovery_sample_count']+EV_FIELDS+EXTENDED_FIELDS+['in_pon','pon_status','germline_overlap','germline_vcf_status','caller_support','origin_label','gnomad_af','gnomad_resource_status','abraom_af','abraom_resource_status']
for caller in CALLERS:
    FIELDS += [caller+'_'+f for f in ('filter','dp','ad_ref','ad_alt','af','record_count','record_conflict','filters_json')]

def atomic_tsv(path,rows,fields):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name('.'+path.name+'.tmp')
    with tmp.open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t',extrasaction='ignore');writer.writeheader();writer.writerows(rows)
    os.replace(tmp,path)

def block_plan(sites,c):
    settings=c.get('recount',{})
    span=settings.get('window_bp',10000);gap=settings.get('gap_bp',1000)
    if not isinstance(span,int) or not 1<=span<=100000 or not isinstance(gap,int) or not 0<=gap<=span:raise ValueError('Invalid recount window settings')
    return list(windows(sites,span,gap))

def populate(sites,c,out,fp):
    """Query each population window once for the cohort, rather than once per sample/site."""
    cache=Path(out)/'population_cache.json';marker=Path(out)/'population_cache.complete.json'
    if marker.exists():
        m=json.loads(marker.read_text())
        if m['fingerprint']!=fp or sha(cache)!=m['sha256']:raise ValueError('Population cache mismatch')
        return json.loads(cache.read_text())
    result={r['variant_id']:{} for r in sites};blocks=block_plan(sites,c)
    for name in ('gnomad','abraom'):
        resource=c['resources'].get(name)
        if not resource:continue
        with pysam.VariantFile(resource['path']) as vcf:
            for rows in blocks:
                chrom=rows[0]['chrom']
                if chrom not in vcf.header.contigs:raise ValueError('Population contig missing: '+chrom)
                desired={r['variant_id'] for r in rows}
                for r in vcf.fetch(chrom,rows[0]['pos']-1,rows[-1]['pos']):
                    af=r.info.get('AF') if 'AF' in r.header.info else None
                    if af is None:continue
                    values=af if isinstance(af,tuple) else (af,)
                    for i,alt in enumerate(r.alts or []):
                        key=f'{r.chrom}:{r.pos}:{r.ref}:{alt}'
                        if key in desired and i<len(values) and values[i] is not None:
                            value=float(values[i]);previous=result[key].get(name)
                            result[key][name]=value if previous is None else max(previous,value)
    write_json(cache,result);write_json(marker,dict(fingerprint=fp,sha256=sha(cache)))
    return result

def integrate(site,ev,sample,records,pon,pop,c):
    key=site['variant_id'];row=dict(site,**ev,sample_id=sample['sample_id'],
        group_curated=sample.get('group',''),genetic_status_curated=sample.get('genetic_status',''),
        clinical_hold_flags=sample.get('hold_reasons',''),carrier_variant_exact=0,
        in_pon=int(key in pon),pon_status='AVAILABLE',germline_overlap='',germline_vcf_status='DEFERRED',
        caller_support=sum(key in records[n] for n in CALLERS),origin_label='NOT_CLASSIFIED_CALLER_EVIDENCE_ONLY')
    for name in ('gnomad','abraom'):
        row[name+'_af']=pop.get(name,'');row[name+'_resource_status']='AVAILABLE' if c['resources'].get(name) else 'NOT_AVAILABLE'
    for name in CALLERS:
        r=records[name].get(key)
        row[name+'_filter']=r['filter'] if r else 'NOT_CALLED'
        for f in ('dp','ad_ref','ad_alt','af'):row[name+'_'+f]=r[f] if r else ''
        row[name+'_record_count']=r['record_count'] if r else 0
        row[name+'_record_conflict']=int(r['record_conflict']) if r else 0
        row[name+'_filters_json']=json.dumps([v['filter'] for v in r['records']],separators=(',',':')) if r else '[]'
    return row

def record_sidecar(path,records,sid):
    tmp=path.with_name('.'+path.name+'.tmp')
    with gzip.open(tmp,'wt') as f:
        w=csv.DictWriter(f,fieldnames=['sample_id','variant_id','caller','record_ordinal','record_json'],delimiter='\t');w.writeheader()
        for caller in CALLERS:
            for key,site in records[caller].items():
                for r in site['records']:w.writerow(dict(sample_id=sid,variant_id=key,caller=caller,record_ordinal=r['ordinal'],record_json=json.dumps(r,sort_keys=True,separators=(',',':'))))
    os.replace(tmp,path)

def recount_sample(sample,sites,c,paths,pop,directory,fp):
    out=Path(directory);out.mkdir(parents=True,exist_ok=True);sid=sample['sample_id'];start=time.monotonic()
    marker=out/'INPUTS.json'
    if marker.exists() and json.loads(marker.read_text())['fingerprint']!=fp:raise ValueError('Recount inputs changed; use a new output')
    write_json(marker,dict(fingerprint=fp))
    chunks=out/'chunks';chunks.mkdir(exist_ok=True)
    records={n:load_records(paths[n],n) for n in CALLERS}
    pon=set(load_records(c['pon_by_assay'][sample['assay_id']],'pon'))
    record_sidecar(out/'caller_records.tsv.gz',records,sid)
    blocks=block_plan(sites,c);done=0;fallbacks=0
    timing=dict(setup_elapsed_seconds=time.monotonic()-start,pileup_elapsed_seconds=0.0,integration_elapsed_seconds=0.0,chunk_io_elapsed_seconds=0.0,cached_windows=0,computed_windows=0)
    def publish(state,row=None):
        write_json(out/'progress.json',dict(sample_id=sid,stage=state,completed_candidates=done,total_candidates=len(sites),
            chromosome=row['chrom'] if row else None,position=row['pos'] if row else None,
            elapsed_seconds=round(time.monotonic()-start,2),timings=timing,saturated_fallback_windows=fallbacks,updated_epoch=time.time()))
    publish('RECOUNT')
    try:
        with pysam.AlignmentFile(sample['bam_path'],'rb',index_filename=sample['bai_path']) as bam,pysam.FastaFile(c['resources']['fasta']['path']) as fasta:
            for i,rows in enumerate(blocks):
                file=chunks/f'{i:07}.tsv';receipt=chunks/f'{i:07}.json'
                cached=False
                if receipt.exists() and file.exists():
                    old=json.loads(receipt.read_text())
                    cached=old['fingerprint']==fp and old['sha256']==sha(file) and old['keys']==[r['variant_id'] for r in rows]
                if not cached:
                    phase_start=time.monotonic()
                    ev,saturated=recount_window(bam,fasta,rows,c['evidence'])
                    pileup_seconds=time.monotonic()-phase_start;timing['pileup_elapsed_seconds']+=pileup_seconds
                    phase_start=time.monotonic()
                    integrated=[integrate(r,ev[r['variant_id']],sample,records,pon,pop[r['variant_id']],c) for r in rows]
                    integration_seconds=time.monotonic()-phase_start;timing['integration_elapsed_seconds']+=integration_seconds
                    phase_start=time.monotonic()
                    atomic_tsv(file,integrated,FIELDS)
                    old=dict(fingerprint=fp,sha256=sha(file),keys=[r['variant_id'] for r in rows],saturated_fallback=saturated,pileup_elapsed_seconds=pileup_seconds,integration_elapsed_seconds=integration_seconds)
                    write_json(receipt,old)
                    timing['chunk_io_elapsed_seconds']+=time.monotonic()-phase_start;timing['computed_windows']+=1
                else:timing['cached_windows']+=1
                fallbacks+=int(old.get('saturated_fallback',False));done+=len(rows);publish('RECOUNT',rows[-1])
        # Ordered bounded-memory assembly. Partial chunks survive failures and can be reused.
        def all_rows():
            for i in range(len(blocks)):
                with (chunks/f'{i:07}.tsv').open() as f:yield from csv.DictReader(f,delimiter='\t')
        phase_start=time.monotonic()
        atomic_tsv(out/'evidence.tsv',all_rows(),FIELDS)
        timing['assembly_elapsed_seconds']=time.monotonic()-phase_start
        timing['elapsed_seconds']=time.monotonic()-start
        timing['scope']='CURRENT_ATTEMPT_WALL_SECONDS; parallel sample durations must not be summed as cohort wall time'
        write_json(out/'timing.json',timing)
        write_json(out/'COMPLETE.json',dict(fingerprint=fp,rows=done,sha256=sha(out/'evidence.tsv'),saturated_fallback_windows=fallbacks))
        publish('COMPLETE',sites[-1] if sites else None)
    except Exception:
        publish('FAILED');raise
    return out/'evidence.tsv'

def run_evidence(samples,c,paths,union,out,fp,workers):
    from .callers import run_jobs
    out=Path(out)
    with pysam.FastaFile(c['resources']['fasta']['path']) as fasta:sites=ordered_sites(union,fasta)
    # Forced sites are not silently added at sample level: any requested site belongs to the global union.
    forced=c['resources'].get('driver_whitelist')
    if forced:
        from .contracts import territory_index,within_territory
        territory=territory_index(c['resources']['common_bed']['path'])
        for row in read_tsv(forced['path']):
            row=dict(row,pos=int(row['pos']));key='{chrom}:{pos}:{ref}:{alt}'.format(**row)
            if not within_territory(territory,row['chrom'],row['pos'],row['ref']):raise ValueError('Forced site outside common territory')
            union.setdefault(key,{k:row[k] for k in ('chrom','pos','ref','alt')})
        with pysam.FastaFile(c['resources']['fasta']['path']) as fasta:sites=ordered_sites(union,fasta)
    atomic_tsv(out/'candidate_union.tsv',sites,['variant_id','chrom','pos','ref','alt','discovered_by','discovery_sample_count'])
    start=time.monotonic();cached=(out/'population_cache.complete.json').exists()
    pop=populate(sites,c,out,fp)
    phase_timing=dict(population_elapsed_seconds=time.monotonic()-start,population_cache_reused=cached,baseline_recount_elapsed_seconds_approx=18000,baseline_phase_timings_available=False,phase_specific_speedup_attribution='NOT_IDENTIFIABLE_FROM_UNINSTRUMENTED_BASELINE')
    write_json(out/'phase_timings.json',phase_timing)
    def run(sample,stage):
        stage('COHORT_WIDE_RECOUNT')
        return recount_sample(sample,sites,c,paths[sample['sample_id']],pop,out/'evidence'/sample['sample_id'],fp)
    start=time.monotonic()
    outputs=run_jobs(samples,workers,out,'COHORT_WIDE_RECOUNT',run)
    phase_timing['sample_recount_wall_elapsed_seconds']=time.monotonic()-start
    phase_timing['samples']={f.parent.name:json.loads((f.parent/'timing.json').read_text()) for f in outputs}
    write_json(out/'phase_timings.json',phase_timing)
    return outputs

def recount_sample_single_pass(sample,sites,c,paths,pop,directory,fp):
    """One full coordinate traversal. Completed samples resume; partial samples restart."""
    from .single_pass import recount_single_pass
    out=Path(directory);out.mkdir(parents=True,exist_ok=True);sid=sample['sample_id'];start=time.monotonic()
    marker=out/'INPUTS.json';complete=out/'COMPLETE.json'
    if marker.exists() and json.loads(marker.read_text())['fingerprint']!=fp:raise ValueError('Single-pass inputs changed; use a new output')
    if complete.exists():
        old=json.loads(complete.read_text())
        if old['fingerprint']==fp and (out/'evidence.tsv').exists() and old['sha256']==sha(out/'evidence.tsv'):
            return out/'evidence.tsv'
    write_json(marker,dict(fingerprint=fp,engine='single_coordinate_pass'))
    records={name:load_records(paths[name],name) for name in CALLERS}
    pon=set(load_records(c['pon_by_assay'][sample['assay_id']],'pon'))
    record_sidecar(out/'caller_records.tsv.gz',records,sid)
    progress_path=out/'progress.json'
    def progress(row,completed,columns):
        write_json(progress_path,dict(sample_id=sid,stage='SINGLE_PASS',completed_candidate_columns=columns,
            integrated_candidates=completed,total_candidates=len(sites),chromosome=row['chrom'],position=row['pos'],
            elapsed_seconds=time.monotonic()-start,updated_epoch=time.time()))
    try:
        scan=time.monotonic()
        with pysam.AlignmentFile(sample['bam_path'],'rb',index_filename=sample['bai_path']) as bam,pysam.FastaFile(c['resources']['fasta']['path']) as fasta:
            evidence,diagnostics=recount_single_pass(bam,fasta,sites,c['evidence']|c.get('single_pass',{}),progress)
        scan_seconds=time.monotonic()-scan;integration=time.monotonic()
        rows=(integrate(site,evidence[site['variant_id']],sample,records,pon,pop[site['variant_id']],c) for site in sites)
        atomic_tsv(out/'evidence.tsv',rows,FIELDS)
        integration_seconds=time.monotonic()-integration
        timing=dict(engine='single_coordinate_pass',scan_elapsed_seconds=scan_seconds,
            integration_and_tsv_elapsed_seconds=integration_seconds,total_elapsed_seconds=time.monotonic()-start,**diagnostics)
        write_json(out/'timing.json',timing)
        write_json(complete,dict(fingerprint=fp,rows=len(sites),sha256=sha(out/'evidence.tsv'),timing=timing))
        write_json(progress_path,dict(sample_id=sid,stage='COMPLETE',integrated_candidates=len(sites),
            total_candidates=len(sites),elapsed_seconds=time.monotonic()-start,updated_epoch=time.time()))
        return out/'evidence.tsv'
    except Exception as error:
        last = json.loads(progress_path.read_text()) if progress_path.exists() else {}
        last.update(sample_id=sid,stage='FAILED',total_candidates=len(sites),
                    error=type(error).__name__,message=str(error),
                    elapsed_seconds=time.monotonic()-start,updated_epoch=time.time())
        write_json(progress_path,last)
        raise

def run_single_pass_evidence(samples,c,paths,union,out,fp,workers):
    from .callers import run_jobs
    out=Path(out)
    with pysam.FastaFile(c['resources']['fasta']['path']) as fasta:sites=ordered_sites(union,fasta)
    atomic_tsv(out/'candidate_union.tsv',sites,['variant_id','chrom','pos','ref','alt','discovered_by','discovery_sample_count'])
    start=time.monotonic();cached=(out/'population_cache.complete.json').exists();pop=populate(sites,c,out,fp)
    phase=dict(engine='single_coordinate_pass',population_elapsed_seconds=time.monotonic()-start,
        population_cache_reused=cached,baseline_total_elapsed_seconds_approx=18000,
        timed_v2_total_elapsed_seconds=11223.317653298378,
        timed_v2_sample_recount_wall_elapsed_seconds=11110.755539759994)
    write_json(out/'phase_timings.json',phase)
    def run(sample,stage):
        stage('SINGLE_PASS_RECOUNT')
        return recount_sample_single_pass(sample,sites,c,paths[sample['sample_id']],pop,out/'evidence'/sample['sample_id'],fp)
    start=time.monotonic();outputs=run_jobs(samples,workers,out,'SINGLE_PASS_RECOUNT',run)
    phase['single_pass_wall_elapsed_seconds']=time.monotonic()-start
    phase['samples']={path.parent.name:json.loads((path.parent/'timing.json').read_text()) for path in outputs}
    write_json(out/'phase_timings.json',phase)
    return outputs
