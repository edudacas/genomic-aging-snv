"""Explicit, resumable paired SAGE BQR benchmark; never prepares cohort config."""
import argparse
from collections import Counter, defaultdict
from decimal import Decimal, InvalidOperation
import fcntl
import gzip
import hashlib
import json
import math
from pathlib import Path

from . import runtime
from .common import load_config, sha, write_json, write_tsv
from .staged_inputs import fingerprint, freeze, selected
from .stages import sample_task, verify_task
from .workbook_input import load_workbook_samples
from .bqr_decision import validate_skip_bqr_decision

PILOTS=('PILOT_A','PILOT_B','PILOT_C')
ARMS=('skip_bqr','enabled')
EDGES=tuple(map(Decimal,('0','0.02','0.03','0.05','0.10','0.20','0.35','0.50','1')))


def validate_comparison(report):
    """Reject placeholders/incomplete pilot reports before accepting a decision."""
    if report.get('status')!='BQR_COMPARISON_COMPLETE':
        raise ValueError('real completed BQR comparison required')
    samples=report.get('samples',[])
    if len(samples)!=3 or {s.get('sample_id') for s in samples}!=set(PILOTS):
        raise ValueError('BQR comparison must cover all three exact pilots')
    for sample in samples:
        for arm in ARMS:
            seconds=sample.get('timing',{}).get(arm,{}).get('sage_wall_elapsed_seconds')
            if type(seconds) not in (int,float) or not math.isfinite(seconds) or seconds<=0:
                raise ValueError('missing real SAGE process timing')
            digest=sample.get('original_vcf_sha256',{}).get(arm,'')
            if len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest):
                raise ValueError('original paired VCF hashes required')
        for arm in ('skip_bqr','bqr_enabled'):
            count=sample.get(arm,{}).get('pass_records')
            if type(count) is not int or count<0:raise ValueError('missing PASS counts')
    return report


def vaf_bucket(value):
    try:
        value=Decimal(str(value))
        if not value.is_finite() or not 0<=value<=1:
            return 'UNKNOWN'
        for low,high in zip(EDGES,EDGES[1:]):
            if low<=value<high or high==1 and value==1:
                return f'[{low},{high}' + (']' if high==1 else ')')
    except (InvalidOperation,TypeError,ValueError):
        pass
    return 'UNKNOWN'


def raw_pass_records(path):
    """Count original records, not normalized alleles. FILTER='.' is not PASS.

    Decimal parsing prevents float32 rounding at prespecified VAF boundaries.
    Multiallelic or absent AF stays UNKNOWN; AD is not silently substituted.
    """
    opener=gzip.open if str(path).endswith('.gz') else open
    counts=defaultdict(Counter);all_pass=Counter();filters=Counter();samples=None;total=0
    with opener(path,'rt') as stream:
        for line in stream:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                samples=line.rstrip('\n').split('\t')[9:]
                if len(samples)!=1:
                    raise ValueError('paired pilot SAGE VCF must have exactly one tumor sample')
                continue
            if line.startswith('#') or not line.strip():
                continue
            fields=line.rstrip('\n').split('\t')
            if samples is None or len(fields)!=10:
                raise ValueError('malformed SAGE VCF')
            chrom,pos,_,ref,alt,_,filt,info,fmt,data=fields
            total+=1;filters[filt]+=1
            if filt!='PASS':
                continue
            identity=(chrom,int(pos),ref,alt)
            tags=dict(item.split('=',1) for item in info.split(';') if '=' in item)
            values=dict(zip(fmt.split(':'),data.split(':')))
            af=values.get('AF','.')
            bucket=vaf_bucket(af) if ',' not in alt and ',' not in af else 'UNKNOWN'
            tier=tags.get('TIER','UNKNOWN')
            if tier in ('','.'):tier='UNKNOWN'
            counts[(bucket,tier)][identity]+=1;all_pass[identity]+=1
    if samples is None:
        raise ValueError('VCF header missing')
    return dict(counts=counts,all_pass=all_pass,filters=dict(filters),total_records=total,
                sample=samples[0])


def compare_sample(sid,skip_path,enabled_path,timing):
    skip,enabled=raw_pass_records(skip_path),raw_pass_records(enabled_path)
    if skip['sample']!=enabled['sample']:
        raise ValueError('paired SAGE sample names differ')
    rows=[]
    for bucket,tier in sorted(set(skip['counts'])|set(enabled['counts'])):
        left=skip['counts'].get((bucket,tier),Counter())
        right=enabled['counts'].get((bucket,tier),Counter())
        rows.append(dict(sample_id=sid,vaf_bin=bucket,sage_tier=tier,
                         skip_bqr_pass=sum(left.values()),bqr_enabled_pass=sum(right.values()),
                         delta_pass=sum(right.values())-sum(left.values()),
                         retained_in_bucket=sum((left&right).values()),
                         lost_from_bucket=sum((left-right).values()),
                         gained_in_bucket=sum((right-left).values())))
    left,right=skip['all_pass'],enabled['all_pass']
    summary=dict(sample_id=sid,timing=timing,
                 skip_bqr=dict(total_records=skip['total_records'],filters=skip['filters'],pass_records=sum(left.values())),
                 bqr_enabled=dict(total_records=enabled['total_records'],filters=enabled['filters'],pass_records=sum(right.values())),
                 retained_pass_records=sum((left&right).values()),
                 lost_pass_records=sum((left-right).values()),gained_pass_records=sum((right-left).values()),
                 original_vcf_sha256=dict(skip_bqr=sha(skip_path),enabled=sha(enabled_path)))
    return summary,rows


def arm_config(base,arm):
    if arm not in ARMS:raise ValueError('unknown BQR arm')
    if base.get('sage_extra_args') not in (['-skip_bqr','-skip_msi_jitter'],['-skip_msi_jitter']):
        raise ValueError('unexpected SAGE extras; preserve frozen protocol')
    if base.get('threads')!=2 or base.get('java_heap_gb')!=4:
        raise ValueError('pilot benchmark requires two threads and 4 GiB heap')
    return dict(base,profile='exploratory',callers=['sage'],
                sage_extra_args=['-skip_msi_jitter'] if arm=='enabled' else ['-skip_bqr','-skip_msi_jitter'])


def arm_fingerprint(fp,sid,arm):
    return hashlib.sha256(json.dumps([fp,sid,arm]).encode()).hexdigest()


def run_benchmark(a):
    base=load_config(a.config)
    for arm in ARMS:arm_config(base,arm)
    rows,input_receipt=load_workbook_samples(a.workbook,a.linkage,expected=1329,check_files=False)
    chosen=[dict(r,germline_vcf='',carrier_variant='') for r in selected(rows,PILOTS)
            if r['sample_id'] in PILOTS]
    if set(r['sample_id'] for r in chosen)!=set(PILOTS) or any(
            r['explore_status']!='SELECTED_NOT_SERVER_VALIDATED' for r in chosen):
        raise ValueError('require exactly the three technically eligible real pilots')
    chosen=sorted(chosen,key=lambda r:PILOTS.index(r['sample_id']))
    out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True)
    with (out/'.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        resources=freeze(arm_config(base,'skip_bqr'),None)
        fp=fingerprint(chosen,base,dict(resources=resources,workbook=input_receipt))
        inputs=out/'BQR_RUN_INPUTS.json'
        if inputs.exists() and load_config(inputs)['fingerprint']!=fp:
            raise ValueError('benchmark inputs changed; use a new directory')
        if not inputs.exists():
            write_json(inputs,dict(fingerprint=fp,selected_samples=list(PILOTS),resources=resources,
                                  config_sha256=sha(a.config),workbook_input=input_receipt,
                                  arm_order='alternating; one SAGE process at a time',
                                  pilot_cache_imported=False,cohort_config_prepared=False))
        summaries=[];buckets=[]
        for i,row in enumerate(chosen):
            sid=row['sample_id'];timing={}
            for arm in (ARMS if i%2==0 else tuple(reversed(ARMS))):
                config=arm_config(base,arm);d=out/arm/sid
                def invoke(directory,config=config,row=row):
                    runtime.call_sample(row,config,directory,None)
                    runtime.run(['bcftools','index','-f','-t',directory/'sage.vcf.gz'],directory/'commands.log')
                    for name in ('sage.vcf.gz','sage.vcf.gz.tbi','sage.norm.vcf.gz',
                                 'sage.norm.vcf.gz.tbi','sage.arguments.txt','sage.timing.json'):
                        if not (directory/name).is_file():raise ValueError('missing benchmark output '+name)
                task,reused=sample_task(d,arm_fingerprint(fp,sid,arm),invoke)
                verify_task(d,arm_fingerprint(fp,sid,arm))
                timing[arm]=dict(sage_wall_elapsed_seconds=load_config(d/'sage.timing.json')['sage_wall_elapsed_seconds'],
                                 sample_total_elapsed_seconds=task['elapsed_seconds'],reused=reused)
            summary,table=compare_sample(sid,out/'skip_bqr'/sid/'sage.vcf.gz',out/'enabled'/sid/'sage.vcf.gz',timing)
            summaries.append(summary);buckets.extend(table)
        write_tsv(out/'pass_by_vaf_tier.tsv',buckets,['sample_id','vaf_bin','sage_tier','skip_bqr_pass',
                  'bqr_enabled_pass','delta_pass','retained_in_bucket','lost_from_bucket','gained_in_bucket'])
        report=dict(status='BQR_COMPARISON_COMPLETE',fingerprint=fp,samples=summaries,
                    selected_bqr_variant=None,cohort_config_prepared=False,
                    strata_table_sha256=sha(out/'pass_by_vaf_tier.tsv'),
                    definitions=dict(pass_definition='exact FILTER=PASS; duplicate original records retained',
                       vaf='FORMAT/AF; decimal edges, left-inclusive/right-exclusive except final bin includes 1',
                       unknown='missing or multiallelic AF is UNKNOWN, never zero',
                       tier='SAGE INFO/TIER, not project calibrated confidence tier',
                       bucket_differences='include changes of VAF/tier for otherwise retained identities',
                       timing='one paired run per sample/arm; cache and order effects remain; no speedup guarantee'))
        validate_comparison(report)
        write_json(out/'bqr_comparison.json',report)
        return 0


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('config','workbook','linkage','out'):parser.add_argument('--'+name,required=True)
    return run_benchmark(parser.parse_args())


if __name__=='__main__':raise SystemExit(main())
