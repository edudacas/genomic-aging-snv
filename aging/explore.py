"""Generate reusable evidence; scientific eligibility and origin labels do not delete rows."""
import collections
import contextlib
import csv
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
from .common import read_tsv, write_tsv, write_json, sha, safe_id, number, load_config
from .contracts import bed_intervals, territory_index, within_territory

PAPER = '10.64898/2026.04.10.26350605'


def add_parser(sub):
    p = sub.add_parser('explore', help='QC, calls and unfiltered evidence; no scientific gate')
    for k in ['config', 'samples', 'out']: p.add_argument('--'+k, required=True)
    p.add_argument('--sample-ids', help='Optional subset of exact sample IDs, one per line')
    p.add_argument('--ibd', help='Optional independently computed shared-haplotype TSV')
    p.add_argument('--preflight-only', action='store_true')
    p = sub.add_parser('label-origins', help='Relabel an evidence TSV without rerunning callers')
    for k in ['evidence', 'out']: p.add_argument('--'+k, required=True)
    p.add_argument('--ibd'); p.add_argument('--tp53-variants', help='Optional exact variant IDs, one per line')
    p.add_argument('--heterozygote-mean', type=float, default=0.46)
    p.add_argument('--mean-source', default='UKB_PAPER_REFERENCE_NOT_LOCAL_CALIBRATION')


def selected(rows, ids=None):
    seen = set(); requested = set(ids) if ids is not None else None
    for r in rows:
        sid = safe_id(r['sample_id'])
        if sid in seen: raise ValueError('duplicate sample ID')
        seen.add(sid)
    if requested is not None and requested-seen: raise ValueError('unknown requested sample ID')
    result=[]
    for original in rows:
        r=dict(original)
        # Missing clinical covariates, pending germline adjudication and batch membership are flags.
        reasons=set(filter(None,r.get('technical_reasons','').split(';')))-{'STRATUM_MISSING'}
        if requested is not None and r['sample_id'] not in requested: state='NOT_REQUESTED'
        elif reasons or not r.get('bam_path') or not r.get('bai_path'): state='TECHNICAL_INPUT_PENDING'
        elif r.get('technical_status')=='HOLD' and not r.get('technical_reasons'): state='TECHNICAL_HOLD_UNSPECIFIED'
        else: state='SELECTED_NOT_SERVER_VALIDATED'
        r['explore_status']=state
        result.append(r)
    return result


def freeze(c, output):
    """Hash actual technical resources; no signed gate or local calibration is needed."""
    res=c['resources']; callers=c.get('callers',['mutect','lofreq','sage'])
    if not callers or len(set(callers))!=len(callers) or set(callers)-{'mutect','lofreq','sage'}:
        raise ValueError('choose one or more distinct callers: mutect, lofreq, sage')
    required={'fasta','fai','dict','common_bed'}
    if 'sage' in callers:
        required.update(['sage_jar','sage_hotspots','sage_high_confidence','sage_driver_panel','sage_ensembl_manifest'])
    if required-set(res): raise ValueError('missing technical resources: '+','.join(sorted(required-set(res))))
    locked={}
    for name, entry in res.items():
        p=Path(entry['path']).resolve(strict=True)
        if not p.is_file() or not p.stat().st_size: raise ValueError('empty resource '+name)
        digest=sha(p)
        if entry.get('sha256') and entry['sha256']!=digest: raise ValueError('checksum mismatch '+name)
        if entry.get('build',c['build'])!=c['build']: raise ValueError('reference build mismatch '+name)
        locked[name]={'path':str(p),'sha256':digest}
    fasta=Path(res['fasta']['path'])
    if Path(str(fasta)+'.fai').resolve()!=Path(res['fai']['path']).resolve(): raise ValueError('FAI must be beside FASTA')
    if fasta.with_suffix('.dict').resolve()!=Path(res['dict']['path']).resolve(): raise ValueError('dictionary must be beside FASTA')
    lengths={x[0]:int(x[1]) for x in (l.split('\t') for l in Path(res['fai']['path']).read_text().splitlines())}
    bed_intervals(res['common_bed']['path'],lengths)
    for name in ['gnomad','abraom','contamination_sites']:
        if name in res:
            idx=res.get(name+'_index',{}).get('path','')
            if not idx or Path(idx).resolve() not in [Path(res[name]['path']+e).resolve() for e in ['.tbi','.csi']]:
                raise ValueError('indexed optional resource required: '+name)
    if 'sage' in callers:
        cache=Path(c['sage_ensembl_dir']).resolve(strict=True)
        entries=read_tsv(res['sage_ensembl_manifest']['path'])
        if not entries: raise ValueError('empty SAGE cache manifest')
        for row in entries:
            p=(cache/row['relative_path']).resolve(strict=True)
            if not p.is_relative_to(cache) or sha(p)!=row['sha256']: raise ValueError('SAGE cache mismatch')
    for assay, entry in c.get('assays',{}).items():
        p=entry.get('capture_bed')
        if p:
            bed_intervals(p,lengths); digest=sha(p)
            if entry.get('capture_sha256') and entry['capture_sha256']!=digest: raise ValueError('capture checksum mismatch')
            locked['capture:'+assay]={'path':p,'sha256':digest}
        cnv=entry.get('cnv_reference')
        if c.get('secondary',{}).get('cnv') and cnv:
            # Optional CNV failures are reported by that module after variant export.
            digest=sha(cnv['path']) if Path(cnv['path']).is_file() else None
            locked['cnv:'+assay]={'path':cnv['path'],'sha256':digest,'expected_sha256':cnv.get('sha256')}
    for assay,p in c.get('pon_by_assay',{}).items():
        locked['pon:'+assay]={'path':p,'sha256':sha(p)}
        indexes=[Path(p+e) for e in ['.tbi','.csi'] if Path(p+e).is_file()]
        if not indexes: raise ValueError('PoN index missing')
        for index in indexes: locked[str(index)]={'sha256':sha(index)}
    if output: write_json(output,{'build':c['build'],'resources':locked,'profile':'exploratory'})
    return locked


def load_ibd(path, build='GRCh38'):
    result=collections.defaultdict(list)
    if not path: return result
    for r in read_tsv(path):
        # Segment coordinates are 1-based inclusive, matching VCF variant_id.
        ch,pos,ref,alt=r['variant_id'].split(':'); pos=int(pos)
        n=number(r.get('shared_cm'))
        if (r.get('build')!=build or ch!=r.get('chrom') or not r.get('source')
            or not r.get('partner_sample_id') or r['partner_sample_id']==r['sample_id']
            or r.get('phased')!='YES' or r.get('partner_same_variant')!='YES'
            or n is None or n<0 or not int(r['start_1based'])<=pos<=int(r['end_1based'])):
            raise ValueError('invalid shared-haplotype evidence')
        result[(r['sample_id'],r['variant_id'])].append(r)
    return result


def labels(r, ibd=None, mean=0.46, mean_source='UKB_PAPER_REFERENCE_NOT_LOCAL_CALIBRATION', tp53=False):
    from scipy.stats import binom
    if not 0<mean<1: raise ValueError('heterozygote mean must be in (0,1)')
    if mean!=0.46 and mean_source=='UKB_PAPER_REFERENCE_NOT_LOCAL_CALIBRATION':
        raise ValueError('custom mean requires its source')
    dp=number(r.get('mutect_dp')); alt=number(r.get('mutect_ad_alt'))
    source='MUTECT_FORMAT_DP_AD'
    if dp is None or alt is None:
        dp=number(r.get('dp')); alt=number(r.get('ad_alt')); source='FORCED_FRAGMENT_COUNTS_ADAPTATION'
    valid=dp is not None and alt is not None and dp==int(dp) and alt==int(alt) and 0<=alt<=dp and dp>0
    p46=float(binom.cdf(alt,dp,.46)) if valid else ''
    p50=float(binom.cdf(alt,dp,.50)) if valid else ''
    p=float(binom.cdf(alt,dp,mean)) if valid else ''
    shares=(ibd or {}).get((r['sample_id'],r['variant_id']),[])
    qualifying=[x for x in shares if float(x['shared_cm'])>=1.5]
    ibd_status='QUALIFYING_SHARED_HAPLOTYPE' if qualifying else ('BELOW_1_5_CM' if shares else 'NOT_AVAILABLE')
    flag='PASS' if valid and dp>=10 and alt>=3 else 'LOW_OR_MISSING_SUPPORT'
    strands=min(number(r.get('alt_forward')) or 0,number(r.get('alt_reverse')) or 0)>=1
    if not valid: hint='COUNTS_UNAVAILABLE'
    elif alt==0: hint='NO_ALT_OBSERVED'
    elif p<.01: hint='SOMATIC_LIKE_VAF'
    else: hint='GERMLINE_COMPATIBLE_VAF'
    if flag!='PASS': origin='INSUFFICIENT_SUPPORT'
    elif hint=='SOMATIC_LIKE_VAF': origin='LIKELY_SOMATIC'
    elif qualifying: origin='LIKELY_GERMLINE'
    else: origin='INCONCLUSIVE'
    # QC is exported independently. Likely labels are model-based, not clinical adjudication.
    return dict(origin_label=origin,vaf_origin_hint=hint,origin_counts_source=source,
        origin_dp=dp if valid else '',origin_ad_alt=alt if valid else '',
        binomial_p_046=p46,binomial_p_050=p50,binomial_p_selected=p,
        heterozygote_mean=mean,heterozygote_mean_source=mean_source,
        paper_depth_alt_flag=flag,fragment_both_strands=int(strands),
        paper_strand_score_exception='NOT_EVALUATED_USE_RAW_MUTECT_VCF',
        ibd_status=ibd_status,ibd_max_cm=max((float(x['shared_cm']) for x in shares),default=''),
        ibd_source=';'.join(sorted({x['source'] for x in shares})),
        classification_scope='TP53_CURATED_ALLELE' if tp53 else 'EXPLORATORY_EXTRAPOLATION',
        origin_method='PAPER_INSPIRED_NOT_CLINICALLY_VALIDATED',origin_reference=PAPER)


def label_table(a):
    ibd=load_ibd(a.ibd)
    tp53=set(Path(a.tp53_variants).read_text().split()) if a.tp53_variants else set()
    rows=read_tsv(a.evidence)
    for r in rows:
        scope=r['variant_id'] in tp53 or (not a.tp53_variants and r.get('classification_scope')=='TP53_CURATED_ALLELE')
        r.update(labels(r,ibd,a.heterozygote_mean,a.mean_source,scope))
    write_tsv(a.out,rows,list(rows[0]) if rows else ['sample_id','variant_id','origin_label'])


def records(path, caller):
    from .caller_records import load_records
    return load_records(path, caller)


def observe(sample, c, paths, output, ibd, tp53, cohort_candidates=None):
    import pysam
    from .evidence import pileup_site,population_af,variants
    callers={name:records(paths.get(name),name) for name in ['mutect','lofreq','sage','germline']}
    forced=records(paths.get('carrier'),'carrier')
    territory=territory_index(c['resources']['common_bed']['path'])
    candidates=dict(cohort_candidates or {})
    for values in list(callers.values())+[forced]:
        for k,r in values.items():
            if within_territory(territory,r['chrom'],r['pos'],r['ref']):
                candidates[k]={x:r[x] for x in ['chrom','pos','ref','alt']}
    if c['resources'].get('driver_whitelist'):
        for r in read_tsv(c['resources']['driver_whitelist']['path']):
            r={x:r[x] for x in ['chrom','pos','ref','alt']}; r['pos']=int(r['pos'])
            if within_territory(territory,r['chrom'],r['pos'],r['ref']):
                candidates['{chrom}:{pos}:{ref}:{alt}'.format(**r)]=r
    ponpath=c.get('pon_by_assay',{}).get(sample['assay_id'])
    pon={k for k,r,alt in variants(ponpath)} if ponpath else None
    model=c.get('origin_model',{}); rows=[]
    with contextlib.ExitStack() as stack:
        bam=stack.enter_context(pysam.AlignmentFile(sample['bam_path'],'rb',index_filename=sample['bai_path']))
        fasta=stack.enter_context(pysam.FastaFile(c['resources']['fasta']['path']))
        populations={n:stack.enter_context(pysam.VariantFile(c['resources'][n]['path']))
                     for n in ['gnomad','abraom'] if c['resources'].get(n)}
        for k,site in sorted(candidates.items()):
            ev=pileup_site(bam,fasta,site,c['evidence']['min_baseq'],c['evidence']['min_mapq'],c['evidence']['max_depth'])
            row=dict(sample_id=sample['sample_id'],variant_id=k,**site,**ev,
                     group_curated=sample.get('group',''),genetic_status_curated=sample.get('genetic_status',''),
                     clinical_hold_flags=sample.get('hold_reasons',''),carrier_variant_exact=int(k in forced),
                     in_pon=int(k in pon) if pon is not None else '',pon_status='AVAILABLE' if pon is not None else 'NOT_AVAILABLE',
                     germline_overlap=int(k in callers['germline']) if paths.get('germline') else '',
                     germline_vcf_status=sample['germline_vcf_status'])
            for name,values in callers.items():
                record=values.get(k)
                row[name+'_filter']=record['filter'] if record else ('NOT_CALLED' if paths.get(name) else 'NOT_AVAILABLE')
                for field in ['dp','ad_alt','ad_ref']: row[name+'_'+field]=record[field] if record else ''
            row['caller_support']=sum(k in callers[n] for n in ['mutect','lofreq','sage'])
            for name in ['gnomad','abraom']:
                row[name+'_af']=population_af(populations[name],site) if name in populations else ''
                row[name+'_resource_status']='AVAILABLE' if name in populations else 'NOT_AVAILABLE'
            if c.get('label_origins',True):
                row.update(labels(row,ibd,model.get('heterozygote_mean',.46),model.get('mean_source','UKB_PAPER_REFERENCE_NOT_LOCAL_CALIBRATION'),k in tp53))
            else:
                row['origin_label']='NOT_CLASSIFIED_CALLER_EVIDENCE_ONLY'
            rows.append(row)
    fields=list(dict.fromkeys(k for r in rows for k in r)) if rows else ['sample_id','variant_id','origin_label']
    write_tsv(output,rows,fields)
    return rows


def export_tables(files, output):
    # Streaming export keeps cohort memory bounded and retains every evidence row.
    fields=[]
    for p in files:
        with Path(p).open() as f:
            for name in csv.DictReader(f,delimiter='\t').fieldnames or []:
                if name not in fields: fields.append(name)
    fields=fields or ['sample_id','variant_id','origin_label']
    tmp=Path(str(output)+'.tmp'); counts=collections.Counter(); total=0
    try:
        with gzip.open(tmp,'wt',newline='') as f:
            w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n'); w.writeheader()
            for p in files:
                with Path(p).open() as src:
                    for row in csv.DictReader(src,delimiter='\t'):
                        w.writerow(row);total+=1;counts[row['origin_label']]+=1
        os.replace(tmp,output)
    finally:
        if tmp.exists(): tmp.unlink()
    return dict(evidence_rows=total,origin_labels=dict(counts),sha256=sha(output))


def explore(a):
    out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True)
    with (out/'.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        return _explore_locked(a,out)


def _explore_locked(a,out):
    from . import runtime
    from .cli import fingerprint,task
    c=load_config(a.config); c['profile']='exploratory'
    if c.get('origin_model',{}).get('alpha',.01)!=.01: raise ValueError('paper profile uses alpha=0.01')
    if not 1<=c['threads']<=int(os.environ.get('AGING_CPU_LIMIT','8')): raise ValueError('thread budget exceeded')
    ids=Path(a.sample_ids).read_text().split() if a.sample_ids else None
    rows=selected(read_tsv(a.samples),ids); tp53={r['carrier_variant'] for r in rows if r.get('carrier_variant')}
    ibd=load_ibd(a.ibd,c['build']); resources=freeze(c,None)
    code={p.name:sha(p) for p in Path(__file__).parent.glob('*.py')}
    run_id=hashlib.sha256(json.dumps(dict(rows=rows,config=c,resources=resources,code=code,ibd_sha256=sha(a.ibd) if a.ibd else None),sort_keys=True).encode()).hexdigest()
    previous=out/'RUN_INPUTS.json'
    if previous.exists() and load_config(previous)['fingerprint']!=run_id: raise ValueError('inputs changed; use a new output directory')
    write_json(out/'resources.lock.json',{'build':c['build'],'resources':resources,'profile':'exploratory'})
    write_json(previous,{'fingerprint':run_id,'profile':'exploratory'})
    files=[]; failures=0; completed=0
    for row in rows:
        if row['explore_status']!='SELECTED_NOT_SERVER_VALIDATED': continue
        sid=row['sample_id']; sample=dict(row)
        v=sample.get('germline_vcf')
        sample['germline_vcf_status']='AVAILABLE' if v and Path(v).is_file() else ('PATH_UNAVAILABLE' if v else 'NOT_LINKED')
        if sample['germline_vcf_status']!='AVAILABLE': sample['germline_vcf']=''
        row['germline_vcf_status']=sample['germline_vcf_status']
        row['contamination_estimation_status']='REQUESTED' if c['resources'].get('contamination_sites') else 'NOT_AVAILABLE'
        try:
            row['bam_sample']=runtime.validate_bam(sample,c)
            if a.preflight_only:
                row['explore_status']='TECHNICAL_PREFLIGHT_PASS'; completed+=1; continue
            fp=fingerprint([sample],c,{'run_id':run_id,'resources':resources})
            task(out/'qc'/sid,fp,lambda d:runtime.qc(sample,c,d))
            d=out/'calls'/sid
            task(d,fp,lambda d:runtime.call_sample(sample,c,d,c.get('pon_by_assay',{}).get(sample['assay_id'])))
            paths={n:str(d/(n+'.norm.vcf.gz')) for n in ['mutect','lofreq','sage','germline','carrier'] if (d/(n+'.norm.vcf.gz')).is_file()}
            # Normalization can change a curated indel key; retain the normalized TP53 identity.
            normalized_tp53=tp53|set(records(paths.get('carrier'),'carrier'))
            task(out/'evidence'/sid,fp,lambda d:observe(sample,c,paths,d/'evidence.tsv',ibd,normalized_tp53))
            files.append(out/'evidence'/sid/'evidence.tsv')
            row['explore_status']='EVIDENCE_COMPLETE'; completed+=1
            row['telseq_status']='NOT_REQUESTED'
            if c.get('secondary',{}).get('telseq'):
                from .secondary import secondary_sample
                sc=dict(c,secondary=dict(c['secondary'],cnv=False))
                try:
                    task(out/'telseq'/sid,fp,lambda td:secondary_sample(sample,sc,out/'calls'/sid,td))
                    row['telseq_status']='RAW_COMPLETE'
                except Exception as e:
                    row['telseq_status']='FAILED';row['telseq_error']=str(e); failures+=1
            row['cnv_status']='NOT_REQUESTED'
            if c.get('secondary',{}).get('cnv'):
                from .secondary import secondary_sample
                sc=dict(c,secondary=dict(c['secondary'],telseq=False))
                try:
                    if not c.get('assays',{}).get(sample['assay_id'],{}).get('cnv_reference'):
                        raise ValueError('CNV requires a reference for this assay')
                    task(out/'cnv'/sid,fp,lambda cd:secondary_sample(sample,sc,out/'calls'/sid,cd))
                    row['cnv_status']='COPY_RATIO_COMPLETE_BAF_SEE_SECONDARY_STATUS'
                except Exception as e:
                    row['cnv_status']='FAILED';row['cnv_error']=str(e);failures+=1
        except Exception as e:
            row['explore_status']='FAILED';row['execution_error']=str(e);failures+=1
        finally:
            write_tsv(out/'samples.all.private.tsv',rows)
    write_tsv(out/'samples.all.private.tsv',rows)
    summary=dict(profile='exploratory',input_samples=len(rows),completed_samples=completed,failed_steps=failures,
        sample_statuses=dict(collections.Counter(r['explore_status'] for r in rows)),preflight_only=a.preflight_only,
        data_scope='PER_SAMPLE_CALLER_UNION_PLUS_AVAILABLE_GERMLINE_AND_FORCED_SITES_WITHIN_CONFIGURED_BED',
        absence_of_row='NOT_ASSESSED_NOT_REFERENCE_GENOTYPE',clinical_filters_applied=False)
    if not a.preflight_only: summary.update(export_tables(files,out/'variants.all.private.tsv.gz'))
    summary['status']='PARTIAL_FAILURE' if failures else ('NO_SAMPLES_PROCESSED' if not completed else ('AVAILABLE_SUBSET_COMPLETE' if completed<len(rows) else 'REQUESTED_INPUTS_COMPLETE'))
    write_json(out/'explore.status.json',summary)
    return 2 if failures or not completed else 0
