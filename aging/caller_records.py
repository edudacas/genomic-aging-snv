"""Lossless per-allele caller records; collisions are explicit, never last-write-wins."""
import json
import pysam

def clean(value):
    if isinstance(value,dict):return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [clean(v) for v in value]
    return value

def load_records(path,caller):
    groups={}
    if not path:return groups
    with pysam.VariantFile(str(path)) as vcf:
        if len(vcf.header.samples)>1:raise ValueError('Explicit sample splitting required')
        for ordinal,r in enumerate(vcf,1):
            sample=next(iter(r.samples.values()),{})
            if caller=='germline' and not any(g and g>0 for g in (sample.get('GT') or [])):continue
            for ai,alt in enumerate(r.alts or [],1):
                if alt.startswith('<') or alt=='*':continue
                key=f'{r.chrom}:{r.pos}:{r.ref}:{alt}'
                ad=sample.get('AD') or ()
                dp=sample.get('DP')
                if caller=='lofreq':dp=r.info.get('DP') if 'DP' in r.header.info else None
                af=r.info.get('AF') if 'AF' in r.header.info else None
                if isinstance(af,tuple):af=af[ai-1] if len(af)>=ai else None
                record=dict(ordinal=ordinal,filter=';'.join(r.filter.keys()) or '.',dp=dp,
                    ad_ref=ad[0] if ad else None,ad_alt=ad[ai] if len(ad)>ai else None,
                    af=af,info=clean(dict(r.info)),samples=clean({s:dict(r.samples[s]) for s in r.samples}),
                    original_alts=list(r.alts),qual=r.qual,id=r.id)
                site=groups.setdefault(key,dict(chrom=r.chrom,pos=r.pos,ref=r.ref,alt=alt,records=[]))
                site['records'].append(record)
    for site in groups.values():
        records=site['records'];site['record_count']=len(records)
        for field in ('filter','dp','ad_ref','ad_alt','af'):
            vals=[r[field] for r in records]
            same=all(v==vals[0] for v in vals)
            site[field]=(vals[0] if vals[0] is not None else '') if same else ('MULTIPLE_RECORDS' if field=='filter' else '')
        site['record_conflict']=len({json.dumps({k:r[k] for k in ('filter','dp','ad_ref','ad_alt','af')},sort_keys=True) for r in records})>1
        site['records_json']=json.dumps(records,sort_keys=True,separators=(',',':'))
    return groups
