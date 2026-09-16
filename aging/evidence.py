"""Allele-specific, fragment-aware reinterrogation; missing depth never becomes negative."""
import collections
import math
from pathlib import Path
from .common import *
from .contracts import bed_intervals

def key(chrom,pos,ref,alt): return f'{chrom}:{pos}:{ref}:{alt}'

def variants(path):
    import pysam
    with pysam.VariantFile(path) as f:
        for r in f:
            for alt in r.alts or []:
                if alt.startswith('<') or alt=='*':continue
                yield key(r.chrom,r.pos,r.ref,alt),r,alt

def union(calls,output,bed,driver_whitelist,calibration_truth=None):
    territory=collections.defaultdict(list)
    for c,s,e in bed_intervals(bed):territory[c].append((s,e))
    data={}
    for sid,paths in calls.items():
        for caller in ['mutect','lofreq','sage']:
            for k,r,alt in variants(paths[caller]):
                if not any(s<=r.pos-1 and r.pos-1+len(r.ref)<=e for s,e in territory[r.chrom]):continue
                v=data.setdefault(k,dict(variant_id=k,chrom=r.chrom,pos=r.pos,ref=r.ref,alt=alt,
                                        discovered_by=set(),discovery_samples=set()))
                v['discovered_by'].add(caller);v['discovery_samples'].add(sid)
    # Frozen forced-driver sites enter even if no caller discovered them.
    for r in read_tsv(driver_whitelist):
        k=key(r['chrom'],int(r['pos']),r['ref'],r['alt'])
        if not any(s<=int(r['pos'])-1 and int(r['pos'])-1+len(r['ref'])<=e for s,e in territory[r['chrom']]):
            raise ValueError('driver whitelist outside common territory')
        data.setdefault(k,dict(variant_id=k,chrom=r['chrom'],pos=int(r['pos']),ref=r['ref'],alt=r['alt'],
                               discovered_by=set(),discovery_samples=set()))
    if calibration_truth:
        for k,r,alt in variants(calibration_truth):
            if not any(s<=r.pos-1 and r.pos-1+len(r.ref)<=e for s,e in territory[r.chrom]):
                raise ValueError('successful calibration truth outside common territory')
            data.setdefault(k,dict(variant_id=k,chrom=r.chrom,pos=r.pos,ref=r.ref,alt=alt,
                                   discovered_by=set(),discovery_samples=set()))
    rows=[]
    for k in sorted(data):
        v=data[k];v['discovered_by']=';'.join(sorted(v['discovered_by']))
        v['discovery_sample_count']=len(v.pop('discovery_samples'));rows.append(v)
    write_tsv(output,rows,['variant_id','chrom','pos','ref','alt','discovered_by','discovery_sample_count'])
    return rows

def allele_from_pileup(p,ref,alt):
    """Exact CIGAR-anchored allele. Complex delins are explicitly unsupported."""
    q=p.query_position
    if q is None or p.is_refskip or p.is_del:return None
    seq=p.alignment.query_sequence
    if len(ref)==len(alt)==1:
        return seq[q].upper() if p.indel==0 else 'OTHER'
    if len(ref)==1 and alt.startswith(ref):
        n=len(alt)-1
        if p.indel==n:return seq[q:q+n+1].upper()
        return ref if p.indel==0 and seq[q].upper()==ref else 'OTHER'
    if len(alt)==1 and ref.startswith(alt):
        if p.indel==-(len(ref)-1) and seq[q].upper()==alt:return alt
        if p.indel==0 and q+len(ref)>len(seq):return None
        return ref if p.indel==0 and seq[q:q+len(ref)].upper()==ref else 'OTHER'
    return None

def pileup_site(bam,fasta,row,min_bq,min_mq,max_depth):
    chrom=row['chrom'];pos=int(row['pos']);ref=row['ref'];alt=row['alt']
    if fasta.fetch(chrom,pos-1,pos-1+len(ref)).upper()!=ref:raise ValueError('candidate REF mismatch')
    supported=(len(ref)==1 or len(alt)==1) and (len(ref)==len(alt) or ref.startswith(alt) or alt.startswith(ref))
    if not supported:return dict(evidence_status='UNSUPPORTED_COMPLEX_ALLELE',dp=0,ad_ref=0,ad_alt=0,vaf='')
    fragments={};saturated=False
    for col in bam.pileup(chrom,pos-1,pos,truncate=True,stepper='samtools',fastafile=fasta,
                          min_base_quality=min_bq,min_mapping_quality=min_mq,max_depth=max_depth,
                          flag_filter=0xF04,ignore_overlaps=False,ignore_orphans=False,compute_baq=False):
        if col.reference_pos!=pos-1:continue
        saturated=col.nsegments>=max_depth
        for p in col.pileups:
            a=allele_from_pileup(p,ref,alt)
            if a is None:continue
            read=p.alignment;q=p.query_position
            n=len(alt) if len(alt)>len(ref) else 1
            qs=read.query_qualities[q:q+n]
            if len(qs)!=n or min(qs)<min_bq:continue
            orient=('F1R2' if (read.is_read1 != read.is_reverse) else 'F2R1') if read.is_paired else 'UNPAIRED'
            item=(a,min(qs),read.mapping_quality,'R' if read.is_reverse else 'F',orient)
            name=read.query_name
            if name not in fragments:fragments[name]=item
            elif fragments[name] is not None:
                prior=fragments[name]
                fragments[name]=max([prior,item],key=lambda x:x[1]) if prior[0]==a else None
    reads=[v for v in fragments.values() if v is not None];dp=len(reads)
    altreads=[v for v in reads if v[0]==alt];nr=sum(v[0]==ref for v in reads);na=len(altreads)
    return dict(evidence_status='DEPTH_SATURATED' if saturated else ('OBSERVED' if dp else 'NO_COVERAGE'),
                dp=dp,ad_ref=nr,ad_alt=na,vaf=na/dp if dp else '',
                alt_forward=sum(v[3]=='F' for v in altreads),alt_reverse=sum(v[3]=='R' for v in altreads),
                alt_f1r2=sum(v[4]=='F1R2' for v in altreads),alt_f2r1=sum(v[4]=='F2R1' for v in altreads),
                alt_mean_bq=sum(v[1] for v in altreads)/na if na else '',
                alt_mean_mq=sum(v[2] for v in altreads)/na if na else '')

def population_af(vcf,row):
    """Exact normalized allele; missing contig or unnormalized resource is not AF=0."""
    chrom=row['chrom'];pos=int(row['pos'])
    if chrom not in vcf.header.contigs:raise ValueError('population resource contig missing')
    afs=[]
    for r in vcf.fetch(chrom,pos-1,pos):
        if r.pos==pos and r.ref==row['ref'] and row['alt'] in (r.alts or []):
            a=r.info.get('AF')
            if a is not None:
                a=a if isinstance(a,tuple) else (a,)
                afs.append(float(a[list(r.alts).index(row['alt'])]))
    return max(afs) if afs else ''

def interrogate(sample,candidates,paths,c,pon_path,output):
    import pysam
    from scipy.stats import binomtest
    carrier_path=Path(paths['germline']).parent/'carrier.norm.vcf.gz'
    carrier_keys={k for k,r,a in variants(carrier_path)} if carrier_path.exists() else set()
    if sample.get('carrier_variant') and not carrier_keys:raise ValueError('normalized carrier identity absent')
    caller_sets={n:{k:';'.join(r.filter.keys()) for k,r,a in variants(paths[n])} for n in ['mutect','lofreq','sage']}
    germ={k for k,r,a in variants(paths['germline']) if any(
        (r.samples[s].get('GT') and any(g and g>0 for g in r.samples[s]['GT'] if g is not None))
        for s in r.samples)}
    with pysam.VariantFile(paths['germline']) as gv:
        if len(gv.header.samples)!=1:raise ValueError('germline VCF requires one sample; split explicitly')
    pon={k for k,r,a in variants(pon_path)}
    whitelist={key(r['chrom'],int(r['pos']),r['ref'],r['alt']) for r in read_tsv(c['resources']['driver_whitelist']['path'])}
    driver=collections.defaultdict(list)
    for ch,s,e in bed_intervals(c['resources']['driver_bed']['path']):driver[ch].append((s,e))
    rows=[]
    with pysam.AlignmentFile(sample['bam_path'],'rb',index_filename=sample['bai_path']) as bam, \
         pysam.FastaFile(c['resources']['fasta']['path']) as fasta, \
         pysam.VariantFile(c['resources']['gnomad']['path']) as gno, \
         pysam.VariantFile(c['resources']['abraom']['path']) as abr:
        for r in candidates:
            k=r['variant_id'];ev=pileup_site(bam,fasta,r,c['evidence']['min_baseq'],c['evidence']['min_mapq'],c['evidence']['max_depth'])
            n=ev['dp'];a=ev['ad_alt']
            p=binomtest(a,n,0.5,alternative='less').pvalue if n else ''
            entry=dict(sample_id=sample['sample_id'],**r,**ev,
                germline_overlap=int(k in germ),germline_binomial_p=p,
                in_pon=int(k in pon),driver_whitelist=int(k in whitelist),
                in_driver_territory=int(any(s<=int(r['pos'])-1 and int(r['pos'])-1+len(r['ref'])<=e for s,e in driver[r['chrom']])),
                carrier_variant_exact=int(k in carrier_keys),gnomad_af=population_af(gno,r),abraom_af=population_af(abr,r))
            for caller in caller_sets:entry[caller+'_filter']=caller_sets[caller].get(k,'NOT_CALLED')
            entry['caller_support']=sum(k in cs for cs in caller_sets.values())
            rows.append(entry)
    fields=list(rows[0]) if rows else ['sample_id','variant_id','dp','ad_alt','vaf','evidence_status']
    write_tsv(output,rows,fields)
    return rows

def classify(r,tier,floor,min_alt,retain_pon=False):
    reasons=[];dp=int(r['dp']);alt=int(r['ad_alt']);vaf=number(r['vaf'])
    if r['evidence_status']!='OBSERVED':reasons.append(r['evidence_status'])
    if str(r.get('carrier_variant_exact'))=='1':reasons.append('CARRIER_ALLELE')
    if dp==0 or vaf is None:reasons.append('NOT_CALLABLE')
    elif not floor<=vaf<0.5:reasons.append('OUTSIDE_VAF_RANGE')
    if alt<min_alt:reasons.append('ALT_SUPPORT')
    if int(r.get('caller_support',0))<tier['min_callers']:reasons.append('CALLER_CONCORDANCE')
    passing=0
    for name in ['mutect','lofreq','sage']:
        flags=set(r.get(name+'_filter','NOT_CALLED').split(';'))
        # Soft PoN status may be waived only by the explicit whitelist or PoN sensitivity.
        if retain_pon or str(r.get('driver_whitelist'))=='1':flags.discard('panel_of_normals')
        if flags<= {'PASS','.'}:passing+=1
    if passing<tier.get('min_pass_callers',0):reasons.append('CALLER_FILTERS')
    if alt and (number(r.get('alt_mean_bq')) or 0)<tier['min_baseq']:reasons.append('BQ')
    if alt and (number(r.get('alt_mean_mq')) or 0)<tier['min_mapq']:reasons.append('MQ')
    if min(int(r.get('alt_forward',0)),int(r.get('alt_reverse',0)))<tier['min_alt_each_strand']:reasons.append('STRAND')
    if min(int(r.get('alt_f1r2',0)),int(r.get('alt_f2r1',0)))<tier['min_alt_each_orientation']:reasons.append('ORIENTATION')
    if str(r.get('germline_overlap'))=='1' and (number(r.get('germline_binomial_p')) is None or float(r['germline_binomial_p'])>=tier['germline_alpha']):
        reasons.append('GERMLINE_COMPATIBLE')
    if any((number(r.get(k)) or 0)>tier['max_population_af'] for k in ['gnomad_af','abraom_af']):reasons.append('COMMON_POPULATION')
    if not retain_pon and str(r.get('in_pon'))=='1' and str(r.get('driver_whitelist'))!='1':reasons.append('PON')
    return 'PASS' if not reasons else ';'.join(reasons)

def callable_threshold(vaf,min_alt,sensitivity,max_depth=100000):
    from scipy.stats import binom
    lo=0;hi=max_depth
    if binom.sf(min_alt-1,hi,vaf)<sensitivity:raise ValueError('no attainable callable depth')
    while lo<hi:
        mid=(lo+hi)//2
        if binom.sf(min_alt-1,mid,vaf)>=sensitivity:hi=mid
        else:lo=mid+1
    return lo
