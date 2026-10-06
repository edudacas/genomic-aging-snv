"""Audited legacy comparator only; exact two functions from aging/evidence.py.
Not a new classifier or cohort launcher. Upstream snapshot 2026-10-06.
"""
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
