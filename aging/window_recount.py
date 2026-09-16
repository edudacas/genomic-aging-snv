"""Bounded pileup windows using the same fragment/allele reducer as the baseline."""
from collections import defaultdict
from .evidence import pileup_site

class ColumnView:
    def __init__(self,col):self.col=col
    def pileup(self,*args,**kwargs):return () if self.col is None else (self.col,)

def ordered_sites(candidates,fasta):
    rank={c:i for i,c in enumerate(fasta.references)}
    rows=[dict(v,variant_id=k) for k,v in candidates.items()]
    for r in rows:
        if r['chrom'] not in rank:raise ValueError('Unknown candidate contig')
    return sorted(rows,key=lambda r:(rank[r['chrom']],int(r['pos']),r['ref'],r['alt']))

def windows(rows,max_span=10000,max_gap=1000):
    chunk=[]
    for row in rows:
        if chunk and (row['chrom']!=chunk[0]['chrom'] or row['pos']-chunk[0]['pos']>=max_span or row['pos']-chunk[-1]['pos']>max_gap):
            yield chunk;chunk=[]
        chunk.append(row)
    if chunk:yield chunk

def recount_window(bam,fasta,rows,settings):
    bq=settings['min_baseq'];mq=settings['min_mapq'];cap=settings['max_depth']
    bypos=defaultdict(list);result={};saturated=False
    for r in rows:bypos[r['pos']-1].append(r)
    for col in bam.pileup(rows[0]['chrom'],rows[0]['pos']-1,rows[-1]['pos'],truncate=True,
                         stepper='samtools',fastafile=fasta,min_base_quality=bq,min_mapping_quality=mq,
                         max_depth=cap,flag_filter=0xF04,ignore_overlaps=False,ignore_orphans=False,compute_baq=False):
        saturated |= col.nsegments>=cap
        for row in bypos.get(col.reference_pos,[]):
            result[row['variant_id']]=pileup_site(ColumnView(col),fasta,row,bq,mq,cap)
    # Depth truncation can depend on pileup start. Preserve baseline behavior on the whole window.
    for row in rows:
        k=row['variant_id']
        if saturated:result[k]=pileup_site(bam,fasta,row,bq,mq,cap)
        elif k not in result:result[k]=pileup_site(ColumnView(None),fasta,row,bq,mq,cap)
        if result[k]['evidence_status']=='UNSUPPORTED_COMPLEX_ALLELE':
            result[k].update(dp='',ad_ref='',ad_alt='',vaf='')
    return result,saturated
