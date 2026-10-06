"""SNV fragment evidence only. No phenotype inputs, endpoint or somatic classifier."""
from collections import Counter, defaultdict
from statistics import median

EXCLUDED_FLAGS=0xF04
COUNT_FIELDS=('usable_fragment_depth','ref_fragments','alt_fragments','other_fragments',
              'alt_forward','alt_reverse','alt_strand_unknown','alt_f1r2','alt_f2r1','alt_unpaired','alt_orientation_unknown')

def strand(o):
    if not o['paired']:return 'R' if o['reverse'] else 'F'
    if o['read_number']==1:return 'R' if o['reverse'] else 'F'
    if o['read_number']==2 and not o['mate_unmapped']:return 'R' if o['mate_reverse'] else 'F'
    return 'UNKNOWN'

def orientation(o):
    if not o['paired']:return 'UNPAIRED'
    if o['read_number'] not in (1,2):return 'UNKNOWN'
    return 'F1R2' if ((o['read_number']==1)!=o['reverse']) else 'F2R1'

def unavailable(status,**diagnostics):
    return dict(evaluation_status=status,**{k:None for k in COUNT_FIELDS},
                vaf_alt_over_ref_alt=None,vaf_alt_over_usable=None,**diagnostics)

def reduce_observations(observations,ref,alt,min_bq=30,min_mq=30,saturated=False):
    """Reduce eligible observations by (RG,QNAME), preserving raw discordance diagnostics.

    Strand anchors to read1 (or read2's mate flag), never the highest-BQ mate.
    BQ/MQ summaries use the minimum across eligible observations per fragment.
    Empty RG is a distinct namespace and explicitly diagnosed, not silently inferred.
    """
    if min_bq<0 or min_mq<0:raise ValueError('Invalid quality thresholds')
    obs=list(observations);raw=defaultdict(set);groups=defaultdict(list);diag=Counter()
    for o in obs:
        key=(o['rg'],o['qname']);raw[key].add(o['base']);diag['raw_read_observations']+=1
        if not o['rg']:diag['reads_missing_rg']+=1
        if o['flag'] & EXCLUDED_FLAGS:diag['excluded_flag_reads']+=1;continue
        if o['mq'] is None or o['mq']==255:diag['unknown_mq_reads']+=1;continue
        if o['mq']<min_mq:diag['low_mq_reads']+=1;continue
        if o['bq'] is None:diag['unknown_bq_reads']+=1;continue
        if o['bq']<min_bq:diag['low_bq_reads']+=1;continue
        groups[key].append(o)
    diag['raw_discordant_fragments']=sum(len(b)>1 for b in raw.values())
    if saturated:return unavailable('DEPTH_SATURATED',saturation_state='CAPPED_NON_EVALUABLE',diagnostics=dict(diag))
    evidence={k:0 for k in COUNT_FIELDS};quality=defaultdict(list);discordant=0;strand_conflicts=0;orientation_conflicts=0
    for values in groups.values():
        bases={o['base'] for o in values}
        if len(bases)!=1:discordant+=1;continue
        base=next(iter(bases));kind='ref' if base==ref else 'alt' if base==alt else 'other'
        evidence['usable_fragment_depth']+=1;evidence[kind+'_fragments']+=1
        quality[kind+'_bq'].append(min(o['bq'] for o in values));quality[kind+'_mq'].append(min(o['mq'] for o in values))
        quality[kind+'_distance_to_end'].append(min(o['distance_to_end'] for o in values))
        if kind=='alt':
            directions={strand(o) for o in values}-{ 'UNKNOWN' }
            direction=next(iter(directions)) if len(directions)==1 else 'UNKNOWN'
            if len(directions)>1:strand_conflicts+=1
            evidence[{'F':'alt_forward','R':'alt_reverse','UNKNOWN':'alt_strand_unknown'}[direction]]+=1
            orientations={orientation(o) for o in values}-{'UNKNOWN'}
            orient=next(iter(orientations)) if len(orientations)==1 else 'UNKNOWN'
            if len(orientations)>1:orientation_conflicts+=1
            evidence[{'F1R2':'alt_f1r2','F2R1':'alt_f2r1','UNPAIRED':'alt_unpaired','UNKNOWN':'alt_orientation_unknown'}[orient]]+=1
    dp=evidence['usable_fragment_depth'];ra=evidence['ref_fragments']+evidence['alt_fragments'];a=evidence['alt_fragments']
    assert sum(evidence[k] for k in ('ref_fragments','alt_fragments','other_fragments'))==dp
    diag.update(eligible_discordant_fragments=discordant,fragment_strand_conflicts=strand_conflicts,
                fragment_orientation_conflicts=orientation_conflicts)
    summaries={}
    for kind in ('ref','alt','other'):
        for field in ('bq','mq','distance_to_end'):
            values=quality[kind+'_'+field]
            summaries[kind+'_'+field+'_min']=min(values) if values else None
            summaries[kind+'_'+field+'_median']=median(values) if values else None
        ends=quality[kind+'_distance_to_end'];summaries[kind+'_near_end_fraction']=sum(x<10 for x in ends)/len(ends) if ends else None
    state='EVALUATED' if dp else 'NO_COVERAGE' if not obs else 'NO_USABLE_FRAGMENTS'
    return dict(evaluation_status=state,saturation_state='NOT_SATURATED',**evidence,
        vaf_alt_over_ref_alt=a/ra if ra else None,vaf_alt_over_usable=a/dp if dp else None,
        diagnostics=dict(diag),quality_summaries=summaries)

def recount_site(bam,fasta,site,min_bq=30,min_mq=30,raw_max_depth=100000):
    """Raw pointwise pysam pileup; complete counts and diagnostics share one traversal.

    Saturation is conservatively non-evaluable, with no inconsistent fallback.
    Only SNVs are supported in this pilot release; indels/MNVs remain explicit.
    """
    if raw_max_depth<1:raise ValueError('Positive raw depth cap required')
    c,pos,ref,alt=site['chrom'],int(site['pos']),site['ref'],site['alt']
    if len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT' or ref==alt:
        return unavailable('UNSUPPORTED_ALLELE',saturation_state='NOT_EVALUATED')
    if pos<1:return unavailable('NOT_EVALUATED',reason='INVALID_COORDINATE')
    try:actual=fasta.fetch(c,pos-1,pos).upper()
    except (KeyError,ValueError):return unavailable('REF_MISMATCH',reason='CONTIG_OR_COORDINATE_ABSENT')
    if actual!=ref:return unavailable('REF_MISMATCH',observed_reference=actual)
    if c not in bam.references:return unavailable('NOT_EVALUATED',reason='BAM_CONTIG_ABSENT')
    observations=[];saturated=False;skipped=Counter()
    for col in bam.pileup(c,pos-1,pos,truncate=True,stepper='nofilter',min_base_quality=0,
            min_mapping_quality=0,max_depth=raw_max_depth,ignore_overlaps=False,ignore_orphans=False,compute_baq=False):
        if col.reference_pos!=pos-1:continue
        saturated=col.nsegments>=raw_max_depth
        for p in col.pileups:
            q=p.query_position
            if q is None or p.is_del or p.is_refskip:skipped['deletion_or_refskip_reads']+=1;continue
            r=p.alignment
            if not r.query_sequence or q>=len(r.query_sequence):skipped['missing_sequence_reads']+=1;continue
            base=r.query_sequence[q].upper() if p.indel==0 else 'OTHER'
            observations.append(dict(rg=r.get_tag('RG') if r.has_tag('RG') else '',qname=r.query_name,
                base=base,flag=r.flag,bq=r.query_qualities[q] if r.query_qualities is not None else None,
                mq=r.mapping_quality,paired=r.is_paired,read_number=1 if r.is_read1 else 2 if r.is_read2 else 0,
                reverse=r.is_reverse,mate_reverse=r.mate_is_reverse,mate_unmapped=r.mate_is_unmapped,
                distance_to_end=min(q,len(r.query_sequence)-1-q)))
    result=reduce_observations(observations,ref,alt,min_bq,min_mq,saturated)
    result.setdefault('diagnostics',{}).update(skipped)
    if result['evaluation_status']=='NO_COVERAGE' and skipped:
        result['evaluation_status']='NO_USABLE_FRAGMENTS'
    return result
