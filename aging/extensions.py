"""Frozen secondary models. Every resource/rule is supplied and hash-pinned."""
import collections
import json
import math
from pathlib import Path
from .common import *

def checked_plan(path):
    p=load_config(path)
    if p.get('status')!='FROZEN' or not p.get('reviewer'):raise ValueError('secondary plan not frozen')
    for entry in p['resources'].values():
        if sha(entry['path'])!=entry['sha256']:raise ValueError('secondary resource mismatch')
    return p

def bh(values):
    order=sorted(range(len(values)),key=lambda i:values[i]);out=[None]*len(values);q=1.
    for rank in range(len(order)-1,-1,-1):
        i=order[rank];q=min(q,values[i]*len(values)/(rank+1));out[i]=q
    return out

def signature_fit(spectrum,plan,output):
    import numpy as np
    from scipy.optimize import nnls
    p=checked_plan(plan);d=load_config(spectrum)
    if d['snvs']<p['min_snvs']:raise ValueError('too few SNVs for frozen signature model')
    matrix=read_tsv(p['resources']['signatures']['path']);names=[k for k in matrix[0] if k!='context']
    if len(matrix)!=96 or len({r['context'] for r in matrix})!=96:raise ValueError('signature matrix must have 96 unique contexts')
    if set(d['sbs96_counts'])-set(r['context'] for r in matrix):raise ValueError('unknown SBS96 context')
    x=np.array([[float(r[k]) for k in names] for r in matrix]);y=np.array([d['sbs96_counts'].get(r['context'],0) for r in matrix])
    if np.any(x<0) or not np.allclose(x.sum(axis=0),1,atol=1e-5):raise ValueError('signature columns must sum to one')
    fit,residual=nnls(x,y);total=fit.sum()
    if total<=0:raise ValueError('empty signature fit')
    rng=np.random.default_rng(p['seed']);boot=[]
    for _ in range(p['bootstrap_reps']):
        b=rng.multinomial(int(y.sum()),y/y.sum());w,_=nnls(x,b);boot.append(w/w.sum() if w.sum() else np.zeros(len(names)))
    boot=np.array(boot);fractions=fit/total
    rows=[dict(signature=k,exposure=float(fit[i]),fraction=float(fractions[i]),
          lower95=float(np.quantile(boot[:,i],.025)),upper95=float(np.quantile(boot[:,i],.975))) for i,k in enumerate(names)]
    if set(p['artifact_signatures'])-set(names):raise ValueError('artifact signature not in basis')
    artifact=sum(fractions[names.index(k)] for k in p['artifact_signatures'])
    write_json(output,dict(status='MODEL_BASED_SIGNATURE_FIT',snvs=int(y.sum()),residual=float(residual),
        exposures=rows,model_artifact_fraction=float(artifact),plan_sha256=sha(plan),
        note='Artifact labels must be validated on technical negatives. Non-clock signatures are not automatically artifacts.'))

def mosaic_candidates(segments,baf,plan,output):
    """Calibrated copy-ratio/BAF screen; candidate detection, not clonal diagnosis."""
    import numpy as np
    p=checked_plan(plan);sites=read_tsv(baf);segments=read_tsv(segments);rows=[]
    from .contracts import bed_intervals
    tp53=bed_intervals(p['resources']['tp53_bed']['path'])
    required=['min_probes','min_het_sites','min_baf_depth','min_segment_bp','log2_loss','log2_gain','baf_deviation']
    if any(k not in p for k in required):raise ValueError('mosaic screen thresholds missing')
    for r in segments:
        ch=r['chromosome'];start=int(r['start']);end=int(r['end']);log2=float(r['log2'])
        hs=[s for s in sites if s['chrom']==ch and start<=int(s['pos'])-1<end and
            int(s['dp'])>=p['min_baf_depth'] and s['evidence_status']=='OBSERVED']
        deviation=float(np.median([abs(float(s['vaf'])-.5) for s in hs])) if hs else None
        state='NO_SCREEN_SIGNAL'
        if end-start<p['min_segment_bp'] or int(r['probes'])<p['min_probes'] or len(hs)<p['min_het_sites']:
            state='NOT_EVALUABLE'
        elif log2<=p['log2_loss']:state='LOSS_WITH_BAF' if deviation>=p['baf_deviation'] else 'LOSS_COPY_RATIO_ONLY'
        elif log2>=p['log2_gain']:state='GAIN_WITH_BAF' if deviation>=p['baf_deviation'] else 'GAIN_COPY_RATIO_ONLY'
        elif deviation>=p['baf_deviation']:state='COPY_NEUTRAL_LOH_CANDIDATE'
        rows.append(dict(chrom=ch,start=start,end=end,log2=log2,n_het_sites=len(hs),baf_deviation=deviation,
                         screen=state,tp53_overlap=int(any(ch==cc and start<b and end>a for cc,a,b in tp53))))
    write_tsv(output,rows,['chrom','start','end','log2','n_het_sites','baf_deviation','screen','tp53_overlap'])
    write_json(str(output)+'.audit.json',dict(status='CALIBRATED_SCREEN_CANDIDATES_REQUIRE_REVIEW',plan_sha256=sha(plan),
        segments=len(rows),candidate_bp=sum(r['end']-r['start'] for r in rows if r['screen'] not in ['NO_SCREEN_SIGNAL','NOT_EVALUABLE']),
        not_evaluable_bp=sum(r['end']-r['start'] for r in rows if r['screen']=='NOT_EVALUABLE'),
        note='TP53 territory is hash-pinned; no validated cell fraction or orthogonal confirmation implied.'))

def coverage_index(intervals):
    index={}
    for ch,s,e in sorted(intervals):
        starts,ends,prefix=index.setdefault(ch,([],[],[0]))
        if ends and s<ends[-1]:raise ValueError('overlapping callable intervals')
        starts.append(s);ends.append(e);prefix.append(prefix[-1]+e-s)
    return index

def covered_bp(index,ch,s,e):
    from bisect import bisect_left,bisect_right
    if ch not in index:return 0
    starts,ends,prefix=index[ch]
    left=bisect_right(ends,s);right=bisect_left(starts,e)
    if right<=left:return 0
    return prefix[right]-prefix[left]-max(0,s-starts[left])-max(0,ends[right-1]-e)

def regional_burdens(feature_dir,regions_bed,plan,output):
    """Gene or pathway territories are merged separately, preventing double-counting exons."""
    from .analysis import model
    p=checked_plan(plan)
    if p['resources']['regions']['path']!=str(regions_bed):raise ValueError('region plan/path differs')
    regions=collections.defaultdict(list)
    for line in Path(regions_bed).read_text().splitlines():
        if not line or line.startswith('#'):continue
        ch,s,e,name=line.split()[:4];s,e=int(s),int(e)
        if s<0 or e<=s:raise ValueError('invalid region')
        safe_id(name);regions[name].append((ch,s,e))
    if not regions:raise ValueError('no gene/pathway territories')
    merged={}
    for name,rs in regions.items():
        m=[]
        for ch,s,e in sorted(rs):
            if m and m[-1][0]==ch and s<=m[-1][2]:m[-1]=(ch,m[-1][1],max(e,m[-1][2]))
            else:m.append((ch,s,e))
        merged[name]=m
    feature_dir=Path(feature_dir);out=Path(output);out.mkdir(parents=True,exist_ok=True)
    samples=read_tsv(feature_dir/'features.private.tsv');byregion=collections.defaultdict(list)
    for sample in samples:
        sid=safe_id(sample['sample_id']);ev=read_tsv(feature_dir/(sid+'.classified.tsv'))
        callable_rows=[]
        for line in (feature_dir/(sid+'.callable.bed')).read_text().splitlines():
            ch,s,e=line.split();callable_rows.append((ch,int(s),int(e)))
        ci=coverage_index(callable_rows)
        for name,rs in merged.items():
            # Regions and callable blocks are merged/non-overlapping within each set.
            bp=sum(covered_bp(ci,ch,s,e) for ch,s,e in rs)
            def inside(v):return any(v['chrom']==ch and s<=int(v['pos'])-1 and int(v['pos'])-1+len(v['ref'])<=e for ch,s,e in rs)
            n=sum(v['confidence']=='PASS' and inside(v) for v in ev)
            pn=sum(v['confidence_retain_pon']=='PASS' and inside(v) for v in ev)
            byregion[name].append(dict(sample,somatic_count=n,retain_pon_count=pn,callable_mb=bp/1e6,
                endpoint_status='EVALUABLE' if bp else 'NO_CALLABLE_TERRITORY'))
    stats=[]
    for name,rows in byregion.items():
        rd=out/name;rd.mkdir(exist_ok=True);write_tsv(rd/'features.private.tsv',rows)
        fits=model(rd/'features.private.tsv',p,rd)
        primary=fits[0];stats.append(dict(region=name,**primary))
    valid=[r for r in stats if 'p' in r]
    for r,q in zip(valid,bh([r['p'] for r in valid])):r['q_bh']=q
    write_json(out/'regional.aggregate.json',dict(status='EXPLORATORY_REGIONAL_BURDEN_NOT_POSITIVE_SELECTION_TEST',
        regions_tested=len(stats),models=stats,plan_sha256=sha(plan)))
