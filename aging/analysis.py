"""Frozen burden definitions and within-stratum age x R337H models."""
import collections
import json
import math
from pathlib import Path
from .common import *
from .contracts import bed_intervals
from .evidence import callable_threshold, classify

def features(samples,c,g,work,output):
    import pysam
    from .evidence import pileup_site
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    rows=[]
    for s in samples:
        sid=s['sample_id'];a=g['assays'][s['assay_id']];floor=a['exome_vaf_floor'] if a['branch']=='A' else a['driver_vaf_floor']
        threshold=callable_threshold(floor,a['min_alt'],a['sensitivity'])
        bed=c['resources']['common_bed' if a['branch']=='A' else 'driver_bed']['path']
        # Identical fragment filter to allele evidence. Enumerate every target base, including absent coverage.
        # Sparse covered-position aggregation avoids a dense per-base output on disk.
        callable_bp=0;total_bp=0;callable_regions=[]
        with pysam.AlignmentFile(s['bam_path'],'rb',index_filename=s['bai_path']) as bam, \
             pysam.FastaFile(c['resources']['fasta']['path']) as fasta:
            for ch,start,end in bed_intervals(bed):
                total_bp+=end-start
                for col in bam.pileup(ch,start,end,truncate=True,stepper='samtools',fastafile=fasta,
                    min_base_quality=c['evidence']['min_baseq'],min_mapping_quality=c['evidence']['min_mapq'],
                    max_depth=c['evidence']['max_depth'],flag_filter=0xF04,ignore_overlaps=False,
                    ignore_orphans=False,compute_baq=False):
                    if col.nsegments>=c['evidence']['max_depth']:continue
                    obs={}
                    for p in col.pileups:
                        if p.is_del or p.is_refskip or p.query_position is None:continue
                        read=p.alignment;q=p.query_position
                        if read.query_qualities[q]<c['evidence']['min_baseq']:continue
                        v=read.query_sequence[q]
                        if read.query_name in obs and obs[read.query_name]!=v:obs[read.query_name]=None
                        else:obs.setdefault(read.query_name,v)
                    if sum(v is not None for v in obs.values())>=threshold:
                        callable_bp+=1;p=col.reference_pos
                        if callable_regions and callable_regions[-1][0]==ch and callable_regions[-1][2]==p:
                            callable_regions[-1]=(ch,callable_regions[-1][1],p+1)
                        else:callable_regions.append((ch,p,p+1))
        atomic(output/(sid+'.callable.bed'),''.join(f'{ch}\t{s}\t{e}\n' for ch,s,e in callable_regions))
        ev=read_tsv(Path(work)/sid/'evidence.tsv');classified=[];counts={'primary':0,'retain_pon':0};vafs=[]
        for r in ev:
            tier=classify(r,g['tier'],floor,a['min_alt'])
            sensitivity=classify(r,g['tier'],floor,a['min_alt'],True)
            if int(r['dp'])<threshold:tier+=';BELOW_CALLABLE_DEPTH';sensitivity+=';BELOW_CALLABLE_DEPTH'
            if a['branch']=='B' and r['in_driver_territory']!='1':tier+=';OUTSIDE_DRIVER';sensitivity+=';OUTSIDE_DRIVER'
            if tier=='PASS':counts['primary']+=1;vafs.append(float(r['vaf']))
            if sensitivity=='PASS':counts['retain_pon']+=1
            classified.append(dict(r,confidence=tier,confidence_retain_pon=sensitivity))
        write_tsv(output/(sid+'.classified.tsv'),classified,list(classified[0]) if classified else ['sample_id','variant_id','confidence'])
        rows.append(dict(s,somatic_count=counts['primary'],retain_pon_count=counts['retain_pon'],
            callable_mb=callable_bp/1e6,target_bp=total_bp,callable_depth=threshold,vaf_floor=floor,branch=a['branch'],
            max_vaf=max(vafs) if vafs else '',vaf_weighted_burden=sum(vafs),
            burden_per_mb=counts['primary']/(callable_bp/1e6) if callable_bp else '',
            endpoint_status='EVALUABLE' if callable_bp else 'NO_CALLABLE_TERRITORY'))
    write_tsv(output/'features.private.tsv',rows)
    return rows

def model(features_path,g,output):
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    df=pd.read_csv(features_path,sep='\t');df=df[df.endpoint_status=='EVALUABLE'].copy()
    df=df[df.group.isin(['LFS_R337H','TP53_negative'])].copy()
    for k in ['age_at_collection','PC1','PC2','PC3','callable_mb','somatic_count']:
        df[k]=pd.to_numeric(df[k],errors='raise')
    df['age_decade']=df.age_at_collection/10;df['r337h']=(df.group=='LFS_R337H').astype(int)
    df['exposed']=(df.treatment_status=='EXPOSED').astype(int)
    results=[]
    scenarios=[('primary',df[df.treatment_status=='NAIVE'],'somatic_count'),
               ('all_treatment',df,'somatic_count'),('age_lt30',df[(df.treatment_status=='NAIVE')&(df.age_at_collection<30)],'somatic_count'),
               ('exclude_lt5',df[(df.treatment_status=='NAIVE')&(df.age_at_collection>=5)],'somatic_count'),
               ('retain_pon',df[df.treatment_status=='NAIVE'],'retain_pon_count')]
    for label,d,outcome in scenarios:
        d=d.copy();keep=[]
        for st,part in d.groupby('stratum'):
            controls=part[part.r337h==0];cases=part[part.r337h==1]
            if len(controls)<g['min_controls_per_stratum'] or len(cases)==0:continue
            # Restrict carrier ages to observed control age support in each exact stratum.
            cases=cases[cases.age_at_collection.between(controls.age_at_collection.min(),controls.age_at_collection.max())]
            if not cases.empty:keep+=list(controls.index)+list(cases.index)
        d=d.loc[keep].copy()
        base=dict(scenario=label,n=len(d),carriers=int(d.r337h.sum()),controls=int((d.r337h==0).sum()),strata=d.stratum.nunique())
        if not len(d) or d[outcome].sum()==0 or d.stratum.nunique()<g['min_model_strata']:
            results.append(dict(base,status='NOT_ESTIMABLE'));continue
        formula=outcome+' ~ age_decade*r337h + C(sex) + PC1 + PC2 + PC3 + C(cancer_history) + C(stratum)'
        if label=='all_treatment':formula+=' + exposed'
        try:
            m=smf.glm(formula,data=d,family=sm.families.Poisson(),offset=np.log(d.callable_mb),missing='raise')
            if np.linalg.matrix_rank(m.exog)<m.exog.shape[1] or len(d)<=m.exog.shape[1]:
                results.append(dict(base,status='NOT_ESTIMABLE_RANK'));continue
            fit=m.fit(cov_type='cluster',cov_kwds={'groups':d.stratum})
            if not fit.converged:raise ValueError('model did not converge')
            term='age_decade:r337h';ci=fit.conf_int().loc[term]
            vals=[float(fit.params[term]),float(fit.pvalues[term]),float(ci.iloc[0]),float(ci.iloc[1])]
            if not all(math.isfinite(v) for v in vals):raise ValueError('nonfinite model result')
            results.append(dict(base,status='ESTIMATED_EXPLORATORY_UNTIL_VALIDATED',log_rate_ratio_per_decade=vals[0],
                rate_ratio_per_decade=math.exp(vals[0]),ci_low=math.exp(vals[2]),ci_high=math.exp(vals[3]),p=vals[1],
                dispersion=float(fit.pearson_chi2/fit.df_resid),formula=formula))
        except (ValueError,np.linalg.LinAlgError) as e:
            results.append(dict(base,status='NOT_ESTIMABLE',reason=type(e).__name__))
    # Secondary sensitivity p-values are descriptive; primary is the sole confirmatory contrast.
    write_json(output/'models.aggregate.json',dict(models=results,
        limitations=['cross-sectional survivorship','tumor-only candidate evidence','model-based specificity unless orthogonally validated']))
    return results

def power(pilot_features,output,seed=1701,reps=1000,design=None):
    """Pilot-based parametric Poisson power grid; variance inflation estimated from pilot."""
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    d=pd.read_csv(pilot_features,sep='\t');d=d[(d.group.isin(['LFS_R337H','TP53_negative'])) & (d.callable_mb>0)].copy()
    if len(d)<10 or d.somatic_count.sum()<5:raise ValueError('pilot insufficient for power simulation')
    d['age_decade']=d.age_at_collection/10;d['r337h']=(d.group=='LFS_R337H').astype(int)
    form='somatic_count ~ age_decade*r337h + C(stratum)'
    m=smf.glm(form,d,family=sm.families.Poisson(),offset=np.log(d.callable_mb))
    if np.linalg.matrix_rank(m.exog)<m.exog.shape[1] or len(d)<=m.exog.shape[1]:raise ValueError('pilot power design not identifiable')
    fit=m.fit();disp=max(1,float(fit.pearson_chi2/fit.df_resid));rng=np.random.default_rng(seed);rows=[]
    null=fit.params.copy();null['age_decade:r337h']=0
    baseline=np.exp(m.exog@null+np.log(d.callable_mb))
    scope='PILOT_DESIGN'
    if design:
        import patsy
        pilot=d.copy()
        target=pd.read_csv(design,sep='\t')
        target=target[(target.group.isin(['LFS_R337H','TP53_negative'])) & (target.treatment_status=='NAIVE')].copy()
        if target.empty or target.person_id.duplicated().any():raise ValueError('invalid final power design')
        target['callable_mb']=pd.to_numeric(target.expected_callable_mb,errors='raise')
        if (target.callable_mb<=0).any():raise ValueError('expected callability must be positive')
        target['age_decade']=target.age_at_collection/10;target['r337h']=(target.group=='LFS_R337H').astype(int)
        basefit=smf.glm('somatic_count ~ age_decade',pilot,family=sm.families.Poisson(),offset=np.log(pilot.callable_mb)).fit()
        baseline=basefit.predict(target,offset=np.log(target.callable_mb))
        exog=patsy.dmatrix('age_decade*r337h + C(stratum)',target,return_type='dataframe')
        if np.linalg.matrix_rank(exog)<exog.shape[1] or len(target)<=exog.shape[1]:raise ValueError('final power design not identifiable')
        from types import SimpleNamespace
        m=SimpleNamespace(exog=np.asarray(exog),exog_names=list(exog.columns))
        d=target;scope='FINAL_DESIGN_WITH_EXPECTED_CALLABILITY'
    for beta in [0,0.1,0.2,0.3,0.5]:
        mu=baseline*np.exp(beta*d.age_decade*d.r337h);detected=0;valid=0
        for _ in range(reps):
            lam=rng.gamma(mu/(disp-1),disp-1) if disp>1 else mu
            y=rng.poisson(lam)
            try:
                f=sm.GLM(y,m.exog,family=sm.families.Poisson(),offset=np.log(d.callable_mb)).fit(scale='X2')
                p=f.pvalues[list(m.exog_names).index('age_decade:r337h')]
                if math.isfinite(p):valid+=1;detected+=int(p<0.05)
            except (ValueError,np.linalg.LinAlgError):pass
        rows.append(dict(beta=beta,valid=valid,reps=reps,power=detected/valid if valid else None))
    write_json(output,dict(status='PILOT_POWER_REVIEW_REQUIRED',scope=scope,seed=seed,dispersion=disp,n=len(d),grid=rows,
        pilot_sha256=sha(pilot_features),design_sha256=sha(design) if design else None,
        note='Parametric planning simulation, pilot-estimated variance. Expected callability and omitted batch-rate heterogeneity require review.'))
