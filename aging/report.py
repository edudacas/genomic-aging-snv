"""Aggregate report contains no patient, sample, BAM, or run identifiers."""
import collections
import html
import json
from pathlib import Path
from .common import *

def report(work,output):
    work=Path(work);sections=[]
    for name in ['preflight.json','QC_COMPLETE.json','PILOT_CALLS_COMPLETE.json','COHORT_CALLS_COMPLETE.json','COMPUTATIONAL_RUN_COMPLETE.json']:
        p=work/name
        if p.exists():
            d=json.loads(p.read_text())
            clean={k:v for k,v in d.items() if k in ['status','samples','candidate_sites','failures','not_claimed']}
            sections.append((name,clean))
    p=work/'features/features.private.tsv'
    if p.exists():
        rows=read_tsv(p);groups=collections.defaultdict(list)
        for r in rows:groups[r['group']].append(r)
        summaries={}
        for name,rs in groups.items():
            summaries[name]=dict(n=len(rs),evaluable=sum(r['endpoint_status']=='EVALUABLE' for r in rs),
                candidate_burden_events=sum(int(r['somatic_count']) for r in rs),
                treatment_naive=sum(r['treatment_status']=='NAIVE' for r in rs))
        sections.append(('Feature denominators',summaries))
    p=work/'analysis/models.aggregate.json'
    if p.exists():sections.append(('Models',json.loads(p.read_text())))
    if not sections:raise ValueError('no completed outputs to report')
    content='<!doctype html><meta charset="utf-8"><title>Genomic aging LFS</title><style>body{max-width:960px;margin:48px auto;font:16px system-ui;color:#172d35}pre{white-space:pre-wrap;background:#f1f4f5;padding:20px}h1{font-size:28px}</style><h1>Genomic aging in LFS</h1><p>Computational evidence; biological conclusions require calibration and scientific review.</p>'
    for label,data in sections:content+='<h2>'+html.escape(label)+'</h2><pre>'+html.escape(json.dumps(data,indent=2))+'</pre>'
    atomic(output,content)
