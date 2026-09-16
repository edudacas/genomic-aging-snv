"""Controlled BAM perturbations and measured calibration evidence ingestion."""
import collections
import json
from pathlib import Path
from .common import *

BATTERY = {'spikein_snv','spikein_indel','downsampling','negative_regions',
           'recurrence','spectrum','het_balance','split_half'}

def resample(bam,output,fraction,seed):
    import pysam
    if not 0<fraction<=1:raise ValueError('fraction must be (0,1]')
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    if output.exists():raise ValueError('refuse overwrite derived BAM')
    # Hashing read name keeps mates and secondary records together, deterministically.
    with pysam.AlignmentFile(bam,'rb') as src,pysam.AlignmentFile(str(output)+'.part','wb',template=src) as dst:
        for r in src:
            u=int(hashlib.sha256(f'{seed}:{r.query_name}'.encode()).hexdigest()[:16],16)/2**64
            if u<fraction:dst.write(r)
    Path(str(output)+'.part').rename(output);pysam.index(str(output));pysam.quickcheck(str(output))
    write_json(str(output)+'.provenance.json',dict(input_sha256=sha(bam),output_sha256=sha(output),
        fraction=fraction,seed=seed,status='COMPUTATIONAL_SUBSAMPLING_NOT_TECHNICAL_REPLICATE'))

def split_half(bam,out,seed):
    import pysam
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    if (out/'half0.bam').exists() or (out/'half1.bam').exists():raise ValueError('refuse overwrite split halves')
    with pysam.AlignmentFile(bam,'rb') as src:
        with pysam.AlignmentFile(str(out/'half0.bam'),'wb',template=src) as a,pysam.AlignmentFile(str(out/'half1.bam'),'wb',template=src) as b:
            for r in src:
                target=int(hashlib.sha256(f'{seed}:{r.query_name}'.encode()).hexdigest()[:16],16)%2
                (a if target==0 else b).write(r)
    for i in [0,1]:pysam.index(str(out/f'half{i}.bam'))
    write_json(out/'split.json',dict(seed=seed,input_sha256=sha(bam),status='STRESS_TEST_NOT_REPRODUCIBILITY'))

def spike(bam,sites,fasta,output,kind,seed,threads):
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    if kind not in ['snv','indel']:raise ValueError('SNV or indel required')
    if output.exists():raise ValueError('refuse overwrite spike BAM')
    run(['addsnv.py' if kind=='snv' else 'addindel.py','-v',sites,'-f',bam,'-r',fasta,
         '-o',output,'--seed',str(seed),'-p',str(threads)],str(output)+'.log')
    run(['samtools','index',output],str(output)+'.log')
    write_json(str(output)+'.provenance.json',dict(input_sha256=sha(bam),sites_sha256=sha(sites),
        output_sha256=sha(output),kind=kind,seed=seed,
        note='Requested spikes are not truth. Use BAMSurgeon successful-spike VCF and verify achieved AD/DP.'))

def summarize(manifest,output):
    """Every row references actual measured evidence, whose bytes are frozen."""
    from scipy.stats import beta
    entries=read_tsv(manifest);seen=set();results=[]
    for r in entries:
        name=r['component'];assay=r['assay_id'];seen.add(name)
        if r['status']!='COMPLETE':raise ValueError('calibration component incomplete '+name)
        if not Path(r['path']).is_file() or sha(r['path'])!=r['sha256']:raise ValueError('calibration evidence mismatch')
        if not r['method'] or not r['reviewer']:raise ValueError('method/reviewer missing')
        out=dict(r)
        if name in ['spikein_snv','spikein_indel','downsampling']:
            measured=read_tsv(r['path'])
            if not measured:raise ValueError('empty sensitivity evidence')
            groups=collections.defaultdict(list)
            for m in measured:
                if m['detected'] not in ['0','1'] or m['successful_truth']!='1':raise ValueError('invalid achieved truth record')
                if not 0<float(m['achieved_vaf'])<0.5 or int(m['dp'])<=0:raise ValueError('invalid achieved VAF/depth')
                groups[(m['vaf_bin'],m['depth_bin'],m['context'])].append(int(m['detected']))
            out['bins']=[]
            for k,vals in groups.items():
                n=len(vals);a=sum(vals)
                out['bins'].append(dict(vaf_bin=k[0],depth_bin=k[1],context=k[2],n=n,detected=a,
                    sensitivity=a/n,lower95=float(beta.ppf(.025,a,n-a+1)) if a else 0,
                    upper95=float(beta.ppf(.975,a+1,n-a)) if a<n else 1))
        results.append(out)
    assays={r['assay_id'] for r in entries}
    for assay in assays:
        have={r['component'] for r in entries if r['assay_id']==assay}
        if not BATTERY<=have:raise ValueError('missing calibration components for '+assay+': '+','.join(sorted(BATTERY-have)))
    if not assays:raise ValueError('empty calibration battery')
    write_json(output,dict(status='CALIBRATION_REVIEW_READY',manifest_sha256=sha(manifest),
        evidence=results,assays=sorted(assays),
        note='No floor automatically approved. Specificity is model-based unless orthogonal validation is supplied.'))

def compare_truth(truth_vcf,evidence,output,vaf_bin,depth_bin,context):
    from .evidence import variants
    truth={k for k,r,a in variants(truth_vcf)}
    obs={r['variant_id']:r for r in read_tsv(evidence)}
    rows=[]
    for k in sorted(truth):
        if k not in obs:raise ValueError('successful spike missing forced evidence; include truth in union')
        r=obs[k]
        if number(r['vaf']) is None or float(r['vaf'])<=0:raise ValueError('spike not verified in BAM')
        rows.append(dict(variant_id=k,successful_truth=1,achieved_vaf=r['vaf'],dp=r['dp'],
            detected=int(r.get('confidence')=='PASS'),vaf_bin=vaf_bin,depth_bin=depth_bin,context=context))
    if not rows:raise ValueError('no successful spike truth')
    write_tsv(output,rows)
