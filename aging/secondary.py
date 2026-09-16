"""Secondary outcomes: TelSeq, copy-ratio/BAF evidence, and spectra. No clonal CNV diagnosis."""
import collections
from pathlib import Path
from .common import *

def secondary_sample(s,c,work,output):
    output=Path(output);output.mkdir(parents=True,exist_ok=True);log=output/'commands.log'
    from .runtime import bam_for_tools
    bam=bam_for_tools(s,output)
    state={}
    if c['secondary']['telseq']:
        run(['telseq','-k',str(c['secondary']['telseq_k']),'-o',output/'telseq.txt',bam],log)
        state['telseq']='RAW_COMPLETE_REQUIRES_RUN_NORMALIZATION'
    if c['secondary']['cnv']:
        r=c['assays'][s['assay_id']]['cnv_reference']
        if (r.get('sha256') or c.get('profile')!='exploratory') and sha(r['path'])!=r['sha256']:
            raise ValueError('CNV reference hash mismatch')
        run(['cnvkit.py','batch',bam,'-r',r['path'],'-d',output/'cnv',
             '-p',str(c['threads']),'--drop-low-coverage'],log)
        state['cnv']='COPY_RATIO_SEGMENTS_ONLY_MOSAIC_SENSITIVITY_NOT_VALIDATED'
        # Reinterrogate heterozygous germline sites from this same blood BAM for BAF.
        import pysam
        from .evidence import variants,pileup_site
        baf=[]
        with pysam.AlignmentFile(s['bam_path'],'rb',index_filename=s['bai_path']) as bam, \
             pysam.FastaFile(c['resources']['fasta']['path']) as fasta:
            germline=Path(work)/'germline.norm.vcf.gz'
            for k,r,a in (variants(germline) if germline.is_file() else []):
                if len(r.ref)!=1 or len(a)!=1:continue
                gt=next(iter(r.samples.values())).get('GT')
                if not gt or set(gt)!={0,1}:continue
                site=dict(variant_id=k,chrom=r.chrom,pos=r.pos,ref=r.ref,alt=a)
                baf.append(dict(site,**pileup_site(bam,fasta,site,c['evidence']['min_baseq'],c['evidence']['min_mapq'],c['evidence']['max_depth'])))
        write_tsv(output/'het_baf.tsv',baf,list(baf[0]) if baf else ['variant_id','chrom','pos','dp','ad_alt','vaf'])
        state['loh']='BAF_EVIDENCE_ONLY_REQUIRES_MOSAIC_CALIBRATION_AND_REVIEW'
        state['baf']='AVAILABLE' if germline.is_file() else 'NOT_AVAILABLE_NO_GERMLINE_VCF'
    write_json(output/'secondary.status.json',state)

def spectra(classified,fasta_path,output):
    import pysam
    counts=collections.Counter();eligible=0
    complement=str.maketrans('ACGT','TGCA')
    with pysam.FastaFile(fasta_path) as f:
        for r in read_tsv(classified):
            if r['confidence']!='PASS' or len(r['ref'])!=1 or len(r['alt'])!=1:continue
            p=int(r['pos']);context=f.fetch(r['chrom'],max(0,p-2),p+1).upper()
            if len(context)!=3 or any(x not in 'ACGT' for x in context):continue
            ref=r['ref'];alt=r['alt']
            if ref in 'AG':context=context.translate(complement)[::-1];ref=ref.translate(complement);alt=alt.translate(complement)
            counts[f'{context[0]}[{ref}>{alt}]{context[2]}']+=1;eligible+=1
    write_json(output,dict(sbs96_counts=dict(counts),snvs=eligible,
        status='DESCRIPTIVE_SPECTRUM_NO_SIGNATURE_OR_ARTIFACT_FRACTION_INFERENCE'))

def normalize_telseq(table,output):
    """Import either validated reused TelSeq or new parsed TelSeq; exact analysis IDs only."""
    import numpy as np
    rows=read_tsv(table);groups=collections.defaultdict(list)
    for r in rows:
        if not r.get('analysis_id') or not r.get('source_sha256') or r.get('status')!='VALIDATED':
            raise ValueError('unvalidated TelSeq import')
        if number(r.get('length_estimate')) is None or float(r['length_estimate'])<=0:raise ValueError('invalid TelSeq')
        groups[(r['assay_id'],r['batch_id'])].append(r)
    if len({r['analysis_id'] for r in rows})!=len(rows):raise ValueError('duplicate TelSeq analyses')
    for _,rs in groups.items():
        ctr=[float(r['length_estimate']) for r in rs if r['group']=='TP53_negative']
        med=float(np.median(ctr)) if ctr else None
        mad=float(np.median(np.abs(np.array(ctr)-med))) if ctr else None
        for r in rs:
            r['control_n']=len(ctr)
            r['run_normalized_z']=(float(r['length_estimate'])-med)/(1.4826*mad) if len(ctr)>=5 and mad and mad>0 else ''
            r['normalization_status']='OK' if r['run_normalized_z']!='' else 'INSUFFICIENT_CONTROLS_OR_ZERO_MAD'
    write_tsv(output,rows)
