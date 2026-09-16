"""Offline exome-wide annotation; no assumption that PARQ180 covers the exome."""
from pathlib import Path
from .common import *
from .extensions import checked_plan

def annotate(candidates,plan,output):
    import pysam
    p=checked_plan(plan);out=Path(output);out.mkdir(parents=True,exist_ok=True)
    cache=Path(p['cache_dir']).resolve(strict=True)
    manifest=read_tsv(p['resources']['cache_manifest']['path'])
    if not manifest:raise ValueError('empty annotation cache manifest')
    for r in manifest:
        f=(cache/r['relative_path']).resolve(strict=True)
        if not f.is_relative_to(cache) or sha(f)!=r['sha256']:raise ValueError('VEP cache changed')
    rows=read_tsv(candidates);header=pysam.VariantHeader();header.add_meta('fileformat','VCFv4.2')
    with pysam.FastaFile(p['resources']['fasta']['path']) as fa:
        for ch,n in zip(fa.references,fa.lengths):header.contigs.add(ch,length=n)
        with pysam.VariantFile(out/'candidates.vcf','w',header=header) as vcf:
            for r in rows:
                pos=int(r['pos']);ref=r['ref']
                if fa.fetch(r['chrom'],pos-1,pos-1+len(ref)).upper()!=ref:raise ValueError('annotation REF mismatch')
                rec=vcf.new_record(contig=r['chrom'],start=pos-1,stop=pos-1+len(ref),alleles=(ref,r['alt']))
                vcf.write(rec)
    run(['vep','--offline','--cache','--dir_cache',cache,'--cache_version',str(p['cache_version']),
         '--assembly',p['build'],'--fasta',p['resources']['fasta']['path'],'--format','vcf','--vcf',
         '--symbol','--canonical','--no_stats','--fork',str(p['threads']),
         '-i',out/'candidates.vcf','-o',out/'annotated.vcf'],out/'commands.log')
    annotations=[]
    with pysam.VariantFile(out/'annotated.vcf') as v:
        fields=v.header.info['CSQ'].description.split('Format: ')[-1].split('|')
        for r in v:
            for csq in r.info.get('CSQ',()):
                item=dict(zip(fields,csq.split('|')))
                annotations.append(dict(variant_id=f'{r.chrom}:{r.pos}:{r.ref}:{r.alts[0]}',
                    gene=item.get('SYMBOL',''),ensembl_gene=item.get('Gene',''),
                    consequence=item.get('Consequence',''),transcript=item.get('Feature',''),canonical=item.get('CANONICAL','')))
    write_tsv(out/'annotations.tsv',annotations,['variant_id','gene','ensembl_gene','consequence','transcript','canonical'])
    write_json(out/'annotation.audit.json',dict(status='OFFLINE_ANNOTATION_COMPLETE',input_variants=len(rows),
        annotation_rows=len(annotations),plan_sha256=sha(plan),input_sha256=sha(candidates)))
