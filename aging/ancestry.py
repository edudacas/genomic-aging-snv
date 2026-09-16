"""Cohort PCs from jointly genotyped common SNPs, with missing genotypes retained."""
from pathlib import Path
from .common import *
from .runtime import validate_bam,resource,bam_for_tools

def ancestry(samples,c,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);log=out/'commands.log'
    if len(samples)<10:raise ValueError('at least 10 samples needed for ancestry PCA')
    sms=[validate_bam(s,c) for s in samples]
    if len(set(sms))!=len(sms):raise ValueError('BAM SM names duplicated; explicit harmonization needed')
    mapping=dict(zip(sms,samples))
    atomic(out/'bams.private.list',''.join(bam_for_tools(s,out/safe_id(s['sample_id']))+'\n' for s in samples))
    run(['bcftools','mpileup','-f',resource(c,'fasta'),'-R',resource(c,'contamination_sites'),
         '-b',out/'bams.private.list','-q','30','-Q','20','-a','FORMAT/DP','-Ob','-o',out/'likelihoods.bcf'],log)
    run(['bcftools','call','-m','-v','-Ob','-o',out/'joint.bcf',out/'likelihoods.bcf'],log)
    run(['bcftools','view','-i','QUAL>=30','-Ob','-o',out/'quality.bcf',out/'joint.bcf'],log)
    run(['bcftools','+setGT',out/'quality.bcf','-Oz','-o',out/'common.vcf.gz','--','-t','q','-n','.','-i','FMT/DP<10'],log)
    run(['plink2','--vcf',out/'common.vcf.gz','--double-id','--snps-only','just-acgt','--max-alleles','2',
         '--autosome','--maf','0.05','--geno','0.02','--mind','0.05','--make-pgen','--out',out/'common',
         '--set-all-var-ids','@:#:$r:$a',
         '--threads',str(c['threads']),'--memory',str(c['java_heap_gb']*1024)],log)
    run(['plink2','--pfile',out/'common','--indep-pairwise','200','50','0.2','--out',out/'prune',
         '--threads',str(c['threads'])],log)
    if len((out/'prune.prune.in').read_text().splitlines())<1000:raise ValueError('fewer than 1000 independent common SNPs; PCA HOLD')
    run(['plink2','--pfile',out/'common','--extract',out/'prune.prune.in','--pca','10','--out',out/'pca',
         '--threads',str(c['threads'])],log)
    rows=[]
    lines=(out/'pca.eigenvec').read_text().splitlines();header=lines[0].split()
    for line in lines[1:]:
        r=dict(zip(header,line.split()));sid=r['IID']
        if sid not in mapping:raise ValueError('PCA sample mapping missing')
        rows.append(dict(cohort_person_id=mapping[sid]['person_id'],analysis_id=mapping[sid]['analysis_id'],
                         PC1=r['PC1'],PC2=r['PC2'],PC3=r['PC3']))
    write_tsv(out/'pcs.private.tsv',rows)
    write_json(out/'ancestry.json',dict(status='COHORT_PCS_REQUIRE_BATCH_CONFOUNDING_REVIEW',input_n=len(samples),output_n=len(rows),
        note='Cohort PCs, not ancestry proportions. Inspect capture/batch confounding and relatedness before clinical-table import.'))
