"""Single-worker Docker runtime. Inputs read-only; resource budget enforced by launcher."""
import collections
import json
import os
import shutil
from pathlib import Path
from .common import *
from .contracts import bed_intervals

def resource(c,name): return c['resources'][name]['path']

def bam_for_tools(row,d):
    """Canonical local links support registries whose index lives in another directory."""
    d=Path(d);d.mkdir(parents=True,exist_ok=True)
    bam=d/'input.bam'
    for link,target in [(bam,row['bam_path']),(d/'input.bam.bai',row['bai_path'])]:
        target=Path(target).resolve(strict=True)
        if link.is_symlink():
            if link.resolve()!=target:raise ValueError('staged input link changed')
        elif link.exists():raise ValueError('staged input already exists as regular file')
        else:link.symlink_to(target)
    return str(bam)

def validate_bam(row,c):
    import pysam
    bam=row['bam_path']; bai=row['bai_path']
    if not Path(bam).is_file() or not Path(bai).is_file():raise ValueError('BAM/BAI absent')
    if row.get('bam_sha256') and sha(bam)!=row['bam_sha256']:raise ValueError('BAM hash mismatch')
    pysam.quickcheck(bam)
    with pysam.AlignmentFile(bam,'rb',index_filename=bai) as b:
        if b.header.to_dict().get('HD',{}).get('SO')!='coordinate':raise ValueError('BAM not coordinate sorted')
        sm={r.get('SM') for r in b.header.to_dict().get('RG',[])}
        if len(sm)!=1 or None in sm:raise ValueError('BAM must contain exactly one RG sample')
        lengths=dict(zip(b.references,b.lengths))
        with pysam.FastaFile(resource(c,'fasta')) as f:
            ref=dict(zip(f.references,f.lengths))
        # Check every BAM contig, not only the exome target chromosomes.
        for chrom,length in lengths.items():
            if ref.get(chrom)!=length:
                raise ValueError('BAM/reference contig mismatch: '+chrom)
        dictionary={}
        for line in Path(resource(c,'dict')).read_text().splitlines():
            if line.startswith('@SQ\t'):
                sq=dict(field.split(':',1) for field in line.split('\t')[1:])
                dictionary[sq['SN']]=sq
        if {k:int(v['LN']) for k,v in dictionary.items()}!=ref:
            raise ValueError('FASTA dictionary differs from FASTA index')
        for sq in b.header.to_dict().get('SQ',[]):
            md5=dictionary.get(sq['SN'],{}).get('M5')
            if sq.get('M5') and md5 and sq['M5'].lower()!=md5.lower():
                raise ValueError('BAM/reference sequence MD5 mismatch: '+sq['SN'])
        for chrom,start,end in bed_intervals(resource(c,'common_bed')):
            if lengths.get(chrom)!=ref.get(chrom):raise ValueError('BAM/reference contig length mismatch')
        # Index validation requires an actual indexed query, not existence alone.
        chrom,start,end=bed_intervals(resource(c,'common_bed'))[0]
        next(b.fetch(chrom,start,end),None)
    return next(iter(sm))

def qc(row,c,d):
    import pysam
    d=Path(d);d.mkdir(parents=True,exist_ok=True);log=d/'commands.log'
    sm=validate_bam(row,c);bam=bam_for_tools(row,d)
    run(['samtools','stats','-@',str(c['threads']-1),bam],log,d/'stats.txt')
    sn={}
    for l in (d/'stats.txt').read_text().splitlines():
        if l.startswith('SN\t'):
            f=l.split('\t');sn[f[1].rstrip(':')]=number(f[2])
    bed=resource(c,'common_bed');total=sum(e-s for _,s,e in bed_intervals(bed))
    # depth's -q is BQ; -Q is MQ. -s avoids counting overlapping paired mates twice.
    run(['samtools','depth','-s','-q',str(c['evidence']['min_baseq']),'-Q',str(c['evidence']['min_mapq']),
         '-G','0xF04','-b',bed,bam],log,d/'depth.tsv')
    compact=c.get('compact_qc')
    coverage_receipt=None
    if compact:
        from .compact_qc import summarize
        mask=compact['mask_path']
        if sha(mask)!=compact['mask_sha256']:raise ValueError('Coverage mask changed')
        filters=dict(min_baseq=c['evidence']['min_baseq'],min_mapq=c['evidence']['min_mapq'],
                     flag_filter='0xF04',overlapping_mates='samtools depth -s',
                     missing_positions='ZERO_WITHIN_COMMON_BED')
        hist,coverage_receipt=summarize(d/'depth.tsv',bed,mask,d,filters,compact['thresholds'])
        if total!=coverage_receipt['territory_bases']:raise ValueError('Coverage territory mismatch')
    else:
        hist=collections.Counter();covered=0
        with (d/'depth.tsv').open() as f:
            for l in f:
                n=int(l.rstrip().split('\t')[2]);hist[n]+=1;covered+=1
        if covered>total:raise ValueError('depth exceeds BED territory')
        hist[0]+=total-covered
        write_tsv(d/'depth_histogram.tsv',[dict(depth=k,bases=v) for k,v in sorted(hist.items())])
    run(['samtools','view','-c','-F','0xF04',bam],log,d/'mapped_count.txt')
    capture=c.get('assays',{}).get(row['assay_id'],{}).get('capture_bed')
    if capture:run(['samtools','view','-c','-F','0xF04','-L',capture,bam],log,d/'on_target_count.txt')
    tags=collections.Counter();n=0
    with pysam.AlignmentFile(bam,'rb',index_filename=row['bai_path']) as b:
        pg=b.header.to_dict().get('PG',[])
        for read in b.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:continue
            n+=1
            for tag in ['RX','MI']:
                if read.has_tag(tag):tags[tag]+=1
            if n>=100000:break
    mapped=int((d/'mapped_count.txt').read_text());ontarget=int((d/'on_target_count.txt').read_text()) if capture else None
    q=dict(sample_id=row['sample_id'],bam_sample=sm,status='TECHNICAL_QC_COMPLETE',
        target_bp=total,mean_depth=sum(k*v for k,v in hist.items())/total,
        mapped_primary_nonduplicate_reads=mapped,on_target_fraction=ontarget/mapped if mapped and ontarget is not None else None,
        duplicate_fraction=(sn.get('reads duplicated',0)/sn['raw total sequences']) if sn.get('raw total sequences') else None,
        read_length=sn.get('average length'),insert_size=sn.get('insert size average'),
        umi_sampled_reads=n,umi_tag_counts=dict(tags),
        umi_status='TAGS_PRESENT_NOT_CONSENSUS_VALIDATED' if tags else 'NO_TAGS_IN_FIRST_100K_NOT_PROOF_OF_ABSENCE',
        program_records=pg,depth_zero_bp=hist[0],bam_path=row['bam_path'],bam_size=Path(row['bam_path']).stat().st_size,
        bam_mtime_ns=Path(row['bam_path']).stat().st_mtime_ns,evidence_filters=c['evidence'],common_bed_sha256=sha(bed))
    if coverage_receipt:
        q.update(coverage_summary_sha256=sha(d/'coverage.summary.json'),
                 x_auto_ratio=coverage_receipt['x_auto_ratio'],y_auto_ratio=coverage_receipt['y_auto_ratio'],
                 sex_inferred='NOT_ASSESSED_GATE_REQUIRED',per_base_depth_retained=False)
    write_json(d/'qc.json',q)
    if coverage_receipt:
        from .compact_qc import discard_verified_depth
        discard_verified_depth(d/'depth.tsv',d)
    return q

def normalize(vcf,out,c,log):
    run(['bcftools','norm','-f',resource(c,'fasta'),'-m','-any','-Oz','-o',out,vcf],log)
    run(['bcftools','index','-f','-t',out],log)

def call_sample(row,c,d,pon,pon_mode=False):
    d=Path(d);d.mkdir(parents=True,exist_ok=True);log=d/'commands.log'
    sm=validate_bam(row,c);bam=bam_for_tools(row,d);fasta=resource(c,'fasta');bed=resource(c,'common_bed')
    java=['gatk','--java-options',f'-Xmx{c["java_heap_gb"]}g -Djava.io.tmpdir={d}']
    exploratory=c.get('profile')=='exploratory'
    callers=c.get('callers',['mutect','lofreq','sage']) if exploratory else ['mutect','lofreq','sage']
    raw=d/'mutect.raw.vcf.gz'
    argv=java+['Mutect2','-R',fasta,'-I',bam,'-tumor',sm,'-L',bed,
        '--native-pair-hmm-threads',str(c['threads']),'--max-mnp-distance','0',
        '--f1r2-tar-gz',d/'f1r2.tar.gz','-O',raw]
    if c['resources'].get('gnomad'):argv+=['--germline-resource',resource(c,'gnomad')]
    if exploratory:argv+=['--genotype-germline-sites','true','--genotype-pon-sites','true']
    if pon:
        argv+=['--panel-of-normals',pon]
        if not exploratory:argv+=['--genotype-pon-sites','true','--genotype-germline-sites','true']
    if 'mutect' in callers:
        run(argv,log)
        if pon_mode:return raw
        run(java+['LearnReadOrientationModel','-I',d/'f1r2.tar.gz','-O',d/'orientation.tar.gz'],log)
        contamination_args=[]
        if c['resources'].get('contamination_sites'):
            run(java+['GetPileupSummaries','-I',bam,'-V',resource(c,'contamination_sites'),
                  '-L',resource(c,'contamination_sites'),'-O',d/'pileups.table'],log)
            run(java+['CalculateContamination','-I',d/'pileups.table','-O',d/'contamination.table',
                  '--tumor-segmentation',d/'segments.table'],log)
            contamination_args=['--contamination-table',d/'contamination.table','--tumor-segmentation',d/'segments.table']
        run(java+['FilterMutectCalls','-R',fasta,'-V',raw,'--stats',str(raw)+'.stats',
                  '--ob-priors',d/'orientation.tar.gz','-O',d/'mutect.filtered.vcf.gz']+contamination_args,log)
        normalize(d/'mutect.filtered.vcf.gz',d/'mutect.norm.vcf.gz',c,log)
    # LoFreq adds indel qualities to a private derived copy; source BAM remains read-only.
    if 'lofreq' in callers:
        run(['lofreq','indelqual','--dindel','-f',fasta,'-o',d/'indelqual.bam',bam],log)
        run(['samtools','index',d/'indelqual.bam'],log)
        run(['lofreq','call','--call-indels','-f',fasta,'-l',bed,'-o',d/'lofreq.vcf',d/'indelqual.bam'],log)
        # LoFreq 2.1.5 omits contig declarations. Keep its original output and
        # supply the validated reference dictionary before bcftools normalization.
        run(['bcftools','reheader','--fai',str(fasta)+'.fai','-o',d/'lofreq.with_contigs.vcf',d/'lofreq.vcf'],log)
        normalize(d/'lofreq.with_contigs.vcf',d/'lofreq.norm.vcf.gz',c,log)
    if 'sage' in callers:
        sage=['java',f'-Xmx{c["java_heap_gb"]}g','-jar',resource(c,'sage_jar'),
              '-tumor',sm,'-tumor_bam',bam,'-ref_genome',fasta,'-ref_genome_version',c['sage_build'],
              '-hotspots',resource(c,'sage_hotspots'),'-high_confidence_bed',resource(c,'sage_high_confidence'),
              '-driver_gene_panel',resource(c,'sage_driver_panel'),'-ensembl_data_dir',c['sage_ensembl_dir'],
              '-threads',str(c['threads']),'-output_vcf',d/'sage.vcf.gz',
              '-max_read_depth',str(c['evidence']['max_depth']),'-max_read_depth_panel',str(c['evidence']['max_depth'])]
        # Restrict discovery, not just the later export. Java @argfiles avoid the
        # OS per-argument limit for an exome containing hundreds of thousands of intervals.
        regions=';'.join(f'{ch}:{start+1}-{end}' for ch,start,end in bed_intervals(bed))
        sage += ['-specific_regions',regions]
        sage += c['sage_extra_args']
        argfile=d/'sage.arguments.txt'
        atomic(argfile,'\n'.join(json.dumps(str(x)) for x in sage[1:])+'\n')
        run(['java','@'+str(argfile)],log)
        normalize(d/'sage.vcf.gz',d/'sage.norm.vcf.gz',c,log)
    if row.get('germline_vcf'):
        normalize(row['germline_vcf'],d/'germline.norm.vcf.gz',c,log)
    # Normalize carrier identity with exactly the same reference and bcftools operation.
    if row.get('carrier_variant'):
        ch,pos,ref,alt=row['carrier_variant'].split(':')
        atomic(d/'carrier.vcf','##fileformat=VCFv4.2\n'+f'##contig=<ID={ch}>\n'+
               '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'+f'{ch}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.\n')
        normalize(d/'carrier.vcf',d/'carrier.norm.vcf.gz',c,log)
    if not c['keep_derived_bam'] and 'lofreq' in callers:
        (d/'indelqual.bam').unlink();(d/'indelqual.bam.bai').unlink()
    return {name:str(d/(name+'.norm.vcf.gz')) for name in callers+['germline'] if (d/(name+'.norm.vcf.gz')).exists()}

def build_pon(rows,c,out):
    # External technical samples are explicitly disjoint from all biological participants.
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    for assay in sorted({r['assay_id'] for r in rows}):
        subset=[r for r in rows if r['assay_id']==assay]
        if len(subset)<c['pon_min_samples']:raise ValueError('too few technical PoN samples')
        sms=[validate_bam(r,c) for r in subset]
        if len(set(sms))!=len(sms):raise ValueError('duplicate BAM sample names within technical PoN')
        if any(number(r.get('age_at_collection')) is None or float(r['age_at_collection'])>c['pon_max_age'] for r in subset):
            raise ValueError('PoN age eligibility missing')
        ad=out/safe_id(assay);ad.mkdir(exist_ok=True);inputs=[]
        for r in subset:inputs.append(call_sample(r,c,ad/safe_id(r['sample_id']),None,True))
        java=['gatk','--java-options',f'-Xmx{c["java_heap_gb"]}g']
        argv=java+['GenomicsDBImport','-R',resource(c,'fasta'),'-L',resource(c,'common_bed'),
                   '--genomicsdb-workspace-path',ad/'db','--reader-threads',str(c['threads'])]
        for v in inputs:argv+=['-V',v]
        run(argv,ad/'commands.log')
        run(java+['CreateSomaticPanelOfNormals','-R',resource(c,'fasta'),
                  '--germline-resource',resource(c,'gnomad'),'-V','gendb://'+str((ad/'db').resolve()),
                  '-O',ad/'pon.vcf.gz'],ad/'commands.log')
    write_json(out/'pon.complete.json',dict(status='COMPLETE',samples=len(rows),
        resource_lock_sha256=sha(out/'resources.lock.json'),
        outputs={str(p.relative_to(out)):sha(p) for p in out.glob('*/pon.vcf.gz*')},
        technical_people=[r['person_id'] for r in rows],technical_analyses=[r['analysis_id'] for r in rows]))
