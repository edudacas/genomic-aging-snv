"""Allele-keyed multi-sample evidence VCF; no invented GT, PASS or somatic classification."""
import csv,json,os,sqlite3,tempfile
from pathlib import Path
import pysam
from .common import sha,write_json

def number(value,typ=int):return None if value in ('',None,'.') else typ(value)

def matrix_rows(files):
    """Transpose one sample at a time on disk; bounded handles even above 1,024 samples."""
    with tempfile.TemporaryDirectory(prefix='.joined-matrix-',dir=files[0].parent.parent.parent) as tmp:
        db=sqlite3.connect(str(Path(tmp)/'matrix.sqlite'))
        try:
            db.execute('PRAGMA temp_store=FILE');db.execute('PRAGMA cache_size=-16384')
            db.execute('CREATE TABLE cells (rownum INTEGER, samplenum INTEGER, data TEXT, PRIMARY KEY(rownum,samplenum)) WITHOUT ROWID')
            count=None
            for sid,file in enumerate(files):
                n=0;batch=[]
                with file.open() as stream:
                    for n,row in enumerate(csv.DictReader(stream,delimiter='\t'),1):
                        batch.append((n,sid,json.dumps(row,separators=(',',':'))))
                        if len(batch)>=1000:db.executemany('INSERT INTO cells VALUES (?,?,?)',batch);batch=[]
                if batch:db.executemany('INSERT INTO cells VALUES (?,?,?)',batch)
                db.commit()
                if count is not None and n!=count:raise ValueError('Incomplete evidence matrix')
                count=n
            group=[]
            for rownum,sid,data in db.execute('SELECT rownum,samplenum,data FROM cells ORDER BY rownum,samplenum'):
                if sid!=len(group):raise ValueError('Missing matrix sample')
                group.append(json.loads(data))
                if len(group)==len(files):yield group;group=[]
            if group:raise ValueError('Incomplete final evidence row')
        finally:db.close()


def export_joined(evidence_files,out,contigs):
    out=Path(out);h=pysam.VariantHeader()
    h.add_meta('source','Genomic_aging_pipeline_evidence_v2')
    h.add_meta('evidenceSemantics','DP and AD are fragment-aware recounts; no GT or biological classification; FILTER dot is not PASS')
    h.add_meta('callerProvenance','Original per-record INFO FORMAT and filters retained in hashed per-sample caller_records.tsv.gz sidecars; no last-record overwrite')
    for chrom,length in contigs:h.contigs.add(chrom,length=length)
    for name,n,typ,desc in [('CALLERS','.', 'String','Union of discovery callers across samples'),('NSDISC',1,'Integer','Number of discovery samples'),('GAF',1,'Float','Exact allele gnomAD AF; missing is not zero'),('AAF',1,'Float','Exact allele ABraOM AF; missing is not zero')]:
        h.add_meta('INFO',items=[('ID',name),('Number',n),('Type',typ),('Description',desc)])
    formats=[('DP',1,'Integer','Recount fragment depth; null for unsupported alleles'),('AD','R','Integer','Recount reference and alternate fragment support; other alleles may make sum less than DP'),('AF',1,'Float','Recount alternate fragment fraction'),('ES',1,'String','Evidence status'),('PON',1,'Integer','Exact allele presence in assay PoN'),('CS',1,'Integer','Number of callers reporting allele in this sample')]
    formats += [('SB',2,'Integer','Alternate fragment forward and reverse counts'),('OR',2,'Integer','Alternate fragment F1R2 and F2R1 counts'),('BQ',1,'Float','Mean alternate base quality'),('MQ',1,'Float','Mean alternate mapping quality')]
    for caller,prefix in [('mutect','M'),('lofreq','L'),('sage','S')]:
        formats += [(prefix+'RC',1,'Integer',caller+' record count; full filters in hashed sidecar'),(prefix+'CF',1,'Integer',caller+' conflicting record summaries'),(prefix+'DP',1,'Integer',caller+' depth only if unambiguous'),(prefix+'AD','R','Integer',caller+' reference and alternate support only if available and unambiguous')]
    for name,n,typ,desc in formats:h.add_meta('FORMAT',items=[('ID',name),('Number',n),('Type',typ),('Description',desc)])
    files=[Path(f) for f in evidence_files];sids=[f.parent.name for f in files]
    if len(set(sids))!=len(sids):raise ValueError('Duplicate sample IDs')
    for sid in sids:h.add_sample(sid)
    matrix=matrix_rows(files)
    tmp=out.with_name('.'+out.stem+'.tmp.vcf.gz');count=0;previous=None;rank={c:i for i,(c,n) in enumerate(contigs)}
    try:
        with pysam.VariantFile(str(tmp),'wz',header=h) as writer:
            for rows in matrix:
                first=rows[0];key=first['variant_id']
                if any(r['variant_id']!=key or r['sample_id']!=sid for r,sid in zip(rows,sids)):raise ValueError('Evidence keys/sample order mismatch')
                order=(rank[first['chrom']],int(first['pos']),first['ref'],first['alt'])
                if previous is not None and order<=previous:raise ValueError('Unsorted/duplicate allele in joined VCF')
                previous=order
                record=writer.new_record(contig=first['chrom'],start=int(first['pos'])-1,alleles=(first['ref'],first['alt']))
                callers=first.get('discovered_by','')
                if callers:record.info['CALLERS']=tuple(callers.split(';'))
                ns=number(first.get('discovery_sample_count'))
                if ns is not None:record.info['NSDISC']=ns
                for field,col in [('GAF','gnomad_af'),('AAF','abraom_af')]:
                    values=[r.get(col,'') for r in rows]
                    if len(set(values))!=1:raise ValueError('Population annotations differ across samples')
                    value=number(values[0],float)
                    if value is not None:record.info[field]=value
                for sid,r in zip(sids,rows):
                    sample=record.samples[sid];sample['DP']=number(r['dp']);sample['AD']=(number(r['ad_ref']),number(r['ad_alt']));sample['AF']=number(r['vaf'],float)
                    sample['ES']=r['evidence_status'];sample['PON']=number(r['in_pon']);sample['CS']=number(r['caller_support'])
                    sample['SB']=(number(r.get('alt_forward')),number(r.get('alt_reverse')))
                    sample['OR']=(number(r.get('alt_f1r2')),number(r.get('alt_f2r1')))
                    sample['BQ']=number(r.get('alt_mean_bq'),float);sample['MQ']=number(r.get('alt_mean_mq'),float)
                    for caller,prefix in [('mutect','M'),('lofreq','L'),('sage','S')]:
                        sample[prefix+'RC']=number(r[caller+'_record_count']);sample[prefix+'CF']=number(r[caller+'_record_conflict'])
                        sample[prefix+'DP']=number(r[caller+'_dp']);sample[prefix+'AD']=(number(r[caller+'_ad_ref']),number(r[caller+'_ad_alt']))
                writer.write(record);count+=1
        pysam.tabix_index(str(tmp),preset='vcf',force=True)
        with pysam.VariantFile(str(tmp)) as check:
            if list(check.header.samples)!=sids:raise ValueError('Joined sample header mismatch')
            observed=sum(1 for _ in check)
        if observed!=count:raise ValueError('Joined VCF row mismatch')
        with pysam.TabixFile(str(tmp)) as tb:
            if sum(sum(1 for _ in tb.fetch(c)) for c in tb.contigs)!=count:raise ValueError('Joined index mismatch')
        os.replace(tmp,out);os.replace(str(tmp)+'.tbi',str(out)+'.tbi')
    finally:
        matrix.close()
    sidecars={sid:{'file':str((file.parent/'caller_records.tsv.gz').relative_to(out.parent)),
                   'sha256':sha(file.parent/'caller_records.tsv.gz')} for sid,file in zip(sids,files)}
    result=dict(file=out.name,records=count,samples=len(sids),sha256=sha(out),index_sha256=sha(Path(str(out)+'.tbi')),no_genotypes_assigned=True,caller_record_sidecars=sidecars)
    write_json(out.with_name('joined.manifest.json'),result)
    return result
