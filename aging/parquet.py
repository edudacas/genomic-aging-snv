"""Typed, bounded-memory analytical exports; blank values remain null, never zero."""
import csv
import json
import os
import tempfile
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq
from .common import sha, write_json

INTS=set('pos dp ad_ref ad_alt alt_forward alt_reverse alt_f1r2 alt_f2r1 carrier_variant_exact in_pon germline_overlap caller_support discovery_sample_count dp_raw_reads alt_raw_reads mq0_raw_reads duplicate_raw_reads raw_depth_saturated filtered_depth_saturated flag_filtered_raw_reads low_mq_raw_reads low_bq_eligible_reads eligible_reads_before_bq eligible_reads_after_bq other_fragments discordant_raw_fragments discordant_eligible_fragments alt_concordant_both_mates ref_forward ref_reverse ref_f1r2 ref_f2r1 saturation_fallback'.split())
FLOATS=set('vaf alt_mean_bq alt_mean_mq gnomad_af abraom_af'.split())
for prefix in ('ref_','alt_ext_'):
    INTS.add(prefix+'unique_fragment_starts')
    FLOATS.update(prefix+name for name in ('min_bq','median_bq','sum_bq','min_mq','median_mq','sum_mq','median_distance_to_end','near_end_fraction','mean_nm','mean_abs_tlen'))
for caller in ('mutect','lofreq','sage','germline'):
    INTS.update(caller+'_'+x for x in ('dp','ad_alt','ad_ref'))
QC_INTS=set('target_bp mapped_primary_nonduplicate_reads umi_sampled_reads depth_zero_bp bam_size bam_mtime_ns'.split())
QC_FLOATS=set('mean_depth on_target_fraction duplicate_fraction read_length insert_size'.split())


def dtype(name):
    if name in INTS or name.removeprefix('qc_') in QC_INTS and name.startswith('qc_'):return pa.int64()
    if name in FLOATS or name.removeprefix('qc_') in QC_FLOATS and name.startswith('qc_'):return pa.float64()
    return pa.string()


def schema(fields):
    return pa.schema([pa.field(k,dtype(k),nullable=True) for k in fields],
                     metadata={b'coordinate_system':b'1-based VCF positions',
                               b'missing_values':b'null is not zero; inspect evidence_status and caller filters',
                               b'schema_version':b'1'})


def convert(row, spec):
    result={}
    for field in spec:
        value=row.get(field.name)
        if value is None or value=='':result[field.name]=None
        elif pa.types.is_integer(field.type):result[field.name]=int(value)
        elif pa.types.is_floating(field.type):result[field.name]=float(value)
        elif isinstance(value,(dict,list)):result[field.name]=json.dumps(value,sort_keys=True)
        else:result[field.name]=str(value)
    return result


def write_rows(path, rows, spec, batch_size=10000):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix='.'+path.name,dir=path.parent);os.close(fd)
    count=0
    try:
        with pq.ParquetWriter(tmp,spec,compression='zstd') as writer:
            batch=[]
            for row in rows:
                batch.append(convert(row,spec))
                if len(batch)>=batch_size:
                    writer.write_table(pa.Table.from_pylist(batch,schema=spec));count+=len(batch);batch=[]
            if batch:writer.write_table(pa.Table.from_pylist(batch,schema=spec));count+=len(batch)
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp):os.unlink(tmp)
    return {'rows':count,'sha256':sha(path)}


def export_run(out, samples, candidates, evidence_files):
    out=Path(out); outputs={}
    fields=['sample_id','variant_id','chrom','pos','ref','alt','evidence_status','dp','ad_ref','ad_alt','vaf']
    for path in evidence_files:
        with open(path) as stream:
            for name in csv.DictReader(stream,delimiter='\t').fieldnames or []:
                if name not in fields:fields.append(name)
    # sample_id is the Hive partition key; omit it from file columns to avoid
    # duplicate string/dictionary schemas when reading the whole dataset.
    evidence_schema=schema([field for field in fields if field!='sample_id'])
    for path in evidence_files:
        path=Path(path);name='evidence_parquet/sample_id='+path.parent.name+'/part-00000.parquet'
        with path.open() as stream:
            outputs[name]=write_rows(out/name,csv.DictReader(stream,delimiter='\t'),evidence_schema)
    candidate_fields=['variant_id','chrom','pos','ref','alt','discovered_by','discovery_sample_count']
    def candidate_rows():
        for key, value in sorted(candidates.items()):
            if value.get('variant_id', key) != key:
                raise ValueError('Candidate key conflicts with variant_id: ' + key)
            yield dict(value, variant_id=key)
    outputs['candidates.parquet']=write_rows(out/'candidates.parquet',
        candidate_rows(),schema(candidate_fields))
    enriched=[];sample_fields=[]
    for sample in samples:
        row=dict(sample)
        qc=out/'qc'/row['sample_id']/'qc.json'
        if qc.is_file():
            row.update({'qc_'+k:v for k,v in json.loads(qc.read_text()).items()})
        enriched.append(row)
        for name in row:
            if name not in sample_fields:sample_fields.append(name)
    outputs['samples.parquet']=write_rows(out/'samples.parquet',enriched,schema(sample_fields))
    write_json(out/'parquet.manifest.json',{'schema_version':2,'canonical_joined_vcf':'joined.evidence.vcf.gz','partitioning':'hive: sample_id','outputs':outputs,
               'evidence_rows':sum(v['rows'] for k,v in outputs.items() if k.startswith('evidence_parquet/'))})
    return {'parquet_manifest':'parquet.manifest.json'}
