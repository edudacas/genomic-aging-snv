"""Cross-validate canonical VCF and typed Parquet against the integrated TSVs."""
import csv
import itertools
import json
import math
from pathlib import Path

import pyarrow.parquet as pq
import pysam

from .common import sha, write_json
from .parquet import convert


def _same_float(actual, expected):
    return math.isclose(float(actual), float(expected), rel_tol=6e-6, abs_tol=1e-8)


def validate_views(out, evidence_files):
    out=Path(out);files=[Path(path) for path in evidence_files]
    streams=[path.open() for path in files];readers=[csv.DictReader(stream,delimiter='\t') for stream in streams]
    pairs=records=0
    try:
        with pysam.VariantFile(str(out/'joined.evidence.vcf.gz')) as vcf:
            if 'GT' in vcf.header.formats or any(name.endswith('FT') for name in vcf.header.formats):
                raise ValueError('Canonical v3 VCF contains GT or verbose filter JSON')
            if list(vcf.header.samples)!=[path.parent.name for path in files]:
                raise ValueError('VCF sample order differs from evidence inputs')
            for record in vcf:
                rows=[next(reader) for reader in readers]
                key=f'{record.chrom}:{record.pos}:{record.ref}:{record.alts[0]}'
                if list(record.filter):raise ValueError('Canonical VCF assigns a global filter')
                for path,row in zip(files,rows):
                    if row['variant_id']!=key:raise ValueError('VCF/TSV key mismatch')
                    sample=record.samples[path.parent.name]
                    for field,column in [('DP','dp'),('AF','vaf'),('PON','in_pon'),('CS','caller_support')]:
                        actual=sample[field];expected=row[column]
                        if expected=='':
                            if actual is not None:raise ValueError('VCF missing-value mismatch: '+field)
                        elif field=='AF':
                            if not _same_float(actual,expected):raise ValueError('VCF AF mismatch')
                        elif actual!=int(expected):raise ValueError('VCF integer mismatch: '+field)
                    ad=sample['AD'];expected=(row['ad_ref'],row['ad_alt'])
                    for actual,value in zip(ad,expected):
                        if (value=='' and actual is not None) or (value!='' and actual!=int(value)):
                            raise ValueError('VCF AD mismatch')
                    pairs+=1
                records+=1
        if any(next(reader,None) is not None for reader in readers):raise ValueError('TSV has extra rows')
    finally:
        for stream in streams:stream.close()
    parquet_rows=0
    for path in files:
        target=out/'evidence_parquet'/('sample_id='+path.parent.name)/'part-00000.parquet'
        parquet=pq.ParquetFile(target);spec=parquet.schema_arrow
        with path.open() as stream:
            expected=(convert(row,spec) for row in csv.DictReader(stream,delimiter='\t'))
            actual=(row for batch in parquet.iter_batches(batch_size=10000) for row in batch.to_pylist())
            for wanted,observed in itertools.zip_longest(expected,actual):
                if wanted!=observed:raise ValueError('Parquet/TSV field mismatch: '+path.parent.name)
                parquet_rows+=1
    result=dict(status='PASS',canonical_vcf_records=records,sample_variant_pairs=pairs,
        parquet_rows=parquet_rows,canonical_vcf_sha256=sha(out/'joined.evidence.vcf.gz'),
        canonical_vcf_index_sha256=sha(out/'joined.evidence.vcf.gz.tbi'),
        parquet_manifest_sha256=sha(out/'parquet.manifest.json'),all_parquet_fields_match_tsv=True,
        vcf_numeric_tolerance='relative 6e-6; absolute 1e-8',no_gt=True,no_global_pass=True,
        verbose_filter_json_in_vcf_format=False)
    write_json(out/'artifact_validation.json',result)
    return result
