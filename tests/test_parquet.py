import tempfile
import unittest
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq
from aging.parquet import schema,write_rows,export_run


class ParquetTests(unittest.TestCase):
    def test_embedded_candidate_id_and_conflict(self):
        with tempfile.TemporaryDirectory() as tmp:
            key = 'chr1:1:A:T'
            value = dict(variant_id=key, chrom='chr1', pos=1, ref='A', alt='T')
            export_run(tmp, [dict(sample_id='S')], {key:value}, [])
            self.assertEqual(pq.read_table(Path(tmp)/'candidates.parquet')['variant_id'].to_pylist(), [key])
            value['variant_id'] = 'conflict'
            with self.assertRaisesRegex(ValueError, 'conflicts'):
                export_run(tmp, [dict(sample_id='S')], {key:value}, [])

    def test_null_zero_and_identifiers_roundtrip(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'part.parquet'
            fields=['sample_id','dp','ad_alt','vaf','mutect_filter','qc_mean_depth']
            rows=[dict(sample_id='001',dp='0',ad_alt='0',vaf='',mutect_filter='NOT_CALLED'),
                  dict(sample_id='002',dp='100',ad_alt='0',vaf='0',mutect_filter='artifact',qc_mean_depth=30.5)]
            result=write_rows(path,rows,schema(fields),batch_size=1)
            table=pq.read_table(path)
            self.assertEqual(result['rows'],2)
            self.assertEqual(table['sample_id'].to_pylist(),['001','002'])
            self.assertEqual(table['vaf'].to_pylist(),[None,0.0])
            self.assertEqual(table['ad_alt'].to_pylist(),[0,0])
            self.assertEqual(table['mutect_filter'].to_pylist(),['NOT_CALLED','artifact'])
            self.assertEqual(table['qc_mean_depth'].to_pylist(),[None,30.5])
            self.assertEqual(table.schema.field('dp').type,pa.int64())

    def test_empty_file_retains_numeric_schema(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'empty.parquet'
            write_rows(p,[],schema(['sample_id','pos','vaf']))
            t=pq.read_table(p)
            self.assertEqual(t.num_rows,0)
            self.assertEqual(t.schema.field('vaf').type,pa.float64())
