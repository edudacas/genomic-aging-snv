import tempfile,unittest,json
from pathlib import Path
from aging.caller_records import load_records

class RecordsTests(unittest.TestCase):
    def test_collision_keeps_pass_and_filtered_no_ad_summing(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'sage.vcf';p.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n##FILTER=<ID=dedupIndel,Description="duplicate">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Alleles">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\nchr1\t10\t.\tA\tC\t.\tPASS\t.\tDP:AD\t20:15,5\nchr1\t10\t.\tA\tC\t.\tdedupIndel\t.\tDP:AD\t18:14,4\n')
            x=load_records(p,'sage')['chr1:10:A:C']
            self.assertEqual(x['record_count'],2);self.assertEqual([r['filter'] for r in x['records']],['PASS','dedupIndel'])
            self.assertEqual(x['filter'],'MULTIPLE_RECORDS');self.assertEqual(x['dp'],'');self.assertEqual(x['ad_alt'],'')
            self.assertEqual(json.loads(x['records_json'])[0]['samples']['S']['AD'],[15,5])
    def test_lofreq_sites_only_keeps_dp_and_af_without_inventing_ad(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'l.vcf';p.write_text('##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100>\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n##INFO=<ID=AF,Number=A,Type=Float,Description="AF">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t10\t.\tA\tC\t.\tPASS\tDP=123;AF=0.04\n')
            x=load_records(p,'lofreq')['chr1:10:A:C'];self.assertEqual(x['dp'],123);self.assertAlmostEqual(x['af'],.04,places=6);self.assertEqual(x['ad_alt'],'')
if __name__=='__main__':unittest.main()
