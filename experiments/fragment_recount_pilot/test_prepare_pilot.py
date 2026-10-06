import tempfile,unittest,gzip,json,subprocess,sys
from pathlib import Path
from prepare_pilot import bed_index,inside,parse_snv,validate,sha,sfp

class PilotTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.bed=self.root/'common.bed';self.bed.write_text('chr1\t0\t10\nchr1\t5\t20\n')
        self.index=bed_index(self.bed)
    def tearDown(self):self.tmp.cleanup()
    def test_bed_boundaries(self):
        self.assertTrue(inside(self.index,'chr1',1));self.assertTrue(inside(self.index,'chr1',20))
        self.assertFalse(inside(self.index,'chr1',21));self.assertFalse(inside(self.index,'1',1))
    def test_filtered_high_vaf_not_discarded(self):
        rows,reason=parse_snv('chr1\t2\t.\tA\tG\t1\tgermline\tAF=0.99\n',self.index)
        self.assertEqual(rows,[('chr1',2,'A','G')]);self.assertIsNone(reason)
    def test_multiallelic_and_indel_preserved_in_skip_audit(self):
        rows,reason=parse_snv('chr1\t2\t.\tA\tG,AT\t1\tPASS\t.\n',self.index)
        self.assertEqual(rows,[('chr1',2,'A','G')]);self.assertEqual(reason,'non_snv_alleles')
    def test_outside_territory(self):
        self.assertEqual(parse_snv('chr1\t21\t.\tA\tG\t.\t.\t.\n',self.index),([], 'outside_common_bed'))
    def test_missing_receipts_fail_closed(self):
        (self.root/'samples.private.tsv').write_text('sample_id\nS1\nS2\nS3\n')
        with self.assertRaises(FileNotFoundError):validate(self.root,['S1','S2','S3'],3)
    def test_complete_synthetic_preparation_and_corruption(self):
        stages=self.root/'stages';stages.mkdir();ids=['S1','S2','S3'];fp='test';previous=None;hashes={}
        (stages/'samples.private.tsv').write_text('sample_id\nS1\nS2\nS3\n')
        for stage in ('qc','sage','mutect','lofreq'):
            fingerprint=sfp(fp,stage,previous);samples={}
            for sid in ids:
                directory=stages/stage/sid;directory.mkdir(parents=True);outputs={}
                if stage!='qc':
                    path=directory/(stage+'.norm.vcf.gz')
                    with gzip.open(path,'wt') as f:
                        f.write('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
                        f.write('chr1\t2\t.\tA\tG\t.\tgermline\tAF=0.99\n')
                    outputs[path.name]=sha(path)
                marker=directory/'COMPLETE.json';marker.write_text(json.dumps(dict(status='COMPLETE',fingerprint=fingerprint,outputs=outputs)))
                samples[sid]={'receipt_sha256':sha(marker)}
            receipt=stages/(stage.upper()+'_COMPLETE.json');receipt.write_text(json.dumps(dict(status=stage.upper()+'_COMPLETE',
                fingerprint=fingerprint,config_fingerprint=fp,samples=samples,failures={})))
            previous=sha(receipt);hashes[stage]=previous
        (stages/'CALLERS_COMPLETE.json').write_text(json.dumps(dict(status='CALLERS_COMPLETE',selected_samples=3,
            config_fingerprint=fp,union_performed=False,recount_performed=False,stages=hashes)))
        command=[sys.executable,str(Path(__file__).parent/'prepare_pilot.py'),'--stages',str(stages),'--bed',str(self.bed),
            '--expected','3','--out',str(self.root/'prepared')]
        for sid in ids:command+=['--pilot-sample',sid]
        subprocess.run(command,check=True,capture_output=True)
        receipt=json.loads((self.root/'prepared/PREPARATION_RECEIPT.PRIVATE.json').read_text())
        self.assertEqual(receipt['candidate_snv_alleles'],1);self.assertEqual(receipt['planned_sample_allele_pairs'],3)
        self.assertFalse(receipt['recount_performed']);self.assertFalse(receipt['scope_gate_frozen'])
        self.assertIn('RECOUNT_PENDING',(self.root/'prepared/candidate_union.pilot.vcf').read_text())
        path=stages/'sage/S1/sage.norm.vcf.gz';path.write_bytes(b'corrupted')
        with self.assertRaisesRegex(ValueError,'Caller VCF hash mismatch'):validate(stages,ids,3)

if __name__=='__main__':unittest.main()
