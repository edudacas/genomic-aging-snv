"""Native pysam BAM/FASTA tests; absence of pysam is explicitly a skip, not validation."""
import tempfile,unittest
from pathlib import Path
try:import pysam
except ImportError:pysam=None
from fragment_recount import recount_site

@unittest.skipIf(pysam is None,'pysam unavailable: native BAM validation pending')
class NativeTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.fa=self.root/'ref.fa';self.fa.write_text('>chr1\n'+'A'*200+'\n');pysam.faidx(str(self.fa))
    def tearDown(self):self.tmp.cleanup()
    def count(self,reads,cap=100000,pos=21,ref='A',alt='G'):
        path=self.root/'reads.bam'
        header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':200}], 'RG':[{'ID':'r1','SM':'synthetic'},{'ID':'r2','SM':'synthetic'}]}
        with pysam.AlignmentFile(str(path),'wb',header=header) as f:
            for name,rg,flag,bq,mq,base in reads:
                r=pysam.AlignedSegment();r.query_name=name;r.query_sequence='A'*20+base+'A'*29
                r.flag=flag;r.reference_id=0;r.reference_start=0;r.mapping_quality=mq;r.cigar=[(0,50)]
                r.next_reference_id=0;r.next_reference_start=0;r.template_length=50 if flag&64 else -50
                r.query_qualities=[bq]*50;r.set_tag('RG',rg);f.write(r)
        pysam.index(str(path))
        with pysam.AlignmentFile(str(path),'rb') as bam,pysam.FastaFile(str(self.fa)) as fasta:
            return recount_site(bam,fasta,dict(chrom='chr1',pos=pos,ref=ref,alt=alt),raw_max_depth=cap)
    def test_overlap_and_rg(self):
        r=self.count([('a','r1',99,40,60,'G'),('a','r1',147,50,60,'G'),('a','r2',99,40,60,'G')])
        self.assertEqual(r['alt_fragments'],2);self.assertEqual(r['alt_forward'],2)
    def test_mate_discordance(self):
        r=self.count([('a','r1',99,40,60,'G'),('a','r1',147,40,60,'A')])
        self.assertEqual(r['usable_fragment_depth'],0);self.assertEqual(r['diagnostics']['eligible_discordant_fragments'],1)
    def test_filters_and_no_coverage(self):
        r=self.count([('a','r1',99|1024,40,60,'G'),('b','r1',99,29,60,'G'),('c','r1',99,40,29,'G')])
        self.assertEqual(r['usable_fragment_depth'],0)
        self.assertEqual(self.count([],pos=100)['evaluation_status'],'NO_COVERAGE')
    def test_ref_mismatch_and_saturation(self):
        reads=[(str(i),'r1',99,40,60,'G') for i in range(5)]
        self.assertEqual(self.count(reads,ref='C')['evaluation_status'],'REF_MISMATCH')
        r=self.count(reads,cap=2);self.assertEqual(r['evaluation_status'],'DEPTH_SATURATED');self.assertIsNone(r['alt_fragments'])

if __name__=='__main__':unittest.main()
