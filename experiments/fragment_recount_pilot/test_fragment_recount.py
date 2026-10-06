import itertools,unittest
from fragment_recount import reduce_observations, recount_site

def obs(name='fragment',base='G',rg='rg',**changes):
    d=dict(rg=rg,qname=name,base=base,flag=0,bq=40,mq=60,paired=True,read_number=1,
           reverse=False,mate_reverse=True,mate_unmapped=False,distance_to_end=20)
    d.update(changes);return d

class FragmentTests(unittest.TestCase):
    def reduce(self,rows,**kw):return reduce_observations(rows,'A','G',**kw)
    def test_concordant_overlap_once(self):
        r=self.reduce([obs(),obs(read_number=2,reverse=True,mate_reverse=False,bq=50)])
        self.assertEqual(r['alt_fragments'],1);self.assertEqual(r['alt_forward'],1);self.assertEqual(r['alt_f1r2'],1)
    def test_discordant_excluded(self):
        r=self.reduce([obs(),obs(base='A',read_number=2)])
        self.assertEqual(r['usable_fragment_depth'],0);self.assertEqual(r['diagnostics']['eligible_discordant_fragments'],1)
    def test_different_rg(self):self.assertEqual(self.reduce([obs(rg='a'),obs(rg='b')])['alt_fragments'],2)
    def test_excluded_flags(self):
        for flag in (0x4,0x100,0x200,0x400,0x800):self.assertEqual(self.reduce([obs(flag=flag)])['alt_fragments'],0)
    def test_low_quality(self):
        for kw in ({'bq':29},{'mq':29},{'mq':255},{'bq':None}):self.assertEqual(self.reduce([obs(**kw)])['alt_fragments'],0)
    def test_low_quality_disagreement_does_not_override(self):
        r=self.reduce([obs(),obs(base='A',bq=10,read_number=2)])
        self.assertEqual(r['alt_fragments'],1);self.assertEqual(r['diagnostics']['raw_discordant_fragments'],1)
    def test_order_and_quality_ties(self):
        reads=[obs(),obs(read_number=2,reverse=True,mate_reverse=False)]
        results=[self.reduce(p) for p in itertools.permutations(reads)]
        self.assertEqual(results[0],results[1])
    def test_single_strand_retained_not_classified(self):
        r=self.reduce([obs(name=str(i)) for i in range(5)])
        self.assertEqual(r['alt_forward'],5);self.assertEqual(r['alt_reverse'],0);self.assertNotIn('PASS',r)
    def test_orientation_and_orphan(self):
        r=self.reduce([obs(),obs(name='b',reverse=True),obs(name='c',paired=False),
                       obs(name='d',read_number=2,mate_unmapped=True)])
        self.assertEqual(r['alt_f1r2'],1);self.assertEqual(r['alt_f2r1'],2);self.assertEqual(r['alt_unpaired'],1)
        self.assertEqual(r['alt_strand_unknown'],1)
    def test_no_coverage(self):
        r=self.reduce([]);self.assertEqual(r['evaluation_status'],'NO_COVERAGE');self.assertIsNone(r['vaf_alt_over_usable'])
    def test_saturation_is_not_zero(self):
        r=self.reduce([obs()],saturated=True);self.assertIsNone(r['alt_fragments']);self.assertEqual(r['evaluation_status'],'DEPTH_SATURATED')
    def test_other_changes_vaf_denominator(self):
        r=self.reduce([obs(),obs(name='ref',base='A'),obs(name='other',base='C')])
        self.assertEqual(r['vaf_alt_over_ref_alt'],.5);self.assertEqual(r['vaf_alt_over_usable'],1/3)
    def test_boundary_not_classified(self):
        for n,alt in ((100,5),(100,35),(500,5),(50,3)):
            r=self.reduce([obs(name=str(i),base='G' if i<alt else 'A') for i in range(n)])
            self.assertEqual(r['vaf_alt_over_usable'],alt/n);self.assertNotIn('countable',r)
    def test_missing_rg_audited(self):self.assertEqual(self.reduce([obs(rg='')])['diagnostics']['reads_missing_rg'],1)
    def test_ref_mismatch_does_not_open_pileup(self):
        class Fasta:
            def fetch(self,*a):return 'C'
        r=recount_site(object(),Fasta(),dict(chrom='chr1',pos=1,ref='A',alt='G'))
        self.assertEqual(r['evaluation_status'],'REF_MISMATCH');self.assertIsNone(r['alt_fragments'])
    def test_unsupported(self):
        r=recount_site(object(),object(),dict(chrom='chr1',pos=1,ref='A',alt='AT'))
        self.assertEqual(r['evaluation_status'],'UNSUPPORTED_ALLELE')

if __name__=='__main__':unittest.main()
