import random,tempfile,unittest
from pathlib import Path
import pysam
from aging.evidence import pileup_site
from aging.window_recount import windows,recount_window,ordered_sites

class WindowTests(unittest.TestCase):
    def test_fragment_recount_matches_baseline(self):
        rng=random.Random(51)
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);fa=p/'ref.fa';fa.write_text('>chr1\n'+'A'*2000+'\n');pysam.faidx(str(fa))
            bam=p/'reads.bam';header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr1','LN':2000}]};reads=[]
            for i in range(250):
                r=pysam.AlignedSegment();r.query_name='fragment'+str(i//2);r.reference_id=0;r.reference_start=rng.randrange(50,500)
                r.flag=1| (64 if i%2==0 else 128) | (16 if i%3==0 else 0)
                if i%19==0:r.flag|=1024
                if i%29==0:r.flag|=256
                r.mapping_quality=10 if i%17==0 else 60
                r.cigartuples=[(0,100)];seq=list('A'*100)
                for pos in [100,150,200,250,300,350,400,450]:
                    q=pos-1-r.reference_start
                    if 0<=q<len(seq) and i%3==0:seq[q]='C'
                r.query_sequence=''.join(seq);r.query_qualities=[10 if i%13==0 else 35]*100
                r.next_reference_id=0;r.next_reference_start=r.reference_start+10;reads.append(r)
            # Anchored insertion/deletion, overlapping mates and conflicting alleles at same fragment.
            for i,cigar,seq in [(300,[(0,20),(1,2),(0,80)],'A'*20+'TT'+'A'*80),(301,[(0,20),(2,2),(0,80)],'A'*100)]:
                r=pysam.AlignedSegment();r.query_name='indel'+str(i);r.reference_id=0;r.reference_start=180;r.mapping_quality=60;r.cigartuples=cigar;r.query_sequence=seq;r.query_qualities=[35]*len(seq);reads.append(r)
            with pysam.AlignmentFile(str(bam),'wb',header=header) as w:
                for r in sorted(reads,key=lambda r:r.reference_start):w.write(r)
            pysam.index(str(bam))
            candidates={}
            for pos in [100,150,200,250,300,350,400,450,900]:
                candidates[f'chr1:{pos}:A:C']=dict(chrom='chr1',pos=pos,ref='A',alt='C')
            for pos,ref,alt in [(200,'A','ATT'),(200,'AAA','A'),(220,'AA','CC')]:
                candidates[f'chr1:{pos}:{ref}:{alt}']=dict(chrom='chr1',pos=pos,ref=ref,alt=alt)
            with pysam.FastaFile(str(fa)) as f,pysam.AlignmentFile(str(bam),'rb') as b:
                sites=ordered_sites(candidates,f)
                self.assertEqual([r['pos'] for r in sites],sorted(r['pos'] for r in sites))
                for cap in [10000,3]:
                    settings={'min_baseq':20,'min_mapq':20,'max_depth':cap};fallback=False
                    for block in windows(sites):
                        actual,sat=recount_window(b,f,block,settings);fallback|=sat
                        for row in block:
                            expected=pileup_site(b,f,row,20,20,cap)
                            if expected['evidence_status']=='UNSUPPORTED_COMPLEX_ALLELE':expected.update(dp='',ad_ref='',ad_alt='',vaf='')
                            self.assertEqual(actual[row['variant_id']],expected,(cap,row['variant_id']))
                    self.assertEqual(fallback,cap==3)

if __name__=='__main__':unittest.main()
