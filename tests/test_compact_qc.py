import csv
import gzip
import io
import json
import tempfile
import unittest
from pathlib import Path
from aging.compact_qc import summarize, discard_verified_depth, split_territory
from aging.common import write_json


class CompactCoverageTest(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory(); self.addCleanup(self.tmp.cleanup)
        self.p=Path(self.tmp.name); self.bed=self.p/'common.bed'; self.mask=self.p/'mask.json'
        self.bed.write_text('chr1\t0\t4\nchrX\t0\t6\nchrY\t0\t8\n')
        write_json(self.mask,dict(build='GRCh38',coordinates='0-based-half-open',version='synthetic',masks=[
            dict(chrom='chrX',start=2,end=4,compartment='X_PAR'),
            dict(chrom='chrY',start=2,end=4,compartment='Y_PAR'),
            dict(chrom='chrY',start=4,end=6,compartment='Y_XTR')]))

    def test_borders_zeros_and_thresholds(self):
        # omitted positions are real zeros, including entire intervals
        depth='chr1\t1\t20\nchr1\t4\t10\nchrX\t2\t10\nchrX\t3\t20\nchrY\t4\t50\nchrY\t5\t80\nchrY\t8\t30\n'
        hist,r=summarize(io.StringIO(depth),self.bed,self.mask,self.p,{'mq':30},[10,20,50])
        self.assertEqual(sum(hist.values()),18)
        self.assertEqual(hist[0],11)
        self.assertAlmostEqual(r['x_auto_ratio'],(10/4)/(30/4))
        self.assertAlmostEqual(r['y_auto_ratio'],(30/4)/(30/4))
        with gzip.open(self.p/'coverage.intervals.tsv.gz','rt') as f: rows=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(len(rows),8)
        self.assertEqual(sum(int(x['depth_sum']) for x in rows),220)
        self.assertEqual(sum(int(x['bases_ge_20']) for x in rows),5)
        self.assertTrue(any(int(x['covered_bases'])==0 for x in rows))

    def test_synthetic_reconstruction_of_pilot_ratios(self):
        # Rounded means from the three-pilot table, not a claim of BAM validation.
        self.bed.write_text('chr1\t0\t100\nchrX\t0\t100\nchrY\t0\t100\n')
        write_json(self.mask,dict(build='GRCh38',coordinates='0-based-half-open',version='synthetic',masks=[]))
        for a,x,y in [(99.91,103.07,0.0),(133.1,69.14,40.55),(153.78,173.62,0.01)]:
            lines=[]
            for chrom,mean in [('chr1',a),('chrX',x),('chrY',y)]:
                total=round(mean*100); q,rem=divmod(total,100)
                lines += [f'{chrom}\t{i+1}\t{q+(i<rem)}\n' for i in range(100) if q+(i<rem)>0]
            _,r=summarize(io.StringIO(''.join(lines)),self.bed,self.mask,self.p,{},[20,50,100])
            self.assertAlmostEqual(r['x_auto_ratio'],x/a,places=12)
            self.assertAlmostEqual(r['y_auto_ratio'],y/a,places=12)

    def test_invalid_input_and_overlaps_fail(self):
        for text in ['chr1\t1\t20\nchr1\t1\t20\n','chr1\t9\t20\n','chr1\t1\t-1\n']:
            with self.assertRaises(ValueError): summarize(io.StringIO(text),self.bed,self.mask,self.p,{},[20])
        self.bed.write_text('chr1\t0\t4\nchr1\t3\t5\n')
        with self.assertRaisesRegex(ValueError,'nonoverlapping'): split_territory(self.bed,self.mask)

    def test_cleanup_only_after_verified_outputs_and_qc(self):
        depth=self.p/'depth.tsv'; depth.write_text('chr1\t1\t20\n')
        summarize(depth,self.bed,self.mask,self.p,{},[20])
        with self.assertRaisesRegex(ValueError,'QC not complete'): discard_verified_depth(depth,self.p)
        self.assertTrue(depth.exists())
        write_json(self.p/'qc.json',{'status':'TECHNICAL_QC_COMPLETE'})
        original=(self.p/'coverage.compartments.tsv').read_bytes()
        (self.p/'coverage.compartments.tsv').write_text('corrupt')
        with self.assertRaisesRegex(ValueError,'changed'): discard_verified_depth(depth,self.p)
        self.assertTrue(depth.exists())
        (self.p/'coverage.compartments.tsv').write_bytes(original)
        discard_verified_depth(depth,self.p)
        self.assertFalse(depth.exists())

    def test_runtime_generates_compact_outputs_without_changing_caller_function(self):
        import pysam
        from unittest.mock import patch
        from aging.runtime import qc
        from aging.common import sha
        bam=self.p/'input-source.bam'; fasta=self.p/'ref.fa'
        fasta.write_text('>chr1\n'+'A'*100+'\n'); pysam.faidx(str(fasta))
        with pysam.AlignmentFile(str(bam),'wb',header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':100}],'RG':[{'ID':'r','SM':'S'}]}) as f: pass
        pysam.index(str(bam))
        c=dict(threads=2,resources={'common_bed':{'path':str(self.bed)}},evidence={'min_baseq':20,'min_mapq':30},
               compact_qc={'mask_path':str(self.mask),'mask_sha256':sha(self.mask),'thresholds':[20,50]})
        def run(argv,log,stdout=None):
            if argv[1]=='stats': Path(stdout).write_text('SN\traw total sequences:\t100\nSN\treads duplicated:\t0\n')
            elif argv[1]=='depth': Path(stdout).write_text('chr1\t1\t20\n')
            elif argv[1]=='view': Path(stdout).write_text('100\n')
        row={'sample_id':'S','bam_path':str(bam),'bai_path':str(bam)+'.bai','assay_id':'A'}
        with patch('aging.runtime.validate_bam',return_value='S'),patch('aging.runtime.run',side_effect=run):
            q=qc(row,c,self.p/'qc')
        self.assertEqual(q['target_bp'],18)
        self.assertEqual(q['depth_zero_bp'],17)
        self.assertFalse((self.p/'qc/depth.tsv').exists())
        self.assertTrue((self.p/'qc/depth_histogram.tsv').exists())


if __name__=='__main__': unittest.main()
