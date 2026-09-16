import tempfile
import unittest
from pathlib import Path

import pysam

from aging.evidence import pileup_site
from aging.single_pass import recount_single_pass
from aging.optimized import recount_sample_single_pass
from aging.common import read_tsv


class SinglePassTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.fasta = self.root / "reference.fa"
        self.fasta.write_text(">chr1\n" + "A" * 1000 + "\n>chr2\n" + "A" * 300 + "\n")
        pysam.faidx(str(self.fasta))
        self.bam = self.root / "reads.bam"
        header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                  "SQ": [{"SN": "chr1", "LN": 1000}, {"SN": "chr2", "LN": 300}]}
        reads = []
        for name, start, allele, flag, mapq, bq in [
            ("alt", 80, "C", 0, 60, 35),
            ("ref", 80, "A", 0, 60, 35),
            ("lowmq", 80, "C", 0, 0, 35),
            ("duplicate", 80, "C", 1024, 60, 35),
            ("lowbq", 80, "C", 0, 60, 10),
        ]:
            read = pysam.AlignedSegment()
            read.query_name = name
            read.flag = flag
            read.reference_id = 0
            read.reference_start = start
            read.mapping_quality = mapq
            read.cigartuples = [(0, 100)]
            sequence = list("A" * 100)
            sequence[100 - start] = allele
            read.query_sequence = "".join(sequence)
            read.query_qualities = [bq] * 100
            read.set_tag("NM", int(allele != "A"))
            reads.append(read)
        with pysam.AlignmentFile(str(self.bam), "wb", header=header) as stream:
            for read in reads:
                stream.write(read)
            read = pysam.AlignedSegment()
            read.query_name = "second_contig_alt"
            read.flag = 0
            read.reference_id = 1
            read.reference_start = 0
            read.mapping_quality = 60
            read.cigartuples = [(0, 50)]
            sequence = list("A" * 50)
            sequence[9] = "C"
            read.query_sequence = "".join(sequence)
            read.query_qualities = [35] * 50
            stream.write(read)
        pysam.index(str(self.bam))

    def test_single_pass_matches_pointwise_filtered_evidence(self):
        sites = [
            dict(variant_id="chr1:101:A:C", chrom="chr1", pos=101, ref="A", alt="C"),
            dict(variant_id="chr1:250:A:C", chrom="chr1", pos=250, ref="A", alt="C"),
            dict(variant_id="chr2:10:A:C", chrom="chr2", pos=10, ref="A", alt="C"),
        ]
        settings = {"min_baseq": 20, "min_mapq": 30, "max_depth": 1000,
                    "raw_max_depth": 10000}
        with pysam.AlignmentFile(str(self.bam), "rb") as bam, pysam.FastaFile(str(self.fasta)) as fasta:
            actual, diagnostics = recount_single_pass(bam, fasta, sites, settings)
        with pysam.AlignmentFile(str(self.bam), "rb") as bam, pysam.FastaFile(str(self.fasta)) as fasta:
            for site in sites:
                expected = pileup_site(bam, fasta, site, 20, 30, 1000)
                for field in ("evidence_status", "dp", "ad_ref", "ad_alt", "vaf",
                              "alt_forward", "alt_reverse", "alt_f1r2", "alt_f2r1",
                              "alt_mean_bq", "alt_mean_mq"):
                    self.assertEqual(actual[site["variant_id"]][field], expected[field], (site, field))
        observed = actual["chr1:101:A:C"]
        self.assertEqual((observed["dp_raw_reads"], observed["alt_raw_reads"]), (5, 4))
        self.assertEqual(observed["mq0_raw_reads"], 1)
        self.assertEqual(observed["duplicate_raw_reads"], 1)
        self.assertEqual(observed["flag_filtered_raw_reads"], 1)
        self.assertEqual(observed["low_mq_raw_reads"], 1)
        self.assertEqual((observed["eligible_reads_before_bq"], observed["low_bq_eligible_reads"],
                          observed["eligible_reads_after_bq"]), (3, 1, 2))
        self.assertEqual((observed["dp"], observed["ad_ref"], observed["ad_alt"]), (2, 1, 1))
        self.assertEqual(actual["chr2:10:A:C"]["ad_alt"], 1)
        self.assertEqual(diagnostics["saturation_fallback_sites"], 0)

    def test_insertion_candidate_at_read_end_matches_baseline(self):
        # Reads cover positions 81..180. At 180 only one quality remains,
        # whereas the insertion candidate requires two, even for REF reads.
        sites = [dict(variant_id=f"chr1:{pos}:A:AT", chrom="chr1", pos=pos,
                      ref="A", alt="AT") for pos in (179, 180)]
        settings = dict(min_baseq=20, min_mapq=30, max_depth=1000)
        with pysam.AlignmentFile(str(self.bam), "rb") as bam, pysam.FastaFile(str(self.fasta)) as fasta:
            actual, _ = recount_single_pass(bam, fasta, sites, settings)
        with pysam.AlignmentFile(str(self.bam), "rb") as bam, pysam.FastaFile(str(self.fasta)) as fasta:
            for site in sites:
                expected = pileup_site(bam, fasta, site, 20, 30, 1000)
                for field, value in expected.items():
                    self.assertEqual(actual[site['variant_id']][field], value, (site, field))
        self.assertEqual(actual['chr1:179:A:AT']['dp'], 2)
        self.assertEqual(actual['chr1:180:A:AT']['dp'], 0)

    def test_saturation_falls_back_to_pointwise_baseline(self):
        site = dict(variant_id="chr1:101:A:C", chrom="chr1", pos=101, ref="A", alt="C")
        settings = {"min_baseq": 20, "min_mapq": 30, "max_depth": 1,
                    "raw_max_depth": 10000}
        with pysam.AlignmentFile(str(self.bam), "rb") as bam, pysam.FastaFile(str(self.fasta)) as fasta:
            actual, diagnostics = recount_single_pass(bam, fasta, [site], settings)
        with pysam.AlignmentFile(str(self.bam), "rb") as bam, pysam.FastaFile(str(self.fasta)) as fasta:
            expected = pileup_site(bam, fasta, site, 20, 30, 1)
        for field in ("evidence_status", "dp", "ad_ref", "ad_alt", "vaf"):
            self.assertEqual(actual[site["variant_id"]][field], expected[field])
        self.assertEqual(actual[site["variant_id"]]["saturation_fallback"], 1)
        self.assertEqual(diagnostics["saturation_fallback_sites"], 1)

    def test_overlapping_mate_disagreement_is_visible_and_excluded(self):
        site = dict(variant_id="chr1:151:A:C", chrom="chr1", pos=151, ref="A", alt="C")
        settings = {"min_baseq": 20, "min_mapq": 30, "max_depth": 1000,
                    "raw_max_depth": 10000}
        for second, expected in [("A", (0, 0, 1, 1)), ("C", (1, 1, 0, 0))]:
            bam_path = self.root / ("mates_" + second + ".bam")
            header = {"HD": {"VN": "1.6", "SO": "coordinate"},
                      "SQ": [{"SN": "chr1", "LN": 1000}]}
            with pysam.AlignmentFile(str(bam_path), "wb", header=header) as stream:
                for flag, start, mate, allele, tlen in [(99, 100, 120, "C", 120),
                                                        (147, 120, 100, second, -120)]:
                    read = pysam.AlignedSegment();read.query_name = "fragment"
                    read.flag = flag;read.reference_id = 0;read.reference_start = start
                    read.mapping_quality = 60;read.cigartuples = [(0, 100)]
                    read.next_reference_id = 0;read.next_reference_start = mate;read.template_length = tlen
                    sequence = list("A" * 100);sequence[150 - start] = allele
                    read.query_sequence = "".join(sequence);read.query_qualities = [35] * 100
                    stream.write(read)
            pysam.index(str(bam_path))
            with pysam.AlignmentFile(str(bam_path), "rb") as bam, pysam.FastaFile(str(self.fasta)) as fasta:
                row, _ = recount_single_pass(bam, fasta, [site], settings)
            row = row[site["variant_id"]]
            self.assertEqual((row["dp"], row["ad_alt"], row["discordant_eligible_fragments"],
                              row["discordant_raw_fragments"]), expected)
            if second == "C":
                self.assertEqual(row["alt_concordant_both_mates"], 1)

    def test_sample_integration_and_completed_sample_resume(self):
        source = self.root / "empty.vcf"
        source.write_text("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
                          "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        vcf = self.root / "empty.vcf.gz"
        pysam.tabix_compress(str(source), str(vcf), force=True)
        pysam.tabix_index(str(vcf), preset="vcf", force=True)
        site = dict(variant_id="chr1:101:A:C", chrom="chr1", pos=101, ref="A", alt="C",
                    discovered_by="mutect", discovery_sample_count=1)
        sample = dict(sample_id="S", assay_id="IDT", bam_path=str(self.bam),
                      bai_path=str(self.bam) + ".bai")
        config = dict(evidence={"min_baseq": 20, "min_mapq": 30, "max_depth": 1000},
                      single_pass={"raw_max_depth": 10000},
                      pon_by_assay={"IDT": str(vcf)},
                      resources={"fasta": {"path": str(self.fasta)}})
        paths = {caller: str(vcf) for caller in ("mutect", "lofreq", "sage")}
        output = recount_sample_single_pass(sample, [site], config, paths,
                                            {site["variant_id"]: {}}, self.root / "integrated", "fp")
        row = read_tsv(output)[0]
        self.assertEqual((row["dp"], row["ad_ref"], row["ad_alt"]), ("2", "1", "1"))
        self.assertEqual(row["dp_raw_reads"], "5")
        first_mtime = output.stat().st_mtime_ns
        again = recount_sample_single_pass(sample, [site], config, paths,
                                           {site["variant_id"]: {}}, self.root / "integrated", "fp")
        self.assertEqual(again.stat().st_mtime_ns, first_mtime)


if __name__ == "__main__":
    unittest.main()
