import unittest

from aging.mnv_groups import build_mnv_groups, components


def record(lps):
    return {"records": [{"info": {"LPS": [lps]}}]}


class MnvGroupTests(unittest.TestCase):
    def test_exact_components_and_shared_phase_are_required(self):
        union = {
            "chr1:10:AC:GT": dict(chrom="chr1", pos=10, ref="AC", alt="GT"),
            "chr1:10:A:G": dict(chrom="chr1", pos=10, ref="A", alt="G"),
            "chr1:11:C:T": dict(chrom="chr1", pos=11, ref="C", alt="T"),
        }
        sage = {key: record(7) for key in union}
        row = build_mnv_groups("S", sage, union)[0]
        self.assertEqual(row["representation_status"], "VERIFIED_PHASED_COMPONENTS")
        self.assertTrue(row["event_group_id"])
        self.assertEqual(components("chr1", 10, "AC", "GT"),
                         ["chr1:10:A:G", "chr1:11:C:T"])
        sage["chr1:11:C:T"] = record(8)
        row = build_mnv_groups("S", sage, union)[0]
        self.assertEqual(row["representation_status"], "UNRESOLVED_NO_SHARED_LPS")
        self.assertEqual(row["event_group_id"], "")

    def test_missing_component_and_overlapping_representation_stay_unresolved(self):
        union = {
            "chr1:10:AC:GT": dict(chrom="chr1", pos=10, ref="AC", alt="GT"),
            "chr1:10:A:G": dict(chrom="chr1", pos=10, ref="A", alt="G"),
        }
        row = build_mnv_groups("S", {}, union)[0]
        self.assertEqual(row["representation_status"], "UNRESOLVED_COMPONENT_NOT_IN_UNION")
        union.update({
            "chr1:11:C:T": dict(chrom="chr1", pos=11, ref="C", alt="T"),
            "chr1:10:ACA:GTA": dict(chrom="chr1", pos=10, ref="ACA", alt="GTA"),
        })
        sage = {key: record(9) for key in union}
        rows = build_mnv_groups("S", sage, union)
        self.assertTrue(all(row["representation_status"] == "UNRESOLVED_OVERLAPPING_MNV_REPRESENTATIONS" for row in rows))


if __name__ == "__main__":
    unittest.main()
