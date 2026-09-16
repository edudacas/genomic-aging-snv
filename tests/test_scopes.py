import copy
import unittest

from aging.scopes import (COHORT_LOW_VAF, DISCOVERY_PLUS_FIXED_CONTROLS,
                          FORCED_ALL, assign_scopes, validate_scope_config)


class ScopeTests(unittest.TestCase):
    def setUp(self):
        self.union = {
            "chr1:10:A:C": dict(chrom="chr1", pos=10, ref="A", alt="C"),
            "chr1:20:A:C": dict(chrom="chr1", pos=20, ref="A", alt="C"),
            "chr1:30:A:C": dict(chrom="chr1", pos=30, ref="A", alt="C"),
            "chr1:40:A:C": dict(chrom="chr1", pos=40, ref="A", alt="C"),
        }
        self.samples = {"case1": "IDT", "case2": "IDT", "control": "IDT"}
        self.discovery = {key: ["case1"] for key in self.union}
        self.rows = [
            dict(sample_id="case1", variant_id="chr1:10:A:C", evidence_status="OBSERVED", ad_alt="9", vaf="0.40", mutect_filter="PASS", gnomad_af="0"),
            dict(sample_id="case1", variant_id="chr1:20:A:C", evidence_status="OBSERVED", ad_alt="2", vaf="0.349", mutect_filter="filtered", gnomad_af="0.2"),
            dict(sample_id="case1", variant_id="chr1:30:A:C", evidence_status="OBSERVED", ad_alt="1", vaf="0.01", mutect_filter="PASS", gnomad_af="0"),
            dict(sample_id="case1", variant_id="chr1:40:A:C", evidence_status="OBSERVED", ad_alt="2", vaf="0.35", mutect_filter="PASS", gnomad_af="0"),
        ]
        self.config = dict(min_alt_fragments=2, max_discovery_vaf=0.35,
                           fixed_controls_by_assay={"IDT": ["control"]},
                           fixed_control_n_by_assay={"IDT": 1}, frozen=True)

    def test_scope_priority_boundaries_and_no_filter_or_af_dependence(self):
        result = assign_scopes(self.union, self.rows, self.discovery,
                               {"drivers": {"chr1:10:A:C"}}, self.config,
                               self.samples, require_frozen=True)
        by_key = {row["variant_id"]: row for row in result}
        self.assertEqual(by_key["chr1:10:A:C"]["recount_scope"], FORCED_ALL)
        self.assertEqual(by_key["chr1:20:A:C"]["recount_scope"], COHORT_LOW_VAF)
        self.assertEqual(by_key["chr1:30:A:C"]["recount_scope"], DISCOVERY_PLUS_FIXED_CONTROLS)
        self.assertEqual(by_key["chr1:40:A:C"]["recount_scope"], DISCOVERY_PLUS_FIXED_CONTROLS)
        self.assertEqual(by_key["chr1:20:A:C"]["planned_recount_sample_ids"], "case1;case2;control")
        self.assertEqual(by_key["chr1:30:A:C"]["planned_recount_sample_ids"], "case1;control")
        changed = copy.deepcopy(self.rows)
        for row in changed:
            row["mutect_filter"] = "opposite"
            row["gnomad_af"] = "0.999"
            row["in_pon"] = "1"
        again = assign_scopes(self.union, changed, self.discovery,
                              {"drivers": {"chr1:10:A:C"}}, self.config,
                              self.samples, require_frozen=True)
        self.assertEqual([(r["variant_id"], r["recount_scope"]) for r in result],
                         [(r["variant_id"], r["recount_scope"]) for r in again])

    def test_joint_condition_not_independent_maxima(self):
        discovery = {"chr1:20:A:C": ["case1", "case2"]}
        rows = [dict(sample_id="case1", variant_id="chr1:20:A:C", evidence_status="OBSERVED", ad_alt="9", vaf="0.40"),
                dict(sample_id="case2", variant_id="chr1:20:A:C", evidence_status="OBSERVED", ad_alt="1", vaf="0.01")]
        result = assign_scopes({"chr1:20:A:C": self.union["chr1:20:A:C"]}, rows,
                               discovery, {}, self.config, self.samples, True)[0]
        self.assertEqual(result["discovery_max_ad_alt"], 9)
        self.assertEqual(result["discovery_max_vaf"], 0.40)
        self.assertEqual(result["discovery_n_samples_scope_eligible"], 0)
        self.assertEqual(result["recount_scope"], DISCOVERY_PLUS_FIXED_CONTROLS)

    def test_external_validated_controls_are_not_run_samples(self):
        config = copy.deepcopy(self.config)
        config["fixed_controls_by_assay"] = {"IDT": ["external"]}
        config["external_control_assays"] = {"external": "IDT"}
        result = assign_scopes(self.union, self.rows, self.discovery, {}, config,
                               self.samples, True)
        row = next(r for r in result if r["variant_id"] == "chr1:30:A:C")
        self.assertEqual(row["planned_recount_sample_ids"], "case1")
        self.assertEqual(row["external_control_reference_ids"], "external")
        config["external_control_assays"] = {"case1": "IDT"}
        with self.assertRaisesRegex(ValueError, "overlap"):
            validate_scope_config(config, self.samples, True)

    def test_frozen_gate_requires_n_and_ids(self):
        broken = copy.deepcopy(self.config)
        broken["fixed_control_n_by_assay"]["IDT"] = 2
        with self.assertRaises(ValueError):
            validate_scope_config(broken, self.samples, True)
        broken = copy.deepcopy(self.config)
        broken["frozen"] = False
        with self.assertRaises(ValueError):
            validate_scope_config(broken, self.samples, True)


if __name__ == "__main__":
    unittest.main()
