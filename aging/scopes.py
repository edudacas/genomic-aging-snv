"""Config-frozen site scopes derived only from discovery fragment evidence."""
import hashlib
import json
from collections import defaultdict


FORCED_ALL = "FORCED_ALL"
COHORT_LOW_VAF = "COHORT_LOW_VAF"
DISCOVERY_PLUS_FIXED_CONTROLS = "DISCOVERY_PLUS_FIXED_CONTROLS"


def validate_scope_config(config, sample_assays, require_frozen=False):
    minimum = config.get("min_alt_fragments")
    maximum = config.get("max_discovery_vaf")
    if not isinstance(minimum, int) or minimum < 1:
        raise ValueError("min_alt_fragments must be a positive integer")
    if not isinstance(maximum, (int, float)) or not 0 < maximum < 0.5:
        raise ValueError("max_discovery_vaf must be between 0 and 0.5")
    controls = config.get("fixed_controls_by_assay", {})
    expected_n = config.get("fixed_control_n_by_assay", {})
    if not isinstance(controls, dict) or not isinstance(expected_n, dict):
        raise ValueError("fixed controls and N must be dictionaries keyed by assay")
    external = config.get("external_control_assays", {})
    if not isinstance(external, dict):
        raise ValueError("external_control_assays must map external IDs to assays")
    if set(external) & set(sample_assays):
        raise ValueError("External controls overlap run sample IDs")
    control_assays = dict(sample_assays, **external)
    all_samples = set(control_assays)
    for assay, ids in controls.items():
        if len(ids) != len(set(ids)):
            raise ValueError("duplicate fixed control: " + assay)
        if expected_n.get(assay) != len(ids):
            raise ValueError("fixed control N does not match IDs: " + assay)
        for sample_id in ids:
            if sample_id not in all_samples or control_assays[sample_id] != assay:
                raise ValueError("fixed control assay/sample mismatch: " + sample_id)
    if require_frozen:
        if config.get("frozen") is not True:
            raise ValueError("cohort scope gate is not frozen")
        assays = set(sample_assays.values())
        if set(controls) != assays or set(expected_n) != assays:
            raise ValueError("every cohort assay requires frozen controls and N")
        if any(expected_n[assay] < 1 for assay in assays):
            raise ValueError("every assay requires at least one fixed control")
    return hashlib.sha256(json.dumps(config, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def discovery_summary(rows, minimum, maximum):
    """Rows must contain only samples that discovered the exact normalized allele."""
    observed = [row for row in rows if row.get("evidence_status") in ("OBSERVED", "DEPTH_SATURATED")]
    def numeric(row, name, cast=float):
        value = row.get(name)
        return None if value in (None, "", ".") else cast(value)
    alt_values = [numeric(row, "ad_alt", int) for row in observed]
    vaf_values = [numeric(row, "vaf") for row in observed]
    alt_values = [value for value in alt_values if value is not None]
    vaf_values = [value for value in vaf_values if value is not None]
    qualifying = []
    alt_ge = []
    for row in observed:
        alt = numeric(row, "ad_alt", int)
        vaf = numeric(row, "vaf")
        if alt is not None and alt >= minimum:
            alt_ge.append(row["sample_id"])
            if vaf is not None and vaf < maximum:
                qualifying.append(row["sample_id"])
    return dict(
        discovery_max_ad_alt=max(alt_values) if alt_values else "",
        discovery_max_vaf=max(vaf_values) if vaf_values else "",
        discovery_n_samples_alt_ge_min=len(set(alt_ge)),
        discovery_n_samples_scope_eligible=len(set(qualifying)),
        discovery_scope_eligible_sample_ids=";".join(sorted(set(qualifying))),
        discovery_evidence_observed_samples=len({row["sample_id"] for row in observed}),
        discovery_evidence_no_coverage_samples=len({row["sample_id"] for row in rows if row.get("evidence_status") == "NO_COVERAGE"}),
        discovery_evidence_unsupported_samples=len({row["sample_id"] for row in rows if row.get("evidence_status") == "UNSUPPORTED_COMPLEX_ALLELE"}),
        discovery_evidence_saturated_samples=len({row["sample_id"] for row in rows if row.get("evidence_status") == "DEPTH_SATURATED"}),
    )


def assign_scopes(union, evidence_rows, discovery_samples, forced_sources, config,
                  sample_assays, require_frozen=False):
    """Return one auditable row per union allele; no caller filter/AF/PoN is read."""
    config_hash = validate_scope_config(config, sample_assays, require_frozen)
    by_key = defaultdict(list)
    for row in evidence_rows:
        if row["sample_id"] in set(discovery_samples.get(row["variant_id"], ())):
            by_key[row["variant_id"]].append(row)
    forced = {}
    for source, keys in forced_sources.items():
        for key in keys:
            forced.setdefault(key, []).append(source)
    result = []
    for key, site in union.items():
        discoverers = sorted(set(discovery_samples.get(key, ())))
        summary = discovery_summary(by_key[key], config["min_alt_fragments"],
                                    config["max_discovery_vaf"])
        if key in forced:
            scope = FORCED_ALL
            reason = "forced_source=" + ";".join(sorted(forced[key]))
        elif summary["discovery_n_samples_scope_eligible"]:
            scope = COHORT_LOW_VAF
            reason = ("same_discovery_sample_ad_alt_ge_%d_and_vaf_lt_%s" %
                      (config["min_alt_fragments"], config["max_discovery_vaf"]))
        else:
            scope = DISCOVERY_PLUS_FIXED_CONTROLS
            reason = "no_discovery_sample_met_joint_low_vaf_support_rule"
        planned = set(discoverers)
        if scope in (FORCED_ALL, COHORT_LOW_VAF):
            planned.update(sample_assays)
        else:
            for ids in config.get("fixed_controls_by_assay", {}).values():
                planned.update(sample_id for sample_id in ids if sample_id in sample_assays)
        external_ids = sorted({
            sample_id for ids in config.get("fixed_controls_by_assay", {}).values()
            for sample_id in ids if sample_id in config.get("external_control_assays", {})
        }) if scope == DISCOVERY_PLUS_FIXED_CONTROLS else []
        result.append(dict(site, variant_id=key, recount_scope=scope,
                           recount_scope_reason=reason,
                           recount_scope_config_sha256=config_hash,
                           discovery_sample_ids=";".join(discoverers),
                           planned_recount_sample_ids=";".join(sorted(planned)),
                           external_control_reference_ids=";".join(external_ids),
                           **summary))
    return result
