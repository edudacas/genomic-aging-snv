"""Audit SAGE-phased MNV/component representations without changing the union."""
import hashlib


def components(chrom, pos, ref, alt):
    if len(ref) != len(alt) or len(ref) < 2:
        return []
    return [f"{chrom}:{pos + index}:{r}:{a}" for index, (r, a) in enumerate(zip(ref, alt)) if r != a]


def _lps(record):
    value = record.get("info", {}).get("LPS", [])
    if value in (None, ""):
        return set()
    if not isinstance(value, list):
        value = [value]
    return {str(item) for item in value}


def build_mnv_groups(sample_id, sage_records, union):
    """One row per union MNV; only exact components with a shared LPS can verify."""
    rows = []
    candidates = []
    for key, site in union.items():
        component_keys = components(site["chrom"], int(site["pos"]), site["ref"], site["alt"])
        if len(component_keys) < 2:
            continue
        base = dict(sample_id=sample_id, mnv_variant_id=key,
                    component_variant_ids=";".join(component_keys),
                    components_in_union=int(all(component in union for component in component_keys)),
                    representation_source="SAGE_LPS_EXACT_COMPONENTS")
        mnv = sage_records.get(key)
        if not base["components_in_union"]:
            rows.append(dict(base, representation_status="UNRESOLVED_COMPONENT_NOT_IN_UNION",
                             lps="", event_group_id=""))
            continue
        if not mnv or any(component not in sage_records for component in component_keys):
            rows.append(dict(base, representation_status="UNRESOLVED_SAGE_RECORD_MISSING",
                             lps="", event_group_id=""))
            continue
        mnv_lps = set().union(*(_lps(record) for record in mnv["records"]))
        shared = set(mnv_lps)
        for component in component_keys:
            shared &= set().union(*(_lps(record) for record in sage_records[component]["records"]))
        if not shared:
            rows.append(dict(base, representation_status="UNRESOLVED_NO_SHARED_LPS",
                             lps="", event_group_id=""))
            continue
        if len(shared) != 1:
            rows.append(dict(base, representation_status="UNRESOLVED_MULTIPLE_SHARED_LPS",
                             lps=";".join(sorted(shared)), event_group_id=""))
            continue
        lps = next(iter(shared))
        candidates.append(dict(base, representation_status="VERIFIED_PHASED_COMPONENTS",
                               lps=lps, event_group_id=hashlib.sha256(
                                   (sample_id + "\0" + key + "\0" + lps).encode()).hexdigest()[:24],
                               _component_set=set(component_keys)))
    # Same local phase set can contain several variants; overlapping MNV
    # representations are ambiguous and are never collapsed automatically.
    overlaps = [any(index != other and row["lps"] == candidate["lps"] and
                    row["_component_set"] & candidate["_component_set"]
                    for other, candidate in enumerate(candidates))
                for index, row in enumerate(candidates)]
    for row, overlap in zip(candidates, overlaps):
        if overlap:
            row.update(representation_status="UNRESOLVED_OVERLAPPING_MNV_REPRESENTATIONS",
                       event_group_id="")
        row.pop("_component_set")
        rows.append(row)
    return sorted(rows, key=lambda row: row["mnv_variant_id"])
