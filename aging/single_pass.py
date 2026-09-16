"""One sequential BAM traversal with baseline-compatible filtered evidence.

The iterator is deliberately raw. Read-level eligibility is applied before the
unchanged v2 fragment reducer sees a synthetic column. This permits additional
diagnostics without a second BAM traversal. Saturated candidate positions are
recounted pointwise with the v2 implementation and are reported as fallbacks.
"""
from collections import defaultdict
from statistics import median

from .evidence import allele_from_pileup, pileup_site

FILTER_FLAGS = 0xF04


class ColumnView:
    def __init__(self, column, pileups=None):
        self.reference_pos = column.reference_pos
        self.pileups = list(column.pileups if pileups is None else pileups)
        self.nsegments = len(self.pileups)


class OneColumnBam:
    def __init__(self, column):
        self.column = column

    def pileup(self, *args, **kwargs):
        return (self.column,)


def _orientation(read):
    if not read.is_paired:
        return "UNPAIRED"
    return "F1R2" if (read.is_read1 != read.is_reverse) else "F2R1"


def _fragment_item(pileup, ref, alt):
    allele = allele_from_pileup(pileup, ref, alt)
    if allele is None:
        return None
    read = pileup.alignment
    qpos = pileup.query_position
    length = len(alt) if len(alt) > len(ref) else 1
    qualities = read.query_qualities[qpos:qpos + length]
    if len(qualities) != length:
        return None
    unclipped_start = read.reference_start - (read.query_alignment_start or 0)
    distance_to_end = min(qpos, len(read.query_sequence) - qpos - 1)
    nm = read.get_tag("NM") if read.has_tag("NM") else None
    return dict(
        allele=allele,
        min_bq=min(qualities),
        mq=read.mapping_quality,
        strand="R" if read.is_reverse else "F",
        orientation=_orientation(read),
        query_name=read.query_name,
        read_group=read.get_tag("RG") if read.has_tag("RG") else "",
        read_number=1 if read.is_read1 else 2 if read.is_read2 else 0,
        unclipped_start=unclipped_start,
        distance_to_end=distance_to_end,
        near_end=int(distance_to_end < 10),
        nm=nm,
        abs_tlen=abs(read.template_length) if read.is_paired and read.template_length else None,
    )


def _eligible_fragments(pileups, ref, alt, min_bq):
    groups = defaultdict(list)
    for pileup in pileups:
        item = _fragment_item(pileup, ref, alt)
        if item is not None and item["min_bq"] >= min_bq:
            # Preserve the v2 reducer exactly for the real-baseline comparison.
            groups[item["query_name"]].append(item)
    retained = []
    discordant = 0
    concordant_alt_both_mates = 0
    for values in groups.values():
        alleles = {value["allele"] for value in values}
        if len(alleles) != 1:
            discordant += 1
            continue
        chosen = max(values, key=lambda value: value["min_bq"])
        retained.append(chosen)
        if chosen["allele"] == alt and {value["read_number"] for value in values} >= {1, 2}:
            concordant_alt_both_mates += 1
    return retained, discordant, concordant_alt_both_mates


def _discordant_fragments(items):
    groups = defaultdict(set)
    for item in items:
        groups[item["query_name"]].add(item["allele"])
    return sum(len(alleles) > 1 for alleles in groups.values())


def _summary(values, prefix, result):
    if not values:
        for name in ("min_bq", "median_bq", "sum_bq", "min_mq", "median_mq", "sum_mq",
                     "median_distance_to_end", "near_end_fraction", "mean_nm", "mean_abs_tlen"):
            result[prefix + name] = ""
        result[prefix + "unique_fragment_starts"] = 0
        return
    for field in ("bq", "mq"):
        source = [value["min_bq" if field == "bq" else "mq"] for value in values]
        result[prefix + "min_" + field] = min(source)
        result[prefix + "median_" + field] = median(source)
        result[prefix + "sum_" + field] = sum(source)
    distances = [value["distance_to_end"] for value in values]
    result[prefix + "median_distance_to_end"] = median(distances)
    result[prefix + "near_end_fraction"] = sum(value["near_end"] for value in values) / len(values)
    result[prefix + "unique_fragment_starts"] = len({(value["read_group"], value["unclipped_start"], value["strand"]) for value in values})
    for field in ("nm", "abs_tlen"):
        observed = [value[field] for value in values if value[field] is not None]
        result[prefix + "mean_" + field] = sum(observed) / len(observed) if observed else ""


def extended_evidence(raw_column, filtered_pileups, fasta, site, settings):
    """Return unchanged v2 evidence plus additive, explicitly named diagnostics."""
    min_bq = settings["min_baseq"]
    filtered_column = ColumnView(raw_column, filtered_pileups)
    baseline = pileup_site(OneColumnBam(filtered_column), fasta, site, min_bq,
                           settings["min_mapq"], settings["max_depth"])
    raw_items = [item for item in (_fragment_item(p, site["ref"], site["alt"])
                                   for p in raw_column.pileups) if item is not None]
    eligible_before_bq = [p for p in filtered_pileups
                          if _fragment_item(p, site["ref"], site["alt"]) is not None]
    eligible_items_before_bq = [_fragment_item(p, site["ref"], site["alt"])
                                for p in eligible_before_bq]
    eligible_items_after_bq = [item for item in eligible_items_before_bq
                               if item["min_bq"] >= min_bq]
    retained, discordant, both_alt = _eligible_fragments(filtered_pileups, site["ref"], site["alt"], min_bq)
    ref = [value for value in retained if value["allele"] == site["ref"]]
    alt = [value for value in retained if value["allele"] == site["alt"]]
    for field, count in (("dp", len(retained)), ("ad_ref", len(ref)), ("ad_alt", len(alt))):
        if baseline[field] != count:
            raise AssertionError(
                f"Fragment count mismatch at {site['variant_id']}: "
                f"{field} baseline={baseline[field]} extended={count}")
    result = dict(baseline)
    result.update(
        dp_raw_reads=len(raw_items),
        alt_raw_reads=sum(value["allele"] == site["alt"] for value in raw_items),
        mq0_raw_reads=sum(value["mq"] == 0 for value in raw_items),
        duplicate_raw_reads=sum(p.alignment.is_duplicate for p in raw_column.pileups),
        raw_depth_saturated=int(raw_column.nsegments >= settings.get("raw_max_depth", 100000)),
        filtered_depth_saturated=int(len(filtered_pileups) >= settings["max_depth"]),
        flag_filtered_raw_reads=sum(bool(p.alignment.flag & FILTER_FLAGS) for p in raw_column.pileups),
        low_mq_raw_reads=sum(p.alignment.mapping_quality < settings["min_mapq"] for p in raw_column.pileups),
        low_bq_eligible_reads=sum(item["min_bq"] < min_bq for item in eligible_items_before_bq),
        eligible_reads_before_bq=len(eligible_items_before_bq),
        eligible_reads_after_bq=len(eligible_items_after_bq),
        other_fragments=len(retained) - len(ref) - len(alt),
        discordant_raw_fragments=_discordant_fragments(raw_items),
        discordant_eligible_fragments=discordant,
        alt_concordant_both_mates=both_alt,
        ref_forward=sum(value["strand"] == "F" for value in ref),
        ref_reverse=sum(value["strand"] == "R" for value in ref),
        ref_f1r2=sum(value["orientation"] == "F1R2" for value in ref),
        ref_f2r1=sum(value["orientation"] == "F2R1" for value in ref),
    )
    _summary(ref, "ref_", result)
    _summary(alt, "alt_ext_", result)
    if baseline["dp"]:
        result["coverage_diagnostic"] = "FILTERED_FRAGMENTS_OBSERVED"
    elif not raw_items:
        result["coverage_diagnostic"] = "NO_ALIGNED_ALLELE_OBSERVATIONS"
    elif not eligible_items_before_bq:
        result["coverage_diagnostic"] = "NO_READS_AFTER_FLAG_AND_MQ"
    elif not eligible_items_after_bq:
        result["coverage_diagnostic"] = "NO_READS_AFTER_BQ"
    else:
        result["coverage_diagnostic"] = "NO_FRAGMENT_AFTER_MATE_RECONCILIATION"
    return result


def recount_single_pass(bam, fasta, sites, settings, progress=None):
    """Traverse the coordinate-sorted BAM once; random access only for saturation fallback."""
    by_position = defaultdict(list)
    for site in sites:
        by_position[(site["chrom"], int(site["pos"]) - 1)].append(site)
    result = {}
    saturated = []
    raw_cap = settings.get("raw_max_depth", max(settings["max_depth"], 100000))
    observed_columns = 0
    for column in bam.pileup(until_eof=True, stepper="nofilter", min_base_quality=0,
                             max_depth=raw_cap, ignore_overlaps=False,
                             ignore_orphans=False, compute_baq=False):
        rows = by_position.get((bam.get_reference_name(column.reference_id), column.reference_pos))
        if not rows:
            continue
        observed_columns += 1
        eligible = [pileup for pileup in column.pileups
                    if not (pileup.alignment.flag & FILTER_FLAGS)
                    and pileup.alignment.mapping_quality >= settings["min_mapq"]]
        filtered_column = ColumnView(column, eligible)
        for site in rows:
            result[site["variant_id"]] = extended_evidence(column, eligible, fasta, site, settings)
            if filtered_column.nsegments >= settings["max_depth"] or column.nsegments >= raw_cap:
                saturated.append(site)
                continue
        if progress:
            progress(rows[-1], len(result), observed_columns)
    # Missing candidate columns retain the v2 absence semantics; diagnostics are explicit zeros.
    empty_fields = dict(dp_raw_reads=0, alt_raw_reads=0, mq0_raw_reads=0,
                        duplicate_raw_reads=0, raw_depth_saturated=0,
                        filtered_depth_saturated=0, flag_filtered_raw_reads=0,
                        low_mq_raw_reads=0, low_bq_eligible_reads=0,
                        eligible_reads_before_bq=0, eligible_reads_after_bq=0,
                        other_fragments=0, discordant_raw_fragments=0,
                        discordant_eligible_fragments=0, alt_concordant_both_mates=0,
                        ref_forward=0, ref_reverse=0, ref_f1r2=0, ref_f2r1=0,
                        coverage_diagnostic="NO_ALIGNED_ALLELE_OBSERVATIONS")
    for site in sites:
        if site["variant_id"] in result:
            continue
        baseline = pileup_site(OneColumnBam(ColumnView(type("Empty", (), {"reference_pos": int(site["pos"]) - 1, "pileups": []})())),
                               fasta, site, settings["min_baseq"], settings["min_mapq"], settings["max_depth"])
        result[site["variant_id"]] = dict(baseline, **empty_fields)
    # Saturation is rare but traversal-start-dependent; retain exact pointwise v2 behavior.
    for site in saturated:
        baseline = pileup_site(bam, fasta, site, settings["min_baseq"], settings["min_mapq"], settings["max_depth"])
        result[site["variant_id"]].update(baseline, saturation_fallback=1)
    for site in sites:
        row = result[site["variant_id"]]
        if row["evidence_status"] == "UNSUPPORTED_COMPLEX_ALLELE":
            row.update(dp="", ad_ref="", ad_alt="", vaf="")
    return result, dict(sequential_columns_with_candidates=observed_columns,
                        saturation_fallback_sites=len(saturated), raw_max_depth=raw_cap)
