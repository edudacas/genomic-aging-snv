# Fragment-aware evidence pilot — review branch only

This experimental, SNV-only release produces BAM evidence, not countable somatic events. It never launches the cohort. The executor requires exactly three samples from a hash-verified prepared pilot receipt. No age, genotype/carrier status, treatment, phenotype, PoN or germline fields are used for recount or acceptance.

## Observation contract

- Fragment namespace `(RG,QNAME)`. Missing RG is an empty namespace, counted in diagnostics; merged BAM libraries without RG require review and cannot establish molecular independence. No UMI deduplication claim.
- Eligible overlapping mates agree → one observation; disagree → no REF/ALT/OTHER contribution. Raw disagreement involving low-quality observations is diagnosed separately; the eligible mate remains usable.
- Fragment strand anchors to read1; read2-only fragments use the mate-strand flag when the mate is mapped. Unresolved strand is UNKNOWN. Unpaired observations retain alignment strand and UNPAIRED orientation. Conflicting mate orientation/strand is UNKNOWN and explicitly counted.
- Paired F1R2/F2R1 follows read-number/alignment orientation, not BQ. BQ ties and iterator ordering do not change the result.
- Exclude flag mask 0xF04 (unmapped, secondary, QC-fail, duplicate, supplementary). Proposed configurable BQ>=30/MQ>=30. MQ255 means unavailable and is excluded with a diagnostic, not treated as excellent mapping quality.
- REF + ALT + OTHER = usable fragment depth. OTHER includes third alleles and insertion-associated SNV observations. Deletions/refskips at the base are diagnosed but cannot enter the usable base denominator.
- Preserve ALT/(REF+ALT) and ALT/usable depth separately. No VAF boundary or ALT cutoff is applied.
- BQ/MQ summaries use the minimum eligible mate value per fragment, then minimum/median over fragments; read-end summary uses closest end among eligible observations.
- REF mismatch/unsupported allele/absent BAM contig are explicit states with missing counts. Genuine evaluated zero ALT is allowed, but is not zero burden. NO_COVERAGE and NO_USABLE_FRAGMENTS have zero observed usable counts and undefined VAF; no endpoint is inferred.
- Raw depth saturation produces DEPTH_SATURATED with missing final counts. There is no pointwise-fallback/capped-diagnostics mixture.

Strand and orientation are complementary diagnostics. Observing only one category does not reject a variant in this release. Thresholds, denominator choice and interpretation are not scientifically frozen.

## Tests

```bash
python3 -m unittest discover -s experiments/fragment_recount_pilot -p 'test_*.py' -v
```

Requires pysam==0.23.3 for four native BAM/FASTA tests. Skipped native tests do not count as native validation; the real pilot executor refuses any skips. Sixteen pure reducer tests cover overlap, discordance, RG collisions, flags, low/unknown qualities, order invariance, strand/orientation, missing coverage, REF mismatch, saturation, OTHER denominators and VAF-boundary examples.

## Prepare and run locally/on the authorized server

Use the bundled `prepare_pilot.py` to produce `candidate_union.pilot.vcf` and `PREPARATION_RECEIPT.PRIVATE.json`. Supply `--stages`, `--bed`, exactly three repeated `--pilot-sample` arguments, and a new `--out` outside production stages. The common BED must be the authorized production territory. The candidate VCF is a locus/allele list only. Confirm the authorized GRCh38 FASTA and its .fai; contig names must match. The preparer preserves all caller FILTER and VAF states, not just a low-VAF compact export. It also verifies stage/task/source-VCF hashes. Six additional preparer tests and two completion-guard tests are included, bringing the isolated suite to 28 tests.

```bash
python3 experiments/fragment_recount_pilot/run_pilot.py \
  --prepared /private/pilot_inputs \
  --fasta /references/GRCh38.fa \
  --out /private/fragment_pilot_v1 \
  --execute-three-pilot-bams
```

This is a pointwise pilot, potentially slow on network storage; it is deliberately NOT a scalable cohort engine. Full BAM, index and FASTA hashing adds I/O before recount. Inputs are immutable and existing output is not overwritten. Interrupted attempts have no completion receipt; use a new output directory after inspection, not a blind restart.

Outputs are private JSONL evidence, per-pair legacy discordances and a compact validation report. The bundled legacy comparator contains exactly the two audited old allele/pileup functions and is never used for classification. Native count/strand/orientation differences require technical explanations, not forced agreement with the legacy algorithm.

`PILOT_EVIDENCE_COMPLETE_VALIDATION_PENDING` is NOT final validation. Non-evaluable candidate pairs block completion and must be resolved in a subsequent run. Each legacy discrepancy requires an explicitly accepted explanation tied to its exact hash and the validation-report SHA. Review JSON schema:

```json
{
  "report_sha256": "sha256-of-PILOT_VALIDATION_PENDING.PRIVATE.json",
  "reviewer": "authorized technical reviewer",
  "discordance_reviews": {
    "discordance-sha256": {"accepted": true, "explanation": "Observed cause and verification"}
  }
}
```

Then `run_pilot.py --out /private/fragment_pilot_v1 --finalize-with-review /private/review.json` rechecks output hashes. `RECOUNT_PILOT_COMPLETE` is emitted only after tests pass without skips, all three BAMs finish, no candidate evaluation failures remain and all differences are reviewed. If no differences or evaluation failures exist, finalization is automatic. The receipt still explicitly denies cohort execution and somatic classification. Human explanation cannot override a non-evaluable state.

No 5–<35% endpoint, >=5 ALT rule, caller agreement filter, PoN/germline exclusion or callability classifier is applied. No joined positive-event VCF is claimed. Patient data and actual validation receipts must never enter this repository.

## Current validation boundary

Implementation and synthetic results are reviewable. Actual three-production-BAM validation remains pending until its private receipt exists. No merge to main or 1,329-sample launch is authorized.
