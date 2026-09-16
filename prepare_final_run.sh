#!/usr/bin/env bash
set -euo pipefail
umask 077
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
[[ $# == 1 ]] || { echo 'Usage: bash prepare_final_run.sh NEW_RUN_DIRECTORY'; exit 2; }
RUN=$(realpath -m "$1")
[[ "$RUN" == /srv/genomic-aging/projects/genomic_aging_lfs/* ]] || exit 2
[[ ! -e "$RUN" ]] || { echo 'Run directory exists; preserving it.' >&2; exit 2; }
PROJECT=/srv/genomic-aging/projects/genomic_aging_lfs
BENCHMARK="$PROJECT/bqr_benchmark_v3_5/benchmark"
WORKBOOK="$PROJECT/inputs/TP53_LFS_WES_truthset.HUMAN_REVIEW.EN.FULL1329.private.simple.xlsx"
BASE_CONFIG="$PROJECT/preflight_v3_3_1/config.json"
ENSEMBL=/srv/genomic-aging/projects/genomic_aging_validation/pilot_workers4_v1/ensembl.manifest.tsv
for FILE in "$WORKBOOK" "$BASE_CONFIG" "$ENSEMBL" "$BENCHMARK/bqr_feasibility.json" \
  "$BENCHMARK/enabled/PILOT_A/FAILURE.json" \
  "$BENCHMARK/enabled/PILOT_A/commands.log" study_linkage.PRIVATE.tsv; do
  [[ -s "$FILE" ]] || { echo "Missing input: $FILE" >&2; exit 2; }
done
python3 finalize_bqr_feasibility.py --benchmark "$BENCHMARK"
install -d -m 700 "$RUN/bqr_evidence"
cp -- "$WORKBOOK" "$RUN/study.xlsx"
cp -- "$BASE_CONFIG" "$RUN/base.config.json"
cp -- "$ENSEMBL" "$RUN/ensembl.manifest.tsv"
cp -- study_linkage.PRIVATE.tsv "$RUN/linkage.tsv"
cp -- docs/ANALYSIS_PLAN_GENOMIC_AGING_LFS.v3.1.md "$RUN/ANALYSIS_PLAN_GENOMIC_AGING_LFS.v3.1.md"
cp -- docs/SPEC_SEX_CHROMOSOME_MOSAIC.md "$RUN/SPEC_SEX_CHROMOSOME_MOSAIC.md"
cp -- "$BENCHMARK/bqr_feasibility.json" "$RUN/bqr_decision.json"
cp -- "$BENCHMARK/enabled/PILOT_A/FAILURE.json" "$RUN/bqr_evidence/FAILURE.json"
cp -- "$BENCHMARK/enabled/PILOT_A/commands.log" "$RUN/bqr_evidence/commands.log"
python3 prepare_cohort_config.py --base "$RUN/base.config.json" --decision "$RUN/bqr_decision.json" \
  --evidence-root "$RUN/bqr_evidence" --release "$PWD" --out "$RUN" \
  --workers 12 --cpus 16 --memory-gib 96
sha256sum "$RUN/study.xlsx" "$RUN/linkage.tsv" "$RUN/ensembl.manifest.tsv" \
  "$RUN/base.config.json" "$RUN/bqr_decision.json" "$RUN/bqr_evidence/FAILURE.json" \
  "$RUN/bqr_evidence/commands.log" "$RUN/config.json" "$RUN/caller_protocol_unchanged.json" \
  "$RUN/ANALYSIS_PLAN_GENOMIC_AGING_LFS.v3.1.md" "$RUN/SPEC_SEX_CHROMOSOME_MOSAIC.md" \
  > "$RUN/PREPARED_INPUTS.sha256"
chmod 600 "$RUN/PREPARED_INPUTS.sha256"
echo "Prepared only; NOT launched: $RUN"
