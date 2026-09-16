#!/usr/bin/env bash
set -euo pipefail
umask 077
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
[[ $# == 3 ]] || { echo 'Usage: bash run_stage.sh qc|sage|primary|mutect|lofreq PREPARED_RUN CONTAINER'; exit 2; }
STAGE=$1
case "$STAGE" in qc|sage|primary|mutect|lofreq) ;; *) exit 2 ;; esac
RUN=$(realpath "$2")
CONTAINER=$3
[[ "$RUN" == /srv/genomic-aging/projects/genomic_aging_lfs/* ]] || exit 2
[[ "$CONTAINER" =~ ^genomic-aging-[a-z0-9-]+$ ]] || exit 2
[[ -z $(docker ps -aq --filter "name=^/${CONTAINER}$") ]] || {
  echo 'Container already exists. Inspect it; this script never removes or restarts it.' >&2; exit 2;
}
sha256sum -c SHA256SUMS
for FILE in config.json study.xlsx linkage.tsv ensembl.manifest.tsv bqr_decision.json \
  bqr_evidence/FAILURE.json bqr_evidence/commands.log caller_protocol_unchanged.json \
  ANALYSIS_PLAN_GENOMIC_AGING_LFS.v3.1.md SPEC_SEX_CHROMOSOME_MOSAIC.md PREPARED_INPUTS.sha256; do
  [[ -s "$RUN/$FILE" ]] || { echo "Missing prepared input: $FILE" >&2; exit 2; }
done
(cd "$RUN" && sha256sum --quiet -c PREPARED_INPUTS.sha256)
read -r CPUS MEMORY WORKERS < <(python3 - "$RUN/config.json" <<'PY'
import json,sys
from aging.stages import resource_budget
c=json.load(open(sys.argv[1]));resource_budget(c)
assert c['execution_mode']=='STAGED_CALLERS_ONLY'
assert c.get('bqr_decision',{}).get('confirmed_by_user') is True
print(c['cpus'],c['memory_gib'],c['sample_workers'])
PY
)
[[ -n "$CPUS" && -n "$MEMORY" && -n "$WORKERS" ]] || exit 2
install -d -m 700 "$RUN/tmp"
docker run -d --name "$CONTAINER" --network none --read-only \
  --cpus "$CPUS" --memory "${MEMORY}g" --memory-swap "${MEMORY}g" --user "$(id -u):$(id -g)" \
  --tmpfs /tmp:rw,size=1g,mode=1777 -e "AGING_CPU_LIMIT=$CPUS" -e TMPDIR=/work/tmp \
  -e JAVA_TOOL_OPTIONS=-XX:ActiveProcessorCount=2 \
  -v "$RUN:/work" -v "$RUN/config.json:/work/config.json:ro" \
  -v "$RUN/study.xlsx:/work/study.xlsx:ro" -v "$RUN/linkage.tsv:/work/linkage.tsv:ro" \
  -v "$RUN/ensembl.manifest.tsv:/work/ensembl.manifest.tsv:ro" \
  -v /srv/genomic-aging/refs:/resources:ro \
  -v /srv/genomic-aging/BAMs:/srv/genomic-aging/BAMs:ro \
  genomic-aging-lfs:3.5.2-staged-final "$STAGE" --workbook /work/study.xlsx --linkage /work/linkage.tsv \
  --config /work/config.json --out /work/stages
