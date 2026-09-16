#!/usr/bin/env bash
set -euo pipefail
# This helper is invoked explicitly; it never starts the production container.
: "${AGING_RUN:?Set AGING_RUN to the existing frozen run directory}"
: "${AGING_EXPORTS:?Set AGING_EXPORTS to a directory outside AGING_RUN}"
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="$SCRIPT_DIR${PYTHONPATH:+:$PYTHONPATH}"
exec python3 -m aging.stage_export "${1:?Usage: bash export_stage.sh sage}" \
  --run "$AGING_RUN" --destination "$AGING_EXPORTS"
