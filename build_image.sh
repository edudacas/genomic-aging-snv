#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")"
sha256sum -c SHA256SUMS
EXPECTED_BASE=sha256:5d5d1d14cd10d77fa38aa82bb5fe2cabb432c6447784d81a6d15fe43dae01f67
[[ $(docker image inspect genomic-aging-lfs:3.0-callers-permfix2 --format '{{.Id}}') == "$EXPECTED_BASE" ]] || {
  echo 'Base image changed; stop for review.' >&2; exit 2;
}
docker build --network none --pull=false -t genomic-aging-lfs:3.5.2-staged-final .
