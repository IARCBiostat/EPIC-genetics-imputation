#!/bin/bash
# Script: src/004_run.sh
# Purpose: Local Nextflow run wrapper.

set -euo pipefail

# Robust environment sourcing
for env_file in ".env" "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../.env" "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../pipeline/.env"; do
  if [ -f "$env_file" ]; then
    set -a; source "$env_file"; set +a
    break
  fi
done

PROFILE=${1:-local}
PIPELINE_DIR="pipeline"

cd "$PIPELINE_DIR"
nextflow run main.nf -profile "$PROFILE" -resume
