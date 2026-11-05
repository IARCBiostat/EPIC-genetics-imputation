#!/bin/bash
# Script: src/misc/export_env_locks.sh
# Purpose: Record the exact software of a finished, known-good run as conda lock files
#   (explicit package URLs + md5), so the repository can pin every environment to it.
#   Exports each Nextflow conda environment of stages 2-4, the host environments
#   (nf_EPIC-genetics incl. Nextflow, renv, Python 2.7), and a SUMMARY.txt of the key
#   tool versions and the dbSNP build.
# Usage:
#   bash src/misc/export_env_locks.sh <run-dir> [<out-dir>]
#     <run-dir>  the run's ${SCRATCH}/${SCRATCH_DATE}, holding stage2/ stage3/ stage4/
#   Optional environment variables for things outside <run-dir>:
#     PY27_PREFIX  the Python 2.7 env used by stage 1 (e.g. <project root>/.conda/py27)
#     TOOLS_DIR    the run's tools directory (records compiled/container tool versions)
#     DBSNP_DIR    the run's reference/dbsnp directory (records the dbSNP build)

set -euo pipefail

RUN_DIR="${1:?Usage: bash src/misc/export_env_locks.sh <run-dir> [<out-dir>]}"
OUT_DIR="${2:-env_locks_$(date +%Y%m%d)}"
[ -d "$RUN_DIR" ] || { echo "ERROR: run directory not found: $RUN_DIR" >&2; exit 1; }
command -v conda >/dev/null 2>&1 || { echo "ERROR: conda is not on PATH" >&2; exit 1; }
mkdir -p "$OUT_DIR"
SUMMARY="$OUT_DIR/SUMMARY.txt"
: > "$SUMMARY"

KEY_PKGS='^(nextflow|openjdk|python|bcftools|htslib|samtools|shapeit5|minimac4|plink|plink2|numpy|pandas|matplotlib-base|seaborn|scipy|r-base|r-haven|r-tidyverse)$'

# Write <out>/<name>.lock.txt for the conda environment at <prefix>, plus a summary line.
export_env() {
    local prefix="$1" name="$2"
    if [ ! -d "$prefix/conda-meta" ]; then
        echo "$name: not found ($prefix)" >> "$SUMMARY"
        return 0
    fi
    conda list --explicit --md5 -p "$prefix" > "$OUT_DIR/$name.lock.txt"
    printf '%s: %s\n' "$name" "$(conda list -p "$prefix" | awk -v re="$KEY_PKGS" '$1 ~ re {printf "%s=%s(%s) ", $1, $2, $4}')" >> "$SUMMARY"
}

echo "== Nextflow conda environments (${RUN_DIR})" >> "$SUMMARY"
for stage in stage2 stage3 stage4; do
    for env in "$RUN_DIR/$stage/conda"/env-*; do
        [ -d "$env" ] && export_env "$env" "${stage}__$(basename "$env")"
    done
done

echo "== Host environments" >> "$SUMMARY"
CONDA_BASE="$(conda info --base)"
export_env "$CONDA_BASE/envs/nf_EPIC-genetics" "host__nf_EPIC-genetics"
export_env "$CONDA_BASE/envs/renv" "host__renv"
[ -n "${PY27_PREFIX:-}" ] && export_env "$PY27_PREFIX" "host__py27"

if [ -n "${TOOLS_DIR:-}" ] && [ -d "$TOOLS_DIR/bin" ]; then
    echo "== Tools in ${TOOLS_DIR}/bin" >> "$SUMMARY"
    for tool in bcftools plink2 plink; do
        [ -x "$TOOLS_DIR/bin/$tool" ] && echo "$tool: $("$TOOLS_DIR/bin/$tool" --version 2>&1 | head -1)" >> "$SUMMARY"
    done
    ls "$TOOLS_DIR/singularity_images" >> "$SUMMARY" 2>/dev/null || true
fi

if [ -n "${DBSNP_DIR:-}" ] && [ -d "$DBSNP_DIR" ]; then
    echo "== dbSNP in ${DBSNP_DIR}" >> "$SUMMARY"
    for vcf in "$DBSNP_DIR"/GCF_000001405.*.gz; do
        [ -f "$vcf" ] || continue
        # Read only the header; gzip is cut off by head, so do not let pipefail see it.
        header="$(set +o pipefail; gzip -dc "$vcf" 2>/dev/null | head -200)"
        build="$(printf '%s\n' "$header" | grep -m1 -E '^##dbSNP_BUILD_ID=' || true)"
        echo "$(basename "$vcf"): ${build:-dbSNP_BUILD_ID not found} ($(du -h "$vcf" | cut -f1))" >> "$SUMMARY"
    done
fi

echo "Wrote $(ls "$OUT_DIR"/*.lock.txt 2>/dev/null | wc -l | tr -d ' ') lock files and SUMMARY.txt to $OUT_DIR"
cat "$SUMMARY"
