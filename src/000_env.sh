#!/bin/bash
# Script: src/000_env.sh
# Purpose: Initialize the Conda environments the pipeline runs from:
#   1. nf_EPIC-genetics: Nextflow + Java, Python 3 (with pandas/matplotlib for the
#      stage-1/3 report figures), and R + haven for 003_data-epic.
#   2. Python 2.7 at ${PYTHON2_ENV} for the legacy stage-1 preprocessing scripts.
# Both are created from lock files (envs/*.lock.txt: exact packages of the June 2026
# run), so every setup gets identical software. Rerunning skips an environment that
# already matches its lock and rebuilds it if the lock has changed.
# Run from the repository root after editing .env: bash src/000_env.sh

set -euo pipefail

ENV_NAME="nf_EPIC-genetics"

# ── Environment ────────────────────────────────────────────────────────────────
ENV_FILE="$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")/.." && pwd)/.env"
if [ ! -f "$ENV_FILE" ]; then
    echo "ERROR: .env not found at ${ENV_FILE}. Run from a clone of the repository, which ships .env." >&2
    exit 1
fi
# shellcheck disable=SC1090
set -a; source "$ENV_FILE"; set +a
REPO_DIR="$(cd "$(dirname "$ENV_FILE")" && pwd -P)"
if [ "$(cd "${GENETICS_PROJECT_ROOT:-}" 2>/dev/null && pwd -P || true)" != "$REPO_DIR" ]; then
    echo "ERROR: GENETICS_PROJECT_ROOT in .env must be this repository's directory." >&2
    echo "       .env has:   ${GENETICS_PROJECT_ROOT:-<unset>}" >&2
    echo "       repository: ${REPO_DIR}" >&2
    exit 1
fi
PYTHON2_ENV="${PYTHON2_ENV:-${GENETICS_PROJECT_ROOT}/.conda/py27}"

echo "=========================================="
echo " Setting up EPIC Genetics Nextflow Environment"
echo "=========================================="

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda is not available in your PATH."
    exit 1
fi

CONDA_CMD="conda"
if command -v mamba >/dev/null 2>&1; then
    CONDA_CMD="mamba"
fi

# Create the environment at <prefix> from <lock> (explicit package URLs: no solving,
# no channels). Skipped if it was already built from this exact lock.
env_from_lock() {
    local prefix="$1" lock="$2" stamp want
    stamp="${prefix}/.epic_lock.md5"
    want="$(md5sum "$lock" | cut -d' ' -f1)"
    if [ -f "$stamp" ] && [ "$(cat "$stamp")" = "$want" ]; then
        echo "  ✓ ${prefix} (matches $(basename "$lock"))"
        return 0
    fi
    if [ -e "$prefix" ]; then
        if [ ! -d "${prefix}/conda-meta" ]; then
            echo "ERROR: ${prefix} exists but is not a conda environment; move it out of the way." >&2
            exit 1
        fi
        echo "  Rebuilding ${prefix} from $(basename "$lock")..."
        rm -rf "$prefix"
    else
        echo "  Creating ${prefix} from $(basename "$lock")..."
    fi
    $CONDA_CMD create -y -p "$prefix" --file "$lock"
    echo "$want" > "$stamp"
}

env_from_lock "$(conda info --base)/envs/${ENV_NAME}" "${REPO_DIR}/envs/nf_EPIC-genetics.lock.txt"
echo ""
env_from_lock "${PYTHON2_ENV}" "${REPO_DIR}/envs/py27.lock.txt"
"${PYTHON2_ENV}/bin/python2.7" --version

echo "=========================================="
echo "Setup complete. Activate with: 'conda activate ${ENV_NAME}'"
echo "Python 2.7: ${PYTHON2_ENV}/bin/python2.7 (PYTHON2_BIN in .env)"
echo "=========================================="
