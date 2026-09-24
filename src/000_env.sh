#!/bin/bash
# Script: src/000_env.sh
# Purpose: Initialize the Conda environments the pipeline runs from:
#   1. nf_EPIC-genetics: Nextflow + Java, Python 3 (with pandas/matplotlib for the
#      stage-1/3 report figures), and R + haven for 003_data-epic.
#   2. Python 2.7 at ${PYTHON2_ENV} for the legacy stage-1 preprocessing scripts.
# Run from the repository root after creating .env: bash src/000_env.sh

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
PYTHON2_ENV="${PYTHON2_ENV:-${GENETICS_PROJECT_ROOT}/.conda/py27}"

echo "=========================================="
echo " Setting up EPIC Genetics Nextflow Environment"
echo "=========================================="

if ! command -v conda >/dev/null 2>&1 && ! command -v mamba >/dev/null 2>&1; then
    echo "ERROR: Neither conda nor mamba is available in your PATH."
    exit 1
fi

CONDA_CMD="conda"
if command -v mamba >/dev/null 2>&1; then
    CONDA_CMD="mamba"
fi

# --override-channels keeps the Anaconda 'defaults' channel (and its terms-of-service
# prompt) out of the solve; everything needed is on conda-forge/bioconda.
echo "Creating environment '${ENV_NAME}'..."
$CONDA_CMD create -n ${ENV_NAME} --override-channels -c conda-forge -c bioconda \
    nextflow openjdk=21 python=3.9 pandas matplotlib-base r-base r-haven -y

echo ""
if [ -x "${PYTHON2_ENV}/bin/python2.7" ]; then
    echo "Python 2.7 environment already present at ${PYTHON2_ENV}"
else
    echo "Creating Python 2.7 environment at ${PYTHON2_ENV}..."
    $CONDA_CMD create -p "${PYTHON2_ENV}" --override-channels -c conda-forge python=2.7 -y
fi
"${PYTHON2_ENV}/bin/python2.7" --version

echo "=========================================="
echo "Setup complete. Activate with: 'conda activate ${ENV_NAME}'"
echo "Python 2.7: ${PYTHON2_ENV}/bin/python2.7 (PYTHON2_BIN in .env)"
echo "=========================================="
