#!/bin/bash
# Script: src/000_env.sh
# Purpose: Initialize a minimal Conda environment to run Nextflow.

set -euo pipefail

ENV_NAME="nf_master"

echo "=========================================="
echo " Setting up Nextflow Master Environment"
echo "=========================================="

if ! command -v conda >/dev/null 2>&1 && ! command -v mamba >/dev/null 2>&1; then
    echo "ERROR: Neither conda nor mamba is available in your PATH."
    exit 1
fi

CONDA_CMD="conda"
if command -v mamba >/dev/null 2>&1; then
    CONDA_CMD="mamba"
fi

echo "Creating environment '${ENV_NAME}'..."
$CONDA_CMD create -n ${ENV_NAME} -c conda-forge -c bioconda nextflow openjdk=21 python=3.9 -y

echo "=========================================="
echo "Setup complete. Activate with: 'conda activate ${ENV_NAME}'"
echo "=========================================="
