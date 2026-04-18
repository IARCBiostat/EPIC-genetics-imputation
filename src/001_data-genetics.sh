#!/bin/bash
# Script: src/001_data-genetics.sh
# Purpose: Sync raw genetics data, manifests, and reference files (Local Sync).

set -euo pipefail
trap 'echo "ERROR: Job failed on line $LINENO" >&2; exit 1' ERR
start_time=$(date +%s)

# ── Environment ────────────────────────────────────────────────────────────────
for env_file in ".env" \
                "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../.env" \
                "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../pipeline_stage2/.env"; do
  if [ -f "$env_file" ]; then
    set -a; source "$env_file"; set +a
    break
  fi
done

SOURCE_ROOT="/data/Epic/subprojects/Genetics/sources/Gwas"
DEST_ROOT="data/genetics"

echo "=========================================="
echo " Synchronizing Genetics Data"
echo " From: $SOURCE_ROOT"
echo " To:   $DEST_ROOT"
echo "=========================================="

mkdir -p "$DEST_ROOT"

# ── Study Definition Mapping ──────────────────────────────────────────────────
# Format: "StudyFolderName"|"StudySubPath"|"DataFolderName" (Optional)
STUDIES=(
    "Brea_01_Erneg|Breast/Brea_01_Erneg"
    "Brea_02_Onco|Breast/Brea_02_Onco"
    "Clrt_01_Gecco|Colonrectum/Clrt_01_Gecco"
    "Ecvd_01|Epic_Cvd/Ecvd_01"
    "Ecvd_02|Epic_Cvd/Ecvd_02"
    "Ecvd_03|Epic_Cvd/Ecvd_03"
    "Glbd_01|Gallbladder/Glbd_01"
    "Inte_01|Interact/Inte_01"
    "Inte_02|Interact/Inte_02"
    "Inte_03|Interact/Inte_03"
    "Kidn_01|Kidney/Kidn_01"
    "Kidn_02|Kidney/Kidn_02"
    "Lung_01|Lung/Lung_01"
    "Lymp_01|Lymphoma/Lymp_01"
    "Neuro_01|Neuro/Neuro_01"
    "Panc_01_PS1|Pancreas/Panc_01_PS1"
    "Panc_02_PS3|Pancreas/Panc_02_PS3"
    "Pros_01_Bpc3|Prostate/Pros_01_Bpc3"
    "Pros_02_Icogs|Prostate/Pros_02_Icogs"
    "Pros_03_Onco|Prostate/Pros_03_Onco"
    "Pros_04_P160555|Prostate/Pros_04_P160555"
    "Ovar_01|Ovar_01|Data_Received_2021"
    "Stom_01|Stomach/Stom_01"
    "Uadt_01|Uadt/Uadt_01"
    "Corpus_Uteri|Corpus_Uteri/GSA2022_922_025_V3/|.|."
)

# ── Sync Loop ──────────────────────────────────────────────────────────────────
for entry in "${STUDIES[@]}"; do
    IFS="|" read -r STUDY_NAME SUB_PATH DATA_OVERRIDE <<< "$entry"
    DATA_FLD="${DATA_OVERRIDE:-Data_Received}"
    
    echo "--- Syncing ${STUDY_NAME} ---"
    
    # Create target paths
    TGT_DIR="${DEST_ROOT}/${STUDY_NAME}"
    mkdir -p "$TGT_DIR"

    # 1. Sync Data folder
    if [ "$DATA_FLD" != "." ]; then
        rsync -avP "${SOURCE_ROOT}/${SUB_PATH}/${DATA_FLD}/" "${TGT_DIR}/${DATA_FLD}/"
    else
        # Flat structure (e.g. Corpus Uteri)
        rsync -avP "${SOURCE_ROOT}/${SUB_PATH}/" "${TGT_DIR}/"
    fi

    # 2. Sync Chip files (Skip if already handled by flat sync or explicitly excluded)
    if [[ "$DATA_FLD" != "." && "$STUDY_NAME" != "Neuro_01" ]]; then
       rsync -avP "${SOURCE_ROOT}/${SUB_PATH}/Chip_files/" "${TGT_DIR}/Chip_files/"
    fi
done

echo "=========================================="
echo " Genetic Data Sync Complete"
echo "=========================================="
