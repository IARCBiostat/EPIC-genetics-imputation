#!/bin/bash
#SBATCH --job-name=001_data_genetics
#SBATCH --output=src/logs/001_data_genetics.out
#SBATCH --error=src/logs/001_data_genetics.err
#SBATCH --time=10-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --partition=low_p

set -euo pipefail
trap 'echo "ERROR: Job failed on line $LINENO" >&2; exit 1' ERR
start_time=$(date +%s)

# ── Environment ────────────────────────────────────────────────────────────────
ENV_FILE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")/.." && pwd)}/.env"
if [ ! -f "$ENV_FILE" ]; then
  echo "ERROR: .env not found at ${ENV_FILE}" >&2; exit 1
fi
# shellcheck disable=SC1090
set -a; source "$ENV_FILE"; set +a
PROJ_ROOT="${GENETICS_PROJECT_ROOT}"

SOURCE_ROOT="${GENETICS_DATA_SOURCE_ROOT}"
DEST_ROOT="${DATA_ROOT}/genetics"
EPIC_REF_DEST_ROOT="${REF_DIR}/Epic"

echo "=========================================="
echo " Synchronizing Genetics Data"
echo " From: $SOURCE_ROOT"
echo " To:   $DEST_ROOT"
echo "=========================================="

mkdir -p "$DEST_ROOT"
mkdir -p "$EPIC_REF_DEST_ROOT"

rsync -avP "${SOURCE_ROOT}/Reference/Epic/Subj_Id_2015.txt" "${EPIC_REF_DEST_ROOT}/Subj_Id_2015.txt"
rsync -avP "${SOURCE_ROOT}/Central_Genetics/genetics_caco.sas7bdat" "${EPIC_REF_DEST_ROOT}/genetics_caco.sas7bdat"
rsync -avP "${SOURCE_ROOT}/Central_Genetics/genetics_id.sas7bdat" "${EPIC_REF_DEST_ROOT}/genetics_id.sas7bdat"
rsync -avP "${SOURCE_ROOT}/Central_Genetics/genetics.sas7bdat" "${EPIC_REF_DEST_ROOT}/genetics.sas7bdat"

# ── Study Definition Mapping ──────────────────────────────────────────────────
# Format:
# "StudyFolderName"|"StudySubPath"|"DataFolderName"|"ChipFolderName"|"ExtraPlinkPrefix"|"ExtraPlinkDest"|"ExtraFilePath"|"ExtraFileDest"
# Optional overrides:
#   DataFolderName: defaults to Data_Received, use "." for flat sync
#   ChipFolderName: defaults to Chip_files, use "." to skip chip sync
#   ExtraPlinkPrefix: source prefix relative to SOURCE_ROOT, synced as .bed/.bim/.fam
#   ExtraPlinkDest: destination directory relative to the study target directory
#   ExtraFilePath: source file relative to SOURCE_ROOT
#   ExtraFileDest: destination file relative to the study target directory
STUDIES=(
    "Brea_01_Erneg|Breast/Brea_01_Erneg"
    "Brea_02_Onco|Breast/Brea_02_Onco"
    "Clrt_01_Gecco|Colonrectum/Clrt_01_Gecco|Data_Received|Chip_files|Extraction/Clrt_Data_Extracted/clrt_gecco_geno|Data_Received/Data_Extracted"
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

    "Pros_03_Onco|Prostate/Pros_03_Onco"
    "Pros_04_P160555|Prostate/Pros_04_P160555"
    "Ovar_01|Ovary/Ovar_01|Data_Received_2022|Chip_files|||Ovary/Ovar_01/Data_Received_2021/Link_Ids_Ovar_01Onco.csv|Data_Received_2021/Link_Ids_Ovar_01Onco.csv"
    "Stom_01|Stomach/Stom_01"
    "Uadt_01|Uadt/Uadt_01"
)

# ── Sync Loop ──────────────────────────────────────────────────────────────────
for entry in "${STUDIES[@]}"; do
    IFS="|" read -r STUDY_NAME SUB_PATH DATA_OVERRIDE CHIP_OVERRIDE EXTRA_PLINK_PREFIX EXTRA_PLINK_DEST EXTRA_FILE_PATH EXTRA_FILE_DEST <<< "$entry"
    DATA_FLD="${DATA_OVERRIDE:-Data_Received}"
    CHIP_FLD="${CHIP_OVERRIDE:-Chip_files}"
    
    echo "--- Syncing ${STUDY_NAME} ---"
    
    # Create target paths
    TGT_DIR="${DEST_ROOT}/${STUDY_NAME}"
    mkdir -p "$TGT_DIR"

    # 1. Sync Data folder
    if [ "$DATA_FLD" != "." ]; then
        rsync -avP "${SOURCE_ROOT}/${SUB_PATH}/${DATA_FLD}/" "${TGT_DIR}/${DATA_FLD}/"
    else
        # Flat structure (DATA_FLD=.)
        rsync -avP "${SOURCE_ROOT}/${SUB_PATH}/" "${TGT_DIR}/"
    fi

    # 2. Sync Chip files (Skip if already handled by flat sync or explicitly excluded)
    if [[ "$CHIP_FLD" != "." && "$STUDY_NAME" != "Neuro_01" ]]; then
       rsync -avP "${SOURCE_ROOT}/${SUB_PATH}/${CHIP_FLD}/" "${TGT_DIR}/${CHIP_FLD}/"
    fi

    # 3. Overlay study-specific raw PLINK prefixes when they live outside the study folder
    if [[ -n "${EXTRA_PLINK_PREFIX:-}" && -n "${EXTRA_PLINK_DEST:-}" ]]; then
        mkdir -p "${TGT_DIR}/${EXTRA_PLINK_DEST}"
        for ext in bed bim fam; do
            rsync -avP "${SOURCE_ROOT}/${EXTRA_PLINK_PREFIX}.${ext}" "${TGT_DIR}/${EXTRA_PLINK_DEST}/"
        done
    fi

    # 4. Sync extra study-specific single files when needed
    if [[ -n "${EXTRA_FILE_PATH:-}" && -n "${EXTRA_FILE_DEST:-}" ]]; then
        mkdir -p "$(dirname "${TGT_DIR}/${EXTRA_FILE_DEST}")"
        rsync -avP "${SOURCE_ROOT}/${EXTRA_FILE_PATH}" "${TGT_DIR}/${EXTRA_FILE_DEST}"
    fi
done

echo "=========================================="
echo " Genetic Data Sync Complete"
echo "=========================================="
