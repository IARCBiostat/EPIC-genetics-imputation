#!/bin/bash
#SBATCH --job-name=001_data_genetics
#SBATCH --output=src/logs/001_data_genetics.out
#SBATCH --error=src/logs/001_data_genetics.err
#SBATCH --time=10-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --partition=low_p

# Script: src/001_data-genetics.sh
# Purpose: Copy the raw EPIC genetics data this pipeline needs into ${DATA_ROOT}:
#   study genotypes + chip manifests -> ${DATA_ROOT}/genetics/<STUDY>/
#   Subj_Id_2015.txt + genetics*.sas7bdat -> ${REF_DIR}/Epic/
# Usage (from the repository root):
#   bash src/001_data-genetics.sh --check   # only verify .env and every source path
#   sbatch src/001_data-genetics.sh         # verify, then copy

set -euo pipefail
trap 'echo "ERROR: Job failed on line $LINENO" >&2; exit 1' ERR
start_time=$(date +%s)

CHECK_ONLY=0
case "${1:-}" in
  --check) CHECK_ONLY=1 ;;
  "") ;;
  *) echo "Usage: $0 [--check]" >&2; exit 1 ;;
esac

# ── Environment ────────────────────────────────────────────────────────────────
ENV_FILE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")/.." && pwd)}/.env"
if [ ! -f "$ENV_FILE" ]; then
  echo "ERROR: .env not found at ${ENV_FILE}" >&2; exit 1
fi
# shellcheck disable=SC1090
set -a; source "$ENV_FILE"; set +a
PROJ_ROOT="${GENETICS_PROJECT_ROOT}"

# Source locations (see .env). The fallbacks cover a single archive that holds
# Reference/Epic/ and Central_Genetics/ next to the study folders.
SOURCE_ROOT="${GENETICS_DATA_SOURCE_ROOT}"
EPIC_REF_SOURCE_DIR="${EPIC_REF_SOURCE_DIR:-${SOURCE_ROOT}/Reference/Epic}"
CENTRAL_GENETICS_DIR="${CENTRAL_GENETICS_DIR:-${SOURCE_ROOT}/Central_Genetics}"
DEST_ROOT="${DATA_ROOT}/genetics"
EPIC_REF_DEST_ROOT="${REF_DIR}/Epic"

echo "=========================================="
echo " Source locations"
echo " Study data:        ${SOURCE_ROOT}"
echo " EPIC ID reference: ${EPIC_REF_SOURCE_DIR}"
echo " Central genetics:  ${CENTRAL_GENETICS_DIR}"
echo "=========================================="

# ── Check .env before copying anything ────────────────────────────────────────
ENV_ERRORS=0
for var in GENETICS_PROJECT_ROOT GENETICS_DATA_SOURCE_ROOT EPIC_REF_SOURCE_DIR CENTRAL_GENETICS_DIR; do
  val="${!var:-}"
  if [ -z "$val" ] || [[ "$val" == /CHANGE* ]]; then
    echo "ERROR: ${var} is unset or still a placeholder ('${val}')." >&2
    ENV_ERRORS=1
    [ "$var" = GENETICS_PROJECT_ROOT ] || SOURCE_UNSET=1
  fi
done
if [ "${SOURCE_UNSET:-0}" -eq 1 ]; then
  echo "       Set GENETICS_SOURCE_ROOT (or the individual source paths) in .env." >&2
fi
REPO_DIR="$(cd "$(dirname "$ENV_FILE")" && pwd -P)"
if [ "$(cd "$PROJ_ROOT" 2>/dev/null && pwd -P || true)" != "$REPO_DIR" ]; then
  echo "ERROR: GENETICS_PROJECT_ROOT must be this repository's directory." >&2
  echo "       .env has:   ${PROJ_ROOT}" >&2
  echo "       repository: ${REPO_DIR}" >&2
  ENV_ERRORS=1
fi
[ "$ENV_ERRORS" -eq 0 ] || exit 1

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

# Split a STUDIES entry into its fields, applying the folder defaults.
parse_study_entry() {
    IFS="|" read -r STUDY_NAME SUB_PATH DATA_OVERRIDE CHIP_OVERRIDE EXTRA_PLINK_PREFIX EXTRA_PLINK_DEST EXTRA_FILE_PATH EXTRA_FILE_DEST <<< "$1"
    DATA_FLD="${DATA_OVERRIDE:-Data_Received}"
    CHIP_FLD="${CHIP_OVERRIDE:-Chip_files}"
}

# ── Check the source layout ───────────────────────────────────────────────────
# Every path the sync reads must exist, so a wrong source setting in .env fails
# here with the full list rather than partway through a multi-hour copy.
MISSING=()
N_CHECKED=0
require() {
    local kind="$1" path="$2"
    N_CHECKED=$((N_CHECKED + 1))
    if [ "$kind" = dir ] && [ ! -d "$path" ]; then MISSING+=("${path}/"); fi
    if [ "$kind" = file ] && [ ! -f "$path" ]; then MISSING+=("${path}"); fi
}

require file "${EPIC_REF_SOURCE_DIR}/Subj_Id_2015.txt"
require file "${CENTRAL_GENETICS_DIR}/genetics_caco.sas7bdat"
require file "${CENTRAL_GENETICS_DIR}/genetics_id.sas7bdat"
require file "${CENTRAL_GENETICS_DIR}/genetics.sas7bdat"
for entry in "${STUDIES[@]}"; do
    parse_study_entry "$entry"
    if [ "$DATA_FLD" != "." ]; then require dir "${SOURCE_ROOT}/${SUB_PATH}/${DATA_FLD}"; else require dir "${SOURCE_ROOT}/${SUB_PATH}"; fi
    if [[ "$CHIP_FLD" != "." && "$STUDY_NAME" != "Neuro_01" ]]; then require dir "${SOURCE_ROOT}/${SUB_PATH}/${CHIP_FLD}"; fi
    if [[ -n "${EXTRA_PLINK_PREFIX:-}" && -n "${EXTRA_PLINK_DEST:-}" ]]; then
        for ext in bed bim fam; do require file "${SOURCE_ROOT}/${EXTRA_PLINK_PREFIX}.${ext}"; done
    fi
    if [[ -n "${EXTRA_FILE_PATH:-}" && -n "${EXTRA_FILE_DEST:-}" ]]; then require file "${SOURCE_ROOT}/${EXTRA_FILE_PATH}"; fi
done

if [ "${#MISSING[@]}" -gt 0 ]; then
    echo "ERROR: ${#MISSING[@]} of ${N_CHECKED} expected source path(s) are missing:" >&2
    printf '         %s\n' "${MISSING[@]}" >&2
    echo "       Check the source settings in .env (GENETICS_SOURCE_ROOT and the three" >&2
    echo "       locations derived from it) against the layout described there." >&2
    exit 1
fi
echo "All ${N_CHECKED} source paths found."
if [ "$CHECK_ONLY" -eq 1 ]; then
    echo "Check only (--check): nothing copied."
    exit 0
fi

echo "=========================================="
echo " Synchronizing Genetics Data"
echo " From: $SOURCE_ROOT"
echo " To:   $DEST_ROOT"
echo "=========================================="

mkdir -p "$DEST_ROOT"
mkdir -p "$EPIC_REF_DEST_ROOT"

rsync -avP "${EPIC_REF_SOURCE_DIR}/Subj_Id_2015.txt" "${EPIC_REF_DEST_ROOT}/Subj_Id_2015.txt"
rsync -avP "${CENTRAL_GENETICS_DIR}/genetics_caco.sas7bdat" "${EPIC_REF_DEST_ROOT}/genetics_caco.sas7bdat"
rsync -avP "${CENTRAL_GENETICS_DIR}/genetics_id.sas7bdat" "${EPIC_REF_DEST_ROOT}/genetics_id.sas7bdat"
rsync -avP "${CENTRAL_GENETICS_DIR}/genetics.sas7bdat" "${EPIC_REF_DEST_ROOT}/genetics.sas7bdat"

# ── Sync Loop ──────────────────────────────────────────────────────────────────
for entry in "${STUDIES[@]}"; do
    parse_study_entry "$entry"

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
