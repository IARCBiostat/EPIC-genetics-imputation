#!/bin/bash
# Script: src/001_data-genetics.sh
# Purpose: Sync raw genetics data, manifests, and reference files.

set -euo pipefail

# Robust environment sourcing
for env_file in ".env" "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../.env" "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../pipeline/.env"; do
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

# --- Breast ---
echo "Syncing Brea_01_Erneg ..."
mkdir -p "${DEST_ROOT}/Brea_01_Erneg"
# rsync -avP "${SOURCE_ROOT}/Breast/Brea_01_Erneg/Data_Received/" "${DEST_ROOT}/Brea_01_Erneg/Data_Received/"
# rsync -avP "${SOURCE_ROOT}/Breast/Brea_01_Erneg/Chip_files/" "${DEST_ROOT}/Brea_01_Erneg/Chip_files/"

echo "Syncing Brea_02_Onco ..."
mkdir -p "${DEST_ROOT}/Brea_02_Onco"
# rsync -avP "${SOURCE_ROOT}/Breast/Brea_02_Onco/Data_Received/" "${DEST_ROOT}/Brea_02_Onco/Data_Received/"
# rsync -avP "${SOURCE_ROOT}/Breast/Brea_02_Onco/Chip_files/" "${DEST_ROOT}/Brea_02_Onco/Chip_files/"

# --- Colonrectum ---
echo "Syncing Clrt_01_Gecco ..."
mkdir -p "${DEST_ROOT}/Clrt_01_Gecco"
# rsync -avP "${SOURCE_ROOT}/Colonrectum/Clrt_01_Gecco/Data_Received/" "${DEST_ROOT}/Clrt_01_Gecco/Data_Received/"
# rsync -avP "${SOURCE_ROOT}/Colonrectum/Clrt_01_Gecco/Chip_files/" "${DEST_ROOT}/Clrt_01_Gecco/Chip_files/"

# --- Epic_Cvd ---
echo "Syncing Ecvd_01 ..."
mkdir -p "${DEST_ROOT}/Ecvd_01"
# rsync -avP "${SOURCE_ROOT}/Epic_Cvd/Ecvd_01/Data_Received/" "${DEST_ROOT}/Ecvd_01/Data_Received/"
# rsync -avP "${SOURCE_ROOT}/Epic_Cvd/Ecvd_01/Chip_files/" "${DEST_ROOT}/Ecvd_01/Chip_files/"

echo "Syncing Ecvd_02 ..."
mkdir -p "${DEST_ROOT}/Ecvd_02"
# rsync -avP "${SOURCE_ROOT}/Epic_Cvd/Ecvd_02/Data_Received/" "${DEST_ROOT}/Ecvd_02/Data_Received/"
# rsync -avP "${SOURCE_ROOT}/Epic_Cvd/Ecvd_02/Chip_files/" "${DEST_ROOT}/Ecvd_02/Chip_files/"

echo "Syncing Ecvd_03 ..."
mkdir -p "${DEST_ROOT}/Ecvd_03"
# rsync -avP "${SOURCE_ROOT}/Epic_Cvd/Ecvd_03/Data_Received/" "${DEST_ROOT}/Ecvd_03/Data_Received/"
# rsync -avP "${SOURCE_ROOT}/Epic_Cvd/Ecvd_03/Chip_files/" "${DEST_ROOT}/Ecvd_03/Chip_files/"

# --- Gallbladder ---
echo "Syncing Glbd_01 ..."
mkdir -p "${DEST_ROOT}/Glbd_01"
# rsync -avP "${SOURCE_ROOT}/Gallbladder/Glbd_01/Data_Received/" "${DEST_ROOT}/Glbd_01/Data_Received/"
# rsync -avP "${SOURCE_ROOT}/Gallbladder/Glbd_01/Chip_files/" "${DEST_ROOT}/Glbd_01/Chip_files/"

# --- Interact ---
echo "Syncing Inte_01 ..."
mkdir -p "${DEST_ROOT}/Inte_01"
# rsync -avP "${SOURCE_ROOT}/Interact/Inte_01/Data_Received/" "${DEST_ROOT}/Inte_01/Data_Received/"
# rsync -avP "${SOURCE_ROOT}/Interact/Inte_01/Chip_files/" "${DEST_ROOT}/Inte_01/Chip_files/"

echo "Syncing Inte_02 ..."
mkdir -p "${DEST_ROOT}/Inte_02"
# rsync -avP "${SOURCE_ROOT}/Interact/Inte_02/Data_Received/" "${DEST_ROOT}/Inte_02/Data_Received/"
# rsync -avP "${SOURCE_ROOT}/Interact/Inte_02/Chip_files/" "${DEST_ROOT}/Inte_02/Chip_files/"

echo "Syncing Inte_03 ..."
mkdir -p "${DEST_ROOT}/Inte_03"
# rsync -avP "${SOURCE_ROOT}/Interact/Inte_03/Data_Received/" "${DEST_ROOT}/Inte_03/Data_Received/"
# rsync -avP "${SOURCE_ROOT}/Interact/Inte_03/Chip_files/" "${DEST_ROOT}/Inte_03/Chip_files/"

# --- Kidney ---
echo "Syncing Kidn_01 ..."
mkdir -p "${DEST_ROOT}/Kidn_01"
# rsync -avP "${SOURCE_ROOT}/Kidney/Kidn_01/Data_Received/" "${DEST_ROOT}/Kidn_01/Data_Received/"
# rsync -avP "${SOURCE_ROOT}/Kidney/Kidn_01/Chip_files/" "${DEST_ROOT}/Kidn_01/Chip_files/"

echo "Syncing Kidn_02 ..."
mkdir -p "${DEST_ROOT}/Kidn_02"
# rsync -avP "${SOURCE_ROOT}/Kidney/Kidn_02/Data_Received/" "${DEST_ROOT}/Kidn_02/Data_Received/"
# rsync -avP "${SOURCE_ROOT}/Kidney/Kidn_02/Chip_files/" "${DEST_ROOT}/Kidn_02/Chip_files/"

# --- Lung ---
echo "Syncing Lung_01 ..."
mkdir -p "${DEST_ROOT}/Lung_01"
# rsync -avP "${SOURCE_ROOT}/Lung/Lung_01/Data_Received/" "${DEST_ROOT}/Lung_01/Data_Received/"
# rsync -avP "${SOURCE_ROOT}/Lung/Lung_01/Chip_files/" "${DEST_ROOT}/Lung_01/Chip_files/"

# --- Lymphoma ---
echo "Syncing Lymp_01 ..."
mkdir -p "${DEST_ROOT}/Lymp_01"
# rsync -avP "${SOURCE_ROOT}/Lymphoma/Lymp_01/Data_Received/" "${DEST_ROOT}/Lymp_01/Data_Received/"
# rsync -avP "${SOURCE_ROOT}/Lymphoma/Lymp_01/Chip_files/" "${DEST_ROOT}/Lymphoma/Lymp_01/Chip_files/"

# --- Neuro ---
echo "Syncing Neuro_01 ..."
mkdir -p "${DEST_ROOT}/Neuro_01"
# rsync -avP "${SOURCE_ROOT}/Neuro/Neuro_01/Data_Received/" "${DEST_ROOT}/Neuro_01/Data_Received/"

# --- Pancreas ---
echo "Syncing Panc_01_PS1 ..."
mkdir -p "${DEST_ROOT}/Panc_01_PS1"
# rsync -avP "${SOURCE_ROOT}/Pancreas/Panc_01_PS1/Data_Received/" "${DEST_ROOT}/Panc_01_PS1/Data_Received/"
# rsync -avP "${SOURCE_ROOT}/Pancreas/Panc_01_PS1/Chip_files/" "${DEST_ROOT}/Panc_01_PS1/Chip_files/"

echo "Syncing Panc_02_PS3 ..."
mkdir -p "${DEST_ROOT}/Panc_02_PS3"
# rsync -avP "${SOURCE_ROOT}/Pancreas/Panc_02_PS3/Data_Received/" "${DEST_ROOT}/Panc_02_PS3/Data_Received/"
# rsync -avP "${SOURCE_ROOT}/Pancreas/Panc_02_PS3/Chip_files/" "${DEST_ROOT}/Panc_02_PS3/Chip_files/"

# --- Prostate ---
echo "Syncing Pros_01_Bpc3 ..."
mkdir -p "${DEST_ROOT}/Pros_01_Bpc3"
# rsync -avP "${SOURCE_ROOT}/Prostate/Pros_01_Bpc3/Data_Received/" "${DEST_ROOT}/Pros_01_Bpc3/Data_Received/"
# rsync -avP "${SOURCE_ROOT}/Prostate/Pros_01_Bpc3/Chip_files/" "${DEST_ROOT}/Pros_01_Bpc3/Chip_files/"

echo "Syncing Pros_02_Icogs ..."
mkdir -p "${DEST_ROOT}/Pros_02_Icogs"
# --- explicit study rsync commands ---

echo "--- Breast ---"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Breast/Brea_01_Erneg/Data_Received/" "${GENETICS_DIR}/Breast/Brea_01_Erneg/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Breast/Brea_01_Erneg/Chip_files/" "${GENETICS_DIR}/Breast/Brea_01_Erneg/Chip_files/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Breast/Brea_02_Onco/Data_Received/" "${GENETICS_DIR}/Breast/Brea_02_Onco/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Breast/Brea_02_Onco/Chip_files/" "${GENETICS_DIR}/Breast/Brea_02_Onco/Chip_files/"

echo "--- InterAct ---"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Inte_01/Data_Received/" "${GENETICS_DIR}/Inte_01/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Inte_01/Chip_files/" "${GENETICS_DIR}/Inte_01/Chip_files/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Inte_02/Data_Received/" "${GENETICS_DIR}/Inte_02/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Inte_02/Chip_files/" "${GENETICS_DIR}/Inte_02/Chip_files/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Inte_03/Data_Received/" "${GENETICS_DIR}/Inte_03/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Inte_03/Chip_files/" "${GENETICS_DIR}/Inte_03/Chip_files/"

echo "--- Lung & Lymphoma ---"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Lung_01/Data_Received/" "${GENETICS_DIR}/Lung_01/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Lung_01/Chip_files/" "${GENETICS_DIR}/Lung_01/Chip_files/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Lymp_01/Data_Received/" "${GENETICS_DIR}/Lymp_01/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Lymp_01/Chip_files/" "${GENETICS_DIR}/Lymp_01/Chip_files/"

echo "--- Kidney ---"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Kidn_01/Data_Received/" "${GENETICS_DIR}/Kidn_01/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Kidn_01/Chip_files/" "${GENETICS_DIR}/Kidn_01/Chip_files/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Kidn_02/Data_Received/" "${GENETICS_DIR}/Kidn_02/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Kidn_02/Chip_files/" "${GENETICS_DIR}/Kidn_02/Chip_files/"

echo "--- Ovary & Pancreas ---"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Ovar_01/Data_Received_2021/" "${GENETICS_DIR}/Ovar_01/Data_Received_2021/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Ovar_01/Chip_files/" "${GENETICS_DIR}/Ovar_01/Chip_files/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Panc_01_PS1/Data_Received/" "${GENETICS_DIR}/Panc_01_PS1/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Panc_01_PS1/Chip_files/" "${GENETICS_DIR}/Panc_01_PS1/Chip_files/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Panc_02_PS3/Data_Received/" "${GENETICS_DIR}/Panc_02_PS3/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Panc_02_PS3/Chip_files/" "${GENETICS_DIR}/Panc_02_PS3/Chip_files/"

echo "--- Epic CVD ---"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Epic_Cvd/Ecvd_01/Data_Received/" "${GENETICS_DIR}/Epic_Cvd/Ecvd_01/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Epic_Cvd/Ecvd_01/Chip_files/" "${GENETICS_DIR}/Epic_Cvd/Ecvd_01/Chip_files/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Epic_Cvd/Ecvd_02/Data_Received/" "${GENETICS_DIR}/Epic_Cvd/Ecvd_02/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Epic_Cvd/Ecvd_02/Chip_files/" "${GENETICS_DIR}/Epic_Cvd/Ecvd_02/Chip_files/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Epic_Cvd/Ecvd_03/Data_Received/" "${GENETICS_DIR}/Epic_Cvd/Ecvd_03/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Epic_Cvd/Ecvd_03/Chip_files/" "${GENETICS_DIR}/Epic_Cvd/Ecvd_03/Chip_files/"

echo "--- Stomach & UADT ---"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Stomach/Stom_01/Data_Received/" "${GENETICS_DIR}/Stomach/Stom_01/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Stomach/Stom_01/Chip_files/" "${GENETICS_DIR}/Stomach/Stom_01/Chip_files/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Uadt/Uadt_01/Data_Received/" "${GENETICS_DIR}/Uadt/Uadt_01/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Uadt/Uadt_01/Chip_files/" "${GENETICS_DIR}/Uadt/Uadt_01/Chip_files/"

echo "--- Prostate ---"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Prostate/Pros_01_Bpc3/Data_Received/" "${GENETICS_DIR}/Prostate/Pros_01_Bpc3/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Prostate/Pros_01_Bpc3/Chip_files/" "${GENETICS_DIR}/Prostate/Pros_01_Bpc3/Chip_files/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Prostate/Pros_02_Icogs/Data_Received/" "${GENETICS_DIR}/Prostate/Pros_02_Icogs/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Prostate/Pros_02_Icogs/Chip_files/" "${GENETICS_DIR}/Prostate/Pros_02_Icogs/Chip_files/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Prostate/Pros_03_Onco/Data_Received/" "${GENETICS_DIR}/Prostate/Pros_03_Onco/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Prostate/Pros_03_Onco/Chip_files/" "${GENETICS_DIR}/Prostate/Pros_03_Onco/Chip_files/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Prostate/Pros_04_P160555/Data_Received/" "${GENETICS_DIR}/Prostate/Pros_04_P160555/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Prostate/Pros_04_P160555/Chip_files/" "${GENETICS_DIR}/Prostate/Pros_04_P160555/Chip_files/"

echo "--- Colonrectum & Gallbladder ---"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Colonrectum/Clrt_01_Gecco/Data_Received/" "${GENETICS_DIR}/Colonrectum/Clrt_01_Gecco/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Colonrectum/Clrt_01_Gecco/Chip_files/" "${GENETICS_DIR}/Colonrectum/Clrt_01_Gecco/Chip_files/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Gallbladder/Glbd_01/Data_Received/" "${GENETICS_DIR}/Gallbladder/Glbd_01/Data_Received/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Gallbladder/Glbd_01/Chip_files/" "${GENETICS_DIR}/Gallbladder/Glbd_01/Chip_files/"

echo "--- Corpus Uteri & Neuro ---"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Corpus_Uteri/GSA2022_922_025_V3/" "${GENETICS_DIR}/Corpus_Uteri/GSA2022_922_025_V3/"
rsync -avP "${REMOTE_HOST}:${SERVER_GWAS}/Neuro/Neuro_01/Data_Received/" "${GENETICS_DIR}/Neuro/Neuro_01/Data_Received/"

echo "=========================================="
echo " Genetic Data Sync Complete"
echo "=========================================="
