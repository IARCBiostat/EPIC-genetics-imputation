#!/bin/bash
#SBATCH --job-name=002_data-reference
#SBATCH --output=src/logs/002_data-reference.out
#SBATCH --error=src/logs/002_data-reference.err
#SBATCH --time=10-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --partition=low_p

# Script: src/002_data-reference.sh
# Purpose: Download the public reference data into ${REF_DIR}:
#   1000G NYGC high-coverage GRCh38 panel, GRCh38 no-alt FASTA, SHAPEIT5 b38 genetic
#   maps, and dbSNP GRCh38. (Stage-1 liftover chain files ship with triple-liftOver,
#   installed by 000_tools.sh.)
# Submit from the repository root: sbatch src/002_data-reference.sh

set -euo pipefail

# Root-only environment sourcing
ENV_FILE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")/.." && pwd)}/.env"
if [ ! -f "$ENV_FILE" ]; then
  echo "ERROR: .env not found at ${ENV_FILE}" >&2; exit 1
fi
# shellcheck disable=SC1090
set -a; source "$ENV_FILE"; set +a
PROJ_ROOT="${GENETICS_PROJECT_ROOT}"

# Defaults (if not in .env)
REF_DIR="${REF_DIR:-${PROJ_ROOT}/data/reference}"

mkdir -p "${REF_DIR}/1000G"
mkdir -p "${REF_DIR}/shapeit5/maps"
mkdir -p "${REF_DIR}/dbsnp"

# Download to <dest>.part and rename on success, so an interrupted or failed download
# never leaves a truncated file (or an HTML error page) that a rerun would skip.
fetch() {
    local url="$1" dest="$2"
    curl -fsSL "$url" -o "${dest}.part"
    mv "${dest}.part" "$dest"
}

echo "=========================================="
echo " Downloading Reference Data (hg38)"
echo " Project Root: ${PROJ_ROOT}"
echo " Reference dir: ${REF_DIR}"
echo "=========================================="

# 1. 1000 Genomes NYGC 2022 high-coverage VCFs (GRCh38)
# Source: EBI FTP
BASE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
echo "Downloading 1000 Genomes NYGC 2022 VCFs (hg38)..."

# ChrX has a slightly different filename pattern (v2)
VCFS=()
for chr in {1..22}; do
    VCFS+=("1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz")
done
VCFS+=("1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz")

for VCF in "${VCFS[@]}"; do
    if [ ! -f "${REF_DIR}/1000G/${VCF}" ]; then
        echo "  Fetching ${VCF}..."
        fetch "${BASE_URL}/${VCF}" "${REF_DIR}/1000G/${VCF}"
    fi
    if [ ! -f "${REF_DIR}/1000G/${VCF}.tbi" ]; then
        fetch "${BASE_URL}/${VCF}.tbi" "${REF_DIR}/1000G/${VCF}.tbi"
    fi
done

# 1b. GRCh38 NO-ALT FASTA (Used by bcftools norm in Nextflow reference_prep module)
FASTA="${REF_DIR}/1000G/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
if [ ! -f "${FASTA}" ]; then
    echo "Downloading GRCh38 FASTA for normalization..."
    curl -fsSL ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gzip -d > "${FASTA}.part"
    mv "${FASTA}.part" "${FASTA}"
fi

# 1c. SHAPEIT5 per-chromosome genetic maps (GRCh38)
# Source: odelaneau/shapeit resources/maps/b38
SHAPEIT5_MAP_DIR="${REF_DIR}/shapeit5/maps"
SHAPEIT5_MAP_TAR_URL="https://github.com/odelaneau/shapeit/raw/main/resources/maps/b38/genetic_maps.b38.tar.gz"
if [ -z "$(ls "${SHAPEIT5_MAP_DIR}"/chr*.b38.gmap.gz 2>/dev/null)" ]; then
    echo "Downloading SHAPEIT5 GRCh38 genetic maps..."
    curl -fsSL "${SHAPEIT5_MAP_TAR_URL}" | tar -xz -C "${SHAPEIT5_MAP_DIR}"
    echo "  ✓ SHAPEIT5 genetic maps"
else
    echo "  ✓ SHAPEIT5 genetic maps (already downloaded)"
fi

# 1d. dbSNP GRCh38 VCF for rsID annotation
# Source: NCBI dbSNP latest_release VCF directory.
# We auto-detect the current GRCh38 assembly file matching GCF_000001405.*.gz.
DBSNP_VCF_DIR_URL="https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/"
DBSNP_FILENAME=""

echo "Resolving latest dbSNP GRCh38 VCF..."
DBSNP_FILENAME="$(curl -fsSL "${DBSNP_VCF_DIR_URL}" 2>/dev/null | tr '"' '\n' | grep -E '^GCF_000001405\.[0-9]+\.gz$' | sort -V | tail -n 1 || true)"
if [ -z "${DBSNP_FILENAME}" ]; then
    echo "  Warning: Could not parse latest dbSNP VCF listing; falling back to GCF_000001405.40.gz"
    DBSNP_FILENAME="GCF_000001405.40.gz"
fi

DBSNP_DEST="${REF_DIR}/dbsnp/${DBSNP_FILENAME}"
if [ ! -f "${DBSNP_DEST}" ]; then
    echo "Downloading dbSNP GRCh38 VCF (${DBSNP_FILENAME})..."
    fetch "${DBSNP_VCF_DIR_URL}/${DBSNP_FILENAME}" "${DBSNP_DEST}"
fi
if [ ! -f "${DBSNP_DEST}.tbi" ]; then
    echo "Downloading dbSNP GRCh38 VCF index..."
    fetch "${DBSNP_VCF_DIR_URL}/${DBSNP_FILENAME}.tbi" "${DBSNP_DEST}.tbi"
fi

echo "Reference data download process complete."
