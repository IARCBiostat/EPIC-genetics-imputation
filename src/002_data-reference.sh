#!/bin/bash
# Script: src/002_data-reference.sh
# Purpose: Download 1000 Genomes reference data into data/reference/

set -euo pipefail

# Root-only environment sourcing
SCRIPT_DIR="$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")" && pwd)"
PROJ_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
export GENETICS_PROJECT_ROOT="${PROJ_ROOT}"
ENV_FILE="${PROJ_ROOT}/.env"
if [ ! -f "$ENV_FILE" ]; then
  echo "ERROR: Root environment file not found: ${ENV_FILE}" >&2
  echo "       Create ${ENV_FILE}; stage-specific .env files are no longer supported." >&2
  exit 1
fi
set -a
# shellcheck disable=SC1090
source "$ENV_FILE"
set +a
export GENETICS_PROJECT_ROOT="${PROJ_ROOT}"

# Defaults (if not in .env)
TOOLS_DIR="${TOOLS_DIR:-${PROJ_ROOT}/tools}"
REF_DIR="${REF_DIR:-${PROJ_ROOT}/data/reference}"

mkdir -p "${REF_DIR}/1000G"
mkdir -p "${REF_DIR}/eagle"
mkdir -p "${REF_DIR}/dbsnp"
mkdir -p "${REF_DIR}/annovar_db"
mkdir -p "${REF_DIR}/Chain_LiftOver"

echo "=========================================="
echo " Downloading Reference Data (hg38)"
echo " Project Root: ${PROJ_ROOT}"
echo "=========================================="

# 1. 1000 Genomes NYGC 2022 high-coverage VCFs (GRCh38)
# Source: EBI FTP
BASE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"
echo "Downloading 1000 Genomes NYGC 2022 VCFs (hg38)..."

for chr in {1..22}; do
    VCF="1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    if [ ! -f "${REF_DIR}/1000G/${VCF}" ]; then
        echo "  Fetching Chromosome ${chr}..."
        curl -sL "${BASE_URL}/${VCF}" -o "${REF_DIR}/1000G/${VCF}"
        curl -sL "${BASE_URL}/${VCF}.tbi" -o "${REF_DIR}/1000G/${VCF}.tbi"
    fi
done

# ChrX has a slightly different filename pattern (v2)
VCF_X="1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
if [ ! -f "${REF_DIR}/1000G/${VCF_X}" ]; then
    echo "  Fetching Chromosome X..."
    curl -sL "${BASE_URL}/${VCF_X}" -o "${REF_DIR}/1000G/${VCF_X}"
    curl -sL "${BASE_URL}/${VCF_X}.tbi" -o "${REF_DIR}/1000G/${VCF_X}.tbi"
fi

# 1b. GRCh38 NO-ALT FASTA (Used by bcftools norm in Nextflow reference_prep module)
if [ ! -f "${REF_DIR}/1000G/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" ]; then
    echo "Downloading GRCh38 FASTA for normalization..."
    curl -sL ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | gzip -d > "${REF_DIR}/1000G/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
fi

# 1c. Eagle hg38 recombination map with chrX
# Source: Broad Eagle download tables
EAGLE_MAP_DEST="${REF_DIR}/eagle/genetic_map_hg38_withX.txt.gz"
if [ ! -f "${EAGLE_MAP_DEST}" ]; then
    echo "Downloading Eagle hg38 genetic map with chrX..."
    EAGLE_MAP_URLS=(
        "https://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz"
        "http://data.broadinstitute.org/alkesgroup/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz"
    )

    downloaded=0
    for url in "${EAGLE_MAP_URLS[@]}"; do
        if curl -fsSL "${url}" -o "${EAGLE_MAP_DEST}"; then
            downloaded=1
            break
        fi
    done

    if [ "${downloaded}" -ne 1 ]; then
        rm -f "${EAGLE_MAP_DEST}"
        echo "ERROR: Failed to download Eagle genetic map from Broad." >&2
        exit 1
    fi
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
    curl -fsSL "${DBSNP_VCF_DIR_URL}/${DBSNP_FILENAME}" -o "${DBSNP_DEST}"
fi
if [ ! -f "${DBSNP_DEST}.tbi" ]; then
    echo "Downloading dbSNP GRCh38 VCF index..."
    curl -fsSL "${DBSNP_VCF_DIR_URL}/${DBSNP_FILENAME}.tbi" -o "${DBSNP_DEST}.tbi"
fi

# 2. LiftOver Chain (hg19 to hg38)
if [ ! -f "${REF_DIR}/Chain_LiftOver/hg19ToHg38.over.chain.gz" ]; then
    echo "Downloading LiftOver Chain..."
    curl -sL http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -o "${REF_DIR}/Chain_LiftOver/hg19ToHg38.over.chain.gz"
fi

# 3. Annovar hg38 Database
if [ -f "${TOOLS_DIR}/annovar/annotate_variation.pl" ]; then
    echo "Downloading Annovar hg38 databases (refGene, avsnp150)..."
    perl "${TOOLS_DIR}/annovar/annotate_variation.pl" -buildver hg38 -downdb -webfrom annovar avsnp150 "${REF_DIR}/annovar_db/"
    perl "${TOOLS_DIR}/annovar/annotate_variation.pl" -buildver hg38 -downdb -webfrom annovar refGene "${REF_DIR}/annovar_db/"
else
    echo "Warning: Annovar tool not found at ${TOOLS_DIR}/annovar/annotate_variation.pl."
fi

echo "Reference data download process complete."
