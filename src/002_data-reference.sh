#!/bin/bash
# Script: src/002_data-reference.sh
# Purpose: Download 1000 Genomes reference data into data/reference/

set -euo pipefail

# Robust environment sourcing
for env_file in ".env" "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../.env" "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../pipeline/.env"; do
  if [ -f "$env_file" ]; then
    set -a; source "$env_file"; set +a
    break
  fi
done

REF_DIR="data/reference"
mkdir -p "${REF_DIR}/1000G"
mkdir -p "${REF_DIR}/annovar_db"
mkdir -p "${REF_DIR}/Chain_LiftOver"

echo "=========================================="
echo " Downloading Reference Data (hg38)"
echo "=========================================="

# 1. 1000 Genomes phase3 VCFs (GRCh38)
# Source: EBI FTP
BASE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL"
echo "Downloading 1000 Genomes Phase 3 VCFs (hg38)..."
for chr in {1..22} X; do
    VCF="ALL.chr${chr}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"
    if [ ! -f "${REF_DIR}/1000G/${VCF}" ]; then
        echo "  Fetching Chromosome ${chr}..."
        wget -q "${BASE_URL}/${VCF}" -P "${REF_DIR}/1000G/"
        wget -q "${BASE_URL}/${VCF}.tbi" -P "${REF_DIR}/1000G/"
    fi
done

# 2. LiftOver Chain (hg19 to hg38)
if [ ! -f "${REF_DIR}/Chain_LiftOver/hg19ToHg38.over.chain.gz" ]; then
    echo "Downloading LiftOver Chain..."
    wget -q http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -P "${REF_DIR}/Chain_LiftOver/"
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
