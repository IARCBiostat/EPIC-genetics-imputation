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
# Source: odelaneau/shapeit resources/maps/b38, pinned to a commit and checksum so the
# maps cannot change under the pipeline.
SHAPEIT5_MAP_DIR="${REF_DIR}/shapeit5/maps"
SHAPEIT5_MAP_COMMIT="f9d726472df3f26120fe738a137fc41cbcc0bbf4"
SHAPEIT5_MAP_SHA256="04f97acc6524d75e1dbc397e72cf6776b2f1e33f72a06ee477ef69359a97c69e"
SHAPEIT5_MAP_TAR_URL="https://raw.githubusercontent.com/odelaneau/shapeit/${SHAPEIT5_MAP_COMMIT}/resources/maps/b38/genetic_maps.b38.tar.gz"
if [ -z "$(ls "${SHAPEIT5_MAP_DIR}"/chr*.b38.gmap.gz 2>/dev/null)" ]; then
    echo "Downloading SHAPEIT5 GRCh38 genetic maps..."
    MAP_TAR="${SHAPEIT5_MAP_DIR}/genetic_maps.b38.tar.gz"
    fetch "${SHAPEIT5_MAP_TAR_URL}" "${MAP_TAR}"
    if [ "$(sha256sum "${MAP_TAR}" | cut -d' ' -f1)" != "${SHAPEIT5_MAP_SHA256}" ]; then
        echo "ERROR: SHAPEIT5 genetic maps checksum mismatch (${MAP_TAR})" >&2
        exit 1
    fi
    tar -xzf "${MAP_TAR}" -C "${SHAPEIT5_MAP_DIR}"
    rm -f "${MAP_TAR}"
    echo "  ✓ SHAPEIT5 genetic maps"
else
    echo "  ✓ SHAPEIT5 genetic maps (already downloaded)"
fi

# 1d. dbSNP GRCh38 VCF for rsID annotation
# Pinned to dbSNP build 157 (the build the June 2026 run annotated with; identical to
# NCBI's latest_release since 15 Jan 2025), from NCBI's permanent archive so rsIDs do
# not change when NCBI publishes a new build. Checked against NCBI's md5.
DBSNP_BUILD="157"
DBSNP_FILENAME="GCF_000001405.40.gz"
DBSNP_VCF_DIR_URL="https://ftp.ncbi.nlm.nih.gov/snp/archive/b${DBSNP_BUILD}/VCF"

DBSNP_DEST="${REF_DIR}/dbsnp/${DBSNP_FILENAME}"
if [ ! -f "${DBSNP_DEST}" ]; then
    echo "Downloading dbSNP build ${DBSNP_BUILD} GRCh38 VCF (${DBSNP_FILENAME}, ~28 GB)..."
    curl -fsSL "${DBSNP_VCF_DIR_URL}/${DBSNP_FILENAME}" -o "${DBSNP_DEST}.part"
    want_md5="$(curl -fsSL "${DBSNP_VCF_DIR_URL}/${DBSNP_FILENAME}.md5" | cut -c1-32)"
    if [ "$(md5sum "${DBSNP_DEST}.part" | cut -d' ' -f1)" != "${want_md5}" ]; then
        echo "ERROR: dbSNP download failed its md5 check; delete ${DBSNP_DEST}.part and rerun." >&2
        exit 1
    fi
    mv "${DBSNP_DEST}.part" "${DBSNP_DEST}"
fi
if [ ! -f "${DBSNP_DEST}.tbi" ]; then
    echo "Downloading dbSNP GRCh38 VCF index..."
    fetch "${DBSNP_VCF_DIR_URL}/${DBSNP_FILENAME}.tbi" "${DBSNP_DEST}.tbi"
fi

echo "Reference data download process complete."
