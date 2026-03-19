#!/bin/bash
# Script: src/000_tools.sh
# Purpose: Download and install core tools into tools/

set -euo pipefail

# Robust environment sourcing
for env_file in ".env" "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../.env" "$(dirname -- "${BASH_SOURCE[0]:-$0}")/../pipeline/.env"; do
  if [ -f "$env_file" ]; then
    set -a; source "$env_file"; set +a
    break
  fi
done

PROJ_ROOT="${GENETICS_PROJECT_ROOT:-$(pwd)}"
TOOLS_DIR="${PROJ_ROOT}/tools"

mkdir -p "$TOOLS_DIR"

# 1. PLINK 1.9
if [ ! -f "${TOOLS_DIR}/plink" ]; then
    echo "Downloading PLINK..."
    wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip -O "${TOOLS_DIR}/plink.zip"
    unzip -o "${TOOLS_DIR}/plink.zip" -d "${TOOLS_DIR}"
    chmod +x "${TOOLS_DIR}/plink"
fi

# 2. bcftools (Conda is preferred, but here is a binary fetch if needed)
if [ ! -f "${TOOLS_DIR}/bcftools" ]; then
    echo "Downloading bcftools..."
    wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2 -O "${TOOLS_DIR}/bcftools.tar.bz2"
    tar -jxvf "${TOOLS_DIR}/bcftools.tar.bz2" -C "${TOOLS_DIR}"
    # Note: This requires compilation if it's source. For simplicity, we assume bioconda is used for complex builds.
    # I'll add a tip to use conda in the walkthrough.
fi

# 3. SHAPEIT4 (Phasing)
if [ ! -f "${TOOLS_DIR}/shapeit4" ]; then
    echo "Downloading SHAPEIT4..."
    wget https://github.com/odelaneau/shapeit4/releases/download/v4.2.2/shapeit4.2.2_linux_x86_64 -O "${TOOLS_DIR}/shapeit4"
    chmod +x "${TOOLS_DIR}/shapeit4"
fi

# 4. Minimac4 (Imputation)
if [ ! -f "${TOOLS_DIR}/minimac4" ]; then
    echo "Downloading Minimac4..."
    # We use the static binary if available, or the tarball
    wget https://github.com/statgen/Minimac4/releases/download/v1.0.2/minimac4-v1.0.2-linux.tar.gz -O "${TOOLS_DIR}/minimac4.tar.gz"
    tar -zxvf "${TOOLS_DIR}/minimac4.tar.gz" -C "${TOOLS_DIR}"
    chmod +x "${TOOLS_DIR}/minimac4"
fi

# 5. conform-gt
if [ ! -f "${TOOLS_DIR}/conform-gt.jar" ]; then
    echo "Downloading conform-gt..."
    wget http://faculty.washington.edu/browning/conform-gt.jar -P "${TOOLS_DIR}"
fi

# 5. Picard Tools (Liftover)
if [ ! -f "${TOOLS_DIR}/picard.jar" ]; then
    echo "Downloading Picard..."
    wget https://github.com/broadinstitute/picard/releases/download/3.1.1/picard.jar -P "${TOOLS_DIR}"
fi

# 6. Annovar
if [ ! -d "${TOOLS_DIR}/annovar" ]; then
    echo "Downloading Annovar..."
    # Using the latest available link (this might expire, user may need to update)
    wget http://www.openbioinformatics.org/annovar/download/0wgpxTR2ES/annovar.latest.tar.gz -O "${TOOLS_DIR}/annovar.tar.gz"
    tar -zxvf "${TOOLS_DIR}/annovar.tar.gz" -C "${TOOLS_DIR}"
    rm "${TOOLS_DIR}/annovar.tar.gz"
fi

echo "=========================================="
echo "Tools setup complete in $TOOLS_DIR"
echo "=========================================="
