#!/bin/bash
#SBATCH --job-name=000_tools
#SBATCH --output=logs/000_tools.out
#SBATCH --error=logs/000_tools.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --partition=low_p

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

PROJ_ROOT="${GENETICS_PROJECT_ROOT:-$(pwd)}"
TOOLS_DIR="${PROJ_ROOT}/tools"
IMG_DIR="${TOOLS_DIR}/singularity_images"
BIN_DIR="${TOOLS_DIR}/bin"
SRC_DIR="${TOOLS_DIR}/src"
RPATH_FLAGS="-Wl,-rpath,${TOOLS_DIR}/lib -Wl,-rpath,${TOOLS_DIR}/lib64"

mkdir -p "$TOOLS_DIR" "$IMG_DIR" "$BIN_DIR" "$SRC_DIR" logs

export PATH="${BIN_DIR}:${PATH}"
export APPTAINER_BINDPATH="/data:/data"

echo "=========================================="
echo " Setting up Genomic Tools Infrastructure"
echo " (Hybrid: Compiled + Apptainer Mode)"
echo " Host: $(hostname) | glibc: $(ldd --version | head -1)"
echo "=========================================="

# ── 1. Compiled Tools (htslib, samtools, bcftools) ────────────────────────────
echo ""
echo "[1/4] Compiling core utilities from source..."

# URLs provided by user
HTSLIB_URL="https://github.com/samtools/htslib/releases/download/1.23.1/htslib-1.23.1.tar.bz2"
SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/1.23.1/samtools-1.23.1.tar.bz2"
BCFTOOLS_URL="https://github.com/samtools/bcftools/releases/download/1.23.1/bcftools-1.23.1.tar.bz2"
MINIMAC4_URL="https://github.com/statgen/Minimac4/releases/latest/download/Minimac4-4.1.6-Linux-static.tgz"

# 1a. HTSLIB
if [ ! -f "${BIN_DIR}/tabix" ]; then
    echo "  Compiling htslib..."
    cd "${SRC_DIR}"
    curl -sL "$HTSLIB_URL" -o htslib.tar.bz2
    tar -xjf htslib.tar.bz2
    cd htslib-1.23.1
    ./configure --prefix="${TOOLS_DIR}" LDFLAGS="${RPATH_FLAGS}"
    make -j4
    make install
    echo "  ✓ htslib"
else
    echo "  ✓ htslib (already installed)"
fi

# 1b. SAMTOOLS
if [ ! -f "${BIN_DIR}/samtools" ]; then
    echo "  Compiling samtools..."
    cd "${SRC_DIR}"
    curl -sL "$SAMTOOLS_URL" -o samtools.tar.bz2
    tar -xjf samtools.tar.bz2
    cd samtools-1.23.1
    ./configure --prefix="${TOOLS_DIR}" --with-htslib="${TOOLS_DIR}" LDFLAGS="${RPATH_FLAGS}"
    make -j4
    make install
    echo "  ✓ samtools"
else
    echo "  ✓ samtools (already installed)"
fi

# 1c. BCFTOOLS
if [ ! -f "${BIN_DIR}/bcftools" ]; then
    echo "  Compiling bcftools..."
    cd "${SRC_DIR}"
    curl -sL "$BCFTOOLS_URL" -o bcftools.tar.bz2
    tar -xjf bcftools.tar.bz2
    cd bcftools-1.23.1
    ./configure --prefix="${TOOLS_DIR}" --with-htslib="${TOOLS_DIR}" LDFLAGS="${RPATH_FLAGS}"
    make -j4
    make install
    echo "  ✓ bcftools"
else
    echo "  ✓ bcftools (already installed)"
fi

# 1d. Minimac4 static binary
if [ ! -x "${BIN_DIR}/minimac4.native" ]; then
    echo "  Installing Minimac4 static binary..."
    cd "${SRC_DIR}"
    rm -rf minimac4-static
    mkdir -p minimac4-static
    curl -sL "$MINIMAC4_URL" -o minimac4-static.tgz
    tar -xzf minimac4-static.tgz -C minimac4-static
    MINIMAC4_NATIVE_PATH="$(find minimac4-static -type f -name minimac4 | head -n 1 || true)"
    if [ -z "${MINIMAC4_NATIVE_PATH}" ]; then
        echo "  ✗ Could not locate minimac4 in extracted archive" >&2
        exit 1
    fi
    cp "${MINIMAC4_NATIVE_PATH}" "${BIN_DIR}/minimac4.native"
    chmod +x "${BIN_DIR}/minimac4.native"
    echo "  ✓ minimac4 static binary"
else
    echo "  ✓ minimac4 static binary (already installed)"
fi

cd "${PROJ_ROOT}"

# ── 2. Container Images ───────────────────────────────────────────────────────
echo ""
echo "[2/4] Pulling Apptainer images for complex tools..."

declare -A CONTAINERS=(
    ["plink"]="quay.io/biocontainers/plink:1.90b7.7--h18e278d_1"
    ["eagle"]="quay.io/biocontainers/eagle:0.9.0--py34_0"
    ["shapeit5"]="quay.io/biocontainers/shapeit5:5.1.1--h34261f4_2"
    ["picard"]="quay.io/biocontainers/picard:3.4.0--hdfd78af_0"
    ["liftover"]="quay.io/biocontainers/ucsc-liftover:469--h9b8f530_0"
)

for tool in "${!CONTAINERS[@]}"; do
    img_path="${IMG_DIR}/${tool}.sif"
    if [ ! -f "$img_path" ]; then
        echo "  Attempting to pull ${tool} from Quay.io..."
        if ! apptainer pull --name "$img_path" "docker://${CONTAINERS[$tool]}"; then
            echo "  ⚠ Pull failed. Attempting direct download from Galaxy Project Depot..."
            # Convert quay tag to galaxy URL format (usually tool:version--build)
            # Example: quay.io/biocontainers/eagle:2.4.1--h9ee0642_1 -> eagle:2.4.1--h9ee0642_1
            TAG_ONLY=$(echo "${CONTAINERS[$tool]}" | sed 's|.*/||')
            GALAXY_URL="https://depot.galaxyproject.org/singularity/${TAG_ONLY}"
            
            if ! curl -sL "$GALAXY_URL" -o "$img_path"; then
                echo "  ✗ FAILED: Could not pull or download ${tool}"
                rm -f "$img_path"
            else
                echo "  ✓ ${tool} (downloaded from Galaxy)"
            fi
        else
            echo "  ✓ ${tool} (pulled from Quay)"
        fi
    else
        echo "  ✓ ${tool} (cached)"
    fi
done

# ── 3. Wrapper Scripts & Manual Downloads ─────────────────────────────────────
echo ""
echo "[3/4] Creating wrappers and downloading manual tools..."

# Apptainer Wrappers (Only for tools not compiled in step 1)

# plink
cat > "${BIN_DIR}/plink" << EOF
#!/bin/bash
apptainer exec "${IMG_DIR}/plink.sif" plink "\$@"
EOF

# eagle
cat > "${BIN_DIR}/eagle" << EOF
#!/bin/bash
apptainer exec "${IMG_DIR}/eagle.sif" Eagle "\$@"
EOF

# minimac4
cat > "${BIN_DIR}/minimac4" << EOF
#!/bin/bash
exec "${BIN_DIR}/minimac4.native" "\$@"
EOF

# shapeit5
cat > "${BIN_DIR}/shapeit5" << EOF
#!/bin/bash
apptainer exec "${IMG_DIR}/shapeit5.sif" phase_common "\$@"
EOF

# picard
cat > "${BIN_DIR}/picard" << EOF
#!/bin/bash
apptainer exec "${IMG_DIR}/picard.sif" picard "\$@"
EOF

# liftOver
cat > "${BIN_DIR}/liftOver" << EOF
#!/bin/bash
apptainer exec "${IMG_DIR}/liftover.sif" liftOver "\$@"
EOF

chmod +x "${BIN_DIR}/plink" "${BIN_DIR}/eagle" "${BIN_DIR}/minimac4" "${BIN_DIR}/shapeit5" "${BIN_DIR}/picard" "${BIN_DIR}/liftOver"

# manual tools logic... (ANNOVAR, conform-gt)

# ANNOVAR
if [ ! -d "${TOOLS_DIR}/annovar" ]; then
    echo "  Downloading ANNOVAR..."
    ANNOVAR_URL="http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz"
    curl -sL "$ANNOVAR_URL" -o "${TOOLS_DIR}/annovar.tar.gz"
    tar -zxf "${TOOLS_DIR}/annovar.tar.gz" -C "${TOOLS_DIR}"
    rm "${TOOLS_DIR}/annovar.tar.gz"
    
    # Create wrappers for main annovar scripts
    cat > "${BIN_DIR}/annotate_variation.pl" << EOF
#!/bin/bash
perl "${TOOLS_DIR}/annovar/annotate_variation.pl" "\$@"
EOF
    cat > "${BIN_DIR}/table_annovar.pl" << EOF
#!/bin/bash
perl "${TOOLS_DIR}/annovar/table_annovar.pl" "\$@"
EOF
    chmod +x "${BIN_DIR}/annotate_variation.pl" "${BIN_DIR}/table_annovar.pl"
    echo "  ✓ ANNOVAR"
else
    echo "  ✓ ANNOVAR (already installed)"
fi

# conform-gt
if [ ! -f "${TOOLS_DIR}/conform-gt.jar" ]; then
    echo "  Downloading conform-gt..."
    curl -sL "https://faculty.washington.edu/browning/conform-gt/conform-gt.24May16.cee.jar" \
         -o "${TOOLS_DIR}/conform-gt.jar"
fi

cat > "${BIN_DIR}/conform-gt" << EOF
#!/bin/bash
java -jar "${TOOLS_DIR}/conform-gt.jar" "\$@"
EOF
chmod +x "${BIN_DIR}/conform-gt"

# triple-liftOver link
if [ -d "${TOOLS_DIR}/triple-liftOver" ]; then
    mkdir -p "${TOOLS_DIR}/triple-liftOver/library"
    ln -sf "${BIN_DIR}/liftOver" "${TOOLS_DIR}/triple-liftOver/library/liftOver"
fi

# ── 4. Verification ────────────────────────────────────────────────────────────
echo ""
echo "[4/4] Verifying tools in ${BIN_DIR}..."
FAILURES=0

# List of tools to check specifically in bin/
CHECK_TOOLS=(
    "bcftools"
    "samtools"
    "bgzip"
    "tabix"
    "plink"
    "eagle"
    "shapeit5"
    "picard"
    "minimac4"
    "table_annovar.pl"
    "conform-gt"
    "liftOver"
)

for tool in "${CHECK_TOOLS[@]}"; do
    if [ -x "${BIN_DIR}/$tool" ]; then
        echo "  ✓ $tool"
    else
        echo "  ✗ FAILED: $tool is missing or not executable in ${BIN_DIR}"
        FAILURES=$((FAILURES + 1))
    fi
done

# Special check for triple-liftOver
if [ -d "${TOOLS_DIR}/triple-liftOver" ] && [ -L "${TOOLS_DIR}/triple-liftOver/library/liftOver" ]; then
    echo "  ✓ triple-liftOver (linked)"
else
    echo "  ✗ FAILED: triple-liftOver not configured correctly"
    FAILURES=$((FAILURES + 1))
fi

# ── Summary ───────────────────────────────────────────────────────────────────
end_time=$(date +%s)
elapsed=$((end_time - start_time))
echo ""
echo "=========================================="
if [ "$FAILURES" -eq 0 ]; then
    echo " ALL TOOLS VERIFIED SUCCESSFULLY"
else
    echo " WARNING: $FAILURES tool(s) failed verification. Check the list above."
fi
echo "=========================================="
echo " Tools directory : ${TOOLS_DIR}"
echo " Wrappers        : ${BIN_DIR}"
echo " Images          : ${IMG_DIR}"
echo " Add to PATH     : export PATH=${BIN_DIR}:\$PATH"
echo " Time taken      : $((elapsed / 60))m $((elapsed % 60))s"
echo "=========================================="
