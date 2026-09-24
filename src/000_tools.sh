#!/bin/bash
#SBATCH --job-name=000_tools
#SBATCH --output=src/logs/000_tools.out
#SBATCH --error=src/logs/000_tools.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --partition=low_p

# Script: src/000_tools.sh
# Purpose: Install the host-side tools the pipeline calls outside Nextflow's
# per-process conda environments, into ${TOOLS_DIR} (wrappers in ${TOOLS_DIR}/bin):
#   - htslib (bgzip, tabix) and bcftools, compiled from source
#   - plink 1.9, SHAPEIT5 and UCSC liftOver, as Apptainer images
#   - plink2 (pinned, in its own conda environment) and R (conda)
#   - triple-liftOver (stage-1 liftover to hg38), from GitHub at a pinned commit
# Submit from the repository root: sbatch src/000_tools.sh

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

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda is not on PATH. Load/initialise conda in the shell you submit from." >&2
  exit 1
fi
CONDA_BASE="$(conda info --base)"

IMG_DIR="${TOOLS_DIR}/singularity_images"
BIN_DIR="${TOOLS_DIR}/bin"
SRC_DIR="${TOOLS_DIR}/src"
RPATH_FLAGS="-Wl,-rpath,${TOOLS_DIR}/lib -Wl,-rpath,${TOOLS_DIR}/lib64"
RENV_NAME="renv"
RENV_PATH="${CONDA_BASE}/envs/${RENV_NAME}"
# plink2 gets its own pinned environment inside TOOLS_DIR: stage-1 and stage-3 QC need
# --king-cutoff-table (added 24 Jun 2024; bioconda 2.0.0a.6.9 is built 29 Jan 2025), and
# a shared or pre-existing env can leave an older plink2 in place.
PLINK2_VERSION="2.0.0a.6.9"
PLINK2_ENV="${TOOLS_DIR}/envs/plink2"

# True when the given plink2 recognises --king-cutoff-table.
plink2_has_king_cutoff_table() {
    local out
    out="$("$1" --king-cutoff-table 2>&1 || true)"
    [[ "$out" == *"--king-cutoff-table requires"* ]]
}
CONDA_BIN="${CONDA_BASE}/bin"
TRIPLE_LIFTOVER_DIR="${TRIPLE_LIFTOVER_DIR:-${TOOLS_DIR}/triple-liftOver}"

mkdir -p "$TOOLS_DIR" "$IMG_DIR" "$BIN_DIR" "$SRC_DIR"

export PATH="${BIN_DIR}:${PATH}"
export APPTAINER_BINDPATH

echo "=========================================="
echo " Setting up Genomic Tools Infrastructure"
echo " (Hybrid: Compiled + Apptainer Mode)"
echo " Host: $(hostname) | glibc: $(ldd --version | head -1)"
echo "=========================================="

# ── 1. Compiled Tools (htslib, bcftools) ──────────────────────────────────────
echo ""
echo "[1/6] Compiling core utilities from source..."

HTSLIB_URL="https://github.com/samtools/htslib/releases/download/1.23.1/htslib-1.23.1.tar.bz2"
BCFTOOLS_URL="https://github.com/samtools/bcftools/releases/download/1.23.1/bcftools-1.23.1.tar.bz2"

# 1a. HTSLIB
if [ ! -f "${BIN_DIR}/tabix" ]; then
    echo "  Compiling htslib..."
    cd "${SRC_DIR}"
    curl -fsSL "$HTSLIB_URL" -o htslib.tar.bz2
    tar -xjf htslib.tar.bz2
    cd htslib-1.23.1
    ./configure --prefix="${TOOLS_DIR}" LDFLAGS="${RPATH_FLAGS}"
    make -j4
    make install
    echo "  ✓ htslib"
else
    echo "  ✓ htslib (already installed)"
fi

# 1b. BCFTOOLS
if [ ! -f "${BIN_DIR}/bcftools" ]; then
    echo "  Compiling bcftools..."
    cd "${SRC_DIR}"
    curl -fsSL "$BCFTOOLS_URL" -o bcftools.tar.bz2
    tar -xjf bcftools.tar.bz2
    cd bcftools-1.23.1
    ./configure --prefix="${TOOLS_DIR}" --with-htslib="${TOOLS_DIR}" LDFLAGS="${RPATH_FLAGS}"
    make -j4
    make install
    echo "  ✓ bcftools"
else
    echo "  ✓ bcftools (already installed)"
fi

# PLINK2 comes from bioconda (its own pinned conda env, step 5);
# remove any previously downloaded invalid .sif file.
rm -f "${IMG_DIR}/plink2.sif"

# ── 2. Container Images ───────────────────────────────────────────────────────
echo ""
echo "[2/6] Pulling Apptainer images for complex tools..."

declare -A CONTAINERS=(
    ["plink"]="quay.io/biocontainers/plink:1.90b7.7--h18e278d_1"
    ["shapeit5"]="quay.io/biocontainers/shapeit5:5.1.1--h34261f4_2"
    ["liftover"]="quay.io/biocontainers/ucsc-liftover:469--h9b8f530_0"
)

for tool in "${!CONTAINERS[@]}"; do
    img_path="${IMG_DIR}/${tool}.sif"
    if [ ! -f "$img_path" ]; then
        echo "  Attempting to pull ${tool} from Quay.io..."
        if ! apptainer pull --name "$img_path" "docker://${CONTAINERS[$tool]}"; then
            echo "  ⚠ Pull failed. Attempting direct download from Galaxy Project Depot..."
            # Convert quay tag to galaxy URL format (usually tool:version--build)
            # Example: quay.io/biocontainers/plink:1.90b7.7--h18e278d_1 -> plink:1.90b7.7--h18e278d_1
            TAG_ONLY=$(echo "${CONTAINERS[$tool]}" | sed 's|.*/||')
            GALAXY_URL="https://depot.galaxyproject.org/singularity/${TAG_ONLY}"

            if ! curl -fsSL "$GALAXY_URL" -o "$img_path"; then
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

# ── 3. Wrapper Scripts ────────────────────────────────────────────────────────
echo ""
echo "[3/6] Creating wrappers..."

# plink
cat > "${BIN_DIR}/plink" << EOF
#!/bin/bash
apptainer exec "${IMG_DIR}/plink.sif" plink "\$@"
EOF

# plink2
cat > "${BIN_DIR}/plink2" << EOF
#!/bin/bash
exec "${PLINK2_ENV}/bin/plink2" "\$@"
EOF

# shapeit5
cat > "${BIN_DIR}/shapeit5" << EOF
#!/bin/bash
apptainer exec "${IMG_DIR}/shapeit5.sif" SHAPEIT5_phase_common "\$@"
EOF

# liftOver
cat > "${BIN_DIR}/liftOver" << EOF
#!/bin/bash
apptainer exec "${IMG_DIR}/liftover.sif" liftOver "\$@"
EOF

chmod +x "${BIN_DIR}/plink" "${BIN_DIR}/plink2" "${BIN_DIR}/shapeit5" "${BIN_DIR}/liftOver"

# ── 4. triple-liftOver ────────────────────────────────────────────────────────
# Stage 1 runs ${TRIPLE_LIFTOVER_DIR}/tripleliftover_v133.pl with the chain files in
# library/chainfiles/ (hg18ToHg38, hg19ToHg38); both ship with the repository.
echo ""
echo "[4/6] Installing triple-liftOver..."

TRIPLE_LIFTOVER_COMMIT="f781ed6c5d016eac75b8be92938003207cdf5f0d"
TRIPLE_LIFTOVER_URL="https://github.com/GraceSheng/triple-liftOver/archive/${TRIPLE_LIFTOVER_COMMIT}.tar.gz"

if [ -f "${TRIPLE_LIFTOVER_DIR}/tripleliftover_v133.pl" ]; then
    echo "  ✓ triple-liftOver (already installed)"
elif [ -e "${TRIPLE_LIFTOVER_DIR}" ]; then
    echo "  ✗ FAILED: ${TRIPLE_LIFTOVER_DIR} exists but has no tripleliftover_v133.pl;"
    echo "    remove it or point TRIPLE_LIFTOVER_DIR in .env at a valid installation."
else
    echo "  Downloading triple-liftOver (${TRIPLE_LIFTOVER_COMMIT:0:7})..."
    tmp_dir="${SRC_DIR}/triple-liftOver.tmp"
    rm -rf "$tmp_dir"
    mkdir -p "$tmp_dir" "$(dirname "${TRIPLE_LIFTOVER_DIR}")"
    curl -fsSL "$TRIPLE_LIFTOVER_URL" | tar -xz --strip-components=1 -C "$tmp_dir"
    mv "$tmp_dir" "${TRIPLE_LIFTOVER_DIR}"
    echo "  ✓ triple-liftOver"
fi

# triple-liftOver calls library/liftOver next to its script; point that at the
# containerised liftOver rather than the bundled binary.
if [ -d "${TRIPLE_LIFTOVER_DIR}/library" ]; then
    ln -sf "${BIN_DIR}/liftOver" "${TRIPLE_LIFTOVER_DIR}/library/liftOver"
fi

# ── 5. Conda environments (plink2, R) ─────────────────────────────────────────
echo ""
echo "[5/6] Installing plink2 and R packages..."
FAILURES=0

# 5a. plink2, pinned, in its own environment (recreated if missing or too old)
if [ -x "${PLINK2_ENV}/bin/plink2" ] && plink2_has_king_cutoff_table "${PLINK2_ENV}/bin/plink2"; then
    echo "  ✓ plink2 (already installed: $("${PLINK2_ENV}/bin/plink2" --version | head -1))"
else
    echo "  Creating plink2 ${PLINK2_VERSION} environment at ${PLINK2_ENV}..."
    mkdir -p "$(dirname "${PLINK2_ENV}")"
    "${CONDA_BIN}/conda" create -y -p "${PLINK2_ENV}" --override-channels -c conda-forge -c bioconda \
        "plink2=${PLINK2_VERSION}"
    if [ -x "${PLINK2_ENV}/bin/plink2" ]; then
        echo "  ✓ plink2 ($("${PLINK2_ENV}/bin/plink2" --version | head -1))"
    else
        echo "  ✗ FAILED: plink2 environment creation failed"
        FAILURES=$((FAILURES + 1))
    fi
fi

# 5b. R. miniconda's R Makeconf hardcodes /opt/rh/devtoolset-8 (CentOS compiler, absent on Ubuntu).
# install.packages() compilation always fails. Fix: dedicated conda env with pre-built
# conda-forge R packages — no compilation involved at all.

if [ ! -d "${RENV_PATH}" ]; then
    echo "  Creating conda R environment (r-base + r-haven + r-tidyverse)..."
    "${CONDA_BIN}/conda" create -y -n "${RENV_NAME}" --override-channels -c conda-forge \
        r-base r-haven r-tidyverse
    if [ ! -x "${RENV_PATH}/bin/Rscript" ]; then
        echo "  ✗ FAILED: conda renv creation failed"
        FAILURES=$((FAILURES + 1))
    else
        echo "  ✓ conda renv (r-base + r-haven + r-tidyverse)"
    fi
else
    if "${RENV_PATH}/bin/Rscript" -e "requireNamespace('haven', quietly=TRUE)" 2>/dev/null; then
        echo "  ✓ conda renv (already set up, haven available)"
    else
        echo "  Adding r-haven to existing renv..."
        "${CONDA_BIN}/conda" install -y -n "${RENV_NAME}" --override-channels -c conda-forge r-haven r-tidyverse
        if ! "${RENV_PATH}/bin/Rscript" -e "requireNamespace('haven', quietly=TRUE)" 2>/dev/null; then
            echo "  ✗ FAILED: r-haven not available in renv after install"
            FAILURES=$((FAILURES + 1))
        else
            echo "  ✓ r-haven added to renv"
        fi
    fi
fi

# Wrapper so scripts using Rscript from BIN_DIR pick up the renv automatically
cat > "${BIN_DIR}/Rscript" << EOF
#!/bin/bash
exec "${RENV_PATH}/bin/Rscript" "\$@"
EOF
chmod +x "${BIN_DIR}/Rscript"

# ── 6. Verification ────────────────────────────────────────────────────────────
echo ""
echo "[6/6] Verifying tools in ${BIN_DIR}..."

# List of tools to check specifically in bin/
CHECK_TOOLS=(
    "bcftools"
    "bgzip"
    "tabix"
    "plink"
    "plink2"
    "shapeit5"
    "liftOver"
    "Rscript"
)

for tool in "${CHECK_TOOLS[@]}"; do
    if [ -x "${BIN_DIR}/$tool" ]; then
        echo "  ✓ $tool"
    else
        echo "  ✗ FAILED: $tool is missing or not executable in ${BIN_DIR}"
        FAILURES=$((FAILURES + 1))
    fi
done

# plink2 must support --king-cutoff-table (stage-1 duplicate removal, stage-3 relatedness)
if [ -x "${BIN_DIR}/plink2" ] && plink2_has_king_cutoff_table "${BIN_DIR}/plink2"; then
    echo "  ✓ plink2 supports --king-cutoff-table"
else
    echo "  ✗ FAILED: ${BIN_DIR}/plink2 does not support --king-cutoff-table"
    FAILURES=$((FAILURES + 1))
fi

# triple-liftOver: script, both chain files, and the liftOver link
if [ -f "${TRIPLE_LIFTOVER_DIR}/tripleliftover_v133.pl" ] \
    && [ -f "${TRIPLE_LIFTOVER_DIR}/library/chainfiles/hg18ToHg38.over.chain.gz" ] \
    && [ -f "${TRIPLE_LIFTOVER_DIR}/library/chainfiles/hg19ToHg38.over.chain.gz" ] \
    && [ -L "${TRIPLE_LIFTOVER_DIR}/library/liftOver" ]; then
    echo "  ✓ triple-liftOver (${TRIPLE_LIFTOVER_DIR})"
else
    echo "  ✗ FAILED: triple-liftOver not configured correctly in ${TRIPLE_LIFTOVER_DIR}"
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
