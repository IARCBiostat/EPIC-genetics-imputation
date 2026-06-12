#!/bin/bash
#SBATCH --job-name=documentation
#SBATCH --output=src/logs/documentation.out
#SBATCH --error=src/logs/documentation.err
#SBATCH --time=1:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --partition=low_p

set -euo pipefail
trap 'echo "ERROR: Documentation script failed on line $LINENO" >&2; exit 1' ERR

SCRIPT_PATH="$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")" && pwd)/$(basename -- "${BASH_SOURCE[0]:-$0}")"

print_help() {
  cat <<'EOF'
Usage: bash src/008_documentation.sh [options]

Copies finalised outputs from scratch to the local repository, regenerates
the sample overlap UpSet plot, and updates the README.md summary table.

Must be run on the HPC login node (or as an sbatch job) where the scratch
filesystem is accessible.

Options:
  --analysis-root <d>  Studies root on scratch (default: from .env)
  --dest-root <d>      Finalised outputs root on scratch (default: from .env)
  --partition <name>   Slurm partition (default: low_p)
  -h, --help           Show this help
EOF
}

ENV_FILE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")/.." && pwd)}/.env"
if [ ! -f "$ENV_FILE" ]; then
  echo "ERROR: .env not found at ${ENV_FILE}" >&2; exit 1
fi
set -a; source "$ENV_FILE"; set +a
PROJ_ROOT="${GENETICS_PROJECT_ROOT}"

# ── Self-submission ────────────────────────────────────────────────────────────
if [ -z "${SLURM_JOB_ID:-}" ]; then
  mkdir -p "${PROJ_ROOT}/src/logs"
  cd "$PROJ_ROOT"
  if ! command -v sbatch >/dev/null 2>&1; then
    # No Slurm: run directly (e.g. local testing)
    exec bash "$SCRIPT_PATH" "$@"
  fi
  echo "Submitting documentation job to Slurm..."
  sbatch \
    --export=ALL \
    --job-name  "epic_documentation" \
    --output    "${PROJ_ROOT}/src/logs/008_documentation_%j.out" \
    --error     "${PROJ_ROOT}/src/logs/008_documentation_%j.err" \
    --time      "1:00:00" \
    --mem       "8G" \
    --cpus-per-task "2" \
    --partition "${SLURM_JOB_PARTITION:-low_p}" \
    "${SCRIPT_PATH}" "$@"
  exit 0
fi

# ── Inside Slurm job ───────────────────────────────────────────────────────────
start_time=$(date +%s)

ANALYSIS_ROOT="${DOC_ANALYSIS_ROOT:-${SCRATCH_RUN}/studies}"
DEST_ROOT="${DOC_DEST_ROOT:-${SCRATCH_RUN}/final}"
CONDA_ENV="${DOC_CONDA_ENV:-epic_documentation}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --analysis-root)  ANALYSIS_ROOT="$2";  shift 2 ;;
    --dest-root)      DEST_ROOT="$2";      shift 2 ;;
    --partition)      shift 2 ;;
    -h|--help)        print_help; exit 0           ;;
    *) echo "ERROR: Unknown argument: $1" >&2; print_help >&2; exit 1 ;;
  esac
done

LOCAL_REPORT_DIR="${PROJ_ROOT}/analysis/final/report"
LOCAL_SUMMARIES_DIR="${PROJ_ROOT}/analysis/final/summaries"
SCRATCH_REPORT_DIR="${DEST_ROOT}/report"
SCRATCH_SUMMARIES_DIR="${DEST_ROOT}/summaries"
DOCS_REPORTS_DIR="${PROJ_ROOT}/src/docs/data/reports"
DOCS_IMG_DIR="${PROJ_ROOT}/src/docs/img"

echo "=========================================="
echo " EPIC Genetics — Documentation"
echo "=========================================="
echo "Slurm job:    ${SLURM_JOB_ID:-interactive}"
echo "Host:         $(hostname)"
echo "Project root: ${PROJ_ROOT}"
echo "Analysis root:${ANALYSIS_ROOT}"
echo "Dest root:    ${DEST_ROOT}"
echo "=========================================="

# ── Activate conda env ────────────────────────────────────────────────────────
source "$(conda info --base)/etc/profile.d/conda.sh" || true

if conda env list | grep -q "^${CONDA_ENV}"; then
  echo "Updating conda environment '${CONDA_ENV}' from envs/documentation.yml ..."
  conda env update -n "${CONDA_ENV}" -f "${PROJ_ROOT}/envs/documentation.yml" --prune
else
  echo "Creating conda environment '${CONDA_ENV}' from envs/documentation.yml ..."
  conda env create -n "${CONDA_ENV}" -f "${PROJ_ROOT}/envs/documentation.yml"
fi

conda activate "${CONDA_ENV}"

PYTHON_BIN="$(which python3)"
echo "Python: ${PYTHON_BIN}"

# ── Copy analysis/final/ from scratch (report/ + summaries/) ─────────────────
echo ""
echo "Copying master report HTMLs: ${SCRATCH_REPORT_DIR}/ -> ${LOCAL_REPORT_DIR}/"
mkdir -p "${LOCAL_REPORT_DIR}"
if compgen -G "${SCRATCH_REPORT_DIR}/*.master-report.html" > /dev/null; then
  cp -v "${SCRATCH_REPORT_DIR}"/*.master-report.html "${LOCAL_REPORT_DIR}/"
else
  echo "WARNING: no master report HTMLs found at ${SCRATCH_REPORT_DIR}" >&2
fi

# ── Copy summaries from scratch to local ──────────────────────────────────────
echo ""
echo "Copying summaries: ${SCRATCH_SUMMARIES_DIR}/*.md -> ${LOCAL_SUMMARIES_DIR}/"
mkdir -p "${LOCAL_SUMMARIES_DIR}"
if compgen -G "${SCRATCH_SUMMARIES_DIR}/*.md" > /dev/null; then
  cp -v "${SCRATCH_SUMMARIES_DIR}"/*.md "${LOCAL_SUMMARIES_DIR}/"
else
  echo "WARNING: no summary markdown files found at ${SCRATCH_SUMMARIES_DIR}" >&2
fi

# ── Sample overlap ─────────────────────────────────────────────────────────────
echo ""
echo "Generating sample overlap UpSet plot ..."
"${PYTHON_BIN}" "${PROJ_ROOT}/src/misc/sample_overlap.py" \
  --analysis-root "${ANALYSIS_ROOT}" \
  --out-dir       "${LOCAL_REPORT_DIR}"

# ── Sync analysis/final/report/ → docs/ (for GitHub Pages) ───────────────────
echo ""
echo "Syncing reports and images to docs/ ..."
mkdir -p "${DOCS_REPORTS_DIR}" "${DOCS_IMG_DIR}"
if compgen -G "${LOCAL_REPORT_DIR}/*.master-report.html" > /dev/null; then
  cp -v "${LOCAL_REPORT_DIR}"/*.master-report.html "${DOCS_REPORTS_DIR}/"
fi
if [ -f "${LOCAL_REPORT_DIR}/sample_overlap_upset.png" ]; then
  cp -v "${LOCAL_REPORT_DIR}/sample_overlap_upset.png" "${DOCS_IMG_DIR}/"
fi

# ── Update summary tables (README.md, docs home page, docs data page) ─────────
echo ""
echo "Updating README.md summary table ..."
"${PYTHON_BIN}" "${PROJ_ROOT}/src/misc/update_summary_table.py" \
  --analysis-root "${ANALYSIS_ROOT}" \
  --summaries-dir "${LOCAL_SUMMARIES_DIR}" \
  --overlap-dir   "${LOCAL_REPORT_DIR}" \
  --readme        "${PROJ_ROOT}/README.md"

echo ""
echo "Updating docs home page and data page summary tables ..."
"${PYTHON_BIN}" "${PROJ_ROOT}/src/misc/update_summary_table.py" \
  --analysis-root "${ANALYSIS_ROOT}" \
  --summaries-dir "${LOCAL_SUMMARIES_DIR}" \
  --overlap-dir   "${LOCAL_REPORT_DIR}" \
  --readme        "${PROJ_ROOT}/src/docs/index.md" \
  --data-page     "${PROJ_ROOT}/src/docs/data/index.md"

# ── Build static site into docs/ (GitHub Pages source) ───────────────────────
echo ""
echo "Building MkDocs site into docs/ ..."
cd "${PROJ_ROOT}"
mkdocs build --config-file "${PROJ_ROOT}/envs/mkdocs.yml" --clean

# ── Done ──────────────────────────────────────────────────────────────────────
end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
echo ""
echo "=========================================="
echo " Documentation complete (${elapsed}s)"
echo "=========================================="
echo ""
echo "Files updated in ${PROJ_ROOT}:"
echo "  analysis/final/report/*.master-report.html"
echo "  analysis/final/report/sample_overlap_upset.png"
echo "  analysis/final/report/sample_overlap_summary.json"
echo "  analysis/final/summaries/*.md"
echo "  src/docs/data/reports/*.master-report.html"
echo "  src/docs/img/sample_overlap_upset.png"
echo "  docs/  (built site — GitHub Pages source)"
echo "  README.md"
echo "  src/docs/index.md"
echo ""
echo "Review changes with: git diff README.md src/docs/index.md"
echo "Then commit with:    git add docs/ src/docs/ README.md && git commit && git push"
