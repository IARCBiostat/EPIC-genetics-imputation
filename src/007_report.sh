#!/bin/bash
#SBATCH --job-name=epic_report
#SBATCH --output=src/logs/007_report.out
#SBATCH --error=src/logs/007_report.err
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=low_p

set -euo pipefail
trap 'echo "ERROR: Report generation failed on line $LINENO" >&2; exit 1' ERR

SCRIPT_DIR="$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")" && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/$(basename -- "${BASH_SOURCE[0]:-$0}")"
DEFAULT_PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

is_repo_root() {
  local candidate="$1"
  [ -d "${candidate}/src" ] && \
  [ -d "${candidate}/analysis" ] && \
  [ -f "${candidate}/src/007_report/master_report.py" ]
}

source_root_env() {
  local root="$1"
  local env_file="${root}/.env"

  if [ ! -f "$env_file" ]; then
    echo "ERROR: Root environment file not found: ${env_file}" >&2
    echo "       Create ${env_file}; stage-specific .env files are no longer supported." >&2
    return 1
  fi

  set -a
  # shellcheck disable=SC1090
  source "$env_file"
  set +a
}

resolve_project_root() {
  local candidates=()
  local candidate
  local resolved

  [ -n "${SLURM_SUBMIT_DIR:-}" ] && candidates+=("${SLURM_SUBMIT_DIR}")
  candidates+=("${SCRIPT_DIR}/.." "$(pwd)")
  [ -n "${GENETICS_PROJECT_ROOT:-}" ] && candidates+=("${GENETICS_PROJECT_ROOT}")
  candidates+=("${DEFAULT_PROJ_ROOT}")

  for candidate in "${candidates[@]}"; do
    [ -z "$candidate" ] && continue
    if resolved="$(cd "$candidate" 2>/dev/null && pwd)"; then
      if is_repo_root "$resolved"; then
        printf '%s\n' "$resolved"
        return 0
      fi
    fi
  done

  echo "ERROR: Could not resolve the EPIC genetics repo root." >&2
  echo "       Tried SLURM_SUBMIT_DIR, script parent, current directory, and GENETICS_PROJECT_ROOT." >&2
  return 1
}

print_help() {
  cat <<'EOF'
Usage: bash src/007_report.sh [options]

Generates cross-stage per-study master reports under analysis/report/.
When run outside Slurm, submits itself as a batch job automatically.

Options:
  --study, --studies <list>    Comma-separated study IDs (default: all)
  --analysis-root <dir>        Analysis root directory (default: <project>/analysis)
  --partition <name>           Slurm partition (default: low_p)
  --time <duration>            Wall time for the job (default: 4:00:00)
  --mem <size>                 Memory for the job (default: 16G)
  --cpus <n>                   CPUs for the job (default: 2)
  -h, --help                   Show this help
EOF
}

PROJ_ROOT="$(resolve_project_root)"
export GENETICS_PROJECT_ROOT="${PROJ_ROOT}"
source_root_env "$PROJ_ROOT"
export GENETICS_PROJECT_ROOT="${PROJ_ROOT}"

if [ -z "${SLURM_JOB_ID:-}" ]; then
  mkdir -p "${PROJ_ROOT}/src/logs"
  cd "$PROJ_ROOT"

  if ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch is not available. Submit this script from the HPC login node." >&2
    exit 1
  fi

  echo "Submitting report generation to Slurm..."
  sbatch \
    --export=ALL \
    --job-name "${REPORT_JOB_NAME:-epic_report}" \
    --output "${REPORT_LOG_OUT:-${PROJ_ROOT}/src/logs/007_report_%j.out}" \
    --error "${REPORT_LOG_ERR:-${PROJ_ROOT}/src/logs/007_report_%j.err}" \
    --time "${REPORT_TIME:-4:00:00}" \
    --mem "${REPORT_MEM:-16G}" \
    --cpus-per-task "${REPORT_CPUS:-2}" \
    --partition "${REPORT_PARTITION:-low_p}" \
    "${SCRIPT_PATH}" "$@"
  exit 0
fi

start_time=$(date +%s)

ANALYSIS_ROOT="${REPORT_ANALYSIS_ROOT:-${PROJ_ROOT}/analysis}"
STUDIES="${REPORT_STUDIES:-all}"
PYTHON3_BIN="${PYTHON3_BIN:-python3}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --study|--studies)
      STUDIES="$2"
      shift 2
      ;;
    --analysis-root)
      ANALYSIS_ROOT="$2"
      shift 2
      ;;
    --partition)
      shift 2
      ;;
    --time)
      shift 2
      ;;
    --mem)
      shift 2
      ;;
    --cpus)
      shift 2
      ;;
    -h|--help)
      print_help
      exit 0
      ;;
    *)
      echo "ERROR: Unknown argument: $1" >&2
      print_help >&2
      exit 1
      ;;
  esac
done

if [ -n "${GENETICS_TOOLS_BIN:-}" ]; then
  export PATH="${GENETICS_TOOLS_BIN}:${PATH}"
elif [ -d "${PROJ_ROOT}/tools/bin" ]; then
  export PATH="${PROJ_ROOT}/tools/bin:${PATH}"
fi

if command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base)"
  # shellcheck disable=SC1091
  source "${CONDA_BASE}/etc/profile.d/conda.sh" || true
  conda activate nf_EPIC-genetics || true
fi

if [ ! -d "${ANALYSIS_ROOT}" ]; then
  echo "ERROR: Analysis root does not exist: ${ANALYSIS_ROOT}" >&2
  exit 1
fi

if [ ! -f "${PROJ_ROOT}/src/007_report/master_report.py" ]; then
  echo "ERROR: Missing report builder: ${PROJ_ROOT}/src/007_report/master_report.py" >&2
  exit 1
fi

echo "============================================"
echo " Running EPIC Genetics Master Report"
echo "============================================"
echo "Slurm job:     ${SLURM_JOB_ID:-interactive}"
echo "Host:          $(hostname)"
echo "Project root:  ${PROJ_ROOT}"
echo "Analysis root: ${ANALYSIS_ROOT}"
echo "Studies:       ${STUDIES}"
echo "Python3:       ${PYTHON3_BIN}"
echo "============================================"
echo ""

"${PYTHON3_BIN}" "${PROJ_ROOT}/src/007_report/master_report.py" \
  --analysis-root "${ANALYSIS_ROOT}" \
  --studies "${STUDIES}"

echo ""
echo "============================================"
echo " Master report generation completed"
echo " Studies: ${STUDIES}"
echo "============================================"

end_time=$(date +%s)
elapsed=$((end_time - start_time))
echo "Total runtime: $(date -u -r "$elapsed" +%H:%M:%S 2>/dev/null || date -u -d "@$elapsed" +%H:%M:%S)"
