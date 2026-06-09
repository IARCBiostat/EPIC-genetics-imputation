#!/bin/bash
#SBATCH --job-name=stage4
#SBATCH --output=src/logs/stage4.out
#SBATCH --error=src/logs/stage4.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --partition=low_p

set -euo pipefail
trap 'echo "ERROR: Stage-4 pipeline failed on line $LINENO" >&2; exit 1' ERR

SCRIPT_PATH="$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")" && pwd)/$(basename -- "${BASH_SOURCE[0]:-$0}")"

abs_path() {
  local target="$1"
  local parent base
  if [[ "$target" = /* ]]; then printf '%s\n' "$target"; return 0; fi
  parent="$(dirname "$target")"
  base="$(basename "$target")"
  [ "$parent" = "." ] && parent="$(pwd)" || parent="$(cd "$parent" 2>/dev/null && pwd)"
  printf '%s/%s\n' "$parent" "$base"
}

print_help() {
  cat <<'EOF'
Usage: sbatch src/007_stage4.sh [options]

Generates per-study master HTML reports and builds the final deliverable
archive for each study under ${STAGE4_DEST_ROOT}/.

Options:
  --profile <name>     Nextflow profile (default: slurm)
  --study <list>       Comma-separated study IDs (default: all)
  --analysis-root <d>  Root containing study outputs (default: from .env)
  --stage2-root <d>    Stage 2 root (default: analysis-root)
  --stage3-root <d>    Stage 3 root (default: analysis-root)
  --dest-root <d>      Final output destination (default: from .env)
  --report-dir <d>     Directory for master HTML copies (default: dest-root/report)
  --partition <name>   Slurm partition (default: low_p)
  --no-resume          Disable Nextflow -resume
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
    echo "ERROR: sbatch is not available. Submit this script from the HPC login node." >&2; exit 1
  fi
  echo "Submitting stage-4 pipeline to Slurm..."
  sbatch \
    --export=ALL \
    --job-name  "${STAGE4_JOB_NAME:-stage4_finalise}" \
    --output    "${STAGE4_LOG_OUT:-${PROJ_ROOT}/src/logs/007_stage4_%j.out}" \
    --error     "${STAGE4_LOG_ERR:-${PROJ_ROOT}/src/logs/007_stage4_%j.err}" \
    --time      "${STAGE4_TIME:-4:00:00}" \
    --mem       "${STAGE4_MEM:-16G}" \
    --cpus-per-task "${STAGE4_CPUS:-2}" \
    --partition "${STAGE4_PARTITION:-low_p}" \
    "${SCRIPT_PATH}" "$@"
  exit 0
fi

# ── Inside Slurm job ───────────────────────────────────────────────────────────
start_time=$(date +%s)

PROFILE="${STAGE4_PROFILE:-slurm}"
STUDY="${STAGE4_STUDY:-all}"
ANALYSIS_ROOT="${STAGE4_ANALYSIS_ROOT:-${SCRATCH_RUN}/studies}"
STAGE2_ROOT="${STAGE4_STAGE2_ROOT:-${ANALYSIS_ROOT}}"
STAGE3_ROOT="${STAGE4_STAGE3_ROOT:-${ANALYSIS_ROOT}}"
DEST_ROOT="${STAGE4_DEST_ROOT:-${SCRATCH_RUN}/final}"
REPORT_DIR="${STAGE4_REPORT_DIR:-${DEST_ROOT}/report}"
SLURM_PARTITION="${STAGE4_PARTITION:-${SLURM_JOB_PARTITION:-low_p}}"
WORKDIR="${STAGE4_WORKDIR:-${SCRATCH_RUN}/stage4/work}"
CONDA_CACHE_DIR="${STAGE4_CONDA_CACHE_DIR:-${SCRATCH_RUN}/stage4/conda}"
CONDA_SOLVER="${STAGE4_CONDA_SOLVER:-classic}"
CONDA_CHANNEL_PRIORITY="${STAGE4_CONDA_CHANNEL_PRIORITY:-strict}"
CACHE_MODE="${STAGE4_CACHE_MODE:-lenient}"
RESUME=1

while [[ $# -gt 0 ]]; do
  case "$1" in
    --profile)        PROFILE="$2";        shift 2 ;;
    --study)          STUDY="$2";          shift 2 ;;
    --analysis-root)  ANALYSIS_ROOT="$2";  shift 2 ;;
    --stage2-root)    STAGE2_ROOT="$2";    shift 2 ;;
    --stage3-root)    STAGE3_ROOT="$2";    shift 2 ;;
    --dest-root)      DEST_ROOT="$2";      shift 2 ;;
    --report-dir)     REPORT_DIR="$2";     shift 2 ;;
    --partition)      SLURM_PARTITION="$2"; shift 2 ;;
    --no-resume)      RESUME=0;            shift   ;;
    -h|--help)        print_help; exit 0           ;;
    *) echo "ERROR: Unknown argument: $1" >&2; print_help >&2; exit 1 ;;
  esac
done

PIPELINE_DIR="${PROJ_ROOT}/pipeline_stage4"
NEXTFLOW_LOG="${SCRATCH_RUN}/stage4/.nextflow.log"

ANALYSIS_ROOT="$(abs_path "${ANALYSIS_ROOT}")"
STAGE2_ROOT="$(abs_path "${STAGE2_ROOT}")"
STAGE3_ROOT="$(abs_path "${STAGE3_ROOT}")"
DEST_ROOT="$(abs_path "${DEST_ROOT}")"
REPORT_DIR="$(abs_path "${REPORT_DIR}")"
WORKDIR="$(abs_path "${WORKDIR}")"
CONDA_CACHE_DIR="$(abs_path "${CONDA_CACHE_DIR}")"
NEXTFLOW_LOG="$(abs_path "${NEXTFLOW_LOG}")"

export NXF_OPTS="${NXF_OPTS:- -Xms2g -Xmx8g}"
export NXF_CONDA_CACHEDIR="${NXF_CONDA_CACHEDIR:-${CONDA_CACHE_DIR}}"
export CONDA_SOLVER
export CONDA_CHANNEL_PRIORITY
export PYTHON3_BIN="${PYTHON3_BIN:-python3}"

source "$(conda info --base)/etc/profile.d/conda.sh" || true
conda activate nf_EPIC-genetics || true
unset R_HOME
if [ -n "${CONDA_PREFIX:-}" ] && [ -x "${CONDA_PREFIX}/bin/java" ]; then
  export JAVA_HOME="${CONDA_PREFIX}"
  export JAVA_CMD="${CONDA_PREFIX}/bin/java"
else
  unset JAVA_CMD
  if [ -n "${JAVA_HOME:-}" ] && [ ! -x "${JAVA_HOME}/bin/java" ]; then unset JAVA_HOME; fi
fi

if command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base)"
  export CONDA_BASE
  export PATH="${CONDA_BASE}/bin:${CONDA_BASE}/condabin:${PATH}"
fi

if [ ! -d "${PIPELINE_DIR}" ]; then
  echo "ERROR: Stage-4 pipeline directory not found: ${PIPELINE_DIR}" >&2; exit 1
fi
if [ ! -d "${STAGE3_ROOT}" ]; then
  echo "ERROR: Stage-3 root not found: ${STAGE3_ROOT}" >&2; exit 1
fi
if [ ! -d "${STAGE2_ROOT}" ]; then
  echo "ERROR: Stage-2 root not found: ${STAGE2_ROOT}" >&2; exit 1
fi

mkdir -p "${WORKDIR}" "${CONDA_CACHE_DIR}" "${SCRATCH_RUN}/stage4" "${DEST_ROOT}" "${REPORT_DIR}"
cd "${SCRATCH_RUN}/stage4"

echo "=========================================="
echo " Running EPIC Genetics Stage-4 Pipeline"
echo "=========================================="
echo "Slurm job:    ${SLURM_JOB_ID:-interactive}"
echo "Host:         $(hostname)"
echo "Project root: ${PROJ_ROOT}"
echo "Pipeline:     ${PIPELINE_DIR}"
echo "Profile:      ${PROFILE}"
echo "Study:        ${STUDY}"
echo "Stage3 root:  ${STAGE3_ROOT}"
echo "Stage2 root:  ${STAGE2_ROOT}"
echo "Destination:  ${DEST_ROOT}"
echo "Report dir:   ${REPORT_DIR}"
echo "Work dir:     ${WORKDIR}"
echo "Cache mode:   ${CACHE_MODE}"
echo "=========================================="

nextflow_cmd=(
  nextflow
  -log "${NEXTFLOW_LOG}"
  run "${PIPELINE_DIR}"
  -params-file "${PIPELINE_DIR}/params.yaml"
  -profile "${PROFILE}"
  -work-dir "${WORKDIR}"
  --pipeline_info_dir "${SCRATCH_RUN}/stage4/pipeline_info"
  --study "${STUDY}"
  --project_root "${PROJ_ROOT}"
  --analysis_root "${ANALYSIS_ROOT}"
  --stage2_root "${STAGE2_ROOT}"
  --stage3_root "${STAGE3_ROOT}"
  --dest_root "${DEST_ROOT}"
  --report_dir "${REPORT_DIR}"
  --slurm_partition "${SLURM_PARTITION}"
  --cache_mode "${CACHE_MODE}"
)

[ "${RESUME}" = "1" ] && nextflow_cmd+=(-resume)

MAX_NF_ATTEMPTS=5
NF_ATTEMPT=0
NF_EXIT=1

while [ "${NF_ATTEMPT}" -lt "${MAX_NF_ATTEMPTS}" ]; do
  NF_ATTEMPT=$(( NF_ATTEMPT + 1 ))
  [ "${NF_ATTEMPT}" -gt 1 ] && \
    echo "Nextflow resume attempt ${NF_ATTEMPT}/${MAX_NF_ATTEMPTS} after NFS IOException..."

  LOG_START=$([ -f "${NEXTFLOW_LOG}" ] && wc -l < "${NEXTFLOW_LOG}" || echo 0)
  LOG_START=$(( LOG_START + 1 ))

  set +e
  "${nextflow_cmd[@]}"
  NF_EXIT=$?
  set -e

  [ "${NF_EXIT}" -eq 0 ] && break

  if ! tail -n +"${LOG_START}" "${NEXTFLOW_LOG}" 2>/dev/null \
          | grep -q "java.io.IOException: Invalid argument"; then
    echo "ERROR: Nextflow failed (exit ${NF_EXIT}) without NFS IOException — not retrying." >&2
    break
  fi

  echo "WARNING: Nextflow session aborted by NFS IOException." >&2
  tail -n +"${LOG_START}" "${NEXTFLOW_LOG}" 2>/dev/null | awk '
    /Handling unexpected condition/ { expect_wdir=1; wdir="" }
    expect_wdir && /work-dir=/ {
      n = split($0, a, "work-dir="); split(a[2], b, ";"); wdir = b[1]; expect_wdir=0; expect_exc=1
    }
    expect_exc && /java\.io\.IOException/ { print wdir; expect_exc=0; wdir="" }
  ' | sort -u | while IFS= read -r wdir; do
    [ -n "${wdir}" ] || continue
    { rm -rf "${wdir}" && echo "  Removed: ${wdir}" >&2; } 2>/dev/null || true
  done

  if ! printf '%s\0' "${nextflow_cmd[@]}" | grep -qz -- '-resume'; then
    nextflow_cmd+=(-resume)
  fi
  sleep 30
done

if [ "${NF_EXIT}" -ne 0 ]; then
  echo "ERROR: Nextflow pipeline failed after ${NF_ATTEMPT} attempt(s) (exit ${NF_EXIT})." >&2
  exit "${NF_EXIT}"
fi

# ── Copy pipeline summaries ────────────────────────────────────────────────────
echo ""
echo "Copying pipeline summaries to ${DEST_ROOT}/summaries/ ..."
mkdir -p "${DEST_ROOT}/summaries"
for summary in stage1-summary.md stage2-summary.md stage3-summary.md; do
  src="${ANALYSIS_ROOT}/${summary}"
  if [ -f "${src}" ]; then
    cp "${src}" "${DEST_ROOT}/summaries/${summary}"
    echo "  Copied: ${summary}"
  else
    echo "  WARNING: Not found, skipping: ${src}" >&2
  fi
done

end_time=$(date +%s)
elapsed=$((end_time - start_time))

echo ""
echo "Stage-4 pipeline complete."
echo "Destination: ${DEST_ROOT}"
echo "Elapsed: $(date -u -r "$elapsed" +%H:%M:%S 2>/dev/null || date -u -d "@$elapsed" +%H:%M:%S)"
