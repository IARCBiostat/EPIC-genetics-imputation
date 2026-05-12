#!/bin/bash
#SBATCH --job-name=stage1
#SBATCH --output=src/logs/stage1.out
#SBATCH --error=src/logs/stage1.err
#SBATCH --time=10-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --partition=low_p

# Script: src/004_stage1.sh
# Purpose: Run stage-1 preprocessing for all bespoke dataset scripts in pipeline_stage1/scripts.

# -----------------------------------------------------------------------------
# Interactive single-study test
# Copy/paste these lines in an HPC shell to test Brea_01_Erneg directly:
#
# cd /data/Epic/subprojects/Genetics/work
# set -a
# source pipeline_stage1/.env
# set +a
# python3 /data/Epic/subprojects/Genetics/work/pipeline_stage1/scripts/process_brea_01_erneg.py \
#   --data-root "${STAGE1_DATA_ROOT}" \
#   --work-root "${STAGE1_WORK_ROOT}" \
#   --plink "${PLINK_BIN}" \
#   --python2 "${PYTHON2_BIN}"
# -----------------------------------------------------------------------------

set -euo pipefail
trap 'echo "ERROR: Stage-1 processing failed on line $LINENO" >&2; exit 1' ERR

SCRIPT_DIR="$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")" && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/$(basename -- "${BASH_SOURCE[0]:-$0}")"
DEFAULT_PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

is_repo_root() {
  local candidate="$1"
  [ -d "${candidate}/src" ] && \
  [ -d "${candidate}/pipeline_stage1/scripts" ] && \
  [ -d "${candidate}/tools" ]
}

resolve_project_root() {
  local candidates=()
  [ -n "${SLURM_SUBMIT_DIR:-}" ] && candidates+=("${SLURM_SUBMIT_DIR}")
  candidates+=("${SCRIPT_DIR}/.." "$(pwd)")
  [ -n "${GENETICS_PROJECT_ROOT:-}" ] && candidates+=("${GENETICS_PROJECT_ROOT}")
  candidates+=("${DEFAULT_PROJ_ROOT}")

  local candidate
  local resolved
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

resolve_data_root() {
  local candidate
  local candidates=()

  [ -n "${STAGE1_DATA_ROOT:-}" ] && candidates+=("${STAGE1_DATA_ROOT}")
  candidates+=("${PROJ_ROOT}/data" "${PROJ_ROOT}")

  for candidate in "${candidates[@]}"; do
    [ -z "$candidate" ] && continue
    if candidate="$(cd "$candidate" 2>/dev/null && pwd)"; then
      if [ -d "${candidate}/genetics" ] && \
         { [ -f "${candidate}/reference/Epic/Subj_Id_2015.txt" ] || \
           [ -f "${candidate}/Reference/Epic/Subj_Id_2015.txt" ]; }; then
        printf '%s\n' "$candidate"
        return 0
      fi
    fi
  done

  echo "ERROR: Could not resolve the stage-1 data root." >&2
  echo "       Set STAGE1_DATA_ROOT to the directory containing genetics/ and reference/Epic/." >&2
  return 1
}

resolve_epic_file() {
  local filename="$1"
  local candidate
  local candidates=()

  candidates+=("${DATA_ROOT}/reference/Epic/${filename}" "${DATA_ROOT}/Reference/Epic/${filename}")

  for candidate in "${candidates[@]}"; do
    if [ -f "$candidate" ]; then
      printf '%s\n' "$candidate"
      return 0
    fi
  done

  echo "ERROR: Could not find EPIC reference file: ${filename}" >&2
  return 1
}

if [ -z "${SLURM_JOB_ID:-}" ]; then
  PROJ_ROOT="$(resolve_project_root)"
  mkdir -p "${PROJ_ROOT}/src/logs"
  cd "$PROJ_ROOT"

  if ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch is not available. Run this script on the HPC login node or submit it with Slurm." >&2
    exit 1
  fi

  echo "Submitting Stage-1 processing to Slurm..."
  sbatch \
    --export=ALL \
    --job-name "${STAGE1_JOB_NAME:-stage1_process}" \
    --output "${STAGE1_LOG_OUT:-${PROJ_ROOT}/src/logs/003_process_%j.out}" \
    --error "${STAGE1_LOG_ERR:-${PROJ_ROOT}/src/logs/003_process_%j.err}" \
    --time "${STAGE1_TIME:-10-00:00:00}" \
    --mem "${STAGE1_MEM:-32G}" \
    --cpus-per-task "${STAGE1_CPUS:-2}" \
    --partition "${STAGE1_PARTITION:-low_p}" \
    "${SCRIPT_PATH}"
  exit 0
fi

start_time=$(date +%s)

# Robust environment sourcing
for env_file in "${DEFAULT_PROJ_ROOT}/pipeline_stage1/.env" \
                "${SCRIPT_DIR}/../pipeline_stage1/.env" \
                ".env" \
                "${DEFAULT_PROJ_ROOT}/.env" \
                "${SCRIPT_DIR}/../.env" \
                "${DEFAULT_PROJ_ROOT}/pipeline_stage2/.env" \
                "${SCRIPT_DIR}/../pipeline_stage2/.env"; do
  if [ -f "$env_file" ]; then
    set -a
    source "$env_file"
    set +a
    break
  fi
done

PROJ_ROOT="$(resolve_project_root)"
export GENETICS_PROJECT_ROOT="${PROJ_ROOT}"

if [ -n "${GENETICS_TOOLS_BIN:-}" ]; then
  export PATH="${GENETICS_TOOLS_BIN}:${PATH}"
elif [ -d "${PROJ_ROOT}/tools/bin" ]; then
  export PATH="${PROJ_ROOT}/tools/bin:${PATH}"
fi

PYTHON3_BIN="${PYTHON3_BIN:-python3}"
PYTHON2_BIN="${PYTHON2_BIN:-python2.7}"
PLINK_BIN="${PLINK_BIN:-plink}"
if [[ "${PLINK_BIN}" == */* ]]; then
  export PATH="$(dirname -- "${PLINK_BIN}"):${PATH}"
fi
export PLINK_BIN
DATA_ROOT="$(resolve_data_root)"
WORK_ROOT="${STAGE1_WORK_ROOT:-$PROJ_ROOT/pipeline_stage1/work}"
SCRIPTS_DIR="${PROJ_ROOT}/pipeline_stage1/scripts"
EPIC_ID_FILE="$(resolve_epic_file "Subj_Id_2015.txt")"
EPIC_CASE_STATUS_FILE="${EPIC_CASE_STATUS_FILE:-$(resolve_epic_file "EPIC_study_case_status.txt")}"
STAGE1_FORCE="${STAGE1_FORCE:-0}"
export EPIC_CASE_STATUS_FILE

if command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base)"
  # shellcheck disable=SC1091
  source "${CONDA_BASE}/etc/profile.d/conda.sh" || true
  conda activate nf_EPIC-genetics || true
  export PATH="${CONDA_BASE}/bin:${CONDA_BASE}/condabin:${PATH}"
fi

if [ -n "${CONDA_PREFIX:-}" ] && [ -x "${CONDA_PREFIX}/bin/java" ]; then
  export JAVA_HOME="${CONDA_PREFIX}"
  export JAVA_CMD="${CONDA_PREFIX}/bin/java"
else
  unset JAVA_CMD
  if [ -n "${JAVA_HOME:-}" ] && [ ! -x "${JAVA_HOME}/bin/java" ]; then
    unset JAVA_HOME
  fi
fi

DATASET_SCRIPTS=(
  "process_brea_01_erneg.py"
  "process_brea_02.py"
  "process_clrt_01.py"
  # "process_corp_01.py" # corpus uteri data isnt available 15/03/2026
  "process_ecvd_01.py"
  "process_ecvd_02.py"
  "process_ecvd_03.py"
  "process_glbd_01.py"
  "process_inte_01.py"
  "process_inte_02.py"
  "process_inte_03.py"
  "process_kidn_01.py"
  "process_kidn_02.py"
  "process_lung_01.py"
  "process_lymp_01.py"
  "process_neuro_01.py"
  "process_ovar_01.py"
  "process_panc_01.py"
  "process_panc_02.py"
  "process_pros_01.py"
  "process_pros_02.py"
  "process_pros_03.py"
  "process_pros_04.py"
  "process_stom_01.py"
  "process_uadt_01.py"
)

if [ -n "${STAGE1_SCRIPTS:-}" ]; then
  IFS=',' read -r -a DATASET_SCRIPTS <<< "${STAGE1_SCRIPTS}"
fi

NORMALIZED_SCRIPTS=()
for script_name in "${DATASET_SCRIPTS[@]}"; do
  script_name="$(echo "$script_name" | xargs)"
  if [[ "$script_name" != process_* ]]; then
    script_name="process_${script_name}"
  fi
  if [[ "$script_name" != *.py ]]; then
    script_name="${script_name}.py"
  fi
  NORMALIZED_SCRIPTS+=("$script_name")
done
DATASET_SCRIPTS=("${NORMALIZED_SCRIPTS[@]}")

mkdir -p "$WORK_ROOT"
mkdir -p "${PROJ_ROOT}/src/logs"
cd "$PROJ_ROOT"

for cmd_name in "${PYTHON3_BIN}" "${PYTHON2_BIN}" "${PLINK_BIN}" perl nextflow; do
  if ! command -v "$cmd_name" >/dev/null 2>&1; then
    echo "ERROR: Required command not found in PATH: ${cmd_name}" >&2
    exit 1
  fi
done

if [ ! -d "$SCRIPTS_DIR" ]; then
  echo "ERROR: Scripts directory not found: ${SCRIPTS_DIR}" >&2
  exit 1
fi

if [ ! -f "$EPIC_ID_FILE" ]; then
  echo "ERROR: Shared EPIC ID reference not found: ${EPIC_ID_FILE}" >&2
  echo "       STAGE1_DATA_ROOT must point to the data root containing genetics/ and reference/Epic/." >&2
  exit 1
fi

if [ ! -f "$EPIC_CASE_STATUS_FILE" ]; then
  echo "ERROR: Shared EPIC case-status reference not found: ${EPIC_CASE_STATUS_FILE}" >&2
  echo "       Run Rscript src/003_data-epic.R first or set EPIC_CASE_STATUS_FILE explicitly." >&2
  exit 1
fi

echo "=========================================="
echo " Running EPIC Genetics Stage-1 Processing"
echo "=========================================="
echo "Slurm job:     ${SLURM_JOB_ID:-interactive}"
echo "Host:          $(hostname)"
echo "Project root: ${PROJ_ROOT}"
echo "Data root:    ${DATA_ROOT}"
echo "Work root:    ${WORK_ROOT}"
echo "Python3:      ${PYTHON3_BIN}"
echo "Python2:      ${PYTHON2_BIN}"
echo "PLINK:        ${PLINK_BIN}"
echo "EPIC status:  ${EPIC_CASE_STATUS_FILE}"
echo "Force rerun:  ${STAGE1_FORCE}"
echo "Datasets:     ${#DATASET_SCRIPTS[@]}"
echo "=========================================="

if [ -n "${STAGE1_STUDY:-}" ]; then
  STUDY_LIST="${STAGE1_STUDY}"
elif [ -n "${STAGE1_SCRIPTS:-}" ]; then
  STUDY_IDS=()
  for script_name in "${DATASET_SCRIPTS[@]}"; do
    script_path="${SCRIPTS_DIR}/${script_name}"
    if [ ! -f "$script_path" ]; then
      echo "ERROR: Dataset script not found: ${script_path}" >&2
      exit 1
    fi
    study_id="$(sed -n "s/^STUDY_ID = '\(.*\)'/\1/p" "$script_path")"
    if [ -z "${study_id}" ]; then
      echo "ERROR: Could not determine STUDY_ID from ${script_path}" >&2
      exit 1
    fi
    STUDY_IDS+=("${study_id}")
  done
  STUDY_LIST="$(IFS=,; echo "${STUDY_IDS[*]}")"
else
  STUDY_LIST="all"
fi

OUTDIR="${STAGE1_OUTDIR:-${PROJ_ROOT}/analysis}"
NEXTFLOW_WORKDIR="${STAGE1_NEXTFLOW_WORKDIR:-${PROJ_ROOT}/pipeline_stage1/.nextflow_work}"
PROFILE="${STAGE1_PROFILE:-slurm}"
PARAMS_FILE="${PROJ_ROOT}/pipeline_stage1/params.yaml"
RESUME="${STAGE1_RESUME:-1}"

mkdir -p "${OUTDIR}" "${WORK_ROOT}" "${NEXTFLOW_WORKDIR}"

NF_CMD=(
  nextflow run "${PROJ_ROOT}/pipeline_stage1/main.nf"
  -profile "${PROFILE}"
  -work-dir "${NEXTFLOW_WORKDIR}"
  -params-file "${PARAMS_FILE}"
  --study "${STUDY_LIST}"
  --data_root "${DATA_ROOT}"
  --outdir "${OUTDIR}"
  --work_root "${WORK_ROOT}"
  --python3_bin "${PYTHON3_BIN}"
  --python2_bin "${PYTHON2_BIN}"
  --plink_bin "${PLINK_BIN}"
)

if [ "${RESUME}" = "1" ] && [ "${STAGE1_FORCE}" != "1" ]; then
  NF_CMD+=(-resume)
fi

echo ""
echo "Launching stage-1 Nextflow pipeline..."
printf '%q ' "${NF_CMD[@]}"
printf '\n'
"${NF_CMD[@]}"

echo ""
echo "Generating consolidated stage-1 summary..."
summary_output="${OUTDIR}/stage1-summary.md"
"${PYTHON3_BIN}" "${PROJ_ROOT}/pipeline_stage1/scripts/summary.py" \
  --analysis-root "${OUTDIR}" \
  --output "${summary_output}"

echo ""
echo "Generating stage-1 figures, tables, and report assets..."
stage1_report_cmd=(
  "${PYTHON3_BIN}" "${PROJ_ROOT}/pipeline_stage1/scripts/run_stage1_reports.py"
  --analysis-root "${OUTDIR}"
  --studies "${STUDY_LIST}"
)
if [ "${STAGE1_FORCE}" = "1" ]; then
  stage1_report_cmd+=(--force)
fi
printf '%q ' "${stage1_report_cmd[@]}"
printf '\n'
"${stage1_report_cmd[@]}"

echo ""
echo "=========================================="
echo " Stage-1 processing completed"
echo " Studies: ${STUDY_LIST}"
echo "=========================================="

end_time=$(date +%s)
elapsed=$((end_time - start_time))
echo "Total runtime: $(date -u -r "$elapsed" +%H:%M:%S 2>/dev/null || date -u -d "@$elapsed" +%H:%M:%S)"
