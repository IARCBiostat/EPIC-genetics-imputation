#!/bin/bash
#SBATCH --job-name=stage1
#SBATCH --output=src/logs/stage1.out
#SBATCH --error=src/logs/stage1.err
#SBATCH --time=10-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --partition=low_p

# Script: src/003_stage1.sh
# Purpose: Run stage-1 preprocessing for all bespoke dataset scripts in pipeline_stage1/scripts.

# -----------------------------------------------------------------------------
# Interactive single-study test
# Copy/paste these lines in an HPC shell to test Brea_01_Erneg directly:
#
# export PATH="/data/Epic/subprojects/Genetics/work/tools/bin:${PATH}"
# python3 /data/Epic/subprojects/Genetics/work/pipeline_stage1/scripts/process_brea_01_erneg.py \
#   --data-root /data/Epic/subprojects/Genetics/sources/Gwas \
#   --work-root /data/Epic/subprojects/Genetics/work/pipeline_stage1/work \
#   --plink plink \
#   --python2 python2.7
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
  candidates+=("${PROJ_ROOT}/../sources/Gwas" "${PROJ_ROOT}/sources/Gwas" "${PROJ_ROOT}")

  for candidate in "${candidates[@]}"; do
    [ -z "$candidate" ] && continue
    if candidate="$(cd "$candidate" 2>/dev/null && pwd)"; then
      if [ -f "${candidate}/Reference/Epic/Subj_Id_2015.txt" ]; then
        printf '%s\n' "$candidate"
        return 0
      fi
    fi
  done

  echo "ERROR: Could not resolve the archive-style stage-1 data root." >&2
  echo "       Set STAGE1_DATA_ROOT explicitly to the directory containing Reference/Epic/ and the study folders." >&2
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
for env_file in ".env" \
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

export PATH="${PROJ_ROOT}/tools/bin:${PATH}"

PYTHON3_BIN="${PYTHON3_BIN:-python3}"
PYTHON2_BIN="${PYTHON2_BIN:-python2.7}"
PLINK_BIN="${PLINK_BIN:-plink}"
DATA_ROOT="$(resolve_data_root)"
WORK_ROOT="${STAGE1_WORK_ROOT:-$PROJ_ROOT/pipeline_stage1/work}"
SCRIPTS_DIR="${PROJ_ROOT}/pipeline_stage1/scripts"
EPIC_ID_FILE="${DATA_ROOT}/Reference/Epic/Subj_Id_2015.txt"
STAGE1_FORCE="${STAGE1_FORCE:-0}"

DATASET_SCRIPTS=(
  "process_brea_01_erneg.py"
  "process_brea_02.py"
  # "process_clrt_01.py" # TO DO
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
  # "process_neur_01.py" # TO DO
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

for cmd_name in "${PYTHON3_BIN}" "${PYTHON2_BIN}" "${PLINK_BIN}" perl; do
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
  echo "       STAGE1_DATA_ROOT must point to the archive-style raw root containing Reference/Epic/." >&2
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
echo "Force rerun:  ${STAGE1_FORCE}"
echo "Datasets:     ${#DATASET_SCRIPTS[@]}"
echo "=========================================="

run_count=0
skip_count=0

for script_name in "${DATASET_SCRIPTS[@]}"; do
  script_path="${SCRIPTS_DIR}/${script_name}"
  if [ ! -f "$script_path" ]; then
    echo "ERROR: Dataset script not found: ${script_path}" >&2
    exit 1
  fi

  study_id="$(sed -n "s/^STUDY_ID = '\(.*\)'/\1/p" "$script_path")"
  raw_rel="$(sed -n "s/^RAW_REL = '\(.*\)'/\1/p" "$script_path")"
  manifest_rel="$(sed -n "s/^MANIFEST_REL = '\(.*\)'/\1/p" "$script_path")"
  id_link_rel="$(sed -n "s/^ID_LINK_REL = '\(.*\)'/\1/p" "$script_path")"

  if [ -z "${study_id}" ]; then
    echo "ERROR: Could not determine STUDY_ID from ${script_path}" >&2
    exit 1
  fi

  final_dir="${PROJ_ROOT}/analysis/${study_id}/stage1"
  final_prefix="${final_dir}/${study_id}"
  final_summary="${final_dir}/summary.txt"

  if [ "${STAGE1_FORCE}" != "1" ] && \
     [ -f "${final_prefix}.bed" ] && \
     [ -f "${final_prefix}.bim" ] && \
     [ -f "${final_prefix}.fam" ] && \
     [ -f "${final_summary}" ]; then
    echo ""
    echo ">>> Skipping ${script_name}"
    echo "    Final stage-1 outputs already exist in analysis/${study_id}/stage1/"
    skip_count=$((skip_count + 1))
    continue
  fi

  for ext in bed bim fam; do
    if [ ! -f "${DATA_ROOT}/${raw_rel}.${ext}" ]; then
      echo "ERROR: Missing raw input file: ${DATA_ROOT}/${raw_rel}.${ext}" >&2
      echo "       STAGE1_DATA_ROOT currently points to: ${DATA_ROOT}" >&2
      exit 1
    fi
  done

  if [ ! -f "${DATA_ROOT}/${manifest_rel}" ]; then
    echo "ERROR: Missing manifest file: ${DATA_ROOT}/${manifest_rel}" >&2
    exit 1
  fi

  if [ -n "${id_link_rel}" ] && [ ! -f "${DATA_ROOT}/${id_link_rel}" ]; then
    echo "ERROR: Missing ID linkage file: ${DATA_ROOT}/${id_link_rel}" >&2
    exit 1
  fi

  echo ""
  echo ">>> Running ${script_name}"
  "${PYTHON3_BIN}" "${script_path}" \
    --data-root "${DATA_ROOT}" \
    --work-root "${WORK_ROOT}" \
    --plink "${PLINK_BIN}" \
    --python2 "${PYTHON2_BIN}"
  run_count=$((run_count + 1))
done

echo ""
echo "Generating consolidated stage-1 summary..."
summary_output="${PROJ_ROOT}/analysis/stage1-summary.md"
"${PYTHON3_BIN}" "${PROJ_ROOT}/pipeline_stage1/scripts/summary.py" \
  --analysis-root "${PROJ_ROOT}/analysis" \
  --output "${summary_output}"

echo ""
echo "=========================================="
echo " Stage-1 processing completed"
echo " Executed: ${run_count}"
echo " Skipped:  ${skip_count}"
echo "=========================================="

end_time=$(date +%s)
elapsed=$((end_time - start_time))
echo "Total runtime: $(date -u -r "$elapsed" +%H:%M:%S 2>/dev/null || date -u -d "@$elapsed" +%H:%M:%S)"
