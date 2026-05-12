#!/bin/bash
#SBATCH --job-name=stage2
#SBATCH --output=src/logs/stage2.out
#SBATCH --error=src/logs/stage2.err
#SBATCH --time=10-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=low_p

set -euo pipefail
trap 'echo "ERROR: Imputation pipeline failed on line $LINENO" >&2; exit 1' ERR

SCRIPT_DIR="$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")" && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/$(basename -- "${BASH_SOURCE[0]:-$0}")"
DEFAULT_PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

is_repo_root() {
  local candidate="$1"
  [ -d "${candidate}/src" ] && \
  [ -d "${candidate}/pipeline_stage2" ] && \
  [ -d "${candidate}/analysis" ]
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
  return 1
}

abs_path() {
  local target="$1"
  local parent
  local base

  if [[ "$target" = /* ]]; then
    printf '%s\n' "$target"
    return 0
  fi

  parent="$(dirname "$target")"
  base="$(basename "$target")"

  if [ "$parent" = "." ]; then
    parent="$(pwd)"
  else
    parent="$(cd "$parent" 2>/dev/null && pwd)"
  fi

  printf '%s/%s\n' "$parent" "$base"
}

print_help() {
  cat <<'EOF'
Usage: bash src/005_stage2.sh [options]

Options:
  --profile <name>       Nextflow profile to use (default: slurm)
  --study <list>         Comma-separated study IDs to process (default: all)
  --chromosomes <list>   Comma-separated chromosomes to process (default: all; e.g. 22)
  --out <dir>            Analysis root where <STUDY>/stage2/ will be written
  --stage1-root <dir>    Root directory containing analysis/<STUDY>/stage1 outputs
  --partition <name>     Slurm partition for internal Nextflow task submissions
  --queue-size <n>       Nextflow Slurm executor queue size (default: 1000)
  --cache-mode <mode>    Nextflow task cache mode (default: deep)
  --no-empirical-validation
                         Disable typed-variant empirical R2 and Dose0 validation
  --empirical-min-samples <n>
                         Minimum matched samples for empirical R2 (default: 20)
  --dose0-min-samples <n>
                         Minimum observed genotype-0 samples for Dose0 (default: 5)
  --genetic-map-file <file>
                         Eagle genetic map file (default: <ref_1000g_dir>/genetic_map_hg38_withX.txt.gz)
  --no-resume            Disable Nextflow -resume
  -h, --help             Show this help

Examples:
  bash src/005_stage2.sh
  bash src/005_stage2.sh --study Brea_01_Erneg
  bash src/005_stage2.sh --study Glbd_01 --chromosomes 22
EOF
}

if [ -z "${SLURM_JOB_ID:-}" ]; then
  PROJ_ROOT="$(resolve_project_root)"
  mkdir -p "${PROJ_ROOT}/src/logs"
  cd "$PROJ_ROOT"

  if ! command -v sbatch >/dev/null 2>&1; then
    echo "ERROR: sbatch is not available. Submit this script from the HPC login node." >&2
    exit 1
  fi

  echo "Submitting imputation pipeline to Slurm..."
  sbatch \
    --export=ALL \
    --job-name "${IMPUTATION_JOB_NAME:-imputation_pipeline}" \
    --output "${IMPUTATION_LOG_OUT:-${PROJ_ROOT}/src/logs/004_pipeline_imputation_%j.out}" \
    --error "${IMPUTATION_LOG_ERR:-${PROJ_ROOT}/src/logs/004_pipeline_imputation_%j.err}" \
    --time "${IMPUTATION_TIME:-10-00:00:00}" \
    --mem "${IMPUTATION_MEM:-32G}" \
    --cpus-per-task "${IMPUTATION_CPUS:-4}" \
    --partition "${IMPUTATION_PARTITION:-low_p}" \
    "${SCRIPT_PATH}" "$@"
  exit 0
fi

start_time=$(date +%s)

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

PROFILE="${IMPUTATION_PROFILE:-slurm}"
STUDY="${IMPUTATION_STUDY:-all}"
CHROMOSOMES="${IMPUTATION_CHROMOSOMES:-all}"
OUTDIR="${IMPUTATION_OUTDIR:-${PROJ_ROOT}/analysis}"
STAGE1_ROOT="${IMPUTATION_STAGE1_ROOT:-${PROJ_ROOT}/analysis}"
SLURM_PARTITION="${IMPUTATION_PARTITION:-${SLURM_JOB_PARTITION:-low_p}}"
REF_1000G_DIR="${IMPUTATION_REF_1000G_DIR:-${PROJ_ROOT}/data/reference/1000G}"
FASTA_REF="${IMPUTATION_FASTA_REF:-${REF_1000G_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna}"
GENETIC_MAP_FILE="${IMPUTATION_GENETIC_MAP_FILE:-${PROJ_ROOT}/data/reference/eagle/genetic_map_hg38_withX.txt.gz}"
WORKDIR="${IMPUTATION_WORKDIR:-${PROJ_ROOT}/pipeline_stage2/work}"
CONDA_CACHE_DIR="${IMPUTATION_CONDA_CACHE_DIR:-${PROJ_ROOT}/pipeline_stage2/conda}"
CONDA_SOLVER="${IMPUTATION_CONDA_SOLVER:-classic}"
CONDA_CHANNEL_PRIORITY="${IMPUTATION_CONDA_CHANNEL_PRIORITY:-strict}"
CACHE_MODE="${IMPUTATION_CACHE_MODE:-deep}"
EXECUTOR_QUEUE_SIZE="${IMPUTATION_EXECUTOR_QUEUE_SIZE:-1000}"
SUMMARY_MIN_R2="${IMPUTATION_SUMMARY_MIN_R2:-0.3}"
SUMMARY_HIGH_QUALITY="${IMPUTATION_SUMMARY_HIGH_QUALITY:-0.8}"
RUN_EMPIRICAL_VALIDATION="${IMPUTATION_RUN_EMPIRICAL_VALIDATION:-true}"
EMPIRICAL_MIN_SAMPLES="${IMPUTATION_EMPIRICAL_MIN_SAMPLES:-20}"
DOSE0_MIN_SAMPLES="${IMPUTATION_DOSE0_MIN_SAMPLES:-5}"
RESUME=1
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --profile)
      PROFILE="$2"
      shift 2
      ;;
    --study)
      STUDY="$2"
      shift 2
      ;;
    --chromosomes|--chromosome|--chr)
      CHROMOSOMES="$2"
      shift 2
      ;;
    --out|--outdir|-o)
      OUTDIR="$2"
      shift 2
      ;;
    --stage1-root)
      STAGE1_ROOT="$2"
      shift 2
      ;;
    --partition)
      SLURM_PARTITION="$2"
      shift 2
      ;;
    --use-stage1-qc-handoff|--no-use-stage1-qc-handoff)
      echo "WARNING: stage-2 now always uses analysis/<STUDY>/stage1/<STUDY> handoff files." >&2
      shift
      ;;
    --cache-mode)
      CACHE_MODE="$2"
      shift 2
      ;;
    --queue-size)
      EXECUTOR_QUEUE_SIZE="$2"
      shift 2
      ;;
    --no-empirical-validation)
      RUN_EMPIRICAL_VALIDATION=false
      shift
      ;;
    --empirical-min-samples)
      EMPIRICAL_MIN_SAMPLES="$2"
      shift 2
      ;;
    --dose0-min-samples)
      DOSE0_MIN_SAMPLES="$2"
      shift 2
      ;;
    --ref-1000g-dir)
      REF_1000G_DIR="$2"
      shift 2
      ;;
    --fasta-ref)
      FASTA_REF="$2"
      shift 2
      ;;
    --genetic-map-file)
      GENETIC_MAP_FILE="$2"
      shift 2
      ;;
    --no-resume)
      RESUME=0
      shift
      ;;
    -h|--help)
      print_help
      exit 0
      ;;
    *)
      EXTRA_ARGS+=("$1")
      shift
      ;;
  esac
done

PIPELINE_DIR="${PROJ_ROOT}/pipeline_stage2"
PARAMS_FILE="${PIPELINE_DIR}/params.yaml"
OUTDIR="$(abs_path "${OUTDIR}")"
STAGE1_ROOT="$(abs_path "${STAGE1_ROOT}")"
REF_1000G_DIR="$(abs_path "${REF_1000G_DIR}")"
FASTA_REF="$(abs_path "${FASTA_REF}")"
GENETIC_MAP_FILE="$(abs_path "${GENETIC_MAP_FILE}")"
WORKDIR="$(abs_path "${WORKDIR}")"
CONDA_CACHE_DIR="$(abs_path "${CONDA_CACHE_DIR}")"

export NXF_OPTS="${NXF_OPTS:- -Xms4g -Xmx24g}"
export NXF_CONDA_CACHEDIR="${NXF_CONDA_CACHEDIR:-${CONDA_CACHE_DIR}}"
export CONDA_SOLVER
export CONDA_CHANNEL_PRIORITY
export PLINK_BIN="${PLINK_BIN:-plink}"
export BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
export EAGLE_BIN="${EAGLE_BIN:-eagle}"
export MINIMAC4_BIN="${MINIMAC4_BIN:-minimac4}"
export PYTHON3_BIN="${PYTHON3_BIN:-python3}"

source "$(conda info --base)/etc/profile.d/conda.sh" || true
conda activate nf_EPIC-genetics || true
unset R_HOME
if [ -n "${CONDA_PREFIX:-}" ] && [ -x "${CONDA_PREFIX}/bin/java" ]; then
  export JAVA_HOME="${CONDA_PREFIX}"
  export JAVA_CMD="${CONDA_PREFIX}/bin/java"
else
  unset JAVA_CMD
  if [ -n "${JAVA_HOME:-}" ] && [ ! -x "${JAVA_HOME}/bin/java" ]; then
    unset JAVA_HOME
  fi
fi

if [ ! -d "${PIPELINE_DIR}" ]; then
  echo "ERROR: Imputation pipeline directory not found: ${PIPELINE_DIR}" >&2
  exit 1
fi

if [ ! -f "${PARAMS_FILE}" ]; then
  echo "ERROR: Missing params file: ${PARAMS_FILE}" >&2
  exit 1
fi

if [ ! -d "${STAGE1_ROOT}" ]; then
  echo "ERROR: Stage-1 root directory not found: ${STAGE1_ROOT}" >&2
  exit 1
fi

if [ ! -d "${REF_1000G_DIR}" ]; then
  echo "ERROR: Reference directory not found: ${REF_1000G_DIR}" >&2
  exit 1
fi

if [ ! -f "${FASTA_REF}" ]; then
  echo "ERROR: Reference FASTA not found: ${FASTA_REF}" >&2
  exit 1
fi

if [ ! -f "${GENETIC_MAP_FILE}" ]; then
  echo "ERROR: Eagle genetic map file not found: ${GENETIC_MAP_FILE}" >&2
  exit 1
fi

if ! [[ "${EXECUTOR_QUEUE_SIZE}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --queue-size must be a positive integer: ${EXECUTOR_QUEUE_SIZE}" >&2
  exit 1
fi

if ! [[ "${EMPIRICAL_MIN_SAMPLES}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --empirical-min-samples must be a positive integer: ${EMPIRICAL_MIN_SAMPLES}" >&2
  exit 1
fi

if ! [[ "${DOSE0_MIN_SAMPLES}" =~ ^[1-9][0-9]*$ ]]; then
  echo "ERROR: --dose0-min-samples must be a positive integer: ${DOSE0_MIN_SAMPLES}" >&2
  exit 1
fi

check_stage1_study() {
  local study_id="$1"
  local prefix
  prefix="${STAGE1_ROOT}/${study_id}/stage1/${study_id}"
  local ext
  for ext in bed bim fam; do
    if [ ! -f "${prefix}.${ext}" ]; then
      echo "ERROR: Missing stage-1 input file: ${prefix}.${ext}" >&2
      return 1
    fi
  done
}

if [ "${STUDY}" = "all" ]; then
  stage1_count="$(find "${STAGE1_ROOT}" -mindepth 2 -maxdepth 2 -type d -name stage1 | wc -l | tr -d ' ')"
  if [ "${stage1_count}" = "0" ]; then
    echo "ERROR: No stage-1 study directories found under ${STAGE1_ROOT}" >&2
    exit 1
  fi
else
  IFS=',' read -r -a requested_studies <<< "${STUDY}"
  for study_id in "${requested_studies[@]}"; do
    study_id="$(echo "$study_id" | xargs)"
    [ -z "${study_id}" ] && continue
    check_stage1_study "${study_id}"
  done
fi

mkdir -p "${OUTDIR}"
mkdir -p "${WORKDIR}"
mkdir -p "${CONDA_CACHE_DIR}"
cd "${PIPELINE_DIR}"

NF_CMD=(nextflow run main.nf -profile "${PROFILE}" -work-dir "${WORKDIR}" -params-file "${PARAMS_FILE}" --outdir "${OUTDIR}" --study "${STUDY}" --chromosomes "${CHROMOSOMES}" --stage1_root "${STAGE1_ROOT}" --slurm_partition "${SLURM_PARTITION}" --executor_queue_size "${EXECUTOR_QUEUE_SIZE}" --ref_1000g_dir "${REF_1000G_DIR}" --fasta_ref "${FASTA_REF}" --genetic_map_file "${GENETIC_MAP_FILE}" --cache_mode "${CACHE_MODE}" --summary_min_r2 "${SUMMARY_MIN_R2}" --summary_high_quality "${SUMMARY_HIGH_QUALITY}" --run_empirical_validation "${RUN_EMPIRICAL_VALIDATION}" --empirical_min_samples "${EMPIRICAL_MIN_SAMPLES}" --dose0_min_samples "${DOSE0_MIN_SAMPLES}")
if [ "${RESUME}" = "1" ]; then
  NF_CMD+=(-resume)
fi
if [ "${#EXTRA_ARGS[@]}" -gt 0 ]; then
  NF_CMD+=("${EXTRA_ARGS[@]}")
fi

echo "=========================================="
echo " Launching EPIC Imputation Pipeline"
echo "=========================================="
echo "Slurm job:    ${SLURM_JOB_ID:-interactive}"
echo "Host:         $(hostname)"
echo "Project root: ${PROJ_ROOT}"
echo "Pipeline:     ${PIPELINE_DIR}"
echo "Work dir:     ${WORKDIR}"
echo "Stage1 root:  ${STAGE1_ROOT}"
echo "Conda cache:  ${NXF_CONDA_CACHEDIR}"
echo "Conda solver: ${CONDA_SOLVER}"
echo "Channel prio: ${CONDA_CHANNEL_PRIORITY}"
echo "Task cache:   ${CACHE_MODE}"
echo "Queue size:   ${EXECUTOR_QUEUE_SIZE}"
echo "Stage1 input: ${STAGE1_ROOT}/<STUDY>/stage1/<STUDY>"
echo "Ref 1000G:    ${REF_1000G_DIR}"
echo "FASTA:        ${FASTA_REF}"
echo "Genetic map:  ${GENETIC_MAP_FILE}"
echo "PLINK bin:    ${PLINK_BIN}"
echo "BCFtools bin: ${BCFTOOLS_BIN}"
echo "Eagle bin:    ${EAGLE_BIN}"
echo "Minimac4 bin: ${MINIMAC4_BIN}"
echo "Study:        ${STUDY}"
echo "Chromosomes:  ${CHROMOSOMES}"
echo "Summary R2:   ${SUMMARY_MIN_R2}"
echo "High R2:      ${SUMMARY_HIGH_QUALITY}"
echo "Empirical validation: ${RUN_EMPIRICAL_VALIDATION}"
echo "Empirical min samples:${EMPIRICAL_MIN_SAMPLES}"
echo "Dose0 min samples:    ${DOSE0_MIN_SAMPLES}"
echo "Profile:      ${PROFILE}"
echo "Partition:    ${SLURM_PARTITION}"
echo "Analysis root:${OUTDIR}"
echo "Resume:       ${RESUME}"
echo "Stage2 summary: ${OUTDIR}/stage2-summary.md"
echo "Study report:   ${OUTDIR}/<STUDY>/stage2/report/report-stage2.html"
echo "=========================================="
printf 'Executing:\n'
printf '%q ' "${NF_CMD[@]}"
printf '\n'

"${NF_CMD[@]}"

end_time=$(date +%s)
elapsed=$((end_time - start_time))
echo "Total runtime: $(date -u -r "$elapsed" +%H:%M:%S 2>/dev/null || date -u -d "@$elapsed" +%H:%M:%S)"
