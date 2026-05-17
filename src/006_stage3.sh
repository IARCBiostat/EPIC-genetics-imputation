#!/bin/bash
#SBATCH --job-name=stage3
#SBATCH --output=src/logs/stage3.out
#SBATCH --error=src/logs/stage3.err
#SBATCH --time=10-00:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=low_p

set -euo pipefail
trap 'echo "ERROR: Stage-3 pipeline failed on line $LINENO" >&2; exit 1' ERR

SCRIPT_DIR="$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")" && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/$(basename -- "${BASH_SOURCE[0]:-$0}")"
DEFAULT_PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}"

is_repo_root() {
  local candidate="$1"
  [ -d "${candidate}/src" ] && \
  [ -d "${candidate}/pipeline_stage3" ] && \
  [ -d "${candidate}/analysis" ]
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

resolve_dbsnp_vcf() {
  local candidate
  local search_dir="${PROJ_ROOT}/data/reference/dbsnp"

  [ -n "${STAGE3_DBSNP_VCF:-}" ] && candidate="${STAGE3_DBSNP_VCF}"
  if [ -n "${candidate:-}" ] && [ -f "${candidate}" ]; then
    printf '%s\n' "$(abs_path "${candidate}")"
    return 0
  fi

  candidate="$(find "${search_dir}" -maxdepth 1 -type f -name 'GCF_000001405.*.gz' ! -name '*.tbi' | sort -V | tail -n 1 || true)"
  if [ -n "${candidate}" ] && [ -f "${candidate}" ]; then
    printf '%s\n' "$(abs_path "${candidate}")"
    return 0
  fi

  echo "ERROR: Could not resolve a dbSNP GRCh38 VCF under ${search_dir}" >&2
  return 1
}

resolve_dbsnp_tbi() {
  local dbsnp_vcf="$1"
  local candidate="${STAGE3_DBSNP_TBI:-${dbsnp_vcf}.tbi}"
  if [ -f "${candidate}" ]; then
    printf '%s\n' "$(abs_path "${candidate}")"
    return 0
  fi

  echo "ERROR: dbSNP index not found: ${candidate}" >&2
  return 1
}

print_help() {
  cat <<'EOF'
Usage: sbatch src/006_stage3.sh [options]

Options:
  --profile <name>                  Nextflow profile to use (default: slurm)
  --study <list>                    Comma-separated study IDs to process (default: all)
  --out <dir>                       Analysis root where <STUDY>/stage3/ will be written
  --stage1-root <dir>               Root directory containing analysis/<STUDY>/stage1 outputs
  --stage2-root <dir>               Root directory containing analysis/<STUDY>/stage2 outputs
  --dbsnp-vcf <file>                dbSNP GRCh38 VCF for rsID annotation
  --dbsnp-tbi <file>                dbSNP GRCh38 VCF index
  --partition <name>                Slurm partition for internal Nextflow task submissions
  --prep-dbsnp-time <duration>            Wall time for PREP_DBSNP_CHROM jobs (default: 72h)
  --filter-chrom-time <duration>          Wall time for FILTER_CHROM jobs (default: 72h)
  --annotate-chrom-time <duration>        Wall time for ANNOTATE_CHROM jobs (default: 72h)
  --import-chrom-time <duration>          Wall time for IMPORT_CHROM jobs (default: 72h)
  --hwe-chrom-time <duration>             Wall time for HWE_CHROM jobs (default: 72h)
  --merge-study-time <duration>           Wall time for MERGE_STUDY jobs (default: 72h)
  --sex-check-time <duration>             Wall time for SEX_CHECK jobs (default: 72h)
  --prune-autosomes-time <duration>       Wall time for PRUNE_AUTOSOMES jobs (default: 72h)
  --king-qc-time <duration>               Wall time for KING_QC jobs (default: 72h)
  --het-pca-qc-time <duration>            Wall time for HET_PCA_QC jobs (default: 72h)
  --sample-review-summary-time <duration> Wall time for SAMPLE_REVIEW_SUMMARY jobs (default: 72h)
  --finalize-study-time <duration>        Wall time for FINALIZE_STUDY jobs (default: 72h)
  --publish-intermediate-plink      Copy intermediate PLINK/PGEN/BED files into analysis output
  --no-publish-intermediate-plink   Do not copy intermediate PLINK/PGEN/BED files (default)
  --min-r2 <value>                  Minimum imputation R2 (default: 0.3)
  --maf <value>                     Minimum MAF (default: 0.01)
  --hwe-p <value>                   HWE p-value threshold (default: 0.000005)
  --run-hwe                         Re-enable stage-3 HWE filtering (disabled by default; stage-1 handoff owns HWE)
  --no-hwe                          Disable HWE filtering
  --exclude-ancestry-outliers       Exclude ancestry outliers from the final dataset
  --no-exclude-ancestry-outliers    Keep ancestry outliers in the final dataset
  --no-resume                       Disable Nextflow -resume
  -h, --help                        Show this help
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

  echo "Submitting stage-3 pipeline to Slurm..."
  sbatch \
    --export=ALL \
    --job-name "${STAGE3_JOB_NAME:-stage3_postimpute}" \
    --output "${STAGE3_LOG_OUT:-${PROJ_ROOT}/src/logs/006_stage3_%j.out}" \
    --error "${STAGE3_LOG_ERR:-${PROJ_ROOT}/src/logs/006_stage3_%j.err}" \
    --time "${STAGE3_TIME:-10-00:00:00}" \
    --mem "${STAGE3_MEM:-32G}" \
    --cpus-per-task "${STAGE3_CPUS:-4}" \
    --partition "${STAGE3_PARTITION:-low_p}" \
    "${SCRIPT_PATH}" "$@"
  exit 0
fi

start_time=$(date +%s)

PROFILE="${STAGE3_PROFILE:-slurm}"
STUDY="${STAGE3_STUDY:-all}"
OUTDIR="${STAGE3_OUTDIR:-${PROJ_ROOT}/analysis}"
STAGE1_ROOT="${STAGE3_STAGE1_ROOT:-${PROJ_ROOT}/analysis}"
STAGE2_ROOT="${STAGE3_STAGE2_ROOT:-${PROJ_ROOT}/analysis}"
SLURM_PARTITION="${STAGE3_PARTITION:-${SLURM_JOB_PARTITION:-low_p}}"
PREP_DBSNP_TIME="${STAGE3_PREP_DBSNP_TIME:-72h}"
FILTER_CHROM_TIME="${STAGE3_FILTER_CHROM_TIME:-72h}"
ANNOTATE_CHROM_TIME="${STAGE3_ANNOTATE_CHROM_TIME:-72h}"
IMPORT_CHROM_TIME="${STAGE3_IMPORT_CHROM_TIME:-72h}"
HWE_CHROM_TIME="${STAGE3_HWE_CHROM_TIME:-72h}"
MERGE_STUDY_TIME="${STAGE3_MERGE_STUDY_TIME:-72h}"
SEX_CHECK_TIME="${STAGE3_SEX_CHECK_TIME:-72h}"
PRUNE_AUTOSOMES_TIME="${STAGE3_PRUNE_AUTOSOMES_TIME:-72h}"
KING_QC_TIME="${STAGE3_KING_QC_TIME:-72h}"
HET_PCA_QC_TIME="${STAGE3_HET_PCA_QC_TIME:-72h}"
SAMPLE_REVIEW_SUMMARY_TIME="${STAGE3_SAMPLE_REVIEW_SUMMARY_TIME:-72h}"
FINALIZE_STUDY_TIME="${STAGE3_FINALIZE_STUDY_TIME:-72h}"
PUBLISH_INTERMEDIATE_PLINK="${STAGE3_PUBLISH_INTERMEDIATE_PLINK:-false}"
WORKDIR="${STAGE3_WORKDIR:-${PROJ_ROOT}/pipeline_stage3/work}"
CONDA_CACHE_DIR="${STAGE3_CONDA_CACHE_DIR:-${PROJ_ROOT}/pipeline_stage3/conda}"
CONDA_SOLVER="${STAGE3_CONDA_SOLVER:-classic}"
CONDA_CHANNEL_PRIORITY="${STAGE3_CONDA_CHANNEL_PRIORITY:-strict}"
CACHE_MODE="${STAGE3_CACHE_MODE:-lenient}"
MIN_R2="${STAGE3_MIN_R2:-0.3}"
MAF="${STAGE3_MAF:-0.01}"
RUN_HWE="${STAGE3_RUN_HWE:-true}"
HWE_P="${STAGE3_HWE_P:-0.000005}"
HWE_K="${STAGE3_HWE_K:-0}"
KING_CUTOFF="${STAGE3_KING_CUTOFF:-0.0884}"
ANCESTRY_PC_COUNT="${STAGE3_ANCESTRY_PC_COUNT:-10}"
ANCESTRY_Z_THRESHOLD="${STAGE3_ANCESTRY_Z_THRESHOLD:-6.0}"
HET_SD_THRESHOLD="${STAGE3_HET_SD_THRESHOLD:-3.0}"
EXCLUDE_ANCESTRY_OUTLIERS="${STAGE3_EXCLUDE_ANCESTRY_OUTLIERS:-true}"
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
    --out|--outdir|-o)
      OUTDIR="$2"
      shift 2
      ;;
    --stage1-root)
      STAGE1_ROOT="$2"
      shift 2
      ;;
    --stage2-root)
      STAGE2_ROOT="$2"
      shift 2
      ;;
    --dbsnp-vcf)
      STAGE3_DBSNP_VCF="$2"
      shift 2
      ;;
    --dbsnp-tbi)
      STAGE3_DBSNP_TBI="$2"
      shift 2
      ;;
    --partition)
      SLURM_PARTITION="$2"
      shift 2
      ;;
    --prep-dbsnp-time)
      PREP_DBSNP_TIME="$2"
      shift 2
      ;;
    --filter-chrom-time)
      FILTER_CHROM_TIME="$2"
      shift 2
      ;;
    --annotate-chrom-time)
      ANNOTATE_CHROM_TIME="$2"
      shift 2
      ;;
    --import-chrom-time)
      IMPORT_CHROM_TIME="$2"
      shift 2
      ;;
    --hwe-chrom-time)
      HWE_CHROM_TIME="$2"
      shift 2
      ;;
    --merge-study-time)
      MERGE_STUDY_TIME="$2"
      shift 2
      ;;
    --sex-check-time)
      SEX_CHECK_TIME="$2"
      shift 2
      ;;
    --prune-autosomes-time)
      PRUNE_AUTOSOMES_TIME="$2"
      shift 2
      ;;
    --king-qc-time)
      KING_QC_TIME="$2"
      shift 2
      ;;
    --het-pca-qc-time)
      HET_PCA_QC_TIME="$2"
      shift 2
      ;;
    --sample-review-summary-time)
      SAMPLE_REVIEW_SUMMARY_TIME="$2"
      shift 2
      ;;
    --finalize-study-time)
      FINALIZE_STUDY_TIME="$2"
      shift 2
      ;;
    --publish-intermediate-plink)
      PUBLISH_INTERMEDIATE_PLINK="true"
      shift
      ;;
    --no-publish-intermediate-plink)
      PUBLISH_INTERMEDIATE_PLINK="false"
      shift
      ;;
    --min-r2)
      MIN_R2="$2"
      shift 2
      ;;
    --maf)
      MAF="$2"
      shift 2
      ;;
    --hwe-p)
      HWE_P="$2"
      shift 2
      ;;
    --run-hwe)
      RUN_HWE="true"
      shift
      ;;
    --no-hwe)
      RUN_HWE="false"
      shift
      ;;
    --exclude-ancestry-outliers)
      EXCLUDE_ANCESTRY_OUTLIERS="true"
      shift
      ;;
    --no-exclude-ancestry-outliers)
      EXCLUDE_ANCESTRY_OUTLIERS="false"
      shift
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

PIPELINE_DIR="${PROJ_ROOT}/pipeline_stage3"
PARAMS_FILE="${PIPELINE_DIR}/params.yaml"
SUMMARY_SCRIPT="${PIPELINE_DIR}/bin/summary.py"
SUMMARY_OUTPUT="${OUTDIR}/stage3-summary.md"
NEXTFLOW_LOG="${PIPELINE_DIR}/.nextflow.log"

OUTDIR="$(abs_path "${OUTDIR}")"
STAGE1_ROOT="$(abs_path "${STAGE1_ROOT}")"
STAGE2_ROOT="$(abs_path "${STAGE2_ROOT}")"
WORKDIR="$(abs_path "${WORKDIR}")"
CONDA_CACHE_DIR="$(abs_path "${CONDA_CACHE_DIR}")"
SUMMARY_OUTPUT="$(abs_path "${SUMMARY_OUTPUT}")"
NEXTFLOW_LOG="$(abs_path "${NEXTFLOW_LOG}")"

DBSNP_VCF="$(resolve_dbsnp_vcf)"
DBSNP_TBI="$(resolve_dbsnp_tbi "${DBSNP_VCF}")"

export NXF_OPTS="${NXF_OPTS:- -Xms4g -Xmx24g}"
export NXF_CONDA_CACHEDIR="${NXF_CONDA_CACHEDIR:-${CONDA_CACHE_DIR}}"
export CONDA_SOLVER
export CONDA_CHANNEL_PRIORITY
export BCFTOOLS_BIN="${BCFTOOLS_BIN:-bcftools}"
export PLINK_BIN="${PLINK_BIN:-plink}"
export PLINK2_BIN="${PLINK2_BIN:-plink2}"
export PYTHON3_BIN="${PYTHON3_BIN:-python3}"
export STAGE3_PREP_DBSNP_TIME="${PREP_DBSNP_TIME}"
export STAGE3_FILTER_CHROM_TIME="${FILTER_CHROM_TIME}"
export STAGE3_ANNOTATE_CHROM_TIME="${ANNOTATE_CHROM_TIME}"
export STAGE3_IMPORT_CHROM_TIME="${IMPORT_CHROM_TIME}"
export STAGE3_HWE_CHROM_TIME="${HWE_CHROM_TIME}"
export STAGE3_MERGE_STUDY_TIME="${MERGE_STUDY_TIME}"
export STAGE3_SEX_CHECK_TIME="${SEX_CHECK_TIME}"
export STAGE3_PRUNE_AUTOSOMES_TIME="${PRUNE_AUTOSOMES_TIME}"
export STAGE3_KING_QC_TIME="${KING_QC_TIME}"
export STAGE3_HET_PCA_QC_TIME="${HET_PCA_QC_TIME}"
export STAGE3_SAMPLE_REVIEW_SUMMARY_TIME="${SAMPLE_REVIEW_SUMMARY_TIME}"
export STAGE3_FINALIZE_STUDY_TIME="${FINALIZE_STUDY_TIME}"
export STAGE3_PUBLISH_INTERMEDIATE_PLINK="${PUBLISH_INTERMEDIATE_PLINK}"

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

if command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base)"
  export CONDA_BASE
  export PATH="${CONDA_BASE}/bin:${CONDA_BASE}/condabin:${PATH}"
fi

if [ ! -d "${PIPELINE_DIR}" ]; then
  echo "ERROR: Stage-3 pipeline directory not found: ${PIPELINE_DIR}" >&2
  exit 1
fi

if [ ! -f "${PARAMS_FILE}" ]; then
  echo "ERROR: Missing params file: ${PARAMS_FILE}" >&2
  exit 1
fi

if [ ! -f "${SUMMARY_SCRIPT}" ]; then
  echo "ERROR: Missing stage-3 summary script: ${SUMMARY_SCRIPT}" >&2
  exit 1
fi

if [ ! -d "${STAGE1_ROOT}" ]; then
  echo "ERROR: Stage-1 root not found: ${STAGE1_ROOT}" >&2
  exit 1
fi

if [ ! -d "${STAGE2_ROOT}" ]; then
  echo "ERROR: Stage-2 root not found: ${STAGE2_ROOT}" >&2
  exit 1
fi

if [ ! -f "${DBSNP_VCF}" ]; then
  echo "ERROR: dbSNP VCF not found: ${DBSNP_VCF}" >&2
  exit 1
fi

if [ ! -f "${DBSNP_TBI}" ]; then
  echo "ERROR: dbSNP VCF index not found: ${DBSNP_TBI}" >&2
  exit 1
fi

if [ "${STUDY}" != "all" ]; then
  IFS=',' read -r -a SELECTED_STUDIES <<< "${STUDY}"
  for study_id in "${SELECTED_STUDIES[@]}"; do
    study_id="$(echo "${study_id}" | xargs)"
    if [ ! -d "${STAGE2_ROOT}/${study_id}/stage2" ]; then
      echo "ERROR: Missing stage-2 directory for ${study_id}: ${STAGE2_ROOT}/${study_id}/stage2" >&2
      exit 1
    fi
  done
fi

mkdir -p "${PROJ_ROOT}/src/logs" "${WORKDIR}" "${CONDA_CACHE_DIR}" "${PIPELINE_DIR}/.nextflow"
cd "${PIPELINE_DIR}"

echo "=========================================="
echo " Running EPIC Genetics Stage-3 Pipeline"
echo "=========================================="
echo "Slurm job:                  ${SLURM_JOB_ID:-interactive}"
echo "Host:                       $(hostname)"
echo "Project root:               ${PROJ_ROOT}"
echo "Pipeline:                   ${PIPELINE_DIR}"
echo "Profile:                    ${PROFILE}"
echo "Study:                      ${STUDY}"
echo "Stage1 root:                ${STAGE1_ROOT}"
echo "Stage2 root:                ${STAGE2_ROOT}"
echo "Output root:                ${OUTDIR}"
echo "dbSNP VCF:                  ${DBSNP_VCF}"
echo "dbSNP TBI:                  ${DBSNP_TBI}"
echo "Work dir:                   ${WORKDIR}"
echo "Conda cache:                ${CONDA_CACHE_DIR}"
echo "PREP_DBSNP_CHROM wall time:      ${PREP_DBSNP_TIME}"
echo "FILTER_CHROM wall time:          ${FILTER_CHROM_TIME}"
echo "ANNOTATE_CHROM wall time:        ${ANNOTATE_CHROM_TIME}"
echo "IMPORT_CHROM wall time:          ${IMPORT_CHROM_TIME}"
echo "HWE_CHROM wall time:             ${HWE_CHROM_TIME}"
echo "MERGE_STUDY wall time:           ${MERGE_STUDY_TIME}"
echo "SEX_CHECK wall time:             ${SEX_CHECK_TIME}"
echo "PRUNE_AUTOSOMES wall time:       ${PRUNE_AUTOSOMES_TIME}"
echo "KING_QC wall time:               ${KING_QC_TIME}"
echo "HET_PCA_QC wall time:            ${HET_PCA_QC_TIME}"
echo "SAMPLE_REVIEW_SUMMARY wall time: ${SAMPLE_REVIEW_SUMMARY_TIME}"
echo "FINALIZE_STUDY wall time:        ${FINALIZE_STUDY_TIME}"
echo "Publish PLINK intermediates:${STAGE3_PUBLISH_INTERMEDIATE_PLINK}"
echo "Min R2:                     ${MIN_R2}"
echo "MAF:                        ${MAF}"
echo "Run HWE:                    ${RUN_HWE}"
echo "HWE p-value:                ${HWE_P}"
echo "Exclude ancestry outliers:  ${EXCLUDE_ANCESTRY_OUTLIERS}"
echo "Cache mode:                 ${CACHE_MODE}"
echo "=========================================="

nextflow_cmd=(
  nextflow
  -log "${NEXTFLOW_LOG}"
  run .
  -params-file "${PARAMS_FILE}"
  -profile "${PROFILE}"
  -work-dir "${WORKDIR}"
  --outdir "${OUTDIR}"
  --stage1_root "${STAGE1_ROOT}"
  --stage2_root "${STAGE2_ROOT}"
  --study "${STUDY}"
  --dbsnp_vcf "${DBSNP_VCF}"
  --dbsnp_tbi "${DBSNP_TBI}"
  --slurm_partition "${SLURM_PARTITION}"
  --cache_mode "${CACHE_MODE}"
  --prep_dbsnp_time "${PREP_DBSNP_TIME}"
  --filter_chrom_time "${FILTER_CHROM_TIME}"
  --annotate_chrom_time "${ANNOTATE_CHROM_TIME}"
  --import_chrom_time "${IMPORT_CHROM_TIME}"
  --hwe_chrom_time "${HWE_CHROM_TIME}"
  --merge_study_time "${MERGE_STUDY_TIME}"
  --sex_check_time "${SEX_CHECK_TIME}"
  --prune_autosomes_time "${PRUNE_AUTOSOMES_TIME}"
  --king_qc_time "${KING_QC_TIME}"
  --het_pca_qc_time "${HET_PCA_QC_TIME}"
  --sample_review_summary_time "${SAMPLE_REVIEW_SUMMARY_TIME}"
  --finalize_study_time "${FINALIZE_STUDY_TIME}"
  --publish_intermediate_plink "${PUBLISH_INTERMEDIATE_PLINK}"
  --min_r2 "${MIN_R2}"
  --maf "${MAF}"
  --run_hwe "${RUN_HWE}"
  --hwe_p "${HWE_P}"
  --hwe_k "${HWE_K}"
  --king_cutoff "${KING_CUTOFF}"
  --ancestry_pc_count "${ANCESTRY_PC_COUNT}"
  --ancestry_z_threshold "${ANCESTRY_Z_THRESHOLD}"
  --het_sd_threshold "${HET_SD_THRESHOLD}"
  --exclude_ancestry_outliers "${EXCLUDE_ANCESTRY_OUTLIERS}"
)

if [ "${RESUME}" = "1" ]; then
  nextflow_cmd+=(-resume)
fi

if [ "${#EXTRA_ARGS[@]}" -gt 0 ]; then
  nextflow_cmd+=("${EXTRA_ARGS[@]}")
fi

"${nextflow_cmd[@]}"

echo ""
echo "Building stage-3 summary..."
"${PYTHON3_BIN}" "${SUMMARY_SCRIPT}" \
  --analysis-root "${OUTDIR}" \
  --stage1-root "${STAGE1_ROOT}" \
  --stage2-root "${STAGE2_ROOT}" \
  --output "${SUMMARY_OUTPUT}"

echo ""
echo "Generating stage-3 tables, figures, flags, and report bundles..."
"${PYTHON3_BIN}" "${PIPELINE_DIR}/bin/run_stage3_reports.py" \
  --analysis-root "${OUTDIR}" \
  --studies "${STUDY}"

end_time=$(date +%s)
elapsed=$((end_time - start_time))

echo ""
echo "Stage-3 pipeline complete."
echo "Elapsed: ${elapsed} seconds"
echo "Summary: ${SUMMARY_OUTPUT}"
