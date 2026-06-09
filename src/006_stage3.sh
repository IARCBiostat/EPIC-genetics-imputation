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

SCRIPT_PATH="$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")" && pwd)/$(basename -- "${BASH_SOURCE[0]:-$0}")"

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

QC thresholds default to values in pipeline_stage3/params.yaml.
Pass a flag only to override the params.yaml value for that run.

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
  --prepare-chrom-time <duration>         Wall time for PREPARE_CHROM jobs (default: 72h)
  --import-chrom-time <duration>          Wall time for IMPORT_CHROM jobs (default: 72h)
  --hwe-chrom-time <duration>             Wall time for HWE_CHROM jobs (default: 72h)
  --merge-study-time <duration>           Wall time for MERGE_FOR_QC jobs (default: 72h)
  --prune-autosomes-time <duration>       Wall time for PRUNE_AUTOSOMES jobs (default: 72h)
  --king-qc-time <duration>               Wall time for KING_QC jobs (default: 72h)
  --het-pca-qc-time <duration>            Wall time for HET_PCA_QC jobs (default: 72h)
  --sample-review-summary-time <duration> Wall time for SAMPLE_REVIEW_SUMMARY jobs (default: 72h)
  --finalize-chrom-time <duration>        Wall time for FINALIZE_CHROM jobs (default: 72h)
  --publish-intermediate-plink      Copy intermediate PLINK/PGEN/BED files into analysis output
  --no-publish-intermediate-plink   Do not copy intermediate PLINK/PGEN/BED files (default)
  --min-r2 <value>                  Minimum imputation R2 (overrides params.yaml)
  --maf <value>                     Minimum MAF (overrides params.yaml)
  --hwe true|false                  Apply HWE filtering on controls only (overrides params.yaml; default: true)
  --hwe-p <value>                   HWE p-value threshold (overrides params.yaml)
  --related true|false              Exclude related samples from the final dataset (overrides params.yaml; default: false)
  --ancestry true|false             Exclude ancestry outliers from the final dataset (overrides params.yaml; default: false)
  --no-resume                       Disable Nextflow -resume
  -h, --help                        Show this help
EOF
}

ENV_FILE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")/.." && pwd)}/.env"
if [ ! -f "$ENV_FILE" ]; then
  echo "ERROR: .env not found at ${ENV_FILE}" >&2; exit 1
fi
# shellcheck disable=SC1090
set -a; source "$ENV_FILE"; set +a
PROJ_ROOT="${GENETICS_PROJECT_ROOT}"

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

# --- Infrastructure settings (environment-specific; shell defaults are appropriate) ---
PROFILE="${STAGE3_PROFILE:-slurm}"
STUDY="${STAGE3_STUDY:-all}"
OUTDIR="${STAGE3_OUTDIR:-${SCRATCH_RUN}/studies}"
STAGE1_ROOT="${STAGE3_STAGE1_ROOT:-${SCRATCH_RUN}/studies}"
STAGE2_ROOT="${STAGE3_STAGE2_ROOT:-${SCRATCH_RUN}/studies}"
SLURM_PARTITION="${STAGE3_PARTITION:-${SLURM_JOB_PARTITION:-low_p}}"
PREP_DBSNP_TIME="${STAGE3_PREP_DBSNP_TIME:-72h}"
PREPARE_CHROM_TIME="${STAGE3_PREPARE_CHROM_TIME:-72h}"
IMPORT_CHROM_TIME="${STAGE3_IMPORT_CHROM_TIME:-72h}"
HWE_CHROM_TIME="${STAGE3_HWE_CHROM_TIME:-72h}"
MERGE_STUDY_TIME="${STAGE3_MERGE_STUDY_TIME:-72h}"
PRUNE_AUTOSOMES_TIME="${STAGE3_PRUNE_AUTOSOMES_TIME:-72h}"
KING_QC_TIME="${STAGE3_KING_QC_TIME:-72h}"
HET_PCA_QC_TIME="${STAGE3_HET_PCA_QC_TIME:-72h}"
SAMPLE_REVIEW_SUMMARY_TIME="${STAGE3_SAMPLE_REVIEW_SUMMARY_TIME:-72h}"
FINALIZE_CHROM_TIME="${STAGE3_FINALIZE_CHROM_TIME:-72h}"
PUBLISH_INTERMEDIATE_PLINK="${STAGE3_PUBLISH_INTERMEDIATE_PLINK:-false}"
WORKDIR="${STAGE3_WORKDIR:-${SCRATCH_RUN}/stage3/work}"
CONDA_CACHE_DIR="${STAGE3_CONDA_CACHE_DIR:-${SCRATCH_RUN}/stage3/conda}"
CONDA_SOLVER="${STAGE3_CONDA_SOLVER:-classic}"
CONDA_CHANNEL_PRIORITY="${STAGE3_CONDA_CHANNEL_PRIORITY:-strict}"
CACHE_MODE="${STAGE3_CACHE_MODE:-lenient}"
RESUME=1

# --- QC parameters: empty means "use params.yaml default"; only set if user passes a flag ---
HWE=""
HWE_P=""
MIN_R2=""
MAF=""
RELATED=""
ANCESTRY=""

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
    --prepare-chrom-time)
      PREPARE_CHROM_TIME="$2"
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
    --finalize-chrom-time)
      FINALIZE_CHROM_TIME="$2"
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
    --hwe)
      if [[ "$2" != "true" && "$2" != "false" ]]; then
        echo "ERROR: --hwe must be 'true' or 'false'" >&2; exit 1
      fi
      HWE="$2"
      shift 2
      ;;
    --hwe-p)
      HWE_P="$2"
      shift 2
      ;;
    --related)
      if [[ "$2" != "true" && "$2" != "false" ]]; then
        echo "ERROR: --related must be 'true' or 'false'" >&2; exit 1
      fi
      RELATED="$2"
      shift 2
      ;;
    --ancestry)
      if [[ "$2" != "true" && "$2" != "false" ]]; then
        echo "ERROR: --ancestry must be 'true' or 'false'" >&2; exit 1
      fi
      ANCESTRY="$2"
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

PIPELINE_DIR="${PROJ_ROOT}/pipeline_stage3"
PARAMS_FILE="${PIPELINE_DIR}/params.yaml"
SUMMARY_SCRIPT="${PIPELINE_DIR}/bin/summary.py"
SUMMARY_OUTPUT="${OUTDIR}/stage3-summary.md"
NEXTFLOW_LOG="${SCRATCH_RUN}/stage3/.nextflow.log"

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
export STAGE3_PREPARE_CHROM_TIME="${PREPARE_CHROM_TIME}"
export STAGE3_IMPORT_CHROM_TIME="${IMPORT_CHROM_TIME}"
export STAGE3_HWE_CHROM_TIME="${HWE_CHROM_TIME}"
export STAGE3_MERGE_STUDY_TIME="${MERGE_STUDY_TIME}"
export STAGE3_PRUNE_AUTOSOMES_TIME="${PRUNE_AUTOSOMES_TIME}"
export STAGE3_KING_QC_TIME="${KING_QC_TIME}"
export STAGE3_HET_PCA_QC_TIME="${HET_PCA_QC_TIME}"
export STAGE3_SAMPLE_REVIEW_SUMMARY_TIME="${SAMPLE_REVIEW_SUMMARY_TIME}"
export STAGE3_FINALIZE_CHROM_TIME="${FINALIZE_CHROM_TIME}"
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

mkdir -p "${PROJ_ROOT}/src/logs" "${WORKDIR}" "${CONDA_CACHE_DIR}" "${SCRATCH_RUN}/stage3"
cd "${SCRATCH_RUN}/stage3"

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
echo "PREPARE_CHROM wall time:         ${PREPARE_CHROM_TIME}"
echo "IMPORT_CHROM wall time:          ${IMPORT_CHROM_TIME}"
echo "HWE_CHROM wall time:             ${HWE_CHROM_TIME}"
echo "MERGE_FOR_QC wall time:          ${MERGE_STUDY_TIME}"
echo "PRUNE_AUTOSOMES wall time:       ${PRUNE_AUTOSOMES_TIME}"
echo "KING_QC wall time:               ${KING_QC_TIME}"
echo "HET_PCA_QC wall time:            ${HET_PCA_QC_TIME}"
echo "SAMPLE_REVIEW_SUMMARY wall time: ${SAMPLE_REVIEW_SUMMARY_TIME}"
echo "FINALIZE_CHROM wall time:        ${FINALIZE_CHROM_TIME}"
echo "Publish PLINK intermediates: ${PUBLISH_INTERMEDIATE_PLINK}"
echo "Min R2:                     ${MIN_R2:-from params.yaml}"
echo "MAF:                        ${MAF:-from params.yaml}"
echo "HWE:                        ${HWE:-from params.yaml}"
echo "HWE p-value:                ${HWE_P:-from params.yaml}"
echo "Exclude related:            ${RELATED:-from params.yaml}"
echo "Exclude ancestry outliers:  ${ANCESTRY:-from params.yaml}"
echo "Cache mode:                 ${CACHE_MODE}"
echo "=========================================="

nextflow_cmd=(
  nextflow
  -log "${NEXTFLOW_LOG}"
  run "${PIPELINE_DIR}"
  -params-file "${PARAMS_FILE}"
  -profile "${PROFILE}"
  -work-dir "${WORKDIR}"
  --outdir "${OUTDIR}"
  --pipeline_info_dir "${SCRATCH_RUN}/stage3/pipeline_info"
  --stage1_root "${STAGE1_ROOT}"
  --stage2_root "${STAGE2_ROOT}"
  --study "${STUDY}"
  --dbsnp_vcf "${DBSNP_VCF}"
  --dbsnp_tbi "${DBSNP_TBI}"
  --slurm_partition "${SLURM_PARTITION}"
  --cache_mode "${CACHE_MODE}"
  --prep_dbsnp_time "${PREP_DBSNP_TIME}"
  --prepare_chrom_time "${PREPARE_CHROM_TIME}"
  --import_chrom_time "${IMPORT_CHROM_TIME}"
  --hwe_chrom_time "${HWE_CHROM_TIME}"
  --merge_study_time "${MERGE_STUDY_TIME}"
  --prune_autosomes_time "${PRUNE_AUTOSOMES_TIME}"
  --king_qc_time "${KING_QC_TIME}"
  --het_pca_qc_time "${HET_PCA_QC_TIME}"
  --sample_review_summary_time "${SAMPLE_REVIEW_SUMMARY_TIME}"
  --finalize_chrom_time "${FINALIZE_CHROM_TIME}"
  --publish_intermediate_plink "${PUBLISH_INTERMEDIATE_PLINK}"
)

# QC overrides: only passed if explicitly set by the user; otherwise params.yaml governs.
[ -n "${MIN_R2}" ]  && nextflow_cmd+=(--min_r2 "${MIN_R2}")
[ -n "${MAF}" ]     && nextflow_cmd+=(--maf "${MAF}")
[ -n "${HWE}" ]     && nextflow_cmd+=(--hwe "${HWE}")
[ -n "${HWE_P}" ]   && nextflow_cmd+=(--hwe_p "${HWE_P}")
[ -n "${RELATED}" ] && nextflow_cmd+=(--related "${RELATED}")
[ -n "${ANCESTRY}" ] && nextflow_cmd+=(--ancestry "${ANCESTRY}")

if [ "${RESUME}" = "1" ]; then
  nextflow_cmd+=(-resume)
fi

if [ "${#EXTRA_ARGS[@]}" -gt 0 ]; then
  nextflow_cmd+=("${EXTRA_ARGS[@]}")
fi

# Auto-retry on transient NFS IOException.
# Nextflow has a resource-leak bug where FileChannelImpl is closed by the JVM
# garbage-collector rather than explicitly; on NFS/Lustre this occasionally
# returns EINVAL on the deferred close(), which Nextflow treats as a
# session-level abort (killing all 150+ running tasks).  The fix is to detect
# that specific abort, remove the empty/partial work dir that triggered it, and
# resume — all completed tasks are preserved by the lenient cache.
MAX_NF_ATTEMPTS=8
NF_ATTEMPT=0
NF_EXIT=1

# Suspend the ERR trap for the retry loop so that Nextflow failures are handled
# by the NFS IOException retry logic below rather than causing an immediate exit.
trap - ERR

while [ "${NF_ATTEMPT}" -lt "${MAX_NF_ATTEMPTS}" ]; do
    NF_ATTEMPT=$(( NF_ATTEMPT + 1 ))
    [ "${NF_ATTEMPT}" -gt 1 ] && \
        echo "Nextflow resume attempt ${NF_ATTEMPT}/${MAX_NF_ATTEMPTS} after NFS IOException..."

    # Record log offset before this attempt so we only scan lines from the current run.
    LOG_START=$([ -f "${NEXTFLOW_LOG}" ] && wc -l < "${NEXTFLOW_LOG}" || echo 0)
    LOG_START=$(( LOG_START + 1 ))

    set +e
    "${nextflow_cmd[@]}"
    NF_EXIT=$?
    set -e

    [ "${NF_EXIT}" -eq 0 ] && break

    # Only retry for the NFS IOException; propagate genuine pipeline failures.
    # Scan only the lines appended during this attempt to avoid false matches from
    # previous runs whose IOException entries persist in the accumulated log file.
    if ! tail -n +"${LOG_START}" "${NEXTFLOW_LOG}" 2>/dev/null \
            | grep -q "java.io.IOException: Invalid argument"; then
        echo "ERROR: Nextflow failed (exit ${NF_EXIT}) without NFS IOException — not retrying." >&2
        break
    fi

    echo "WARNING: Nextflow session aborted by NFS IOException." >&2
    echo "  Identifying and removing partial work directories..." >&2

    tail -n +"${LOG_START}" "${NEXTFLOW_LOG}" 2>/dev/null | awk '
        /Handling unexpected condition/ { expect_wdir=1; wdir="" }
        expect_wdir && /work-dir=/ {
            n = split($0, a, "work-dir=")
            split(a[2], b, ";")
            wdir = b[1]; expect_wdir=0; expect_exc=1
        }
        expect_exc && /java\.io\.IOException/ { print wdir; expect_exc=0; wdir="" }
    ' | sort -u | \
    while IFS= read -r wdir; do
        [ -n "${wdir}" ] || continue
        { rm -rf "${wdir}" && echo "  Removed: ${wdir}" >&2; } 2>/dev/null || true
    done

    # Ensure -resume is active for all subsequent attempts.
    if ! printf '%s\0' "${nextflow_cmd[@]}" | grep -qz -- '-resume'; then
        nextflow_cmd+=(-resume)
    fi

    sleep 30
done

# Restore ERR trap for the remaining steps.
trap 'echo "ERROR: Stage-3 pipeline failed on line $LINENO" >&2; exit 1' ERR

if [ "${NF_EXIT}" -ne 0 ]; then
    echo "ERROR: Nextflow pipeline failed after ${NF_ATTEMPT} attempt(s) (exit ${NF_EXIT})." >&2
    exit "${NF_EXIT}"
fi

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
