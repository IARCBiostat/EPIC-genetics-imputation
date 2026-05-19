#!/bin/bash
#SBATCH --job-name=epic_finalise
#SBATCH --output=src/logs/008_finalise.out
#SBATCH --error=src/logs/008_finalise.err
#SBATCH --time=2:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --partition=low_p

set -euo pipefail
trap 'echo "ERROR: Finalise step failed on line $LINENO" >&2; exit 1' ERR

SCRIPT_PATH="$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")" && pwd)/$(basename -- "${BASH_SOURCE[0]:-$0}")"

print_help() {
  cat <<'EOF'
Usage: bash src/008_finalise.sh [options]

Collects final PLINK2 pfiles and report artifacts into a date-labelled
output directory. When run outside Slurm, submits itself as a batch job.

Options:
  --study, --studies <list>    Comma-separated study IDs (default: all)
  --analysis-root <dir>        Analysis root directory (default: <project>/analysis)
  --outdir <dir>               Destination directory (default: <project>/final_YYYY-MM-DD)
  --dry-run                    Print what would be copied without copying
  --partition <name>           Slurm partition (default: low_p)
  --time <duration>            Wall time for the job (default: 2:00:00)
  --mem <size>                 Memory for the job (default: 4G)
  --cpus <n>                   CPUs for the job (default: 1)
  -h, --help                   Show this help
EOF
}

ENV_FILE="$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")/.." && pwd)/.env"
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

  echo "Submitting finalise step to Slurm..."
  sbatch \
    --export=ALL \
    --job-name "${FINALISE_JOB_NAME:-epic_finalise}" \
    --output "${FINALISE_LOG_OUT:-${PROJ_ROOT}/src/logs/008_finalise_%j.out}" \
    --error "${FINALISE_LOG_ERR:-${PROJ_ROOT}/src/logs/008_finalise_%j.err}" \
    --time "${FINALISE_TIME:-2:00:00}" \
    --mem "${FINALISE_MEM:-4G}" \
    --cpus-per-task "${FINALISE_CPUS:-1}" \
    --partition "${FINALISE_PARTITION:-low_p}" \
    "${SCRIPT_PATH}" "$@"
  exit 0
fi

start_time=$(date +%s)

ANALYSIS_ROOT="${FINALISE_ANALYSIS_ROOT:-${PROJ_ROOT}/analysis}"
DATE_LABEL="${FINALISE_DATE:-$(date +%Y-%m-%d)}"
DEST_ROOT="${FINALISE_OUTDIR:-${PROJ_ROOT}/final_${DATE_LABEL}}"
STUDIES="${FINALISE_STUDIES:-all}"
DRY_RUN=0

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
    --outdir)
      DEST_ROOT="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
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

if ! command -v rsync >/dev/null 2>&1; then
  echo "ERROR: rsync is required but was not found in PATH." >&2
  exit 1
fi

if [ ! -d "${ANALYSIS_ROOT}" ]; then
  echo "ERROR: Analysis root does not exist: ${ANALYSIS_ROOT}" >&2
  exit 1
fi

mkdir -p "${DEST_ROOT}"
MANIFEST="${DEST_ROOT}/finalise_manifest.tsv"
printf "study\tartifact_type\tsource\tdestination\tstatus\n" > "${MANIFEST}"

rsync_cmd=(rsync -a)
if [ "${DRY_RUN}" = "1" ]; then
  rsync_cmd+=(--dry-run)
fi

log_manifest() {
  local study="$1"
  local artifact_type="$2"
  local source="$3"
  local destination="$4"
  local status="$5"
  printf "%s\t%s\t%s\t%s\t%s\n" "${study}" "${artifact_type}" "${source}" "${destination}" "${status}" >> "${MANIFEST}"
}

sync_dir() {
  local study="$1"
  local artifact_type="$2"
  local source="$3"
  local destination="$4"

  if [ -d "${source}" ]; then
    mkdir -p "${destination}"
    "${rsync_cmd[@]}" "${source}/" "${destination}/"
    log_manifest "${study}" "${artifact_type}" "${source}" "${destination}" "copied"
  else
    log_manifest "${study}" "${artifact_type}" "${source}" "${destination}" "missing"
  fi
}

sync_file() {
  local study="$1"
  local artifact_type="$2"
  local source="$3"
  local destination_dir="$4"

  if [ -f "${source}" ]; then
    mkdir -p "${destination_dir}"
    "${rsync_cmd[@]}" "${source}" "${destination_dir}/"
    log_manifest "${study}" "${artifact_type}" "${source}" "${destination_dir}/$(basename "${source}")" "copied"
  else
    log_manifest "${study}" "${artifact_type}" "${source}" "${destination_dir}/$(basename "${source}")" "missing"
  fi
}

discover_studies() {
  if [ "${STUDIES}" != "all" ]; then
    printf "%s\n" "${STUDIES}" | tr ',' '\n' | awk '{$1=$1; print}' | sed '/^$/d' | sort -u
    return
  fi

  find "${ANALYSIS_ROOT}" -mindepth 1 -maxdepth 1 -type d \
    ! -name cohort \
    ! -name deprecated \
    ! -name report \
    ! -name reports \
    ! -name 'stage*-pipeline_info' \
    -exec test -d '{}/stage3/final' ';' -print \
    | awk -F/ '{print $NF}' \
    | sort -u
}

mapfile -t STUDY_LIST < <(discover_studies)
if [ "${#STUDY_LIST[@]}" -eq 0 ]; then
  echo "ERROR: No studies with Stage 3 final outputs were found under ${ANALYSIS_ROOT}." >&2
  exit 1
fi

echo "============================================"
echo " Running EPIC Genetics Finalise"
echo "============================================"
echo "Slurm job:     ${SLURM_JOB_ID:-interactive}"
echo "Host:          $(hostname)"
echo "Project root:  ${PROJ_ROOT}"
echo "Analysis root: ${ANALYSIS_ROOT}"
echo "Destination:   ${DEST_ROOT}"
echo "Studies:       ${#STUDY_LIST[@]} found"
echo "Dry run:       ${DRY_RUN}"
echo "============================================"
echo ""

for study in "${STUDY_LIST[@]}"; do
  study_root="${ANALYSIS_ROOT}/${study}"
  final_dir="${study_root}/stage3/final"

  if [ ! -f "${final_dir}/${study}.pgen" ] || \
     [ ! -f "${final_dir}/${study}.pvar" ] || \
     [ ! -f "${final_dir}/${study}.psam" ]; then
    echo "WARN: Missing final pfiles for ${study} — skipping pfile copy"
    log_manifest "${study}" "stage3_final_pfiles" "${final_dir}/${study}.{pgen,pvar,psam}" "${DEST_ROOT}/studies/${study}/stage3/final" "missing"
    continue
  fi

  echo "  ${study}"
  sync_dir "${study}" "stage3_final_pfiles" "${final_dir}" "${DEST_ROOT}/studies/${study}/stage3/final"

  for stage in stage1 stage2 stage3; do
    sync_dir "${study}" "${stage}_report" "${study_root}/${stage}/report" "${DEST_ROOT}/studies/${study}/${stage}/report"
  done

  sync_file "${study}" "master_report" "${ANALYSIS_ROOT}/report/${study}.master-report.html" "${DEST_ROOT}/reports/master"
  sync_file "${study}" "stage2_report_copy" "${ANALYSIS_ROOT}/report/stage2/${study}.report-stage2.html" "${DEST_ROOT}/reports/stage2"
  sync_file "${study}" "stage3_report_copy" "${ANALYSIS_ROOT}/report/stage3/${study}.report-stage3.html" "${DEST_ROOT}/reports/stage3"
done

for summary in stage1-summary.md stage2-summary.md stage3-summary.md; do
  sync_file "all" "summary" "${ANALYSIS_ROOT}/${summary}" "${DEST_ROOT}/summaries"
done

echo ""
echo "============================================"
echo " Finalise completed"
echo " Destination: ${DEST_ROOT}"
echo " Manifest:    ${MANIFEST}"
echo "============================================"

end_time=$(date +%s)
elapsed=$((end_time - start_time))
echo "Total runtime: $(date -u -r "$elapsed" +%H:%M:%S 2>/dev/null || date -u -d "@$elapsed" +%H:%M:%S)"
