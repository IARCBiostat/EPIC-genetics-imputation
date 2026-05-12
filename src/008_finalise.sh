#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

ANALYSIS_ROOT="${ANALYSIS_ROOT:-${ROOT}/analysis}"
DATE_LABEL="${FINALISE_DATE:-$(date +%Y-%m-%d)}"
DEST_ROOT="${FINALISE_OUTDIR:-${ROOT}/final_${DATE_LABEL}}"
STUDIES="all"
DRY_RUN=0

usage() {
  cat <<'EOF'
Usage: bash src/008_finalise.sh [--study STUDY[,STUDY2]] [--analysis-root PATH] [--outdir PATH] [--dry-run]

Collect final study PLINK2 pfiles and report artifacts into a date-labelled
directory at the project root. By default the destination is:

  final_YYYY-MM-DD/
EOF
}

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
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "ERROR: Unknown argument: $1" >&2
      usage >&2
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
  echo "ERROR: No studies with Stage 3 final outputs were found." >&2
  exit 1
fi

for study in "${STUDY_LIST[@]}"; do
  study_root="${ANALYSIS_ROOT}/${study}"
  final_dir="${study_root}/stage3/final"

  if [ ! -f "${final_dir}/${study}.pgen" ] || [ ! -f "${final_dir}/${study}.pvar" ] || [ ! -f "${final_dir}/${study}.psam" ]; then
    log_manifest "${study}" "stage3_final_pfiles" "${final_dir}/${study}.{pgen,pvar,psam}" "${DEST_ROOT}/studies/${study}/stage3/final" "missing"
    continue
  fi

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

echo "Finalised outputs: ${DEST_ROOT}"
echo "Manifest: ${MANIFEST}"
