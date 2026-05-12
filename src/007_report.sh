#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

ANALYSIS_ROOT="${ANALYSIS_ROOT:-${ROOT}/analysis}"
STUDIES="all"
PYTHON3_BIN="${PYTHON3_BIN:-python3}"

usage() {
  cat <<USAGE
Usage: bash src/007_report.sh [--study STUDY[,STUDY2]] [--analysis-root PATH]

Creates cross-stage per-study master reports under analysis/report/.

Options:
  --study, --studies    Study ID or comma-separated study IDs. Default: all.
  --analysis-root       Analysis root. Default: ${ROOT}/analysis.
  -h, --help            Show this help.
USAGE
}

while [ "$#" -gt 0 ]; do
  case "$1" in
    --study|--studies)
      STUDIES="$2"
      shift 2
      ;;
    --analysis-root)
      ANALYSIS_ROOT="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [ ! -d "${ANALYSIS_ROOT}" ]; then
  echo "ERROR: Analysis root does not exist: ${ANALYSIS_ROOT}" >&2
  exit 1
fi

if [ ! -f "${ROOT}/src/007_report/master_report.py" ]; then
  echo "ERROR: Missing report builder: ${ROOT}/src/007_report/master_report.py" >&2
  exit 1
fi

"${PYTHON3_BIN}" "${ROOT}/src/007_report/master_report.py" \
  --analysis-root "${ANALYSIS_ROOT}" \
  --studies "${STUDIES}"
