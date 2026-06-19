#!/bin/bash
#SBATCH --job-name=copy_final
#SBATCH --output=src/logs/copy.out
#SBATCH --error=src/logs/copy.err
#SBATCH --time=8:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=25
#SBATCH --partition=low_p

# Move the finalised study outputs (final/) from scratch to the project root on
# /data, naming the destination directory after the original run date. The copy
# is parallelised across the top-level entries of final/ (one rsync per entry,
# up to N in parallel), then verified by file count before the scratch source is
# deleted — a safe move.
#
# Usage:
#   sbatch src/copy.sh                  # copy + verify only (leave scratch intact)
#   sbatch src/copy.sh --delete-source  # copy + verify + delete scratch source

set -euo pipefail
trap 'echo "ERROR: copy.sh failed on line $LINENO" >&2; exit 1' ERR

# ── Load environment ───────────────────────────────────────────────────────────
ENV_FILE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")/.." && pwd)}/.env"
if [ ! -f "$ENV_FILE" ]; then
  echo "ERROR: .env not found at ${ENV_FILE}" >&2; exit 1
fi
set -a; source "$ENV_FILE"; set +a

DELETE_SOURCE=0
[ "${1:-}" = "--delete-source" ] && DELETE_SOURCE=1

# Original run date: prefer SCRATCH_DATE, else the basename of SCRATCH_RUN.
RUN_DATE="${SCRATCH_DATE:-$(basename "$SCRATCH_RUN")}"
SRC="${SCRATCH_RUN}/final"
DEST="${GENETICS_PROJECT_ROOT}/${RUN_DATE}"

# One parallel rsync per top-level entry of final/ — i.e. one job per study
# (plus report/ and summaries/). Counted dynamically so it scales to the run.
if [ -d "$SRC" ]; then
  NPROC=$(find "$SRC" -mindepth 1 -maxdepth 1 | wc -l)
fi
[ "${NPROC:-0}" -lt 1 ] && NPROC=1

echo "=========================================="
echo " EPIC Genetics — Move final/ to /data"
echo "=========================================="
echo "Slurm job:   ${SLURM_JOB_ID:-interactive}"
echo "Host:        $(hostname)"
echo "Source:      ${SRC}"
echo "Dest:        ${DEST}"
echo "Parallelism: ${NPROC}"
echo "Delete src:  ${DELETE_SOURCE}"
echo "=========================================="

if [ ! -d "$SRC" ]; then
  echo "ERROR: source directory not found: ${SRC}" >&2; exit 1
fi
mkdir -p "$DEST"

# ── Parallel copy ──────────────────────────────────────────────────────────────
# One rsync per top-level entry of final/, up to NPROC concurrently. Each entry
# (dir or file) is copied into DEST preserving its name. xargs exits non-zero if
# any rsync fails, which trips the ERR trap above.
#
# --no-group/--no-owner: do NOT preserve the scratch group/owner. /data enforces
#   its own setgid group, so chgrp would fail with "Permission denied"; instead
#   new files inherit the destination directory's group.
# --delete: remove anything in the dest entry not present in the source, so
#   orphaned partial temp files from an interrupted run don't linger (and don't
#   inflate the verification file count).
echo ""
echo "Copying ${SRC}/ -> ${DEST}/ with ${NPROC} parallel rsync workers ..."
find "$SRC" -mindepth 1 -maxdepth 1 -print0 \
  | xargs -0 -P "$NPROC" -I{} rsync -rlptD --no-group --no-owner --delete "{}" "$DEST/"

# ── Verify (file counts must match) ────────────────────────────────────────────
echo ""
echo "Verifying ..."
SRC_COUNT=$(find "$SRC"  -type f | wc -l)
DST_COUNT=$(find "$DEST" -type f | wc -l)
echo "  Source files: ${SRC_COUNT}"
echo "  Dest files:   ${DST_COUNT}"
echo "  Source size:  $(du -sh "$SRC"  | cut -f1)"
echo "  Dest size:    $(du -sh "$DEST" | cut -f1)"

if [ "$SRC_COUNT" -ne "$DST_COUNT" ]; then
  echo "ERROR: file counts differ — leaving scratch source intact." >&2
  exit 1
fi
echo "  File counts match."

# ── Remove scratch source (only with --delete-source) ──────────────────────────
if [ "$DELETE_SOURCE" -eq 1 ]; then
  echo ""
  echo "Removing scratch source: ${SRC}"
  rm -rf "$SRC"
else
  echo ""
  echo "Leaving scratch source in place: ${SRC}"
  echo "(pass --delete-source to remove it after a verified copy)"
fi

echo ""
echo "=========================================="
echo " Done. final/ now at ${DEST}"
echo "=========================================="
