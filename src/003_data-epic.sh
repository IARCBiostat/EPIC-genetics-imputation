#!/bin/bash
#SBATCH --job-name=003_data_epic
#SBATCH --output=src/logs/003_data_epic.out
#SBATCH --error=src/logs/003_data_epic.err
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=low_p

set -euo pipefail
trap 'echo "ERROR: Job failed on line $LINENO" >&2; exit 1' ERR
start_time=$(date +%s)

# ── Environment ────────────────────────────────────────────────────────────────
ENV_FILE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")/.." && pwd)}/.env"
if [ ! -f "$ENV_FILE" ]; then
  echo "ERROR: .env not found at ${ENV_FILE}" >&2; exit 1
fi
# shellcheck disable=SC1090
set -a; source "$ENV_FILE"; set +a

PROJ_ROOT="${SLURM_SUBMIT_DIR:-$(cd "$(dirname -- "${BASH_SOURCE[0]:-$0}")/.." && pwd)}"

# Activate the nf_EPIC-genetics conda environment (provides r-base + r-haven).
if command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base)"
  # shellcheck disable=SC1091
  source "${CONDA_BASE}/etc/profile.d/conda.sh" || true
  conda activate nf_EPIC-genetics || true
fi

RSCRIPT="${RSCRIPT_BIN:-$(command -v Rscript)}"

echo "=========================================="
echo " Building EPIC Study Case Status File"
echo " Host: $(hostname)"
echo "=========================================="

"${RSCRIPT}" "${PROJ_ROOT}/src/misc/003_data-epic.R" "$@"

end_time=$(date +%s)
elapsed=$(( end_time - start_time ))
echo "=========================================="
echo " Done. Time taken: $((elapsed / 60))m $((elapsed % 60))s"
echo "=========================================="
