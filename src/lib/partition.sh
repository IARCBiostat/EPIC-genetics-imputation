# src/lib/partition.sh — sourced by the src/ job scripts right after they load .env.
#
# PARTITION in .env is the one Slurm partition for every job and every Nextflow task.
#
# on_partition <script> [sbatch options...] [-- script arguments...]
#   Makes the calling script run as a Slurm job on ${PARTITION}:
#   - started with `bash src/<script>`: submits it with sbatch on ${PARTITION}, exits;
#   - started with `sbatch src/<script>`: Slurm has placed it on the partition in the
#     script's #SBATCH header; if that is not ${PARTITION}, resubmits it once on
#     ${PARTITION} and exits;
#   - already running on ${PARTITION}: returns, and the script carries on.
#   <script> is the script's path relative to the repository root (a job's own copy
#   runs from Slurm's spool directory, so it cannot find itself).

PARTITION="${PARTITION:-low_p}"
REPO_DIR="$(cd "$(dirname "${ENV_FILE:?load .env before sourcing partition.sh}")" && pwd)"

on_partition() {
    local script="$1"
    shift
    local opts=() args=()
    while [ $# -gt 0 ]; do
        if [ "$1" = "--" ]; then
            shift
            args=("$@")
            break
        fi
        opts+=("$1")
        shift
    done

    if [ -n "${SLURM_JOB_ID:-}" ]; then
        if [ "${SLURM_JOB_PARTITION:-}" = "${PARTITION}" ] || [ -n "${EPIC_ON_PARTITION:-}" ]; then
            return 0
        fi
        echo "Job ${SLURM_JOB_ID} started on partition '${SLURM_JOB_PARTITION:-unknown}';" \
             "resubmitting on '${PARTITION}' (PARTITION in .env)."
    fi
    if ! command -v sbatch >/dev/null 2>&1; then
        echo "ERROR: sbatch is not available. Run this on the HPC login node." >&2
        exit 1
    fi
    mkdir -p "${REPO_DIR}/src/logs"
    cd "${REPO_DIR}"
    # EPIC_ON_PARTITION marks the new job as placed, so it runs instead of moving again.
    EPIC_ON_PARTITION=1 sbatch --export=ALL --partition "${PARTITION}" \
        ${opts[@]+"${opts[@]}"} "${REPO_DIR}/${script}" ${args[@]+"${args[@]}"}
    exit 0
}
