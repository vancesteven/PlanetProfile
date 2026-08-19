#!/bin/bash
#SBATCH --job-name=sbc_corr
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=slurm_logs/sbc_corr_%j.out
#SBATCH --error=slurm_logs/sbc_corr_%j.err
# ---- EDIT FOR YOUR ALLOCATION -----------------------------------------
##SBATCH --account=CHANGEME
##SBATCH --partition=CHANGEME
# -----------------------------------------------------------------------
#
# Corrected-pipeline SBC (§0.28 protocol) on one composition.
#
# Single node, many cores: the driver already forks a ProcessPoolExecutor
# over pairs and each pair is fully independent, so --n-jobs
# $SLURM_CPUS_PER_TASK is all that is needed -- no array sharding, no new
# code, and the pair->rank bookkeeping stays in one process where the
# BH-FDR step can see every surviving pair at once.
#
# Measured on the completed NH3 run: 613 pairs, ~6.2 min/pair
# serial-equivalent = ~63 CPU-h, which finished in ~14h wall at 8 workers.
# At 32 cores that is ~2h; --time=12:00:00 leaves room for a slower node
# and for the salts (whose caches are larger).
#
# --mem=64G: the NH3 run held ~1.4GB RSS per worker (N=20000 posterior
# draws plus the structure cache), so 32 workers need ~45GB plus the
# parent. Scale this WITH --cpus-per-task, or workers will be OOM-killed
# mid-sweep and their pairs recorded as errors.
#
# Submit:
#   sbatch plans/cluster/nh3_corrected_sbc.sbatch                  # NH3
#   DRIVER_ARGS="--n-post 20000" sbatch plans/cluster/nh3_corrected_sbc.sbatch
#
# ALWAYS run the fiducial reproduction check first (cheap, ~12 min, and it
# is the preregistered §5 precondition for the real run):
#   VERIFY_ONLY=1 sbatch plans/cluster/nh3_corrected_sbc.sbatch

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/pp_env.sh"

DRIVER="${DRIVER:-plans/scripts/nh3_corrected_sbc.py}"
N_JOBS="${N_JOBS:-${SLURM_CPUS_PER_TASK:-1}}"

if [[ -n "${VERIFY_ONLY:-}" ]]; then
  echo "preregistration §5 check only (fiducial reproduction + rank smoke)"
  srun python "$DRIVER" --verify-only
  exit 0
fi

echo "running with --n-jobs $N_JOBS"
srun python "$DRIVER" --n-jobs "$N_JOBS" ${DRIVER_ARGS:-}

echo
echo "NEXT: the driver writes its own report + verdict. Do NOT reinterpret"
echo "a FAIL or INCONCLUSIVE-BY-COVERAGE locally -- both are preregistered"
echo "outcomes routed to the scientific-reviewer."
