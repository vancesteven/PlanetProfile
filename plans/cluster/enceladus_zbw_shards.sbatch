#!/bin/bash
#SBATCH --job-name=enc_zbw
#SBATCH --array=0-10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=slurm_logs/enc_zbw_%A_%a.out
#SBATCH --error=slurm_logs/enc_zbw_%A_%a.err
# ---- EDIT FOR YOUR ALLOCATION -----------------------------------------
##SBATCH --account=CHANGEME
##SBATCH --partition=CHANGEME
# -----------------------------------------------------------------------
#
# Enceladus (zb, w) production cache -- ocean grid, sharded by zb.
#
# 11 array tasks x 8 zb nodes x 40 w = 3520 ocean nodes, plus the 39-node
# frozen axis on task 0. Measured single-node rate on the reference run was
# ~350 ocean nodes/h, so a 320-node shard is ~1h; the 6h wall request is
# ~6x headroom (per-node cost varies with shell thickness, and the frozen
# axis adds ~35min to task 0). Whole-grid serial equivalent was ~10h.
#
# --cpus-per-task=1 is deliberate: the builder is single-threaded per node
# and pp_env.sh pins BLAS to 1 thread, so parallelism comes from the array,
# not from threads. Asking for more cores would idle them.
#
# --mem=8G: the reference single-node run peaked ~9GB RSS holding the WHOLE
# 3520-node structure list; one 320-node shard holds ~1/11 of that, so 8G
# is generous. Task 0 also holds the frozen axis (39 small nodes).
#
# Submit:
#   mkdir -p slurm_logs
#   sbatch plans/cluster/enceladus_zbw_shards.sbatch
#   # then, after ALL tasks succeed:
#   sbatch --dependency=afterok:<jobid> plans/cluster/enceladus_zbw_merge.sbatch
#
# Or in one shot (merge auto-waits on the array):
#   JID=$(sbatch --parsable plans/cluster/enceladus_zbw_shards.sbatch)
#   sbatch --dependency=afterok:$JID plans/cluster/enceladus_zbw_merge.sbatch
#
# Check the decomposition WITHOUT burning an allocation first:
#   python plans/scripts/enceladus_zbw_shard_build.py --n-shards 11 --plan
#
# IMPORTANT: keep --array and N_SHARDS consistent. The merge step refuses
# to assemble a partial array (it checks index coverage against the
# canonical 88-node axis), so a truncated array fails loudly rather than
# silently producing a cache with a hole in it.

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/pp_env.sh"

N_SHARDS="${N_SHARDS:-11}"
SHARD_DIR="${SHARD_DIR:-$PP_SCRATCH/enc_zbw_shards}"
# Shards must OUTLIVE the array tasks for the merge job to read them, so
# they cannot live in a per-task SLURM_TMPDIR. Default to a shared path.
if [[ -n "${SLURM_TMPDIR:-}" && "$SHARD_DIR" == "$SLURM_TMPDIR"* ]]; then
  SHARD_DIR="${SCRATCH:-$PP_REPO}/enc_zbw_shards_${SLURM_ARRAY_JOB_ID:-manual}"
fi
mkdir -p "$SHARD_DIR"
echo "shard output dir: $SHARD_DIR"
echo "  (pass the SAME SHARD_DIR to the merge job)"

# Exactly ONE task builds the frozen axis; the merge asserts exactly one.
FROZEN_FLAG=""
if [[ "${SLURM_ARRAY_TASK_ID:-0}" == "0" ]]; then
  FROZEN_FLAG="--with-frozen"
  echo "this task ALSO builds the 39-node frozen axis"
fi

srun python plans/scripts/enceladus_zbw_shard_build.py \
  --shard-index "${SLURM_ARRAY_TASK_ID:-0}" \
  --n-shards "$N_SHARDS" \
  --out-dir "$SHARD_DIR" \
  $FROZEN_FLAG

echo "shard ${SLURM_ARRAY_TASK_ID:-0} done"
