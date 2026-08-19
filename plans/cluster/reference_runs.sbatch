#!/bin/bash
#SBATCH --job-name=pp_ref
#SBATCH --array=0-5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_logs/pp_ref_%A_%a.out
#SBATCH --error=slurm_logs/pp_ref_%A_%a.err
# ---- EDIT FOR YOUR ALLOCATION -----------------------------------------
##SBATCH --account=CHANGEME
##SBATCH --partition=CHANGEME
# -----------------------------------------------------------------------
#
# Multi-seed reference sampler runs, one array task per (config, seed).
#
# Built for r7's post-build step -- "two nested-sampling runs (ocean
# 5-param, frozen 3-param), n_eff=8000 n_active=4096, 3 seeds each, no
# branch coordinate" -- i.e. 2 configs x 3 seeds = 6 tasks, the default
# --array=0-5. Also serves the B3 reference-wander protocol (1 config x
# >=3 seeds: set CONFIGS to one path and --array=0-2).
#
# Task index maps as: config = CONFIGS[idx / n_seeds], seed = SEEDS[idx % n_seeds]
#
# Submit (set CONFIGS to the ratified branch configs once they exist):
#   CONFIGS="path/to/ocean_5param.json path/to/frozen_3param.json" \
#   SEEDS="8001 8002 8003" \
#     sbatch plans/cluster/reference_runs.sbatch
#
# The seed reaches the sampler through a scratch-only patched copy of the
# committed config (cluster_seeded_config.py) because run_inference_cli
# takes its seed from config.random_state, not a flag. The committed
# config stays the source of truth; the patch is recorded inside each
# derived config's metadata so a result pickle traces back to its parent.
#
# --cpus-per-task=16: nested samplers parallelize their likelihood calls,
# and unlike the cache build these ARE worth threading. Adjust to your
# node size. NOTE pp_env.sh pins BLAS to 1 thread on purpose -- the
# sampler's own pool provides the parallelism; nested BLAS threading on
# top of it oversubscribes.
#
# --time=24:00:00 is a guess scaled from the n_eff=2000 reference runs
# this project has already done; n_eff=8000 with n_active=4096 is
# substantially heavier. CHECK the first task's rate early and resubmit
# with a longer wall if it is going to clip -- a sampler killed at the
# wall leaves a partial chain, which must NOT be pooled as if converged.

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/pp_env.sh"

CONFIGS="${CONFIGS:?set CONFIGS to a space-separated list of config JSON paths}"
SEEDS="${SEEDS:-8001 8002 8003}"
OUT_DIR="${OUT_DIR:-$PP_REPO/validation_reports/cluster_reference_runs}"

read -ra CFG_ARR <<< "$CONFIGS"
read -ra SEED_ARR <<< "$SEEDS"
N_SEEDS=${#SEED_ARR[@]}
N_CFG=${#CFG_ARR[@]}
N_TASKS=$(( N_CFG * N_SEEDS ))

IDX="${SLURM_ARRAY_TASK_ID:-0}"
if (( IDX >= N_TASKS )); then
  echo "ERROR: array index $IDX exceeds ${N_CFG} configs x ${N_SEEDS} " \
       "seeds = $N_TASKS tasks. Fix --array to 0-$(( N_TASKS - 1 ))." >&2
  exit 1
fi

CFG="${CFG_ARR[$(( IDX / N_SEEDS ))]}"
SEED="${SEED_ARR[$(( IDX % N_SEEDS ))]}"
TAG="$(basename "$CFG" .json)_seed${SEED}"

echo "task $IDX -> config=$CFG seed=$SEED"

PATCHED="$PP_SCRATCH/${TAG}.json"
python plans/scripts/cluster_seeded_config.py \
  --config "$CFG" --seed "$SEED" --out "$PATCHED" \
  ${EXTRA_SET:+$EXTRA_SET}

mkdir -p "$OUT_DIR"
# Stage the result in scratch, then move -- see the merge job's note about
# writing large pickles straight into a synced repo path.
STAGED="$PP_SCRATCH/${TAG}_result.pkl"

srun python -m PlanetProfile.Inference.run_inference_cli \
  --config "$PATCHED" \
  --output "$STAGED" \
  --progress "$PP_SCRATCH/${TAG}_progress.json" \
  --verbose

mv "$STAGED" "$OUT_DIR/${TAG}_result.pkl"
cp "$PATCHED" "$OUT_DIR/${TAG}_config.json"   # provenance beside the result
echo "done -> $OUT_DIR/${TAG}_result.pkl"
