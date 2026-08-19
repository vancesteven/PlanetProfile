#!/bin/bash
#SBATCH --job-name=enc_zbw_merge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=slurm_logs/enc_zbw_merge_%j.out
#SBATCH --error=slurm_logs/enc_zbw_merge_%j.err
# ---- EDIT FOR YOUR ALLOCATION -----------------------------------------
##SBATCH --account=CHANGEME
##SBATCH --partition=CHANGEME
# -----------------------------------------------------------------------
#
# Assemble the sharded Enceladus (zb, w) cache and emit the six r7 §0.29
# acceptance gates. Pure I/O + assertions -- no physics is recomputed.
#
# --mem=32G: this is the one step that holds ALL 3520 structures plus the
# 39 frozen nodes in memory at once (the reference single-node build peaked
# ~9GB; 32G leaves room for the shard-by-shard unpickle overlap).
#
# Run AFTER the shard array succeeds, with the SAME SHARD_DIR:
#   SHARD_DIR=/scratch/me/enc_zbw_shards_12345 \
#     sbatch --dependency=afterok:12345 \
#     plans/cluster/enceladus_zbw_merge.sbatch
#
# `afterok` matters: it will not start if ANY array task failed, which is
# what you want -- the merge refuses partial coverage anyway, so this just
# turns a confusing merge error into a clean "never ran".

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/pp_env.sh"

SHARD_DIR="${SHARD_DIR:?set SHARD_DIR to the directory the shard array wrote to}"
OUT_PATH="${OUT_PATH:-$PP_REPO/PlanetProfile/Test/mcmc_results/Enceladus/Cassini_isostasy/enceladus_zbw_production_cache.pkl}"

echo "merging shards from : $SHARD_DIR"
echo "writing cache to    : $OUT_PATH"

# Build into scratch first, then move -- the repo lives in Dropbox on the
# workstation side, and writing a multi-GB pickle directly into a synced
# folder invites a partial-file sync. (Harmless on a pure-cluster path.)
STAGED="$PP_SCRATCH/$(basename "$OUT_PATH")"

srun python plans/scripts/enceladus_zbw_merge_shards.py \
  --shard-dir "$SHARD_DIR" \
  --out "$STAGED"

mkdir -p "$(dirname "$OUT_PATH")"
mv "$STAGED" "$OUT_PATH"
echo "cache in place: $OUT_PATH"
ls -la "$OUT_PATH"

echo
echo "NEXT: report the six acceptance gates from"
echo "  validation_reports/enceladus_isostasy/r7_production_build_report.json"
echo "to Machine A for verification (r7 §0.29 sequence: gates -> A verifies"
echo "-> two-run reference protocol). Gates are interpreted, never tuned."
