#!/usr/bin/env bash
# Merge a finished Enceladus (zb, w) shard build and run the acceptance
# gates (§0.31 A1-A4, r7 ACCEPTANCE, §0.30 axes). Run from anywhere.
#
#   plans/cluster/enceladus_zbw_finish.sh [SHARD_DIR]
#
# The merged cache is staged in PP_SCRATCH (default ~/scratch) and is NOT
# copied into the Dropbox-synced repo: that waits until the gates are read
# and adjudicated. An existing merged cache is reused, never overwritten
# (delete it yourself to re-merge). Output is also saved to
# $PP_SCRATCH/enceladus_zbw_finish.<timestamp>.log.

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PP_SCRATCH="${PP_SCRATCH:-$HOME/scratch}"
SHARD_DIR="${1:-${SHARD_DIR:-$PP_SCRATCH/enc_zbw_shards}}"
N_SHARDS="${N_SHARDS:-11}"
CACHE="${CACHE:-$PP_SCRATCH/enceladus_zbw_production_cache.pkl}"
PPHOME="${PPHOME:-$PP_SCRATCH/pphome}"
PP_CONDA_ENV="${PP_CONDA_ENV:-PPcl}"
PP_CONDA_SH="${PP_CONDA_SH:-$HOME/mamba/etc/profile.d/conda.sh}"
LOG="$PP_SCRATCH/enceladus_zbw_finish.$(date +%Y%m%d_%H%M%S).log"

exec > >(tee "$LOG") 2>&1

if [[ -f "$PP_CONDA_SH" ]]; then
  # shellcheck disable=SC1090
  source "$PP_CONDA_SH"
  conda activate "$PP_CONDA_ENV"
fi
cd "$REPO"
echo "repo $REPO @ $(git rev-parse --short HEAD); shards $SHARD_DIR; $(date)"

# Preconditions: refuse to merge an unfinished or C2-failed build.
n_meta=$(ls "$SHARD_DIR"/shard_*_of_*.meta.json 2>/dev/null | wc -l | tr -d ' ')
[[ -f "$SHARD_DIR/DONE" ]] || { echo "ERROR: $SHARD_DIR/DONE missing" >&2; exit 1; }
(( n_meta == N_SHARDS )) || { echo "ERROR: $n_meta/$N_SHARDS .meta.json" >&2; exit 1; }
if ls "$SHARD_DIR"/*QUARANTINED* >/dev/null 2>&1; then
  echo "ERROR: C2 quarantine dump present; do not merge:" >&2
  ls "$SHARD_DIR"/*QUARANTINED* >&2
  exit 1
fi
if grep -l Traceback "$SHARD_DIR"/shard_*.log >/dev/null 2>&1; then
  echo "ERROR: Traceback in:" >&2
  grep -l Traceback "$SHARD_DIR"/shard_*.log >&2
  exit 1
fi
echo "build commit: $(cat "$SHARD_DIR/BUILD_COMMIT" 2>/dev/null || echo missing)"

if [[ -f "$CACHE" ]]; then
  echo "=== merged cache exists, reusing: $CACHE"
else
  echo "=== merge -> $CACHE"
  HOME="$PPHOME" python plans/scripts/enceladus_zbw_merge_shards.py \
    --shard-dir "$SHARD_DIR" --out "$CACHE" \
    || { echo "ERROR: merge failed; no cache written" >&2; exit 1; }
fi

echo "=== acceptance gates"
HOME="$PPHOME" python plans/scripts/enceladus_zbw_acceptance.py \
  --cache "$CACHE" --shard-dir "$SHARD_DIR"
rc=$?

echo "=== acceptance exit $rc  (0 = no FAIL; UNCHECKABLE gates still need review)"
echo "json: $CACHE.acceptance.json"
echo "log : $LOG"
echo "Do NOT move the cache into PlanetProfile/Test/mcmc_results/ until the gates are adjudicated."
exit $rc
