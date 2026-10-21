#!/usr/bin/env bash
# Run the r7 rebuild adjudication conditions C1-C3 (see
# validation_reports/enceladus_isostasy/r7_rebuild_acceptance_2026-10-21.md)
# on the machine that holds the shards and the merged cache:
#   C1  tally the logged exception behind every None node
#   C2  re-merge so r7_production_build_report.json carries true provenance
#       (rewrites the merged cache with identical content)
#   C3  mass-deficit study (~20-30 min)
# Each step runs even if an earlier one fails; exit codes are summarised.
# Outputs land in validation_reports/enceladus_isostasy/ for review.
#
#   plans/cluster/enceladus_zbw_conditions.sh [SHARD_DIR]

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PP_SCRATCH="${PP_SCRATCH:-$HOME/scratch}"
SHARD_DIR="${1:-${SHARD_DIR:-$PP_SCRATCH/enc_zbw_shards}}"
CACHE="${CACHE:-$PP_SCRATCH/enceladus_zbw_production_cache.pkl}"
PPHOME="${PPHOME:-$PP_SCRATCH/pphome}"
PP_CONDA_ENV="${PP_CONDA_ENV:-PPcl}"
PP_CONDA_SH="${PP_CONDA_SH:-$HOME/mamba/etc/profile.d/conda.sh}"
VR="$REPO/validation_reports/enceladus_isostasy"
LOG="$VR/r7_rebuild_conditions_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee "$LOG") 2>&1

if [[ -f "$PP_CONDA_SH" ]]; then
  # shellcheck disable=SC1090
  source "$PP_CONDA_SH"
  conda activate "$PP_CONDA_ENV"
fi
cd "$REPO"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
echo "repo @ $(git rev-parse --short HEAD); shards $SHARD_DIR; cache $CACHE; $(date)"
[[ -f "$CACHE" ]] || { echo "ERROR: no merged cache at $CACHE (run enceladus_zbw_finish.sh first)" >&2; exit 1; }

echo; echo "=== C1: None reasons"
plans/cluster/enceladus_zbw_none_reasons.sh "$SHARD_DIR"; rc1=$?

echo; echo "=== C2: re-merge with provenance"
HOME="$PPHOME" python plans/scripts/enceladus_zbw_merge_shards.py \
  --shard-dir "$SHARD_DIR" --out "$CACHE"; rc2=$?
grep -E '"(build_mode|build_commit|merge_commit|merge_host)"' \
  "$VR/r7_production_build_report.json"

echo; echo "=== C3: mass-deficit study ($(date))"
HOME="$PPHOME" python plans/scripts/enceladus_zbw_mass_deficit_study.py \
  --cache "$CACHE" --out-dir "$VR/c3_mass_deficit"; rc3=$?

echo
echo "=== SUMMARY $(date)"
echo "C1 none reasons : exit $rc1  -> $VR/r7_rebuild_none_reasons.txt"
echo "C2 re-merge     : exit $rc2  -> $VR/r7_production_build_report.json"
echo "C3 mass deficit : exit $rc3 (0 = ran; PASS/FAIL is in the report) -> $VR/c3_mass_deficit/"
echo "log             : $LOG"
echo "Nothing is adjudicated here; report results as-is."
(( rc1 == 0 && rc2 == 0 && rc3 == 0 ))
