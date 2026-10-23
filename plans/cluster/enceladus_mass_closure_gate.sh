#!/usr/bin/env bash
# Run the Enceladus mass-closure jitter gate (reviewer ruling 2026-10-21,
# validation_reports/enceladus_isostasy/c3_root_cause_2026-10-21/README.md)
# on the machine that holds the merged (zb, w) cache: 10 nodes x
# deltaP/{1,2,4} x nPoros {8,12,16}, closure ON and OFF, plus the frozen
# two-layer stack mass check, then (mass_moi only) the no-rebuild feasibility
# scan of every cached ocean node. ~90 builds of ~75 s each; JOBS parallel.
# Results land in validation_reports/enceladus_isostasy/mass_closure_gate_<MODE>/
# (the 2026-10-22 rock_density run is in mass_closure_gate/).
#
#   plans/cluster/enceladus_mass_closure_gate.sh [JOBS] [MODE]   MODE: mass_moi (default) | rock_density

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PP_SCRATCH="${PP_SCRATCH:-$HOME/scratch}"
JOBS="${1:-${JOBS:-6}}"
MODE="${2:-${MODE:-mass_moi}}"
CACHE="${CACHE:-$PP_SCRATCH/enceladus_zbw_production_cache.pkl}"
PPHOME="${PPHOME:-$PP_SCRATCH/pphome}"
PP_CONDA_ENV="${PP_CONDA_ENV:-PPcl}"
PP_CONDA_SH="${PP_CONDA_SH:-$HOME/mamba/etc/profile.d/conda.sh}"
OUT="$REPO/validation_reports/enceladus_isostasy/mass_closure_gate_$MODE"
mkdir -p "$OUT"
LOG="$OUT/gate_$(date +%Y%m%d_%H%M%S).log"

exec > >(tee "$LOG") 2>&1

if [[ -f "$PP_CONDA_SH" ]]; then
  # shellcheck disable=SC1090
  source "$PP_CONDA_SH"
  conda activate "$PP_CONDA_ENV"
fi
cd "$REPO"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
echo "repo @ $(git rev-parse --short HEAD); cache $CACHE; jobs $JOBS; mode $MODE; $(date)"
[[ -f "$CACHE" ]] || { echo "ERROR: no merged cache at $CACHE" >&2; exit 1; }

HOME="$PPHOME" python plans/scripts/enceladus_mass_closure_jitter_gate.py \
  --cache "$CACHE" --out-dir "$OUT" --jobs "$JOBS" --frozen-check \
  --closure-mode "$MODE"
rc=$?
rc_scan=0
if [[ "$MODE" == mass_moi ]]; then
  echo; echo "=== cache-wide mass_moi feasibility scan $(date)"
  HOME="$PPHOME" python plans/scripts/enceladus_mass_closure_jitter_gate.py \
    --cache "$CACHE" --out-dir "$OUT" --scan-cache
  rc_scan=$?
fi
echo; echo "=== gate exit $rc, scan exit $rc_scan $(date); outputs in $OUT; log $LOG"
echo "Nothing is adjudicated here; report results as-is."
(( rc == 0 && rc_scan == 0 ))
