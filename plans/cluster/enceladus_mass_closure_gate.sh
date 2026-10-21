#!/usr/bin/env bash
# Run the Enceladus mass-closure jitter gate (reviewer ruling 2026-10-21,
# validation_reports/enceladus_isostasy/c3_root_cause_2026-10-21/README.md)
# on the machine that holds the merged (zb, w) cache: 10 nodes x
# deltaP/{1,2,4} x nPoros {8,12,16}, closure ON and OFF, plus the frozen
# two-layer stack mass check. ~90 builds of ~45 s each; JOBS parallel.
# Results land in validation_reports/enceladus_isostasy/mass_closure_gate/.
#
#   plans/cluster/enceladus_mass_closure_gate.sh [JOBS]

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PP_SCRATCH="${PP_SCRATCH:-$HOME/scratch}"
JOBS="${1:-${JOBS:-6}}"
CACHE="${CACHE:-$PP_SCRATCH/enceladus_zbw_production_cache.pkl}"
PPHOME="${PPHOME:-$PP_SCRATCH/pphome}"
PP_CONDA_ENV="${PP_CONDA_ENV:-PPcl}"
PP_CONDA_SH="${PP_CONDA_SH:-$HOME/mamba/etc/profile.d/conda.sh}"
OUT="$REPO/validation_reports/enceladus_isostasy/mass_closure_gate"
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
echo "repo @ $(git rev-parse --short HEAD); cache $CACHE; jobs $JOBS; $(date)"
[[ -f "$CACHE" ]] || { echo "ERROR: no merged cache at $CACHE" >&2; exit 1; }

HOME="$PPHOME" python plans/scripts/enceladus_mass_closure_jitter_gate.py \
  --cache "$CACHE" --out-dir "$OUT" --jobs "$JOBS" --frozen-check
rc=$?
echo; echo "=== exit $rc $(date); outputs in $OUT; log $LOG"
echo "Nothing is adjudicated here; report results as-is."
exit $rc
