#!/usr/bin/env bash
# No-rebuild per-node mass_moi closure scan of the merged Enceladus (zb, w)
# cache at rho_ice = none, 915, 935 kg/m3 (reporting + the runner's loose
# rejects rho_s_above_grain_max / dR_s_gt_10km). Seconds-to-minutes; no
# builds. Output: validation_reports/enceladus_isostasy/mass_moi_scan/
# mass_moi_cache_scan.{json,md} (per-node rows, 20 most extreme nodes by
# |dR_s| and rho_s, rejects with reasons, per rho_ice).
#
#   plans/cluster/enceladus_mass_moi_scan.sh [RHO_ICE_LIST]   # default none,915,935
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PP_SCRATCH="${PP_SCRATCH:-$HOME/scratch}"
RHO_ICE="${1:-${RHO_ICE:-none,915,935}}"
CACHE="${CACHE:-$PP_SCRATCH/enceladus_zbw_production_cache.pkl}"
PPHOME="${PPHOME:-$PP_SCRATCH/pphome}"
PP_CONDA_ENV="${PP_CONDA_ENV:-PPcl}"
PP_CONDA_SH="${PP_CONDA_SH:-$HOME/mamba/etc/profile.d/conda.sh}"
OUT="$REPO/validation_reports/enceladus_isostasy/mass_moi_scan"
mkdir -p "$OUT"
if [[ -f "$PP_CONDA_SH" ]]; then
  # shellcheck disable=SC1090
  source "$PP_CONDA_SH"; conda activate "$PP_CONDA_ENV"
fi
cd "$REPO"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
[[ -f "$CACHE" ]] || { echo "ERROR: no merged cache at $CACHE" >&2; exit 1; }
echo "repo @ $(git rev-parse --short HEAD); cache $CACHE; rho_ice $RHO_ICE"
HOME="$PPHOME" python plans/scripts/enceladus_mass_closure_jitter_gate.py \
  --cache "$CACHE" --out-dir "$OUT" --scan-cache --rho-ice "$RHO_ICE"
echo "outputs in $OUT"
