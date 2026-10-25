#!/usr/bin/env bash
# No-rebuild per-node mass_rho_rock closure scan of the merged Enceladus
# (zb, w) cache at rho_ice = none, 915, 935 kg/m3 and rho_s = 2200 and 2600
# (the rho_rock_kgm3 prior-box edges): feasibility, remaining ocean thickness
# R_o - R_s, derived C/MR^2, dR_s vs PP, binding nodes, and the residual
# ocean-density double-count check (finite difference in rho_o). Minutes; no
# builds. Output: validation_reports/enceladus_isostasy/mass_rho_rock_scan/
# mass_rho_rock_cache_scan.{json,md}.
#
#   plans/cluster/enceladus_mass_rho_rock_scan.sh [RHO_ICE_LIST]   # default none,915,935
#   EXTRA_ARGS="--no-dc-check" plans/cluster/enceladus_mass_rho_rock_scan.sh   # skip FD check
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PP_SCRATCH="${PP_SCRATCH:-$HOME/scratch}"
RHO_ICE="${1:-${RHO_ICE:-none,915,935}}"
CACHE="${CACHE:-$PP_SCRATCH/enceladus_zbw_production_cache.pkl}"
PPHOME="${PPHOME:-$PP_SCRATCH/pphome}"
PP_CONDA_ENV="${PP_CONDA_ENV:-PPcl}"
PP_CONDA_SH="${PP_CONDA_SH:-$HOME/mamba/etc/profile.d/conda.sh}"
OUT="$REPO/validation_reports/enceladus_isostasy/mass_rho_rock_scan"
mkdir -p "$OUT"
if [[ -f "$PP_CONDA_SH" ]]; then
  # shellcheck disable=SC1090
  source "$PP_CONDA_SH"; conda activate "$PP_CONDA_ENV"
fi
cd "$REPO"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
[[ -f "$CACHE" ]] || { echo "ERROR: no merged cache at $CACHE" >&2; exit 1; }
echo "repo @ $(git rev-parse --short HEAD); cache $CACHE; rho_ice $RHO_ICE"
# shellcheck disable=SC2086
HOME="$PPHOME" python plans/scripts/enceladus_mass_closure_jitter_gate.py \
  --cache "$CACHE" --out-dir "$OUT" --scan-cache --closure-mode mass_rho_rock \
  --rho-ice "$RHO_ICE" ${EXTRA_ARGS:-}
echo "outputs in $OUT"
