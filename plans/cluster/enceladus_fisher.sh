#!/usr/bin/env bash
# Fisher identifiability of the Enceladus ocean-branch (zb, log10 w, rho_ice,
# rho_s, C2, libration_sys_frac) inference at a few reference points, through
# the runner's own forward model + full data covariance, BEFORE any MCMC.
# Minutes; no builds. Output: validation_reports/enceladus_isostasy/
# fisher_identifiability/fisher_identifiability.{md,json}.
#
#   plans/cluster/enceladus_fisher.sh                 # default reference points
#   POINTS="a:zb=25,logw=mid,rho_ice=925,rho_s=2400,C2=0.75,sys=0" plans/cluster/enceladus_fisher.sh
#   EXTRA_ARGS="--steps 1,2,3" plans/cluster/enceladus_fisher.sh
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PP_SCRATCH="${PP_SCRATCH:-$HOME/scratch}"
CACHE="${CACHE:-$PP_SCRATCH/enceladus_zbw_production_cache.pkl}"
PPHOME="${PPHOME:-$PP_SCRATCH/pphome}"
PP_CONDA_ENV="${PP_CONDA_ENV:-PPcl}"
PP_CONDA_SH="${PP_CONDA_SH:-$HOME/mamba/etc/profile.d/conda.sh}"
OUT="$REPO/validation_reports/enceladus_isostasy/fisher_identifiability"
mkdir -p "$OUT"
if [[ -f "$PP_CONDA_SH" ]]; then
  # shellcheck disable=SC1090
  source "$PP_CONDA_SH"; conda activate "$PP_CONDA_ENV"
fi
cd "$REPO"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
[[ -f "$CACHE" ]] || { echo "ERROR: no merged cache at $CACHE" >&2; exit 1; }
echo "repo @ $(git rev-parse --short HEAD); cache $CACHE"
PTS=()
[[ -n "${POINTS:-}" ]] && PTS=(--points "$POINTS")
# shellcheck disable=SC2086
HOME="$PPHOME" python plans/scripts/enceladus_fisher_identifiability.py \
  --cache "$CACHE" --out-dir "$OUT" ${PTS[@]+"${PTS[@]}"} ${EXTRA_ARGS:-}
echo "outputs in $OUT"
