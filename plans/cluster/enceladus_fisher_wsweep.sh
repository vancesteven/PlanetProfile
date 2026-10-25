#!/usr/bin/env bash
# Fisher identifiability at four salinities (log10 w = -0.5, 1.0, 1.5, 1.9;
# zb 25, rho_ice 925, rho_s 2400) -- reviewer condition 1, 2026-10-25.
# Also saves the rho_grain_max check to JSON.
#   plans/cluster/enceladus_fisher_wsweep.sh
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
VR="$REPO/validation_reports/enceladus_isostasy"
P="zb=25,rho_ice=925,rho_s=2400,C2=0.75,sys=0"
POINTS="lw-0.5:logw=-0.5,$P;lw1.0:logw=1.0,$P;lw1.5:logw=1.5,$P;lw1.9:logw=1.9,$P" \
EXTRA_ARGS="--out-dir $VR/fisher_identifiability_wsweep" \
  "$REPO/plans/cluster/enceladus_fisher.sh"; rc=$?
PP_CONDA_SH="${PP_CONDA_SH:-$HOME/mamba/etc/profile.d/conda.sh}"
# shellcheck disable=SC1090
[[ -f "$PP_CONDA_SH" ]] && { source "$PP_CONDA_SH"; conda activate "${PP_CONDA_ENV:-PPcl}"; }
cd "$REPO" && HOME="${PPHOME:-$HOME/scratch/pphome}" python plans/scripts/enceladus_rho_grain_max.py \
  --json-out "$VR/rho_grain_max.json" >/dev/null && echo "rho_grain_max -> $VR/rho_grain_max.json"
exit $rc
