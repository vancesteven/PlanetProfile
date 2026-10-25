#!/usr/bin/env bash
# Enceladus ocean6D rho_rock headline posterior, run (a): the committed config
# run TWICE with the same seed (reproducibility check), against the Mac
# production cache. Reviewer-cleared 2026-10-25 (Fisher checkpoint).
#
#   PILOT=1 plans/cluster/enceladus_posterior_a.sh   # small timing pilot (n_eff 500)
#   plans/cluster/enceladus_posterior_a.sh           # full run, both repeats
#
# Output: validation_reports/enceladus_isostasy/posterior_a[_pilot]/
#   rep{1,2}_result.pkl, rep{1,2}_config.json, rep{1,2}_progress.json, run log
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PP_SCRATCH="${PP_SCRATCH:-$HOME/scratch}"
CACHE="${CACHE:-$PP_SCRATCH/enceladus_zbw_production_cache.pkl}"
PPHOME="${PPHOME:-$PP_SCRATCH/pphome}"
PP_CONDA_SH="${PP_CONDA_SH:-$HOME/mamba/etc/profile.d/conda.sh}"
CFG="PlanetProfile/Inference/configs/enceladus_cassini_isostasy_ocean6D_rhorock.json"
SEED="${SEED:-75}"
NAME=posterior_a; SETS=(--set "structure_cache_path=\"$CACHE\"")
if [[ -n "${PILOT:-}" ]]; then
  NAME=posterior_a_pilot
  SETS+=(--set sampler_settings.n_effective=500 --set sampler_settings.n_active=256
         --set sampler_settings.n_posterior_samples=2000)
fi
OUT="$REPO/validation_reports/enceladus_isostasy/$NAME"; mkdir -p "$OUT"
exec > >(tee "$OUT/run_$(date +%Y%m%d_%H%M%S).log") 2>&1
# shellcheck disable=SC1090
[[ -f "$PP_CONDA_SH" ]] && { source "$PP_CONDA_SH"; conda activate "${PP_CONDA_ENV:-PPcl}"; }
cd "$REPO"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
[[ -f "$CACHE" ]] || { echo "ERROR: no cache at $CACHE" >&2; exit 1; }
echo "repo @ $(git rev-parse --short HEAD); cache $CACHE; seed $SEED; $NAME; $(date)"
rc=0
for rep in 1 2; do
  P="$OUT/rep${rep}_config.json"
  python plans/scripts/cluster_seeded_config.py --config "$CFG" --seed "$SEED" \
    --out "$P" "${SETS[@]}" || exit 1
  echo "=== rep $rep start $(date)"; t0=$(date +%s)
  HOME="$PPHOME" python -m PlanetProfile.Inference.run_inference_cli --config "$P" \
    --output "$PP_SCRATCH/${NAME}_rep${rep}_result.pkl" \
    --progress "$OUT/rep${rep}_progress.json" --verbose || rc=1
  [[ -f "$PP_SCRATCH/${NAME}_rep${rep}_result.pkl" ]] && \
    mv "$PP_SCRATCH/${NAME}_rep${rep}_result.pkl" "$OUT/rep${rep}_result.pkl"
  echo "=== rep $rep end $(date)  ($(( $(date +%s) - t0 )) s)"
done
echo "outputs in $OUT"; exit $rc
