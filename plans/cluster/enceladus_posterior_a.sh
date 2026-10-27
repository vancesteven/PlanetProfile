#!/usr/bin/env bash
# Enceladus ocean6D rho_rock headline posterior, run (a): the committed config
# run TWICE with the same seed (reproducibility check), against the Mac
# production cache. Reviewer-cleared 2026-10-25 (Fisher checkpoint).
#
#   PILOT=1 plans/cluster/enceladus_posterior_a.sh   # small timing pilot (n_eff 500)
#   PREFLIGHT_ONLY=1 plans/cluster/enceladus_posterior_a.sh   # loglike preflight only
#   plans/cluster/enceladus_posterior_a.sh           # full run, both repeats
#   SEED=2027 REPS=1 plans/cluster/enceladus_posterior_a.sh  # different-seed replicate
#                                  # -> posterior_a_seed2027/ (does not overwrite seed 75)
#   VARIANT=g2 REPS=1 plans/cluster/enceladus_posterior_a.sh  # run (g2): log10 w in [-1,1.3]
#   VARIANT=g1 REPS=1 plans/cluster/enceladus_posterior_a.sh  # run (g1): w uniform in ppt
#   THIN_CHECK=1 plans/cluster/enceladus_posterior_a.sh  # max loglike over valid cells at
#                                  # zb 13-15 km (reviewer 2026-10-27) -> posterior_a/thin_zb13_15.json
#
# Output: validation_reports/enceladus_isostasy/posterior_a[_pilot]/
#   preflight.json, rep{1,2}_result.pkl, rep{1,2}_config.json, rep{1,2}_progress.json, run log
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PP_SCRATCH="${PP_SCRATCH:-$HOME/scratch}"
CACHE="${CACHE:-$PP_SCRATCH/enceladus_zbw_production_cache.pkl}"
PPHOME="${PPHOME:-$PP_SCRATCH/pphome}"
PP_CONDA_SH="${PP_CONDA_SH:-$HOME/mamba/etc/profile.d/conda.sh}"
CFG="PlanetProfile/Inference/configs/enceladus_cassini_isostasy_ocean6D_rhorock.json"
SEED="${SEED:-75}"
NAME=posterior_a; [[ "$SEED" != 75 ]] && NAME="posterior_a_seed$SEED"
SETS=(--set "structure_cache_path=\"$CACHE\"")
# Sensitivity variants (reviewer Fisher-checkpoint condition 2, run (g)).
#   VARIANT=g2 : log10 w ~ U[-1, 1.3]  (bounds override only)
#   VARIANT=g1 : w ~ U[0.1, 100] ppt, i.e. log10 w with density ∝ 10^x
#                (prior_type 'uniform-in-linear'; mcmc_runner)
case "${VARIANT:-}" in
  "") ;;
  g2) NAME="posterior_g2"; SETS+=(--set 'param_space.log10_wOcean_ppt.bounds=[-1.0,1.3]') ;;
  g1) NAME="posterior_g1"; SETS+=(--set 'param_space.log10_wOcean_ppt.prior_type="uniform-in-linear"') ;;
  *) echo "unknown VARIANT=$VARIANT" >&2; exit 1 ;;
esac
[[ -n "${VARIANT:-}" && "$SEED" != 75 ]] && NAME="${NAME}_seed$SEED"
if [[ -n "${PILOT:-}" ]]; then
  NAME=posterior_a_pilot
  SETS+=(--set sampler_settings.n_effective=500 --set sampler_settings.n_active=256
         --set sampler_settings.n_posterior_samples=2000)
fi
OUT="${OUT_ROOT:-$REPO/validation_reports/enceladus_isostasy}/$NAME"; mkdir -p "$OUT"
exec > >(tee "$OUT/run_$(date +%Y%m%d_%H%M%S).log") 2>&1
# shellcheck disable=SC1090
[[ -f "$PP_CONDA_SH" ]] && { source "$PP_CONDA_SH"; conda activate "${PP_CONDA_ENV:-PPcl}"; }
cd "$REPO"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
[[ -f "$CACHE" ]] || { echo "ERROR: no cache at $CACHE" >&2; exit 1; }
echo "repo @ $(git rev-parse --short HEAD); cache $CACHE; seed $SEED; $NAME; $(date)"
if [[ -n "${THIN_CHECK:-}" ]]; then
  # Prior draws restricted to zb 13-15 km (just above the None-corner hole at
  # zb <= 13): the best loglike found there should be far below the posterior's.
  T="$OUT/thin_config.json"
  python plans/scripts/cluster_seeded_config.py --config "$CFG" --seed "$SEED" \
    --out "$T" "${SETS[@]}" --set 'param_space.zb_km.bounds=[13.0,15.0]' >/dev/null || exit 1
  HOME="$PPHOME" python plans/scripts/enceladus_loglike_preflight.py --config "$T" \
    --n "${PREFLIGHT_N:-5000}" --json-out "$OUT/thin_zb13_15.json"
  echo "THIN_CHECK written to $OUT/thin_zb13_15.json"; exit 0
fi
# Preflight on the real cache (reviewer 2026-10-25): one central loglike call
# + 1000 prior draws (accept/reject reasons, ms/call). Abort if it fails.
PRE="$OUT/preflight_config.json"
python plans/scripts/cluster_seeded_config.py --config "$CFG" --seed "$SEED" \
  --out "$PRE" "${SETS[@]}" >/dev/null || exit 1
echo "=== preflight $(date)"
HOME="$PPHOME" python plans/scripts/enceladus_loglike_preflight.py --config "$PRE" \
  --n "${PREFLIGHT_N:-1000}" --json-out "$OUT/preflight.json" || {
  echo "PREFLIGHT FAILED -- not starting the sampler" >&2; exit 1; }
[[ -n "${PREFLIGHT_ONLY:-}" ]] && { echo "PREFLIGHT_ONLY set; stopping"; exit 0; }
rc=0
for rep in $(seq 1 "${REPS:-2}"); do
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
