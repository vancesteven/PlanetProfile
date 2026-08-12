"""§0.20 R1 — SNIS-side ocean-fraction weight-treatment cross-check.

Re-conditions the NH3 flow at the fiducial (seed 72, deterministic), asserts
the recomputed log-weights reproduce the committed log_weights.npy byte-for-
byte (determinism check), then reports the corrected ocean fraction under:
  raw       : normalised raw weights (what the committed report used; k<0.5)
  psis      : PSIS-smoothed tail weights (forced, to see tail sensitivity)
  ablate01  : top-0.1% weights zeroed then renormalised

Preregistered reading (manager): if the +0.0149 residual moves materially
under any treatment, the corrected (SNIS) side reopens. Diagnostic only —
does NOT consume the corrected ocean fraction as a gate.

Output: /tmp/r1_snis_audit.json
"""
import json
import time
from pathlib import Path

import numpy as np

REPO = Path('/Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai')
SAVED_LW = REPO / ('validation_reports/titan_freegrav_nh3_1m/is_correction/'
                   'log_weights.npy')
CONFIG = REPO / ('PlanetProfile/Inference/configs/'
                 'test54_titan_nh3_freegrav.json')
N_POST = 20000


def main():
    from PlanetProfile.Inference.inference_core import InferenceConfig
    from PlanetProfile.Inference.sbi_runner import SBIRunner
    from PlanetProfile.Inference.is_correction import (
        compute_is_correction, pareto_k_fit, _psis_smooth,
        PARETO_K_CLEAN, PARETO_K_MAX)

    cfg = InferenceConfig.from_json(str(CONFIG))
    cfg.mode = 'sbi'
    run_seed = int(cfg.random_state)  # seed offset 0 == fiducial

    t0 = time.time()
    runner = SBIRunner(cfg)
    base_hash = cfg.generate_hash()
    art = runner._find_matching_artifact(base_hash)
    if art is None:
        raise SystemExit(f'no artifact matches base hash {base_hash}')
    artifact = runner._load_artifact_dict(art)
    runner._posterior = artifact['posterior']
    runner._train_info = {k: artifact.get(k) for k in
                          ('density_estimator', 'seed', 'n_train_effective',
                           'rejection_stats', 'theta_norm', 'x_norm')}
    x_obs = {n: s[0] for n, s in cfg.observables.items()}
    res = runner._condition_and_package(
        x_obs=x_obs, n_posterior_samples=N_POST, seed=run_seed,
        n_reeval=int(cfg.sampler_settings.get('n_reeval', 500)),
        artifact_path=art, reused=True,
        rejection_stats=artifact.get('rejection_stats'),
        config_hash=base_hash, t0=t0, n_derived=None)
    wall = time.time() - t0

    k2_bounds = cfg.sampler_settings.get('k2_support_bounds') or {}
    corr = compute_is_correction(res, k2_support_bounds=k2_bounds)

    # determinism check vs committed log_weights.npy
    saved_lw = np.load(SAVED_LW)
    now_lw = corr.log_w
    both_finite = np.isfinite(saved_lw) & np.isfinite(now_lw)
    both_inf = np.isinf(saved_lw) & np.isinf(now_lw)
    max_abs = float(np.max(np.abs(saved_lw[both_finite] - now_lw[both_finite])))
    reproduces = bool(both_finite.sum() + both_inf.sum() == len(saved_lw)
                      and max_abs < 1e-9)

    # ocean mask over flow draws
    d_ocean = np.nan_to_num(np.asarray(res.D_ocean_results, float), nan=0.0)
    has_ocean = d_ocean > 0.5

    finite = np.isfinite(now_lw)
    lw = now_lw[finite] - now_lw[finite].max()
    w_raw_f = np.exp(lw)
    N = len(now_lw)

    def frac(wfull):
        wn = wfull / wfull.sum()
        return float(wn[has_ocean].sum())

    # raw
    w_raw = np.zeros(N)
    w_raw[finite] = w_raw_f
    f_raw = frac(w_raw)

    # PSIS-smoothed (force, regardless of k gate, to probe tail sensitivity)
    k_hat, tail_idx, gpd = pareto_k_fit(w_raw_f)
    w_psis_f = _psis_smooth(w_raw_f.copy(), tail_idx, gpd, k_hat)
    w_psis = np.zeros(N)
    w_psis[finite] = w_psis_f
    f_psis = frac(w_psis)

    # top-0.1% ablation (zero the largest 0.1% of finite weights)
    w_ab = w_raw.copy()
    n_ab = max(1, int(round(0.001 * finite.sum())))
    fin_idx = np.where(finite)[0]
    top = fin_idx[np.argsort(w_raw[fin_idx])[-n_ab:]]
    w_ab[top] = 0.0
    f_ab = frac(w_ab)

    ref_A = 0.9172547835733885  # R1 correct-pooled reference
    treatments = {
        'raw': f_raw,
        'psis_forced': f_psis,
        'ablate_top0.1pct': f_ab,
    }
    residuals = {k: float(v - ref_A) for k, v in treatments.items()}
    swing = float(max(treatments.values()) - min(treatments.values()))

    report = {
        'purpose': 'R1 SNIS-side ocean-fraction weight-treatment cross-check',
        'seed': run_seed,
        'wall_condition_s': wall,
        'reproduces_saved_log_weights': reproduces,
        'max_abs_logw_diff': max_abs,
        'n_rejected': int((~finite).sum()),
        'pareto_k': float(k_hat),
        'psis_gate_would_smooth': bool(PARETO_K_CLEAN < k_hat <= PARETO_K_MAX),
        'n_ablated': int(n_ab),
        'ocean_fraction_by_treatment': treatments,
        'residual_vs_referenceA': residuals,
        'treatment_swing': swing,
        'committed_residual': 0.0149,
        'reading': (
            'psis_forced and ablate treatments probe whether a few heavy tail '
            'weights drive the ocean mass. If all three fractions stay within '
            'a small fraction of the +0.0149 residual, the corrected side is '
            'weight-treatment robust and the residual is NOT a heavy-tail / '
            'SNIS artifact -> reference side (R3) is the remaining suspect.')
    }
    out = Path('/tmp/r1_snis_audit.json')
    out.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    print(f'\nreport -> {out}')


if __name__ == '__main__':
    main()
