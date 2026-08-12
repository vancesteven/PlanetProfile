"""§0.20 R2 — N=100k corrected conditioning at the fiducial (tighter corrected SE).

Reviewer-cleared 2026-08-11 (PASS-WITH-CONCERNS on R1). Re-conditions the
deployed NH3 flow at the fiducial x with N_posterior=100k (deterministic seed
72), computes the SNIS-corrected ocean fraction and a bootstrap SE on it.

Preregistered reading (Machine A §0.20 — NOT self-adjudicated; the
scientific-reviewer interprets):
  * residual PERSISTS at the tighter corrected SE  -> tension is NOT finite-N
    corrected-side bias; proceed to R3 (decisive reference recompute).
  * residual SHRINKS materially                    -> 1/ESS-ceiling argument
    fails, the corrected/SNIS side reopens, R3 pauses.

Diagnostic only. Consumes no gate. Does NOT re-ratify C16 (R3 is decisive).
Reuses the exact R1 SNIS conditioning machinery, scaled 20k -> 100k, plus a
bootstrap SE on the weighted corrected fraction.

Output: /tmp/r2_corrected_100k.json
"""
import json
import time
from pathlib import Path

import numpy as np

REPO = Path('/Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai')
CONFIG = REPO / ('PlanetProfile/Inference/configs/'
                 'test54_titan_nh3_freegrav.json')
N_POST = 100000
N_BOOT = 2000

# committed references (from R1 / C13)
REF_A = 0.9172547835733885            # R1 correct-pooled reference (fiducial cmp)
FIDUCIAL_R1 = 0.9292539223381953      # R1 N=20k corrected fiducial
POOLED_RESIDUAL = 0.0149              # committed pooled 3-seed residual (context)


def main():
    from PlanetProfile.Inference.inference_core import InferenceConfig
    from PlanetProfile.Inference.sbi_runner import SBIRunner
    from PlanetProfile.Inference.is_correction import (
        compute_is_correction, pareto_k_fit,
        PARETO_K_CLEAN, PARETO_K_MAX)

    cfg = InferenceConfig.from_json(str(CONFIG))
    cfg.mode = 'sbi'
    run_seed = int(cfg.random_state)  # seed offset 0 == fiducial (seed 72)

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
    lw = np.asarray(corr.log_w, float)

    # ocean mask over flow draws (same convention as is_correction.py:251 / R1)
    d_ocean = np.nan_to_num(np.asarray(res.D_ocean_results, float), nan=0.0)
    has_ocean = d_ocean > 0.5

    finite = np.isfinite(lw)
    n_rejected = int((~finite).sum())
    lwf = lw[finite] - lw[finite].max()
    w = np.exp(lwf)                       # unnormalised finite weights
    o = has_ocean[finite].astype(float)   # ocean indicator on finite draws

    # point estimate: weighted corrected ocean fraction
    frac = float((w * o).sum() / w.sum())

    # ESS (Kish) on the corrected weights
    ess = float((w.sum() ** 2) / np.sum(w ** 2))
    ess_over_n = ess / w.size

    # Pareto-k tail diagnostic (deployed convention)
    k_hat, _, _ = pareto_k_fit(w)
    psis_would_smooth = bool(PARETO_K_CLEAN < k_hat <= PARETO_K_MAX)

    # bootstrap SE of the weighted fraction: resample finite draws w/ replacement,
    # recompute weighted fraction each rep. Seeded (fixed offset) for determinism.
    rng = np.random.default_rng(run_seed)
    nf = w.size
    boot = np.empty(N_BOOT)
    for b in range(N_BOOT):
        idx = rng.integers(0, nf, nf)
        wb = w[idx]
        boot[b] = (wb * o[idx]).sum() / wb.sum()
    boot_se = float(boot.std(ddof=1))
    boot_mean = float(boot.mean())
    ci_lo, ci_hi = (float(np.percentile(boot, 2.5)),
                    float(np.percentile(boot, 97.5)))

    residual_vs_refA = float(frac - REF_A)
    # preregistered reading helpers (interpreted by reviewer, not here)
    resid_over_se = (residual_vs_refA / boot_se) if boot_se > 0 else float('inf')

    report = {
        'purpose': 'R2 N=100k corrected conditioning at fiducial (tighter SE)',
        'seed': run_seed,
        'n_posterior_requested': N_POST,
        'n_finite': int(nf),
        'n_rejected': n_rejected,
        'wall_condition_s': wall,
        'pareto_k': float(k_hat),
        'psis_gate_would_smooth': psis_would_smooth,
        'kish_ess': ess,
        'ess_over_n': ess_over_n,
        'corrected_ocean_fraction': frac,
        'bootstrap': {
            'n_boot': N_BOOT,
            'se': boot_se,
            'mean': boot_mean,
            'ci95': [ci_lo, ci_hi],
        },
        'reference_A_fiducial': REF_A,
        'residual_vs_referenceA': residual_vs_refA,
        'residual_over_bootSE': resid_over_se,
        'context': {
            'r1_fiducial_N20k': FIDUCIAL_R1,
            'committed_pooled_residual': POOLED_RESIDUAL,
            'note': ('residual here is fiducial (corrected - refA). Pooled '
                     '3-seed residual +0.0149 is R3 remit. R2 tightens the '
                     'CORRECTED-side SE at the fiducial only.'),
        },
        'preregistered_reading': (
            'If residual_vs_referenceA is many bootSE above 0 and comparable to '
            'the R1 N=20k fiducial residual (+0.0120), the corrected side is not '
            'finite-N biased -> proceed to R3. If it shrinks toward 0 relative '
            'to R1, the 1/ESS-ceiling argument fails and the corrected side '
            'reopens. Reviewer adjudicates.'),
    }
    out = Path('/tmp/r2_corrected_100k.json')
    out.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    print(f'\nreport -> {out}')


if __name__ == '__main__':
    main()
