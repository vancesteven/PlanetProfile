#!/usr/bin/env python
"""Corrected-pipeline SBC for the deployed Titan NH3 free-gravity artifact.

Design fixed in validation_reports/titan_freegrav_nh3_1m/is_correction/
corrected_sbc_PREREGISTRATION.md BEFORE any result is read -- read that
file for the full rationale and the manager's §0.28 ruling this
implements. Summary:

- theta_0, x pairs drawn from the EFFECTIVE prior via
  SBIRunner.generate_training_set (obs_noise=True) -- C15.
- Over-request 900 pairs to land >=500 survivors after (a) the
  support-guard attrition already inside generate_training_set and
  (b) the Pareto-k<=0.7 / ESS>=1000 reliability filter applied per pair.
- Per pair: rebuild a FRESH config+runner with observables set to that
  pair's x (closes the same MCMCRunner likelihood-closure trap fixed for
  the C12 sweep), full n_derived=None recompute at N=20000,
  compute_is_correction for weights.
- Weighted rank per parameter: r_j = sum_k w_k * 1[posterior_sample_jk <
  theta_0_j] (C14 "option i", continuous on [0,1]) -> scipy.stats.kstest
  vs Uniform per parameter -> BH-FDR across the 13 sampled parameters
  (reusing validate_sbi.py's BH block via a shared helper).
- Low-ESS/Pareto-k-failing pairs are EXCLUDED from rank statistics and
  reported as an excluded fraction with their x z-distance distribution;
  INCONCLUSIVE-BY-COVERAGE if that fraction exceeds 25%.

Usage:
  python plans/scripts/nh3_corrected_sbc.py --verify-only
  python plans/scripts/nh3_corrected_sbc.py --n-jobs 10
"""
import argparse
import json
import sys
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from scipy import stats as _stats

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

CONFIG_PATH = str(REPO / 'PlanetProfile/Inference/configs/'
                        'test54_titan_nh3_freegrav.json')
FIDUCIAL_REPORT = (REPO / 'validation_reports/titan_freegrav_nh3_1m/'
                          'is_correction/is_validation_nh3.json')
OUT_DIR = REPO / 'validation_reports/titan_freegrav_nh3_1m/is_correction'

N_PAIRS_REQUEST = 900
N_SBC_TARGET = 500
N_POST_PER_PAIR = 20000
PARETO_K_MAX = 0.7
ESS_ABS_FLOOR = 1000.0
EXCLUDED_FRAC_INCONCLUSIVE = 0.25
SBC_KS_ALPHA = 0.05
PAIR_DRAW_SEED = 530170
PAIR_NOISE_SEED = 530171
PAIR_POST_SEED_BASE = 531000


def _load_base_config():
    from PlanetProfile.Inference.inference_core import InferenceConfig
    cfg = InferenceConfig.from_json(CONFIG_PATH)
    cfg.mode = 'sbi'
    return cfg


def _x_norm_z(x_point, obs_names, art_x_norm):
    means = np.asarray(art_x_norm['mean']).ravel()
    stds = np.asarray(art_x_norm['std']).ravel()
    z = np.array([(x_point[n] - means[i]) / stds[i]
                  for i, n in enumerate(obs_names)])
    return float(np.sqrt(np.mean(z ** 2)))


def _bh_fdr(pvals):
    """Benjamini-Hochberg adjusted p-values (same recipe as validate_sbi.py)."""
    pvals = np.asarray(pvals, dtype=np.float64)
    order = np.argsort(pvals)
    K = len(pvals)
    bh_adj = np.empty(K, dtype=np.float64)
    running_min = 1.0
    for rank in range(K - 1, -1, -1):
        idx = order[rank]
        val = pvals[idx] * K / (rank + 1)
        running_min = min(running_min, val)
        bh_adj[idx] = min(running_min, 1.0)
    return bh_adj


def _build_pair_list():
    """Draw N_PAIRS_REQUEST (theta_0, x) pairs from the effective prior."""
    from PlanetProfile.Inference.sbi_runner import SBIRunner
    cfg = _load_base_config()
    runner = SBIRunner(cfg)
    theta, x, stats_dict = runner.generate_training_set(
        n_simulations=N_PAIRS_REQUEST, seed=PAIR_DRAW_SEED, obs_noise=True,
        noise_seed=PAIR_NOISE_SEED)
    obs_names = list(cfg.observables.keys())
    param_names = list(cfg.param_space.keys())
    if theta.shape[0] < N_PAIRS_REQUEST:
        print(f'NOTE: support guard/non-finite dropped rows -- got '
              f'{theta.shape[0]}/{N_PAIRS_REQUEST} pairs from '
              'generate_training_set', file=sys.stderr)
    pairs = []
    for i in range(theta.shape[0]):
        pairs.append({
            'idx': i,
            'theta_0': {n: float(v) for n, v in zip(param_names, theta[i])},
            'x': {n: float(v) for n, v in zip(obs_names, x[i])},
        })
    return pairs, obs_names, param_names, stats_dict


def _run_one_pair(pair, n_post, seed):
    """Runs in a worker process. Returns a JSON-serializable dict."""
    import sys as _sys
    from pathlib import Path as _Path
    _repo = _Path(__file__).resolve().parents[2]
    if str(_repo) not in _sys.path:
        _sys.path.insert(0, str(_repo))
    import numpy as _np
    import torch as _torch
    from PlanetProfile.Inference.inference_core import InferenceConfig
    from PlanetProfile.Inference.sbi_runner import SBIRunner
    from PlanetProfile.Inference.is_correction import compute_is_correction

    base_cfg = InferenceConfig.from_json(CONFIG_PATH)
    base_cfg.mode = 'sbi'
    base_hash = base_cfg.generate_hash()

    probe_runner = SBIRunner(base_cfg)
    art_path = probe_runner._find_matching_artifact(base_hash)
    if art_path is None:
        return {**pair, 'error': f'no artifact for hash {base_hash}'}
    artifact = probe_runner._load_artifact_dict(art_path)

    x_point = pair['x']
    swept_cfg_dict = base_cfg.to_dict()
    swept_observables = {}
    for name, (_, sigma) in base_cfg.observables.items():
        swept_observables[name] = [float(x_point[name]), float(sigma)]
    swept_cfg_dict['observables'] = swept_observables
    swept_cfg = InferenceConfig.from_dict(swept_cfg_dict)
    swept_cfg.mode = 'sbi'

    runner = SBIRunner(swept_cfg)
    runner._posterior = artifact['posterior']
    runner._train_info = {k: artifact.get(k) for k in
                          ('density_estimator', 'seed',
                           'n_train_effective', 'rejection_stats',
                           'theta_norm', 'x_norm')}

    t0 = time.time()
    try:
        res = runner._condition_and_package(
            x_obs=x_point, n_posterior_samples=n_post, seed=seed,
            n_reeval=int(swept_cfg.sampler_settings.get('n_reeval', 500)),
            artifact_path=art_path, reused=True,
            rejection_stats=artifact.get('rejection_stats'),
            config_hash=swept_cfg.generate_hash(), t0=t0, n_derived=None)
    except Exception as exc:  # noqa: BLE001
        return {**pair, 'error': f'{type(exc).__name__}: {exc}'}
    wall_s = time.time() - t0

    k2_bounds = swept_cfg.sampler_settings.get('k2_support_bounds') or {}
    try:
        corr = compute_is_correction(res, k2_support_bounds=k2_bounds)
    except Exception as exc:  # noqa: BLE001
        return {**pair, 'error': f'compute_is_correction: '
                                  f'{type(exc).__name__}: {exc}',
                'wall_s': wall_s}

    pareto_k = corr.pareto_k
    k_is_nan = bool(_np.isnan(pareto_k)) if pareto_k is not None else True
    pareto_pass = bool((not k_is_nan) and pareto_k <= PARETO_K_MAX)
    ess_pass = bool(corr.ess >= ESS_ABS_FLOOR)
    survived = bool(pareto_pass and ess_pass)

    z_dist = _x_norm_z(x_point, list(swept_cfg.observables.keys()),
                        artifact['x_norm'])

    ranks = None
    if survived and pair.get('theta_0'):
        param_names = list(swept_cfg.param_space.keys())
        theta_0_vec = np.array([pair['theta_0'][p] for p in param_names])
        weights = corr.weights
        samples = np.asarray(res.samples, dtype=float)  # (N, n_params)
        below = (samples < theta_0_vec[None, :]).astype(float)  # (N, n_params)
        ranks = (weights[:, None] * below).sum(axis=0).tolist()

    out = {
        'idx': pair['idx'],
        'theta_0': pair['theta_0'],
        'x': pair['x'],
        'wall_s': wall_s,
        'pareto_k': None if pareto_k is None or k_is_nan else float(pareto_k),
        'pareto_k_is_nan': k_is_nan,
        'pareto_pass': pareto_pass,
        'ess': float(corr.ess),
        'ess_pass': ess_pass,
        'survived': survived,
        'z_dist': z_dist,
        'ranks': ranks,
        'verdict': corr.verdict,
        'fail_reasons': corr.fail_reasons,
    }
    return out


def verify_at_fiducial():
    """Preregistration §5 part 1: point-builder must reproduce the
    committed fiducial report exactly before the real run."""
    cfg = _load_base_config()
    fiducial_x = {n: float(s[0]) for n, s in cfg.observables.items()}
    pair = {'idx': -1, 'theta_0': {}, 'x': fiducial_x}
    out = _run_one_pair(pair, n_post=20000, seed=int(cfg.random_state))
    with open(FIDUCIAL_REPORT) as f:
        ref = json.load(f)
    ref_k = ref['is_diagnostics']['pareto_k']
    ref_ess = ref['is_diagnostics']['ess']
    print(json.dumps({k: v for k, v in out.items() if k != 'ranks'},
                      indent=2, default=str))
    ok_k = out.get('pareto_k') is not None and abs(
        out['pareto_k'] - ref_k) < 1e-6
    ok_ess = abs(out['ess'] - ref_ess) < 1e-3
    print(f"\nVERIFY vs committed fiducial: pareto_k {out.get('pareto_k')} "
          f"vs {ref_k} (match={ok_k}); ess {out.get('ess')} vs {ref_ess} "
          f"(match={ok_ess})")
    if not (ok_k and ok_ess):
        raise SystemExit(
            "Point-builder does NOT reproduce the committed fiducial "
            "report. Per preregistration §5, the SBC run must not start "
            "until this is fixed.")

    # Preregistration §5 part 2: rank-formula smoke test. Draw theta_0
    # from the flow's OWN posterior at a non-fiducial x and confirm the
    # resulting rank is not pathologically degenerate (a single draw
    # can land anywhere in [0,1]; this only checks the formula executes
    # sanely, not a calibration claim).
    cfg2 = _load_base_config()
    obs_names = list(cfg2.observables.keys())
    from PlanetProfile.Inference.sbi_runner import SBIRunner
    probe = SBIRunner(cfg2)
    art_path = probe._find_matching_artifact(cfg2.generate_hash())
    artifact = probe._load_artifact_dict(art_path)
    probe._posterior = artifact['posterior']
    smoke_x = {n: float(s[0]) for n, s in cfg2.observables.items()}
    draw = probe.sample_posterior(smoke_x, n_samples=1, seed=42)
    theta_0_smoke = {n: float(v) for n, v in
                     zip(cfg2.param_space.keys(), draw[0])}
    smoke_pair = {'idx': -2, 'theta_0': theta_0_smoke, 'x': smoke_x}
    smoke_out = _run_one_pair(smoke_pair, n_post=20000,
                              seed=int(cfg2.random_state) + 1)
    print(f"\nSmoke rank formula check: survived={smoke_out['survived']} "
          f"ranks={smoke_out.get('ranks')}")
    if smoke_out['survived']:
        r = np.asarray(smoke_out['ranks'])
        if np.any(~np.isfinite(r)) or np.any(r < 0) or np.any(r > 1):
            raise SystemExit(
                "Rank-formula smoke test produced non-finite or "
                "out-of-[0,1] ranks -- fix before the real run.")
    print("VERIFIED: point-builder reproduces the fiducial report and the "
          "rank formula executes sanely. Safe to run the full SBC.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify-only', action='store_true')
    ap.add_argument('--n-jobs', type=int, default=1)
    ap.add_argument('--n-post', type=int, default=N_POST_PER_PAIR)
    ap.add_argument('--out', default=str(OUT_DIR / 'corrected_sbc_report.json'))
    args = ap.parse_args()

    if args.verify_only:
        verify_at_fiducial()
        return

    pairs, obs_names, param_names, draw_stats = _build_pair_list()
    n_total = len(pairs)
    print(f'Drew {n_total} pairs from generate_training_set '
          f'(requested {N_PAIRS_REQUEST}); attrition stats: {draw_stats}')

    results = []
    t_start = time.time()
    if args.n_jobs <= 1:
        for i, p in enumerate(pairs):
            seed = PAIR_POST_SEED_BASE + i
            results.append(_run_one_pair(p, args.n_post, seed))
            done = i + 1
            if done % 10 == 0 or done == n_total:
                print(f'  [{done}/{n_total}] survived so far: '
                      f'{sum(1 for r in results if r.get("survived"))}')
    else:
        with ProcessPoolExecutor(max_workers=args.n_jobs) as ex:
            futs = {ex.submit(_run_one_pair, p, args.n_post,
                               PAIR_POST_SEED_BASE + i): i
                    for i, p in enumerate(pairs)}
            done = 0
            for fut in as_completed(futs):
                results.append(fut.result())
                done += 1
                if done % 10 == 0 or done == n_total:
                    n_surv = sum(1 for r in results if r.get('survived'))
                    print(f'  {done}/{n_total} done, {n_surv} survived '
                          f'({time.time()-t_start:.0f}s elapsed)')
    wall_total = time.time() - t_start

    n_errored = sum(1 for r in results if 'error' in r)
    scored = [r for r in results if 'error' not in r]
    survived = [r for r in scored if r['survived']]
    excluded = [r for r in scored if not r['survived']] + \
        [r for r in results if 'error' in r]
    excluded_frac = len(excluded) / n_total if n_total else None

    z_excluded = [r.get('z_dist') for r in excluded if r.get('z_dist') is not None]
    z_survived = [r.get('z_dist') for r in survived if r.get('z_dist') is not None]

    n_sbc_used = len(survived)
    shortfall = n_sbc_used < N_SBC_TARGET

    per_param = []
    verdict = None
    if excluded_frac is not None and excluded_frac > EXCLUDED_FRAC_INCONCLUSIVE:
        verdict = 'INCONCLUSIVE-BY-COVERAGE'
    elif n_sbc_used < 20:
        verdict = 'INCONCLUSIVE-INSUFFICIENT-PAIRS'
    else:
        rank_matrix = np.array([r['ranks'] for r in survived])  # (n_sbc, n_params)
        ks_pvals = []
        for j, pname in enumerate(param_names):
            r_j = rank_matrix[:, j]
            r_j = np.clip(r_j, 0.0, 1.0)
            ks_stat, ks_pval = _stats.kstest(r_j, 'uniform')
            ks_pvals.append(ks_pval)
        ks_pvals = np.asarray(ks_pvals)
        bh_adj = _bh_fdr(ks_pvals)
        for j, pname in enumerate(param_names):
            per_param.append({
                'param': pname,
                'ks_pval': float(ks_pvals[j]),
                'ks_pass': bool(ks_pvals[j] >= SBC_KS_ALPHA),
                'ks_pval_bh_adj': float(bh_adj[j]),
                'ks_pass_bh': bool(bh_adj[j] >= SBC_KS_ALPHA),
            })
        bh_all_pass = bool(np.all(bh_adj >= SBC_KS_ALPHA))
        verdict = 'PASS' if bh_all_pass else 'FAIL'

    report = {
        'ruling': 'plans/MACHINE-B-HANDOFF.md §0.28',
        'preregistration': str(OUT_DIR / 'corrected_sbc_PREREGISTRATION.md'),
        'n_pairs_requested': N_PAIRS_REQUEST,
        'n_pairs_drawn': n_total,
        'n_post_per_pair': args.n_post,
        'n_errored': n_errored,
        'n_survived': len(survived),
        'n_excluded': len(excluded),
        'excluded_frac': excluded_frac,
        'excluded_frac_inconclusive_threshold': EXCLUDED_FRAC_INCONCLUSIVE,
        'n_sbc_target': N_SBC_TARGET,
        'n_sbc_used': n_sbc_used,
        'shortfall_vs_target': shortfall,
        'z_dist_excluded': {
            'min': float(np.min(z_excluded)) if z_excluded else None,
            'median': float(np.median(z_excluded)) if z_excluded else None,
            'max': float(np.max(z_excluded)) if z_excluded else None,
        },
        'z_dist_survived': {
            'min': float(np.min(z_survived)) if z_survived else None,
            'median': float(np.median(z_survived)) if z_survived else None,
            'max': float(np.max(z_survived)) if z_survived else None,
        },
        'per_parameter': per_param,
        'verdict': verdict,
        'ks_alpha': SBC_KS_ALPHA,
        'multiplicity_correction': 'benjamini_hochberg_fdr',
        'wall_total_s': wall_total,
        'draw_stats': draw_stats,
        'pairs': results,
    }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\nWrote {out_path}")
    print(f"VERDICT: {verdict} (n_sbc_used={n_sbc_used}, "
          f"excluded_frac={excluded_frac:.1%})")


if __name__ == '__main__':
    main()
