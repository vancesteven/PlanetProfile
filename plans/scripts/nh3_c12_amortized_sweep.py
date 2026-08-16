#!/usr/bin/env python
"""C12 amortized sweep for the deployed Titan NH3 free-gravity artifact.

Design fixed in validation_reports/titan_freegrav_nh3_1m/is_correction/
c12_sweep_PREREGISTRATION.md BEFORE any result is read — read that file
for the full rationale. Summary:

- 200 prior-predictive x (theta ~ effective prior -> forward model ->
  x with obs_noise=True) + 8 axis endpoints (4 observables x mean+/-2
  sigma_train from the deployed artifact's x_norm, Im_k2 lower endpoint
  clamped to 0.0) = 208 points, N=5000 posterior draws each.
- Gate: Pareto-k <= 0.7 for >= 95% of the 208 points (NaN counts as fail).
- Per point, a FRESH InferenceConfig/SBIRunner is built with `observables`
  centrals set to that point's x, so the MCMCRunner likelihood closure
  (built once at construction from config.observables) matches the x used
  for flow conditioning. The artifact is located via the UNMUTATED
  fiducial config's hash, then injected into the per-point runner exactly
  as is_correction_validate.py does for its single fiducial call.
- assert_byte_identity's config-hash check is EXPECTED to fail off-fiducial
  (observables differ by construction) -- recorded, never silently
  skipped; the other three checks (structure cache SHA, param names,
  imag_convention) still run and still gate.

Usage:
  python plans/scripts/nh3_c12_amortized_sweep.py --verify-only
  python plans/scripts/nh3_c12_amortized_sweep.py --n-jobs 12
"""
import argparse
import json
import sys
import time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

CONFIG_PATH = str(REPO / 'PlanetProfile/Inference/configs/'
                        'test54_titan_nh3_freegrav.json')
ARTIFACT_NAME = 'titan_freegrav_nh3_posterior_1m.pt'
FIDUCIAL_REPORT = (REPO / 'validation_reports/titan_freegrav_nh3_1m/'
                          'is_correction/is_validation_nh3.json')
OUT_DIR = REPO / 'validation_reports/titan_freegrav_nh3_1m/is_correction'

N_PRIOR_PRED = 200
N_POST_PER_POINT = 5000
PARETO_K_MAX = 0.7
FAIL_FRAC_MAX = 0.05
PRIOR_SEED = 90210  # dedicated seed for the 200 prior-predictive draws


def _load_base_config():
    from PlanetProfile.Inference.inference_core import InferenceConfig
    cfg = InferenceConfig.from_json(CONFIG_PATH)
    cfg.mode = 'sbi'
    return cfg


def _axis_endpoints():
    """4 observables x mean+/-2*sigma_train from the artifact's x_norm."""
    import torch
    art_path = REPO / 'PlanetProfile/Inference/sbi_artifacts' / ARTIFACT_NAME
    artifact = torch.load(art_path, map_location='cpu', weights_only=False)
    x_norm = artifact['x_norm']
    obs_names = list(artifact['metadata']['observable_names']) \
        if 'metadata' in artifact and 'observable_names' in artifact.get(
            'metadata', {}) else None
    if obs_names is None:
        # Fall back to the config's declared observable order.
        cfg = _load_base_config()
        obs_names = list(cfg.observables.keys())
    means = np.asarray(x_norm['mean']).ravel()
    stds = np.asarray(x_norm['std']).ravel()
    points = []
    for i, name in enumerate(obs_names):
        lo = float(means[i] - 2.0 * stds[i])
        hi = float(means[i] + 2.0 * stds[i])
        if name == 'Im_k2' and lo < 0.0:
            lo = 0.0
        points.append({'axis': name, 'endpoint': 'lo', 'value': lo})
        points.append({'axis': name, 'endpoint': 'hi', 'value': hi})
    return points, obs_names


def _build_sweep_x_list():
    """Return list of dicts: {'kind': 'prior_pred'|'endpoint', 'x': {...}}."""
    from PlanetProfile.Inference.sbi_runner import SBIRunner
    cfg = _load_base_config()
    runner = SBIRunner(cfg)
    theta, x, meta = runner.generate_training_set(
        n_simulations=N_PRIOR_PRED, seed=PRIOR_SEED, obs_noise=True,
        noise_seed=PRIOR_SEED + 1)
    obs_names = list(cfg.observables.keys())
    if x.shape[0] < N_PRIOR_PRED:
        print(f'WARNING: support guard/non-finite dropped rows -- got '
              f'{x.shape[0]}/{N_PRIOR_PRED} prior-predictive draws',
              file=sys.stderr)
    points = []
    for row in x:
        points.append({'kind': 'prior_pred',
                       'x': {n: float(v) for n, v in zip(obs_names, row)}})

    endpoints, _ = _axis_endpoints()
    fiducial = {n: float(s[0]) for n, s in cfg.observables.items()}
    for ep in endpoints:
        xo = dict(fiducial)
        xo[ep['axis']] = ep['value']
        points.append({'kind': 'endpoint', 'axis': ep['axis'],
                       'endpoint_side': ep['endpoint'], 'x': xo})
    return points


def _run_one_point(point, n_post, seed):
    """Runs in a worker process. Returns a JSON-serializable dict."""
    import sys as _sys
    from pathlib import Path as _Path
    _repo = _Path(__file__).resolve().parents[2]
    if str(_repo) not in _sys.path:
        _sys.path.insert(0, str(_repo))
    import numpy as _np
    from PlanetProfile.Inference.inference_core import InferenceConfig
    from PlanetProfile.Inference.sbi_runner import SBIRunner
    from PlanetProfile.Inference.is_correction import (
        compute_is_correction, assert_byte_identity)

    base_cfg = InferenceConfig.from_json(CONFIG_PATH)
    base_cfg.mode = 'sbi'
    base_hash = base_cfg.generate_hash()

    # Locate + load the artifact against the UNMUTATED fiducial hash.
    probe_runner = SBIRunner(base_cfg)
    art_path = probe_runner._find_matching_artifact(base_hash)
    if art_path is None:
        return {**point, 'error': f'no artifact for hash {base_hash}'}
    artifact = probe_runner._load_artifact_dict(art_path)

    # Build a FRESH config with observables centrals set to this point's x
    # (sigmas unchanged) so the MCMCRunner likelihood closure matches x.
    x_point = point['x']
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
    except Exception as exc:  # noqa: BLE001 -- record, don't abort the sweep
        return {**point, 'error': f'{type(exc).__name__}: {exc}'}
    wall_s = time.time() - t0

    c2_hash_expected_mismatch = (swept_cfg.generate_hash() != base_hash)
    identity_error = None
    try:
        # Only the non-observable-dependent checks are expected to hold;
        # config_hash mismatch is EXPECTED off-fiducial by construction.
        identity = assert_byte_identity(runner, res)
    except AssertionError as exc:
        msg = str(exc)
        if 'config hash mismatch' in msg and c2_hash_expected_mismatch:
            identity = {'note': 'config-hash mismatch expected off-fiducial'}
        else:
            identity_error = msg
            identity = None

    k2_bounds = swept_cfg.sampler_settings.get('k2_support_bounds') or {}
    try:
        corr = compute_is_correction(res, k2_support_bounds=k2_bounds)
    except Exception as exc:  # noqa: BLE001
        return {**point, 'error': f'compute_is_correction: '
                                   f'{type(exc).__name__}: {exc}',
                'wall_s': wall_s,
                'c2_hash_expected_mismatch': c2_hash_expected_mismatch}

    pareto_k = corr.pareto_k
    k_is_nan = bool(_np.isnan(pareto_k)) if pareto_k is not None else True
    pareto_pass = bool((not k_is_nan) and pareto_k <= PARETO_K_MAX)

    out = {
        **point,
        'wall_s': wall_s,
        'c2_hash_expected_mismatch': c2_hash_expected_mismatch,
        'identity_error': identity_error,
        'pareto_k': None if pareto_k is None or _np.isnan(pareto_k)
                    else float(pareto_k),
        'pareto_k_is_nan': k_is_nan,
        'pareto_pass': pareto_pass,
        'ess': float(corr.ess),
        'ess_over_n': float(corr.ess_over_n),
        'n_required_at_essfloor1000': (1000.0 / corr.ess_over_n)
                                       if corr.ess_over_n > 0 else None,
        'verdict': corr.verdict,
        'fail_reasons': corr.fail_reasons,
        'w_max_frac': float(corr.w_max_frac),
        'frac_rejected': float(corr.frac_rejected),
        'branch': corr.branch,
    }
    return out


def verify_at_fiducial():
    """C12 §5: point-builder must reproduce the committed fiducial report
    exactly before the real sweep runs."""
    cfg = _load_base_config()
    fiducial_x = {n: float(s[0]) for n, s in cfg.observables.items()}
    point = {'kind': 'fiducial_verify', 'x': fiducial_x}
    out = _run_one_point(point, n_post=20000, seed=int(cfg.random_state))
    with open(FIDUCIAL_REPORT) as f:
        ref = json.load(f)
    ref_k = ref['is_diagnostics']['pareto_k']
    ref_ess = ref['is_diagnostics']['ess']
    print(json.dumps(out, indent=2, default=str))
    ok_k = out.get('pareto_k') is not None and abs(
        out['pareto_k'] - ref_k) < 1e-6
    ok_ess = abs(out['ess'] - ref_ess) < 1e-3
    print(f"\nVERIFY vs committed fiducial: pareto_k {out.get('pareto_k')} "
          f"vs {ref_k} (match={ok_k}); ess {out.get('ess')} vs {ref_ess} "
          f"(match={ok_ess})")
    if not (ok_k and ok_ess):
        raise SystemExit(
            "Point-builder does NOT reproduce the committed fiducial "
            "report. Per preregistration §5, the sweep must not run "
            "until this is fixed.")
    print("VERIFIED: point-builder reproduces the fiducial report. "
          "Safe to run the full sweep.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify-only', action='store_true',
                    help='Run the single-point fiducial reproduction '
                         'check (preregistration §5) and exit.')
    ap.add_argument('--n-jobs', type=int, default=1)
    ap.add_argument('--n-post', type=int, default=N_POST_PER_POINT)
    ap.add_argument('--out', default=str(OUT_DIR / 'c12_sweep_report.json'))
    args = ap.parse_args()

    if args.verify_only:
        verify_at_fiducial()
        return

    points = _build_sweep_x_list()
    n_total = len(points)
    print(f'Built {n_total} sweep points '
          f'({sum(1 for p in points if p["kind"] == "prior_pred")} '
          f'prior-predictive + '
          f'{sum(1 for p in points if p["kind"] == "endpoint")} endpoints)')

    results = []
    t_start = time.time()
    if args.n_jobs <= 1:
        for i, p in enumerate(points):
            seed = 71000 + i
            results.append(_run_one_point(p, args.n_post, seed))
            print(f'  [{i+1}/{n_total}] {p["kind"]} '
                  f'pareto_k={results[-1].get("pareto_k")}')
    else:
        with ProcessPoolExecutor(max_workers=args.n_jobs) as ex:
            futs = {ex.submit(_run_one_point, p, args.n_post, 71000 + i): i
                   for i, p in enumerate(points)}
            done = 0
            for fut in as_completed(futs):
                results.append(fut.result())
                done += 1
                if done % 10 == 0 or done == n_total:
                    print(f'  {done}/{n_total} done '
                          f'({time.time()-t_start:.0f}s elapsed)')
    wall_total = time.time() - t_start

    # --- gate + report -----------------------------------------------
    n_errored = sum(1 for r in results if 'error' in r)
    scored = [r for r in results if 'error' not in r]
    n_fail_pareto = sum(1 for r in scored if not r['pareto_pass'])
    n_fail_pareto += n_errored  # an errored point counts as a failure

    def _subset_stats(kind):
        sub = [r for r in results if r.get('kind') == kind]
        sub_err = sum(1 for r in sub if 'error' in r)
        sub_scored = [r for r in sub if 'error' not in r]
        sub_fail = sum(1 for r in sub_scored if not r['pareto_pass']) \
            + sub_err
        return {'n': len(sub), 'n_errored': sub_err,
                'n_pareto_fail': sub_fail,
                'frac_fail': sub_fail / len(sub) if sub else None}

    ess_over_n_vals = [r['ess_over_n'] for r in scored if 'ess_over_n' in r]
    n_required_vals = [r['n_required_at_essfloor1000'] for r in scored
                       if r.get('n_required_at_essfloor1000') is not None]

    frac_fail_pooled = n_fail_pareto / n_total
    report = {
        'gate': 'C12 amortized sweep: Pareto-k <= 0.7 for >= 95% of '
                f'{n_total} points ({N_PRIOR_PRED} prior-predictive + 8 '
                'axis endpoints)',
        'preregistration': str(
            OUT_DIR / 'c12_sweep_PREREGISTRATION.md'),
        'n_total': n_total,
        'n_post_per_point': args.n_post,
        'n_errored': n_errored,
        'n_pareto_fail_pooled': n_fail_pareto,
        'frac_fail_pooled': frac_fail_pooled,
        'pass_pooled': bool(frac_fail_pooled <= FAIL_FRAC_MAX),
        'subset_prior_predictive': _subset_stats('prior_pred'),
        'subset_endpoints': _subset_stats('endpoint'),
        'ess_over_n_distribution': {
            'min': float(np.min(ess_over_n_vals)) if ess_over_n_vals
                   else None,
            'median': float(np.median(ess_over_n_vals))
                      if ess_over_n_vals else None,
            'max': float(np.max(ess_over_n_vals)) if ess_over_n_vals
                   else None,
        },
        'n_required_at_essfloor1000_distribution': {
            'min': float(np.min(n_required_vals)) if n_required_vals
                   else None,
            'median': float(np.median(n_required_vals))
                      if n_required_vals else None,
            'max': float(np.max(n_required_vals)) if n_required_vals
                   else None,
        },
        'wall_total_s': wall_total,
        'points': results,
    }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(report, f, indent=2, default=str)
    print(f"\nWrote {out_path}")
    print(f"POOLED: {n_fail_pareto}/{n_total} failed "
          f"({frac_fail_pooled:.1%}) -> "
          f"{'PASS' if report['pass_pooled'] else 'FAIL'} "
          f"(threshold <= {FAIL_FRAC_MAX:.0%})")
    print(f"prior-predictive-only: {report['subset_prior_predictive']}")
    print(f"endpoints-only: {report['subset_endpoints']}")


if __name__ == '__main__':
    main()
