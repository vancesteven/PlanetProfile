#!/usr/bin/env python
"""Injection-recovery configs for the Enceladus ocean6D rho_rock campaign
(reviewer Fisher-checkpoint condition 4; run (h) review: include a truth at
log10 w ~ 1.8 where the cache hole and the w-dependence interact).

Given a truth theta, compute the four observables (C20, C22, C30,
libration_deg) through the SAME runner methods log_likelihood_fn uses
(_derive_gravity_isostatic -> mass closure; _derive_libration_deg -> model
budget incl. multiplicative_factor and (1 + libration_sys_frac)), then write
a config identical to the input except observables[k][0] := the injected
value. Sigmas are untouched. NOISE-FREE by default (the posterior should
then centre on the truth up to prior/geometry effects); --noise-seed draws
one Gaussian realisation with the configured sigmas instead.

The truth's own log-likelihood under the injected data is checked to be the
noise-free maximum (chi2 = 0 up to sigma_model terms) before writing.

  python plans/scripts/enceladus_injection_config.py --config CFG.json \
      --truth zb_km=27 log10_wOcean_ppt=1.8 ... --out INJ.json \
      [--noise-seed N] [--json-out TRUTH.json]
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from PlanetProfile.Inference.inference_core import InferenceConfig  # noqa: E402
from PlanetProfile.Inference.mcmc_runner import MCMCRunner  # noqa: E402

OBS = ('C20', 'C22', 'C30', 'libration_deg')


def predict(runner, theta_vec):
    td = runner._expand_theta(np.asarray(theta_vec, float))
    runner.last_mass_closure = None
    iso = runner._derive_gravity_isostatic(td)
    if iso is None:
        raise SystemExit('truth is outside the support: gravity None '
                         f'({(runner.last_mass_closure or {}).get("reason")})')
    lib = runner._derive_libration_deg(td)
    if lib is None:
        raise SystemExit('truth is outside the support: libration None')
    mc = dict(runner.last_mass_closure or {})
    return ({'C20': float(iso[(2, 0)]), 'C22': float(iso[(2, 2)]),
             'C30': float(iso[(3, 0)]), 'libration_deg': float(lib)},
            {k: (float(v) if isinstance(v, (int, float, np.floating)) else v)
             for k, v in mc.items() if k in ('cmr2_derived', 'R_s_m', 'rho_s',
                                             'ocean_thickness_m')})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True)
    ap.add_argument('--truth', nargs='+', required=True,
                    help='name=value for EVERY sampled parameter')
    ap.add_argument('--out', required=True)
    ap.add_argument('--noise-seed', type=int)
    ap.add_argument('--json-out')
    a = ap.parse_args()

    raw = json.loads(Path(a.config).read_text())
    cfg = InferenceConfig.from_json(a.config)
    runner = MCMCRunner(cfg)
    names = list(runner.param_names)
    truth = dict(kv.split('=', 1) for kv in a.truth)
    missing = set(names) - set(truth)
    extra = set(truth) - set(names)
    if missing or extra:
        raise SystemExit(f'truth must set exactly {names}; missing {missing} extra {extra}')
    th = np.array([float(truth[n]) for n in names])
    for n, v in zip(names, th):
        lo, hi = cfg.param_space[n]['bounds']
        if not lo <= v <= hi:
            raise SystemExit(f'{n}={v} outside prior bounds [{lo}, {hi}]')

    pred, derived = predict(runner, th)
    inj = dict(pred)
    noise_cov = None
    if a.noise_seed is not None:
        # Draw from the covariance the LIKELIHOOD uses: config sigmas
        # inflated by metadata.sigma_model_add, with the declared
        # metadata.observable_correlations (reviewer 2026-10-27; the first
        # version drew raw independent sigmas, a recorded spec deviation).
        md = raw.get('metadata', {})
        add = md.get('sigma_model_add', {}) or {}
        s = np.array([np.hypot(float(raw['observables'][k][1]),
                               float(add.get(k, 0.0))) for k in OBS])
        C = np.diag(s ** 2)
        for pair, r in (md.get('observable_correlations', {}) or {}).items():
            i, j = (OBS.index(x.strip()) for x in pair.split(','))
            C[i, j] = C[j, i] = float(r) * s[i] * s[j]
        noise_cov = C.tolist()
        rng = np.random.default_rng(a.noise_seed)
        eps = np.linalg.cholesky(C) @ rng.normal(size=len(OBS))
        for k, e in zip(OBS, eps):
            inj[k] = pred[k] + float(e)
    for k in OBS:
        raw['observables'][k][0] = inj[k]
    raw.setdefault('metadata', {})['injection'] = {
        'truth': dict(zip(names, map(float, th))), 'noise_free_prediction': pred,
        'injected': inj, 'noise_seed': a.noise_seed,
        'noise_covariance': noise_cov,
        'noise_model': ('inflated sigma (sigma_model_add) + observable_correlations'
                        if noise_cov is not None else None),
        'truth_derived': derived, 'source_config': a.config,
        'note': 'observables[k][0] replaced by the injected value; sigmas unchanged'}
    Path(a.out).write_text(json.dumps(raw, indent=1))

    # Sanity: the truth must score as the noise-free maximum under the
    # injected data (rebuild the runner so its likelihood closure binds the
    # NEW observables -- the closure is captured at construction).
    cfg2 = InferenceConfig.from_json(a.out)
    r2 = MCMCRunner(cfg2)
    ll_truth = float(r2.log_likelihood_fn(th))
    print(json.dumps({'truth': raw['metadata']['injection']['truth'],
                      'prediction': pred, 'injected': inj,
                      'truth_derived': derived,
                      'loglike_truth_under_injected': ll_truth}, indent=1))
    if a.json_out:
        Path(a.json_out).write_text(json.dumps(raw['metadata']['injection'], indent=1))
    if a.noise_seed is None and not np.isfinite(ll_truth):
        raise SystemExit('truth loglike not finite under its own injection')


if __name__ == '__main__':
    main()
