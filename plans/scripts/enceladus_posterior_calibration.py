#!/usr/bin/env python
"""Local calibration check of the Enceladus run (a) posterior (final-acceptance
review 2026-10-27): noise-free injections test BIAS only; this tests
CALIBRATION in the data neighbourhood.

1. Rebuild the exact per-sample predictions P = (C20, C22, C30, libration)
   for every run (a) posterior sample: C20/C22/C30 from the stored
   isostatic_derived arrays, libration recomputed through
   MCMCRunner._derive_libration_deg (the likelihood's own path).
   Assert the reconstructed log-likelihood -- config sigmas inflated by
   metadata.sigma_model_add, rho(C20,C22) from
   metadata.observable_correlations -- matches the stored log_likelihoods.
2. N replicates: draw a truth sample j from posterior (a) by weight, draw a
   data realisation x = P_j + Gaussian noise (full inflated covariance),
   form the exact importance-reweighted posterior w_i ∝ w_i^(a) L(x|i)/L(x_obs|i)
   ... equivalently log w_i = log w_i^(a) - loglike_obs_i + loglike_x_i, with
   the truth sample excluded (leave-one-out), and record the truth's rank
   (posterior CDF at the truth) for zb, rho_s, C2.
   Calibrated => ranks ~ U(0,1): report 68% coverage and KS distance.
3. Control: the same with P replaced by its linear least-squares fit in
   theta (an exactly Gaussian forward model on the same samples).

  python plans/scripts/enceladus_posterior_calibration.py [--n 1000] \
      [--seed 4] [--json-out FILE]
"""
import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from PlanetProfile.Inference.inference_core import InferenceConfig  # noqa: E402
from PlanetProfile.Inference.mcmc_runner import MCMCRunner  # noqa: E402

D = REPO / 'validation_reports/enceladus_isostasy/posterior_a'
OBS = ('C20', 'C22', 'C30', 'libration_deg')
GATED = ('zb_km', 'rho_rock_kgm3', 'compensation_C2')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=1000)
    ap.add_argument('--seed', type=int, default=4)
    ap.add_argument('--json-out')
    a = ap.parse_args()

    r = pickle.load(open(D / 'rep1_result.pkl', 'rb'))
    w = np.asarray(r.weights, float)
    w /= w.sum()
    raw = json.loads((D / 'rep1_config.json').read_text())
    run = MCMCRunner(InferenceConfig.from_json(str(D / 'rep1_config.json')))
    S = np.asarray(r.samples, float)
    names = list(r.param_names)
    lib = np.array([np.nan if (v := run._derive_libration_deg(run._expand_theta(s)))
                    is None else v for s in S])
    d = r.metadata['isostatic_derived']
    P = np.c_[d['C20'], d['C22'], d['C30'], lib]
    assert np.all(np.isfinite(P)), 'non-finite prediction in posterior (a)'

    md = raw['metadata']
    obs = np.array([raw['observables'][k][0] for k in OBS])
    sig = np.array([raw['observables'][k][1] for k in OBS])
    sm = np.array([md['sigma_model_add'][k] for k in OBS])
    s = np.sqrt(sig ** 2 + sm ** 2)
    rho = float(md['observable_correlations']['C20,C22'])
    C = np.diag(s ** 2)
    C[0, 1] = C[1, 0] = rho * s[0] * s[1]
    Ci = np.linalg.inv(C)
    Lc = np.linalg.cholesky(C)

    def ll(PP, o):
        r_ = PP - o
        return -0.5 * np.einsum('ni,ij,nj->n', r_, Ci, r_)

    ll_obs = ll(P, obs)
    recon = float(np.max(np.abs(ll_obs - np.asarray(r.log_likelihoods, float))))
    print(f'reconstructed loglike vs stored: max abs diff {recon:.3e}')
    assert recon < 1e-6, 'likelihood reconstruction does not match the sampler'

    X = np.c_[np.ones(len(S)), S]
    B, *_ = np.linalg.lstsq(X, P, rcond=None)
    PL = X @ B
    nonlin = (np.sqrt(np.mean((P - PL) ** 2, 0)) / s).tolist()
    gi = [names.index(n) for n in GATED]
    out = {'n': a.n, 'seed': a.seed, 'loglike_recon_max_abs_diff': recon,
           'inflated_sigma': dict(zip(OBS, s.tolist())), 'rho_C20_C22': rho,
           'nonlinearity_rms_over_sigma': dict(zip(OBS, nonlin)), 'runs': {}}

    for tag, PP in (('nonlinear', P), ('linear_control', PL)):
        rng = np.random.default_rng(a.seed)
        base = np.log(w) - ll(PP, obs)
        R, E = [], []
        for _ in range(a.n):
            j = rng.choice(len(S), p=w)
            x = PP[j] + Lc @ rng.normal(size=4)
            lw = base + ll(PP, x)
            lw[j] = -np.inf
            lw -= lw.max()
            ww = np.exp(lw)
            ww /= ww.sum()
            E.append(1.0 / np.sum(ww ** 2))
            R.append([np.sum(ww[S[:, i] <= S[j, i]]) for i in gi])
        R = np.array(R)
        res = {'ess_p5': float(np.percentile(E, 5)), 'ess_min': float(np.min(E))}
        for k, n in enumerate(GATED):
            rk = np.sort(R[:, k])
            res[n] = {'cov68': float(np.mean((rk > .16) & (rk < .84))),
                      'cov95': float(np.mean((rk > .025) & (rk < .975))),
                      'rank_mean': float(rk.mean()),
                      'ks': float(np.max(np.abs(rk - np.arange(1, a.n + 1) / a.n)))}
        out['runs'][tag] = res
        print(tag, json.dumps(res))
    out['ks_crit_5pct'] = 1.36 / np.sqrt(a.n)
    print('KS 5% critical', round(out['ks_crit_5pct'], 4),
          '; nonlinearity/sigma', np.round(nonlin, 3))
    if a.json_out:
        Path(a.json_out).write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
