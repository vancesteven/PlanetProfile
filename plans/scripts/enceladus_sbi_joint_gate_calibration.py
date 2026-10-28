#!/usr/bin/env python
"""Joint-gate calibration for the Enceladus ocean-branch SBI preregistration
(reviewer APPROVE WITH EDITS 2026-10-27, MAJOR-3). REFERENCE DRAWS ONLY --
no flow exists yet, so nothing here can be tuned against a flow result.

Joint gate (G2-ii): energy distance between flow and reference in whitened
(zb, rho_s, C2) space, gated at k x the reference self-distance p99 floor at
matched n (N_flow = 40,000 vs reference ESS 32,273), plus
|corr_flow(rho_s, C2) - corr_ref| <= 0.05.

This script:
  1. builds the pooled run (a) reference (seeds 75 + 2027, treatment A);
  2. computes the self-distance null: energy distance between two INDEPENDENT
     weight-resamples of the reference at the matched sizes, R times;
  3. computes the corr null the same way;
  4. MUTATION sanity checks (required by the reviewer): the gate must FAIL on
     a reference sample whose rho_s and C2 columns are independently
     permuted, and is reported for a subtler partial (20%) permutation.

Energy distance uses an m-point subsample per side (O(m^2) cost); m and R
are recorded. Output: JSON with floor, k, and mutation results.

  python plans/scripts/enceladus_sbi_joint_gate_calibration.py --out FILE
"""
import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
D = REPO / 'validation_reports/enceladus_isostasy'
JOINT = ('zb_km', 'rho_rock_kgm3', 'compensation_C2')


def pooled_reference():
    S, W = [], []
    for p in ('posterior_a/rep1_result.pkl', 'posterior_a_seed2027/rep1_result.pkl'):
        r = pickle.load(open(D / p, 'rb'))
        w = np.asarray(r.weights, float)
        S.append(np.asarray(r.samples, float))
        W.append(w / w.sum() / 2.0)
    W = np.concatenate(W)
    return np.concatenate(S), W / W.sum(), list(r.param_names)


def wcorr(x, y, w):
    mx, my = np.average(x, weights=w), np.average(y, weights=w)
    return float(np.average((x - mx) * (y - my), weights=w)
                 / np.sqrt(np.average((x - mx) ** 2, weights=w)
                           * np.average((y - my) ** 2, weights=w)))


def energy_distance(a, b, rng, m):
    a = a[rng.choice(len(a), min(m, len(a)), replace=False)]
    b = b[rng.choice(len(b), min(m, len(b)), replace=False)]

    def d(u, v):
        return float(np.mean(np.linalg.norm(u[:, None, :] - v[None, :, :], axis=2)))
    return 2.0 * d(a, b) - d(a, a) - d(b, b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out', required=True)
    ap.add_argument('--n-flow', type=int, default=40000)
    ap.add_argument('--reps', type=int, default=100)
    ap.add_argument('--m', type=int, default=2000)
    ap.add_argument('--seed', type=int, default=20261027)
    a = ap.parse_args()

    S, W, names = pooled_reference()
    ess = 1.0 / np.sum(W ** 2)
    n_ref = int(np.floor(ess))
    ix = [names.index(n) for n in JOINT]
    X = S[:, ix]
    mu = np.average(X, axis=0, weights=W)
    C = np.cov(X.T, aweights=W)
    Lw = np.linalg.cholesky(np.linalg.inv(C))
    Z = (X - mu) @ Lw
    Linv = np.linalg.inv(Lw)
    i_rs, i_c2 = JOINT.index('rho_rock_kgm3'), JOINT.index('compensation_C2')
    c_ref = wcorr(X[:, i_rs], X[:, i_c2], W)
    rng = np.random.default_rng(a.seed)

    def draw(n):
        return Z[rng.choice(len(Z), size=n, p=W)]

    def corr_of(z):
        x = z @ Linv + mu
        return float(np.corrcoef(x[:, i_rs], x[:, i_c2])[0, 1])

    null_e, null_c = [], []
    for _ in range(a.reps):
        f = draw(a.n_flow)
        null_e.append(energy_distance(f, draw(n_ref), rng, a.m))
        null_c.append(abs(corr_of(f) - c_ref))
    p99_e = float(np.percentile(null_e, 99))
    p99_c = float(np.percentile(null_c, 99))
    K = 2.0
    gate_e = K * p99_e

    def mutate(frac):
        # Permute in PHYSICAL space (the whitening already decorrelates
        # rho_s and C2, so permuting whitened columns is a near no-op --
        # the first version of this check did exactly that and could not
        # fail). Then re-whiten with the REFERENCE transform.
        x = draw(a.n_flow) @ Linv + mu
        k = rng.random(len(x)) < frac
        x[k, i_rs] = rng.permutation(x[k, i_rs])
        x[k, i_c2] = rng.permutation(x[k, i_c2])
        return (x - mu) @ Lw

    mut = {}
    for frac in (1.0, 0.2):
        e, c = [], []
        for _ in range(20):
            z = mutate(frac)
            e.append(energy_distance(z, draw(n_ref), rng, a.m))
            c.append(abs(corr_of(z) - c_ref))
        mut[str(frac)] = {
            'energy_min': float(min(e)), 'energy_median': float(np.median(e)),
            'energy_fail_rate': float(np.mean(np.array(e) > gate_e)),
            'dcorr_min': float(min(c)),
            'dcorr_fail_rate': float(np.mean(np.array(c) > 0.05)),
            'joint_fail_rate': float(np.mean((np.array(e) > gate_e)
                                             | (np.array(c) > 0.05)))}

    out = {
        'kind': 'enceladus_sbi_joint_gate_calibration',
        'note': 'REFERENCE draws only; no flow. Fixes k and the floor recipe before any training.',
        'reference': 'pooled run (a) seeds 75 + 2027, treatment A',
        'ess_ref': float(ess), 'n_ref_draws': n_ref, 'n_flow': a.n_flow,
        'whitened_params': list(JOINT), 'energy_subsample_m': a.m,
        'null_reps': a.reps, 'rng_seed': a.seed,
        'corr_ref_rho_s_C2': c_ref,
        'energy_null': {'median': float(np.median(null_e)), 'p99': p99_e},
        'dcorr_null': {'median': float(np.median(null_c)), 'p99': p99_c},
        'k': K, 'energy_gate': gate_e, 'dcorr_gate': 0.05,
        'null_false_fail_rate': {
            'energy': float(np.mean(np.array(null_e) > gate_e)),
            'dcorr': float(np.mean(np.array(null_c) > 0.05))},
        'mutation': mut,
        'mutation_sanity_pass': bool(mut['1.0']['joint_fail_rate'] == 1.0),
    }
    Path(a.out).write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
