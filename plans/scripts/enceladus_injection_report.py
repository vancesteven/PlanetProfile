#!/usr/bin/env python
"""Injection-recovery + rho_ice U[880,940] report for the Enceladus ocean6D
rho_rock campaign.

Recovery criterion (stated in plans/STATUS.md BEFORE results, 2026-10-27):
each truth lies inside the 16-84% interval for zb, rho_s and C2, or the
reason is understood (w and rho_ice are prior-dominated, so their
posteriors revert toward the prior instead of recovering the truth).

For each injection also reports the truth's position as a posterior
percentile and its offset in posterior sd, for every parameter and for
derived C/MR^2. For ri880 reports the shift vs run (a) for every parameter.

  python plans/scripts/enceladus_injection_report.py [--json-out FILE]
"""
import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

D = Path('validation_reports/enceladus_isostasy')
INJ = ('inj_w18', 'inj_ri918', 'inj_ri933')
GATED = ('zb_km', 'rho_rock_kgm3', 'compensation_C2')


def load(name):
    r = pickle.load(open(D / name / 'rep1_result.pkl', 'rb'))
    w = np.asarray(r.weights, float)
    return r, w / w.sum()


def wq(x, w, q):
    i = np.argsort(x)
    c = np.cumsum(w[i])
    c /= c[-1]
    return np.interp(q, c, x[i])


def wpct(x, w, v):
    return float(np.sum(w[x <= v]) / np.sum(w))


def summ(x, w):
    m = float(np.average(x, weights=w))
    s = float(np.sqrt(np.average((x - m) ** 2, weights=w)))
    return dict(q16=float(wq(x, w, .16)), q50=float(wq(x, w, .5)),
                q84=float(wq(x, w, .84)), mean=m, sd=s)


def cmr2(r, w):
    c = np.asarray(r.metadata['isostatic_derived']['cmr2_derived'], float)
    f = np.isfinite(c)
    return c[f], w[f] / w[f].sum()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--json-out')
    a = ap.parse_args()
    out = {'criterion': ('truth inside 16-84% for zb, rho_s, C2 (or reason '
                         'understood); stated before results 2026-10-27')}
    A, wa = load('posterior_a')
    names = list(A.param_names)

    for v in INJ:
        r, w = load(f'posterior_{v}')
        assert list(r.param_names) == names
        cfg = json.loads((D / f'posterior_{v}' / 'rep1_config.json').read_text())
        inj = cfg['metadata']['injection']
        truth = inj['truth']
        rows = {}
        for j, n in enumerate(names):
            x = r.samples[:, j]
            s = summ(x, w)
            t = float(truth[n])
            s.update(truth=t, truth_pct=wpct(x, w, t),
                     truth_offset_sd=(t - s['mean']) / s['sd'],
                     inside_16_84=bool(s['q16'] <= t <= s['q84']),
                     prior_sd_ratio=None)
            rows[n] = s
        c, wc = cmr2(r, w)
        tc = float(inj['truth_derived']['cmr2_derived'])
        sc = summ(c, wc)
        sc.update(truth=tc, truth_pct=wpct(c, wc, tc),
                  truth_offset_sd=(tc - sc['mean']) / sc['sd'],
                  inside_16_84=bool(sc['q16'] <= tc <= sc['q84']))
        rows['cmr2_derived'] = sc
        gated_pass = all(rows[n]['inside_16_84'] for n in GATED)
        out[v] = dict(truth=truth, ess=float(1 / np.sum(w ** 2)),
                      log_Z=r.metadata.get('log_Z'),
                      max_loglike=float(np.max(r.log_likelihoods)),
                      params=rows, gated_inside_16_84=gated_pass)

    r, w = load('posterior_ri880')
    rows = {}
    for j, n in enumerate(names):
        sa, sb = summ(A.samples[:, j], wa), summ(r.samples[:, j], w)
        rows[n] = dict(a=sa, ri880=sb,
                       shift_median=sb['q50'] - sa['q50'],
                       shift_over_sd_a=(sb['mean'] - sa['mean']) / sa['sd'],
                       sd_ratio=sb['sd'] / sa['sd'])
    ca, wca = cmr2(A, wa)
    cb, wcb = cmr2(r, w)
    sa, sb = summ(ca, wca), summ(cb, wcb)
    rows['cmr2_derived'] = dict(a=sa, ri880=sb, shift_median=sb['q50'] - sa['q50'],
                                shift_over_sd_a=(sb['mean'] - sa['mean']) / sa['sd'],
                                sd_ratio=sb['sd'] / sa['sd'])
    out['ri880'] = dict(ess=float(1 / np.sum(w ** 2)), log_Z=r.metadata.get('log_Z'),
                        log_Z_a=A.metadata.get('log_Z'), params=rows)

    for v in INJ:
        o = out[v]
        print(f"== {v}: ESS {o['ess']:.0f}  max logL {o['max_loglike']:.3f}  "
              f"gated inside 16-84%: {o['gated_inside_16_84']}")
        for n, s in o['params'].items():
            print('   %-20s truth %10.5g  post %s  pct %.2f  off %+.2f sd  %s' % (
                n, s['truth'], np.round([s['q16'], s['q50'], s['q84']], 5),
                s['truth_pct'], s['truth_offset_sd'],
                'IN' if s['inside_16_84'] else 'out'))
    o = out['ri880']
    print(f"== ri880 vs (a): ESS {o['ess']:.0f}  logZ {o['log_Z']:.3f} vs {o['log_Z_a']:.3f}")
    for n, s in o['params'].items():
        print('   %-20s median %10.5g -> %10.5g  shift/sd %+.3f  sd ratio %.3f' % (
            n, s['a']['q50'], s['ri880']['q50'], s['shift_over_sd_a'], s['sd_ratio']))
    if a.json_out:
        Path(a.json_out).write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
