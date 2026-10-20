"""R1/R2 prereg reference-only calibration (reviewer 2026-10-20 edits 1-2).
Reads ONLY the n_eff=8000 reference seed pkls. No corrected-side quantity."""
import json, pickle, sys
import numpy as np
sys.path.insert(0, '.')
from PlanetProfile.Inference.validate_sbi import _run_crosscheck, _normalize_weights, _weighted_median_1d
SRC = 'validation_reports/nh3_diagnosis/matched_reference_neff8000'
L = 1270
C9 = ["log10_eta_Ih", "log10_eta_III", "log10_eta_V", "log10_eta_VI", "log10_eta_sil", "log10_zeta", "alpha"]
seeds, S, W = [], [], []
for sd in (72, 172, 272):
    o = pickle.load(open(f'{SRC}/nh3_reference_seed{sd}_neff8000.pkl', 'rb'))
    names = list(o.param_names)
    w = np.asarray(o.weights, float)
    S.append(np.asarray(o.samples, float)); W.append(w / w.sum() / 3.0)
    seeds.append((sd, np.asarray(o.samples, float), w / w.sum()))
S = np.concatenate(S); W = np.concatenate(W); W = W / W.sum()
ess_ref = 1.0 / np.sum(W ** 2)

# Per-seed medians -> s_seed (as before) and pooled median
s_seed = {}
for k in C9:
    j = names.index(k)
    meds = [_weighted_median_1d(s_[:, j], w_) for _, s_, w_ in seeds]
    s_seed[k] = float(np.std(meds, ddof=1))
m_ref = {k: _weighted_median_1d(S[:, names.index(k)], W) for k in C9}

# SE_corr: std of the median of L weight-resampled pooled-reference draws, B=1000, RandomState(72)
rng = np.random.RandomState(72)
B = 1000
meds = np.empty((B, len(C9)))
cdf = np.cumsum(W)
for b in range(B):
    idx = np.searchsorted(cdf, rng.random_sample(L) * cdf[-1])
    meds[b] = [np.median(S[idx, names.index(k)]) for k in C9]
se_corr = {k: float(meds[:, i].std(ddof=1)) for i, k in enumerate(C9)}
tol = {k: 2.0 * float(np.hypot(s_seed[k], se_corr[k])) for k in C9}
null_c9_pass = {k: float(np.mean(np.abs(meds[:, i] - m_ref[k]) <= tol[k])) for i, k in enumerate(C9)}

# Null-fail table: 10 reference self-draw replicates of size L through _run_crosscheck, stock vs matched floor
rng2 = np.random.RandomState(1072)
table = {'stock': {n: 0 for n in names}, 'matched': {n: 0 for n in names}}
for r in range(10):
    idx = np.searchsorted(cdf, rng2.random_sample(L) * cdf[-1])
    null = S[idx]
    for mode, flag in (('stock', False), ('matched', True)):
        res = _run_crosscheck(names, S, W, null, seed=72, matched_n_floor=flag)
        for p in res['per_parameter']:
            if not (p['mean_pass'] and p['sigma_pass'] and p['shape_pass']):
                table[mode][p['param']] += 1
    print('replicate', r, 'done', flush=True)

out = dict(
    kind='r1_r2_c9_prereg_bands',
    note=('Computed from the n_eff=8000 REFERENCE seed pkls ONLY; no corrected quantity read. '
          'C9 rule (reviewer-edited 2026-10-20): PASS-by-C9 iff sigma_pass AND '
          '|m_corr - m_ref_pooled| <= tol = 2*sqrt(s_seed^2 + SE_corr^2).'),
    L=L, ess_ref_pooled=float(ess_ref),
    c9={k: dict(s_seed_ddof1=s_seed[k], se_corr_L1270=se_corr[k], tol=tol[k],
                m_ref_pooled=float(m_ref[k]), null_pass_rate_B1000=null_c9_pass[k]) for k in C9},
    se_corr_method='std of median of L=1270 weight-resampled pooled-reference draws, B=1000, RandomState(72)',
    null_fail_table=dict(replicates=10, sample_size=L, rng='RandomState(1072)',
                         crosscheck_seed=72, fails_per_param=table),
)
json.dump(out, open('validation_reports/titan_freegrav_nh3_1m/is_correction/r1_r2_prereg_bands.json', 'w'), indent=1)
print(json.dumps({k: {kk: round(vv, 4) for kk, vv in v.items()} for k, v in out['c9'].items()}, indent=0))
print('ESS_ref', round(ess_ref)); print(json.dumps(table))
