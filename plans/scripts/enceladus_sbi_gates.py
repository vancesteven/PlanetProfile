#!/usr/bin/env python
"""Gate evaluation (G1-G5 + query-time guard) for the Enceladus ocean-branch
SBI artifact, EXACTLY as preregistered in
validation_reports/enceladus_isostasy/sbi/DESIGN_PREREGISTRATION.md
(commit 08659c67). No threshold is a CLI option; --smoke only reduces SIZES
(and its verdict is meaningless). Never tune after a result: a FAIL is
recorded and escalated.

  KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/enceladus_sbi_gates.py \
      --artifact ~/scratch/enc_sbi/enceladus_ocean6D_rhorock_posterior_1m.pt \
      --out validation_reports/enceladus_isostasy/sbi/gates_report.json

libomp: torch (sbi) and the numba forward model share ONE process here (as in
is_correction_validate.py / nh3_r1r2_gates.py); only KMP_DUPLICATE_LIB_OK=TRUE
is required. MCMCRunner (numba) is imported/built BEFORE torch is touched.
"""
import os
os.environ.setdefault('KMP_DUPLICATE_LIB_OK', 'TRUE')
os.environ.setdefault('OMP_NUM_THREADS', '2')

import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from PlanetProfile.Inference.inference_core import InferenceConfig  # noqa: E402
from PlanetProfile.Inference.mcmc_runner import MCMCRunner  # noqa: E402
from PlanetProfile.Inference.forward_models import (  # noqa: E402
    apply_parameters, UnservableSampleError)
from PlanetProfile.Inference.validate_sbi import (  # noqa: E402
    _run_crosscheck, _weighted_median_1d, _weighted_mean_std)
from enceladus_sbi_joint_gate_calibration import (  # noqa: E402
    pooled_reference, wcorr, energy_distance, JOINT)

D = REPO / 'validation_reports/enceladus_isostasy'
CFG = REPO / 'PlanetProfile/Inference/configs/enceladus_cassini_isostasy_ocean6D_rhorock.json'
CACHE = '/Users/svance/scratch/enceladus_zbw_production_cache.pkl'
SCR = Path('/Users/svance/scratch/enc_sbi')
CALIB = D / 'sbi/joint_gate_calibration.json'

# ---- preregistered constants (DESIGN_PREREGISTRATION.md) -------------------
G1_SEED = 75
G2_MEDIAN_SD = 0.25
G2_KS = 0.10
G2_OBS_SIGMA = 0.1
G2_DCORR = 0.05
G2_C2_DP = 0.02
G2_C2_THRESH = 0.95
G3_PVAL = 0.05
G3_COV_LO, G3_COV_HI = 0.65, 0.71
G4_MED_SD = 0.3
G4_SD_RATIO = (0.7, 1.4)
G5_MEAN = 0.01
G5_P95 = 0.05
GUARD_FRAC = 0.01
ZB_THIN_KM = 15.0
ENERGY_M = 2000
ENERGY_RNG = 20261027
INJ = ['inj_w18', 'inj_ri918', 'inj_ri933', 'inj_noise1', 'inj_noise2', 'inj_noise3']
G4_PARAMS = ['zb_km', 'rho_rock_kgm3', 'compensation_C2']

T0 = time.time()
TIMES = {}


def log(*a):
    print(f'[{time.time() - T0:7.0f}s]', *a, flush=True)


class Timer:
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.t = time.time()

    def __exit__(self, *e):
        TIMES[self.name] = time.time() - self.t
        log(f'{self.name} done in {TIMES[self.name]:.1f}s')


# ---------------------------------------------------------------- helpers ---
def pool(paths):
    """Pooled treatment A: each seed's weights w/sum(w)/n_seeds."""
    S, W = [], []
    for p in paths:
        r = pickle.load(open(D / p / 'rep1_result.pkl', 'rb'))
        w = np.asarray(r.weights, float)
        S.append(np.asarray(r.samples, float))
        W.append(w / w.sum() / len(paths))
        names = list(r.param_names)
    W = np.concatenate(W)
    return np.concatenate(S), W / W.sum(), names


def wquantile(v, w, q):
    o = np.argsort(v)
    c = np.cumsum(w[o])
    c /= c[-1]
    return float(v[o][min(np.searchsorted(c, q), len(v) - 1)])


def wsd(v, w):
    m = np.average(v, weights=w)
    return float(np.sqrt(np.average((v - m) ** 2, weights=w)))


def ks_weighted(a, wa, b, wb):
    """Two-sample KS between weighted samples (exact on the pooled grid)."""
    wa = np.asarray(wa, float) / np.sum(wa)
    wb = np.asarray(wb, float) / np.sum(wb)
    oa, ob = np.argsort(a), np.argsort(b)
    a_s, b_s = a[oa], b[ob]
    ca, cb = np.cumsum(wa[oa]), np.cumsum(wb[ob])
    grid = np.concatenate([a_s, b_s])
    fa = np.where(np.searchsorted(a_s, grid, side='right') > 0,
                  ca[np.maximum(np.searchsorted(a_s, grid, side='right') - 1, 0)], 0.0)
    fb = np.where(np.searchsorted(b_s, grid, side='right') > 0,
                  cb[np.maximum(np.searchsorted(b_s, grid, side='right') - 1, 0)], 0.0)
    return float(np.max(np.abs(fa - fb)))


def bh_adjust(p):
    p = np.asarray(p, float)
    m = len(p)
    o = np.argsort(p)
    adj = np.empty(m)
    run = 1.0
    for rank in range(m, 0, -1):
        i = o[rank - 1]
        run = min(run, p[i] * m / rank)
        adj[i] = run
    return adj


def is_hole(runner, th):
    """Hole cell = apply_parameters raises UnservableSampleError OR the
    ocean zb x w struct is None. Same two tests, same order, as
    enceladus_prior_pushforward.classify_and_derive (without the expensive
    gravity derivation, which does not affect the hole label)."""
    td = runner._expand_theta(th)
    try:
        apply_parameters(td, runner.structure_data)
    except UnservableSampleError:
        return True
    if 'zb_km' in td and runner._ocean_zbw_struct_for_theta(td) is None:
        return True
    return False


def hole_mask(runner, TH):
    return np.array([is_hole(runner, t) for t in TH], bool)


def push_noise_free(runner, TH, ref_names=None):
    """Noise-free observables via generate_sbi_dataset(theta_override=...)
    (same path as the training x0) PLUS C/MR^2 via _expand_theta +
    _derive_gravity_isostatic (runner.last_mass_closure['cmr2_derived']).
    Returns dict name -> array (NaN where unservable)."""
    obs = list(runner.config.observables.keys())
    TH = np.asarray(TH, float)
    _, X = runner.generate_sbi_dataset(theta_override=TH, drop_nonfinite=False,
                                       imag_convention='abs')
    X = np.asarray(X, float)
    if X.shape != (len(TH), len(obs)):
        raise RuntimeError(f'pushforward shape {X.shape} != {(len(TH), len(obs))}: missing column')
    cm = np.full(len(TH), np.nan)
    for i, th in enumerate(TH):
        td = runner._expand_theta(th)
        runner.last_mass_closure = None
        try:
            iso = runner._derive_gravity_isostatic(td)
        except UnservableSampleError:
            iso = None
        if iso is not None and runner.last_mass_closure:
            cm[i] = float(runner.last_mass_closure.get('cmr2_derived', np.nan))
    out = {o: X[:, j] for j, o in enumerate(obs)}
    out['cmr2'] = cm
    return out


def draw(srun, xdict, n, seed, reject=True):
    return srun.sample_posterior(xdict, n_samples=n, seed=seed, reject_outside_prior=reject)


def reorder(arr, from_names, to_names):
    if list(from_names) == list(to_names):
        return arr
    assert set(from_names) == set(to_names), (from_names, to_names)
    return arr[:, [list(from_names).index(n) for n in to_names]]


# ------------------------------------------------------------------- main ---
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--artifact', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--smoke', action='store_true',
                    help='reduced SIZES only (verdict meaningless); refuses to write into validation_reports/')
    ap.add_argument('--cache', default=CACHE)
    a = ap.parse_args()
    outp = Path(a.out).resolve()
    if a.smoke and str(REPO / 'validation_reports') in str(outp):
        raise SystemExit('--smoke must not write into validation_reports/')

    S = dict(n_post=40000, n_push=6000, g3_n=2000, g3_draws=1000, g4_n=40000,
             guard_n=10000, thin_n=None, leak_n=10000)
    if a.smoke:
        S = dict(n_post=2000, n_push=300, g3_n=50, g3_draws=200, g4_n=2000,
                 guard_n=1000, thin_n=30, leak_n=1000)

    calib = json.load(open(CALIB))
    K_E = calib['k']
    E_GATE = calib['energy_gate']
    CORR_REF_JSON = calib['corr_ref_rho_s_C2']

    raw = json.loads(CFG.read_text())
    raw['structure_cache_path'] = a.cache
    cfg = InferenceConfig(**raw)
    runner = MCMCRunner(cfg)                      # numba side first
    names = list(runner.param_names)
    obs = list(cfg.observables.keys())
    xobs = {k: float(v[0]) for k, v in cfg.observables.items()}
    sig = np.array([np.hypot(raw['observables'][k][1],
                             (raw['metadata'].get('sigma_model_add') or {}).get(k, 0.0))
                    for k in obs])
    bounds = {k: v['bounds'] for k, v in raw['param_space'].items()}

    from PlanetProfile.Inference.sbi_runner import SBIRunner   # torch after numba
    import torch
    torch.set_num_threads(int(os.environ.get('OMP_NUM_THREADS', '2')))
    srun = SBIRunner.load_artifact(a.artifact)
    assert list(srun.param_names) == names, (srun.param_names, names)
    assert list(srun.obs_names) == obs, (srun.obs_names, obs)
    rep = {'kind': 'enceladus_sbi_gates', 'smoke': bool(a.smoke), 'artifact': str(a.artifact),
           'config': str(CFG), 'cache': a.cache, 'sizes': S, 'preregistration_commit': '08659c67',
           'artifact_config_hash': (srun.artifact_meta or {}).get('config_hash'),
           'param_names': names, 'obs_names': obs, 'sigma_infl': dict(zip(obs, sig.tolist())),
           'x_obs': xobs}

    # -------- references (pooled treatment A) --------
    RS, RW, rnames = pooled_reference()
    assert set(rnames) == set(names)
    RS = reorder(RS, rnames, names)
    ess_ref = float(1.0 / np.sum(RW ** 2))
    n_ref = int(np.floor(ess_ref))
    ix = {n: names.index(n) for n in names}

    # =========================== G1 =====================================
    with Timer('G1'):
        FL = draw(srun, xobs, S['n_post'], seed=G1_SEED)
        cc = _run_crosscheck(names, RS, RW, FL, seed=G1_SEED, matched_n_floor=False)
        g1_pass = bool(cc['all_pass'])
        rep['G1'] = {'pass': g1_pass, 'n_ks': [p['ks_n'] for p in cc['per_parameter']],
                     'n_ks_unique': sorted({int(p['ks_n']) for p in cc['per_parameter']}),
                     'floor_ess': cc.get('floor_ess', cc['ess_mcmc']),
                     'ess_ref': cc['ess_mcmc'], 'n_flow': cc['n_sbi'],
                     'matched_n_floor': False, 'crosscheck': cc,
                     'note': 'n_ks=floor(ESS_ref)=32,273 only when n_post=40,000 (full run)'}

    # =========================== G2 =====================================
    with Timer('G2'):
        g2 = {'checks': {}}
        # --- parameter-level (full weighted reference vs full flow) ---
        jn = list(JOINT)
        ji = [ix[n] for n in jn]
        X = RS[:, ji]
        mu = np.average(X, axis=0, weights=RW)
        Cm = np.cov(X.T, aweights=RW)
        Lw = np.linalg.cholesky(np.linalg.inv(Cm))
        rng = np.random.default_rng(ENERGY_RNG)
        ref_draw = X[rng.choice(len(X), size=n_ref, p=RW)]          # matched n_ref
        Zf, Zr = (FL[:, ji] - mu) @ Lw, (ref_draw - mu) @ Lw
        energy = energy_distance(Zf, Zr, rng, ENERGY_M)
        i_rs, i_c2 = jn.index('rho_rock_kgm3'), jn.index('compensation_C2')
        corr_ref = wcorr(X[:, i_rs], X[:, i_c2], RW)
        corr_flow = float(np.corrcoef(FL[:, ji][:, i_rs], FL[:, ji][:, i_c2])[0, 1])
        dcorr = abs(corr_flow - corr_ref)
        joint_pass = bool(dcorr <= G2_DCORR and energy <= E_GATE)
        g2['joint'] = {'corr_flow': corr_flow, 'corr_ref': corr_ref,
                       'corr_ref_json': CORR_REF_JSON, 'dcorr': dcorr, 'dcorr_gate': G2_DCORR,
                       'energy': energy, 'energy_gate': E_GATE, 'k': K_E,
                       'n_flow': int(len(FL)), 'n_ref_matched': n_ref, 'm': ENERGY_M,
                       'rng': ENERGY_RNG, 'pass': joint_pass,
                       'note': 'ref transform whitening, ref drawn with replacement ~ W at n=floor(ESS)'}
        c2 = FL[:, ix['compensation_C2']]
        p_flow = float(np.mean(c2 > G2_C2_THRESH))
        p_ref = float(np.sum(RW[RS[:, ix['compensation_C2']] > G2_C2_THRESH]))
        g2['C2_gt_0.95'] = {'p_flow': p_flow, 'p_ref': p_ref, 'abs_diff': abs(p_flow - p_ref),
                            'gate': G2_C2_DP, 'pass': bool(abs(p_flow - p_ref) <= G2_C2_DP)}
        # out-of-box leakage (diagnostic draw without the prior-box rejection)
        try:
            FLo = draw(srun, xobs, S['leak_n'], seed=G1_SEED + 1, reject=False)
            out = np.zeros(len(FLo), bool)
            for n_, j in ix.items():
                lo, hi = bounds[n_]
                out |= (FLo[:, j] < lo) | (FLo[:, j] > hi)
            g2['out_of_box_leakage'] = {'fraction': float(out.mean()), 'n': int(len(FLo)),
                                        'per_param': {n_: float(np.mean((FLo[:, j] < bounds[n_][0]) | (FLo[:, j] > bounds[n_][1])))
                                                      for n_, j in ix.items()}}
        except Exception as e:
            g2['out_of_box_leakage'] = {'error': f'{type(e).__name__}: {e}'}

        # --- pushforward (noise-free, same theta_override path both sides) ---
        r2 = np.random.default_rng(ENERGY_RNG + 2)
        n_push = S['n_push']
        fi = r2.choice(len(FL), size=min(n_push, len(FL)), replace=False)
        ri = r2.choice(len(RS), size=min(n_push, len(RS)), replace=False)  # uniform subset, weights carried
        TF, TR, WR = FL[fi], RS[ri], RW[ri] / RW[ri].sum()
        t1 = time.time()
        PF = push_noise_free(runner, TF)
        PR = push_noise_free(runner, TR)
        g2['pushforward_seconds'] = time.time() - t1
        g2['n_push_flow'], g2['n_push_ref'] = int(len(TF)), int(len(TR))
        g2['ref_push_ess'] = float(1.0 / np.sum(WR ** 2))
        for k, v in PR.items():
            if not np.all(np.isfinite(v)):
                raise RuntimeError(f'G2: reference pushforward column {k} has non-finite entries '
                                   f'({np.sum(~np.isfinite(v))}); failing loudly')
        for k in obs + ['cmr2']:
            if k not in PF or k not in PR:
                raise RuntimeError(f'G2: missing pushforward column {k}')
        flow_fin = {k: np.isfinite(PF[k]) for k in PF}
        g2['flow_push_nonfinite_fraction'] = {k: float(1 - flow_fin[k].mean()) for k in PF}
        chk = {}
        # medians (<=0.25 posterior sd) for zb, rho_s (=rho_rock_kgm3), C2, C/MR2
        for nm, fv, rv, rw_, sd_ in (
                [(n_, FL[:, ix[n_]], RS[:, ix[n_]], RW, wsd(RS[:, ix[n_]], RW))
                 for n_ in ('zb_km', 'rho_rock_kgm3', 'compensation_C2')]
                + [('cmr2', PF['cmr2'][flow_fin['cmr2']], PR['cmr2'], WR, wsd(PR['cmr2'], WR))]):
            gap = abs(float(np.median(fv)) - wquantile(rv, rw_, 0.5)) / sd_
            chk[f'median_{nm}'] = {'gap_over_ref_sd': gap, 'gate': G2_MEDIAN_SD,
                                   'flow_median': float(np.median(fv)),
                                   'ref_median': wquantile(rv, rw_, 0.5), 'ref_sd': sd_,
                                   'pass': bool(gap <= G2_MEDIAN_SD)}
        # KS: C/MR2 and each observable pushforward
        for k in ['cmr2'] + obs:
            fv = PF[k][flow_fin[k]]
            d = ks_weighted(fv, np.ones(len(fv)), PR[k], WR)
            chk[f'ks_{k}'] = {'ks': d, 'gate': G2_KS, 'pass': bool(d <= G2_KS)}
        # per-observable median gap <= 0.1 sigma_infl
        for j, k in enumerate(obs):
            fv = PF[k][flow_fin[k]]
            gap = abs(float(np.median(fv)) - wquantile(PR[k], WR, 0.5)) / sig[j]
            chk[f'obs_median_{k}'] = {'gap_over_sigma_infl': gap, 'gate': G2_OBS_SIGMA,
                                      'pass': bool(gap <= G2_OBS_SIGMA)}
        g2['checks'] = chk
        g2['pass'] = bool(all(v['pass'] for v in chk.values()) and joint_pass
                          and g2['C2_gt_0.95']['pass'])
        rep['G2'] = g2

    # ===== helper for G3/G5: per-x flow draws + log q =====
    def lp(theta_np, x_vec):
        with torch.no_grad():
            t = torch.as_tensor(np.asarray(theta_np), dtype=torch.float32)
            xv = torch.as_tensor(x_vec, dtype=torch.float32)
            return srun._posterior.log_prob(t, x=xv, norm_posterior=False).cpu().numpy()

    def xdict(xv):
        return dict(zip(obs, [float(v) for v in xv]))

    # =========================== G3 =====================================
    with Timer('G3'):
        hs = sorted(SCR.glob('heldout_shard_*.npz'))
        assert hs, 'no held-out shards'
        TH0 = np.concatenate([np.load(f)['theta'] for f in hs])[:S['g3_n']]
        X0 = np.concatenate([np.load(f)['x'] for f in hs])[:S['g3_n']]
        ndr = S['g3_draws']
        frac = np.zeros((len(TH0), len(names)))
        cov68 = np.zeros((len(TH0), len(names)), bool)
        hpd = np.zeros(len(TH0))
        holefrac = np.zeros(len(TH0))
        for i in range(len(TH0)):
            dr = draw(srun, xdict(X0[i]), ndr, seed=900000 + i)
            frac[i] = np.mean(dr < TH0[i], axis=0)
            lo, hi = np.percentile(dr, [16, 84], axis=0)
            cov68[i] = (TH0[i] >= lo) & (TH0[i] <= hi)
            lq_d = lp(dr, X0[i])
            lq_0 = float(lp(TH0[i][None, :], X0[i])[0])
            hpd[i] = np.mean(lq_d < lq_0)
            if TH0[i][ix['zb_km']] < ZB_THIN_KM:
                holefrac[i] = hole_mask(runner, dr).mean()
            if (i + 1) % 200 == 0:
                log(f'G3 {i + 1}/{len(TH0)}')
        from scipy.stats import kstest

        def sbc_stats(sel):
            n = int(sel.sum())
            if n < 10:
                return {'n': n, 'note': 'too few pairs'}
            p = np.array([kstest(frac[sel, j], 'uniform').pvalue for j in range(len(names))])
            padj = bh_adjust(p)
            cov = cov68[sel].mean(axis=0)
            hp = kstest(hpd[sel], 'uniform')
            return {'n': n,
                    'per_param': {n_: {'ks_p': float(p[j]), 'ks_p_bh': float(padj[j]),
                                       'coverage68': float(cov[j]),
                                       'ks_pass': bool(padj[j] >= G3_PVAL),
                                       'cov_pass': bool(G3_COV_LO <= cov[j] <= G3_COV_HI)}
                                  for j, n_ in enumerate(names)},
                    'hpd_ks_stat': float(hp.statistic), 'hpd_ks_p': float(hp.pvalue),
                    'hpd_pass': bool(hp.pvalue >= G3_PVAL)}
        allsel = np.ones(len(TH0), bool)
        g3 = sbc_stats(allsel)
        g3['pass'] = bool(all(v['ks_pass'] and v['cov_pass'] for v in g3['per_param'].values())
                          and g3['hpd_pass'])
        g3['n_draws_per_pair'] = ndr
        g3['heldout_files'] = [str(f) for f in hs]
        thin_sel = TH0[:, ix['zb_km']] < ZB_THIN_KM
        g3['stratum_zb0_lt_15_nogate'] = sbc_stats(thin_sel)
        rep['G3'] = g3
        g3_holefrac = holefrac[thin_sel]

    # =========================== G4 =====================================
    with Timer('G4'):
        g4 = {'datasets': {}}
        ok_all = True
        for nm in INJ:
            cfgj = json.load(open(D / f'posterior_{nm}' / 'rep1_config.json'))
            xv = {k: float(cfgj['observables'][k][0]) for k in obs}
            dirs = [f'posterior_{nm}'] + (['posterior_inj_noise3s2027_seed2027'] if nm == 'inj_noise3' else [])
            PS, PW, pn = pool(dirs)
            PS = reorder(PS, pn, names)
            fl = draw(srun, xv, S['g4_n'], seed=4000 + INJ.index(nm))
            res = {'x': xv, 'ref_dirs': dirs, 'ref_ess': float(1 / np.sum(PW ** 2)), 'params': {}}
            for p_ in G4_PARAMS:
                j = ix[p_]
                med_r = wquantile(PS[:, j], PW, 0.5)
                sd_r = wsd(PS[:, j], PW)
                med_f, sd_f = float(np.median(fl[:, j])), float(fl[:, j].std())
                gap, ratio = abs(med_f - med_r) / sd_r, sd_f / sd_r
                ok = bool(gap <= G4_MED_SD and G4_SD_RATIO[0] <= ratio <= G4_SD_RATIO[1])
                ok_all &= ok
                res['params'][p_] = {'flow_median': med_f, 'ref_median': med_r, 'flow_sd': sd_f,
                                     'ref_sd': sd_r, 'gap_over_ref_sd': gap, 'sd_ratio': ratio,
                                     'gap_gate': G4_MED_SD, 'sd_ratio_bounds': list(G4_SD_RATIO),
                                     'pass': ok}
            res['pass'] = all(v['pass'] for v in res['params'].values())
            g4['datasets'][nm] = res
        g4['pass'] = bool(ok_all)
        rep['G4'] = g4

    # =========================== G5 =====================================
    with Timer('G5'):
        hf = []
        thin = np.load(SCR / 'thin_shard_00.npz')
        THt, Xt = thin['theta'], thin['x']
        if S['thin_n']:
            THt, Xt = THt[:S['thin_n']], Xt[:S['thin_n']]
        thin_hf = np.zeros(len(THt))
        for i in range(len(THt)):
            dr = draw(srun, xdict(Xt[i]), ndr, seed=950000 + i)
            thin_hf[i] = hole_mask(runner, dr).mean()
        pooled = np.concatenate([g3_holefrac, thin_hf])
        stat = lambda v: ({'n_pairs': int(len(v)), 'mean': float(v.mean()),
                           'p95': float(np.percentile(v, 95))} if len(v) else {'n_pairs': 0})
        g5 = {'g3_zb0_lt_15': stat(g3_holefrac), 'thin_shell': stat(thin_hf), 'pooled': stat(pooled),
              'gates': {'mean': G5_MEAN, 'p95': G5_P95}, 'n_draws_per_pair': ndr}
        g5['pass'] = bool(len(pooled) > 0 and pooled.mean() <= G5_MEAN
                          and np.percentile(pooled, 95) <= G5_P95)
        # x_obs: report only
        fo = draw(srun, xobs, S['guard_n'], seed=G1_SEED + 5)
        h_o = hole_mask(runner, fo)
        thin_o = fo[:, ix['zb_km']] < ZB_THIN_KM
        g5['at_x_obs_report_only'] = {'hole_fraction': float(h_o.mean()), 'n': int(len(fo))}
        rep['G5'] = g5
        rep['query_guard_at_x_obs'] = {
            'n': int(len(fo)), 'frac_zb_lt_15km': float(thin_o.mean()),
            'frac_hole': float(h_o.mean()), 'frac_union': float((thin_o | h_o).mean()),
            'threshold': GUARD_FRAC, 'label_required': bool((thin_o | h_o).mean() > GUARD_FRAC),
            'label': 'effective-prior posterior; servable cells only -- the cache hole is excluded '
                     'and the physical posterior there is unavailable to either SBI or pocoMC'}

    gates = {g: bool(rep[g]['pass']) for g in ('G1', 'G2', 'G3', 'G4', 'G5')}
    rep['gate_pass'] = gates
    rep['timings_s'] = TIMES
    rep['verdict'] = ('PASS: deploy as split-status-free "SBI, validated against pocoMC at x_obs and '
                      'globally calibrated (SBC)"' if all(gates.values()) else
                      'FAIL (' + ','.join(g for g, v in gates.items() if not v) +
                      '): record, do not deploy as validated, escalate; pocoMC run (a) stays authoritative')
    if a.smoke:
        rep['verdict'] = '[SMOKE - MEANINGLESS] ' + rep['verdict']
    outp.parent.mkdir(parents=True, exist_ok=True)
    outp.write_text(json.dumps(rep, indent=1, default=lambda o: o.tolist() if hasattr(o, 'tolist') else str(o)))
    log('wrote', outp, rep['verdict'])


if __name__ == '__main__':
    main()
