#!/usr/bin/env python
"""R1/R2 driver for the NH3 tidal-sector quarantine lift.

Binding spec: validation_reports/titan_freegrav_nh3_1m/is_correction/
r1_r2_PREREGISTRATION.md (PREREGISTERED, commit 634ea0cd) and
r1_r2_prereg_bands.json. No threshold, seed or rule here is the driver's own.

  --dry-run : ALL run-time provenance asserts + the two 200-sample plumbing
              recomputes (each side compared ONLY with ITSELF). Computes NO
              corrected-vs-reference statistic (no crosscheck, no m_corr, no
              KS, no median shift).
  (default) : dry-run checks, then R1a, R1b, R2, outcome; writes
              r1_r2_report.json. Refuses to overwrite an existing report.

Run with HOME=/tmp/pphome (TidalPy config) from the repo root.
"""
import argparse
import hashlib
import json
import pickle
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'plans/scripts'))

CONFIG = 'PlanetProfile/Inference/configs/test54_titan_nh3_freegrav.json'
ISC = REPO / 'validation_reports/titan_freegrav_nh3_1m/is_correction'
BANDS = ISC / 'r1_r2_prereg_bands.json'
LOGW = ISC / 'log_weights.npy'
RIDX = ISC / 'resampled_indices.npy'
REPORT = ISC / 'r1_r2_report.json'
BUNDLE_NPZ = REPO / 'PlanetProfile/Inference/sbi_artifacts/titan_freegrav_nh3_is_fiducial.npz'
BUNDLE_JSON = BUNDLE_NPZ.with_suffix('.json')
REFDIR = REPO / 'validation_reports/nh3_diagnosis/matched_reference_neff8000'
REF_REPORT = REFDIR / 'matched_reference_report.json'
SEEDS = (72, 172, 272)

# ---- preregistered constants (spec section in comments) -------------------
L_PREREG = 1270                 # Inputs
RESAMPLE_SEED = 72              # Inputs / Seeds
CROSSCHECK_SEED = 72            # Seeds
SE_BOOT_B, SE_BOOT_SEED = 1000, 72   # Seeds (C9 SE_corr)
PLUMB_N, PLUMB_SEED, PLUMB_RTOL = 200, 72, 1e-9   # R2 plumbing
CONFIG_HASH_PREFIX = '1611b65f'
NPZ_SHA_PREFIX = '0818b4e5'
CACHE_SHA_PREFIX = '3d837cf8'
SIGMA_C20, SIGMA_C22 = 4.71e-7, 6.2e-8      # R2 gate
SIGMA_RE, SIGMA_IM = 0.048, 0.035           # R1b
C9_NUISANCE = ['log10_eta_Ih', 'log10_eta_III', 'log10_eta_V',
               'log10_eta_VI', 'log10_eta_sil', 'log10_zeta', 'alpha']
NON_C9 = ['Tb_K', 'log10_wOcean_ppt', 'R_core_km', 'rho_core_kgm3',
          'dC20_nh', 'dC22_nh']
KS_MAX = 0.15                   # R1b
MEDIAN_GAP_MAX = 0.5            # R1b, G2, in sigma_obs
GRAV_MAX = 0.1                  # R2, in sigma_obs


# ===========================================================================
# Pure helpers (unit-tested in tests/nh3_r1r2_gates_test.py)
# ===========================================================================

def c9_decision(stock_pass, sigma_pass, m_corr, m_ref, tol):
    """C9 rule for an unidentified nuisance. Returns 'PASS', 'PASS_BY_C9'
    or 'FAIL'. Stock PASS => 'PASS'. Stock FAIL => PASS_BY_C9 iff sigma_pass
    AND |m_corr - m_ref| <= tol; otherwise FAIL (prereg R1a)."""
    if stock_pass:
        return 'PASS'
    if sigma_pass and abs(m_corr - m_ref) <= tol:
        return 'PASS_BY_C9'
    return 'FAIL'


def outcome_mapping(param_status, r1b_pass, r2_pass, c9_names=C9_NUISANCE):
    """Outcome rule + lift mapping.

    param_status: {name: 'PASS'|'PASS_BY_C9'|'FAIL'} for all 13 params.
    Lift takes effect iff all 6 non-C9 params PASS, R1b PASS, R2 PASS.
    C9 FAILs stay quarantined by name (lift mapping), and do not by
    themselves block the lift. 'all_13_pass' is reported separately so the
    manager can see the stricter reading of the outcome-rule table."""
    non_c9 = [n for n in param_status if n not in c9_names]
    non_c9_fail = sorted(n for n in non_c9 if param_status[n] == 'FAIL')
    quarantined = sorted(n for n in c9_names
                         if param_status.get(n) == 'FAIL')
    lift = bool(not non_c9_fail and r1b_pass and r2_pass)
    return {
        'verdict': 'LIFT_TAKES_EFFECT' if lift else 'NO_LIFT_REESCALATE',
        'non_c9_fails': non_c9_fail,
        'quarantined_by_name': quarantined,
        'c9_pass_by_c9': sorted(n for n in c9_names
                                if param_status.get(n) == 'PASS_BY_C9'),
        'r1a_structure_sector_pass': not non_c9_fail,
        'r1b_pass': bool(r1b_pass), 'r2_pass': bool(r2_pass),
        'all_13_pass': bool(all(v != 'FAIL' for v in param_status.values())),
    }


def max_rel_diff(new, old):
    new, old = np.asarray(new, float), np.asarray(old, float)
    bad = ~np.isfinite(new) | ~np.isfinite(old)
    if bad.any():
        return float('inf')
    return float(np.max(np.abs(new - old) / np.maximum(np.abs(old), 1e-300)))


def weights_from_log_w(log_w):
    """Same construction as is_correction.compute_is_correction (max-shifted
    exp, -inf -> 0, normalised)."""
    lw = np.asarray(log_w, float)
    fin = np.isfinite(lw)
    w = np.zeros_like(lw)
    w[fin] = np.exp(lw[fin] - lw[fin].max())
    return w / w.sum()


def pool_treatment_A(per_seed_weights):
    """Treatment A (r3_pool_and_gate.py): each seed w/sum(w)/n_seeds, then the
    pooled weights renormalised to sum 1."""
    n = len(per_seed_weights)
    W = np.concatenate([np.asarray(w, float) / np.sum(w) / n
                        for w in per_seed_weights])
    return W / W.sum()


# ===========================================================================
# Provenance + plumbing
# ===========================================================================

def sha256_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()


def git_sha():
    try:
        return subprocess.check_output(
            ['git', '-C', str(REPO), 'rev-parse', 'HEAD'], text=True).strip()
    except Exception:
        return None


class ProvenanceError(AssertionError):
    pass


def need(cond, msg):
    if not cond:
        raise ProvenanceError(msg)
    print(f'  [ok] {msg}')


def load_inputs():
    """Load bundle, reference seeds, bands, config; run every provenance
    assert that needs no runner. Returns a dict of loaded objects."""
    from PlanetProfile.Inference.inference_core import InferenceConfig
    from PlanetProfile.Inference.is_correction import systematic_resample
    from PlanetProfile.Inference.validate_sbi import _weighted_median_1d

    print('== provenance asserts ==')
    man = json.load(open(BUNDLE_JSON))
    bands = json.load(open(BANDS))
    npz_sha = sha256_file(BUNDLE_NPZ)
    need(npz_sha == man['npz_sha256'],
         f'bundle npz sha256 == manifest npz_sha256 ({npz_sha[:12]})')
    need(npz_sha.startswith(NPZ_SHA_PREFIX),
         f'bundle sha prefix == {NPZ_SHA_PREFIX} (prereg)')
    bz = np.load(BUNDLE_NPZ, allow_pickle=False)
    B = {k: bz[k] for k in bz.files}
    N = len(B['log_w'])
    need(N == 20000 == man['n_samples'], 'bundle N == 20000 == manifest')
    need(man['effective_seed'] == 72, 'bundle effective_seed == 72')

    lw_c = np.load(LOGW)
    lw = B['log_w']
    need(lw.shape == lw_c.shape, f'log_w shape {lw.shape} == committed')
    need(np.array_equal(np.isneginf(lw), np.isneginf(lw_c)),
         'log_w -inf mask identical to committed log_weights.npy')
    need(np.array_equal(np.isfinite(lw), np.isfinite(lw_c)),
         'log_w finite mask identical (no NaN/+inf differences)')
    fin = np.isfinite(lw_c)
    need(np.array_equal(lw[fin], lw_c[fin]),
         f'log_w finite entries bit-exact (n_finite={int(fin.sum())}, '
         f'n_-inf={int((~fin).sum())})')

    w = weights_from_log_w(lw)
    idx_c = np.load(RIDX)
    idx = systematic_resample(w, L_PREREG, seed=RESAMPLE_SEED)
    need(len(idx_c) == L_PREREG, f'committed resampled_indices L == {L_PREREG}')
    need(np.array_equal(idx, idx_c),
         'systematic_resample(w, 1270, seed=72) == committed '
         'resampled_indices.npy exactly')
    need(np.all(np.isfinite(lw[idx_c])),
         'no resampled index lands on a -inf (rejected) draw')

    cfg = InferenceConfig.from_json(str(REPO / CONFIG))
    # reference seeds
    ref = []
    hashes = []
    for sd in SEEDS:
        o = pickle.load(open(REFDIR / f'nh3_reference_seed{sd}_neff8000.pkl',
                             'rb'))
        hashes.append(o.config.generate_hash())
        ref.append(o)
    need(len(set(hashes)) == 1 and hashes[0].startswith(CONFIG_HASH_PREFIX),
         f'config_hash identical across 3 seeds == {hashes[0]}')
    need(cfg.generate_hash() == hashes[0],
         'current config file hash == reference config hash')
    names = list(ref[0].param_names)
    need(all(list(o.param_names) == names for o in ref),
         'param_names identical across 3 seeds')
    need(names == [str(x) for x in B['param_names']] == list(man['param_names']),
         'reference param_names == bundle param_names (same order)')
    for o in ref:
        need(all(np.size(getattr(o, a, None)) for a in
                 ('c20_results', 'c22_results', 'k2_results')),
             'reference seed carries c20/c22/k2 arrays')

    rr = json.load(open(REF_REPORT))
    cache_ref = rr['provenance']['structure_cache_sha256']
    need(cache_ref.startswith(CACHE_SHA_PREFIX) and
         cache_ref == man['structure_cache_sha256'],
         f'matched_reference_report cache sha == bundle manifest cache sha '
         f'({cache_ref[:12]})')
    cpath = REPO / cfg.structure_cache_path
    cache_disk = sha256_file(cpath)
    need(cache_disk == cache_ref,
         'structure-cache file on disk sha == both recorded shas')

    obs = {k: float(v[1]) for k, v in cfg.observables.items()}
    need(obs['C20'] == SIGMA_C20 and obs['C22'] == SIGMA_C22,
         f"config sigmas C20={obs['C20']} C22={obs['C22']} match prereg")
    need(obs['Re_k2'] == SIGMA_RE and obs['Im_k2'] == SIGMA_IM,
         f"config sigmas Re_k2={obs['Re_k2']} Im_k2={obs['Im_k2']} match prereg")
    need(man['obs_sigma'] == obs, 'bundle manifest obs_sigma == config sigmas')

    # c20/c22 availability (R2: never skip)
    for k in ('c20_results', 'c22_results', 'k2_results'):
        need(k in B, f'bundle carries {k}')
    nonrej = np.isfinite(lw)
    for k in ('c20_results', 'c22_results'):
        need(np.isfinite(B[k][nonrej]).all(),
             f'bundle {k} finite on all {int(nonrej.sum())} non-rejected draws')
        need(np.all(~np.isfinite(lw[~np.isfinite(B[k])])),
             f'bundle {k} non-finite entries all coincide with non-finite log_w')
    need(np.isfinite(B['k2_results'][nonrej]).all(),
         'bundle k2_results finite on non-rejected draws')
    for o in ref:
        need(np.isfinite(o.c20_results).all() and np.isfinite(o.c22_results).all(),
             f'reference seed c20/c22 100% finite (n={len(o.c20_results)})')

    # pooled reference, treatment A
    S = np.concatenate([np.asarray(o.samples, float) for o in ref])
    W = pool_treatment_A([np.asarray(o.weights, float) for o in ref])
    need(abs(W.sum() - 1) < 1e-12, 'pooled treatment-A weights sum to 1')
    ess = float(1.0 / np.sum(W ** 2))
    need(abs(ess - bands['ess_ref_pooled']) / ess < 1e-9,
         f"pooled ESS {ess:.3f} == bands ess_ref_pooled")

    # reference-only band re-derivation (guards the bands JSON against drift)
    per_seed = [(np.asarray(o.samples, float),
                 np.asarray(o.weights, float) / np.sum(o.weights))
                for o in ref]
    rng = np.random.RandomState(SE_BOOT_SEED)
    cdf = np.cumsum(W)
    meds = np.empty((SE_BOOT_B, len(C9_NUISANCE)))
    cols = [names.index(k) for k in C9_NUISANCE]
    for b in range(SE_BOOT_B):
        ii = np.searchsorted(cdf, rng.random_sample(L_PREREG) * cdf[-1])
        meds[b] = [np.median(S[ii, c]) for c in cols]
    for i, k in enumerate(C9_NUISANCE):
        bk = bands['c9'][k]
        c = names.index(k)
        s_seed = float(np.std([_weighted_median_1d(s_[:, c], w_)
                               for s_, w_ in per_seed], ddof=1))
        m_ref = _weighted_median_1d(S[:, c], W)
        se = float(meds[:, i].std(ddof=1))
        tol = 2.0 * float(np.hypot(s_seed, se))
        for nm, a, b_ in (('s_seed', s_seed, bk['s_seed_ddof1']),
                          ('m_ref', m_ref, bk['m_ref_pooled']),
                          ('SE_corr', se, bk['se_corr_L1270']),
                          ('tol', tol, bk['tol'])):
            need(abs(a - b_) <= 1e-9 * max(abs(b_), 1.0),
                 f'bands {k}.{nm} reproduced from reference pkls ({a:.6g})')

    return dict(cfg=cfg, man=man, bands=bands, B=B, lw=lw, w=w, idx=idx_c,
                ref=ref, names=names, S=S, W=W, ess=ess, obs=obs,
                npz_sha=npz_sha, cache_sha=cache_ref, cache_path=str(cpath),
                config_hash=hashes[0])


def plumbing(I):
    """Blocking R2 plumbing: each side recomputed and compared with ITSELF."""
    from PlanetProfile.Inference.mcmc_runner import MCMCRunner
    from PlanetProfile.Inference.sbi_runner import SBIRunner
    from PlanetProfile.Inference.inference_core import InferenceConfig

    print('== R2 plumbing recompute (self vs self only) ==')
    out = {}
    # Reference side: plain MCMCRunner on the (mode='mcmc') config, as the
    # reference campaign used. Corrected side: the SBI path's runner
    # (SBIRunner._get_mcmc_runner()), as the bundle builder used.
    t0 = time.time()
    cfg_ref = InferenceConfig.from_json(str(REPO / CONFIG))
    run_ref = MCMCRunner(cfg_ref)
    print(f'  MCMCRunner built ({time.time()-t0:.1f}s)')
    cfg_sbi = InferenceConfig.from_json(str(REPO / CONFIG))
    cfg_sbi.mode = 'sbi'
    t0 = time.time()
    run_sbi = SBIRunner(cfg_sbi)._get_mcmc_runner()
    print(f'  SBIRunner._get_mcmc_runner() built ({time.time()-t0:.1f}s)')

    def recompute(runner, samples):
        c20 = np.full(len(samples), np.nan)
        c22 = np.full(len(samples), np.nan)
        for i, th in enumerate(samples):
            pair = runner._derive_gravity_pair(runner._expand_theta(th))
            if pair is not None:
                c20[i], c22[i] = pair
        return c20, c22

    # reference: 200 uniform without replacement from the pooled reference
    rs = np.random.RandomState(PLUMB_SEED)
    ref_c20 = np.concatenate([np.asarray(o.c20_results, float) for o in I['ref']])
    ref_c22 = np.concatenate([np.asarray(o.c22_results, float) for o in I['ref']])
    ir = np.sort(rs.choice(len(I['S']), PLUMB_N, replace=False))
    t0 = time.time()
    n20, n22 = recompute(run_ref, I['S'][ir])
    d20, d22 = max_rel_diff(n20, ref_c20[ir]), max_rel_diff(n22, ref_c22[ir])
    out['reference'] = dict(n=PLUMB_N, seed=PLUMB_SEED, max_rel_c20=d20,
                            max_rel_c22=d22, elapsed_s=time.time() - t0,
                            index_sha=hashlib.sha256(ir.tobytes()).hexdigest()[:12],
                            ok=bool(d20 <= PLUMB_RTOL and d22 <= PLUMB_RTOL))
    print(f"  reference: max rel diff C20={d20:.3e} C22={d22:.3e} "
          f"({out['reference']['elapsed_s']:.1f}s) ok={out['reference']['ok']}")

    # corrected: 200 uniform without replacement from the NON-rejected draws
    # (rejected draws carry NaN/sentinel stored values by construction)
    rs = np.random.RandomState(PLUMB_SEED)
    pool = np.flatnonzero(np.isfinite(I['lw']))
    ic = np.sort(rs.choice(pool, PLUMB_N, replace=False))
    t0 = time.time()
    n20, n22 = recompute(run_sbi, np.asarray(I['B']['samples'], float)[ic])
    d20 = max_rel_diff(n20, I['B']['c20_results'][ic])
    d22 = max_rel_diff(n22, I['B']['c22_results'][ic])
    out['corrected'] = dict(n=PLUMB_N, seed=PLUMB_SEED, max_rel_c20=d20,
                            max_rel_c22=d22, elapsed_s=time.time() - t0,
                            index_sha=hashlib.sha256(ic.tobytes()).hexdigest()[:12],
                            ok=bool(d20 <= PLUMB_RTOL and d22 <= PLUMB_RTOL))
    print(f"  corrected: max rel diff C20={d20:.3e} C22={d22:.3e} "
          f"({out['corrected']['elapsed_s']:.1f}s) ok={out['corrected']['ok']}")
    out['ok'] = bool(out['reference']['ok'] and out['corrected']['ok'])
    return out


# ===========================================================================
# Real mode
# ===========================================================================

def run_real(I, plumb):
    from PlanetProfile.Inference.validate_sbi import (
        _run_crosscheck, _weighted_median_1d)
    from PlanetProfile.Inference.is_correction import weighted_quantile
    from is_correction_validate import (   # original C10 / C11 implementation
        weighted_ks, MEDIAN_GAP_SIGMA_MAX, WKS_D_MAX, GRAV_NOREG_SIGMA)
    assert (MEDIAN_GAP_SIGMA_MAX, WKS_D_MAX, GRAV_NOREG_SIGMA) == \
        (MEDIAN_GAP_MAX, KS_MAX, GRAV_MAX), 'prereg thresholds != C10/C11 consts'

    names, B, S, W, bands = I['names'], I['B'], I['S'], I['W'], I['bands']
    w = I['w']
    corr = np.asarray(B['samples'], float)[I['idx']]       # unweighted, L=1270

    # ---------------- R1a ----------------
    cc = _run_crosscheck(names, S, W, corr, seed=CROSSCHECK_SEED,
                         empirical_floor=None, matched_n_floor=True)
    per, status = [], {}
    for p in cc['per_parameter']:
        n = p['param']
        stock = bool(p['param_pass'])
        row = dict(p, stock_pass=stock)
        if n in C9_NUISANCE:
            bk = bands['c9'][n]
            m_corr = float(np.median(corr[:, names.index(n)]))
            st = c9_decision(stock, p['sigma_pass'], m_corr,
                             bk['m_ref_pooled'], bk['tol'])
            row.update(c9=True, m_corr=m_corr, m_ref=bk['m_ref_pooled'],
                       c9_tol=bk['tol'],
                       c9_abs_median_diff=abs(m_corr - bk['m_ref_pooled']),
                       c9_s_seed=bk['s_seed_ddof1'],
                       c9_se_corr=bk['se_corr_L1270'])
        else:
            st = 'PASS' if stock else 'FAIL'
            row.update(c9=False)
        row['verdict'] = st
        status[n] = st
        per.append(row)
    r1a = dict(n_mcmc=cc['n_mcmc'], ess_mcmc=cc['ess_mcmc'],
               n_sbi=cc['n_sbi'], floor_ess=cc.get('floor_ess'),
               matched_n_floor=cc.get('matched_n_floor'),
               stock_all_pass=cc['all_pass'], per_parameter=per,
               verdict_by_param=status)

    # ---------------- R1b (C10) ----------------
    kref = np.concatenate([np.asarray(o.k2_results, float) for o in I['ref']])
    kc = np.asarray(B['k2_results'], float)
    re_c, im_c = kc[:, 0], np.abs(kc[:, 1])
    re_r, im_r = kref[:, 0], np.abs(kref[:, 1])

    def pp(v, wt):
        return float(weighted_quantile(v, wt, 0.5)[0])

    okc, okr = np.isfinite(im_c), np.isfinite(im_r)
    ks_im = weighted_ks(im_c[okc], w[okc], im_r[okr], W[okr])
    ks_re = weighted_ks(re_c[okc], w[okc], re_r[okr], W[okr])
    gap_im = abs(pp(im_c, w) - pp(im_r, W)) / SIGMA_IM
    gap_re = abs(pp(re_c, w) - pp(re_r, W)) / SIGMA_RE
    r1b = dict(wks_re=ks_re, wks_im=ks_im, wks_max=KS_MAX,
               gap_re_sigma=gap_re, gap_im_sigma=gap_im,
               median_gap_max=MEDIAN_GAP_MAX,
               median_corr=dict(re=pp(re_c, w), im=pp(im_c, w)),
               median_ref=dict(re=pp(re_r, W), im=pp(im_r, W)),
               pass_=bool(ks_re <= KS_MAX and ks_im <= KS_MAX and
                          gap_re <= MEDIAN_GAP_MAX and gap_im <= MEDIAN_GAP_MAX))

    # ---------------- R2 (C11) ----------------
    r2 = dict(plumbing_ok=plumb['ok'], gate_sigma=GRAV_MAX, params={})
    if plumb['ok']:
        fin = np.isfinite(I['lw'])
        for nm, key, sig in (('C20', 'c20_results', SIGMA_C20),
                             ('C22', 'c22_results', SIGMA_C22)):
            vr = np.concatenate([np.asarray(getattr(o, key), float)
                                 for o in I['ref']])
            vc = np.asarray(B[key], float)
            mc, mr = pp(vc[fin], w[fin]), pp(vr, W)
            shift = abs(mc - mr) / sig
            ksv = weighted_ks(vc[fin], w[fin], vr, W)
            r2['params'][nm] = dict(median_corr=mc, median_ref=mr,
                                    shift_sigma_obs=shift, sigma_obs=sig,
                                    pass_=bool(shift <= GRAV_MAX),
                                    ks_informational=ksv)
        r2['pass_'] = bool(all(v['pass_'] for v in r2['params'].values()))
    else:
        r2['pass_'] = False
        r2['note'] = 'plumbing failed => R2 FAIL (prereg); stats not computed'

    outcome = outcome_mapping(status, r1b['pass_'], r2['pass_'])
    return dict(r1a=r1a, r1b=r1b, r2=r2, outcome=outcome)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()
    t0 = time.time()
    if not args.dry_run and REPORT.exists():
        raise SystemExit(f'{REPORT} exists; refusing to overwrite')
    I = load_inputs()
    plumb = plumbing(I)
    if args.dry_run:
        print('== dry-run summary ==')
        print(f"  N={len(I['lw'])} finite log_w={int(np.isfinite(I['lw']).sum())} "
              f"L={len(I['idx'])}  ESS_ref_pooled={I['ess']:.1f}")
        print(f"  config_hash={I['config_hash']}  npz_sha={I['npz_sha'][:12]}  "
              f"cache_sha={I['cache_sha'][:12]}")
        print(f"  plumbing ref   max_rel C20/C22 = "
              f"{plumb['reference']['max_rel_c20']:.3e}/"
              f"{plumb['reference']['max_rel_c22']:.3e}")
        print(f"  plumbing corr  max_rel C20/C22 = "
              f"{plumb['corrected']['max_rel_c20']:.3e}/"
              f"{plumb['corrected']['max_rel_c22']:.3e}")
        print(f"  plumbing ok = {plumb['ok']}  (rtol {PLUMB_RTOL})")
        print(f'  total runtime {time.time()-t0:.1f}s. NO R1/R2 statistic computed.')
        return 0 if plumb['ok'] else 1
    res = run_real(I, plumb)
    report = dict(
        kind='r1_r2_report', git_sha=git_sha(), driver='plans/scripts/nh3_r1r2_gates.py',
        prereg='r1_r2_PREREGISTRATION.md (634ea0cd)',
        inputs=dict(
            bundle_npz_sha256=I['npz_sha'], bundle_json_sha256=sha256_file(BUNDLE_JSON),
            log_weights_sha256=sha256_file(LOGW), resampled_indices_sha256=sha256_file(RIDX),
            bands_sha256=sha256_file(BANDS), config_sha256=sha256_file(REPO / CONFIG),
            config_hash=I['config_hash'], structure_cache_sha256=I['cache_sha'],
            reference_pkl_sha256={str(sd): sha256_file(REFDIR / f'nh3_reference_seed{sd}_neff8000.pkl')
                                  for sd in SEEDS}),
        seeds=dict(crosscheck=CROSSCHECK_SEED, resample=RESAMPLE_SEED,
                   se_boot=SE_BOOT_SEED, plumbing=PLUMB_SEED),
        plumbing=plumb, runtime_s=time.time() - t0, **res)
    with open(REPORT, 'w') as f:
        json.dump(report, f, indent=1, default=float)
    print(json.dumps(res['outcome'], indent=1))
    print(f'report -> {REPORT}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
