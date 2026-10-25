#!/usr/bin/env python
"""Fisher identifiability of the Enceladus ocean-branch (6D + rho_s) inference.

Reviewer requirement: BEFORE any MCMC, compute the Fisher matrix at reference
points in (zb_km, log10_wOcean_ppt, rho_ice_kgm3, rho_rock_kgm3 [= rho_s],
compensation_C2, libration_sys_frac) with the FULL data covariance, and report
corr(rho_s, C2), sigma(C2) and sigma(zb) with rho_s free vs fixed, sigma(rho_s),
the derived sigma(C/MR^2), the full correlation matrix and eigenvalues. If
|corr(rho_s, C2)| > 0.98 the point is flagged "C2 not separately constrained".

METHOD
  * Forward model = the runner's OWN functions, nothing re-implemented:
      MCMCRunner._derive_gravity_isostatic(theta)   -> C20, C22, C30
      MCMCRunner._derive_libration_deg(theta)       -> libration [deg]
      runner.last_mass_closure['cmr2_derived']      -> derived C/MR^2
    with metadata.isostasy.stack_mass_closure taken from the config
    ('mass_rho_rock' for the ocean6D_rhorock config), the theta dict carrying
    zb_km, log10_wOcean_ppt, rho_ice_kgm3, rho_rock_kgm3, compensation_C2,
    libration_sys_frac exactly as log_likelihood's _expand_theta would.
  * Data covariance = the likelihood's: sigma_obs inflated in quadrature with
    metadata.sigma_model_add (mcmc_runner._inflate_sigma), rho(C20,C22) =
    metadata.observable_correlations['C20,C22'] applied to the INFLATED
    sigmas, C30 and libration independent.
  * Jacobian by central finite differences (one-sided where a bound/edge is
    hit, flagged), at two step scales (--steps, default 1,2); the reported J
    is the smaller-step one and the per-column relative difference of the
    covariance-whitened columns is reported as the stability check.
  * zb / w DERIVATIVE HANDLING. The runner does NOT use nearest-node lookup:
    _ocean_zbw_struct_for_theta interpolates BILINEARLY in (zb, log10 w) via
    grid_interp_2d.bilinear_weights + blend_bilinear_corners. A bilinear
    interpolant is piecewise linear, so its derivative is piecewise constant
    with jumps at nodes, and the closure's discrete porosity / R_sil
    candidate choice adds node-to-node jitter. A step far below a node
    spacing would return one cell's slope; so the zb and log10 w steps are
    1x and 2x the LARGEST local node spacing of the cache axis (the central
    difference then spans >= one node and averages the kink). If the stencil
    leaves the cache span, the difference is one-sided and flagged.
    Smooth parameters (rho_ice, rho_s, C2, libration_sys_frac) use fixed
    steps (rho_ice 2.5, rho_s 10, C2 0.05, sys_frac 0.001) x multiplier.
  * F = J^T S^-1 J. Variants reported:
      A  likelihood only (all priors flat)
      B  A + Gaussian prior on libration_sys_frac (the only Gaussian prior;
         the truncation at +/-3 sigma is ignored); uniform priors flat.
         HEADLINE.
      C  B + uniform-prior widths as Gaussian proxies, sigma_u = width/sqrt(12)
         (a regularisation proxy for how much a uniform box bounds a weakly
         identified direction; NOT a proper treatment of a truncated prior).
    Marginal sigmas = sqrt(diag(F^-1)); "rho_s fixed" = delete the rho_s
    row/column from F, then invert (conditional sigmas). Columns with ~zero
    sensitivity are reported "unconstrained" (sigma = inf) and dropped from
    the inversion instead of poisoning it.
  * sigma(C/MR^2) = sqrt(g^T F^-1 g), g = gradient of the derived C/MR^2
    (same finite-difference stencil), linearised.

Usage (production cache is on the Mac only):
  python plans/scripts/enceladus_fisher_identifiability.py \
      --cache ~/scratch/enceladus_zbw_production_cache.pkl --out-dir DIR
  python plans/scripts/enceladus_fisher_identifiability.py --self-test
  --points "name:zb=25,logw=mid,rho_ice=925,rho_s=2400,C2=0.75,sys=0;..."
  --steps 1,2     multipliers on the base steps (first = reported J)
"""
import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

CONFIG_DEFAULT = str(REPO / 'PlanetProfile/Inference/configs/'
                     'enceladus_cassini_isostasy_ocean6D_rhorock.json')
SMOKE_CACHE = str(REPO / 'PlanetProfile/Test/mcmc_results/Enceladus/'
                  'Cassini_smoke/enceladus_seawater_smoke_grid.pkl')

PARAMS = ['zb_km', 'log10_wOcean_ppt', 'rho_ice_kgm3', 'rho_rock_kgm3',
          'compensation_C2', 'libration_sys_frac']
SHORT = {'zb_km': 'zb', 'log10_wOcean_ppt': 'log10w', 'rho_ice_kgm3': 'rho_ice',
         'rho_rock_kgm3': 'rho_s', 'compensation_C2': 'C2',
         'libration_sys_frac': 'sys'}
OBS = ['C20', 'C22', 'C30', 'libration_deg']
BASE_STEP = {'rho_ice_kgm3': 2.5, 'rho_rock_kgm3': 10.0,
             'compensation_C2': 0.05, 'libration_sys_frac': 0.001}
# prior bounds (fallback if the config lacks them); stencils never cross these
PRIOR_BOUNDS = {'zb_km': (5.0, 45.0), 'log10_wOcean_ppt': (-1.0, 2.0),
                'rho_ice_kgm3': (915.0, 935.0),
                'rho_rock_kgm3': (2200.0, 2600.0),
                'compensation_C2': (0.0, 1.0),
                'libration_sys_frac': (-0.012, 0.012)}
CORR_FLAG = 0.98
NULL_TOL = 1e-8     # whitened-column norm below this -> no sensitivity


# ------------------------------------------------------------------ pure math
def fisher_from_jacobian(J, cov, prior_prec=None):
    """F = J^T cov^-1 J (+ diag(prior_prec)). J: (n_obs, n_par)."""
    J = np.asarray(J, float)
    F = J.T @ np.linalg.solve(np.asarray(cov, float), J)
    F = 0.5 * (F + F.T)
    if prior_prec is not None:
        F = F + np.diag(np.asarray(prior_prec, float))
    return F


RANK_TOL = 1e-9     # scaled-Fisher eigenvalue / max below this -> null direction


def invert_fisher(F, tol=NULL_TOL, rank_tol=RANK_TOL):
    """Rank-aware covariance of a (possibly singular) Fisher matrix.

    F is scaled to unit diagonal, eigendecomposed; eigenvalues below
    rank_tol*max are NULL directions (infinite variance). Parameters with
    ~zero information (diag(F) <= tol^2 max diag) are dropped as inactive.
    Returns dict: C (n,n; inf on diagonal of params with a null-direction
    component, nan for pairs involving them), active, null_vecs (columns,
    in parameter space, unit-diagonal-scaled), rank, inf_mask, Ci (scaled
    covariance factor pieces used for linear functionals)."""
    F = np.asarray(F, float)
    n = F.shape[0]
    d = np.diag(F)
    active = d > (tol ** 2) * max(d.max(), 1e-300)
    idx = np.where(active)[0]
    C = np.full((n, n), np.nan)
    C[np.arange(n), np.arange(n)] = np.inf
    out = {'C': C, 'active': active, 'idx': idx, 'n': n}
    if idx.size == 0:
        out.update(rank=0, nullity=0, inf_mask=np.ones(n, bool),
                   null_scaled=np.zeros((n, 0)),
                   null_vecs=np.zeros((n, 0)), D=np.ones(n), lam=np.array([]),
                   V=np.zeros((0, 0)))
        return out
    sub = F[np.ix_(idx, idx)]
    D = np.sqrt(np.diag(sub))
    Fs = sub / np.outer(D, D)
    lam, V = np.linalg.eigh(Fs)
    good = lam > rank_tol * lam.max()
    Cs = (V[:, good] / lam[good]) @ V[:, good].T
    Csub = Cs / np.outer(D, D)
    nullv = V[:, ~good]
    inf_loc = np.any(np.abs(nullv) > 1e-6, axis=1) if nullv.shape[1] else \
        np.zeros(idx.size, bool)
    inf_mask = np.ones(n, bool)       # inactive params are unconstrained too
    inf_mask[idx] = inf_loc
    for a, ia in enumerate(idx):
        for b, ib in enumerate(idx):
            if not (inf_loc[a] or inf_loc[b]):
                C[ia, ib] = Csub[a, b]
    full_null = np.zeros((n, nullv.shape[1]))
    full_null[idx] = nullv / D[:, None]
    null_scaled = np.zeros((n, nullv.shape[1]))
    null_scaled[idx] = nullv
    out.update(C=C, rank=int(good.sum()), nullity=int((~good).sum()),
               inf_mask=inf_mask, null_vecs=full_null, null_scaled=null_scaled, D=D, lam=lam, V=V,
               good=good, Cs=Cs, inf_loc=inf_loc)
    return out


def corr_from_cov(inv):
    """Posterior correlation. Pairs inside an infinite (null) direction get
    the LIMITING correlation +/-1 when the nullity is 1 (sign of the null
    vector components), 0 when only one member is infinite, nan otherwise."""
    C, n = inv['C'], inv['n']
    R = np.full((n, n), np.nan)
    sd = np.sqrt(np.abs(np.diag(C)))
    for a in range(n):
        for b in range(n):
            if a == b:
                R[a, b] = 1.0
            elif not inv['inf_mask'][a] and not inv['inf_mask'][b]:
                R[a, b] = C[a, b] / (sd[a] * sd[b])
            elif inv['inf_mask'][a] != inv['inf_mask'][b]:
                R[a, b] = 0.0
            elif inv['nullity'] == 1:
                # limiting corr is +/-1 mathematically; only meaningful when
                # both members carry real weight (>0.1 of the unit null
                # vector in unit-diagonal-scaled coordinates)
                v = inv['null_scaled'][:, 0]
                R[a, b] = (float(np.sign(v[a] * v[b]))
                           if min(abs(v[a]), abs(v[b])) > 0.1 else np.nan)
    return R


def functional_sigma(inv, g):
    """sqrt(g^T C g) for a linear functional g; inf if g has a null-space
    component (the functional is not constrained by the data+priors)."""
    g = np.asarray(g, float)
    if not np.all(np.isfinite(g)) or inv['idx'].size == 0:
        return float('nan')
    idx = inv['idx']
    gmax = max(np.abs(g).max(), 1e-300)
    if np.any(np.abs(g[~inv['active']]) > 1e-12 * gmax):
        return float('inf')     # sensitive to a parameter with no information
    Vn = inv['V'][:, ~inv['good']]
    gD = g[idx] / inv['D']
    if Vn.shape[1] and np.linalg.norm(Vn.T @ gD) > 1e-6 * max(np.linalg.norm(gD), 1e-300):
        return float('inf')
    return float(np.sqrt(max(gD @ inv['Cs'] @ gD, 0.0)))


def analyze_fisher(F, names, g_cmr2=None, i_rho_s=None, i_C2=None, i_zb=None,
                   corr_flag=CORR_FLAG):
    """Identifiability summary of one Fisher matrix F (names = param order)."""
    n = len(names)
    inv = invert_fisher(F)
    C = inv['C']
    sig = np.sqrt(np.abs(np.diag(C)))
    R = corr_from_cov(inv)
    out = {'rank': inv['rank'], 'nullity': inv['nullity'],
           'sigma': {nm: _f(sig[k]) for k, nm in enumerate(names)},
           'unconstrained': {nm: bool(inv['inf_mask'][k])
                             for k, nm in enumerate(names)},
           'null_directions': [
               {nm: float(inv['null_vecs'][k, j] * np.sqrt(F[k, k]))
                for k, nm in enumerate(names)}
               for j in range(inv['null_vecs'].shape[1])],
           # same directions in raw parameter units (km, dex, kg/m^3, ...),
           # unit-normalised; the scaled form above is what the eigvecs show
           'null_directions_raw': [
               {nm: float(inv['null_vecs'][k, j]
                          / np.linalg.norm(inv['null_vecs'][:, j]))
                for k, nm in enumerate(names)}
               for j in range(inv['null_vecs'].shape[1])],
           'corr_posterior': _mat(R)}
    d = np.sqrt(np.maximum(np.diag(F), 1e-300))
    Fs = F / np.outer(d, d)
    w, v = np.linalg.eigh(Fs)
    out['eigvals_scaled_fisher'] = [float(x) for x in w]
    out['weakest_eigvec_scaled'] = {nm: float(v[k, 0])
                                    for k, nm in enumerate(names)}
    out['eigvals_fisher_raw'] = [float(x) for x in np.linalg.eigvalsh(F)]
    out['cond_scaled'] = float(w[-1] / max(w[0], 1e-300))
    if i_rho_s is not None and i_C2 is not None:
        c = R[i_rho_s, i_C2]
        out['corr_rho_s_C2'] = _f(c)
        out['corr_rho_s_C2_is_limit_of_singular_F'] = bool(
            inv['inf_mask'][i_rho_s] and inv['inf_mask'][i_C2])
        lim = out['corr_rho_s_C2_is_limit_of_singular_F']
        # the flag is set from a FINITE correlation only; a singular F is
        # reported as such (rank/nullity/null directions), not as a corr.
        out['flag_C2_not_separately_constrained'] = bool(
            np.isfinite(c) and not lim and abs(c) > corr_flag)
        out['singular_rho_s_C2_unconstrained'] = bool(lim or (
            inv['inf_mask'][i_rho_s] and inv['inf_mask'][i_C2]))
    if i_rho_s is not None:
        keep = [k for k in range(n) if k != i_rho_s]
        invf = invert_fisher(F[np.ix_(keep, keep)])
        sf = np.sqrt(np.abs(np.diag(invf['C'])))
        out['rank_rho_s_fixed'] = invf['rank']
        out['sigma_rho_s_fixed'] = {names[k]: _f(sf[j])
                                    for j, k in enumerate(keep)}
    if g_cmr2 is not None:
        out['sigma_cmr2_derived'] = _f(functional_sigma(inv, g_cmr2))
        if i_rho_s is not None:
            gk = np.asarray(g_cmr2, float)[keep]
            out['sigma_cmr2_derived_rho_s_fixed'] = _f(
                functional_sigma(invf, gk))
    return out



def _f(x):
    x = float(x)
    return x if np.isfinite(x) else (None if np.isnan(x) else 'inf')


def _mat(M):
    return [[_f(x) for x in row] for row in np.asarray(M, float)]


def build_covariance(obs, sigma_model_add, rho_c20_c22):
    """4x4 covariance in order OBS. Mirrors log_likelihood: inflate each
    sigma in quadrature with sigma_model_add (mcmc_runner._inflate_sigma),
    then apply rho(C20,C22) to the INFLATED sigmas; C30/libration
    independent."""
    s = []
    for nm in OBS:
        so = float(obs[nm][1])
        ex = sigma_model_add.get(nm)
        s.append(float(np.sqrt(so ** 2 + float(ex) ** 2)) if ex else so)
    s = np.array(s)
    cov = np.diag(s ** 2)
    if rho_c20_c22 is not None:
        cov[0, 1] = cov[1, 0] = float(rho_c20_c22) * s[0] * s[1]
    return cov, s


def fd_column(f, x0, i, h, lo, hi):
    """Derivative of vector function f wrt x[i]: central, one-sided at a
    bound. Returns (d, kind)."""
    xp, xm = x0.copy(), x0.copy()
    if x0[i] - h < lo:
        xp[i] = x0[i] + h
        return (f(xp) - f(x0)) / h, 'forward'
    if x0[i] + h > hi:
        xm[i] = x0[i] - h
        return (f(x0) - f(xm)) / h, 'backward'
    xp[i], xm[i] = x0[i] + h, x0[i] - h
    return (f(xp) - f(xm)) / (2.0 * h), 'central'


def jacobian(f, x0, steps, bounds):
    """f: x -> vector (obs..., cmr2). steps/bounds per parameter."""
    cols, kinds = [], []
    for i in range(len(x0)):
        d, k = fd_column(f, x0, i, steps[i], *bounds[i])
        cols.append(d)
        kinds.append(k)
    return np.array(cols).T, kinds


def whiten(J, cov):
    L = np.linalg.cholesky(cov)
    return np.linalg.solve(L, J)


# ------------------------------------------------------------- runner access
def make_runner(cfg_path):
    from PlanetProfile.Inference.inference_core import InferenceConfig
    from PlanetProfile.Inference.mcmc_runner import MCMCRunner
    config = InferenceConfig.from_json(cfg_path)
    runner = object.__new__(MCMCRunner)
    runner.config = config
    return runner


def make_forward(runner):
    """x (6 params, PARAMS order) -> [C20, C22, C30, libration_deg, cmr2]."""
    def fwd(x):
        th = dict(zip(PARAMS, [float(v) for v in x]))
        runner.last_mass_closure = None
        g = runner._derive_gravity_isostatic(th)
        lib = runner._derive_libration_deg(th)
        if g is None or lib is None:
            raise RuntimeError('forward model rejected theta=%s' % th)
        d = dict(getattr(runner, 'last_mass_closure', None) or {})
        cm = d.get('cmr2_derived', np.nan)
        return np.array([g[(2, 0)], g[(2, 2)], g[(3, 0)], lib, cm], float)
    return fwd


def axis_info(cache):
    zb = np.asarray(cache['zb_km_grid'], float)
    w = np.asarray(cache['wOcean_ppt_grid'], float)
    lw = np.log10(w)
    return {'zb': zb, 'logw': lw,
            'zb_h': float(np.max(np.diff(zb))), 'logw_h': float(np.max(np.diff(lw))),
            'zb_span': (float(zb[0]), float(zb[-1])),
            'logw_span': (float(lw[0]), float(lw[-1]))}


DEFAULT_POINTS = (
    'rho_s2300:zb=25,logw=mid,rho_ice=925,rho_s=2300,C2=0.75,sys=0;'
    'best_fit_rho_s2400:zb=25,logw=mid,rho_ice=925,rho_s=2400,C2=0.75,sys=0;'
    'rho_s2500:zb=25,logw=mid,rho_ice=925,rho_s=2500,C2=0.75,sys=0;'
    'zb20:zb=20,logw=mid,rho_ice=925,rho_s=2400,C2=0.75,sys=0;'
    'zb30:zb=30,logw=mid,rho_ice=925,rho_s=2400,C2=0.75,sys=0')


def parse_points(s, ax):
    pts = []
    for item in [t for t in s.split(';') if t.strip()]:
        name, _, body = item.partition(':')
        kv = dict(p.split('=') for p in body.split(','))
        logw = kv['logw']
        logw = (0.5 * (ax['logw_span'][0] + ax['logw_span'][1])
                if logw == 'mid' else float(logw))
        pts.append({'name': name.strip(), 'x': np.array([
            float(kv['zb']), logw, float(kv['rho_ice']), float(kv['rho_s']),
            float(kv['C2']), float(kv['sys'])])})
    return pts


def prior_info(cfg):
    ps = cfg.param_space
    bounds, gauss_prec, unif_prec = {}, np.zeros(6), np.zeros(6)
    for k, nm in enumerate(PARAMS):
        p = ps.get(nm, {}) or {}
        b = p.get('bounds', PRIOR_BOUNDS[nm])
        bounds[nm] = (float(b[0]), float(b[1]))
        if 'gauss' in str(p.get('prior_type', '')) or \
                p.get('prior_type') in ('normal',):
            gauss_prec[k] = 1.0 / float(p['std']) ** 2
        else:
            unif_prec[k] = 12.0 / (bounds[nm][1] - bounds[nm][0]) ** 2
    return bounds, gauss_prec, unif_prec


def run_point(fwd, pt, cov, ax, bounds, gauss_prec, unif_prec, mults):
    x0 = pt['x']
    lo_hi = [bounds[nm] for nm in PARAMS]
    # zb / logw stencils must also stay inside the cache span (else None)
    lo_hi[0] = (max(lo_hi[0][0], ax['zb_span'][0]),
                min(lo_hi[0][1], ax['zb_span'][1]))
    lo_hi[1] = (max(lo_hi[1][0], ax['logw_span'][0]),
                min(lo_hi[1][1], ax['logw_span'][1]))
    base = [ax['zb_h'], ax['logw_h']] + [BASE_STEP[nm] for nm in PARAMS[2:]]
    Js, kinds, steps_used = [], None, []
    f0 = fwd(x0)
    for m in mults:
        st = [b * m for b in base]
        J, kinds = jacobian(fwd, x0, st, lo_hi)
        Js.append(J)
        steps_used.append(dict(zip(PARAMS, st)))
    J = Js[0]
    A1, A2 = whiten(J[:4], cov), whiten(Js[-1][:4], cov)
    # relative difference of whitened columns, small vs large step
    rel = []
    for k in range(6):
        n2 = np.linalg.norm(A2[:, k])
        rel.append(float(np.linalg.norm(A1[:, k] - A2[:, k]) / n2)
                   if n2 > NULL_TOL else None)
    g = J[4]
    Fl = fisher_from_jacobian(J[:4], cov)
    variants = {'A_likelihood_only': Fl,
                'B_gauss_prior_sys (headline)': Fl + np.diag(gauss_prec),
                'C_B_plus_uniform_width_proxy': Fl + np.diag(gauss_prec)
                + np.diag(unif_prec)}
    res = {'name': pt['name'], 'x': dict(zip(PARAMS, map(float, x0))),
           'pred': dict(zip(OBS + ['cmr2_derived'], map(float, f0))),
           'fd_kind': dict(zip(PARAMS, kinds)),
           'steps': steps_used,
           'step_stability_rel_diff_whitened_col': dict(zip(PARAMS, rel)),
           'jacobian_per_unit': {OBS[i] if i < 4 else 'cmr2_derived':
                                 dict(zip(PARAMS, map(float, J[i])))
                                 for i in range(5)},
           'jacobian_whitened_norm': dict(zip(PARAMS, map(float, np.linalg.norm(A1, axis=0)))),
           'variants': {}}
    for vn, F in variants.items():
        res['variants'][vn] = analyze_fisher(
            F, PARAMS, g_cmr2=g, i_rho_s=3, i_C2=4, i_zb=0)
        res['variants'][vn]['fisher'] = _mat(F)
    return res


# ------------------------------------------------------------------- report
def _fmt(x, p=4):
    if x is None:
        return 'n/a'
    if x == 'inf':
        return 'inf'
    return f'{x:.{p}g}'


def render_md(rep):
    L = ['# Enceladus ocean-branch Fisher identifiability', '',
         f"Config: `{rep['config']}`  ", f"Cache: `{rep['cache']}`  ",
         f"Closure: `{rep['closure']}`; covariance sigmas (inflated): "
         f"{rep['sigma_inflated']}, rho(C20,C22)={rep['rho_c20_c22']}  ",
         f"Steps multipliers: {rep['mults']}; zb/log10w base step = largest "
         f"node spacing ({rep['axis']['zb_h']:.4g} km, "
         f"{rep['axis']['logw_h']:.4g} dex).  ",
         f"Self-test: {rep['self_test']}", '',
         '**Linearised (Fisher/Cramer-Rao) widths at a point; not posterior '
         'widths. Headline variant = B (Gaussian prior on libration_sys_frac '
         'only).** Four observables cannot identify six parameters: the '
         'data Fisher has rank <= 4 (variant A) and <= 5 with the sys_frac '
         'prior (B), so A/B are SINGULAR with all six free; their sigmas are '
         'inf along the null directions listed and correlations on those '
         'directions are the limiting +/-1. Finite widths in the free case '
         'come only from the prior box (variant C proxy). The rho_s-fixed '
         'columns are well posed in B.', '',
         ]
    for vkey in ('B_gauss_prior_sys (headline)', 'C_B_plus_uniform_width_proxy'):
        L += [f'## Summary (variant {vkey})', '',
              '| point | rank free/fixed | corr(rho_s,C2) | flag | sig(C2) free '
              '| sig(C2) rho_s fixed | sig(zb) free | sig(zb) rho_s fixed '
              '| sig(rho_s) | sig(C/MR2) free | sig(C/MR2) rho_s fixed |',
              '|' + '---|' * 11]
        for p in rep['points']:
            v = p['variants'][vkey]
            L.append('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |' % (
                p['name'], '%d/6, %d/5' % (v['rank'], v['rank_rho_s_fixed']),
                _fmt(v.get('corr_rho_s_C2')) + (
                    ' (limit, singular F)' if v.get(
                        'corr_rho_s_C2_is_limit_of_singular_F') else ''),
                'C2 NOT SEPARATELY CONSTRAINED' if v.get(
                    'flag_C2_not_separately_constrained') else 'ok',
                _fmt(v['sigma']['compensation_C2']),
                _fmt(v['sigma_rho_s_fixed']['compensation_C2']),
                _fmt(v['sigma']['zb_km']), _fmt(v['sigma_rho_s_fixed']['zb_km']),
                _fmt(v['sigma']['rho_rock_kgm3']),
                _fmt(v['sigma_cmr2_derived']),
                _fmt(v['sigma_cmr2_derived_rho_s_fixed'])))
        L.append('')
    for p in rep['points']:
        L += ['', f"## Point `{p['name']}`", '',
              'x = ' + ', '.join(f'{SHORT[k]}={v:.6g}' for k, v in p['x'].items()),
              '', 'pred = ' + ', '.join(f'{k}={v:.6g}' for k, v in p['pred'].items()),
              '', 'FD stencil: ' + ', '.join(f'{SHORT[k]}:{v}' for k, v in p['fd_kind'].items()),
              '', 'Step stability (relative difference of covariance-whitened '
              'Jacobian columns, step x%s vs x%s): ' % (rep['mults'][0], rep['mults'][-1])
              + ', '.join(f"{SHORT[k]}={_fmt(v, 3)}" for k, v in
                          p['step_stability_rel_diff_whitened_col'].items()),
              '', '||whitened J col||: ' + ', '.join(
                  f"{SHORT[k]}={v:.3g}" for k, v in p['jacobian_whitened_norm'].items()),
              '']
        for vn, v in p['variants'].items():
            L += [f'### {vn}', '',
                  'sigma (marginal): ' + ', '.join(
                      f"{SHORT[k]}={_fmt(s)}" for k, s in v['sigma'].items()),
                  '', 'sigma (rho_s fixed): ' + ', '.join(
                      f"{SHORT[k]}={_fmt(s)}" for k, s in v['sigma_rho_s_fixed'].items()),
                  '', f"corr(rho_s,C2) = {_fmt(v.get('corr_rho_s_C2'))}; "
                  f"sigma(C/MR2) = {_fmt(v['sigma_cmr2_derived'])} "
                  f"(rho_s fixed {_fmt(v['sigma_cmr2_derived_rho_s_fixed'])})",
                  '', f"rank {v['rank']}/6, nullity {v['nullity']}; "
                  'null directions (unit-diagonal-scaled Fisher '
                  'coordinates, i.e. raw component x sqrt(F_kk)): ' + (json.dumps(
                      [{SHORT[k]: round(x, 3) for k, x in d.items()}
                       for d in v['null_directions']]) or '[]'),
                  '', 'null directions (raw parameter units, unit-normalised): '
                  + (json.dumps(
                      [{SHORT[k]: float(f'{x:.3g}') for k, x in d.items()}
                       for d in v.get('null_directions_raw', [])]) or '[]'),
                  '', 'Posterior correlation matrix (order: '
                  + ', '.join(SHORT[k] for k in PARAMS) + '):', '', '```']
            for row in v['corr_posterior']:
                L.append(' '.join('%8s' % ('n/a' if x is None else
                                           ('inf' if x == 'inf' else '%.3f' % x))
                                  for x in row))
            L += ['```', '', 'Eigenvalues of unit-diagonal-scaled Fisher '
                  '(ascending): ' + ', '.join(f'{e:.4g}' for e in v['eigvals_scaled_fisher'])
                  + f"; cond = {v['cond_scaled']:.3g}", '',
                  'Weakest eigenvector: ' + ', '.join(
                      f"{SHORT[k]}={x:+.3f}" for k, x in v['weakest_eigvec_scaled'].items()),
                  '']
    return '\n'.join(L) + '\n'


# --------------------------------------------------------------- self-test
def synthetic_zbw_cache(zb_grid, w_grid):
    """Mass-conserved 3-layer Enceladus-like nodes in a v3.1-zbw container
    (same recipe as tests/enceladus_ocean_zbw_lookup_test._b2prime_struct)."""
    R_T, M, OM, ECC = 252.22e3, 1.08e20, 5.307e-5, 0.0047
    V = 4.0 * np.pi / 3.0
    structs = []
    for zb in zb_grid:
        for w in w_grid:
            R_b = R_T - zb * 1e3
            R_c = R_b - 38e3 * (1.0 + 0.02 * np.log10(w))
            rho_o = 1000.0 + 0.8 * w
            m_s = 925.0 * V * (R_T ** 3 - R_b ** 3)
            m_o = rho_o * V * (R_b ** 3 - R_c ** 3)
            rho_c = (M - m_s - m_o) / (V * R_c ** 3)
            structs.append({
                'r_m': np.array([R_c, R_b, R_T]),
                'rho': np.array([rho_c, rho_o, 925.0]),
                'phases': np.array([50, 0, 1]), 'omega': OM,
                'eccentricity': ECC, 'Mtot_kg': M,
                'changeIndices': np.array([0, 1, 2, 3]), 'n_layers': 3,
                'region_phases': ['Sil', 'water', 'Ih'],
                'layer_types': ['solid', 'liquid', 'solid'],
                'layer_upper_radii': [R_c, R_b, R_T], 'R_body_m': R_T})
    return {'zb_km_grid': np.asarray(zb_grid, float),
            'wOcean_ppt_grid': np.asarray(w_grid, float),
            'structures': structs, 'schema_version': 'v3.1-zbw',
            'ocean_comp': 'Seawater', 'frozen_branch_supported': False}


def smoke_exercise(runner, out):
    """The shipped smoke cache is 1-D (Tb) with no zb/w axes: pin its first
    ocean structure via a monkeypatched lookup so the rho_s / C2 / rho_ice /
    sys_frac block (zb, w columns are exactly zero) goes through the real
    gravity + libration + closure code."""
    with open(SMOKE_CACHE, 'rb') as fh:
        cache = pickle.load(fh)
    st = None
    for s in cache['structures']:
        if s is not None and np.any(np.asarray(s['phases']) == 0):
            st = s
            break
    if st is None:
        out['smoke'] = 'no ocean structure in smoke cache'
        return
    runner.structure_data = st
    runner._struct_for_hydrosphere = lambda theta, _s=st: _s
    fwd = make_forward(runner)
    x0 = np.array([25., 0.5, 925., 2400., 0.75, 0.0])
    try:
        f0 = fwd(x0)
    except Exception as e:  # smoke structure may not satisfy the closure
        out['smoke'] = 'forward failed on smoke structure: %r' % (e,)
        return
    out['smoke_pred'] = dict(zip(OBS + ['cmr2_derived'], map(float, f0)))
    out['smoke'] = 'forward OK (pinned structure; zb, log10w columns are zero)'
    cov, _ = build_covariance(_obs_from(runner), _sm_add(runner), _rho(runner))
    bounds, gp, up = prior_info(runner.config)
    ax = {'zb_span': (-1e9, 1e9), 'logw_span': (-1e9, 1e9),
          'zb_h': 1.0, 'logw_h': 0.1}
    res = run_point(fwd, {'name': 'smoke_pinned', 'x': x0}, cov, ax, bounds, gp,
                    up, [1.0, 2.0])
    out['smoke_point'] = res


def _obs_from(runner):
    return {k: tuple(v) for k, v in runner.config.observables.items()}


def _sm_add(runner):
    return (runner.config.metadata or {}).get('sigma_model_add', {}) or {}


def _rho(runner):
    return ((runner.config.metadata or {}).get('observable_correlations', {})
            or {}).get('C20,C22')


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--cache')
    ap.add_argument('--config', default=CONFIG_DEFAULT)
    ap.add_argument('--out-dir', default=None)
    ap.add_argument('--points', default=DEFAULT_POINTS)
    ap.add_argument('--steps', default='1,2',
                    help='step multipliers (first = reported Jacobian)')
    ap.add_argument('--self-test', action='store_true')
    a = ap.parse_args(argv)
    mults = [float(s) for s in a.steps.split(',')]
    runner = make_runner(a.config)
    st_note = 'no'
    extra = {}
    if a.self_test:
        zb_grid = np.arange(14.0, 36.0001, 1.0)
        w_grid = 10.0 ** np.linspace(-1.0, 2.0, 13)
        cache = synthetic_zbw_cache(zb_grid, w_grid)
        cache_name = 'SYNTHETIC zbw cache (%d zb x %d w nodes)' % (
            zb_grid.size, w_grid.size)
        st_note = ('YES -- synthetic mass-conserved 3-layer zbw cache through '
                   'the real bilinear lookup, closure, gravity and libration '
                   'code; numbers test code paths only, NOT science. Plus a '
                   'pinned-smoke-structure exercise (see JSON smoke_*).')
        smoke_runner = make_runner(a.config)
        smoke_exercise(smoke_runner, extra)
    else:
        if not a.cache:
            raise SystemExit('--cache required (or --self-test)')
        with open(Path(a.cache).expanduser(), 'rb') as fh:
            cache = pickle.load(fh)
        cache_name = str(a.cache)
    runner.structure_data = cache
    ax = axis_info(cache)
    cfg = runner.config
    cov, s = build_covariance(_obs_from(runner), _sm_add(runner), _rho(runner))
    bounds, gp, up = prior_info(cfg)
    fwd = make_forward(runner)
    pts = parse_points(a.points, ax)
    results = []
    for pt in pts:
        try:
            results.append(run_point(fwd, pt, cov, ax, bounds, gp, up, mults))
        except Exception as e:
            results.append({'name': pt['name'], 'error': repr(e)})
            print('POINT FAILED', pt['name'], repr(e), file=sys.stderr)
    good = [r for r in results if 'error' not in r]
    rep = {'config': str(a.config), 'cache': cache_name, 'self_test': st_note,
           'closure': (cfg.metadata.get('isostasy') or {}).get(
               'stack_mass_closure'),
           'sigma_inflated': dict(zip(OBS, map(float, s))),
           'rho_c20_c22': _rho(runner), 'mults': mults,
           'axis': {'zb_h': ax['zb_h'], 'logw_h': ax['logw_h'],
                    'zb_span': ax['zb_span'], 'logw_span': ax['logw_span']},
           'covariance': _mat(cov), 'obs_order': OBS, 'param_order': PARAMS,
           'priors': {'gaussian_precision': list(map(float, gp)),
                      'uniform_proxy_precision': list(map(float, up))},
           'points': good, 'failed_points': [r for r in results if 'error' in r],
           **extra}
    md = render_md(rep)
    if a.out_dir:
        od = Path(a.out_dir)
        od.mkdir(parents=True, exist_ok=True)
        (od / 'fisher_identifiability.json').write_text(
            json.dumps(rep, indent=1, default=str))
        (od / 'fisher_identifiability.md').write_text(md)
        print('wrote', od)
    print(md[:6000])
    return 0 if good else 1


if __name__ == '__main__':
    sys.exit(main())
