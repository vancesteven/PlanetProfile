#!/usr/bin/env python
"""Mass-closure jitter gate for the Enceladus (zb, w) cache.

Replaces the C3 acceptance condition (scientific-reviewer ruling
2026-10-21, validation_reports/enceladus_isostasy/c3_root_cause_2026-10-21/
README.md). The C3 extrapolate-to-zero-deficit study is invalid because C20
and C22 are a jagged function of (zb, w) through the discrete porosity /
R_sil candidate choice, not just the mass deficit.

PROCEDURE
  * 8-10 fixed (zb, w) nodes (default includes cache indices (70,18) and
    (77,10) plus a spread over the grid).
  * Rebuild each node over Ocean.deltaP/{1,2,4} (0.1, 0.05, 0.025 MPa)
    x Steps.nPoros {8, 12, 16} (9 builds/node) through the builder's own
    override mechanisms (ocean_overrides['deltaP'], planet_overrides
    ['Steps'] with nPoros) -- same per-node call as production
    (enceladus_zbw_mass_deficit_study._production_kwargs).
  * C20/C22/C30/libration through the campaign's own forward model
    (MCMCRunner._derive_gravity_isostatic / _derive_libration_deg) with
    metadata.isostasy.stack_mass_closure = 'rock_density' (ON, the gated
    quantity) and without it (OFF, reference). Worst case over
    compensation_C2 in {0, 0.5, 1}.

PASS (all required)
  1. C20 and C22: spread (max-min over the 9 builds) <= 0.5 sigma_obs at
     EVERY node. Fallback: if spread fails, sigma_model := max(std over
     builds (ddof=1), |bias|) (max over nodes), where bias = production
     (deltaP/1, nPoros 8) minus finest (deltaP/4, nPoros 16); PASS if
     sigma_model <= 0.5 sigma_obs. Any failed build -> INCONCLUSIVE (INTERPRETATION
     of "in quadrature": see open question in the report; the value that
     must then be added to metadata.sigma_model_add.C20/C22 is reported, and
     the likelihood combines it as hypot(sigma_obs_inflated, sigma_model)).
  2. libration spread <= 0.1 sigma_obs at every node.
  3. Local-quadratic check: the production-setting (cached) node value
     (closure applied to the cache's neighbouring structures) must not depart
     from a quadratic fit in (zb, log10 w) over its 5x5 neighbourhood
     (node itself excluded) by > 3x the node's jitter std. Needs --cache;
     otherwise reported SKIPPED and the verdict is INCONCLUSIVE.
  4. TODO (later step): posterior-median shift < 0.2 posterior sigma
     between sigma_model on and off. Not computed here.

MODE mass_moi (default; Option 4, reviewer ruling 2026-10-22, README in
validation_reports/enceladus_isostasy/mass_closure_gate/): re-run on the same
nodes/builds. PASS needs: C20/C22 spread <= 0.1 sigma_obs at EVERY node;
libration spread <= 0.05 sigma; closed-stack mass residual ~1e-15 (limit
1e-12); every node/build feasible; offset (mass_moi - rock_density, finest deltaP) within
0.1 sigma_obs of the predicted (moi_target - MoI_rd)*dX/dC at every node; local-quadratic check
(needs --cache). Reported per node: R_s - R_sil(PP), rho_s, stack MoI before/
after vs PP CMR2. Posterior-median check: TODO (not computed).
--closure-mode rock_density reproduces the original 2026-10-21 gate.

Usage:
  HOME=~/scratch/pphome OMP_NUM_THREADS=1 python \\
     plans/scripts/enceladus_mass_closure_jitter_gate.py \\
     --cache PlanetProfile/Test/mcmc_results/Enceladus/Cassini_isostasy/enceladus_zbw_production_cache.pkl \\
     --out-dir ~/scratch/jitter_gate [--nodes "70,18;77,10;..."] [--jobs 4]
  Cache scan (no rebuilds, reporting only):  --cache ... --scan-cache
  Reduced exercise, no cache:  --self-test   (1 node, 2 builds by default)
"""
import argparse
import copy
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from enceladus_zbw_production_build import (  # noqa: E402
    CONFIG_PATH, _load_config, _zb_axis, _w_axis)
from enceladus_zbw_mass_deficit_study import (  # noqa: E402
    _production_kwargs, make_runner, sigma_obs, OBS_NAMES, C2_GRID)

DEFAULT_NODES = [(70, 18), (77, 10), (12, 6), (25, 30), (38, 18), (50, 34),
                 (58, 12), (64, 26), (84, 4), (86, 22)]
DP_FACTORS = (1, 2, 4)          # deltaP 0.1, 0.05, 0.025
NPOROS = (8, 12, 16)
SPREAD_C = 0.5                  # rock_density mode: C20/C22 spread <= 0.5
SPREAD_LIB = 0.1                # rock_density mode: libration spread <= 0.1
# mass_moi (Option 4, reviewer ruling 2026-10-22) limits:
MM_SPREAD_C = 0.1               # C20/C22 spread <= 0.1 sigma_obs, EVERY node
MM_SPREAD_LIB = 0.05            # libration spread <= 0.05 sigma_obs
MM_RESIDUAL = 1e-12             # closed-stack |M/Mtot - 1| (expected ~1e-15)
MM_OFFSET = 0.1                 # |measured - predicted| offset vs rock_density
DC_STEP = 5e-4                  # moi_target +- step for dX/dC
LOCALFIT_K = 3.0                # departure <= 3 x jitter std
GATED = ('C20', 'C22')


def parse_nodes(s):
    return [tuple(int(x) for x in p.split(',')) for p in s.split(';') if p]


def build_task(task):
    """Top-level (picklable) single build. task = (i, j, zb, w, dpf, npor)."""
    import importlib
    from PlanetProfile.Inference.cache_builder import (
        build_single_structure, _ocean_node_zb_km)
    i, j, zb_km, w_ppt, dpf, npor = task
    cfg = _load_config()
    kw = _production_kwargs(cfg, zb_km, w_ppt)
    mod = importlib.import_module(kw['planet_template_module'])
    dP0 = float(mod.Planet.Ocean.deltaP)
    dP1 = round(dP0 / dpf, 3)
    if abs(dP0 / dpf - dP1) > 1e-12:
        raise ValueError(f'deltaP {dP0}/{dpf} needs >3 decimals')
    if dpf != 1:
        kw['ocean_overrides']['deltaP'] = dP1
    if int(npor) != int(mod.Planet.Steps.nPoros):
        steps = copy.deepcopy(mod.Planet.Steps)
        steps.nPoros = int(npor)
        kw['planet_overrides'] = {'Steps': steps}
    t0 = time.time()
    try:
        s = build_single_structure(**kw)
        s['zb_km_actual'] = _ocean_node_zb_km(s)
        err = None
    except Exception as exc:  # recorded, not fatal
        s, err = None, f'{type(exc).__name__}: {exc}'
    return {'task': task, 'deltaP': dP1, 'nPoros': int(npor), 'struct': s,
            'wall_s': time.time() - t0, 'error': err}


def make_runners(mode='mass_moi'):
    """(on, off, ref): on = gated closure `mode`; off = no closure;
    ref = 'rock_density' (the legacy closure; reference for the mass_moi
    offset check)."""
    on = make_runner(CONFIG_PATH)
    on.config.metadata['isostasy']['stack_mass_closure'] = mode
    ref = make_runner(CONFIG_PATH)
    ref.config.metadata['isostasy']['stack_mass_closure'] = 'rock_density'
    # capture the rock_density-closed reduced stack (for its MoI)
    ref._captured = {}
    _orig = ref._close_reduced_stack_mass

    def _wrap(r, rho, st, _orig=_orig, _ref=ref):
        out = _orig(r, rho, st)
        if out is not None:
            _ref._captured['v'] = (np.array(r, float), np.array(out, float))
        return out
    ref._close_reduced_stack_mass = _wrap
    # +-DC_STEP moi_target runners (derivative dX/dC). The consistency guard
    # is bypassed on purpose: these probe off-target closures.
    pm = {}
    for sgn in (+1, -1):
        rn = make_runner(CONFIG_PATH)
        rn.config.metadata['isostasy']['stack_mass_closure'] = 'mass_moi'
        rn.config.metadata['isostasy']['moi_target'] = float(
            rn.config.metadata['isostasy']['moi_target']) + sgn * DC_STEP
        rn.config.planet_template_module = None
        rn.config.metadata.pop('bulk_overrides', None)
        pm[sgn] = rn
    ref._pm = pm
    off = make_runner(CONFIG_PATH)
    off.config.metadata['isostasy'].pop('stack_mass_closure', None)
    return on, off, ref


def observables(runner, struct):
    """{obs: max-over-C2 handled by caller -> {C2: value}}; also the closure
    diagnostic of the gravity call."""
    runner.structure_data = struct
    out = {k: {} for k in OBS_NAMES}
    diag = None
    for c2 in C2_GRID:
        th = {'compensation_C2': c2}
        g = runner._derive_gravity_isostatic(th)
        if c2 == 0.5:
            diag = getattr(runner, 'last_mass_closure', None)
        lib = runner._derive_libration_deg(th)
        for k, lm in (('C20', (2, 0)), ('C22', (2, 2)), ('C30', (3, 0))):
            out[k][c2] = np.nan if g is None else float(g[lm])
        out['libration_deg'][c2] = np.nan if lib is None else float(lib)
    return out, diag


def spread_std(vals_by_build, sig):
    """vals_by_build: list of {C2: v}. Worst case over C2 of spread/std,
    in sigma_obs units."""
    sp, sd = [], []
    for c2 in C2_GRID:
        v = np.array([b[c2] for b in vals_by_build], float)
        if not np.all(np.isfinite(v)):
            return None, None
        sp.append((v.max() - v.min()) / sig)
        sd.append(v.std(ddof=1) / sig if v.size > 1 else 0.0)
    return float(max(sp)), float(max(sd))


def local_quadratic_check(cache, runner_on, i, j, std_sig_by_obs, sig,
                          half=2):
    """Departure of cached node (i,j) (closure ON) from a quadratic fit in
    (zb, log10 w) over its (2*half+1)^2 neighbourhood excluding itself.
    Returns {obs: worst-over-C2 |dep| / (3-sigma-limit)} style dict:
    {'dep_sigma', 'jitter_std_sigma', 'ratio'}."""
    zb = np.asarray(cache['zb_km_grid'], float)
    w = np.log10(np.asarray(cache['wOcean_ppt_grid'], float))
    nw = w.size
    pts = []
    for a in range(max(0, i - half), min(zb.size, i + half + 1)):
        for b in range(max(0, j - half), min(nw, j + half + 1)):
            s = cache['structures'][a * nw + b]
            if s is None:
                continue
            o, _ = observables(runner_on, s)
            pts.append((a, b, zb[a], w[b], o))
    center = [p for p in pts if (p[0], p[1]) == (i, j)]
    others = [p for p in pts if (p[0], p[1]) != (i, j)]
    if not center or len(others) < 8:
        return None
    x0, y0 = zb[i], w[j]
    sx = max(1e-9, np.ptp([p[2] for p in pts]))
    sy = max(1e-9, np.ptp([p[3] for p in pts]))

    def design(x, y):
        u, v = (x - x0) / sx, (y - y0) / sy
        return np.stack([np.ones_like(u), u, v, u * u, u * v, v * v], -1)
    X = design(np.array([p[2] for p in others]),
               np.array([p[3] for p in others]))
    res = {}
    for k in OBS_NAMES:
        worst, jit = 0.0, std_sig_by_obs.get(k)
        for c2 in C2_GRID:
            y = np.array([p[4][k][c2] for p in others])
            if not np.all(np.isfinite(y)):
                continue
            beta, *_ = np.linalg.lstsq(X, y, rcond=None)
            pred = beta[0]       # design at the node itself is [1,0,0,0,0,0]
            dep = abs(center[0][4][k][c2] - pred) / sig[k]
            worst = max(worst, dep)
        res[k] = {'dep_sigma_obs': float(worst),
                  'jitter_std_sigma_obs': jit,
                  'ratio_dep_over_jitter_std': (float(worst / jit)
                                                if jit else None),
                  'pass': bool(jit is not None and worst <= LOCALFIT_K * jit)}
    return res



def _pair(ons):
    return (np.array([o['C20'][0.5] for o in ons], float),
            np.array([o['C22'][0.5] for o in ons], float))


def _corr_pair(ons):
    if len(ons) < 3:
        return None
    a, b = _pair(ons)
    if a.std() == 0 or b.std() == 0:
        return None
    return float(np.corrcoef(a, b)[0, 1])


def _dev_pair(ons):
    if len(ons) < 3:
        return None
    a, b = _pair(ons)
    return (a - a.mean()).tolist(), (b - b.mean()).tolist()


def bias_prod_minus_finest(blds, runner_on, sig):
    """Systematic bias = production setting (deltaP/1, nPoros 8) minus the
    finest available build (largest deltaP divisor, then largest nPoros;
    full grid: deltaP/4, nPoros 16), closure ON, per observable, in
    sigma_obs, worst case over compensation_C2 (signed value of the
    largest |.|)."""
    ok = [b for b in blds if b['struct'] is not None]
    prod = [b for b in ok if b['task'][4] == 1 and b['nPoros'] == 8]
    if not ok or not prod:
        return None
    fine = max(ok, key=lambda b: (b['task'][4], b['nPoros']))
    if fine is prod[0]:
        return None
    op, _ = observables(runner_on, prod[0]['struct'])
    of, _ = observables(runner_on, fine['struct'])
    out = {'finest': {'deltaP': fine['deltaP'], 'nPoros': fine['nPoros']}}
    for k in OBS_NAMES:
        d = [(op[k][c] - of[k][c]) / sig[k] for c in C2_GRID]
        out[k] = (None if not np.all(np.isfinite(d))
                  else float(max(d, key=abs)))
    return out


def run(args):
    cfg = _load_config()
    sig, sig_infl = sigma_obs(cfg)
    zb_ax, w_ax = _zb_axis(), _w_axis()
    cache = None
    if args.cache:
        with open(Path(args.cache).expanduser(), 'rb') as f:
            cache = pickle.load(f)
        zb_ax = np.asarray(cache['zb_km_grid'], float)
        w_ax = np.asarray(cache['wOcean_ppt_grid'], float)
    if args.self_test:
        nodes = [(70, 18)]
        dps = tuple(int(x) for x in args.dp_factors.split(',')) \
            if args.dp_factors else (1, 2)
        nps = tuple(int(x) for x in args.nporos.split(',')) \
            if args.nporos else (8,)
    else:
        nodes = parse_nodes(args.nodes) if args.nodes else DEFAULT_NODES
        dps = tuple(int(x) for x in args.dp_factors.split(',')) \
            if args.dp_factors else DP_FACTORS
        nps = tuple(int(x) for x in args.nporos.split(',')) \
            if args.nporos else NPOROS
    tasks = [(i, j, float(zb_ax[i]), float(w_ax[j]), d, n)
             for (i, j) in nodes for d in dps for n in nps]
    print(f'{len(nodes)} nodes x {len(dps)} deltaP x {len(nps)} nPoros = '
          f'{len(tasks)} builds (jobs={args.jobs})', flush=True)
    t0 = time.time()
    if args.jobs > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(args.jobs) as ex:
            built = list(ex.map(build_task, tasks))
    else:
        built = []
        for t in tasks:
            b = build_task(t)
            print(f"  node ({t[0]},{t[1]}) dP/{t[4]} nPoros={t[5]}: "
                  f"{b['wall_s']:.1f}s {b['error'] or 'ok'}", flush=True)
            built.append(b)
    total_wall = time.time() - t0

    mode = args.closure_mode
    run_on, run_off, run_ref = make_runners(mode)
    per_node, nodes_out = {}, []
    for b in built:
        per_node.setdefault(b['task'][:2], []).append(b)
    for (i, j) in nodes:
        rec = {'i_zb': i, 'i_w': j, 'zb_km': float(zb_ax[i]),
               'wOcean_ppt': float(w_ax[j]), 'builds': []}
        ons, offs = [], []
        for b in per_node[(i, j)]:
            br = {'deltaP': b['deltaP'], 'nPoros': b['nPoros'],
                  'wall_s': round(b['wall_s'], 2), 'error': b['error']}
            if b['struct'] is not None:
                oo, diag = observables(run_on, b['struct'])
                of, _ = observables(run_off, b['struct'])
                ref_o, ref_moi, dxdc = None, None, None
                if mode == 'mass_moi':
                    ref_o, _ = observables(run_ref, b['struct'])
                    rr_, rrho_ = run_ref._captured['v']
                    ref_moi = run_ref._stack_cmr2(rr_, rrho_)
                    op_, _ = observables(run_ref._pm[+1], b['struct'])
                    om_, _ = observables(run_ref._pm[-1], b['struct'])
                    dxdc = {k: (op_[k][0.5] - om_[k][0.5]) / (2 * DC_STEP)
                            for k in OBS_NAMES}
                s = b['struct']
                br.update({
                    'on': {k: oo[k][0.5] for k in OBS_NAMES},
                    'off': {k: of[k][0.5] for k in OBS_NAMES},
                    'closure_diag': diag,
                    'ref_rock_density': (None if ref_o is None else
                                         {k: ref_o[k][0.5]
                                          for k in OBS_NAMES}),
                    'pp_CMR2': s.get('CMR2'),
                    'ref_rock_density_moi': ref_moi,
                    'dXdC': dxdc,
                    'mass_residual_frac': s.get('mass_residual_frac'),
                    'zb_km_actual': s.get('zb_km_actual')})
                ons.append(oo)
                offs.append(of)
            rec['builds'].append(br)
        rec['n_ok'] = len(ons)
        stats = {}
        for label, lst in (('on', ons), ('off', offs)):
            stats[label] = {}
            for k in OBS_NAMES:
                if len(lst) < 2:
                    stats[label][k] = {'spread_sigma_obs': None,
                                       'std_sigma_obs': None}
                    continue
                sp, sd = spread_std([o[k] for o in lst], sig[k])
                stats[label][k] = {'spread_sigma_obs': sp,
                                   'std_sigma_obs': sd}
        rec['stats'] = stats
        rec['failed_builds'] = [(b['deltaP'], b['nPoros'])
                                for b in rec['builds'] if b['error']
                                or 'on' not in b]
        rec['bias'] = bias_prod_minus_finest(per_node[(i, j)], run_on, sig)
        rec['c20_c22_corr'] = _corr_pair(ons)
        rec['_dev'] = _dev_pair(ons)
        rec['local_quadratic'] = None
        if cache is not None:
            std_s = {k: stats['on'][k]['std_sigma_obs'] for k in OBS_NAMES}
            try:
                rec['local_quadratic'] = local_quadratic_check(
                    cache, run_on, i, j, std_s, sig)
            except Exception as exc:
                rec['local_quadratic'] = {'error': f'{type(exc).__name__}: {exc}'}
        nodes_out.append(rec)
    verdict = adjudicate(nodes_out, sig, cache is not None, cfg, mode)
    frozen = frozen_check(cache) if (args.frozen_check and cache) else None
    for n in nodes_out:
        n.pop('_dev', None)
    meta = {'closure_mode': mode, 'cache': args.cache, 'self_test': args.self_test,
            'dp_factors': list(dps), 'nporos': list(nps),
            'n_builds': len(tasks), 'total_wall_s': round(total_wall, 1),
            'mean_build_wall_s': round(float(np.mean(
                [b['wall_s'] for b in built])), 2),
            'sigma_obs': sig, 'sigma_obs_inflated': sig_infl,
            'config': CONFIG_PATH}
    out = {'meta': meta, 'verdict': verdict, 'nodes': nodes_out,
           'frozen_check': frozen}
    out_dir = Path(args.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / 'jitter_gate.json').write_text(
        json.dumps(out, indent=1, default=float))
    (out_dir / 'jitter_gate.md').write_text(render_md(out))
    print(render_md(out))
    return out


def adjudicate(nodes_out, sig, have_cache, cfg=None, mode='rock_density'):
    """sigma_model per observable = max over nodes of max(std(ddof=1),
    |bias|) in raw sigma_obs units. PASS needs C20/C22 sigma_model <= 0.5
    sigma_obs and libration spread <= 0.1 sigma_obs at every node."""
    v = {'criteria': {}, 'sigma_model': {}}
    full = [n for n in nodes_out if n['n_ok'] >= 2]
    failed = {f"({n['i_zb']},{n['i_w']})": n['failed_builds']
              for n in nodes_out if n['failed_builds']}
    v['failed_builds_by_node'] = failed
    add = ((cfg or {}).get('metadata', {}).get('sigma_model_add', {}) or {})
    for k in OBS_NAMES:
        sds = [n['stats']['on'][k]['std_sigma_obs'] for n in full]
        bs = [abs(n['bias'][k]) for n in full
              if n.get('bias') and n['bias'].get(k) is not None]
        std_m = max([x for x in sds if x is not None], default=None)
        bias_m = max(bs, default=None)
        cand = [x for x in (std_m, bias_m) if x is not None]
        sm = max(cand) if cand else None
        v['sigma_model'][k] = {
            'max_std_over_sigma_obs': std_m,
            'max_abs_bias_over_sigma_obs': bias_m,
            'sigma_model_over_sigma_obs': sm,
            'sigma_model_abs': None if sm is None else float(sm * sig[k]),
            'recommended_sigma_model_add_abs': (
                None if sm is None or k not in add else
                float(np.hypot(add[k], sm * sig[k]))),
            'existing_sigma_model_add': add.get(k)}
    ok1 = True
    mm = mode == 'mass_moi'
    lim_c, lim_l = (MM_SPREAD_C, MM_SPREAD_LIB) if mm else (SPREAD_C,
                                                            SPREAD_LIB)
    if mm:
        for k in GATED:
            sp = [n['stats']['on'][k]['spread_sigma_obs'] for n in full]
            good = bool(sp) and len(sp) == len(nodes_out) and all(
                x is not None and x <= lim_c for x in sp)
            v['criteria'][f'{k}_spread<={lim_c}sigma_obs_every_node'] = good
            ok1 = ok1 and good
        feas = [bool(b.get('closure_diag') and b['closure_diag'].get(
                    'feasible'))
                for n in nodes_out for b in n['builds'] if 'on' in b]
        v['criteria']['every_node_build_feasible'] = bool(feas) and all(feas)
        ok1 = ok1 and v['criteria']['every_node_build_feasible']
        res = [abs(b['closure_diag']['mass_residual_after'])
               for n in nodes_out for b in n['builds']
               if b.get('closure_diag') and 'mass_residual_after'
               in b['closure_diag']]
        v['closed_stack_max_abs_mass_residual'] = max(res, default=None)
        v['criteria'][f'closed_stack_residual<={MM_RESIDUAL:g}'] = (
            bool(res) and max(res) <= MM_RESIDUAL)
        ok1 = ok1 and v['criteria'][f'closed_stack_residual<={MM_RESIDUAL:g}']
        off = offset_vs_rock_density(nodes_out, sig, cfg)
        v['offset_vs_rock_density'] = off
        good = (off is not None and all(
            off[k]['max_abs_residual_sigma_obs'] is not None
            and off[k]['max_abs_residual_sigma_obs'] <= MM_OFFSET
            for k in GATED))
        v['criteria'][
            f'|offset - predicted|<={MM_OFFSET}sigma_every_node'] = good
        ok1 = ok1 and good
    else:
        for k in GATED:
            sm = v['sigma_model'][k]['sigma_model_over_sigma_obs']
            good = sm is not None and sm <= SPREAD_C
            v['criteria'][f'{k}_sigma_model<=0.5sigma_obs'] = good
            ok1 = ok1 and good
    sp_l = [n['stats']['on']['libration_deg']['spread_sigma_obs']
            for n in full]
    lib_ok = bool(sp_l) and all(x is not None and x <= lim_l
                                for x in sp_l)
    v['criteria'][f'libration_spread<={lim_l}sigma_all_nodes'] = lib_ok
    if have_cache:
        lq = []
        for n in nodes_out:
            q = n.get('local_quadratic')
            if not q or 'error' in q:
                lq.append(None)
            else:
                lq.append(all(q[k]['pass'] for k in GATED
                              if q[k]['jitter_std_sigma_obs'] is not None))
        v['criteria']['local_quadratic<=3x_jitter_std'] = (
            None if any(x is None for x in lq) else all(lq))
    else:
        v['criteria']['local_quadratic<=3x_jitter_std'] = 'SKIPPED (no --cache)'
    v['criteria']['posterior_median_shift<0.2_sigma'] = (
        'TODO: later step (sigma_model on vs off posterior runs)')
    # reporting only: C20-C22 jitter correlation vs rho_obs
    rho_obs = ((cfg or {}).get('metadata', {}).get(
        'observable_correlations', {}) or {}).get('C20,C22')
    da, db = [], []
    for n in nodes_out:
        if n.get('_dev'):
            da += n['_dev'][0]
            db += n['_dev'][1]
    pooled = (float(np.corrcoef(da, db)[0, 1])
              if len(da) >= 3 and np.std(da) > 0 and np.std(db) > 0 else None)
    v['c20_c22_jitter_corr'] = {
        'rho_obs_in_config': rho_obs, 'pooled': pooled,
        'per_node': {f"({n['i_zb']},{n['i_w']})": n['c20_c22_corr']
                     for n in nodes_out},
        'flag_differs_strongly': (None if pooled is None or rho_obs is None
                                  else bool(abs(pooled - float(rho_obs))
                                            > 0.3)),
        'note': 'runner applies rho_obs to the INFLATED sigmas; flag = '
                '|pooled - rho_obs| > 0.3 (reporting only, not gating)'}
    lq_c = v['criteria']['local_quadratic<=3x_jitter_std']
    if failed or not full:
        v['overall'] = ('INCONCLUSIVE (build failures: see '
                        'failed_builds_by_node)')
    elif not (ok1 and lib_ok) or lq_c is False:
        v['overall'] = 'FAIL'
    elif lq_c is not True:
        v['overall'] = 'INCONCLUSIVE (local-fit check not run)'
    else:
        v['overall'] = 'PASS (all but posterior-median check, TODO)'
    return v


def offset_vs_rock_density(nodes_out, sig, cfg):
    """Criterion 3 (reviewer 2026-10-22). Per node, at the finest deltaP:
    measured offset = mean_nPoros[mass_moi - rock_density] (C2=0.5, in
    sigma_obs); predicted offset = (moi_target - MoI of the rock_density-closed
    stack) * dX/dC  [sign: mass_moi at target C_rd reproduces the
    rock_density stack exactly, so mass_moi(Ct) - rock_density = X(Ct) -
    X(C_rd) = (Ct - C_rd) dX/dC; the reviewer text wrote (C_rd - Ct)], dX/dC from mass_moi at moi_target +- 5e-4 on the same
    builds. PASS if |measured - predicted| <= 0.1 sigma at every node
    (C20, C22). Reports the signed and largest-absolute offset per node."""
    ct = float(cfg['metadata']['isostasy']['moi_target'])
    res = {}
    for k in OBS_NAMES:
        per = []
        for n in nodes_out:
            ok = [b for b in n['builds'] if 'on' in b
                  and b.get('ref_rock_density')]
            if not ok:
                continue
            fine = min(b['deltaP'] for b in ok)
            fb = [b for b in ok if b['deltaP'] == fine]
            meas = np.mean([(b['on'][k] - b['ref_rock_density'][k])
                            for b in fb]) / sig[k]
            pred = np.mean([(ct - b['ref_rock_density_moi']) * b['dXdC'][k]
                            for b in fb]) / sig[k]
            per.append({'node': f"({n['i_zb']},{n['i_w']})",
                        'measured_sigma_obs': float(meas),
                        'predicted_sigma_obs': float(pred),
                        'residual_sigma_obs': float(meas - pred)})
        meas = [p['measured_sigma_obs'] for p in per]
        resid = [abs(p['residual_sigma_obs']) for p in per]
        res[k] = {'per_node': per,
                  'largest_abs_offset_sigma_obs': (
                      float(max(abs(x) for x in meas)) if meas else None),
                  'max_abs_residual_sigma_obs': (
                      float(max(resid)) if resid else None)}
    return res


def frozen_check(cache):
    """(M_stack - M)/M of the frozen branch's reduced two-layer stack (the
    same reduction as MCMCRunner._frozen_two_layer_stack, taken BEFORE its
    mass-neutral step). I-F2 permits up to ~one outer shell's deficit
    (~0.1-0.7% M); not previously measured."""
    res = []
    for s in cache.get('frozen_structures') or []:
        if s is None:
            continue
        r = np.asarray(s['r_m'], float)
        rho = np.asarray(s['rho'], float)
        ph = np.asarray(s['phases'])
        o = np.argsort(r)
        r, rho, ph = r[o], rho[o], ph[o]
        k = np.isfinite(r) & np.isfinite(rho) & (r > 0)
        r, rho, ph = r[k], rho[k], ph[k]
        it = np.where(ph >= 50)[0]
        if it.size == 0 or int(it[-1]) >= r.size - 1:
            continue
        li = int(it[-1])
        rin = np.concatenate(([0.0], r))
        vol = r ** 3 - rin[:-1] ** 3

        def vw(a, b):
            g = vol[a:b + 1] > 0
            return np.sum(rho[a:b + 1][g] * vol[a:b + 1][g]) / np.sum(
                vol[a:b + 1][g])
        rr = np.array([r[li], r[-1]])
        rh = np.array([vw(0, li), vw(li + 1, r.size - 1)])
        M = float(np.sum(rh * (4 * np.pi / 3) * (rr ** 3 - np.array(
            [0.0, rr[0]]) ** 3)))
        res.append((M - float(s['Mtot_kg'])) / float(s['Mtot_kg']))
    if not res:
        return {'n': 0}
    a = np.array(res)
    return {'n': int(a.size), 'min': float(a.min()), 'max': float(a.max()),
            'median': float(np.median(a)),
            'min_abs': float(np.abs(a).min()),
            'max_abs': float(np.abs(a).max())}


def render_md(out):
    m, v = out['meta'], out['verdict']
    L = ['# Enceladus mass-closure jitter gate', '',
         f"Overall: **{v['overall']}**", '',
         f"builds={m['n_builds']} deltaP/{m['dp_factors']} "
         f"nPoros={m['nporos']} mean build {m['mean_build_wall_s']} s, "
         f"total {m['total_wall_s']} s; cache={m['cache']}", '',
         '## Criteria', '']
    for k, x in v['criteria'].items():
        L.append(f'- {k}: {x}')
    L += ['', '## sigma_model = max(std, |bias|)', '',
          json.dumps(v['sigma_model'], indent=1, default=float), '',
          '## C20-C22 jitter correlation', '',
          json.dumps(v['c20_c22_jitter_corr'], default=float), '',
          '## Failed builds', '', json.dumps(v['failed_builds_by_node']), '',
          '## Bias (production - finest), sigma_obs', '']
    for n in out['nodes']:
        L.append(f"- ({n['i_zb']},{n['i_w']}): "
                 f"{json.dumps(n['bias'], default=float)}")
    if m.get('closure_mode') == 'mass_moi':
        L += ['', '## mass_moi closure per node (production build = first, '
              'dP/1 nPoros 8 if present; range over builds in brackets)', '',
              '| node | feasible | R_s-R_sil(PP) km | rho_s kg/m3 | '
              'rho_s-rho_sil(PP) | MoI stack before | MoI after (target) | '
              'PP CMR2 | max abs M resid |', '|' + '---|' * 9]
        for n in out['nodes']:
            ds = [b for b in n['builds'] if b.get('closure_diag')]
            if not ds:
                continue
            prod = next((b for b in ds if b['deltaP'] == max(
                x['deltaP'] for x in ds) and b['nPoros'] == min(
                x['nPoros'] for x in ds)), ds[0])
            d = prod['closure_diag']

            def rg(key, sc=1.0):
                x = [b['closure_diag'].get(key) for b in ds
                     if b['closure_diag'].get(key) is not None]
                return (f"[{min(x) * sc:.3f}, {max(x) * sc:.3f}]"
                        if x else 'n/a')
            feas = all(b['closure_diag'].get('feasible') for b in ds)
            L.append(
                f"| ({n['i_zb']},{n['i_w']}) | {feas} | "
                f"{d.get('dR_s_vs_pp_m', float('nan')) / 1e3:.3f} "
                f"{rg('dR_s_vs_pp_m', 1e-3)} | "
                f"{d.get('rho_s', float('nan')):.1f} {rg('rho_s')} | "
                f"{d.get('drho_s_vs_pp', float('nan')):.1f} | "
                f"{d.get('moi_stack_before', float('nan')):.5f} "
                f"{rg('moi_stack_before')} | "
                f"{d.get('moi_stack_after', float('nan')):.5f} | "
                f"{prod.get('pp_CMR2')} | "
                f"{max((abs(b['closure_diag'].get('mass_residual_after', 0.0)) for b in ds), default=0):.2e} |")
        L += ['', 'closure infeasible reasons: ' + json.dumps({
            f"({n['i_zb']},{n['i_w']})": sorted({
                b['closure_diag'].get('reason') for b in n['builds']
                if b.get('closure_diag')
                and not b['closure_diag'].get('feasible')})
            for n in out['nodes']}), '',
            '## Offset mass_moi - rock_density (finest deltaP), measured vs '
            'predicted = (moi_target - MoI_rd) dX/dC, sigma_obs', '']
        for k in GATED:
            o = (out['verdict'].get('offset_vs_rock_density') or {}).get(k)
            if not o:
                continue
            L.append(f"- {k}: largest |offset| "
                     f"{o['largest_abs_offset_sigma_obs']}, max |meas-pred| "
                     f"{o['max_abs_residual_sigma_obs']}")
            for p_ in o['per_node']:
                L.append(f"    - {p_['node']}: signed offset "
                         f"{p_['measured_sigma_obs']:+.4f}, predicted "
                         f"{p_['predicted_sigma_obs']:+.4f}, residual "
                         f"{p_['residual_sigma_obs']:+.4f}")
        L.append('')
    L += ['', '## Frozen two-layer stack (M_stack-M)/M',
          json.dumps(out.get('frozen_check')), '',
          '', '## Per node: spread / std in sigma_obs (closure ON | OFF)', '',
          '| node (i,j) | zb | w | n | C20 sp | C22 sp | C30 sp | lib sp | '
          'C20 sp OFF | C22 sp OFF | lib sp OFF |', '|' + '---|' * 11]

    def f(x):
        return 'n/a' if x is None else f'{x:.3f}'
    for n in out['nodes']:
        s, o = n['stats']['on'], n['stats']['off']
        L.append(f"| ({n['i_zb']},{n['i_w']}) | {n['zb_km']:.2f} | "
                 f"{n['wOcean_ppt']:.3f} | {n['n_ok']} | "
                 f"{f(s['C20']['spread_sigma_obs'])} | "
                 f"{f(s['C22']['spread_sigma_obs'])} | "
                 f"{f(s['C30']['spread_sigma_obs'])} | "
                 f"{f(s['libration_deg']['spread_sigma_obs'])} | "
                 f"{f(o['C20']['spread_sigma_obs'])} | "
                 f"{f(o['C22']['spread_sigma_obs'])} | "
                 f"{f(o['libration_deg']['spread_sigma_obs'])} |")
    L += ['', '## Local quadratic check (cache neighbourhood, closure ON)', '']
    for n in out['nodes']:
        L.append(f"- ({n['i_zb']},{n['i_w']}): "
                 f"{json.dumps(n['local_quadratic'], default=float)}")
    L += ['', 'TODO (later step): posterior-median check, sigma_model on vs '
          'off (<0.2 posterior sigma).', '']
    return '\n'.join(L)


def scan_cache(args):
    """--scan-cache: run EVERY cached ocean node through the mass_moi closure
    (no rebuilds) and report: feasible fraction, rho_s range, R_s - R_sil(PP)
    range, and the rejects with reasons. REPORTING ONLY: no loose physical
    rejects (rho_s window etc.) are applied -- those come after this scan."""
    if not args.cache:
        raise SystemExit('--scan-cache needs --cache')
    with open(Path(args.cache).expanduser(), 'rb') as f:
        cache = pickle.load(f)
    if 'zb_km_grid' in cache:
        zb = np.asarray(cache['zb_km_grid'], float)
        w = np.asarray(cache['wOcean_ppt_grid'], float)
    else:       # 1-D (e.g. Tb-grid smoke) cache: index only
        n_ = len(cache['structures'])
        zb, w = np.full(n_, np.nan), np.array([np.nan])
    nw = w.size
    runner, _, _ = make_runners('mass_moi')
    rows, rejects, n_total, n_none, n_nonocean = [], [], 0, 0, 0
    for idx, st in enumerate(cache['structures']):
        n_total += 1
        i, j = divmod(idx, nw)
        if st is None:
            n_none += 1
            continue
        if not np.any(np.asarray(st['phases']) == 0):
            n_nonocean += 1
            continue
        runner.structure_data = st
        runner.last_mass_closure = None
        g = runner._derive_gravity_isostatic({'compensation_C2': 0.5})
        d = runner.last_mass_closure
        if d is None:
            reason = 'not_reached_closure (reduction/geometry reject)'
        elif not d.get('feasible'):
            reason = d.get('reason')
        elif g is None:
            reason = 'closure_feasible_but_downstream_none'
        else:
            reason = None
        if reason:
            rejects.append({'i_zb': i, 'i_w': j, 'zb_km': float(zb[i]),
                            'wOcean_ppt': float(w[j]), 'reason': reason,
                            'a': (d or {}).get('a'), 'b': (d or {}).get('b'),
                            'pre_closure_residual_frac':
                                (d or {}).get('pre_closure_residual_frac')})
        else:
            rows.append({'i_zb': i, 'i_w': j, 'rho_s': d['rho_s'],
                         'dR_s_km': d['dR_s_vs_pp_m'] / 1e3,
                         'rho_s_pp': d['rho_sil_pp'],
                         'moi_before': d['moi_stack_before']})
    n_oc = len(rows) + len(rejects)

    def rng(key):
        x = np.array([r[key] for r in rows], float)
        return (None if x.size == 0 else
                {'min': float(x.min()), 'median': float(np.median(x)),
                 'max': float(x.max())})
    reasons = {}
    for r in rejects:
        reasons[r['reason']] = reasons.get(r['reason'], 0) + 1
    out = {'cache': args.cache, 'n_slots': n_total, 'n_none': n_none,
           'n_frozen_or_no_ocean_skipped': n_nonocean, 'n_ocean': n_oc,
           'n_feasible': len(rows),
           'feasible_fraction': (len(rows) / n_oc if n_oc else None),
           'rho_s_kgm3': rng('rho_s'), 'rho_sil_pp_kgm3': rng('rho_s_pp'),
           'R_s_minus_R_sil_pp_km': rng('dR_s_km'),
           'moi_stack_before': rng('moi_before'),
           'reject_reasons': reasons, 'rejects': rejects}
    out_dir = Path(args.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / 'mass_moi_cache_scan.json').write_text(
        json.dumps(out, indent=1, default=float))
    L = ['# mass_moi cache scan (reporting only; no loose rejects applied)',
         '', f"cache: {args.cache}",
         f"slots {n_total} (None {n_none}, frozen/no-ocean skipped "
         f"{n_nonocean}), ocean nodes {n_oc}, feasible {len(rows)} "
         f"({out['feasible_fraction']})", '',
         f"rho_s: {out['rho_s_kgm3']}", f"rho_sil(PP): "
         f"{out['rho_sil_pp_kgm3']}",
         f"R_s - R_sil(PP) km: {out['R_s_minus_R_sil_pp_km']}",
         f"MoI stack before: {out['moi_stack_before']}", '',
         f"reject reasons: {reasons}", '', '## Rejects', '']
    L += [f"- ({r['i_zb']},{r['i_w']}) zb={r['zb_km']:.2f} "
          f"w={r['wOcean_ppt']:.3f}: {r['reason']}" for r in rejects]
    (out_dir / 'mass_moi_cache_scan.md').write_text('\n'.join(L))
    print('\n'.join(L[:14]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--cache')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--nodes', help='"i,j;i,j;..." cache indices '
                    '(default: 10 fixed nodes incl. (70,18),(77,10))')
    ap.add_argument('--dp-factors', help='comma list, default 1,2,4')
    ap.add_argument('--nporos', help='comma list, default 8,12,16')
    ap.add_argument('--frozen-check', action='store_true',
                    help='with --cache: report frozen two-layer stack mass '
                         'residual range')
    ap.add_argument('--closure-mode', choices=('rock_density', 'mass_moi'),
                    default='mass_moi',
                    help='closure under test (default mass_moi; limits: '
                         'C20/C22 spread <=0.1, libration <=0.05 sigma, '
                         'residual ~1e-15, all feasible, offset <=0.3)')
    ap.add_argument('--jobs', type=int, default=1)
    ap.add_argument('--self-test', action='store_true',
                    help='1 node (70,18), deltaP/{1,2} x nPoros {8} unless '
                         'overridden; no cache needed')
    ap.add_argument('--scan-cache', action='store_true',
                    help='with --cache: run every cached ocean node through '
                         'the mass_moi closure (no rebuilds); report '
                         'feasible fraction, rho_s and R_s-R_sil ranges, '
                         'rejects+reasons; exits without running the gate')
    args = ap.parse_args()
    if args.scan_cache:
        scan_cache(args)
        return
    if not args.self_test and not args.cache:
        print('NOTE: no --cache: local-quadratic check will be SKIPPED '
              '(verdict INCONCLUSIVE at best).', file=sys.stderr)
    run(args)


if __name__ == '__main__':
    main()
