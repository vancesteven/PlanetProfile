#!/usr/bin/env python
"""C3 mass-deficit study for the Enceladus (zb, w) production cache.

Scientific-reviewer condition C3 (r7 follow-up), verbatim:
  "Before quoting evidence: rebuild the ~5 nodes with the largest mass
   deficit using a finer silicate radial step; compare C20/C22/libration
   from the existing cache with the finer-step rebuilds. The shift must be
   <=0.1 sigma_obs; otherwise set a bound or correct the deficit. Also map
   the deficit over (zb, w) to check for jagged noise that would show up in
   interpolation."

WHAT THE "SILICATE RADIAL STEP" ACTUALLY IS (read before trusting a PASS)
-------------------------------------------------------------------------
The one-signed mass deficit comes from Silicates.py:~98
(`indsSilValid = np.where(Mtot_kg <= Planet.Bulk.M_kg)[0]`, then `[0]`):
the solver builds one candidate silicate/core profile per HYDROSPHERE layer
(nProfiles = nSurfIce - iSilStart + nOceanMax; the silicate top radius is a
hydrosphere-layer boundary) and takes the first whose total mass does not
exceed the target. The attainable masses are therefore quantised by the
radial thickness of ONE hydrosphere/ocean layer at the rock-ocean boundary,
which is set by `Planet.Ocean.deltaP` (pressure step of the ocean grid,
Enceladus template = 0.1 MPa, PPEnceladus.py:35). `Planet.Steps.nSilMax`
(PPEnceladus.py:28, 120) only sets the number of sub-layers INSIDE a given
candidate silicate body (Geophysical.py:791), so it does NOT change which
candidate is selected and cannot, by itself, reduce the deficit.

So the default refinement knob here is `Ocean.deltaP / refine_factor`
(applied through the builder's existing `ocean_overrides`, nothing else
touched). `--refine-knob nSilMax` is also supported (via the builder's
`planet_overrides={'Steps': <copy with nSilMax*factor>}`) for completeness;
`--refine-knob both` applies both. Note a deltaP change also refines the
ocean EOS pressure grid for the node, which is part of what "finer" means
here; the reproduction control and the cache differ ONLY by that knob.

RULED GATE (reviewer C3): refinement only SHRINKS the deficit, so the gate
is on the shift extrapolated to zero deficit,
  dobs_full = dobs * def_repro / (def_repro - def_refined)   (factor 4),
|dobs_full| <= 0.1 sigma_obs, with a linearity check against factor 2 (flag
if the two extrapolations differ >20%). Asserted: Tb and zb identical
between repro and refined builds (Tb root-find uses only IceLayers), refined
deficit >= repro deficit (nested candidates), |deficit| <= the one-layer
mass bound, and refined residual well below repro (else a second error
source). Factors must keep deltaP at <=3 decimals (HydroEOS.py:1017).
nSilMax is run on the single worst node only, as an ungated diagnostic.

Same builder entry point as the production cache: this script reuses the
exact per-node call in cache_builder.build_zbw_grid_cache
(cache_builder.py ~2047) -- build_single_structure with
  ocean_overrides = config ocean_overrides + {wOcean_ppt: w}
  bulk_overrides  = bulk_overrides_from_config(cfg) + {zb_approximate_km,
                    zbTol_km = 0.4 * ZB_TOL_KM}
  do_overrides    = {ICEIh_THICKNESS: True}
  extrap_ocean=False, Tb placeholder 272.0 K
and the production script's config loader / constants
(enceladus_zbw_production_build.py).

Observables use the campaign's own forward model, never reimplemented:
  MCMCRunner._derive_gravity_isostatic  (mcmc_runner.py:2240) -> C20, C22, C30
  MCMCRunner._derive_libration_deg      (mcmc_runner.py:2620) -> libration_deg
evaluated at compensation_C2 in {0, 0.5, 1} (the one nuisance they consume;
the worst case over C2 is what is gated). sigma_obs = the config's
`observables` second entry (raw observational sigma; the model-inflated
sigma is also reported).

Usage on the Mac:
  HOME=~/scratch/pphome python plans/scripts/enceladus_zbw_mass_deficit_study.py \
     --cache ~/scratch/enceladus_zbw_production_cache.pkl \
     --out-dir ~/scratch/c3_mass_deficit [--n-nodes 5]
Map only (no rebuilds):  add --map-only
Exercise the code path with no cache (builds 1 node):  --self-test
"""
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

from enceladus_zbw_production_build import (  # noqa: E402
    CONFIG_PATH, ZB_TOL_KM, _load_config)

TB_PLACEHOLDER_K = 272.0
C2_GRID = (0.0, 0.5, 1.0)
PASS_SIGMA = 0.1
OBS_NAMES = ('C20', 'C22', 'C30', 'libration_deg')


# ---------------------------------------------------------------- cache I/O

def _mass_res(s):
    if s is None:
        return np.nan
    v = s.get('mass_residual_frac')
    if v is None and s.get('Mtot_achieved_kg') is not None:
        v = s['Mtot_achieved_kg'] / s['Mtot_kg'] - 1.0
    return np.nan if v is None else float(v)


def _grid(cache):
    zb = np.asarray(cache['zb_km_grid'], float)
    w = np.asarray(cache['wOcean_ppt_grid'], float)
    Z = np.full((zb.size, w.size), np.nan)
    for i in range(zb.size):
        for j in range(w.size):
            Z[i, j] = _mass_res(cache['structures'][i * w.size + j])
    return zb, w, Z


# --------------------------------------------------------------- 1. the map

def _interp_dev(x, f, axis):
    """Deviation of each interior point from the linear interpolation of its
    two neighbours along `axis`, respecting non-uniform spacing x."""
    f = np.moveaxis(f, axis, 0)
    xi = x.reshape((-1,) + (1,) * (f.ndim - 1))
    dev = np.full(f.shape, np.nan)
    h0 = xi[1:-1] - xi[:-2]
    h1 = xi[2:] - xi[1:-1]
    lin = (f[:-2] * h1 + f[2:] * h0) / (h0 + h1)
    dev[1:-1] = f[1:-1] - lin
    return np.moveaxis(dev, 0, axis)


def _sign_changes(x, f, axis):
    """Count slope sign changes along `axis`, over triples with all 3 finite.
    Returns (n_changes, n_triples)."""
    f = np.moveaxis(f, axis, 0)
    s = np.diff(f, axis=0) / np.diff(x).reshape((-1,) + (1,) * (f.ndim - 1))
    a, b = s[:-1], s[1:]
    ok = np.isfinite(a) & np.isfinite(b) & (a != 0) & (b != 0)
    return int(np.sum(ok & (np.sign(a) != np.sign(b)))), int(np.sum(ok))


def map_and_roughness(zb, w, Z, out_dir):
    lw = np.log10(w)
    rows = []
    for i in range(zb.size):
        for j in range(w.size):
            rows.append((i, j, zb[i], w[j], Z[i, j]))
    with open(out_dir / 'mass_residual_map.csv', 'w') as f:
        f.write('i_zb,i_w,zb_km,wOcean_ppt,mass_residual_frac\n')
        for i, j, a, b, v in rows:
            f.write(f'{i},{j},{a:.6f},{b:.8g},'
                    f'{"" if not np.isfinite(v) else repr(float(v))}\n')

    rough = {}
    fin = Z[np.isfinite(Z)]
    rng = float(fin.max() - fin.min()) if fin.size else float('nan')
    for name, x, axis in (('along_zb', zb, 0), ('along_log10w', lw, 1)):
        dev = _interp_dev(x, Z, axis)
        absdev = np.abs(dev[np.isfinite(dev)])
        with np.errstate(invalid='ignore', divide='ignore'):
            relmap = np.abs(dev) / np.abs(Z)
        rel = relmap[np.isfinite(relmap)]
        nsc, ntr = _sign_changes(x, Z, axis)
        rough[name] = {
            'n_interior_points': int(absdev.size),
            'abs_dev_median': float(np.median(absdev)) if absdev.size else None,
            'abs_dev_p90': float(np.percentile(absdev, 90)) if absdev.size else None,
            'abs_dev_max': float(absdev.max()) if absdev.size else None,
            'abs_dev_max_over_range': (float(absdev.max() / rng)
                                       if absdev.size and rng > 0 else None),
            'rel_dev_median': float(np.median(rel)) if rel.size else None,
            'rel_dev_p90': float(np.percentile(rel, 90)) if rel.size else None,
            'rel_dev_max': float(rel.max()) if rel.size else None,
            'slope_sign_changes': nsc,
            'slope_triples': ntr,
            # smooth field -> ~0; independent noise -> ~0.67
            'slope_sign_change_fraction': (nsc / ntr if ntr else None),
        }
    fr = [rough[k]['slope_sign_change_fraction'] for k in rough
          if rough[k]['slope_sign_change_fraction'] is not None]
    rough['jagged_flag'] = bool(fr and max(fr) > 0.4)
    rough['jagged_flag_rule'] = ('max slope-sign-change fraction > 0.4 '
                                 '(smooth ~0, white noise ~0.67). A '
                                 'quantised-step deficit is sawtooth by '
                                 'construction, so judge this against the '
                                 'observable shifts, not in isolation.')
    rough['range'] = {'min': float(fin.min()) if fin.size else None,
                      'max': float(fin.max()) if fin.size else None,
                      'n_finite': int(fin.size), 'n_masked': int(Z.size - fin.size)}

    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 6.5))
        M = np.ma.masked_invalid(Z)
        pc = ax.pcolormesh(lw, zb, M, shading='nearest', cmap='viridis')
        ax.set_xlabel('log10 wOcean_ppt')
        ax.set_ylabel('zb  [km]')
        ax.set_title('Enceladus ocean-node mass_residual_frac '
                     '(Mtot_achieved/Mtot_target - 1); masked = None node')
        fig.colorbar(pc, ax=ax, label='mass_residual_frac')
        fig.tight_layout()
        fig.savefig(out_dir / 'mass_residual_heatmap.png', dpi=150)
        plt.close(fig)
    except Exception as exc:  # plotting must never block the study
        rough['png_error'] = f'{type(exc).__name__}: {exc}'
    return rough


# ------------------------------------------------------- 2. node selection

def pick_nodes(zb, w, Z, n_nodes, n_control):
    flat = [(Z[i, j], i, j) for i in range(zb.size) for j in range(w.size)
            if np.isfinite(Z[i, j])]
    flat.sort()
    worst = [(i, j, 'worst') for _, i, j in flat[:n_nodes]]
    med = float(np.median([v for v, _, _ in flat]))
    taken = {(i, j) for i, j, _ in worst}
    ctrl = []
    for v, i, j in sorted(flat, key=lambda t: abs(t[0] - med)):
        if (i, j) in taken:
            continue
        ctrl.append((i, j, 'control_median'))
        taken.add((i, j))
        if len(ctrl) >= n_control:
            break
    return worst + ctrl


# ------------------------------------------------------ 3. builds (rebuild)

def _production_kwargs(cfg, zb_km, w_ppt):
    """The exact per-node build_single_structure arguments of
    build_zbw_grid_cache (cache_builder.py ~2010-2060)."""
    from PlanetProfile.Inference.cache_builder import bulk_overrides_from_config
    ocean = dict(cfg['ocean_overrides'])
    ocean.pop('wOcean_ppt', None)
    ocean['wOcean_ppt'] = float(w_ppt)
    bulk = {**bulk_overrides_from_config(cfg),
            'zb_approximate_km': float(zb_km),
            'zbTol_km': 0.4 * ZB_TOL_KM}
    return dict(planet_template_module=cfg['planet_template_module'],
                Tb_K=float(TB_PLACEHOLDER_K),
                ocean_overrides=ocean, bulk_overrides=bulk,
                planet_overrides=None,
                do_overrides={'ICEIh_THICKNESS': True},
                extrap_ocean=False)


def build_node(cfg, zb_km, w_ppt, dp_factor=1.0, nsil_factor=1.0):
    """Build one node. dp_factor divides Ocean.deltaP (ruled C3 knob);
    nsil_factor multiplies Steps.nSilMax (diagnostic only). The refined
    deltaP must be representable at 3 dp (HydroEOS.py:1017 rounds the EOS
    deltaP to 3 dp), otherwise the request is rejected."""
    from copy import deepcopy
    import importlib
    from PlanetProfile.Inference.cache_builder import (
        build_single_structure, _ocean_node_zb_km)
    kw = _production_kwargs(cfg, zb_km, w_ppt)
    applied = {}
    if dp_factor != 1.0:
        mod = importlib.import_module(kw['planet_template_module'])
        dP0 = float(mod.Planet.Ocean.deltaP)
        dP1 = dP0 / dp_factor
        if abs(round(dP1, 3) - dP1) > 1e-12:
            raise SystemExit(
                f'deltaP {dP0}/{dp_factor} = {dP1} has more than 3 decimals; '
                f'HydroEOS.py:1017 rounds deltaP to 3 dp -- rejected.')
        kw['ocean_overrides']['deltaP'] = round(dP1, 3)
        applied['Ocean.deltaP'] = (dP0, round(dP1, 3))
    if nsil_factor != 1.0:
        mod = importlib.import_module(kw['planet_template_module'])
        steps = deepcopy(mod.Planet.Steps)
        n0 = int(steps.nSilMax)
        steps.nSilMax = int(round(n0 * nsil_factor))
        kw['planet_overrides'] = {'Steps': steps}
        applied['Steps.nSilMax'] = (n0, int(steps.nSilMax))
    t0 = time.time()
    s = build_single_structure(**kw)
    s['zb_km_actual'] = _ocean_node_zb_km(s)
    s['_build_wall_s'] = time.time() - t0
    s['_applied_knob'] = applied
    return s


def one_layer_bound(s):
    """4B: predicted one-layer mass bound
    4 pi r_s^2 dr (rho_rock_top - rho_ocean_bot) / M, with r_s the rock top
    radius and dr the bottom ocean layer thickness read from the structure
    (PP outer-edge convention: layer i spans [r[i-1], r[i]])."""
    r = np.asarray(s['r_m'], float)
    rho = np.asarray(s['rho'], float)
    ph = np.asarray(s['phases'])
    o = np.argsort(r)
    r, rho, ph = r[o], rho[o], ph[o]
    oc = np.where(ph == 0)[0]
    if oc.size == 0 or oc[0] == 0:
        return float('nan')
    f = int(oc[0])
    dr = r[f] - r[f - 1]
    return float(4 * np.pi * r[f - 1] ** 2 * dr * (rho[f - 1] - rho[f])
                 / s['Mtot_kg'])


# --------------------------------------------------------- 4. observables

def make_runner(cfg_path):
    from PlanetProfile.Inference.inference_core import InferenceConfig
    from PlanetProfile.Inference.mcmc_runner import MCMCRunner
    config = InferenceConfig.from_json(cfg_path)
    runner = object.__new__(MCMCRunner)
    runner.config = config
    return runner


def observables(runner, struct):
    """{obs: {C2: value}} through the production methods."""
    runner.structure_data = struct
    out = {k: {} for k in OBS_NAMES}
    for c2 in C2_GRID:
        theta = {'compensation_C2': c2}
        g = runner._derive_gravity_isostatic(theta)
        lib = runner._derive_libration_deg(theta)
        out['C20'][c2] = np.nan if g is None else float(g[(2, 0)])
        out['C22'][c2] = np.nan if g is None else float(g[(2, 2)])
        out['C30'][c2] = np.nan if g is None else float(g[(3, 0)])
        out['libration_deg'][c2] = np.nan if lib is None else float(lib)
    return out


def sigma_obs(cfg):
    obs = cfg['observables']
    add = cfg['metadata'].get('sigma_model_add', {})
    sig = {k: float(obs[k][1]) for k in OBS_NAMES}
    sig_infl = {k: float(np.hypot(sig[k], add.get(k, 0.0))) for k in OBS_NAMES}
    return sig, sig_infl


def dsigma(a, b, sig):
    """max over C2 of |a-b|/sigma per observable."""
    return {k: float(np.nanmax([abs(a[k][c] - b[k][c]) for c in C2_GRID])
                     / sig[k]) if all(np.isfinite([a[k][c] - b[k][c]
                                                  for c in C2_GRID])) else None
            for k in OBS_NAMES}


# ------------------------------------------------------------------- main

FACTORS = (2.0, 4.0)   # ruled: gate at 4, linearity check vs 2
GATE_FACTOR = 4.0


def _extrapolate(o_ref, o_rep, d_rep, d_ref, sig):
    """Delta_obs_full = Delta_obs * def_rep / (def_rep - def_ref), per
    observable, signed, worst |.| over compensation_C2, in sigma_obs."""
    scale = d_rep / (d_rep - d_ref)
    out = {}
    for k in OBS_NAMES:
        v = [(o_ref[k][c] - o_rep[k][c]) * scale / sig[k] for c in C2_GRID]
        out[k] = (None if not np.all(np.isfinite(v))
                  else float(max(v, key=abs)))
    return out


def study_nodes(cfg, cache, picks, zb, w, args, out_dir, sig, sig_infl):
    runner = make_runner(CONFIG_PATH)
    results = []
    for n_done, (i, j, role) in enumerate(picks):
        zb_km, w_ppt = float(zb[i]), float(w[j])
        cached = cache['structures'][i * w.size + j]
        print(f'\n=== node ({i},{j}) zb={zb_km:.4f} w={w_ppt:.4f} [{role}] '
              f'cached mass_residual={_mass_res(cached):+.6f}', flush=True)
        rec = {'i_zb': i, 'i_w': j, 'zb_km': zb_km, 'wOcean_ppt': w_ppt,
               'role': role, 'findings': []}
        try:
            repro = build_node(cfg, zb_km, w_ppt)
            ref = {f: build_node(cfg, zb_km, w_ppt, dp_factor=f)
                   for f in FACTORS}
            diag = None
            if n_done == 0 and role == 'worst':
                diag = build_node(cfg, zb_km, w_ppt, dp_factor=GATE_FACTOR,
                                  nsil_factor=2.0)   # "both", diagnostic
        except SystemExit:
            raise
        except Exception as exc:
            rec['error'] = f'{type(exc).__name__}: {exc}'
            print('  BUILD FAILED:', rec['error'])
            results.append(rec)
            continue

        oc, orp = observables(runner, cached), observables(runner, repro)
        orf = {f: observables(runner, ref[f]) for f in FACTORS}
        d_c, d_r = _mass_res(cached), _mass_res(repro)
        d_f = {f: _mass_res(ref[f]) for f in FACTORS}
        rec.update({
            'mass_residual_cached': d_c, 'mass_residual_reproduced': d_r,
            'mass_residual_refined': {str(f): d_f[f] for f in FACTORS},
            'deltaP_refined': {str(f): ref[f]['_applied_knob'] for f in FACTORS},
            'zb_actual': {'cached': float(cached.get('zb_km_actual', np.nan)),
                          'repro': repro['zb_km_actual'],
                          **{f'f{int(f)}': ref[f]['zb_km_actual']
                             for f in FACTORS}},
            'Tb_K': {'cached': float(cached.get('Tb_K', np.nan)),
                     'repro': float(repro.get('Tb_K', np.nan)),
                     **{f'f{int(f)}': float(ref[f].get('Tb_K', np.nan))
                        for f in FACTORS}},
            'obs_cached': oc, 'obs_reproduced': orp,
            'obs_refined': {str(f): orf[f] for f in FACTORS},
            'dsigma_repro_minus_cached': dsigma(orp, oc, sig),
            'build_wall_s': {'repro': repro['_build_wall_s'],
                             **{f'f{int(f)}': ref[f]['_build_wall_s']
                                for f in FACTORS}},
        })
        rec['reproduces_cache'] = bool(
            abs(d_r - d_c) < 1e-8 and all(
                v is not None and v < 1e-3
                for v in rec['dsigma_repro_minus_cached'].values()))

        # (3) root-find identity assertion -- a finding if it fails.
        for f in FACTORS:
            dTb = ref[f].get('Tb_K', np.nan) - repro.get('Tb_K', np.nan)
            dzb = ref[f]['zb_km_actual'] - repro['zb_km_actual']
            rec[f'dTb_f{int(f)}'], rec[f'dzb_f{int(f)}_km'] = float(dTb), float(dzb)
            if not (dTb == 0.0 and abs(dzb) <= 1e-9):
                rec['findings'].append(
                    f'ROOT-FIND NOT IDENTICAL at factor {f:g}: dTb={dTb:.3e} K, '
                    f'dzb={dzb:.3e} km (expected 0: Tb root-find uses only '
                    f'IceLayers, LayerPropagators.py:622-640). Observable '
                    f'shifts at this factor include a zb change; NOT corrected.')
        # (4a) nesting: refined deficit must be >= repro deficit.
        for f in FACTORS:
            if not d_f[f] >= d_r - 1e-12:
                rec['findings'].append(
                    f'NESTING VIOLATED at factor {f:g}: refined deficit '
                    f'{d_f[f]:+.3e} < repro {d_r:+.3e} (candidate set should '
                    f'be nested).')
        if any('NESTING' in x for x in rec['findings']) or any(
                abs(d_r - d_f[f]) < 1e-15 for f in FACTORS):
            rec['findings'].append('Extrapolation undefined (no/negative '
                                   'deficit change); gate FAILS.')
            ext = {f: {k: None for k in OBS_NAMES} for f in FACTORS}
        else:
            ext = {f: _extrapolate(orf[f], orp, d_r, d_f[f], sig)
                   for f in FACTORS}
        rec['extrapolated_dsigma_full'] = {str(f): ext[f] for f in FACTORS}
        # linearity between factor 2 and 4
        lin = {}
        for k in OBS_NAMES:
            a2, a4 = ext[2.0][k], ext[4.0][k]
            if a2 is None or a4 is None:
                lin[k] = None
            elif max(abs(a2), abs(a4)) < 0.01:
                lin[k] = {'ratio': None, 'ok': True,
                          'note': 'both < 0.01 sigma; ratio meaningless'}
            else:
                ratio = a2 / a4 if a4 != 0 else float('inf')
                lin[k] = {'ratio': float(ratio),
                          'ok': bool(abs(ratio - 1.0) <= 0.2)}
        rec['linearity_f2_over_f4'] = lin
        rec['linearity_ok'] = all(v is None or v['ok'] for v in lin.values())
        # 4b one-layer bound (repro structure, outer-edge convention)
        bound = one_layer_bound(repro)
        rec['one_layer_bound'] = bound
        rec['deficit_within_bound'] = bool(np.isfinite(bound)
                                           and abs(d_r) <= bound * 1.01)
        if not rec['deficit_within_bound']:
            rec['findings'].append(
                f'|deficit| {abs(d_r):.3e} exceeds predicted one-layer bound '
                f'{bound:.3e} (1% slack): deficit is NOT explained by '
                f'quantisation alone.')
        # 4d floor
        rec['floor_flag'] = bool(abs(d_f[GATE_FACTOR]) >= 3e-4
                                 or abs(d_f[GATE_FACTOR]) > 0.5 * abs(d_r))
        if rec['floor_flag']:
            rec['findings'].append(
                f'FLOOR: refined (x{GATE_FACTOR:g}) |residual| '
                f'{abs(d_f[GATE_FACTOR]):.3e} did not fall well below the '
                f'repro {abs(d_r):.3e} (or stays >= 3e-4): a second, '
                f'non-quantisation mass error source exists.')
        # diagnostic nSilMax (worst node #1 only, not gated)
        if diag is not None:
            od = observables(runner, diag)
            rec['diagnostic_deltaP4_plus_nSilMax2'] = {
                'knob': diag['_applied_knob'],
                'mass_residual': _mass_res(diag),
                'dsigma_vs_deltaP4_alone': dsigma(od, orf[GATE_FACTOR], sig),
                'note': 'DIAGNOSTIC, not gated'}
        g = ext[GATE_FACTOR]
        if not rec['linearity_ok']:
            rec['findings'].append(
                'LINEARITY: extrapolations from factors 2 and 4 disagree by '
                '>20% (see linearity_f2_over_f4); extrapolation unreliable.')
        rec['pass'] = bool(
            all(v is not None and abs(v) <= PASS_SIGMA for v in g.values())
            and not rec['findings'])
        print('  mass_res cached/repro/f2/f4: '
              f'{d_c:+.6f} / {d_r:+.6f} / {d_f[2.0]:+.6f} / {d_f[4.0]:+.6f}')
        print('  dTb f2/f4 = %.3e/%.3e K; dzb f2/f4 = %.3e/%.3e km; '
              'reproduces_cache=%s' % (rec['dTb_f2'], rec['dTb_f4'],
                                      rec['dzb_f2_km'], rec['dzb_f4_km'],
                                      rec['reproduces_cache']))
        print('  one-layer bound %.4e; |deficit| %.4e; within=%s'
              % (bound, abs(d_r), rec['deficit_within_bound']))
        for f in FACTORS:
            print(f'  extrapolated |dobs_full|/sigma (x{f:g}):',
                  {k: (None if v is None else round(v, 4))
                   for k, v in ext[f].items()})
        print('  linearity f2/f4 ok:', rec['linearity_ok'], '| PASS' if
              rec['pass'] else '| FAIL')
        print('  build wall s:', rec['build_wall_s'])
        for x in rec['findings']:
            print('  *** FINDING:', x)
        results.append(rec)
    return results


def write_reports(out_dir, meta, rough, results, sig, sig_infl):
    gate = [r for r in results if r.get('role') == 'worst']
    overall = bool(gate) and all(r.get('pass') for r in gate) \
        and not any('error' in r for r in gate)
    doc = {'meta': meta, 'sigma_obs': sig, 'sigma_obs_model_inflated': sig_infl,
           'pass_threshold_sigma': PASS_SIGMA, 'roughness': rough,
           'nodes': results, 'overall_pass_worst_nodes': overall}
    with open(out_dir / 'c3_mass_deficit_study.json', 'w') as f:
        json.dump(doc, f, indent=2, default=str)

    L = ['# C3 mass-deficit study', '',
         f"Status vocabulary: results below are `implemented, unverified` "
         f"until the user runs them on the Mac cache.", '',
         f"- cache: `{meta.get('cache')}`", f"- refine: {meta.get('refine')}",
         f"- sigma_obs (raw): {sig}", f"- sigma (model-inflated): {sig_infl}",
         f"- gate: shift extrapolated to zero deficit <= {PASS_SIGMA} sigma_obs "
         f"(max over compensation_C2 in {list(C2_GRID)})", '',
         '## Map / roughness', '']
    if rough:
        L += [f"- finite nodes {rough['range']['n_finite']}, masked "
              f"{rough['range']['n_masked']}, range "
              f"[{rough['range']['min']}, {rough['range']['max']}]",
              f"- jagged_flag: **{rough['jagged_flag']}** ({rough['jagged_flag_rule']})"]
        for k in ('along_zb', 'along_log10w'):
            r = rough[k]
            L.append(f"- {k}: abs dev median {r['abs_dev_median']}, "
                     f"p90 {r['abs_dev_p90']}, max {r['abs_dev_max']} "
                     f"({r['abs_dev_max_over_range']} of range); rel median "
                     f"{r['rel_dev_median']}, max {r['rel_dev_max']}; slope sign "
                     f"changes {r['slope_sign_changes']}/{r['slope_triples']}")
    L += ['', '## Nodes (gate: extrapolated-to-zero-deficit shift at x4, '
           'worst over C2, in sigma_obs)', '',
          '| node | role | zb | w | def cached | repro | x2 | x4 | bound | '
          + ' | '.join(f'{k}' for k in OBS_NAMES) + ' | PASS |',
          '|' + '---|' * (10 + len(OBS_NAMES))]
    findings = []
    for r in results:
        if 'error' in r:
            L.append(f"| ({r['i_zb']},{r['i_w']}) | {r['role']} | "
                     f"{r['zb_km']:.3f} | {r['wOcean_ppt']:.3f} | ERROR: "
                     f"{r['error']} |")
            continue
        e = r['extrapolated_dsigma_full'][str(GATE_FACTOR)]
        m = r['mass_residual_refined']
        L.append(
            f"| ({r['i_zb']},{r['i_w']}) | {r['role']} | {r['zb_km']:.3f} | "
            f"{r['wOcean_ppt']:.3f} | {r['mass_residual_cached']:+.5f} | "
            f"{r['mass_residual_reproduced']:+.5f} | {m['2.0']:+.5f} | "
            f"{m['4.0']:+.5f} | {r['one_layer_bound']:.5f} | "
            + ' | '.join('n/a' if e[k] is None else f'{e[k]:+.4f}'
                         for k in OBS_NAMES)
            + f" | {'PASS' if r['pass'] else 'FAIL'} |")
        findings += [f"({r['i_zb']},{r['i_w']}): {x}" for x in r['findings']]
    L += ['', '## Findings', ''] + ([f'- {x}' for x in findings]
                                     or ['- none'])
    L += ['', '## Diagnostics', '']
    for r in results:
        if 'diagnostic_deltaP4_plus_nSilMax2' in r:
            L.append(f"- nSilMax x2 on top of deltaP/4, node "
                     f"({r['i_zb']},{r['i_w']}) (not gated): "
                     f"{r['diagnostic_deltaP4_plus_nSilMax2']}")
    L += ['', f"**Overall (worst-deficit nodes only): "
          f"{'PASS' if overall else 'FAIL / incomplete'}**"]
    (out_dir / 'c3_mass_deficit_study.md').write_text('\n'.join(L) + '\n')
    return overall


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--cache')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--n-nodes', type=int, default=5)
    ap.add_argument('--n-control', type=int, default=2)
    ap.add_argument('--map-only', action='store_true')
    ap.add_argument('--self-test', action='store_true',
                    help='no cache: build ONE node (zb=--self-test-zb, '
                         'w=--self-test-w) as the "cached" node, then run '
                         'the full repro+refined+observables path on it.')
    ap.add_argument('--self-test-zb', type=float, default=25.0)
    ap.add_argument('--self-test-w', type=float, default=10.0)
    args = ap.parse_args()

    out_dir = Path(args.out_dir).expanduser()
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = _load_config()
    sig, sig_infl = sigma_obs(cfg)

    if args.self_test:
        print('SELF-TEST: building the stand-in "cached" node ...', flush=True)
        s0 = build_node(cfg, args.self_test_zb, args.self_test_w)
        cache = {'zb_km_grid': np.array([args.self_test_zb]),
                 'wOcean_ppt_grid': np.array([args.self_test_w]),
                 'structures': [s0]}
        # NOTE: the stand-in is itself a production-setting build, so the
        # reproduction control must match it exactly (a real determinism
        # check on this machine).
    else:
        if not args.cache:
            ap.error('--cache is required (or use --self-test)')
        with open(Path(args.cache).expanduser(), 'rb') as f:
            cache = pickle.load(f)

    zb, w, Z = _grid(cache)
    rough = map_and_roughness(zb, w, Z, out_dir)
    print(f"Map: {rough['range']}")
    for k in ('along_zb', 'along_log10w'):
        r = rough[k]
        print(f"  {k}: abs_dev max={r['abs_dev_max']} "
              f"(frac of range {r['abs_dev_max_over_range']}), slope sign "
              f"change fraction={r['slope_sign_change_fraction']}")
    print(f"  jagged_flag={rough['jagged_flag']}")

    meta = {'cache': args.cache or 'SELF-TEST', 'self_test': args.self_test,
            'refine': f'Ocean.deltaP /{list(FACTORS)} (gate x{GATE_FACTOR:g}); nSilMax x2 diagnostic',
            'n_nodes': args.n_nodes, 'zb_tol_km': ZB_TOL_KM,
            'config': CONFIG_PATH}
    if args.map_only:
        write_reports(out_dir, meta, rough, [], sig, sig_infl)
        return
    picks = pick_nodes(zb, w, Z, args.n_nodes, args.n_control)
    if args.self_test:
        picks = [(0, 0, 'worst')]
    results = study_nodes(cfg, cache, picks, zb, w, args, out_dir, sig,
                          sig_infl)
    overall = write_reports(out_dir, meta, rough, results, sig, sig_infl)
    print(f'\nOVERALL (worst nodes): {"PASS" if overall else "FAIL/incomplete"}'
          f'  -> {out_dir}/c3_mass_deficit_study.{{json,md}}')


if __name__ == '__main__':
    main()
