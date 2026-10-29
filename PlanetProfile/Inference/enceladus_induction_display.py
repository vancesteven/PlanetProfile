"""Enceladus induction-response DISPLAY reference (not an observable).

Config ruling (enceladus_cassini_isostasy_ocean6D_rhorock.json
metadata.induction_provenance; plans/active/enceladus-config-freeze.md):
induction is a display reference ONLY -- no conditioning and no support cut.
The inducing-field model scatter (2.4 nT in Bz, Saur et al. 2024 p. 24)
exceeds the 1-3 nT candidate signal, and no published amplitude plus
uncertainty exists.

For each posterior sample (zb_km, log10_wOcean_ppt) this maps the
precomputed per-node complex response Ae (built by
plans/scripts/enceladus_ae_sidecar.py) to the sample with the SAME bilinear
(zb, log10 w) weights as the ocean-branch structure lookup
(``grid_interp_2d.bilinear_weights``) and the same rule: any unbuilt corner
gives NaN, never a nearest-node fallback.

Two core-conductivity curves are carried: PlanetProfile's self-consistent
porous-rock sigma, and an IMPOSED 25 S/m core (Saur et al. 2024's
conductive-core scenario). The 25 S/m curve is an overlay, not a sampled
nuisance.

Field amplitudes use the convention of the existing synodic band
(``forward_models.induced_amplitude_band_nT``): |B_ind| = |Ae| x B_exc at the
model surface (R = 252.1 km).
"""
from __future__ import annotations

import pickle
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np

from .grid_interp_2d import bilinear_weights

# Saur et al. 2024 (PSJ 5:245) inducing-field amplitudes.
# Orbital (eccentricity) period: delta B_Theta = +/-3 e B0 = +/-5.0 nT (Sec 2.1).
ORBITAL_BEXC_NT = 5.0
# Synodic (PPO): "amplitudes of up to 1-2 nT" (Sec 2.2); non-stationary.
SYNODIC_BEXC_NT = (1.0, 2.0)
# Candidate induction signal, abstract: "perturbations of 1-3 nT".
SAUR_CANDIDATE_NT = (1.0, 3.0)
LABELS = {'orbital': 'true anomaly', 'synodic': 'synodic'}

CAPTION = (
    'Display reference only. Induction is NOT conditioned on and imposes no '
    'support cut: the inducing-field model scatter (2.4 nT in Bz, Saur et al. '
    '2024) exceeds the 1-3 nT candidate signal, and no published amplitude '
    'with uncertainty exists. '
    '(1) The orbital (eccentricity) period, 32.9 h, is the primary channel: '
    'Saturn drives ±5.0 nT in B_Θ there. The synodic (PPO) period, 15.56 h, '
    'has a non-stationary driver (time-varying amplitude, drifting phase), so '
    'only the band |Ae| × [1, 2] nT is meaningful. '
    '(2) Two core-conductivity curves: PlanetProfile\'s self-consistent '
    'porous-rock σ, and an IMPOSED 25 S/m core (Saur\'s conductive-core '
    'scenario) -- an overlay, not a sampled parameter. '
    '(3) Amplitude convention: Saur\'s A is the response at the surface '
    'through a 5-10 km insulating polar crust; ours is through the '
    'posterior mean shell (~27 km), so our A is roughly 25 % lower for the '
    'same ocean. '
    '(4) Radius conventions: Saur uses R_E = 256.6 km (also the gravity '
    'reference radius here); the interior model uses R = 252.1 km. '
    '(5) |B_ind| = |Ae| × B_exc at the model surface, the peak (polar) '
    'equivalent amplitude. '
    '(6) Upgrade to a conditioned observable only if a published amplitude '
    'and uncertainty exist at a named period, σ(B_ind) < 0.5 |B_ind|, and '
    'the plasma systematic is marginalized.')


def sidecar_path_for(cache_path) -> Path:
    return Path(str(cache_path) + '.ae_sidecar.pkl')


def _node_ocean_sigma(structs) -> np.ndarray:
    """Thickness-weighted mean ocean conductivity (S/m) per node; NaN if none."""
    out = np.full(len(structs), np.nan)
    for i, s in enumerate(structs):
        if not s or not s.get('D_ocean_km'):
            continue
        r = np.asarray(s['rSigChange_m'], float)
        sig = np.asarray(s['sigmaLayers_Sm'], float)
        r_top = s['R_body_m'] - 1e3 * s['zb_km_actual']
        r_floor = r_top - 1e3 * s['D_ocean_km']
        r_in = np.concatenate([[0.0], r[:-1]])
        dr = np.clip(np.minimum(r, r_top) - np.maximum(r_in, r_floor), 0, None)
        if dr.sum() > 0:
            out[i] = float(np.sum(sig * dr) / dr.sum())
    return out


def _node_array(grid: Dict[int, Optional[Dict[str, complex]]], n: int,
                label: str) -> np.ndarray:
    a = np.full(n, np.nan + 1j * np.nan, dtype=complex)
    for i, d in grid.items():
        if d and d.get(label) is not None:
            a[int(i)] = complex(d[label])
    return a


def _blend(node_vals: np.ndarray, zb_grid, w_grid, zb, w_ppt):
    out = np.full(len(zb), np.nan + 1j * np.nan, dtype=node_vals.dtype) \
        if np.iscomplexobj(node_vals) else np.full(len(zb), np.nan)
    for k, (z, w) in enumerate(zip(zb, w_ppt)):
        if not (zb_grid[0] <= z <= zb_grid[-1]
                and w_grid[0] <= w <= w_grid[-1]):
            continue
        corners, weights = bilinear_weights(zb_grid, w_grid, float(z), float(w))
        vals = node_vals[corners]
        if np.any(~np.isfinite(vals)):
            continue          # an unbuilt corner: no fallback
        out[k] = np.dot(weights, vals)
    return out


def induction_display(samples, param_names, structure_data,
                      sidecar_path) -> Optional[Dict[str, Any]]:
    """Per-sample display quantities, or None if the sidecar is absent.

    Returns ``{'sigma_ocean_Sm': [...], 'curves': {curve: {period:
    {'Ae_re', 'Ae_im'}}}, 'bexc_nT': {...}, 'caption': CAPTION}`` with
    NaN entries where a bracketing node is unbuilt.
    """
    p = Path(sidecar_path)
    if not p.exists():
        return None
    with open(p, 'rb') as f:
        payload = pickle.load(f)
    zb_grid = np.asarray(structure_data['zb_km_grid'], float)
    w_grid = np.asarray(structure_data['wOcean_ppt_grid'], float)
    structs = structure_data['structures']
    n = len(structs)
    if payload.get('key', {}).get('n') != n:
        raise ValueError(f'Ae display sidecar {p} does not match the cache '
                         f'({payload.get("key", {}).get("n")} vs {n} nodes)')
    th = np.asarray(samples, float)
    names = list(param_names)
    zb = th[:, names.index('zb_km')]
    w_ppt = 10.0 ** th[:, names.index('log10_wOcean_ppt')]
    out = {'sigma_ocean_Sm': _blend(_node_ocean_sigma(structs), zb_grid,
                                    w_grid, zb, w_ppt).tolist(),
           'curves': {}, 'caption': CAPTION,
           'bexc_nT': {'orbital': ORBITAL_BEXC_NT,
                       'synodic': list(SYNODIC_BEXC_NT)},
           'saur_candidate_nT': list(SAUR_CANDIDATE_NT),
           'periods_hr': payload.get('periods_hr', {}),
           'core25_Sm': payload.get('core25_Sm')}
    for curve, key in (('pp_core', 'ae_grid'), ('core_25Sm', 'ae_grid_core25')):
        grid = payload.get(key)
        if grid is None:
            continue
        out['curves'][curve] = {}
        for per, lab in LABELS.items():
            ae = _blend(_node_array(grid, n, lab), zb_grid, w_grid, zb, w_ppt)
            out['curves'][curve][per] = {'Ae_re': ae.real.tolist(),
                                         'Ae_im': ae.imag.tolist()}
    return out


def summarize(disp: Dict[str, Any]) -> Dict[str, Dict[str, Dict[str, float]]]:
    """Median and 16/84 % of |Ae|, phase lag Φ = -arg(Ae) (deg) and |B_ind| (nT)
    per curve and period (NaN rows excluded)."""
    def q(v):
        v = np.asarray(v, float)
        v = v[np.isfinite(v)]
        if v.size == 0:
            return {'p16': np.nan, 'p50': np.nan, 'p84': np.nan, 'n': 0}
        p16, p50, p84 = np.percentile(v, [16, 50, 84])
        return {'p16': float(p16), 'p50': float(p50), 'p84': float(p84),
                'n': int(v.size)}
    res = {}
    for curve, per_d in disp['curves'].items():
        res[curve] = {}
        for per, d in per_d.items():
            ae = np.asarray(d['Ae_re']) + 1j * np.asarray(d['Ae_im'])
            amp = np.abs(ae)
            row = {'A': q(amp), 'phi_deg': q(-np.degrees(np.angle(ae)))}
            if per == 'orbital':
                row['B_nT'] = q(amp * disp['bexc_nT']['orbital'])
            else:
                lo, hi = disp['bexc_nT']['synodic']
                row['B_lo_nT'], row['B_hi_nT'] = q(amp * lo), q(amp * hi)
            res[curve][per] = row
    return res
