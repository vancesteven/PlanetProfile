"""Per-sample radial interior profiles for the inference globe panel.

The sample picker selects one posterior draw; this module rebuilds that
draw's full radial structure by pushing its parameter vector through the
SAME bilinear (Tb x w) cache interpolation the likelihood used
(`forward_models.apply_parameters`), then renders standard
PlanetProfile-style property profiles (rho, T, P, sigma) and a
proportional layer-stack figure. All figures are matplotlib and are meant
to be displayed through Utilities.crisp_figs (vector SVG + cached
exports).

The structure cache is the file named by
``result.config.structure_cache_path`` (shipped with the app for the
amortized slots). When it is missing the profile view degrades to the
layer-stack figure, which needs only the per-sample derived arrays every
result already carries.
"""
from __future__ import annotations

import pickle
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import streamlit as st

from Utilities.globe_view import LAYER_COLORS, LAYER_LABELS

# PlanetProfile phase codes -> globe layer kinds (colors/labels shared
# with the 3D view). Codes: ocean=0, ice I..VI = 1..6, clathrate=30,
# silicate=50 (+ice-in-pore variants 5x), iron core=100.
_PHASE_KINDS = (
    (lambda p: p == 0, 'ocean'),
    (lambda p: p == 1, 'ice_Ih'),
    (lambda p: p in (2, 3), 'ice_III'),   # II folded into III color band
    (lambda p: p == 5, 'ice_V'),
    (lambda p: p == 6, 'ice_VI'),
    (lambda p: 30 <= p < 40, 'ice_Ih'),   # clathrate: show with Ih band
    (lambda p: 50 <= p < 100, 'silicate'),
    (lambda p: p >= 100, 'core'),
)


def phase_kind(p) -> str:
    p = int(p)
    for test, kind in _PHASE_KINDS:
        if test(p):
            return kind
    return 'silicate'


@st.cache_resource(show_spinner=False)
def _load_structure_cache(path_str: str):
    with open(path_str, 'rb') as f:
        return pickle.load(f)


def resolve_cache_path(result, parent_directory) -> Optional[Path]:
    p = getattr(result.config, 'structure_cache_path', '') or ''
    if not p:
        return None
    cand = Path(p)
    if not cand.is_absolute():
        cand = Path(parent_directory) / p
    return cand if cand.exists() else None


def _theta_for(result, sel_idx: Optional[int]) -> Dict[str, float]:
    samples = np.asarray(result.samples, float)
    names = list(result.param_names)
    row = (samples[sel_idx] if sel_idx is not None
           else np.median(samples, axis=0))
    return {n: float(row[j]) for j, n in enumerate(names)}


def profile_for_sample(result, sel_idx: Optional[int], parent_directory
                       ) -> Tuple[Optional[Dict], str]:
    """Radial structure for one posterior draw (or the marginal median
    when sel_idx is None). Returns (profile dict, note) — profile None
    with the reason in note when the cache can't serve it."""
    cache_path = resolve_cache_path(result, parent_directory)
    if cache_path is None:
        return None, ("structure cache not found for this result "
                      "(config.structure_cache_path missing or not "
                      "shipped) — layer stack only")
    try:
        structure_data = _load_structure_cache(str(cache_path))
    except Exception as e:
        return None, f"structure cache unreadable: {e}"

    theta = _theta_for(result, sel_idx)

    # v5/v6 thick-ice reparameterization: Tb_K is not sampled; derive it
    # by inverting the cached D_iceIh(Tb, w) field — the same inversion
    # the runner uses (grid_interp_2d).
    if 'D_iceIh_km' in theta and 'Tb_K' not in theta:
        try:
            from PlanetProfile.Inference.grid_interp_2d import (
                invert_d_iceIh_to_Tb, wOcean_ppt_from_theta)
            d_flat = [(None if s is None else s.get('D_iceIh_km'))
                      for s in structure_data.get('structures', [])]
            Tb = invert_d_iceIh_to_Tb(
                np.asarray(structure_data['Tb_K_grid'], float),
                np.asarray(structure_data['wOcean_ppt_grid'], float),
                d_flat, float(theta['D_iceIh_km']),
                wOcean_ppt_from_theta(theta))
            if Tb is None:
                return None, ("selected draw's ice thickness is outside "
                              "the cached (Tb, w) band — no structure")
            theta['Tb_K'] = float(Tb)
        except Exception as e:
            return None, f"D_iceIh → Tb inversion unavailable: {e}"

    try:
        from PlanetProfile.Inference.forward_models import (
            apply_parameters, UnservableSampleError)
        modified = apply_parameters(theta, structure_data)
    except UnservableSampleError:
        return None, ("selected draw falls in an unbuilt corner of the "
                      "(Tb, w) cache — no structure")
    except Exception as e:
        return None, f"structure interpolation failed: {e}"

    r_m = modified.get('r_m')
    if r_m is None:
        return None, ("cache did not yield a radial structure for this "
                      "parameterization")

    order = np.argsort(np.asarray(r_m, float))  # plot outward
    prof = {'r_km': np.asarray(r_m, float)[order] / 1e3}
    for key in ('rho', 'T_K', 'P_MPa', 'sigma_Sm', 'eta_Pa_base',
                'mu_Pa', 'phases'):
        v = modified.get(key)
        prof[key] = (np.asarray(v, float)[order]
                     if v is not None and np.size(v) == order.size else None)
    prof['Tb_K'] = theta.get('Tb_K')
    return prof, ''


def _phase_bands(r_km, phases) -> List[Tuple[float, float, str]]:
    """Contiguous (r_lo, r_hi, kind) runs of one material, outward."""
    bands = []
    if phases is None:
        return bands
    kinds = [phase_kind(p) for p in phases]
    lo = r_km[0]
    for i in range(1, len(kinds) + 1):
        if i == len(kinds) or kinds[i] != kinds[i - 1]:
            bands.append((lo, r_km[i - 1] if i == len(kinds) else r_km[i],
                          kinds[i - 1]))
            if i < len(kinds):
                lo = r_km[i]
    return bands


def build_profile_figure(prof: Dict, title: str):
    """2x2 standard-property profiles vs radius, phase-banded like the
    CLI hydrosphere plots: density, temperature, pressure, electrical
    conductivity."""
    import matplotlib.pyplot as plt

    r = prof['r_km']
    panels = [
        ('rho', r'Density $\rho$ (kg m$^{-3}$)', 'linear'),
        ('T_K', r'Temperature $T$ (K)', 'linear'),
        ('P_MPa', r'Pressure $P$ (MPa)', 'linear'),
        ('sigma_Sm', r'Conductivity $\sigma$ (S m$^{-1}$)', 'log'),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7), sharex=True)
    bands = _phase_bands(r, prof.get('phases'))
    seen_kinds: List[str] = []
    for (key, ylab, yscale), ax in zip(panels, axes.ravel()):
        v = prof.get(key)
        for lo, hi, kind in bands:
            ax.axvspan(lo, hi, color=_rgb01(LAYER_COLORS[kind]),
                       alpha=0.25, linewidth=0)
            if kind not in seen_kinds:
                seen_kinds.append(kind)
        if v is None:
            ax.text(0.5, 0.5, 'not in cache', transform=ax.transAxes,
                    ha='center', color='0.5')
        else:
            vv = np.where(np.asarray(v) > 0, v, np.nan) \
                if yscale == 'log' else v
            ax.plot(r, vv, color='k', linewidth=1.4)
            ax.set_yscale(yscale)
        ax.set_ylabel(ylab)
    for ax in axes[1]:
        ax.set_xlabel('Radius $r$ (km)')
    handles = [plt.Rectangle((0, 0), 1, 1,
                             color=_rgb01(LAYER_COLORS[k]), alpha=0.45)
               for k in seen_kinds]
    fig.legend(handles, [LAYER_LABELS.get(k, k) for k in seen_kinds],
               loc='lower center', ncol=min(4, len(seen_kinds) or 1),
               fontsize=8, frameon=False)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    return fig


def build_stack_figure(R_km: float, thicknesses: List[Tuple[str, float]],
                       title: str):
    """One proportional radial column, surface at top: each present
    material drawn to scale with its thickness labeled (globe colors)."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(4.2, 6.4))
    top = float(R_km)
    for kind, d in thicknesses:
        if d <= 0:
            continue
        ax.bar(0, d, width=0.9, bottom=top - d,
               color=_rgb01(LAYER_COLORS[kind]), edgecolor='white',
               linewidth=0.6)
        label = f"{LAYER_LABELS.get(kind, kind)}  {d:.0f} km"
        # annotate inside when the band is wide enough, else beside it
        if d > 0.045 * R_km:
            ax.text(0, top - d / 2, label, ha='center', va='center',
                    fontsize=8)
        else:
            y = top - d / 2
            ax.plot([0.46, 0.56], [y, y], color='0.4', linewidth=0.6)
            ax.text(0.58, y, label, ha='left', va='center', fontsize=7)
        top -= d
    ax.set_xlim(-0.6, 1.4)
    ax.set_ylim(0, R_km * 1.02)
    ax.set_xticks([])
    ax.set_ylabel('Radius (km)')
    ax.set_title(title, fontsize=10)
    for s in ('top', 'right', 'bottom'):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    return fig


def _rgb01(rgb_str: str):
    """'rgb(r, g, b)' (plotly convention in globe_view) -> matplotlib tuple."""
    vals = rgb_str[rgb_str.index('(') + 1:rgb_str.index(')')].split(',')
    return tuple(int(v) / 255.0 for v in vals)
