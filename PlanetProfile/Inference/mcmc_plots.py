"""
Body-agnostic MCMC visualisation helpers for PlanetProfile Bayesian inference.

Extracted from Test48_mcmc_andrade_yao2014.py so the same plot code can be
reused for any ocean world (Europa, Ganymede, Titan, ...) by passing
body-specific scalars as explicit arguments rather than reading module-level
constants.

All functions:
  - Accept all body-specific values as explicit arguments.
  - Accept an ``output_path`` argument (absolute path to the output PNG).
  - Save to ``output_path`` at dpi=150 (or as documented), log at INFO level,
    and close the figure before returning.
  - Import heavy libraries (matplotlib, corner, seaborn) inside the function
    body to keep module-level import fast.

Public functions
----------------
plot_corner
    Corner plot of the full posterior sample matrix.

plot_k2_scatter_by
    k2 Re/Im scatter coloured by an arbitrary array (D_hydro, f_sil, etc.)
    with Petricca-style 1-sigma / 2-sigma observation ellipses.

plot_ice_comparison
    Two-panel: (a) k2 scatter coloured by Ice Ih fraction of ice heating,
    (b) log10 Ice Ih vs Ice V heating scatter coloured by D_ocean.

plot_heating_vs_parameters
    2x5 scatter grid of log10 tidal heating power vs each parameter, plus an
    optional full-width cumulative heating-fraction bar below.

plot_mass_cmr2_diagnostics
    Two-panel: total mass vs D_hydro (left), CMR2 vs D_hydro (right) with
    observation bands.

plot_cmr2_surface
    pcolormesh of CMR2 across the (Tb, D_hydro) grid.

plot_tb_structure
    Ice Ih shell thickness vs Tb from the grid with posterior Tb histogram
    overlay.

plot_layers_vs_docean
    3-panel: (top) D_ocean posterior histogram, (middle) cumulative-thickness
    stackplot sorted by (f_sil ASC, D_ocean ASC), (bottom) per-phase tidal
    heating on the same x-ordering.

plot_structure_wedge
    Wedge diagram of the posterior interior structure with labels showing
    median layer thicknesses and 5/95-percentile uncertainty arcs.
"""

from __future__ import annotations

import importlib
import logging
import os
from copy import deepcopy
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helper — shared colour / label dictionaries used by multiple functions
# ---------------------------------------------------------------------------

_PHASE_COLORS_STACK: Dict[str, str] = {
    'Clath': '#D4F1F9',
    'Ih':    '#AEE1F8',
    '0':     '#1E90FF',
    'III':   '#C97BAE',
    'II':    '#B0E0E6',
    'V':     '#9B59B6',
    'VI':    '#6C3483',
    'Sil':   '#C8A96E',
    'Rock':  '#C8A96E',  # Alias for Silicate
    'Core':  '#8B5A2B',
}

_PHASE_LABELS_STACK: Dict[str, str] = {
    'Clath': 'Clathrate',
    'Ih':    'Ice Ih',
    '0':     'Ocean',
    'III':   'Ice III',
    'II':    'Ice II',
    'V':     'Ice V',
    'VI':    'Ice VI',
    'Sil':   'Silicate',
    'Rock':  'Rock',
    'Core':  'Core',
}

_FULL_PHASE_COLORS: Dict[str, str] = {
    'Clath': '#D4F1F9',
    'Ih':    '#AEE1F8',
    '0':     '#1E90FF',
    'III':   '#C97BAE',
    'II':    '#B0E0E6',
    'V':     '#9B59B6',
    'VI':    '#6C3483',
    'Sil':   '#C8A96E',
    'Rock':  '#C8A96E',
    'Core':  '#8B5A2B',
}

_FULL_PHASE_LABELS: Dict[str, str] = {
    'Clath': 'Clathrate',
    'Ih':    'Ice Ih',
    '0':     'Ocean',
    'III':   'Ice III',
    'II':    'Ice II',
    'V':     'Ice V',
    'VI':    'Ice VI',
    'Sil':   'Silicate',
    'Rock':  'Rock',
    'Core':  'Core',
}

_WEDGE_COLORS: Dict[str, str] = {
    'Ice Ih':     '#AEE1F8',
    'Ocean':      '#1E90FF',
    'Ice III':    '#C97BAE',
    'Ice II':     '#B0E0E6',
    'Ice V':      '#9B59B6',
    'Ice VI':     '#6C3483',
    'Silicate':   '#C8A96E',
    'Rock':       '#C8A96E',
    'Dense core': '#8B5A2B',
    'Core':       '#8B5A2B',
}


def _wedge_color_map() -> Dict[str, str]:
    """Return a color dict keyed by layer name using PP's canonical Color config.

    Falls back to _WEDGE_COLORS hex values if Color is unavailable or uninitialized.
    The ocean and silicate entries are sampled at the colormap midpoint (0.5).
    """
    try:
        from PlanetProfile.GetConfig import Color
        def _get(attr, fallback, cmap=False):
            val = getattr(Color, attr, None)
            if val is None:
                return fallback
            if cmap:
                return val(0.5) if callable(val) else fallback
            return val
        sil = _get('silCondCmap', _WEDGE_COLORS['Silicate'], cmap=True)
        core = _get('Fe', _WEDGE_COLORS['Dense core'])
        return {
            'Ice Ih':     _get('iceIcond', _WEDGE_COLORS['Ice Ih']),
            'Ocean':      _get('oceanCmap', _WEDGE_COLORS['Ocean'], cmap=True),
            'Ice II':     _get('iceII',  _WEDGE_COLORS['Ice II']),
            'Ice III':    _get('iceIII', _WEDGE_COLORS['Ice III']),
            'Ice V':      _get('iceV',   _WEDGE_COLORS['Ice V']),
            'Ice VI':     _get('iceVI',  _WEDGE_COLORS['Ice VI']),
            'Silicate':   sil,
            'Rock':       sil,
            'Dense core': core,
            'Core':       core,
        }
    except Exception:
        return dict(_WEDGE_COLORS)


# ---------------------------------------------------------------------------
# Default column-index mapping (Test48 / 10-D parameter layout)
# ---------------------------------------------------------------------------

_DEFAULT_PARAM_INDICES: Dict[str, int] = {
    'Tb':      5,
    'D_hydro': 6,
    'f_core':  9,
}


# ===========================================================================
# 1. Corner plot
# ===========================================================================

def plot_corner(
    samples: np.ndarray,
    labels: Sequence[str],
    title: str,
    output_path: str,
    quantiles: Tuple[float, ...] = (0.16, 0.5, 0.84),
) -> None:
    """Corner plot of the full posterior sample matrix.

    Handles degenerate columns (zero variance) by expanding the axis range by
    ±1 in native units so the ``corner`` library does not crash.

    Args:
        samples:     Array of shape (N_samples, N_params).
        labels:      Axis label for each column of ``samples``.
        title:       Figure suptitle string.
        output_path: Absolute path to the output PNG file.
        quantiles:   Quantile values to show as vertical lines on 1-D histograms.
    """
    import corner
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style='white', font_scale=0.9)

    cs_std = np.nanstd(samples, axis=0)
    cs_med = np.nanmedian(samples, axis=0)
    corner_range = [
        (
            float(np.nanmin(samples[:, i])) - 0.1 * abs(cs_std[i]),
            float(np.nanmax(samples[:, i])) + 0.1 * abs(cs_std[i]),
        )
        if cs_std[i] > 0
        else (float(cs_med[i]) - 1.0, float(cs_med[i]) + 1.0)
        for i in range(samples.shape[1])
    ]

    fig = corner.corner(
        samples,
        labels=list(labels),
        color='steelblue',
        range=corner_range,
        quantiles=list(quantiles),
        show_titles=True,
        title_fmt='.3f',
        title_kwargs={'fontsize': 10},
    )
    fig.suptitle(title, fontsize=14, y=1.02)

    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    log.info(f'Saved {output_path}')
    plt.close(fig)


# ===========================================================================
# 2. k2 scatter coloured by an arbitrary array
# ===========================================================================

def plot_k2_scatter_by(
    k2_results: np.ndarray,
    color_values: np.ndarray,
    colorbar_label: str,
    obs_re: float,
    obs_im: float,
    obs_re_err: float,
    obs_im_err: float,
    title: str,
    output_path: str,
    cmap: str = 'plasma_r',
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
) -> None:
    """k2 Re/Im scatter coloured by an arbitrary per-sample array.

    Draws Petricca-style 1-sigma (dashed) and 2-sigma (dotted) observation
    ellipses in red.

    Args:
        k2_results:      Array of shape (N, 2) with columns [Re(k2), Im(k2)].
                         Im(k2) is taken as abs() internally so negative values
                         are handled correctly.
        color_values:    Per-sample colour array, shape (N,).
        colorbar_label:  Label for the colour bar.
        obs_re:          Observed Re(k2) central value.
        obs_im:          Observed |Im(k2)| central value.
        obs_re_err:      1-sigma uncertainty on obs_re.
        obs_im_err:      1-sigma uncertainty on obs_im.
        title:           Axes title string.
        output_path:     Absolute path to the output PNG file.
        cmap:            Matplotlib colour map name (default ``'plasma_r'``).
        vmin:            Lower clamp for colour scale (None = data min).
        vmax:            Upper clamp for colour scale (None = data max).
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.patches import Ellipse

    sns.set_theme(style='white', font_scale=0.9)

    Re_arr = k2_results[:, 0]
    Im_arr = np.abs(k2_results[:, 1])

    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(Re_arr, Im_arr, c=color_values, cmap=cmap, s=8, alpha=0.6,
                    vmin=vmin, vmax=vmax)
    plt.colorbar(sc, ax=ax, label=colorbar_label)

    for n_sigma, ls, lw, lbl in [
        (1, '--', 2, r'obs 1$\sigma$'),
        (2, ':',  1, r'obs 2$\sigma$'),
    ]:
        ax.add_patch(Ellipse(
            (obs_re, obs_im),
            2 * n_sigma * obs_re_err,
            2 * n_sigma * obs_im_err,
            fill=False, color='red', linewidth=lw, linestyle=ls, label=lbl,
        ))

    ax.set_xlabel(r'$\mathrm{Re}(k_2)$')
    ax.set_ylabel(r'$|\mathrm{Im}(k_2)|$')
    ax.set_title(title)
    ax.legend()

    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    log.info(f'Saved {output_path}')
    plt.close(fig)


# ===========================================================================
# 3. Ice-phase comparison (two-panel)
# ===========================================================================

def plot_ice_comparison(
    k2_results: np.ndarray,
    heating_results: List[Dict[str, float]],
    d_ocean_eval: np.ndarray,
    obs_re: float,
    obs_im: float,
    obs_re_err: float,
    obs_im_err: float,
    output_path: str,
    w_floor: float = 1e8,
) -> None:
    """Two-panel ice-phase comparison figure.

    Panel (a): k2 scatter coloured by Ice Ih / (total ice heating) ratio.
               Green = Ih dominated, red = HP ice (V/VI) dominated.
    Panel (b): log10 Ice Ih vs log10 Ice V heating scatter, coloured by
               D_ocean.  Includes a 1:1 equal-power reference line.

    Args:
        k2_results:      Array (N, 2) of [Re(k2), Im(k2)] for eval samples.
        heating_results: List of per-sample dicts mapping phase name to watts.
        d_ocean_eval:    Per-sample ocean thickness in km, shape (N,).
        obs_re:          Observed Re(k2).
        obs_im:          Observed |Im(k2)|.
        obs_re_err:      1-sigma uncertainty on obs_re.
        obs_im_err:      1-sigma uncertainty on obs_im.
        output_path:     Absolute path to the output PNG file.
        w_floor:         Heating-power floor in Watts (default 1e8 W).
                         Samples below this threshold are excluded from panel (b).
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.patches import Ellipse

    sns.set_theme(style='white', font_scale=0.9)

    Re_arr = k2_results[:, 0]
    Im_arr = np.abs(k2_results[:, 1])

    Ih_pw  = np.array([h.get('Ih',  0.0) for h in heating_results])
    V_pw   = np.array([h.get('V',   0.0) for h in heating_results])
    III_pw = np.array([h.get('III', 0.0) for h in heating_results])
    VI_pw  = np.array([h.get('VI',  0.0) for h in heating_results])
    ice_tot = Ih_pw + III_pw + V_pw + VI_pw
    f_ih = np.where(ice_tot > w_floor, Ih_pw / (ice_tot + 1e-30), np.nan)
    valid_ice = np.isfinite(f_ih) & (ice_tot > w_floor)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # --- Panel (a) ---
    ax = axes[0]
    if np.any(valid_ice):
        sc = ax.scatter(Re_arr[valid_ice], Im_arr[valid_ice],
                        c=f_ih[valid_ice], cmap='RdYlGn', s=8, alpha=0.7,
                        vmin=0, vmax=1)
        plt.colorbar(sc, ax=ax, label='Ice Ih fraction of ice heating')
    else:
        ax.scatter(Re_arr, Im_arr, s=8, alpha=0.5, c='gray')

    for n_sigma, ls, lw, lbl in [
        (1, '--', 2, r'obs 1$\sigma$'),
        (2, ':',  1, r'obs 2$\sigma$'),
    ]:
        ax.add_patch(Ellipse(
            (obs_re, obs_im),
            2 * n_sigma * obs_re_err,
            2 * n_sigma * obs_im_err,
            fill=False, color='red', linewidth=lw, linestyle=ls, label=lbl,
        ))
    ax.set_xlabel(r'$\mathrm{Re}(k_2)$')
    ax.set_ylabel(r'$|\mathrm{Im}(k_2)|$')
    ax.set_title(r'$k_2$ posterior: Ice Ih vs HP ice dominance (green=Ih, red=V/VI)')
    ax.legend(fontsize=8)

    # --- Panel (b) ---
    ax = axes[1]
    ih_mask = valid_ice & (Ih_pw > w_floor) & (V_pw > w_floor)
    if np.any(ih_mask):
        sc = ax.scatter(
            np.log10(Ih_pw[ih_mask]),
            np.log10(V_pw[ih_mask]),
            c=d_ocean_eval[ih_mask],
            cmap='plasma', s=8, alpha=0.7,
        )
        plt.colorbar(sc, ax=ax, label=r'$D_\mathrm{ocean}$ (km)')
        mn = min(np.log10(Ih_pw[ih_mask]).min(), np.log10(V_pw[ih_mask]).min())
        mx = max(np.log10(Ih_pw[ih_mask]).max(), np.log10(V_pw[ih_mask]).max())
        ax.plot([mn, mx], [mn, mx], 'k--', lw=1, label='equal power')
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5,
                'No samples with both\nIce Ih and Ice V heating',
                ha='center', va='center', transform=ax.transAxes)
    ax.set_xlabel(r'$\log_{10}$ Ice Ih heating (W)')
    ax.set_ylabel(r'$\log_{10}$ Ice V heating (W)')
    ax.set_title('Ice Ih vs Ice V heating power')

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    log.info(f'Saved {output_path}')
    plt.close(fig)


# ===========================================================================
# 4. Heating scatter vs parameters (2x5 grid + optional cumulative bar)
# ===========================================================================

def plot_heating_vs_parameters(
    eval_samples: np.ndarray,
    heating_results: List[Dict[str, float]],
    param_labels: Sequence[str],
    extra_xvals: List[np.ndarray],
    extra_xlabels: List[str],
    output_path: str,
    w_floor: float = 1e8,
    cumulative_bar: bool = True,
    eval_d_ocean: Optional[np.ndarray] = None,
    title: Optional[str] = None,
) -> None:
    """2x5 scatter grid of log10 tidal heating power vs each parameter.

    Optionally adds a full-width cumulative stacked-fraction bar below.

    The scatter panels show one point per eval sample per active phase.
    An "active" phase is one that carries > ``w_floor`` W in at least 5 samples.

    The cumulative bar panel shows per-model heating fractions sorted by
    (f_sil ASC, D_ocean ASC) so cross-panel ordering matches
    ``plot_layers_vs_docean``.

    Args:
        eval_samples:   Array (N_eval, N_params) of parameter values for the
                        evaluated subset.
        heating_results: List (N_eval,) of dicts mapping phase name to watts.
        param_labels:   Labels for the first N columns of eval_samples.
                        Must have length >= 5 (columns 0..4 are first panel row).
                        Column 6 onward are second row.  See ``extra_xvals`` for
                        how to insert derived quantities (e.g. D_iceIh) into the
                        layout.
        extra_xvals:    List of additional x-value arrays to append after the
                        first 5 columns of eval_samples.  Typical usage: pass
                        ``[d_iceIh_eval]`` to insert the Ice Ih thickness
                        derived from the grid lookup.
        extra_xlabels:  Labels matching ``extra_xvals``.
        output_path:    Absolute path to the output PNG file.
        w_floor:        Heating floor in Watts; phases below this are omitted
                        from scatter panels (default 1e8 W).
        cumulative_bar: If True (default), add the full-width cumulative
                        heating-fraction bar as a third row.
        eval_d_ocean:   Per-eval-sample D_ocean (km) used for the secondary
                        sort key in the cumulative bar.  If None and
                        ``cumulative_bar`` is True, the sort uses only f_sil.
        title:          Figure suptitle.  If None, no suptitle is added.
    """
    import matplotlib.gridspec as gridspec
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style='white', font_scale=0.9)

    ALL_PHASES = ['Ih', 'III', 'V', 'VI', 'Sil', 'Core', 'Clath']
    phase_colors = {ph: _FULL_PHASE_COLORS.get(ph, '#cccccc') for ph in ALL_PHASES}

    heating_power = {
        ph: np.array([h.get(ph, 0.0) for h in heating_results])
        for ph in ALL_PHASES
    }
    active_phases = [ph for ph in ALL_PHASES
                     if np.sum(heating_power[ph] > 1e6) > 5]

    # Build x-value / label sequence for the 2×5 scatter grid
    # First 5 params come from eval_samples columns 0..4, then extra_xvals,
    # then eval_samples columns 6..9 (columns beyond 5).
    plot_xvals = (
        [eval_samples[:, i] for i in range(5)]
        + list(extra_xvals)
        + [eval_samples[:, i] for i in range(6, eval_samples.shape[1])]
    )
    plot_xlabels = (
        list(param_labels[:5])
        + list(extra_xlabels)
        + list(param_labels[6:])
    )

    n_scatter = 10  # always 2 rows × 5 cols
    if len(plot_xvals) < n_scatter:
        # Pad if fewer than 10 panels are available
        plot_xvals   = plot_xvals   + [None] * (n_scatter - len(plot_xvals))
        plot_xlabels = plot_xlabels + ['']   * (n_scatter - len(plot_xlabels))
    else:
        plot_xvals   = plot_xvals[:n_scatter]
        plot_xlabels = plot_xlabels[:n_scatter]

    # --- Layout ---
    if cumulative_bar:
        fig = plt.figure(figsize=(25, 13))
        gs = gridspec.GridSpec(3, 5, height_ratios=[1, 1, 1.1],
                               hspace=0.35, wspace=0.3)
        scatter_axes = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(5)]
        ax_cum = fig.add_subplot(gs[2, :])
    else:
        fig = plt.figure(figsize=(25, 9))
        gs = gridspec.GridSpec(2, 5, hspace=0.35, wspace=0.3)
        scatter_axes = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(5)]
        ax_cum = None

    # --- Scatter panels ---
    for ip, (xvals, xlabel) in enumerate(zip(plot_xvals, plot_xlabels)):
        ax = scatter_axes[ip]
        if xvals is None:
            ax.set_visible(False)
            continue
        for ph in active_phases:
            pw = heating_power[ph]
            mask = pw > w_floor
            if np.any(mask):
                ax.scatter(xvals[mask], np.log10(pw[mask]),
                           s=4, alpha=0.4, color=phase_colors[ph], label=ph)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'$\log_{10}$ Heating (W)')
        if ip == 0 and active_phases:
            ax.legend(fontsize=8, loc='best')

    # --- Cumulative heating-fraction bar ---
    if cumulative_bar and ax_cum is not None:
        stack_phases = ['Core', 'Sil', 'VI', 'V', 'III', '0', 'Ih', 'Clath']

        # Only sum individual phases, excluding the new aggregate '_W' keys 
        # (Silicate_W, HP_Ice_W, Ice_Ih_W) to avoid double-counting.
        agg_keys = {'Silicate_W', 'HP_Ice_W', 'Ice_Ih_W'}
        total_W = np.array([sum(v for k, v in h.items() if k not in agg_keys) 
                            for h in heating_results])
        safe_tot = np.where(total_W > 1e-30, total_W, 1e-30)
        fracs: Dict[str, np.ndarray] = {}
        for ph in stack_phases:
            fracs[ph] = np.array([h.get(ph, 0.0) for h in heating_results]) / safe_tot

        f_sil_per = fracs['Sil'].copy()

        if eval_d_ocean is not None:
            d_ocean_per = eval_d_ocean
        else:
            d_ocean_per = np.zeros(len(heating_results))

        order = np.lexsort((d_ocean_per, f_sil_per))
        n_models = len(order)
        x = np.arange(n_models)

        bottom = np.zeros(n_models)
        for ph in stack_phases:
            h_vals = fracs[ph][order]
            ax_cum.bar(x, h_vals, bottom=bottom, width=1.0,
                       color=_PHASE_COLORS_STACK[ph], edgecolor='none',
                       label=_PHASE_LABELS_STACK[ph])
            bottom += h_vals

        ax_cum.set_ylim(0, 1.0)
        ax_cum.set_xlim(0, n_models - 1)
        ax_cum.set_ylabel('Cumulative heating fraction', fontsize=11)
        ax_cum.set_title(
            r'Per-model heating fractions  (models sorted by $f_\mathrm{sil}$ ASC, '
            r'then $D_\mathrm{ocean}$ ASC — matches layers_vs_docean ordering)',
            fontsize=11,
        )

        n_ticks = 8
        tick_positions = np.linspace(0, n_models - 1, n_ticks, dtype=int)
        tick_labels = [
            f'{f_sil_per[order][i]:.2f}\n({d_ocean_per[order][i]:.0f} km)'
            for i in tick_positions
        ]
        ax_cum.set_xticks(tick_positions)
        ax_cum.set_xticklabels(tick_labels, fontsize=8)
        ax_cum.set_xlabel(
            r'$f_\mathrm{sil}$  ($D_\mathrm{ocean}$)   [models sorted ascending]',
            fontsize=11,
        )
        handles, labels = ax_cum.get_legend_handles_labels()
        ax_cum.legend(handles[::-1], labels[::-1],
                      loc='center left', bbox_to_anchor=(1.005, 0.5), fontsize=9)

    if title is not None:
        fig.suptitle(title, fontsize=14)

    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    log.info(f'Saved {output_path}')
    plt.close(fig)


# ===========================================================================
# 5. Mass and CMR2 diagnostics vs D_hydro
# ===========================================================================

def plot_mass_cmr2_diagnostics(
    d_hydro_values: np.ndarray,
    mtot_results: np.ndarray,
    cmr2_results: np.ndarray,
    obs_mtot: float,
    obs_mtot_err: float,
    obs_cmr2: float,
    obs_cmr2_err: float,
    output_path: str,
) -> None:
    """Two-panel mass and CMR2 diagnostics vs hydrosphere thickness.

    Panel (left):  Total mass vs D_hydro with observed mass band.
    Panel (right): CMR2 vs D_hydro with observed CMR2 band.

    The mass axis is scaled to 1e23 kg; the legend for the mass panel
    formats the uncertainty in units of 1e20 kg.

    Args:
        d_hydro_values: Per-sample D_hydro in km, shape (N,).
        mtot_results:   Per-sample total mass in kg, shape (N,).
        cmr2_results:   Per-sample CMR2 (dimensionless), shape (N,).
        obs_mtot:       Observed total mass in kg.
        obs_mtot_err:   1-sigma uncertainty on obs_mtot in kg.
        obs_cmr2:       Observed CMR2.
        obs_cmr2_err:   1-sigma uncertainty on obs_cmr2.
        output_path:    Absolute path to the output PNG file.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style='white', font_scale=0.9)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: total mass
    ax = axes[0]
    ax.scatter(d_hydro_values, mtot_results / 1e23, s=6, alpha=0.5, c='steelblue')
    ax.axhline(obs_mtot / 1e23, color='red', ls='--',
               label=f'obs ± {obs_mtot_err / 1e20:.1f}×10²⁰ kg')
    ax.axhspan(
        (obs_mtot - obs_mtot_err) / 1e23,
        (obs_mtot + obs_mtot_err) / 1e23,
        alpha=0.15, color='red',
    )
    ax.set_xlabel(r'$D_\mathrm{hydro}$ (km)')
    ax.set_ylabel(r'$M_\mathrm{tot}$ (×10²³ kg)')
    ax.set_title('Total Mass vs Hydrosphere Thickness')
    ax.legend()

    # Right: CMR2
    ax = axes[1]
    ax.scatter(d_hydro_values, cmr2_results, s=6, alpha=0.5, c='darkorange')
    ax.axhline(obs_cmr2, color='red', ls='--',
               label=f'obs={obs_cmr2}+/-{obs_cmr2_err}')
    ax.axhspan(obs_cmr2 - obs_cmr2_err, obs_cmr2 + obs_cmr2_err,
               alpha=0.15, color='red')
    ax.set_xlabel(r'$D_\mathrm{hydro}$ (km)')
    ax.set_ylabel('C/MR²')
    ax.set_title('CMR2 vs Hydrosphere Thickness')
    ax.legend()

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    log.info(f'Saved {output_path}')
    plt.close(fig)


# ===========================================================================
# 6. CMR2 surface over the (Tb, D_hydro) grid
# ===========================================================================

def plot_cmr2_surface(
    tb_vals: np.ndarray,
    d_vals: np.ndarray,
    grid_cache: Dict[Tuple[float, float], Dict],
    output_path: str,
) -> None:
    """pcolormesh of CMR2 across the (Tb, D_hydro) grid.

    Args:
        tb_vals:     1-D array of basal temperature grid values in K.
        d_vals:      1-D array of D_hydro grid values in km.
        grid_cache:  Dict keyed by (tb, d) tuples mapping to structure dicts.
                     Each structure dict should contain a ``'CMR2'`` entry.
        output_path: Absolute path to the output PNG file.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style='white', font_scale=0.9)

    CMR2_grid = np.full((len(d_vals), len(tb_vals)), np.nan)
    for j, d in enumerate(d_vals):
        for i, tb in enumerate(tb_vals):
            key = (float(tb), float(d))
            if key in grid_cache:
                CMR2_grid[j, i] = grid_cache[key].get('CMR2', np.nan)

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.pcolormesh(tb_vals, d_vals, CMR2_grid, cmap='viridis', shading='nearest')
    plt.colorbar(im, ax=ax, label='C/MR²')
    ax.set_xlabel(r'$T_b$ (K)')
    ax.set_ylabel(r'$D_\mathrm{hydro}$ (km)')
    ax.set_title('CMR2 across the grid')

    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    log.info(f'Saved {output_path}')
    plt.close(fig)


# ===========================================================================
# 7. Ice Ih shell thickness vs Tb with posterior histogram overlay
# ===========================================================================

def plot_tb_structure(
    tb_vals: np.ndarray,
    d_vals: np.ndarray,
    grid_cache: Dict[Tuple[float, float], Dict],
    samples: np.ndarray,
    output_path: str,
) -> None:
    """Ice Ih shell thickness vs Tb from the grid with posterior Tb overlay.

    The main line shows the maximum D_iceIh_km value seen at each Tb across
    all D_hydro entries (truncated no-ocean entries have D_iceIh < true ice
    thickness; taking the max recovers the physical relationship).

    The twin y-axis shows a normalised histogram of the posterior Tb samples.

    Args:
        tb_vals:     1-D array of basal temperature grid values in K.
        d_vals:      1-D array of D_hydro grid values in km (used for the
                     grid lookup; not plotted directly).
        grid_cache:  Dict keyed by (tb, d) tuples.
        samples:     Full posterior sample array (N_samples, N_params).
                     Column 5 is assumed to be Tb_K (the default Test48 layout).
        output_path: Absolute path to the output PNG file.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style='white', font_scale=0.9)

    # Build Tb -> max D_iceIh lookup
    tb_to_diceIh: Dict[float, float] = {}
    for (tb, _d), struct in grid_cache.items():
        val = struct.get('D_iceIh_km', np.nan)
        if np.isfinite(val) and (tb not in tb_to_diceIh or val > tb_to_diceIh[tb]):
            tb_to_diceIh[tb] = val

    tb_sorted   = sorted(tb_to_diceIh)
    diceIh_line = [tb_to_diceIh[tb] for tb in tb_sorted]

    fig, ax = plt.subplots(figsize=(8, 5))
    if any(np.isfinite(v) for v in diceIh_line):
        ax.plot(tb_sorted, diceIh_line, 'k-o', markersize=4)
        ax2 = ax.twinx()
        ax2.hist(samples[:, 5], bins=30, alpha=0.3, color='blue',
                 density=True, label=r'Posterior $T_b$')
        ax2.set_ylabel('Posterior density', color='blue')
        ax.set_xlabel(r'$T_b$ (K)')
        ax.set_ylabel('Ice Ih shell thickness (km)')
        ax.set_title('Ice Shell Thickness vs Basal Temperature')
        fig.tight_layout()

    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    log.info(f'Saved {output_path}')
    plt.close(fig)


# ===========================================================================
# 8. Layer structure vs D_ocean (3-panel)
# ===========================================================================

def plot_layers_vs_docean(
    samples: np.ndarray,
    eval_idx: np.ndarray,
    grid_cache: Dict[Tuple[float, float], Dict],
    heating_results: List[Dict[str, float]],
    output_path: str,
    R_body_km: float,
    body_name: str = 'Titan',
    param_indices: Optional[Dict[str, int]] = None,
    equil_heating_GW: Optional[float] = None,
    equil_heating_label: Optional[str] = None,
) -> None:
    """3-panel interior-structure figure sorted by (f_sil ASC, D_ocean ASC).

    Panel layout (shared x-axis for middle and bottom):
      Top:    D_ocean posterior histogram (independent x-axis).
      Middle: Cumulative-thickness stackplot — one column per eval model.
      Bottom: Per-phase tidal heating power on the same x-ordering.

    Args:
        samples:         Full posterior array (N_samples, N_params).
        eval_idx:        Integer indices into ``samples`` for the evaluated
                         subset, shape (N_eval,).
        grid_cache:      Dict keyed by (tb, d) tuples.
        heating_results: List (N_eval,) of per-phase heating dicts in W.
        output_path:     Absolute path to the output PNG file.
        R_body_km:       Body surface radius in km.
        body_name:       Body name string used in the title (default ``'Titan'``).
        param_indices:   Dict mapping semantic keys to column indices in
                         ``samples``.  Required keys: ``'Tb'``, ``'D_hydro'``,
                         ``'f_core'``.  If None, uses the Test48 defaults
                         ``{'Tb': 5, 'D_hydro': 6, 'f_core': 9}``.
    """
    import matplotlib.gridspec as gridspec
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_theme(style='white', font_scale=0.9)

    if param_indices is None:
        param_indices = _DEFAULT_PARAM_INDICES

    i_Tb      = param_indices['Tb']
    i_D_hydro = param_indices['D_hydro']
    i_f_core  = param_indices['f_core']

    PHASE_LIST = ['Core', 'Sil', 'VI', 'V', 'III', '0', 'Ih', 'Clath']

    # Extract per-model structural data
    model_data = []
    for ii, si in enumerate(eval_idx):
        tb       = samples[si, i_Tb]
        d_hydro  = samples[si, i_D_hydro]
        f_core_i = samples[si, i_f_core]

        # Find closest structure in grid (robust to correlated grids)
        best_pt = None
        min_dist = float('inf')
        for (ktb, kd) in grid_cache.keys():
            dist = ((ktb - tb)/5.0)**2 + ((kd - d_hydro)/50.0)**2
            if dist < min_dist:
                min_dist = dist
                best_pt = (ktb, kd)
        
        pt = grid_cache.get(best_pt)
        if pt is None:
            continue

        d_ih  = pt.get('D_iceIh_km',  0.0)
        d_oc  = pt.get('D_ocean_km',  0.0)
        d_iii = pt.get('D_iceIII_km', 0.0)
        d_v   = pt.get('D_iceV_km',   0.0)
        d_vi  = pt.get('D_iceVI_km',  0.0)
        d_hp  = pt.get('D_hp_ice_km', 0.0)
        if d_iii == 0 and d_v == 0 and d_vi == 0 and d_hp > 0:
            d_v = d_hp

        d_hydro_actual = d_ih + d_oc + d_iii + d_v + d_vi
        r_sil_km  = R_body_km - d_hydro_actual
        r_core_km = f_core_i * r_sil_km
        d_sil  = r_sil_km - r_core_km
        d_core = r_core_km

        h = heating_results[ii]
        total = sum(h.values())
        f_sil = h.get('Sil', 0.0) / max(total, 1e-30) if total > 1e-6 else 0.0
        heat_by_phase = {ph: h.get(ph, 0.0) for ph in PHASE_LIST}

        model_data.append({
            'D_ocean': d_oc,
            'Ih':  d_ih,
            '0':   d_oc,
            'III': d_iii,
            'V':   d_v,
            'VI':  d_vi,
            'Sil': d_sil,
            'Core': d_core,
            'idx': ii,
            'f_sil': f_sil,
            'heat': heat_by_phase,
            'total_heat': total,
        })

    if not model_data:
        log.warning('plot_layers_vs_docean: no valid samples, skipping.')
        return

    # Sort by ascending Ice Ih thickness — the most informative ordering for
    # ocean-world posteriors (warm-Tb thin-shell models on the left, cold-Tb
    # thick-shell models on the right). Test48's original f_sil ordering made
    # sense for the no-ocean Titan case where silicate dissipation dominated;
    # for Europa/Ganymede/Callisto with active oceans, ice-shell thickness
    # is the more interpretable abscissa.
    model_data.sort(key=lambda r: r.get('Ih', 0.0))
    d_ih_vals    = np.array([r.get('Ih', 0.0) for r in model_data])
    d_ocean_vals = np.array([r['D_ocean'] for r in model_data])
    f_sil_vals   = np.array([r['f_sil']  for r in model_data])
    no_ocean_frac = np.mean(d_ocean_vals < 0.5)

    fig = plt.figure(figsize=(10, 12))
    gs  = gridspec.GridSpec(3, 1, height_ratios=[1, 4, 3], hspace=0.15)
    ax_dens   = fig.add_subplot(gs[0])
    ax_struct = fig.add_subplot(gs[1])
    ax_heat   = fig.add_subplot(gs[2], sharex=ax_struct)

    # Top: D_ocean histogram
    ax_dens.hist(d_ocean_vals, bins=40, alpha=0.5, color='steelblue',
                 density=True, edgecolor='none')
    ax_dens.axvline(0.5, color='red', ls=':', lw=1.5, label='no-ocean boundary')
    ax_dens.set_xlabel(r'$D_\mathrm{ocean}$ (km)', fontsize=10)
    ax_dens.set_ylabel('Density', fontsize=9)
    ax_dens.set_title(
        f'{body_name} Interior Structure — All Posterior Models'
        f'  (no-ocean: {no_ocean_frac:.0%};  sorted by ascending $D_\\mathrm{{Ih}}$)',
        fontsize=12,
    )
    ax_dens.tick_params(labelsize=8)
    ax_dens.legend(fontsize=8)

    # Middle: stackplot. Drop phases that contribute zero across all models
    # so the legend only lists phases physically present (e.g., Europa under
    # Seawater has no HP ices; Callisto with NaCl 100 ppt has no HP ices and
    # may have no Fe core depending on the prior). Threshold 0.1 km handles
    # numerical noise without dropping genuinely-thin layers.
    _all_stack_phases = ['Core', 'Sil', 'VI', 'V', 'III', '0', 'Ih', 'Clath']
    stack_phases = [p for p in _all_stack_phases
                    if max(r.get(p, 0.0) for r in model_data) > 0.1]
    x      = np.arange(len(model_data))
    stacks = [np.array([r.get(p, 0.0) for r in model_data]) for p in stack_phases]
    colors = [_FULL_PHASE_COLORS.get(p, '#cccccc') for p in stack_phases]
    labels = [_FULL_PHASE_LABELS.get(p, p) for p in stack_phases]

    polys = ax_struct.stackplot(x, *stacks, colors=colors, labels=labels)
    for poly in polys:
        poly.set_linewidth(0)
        poly.set_edgecolor('none')

    ax_struct.axhline(R_body_km, color='k', ls='-', lw=0.5, alpha=0.3)
    ax_struct.set_ylabel('Cumulative thickness (km)', fontsize=12)
    ax_struct.set_ylim(0, R_body_km * 1.02)
    ax_struct.set_xlim(0, len(model_data) - 1)
    ax_struct.tick_params(labelsize=9)
    import matplotlib.pyplot as _plt
    _plt.setp(ax_struct.get_xticklabels(), visible=False)

    handles, lbls = ax_struct.get_legend_handles_labels()
    ax_struct.legend(handles[::-1], lbls[::-1],
                     loc='center left', bbox_to_anchor=(1.02, 0.5),
                     fontsize=9, frameon=True)

    # Bottom: per-phase heating
    for ph in PHASE_LIST:
        vals = np.array([r['heat'][ph] / 1e9 for r in model_data])
        mask = vals > 1e-3
        if np.any(mask):
            ax_heat.scatter(x[mask], vals[mask], s=8, alpha=0.5,
                            color=_FULL_PHASE_COLORS.get(ph, '#C8A96E'),
                            label=_FULL_PHASE_LABELS.get(ph, ph))

    ax_heat.set_yscale('log')
    ax_heat.set_ylabel('Phase heating power (GW)', fontsize=12)
    if equil_heating_GW is not None:
        lbl = equil_heating_label or f'equil. (~{equil_heating_GW:.0f} GW)'
        ax_heat.axhline(equil_heating_GW, color='gray', ls='--', lw=1, alpha=0.6, label=lbl)
    ax_heat.legend(fontsize=8, loc='lower left', ncol=2)
    ax_heat.tick_params(labelsize=9)
    ax_heat.set_title(
        'Per-Phase Tidal Heating  (x-axis shared with structure panel)',
        fontsize=11,
    )

    n_ticks = 8
    tick_positions = np.linspace(0, len(model_data) - 1, n_ticks, dtype=int)
    tick_labels = [
        f'{d_ih_vals[i]:.0f}\n({d_ocean_vals[i]:.0f} km)'
        for i in tick_positions
    ]
    ax_heat.set_xticks(tick_positions)
    ax_heat.set_xticklabels(tick_labels, fontsize=8)
    ax_heat.set_xlabel(
        r'$D_\mathrm{Ih}$ (km)   ($D_\mathrm{ocean}$)   '
        '[models sorted ascending in $D_\mathrm{Ih}$]',
        fontsize=11,
    )

    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    log.info(f'Saved {output_path}')
    _plt.close(fig)


# ===========================================================================
# 9. Wedge diagram of posterior interior structure
# ===========================================================================

def plot_structure_wedge(
    samples: np.ndarray,
    eval_idx: np.ndarray,
    grid_cache: Dict[Tuple[float, float], Dict],
    output_path: str,
    R_body_km: float,
    body_name: str = 'Titan',
    param_indices: Optional[Dict[str, int]] = None,
) -> None:
    """Wedge diagram of the posterior interior structure.

    Shows median layer radii as filled wedge patches with 5th/95th-percentile
    uncertainty arcs.  Matches PlanetProfile's ``PlotWedge()`` visual style.

    Args:
        samples:       Full posterior array (N_samples, N_params).
        eval_idx:      Integer indices into ``samples`` for the evaluated subset.
        grid_cache:    Dict keyed by (tb, d) tuples.
        output_path:   Absolute path to the output PNG file.
        R_body_km:     Body surface radius in km.
        body_name:     Body name string used in the title (default ``'Titan'``).
        param_indices: Dict mapping semantic keys to column indices in
                       ``samples``.  Required keys: ``'Tb'``, ``'D_hydro'``,
                       ``'f_core'``.  If None, uses the Test48 defaults
                       ``{'Tb': 5, 'D_hydro': 6, 'f_core': 9}``.
    """
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    from matplotlib.patches import Wedge as MplWedge

    if param_indices is None:
        param_indices = _DEFAULT_PARAM_INDICES

    i_Tb      = param_indices['Tb']
    i_D_hydro = param_indices['D_hydro']
    i_f_core  = param_indices['f_core']

    ANG1, ANG2 = 55, 125  # wedge angular extent (degrees)

    # Collect per-sample layer boundary radii
    r_iceIh_bot:  List[float] = []
    r_ocean_bot:  List[float] = []
    r_iceII_bot:  List[float] = []
    r_iceIII_bot: List[float] = []
    r_iceV_bot:   List[float] = []
    r_iceVI_bot:  List[float] = []
    r_sil_bot:    List[float] = []
    r_core_bot:   List[float] = []

    for i in eval_idx:
        tb       = samples[i, i_Tb]
        d_hydro  = samples[i, i_D_hydro]
        f_core_i = samples[i, i_f_core]

        # Find closest structure in grid (robust to correlated grids)
        best_pt = None
        min_dist = float('inf')
        for (ktb, kd) in grid_cache.keys():
            # Normalized distance (Tb ~ 10K range, D ~ 100km range)
            dist = ((ktb - tb)/5.0)**2 + ((kd - d_hydro)/50.0)**2
            if dist < min_dist:
                min_dist = dist
                best_pt = (ktb, kd)
        
        pt = grid_cache.get(best_pt)
        if pt is None:
            continue

        d_ih  = pt.get('D_iceIh_km',  0.0)
        d_oc  = pt.get('D_ocean_km',  0.0)
        d_ii  = pt.get('D_iceII_km',  0.0)
        d_iii = pt.get('D_iceIII_km', 0.0)
        d_v   = pt.get('D_iceV_km',   0.0)
        d_vi  = pt.get('D_iceVI_km',  0.0)
        d_hp  = pt.get('D_hp_ice_km', 0.0)
        if d_ii == 0 and d_iii == 0 and d_v == 0 and d_vi == 0 and d_hp > 0:
            d_v = d_hp

        r_ih_b  = R_body_km - d_ih
        r_oc_b  = r_ih_b - d_oc
        r_ii_b  = r_oc_b - d_ii
        r_iii_b = r_ii_b - d_iii
        r_v_b   = r_iii_b - d_v
        r_vi_b  = r_v_b - d_vi
        r_sil_b = r_vi_b           # top of silicate = bottom of last ice
        r_core_top = f_core_i * r_sil_b

        r_iceIh_bot.append(r_ih_b)
        r_ocean_bot.append(r_oc_b)
        r_iceII_bot.append(r_ii_b)
        r_iceIII_bot.append(r_iii_b)
        r_iceV_bot.append(r_v_b)
        r_iceVI_bot.append(r_vi_b)
        r_sil_bot.append(r_core_top)
        r_core_bot.append(0.0)

    if not r_iceIh_bot:
        log.warning('plot_structure_wedge: no valid samples, skipping.')
        return

    def pct(arr: List[float]) -> np.ndarray:
        return np.percentile(arr, [5, 50, 95])

    p_ih  = pct(r_iceIh_bot)
    p_oc  = pct(r_ocean_bot)
    p_ii  = pct(r_iceII_bot)
    p_iii = pct(r_iceIII_bot)
    p_v   = pct(r_iceV_bot)
    p_vi  = pct(r_iceVI_bot)
    p_sil = pct(r_sil_bot)

    # Use median radii for the wedge patches
    layers = [
        ('Ice Ih',     R_body_km, p_ih[1]),
        ('Ocean',      p_ih[1],   p_oc[1]),
        ('Ice II',     p_oc[1],   p_ii[1]),
        ('Ice III',    p_ii[1],   p_iii[1]),
        ('Ice V',      p_iii[1],  p_v[1]),
        ('Ice VI',     p_v[1],    p_vi[1]),
        ('Silicate',   p_vi[1],   p_sil[1]),
        ('Dense core', p_sil[1],  0.0),
    ]

    wedge_colors = _wedge_color_map()

    fig, ax = plt.subplots(figsize=(6, 8))
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.05, 1.05)
    ax.set_aspect('equal')
    ax.axis('off')

    cx, cy = 0.5, 0.0  # wedge centre at bottom

    for name, r_top, r_bot in layers:
        thickness_km = r_top - r_bot
        if thickness_km < 0.5:
            continue
        r_outer_norm = r_top / R_body_km
        width_norm   = thickness_km / R_body_km
        wedge = MplWedge(
            (cx, cy), r_outer_norm, ANG1, ANG2,
            width=width_norm,
            fc=wedge_colors[name], ec='#333333', lw=0.8,
        )
        ax.add_patch(wedge)

    # Labels on the left side with leader arrows
    tick_angle_rad = np.radians(ANG1 - 3)
    label_entries = []
    for name, r_top, r_bot in layers:
        if r_top - r_bot < 5:
            continue
        r_mid  = (r_top + r_bot) / 2
        r_norm = r_mid / R_body_km
        x = cx + r_norm * np.cos(tick_angle_rad)
        y = cy + r_norm * np.sin(tick_angle_rad)
        label_entries.append((x, y, name, r_top - r_bot))

    # Space labels vertically (min 0.06 apart in normalised coords)
    label_y_positions = [e[1] for e in label_entries]
    spaced_y: List[float] = []
    for i, y in enumerate(sorted(label_y_positions)):
        if i > 0 and y - spaced_y[-1] < 0.06:
            y = spaced_y[-1] + 0.06
        spaced_y.append(y)
    sorted_idx = sorted(range(len(label_y_positions)), key=lambda i: label_y_positions[i])
    final_y = [0.0] * len(label_entries)
    for rank, orig_i in enumerate(sorted_idx):
        final_y[orig_i] = spaced_y[rank]

    for i, (x, y, name, thick) in enumerate(label_entries):
        label_x = -0.08
        label_y = final_y[i]
        ax.annotate(
            f'{name} ({thick:.0f} km)',
            xy=(x, y), xytext=(label_x, label_y),
            fontsize=8, ha='right', va='center', color='#333333',
            arrowprops=dict(arrowstyle='-', color='#666666', lw=0.6),
        )

    ax.text(
        cx, R_body_km / R_body_km + 0.02 + cy,
        f'R = {R_body_km:.0f} km',
        ha='center', va='bottom', fontsize=9,
    )

    ocean_thick  = p_ih[1] - p_oc[1]
    ice_ih_thick = R_body_km - p_ih[1]
    hp_thick     = p_oc[1] - p_vi[1]  # All HP ices (II, III, V, VI)
    ax.set_title(
        f'{body_name} Interior Structure (Posterior Median)\n'
        f'Ice Ih: {ice_ih_thick:.0f} km | Ocean: {ocean_thick:.0f} km | '
        f'HP ices: {hp_thick:.0f} km',
        fontsize=11, pad=10,
    )

    handles = [
        mpatches.Patch(color=c, label=l)
        for l, c in wedge_colors.items()
        if any(l == name and r_top - r_bot > 0.5 for name, r_top, r_bot in layers)
    ]
    ax.legend(handles=handles, loc='lower left', fontsize=8, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    log.info(f'Saved {output_path}')
    plt.close(fig)


# ===========================================================================
# 10. PP-canonical wedge helpers
# ===========================================================================

def _run_pp_with_overrides(
    planet_template_module: str,
    overrides: Dict[str, Any],
    ocean_overrides: Optional[Dict[str, Any]] = None,
) -> Tuple[Any, Any]:
    """Deep-copy Planet from planet_template_module, apply scalar overrides,
    run PlanetProfile with CALC_NEW=True/NO_SAVEFILE=True/SKIP_PLOTS=True.
    Returns (Planet, Params).

    Parameters
    ----------
    planet_template_module
        Importable dotted module path that exposes a ``Planet`` object at
        module top-level (e.g. ``'PlanetProfile.Test.PPTest48'``).
    overrides
        Dict mapping dotted attribute paths (relative to ``Planet``) to values.
        E.g. ``{'Bulk.Tb_K': 250.5, 'Ocean.wOcean_ppt': 10.0}``.
    ocean_overrides
        Dict mapping ``Planet.Ocean.<key>`` attribute names to values.
        Applied after ``overrides``.  Use for composition switches:
        ``{'comp': 'NaCl', 'wOcean_ppt': 100.0}``.
    """
    import sys as _sys

    from PlanetProfile.GetConfig import Params as configParams
    from PlanetProfile.Main import PlanetProfile as RunPP

    # Reload module to pick up any in-process edits, then deep-copy Planet.
    if planet_template_module in _sys.modules:
        importlib.reload(_sys.modules[planet_template_module])
    else:
        importlib.import_module(planet_template_module)
    mod = _sys.modules[planet_template_module]
    Planet = deepcopy(mod.Planet)

    # Apply dotted-path overrides (e.g. 'Bulk.Tb_K' -> Planet.Bulk.Tb_K)
    for dotted_key, value in overrides.items():
        parts = dotted_key.split('.')
        obj = Planet
        for part in parts[:-1]:
            obj = getattr(obj, part)
        setattr(obj, parts[-1], value)

    # Apply ocean_overrides to Planet.Ocean
    if ocean_overrides:
        for key, value in ocean_overrides.items():
            if not hasattr(Planet.Ocean, key):
                log.warning(
                    f'_run_pp_with_overrides: Planet.Ocean has no attribute {key!r}; '
                    'setting it anyway.'
                )
            setattr(Planet.Ocean, key, value)

    tb_val = overrides.get('Bulk.Tb_K', getattr(getattr(Planet, 'Bulk', None), 'Tb_K', None))
    log.info(
        f'_run_pp_with_overrides: template={planet_template_module!r}  Tb_K={tb_val}'
    )

    configParams.CALC_NEW = True
    configParams.NO_SAVEFILE = True
    configParams.SKIP_PLOTS = True

    Planet, Params = RunPP(Planet, configParams)
    return Planet, Params


def _validate_planet_against_cache(
    Planet: Any,
    cache_pt: Dict[str, Any],
    strict: bool = False,
) -> Dict[str, float]:
    """Compare key structural arrays between a freshly-run Planet and a cached grid point.

    Both the fresh Planet and the cache point expose the following arrays
    (same keys used in ``cache_builder.build_single_structure``):

    * ``r_m``   — radius grid (m)
    * ``rho``   — density (kg/m^3)
    * ``T_K``   — temperature (K)
    * ``P_MPa`` — pressure (MPa)
    * ``mu_Pa`` — shear modulus (Pa)
    * ``K_Pa``  — bulk modulus (Pa)

    The fresh Planet arrays are extracted from
    ``Planet.Gravity.ALMAModel['model']`` using the same column-index
    logic as ``cache_builder.build_single_structure``.  If the ALMA model
    is unavailable, the function falls back to ``Planet.r_m`` / ``Planet.T_K``
    / ``Planet.P_MPa`` for the scalar profiles and skips elastic moduli.

    Arrays of different lengths are reconciled by interpolating both onto
    a common radius grid before comparison.

    Parameters
    ----------
    Planet
        A ``Planet`` object returned by ``PlanetProfile.Main.PlanetProfile``.
    cache_pt
        A single structure dict from the grid cache (v2.0 or v2.1 format).
    strict
        If ``True``, raise ``ValueError`` on the first tolerance violation.

    Returns
    -------
    Dict mapping field name to its maximum relative (or absolute) error.
    """
    # Tolerances: (rtol, atol).  None means "not used for that mode".
    _TOLERANCES: Dict[str, Tuple[Optional[float], Optional[float]]] = {
        'r_m':   (1e-4, 1.0),
        'rho':   (1e-3, None),
        'T_K':   (None, 0.5),
        'P_MPa': (1e-2, None),
        'mu_Pa': (5e-2, None),
        'K_Pa':  (5e-2, None),
    }

    # --- Extract fresh arrays from Planet ---
    fresh: Dict[str, Optional[np.ndarray]] = {k: None for k in _TOLERANCES}
    r_m_fresh: Optional[np.ndarray] = None
    r_m_primary: Optional[np.ndarray] = None
    sort_idx: Optional[np.ndarray] = None

    try:
        model = Planet.Gravity.ALMAModel['model']
        cols = Planet.Gravity.columns
        r_col   = cols.index('r')
        rho_col = cols.index('rho')
        VP_col  = cols.index('VP')
        GS_col  = cols.index('GS')

        r_m_fresh   = model[:, r_col].astype(np.float64)
        rho_fresh   = model[:, rho_col].astype(np.float64)
        mu_fresh    = model[:, GS_col].astype(np.float64)
        VP_fresh    = model[:, VP_col].astype(np.float64)
        K_fresh_raw = rho_fresh * VP_fresh ** 2 - (4.0 / 3.0) * mu_fresh
        K_fresh     = np.maximum(K_fresh_raw, 1e6)

        fresh['r_m']   = r_m_fresh
        fresh['rho']   = rho_fresh
        fresh['mu_Pa'] = mu_fresh
        fresh['K_Pa']  = K_fresh
    except Exception as exc:
        log.warning(
            f'_validate_planet_against_cache: ALMA model extraction failed '
            f'({exc}); will skip elastic-moduli comparison.'
        )

    # Temperature and pressure from the primary PP profile
    try:
        T_K_arr = np.asarray(Planet.T_K, dtype=np.float64)
        r_m_primary = np.asarray(Planet.r_m[: T_K_arr.size], dtype=np.float64)
        sort_idx    = np.argsort(r_m_primary)
        r_ref       = r_m_fresh if r_m_fresh is not None else r_m_primary[sort_idx]
        fresh['T_K'] = np.interp(r_ref, r_m_primary[sort_idx], T_K_arr[sort_idx])
        if fresh['r_m'] is None:
            fresh['r_m'] = r_ref
    except Exception as exc:
        log.warning(f'_validate_planet_against_cache: T_K extraction failed ({exc}).')

    try:
        P_arr   = np.asarray(Planet.P_MPa, dtype=np.float64)
        r_ref2  = r_m_fresh if r_m_fresh is not None else (
            r_m_primary[sort_idx] if r_m_primary is not None and sort_idx is not None
            else None
        )
        if r_ref2 is not None and r_m_primary is not None and sort_idx is not None:
            P_trimmed = P_arr[: r_m_primary.size]
            fresh['P_MPa'] = np.interp(r_ref2, r_m_primary[sort_idx], P_trimmed[sort_idx])
    except Exception as exc:
        log.warning(f'_validate_planet_against_cache: P_MPa extraction failed ({exc}).')

    # --- Compare field by field ---
    errors: Dict[str, float] = {}
    violations: List[str] = []

    for field_name, (rtol, atol) in _TOLERANCES.items():
        arr_fresh = fresh.get(field_name)
        arr_cache = cache_pt.get(field_name)

        if arr_fresh is None or arr_cache is None:
            log.debug(f'_validate_planet_against_cache: skipping {field_name} (unavailable).')
            continue

        arr_fresh = np.asarray(arr_fresh, dtype=np.float64)
        arr_cache = np.asarray(arr_cache, dtype=np.float64)

        # Interpolate to a common radius grid if lengths differ
        if arr_fresh.size != arr_cache.size:
            r_fresh_ref = fresh.get('r_m')
            r_cache_ref = cache_pt.get('r_m')
            if r_fresh_ref is not None and r_cache_ref is not None and field_name != 'r_m':
                r_f = np.asarray(r_fresh_ref, dtype=np.float64)
                r_c = np.asarray(r_cache_ref, dtype=np.float64)
                r_common = np.linspace(
                    max(float(r_f[0]), float(r_c[0])),
                    min(float(r_f[-1]), float(r_c[-1])),
                    min(arr_fresh.size, arr_cache.size),
                )
                arr_fresh = np.interp(r_common, r_f, arr_fresh)
                arr_cache = np.interp(r_common, r_c, arr_cache)
            else:
                # Last resort: trim to shorter length
                n = min(arr_fresh.size, arr_cache.size)
                arr_fresh = arr_fresh[:n]
                arr_cache = arr_cache[:n]

        diff = np.abs(arr_fresh - arr_cache)

        if rtol is not None:
            scale = np.maximum(np.abs(arr_cache), 1e-30)
            err_val = float(np.nanmax(diff / scale))
            exceeded = err_val > rtol
            tol_desc = f'rtol={rtol}'
        else:
            err_val = float(np.nanmax(diff))
            exceeded = err_val > (atol if atol is not None else 0.0)
            tol_desc = f'atol={atol}'

        errors[field_name] = err_val
        if exceeded:
            msg = (
                f'_validate_planet_against_cache: {field_name} exceeds tolerance '
                f'({tol_desc}): max_err={err_val:.4g}'
            )
            log.warning(msg)
            violations.append(msg)

    if strict and violations:
        raise ValueError(
            'Cache validation failed (strict=True). Violations:\n'
            + '\n'.join(violations)
        )

    return errors


def plot_structure_wedge_pp(
    result: Any,
    grid_cache: Dict[str, Any],
    output_path: str,
    *,
    planet_template_module: Optional[str] = None,
    param_overrides: Optional[Dict[str, str]] = None,
    use: str = 'median',
    sample_index: Optional[int] = None,
    strict_validate: bool = False,
    fig_format: Optional[str] = None,
) -> None:
    """Re-run PlanetProfile at the posterior point and produce the canonical PP wedge.

    Parameters
    ----------
    result
        An ``InferenceResult`` object with attributes ``samples``,
        ``param_names``, ``log_likelihoods``, and ``config``.
    grid_cache
        Structure cache dict.  Supports both v2.1 list format
        (``{'Tb_K_grid': [...], 'structures': [...]}`` ) and v2.0 dict
        format (``{'grid_cache': {(Tb, D): {...}}, 'grid_Tb_values': [...]}``).
    output_path
        Absolute path for the output figure (including extension).
    planet_template_module
        Importable module path (e.g. ``'PlanetProfile.Test.PPTest48'``).
        Defaults to ``result.config.planet_template_module``.
        Raises ``ValueError`` if still ``None`` after resolution.
    param_overrides
        Extra mappings from MCMC param name to PP dotted-attribute path,
        extending the built-in defaults
        ``{'Tb_K': 'Bulk.Tb_K', 'wOcean_ppt': 'Ocean.wOcean_ppt'}``.
    use
        Which posterior point to use: ``'median'`` (default),
        ``'best_fit'`` (highest log-likelihood), or ``'sample'``
        (requires ``sample_index``).
    sample_index
        Row index into ``result.samples`` used when ``use='sample'``.
    strict_validate
        If ``True``, raise on any cache-vs-fresh-Planet tolerance violation.
    fig_format
        Output figure format (one of ``'png'``, ``'pdf'``, ``'eps'``, ``'svg'``,
        ``'jpg'``, ``'jpeg'``, ``'tif'``, ``'tiff'``).  When ``None`` (default)
        the function honors whatever ``FigMisc.figFormat`` is currently set to
        on the PP-wide config.  When given, the function temporarily sets
        ``FigMisc.figFormat`` (and ``FigMisc.xtn``) to this value for the
        duration of the call and restores the prior values in a ``finally``
        block.  Either way, the extension on ``output_path`` is rewritten to
        match the resolved format so the on-disk filename always agrees with
        the byte stream PP writes.
    """
    from PlanetProfile.Plotting.ProfilePlots import PlotWedge

    # 1. Resolve planet_template_module
    if planet_template_module is None:
        planet_template_module = getattr(result.config, 'planet_template_module', None)
    if planet_template_module is None:
        raise ValueError(
            'plot_structure_wedge_pp: planet_template_module was not provided and '
            'result.config.planet_template_module is None. '
            'Pass it explicitly or set it on result.config.'
        )

    # 2. Choose posterior point
    samples = np.asarray(result.samples)
    param_names = list(result.param_names)

    if use == 'median':
        theta = np.median(samples, axis=0)
    elif use == 'best_fit':
        lls = np.asarray(result.log_likelihoods, dtype=float)
        if np.all(~np.isfinite(lls)):
            log.warning(
                'plot_structure_wedge_pp: all log_likelihoods are non-finite; '
                'falling back to median.'
            )
            theta = np.median(samples, axis=0)
        else:
            theta = samples[int(np.nanargmax(lls))]
    elif use == 'sample':
        if sample_index is None:
            raise ValueError("use='sample' requires sample_index to be set.")
        theta = samples[sample_index]
    else:
        raise ValueError(f"use must be 'median', 'best_fit', or 'sample'; got {use!r}.")

    theta_dict: Dict[str, float] = dict(zip(param_names, theta.tolist()))

    # 3. Build PP override dict
    _DEFAULT_PARAM_MAP: Dict[str, str] = {
        'Tb_K':       'Bulk.Tb_K',
        'wOcean_ppt': 'Ocean.wOcean_ppt',
    }
    param_map = dict(_DEFAULT_PARAM_MAP)
    if param_overrides:
        param_map.update(param_overrides)

    pp_overrides: Dict[str, Any] = {}

    # From sampled params
    for mcmc_key, pp_path in param_map.items():
        if mcmc_key in theta_dict:
            pp_overrides[pp_path] = theta_dict[mcmc_key]

    # From fixed_params (e.g. Tb_K fixed in no-ocean runs)
    fixed: Dict[str, Any] = getattr(result.config, 'fixed_params', {}) or {}
    for mcmc_key, pp_path in param_map.items():
        if mcmc_key in fixed and pp_path not in pp_overrides:
            pp_overrides[pp_path] = fixed[mcmc_key]

    log.info(f'plot_structure_wedge_pp: PP overrides = {pp_overrides}  (use={use!r})')

    # 4. Get ocean_overrides from config
    ocean_overrides: Dict[str, Any] = getattr(result.config, 'ocean_overrides', {}) or {}

    # 5. Run PlanetProfile
    Planet, Params = _run_pp_with_overrides(
        planet_template_module, pp_overrides, ocean_overrides
    )

    # 5b. Sanitize wedge-rendering inputs that PP may leave as NaN/None.
    #
    # ProfilePlots.PlotWedge (core PP, not modified here) uses the following on
    # the silicate shell:
    #   silOuter radius = Planet.Sil.Rmean_m / 1e3 / rMax_km
    #   silOuter width  = (Planet.Sil.Rmean_m - Planet.Core.Rmean_m) / rMax_km
    #   dzSilCond_km    = (Planet.Sil.Rmean_m - Planet.Core.Rmean_m) / 1e3
    #                     - Planet.dzSilPorous_km
    #   if dzSilCond_km > 0: draw conductive gradient patches
    #
    # When inference re-runs PP via _run_pp_with_overrides at a posterior
    # median (especially with CONSTANT_INNER_DENSITY=True or POROUS_ROCK=True
    # variants), Planet.dzSilPorous_km can be NaN and Planet.Core.Rmean_m can
    # be None. NaN propagates through dzSilCond_km, making the conductive-
    # gradient guard fail (NaN > 0 is False) so no silicate patches are drawn.
    # Per-patch widths in the porous block also become NaN, silently rendering
    # nothing. This patch leaves PP physics untouched and only repairs the
    # plotting attributes; nothing scientific is altered.
    def _is_missing(v):
        """True if v is None, NaN, or non-finite; False otherwise."""
        if v is None:
            return True
        try:
            return not np.isfinite(v)
        except (TypeError, ValueError):
            return True

    try:
        # mantleEOS sanitization: PlotWedge (ProfilePlots.py) does
        #   ``if 'Comet' in Planet.Sil.mantleEOS:``
        # which raises TypeError when mantleEOS is None, and produces a
        # misleading label ('no chondrite') when PP has overwritten it with
        # the sentinel value 'none' (set by SetupInit when
        # CONSTANT_INNER_DENSITY=True).  Restore the original template value
        # from the already-imported module so the legend shows the real EOS
        # name (e.g. 'CV3hy1wt_678_1.tab' → 'CV chondrite').
        _cur_meos = getattr(getattr(Planet, 'Sil', None), 'mantleEOS', None)
        if not _cur_meos or _cur_meos == 'none':
            import sys as _sys_meos
            _tmpl_mod = _sys_meos.modules.get(planet_template_module)
            _tmpl_meos = (
                getattr(getattr(getattr(_tmpl_mod, 'Planet', None), 'Sil', None),
                        'mantleEOS', None)
                if _tmpl_mod is not None else None
            )
            if _tmpl_meos and _tmpl_meos != 'none':
                Planet.Sil.mantleEOS = _tmpl_meos
                log.debug(
                    f'plot_structure_wedge_pp: restored mantleEOS={_tmpl_meos!r} '
                    'from template module (PP had overwritten it with '
                    f'{_cur_meos!r}).'
                )
            else:
                # Last resort: keep 'none' if we truly have no better value,
                # but replace None with the sentinel so PlotWedge does not crash.
                if _cur_meos is None:
                    Planet.Sil.mantleEOS = 'none'
                    log.warning(
                        'plot_structure_wedge_pp: mantleEOS is None and could not '
                        'be recovered from template module; set to "none". '
                        'The wedge legend will show "no chondrite" — check '
                        f'planet_template_module={planet_template_module!r}.'
                    )

        if getattr(Planet.Core, 'Rmean_m', None) is None or not np.isfinite(
            Planet.Core.Rmean_m
        ):
            Planet.Core.Rmean_m = 0.0
        if getattr(Planet.Sil, 'Rmean_m', None) is None or not np.isfinite(
            Planet.Sil.Rmean_m
        ) or Planet.Sil.Rmean_m == 0:
            sil_mask = np.asarray(Planet.phase) == 50  # 50 = silicate in PP
            if sil_mask.any():
                Planet.Sil.Rmean_m = float(np.max(np.asarray(Planet.r_m)[sil_mask]))

        # Ice / clathrate layer geometry: PlotWedge guards each HP-ice patch
        # with ``if Planet.dzIce*_km > 0`` and then computes the wedge radius
        # from ``Planet.zIce*_m`` (or ``zClath_km``).  When PP runs at a
        # posterior median via _run_pp_with_overrides, those scalar layer
        # summaries are sometimes left as NaN even though Planet.r_m and
        # Planet.phase fully describe the radial structure.  NaN > 0 is False,
        # so the layer is silently suppressed and the wedge ends up empty
        # outside the silicate shell.  Reconstruct any missing dz/z values
        # from (r_m, phase) so that the layers render.  PP's own values are
        # left untouched whenever they are finite.
        #
        # PP phase codes (PlanetProfile.Utilities.defineStructs.Constants):
        #   1=IceIh, 2=IceII, 3=IceIII, 5=IceV, 6=IceVI,
        #   30=clathrate (phaseClath), 50=silicate, 100=Fe.
        _phase_to_attrs = {
            1:  ('dzIceI_km',   'zIceI_m'),
            2:  ('dzIceII_km',  'zIceII_m'),
            3:  ('dzIceIII_km', 'zIceIII_m'),
            5:  ('dzIceV_km',   'zIceV_m'),
            6:  ('dzIceVI_km',  'zIceVI_m'),
            30: ('dzClath_km',  'zClath_km'),  # zClath stored in km, not m
        }

        r_m_arr = np.asarray(getattr(Planet, 'r_m', []), dtype=float)
        phase_arr_raw = getattr(Planet, 'phase', None)
        phase_arr = (
            np.asarray(phase_arr_raw, dtype=int)
            if phase_arr_raw is not None else np.array([], dtype=int)
        )
        bulk_R_m = getattr(Planet.Bulk, 'R_m', None)
        if bulk_R_m is not None and np.isfinite(bulk_R_m):
            surface_r_m = float(bulk_R_m)
        elif r_m_arr.size:
            surface_r_m = float(r_m_arr[0])
        else:
            surface_r_m = 0.0

        # PP uses a cell-centered radial grid: ``Planet.r_m`` holds N+1
        # boundary radii (descending from surface to centre) and
        # ``Planet.phase`` holds N cell phases.  Cell ``i`` spans the radii
        # ``r_m[i]`` (top) to ``r_m[i+1]`` (bottom) and carries phase
        # ``phase[i]``.  Some PP runs return r_m with the same length as
        # phase (boundaries-only convention); handle both.
        n_phase = phase_arr.size
        n_r = r_m_arr.size
        if n_phase and n_r and n_r in (n_phase, n_phase + 1):
            for ph_code, (dz_attr, z_attr) in _phase_to_attrs.items():
                mask = phase_arr == ph_code
                dz_curr = getattr(Planet, dz_attr, None)
                z_curr  = getattr(Planet, z_attr,  None)
                if mask.any():
                    idx = np.where(mask)[0]
                    first = int(idx[0])
                    last = int(idx[-1])
                    # Top of this phase is the upper boundary of the first
                    # cell of this phase = r_m[first].
                    r_top = float(r_m_arr[first])
                    # Bottom is the lower boundary of the last cell of this
                    # phase = r_m[last + 1] when boundaries-array is N+1
                    # long, else r_m[last].
                    if n_r == n_phase + 1:
                        r_bot = float(r_m_arr[last + 1])
                    elif last + 1 < n_r:
                        r_bot = float(r_m_arr[last + 1])
                    else:
                        r_bot = float(r_m_arr[last])
                    z_top_m = surface_r_m - r_top
                    dz_km = max(r_top - r_bot, 0.0) / 1e3
                    if _is_missing(dz_curr):
                        setattr(Planet, dz_attr, dz_km)
                    if _is_missing(z_curr):
                        if z_attr.endswith('_km'):
                            setattr(Planet, z_attr, z_top_m / 1e3)
                        else:
                            setattr(Planet, z_attr, z_top_m)
                else:
                    # Phase absent from the radial structure: 0 is the
                    # PlotWedge "draw nothing" value (guards use ``> 0``).
                    if _is_missing(dz_curr):
                        setattr(Planet, dz_attr, 0.0)
                    if _is_missing(z_curr):
                        setattr(Planet, z_attr, 0.0)

        # Layer attributes that PlotWedge guards with ``> 0``.  When PP did
        # not exercise the corresponding mode (POROUS_ROCK, FeS layer, surface
        # HP undifferentiated ices, wet HP ices), the attribute is left NaN.
        # Coerce to 0 so the guard correctly suppresses that wedge patch
        # instead of corrupting the figure with NaN geometry.
        for _dz_attr in (
            'dzSilPorous_km', 'dzFeS_km',
            'dzIceIIIund_km', 'dzIceVund_km',
            'dzWetHPs_km',
        ):
            _v = getattr(Planet, _dz_attr, None)
            if _is_missing(_v):
                setattr(Planet, _dz_attr, 0.0)
        for _z_attr in ('zIceIIIund_m', 'zIceVund_m'):
            _v = getattr(Planet, _z_attr, None)
            if _is_missing(_v):
                setattr(Planet, _z_attr, 0.0)

        # Convective-thermal attributes used by ice Ih shell rendering
        # (lines 1199-1206 of ProfilePlots.PlotWedge).  NaN here causes the
        # surface ice shell to render with NaN width and become invisible.
        # 0 is safe: no convective layer is drawn, the conductive lid covers
        # the full ice I thickness via dzIceI_km.
        for _conv_attr in ('eLid_m', 'Dconv_m', 'deltaTBL_m'):
            _v = getattr(Planet, _conv_attr, None)
            if _is_missing(_v):
                setattr(Planet, _conv_attr, 0.0)

        # If eLid_m is set but Dconv_m + deltaTBL_m was zeroed from NaN,
        # and the total ice I thickness exceeds the conductive lid, there
        # is a convective zone that PlotWedge would leave unrendered
        # (the guard at ProfilePlots line 1203 / 1171 is
        # ``(Planet.Dconv_m + Planet.deltaTBL_m) > 0``).  Reconstruct the
        # missing convective thickness from the difference between the
        # full ice I shell (dzIceI_km, already reconstructed from r_m /
        # phase above) and the conductive stagnant lid.
        #
        # Safety:
        #   - If PP itself populated Dconv_m > 0, this branch is skipped
        #     (guard ``_dconv + _dtbl <= 0``) and PP's value is preserved.
        #   - If the run is fully conductive (no convection in PP),
        #     eLid_m == dzIceI_km*1e3 so _conv_thickness_m == 0 and
        #     nothing changes.
        #   - deltaTBL_m is intentionally left at 0; PlotWedge only uses
        #     the sum (Dconv_m + deltaTBL_m) for the convective patch
        #     width, so attributing the recovered thickness to Dconv_m
        #     alone yields the same wedge geometry without inventing a
        #     thermal-boundary-layer thickness PP did not compute.
        _eLid_m = float(getattr(Planet, 'eLid_m', 0.0) or 0.0)
        _dzIceI_m = float(getattr(Planet, 'dzIceI_km', 0.0) or 0.0) * 1e3
        _dconv = float(getattr(Planet, 'Dconv_m', 0.0) or 0.0)
        _dtbl = float(getattr(Planet, 'deltaTBL_m', 0.0) or 0.0)
        if _dzIceI_m > 0 and (_dconv + _dtbl) <= 0:
            _conv_thickness_m = max(0.0, _dzIceI_m - _eLid_m)
            if _conv_thickness_m > 0:
                Planet.Dconv_m = _conv_thickness_m
                # deltaTBL_m stays 0; PlotWedge uses the sum
                # Dconv_m + deltaTBL_m, so this yields a single
                # iceIconv-colored convective wedge of the correct
                # thickness without inventing a TBL value.

        # Ocean attributes.  In no-ocean runs PP may leave D_km / zb_km as
        # NaN; PlotWedge guards ocean rendering with ``Planet.D_km > 0``.
        for _ocean_attr in ('D_km', 'zb_km'):
            _v = getattr(Planet, _ocean_attr, None)
            if _is_missing(_v):
                setattr(Planet, _ocean_attr, 0.0)
    except Exception as exc:  # noqa: BLE001
        log.warning(
            f'plot_structure_wedge_pp: wedge attribute sanitization failed: {exc!r}'
        )

    # 6. Find nearest cache point for the resolved Tb_K
    resolved_Tb: Optional[float] = None
    if 'Bulk.Tb_K' in pp_overrides:
        resolved_Tb = float(pp_overrides['Bulk.Tb_K'])
    elif 'Tb_K' in fixed:
        resolved_Tb = float(fixed['Tb_K'])

    cache_pt: Optional[Dict[str, Any]] = None

    # v2.1 list format: {'Tb_K_grid': [...], 'structures': [...]}
    if 'Tb_K_grid' in grid_cache and 'structures' in grid_cache:
        tb_grid = np.asarray(grid_cache['Tb_K_grid'], dtype=np.float64)
        if resolved_Tb is not None:
            idx = int(np.argmin(np.abs(tb_grid - resolved_Tb)))
        else:
            idx = 0
        structs = grid_cache['structures']
        if 0 <= idx < len(structs):
            cache_pt = structs[idx]

    # v2.0 dict format: {'grid_cache': {(Tb, D): {...}}, 'grid_Tb_values': [...]}
    elif 'grid_cache' in grid_cache:
        inner: Dict[Any, Dict[str, Any]] = grid_cache['grid_cache']
        if resolved_Tb is not None and inner:
            tb_arr = np.asarray([k[0] for k in inner.keys()], dtype=np.float64)
            nearest_tb = float(tb_arr[np.argmin(np.abs(tb_arr - resolved_Tb))])
            for key, val in inner.items():
                if abs(key[0] - nearest_tb) < 0.01:
                    cache_pt = val
                    break
        elif inner:
            cache_pt = next(iter(inner.values()))

    if cache_pt is None:
        log.warning(
            'plot_structure_wedge_pp: could not locate a cache point for '
            f'Tb_K={resolved_Tb}; skipping validation.'
        )
    else:
        # 7. Validate fresh Planet against the cache point
        _validate_planet_against_cache(Planet, cache_pt, strict=strict_validate)

    # 8. Resolve fig_format and rewrite output_path extension to match.
    #
    # PlanetProfile's PlotWedge delegates to
    # ``fig.savefig(path, format=FigMisc.figFormat, ...)``, so the byte stream
    # is governed by ``FigMisc.figFormat`` regardless of the extension on the
    # filename.  Resolve the format first, then force the on-disk filename to
    # agree with it — the extension on the caller-supplied ``output_path`` is
    # never authoritative.
    from PlanetProfile.GetConfig import FigMisc as _FigMisc
    _supported = {'png', 'pdf', 'eps', 'svg', 'jpg', 'jpeg', 'tif', 'tiff'}
    _prev_figFormat = getattr(_FigMisc, 'figFormat', None)
    _prev_xtn = getattr(_FigMisc, 'xtn', None)

    if fig_format is None:
        _resolved_format = (_prev_figFormat or 'pdf').lower()
    else:
        _resolved_format = str(fig_format).lower().lstrip('.')

    if _resolved_format not in _supported:
        raise ValueError(
            f"plot_structure_wedge_pp: fig_format={_resolved_format!r} is not "
            f"supported.  Choose one of {sorted(_supported)}."
        )

    # Temporarily set the PP-wide format so PlotWedge writes the right bytes.
    _FigMisc.figFormat = _resolved_format
    _FigMisc.xtn = '.' + _resolved_format

    # Rewrite output_path's extension to match the resolved format so the
    # filename always agrees with the byte stream PlotWedge produces.
    _path_root, _ = os.path.splitext(output_path)
    output_path = _path_root + '.' + _resolved_format

    # 8b. Configure output path and flags on the Params object.
    Params.FigureFiles.vwedg = output_path
    Params.ALL_ONE_BODY = True
    # PlotWedge does not gate on SKIP_PLOTS, but make intent explicit for any
    # future call sites that might.
    Params.SKIP_PLOTS = False

    # 8c. Self-contained guards around PlotWedge so callers do not need to
    # monkey-patch matplotlib at the driver level:
    #   - ``Figure.tight_layout`` raises
    #     ``ValueError("need at least one array to concatenate")`` when any
    #     ``Wedge`` patch has zero radial extent (e.g. ocean depth = 0,
    #     no HP ice).  Catch that specific ValueError and continue; the
    #     figure is still valid.
    #   - When ``tight_layout`` is skipped (or simply produces tight axes
    #     against the figure bbox), labels at the figure edge are clipped on
    #     save.  Inject ``bbox_inches='tight'`` and ``pad_inches`` into the
    #     ``savefig`` call so labels are never clipped, regardless of which
    #     path tight_layout took.
    import matplotlib.figure as _mpl_fig
    _orig_tight_layout = _mpl_fig.Figure.tight_layout
    _orig_savefig = _mpl_fig.Figure.savefig

    def _safe_tight_layout(self, *args, **kwargs):
        try:
            return _orig_tight_layout(self, *args, **kwargs)
        except ValueError as e:
            if 'concatenate' in str(e):
                log.warning(
                    'plot_structure_wedge_pp: tight_layout skipped due to '
                    f'empty Wedge patch ({e}); relying on bbox_inches="tight".'
                )
                return None
            raise

    def _safe_savefig(self, fname, *args, **kwargs):
        # Try with bbox_inches='tight' first so labels at the figure edge
        # are never clipped.  If matplotlib's tight-bbox computation raises
        # the empty-Wedge ``concatenate`` ValueError (same root cause as
        # the tight_layout crash), strip the tight bbox kwargs, widen the
        # figure margins instead, and retry so we still keep labels visible
        # on figures that contain a zero-extent ``Wedge`` patch.
        kwargs.setdefault('bbox_inches', 'tight')
        kwargs.setdefault('pad_inches', 0.1)
        try:
            return _orig_savefig(self, fname, *args, **kwargs)
        except ValueError as e:
            if 'concatenate' not in str(e):
                raise
            log.warning(
                'plot_structure_wedge_pp: bbox_inches="tight" failed due to '
                f'empty Wedge patch ({e}); retrying with subplots_adjust '
                'margins so labels stay visible.'
            )
            kwargs.pop('bbox_inches', None)
            kwargs.pop('pad_inches', None)
            try:
                self.subplots_adjust(left=0.12, right=0.95, top=0.93, bottom=0.10)
            except Exception:  # pragma: no cover - defensive
                pass
            return _orig_savefig(self, fname, *args, **kwargs)

    _mpl_fig.Figure.tight_layout = _safe_tight_layout
    _mpl_fig.Figure.savefig = _safe_savefig

    try:
        # 9. Produce the PP wedge figure
        PlotWedge([Planet], Params)
        log.info(f'plot_structure_wedge_pp: saved PP wedge to {output_path}')
    finally:
        # Restore matplotlib and FigMisc state regardless of plot outcome.
        _mpl_fig.Figure.tight_layout = _orig_tight_layout
        _mpl_fig.Figure.savefig = _orig_savefig
        if _prev_figFormat is not None:
            _FigMisc.figFormat = _prev_figFormat
        if _prev_xtn is not None:
            _FigMisc.xtn = _prev_xtn
