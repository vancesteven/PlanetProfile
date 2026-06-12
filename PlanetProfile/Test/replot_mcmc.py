"""
Regenerate MCMC plots from saved pickle data.

Produces Petricca et al. (2025)-style filled KDE corner plots,
k2 scatter plots, and heating distribution plots for:
  - Andrade no-ocean (PPTest41)
  - Maxwell ocean (PPTest42)
  - Andrade Arrhenius no-ocean (PPTest43)
  - Maxwell Arrhenius ocean (PPTest44)

Usage:
  python PlanetProfile/Test/replot_mcmc.py
"""
import numpy as np
import os
import pickle
import logging

logging.basicConfig(level=logging.INFO, format='%(name)s - %(message)s')
log = logging.getLogger('replot_mcmc')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse
import seaborn as sns

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_BASE = os.path.join(_THIS_DIR, 'mcmc_results')
DIR41 = os.path.join(_BASE, 'Titan', 'Test41_andrade_no_ocean')
DIR42 = os.path.join(_BASE, 'Titan', 'Test42_maxwell_ocean')
DIR43 = os.path.join(_BASE, 'Titan', 'Test43_andrade_arrhenius_no_ocean')
DIR44 = os.path.join(_BASE, 'Titan', 'Test44_maxwell_arrhenius_ocean')

# Observational constraints
RE_K2_OBS, RE_K2_ERR = 0.608, 0.048
IM_K2_OBS, IM_K2_ERR = 0.135, 0.035

# Consistent phase color scheme used across all plots (matches mcmc_plots.py)
PHASE_COLORS = {
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
    'Fe':    '#8B5A2B',  # Alias for Core
}
PHASE_LABELS = {
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
    'Fe':    'Iron Core',
}
# Canonical radial order (center to surface)
PHASE_ORDER_CANONICAL = ['Core', 'Fe', 'Sil', 'Rock', 'VI', 'V', 'III', 'II', '0', 'Ih', 'Clath']


# ============================================================
# Colored KDE corner plot
# ============================================================
def kde_corner(samples, labels, title, output_path, cmap='Blues',
               color='steelblue', truths=None):
    """Petricca et al.-style filled KDE corner plot.

    Args:
        samples: (n_samples, n_dim) array
        labels: list of axis labels (LaTeX OK)
        title: figure title
        output_path: path for saved PNG
        cmap: colormap for 2D KDE panels
        color: color for 1D marginal fills
        truths: optional list of true values to mark
    """
    n_dim = samples.shape[1]
    fig = plt.figure(figsize=(2.5 * n_dim, 2.5 * n_dim))
    gs = gridspec.GridSpec(n_dim, n_dim, hspace=0.08, wspace=0.08)

    for i in range(n_dim):
        for j in range(n_dim):
            if j > i:
                continue

            ax = fig.add_subplot(gs[i, j])

            if i == j:
                # 1D marginal — filled KDE
                sns.kdeplot(samples[:, i], fill=True, color=color,
                            alpha=0.5, ax=ax, linewidth=1.2)
                # Percentile lines
                q16, q50, q84 = np.percentile(samples[:, i], [16, 50, 84])
                for q, ls in [(q16, ':'), (q50, '-'), (q84, ':')]:
                    ax.axvline(q, color='k', linestyle=ls, linewidth=0.8)
                # Title with median +/- uncertainties
                minus = q50 - q16
                plus = q84 - q50
                ax.set_title(f'{q50:.2f}$_{{-{minus:.2f}}}^{{+{plus:.2f}}}$',
                             fontsize=9)
                ax.set_yticks([])
                if truths is not None and truths[i] is not None:
                    ax.axvline(truths[i], color='red', linewidth=1, linestyle='--')
            else:
                # 2D KDE — filled contours
                try:
                    sns.kdeplot(x=samples[:, j], y=samples[:, i],
                                fill=True, cmap=cmap, levels=8,
                                thresh=0.05, ax=ax)
                    sns.kdeplot(x=samples[:, j], y=samples[:, i],
                                levels=[0.393, 0.865], colors=['white', 'white'],
                                linewidths=[1.0, 0.6], ax=ax)
                except Exception:
                    # Fallback to scatter if KDE fails
                    ax.scatter(samples[:, j], samples[:, i], s=1, alpha=0.2,
                               color=color)

                if truths is not None:
                    if truths[j] is not None:
                        ax.axvline(truths[j], color='red', linewidth=0.8,
                                   linestyle='--')
                    if truths[i] is not None:
                        ax.axhline(truths[i], color='red', linewidth=0.8,
                                   linestyle='--')

            # Labels
            if i == n_dim - 1:
                ax.set_xlabel(labels[j], fontsize=11)
            else:
                ax.set_xticklabels([])
            if j == 0 and i != 0:
                ax.set_ylabel(labels[i], fontsize=11)
            elif j != 0:
                ax.set_yticklabels([])

            ax.tick_params(labelsize=8)

    fig.suptitle(title, fontsize=14, y=0.98)
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    log.info(f'Corner plot saved to {output_path}')
    plt.close(fig)


# ============================================================
# k2 scatter with error ellipse
# ============================================================
def plot_k2_scatter(k2_results, heating_results, title, output_path):
    fig, ax = plt.subplots(figsize=(8, 6))
    Re_arr = k2_results[:, 0]
    Im_arr = np.abs(k2_results[:, 1])

    f_sil = []
    for h in heating_results:
        total = sum(h.values()) if h else 1e-30
        f_sil.append(h.get('Sil', 0) / total if total > 0 else 0)
    f_sil = np.array(f_sil)

    sc = ax.scatter(Re_arr, Im_arr, c=f_sil, cmap='RdYlBu_r', s=8,
                    alpha=0.6, vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax, label='Silicate heating fraction')

    ell = Ellipse((RE_K2_OBS, IM_K2_OBS), 2*RE_K2_ERR, 2*IM_K2_ERR,
                  fill=False, color='red', linewidth=2, linestyle='--',
                  label=r'Petricca 1$\sigma$')
    ax.add_patch(ell)
    ell2 = Ellipse((RE_K2_OBS, IM_K2_OBS), 4*RE_K2_ERR, 4*IM_K2_ERR,
                   fill=False, color='red', linewidth=1, linestyle=':',
                   label=r'Petricca 2$\sigma$')
    ax.add_patch(ell2)
    ax.set_xlabel(r'Re(k$_2$)')
    ax.set_ylabel(r'|Im(k$_2$)|')
    ax.set_title(title)
    ax.legend()
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    log.info(f'k2 scatter saved to {output_path}')
    plt.close(fig)


# ============================================================
# Heating distribution plots
# ============================================================
def plot_heating(samples, eval_idx, heating_results, param_names,
                 param_labels, title_prefix, output_path):
    """Heating fraction vs each parameter, plus stacked bar sorted by f_sil."""
    eval_samples = samples[eval_idx]
    n_params = len(param_names)
    n_cols = min(n_params + 1, 3)
    n_rows = (n_params + 1 + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6 * n_cols, 5 * n_rows))
    axes = np.atleast_2d(axes).flatten()

    # Compute heating fractions per phase
    heating_fracs = {}
    for phase in ['Ih', 'III', 'V', 'VI', 'Sil']:
        fracs = []
        for h in heating_results:
            total = sum(h.values()) if h else 1e-30
            fracs.append(h.get(phase, 0) / total if total > 0 else 0)
        heating_fracs[phase] = np.array(fracs)

    # Scatter panels: heating fraction vs each parameter
    for ip, (pname, plabel) in enumerate(zip(param_names, param_labels)):
        ax = axes[ip]
        x = eval_samples[:, ip]
        for phase in ['Ih', 'III', 'V', 'VI', 'Sil']:
            if phase in heating_fracs:
                ax.scatter(x, heating_fracs[phase], s=4, alpha=0.3,
                          color=PHASE_COLORS.get(phase, 'gray'),
                          label=PHASE_LABELS.get(phase, phase))
        ax.set_xlabel(plabel)
        ax.set_ylabel('Heating fraction')
        ax.set_ylim(-0.05, 1.05)
        if ip == 0:
            ax.legend(fontsize=8, loc='best')

    # Stacked bar sorted by silicate heating fraction
    ax = axes[n_params]
    f_sil = heating_fracs.get('Sil', np.zeros(len(heating_results)))
    sort_idx = np.argsort(f_sil)
    x_vals = f_sil[sort_idx]

    phases_to_plot = ['Sil', 'VI', 'V', 'III', 'Ih']
    bottom = np.zeros(len(heating_results))
    for phase in phases_to_plot:
        if phase in heating_fracs:
            sorted_fracs = heating_fracs[phase][sort_idx]
            ax.bar(range(len(sorted_fracs)), sorted_fracs, bottom=bottom,
                   color=PHASE_COLORS.get(phase, 'gray'),
                   label=PHASE_LABELS.get(phase, phase),
                   width=1.0, linewidth=0)
            bottom += sorted_fracs

    # Label some x-axis ticks with the actual silicate fraction
    n_samples = len(heating_results)
    tick_positions = np.linspace(0, n_samples - 1, 6, dtype=int)
    tick_labels = [f'{x_vals[t]:.2f}' for t in tick_positions]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels)
    ax.set_xlabel('Silicate heating fraction')
    ax.set_ylabel('Heating fraction')
    ax.legend(fontsize=8)
    ax.set_title('Per-phase heating across posterior')

    # Hide unused panels
    for ip in range(n_params + 1, len(axes)):
        axes[ip].set_visible(False)

    fig.suptitle(f'{title_prefix}: Heating Distribution', fontsize=14)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    log.info(f'Heating plot saved to {output_path}')
    plt.close(fig)


# ============================================================
# D_ocean structure plot (Maxwell only)
# ============================================================
def plot_ocean_structure(cache, samples, output_path):
    """Ocean thickness vs Tb_K with posterior density overlay."""
    fig, ax = plt.subplots(figsize=(8, 5))
    Tb_vals = sorted(cache.keys())
    D_ocean_vals = []
    for Tb in Tb_vals:
        d = cache[Tb]
        phases = d['phases']
        ci = d['changeIndices']
        n = d['n_layers']
        ocean_top = ocean_bot = None
        for i in range(n):
            s, e = ci[i], ci[i+1]
            if int(phases[s]) == 0:
                ocean_bot = d['r_m'][s]
                ocean_top = d['r_m'][e-1]
        if ocean_top is not None:
            D_ocean_vals.append((ocean_top - ocean_bot) / 1e3)
        else:
            D_ocean_vals.append(0)

    ax.plot(Tb_vals, D_ocean_vals, 'k-o', markersize=3)
    ax2 = ax.twinx()
    ax2.hist(samples[:, 3], bins=30, alpha=0.3, color='steelblue',
             density=True, label='Posterior $T_b$')
    ax2.set_ylabel('Posterior density', color='steelblue')
    ax.set_xlabel('$T_b$ (K)')
    ax.set_ylabel('$D_\\mathrm{ocean}$ (km)')
    ax.set_title('Ocean Thickness vs Basal Temperature')
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    log.info(f'Ocean structure plot saved to {output_path}')
    plt.close(fig)


# ============================================================
# Tb_K -> D_ocean conversion
# ============================================================
def build_Tb_to_Docean(cache):
    """Build interpolation from Tb_K to D_ocean (km)."""
    from scipy.interpolate import interp1d
    Tb_vals = sorted(cache.keys())
    D_ocean_vals = []
    for Tb in Tb_vals:
        d = cache[Tb]
        phases = d['phases']
        ci = d['changeIndices']
        n = d['n_layers']
        ocean_top = ocean_bot = None
        for i in range(n):
            s, e = ci[i], ci[i+1]
            if int(phases[s]) == 0:
                ocean_bot = d['r_m'][s]
                ocean_top = d['r_m'][e-1]
        if ocean_top is not None:
            D_ocean_vals.append((ocean_top - ocean_bot) / 1e3)
        else:
            D_ocean_vals.append(0)
    return interp1d(Tb_vals, D_ocean_vals, kind='linear',
                    bounds_error=False, fill_value='extrapolate')


# ============================================================
# Layer structure helpers
# ============================================================
def _extract_layer_boundaries(cache):
    """Extract per-phase layer boundaries from the structure cache.

    Returns dict: Tb_K -> list of (phase_str, r_bottom_km, r_top_km).
    """
    result = {}
    for Tb in sorted(cache.keys()):
        d = cache[Tb]
        ci = d['changeIndices']
        r_m = d['r_m']
        n = d['n_layers']
        rp = d['region_phases']
        layers = []
        for i in range(n):
            s, e = ci[i], ci[i + 1]
            r_bot = r_m[s] / 1e3      # km
            r_top = r_m[e - 1] / 1e3   # km
            phase = rp[i].replace('_conv', '')  # merge convective sublayers
            layers.append((phase, r_bot, r_top))
        result[Tb] = layers
    return result


def _phase_thicknesses(layers):
    """Sum layer thicknesses per phase group.

    Returns dict: phase_group -> thickness_km.
    Groups: 'Ih', 'HP' (III+V+VI), 'Sil', 'Ocean', 'Fe', 'Clath'.
    """
    HP_PHASES = {'III', 'V', 'VI'}
    groups = {}
    for phase, r_bot, r_top in layers:
        thick = r_top - r_bot
        if phase in HP_PHASES:
            key = 'HP'
        elif phase == '0':
            key = 'Ocean'
        else:
            key = phase
        groups[key] = groups.get(key, 0) + thick
    return groups


# ============================================================
# Interior structure cross-section vs D_ocean
# ============================================================
def plot_structure_profile(cache, samples, eval_idx, heating_results,
                           output_path):
    """Two-panel: layer boundaries vs D_ocean + per-phase heating vs D_ocean."""
    layer_data = _extract_layer_boundaries(cache)
    Tb_to_Docean = build_Tb_to_Docean(cache)
    Tb_sorted = sorted(layer_data.keys())
    D_ocean_arr = np.array([float(Tb_to_Docean(Tb)) for Tb in Tb_sorted])

    # Collect unique phases from ALL cache entries (not just median)
    all_phases = set()
    for Tb in Tb_sorted:
        for phase, _, _ in layer_data[Tb]:
            all_phases.add(phase.replace('_conv', ''))
    phase_order = [p for p in PHASE_ORDER_CANONICAL if p in all_phases]

    fig = plt.figure(figsize=(10, 11))
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 4, 3], hspace=0.12)
    ax_dens = fig.add_subplot(gs[0])
    ax_struct = fig.add_subplot(gs[1], sharex=ax_dens)
    ax_heat = fig.add_subplot(gs[2], sharex=ax_struct)

    # --- Top strip: posterior D_ocean density ---
    D_ocean_posterior = Tb_to_Docean(samples[:, 3])
    ax_dens.hist(D_ocean_posterior, bins=40, alpha=0.5, color='steelblue',
                 density=True, edgecolor='none')
    ax_dens.set_ylabel('Density', fontsize=9)
    ax_dens.set_title('Titan Interior Structure vs Ocean Thickness',
                      fontsize=13)
    ax_dens.tick_params(labelsize=8)
    plt.setp(ax_dens.get_xticklabels(), visible=False)

    # --- Middle panel: stackplot of phase thicknesses ---
    # Build thickness (km) arrays for each phase at each D_ocean
    phase_thick_km = {}
    for phase in phase_order:
        thicks = []
        for Tb in Tb_sorted:
            layers = layer_data[Tb]
            total = 0
            for p, rb, rt in layers:
                if p.replace('_conv', '') == phase:
                    total += (rt - rb)
            thicks.append(total)
        phase_thick_km[phase] = thicks

    stacks = [phase_thick_km[p] for p in phase_order]
    colors = [PHASE_COLORS.get(p, '#cccccc') for p in phase_order]
    labels = [PHASE_LABELS.get(p, p) for p in phase_order]

    polys = ax_struct.stackplot(D_ocean_arr, *stacks, colors=colors, labels=labels)
    for poly in polys:
        poly.set_linewidth(0)
        poly.set_edgecolor('none')
    ax_struct.set_ylabel('Cumulative thickness (km)', fontsize=12)
    # Reverse legend order so surface layers are on top in legend
    handles, lbls = ax_struct.get_legend_handles_labels()
    ax_struct.legend(handles[::-1], lbls[::-1],
                     loc='center left', bbox_to_anchor=(1.02, 0.5),
                     fontsize=9, frameon=True)
    ax_struct.tick_params(labelsize=9)
    plt.setp(ax_struct.get_xticklabels(), visible=False)

    # --- Bottom panel: per-phase heating vs D_ocean ---
    eval_Tb = samples[eval_idx, 3]
    eval_D_ocean = np.array([float(Tb_to_Docean(t)) for t in eval_Tb])

    heat_phases = {}
    for phase in ['Ih', 'III', 'V', 'VI', 'Sil']:
        vals = np.array([h.get(phase, 0) / 1e9 for h in heating_results])
        heat_phases[phase] = vals

    for phase in ['Sil', 'VI', 'V', 'III', 'Ih']:
        ax_heat.scatter(eval_D_ocean, heat_phases[phase], s=6, alpha=0.4,
                        color=PHASE_COLORS[phase],
                        label=PHASE_LABELS[phase])

    ax_heat.set_xlabel('$D_\\mathrm{ocean}$ (km)', fontsize=12)
    ax_heat.set_ylabel('Heating power (GW)', fontsize=12)
    ax_heat.legend(fontsize=9)
    ax_heat.tick_params(labelsize=9)
    # Place title inside the plot area to avoid overlap with panel above
    ax_heat.text(0.5, 0.97, 'Per-Phase Tidal Heating vs Ocean Thickness',
                 transform=ax_heat.transAxes, fontsize=13,
                 ha='center', va='top')

    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    log.info(f'Structure profile saved to {output_path}')
    plt.close(fig)


# ============================================================
# Per-phase thickness vs per-phase heating scatter
# ============================================================
def plot_thickness_vs_heating(cache, samples, eval_idx, heating_results,
                              output_path):
    """Scatter: phase thickness (km) vs phase heating (GW), colored by D_ocean."""
    layer_data = _extract_layer_boundaries(cache)
    Tb_to_Docean = build_Tb_to_Docean(cache)
    Tb_values = np.array(sorted(cache.keys()))

    # For each eval sample, compute per-phase thickness and heating
    thicks_Ih, thicks_HP, thicks_Sil = [], [], []
    heat_Ih, heat_HP, heat_Sil = [], [], []
    d_ocean_vals = []
    HP_PHASES = {'III', 'V', 'VI'}

    for j, si in enumerate(eval_idx):
        Tb_K = samples[si, 3]
        idx_near = np.argmin(np.abs(Tb_values - Tb_K))
        Tb_near = Tb_values[idx_near]
        layers = layer_data[Tb_near]
        pt = _phase_thicknesses(layers)

        thicks_Ih.append(pt.get('Ih', 0))
        thicks_HP.append(pt.get('HP', 0))
        thicks_Sil.append(pt.get('Sil', 0))

        h = heating_results[j]
        heat_Ih.append(h.get('Ih', 0) / 1e9)
        heat_HP.append(sum(h.get(p, 0) for p in HP_PHASES) / 1e9)
        heat_Sil.append(h.get('Sil', 0) / 1e9)

        d_ocean_vals.append(float(Tb_to_Docean(Tb_K)))

    thicks_Ih = np.array(thicks_Ih)
    thicks_HP = np.array(thicks_HP)
    thicks_Sil = np.array(thicks_Sil)
    heat_Ih = np.array(heat_Ih)
    heat_HP = np.array(heat_HP)
    heat_Sil = np.array(heat_Sil)
    d_ocean_vals = np.array(d_ocean_vals)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    vmin, vmax = d_ocean_vals.min(), d_ocean_vals.max()

    # Ice Ih: thickness vs heating
    sc = axes[0].scatter(thicks_Ih, heat_Ih, c=d_ocean_vals, cmap='viridis',
                         s=10, alpha=0.6, vmin=vmin, vmax=vmax)
    axes[0].set_xlabel('Ice Ih thickness (km)', fontsize=11)
    axes[0].set_ylabel('Heating (GW)', fontsize=11)
    axes[0].set_title(PHASE_LABELS['Ih'], fontsize=12)
    axes[0].tick_params(labelsize=9)

    # HP Ice: thickness vs heating
    axes[1].scatter(thicks_HP, heat_HP, c=d_ocean_vals, cmap='viridis',
                    s=10, alpha=0.6, vmin=vmin, vmax=vmax)
    axes[1].set_xlabel('HP ice (III+V+VI) thickness (km)', fontsize=11)
    axes[1].set_ylabel('Heating (GW)', fontsize=11)
    axes[1].set_title('HP Ice (III+V+VI)', fontsize=12)
    axes[1].tick_params(labelsize=9)

    # Silicate: HP ice thickness vs silicate heating
    axes[2].scatter(thicks_HP, heat_Sil, c=d_ocean_vals, cmap='viridis',
                    s=10, alpha=0.6, vmin=vmin, vmax=vmax)
    axes[2].set_xlabel('HP ice (III+V+VI) thickness (km)', fontsize=11)
    axes[2].set_ylabel('Heating (GW)', fontsize=11)
    axes[2].set_title(PHASE_LABELS['Sil'], fontsize=12)
    axes[2].tick_params(labelsize=9)

    # Colorbar on the far right, outside the axes
    cbar_ax = fig.add_axes([0.93, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(sc, cax=cbar_ax)
    cbar.set_label('$D_\\mathrm{ocean}$ (km)', fontsize=11)

    fig.suptitle('Maxwell Ocean: Layer Thickness vs Tidal Heating',
                 fontsize=14, y=1.02)
    fig.subplots_adjust(wspace=0.35, right=0.91)
    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    log.info(f'Thickness vs heating saved to {output_path}')
    plt.close(fig)


# ============================================================
# Main: load pickles and regenerate
# ============================================================
if __name__ == '__main__':
    # ---- Andrade no-ocean ----
    andrade_pkl = os.path.join(DIR41, 'andrade_no_ocean_mcmc.pkl')
    if os.path.exists(andrade_pkl):
        log.info('Loading Andrade no-ocean results...')
        with open(andrade_pkl, 'rb') as f:
            r = pickle.load(f)

        samples = r['samples']
        param_labels = [r'$\alpha$', r'$\log_{10}\zeta$',
                        r'$\log_{10}\eta_\mathrm{Ih}$',
                        r'$\log_{10}\eta_\mathrm{HP}$',
                        r'$\log_{10}\eta_\mathrm{sil}$']
        param_names = r['param_names']

        kde_corner(samples, param_labels,
                   'Andrade No-Ocean Titan: Posterior',
                   os.path.join(DIR41, 'andrade_no_ocean_corner.png'))

        plot_k2_scatter(r['k2_results'], r['heating_results'],
                        r'Andrade No-Ocean: k$_2$ Posterior',
                        os.path.join(DIR41, 'andrade_no_ocean_k2_scatter.png'))

        plot_heating(samples, r['eval_idx'], r['heating_results'],
                     param_names, param_labels,
                     'Andrade No-Ocean',
                     os.path.join(DIR41, 'andrade_no_ocean_heating.png'))
    else:
        log.warning(f'Not found: {andrade_pkl}')

    # ---- Maxwell ocean ----
    maxwell_pkl = os.path.join(DIR42, 'maxwell_ocean_mcmc.pkl')
    cache_pkl = os.path.join(DIR42, 'maxwell_ocean_structure_cache.pkl')
    if os.path.exists(maxwell_pkl) and os.path.exists(cache_pkl):
        log.info('Loading Maxwell ocean results...')
        with open(maxwell_pkl, 'rb') as f:
            r = pickle.load(f)
        with open(cache_pkl, 'rb') as f:
            cache = pickle.load(f)

        samples = r['samples']
        param_names = r['param_names']

        # Convert Tb_K column to D_ocean
        Tb_to_Docean = build_Tb_to_Docean(cache)
        samples_plot = samples.copy()
        samples_plot[:, 3] = Tb_to_Docean(samples[:, 3])

        param_labels_ocean = [
            r'$\log_{10}\eta_\mathrm{Ih}$',
            r'$\log_{10}\eta_\mathrm{HP}$',
            r'$\log_{10}\eta_\mathrm{sil}$',
            r'$D_\mathrm{ocean}$ (km)',
        ]
        param_names_ocean = list(param_names)
        param_names_ocean[3] = 'D_ocean_km'

        kde_corner(samples_plot, param_labels_ocean,
                   'Maxwell Ocean Titan: Posterior',
                   os.path.join(DIR42, 'maxwell_ocean_corner.png'))

        plot_k2_scatter(r['k2_results'], r['heating_results'],
                        r'Maxwell Ocean: k$_2$ Posterior',
                        os.path.join(DIR42, 'maxwell_ocean_k2_scatter.png'))

        plot_heating(samples_plot, r['eval_idx'], r['heating_results'],
                     param_names_ocean, param_labels_ocean,
                     'Maxwell Ocean',
                     os.path.join(DIR42, 'maxwell_ocean_heating.png'))

        plot_ocean_structure(cache, samples,
                            os.path.join(DIR42, 'maxwell_ocean_Tb_structure.png'))

        plot_structure_profile(cache, samples, r['eval_idx'], r['heating_results'],
                               os.path.join(DIR42, 'maxwell_ocean_structure_profile.png'))

        plot_thickness_vs_heating(cache, samples, r['eval_idx'], r['heating_results'],
                                  os.path.join(DIR42, 'maxwell_ocean_thickness_vs_heating.png'))
    else:
        log.warning(f'Not found: {maxwell_pkl} or {cache_pkl}')

    # ---- Andrade Arrhenius no-ocean (PPTest43) ----
    arr_andrade_pkl = os.path.join(DIR43, 'andrade_arrhenius_no_ocean_mcmc.pkl')
    if os.path.exists(arr_andrade_pkl):
        log.info('Loading Andrade Arrhenius no-ocean results...')
        with open(arr_andrade_pkl, 'rb') as f:
            r = pickle.load(f)

        samples = r['samples']
        param_labels = [r'$\alpha$', r'$\log_{10}\zeta$',
                        r'$\log_{10}\eta_\mathrm{melt,Ih}$',
                        r'$\log_{10}\eta_\mathrm{melt,HP}$',
                        r'$\log_{10}\eta_\mathrm{melt,sil}$']
        param_names = r['param_names']

        kde_corner(samples, param_labels,
                   'Andrade Arrhenius No-Ocean Titan: Posterior',
                   os.path.join(DIR43, 'andrade_arrhenius_no_ocean_corner.png'))

        plot_k2_scatter(r['k2_results'], r['heating_results'],
                        r'Andrade Arrhenius No-Ocean: k$_2$ Posterior',
                        os.path.join(DIR43, 'andrade_arrhenius_no_ocean_k2_scatter.png'))

        plot_heating(samples, r['eval_idx'], r['heating_results'],
                     param_names, param_labels,
                     'Andrade Arrhenius No-Ocean',
                     os.path.join(DIR43, 'andrade_arrhenius_no_ocean_heating.png'))
    else:
        log.warning(f'Not found: {arr_andrade_pkl}')

    # ---- Maxwell Arrhenius ocean (PPTest44) ----
    arr_maxwell_pkl = os.path.join(DIR44, 'maxwell_arrhenius_ocean_mcmc.pkl')
    arr_cache_pkl = os.path.join(DIR44, 'maxwell_arrhenius_ocean_structure_cache.pkl')
    if os.path.exists(arr_maxwell_pkl) and os.path.exists(arr_cache_pkl):
        log.info('Loading Maxwell Arrhenius ocean results...')
        with open(arr_maxwell_pkl, 'rb') as f:
            r = pickle.load(f)
        with open(arr_cache_pkl, 'rb') as f:
            cache = pickle.load(f)

        samples = r['samples']
        param_names = r['param_names']

        # Convert Tb_K column to D_ocean
        Tb_to_Docean = build_Tb_to_Docean(cache)
        samples_plot = samples.copy()
        samples_plot[:, 3] = Tb_to_Docean(samples[:, 3])

        param_labels_ocean = [
            r'$\log_{10}\eta_\mathrm{melt,Ih}$',
            r'$\log_{10}\eta_\mathrm{melt,HP}$',
            r'$\log_{10}\eta_\mathrm{melt,sil}$',
            r'$D_\mathrm{ocean}$ (km)',
        ]
        param_names_ocean = list(param_names)
        param_names_ocean[3] = 'D_ocean_km'

        kde_corner(samples_plot, param_labels_ocean,
                   'Maxwell Arrhenius Ocean Titan: Posterior',
                   os.path.join(DIR44, 'maxwell_arrhenius_ocean_corner.png'))

        plot_k2_scatter(r['k2_results'], r['heating_results'],
                        r'Maxwell Arrhenius Ocean: k$_2$ Posterior',
                        os.path.join(DIR44, 'maxwell_arrhenius_ocean_k2_scatter.png'))

        plot_heating(samples_plot, r['eval_idx'], r['heating_results'],
                     param_names_ocean, param_labels_ocean,
                     'Maxwell Arrhenius Ocean',
                     os.path.join(DIR44, 'maxwell_arrhenius_ocean_heating.png'))

        plot_ocean_structure(cache, samples,
                            os.path.join(DIR44, 'maxwell_arrhenius_ocean_Tb_structure.png'))

        plot_structure_profile(cache, samples, r['eval_idx'], r['heating_results'],
                               os.path.join(DIR44, 'maxwell_arrhenius_ocean_structure_profile.png'))

        plot_thickness_vs_heating(cache, samples, r['eval_idx'], r['heating_results'],
                                  os.path.join(DIR44, 'maxwell_arrhenius_ocean_thickness_vs_heating.png'))
    else:
        log.warning(f'Not found: {arr_maxwell_pkl} or {arr_cache_pkl}')

    log.info('Done.')
