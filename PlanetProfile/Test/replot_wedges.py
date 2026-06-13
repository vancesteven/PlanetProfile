"""
Regenerate structure-wedge plots for all MCMC tests that have existing pkl results.
Does NOT import TidalPy — safe to run in any environment.

Usage:
    python PlanetProfile/Test/replot_wedges.py
"""
import logging
import os
import pickle
import sys

import matplotlib
matplotlib.use('Agg')

_pp_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _pp_root)

logging.basicConfig(level=logging.INFO, format='%(name)s - %(message)s')
log = logging.getLogger('replot_wedges')

import importlib.util as _ilu
_spec = _ilu.spec_from_file_location(
    'mcmc_plots',
    os.path.join(_pp_root, 'PlanetProfile', 'Inference', 'mcmc_plots.py')
)
mp = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(mp)

_BASE = os.path.join(_pp_root, 'PlanetProfile', 'Test', 'mcmc_results')

TITAN_R_KM = 2575.5
EUROPA_R_KM = 1560.8
GANYMEDE_R_KM = 2634.1
CALLISTO_R_KM = 2410.3


def _load(path):
    with open(path, 'rb') as f:
        return pickle.load(f)


def _load_grid(path):
    """Load grid cache pkl, unwrapping the 'grid_cache' envelope if present."""
    d = _load(path)
    if isinstance(d, dict) and 'grid_cache' in d and isinstance(d['grid_cache'], dict):
        return d['grid_cache']
    return d


def replot_titan_test46():
    pkl    = os.path.join(_BASE, 'Titan', 'Test46_andrade_hybrid_hydro', 'hybrid_hydro_andrade_mcmc.pkl')
    grid   = os.path.join(_BASE, 'Titan', 'Test46_andrade_hybrid_hydro', 'titan_maxwell_hybrid_hydro_grid_cache.pkl')
    outdir = os.path.join(_BASE, 'Titan', 'Test46_andrade_hybrid_hydro')
    if not (os.path.exists(pkl) and os.path.exists(grid)):
        log.warning('Test46 andrade hybrid: pkls not found, skipping.')
        return
    r = _load(pkl)
    grid_cache = _load_grid(grid)
    mp.plot_structure_wedge(
        r['samples'], r['eval_idx'], grid_cache,
        output_path=os.path.join(outdir, 'hybrid_hydro_andrade_structure_profile.png'),
        R_body_km=TITAN_R_KM,
        body_name='Titan',
        param_indices={'Tb': 5, 'D_hydro': 6, 'f_core': 9},
    )


def replot_titan_test48():
    pkl    = os.path.join(_BASE, 'Titan', 'Test48_andrade_yao2014', 'hybrid_hydro_andrade_yao2014_mcmc.pkl')
    grid   = os.path.join(_BASE, 'Titan', 'Test48_andrade_yao2014', 'titan_yao2014_hybrid_hydro_grid_cache.pkl')
    outdir = os.path.join(_BASE, 'Titan', 'Test48_andrade_yao2014')
    if not (os.path.exists(pkl) and os.path.exists(grid)):
        log.warning('Test48 yao2014: pkls not found, skipping.')
        return
    r = _load(pkl)
    grid_cache = _load_grid(grid)
    mp.plot_structure_wedge(
        r['samples'], r['eval_idx'], grid_cache,
        output_path=os.path.join(outdir, 'hybrid_hydro_andrade_yao2014_structure_profile.png'),
        R_body_km=TITAN_R_KM,
        body_name='Titan',
        param_indices={'Tb': 5, 'D_hydro': 6, 'f_core': 9},
    )


def replot_titan_test49(tag, pkl_name, grid_name, outname):
    subdir = 'Test49_clathrate4km' if '4km' in tag else 'Test49_clathrate2km'
    dirpath = os.path.join(_BASE, 'Titan', subdir)
    pkl  = os.path.join(dirpath, pkl_name)
    grid = os.path.join(dirpath, grid_name)
    if not (os.path.exists(pkl) and os.path.exists(grid)):
        log.warning(f'{tag}: pkls not found, skipping.')
        return
    import numpy as np
    r = _load(pkl)
    grid_cache = _load_grid(grid)
    samples = r['samples']
    eval_idx = r.get('eval_idx', np.arange(len(samples)))
    mp.plot_structure_wedge(
        samples, eval_idx, grid_cache,
        output_path=os.path.join(dirpath, outname),
        R_body_km=TITAN_R_KM,
        body_name='Titan',
        param_indices={'Tb': 5, 'D_hydro': 6, 'f_core': 9},
    )


def replot_europa_test51():
    """Generate the Test51 Europa seawater structure-wedge plot.

    Test51 uses the v2.1 cache format (``{'Tb_K_grid': ..., 'structures': [...]}``)
    and an ``InferenceResult`` pkl whose posterior parameter space is
    ``['alpha', 'log10_zeta', 'log10_eta_Ih', 'log10_eta_sil', 'Tb_K',
    'R_core_km', 'rho_core_kgm3']``.  This is incompatible with the legacy
    ``plot_structure_wedge`` signature which expects ``D_hydro`` and ``f_core``
    as direct sample columns and a ``(Tb, D)``-keyed grid.

    Instead, this function builds layer radii per-sample using:
      - ``InferenceResult.D_iceIh_results`` (already computed by the runner)
      - ``InferenceResult.D_ocean_results`` (already computed by the runner)
      - HP-ice thicknesses looked up from the v2.1 grid by nearest ``Tb_K``
      - Core radius from the sampled ``R_core_km`` (param index 5)
    then calls ``mp.plot_structure_wedge`` with a pre-keyed ``(Tb, D)``-style
    grid that it can consume, or draws the wedge directly via the underlying
    matplotlib patches.
    """
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    from matplotlib.patches import Wedge as MplWedge

    dirpath = os.path.join(_BASE, 'Europa', 'Test51_seawater')
    pkl  = os.path.join(dirpath, 'test51_europa_seawater_results.pkl')
    grid = os.path.join(dirpath, 'europa_seawater_structure_grid.pkl')
    if not os.path.exists(grid):
        log.warning('Test51 Europa: structure grid not found, skipping.')
        return
    if not os.path.exists(pkl):
        log.warning('Test51 Europa: results pkl not found, skipping.')
        return
    try:
        r = _load(pkl)
        raw_grid = _load(grid)
    except Exception as e:
        log.warning(f'Test51 Europa: load failed ({e}), skipping.')
        return

    # --- Extract posterior arrays -------------------------------------------
    # InferenceResult: param order is
    #   [alpha, log10_zeta, log10_eta_Ih, log10_eta_sil, Tb_K, R_core_km, rho_core_kgm3]
    if isinstance(r, dict):
        samples       = np.asarray(r['samples'])
        D_iceIh_arr   = r.get('D_iceIh_results')
        D_ocean_arr   = r.get('D_ocean_results')
        i_Tb          = r.get('i_Tb', 4)
        i_Rcore       = r.get('i_Rcore', 5)
    else:
        samples       = np.asarray(r.samples)
        # Try runner-computed arrays first (most reliable)
        D_iceIh_arr   = getattr(r, 'D_iceIh_results', None)
        D_ocean_arr   = getattr(r, 'D_ocean_results', None)
        pnames        = list(getattr(r, 'param_names', []))
        i_Tb          = pnames.index('Tb_K')    if 'Tb_K'       in pnames else 4
        i_Rcore       = pnames.index('R_core_km') if 'R_core_km' in pnames else 5

    n_samples = len(samples)
    if n_samples == 0:
        log.warning('Test51 Europa: no posterior samples, skipping.')
        return

    # --- Parse v2.1 grid for HP-ice thicknesses by Tb ------------------------
    Tb_grid_arr   = None
    struct_list   = None
    if isinstance(raw_grid, dict):
        if 'Tb_K_grid' in raw_grid and 'structures' in raw_grid:
            Tb_grid_arr = np.asarray(raw_grid['Tb_K_grid'], dtype=float)
            struct_list = raw_grid['structures']
        elif 'grid_cache' in raw_grid and isinstance(raw_grid['grid_cache'], dict):
            # v2.0 envelope — convert to Tb-keyed list
            inner = raw_grid['grid_cache']
            Tb_grid_arr = np.array(sorted(inner.keys()), dtype=float)
            struct_list = [inner[tb] for tb in Tb_grid_arr]

    def _hp_ice_from_grid(Tb_K: float):
        """Return (D_iceII, D_iceIII, D_iceV, D_iceVI) km from nearest grid point."""
        if Tb_grid_arr is None or struct_list is None:
            return 0.0, 0.0, 0.0, 0.0
        idx = int(np.argmin(np.abs(Tb_grid_arr - Tb_K)))
        pt  = struct_list[idx]
        return (
            float(pt.get('D_iceII_km',  pt.get('D_iceIII_km', 0.0))),
            float(pt.get('D_iceIII_km', 0.0)),
            float(pt.get('D_iceV_km',   0.0)),
            float(pt.get('D_iceVI_km',  0.0)),
        )

    # --- Fall back to grid-derived D_iceIh / D_ocean when not in result ------
    if D_iceIh_arr is None or len(D_iceIh_arr) != n_samples:
        if Tb_grid_arr is not None and struct_list is not None:
            D_iceIh_arr = np.array([
                float(struct_list[int(np.argmin(np.abs(Tb_grid_arr - samples[i, i_Tb])))].get(
                    'D_iceIh_km', 0.0))
                for i in range(n_samples)
            ])
        else:
            D_iceIh_arr = np.zeros(n_samples)
    else:
        D_iceIh_arr = np.asarray(D_iceIh_arr, dtype=float)

    if D_ocean_arr is None or len(D_ocean_arr) != n_samples:
        if Tb_grid_arr is not None and struct_list is not None:
            D_ocean_arr = np.array([
                float(struct_list[int(np.argmin(np.abs(Tb_grid_arr - samples[i, i_Tb])))].get(
                    'D_ocean_km', 0.0))
                for i in range(n_samples)
            ])
        else:
            D_ocean_arr = np.zeros(n_samples)
    else:
        D_ocean_arr = np.asarray(D_ocean_arr, dtype=float)

    # --- Per-sample layer radii (km) -----------------------------------------
    r_iceIh_bot  = []
    r_ocean_bot  = []
    r_iceIII_bot = []
    r_iceV_bot   = []
    r_iceVI_bot  = []
    r_sil_bot    = []
    r_core_bot   = []

    for i in range(n_samples):
        Tb_K_i   = float(samples[i, i_Tb])
        R_core_i = float(samples[i, i_Rcore])  # km
        D_ih     = float(D_iceIh_arr[i])
        D_oc     = float(D_ocean_arr[i])
        _d2, _d3, _d5, _d6 = _hp_ice_from_grid(Tb_K_i)

        r_ih  = EUROPA_R_KM - D_ih
        r_oc  = r_ih - D_oc
        r_iii = r_oc - _d3
        r_v   = r_iii - _d5
        r_vi  = r_v   - _d6
        # silicate base = core top
        r_core_top = R_core_i  # sampled directly as km

        r_iceIh_bot.append(r_ih)
        r_ocean_bot.append(r_oc)
        r_iceIII_bot.append(r_iii)
        r_iceV_bot.append(r_v)
        r_iceVI_bot.append(r_vi)
        r_sil_bot.append(r_core_top)
        r_core_bot.append(0.0)

    if not r_iceIh_bot:
        log.warning('Test51 Europa: no valid samples after layer computation, skipping.')
        return

    def pct(arr):
        return np.percentile(arr, [5, 50, 95])

    p_ih  = pct(r_iceIh_bot)
    p_oc  = pct(r_ocean_bot)
    p_iii = pct(r_iceIII_bot)
    p_v   = pct(r_iceV_bot)
    p_vi  = pct(r_iceVI_bot)
    p_sil = pct(r_sil_bot)

    layers = [
        ('Ice Ih',     EUROPA_R_KM, p_ih[1]),
        ('Ocean',      p_ih[1],     p_oc[1]),
        ('Ice III',    p_oc[1],     p_iii[1]),
        ('Ice V',      p_iii[1],    p_v[1]),
        ('Ice VI',     p_v[1],      p_vi[1]),
        ('Silicate',   p_vi[1],     p_sil[1]),
        ('Dense core', p_sil[1],    0.0),
    ]

    wedge_colors = mp._wedge_color_map()

    ANG1, ANG2 = 55, 125
    fig, ax = plt.subplots(figsize=(6, 8))
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.05, 1.05)
    ax.set_aspect('equal')
    ax.axis('off')
    cx, cy = 0.5, 0.0

    for name, r_top, r_bot in layers:
        thickness_km = r_top - r_bot
        if thickness_km < 0.5:
            continue
        r_outer_norm = r_top / EUROPA_R_KM
        width_norm   = thickness_km / EUROPA_R_KM
        wedge = MplWedge(
            (cx, cy), r_outer_norm, ANG1, ANG2,
            width=width_norm,
            fc=wedge_colors.get(name, '#888888'), ec='#333333', lw=0.8,
        )
        ax.add_patch(wedge)

    tick_angle_rad = np.radians(ANG1 - 3)
    label_entries  = []
    for name, r_top, r_bot in layers:
        if r_top - r_bot < 5:
            continue
        r_mid  = (r_top + r_bot) / 2
        r_norm = r_mid / EUROPA_R_KM
        x = cx + r_norm * np.cos(tick_angle_rad)
        y = cy + r_norm * np.sin(tick_angle_rad)
        label_entries.append((x, y, name, r_top - r_bot))

    label_y_positions = [e[1] for e in label_entries]
    spaced_y = []
    for i, y in enumerate(sorted(label_y_positions)):
        if i > 0 and y - spaced_y[-1] < 0.06:
            y = spaced_y[-1] + 0.06
        spaced_y.append(y)
    sorted_idx = sorted(range(len(label_y_positions)), key=lambda k: label_y_positions[k])
    final_y    = [0.0] * len(label_entries)
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
        cx, EUROPA_R_KM / EUROPA_R_KM + 0.02 + cy,
        f'R = {EUROPA_R_KM:.0f} km',
        ha='center', va='bottom', fontsize=9,
    )

    ice_ih_thick = EUROPA_R_KM - p_ih[1]
    ocean_thick  = p_ih[1] - p_oc[1]
    hp_thick     = p_oc[1] - p_vi[1]
    ax.set_title(
        'Europa Interior Structure (Posterior Median)\n'
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

    outpath = os.path.join(dirpath, 'europa_test51_structure_wedge.png')
    fig.tight_layout()
    fig.savefig(outpath, dpi=150, bbox_inches='tight')
    log.info(f'Saved {outpath}')
    plt.close(fig)


if __name__ == '__main__':
    log.info('Regenerating structure-wedge plots with updated PP colors...')

    replot_titan_test46()
    replot_titan_test48()

    replot_titan_test49(
        tag='Test49_4km',
        pkl_name='test49_clathrate4km_mcmc_results.pkl',
        grid_name='titan_yao2014_clathrate4km_hybrid_hydro_grid_cache.pkl',
        outname='test49_clathrate4km_structure_wedge.png',
    )
    replot_titan_test49(
        tag='Test49_2km',
        pkl_name='test49_clathrate2km_mcmc_results.pkl',
        grid_name='titan_yao2014_clathrate2km_hybrid_hydro_grid_cache.pkl',
        outname='test49_clathrate2km_structure_wedge.png',
    )

    replot_europa_test51()

    log.info('Done.')
