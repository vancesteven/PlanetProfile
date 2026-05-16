"""
MCMC exploration: Andrade rheology, all-ice Titan with Yao 2014 spherical Ih convection (PPTest50).

Explores Titan's interior structure under the hypothesis of NO liquid ocean —
the entire hydrosphere is ice (Ih + III + V + VI), as in Petricca et al. (2025).

With NO_OCEAN_EXCEPT_INNER_ICES=True and Fe_CORE=False, the ice structure is
self-consistently determined by PP and barely varies with Tb_K (< 1 km over the
entire valid range 240–251 K). The all-ice hypothesis produces a UNIQUE
structure: D_hsphere ≈ 494 km, rho_sil ≈ 2555 kg/m³, CMR2 ≈ 0.343.

Default config: PlanetProfile/Inference/configs/test50_titan_noocean_andrade_8D.json
5D variant   : PlanetProfile/Inference/configs/test50_titan_noocean_andrade_5D.json

Usage:
  mamba activate PPcl
  # Run default 8D inference
  python PlanetProfile/Test/Test50_mcmc_andrade_noocean_yao2014.py

  # Run with a specific config (e.g. 5D variant)
  python PlanetProfile/Test/Test50_mcmc_andrade_noocean_yao2014.py \\
      --config PlanetProfile/Inference/configs/test50_titan_noocean_andrade_5D.json

  # Rebuild structure cache and run
  python PlanetProfile/Test/Test50_mcmc_andrade_noocean_yao2014.py --rebuild-structure

  # Build structure only, no MCMC
  python PlanetProfile/Test/Test50_mcmc_andrade_noocean_yao2014.py --build-structure

  # Regenerate plots from saved results
  python PlanetProfile/Test/Test50_mcmc_andrade_noocean_yao2014.py --replot
"""
import argparse
import logging
import os
import pickle
import sys
import tempfile
import time
import importlib

import numpy as np

# --- Environment setup ---
import platformdirs
_pp_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
platformdirs.user_documents_dir = lambda: _pp_root
sys.path.insert(0, _pp_root)

logging.basicConfig(level=logging.WARNING, format='%(name)s - %(message)s')
log = logging.getLogger('PPTest50_MCMC_Yao2014')
log.setLevel(logging.INFO)

from PlanetProfile.Utilities.defineStructs import EOSlist, Constants
from PlanetProfile.Utilities.Indexing import PhaseConv
from PlanetProfile.Main import PlanetProfile as RunPP
from PlanetProfile.Gravity.Gravity import SetupGravity

# ============================================================
# Configuration
# ============================================================
OUTPUT_DIR = os.path.join(_pp_root, 'PlanetProfile', 'Test', 'mcmc_results')
os.makedirs(OUTPUT_DIR, exist_ok=True)

STRUCTURE_CACHE_PATH = os.path.join(OUTPUT_DIR, 'titan_allice_yao2014_structure_grid.pkl')

# Default config path (8D)
_CONFIGS_DIR = os.path.join(_pp_root, 'PlanetProfile', 'Inference', 'configs')
DEFAULT_CONFIG_JSON = os.path.join(_CONFIGS_DIR, 'test50_titan_noocean_andrade_8D.json')

# Observational constraints (Petricca et al. 2025, Cassini radio science)
OBS = {
    'Re_k2': 0.608,  'Re_k2_err': 0.048,
    'Im_k2': 0.135,  'Im_k2_err': 0.035,
    'CMR2':  0.343,  'CMR2_err':  0.001,
}

# Titan reference
TITAN_R_M = 2574.73e3

# Tb_K grid for structural interpolation.  Upper bound = Ih-III-L triple
# point minus epsilon (grid-resolution-safe for nIceI=200).  Lower bound
# simulates solute-driven triple-point lowering.
STRUCTURE_TB_GRID = np.linspace(249.0, 250.965, 9)
STRUCTURE_TB_K = float(STRUCTURE_TB_GRID[-1])  # back-compat alias for logging

# MCMC / output settings
RANDOM_STATE = 42
N_REEVAL = 500

# Canonical output pkl name (8D default)
PKL_PATH = os.path.join(OUTPUT_DIR, 'allice_yao2014_andrade_mcmc_results.pkl')


# ============================================================
# Step 1: Build / load structural grid
# ============================================================

def _save_cache(data, filepath):
    """Atomic pickle write."""
    filepath = str(filepath)
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    fd, tmp_path = tempfile.mkstemp(dir=os.path.dirname(filepath), suffix='.pkl.tmp')
    try:
        with os.fdopen(fd, 'wb') as f:
            pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp_path, filepath)
    except BaseException:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)
        raise


def _build_single_structure(Tb_K):
    """Run PlanetProfile for one Tb_K point and extract arrays."""
    EOSlist.loaded.clear()
    from PlanetProfile.GetConfig import Params as configParams

    test_module = 'PlanetProfile.Test.PPTest50'
    if test_module in sys.modules:
        importlib.reload(sys.modules[test_module])
    else:
        importlib.import_module(test_module)
    mod = sys.modules[test_module]
    from copy import deepcopy
    Planet = deepcopy(mod.Planet)

    Planet.Bulk.Tb_K = Tb_K

    configParams.Gravity.backend = 'tidalpy'
    configParams.CALC_NEW = True
    configParams.CALC_NEW_GRAVITY = False
    configParams.NO_SAVEFILE = True
    configParams.SKIP_PLOTS = True

    Planet, Params = RunPP(Planet, configParams)
    Params.CALC_NEW_GRAVITY = True
    Planet, Params = SetupGravity(Planet, Params)

    model = Planet.Gravity.ALMAModel['model']
    cols = Planet.Gravity.columns
    rIndex = cols.index('r')
    rhoIndex = cols.index('rho')
    VPIndex = cols.index('VP')
    GSIndex = cols.index('GS')
    etaIndex = cols.index('eta')
    pIndex = cols.index('phase')

    r_m = model[:, rIndex].astype(np.float64)
    rho = model[:, rhoIndex].astype(np.float64)
    mu_Pa = model[:, GSIndex].astype(np.float64)
    VP_ms = model[:, VPIndex].astype(np.float64)
    eta_Pa_base = model[:, etaIndex].astype(np.float64)
    phases = model[:, pIndex]

    try:
        T_K_primary = np.asarray(Planet.T_K, dtype=np.float64)
        r_m_primary = np.asarray(Planet.r_m[:T_K_primary.size], dtype=np.float64)
        sort_idx = np.argsort(r_m_primary)
        T_K = np.interp(r_m, r_m_primary[sort_idx], T_K_primary[sort_idx])
    except (AttributeError, ValueError) as _exc:
        T_K = np.full(r_m.size, np.nan)
        log.warning(f'Planet.T_K extraction failed ({_exc}) — Arrhenius Ih viscosity will be skipped.')

    try:
        P_MPa_primary = np.asarray(Planet.P_MPa[:T_K_primary.size], dtype=np.float64)
        P_MPa = np.interp(r_m, r_m_primary[sort_idx], P_MPa_primary[sort_idx])
    except (AttributeError, ValueError, NameError) as _exc:
        P_MPa = np.full(r_m.size, np.nan)
        log.warning(f'Planet.P_MPa extraction failed ({_exc}) — no-ocean safeguard will be disabled.')

    K_Pa = rho * VP_ms**2 - (4.0 / 3.0) * mu_Pa
    nan_mask = ~np.isfinite(K_Pa) | (K_Pa <= 0)
    if np.any(nan_mask):
        for i in np.where(nan_mask)[0]:
            ph = int(phases[i])
            if 50 <= ph < 100:
                nu = 0.25
            elif ph >= 100:
                nu = 0.29
            else:
                nu = 0.33
            K_Pa[i] = 2.0 * mu_Pa[i] * (1.0 + nu) / (3.0 * (1.0 - 2.0 * nu))
    K_Pa = np.maximum(K_Pa, 1e6)

    changeIndices = np.max(Planet.Reduced.changeIndices) - np.flipud(Planet.Reduced.changeIndices)
    n_layers = len(changeIndices) - 1

    MIN_POINTS = 5
    needs_padding = any(
        changeIndices[i+1] - changeIndices[i] < MIN_POINTS
        for i in range(n_layers)
    )

    _orig_iConv = np.flipud(Planet.Reduced.iConv)
    region_phases = []
    for i_layer in range(n_layers):
        start = changeIndices[i_layer]
        phase = phases[start]
        if phase >= Constants.phaseClath and phase < Constants.phaseClath + 10:
            phase = Constants.phaseClath
        convection = _orig_iConv[start]
        phase_str = PhaseConv(phase, liq='0')
        if convection:
            phase_str += '_conv'
        region_phases.append(phase_str)

    bulk_visc = np.zeros_like(eta_Pa_base)

    if needs_padding:
        new_r, new_rho, new_K, new_mu, new_eta, new_phases, new_bv, new_T = \
            [], [], [], [], [], [], [], []
        new_ci = [0]
        for i_layer in range(n_layers):
            s, e = changeIndices[i_layer], changeIndices[i_layer + 1]
            n_pts = e - s
            if n_pts < MIN_POINTS and n_pts >= 2:
                r_layer = r_m[s:e]
                r_interp = np.linspace(r_layer[0], r_layer[-1], MIN_POINTS)
                new_r.append(r_interp)
                new_rho.append(np.interp(r_interp, r_layer, rho[s:e]))
                new_K.append(np.interp(r_interp, r_layer, K_Pa[s:e]))
                new_mu.append(np.interp(r_interp, r_layer, mu_Pa[s:e]))
                new_eta.append(np.interp(r_interp, r_layer, eta_Pa_base[s:e]))
                new_phases.append(np.full(MIN_POINTS, phases[s]))
                new_bv.append(np.zeros(MIN_POINTS))
                new_T.append(np.interp(r_interp, r_layer, T_K[s:e]))
                new_ci.append(new_ci[-1] + MIN_POINTS)
            elif n_pts == 1:
                r_top = r_m[s]
                r_bot = r_m[s - 1] if s > 0 else r_top - 1.0
                if r_bot >= r_top:
                    r_bot = r_top - 1.0
                r_interp = np.linspace(r_bot, r_top, MIN_POINTS)
                new_r.append(r_interp)
                new_rho.append(np.full(MIN_POINTS, rho[s]))
                new_K.append(np.full(MIN_POINTS, K_Pa[s]))
                new_mu.append(np.full(MIN_POINTS, mu_Pa[s]))
                new_eta.append(np.full(MIN_POINTS, eta_Pa_base[s]))
                new_phases.append(np.full(MIN_POINTS, phases[s]))
                new_bv.append(np.zeros(MIN_POINTS))
                T_lo = T_K[s - 1] if s > 0 and np.isfinite(T_K[s - 1]) else T_K[s]
                T_hi = T_K[s] if np.isfinite(T_K[s]) else T_lo
                new_T.append(np.linspace(T_lo, T_hi, MIN_POINTS))
                new_ci.append(new_ci[-1] + MIN_POINTS)
            else:
                new_r.append(r_m[s:e])
                new_rho.append(rho[s:e])
                new_K.append(K_Pa[s:e])
                new_mu.append(mu_Pa[s:e])
                new_eta.append(eta_Pa_base[s:e])
                new_phases.append(phases[s:e])
                new_bv.append(bulk_visc[s:e])
                new_T.append(T_K[s:e])
                new_ci.append(new_ci[-1] + (e - s))

        r_m = np.concatenate(new_r)
        rho = np.concatenate(new_rho)
        K_Pa = np.concatenate(new_K)
        mu_Pa = np.concatenate(new_mu)
        eta_Pa_base = np.concatenate(new_eta)
        phases = np.concatenate(new_phases)
        bulk_visc = np.concatenate(new_bv)
        T_K = np.concatenate(new_T)
        changeIndices = np.array(new_ci)

    layer_upper_radii = []
    layer_types = []
    for i_layer in range(n_layers):
        end = changeIndices[i_layer + 1]
        layer_upper_radii.append(r_m[end - 1])
        layer_types.append('liquid' if phases[changeIndices[i_layer]] == 0 else 'solid')

    omega = Planet.Bulk.meanMotion_radps
    ecc = Planet.Bulk.eccentricity
    host_mass = Constants.parentMass_kg[Planet.parent]
    a_m = (Constants.G * host_mass / omega**2) ** (1.0 / 3.0)

    D_iceIh_km = D_iceIII_km = D_iceV_km = D_iceVI_km = 0.0
    for i_layer in range(n_layers):
        s = changeIndices[i_layer]
        e = changeIndices[i_layer + 1]
        ph = int(phases[s])
        thick_km = (r_m[e-1] - r_m[s]) / 1e3
        if ph == 1:
            D_iceIh_km += thick_km
        elif ph == 3:
            D_iceIII_km += thick_km
        elif ph == 5:
            D_iceV_km += thick_km
        elif ph == 6:
            D_iceVI_km += thick_km

    CMR2_pp = Planet.Bulk.Cmeasured if hasattr(Planet.Bulk, 'Cmeasured') else np.nan
    try:
        CMR2_pp = float(Planet.CMR2mean)
    except (AttributeError, TypeError):
        pass

    R_sil = getattr(Planet.Sil, 'Rmean_m', r_m[0])
    D_hsphere_km = (TITAN_R_M - R_sil) / 1e3

    return {
        'r_m': np.ascontiguousarray(r_m),
        'rho': np.ascontiguousarray(rho),
        'K_Pa': np.ascontiguousarray(K_Pa),
        'mu_Pa': np.ascontiguousarray(mu_Pa),
        'eta_Pa_base': eta_Pa_base,
        'phases': phases,
        'T_K': np.ascontiguousarray(T_K),
        'P_MPa': np.ascontiguousarray(P_MPa),
        'bulk_visc': np.ascontiguousarray(bulk_visc),
        'changeIndices': changeIndices,
        'n_layers': n_layers,
        'layer_upper_radii': tuple(layer_upper_radii),
        'layer_types': tuple(layer_types),
        'region_phases': region_phases,
        'omega': omega,
        'eccentricity': ecc,
        'host_mass': host_mass,
        'a_m': a_m,
        'R_body_m': TITAN_R_M,
        'Mtot_kg': Planet.Bulk.M_kg,
        'CMR2': CMR2_pp,
        'Tb_K': Tb_K,
        'rhoSil_kgm3': getattr(Planet.Sil, 'rhoMean_kgm3', np.nan),
        'D_hsphere_km': D_hsphere_km,
        'D_iceIh_km': D_iceIh_km,
        'D_iceIII_km': D_iceIII_km,
        'D_iceV_km': D_iceV_km,
        'D_iceVI_km': D_iceVI_km,
    }


def build_or_load_structure_grid(force_rebuild=False):
    """Build or load a grid of all-ice structures across STRUCTURE_TB_GRID.

    Returns dict with 'Tb_K_grid' and 'structures' keys, consumed by MCMCRunner
    via the Tb_K parameter hook in forward_models.py.
    """
    if not force_rebuild and os.path.exists(STRUCTURE_CACHE_PATH):
        try:
            with open(STRUCTURE_CACHE_PATH, 'rb') as f:
                grid = pickle.load(f)
            if (isinstance(grid, dict)
                    and 'Tb_K_grid' in grid and 'structures' in grid
                    and np.allclose(grid['Tb_K_grid'], STRUCTURE_TB_GRID)):
                log.info(f'Loaded cached Tb grid ({len(grid["structures"])} structures) '
                         f'from {STRUCTURE_CACHE_PATH}')
                mid = grid['structures'][len(grid['structures']) // 2]
                log.info(f'  Mid-grid (Tb={mid["Tb_K"]:.3f} K): '
                         f'D_hsphere={mid["D_hsphere_km"]:.1f} km, '
                         f'CMR2={mid["CMR2"]:.4f}, '
                         f'rho_sil={mid["rhoSil_kgm3"]:.0f} kg/m³')
                return grid
            log.warning('Cached grid schema or Tb array mismatch — rebuilding.')
        except (EOFError, pickle.UnpicklingError, KeyError):
            log.warning('Cache corrupted — rebuilding.')

    log.info(f'Building Tb grid: {len(STRUCTURE_TB_GRID)} structures over '
             f'[{STRUCTURE_TB_GRID[0]:.3f}, {STRUCTURE_TB_GRID[-1]:.3f}] K')
    structures = []
    ref_ci = None
    ref_layer_types = None
    ref_region_phases = None
    for i, Tb in enumerate(STRUCTURE_TB_GRID):
        log.info(f'  [{i+1}/{len(STRUCTURE_TB_GRID)}] Tb = {Tb:.4f} K ...')
        s = _build_single_structure(float(Tb))
        if ref_ci is None:
            ref_ci = s['changeIndices']
            ref_layer_types = s['layer_types']
            ref_region_phases = s['region_phases']
        else:
            if not np.array_equal(s['changeIndices'], ref_ci):
                raise RuntimeError(
                    f'Grid build at Tb={Tb:.3f} produced different changeIndices '
                    f'than the first grid point.'
                )
            if s['layer_types'] != ref_layer_types:
                raise RuntimeError(f'layer_types changed at Tb={Tb:.3f}')
            if s['region_phases'] != ref_region_phases:
                raise RuntimeError(f'region_phases changed at Tb={Tb:.3f}')
        structures.append(s)
        log.info(f'    D_Ih={s["D_iceIh_km"]:.1f}, D_III={s["D_iceIII_km"]:.1f}, '
                 f'D_V={s["D_iceV_km"]:.1f}, D_VI={s["D_iceVI_km"]:.1f} km, '
                 f'CMR2={s["CMR2"]:.4f}')

    grid = {'Tb_K_grid': np.asarray(STRUCTURE_TB_GRID, dtype=np.float64),
            'structures': structures}
    _save_cache(grid, STRUCTURE_CACHE_PATH)
    log.info(f'  Saved grid ({len(structures)} structures) to {STRUCTURE_CACHE_PATH}')
    return grid


# ============================================================
# Plotting (delegates to shared mcmc_plots toolkit)
# ============================================================

def make_plots(result, structure_grid, output_tag='allice_yao2014_andrade'):
    """Generate diagnostic plots from an InferenceResult."""
    import matplotlib
    matplotlib.use('Agg')
    from PlanetProfile.Inference import mcmc_plots as mp

    samples = result.samples
    param_labels = result.param_labels

    # Representative mid-grid structure for headline diagnostic numbers
    if isinstance(structure_grid, dict) and 'structures' in structure_grid:
        rep = structure_grid['structures'][len(structure_grid['structures']) // 2]
    else:
        rep = structure_grid

    out = lambda name: os.path.join(OUTPUT_DIR, f'{output_tag}_{name}.png')
    title_prefix = (f'All-Ice Titan: D_hsphere={rep["D_hsphere_km"]:.0f} km, '
                    f'CMR2={rep["CMR2"]:.4f}')

    # k2 results and heating from the InferenceResult
    k2_arr = result.k2_results  # (n_samples, 2)
    heating_list = result.heating_results or []

    n_eval = len(heating_list)
    heating_idx = result.metadata.get('heating_indices', np.arange(n_eval))
    eval_samples = samples[heating_idx] if n_eval > 0 else samples[:0]

    totals = np.array([sum(h.values()) if h else 1e-30 for h in heating_list])
    safe_tot = np.where(totals > 1e-30, totals, 1e-30)
    f_sil = np.array([h.get('Sil', 0.0) for h in heating_list]) / safe_tot

    # k2 evaluated for heating subset
    k2_eval = k2_arr[heating_idx] if n_eval > 0 and k2_arr is not None else k2_arr

    # Corner plot
    mp.plot_corner(
        samples, param_labels,
        title=title_prefix + ' — Posterior',
        output_path=out('corner'),
    )

    # k2 scatter
    if k2_eval is not None and len(k2_eval):
        k2_plot = np.column_stack([k2_eval[:, 0], np.abs(k2_eval[:, 1])])
        mp.plot_k2_scatter_by(
            k2_plot, color_values=f_sil,
            colorbar_label='Silicate heating fraction',
            obs_re=OBS['Re_k2'], obs_im=OBS['Im_k2'],
            obs_re_err=OBS['Re_k2_err'], obs_im_err=OBS['Im_k2_err'],
            title=title_prefix + r' — $k_2$ Posterior',
            output_path=out('k2_scatter'),
            cmap='RdYlBu_r', vmin=0, vmax=1,
        )

    # Heating vs parameters
    if n_eval > 0:
        mp.plot_heating_vs_parameters(
            eval_samples, heating_list, param_labels,
            extra_xvals=[], extra_xlabels=[],
            output_path=out('heating'),
            cumulative_bar=True,
            eval_d_ocean=np.zeros(n_eval),
            title=(title_prefix + ' — Tidal Heating vs Parameters (top) '
                   'and Per-Model Fractions (bottom)'),
        )

    # Summary statistics
    log.info(f'=== Posterior Summary ({output_tag}) ===')
    log.info(f'  Representative structure (Tb={rep["Tb_K"]:.3f} K): '
             f'D_hsphere={rep["D_hsphere_km"]:.1f} km, '
             f'CMR2={rep["CMR2"]:.4f}, rho_sil={rep["rhoSil_kgm3"]:.0f} kg/m3')
    stats = result.get_summary_stats()
    for i, name in enumerate(result.param_names):
        log.info(f'  {name}: {stats["median"][i]:.3f} '
                 f'[{stats["q16"][i]:.3f}, {stats["q84"][i]:.3f}]')


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='PPTest50 Andrade all-ice Titan MCMC (Yao 2014 convection)')
    parser.add_argument('--config', type=str, default=DEFAULT_CONFIG_JSON,
                        help='Path to InferenceConfig JSON (default: 8D config)')
    parser.add_argument('--rebuild-structure', action='store_true',
                        help='Force rebuild of structural cache')
    parser.add_argument('--build-structure', action='store_true',
                        help='Build structure then exit without running MCMC')
    parser.add_argument('--no-plots', action='store_true',
                        help='Skip plot generation')
    parser.add_argument('--replot', action='store_true',
                        help='Load existing InferenceResult pkl and regenerate plots')
    args = parser.parse_args()

    # Build / load Tb grid of all-ice structures
    structure_grid = build_or_load_structure_grid(force_rebuild=args.rebuild_structure)

    if args.build_structure:
        log.info('--build-structure: exiting after structure build.')
        sys.exit(0)

    # Derive output pkl path from config filename
    config_stem = os.path.splitext(os.path.basename(args.config))[0]
    pkl_path = os.path.join(OUTPUT_DIR, f'{config_stem}_result.pkl')
    output_tag = config_stem  # used for plot filenames

    # Load InferenceConfig from JSON
    from PlanetProfile.Inference.inference_core import InferenceConfig
    from PlanetProfile.Inference.mcmc_runner import MCMCRunner

    config = InferenceConfig.from_json(args.config)
    # Resolve structure_cache_path relative to repo root if not absolute
    if not os.path.isabs(config.structure_cache_path):
        config.structure_cache_path = os.path.join(_pp_root, config.structure_cache_path)

    if args.replot:
        log.info(f'Loading results from {pkl_path}')
        from PlanetProfile.Inference.inference_core import InferenceResult
        result = InferenceResult.load(pkl_path)
    else:
        runner = MCMCRunner(config)
        result = runner.run()
        result.save(pkl_path)
        log.info(f'MCMC results saved to {pkl_path}')

    if not args.no_plots:
        make_plots(result, structure_grid, output_tag=output_tag)

    log.info('Done.')
