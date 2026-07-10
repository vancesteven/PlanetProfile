"""
MCMC exploration: Maxwell rheology, hybrid hydrosphere-thickness Titan (PPTest45).

Extends PPTest42 (Maxwell ocean) by treating total hydrosphere thickness D_hydro_km
as a free parameter instead of letting PlanetProfile close on MoI.  The hydrosphere
PT remains self-consistent (driven by Tb_K), but the silicate top is placed at
R - D_hydro_km, so Mtot_kg and CMR2 become model outputs.

This allows exploration of interior structure space without requiring PlanetProfile
to satisfy the CMR2 closure condition at each MCMC step.

Parameter space (6D):
  log10(eta_Ih):   Ice Ih viscosity (Pa s)         [12, 16]
  log10(eta_HP):   HP ice viscosity (Pa s)          [10, 18]
  log10(eta_sil):  Silicate viscosity (Pa s)        [12, 22]
  Tb_K:            Basal temperature (K)            [252, 270]
  D_hydro_km:      Total hydrosphere thickness (km) [50, 800]
  rho_sil:         Silicate density (kg/m³)         [2000, 6000]

Observational constraints (Petricca et al. 2025 + Titan bulk properties):
  Re(k2)   = 0.608 +/- 0.048
  |Im(k2)| = 0.135 +/- 0.035
  CMR2     = 0.343 +/- 0.001   (Petricca et al. 2025)
  Mtot derived from rho=1.880±0.004 g/cm³, R=2575.5±2.0 km

Grid pre-computation:
  - Tb_K grid:    [252, 270] K at 0.1 K spacing → 181 points
                  252 K → D_iceIh ~ 157 km; 270 K → D_iceIh ~ 33 km
  - D_hydro_km:   [50, 800] km at 50 km spacing → 16 points
                  D_hydro_km < D_iceIh_km → no-ocean (all-ice) structure
  - Total:        181 × 16 = 2896 grid points (some no-ocean pts may be skipped)
  Estimated build time: ~10 hours (incremental saves; safe to interrupt/resume).
  The hybrid grid is generated via PlanetProfile/Inference/hybrid_structure_cache.py.

Usage:
  conda activate ./venvPP
  python PlanetProfile/Test/Test45_mcmc_maxwell_hybrid_hydro.py

  # To force grid rebuild:
  python PlanetProfile/Test/Test45_mcmc_maxwell_hybrid_hydro.py --rebuild-grid
"""
import argparse
import logging
import os
import pickle
import sys
import time

import numpy as np

# --- Environment setup ---
_pp_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _pp_root)

logging.basicConfig(level=logging.WARNING, format='%(name)s - %(message)s')
log = logging.getLogger('PPTest45_MCMC')
log.setLevel(logging.INFO)

# TidalPy imports
from TidalPy.RadialSolver import build_rs_input_from_data, radial_solver
from TidalPy.rheology import Elastic, Maxwell
from TidalPy.tides.multilayer.heating import (
    calc_radial_volumetric_tidal_heating_from_rs_solution,
)

# ============================================================
# Configuration
# ============================================================
OUTPUT_DIR = os.path.join(_pp_root, 'PlanetProfile', 'Test', 'mcmc_results', 'Titan', 'Test45_maxwell_hybrid_hydro')
os.makedirs(OUTPUT_DIR, exist_ok=True)

GRID_CACHE_PATH = os.path.join(OUTPUT_DIR, 'titan_maxwell_hybrid_hydro_grid_cache.pkl')

# Observational constraints
RE_K2_OBS  = 0.608;  RE_K2_ERR  = 0.048
IM_K2_OBS  = 0.135;  IM_K2_ERR  = 0.035
CMR2_OBS   = 0.343;  CMR2_ERR   = 0.001
# Mass derived from density × volume; σ propagated from σ_ρ and σ_R
TITAN_RHO_GCM3     = 1.880;     TITAN_RHO_ERR_GCM3 = 0.004   # g/cm³
TITAN_R_KM         = 2575.5;    TITAN_R_ERR_KM     = 2.0      # km
TITAN_R_M          = TITAN_R_KM * 1e3
MTOT_OBS = (4.0/3.0) * np.pi * TITAN_R_M**3 * (TITAN_RHO_GCM3 * 1e3)   # ≈ 1.34534e23 kg
MTOT_ERR = MTOT_OBS * np.sqrt((TITAN_RHO_ERR_GCM3/TITAN_RHO_GCM3)**2
                               + (3*TITAN_R_ERR_KM/TITAN_R_KM)**2)        # ≈ 4.25e20 kg
# Reference mass used in CMR2 denominator during grid build (Planet.Bulk.M_kg from PPTest45)
TITAN_M_REF_KG = 1.3452e23

# Grid extents — 0.1 K Tb_K resolution covers D_iceIh ~157 km (252 K) to ~33 km (270 K)
TB_MIN, TB_MAX, TB_STEP     = 252.0, 270.0, 0.1    # 181 Tb_K pts; ~0.6 km D_iceIh resolution
D_MIN,  D_MAX,  D_STEP      = 50.0,  800.0, 50.0   # 16 D_hydro pts; 50-160 km = no-ocean

# MCMC settings
N_EFF        = 500
RANDOM_STATE = 42
N_REEVAL     = 500

PARAM_NAMES  = ['log10_eta_Ih', 'log10_eta_HP', 'log10_eta_sil', 'Tb_K', 'D_hydro_km', 'rho_sil']
PARAM_LABELS = [
    r'$\log_{10}\eta_\mathrm{Ih}$',
    r'$\log_{10}\eta_\mathrm{HP}$',
    r'$\log_{10}\eta_\mathrm{sil}$',
    r'$T_b$ (K)',
    r'$D_\mathrm{hydro}$ (km)',
    r'$\rho_\mathrm{sil}$ (kg m$^{-3}$)',
]
N_DIM = 6


# ============================================================
# Step 1: Build / load hybrid grid
# ============================================================

def build_or_load_grid(force_rebuild: bool = False):
    """Build the Tb_K × D_hydro_km hybrid grid or load from cache.

    Incremental resume: if the cache exists but is incomplete, missing points
    are computed and appended automatically.
    """
    tb_grid = list(np.arange(TB_MIN, TB_MAX + TB_STEP / 2, TB_STEP))
    d_grid  = list(np.arange(D_MIN,  D_MAX  + D_STEP  / 2, D_STEP))
    n_expected = len(tb_grid) * len(d_grid)

    from PlanetProfile.Inference.hybrid_structure_cache import build_hybrid_hydrosphere_grid

    if os.path.exists(GRID_CACHE_PATH) and not force_rebuild:
        with open(GRID_CACHE_PATH, 'rb') as f:
            cache_data = pickle.load(f)
        n_pts = len(cache_data.get('grid_cache', {}))
        build_complete = cache_data.get('grid_metadata', {}).get('build_complete', False)
        log.info(f'Loaded {n_pts} cached grid points from {GRID_CACHE_PATH}'
                 f'{" (build complete)" if build_complete else ""}')
        if n_pts < n_expected and not build_complete:
            log.info(f'Grid incomplete ({n_pts}/{n_expected}). Resuming build for missing points...')
            cache_data = build_hybrid_hydrosphere_grid(
                'PlanetProfile.Test.PPTest45',
                tb_grid,
                d_grid,
                GRID_CACHE_PATH,
                rheology='maxwell',
                force_rebuild=False,
            )
    else:
        log.info(f'Building {len(tb_grid)} × {len(d_grid)} hybrid hydrosphere grid '
                 f'({n_expected} points)...')
        cache_data = build_hybrid_hydrosphere_grid(
            'PlanetProfile.Test.PPTest45',
            tb_grid,
            d_grid,
            GRID_CACHE_PATH,
            rheology='maxwell',
            force_rebuild=force_rebuild,
        )

    grid_cache = cache_data['grid_cache']
    tb_vals = np.array(sorted(set(k[0] for k in grid_cache)))
    d_vals  = np.array(sorted(set(k[1] for k in grid_cache)))
    log.info(f'Grid: Tb_K [{tb_vals[0]:.1f}–{tb_vals[-1]:.1f}] K  '
             f'D_hydro [{d_vals[0]:.0f}–{d_vals[-1]:.0f}] km  '
             f'({len(tb_vals)}×{len(d_vals)} = {len(grid_cache)} points)')
    return grid_cache, tb_vals, d_vals


# ============================================================
# Step 2: Forward model (2D nearest-neighbour lookup)
# ============================================================

def forward_model(theta, grid_cache, tb_vals, d_vals, return_heating=False):
    """Compute k2 and optionally per-phase heating for hybrid structure.

    Args:
        theta: [log10_eta_Ih, log10_eta_HP, log10_eta_sil, Tb_K, D_hydro_km]

    Returns:
        (Re_k2, Im_k2, Mtot_kg, CMR2, perPhase_W)
        perPhase_W is {} when return_heating=False or eccentricity=0.
    """
    log10_eta_Ih, log10_eta_HP, log10_eta_sil, Tb_K, D_hydro_km, rho_sil = theta
    eta_Ih  = 10 ** log10_eta_Ih
    eta_HP  = 10 ** log10_eta_HP
    eta_sil = 10 ** log10_eta_sil

    # 2D nearest-neighbour lookup (some points may be absent if grid build skipped them)
    idx_tb = int(np.argmin(np.abs(tb_vals - Tb_K)))
    idx_d  = int(np.argmin(np.abs(d_vals  - D_hydro_km)))
    data = grid_cache.get((float(tb_vals[idx_tb]), float(d_vals[idx_d])))
    if data is None:
        return np.nan, np.nan, np.nan, np.nan, {}

    # Analytically recompute Mtot and CMR2 for sampled rho_sil.
    # Grid stores values with rho_sil_cached (PPTest42 MoI-closed density, ~3595 kg/m³).
    # Decompose into hydrosphere + silicate sphere, substitute new rho_sil.
    r_sil_top    = data['R_body_m'] - data['D_hydro_km'] * 1e3
    rho_sil_c    = data['rhoSil_kgm3']
    V_sil        = (4.0/3.0) * np.pi * r_sil_top**3
    I_sil_factor = (8.0*np.pi/15.0) * r_sil_top**5
    M_hydro      = data['Mtot_kg'] - rho_sil_c * V_sil
    I_hydro      = data['CMR2'] * TITAN_M_REF_KG * data['R_body_m']**2 - rho_sil_c * I_sil_factor
    Mtot_kg      = M_hydro + rho_sil * V_sil
    I_total      = I_hydro + rho_sil * I_sil_factor
    CMR2         = I_total / (Mtot_kg * data['R_body_m']**2)

    # Apply viscosity overrides
    eta_mod   = data['eta_Pa_base'].copy()
    phases    = data['phases']
    ci        = data['changeIndices']
    n_layers  = data['n_layers']

    for i in range(n_layers):
        s, e = ci[i], ci[i + 1]
        ph = int(phases[s])
        if ph == 1:              eta_mod[s:e] = eta_Ih
        elif ph in (3, 5, 6):   eta_mod[s:e] = eta_HP
        elif 50 <= ph < 100:    eta_mod[s:e] = eta_sil

    # Build Maxwell rheology models per layer
    shear = []
    for rp in data['region_phases']:
        base = rp.replace('_conv', '')
        shear.append(Elastic() if base in ('0', 'Clath') else Maxwell())
    bulk = [Elastic() for _ in shear]

    try:
        bd = build_rs_input_from_data(
            data['omega'],
            data['r_m'],
            data['rho'],
            data['K_Pa'],
            data['mu_Pa'],
            data['bulk_visc'],
            np.ascontiguousarray(eta_mod),
            data['layer_upper_radii'],
            data['layer_types'],
            tuple([False] * n_layers),
            tuple([False] * n_layers),
            tuple(shear),
            tuple(bulk),
            perform_checks=False,
            warnings=False,
        )
        result = radial_solver(
            bd.radius_array, bd.density_array,
            bd.complex_bulk_modulus_array, bd.complex_shear_modulus_array,
            bd.frequency, bd.planet_bulk_density,
            bd.layer_types, bd.is_static_bylayer, bd.is_incompressible_bylayer,
            bd.upper_radius_bylayer_array,
            degree_l=2, solve_for=('tidal',), verbose=False, raise_on_fail=False,
        )
        if not result.success:
            return np.nan, np.nan, Mtot_kg, CMR2, {}

        k2    = complex(result.k)
        Re_k2 = k2.real
        Im_k2 = k2.imag

        perPhase_W = {}
        if return_heating and data['eccentricity'] > 0:
            perPhase_W = _compute_heating(result, data)

        return Re_k2, Im_k2, Mtot_kg, CMR2, perPhase_W

    except Exception as exc:
        log.debug(f'TidalPy failed: {exc}')
        return np.nan, np.nan, Mtot_kg, CMR2, {}


def _compute_heating(result, data):
    from PlanetProfile.Utilities.Indexing import PhaseConv
    hp = calc_radial_volumetric_tidal_heating_from_rs_solution(
        data['eccentricity'], data['omega'], data['a_m'],
        data['host_mass'], result, perform_checks=False,
    )
    rr    = np.asarray(result.radius_array)
    r_m   = data['r_m']
    ci    = data['changeIndices']
    n     = data['n_layers']
    pmap  = data.get('phase_map', {0: '0', 1: 'Ih', 2: 'II', 3: 'III', 5: 'V', 6: 'VI'})
    out   = {}
    for i in range(n):
        s, e  = ci[i], ci[i + 1]
        ph    = int(data['phases'][s])
        pname = pmap.get(ph, PhaseConv(ph, liq='0'))
        lo, hi = r_m[s], r_m[e - 1]
        mask  = (rr >= lo - 1.0) & (rr <= hi + 1.0)
        if np.any(mask):
            lr, lh = rr[mask], hp[mask]
            power  = (np.trapezoid(lh * 4 * np.pi * lr**2, lr)
                      if len(lr) > 1 else lh[0] * (4/3) * np.pi * (hi**3 - lo**3))
            out[pname] = out.get(pname, 0.0) + power
    return out


# ============================================================
# Step 3: Log-likelihood
# ============================================================

def log_likelihood(theta, grid_cache, tb_vals, d_vals):
    """Gaussian log-likelihood on k2, CMR2, and Titan mass."""
    Re_k2, Im_k2, Mtot_kg, CMR2, _ = forward_model(
        theta, grid_cache, tb_vals, d_vals, return_heating=False
    )
    if np.isnan(Re_k2):
        return -1e30
    chi2 = (
        ((Re_k2 - RE_K2_OBS)     / RE_K2_ERR) ** 2 +
        ((abs(Im_k2) - IM_K2_OBS) / IM_K2_ERR) ** 2
    )
    if np.isfinite(CMR2):
        chi2 += ((CMR2 - CMR2_OBS) / CMR2_ERR) ** 2
    if np.isfinite(Mtot_kg):
        chi2 += ((Mtot_kg - MTOT_OBS) / MTOT_ERR) ** 2
    return -0.5 * chi2


# ============================================================
# Step 4: MCMC
# ============================================================

def run_mcmc(grid_cache, tb_vals, d_vals):
    import pocomc as pc
    from scipy.stats import uniform

    log.info(f'Starting pocoMC MCMC ({N_DIM}D, n_eff={N_EFF})')

    prior = pc.Prior([
        uniform(loc=12.0, scale=4.0),                   # log10_eta_Ih: [12, 16]
        uniform(loc=10.0, scale=8.0),                   # log10_eta_HP: [10, 18]
        uniform(loc=12.0, scale=10.0),                  # log10_eta_sil: [12, 22]
        uniform(loc=TB_MIN, scale=TB_MAX - TB_MIN),     # Tb_K: [252, 270] K
        uniform(loc=D_MIN,  scale=D_MAX  - D_MIN),      # D_hydro_km: [50, 800] km
        uniform(loc=2000.0, scale=4000.0),               # rho_sil: [2000, 6000] kg/m³
    ])

    def _log_like(theta):
        return log_likelihood(theta, grid_cache, tb_vals, d_vals)

    t0 = time.time()
    sampler = pc.Sampler(
        prior=prior,
        likelihood=_log_like,
        n_effective=N_EFF,
        random_state=RANDOM_STATE,
    )
    sampler.run()
    elapsed = time.time() - t0

    samples, log_likes, logp, logw = sampler.posterior()
    log.info(f'MCMC done in {elapsed/60:.1f} min, {len(samples)} samples')
    log.info(f'  Best log-like: {np.max(log_likes):.2f}')
    return samples, log_likes, sampler


# ============================================================
# Step 5: Re-evaluate heating on posterior
# ============================================================

def evaluate_heating(samples, grid_cache, tb_vals, d_vals, n_eval=None):
    if n_eval is None:
        n_eval = min(len(samples), N_REEVAL)
    idx = np.random.choice(len(samples), n_eval, replace=False)
    idx.sort()

    k2_results      = []
    mtot_results    = []
    cmr2_results    = []
    heating_results = []

    log.info(f'Re-evaluating {n_eval} posterior samples...')
    t0 = time.time()
    for i, si in enumerate(idx):
        Re, Im, Mtot, CMR2, phW = forward_model(
            samples[si], grid_cache, tb_vals, d_vals, return_heating=True
        )
        k2_results.append((Re, Im))
        mtot_results.append(Mtot)
        cmr2_results.append(CMR2)
        heating_results.append(phW)
        if (i + 1) % 100 == 0:
            log.info(f'  {i+1}/{n_eval} ({time.time()-t0:.0f}s)')

    log.info(f'Done in {time.time()-t0:.0f}s')
    return idx, np.array(k2_results), np.array(mtot_results), np.array(cmr2_results), heating_results


# ============================================================
# Step 6: Plots
# ============================================================

def make_plots(samples, log_likes, k2_results, mtot_results, cmr2_results,
               heating_results, eval_idx, grid_cache, tb_vals, d_vals):
    import corner
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.patches import Ellipse

    sns.set_theme(style='white', font_scale=0.9)

    eval_samples = samples[eval_idx]
    Re_arr = k2_results[:, 0]
    Im_arr = np.abs(k2_results[:, 1])
    D_arr  = eval_samples[:, 4]  # D_hydro_km

    # Build Tb_K → D_iceIh_km lookup (use max across D_hydro values to get true ice thickness;
    # truncated no-ocean entries have D_iceIh_km = D_hydro_km < true thickness)
    tb_to_diceIh = {}
    for (tb, _d), struct in grid_cache.items():
        val = struct.get('D_iceIh_km', np.nan)
        if np.isfinite(val) and (tb not in tb_to_diceIh or val > tb_to_diceIh[tb]):
            tb_to_diceIh[tb] = val

    def _tb_to_diceIh(tb_k):
        idx = int(np.argmin(np.abs(tb_vals - tb_k)))
        return tb_to_diceIh.get(float(tb_vals[idx]), np.nan)

    d_iceIh_all  = np.array([_tb_to_diceIh(theta[3]) for theta in samples])
    d_iceIh_eval = d_iceIh_all[eval_idx]

    # Silicate heating fraction per eval sample
    f_sil = np.array([
        h.get('Sil', 0) / max(sum(h.values()), 1e-30)
        for h in heating_results
    ])

    # --- 1. Corner plot: D_iceIh in place of Tb_K ---
    corner_samples = np.column_stack([
        samples[:, :3], d_iceIh_all, samples[:, 4], samples[:, 5],
    ])
    corner_labels = (list(PARAM_LABELS[:3])
                     + [r'$D_\mathrm{IceIh}$ (km)']
                     + list(PARAM_LABELS[4:]))
    # Compute per-column range; expand degenerate columns by ±1 to avoid corner crash
    cs_std = np.nanstd(corner_samples, axis=0)
    cs_med = np.nanmedian(corner_samples, axis=0)
    corner_range = [
        (float(np.nanmin(corner_samples[:, i])) - 0.1 * abs(cs_std[i]),
         float(np.nanmax(corner_samples[:, i])) + 0.1 * abs(cs_std[i]))
        if cs_std[i] > 0
        else (float(cs_med[i]) - 1.0, float(cs_med[i]) + 1.0)
        for i in range(corner_samples.shape[1])
    ]
    fig = corner.corner(
        corner_samples, labels=corner_labels,
        color='steelblue', range=corner_range,
        quantiles=[0.16, 0.5, 0.84], show_titles=True,
        title_fmt='.3f', title_kwargs={'fontsize': 10},
    )
    fig.suptitle('Hybrid Hydrosphere Titan: Posterior (Maxwell)', fontsize=14, y=1.02)
    _save(fig, 'hybrid_hydro_maxwell_corner.png')

    # --- 2. k2 scatter coloured by D_hydro_km (1σ + 2σ ellipses) ---
    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(Re_arr, Im_arr, c=D_arr, cmap='plasma_r', s=8, alpha=0.6,
                    vmin=D_MIN, vmax=D_MAX)
    plt.colorbar(sc, ax=ax, label=r'$D_\mathrm{hydro}$ (km)')
    for nσ, ls, lw, lbl in [(1, '--', 2, r'Petricca 1$\sigma$'), (2, ':', 1, r'Petricca 2$\sigma$')]:
        ax.add_patch(Ellipse((RE_K2_OBS, IM_K2_OBS),
                             2*nσ*RE_K2_ERR, 2*nσ*IM_K2_ERR,
                             fill=False, color='red', linewidth=lw, linestyle=ls, label=lbl))
    ax.set_xlabel(r'$\mathrm{Re}(k_2)$');  ax.set_ylabel(r'$|\mathrm{Im}(k_2)|$')
    ax.set_title(r'Hybrid Hydrosphere: $k_2$ Posterior (by $D_\mathrm{hydro}$)')
    ax.legend()
    _save(fig, 'hybrid_hydro_maxwell_k2_scatter.png')

    # --- 3. k2 scatter coloured by silicate heating fraction (Test42-style) ---
    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(Re_arr, Im_arr, c=f_sil, cmap='RdYlBu_r', s=8, alpha=0.6, vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax, label='Silicate heating fraction')
    for nσ, ls, lw, lbl in [(1, '--', 2, r'Petricca 1$\sigma$'), (2, ':', 1, r'Petricca 2$\sigma$')]:
        ax.add_patch(Ellipse((RE_K2_OBS, IM_K2_OBS),
                             2*nσ*RE_K2_ERR, 2*nσ*IM_K2_ERR,
                             fill=False, color='red', linewidth=lw, linestyle=ls, label=lbl))
    ax.set_xlabel(r'$\mathrm{Re}(k_2)$');  ax.set_ylabel(r'$|\mathrm{Im}(k_2)|$')
    ax.set_title(r'Hybrid Hydrosphere: $k_2$ Posterior (silicate heating)')
    ax.legend()
    _save(fig, 'hybrid_hydro_maxwell_k2_scatter_heating.png')

    # --- 4. Heating distribution vs parameters (log10 power, 2×3 for 6 params) ---
    # Heating fractions are bimodal (0 or 1) because the posterior favours thin ice shells
    # (Tb_K → 269 K → no HP phases) and high silicate viscosity (near-elastic silicate).
    # Log10 power reveals the actual heating magnitudes and their parameter dependence.
    ALL_PHASES = ['Ih', 'III', 'V', 'VI', 'Sil']
    phase_colors = {'Ih': 'C0', 'III': 'C1', 'V': 'C2', 'VI': 'C3', 'Sil': 'C4'}
    heating_power = {
        ph: np.array([h.get(ph, 0.0) for h in heating_results])
        for ph in ALL_PHASES
    }
    # Only show phases that carry non-trivial power in at least a few samples
    active_phases = [ph for ph in ALL_PHASES
                     if np.sum(heating_power[ph] > 1e6) > 5]
    W_FLOOR = 1e8   # plot floor: 1e8 W (below this → not shown)

    plot_xvals   = [eval_samples[:, 0], eval_samples[:, 1], eval_samples[:, 2],
                    d_iceIh_eval, eval_samples[:, 4], eval_samples[:, 5]]
    plot_xlabels = (list(PARAM_LABELS[:3])
                    + [r'$D_\mathrm{IceIh}$ (km)']
                    + list(PARAM_LABELS[4:]))

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for ip, (xvals, xlabel) in enumerate(zip(plot_xvals, plot_xlabels)):
        ax = axes.flat[ip]
        for ph in active_phases:
            pw = heating_power[ph]
            mask = pw > W_FLOOR
            if np.any(mask):
                ax.scatter(xvals[mask], np.log10(pw[mask]), s=4, alpha=0.4,
                           color=phase_colors[ph], label=ph)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(r'$\log_{10}$ Heating (W)')
        if ip == 0 and active_phases:
            ax.legend(fontsize=8, loc='best')
    fig.suptitle('Hybrid Hydrosphere: Tidal Heating Power', fontsize=14)
    fig.tight_layout()
    _save(fig, 'hybrid_hydro_maxwell_heating.png')

    # --- 5. Mtot and CMR2 vs D_hydro_km ---
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax = axes[0]
    ax.scatter(D_arr, mtot_results / 1e23, s=6, alpha=0.5, c='steelblue')
    ax.axhline(MTOT_OBS / 1e23, color='red', ls='--',
               label=f'obs ± {MTOT_ERR/1e20:.1f}×10²⁰ kg')
    ax.axhspan((MTOT_OBS - MTOT_ERR) / 1e23, (MTOT_OBS + MTOT_ERR) / 1e23,
               alpha=0.15, color='red')
    ax.set_xlabel(r'$D_\mathrm{hydro}$ (km)');  ax.set_ylabel(r'$M_\mathrm{tot}$ (×10²³ kg)')
    ax.set_title('Total Mass vs Hydrosphere Thickness');  ax.legend()
    ax = axes[1]
    ax.scatter(D_arr, cmr2_results, s=6, alpha=0.5, c='darkorange')
    ax.axhline(CMR2_OBS, color='red', ls='--', label=f'obs={CMR2_OBS}±{CMR2_ERR}')
    ax.axhspan(CMR2_OBS - CMR2_ERR, CMR2_OBS + CMR2_ERR, alpha=0.15, color='red')
    ax.set_xlabel(r'$D_\mathrm{hydro}$ (km)');  ax.set_ylabel('C/MR²')
    ax.set_title('CMR2 vs Hydrosphere Thickness');  ax.legend()
    fig.tight_layout()
    _save(fig, 'hybrid_hydro_maxwell_mass_cmr2.png')

    # --- 6. CMR2 surface over grid (Tb_K × D_hydro_km) ---
    fig, ax = plt.subplots(figsize=(8, 6))
    CMR2_grid = np.full((len(d_vals), len(tb_vals)), np.nan)
    for j, d in enumerate(d_vals):
        for i, tb in enumerate(tb_vals):
            key = (float(tb), float(d))
            if key in grid_cache:
                CMR2_grid[j, i] = grid_cache[key].get('CMR2', np.nan)
    im = ax.pcolormesh(tb_vals, d_vals, CMR2_grid, cmap='viridis', shading='nearest')
    plt.colorbar(im, ax=ax, label='C/MR²')
    ax.set_xlabel(r'$T_b$ (K)');  ax.set_ylabel(r'$D_\mathrm{hydro}$ (km)')
    ax.set_title('CMR2 across the grid')
    _save(fig, 'hybrid_hydro_maxwell_cmr2_surface.png')

    # --- 7. Ice shell thickness vs Tb_K from grid (Test42-style), posterior overlay ---
    fig, ax = plt.subplots(figsize=(8, 5))
    tb_sorted   = sorted(tb_to_diceIh)
    diceIh_line = [tb_to_diceIh[tb] for tb in tb_sorted]
    if any(np.isfinite(v) for v in diceIh_line):
        ax.plot(tb_sorted, diceIh_line, 'k-o', markersize=4)
        ax2 = ax.twinx()
        ax2.hist(samples[:, 3], bins=30, alpha=0.3, color='blue',
                 density=True, label=r'Posterior $T_b$')
        ax2.set_ylabel('Posterior density', color='blue')
        ax.set_xlabel(r'$T_b$ (K)')
        ax.set_ylabel(r'Ice Ih shell thickness (km)')
        ax.set_title('Ice Shell Thickness vs Basal Temperature')
        fig.tight_layout()
    _save(fig, 'hybrid_hydro_maxwell_Tb_structure.png')

    plt.close('all')


def _save(fig, fname):
    import matplotlib.pyplot as plt
    path = os.path.join(OUTPUT_DIR, fname)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    log.info(f'Saved {path}')
    plt.close(fig)


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PPTest45 hybrid MCMC')
    parser.add_argument('--rebuild-grid', action='store_true',
                        help='Force rebuild of hybrid structure grid from scratch')
    parser.add_argument('--grid-only', action='store_true',
                        help='Build/resume grid then exit without running MCMC')
    parser.add_argument('--no-plots', action='store_true',
                        help='Skip plot generation')
    parser.add_argument('--replot', action='store_true',
                        help='Load existing pkl and regenerate plots without re-running MCMC')
    args = parser.parse_args()

    # Build / load grid
    grid_cache, tb_vals, d_vals = build_or_load_grid(force_rebuild=args.rebuild_grid)

    if args.grid_only:
        log.info('--grid-only: grid complete, exiting.')
        sys.exit(0)

    pkl_path = os.path.join(OUTPUT_DIR, 'hybrid_hydro_maxwell_mcmc.pkl')

    if args.replot:
        log.info(f'Replot mode: loading existing results from {pkl_path}')
        with open(pkl_path, 'rb') as f:
            results = pickle.load(f)
        samples       = results['samples']
        log_likes     = results['log_likes']
        eval_idx      = results['eval_idx']
        k2_results    = results['k2_results']
        mtot_results  = results['mtot_results']
        cmr2_results  = results['cmr2_results']
        heating_results = results['heating_results']
        if not args.no_plots:
            make_plots(samples, log_likes, k2_results, mtot_results, cmr2_results,
                       heating_results, eval_idx, grid_cache, tb_vals, d_vals)
        sys.exit(0)

    # Sanity check forward model
    theta_test = [14.0, 13.0, 15.0, 260.0, 500.0, 3300.0]
    Re, Im, Mtot, CMR2, _ = forward_model(theta_test, grid_cache, tb_vals, d_vals)
    log.info(f'Sanity check: Re(k2)={Re:.4f}, Im(k2)={Im:.4f}, '
             f'Mtot={Mtot:.4e} kg, CMR2={CMR2:.5f}')
    if np.isnan(Re):
        log.error('Forward model returned NaN — check grid or TidalPy installation.')
        sys.exit(1)

    # Run MCMC
    samples, log_likes, sampler = run_mcmc(grid_cache, tb_vals, d_vals)

    # Re-evaluate heating on posterior subset
    eval_idx, k2_results, mtot_results, cmr2_results, heating_results = evaluate_heating(
        samples, grid_cache, tb_vals, d_vals
    )

    # Save results
    results = {
        'samples':          samples,
        'log_likes':        log_likes,
        'param_names':      PARAM_NAMES,
        'param_labels':     PARAM_LABELS,
        'eval_idx':         eval_idx,
        'k2_results':       k2_results,
        'mtot_results':     mtot_results,
        'cmr2_results':     cmr2_results,
        'heating_results':  heating_results,
        'tb_vals':          tb_vals,
        'd_vals':           d_vals,
        'observational_constraints': {
            'Re_k2':     (RE_K2_OBS,  RE_K2_ERR),
            'abs_Im_k2': (IM_K2_OBS,  IM_K2_ERR),
            'CMR2':      (CMR2_OBS,   CMR2_ERR),
            'Mtot_kg':   (MTOT_OBS,   MTOT_ERR),
        },
    }
    pkl_path = os.path.join(OUTPUT_DIR, 'hybrid_hydro_maxwell_mcmc.pkl')
    with open(pkl_path, 'wb') as f:
        pickle.dump(results, f)
    log.info(f'Results saved to {pkl_path}')

    if not args.no_plots:
        make_plots(samples, log_likes, k2_results, mtot_results, cmr2_results,
                   heating_results, eval_idx, grid_cache, tb_vals, d_vals)

    # Summary
    print('\n' + '=' * 65)
    print('HYBRID HYDROSPHERE TITAN MCMC SUMMARY')
    print('=' * 65)
    for ip, name in enumerate(PARAM_NAMES):
        med = np.median(samples[:, ip])
        lo  = np.percentile(samples[:, ip], 16)
        hi  = np.percentile(samples[:, ip], 84)
        print(f'  {name:20s}: {med:.3f}  [{lo:.3f}, {hi:.3f}]')

    eval_s = samples[eval_idx]
    print(f'\n  Re(k2) median:    {np.nanmedian(k2_results[:,0]):.4f}  '
          f'(obs {RE_K2_OBS:.3f} ± {RE_K2_ERR:.3f})')
    print(f'  |Im(k2)| median:  {np.nanmedian(np.abs(k2_results[:,1])):.4f}  '
          f'(obs {IM_K2_OBS:.3f} ± {IM_K2_ERR:.3f})')
    print(f'  CMR2 median:      {np.nanmedian(cmr2_results):.5f}  '
          f'(obs {CMR2_OBS:.3f} ± {CMR2_ERR:.3f})')
    print(f'  Mtot median:      {np.nanmedian(mtot_results):.4e} kg  '
          f'(obs {MTOT_OBS:.4e} ± {MTOT_ERR:.1e})')
    print(f'\n  Output: {OUTPUT_DIR}')
