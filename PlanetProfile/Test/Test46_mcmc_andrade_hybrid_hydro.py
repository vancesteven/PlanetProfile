"""
MCMC exploration: Andrade rheology, hybrid hydrosphere-thickness Titan (PPTest46).

Extends PPTest42 (Maxwell ocean) by treating total hydrosphere thickness D_hydro_km
as a free parameter instead of letting PlanetProfile close on MoI.  The hydrosphere
PT remains self-consistent (driven by Tb_K), but the silicate top is placed at
R - D_hydro_km, so Mtot_kg and CMR2 become model outputs.

This allows exploration of interior structure space without requiring PlanetProfile
to satisfy the CMR2 closure condition at each MCMC step.

Parameter space (10D):
  alpha:           Andrade alpha exponent                [0.2, 0.4]
  log10(zeta):     Andrade zeta (Pa s^alpha)             [-2, 2]
  log10(eta_Ih):   Ice Ih viscosity (Pa s)               [12, 16]
  log10(eta_HP):   HP ice viscosity (Pa s)               [10, 18]
  log10(eta_sil):  Silicate viscosity (Pa s)             [18, 22]
  Tb_K:            Basal temperature (K)                 [252, 270]
  D_hydro_km:      Total hydrosphere thickness (km)      [50, 800]
  rho_sil:         Outer silicate mantle density (kg/m³) [2000, 5000]
  rho_core:        Dense inner silicate/rock density     [2000, 5000]
  f_core:          R_core / r_sil_top (radius fraction)  [0, 0.8]
  Constraint: rho_core >= rho_sil (stable stratification; equality = single layer)

Observational constraints (Petricca et al. 2025, Cassini radio science):
  Re(k2)   = 0.608 +/- 0.048
  |Im(k2)| = 0.135 +/- 0.035
  CMR2     = 0.343 +/- 0.001
  Mtot derived from rho=1.880+/-0.004 g/cm3, R=2575.5+/-2.0 km

Grid pre-computation:
  - Tb_K grid:    [252, 253) K at 0.02 K (fine) + [253, 270] at 0.1 K → ~221 points
                  Fine spacing near cold boundary where posterior lives.
                  252 K → D_iceIh ~ 157 km; 270 K → D_iceIh ~ 33 km
  - D_hydro_km:   [50, 800] km at 25 km spacing → 31 points
                  D_hydro_km < D_iceIh_km → no-ocean (all-ice) structure
  - Total:        ~221 × 31 = 6851 grid points (some no-ocean pts may be skipped)
  Estimated build time: ~10 hours (incremental saves; safe to interrupt/resume).
  The hybrid grid is generated via PlanetProfile/Inference/hybrid_structure_cache.py.
  Grid cache is shared with PPTest45 (maxwell).

Usage:
  conda activate ./venvPP
  python PlanetProfile/Test/Test46_mcmc_andrade_hybrid_hydro.py

  # To force grid rebuild:
  python PlanetProfile/Test/Test46_mcmc_andrade_hybrid_hydro.py --rebuild-grid
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
log = logging.getLogger('PPTest46_MCMC')
log.setLevel(logging.INFO)

# TidalPy imports
from TidalPy.RadialSolver import build_rs_input_from_data, radial_solver
from TidalPy.rheology import Andrade, Elastic
from TidalPy.tides.multilayer.heating import (
    calc_radial_volumetric_tidal_heating_from_rs_solution,
)

# ============================================================
# Configuration
# ============================================================
OUTPUT_DIR = os.path.join(_pp_root, 'PlanetProfile', 'Test', 'mcmc_results', 'Titan', 'Test46_andrade_hybrid_hydro')
os.makedirs(OUTPUT_DIR, exist_ok=True)

GRID_CACHE_PATH = os.path.join(OUTPUT_DIR, 'titan_maxwell_hybrid_hydro_grid_cache.pkl')

# Observational constraints (Petricca et al. 2025, Cassini radio science)
OBS = {
    'Re_k2': 0.608,  'Re_k2_err': 0.048,
    'Im_k2': 0.135,  'Im_k2_err': 0.035,
    'CMR2':  0.343,  'CMR2_err':  0.001,
}
# Mass derived from density × volume; σ propagated from σ_ρ and σ_R
TITAN_RHO_GCM3     = 1.880;     TITAN_RHO_ERR_GCM3 = 0.004   # g/cm³
TITAN_R_KM         = 2575.5;    TITAN_R_ERR_KM     = 2.0      # km
TITAN_R_M          = TITAN_R_KM * 1e3
MTOT_OBS = (4.0/3.0) * np.pi * TITAN_R_M**3 * (TITAN_RHO_GCM3 * 1e3)   # ≈ 1.34534e23 kg
MTOT_ERR = MTOT_OBS * np.sqrt((TITAN_RHO_ERR_GCM3/TITAN_RHO_GCM3)**2
                               + (3*TITAN_R_ERR_KM/TITAN_R_KM)**2)        # ≈ 4.25e20 kg
# Reference mass used in CMR2 denominator during grid build (Planet.Bulk.M_kg from PPTest46)
TITAN_M_REF_KG = 1.3452e23

# Grid extents — fine spacing near cold boundary where posterior lives
TB_MIN, TB_MAX              = 252.0, 270.0
TB_FINE_STEP, TB_FINE_MAX   = 0.02,  253.0   # 50 pts in [252, 253); D_iceIh ~0.1 km resolution
TB_COARSE_STEP              = 0.1             # 171 pts in [253, 270]; ~0.6 km resolution
D_MIN,  D_MAX,  D_STEP      = 50.0,  800.0, 25.0   # 31 D_hydro pts; 50-160 km = no-ocean

# MCMC settings
N_EFF        = 500
RANDOM_STATE = 42
N_REEVAL     = 500

PARAM_NAMES  = ['alpha', 'log10_zeta', 'log10_eta_Ih', 'log10_eta_HP', 'log10_eta_sil',
                'Tb_K', 'D_hydro_km', 'rho_sil', 'rho_core', 'f_core']
PARAM_LABELS = [
    r'$\alpha$',
    r'$\log_{10}\zeta$',
    r'$\log_{10}\eta_\mathrm{Ih}$',
    r'$\log_{10}\eta_\mathrm{HP}$',
    r'$\log_{10}\eta_\mathrm{sil}$',
    r'$T_b$ (K)',
    r'$D_\mathrm{hydro}$ (km)',
    r'$\rho_\mathrm{sil}$ (kg m$^{-3}$)',
    r'$\rho_\mathrm{core}$ (kg m$^{-3}$)',
    r'$f_\mathrm{core}$',
]
N_DIM = 10


# ============================================================
# Step 1: Build / load hybrid grid
# ============================================================

def build_or_load_grid(force_rebuild: bool = False):
    """Build the Tb_K × D_hydro_km hybrid grid or load from cache.

    Incremental resume: if the cache exists but is incomplete, missing points
    are computed and appended automatically.
    """
    tb_fine   = np.arange(TB_MIN, TB_FINE_MAX, TB_FINE_STEP)
    tb_coarse = np.arange(TB_FINE_MAX, TB_MAX + TB_COARSE_STEP / 2, TB_COARSE_STEP)
    tb_grid   = list(np.unique(np.concatenate([tb_fine, tb_coarse])))
    d_grid    = list(np.arange(D_MIN, D_MAX + D_STEP / 2, D_STEP))
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
                'PlanetProfile.Test.PPTest46',
                tb_grid,
                d_grid,
                GRID_CACHE_PATH,
                rheology='maxwell',
                force_rebuild=False,
            )
        elif n_pts < n_expected and build_complete:
            log.info(f'Grid marked complete with {n_pts} points '
                     f'(expected {n_expected}; float rounding in Tb grid). Proceeding.')
    else:
        log.info(f'Building {len(tb_grid)} × {len(d_grid)} hybrid hydrosphere grid '
                 f'({n_expected} points)...')
        cache_data = build_hybrid_hydrosphere_grid(
            'PlanetProfile.Test.PPTest46',
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
        theta: [alpha, log10_zeta, log10_eta_Ih, log10_eta_HP, log10_eta_sil,
                Tb_K, D_hydro_km, rho_sil, rho_core, f_core]
               rho_core >= rho_sil required (stable stratification).
               f_core = R_core / r_sil_top in [0, 0.8].

    Returns:
        (Re_k2, Im_k2, Mtot_kg, CMR2, perPhase_W)
        perPhase_W is {} when return_heating=False or eccentricity=0.
    """
    alpha, log10_zeta, log10_eta_Ih, log10_eta_HP, log10_eta_sil, Tb_K, D_hydro_km, rho_sil, rho_core, f_core = theta
    eta_Ih  = 10 ** log10_eta_Ih
    eta_HP  = 10 ** log10_eta_HP
    eta_sil = 10 ** log10_eta_sil

    # 2D nearest-neighbour lookup (some points may be absent if grid build skipped them)
    idx_tb = int(np.argmin(np.abs(tb_vals - Tb_K)))
    idx_d  = int(np.argmin(np.abs(d_vals  - D_hydro_km)))
    data = grid_cache.get((float(tb_vals[idx_tb]), float(d_vals[idx_d])))
    if data is None:
        return np.nan, np.nan, np.nan, np.nan, {}

    # Analytically recompute Mtot and CMR2 for two-layer interior.
    # Grid stores values with rho_sil_cached (PPTest42 MoI-closed density, ~3595 kg/m³).
    # Decompose into hydrosphere + full silicate sphere, then split silicate into
    # outer mantle (rho_sil) + dense inner core (rho_core, radius R_core = f_core * r_sil_top).
    r_sil_top    = data['R_body_m'] - data['D_hydro_km'] * 1e3
    rho_sil_c    = data['rhoSil_kgm3']
    R_core       = f_core * r_sil_top
    V_sil_full   = (4.0/3.0) * np.pi * r_sil_top**3
    V_core       = (4.0/3.0) * np.pi * R_core**3
    I_sil_full   = (8.0*np.pi/15.0) * r_sil_top**5
    I_core_fac   = (8.0*np.pi/15.0) * R_core**5
    M_hydro      = data['Mtot_kg'] - rho_sil_c * V_sil_full
    I_hydro      = data['CMR2'] * TITAN_M_REF_KG * data['R_body_m']**2 - rho_sil_c * I_sil_full
    Mtot_kg      = M_hydro + rho_sil * (V_sil_full - V_core) + rho_core * V_core
    I_total      = I_hydro + rho_sil * (I_sil_full - I_core_fac) + rho_core * I_core_fac
    CMR2         = I_total / (Mtot_kg * data['R_body_m']**2)

    # No-ocean safeguard: when D_ocean_km == 0 no layer should be 'liquid'.
    layer_types_use   = list(data['layer_types'])
    region_phases_use = list(data['region_phases'])
    if data.get('D_ocean_km', 1.0) < 0.5:
        layer_types_use   = ['solid' if lt == 'liquid' else lt for lt in layer_types_use]
        region_phases_use = ['Ih' if rp == '0' else rp for rp in region_phases_use]

    # Early mass guard: reject wildly unphysical samples before radial solver
    if abs(Mtot_kg - MTOT_OBS) > 5 * MTOT_ERR:
        return np.nan, np.nan, Mtot_kg, CMR2, {}

    # --- A1: Rebuild radial arrays with MCMC-sampled interior densities ---
    # The cached structure has a single silicate layer at rho_sil_cached.
    # Split it into outer mantle (rho_sil) + inner core (rho_core) at R_core.
    r_m       = data['r_m'].copy()
    rho_mod   = data['rho'].copy()
    K_Pa_mod  = data['K_Pa'].copy()
    mu_Pa_mod = data['mu_Pa'].copy()
    bv_mod    = data['bulk_visc'].copy()
    eta_mod   = data['eta_Pa_base'].copy()
    phases    = data['phases'].copy()
    ci        = list(data['changeIndices'])
    n_layers  = data['n_layers']
    layer_upper_radii = list(data['layer_upper_radii'])

    # Find silicate layer(s) and apply two-layer density
    sil_layers = []
    for i in range(n_layers):
        s, e = ci[i], ci[i + 1]
        if 50 <= int(phases[s]) < 100:
            sil_layers.append(i)

    if sil_layers and f_core > 0.001:
        # Split the innermost silicate layer at R_core
        i_sil = sil_layers[0]
        s_sil, e_sil = ci[i_sil], ci[i_sil + 1]
        r_sil_layer = r_m[s_sil:e_sil]

        # Find split point closest to R_core
        idx_split = int(np.searchsorted(r_sil_layer, R_core))
        idx_split = max(2, min(idx_split, len(r_sil_layer) - 2))
        abs_split = s_sil + idx_split

        # Set densities and scale elastic moduli (K ∝ ρ, μ ∝ ρ at constant Vp, Vs)
        rho_cached = rho_mod[s_sil:e_sil].copy()
        rho_mod[s_sil:abs_split] = rho_core
        rho_mod[abs_split:e_sil] = rho_sil
        scale_core = np.where(rho_cached[: idx_split] > 0,
                              rho_core / rho_cached[: idx_split], 1.0)
        scale_mantle = np.where(rho_cached[idx_split:] > 0,
                                rho_sil / rho_cached[idx_split:], 1.0)
        K_Pa_mod[s_sil:abs_split] *= scale_core
        mu_Pa_mod[s_sil:abs_split] *= scale_core
        K_Pa_mod[abs_split:e_sil] *= scale_mantle
        mu_Pa_mod[abs_split:e_sil] *= scale_mantle

        # For any additional silicate layers above the first, use rho_sil
        for i_extra in sil_layers[1:]:
            s_e, e_e = ci[i_extra], ci[i_extra + 1]
            rho_cached_e = rho_mod[s_e:e_e].copy()
            rho_mod[s_e:e_e] = rho_sil
            scale_e = np.where(rho_cached_e > 0, rho_sil / rho_cached_e, 1.0)
            K_Pa_mod[s_e:e_e] *= scale_e
            mu_Pa_mod[s_e:e_e] *= scale_e

        # Insert a new layer boundary at the split point
        R_core_actual = float(r_m[abs_split])
        # Shift all layer indices after the split
        new_ci = ci[:i_sil + 1] + [abs_split] + ci[i_sil + 1:]
        # Insert layer metadata for the new inner core layer
        layer_upper_radii.insert(i_sil, R_core_actual)
        layer_types_use.insert(i_sil, 'solid')
        region_phases_use.insert(i_sil, 'Sil')
        n_layers += 1
        ci = new_ci
    else:
        # f_core ~ 0: single uniform density throughout silicate
        for i_s in sil_layers:
            s_s, e_s = ci[i_s], ci[i_s + 1]
            rho_cached_s = rho_mod[s_s:e_s].copy()
            rho_mod[s_s:e_s] = rho_sil
            scale_s = np.where(rho_cached_s > 0, rho_sil / rho_cached_s, 1.0)
            K_Pa_mod[s_s:e_s] *= scale_s
            mu_Pa_mod[s_s:e_s] *= scale_s

    # Apply viscosity overrides
    for i in range(n_layers):
        s, e = ci[i], ci[i + 1]
        ph = int(phases[min(s, len(phases) - 1)])
        if ph == 1:              eta_mod[s:e] = eta_Ih
        elif ph in (3, 5, 6):   eta_mod[s:e] = eta_HP
        elif 50 <= ph < 100:    eta_mod[s:e] = eta_sil

    # Build Andrade rheology models per layer
    zeta_pa = 10 ** log10_zeta
    zeta_tp = zeta_pa ** (1.0 / alpha)
    shear = []
    for rp in region_phases_use:
        base = rp.replace('_conv', '')
        shear.append(Elastic() if base in ('0', 'Clath') else Andrade(args=(alpha, zeta_tp)))
    bulk = [Elastic() for _ in shear]

    try:
        bd = build_rs_input_from_data(
            data['omega'],
            r_m,
            rho_mod,
            K_Pa_mod,
            mu_Pa_mod,
            bv_mod,
            np.ascontiguousarray(eta_mod),
            tuple(layer_upper_radii),
            tuple(layer_types_use),
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
    rho_sil, rho_core = theta[7], theta[8]
    if rho_core < rho_sil:   # unstable stratification — reject
        return -1e30
    Re_k2, Im_k2, Mtot_kg, CMR2, _ = forward_model(
        theta, grid_cache, tb_vals, d_vals, return_heating=False
    )
    if np.isnan(Re_k2):
        return -1e30
    chi2 = (
        ((Re_k2 - OBS['Re_k2'])     / OBS['Re_k2_err']) ** 2 +
        ((abs(Im_k2) - OBS['Im_k2']) / OBS['Im_k2_err']) ** 2
    )
    if np.isfinite(CMR2):
        chi2 += ((CMR2 - OBS['CMR2']) / OBS['CMR2_err']) ** 2
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
        uniform(loc=0.2,   scale=0.2),                   # alpha: [0.2, 0.4]
        uniform(loc=-2.0,  scale=4.0),                   # log10_zeta: [-2, 2]
        uniform(loc=12.0,  scale=4.0),                   # log10_eta_Ih: [12, 16]
        uniform(loc=10.0,  scale=8.0),                   # log10_eta_HP: [10, 18]
        uniform(loc=18.0,  scale=4.0),                   # log10_eta_sil: [18, 22]
        uniform(loc=TB_MIN, scale=TB_MAX - TB_MIN),      # Tb_K: [252, 270] K
        uniform(loc=D_MIN,  scale=D_MAX  - D_MIN),       # D_hydro_km: [50, 800] km
        uniform(loc=2000.0, scale=3000.0),               # rho_sil (outer): [2000, 5000] kg/m³
        uniform(loc=2000.0, scale=3000.0),               # rho_core (inner): [2000, 5000] kg/m³
        uniform(loc=0.0,    scale=0.8),                  # f_core = R_core/r_sil_top: [0, 0.8]
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
    D_arr  = eval_samples[:, 6]  # D_hydro_km

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

    d_iceIh_all  = np.array([_tb_to_diceIh(theta[5]) for theta in samples])
    d_iceIh_eval = d_iceIh_all[eval_idx]

    # D_ocean for each eval sample (specific to the (Tb_K, D_hydro_km) grid point)
    def _get_docean(tb_k, d_hydro):
        _i_tb = int(np.argmin(np.abs(tb_vals - tb_k)))
        _i_d  = int(np.argmin(np.abs(d_vals  - d_hydro)))
        _pt   = grid_cache.get((float(tb_vals[_i_tb]), float(d_vals[_i_d])))
        return _pt.get('D_ocean_km', 0.0) if _pt is not None else np.nan

    d_ocean_eval = np.array([_get_docean(eval_samples[i, 5], eval_samples[i, 6])
                              for i in range(len(eval_idx))])

    # Silicate heating fraction per eval sample
    f_sil = np.array([
        h.get('Sil', 0) / max(sum(h.values()), 1e-30)
        for h in heating_results
    ])

    # --- 1. Corner plot: D_iceIh in place of Tb_K (all 10 params) ---
    corner_samples = np.column_stack([
        samples[:, :5], d_iceIh_all, samples[:, 6:],
    ])
    corner_labels = (list(PARAM_LABELS[:5])
                     + [r'$D_\mathrm{IceIh}$ (km)']
                     + list(PARAM_LABELS[6:]))
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
    fig.suptitle('Hybrid Hydrosphere Titan: Posterior (Andrade)', fontsize=14, y=1.02)
    _save(fig, 'hybrid_hydro_andrade_corner.png')

    # --- 2. k2 scatter coloured by D_hydro_km (1σ + 2σ ellipses) ---
    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(Re_arr, Im_arr, c=D_arr, cmap='plasma_r', s=8, alpha=0.6,
                    vmin=D_MIN, vmax=D_MAX)
    plt.colorbar(sc, ax=ax, label=r'$D_\mathrm{hydro}$ (km)')
    for nσ, ls, lw, lbl in [(1, '--', 2, r'obs 1$\sigma$'), (2, ':', 1, r'obs 2$\sigma$')]:
        ax.add_patch(Ellipse((OBS['Re_k2'], OBS['Im_k2']),
                             2*nσ*OBS['Re_k2_err'], 2*nσ*OBS['Im_k2_err'],
                             fill=False, color='red', linewidth=lw, linestyle=ls, label=lbl))
    ax.set_xlabel(r'$\mathrm{Re}(k_2)$');  ax.set_ylabel(r'$|\mathrm{Im}(k_2)|$')
    ax.set_title(r'Hybrid Hydrosphere: $k_2$ Posterior (by $D_\mathrm{hydro}$)')
    ax.legend()
    _save(fig, 'hybrid_hydro_andrade_k2_scatter.png')

    # --- 3. k2 scatter coloured by silicate heating fraction (Test42-style) ---
    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(Re_arr, Im_arr, c=f_sil, cmap='RdYlBu_r', s=8, alpha=0.6, vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax, label='Silicate heating fraction')
    for nσ, ls, lw, lbl in [(1, '--', 2, r'obs 1$\sigma$'), (2, ':', 1, r'obs 2$\sigma$')]:
        ax.add_patch(Ellipse((OBS['Re_k2'], OBS['Im_k2']),
                             2*nσ*OBS['Re_k2_err'], 2*nσ*OBS['Im_k2_err'],
                             fill=False, color='red', linewidth=lw, linestyle=ls, label=lbl))
    ax.set_xlabel(r'$\mathrm{Re}(k_2)$');  ax.set_ylabel(r'$|\mathrm{Im}(k_2)|$')
    ax.set_title(r'Hybrid Hydrosphere: $k_2$ Posterior (silicate heating)')
    ax.legend()
    _save(fig, 'hybrid_hydro_andrade_k2_scatter_heating.png')

    # --- 4. Heating distribution vs parameters (log10 power, 2×4 for 8 params) ---
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

    plot_xvals   = ([eval_samples[:, i] for i in range(5)]
                    + [d_iceIh_eval]
                    + [eval_samples[:, i] for i in range(6, 10)])
    plot_xlabels = (list(PARAM_LABELS[:5])
                    + [r'$D_\mathrm{IceIh}$ (km)']
                    + list(PARAM_LABELS[6:]))

    fig, axes = plt.subplots(2, 5, figsize=(25, 8))
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
    _save(fig, 'hybrid_hydro_andrade_heating.png')

    # --- 4b. Ice phase comparison: k2 by Ih dominance + Ih vs V scatter ---
    Ih_pw  = heating_power['Ih']
    V_pw   = heating_power['V']
    III_pw = heating_power['III']
    VI_pw  = heating_power['VI']
    ice_tot = Ih_pw + III_pw + V_pw + VI_pw
    f_ih = np.where(ice_tot > W_FLOOR, Ih_pw / (ice_tot + 1e-30), np.nan)
    valid_ice = np.isfinite(f_ih) & (ice_tot > W_FLOOR)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    ax = axes[0]
    if np.any(valid_ice):
        sc = ax.scatter(Re_arr[valid_ice], Im_arr[valid_ice], c=f_ih[valid_ice],
                        cmap='RdYlGn', s=8, alpha=0.7, vmin=0, vmax=1)
        plt.colorbar(sc, ax=ax, label='Ice Ih fraction of ice heating')
    else:
        ax.scatter(Re_arr, Im_arr, s=8, alpha=0.5, c='gray')
    for nσ, ls, lw, lbl in [(1, '--', 2, r'obs 1$\sigma$'), (2, ':', 1, r'obs 2$\sigma$')]:
        ax.add_patch(Ellipse((OBS['Re_k2'], OBS['Im_k2']),
                             2*nσ*OBS['Re_k2_err'], 2*nσ*OBS['Im_k2_err'],
                             fill=False, color='red', linewidth=lw, linestyle=ls, label=lbl))
    ax.set_xlabel(r'$\mathrm{Re}(k_2)$'); ax.set_ylabel(r'$|\mathrm{Im}(k_2)|$')
    ax.set_title(r'$k_2$ posterior: Ice Ih vs HP ice dominance (green=Ih, red=V/VI)')
    ax.legend(fontsize=8)

    ax = axes[1]
    ih_mask = valid_ice & (Ih_pw > W_FLOOR) & (V_pw > W_FLOOR)
    if np.any(ih_mask):
        sc = ax.scatter(np.log10(Ih_pw[ih_mask]), np.log10(V_pw[ih_mask]),
                        c=d_ocean_eval[ih_mask], cmap='plasma', s=8, alpha=0.7)
        plt.colorbar(sc, ax=ax, label=r'$D_\mathrm{ocean}$ (km)')
        mn = min(np.log10(Ih_pw[ih_mask]).min(), np.log10(V_pw[ih_mask]).min())
        mx = max(np.log10(Ih_pw[ih_mask]).max(), np.log10(V_pw[ih_mask]).max())
        ax.plot([mn, mx], [mn, mx], 'k--', lw=1, label='equal power')
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, 'No samples with both\nIce Ih and Ice V heating',
                ha='center', va='center', transform=ax.transAxes)
    ax.set_xlabel(r'$\log_{10}$ Ice Ih heating (W)')
    ax.set_ylabel(r'$\log_{10}$ Ice V heating (W)')
    ax.set_title('Ice Ih vs Ice V heating power')
    fig.tight_layout()
    _save(fig, 'hybrid_hydro_andrade_ice_comparison.png')

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
    ax.axhline(OBS['CMR2'], color='red', ls='--',
               label=f"obs={OBS['CMR2']}+/-{OBS['CMR2_err']}")
    ax.axhspan(OBS['CMR2'] - OBS['CMR2_err'], OBS['CMR2'] + OBS['CMR2_err'], alpha=0.15, color='red')
    ax.set_xlabel(r'$D_\mathrm{hydro}$ (km)');  ax.set_ylabel('C/MR²')
    ax.set_title('CMR2 vs Hydrosphere Thickness');  ax.legend()
    fig.tight_layout()
    _save(fig, 'hybrid_hydro_andrade_mass_cmr2.png')

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
    _save(fig, 'hybrid_hydro_andrade_cmr2_surface.png')

    # --- 7. Ice shell thickness vs Tb_K from grid (Test42-style), posterior overlay ---
    fig, ax = plt.subplots(figsize=(8, 5))
    tb_sorted   = sorted(tb_to_diceIh)
    diceIh_line = [tb_to_diceIh[tb] for tb in tb_sorted]
    if any(np.isfinite(v) for v in diceIh_line):
        ax.plot(tb_sorted, diceIh_line, 'k-o', markersize=4)
        ax2 = ax.twinx()
        ax2.hist(samples[:, 5], bins=30, alpha=0.3, color='blue',
                 density=True, label=r'Posterior $T_b$')
        ax2.set_ylabel('Posterior density', color='blue')
        ax.set_xlabel(r'$T_b$ (K)')
        ax.set_ylabel(r'Ice Ih shell thickness (km)')
        ax.set_title('Ice Shell Thickness vs Basal Temperature')
        fig.tight_layout()
    _save(fig, 'hybrid_hydro_andrade_Tb_structure.png')

    # --- 8. Interior structure profile (posterior layer boundaries) ---
    _plot_structure_profile(
        samples, grid_cache, tb_vals, d_vals, eval_idx,
        os.path.join(OUTPUT_DIR, 'hybrid_hydro_andrade_structure_profile.png'),
    )

    # --- 9. Layer structure vs D_ocean (3-panel: histogram / thicknesses / heating) ---
    _plot_layer_structure_vs_docean(
        samples, grid_cache, tb_vals, d_vals, eval_idx, heating_results,
        os.path.join(OUTPUT_DIR, 'hybrid_hydro_andrade_layers_vs_docean.png'),
    )

    plt.close('all')


def _save(fig, fname):
    import matplotlib.pyplot as plt
    path = os.path.join(OUTPUT_DIR, fname)
    fig.savefig(path, dpi=150, bbox_inches='tight')
    log.info(f'Saved {path}')
    plt.close(fig)


def _plot_layer_structure_vs_docean(samples, grid_cache, tb_vals, d_vals,
                                     eval_idx, heating_results, output_path):
    """3-panel: all posterior models plotted individually.

    Top:    D_ocean posterior histogram
    Middle: cumulative layer thickness scatter (each eval sample = one column)
            sorted by D_ocean — shows full distribution of models
    Bottom: per-phase heating (GW) vs silicate heating fraction
    """
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec

    PHASE_COLORS = {
        'Clath': '#D4F1F9', 'Ih': '#AEE1F8', '0': '#1E90FF',
        'III': '#C97BAE', 'V': '#9B59B6', 'VI': '#6C3483',
        'Sil': '#C8A96E', 'Core': '#8B5A2B',
    }
    PHASE_LABELS = {
        'Clath': 'Clathrate', 'Ih': 'Ice Ih', '0': 'Ocean',
        'III': 'Ice III', 'V': 'Ice V', 'VI': 'Ice VI',
        'Sil': 'Silicate', 'Core': 'Core',
    }
    R_T_km = TITAN_R_M / 1e3

    # Extract layer thicknesses for EACH evaluated posterior sample
    model_data = []
    for ii, si in enumerate(eval_idx):
        tb, d_hydro, f_core_i = samples[si, 5], samples[si, 6], samples[si, 9]
        idx_tb = int(np.argmin(np.abs(tb_vals - tb)))
        idx_d  = int(np.argmin(np.abs(d_vals  - d_hydro)))
        pt = grid_cache.get((float(tb_vals[idx_tb]), float(d_vals[idx_d])))
        if pt is None:
            continue
        d_ih  = pt.get('D_iceIh_km', 0.0)
        d_oc  = pt.get('D_ocean_km', 0.0)
        d_iii = pt.get('D_iceIII_km', 0.0)
        d_v   = pt.get('D_iceV_km', 0.0)
        d_vi  = pt.get('D_iceVI_km', 0.0)
        d_hp  = pt.get('D_hp_ice_km', 0.0)
        if d_iii == 0 and d_v == 0 and d_vi == 0 and d_hp > 0:
            d_v = d_hp
        d_hydro_actual = d_ih + d_oc + d_iii + d_v + d_vi
        r_sil_km  = R_T_km - d_hydro_actual
        r_core_km = f_core_i * r_sil_km
        d_sil = r_sil_km - r_core_km
        d_core = r_core_km
        model_data.append({
            'D_ocean': d_oc, 'Ih': d_ih, '0': d_oc,
            'III': d_iii, 'V': d_v, 'VI': d_vi,
            'Sil': d_sil, 'Core': d_core, 'idx': ii,
        })

    if not model_data:
        log.warning('_plot_layer_structure_vs_docean: no valid samples, skipping.')
        return

    model_data.sort(key=lambda r: r['D_ocean'])
    d_ocean_sorted = np.array([r['D_ocean'] for r in model_data])
    no_ocean_frac = np.mean(d_ocean_sorted < 0.5)

    fig = plt.figure(figsize=(10, 12))
    gs  = gridspec.GridSpec(3, 1, height_ratios=[1, 4, 3], hspace=0.12)
    ax_dens   = fig.add_subplot(gs[0])
    ax_struct = fig.add_subplot(gs[1], sharex=ax_dens)
    ax_heat   = fig.add_subplot(gs[2])

    # --- Top: posterior D_ocean density ---
    ax_dens.hist(d_ocean_sorted, bins=40, alpha=0.5, color='steelblue',
                 density=True, edgecolor='none')
    ax_dens.axvline(0.5, color='red', ls=':', lw=1.5, label='no-ocean boundary')
    ax_dens.set_ylabel('Density', fontsize=9)
    ax_dens.set_title(f'Titan Interior Structure — All Posterior Models'
                      f'  (no-ocean: {no_ocean_frac:.0%})', fontsize=13)
    ax_dens.tick_params(labelsize=8)
    ax_dens.legend(fontsize=8)
    plt.setp(ax_dens.get_xticklabels(), visible=False)

    # --- Middle: stackplot of ALL models sorted by D_ocean ---
    # Stack order (bottom to top): Core, Sil, VI, V, III, Ocean, Ih
    stack_phases = ['Core', 'Sil', 'VI', 'V', 'III', '0', 'Ih']
    x = np.arange(len(model_data))
    stacks = [np.array([r.get(p, 0.0) for r in model_data]) for p in stack_phases]
    colors = [PHASE_COLORS.get(p, '#cccccc') for p in stack_phases]
    labels = [PHASE_LABELS.get(p, p) for p in stack_phases]

    polys = ax_struct.stackplot(x, *stacks, colors=colors, labels=labels)
    for poly in polys:
        poly.set_linewidth(0)
        poly.set_edgecolor('none')

    ax_struct.axhline(R_T_km, color='k', ls='-', lw=0.5, alpha=0.3)
    ax_struct.set_ylabel('Cumulative thickness (km)', fontsize=12)
    ax_struct.set_ylim(0, R_T_km * 1.02)

    # Add a secondary x-axis showing D_ocean values
    n_ticks = 8
    tick_positions = np.linspace(0, len(model_data) - 1, n_ticks, dtype=int)
    ax_struct.set_xticks(tick_positions)
    ax_struct.set_xticklabels([f'{d_ocean_sorted[i]:.0f}' for i in tick_positions])
    ax_struct.set_xlabel(r'$D_\mathrm{ocean}$ (km)', fontsize=11)

    handles, lbls = ax_struct.get_legend_handles_labels()
    ax_struct.legend(handles[::-1], lbls[::-1],
                     loc='center left', bbox_to_anchor=(1.02, 0.5),
                     fontsize=9, frameon=True)
    ax_struct.tick_params(labelsize=9)

    # --- Bottom: per-phase heating vs silicate heating fraction ---
    f_sil = np.array([
        h.get('Sil', 0) / max(sum(h.values()), 1e-30)
        for h in heating_results
    ])
    total_heat = np.array([sum(h.values()) for h in heating_results])
    valid_heat = total_heat > 1e6

    for phase in ['Ih', 'III', 'V', 'VI', 'Sil']:
        vals = np.array([h.get(phase, 0) / 1e9 for h in heating_results])
        mask = valid_heat & (vals > 1e-3)
        if np.any(mask):
            ax_heat.scatter(f_sil[mask], vals[mask], s=8, alpha=0.4,
                            color=PHASE_COLORS.get(phase, '#C8A96E'),
                            label=PHASE_LABELS.get(phase, phase))
    ax_heat.set_yscale('log')
    ax_heat.set_xlabel('Silicate heating fraction', fontsize=12)
    ax_heat.set_ylabel('Phase heating power (GW)', fontsize=12)
    ax_heat.set_xlim(-0.02, 1.02)
    ax_heat.axhline(400, color='gray', ls='--', lw=1, alpha=0.6, label='Titan equil. (~400 GW)')
    ax_heat.legend(fontsize=9, loc='upper right')
    ax_heat.tick_params(labelsize=9)
    ax_heat.set_title('Per-Phase Tidal Heating vs Silicate Fraction', fontsize=12)

    fig.savefig(output_path, dpi=200, bbox_inches='tight')
    log.info(f'Saved {output_path}')
    plt.close(fig)


def _plot_structure_profile(samples, grid_cache, tb_vals, d_vals, eval_idx, output_path):
    """Posterior interior structure as a PP-style wedge diagram.

    Shows median layer radii with 5/95 percentile uncertainty arcs.
    Matches PlanetProfile's PlotWedge() visual style.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Wedge as MplWedge, FancyArrowPatch
    import matplotlib.patches as mpatches
    from PlanetProfile.Inference.mcmc_plots import _wedge_color_map

    COLORS = _wedge_color_map()
    R_T_km = TITAN_R_M / 1e3
    ANG1, ANG2 = 55, 125  # wedge angular extent

    # Collect radii for each layer boundary across posterior
    r_iceIh_bot, r_ocean_bot = [], []
    r_iceIII_bot, r_iceV_bot, r_iceVI_bot = [], [], []
    r_sil_bot, r_core_bot = [], []

    for i in eval_idx:
        tb, d_hydro, f_core_i = samples[i, 5], samples[i, 6], samples[i, 9]
        idx_tb = int(np.argmin(np.abs(tb_vals - tb)))
        idx_d  = int(np.argmin(np.abs(d_vals  - d_hydro)))
        pt = grid_cache.get((float(tb_vals[idx_tb]), float(d_vals[idx_d])))
        if pt is None:
            continue
        d_ih  = pt.get('D_iceIh_km', 0.0)
        d_oc  = pt.get('D_ocean_km', 0.0)
        d_iii = pt.get('D_iceIII_km', 0.0)
        d_v   = pt.get('D_iceV_km', 0.0)
        d_vi  = pt.get('D_iceVI_km', 0.0)
        d_hp  = pt.get('D_hp_ice_km', 0.0)
        if d_iii == 0 and d_v == 0 and d_vi == 0 and d_hp > 0:
            d_v = d_hp

        # Radii (from center): surface = R_T_km
        r_ih_b  = R_T_km - d_ih
        r_oc_b  = r_ih_b - d_oc
        r_iii_b = r_oc_b - d_iii
        r_v_b   = r_iii_b - d_v
        r_vi_b  = r_v_b - d_vi
        r_sil_b = r_vi_b  # top of silicate = bottom of last ice
        r_core_top = f_core_i * r_sil_b  # core radius

        r_iceIh_bot.append(r_ih_b)
        r_ocean_bot.append(r_oc_b)
        r_iceIII_bot.append(r_iii_b)
        r_iceV_bot.append(r_v_b)
        r_iceVI_bot.append(r_vi_b)
        r_sil_bot.append(r_core_top)
        r_core_bot.append(0.0)

    if not r_iceIh_bot:
        log.warning('_plot_structure_profile: no valid samples, skipping.')
        return

    # Compute percentiles for each boundary
    def pct(arr):
        return np.percentile(arr, [5, 50, 95])

    p_ih  = pct(r_iceIh_bot)
    p_oc  = pct(r_ocean_bot)
    p_iii = pct(r_iceIII_bot)
    p_v   = pct(r_iceV_bot)
    p_vi  = pct(r_iceVI_bot)
    p_sil = pct(r_sil_bot)

    # Use median radii for the wedge
    layers = [
        ('Ice Ih',     R_T_km,   p_ih[1]),
        ('Ocean',      p_ih[1],  p_oc[1]),
        ('Ice III',    p_oc[1],  p_iii[1]),
        ('Ice V',      p_iii[1], p_v[1]),
        ('Ice VI',     p_v[1],   p_vi[1]),
        ('Silicate',   p_vi[1],  p_sil[1]),
        ('Dense core', p_sil[1], 0.0),
    ]

    fig, ax = plt.subplots(figsize=(6, 8))
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.05, 1.05)
    ax.set_aspect('equal')
    ax.axis('off')

    cx, cy = 0.5, 0.0  # wedge center at bottom

    # Draw layers as wedge patches (outer to inner)
    for name, r_top, r_bot in layers:
        thickness_km = r_top - r_bot
        if thickness_km < 0.5:
            continue
        r_outer_norm = r_top / R_T_km
        width_norm = thickness_km / R_T_km
        wedge = MplWedge((cx, cy), r_outer_norm, ANG1, ANG2,
                         width=width_norm,
                         fc=COLORS[name], ec='#333333', lw=0.8)
        ax.add_patch(wedge)

    # Add labels on the left side with arrows to avoid overlap
    tick_angle_rad = np.radians(ANG1 - 3)
    label_entries = []
    for name, r_top, r_bot in layers:
        if r_top - r_bot < 5:
            continue
        r_mid = (r_top + r_bot) / 2
        r_norm = r_mid / R_T_km
        x = cx + r_norm * np.cos(tick_angle_rad)
        y = cy + r_norm * np.sin(tick_angle_rad)
        label_entries.append((x, y, name, r_top - r_bot))

    # Space labels vertically to avoid overlap (min 0.06 apart in norm coords)
    label_y_positions = [e[1] for e in label_entries]
    spaced_y = []
    for i, y in enumerate(sorted(label_y_positions)):
        if i > 0 and y - spaced_y[-1] < 0.06:
            y = spaced_y[-1] + 0.06
        spaced_y.append(y)
    # Map back in original order
    sorted_idx = sorted(range(len(label_y_positions)), key=lambda i: label_y_positions[i])
    final_y = [0.0] * len(label_entries)
    for rank, orig_i in enumerate(sorted_idx):
        final_y[orig_i] = spaced_y[rank]

    for i, (x, y, name, thick) in enumerate(label_entries):
        label_x = -0.08
        label_y = final_y[i]
        ax.annotate(f'{name} ({thick:.0f} km)',
                    xy=(x, y), xytext=(label_x, label_y),
                    fontsize=8, ha='right', va='center', color='#333333',
                    arrowprops=dict(arrowstyle='-', color='#666666', lw=0.6))

    # Surface radius label
    ax.text(cx, R_T_km / R_T_km + 0.02 + cy,
            f'R = {R_T_km:.0f} km', ha='center', va='bottom', fontsize=9)

    # Title with posterior info
    d_ocean_med = R_T_km - p_ih[1] - (p_ih[1] - p_oc[1])
    # Actually: ocean thickness = p_ih[1] - p_oc[1]
    ocean_thick = p_ih[1] - p_oc[1]
    ice_ih_thick = R_T_km - p_ih[1]
    hp_thick = p_oc[1] - p_vi[1]
    ax.set_title(f'Titan Interior Structure (Posterior Median)\n'
                 f'Ice Ih: {ice_ih_thick:.0f} km | Ocean: {ocean_thick:.0f} km | '
                 f'HP ice: {hp_thick:.0f} km',
                 fontsize=11, pad=10)

    # Legend
    handles = [mpatches.Patch(color=c, label=l) for l, c in COLORS.items()
               if any(l == name and r_top - r_bot > 0.5 for name, r_top, r_bot in layers)]
    ax.legend(handles=handles, loc='lower left', fontsize=8, framealpha=0.9)

    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    log.info(f'Saved {output_path}')
    plt.close(fig)


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PPTest46 Andrade hybrid MCMC')
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

    pkl_path = os.path.join(OUTPUT_DIR, 'hybrid_hydro_andrade_mcmc.pkl')

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

    # Sanity check: probe a spread of grid points with mass-consistent densities.
    rng_sc = np.random.default_rng(0)
    _sc_tb = rng_sc.choice(tb_vals, 20)
    _sc_d  = rng_sc.choice(d_vals,  20)
    _n_ok = 0
    for _tb, _d in zip(_sc_tb, _sc_d):
        # Derive mass-consistent rho_sil (uniform, no core) for this grid point
        _key = (float(tb_vals[int(np.argmin(np.abs(tb_vals - _tb)))]),
                float(d_vals[int(np.argmin(np.abs(d_vals - _d)))]))
        _data = grid_cache.get(_key)
        if _data is None:
            log.error(f'Grid lookup failed for Tb={_tb:.1f} K, D={_d:.0f} km — cache may be corrupt.')
            sys.exit(1)
        _r_sil = _data['R_body_m'] - _data['D_hydro_km'] * 1e3
        _V_sil = (4.0/3.0) * np.pi * _r_sil**3
        _M_hydro = _data['Mtot_kg'] - _data['rhoSil_kgm3'] * _V_sil
        _rho_sc = (MTOT_OBS - _M_hydro) / _V_sil
        if _rho_sc < 2000 or _rho_sc > 5000:
            continue  # Skip grid points where no physical density matches Titan mass
        _theta = [0.3, -1.0, 14.0, 13.0, 20.0, float(_tb), float(_d), _rho_sc, _rho_sc, 0.0]
        _Re, _Im, _M, _C, _ = forward_model(_theta, grid_cache, tb_vals, d_vals)
        if not np.isnan(_Re):
            _n_ok += 1
    log.info(f'Grid health check: {_n_ok}/20 sampled points returned valid k2 '
             f'(remainder are mass-inconsistent or TidalPy solver failures)')
    if _n_ok == 0:
        log.error('All 20 sampled grid points returned NaN k2 — TidalPy may be broken.')
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
            'Re_k2':     (OBS['Re_k2'],  OBS['Re_k2_err']),
            'abs_Im_k2': (OBS['Im_k2'],  OBS['Im_k2_err']),
            'CMR2':      (OBS['CMR2'],   OBS['CMR2_err']),
            'Mtot_kg':   (MTOT_OBS,   MTOT_ERR),
        },
    }
    pkl_path = os.path.join(OUTPUT_DIR, 'hybrid_hydro_andrade_mcmc.pkl')
    with open(pkl_path, 'wb') as f:
        pickle.dump(results, f)
    log.info(f'Results saved to {pkl_path}')

    if not args.no_plots:
        make_plots(samples, log_likes, k2_results, mtot_results, cmr2_results,
                   heating_results, eval_idx, grid_cache, tb_vals, d_vals)

    # Summary
    print('\n' + '=' * 65)
    print('HYBRID HYDROSPHERE TITAN MCMC SUMMARY (ANDRADE)')
    print('=' * 65)
    for ip, name in enumerate(PARAM_NAMES):
        med = np.median(samples[:, ip])
        lo  = np.percentile(samples[:, ip], 16)
        hi  = np.percentile(samples[:, ip], 84)
        print(f'  {name:20s}: {med:.3f}  [{lo:.3f}, {hi:.3f}]')

    eval_s = samples[eval_idx]

    # No-ocean fraction from posterior
    def _docean_for_sample(s):
        _i_tb = int(np.argmin(np.abs(tb_vals - s[5])))
        _i_d  = int(np.argmin(np.abs(d_vals  - s[6])))
        _pt   = grid_cache.get((float(tb_vals[_i_tb]), float(d_vals[_i_d])))
        return _pt.get('D_ocean_km', 0.0) if _pt is not None else np.nan
    d_ocean_all = np.array([_docean_for_sample(samples[i]) for i in eval_idx])
    no_ocean_frac = np.mean(d_ocean_all < 0.5) if len(d_ocean_all) > 0 else 0

    _re = OBS['Re_k2'];  _re_e = OBS['Re_k2_err']
    _im = OBS['Im_k2'];  _im_e = OBS['Im_k2_err']
    _cm = OBS['CMR2'];   _cm_e = OBS['CMR2_err']
    print(f"\n  Re(k2) median:    {np.nanmedian(k2_results[:,0]):.4f}  "
          f"(obs {_re:.3f} +/- {_re_e:.3f})")
    print(f"  |Im(k2)| median:  {np.nanmedian(np.abs(k2_results[:,1])):.4f}  "
          f"(obs {_im:.3f} +/- {_im_e:.3f})")
    print(f"  CMR2 median:      {np.nanmedian(cmr2_results):.5f}  "
          f"(obs {_cm:.3f} +/- {_cm_e:.3f})")
    print(f'  Mtot median:      {np.nanmedian(mtot_results):.4e} kg  '
          f'(obs {MTOT_OBS:.4e} ± {MTOT_ERR:.1e})')
    print(f'  No-ocean fraction: {no_ocean_frac:.1%}  (D_ocean < 0.5 km)')
    print(f'\n  Output: {OUTPUT_DIR}')
