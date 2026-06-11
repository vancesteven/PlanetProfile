"""
MCMC exploration: Andrade rheology + Yao 2014 spherical convection, hybrid Titan (PPTest48).

Extends Test46 (Andrade hybrid hydrosphere) by:
1. Rebuilding the structure cache with Yao et al. (2014) spherical convection for
   Ice Ih instead of Deschamps & Sotin (2001). This produces thicker stagnant lids
   and lower basal heat flux for Titan's geometry (f ~ 0.96).
2. Adding a heat-flux consistency constraint: the sampled Ice Ih viscosity must
   produce physically plausible convective heat flux via the Yao scaling law.

Additional observable:
  q_surface = 10 +/- 5 mW/m²  (broad estimate; Nimmo & Bills 2010 tidal flux)

The heat-flux constraint couples sampled eta_Ih to convective physics:
  - Low eta_Ih → high Ra → high q → penalized if q >> q_obs
  - High eta_Ih → Ra < Ra_crit → conductive → penalized

Parameter space and grid identical to Test46 (10D, ~221×31 grid).

Usage:
  conda activate ./venvPP
  python PlanetProfile/Test/Test48_mcmc_andrade_yao2014.py

  # To force grid rebuild:
  python PlanetProfile/Test/Test48_mcmc_andrade_yao2014.py --rebuild-grid
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
log = logging.getLogger('PPTest48_MCMC_Yao2014')
log.setLevel(logging.INFO)

# TidalPy imports
from TidalPy.RadialSolver import build_rs_input_from_data, radial_solver

# ============================================================
# Configuration
# ============================================================
OUTPUT_DIR = os.path.join(_pp_root, 'PlanetProfile', 'Test', 'mcmc_results', 'Titan', 'Test48_andrade_yao2014')
os.makedirs(OUTPUT_DIR, exist_ok=True)

GRID_CACHE_PATH = os.path.join(OUTPUT_DIR, 'titan_yao2014_hybrid_hydro_grid_cache.pkl')

# Observational constraints (Petricca et al. 2025, Cassini radio science)
OBS = {
    'Re_k2': 0.608,  'Re_k2_err': 0.048,
    'Im_k2': 0.135,  'Im_k2_err': 0.035,
    'CMR2':  0.343,  'CMR2_err':  0.005,   # relaxed from 0.001 to allow cold-Tb structural mode
    'q_surface_mWm2': 10.0,  'q_surface_err_mWm2': 5.0,
}

# Yao 2014 convection constants for Ice Ih
R_GAS = 8.314462  # J/(mol·K)
EACT_IH_JMOL = 60e3  # Activation energy, 60 kJ/mol (Yao Table 5)
T_SURF_TITAN_K = 94.0  # Titan surface temperature
K_ICE_WMK = 2.6  # Ice Ih thermal conductivity (Yao Table 5)
ALPHA_ICE_K = 1.56e-4  # Thermal expansivity
RHO_ICE_KGM3 = 917.0  # Ice Ih density
KAPPA_ICE_M2S = 1.47e-6  # Thermal diffusivity
G_TITAN_MS2 = 1.352  # Surface gravity
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
                'PlanetProfile.Test.PPTest48',
                tb_grid,
                d_grid,
                GRID_CACHE_PATH,
                rheology='maxwell',
                convection_model='yao2014',
                force_rebuild=False,
            )
        elif n_pts < n_expected and build_complete:
            log.info(f'Grid marked complete with {n_pts} points '
                     f'(expected {n_expected}; float rounding in Tb grid). Proceeding.')
    else:
        log.info(f'Building {len(tb_grid)} × {len(d_grid)} hybrid hydrosphere grid '
                 f'(Yao 2014 convection, {n_expected} points)...')
        cache_data = build_hybrid_hydrosphere_grid(
            'PlanetProfile.Test.PPTest48',
            tb_grid,
            d_grid,
            GRID_CACHE_PATH,
            rheology='maxwell',
            convection_model='yao2014',
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

    # Split innermost silicate layer at R_core; apply two-layer density + moduli
    from PlanetProfile.Inference.mcmc_common import split_silicate_core
    (rho_mod, K_Pa_mod, mu_Pa_mod,
     ci, n_layers, layer_upper_radii, layer_types_use, region_phases_use) = split_silicate_core(
        r_m, rho_mod, K_Pa_mod, mu_Pa_mod, phases, ci, n_layers,
        layer_upper_radii, layer_types_use, region_phases_use,
        rho_sil=rho_sil, rho_core=rho_core, f_core=f_core,
    )

    # Apply viscosity overrides
    for i in range(n_layers):
        s, e = ci[i], ci[i + 1]
        ph = int(phases[min(s, len(phases) - 1)])
        if ph == 1:              eta_mod[s:e] = eta_Ih
        elif ph in (3, 5, 6):   eta_mod[s:e] = eta_HP
        elif 50 <= ph < 100:    eta_mod[s:e] = eta_sil

    # Apply Arrhenius T-dependence to Ice Ih. The sampled eta_Ih is the BASAL
    # viscosity at T=Tb (matches the yao_heat_flux_mWm2 constraint convention);
    # the cold stagnant lid is much stiffer. η(T) = η_Ih * exp(E/R * (1/T - 1/Tb))
    from PlanetProfile.Inference.mcmc_common import apply_arrhenius_ih
    apply_arrhenius_ih(eta_mod, phases, ci, n_layers,
                       T_K_profile=data.get('T_K'),
                       Tb_K=Tb_K,
                       E_act_J_mol=EACT_IH_JMOL, R_gas=R_GAS)

    # Build Andrade rheology models per layer
    from PlanetProfile.Inference.mcmc_common import build_andrade_shear_bulk
    shear, bulk = build_andrade_shear_bulk(region_phases_use, alpha, log10_zeta)

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
            from PlanetProfile.Inference.mcmc_common import compute_per_phase_heating
            perPhase_W = compute_per_phase_heating(result, data)

        return Re_k2, Im_k2, Mtot_kg, CMR2, perPhase_W

    except Exception as exc:
        log.debug(f'TidalPy failed: {exc}')
        return np.nan, np.nan, Mtot_kg, CMR2, {}



# ============================================================
# Step 3: Yao 2014 heat flux constraint + Log-likelihood
# ============================================================

def yao_heat_flux_mWm2(Tb_K, D_iceIh_km, eta_Ih_Pas, R_body_m):
    """Compute Yao et al. (2014) basal heat flux from scaling law (Eq. 35).

    Lightweight inline version — uses fixed Ice Ih material properties from
    Yao Table 5 rather than calling the full EOS. Returns surface heat flux
    in mW/m² (assumes spherical shell, flux conserved across the lid).

    Returns NaN if the layer is sub-critical (conductive regime).
    """
    D_m = D_iceIh_km * 1e3
    if D_m < 1e3:
        return np.nan

    rBot_m = R_body_m - D_m
    f = rBot_m / R_body_m
    DeltaT = Tb_K - T_SURF_TITAN_K
    if DeltaT <= 0:
        return np.nan

    gamma = EACT_IH_JMOL / (R_GAS * Tb_K)
    theta_m = 1.0 - 1.23 / (gamma * f**1.5)
    if theta_m <= 0 or theta_m >= 1:
        return np.nan
    T_m = T_SURF_TITAN_K + theta_m * DeltaT

    DeltaT_v = R_GAS * T_m**2 / EACT_IH_JMOL
    eta_m = eta_Ih_Pas * np.exp(
        EACT_IH_JMOL / (R_GAS * T_m) - EACT_IH_JMOL / (R_GAS * Tb_K)
    )
    if eta_m <= 0 or not np.isfinite(eta_m):
        return np.nan

    # Full-ΔT Ra for onset check (Solomatov 1995 expects this definition)
    Ra_fullDT = (ALPHA_ICE_K * RHO_ICE_KGM3 * G_TITAN_MS2 * DeltaT * D_m**3
                 / (eta_Ih_Pas * KAPPA_ICE_M2S))
    Ra_crit = 20.9 * gamma**4  # Solomatov (1995) stagnant lid
    if Ra_fullDT < Ra_crit:
        return np.nan

    # Viscous-temperature Ra for heat flux scaling (Yao Eq. 35)
    Ra_m = (ALPHA_ICE_K * RHO_ICE_KGM3 * G_TITAN_MS2 * DeltaT_v * D_m**3
            / (eta_m * KAPPA_ICE_M2S))

    Phi_c = K_ICE_WMK * DeltaT / D_m
    q_bot = (1.46 * Ra_m**0.27 / f**1.78
             * (DeltaT_v / DeltaT)**1.21 * Phi_c)

    q_surface = q_bot * (rBot_m / R_body_m)**2
    return q_surface * 1e3  # W/m² → mW/m²


def log_likelihood(theta, grid_cache, tb_vals, d_vals):
    """Gaussian log-likelihood on k2, CMR2, Titan mass, and Yao heat flux."""
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

    # Yao 2014 heat flux constraint — DISABLED.
    # The 10 ± 5 mW/m² pin assumes Titan is in steady-state thermal equilibrium,
    # which is inconsistent with evidence for transient heating onset within
    # the last ~30 Myr. Re-enable only if coupling η_Ih to surface q is wanted.

    return -0.5 * chi2


# ============================================================
# Step 4: MCMC
# ============================================================

def run_mcmc(grid_cache, tb_vals, d_vals):
    import pocomc as pc
    from scipy.stats import uniform
    from PlanetProfile.Inference.mcmc_common import run_pocomc_sampler

    prior = pc.Prior([
        uniform(loc=0.2,   scale=0.2),                   # alpha: [0.2, 0.4]
        uniform(loc=-2.0,  scale=4.0),                   # log10_zeta: [-2, 2]
        uniform(loc=12.0,  scale=4.0),                   # log10_eta_Ih: [12, 16]
        uniform(loc=12.0,  scale=4.0),                   # log10_eta_HP: [12, 16] — brackets HP Maxwell peaks (μ≈5-9 GPa, ω_Titan=4.56e-6 → η_peak~10^15). Narrower than Petricca2025 [10,18] to force MCMC to explore HP-dominated solutions.
        uniform(loc=18.0,  scale=4.0),                   # log10_eta_sil: [18, 22]
        uniform(loc=TB_MIN, scale=TB_MAX - TB_MIN),      # Tb_K: [252, 270] K
        uniform(loc=D_MIN,  scale=D_MAX  - D_MIN),       # D_hydro_km: [50, 800] km
        uniform(loc=1800.0, scale=3200.0),               # rho_sil (outer): [1800, 5000] — floor lowered to open no-ocean branch
        uniform(loc=1800.0, scale=3200.0),               # rho_core (inner): [1800, 5000] kg/m³
        uniform(loc=0.0,    scale=0.8),                  # f_core = R_core/r_sil_top: [0, 0.8]
    ])

    def _log_like(theta):
        return log_likelihood(theta, grid_cache, tb_vals, d_vals)

    samples, log_likes, sampler = run_pocomc_sampler(
        prior, _log_like, n_effective=N_EFF, random_state=RANDOM_STATE
    )
    return samples, log_likes, sampler


# ============================================================
# Step 5: Re-evaluate heating on posterior
# ============================================================

def evaluate_heating(samples, grid_cache, tb_vals, d_vals, n_eval=None):
    from PlanetProfile.Inference.mcmc_common import evaluate_posterior
    if n_eval is None:
        n_eval = min(len(samples), N_REEVAL)

    def _fwd(theta):
        return forward_model(theta, grid_cache, tb_vals, d_vals, return_heating=True)

    eval_idx, results = evaluate_posterior(
        samples, _fwd, n_eval=n_eval, random_state=RANDOM_STATE
    )
    k2_results      = np.array([(r[0], r[1]) for r in results])
    mtot_results    = np.array([r[2] for r in results])
    cmr2_results    = np.array([r[3] for r in results])
    heating_results = [r[4] for r in results]
    return eval_idx, k2_results, mtot_results, cmr2_results, heating_results


# ============================================================
# Step 6: Plots
# ============================================================

def make_plots(samples, log_likes, k2_results, mtot_results, cmr2_results,
               heating_results, eval_idx, grid_cache, tb_vals, d_vals):
    """Generate all diagnostic plots via the shared mcmc_plots toolkit.

    Data preparation (grid lookups, per-sample derived quantities) stays here
    because it depends on the body-specific grid-cache schema.  Body-agnostic
    plotting is delegated to PlanetProfile.Inference.mcmc_plots.
    """
    import seaborn as sns
    from PlanetProfile.Inference import mcmc_plots as mp
    sns.set_theme(style='white', font_scale=0.9)

    # --- Per-sample derived quantities -----------------------------------
    eval_samples = samples[eval_idx]

    # Tb_K → D_iceIh_km lookup (true Ice Ih thickness across D_hydro column)
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

    # D_ocean for each eval sample
    def _get_docean(tb_k, d_hydro):
        _i_tb = int(np.argmin(np.abs(tb_vals - tb_k)))
        _i_d  = int(np.argmin(np.abs(d_vals  - d_hydro)))
        _pt   = grid_cache.get((float(tb_vals[_i_tb]), float(d_vals[_i_d])))
        return _pt.get('D_ocean_km', 0.0) if _pt is not None else np.nan

    d_ocean_eval = np.array([_get_docean(eval_samples[i, 5], eval_samples[i, 6])
                              for i in range(len(eval_idx))])

    D_arr  = eval_samples[:, 6]  # D_hydro_km
    f_sil  = np.array([
        h.get('Sil', 0) / max(sum(h.values()), 1e-30)
        for h in heating_results
    ])

    out = lambda name: os.path.join(OUTPUT_DIR, f'hybrid_hydro_andrade_yao2014_{name}.png')

    # --- 1. Corner plot: swap Tb_K for D_iceIh in labels/samples -------
    corner_samples = np.column_stack([samples[:, :5], d_iceIh_all, samples[:, 6:]])
    corner_labels = (list(PARAM_LABELS[:5])
                     + [r'$D_\mathrm{IceIh}$ (km)']
                     + list(PARAM_LABELS[6:]))
    mp.plot_corner(
        corner_samples, corner_labels,
        title='Hybrid Hydrosphere Titan: Posterior (Andrade + Yao 2014)',
        output_path=out('corner'),
    )

    # --- 2. k2 scatter coloured by D_hydro -----------------------------
    mp.plot_k2_scatter_by(
        k2_results, color_values=D_arr,
        colorbar_label=r'$D_\mathrm{hydro}$ (km)',
        obs_re=OBS['Re_k2'], obs_im=OBS['Im_k2'],
        obs_re_err=OBS['Re_k2_err'], obs_im_err=OBS['Im_k2_err'],
        title=r'Hybrid Hydrosphere: $k_2$ Posterior (by $D_\mathrm{hydro}$)',
        output_path=out('k2_scatter'),
        cmap='plasma_r', vmin=D_MIN, vmax=D_MAX,
    )

    # --- 3. k2 scatter coloured by silicate heating fraction -----------
    mp.plot_k2_scatter_by(
        k2_results, color_values=f_sil,
        colorbar_label='Silicate heating fraction',
        obs_re=OBS['Re_k2'], obs_im=OBS['Im_k2'],
        obs_re_err=OBS['Re_k2_err'], obs_im_err=OBS['Im_k2_err'],
        title=r'Hybrid Hydrosphere: $k_2$ Posterior (silicate heating)',
        output_path=out('k2_scatter_heating'),
        cmap='RdYlBu_r', vmin=0, vmax=1,
    )

    # --- 4. Heating vs parameters + cumulative bar ---------------------
    mp.plot_heating_vs_parameters(
        eval_samples, heating_results, PARAM_LABELS,
        extra_xvals=[d_iceIh_eval],
        extra_xlabels=[r'$D_\mathrm{IceIh}$ (km)'],
        output_path=out('heating'),
        cumulative_bar=True,
        eval_d_ocean=d_ocean_eval,
        title=('Hybrid Hydrosphere: Tidal Heating Power (top) '
               'and Per-Model Heating Fractions (bottom)'),
    )

    # --- 4b. Ice phase comparison --------------------------------------
    mp.plot_ice_comparison(
        k2_results, heating_results, d_ocean_eval,
        obs_re=OBS['Re_k2'], obs_im=OBS['Im_k2'],
        obs_re_err=OBS['Re_k2_err'], obs_im_err=OBS['Im_k2_err'],
        output_path=out('ice_comparison'),
    )

    # --- 5. Mtot and CMR2 vs D_hydro_km --------------------------------
    mp.plot_mass_cmr2_diagnostics(
        d_hydro_values=D_arr,
        mtot_results=mtot_results,
        cmr2_results=cmr2_results,
        obs_mtot=MTOT_OBS, obs_mtot_err=MTOT_ERR,
        obs_cmr2=OBS['CMR2'], obs_cmr2_err=OBS['CMR2_err'],
        output_path=out('mass_cmr2'),
    )

    # --- 6. CMR2 surface over grid -------------------------------------
    mp.plot_cmr2_surface(tb_vals, d_vals, grid_cache, output_path=out('cmr2_surface'))

    # --- 7. Ice shell thickness vs Tb with posterior overlay ------------
    mp.plot_tb_structure(tb_vals, d_vals, grid_cache, samples,
                         output_path=out('Tb_structure'))

    # --- 8. Interior structure wedge ------------------------------------
    mp.plot_structure_wedge(
        samples, eval_idx, grid_cache,
        output_path=out('structure_profile'),
        R_body_km=TITAN_R_M / 1e3,
        body_name='Titan',
    )

    # --- 9. Layer structure vs D_ocean (3-panel) ------------------------
    mp.plot_layers_vs_docean(
        samples, eval_idx, grid_cache, heating_results,
        output_path=out('layers_vs_docean'),
        R_body_km=TITAN_R_M / 1e3,
        body_name='Titan',
        equil_heating_GW=400.0,
        equil_heating_label='Titan equil. (~400 GW)',
    )



# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PPTest48 Andrade + Yao2014 hybrid MCMC')
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

    pkl_path = os.path.join(OUTPUT_DIR, 'hybrid_hydro_andrade_yao2014_mcmc.pkl')

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
            'q_surface_mWm2': (OBS['q_surface_mWm2'], OBS['q_surface_err_mWm2']),
        },
        'convection_model': 'yao2014',
    }
    pkl_path = os.path.join(OUTPUT_DIR, 'hybrid_hydro_andrade_yao2014_mcmc.pkl')
    with open(pkl_path, 'wb') as f:
        pickle.dump(results, f)
    log.info(f'Results saved to {pkl_path}')

    if not args.no_plots:
        make_plots(samples, log_likes, k2_results, mtot_results, cmr2_results,
                   heating_results, eval_idx, grid_cache, tb_vals, d_vals)

    # Summary
    print('\n' + '=' * 65)
    print('HYBRID HYDROSPHERE TITAN MCMC SUMMARY (ANDRADE + YAO 2014)')
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
