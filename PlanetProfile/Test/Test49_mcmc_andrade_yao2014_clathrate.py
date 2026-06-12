"""
MCMC exploration: Andrade rheology + Yao 2014 spherical convection, 4km Clathrate Titan (PPTest49).

Perturbation check of Test48 (Path B):
1. Uses PPTest49 (4 km clathrate cap instead of 5 km).
2. Uses Yao et al. (2014) spherical convection for Ice Ih.
3. Goal: evaluate sensitivity of the Path B result (100% Ih heating) to clathrate thickness.

Usage:
  python PlanetProfile/Test/Test49_mcmc_andrade_yao2014_clathrate.py
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
log = logging.getLogger('PPTest49_MCMC_Clath4km')
log.setLevel(logging.INFO)

# TidalPy imports
from TidalPy.RadialSolver import build_rs_input_from_data, radial_solver

# ============================================================
# Configuration
# ============================================================
OUTPUT_DIR = os.path.join(_pp_root, 'PlanetProfile', 'Test', 'mcmc_results', 'Titan', 'Test49_clathrate4km')
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Unique cache for the 4km clathrate case
GRID_CACHE_PATH = os.path.join(OUTPUT_DIR, 'titan_yao2014_clathrate4km_hybrid_hydro_grid_cache.pkl')

# Observational constraints (Petricca et al. 2025, Cassini radio science)
OBS = {
    'Re_k2': 0.608,  'Re_k2_err': 0.048,
    'Im_k2': 0.135,  'Im_k2_err': 0.035,
    'CMR2':  0.343,  'CMR2_err':  0.005,
    'q_surface_mWm2': 10.0,  'q_surface_err_mWm2': 5.0,
}

# Yao 2014 convection constants for Ice Ih
R_GAS = 8.314462  
EACT_IH_JMOL = 60e3  
T_SURF_TITAN_K = 94.0  
K_ICE_WMK = 2.6  
ALPHA_ICE_K = 1.56e-4  
RHO_ICE_KGM3 = 917.0  
KAPPA_ICE_M2S = 1.47e-6  
G_TITAN_MS2 = 1.352  

TITAN_RHO_GCM3     = 1.880;     TITAN_RHO_ERR_GCM3 = 0.004   
TITAN_R_KM         = 2575.5;    TITAN_R_ERR_KM     = 2.0      
TITAN_R_M          = TITAN_R_KM * 1e3
MTOT_OBS = (4.0/3.0) * np.pi * TITAN_R_M**3 * (TITAN_RHO_GCM3 * 1e3)   
MTOT_ERR = MTOT_OBS * np.sqrt((TITAN_RHO_ERR_GCM3/TITAN_RHO_GCM3)**2
                               + (3*TITAN_R_ERR_KM/TITAN_R_KM)**2)        
TITAN_M_REF_KG = 1.3452e23

# Grid extents (matching Test48)
TB_MIN, TB_MAX              = 252.0, 270.0
TB_FINE_STEP, TB_FINE_MAX   = 0.02,  253.0   
TB_COARSE_STEP              = 0.1             
D_MIN,  D_MAX,  D_STEP      = 50.0,  800.0, 25.0   

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
        log.info(f'Loaded {n_pts} cached grid points from {GRID_CACHE_PATH}')
        if n_pts < n_expected and not build_complete:
            log.info(f'Grid incomplete ({n_pts}/{n_expected}). Resuming build...')
            cache_data = build_hybrid_hydrosphere_grid(
                'PlanetProfile.Test.PPTest49',
                tb_grid,
                d_grid,
                GRID_CACHE_PATH,
                rheology='maxwell',
                convection_model='yao2014',
                force_rebuild=False,
            )
    else:
        log.info(f'Building grid for PPTest49 (4km Clathrate, {n_expected} points)...')
        cache_data = build_hybrid_hydrosphere_grid(
            'PlanetProfile.Test.PPTest49',
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
    return grid_cache, tb_vals, d_vals


# ============================================================
# Step 2: Forward model
# ============================================================

def forward_model(theta, grid_cache, tb_vals, d_vals, return_heating=False):
    alpha, log10_zeta, log10_eta_Ih, log10_eta_HP, log10_eta_sil, Tb_K, D_hydro_km, rho_sil, rho_core, f_core = theta
    eta_Ih  = 10 ** log10_eta_Ih
    eta_HP  = 10 ** log10_eta_HP
    eta_sil = 10 ** log10_eta_sil

    idx_tb = int(np.argmin(np.abs(tb_vals - Tb_K)))
    idx_d  = int(np.argmin(np.abs(d_vals  - D_hydro_km)))
    data = grid_cache.get((float(tb_vals[idx_tb]), float(d_vals[idx_d])))
    if data is None:
        return np.nan, np.nan, np.nan, np.nan, {}

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

    layer_types_use   = list(data['layer_types'])
    region_phases_use = list(data['region_phases'])
    if data.get('D_ocean_km', 1.0) < 0.5:
        layer_types_use   = ['solid' if lt == 'liquid' else lt for lt in layer_types_use]
        region_phases_use = ['Ih' if rp == '0' else rp for rp in region_phases_use]

    if abs(Mtot_kg - MTOT_OBS) > 5 * MTOT_ERR:
        return np.nan, np.nan, Mtot_kg, CMR2, {}

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

    from PlanetProfile.Inference.mcmc_common import apply_arrhenius_ih
    apply_arrhenius_ih(eta_mod, phases, ci, n_layers,
                       T_K_profile=data.get('T_K'),
                       Tb_K=Tb_K,
                       E_act_J_mol=EACT_IH_JMOL, R_gas=R_GAS)

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
        Im_k2 = np.abs(k2.imag)

        perPhase_W = {}
        if return_heating and data['eccentricity'] > 0:
            from PlanetProfile.Inference.mcmc_common import compute_per_phase_heating
            perPhase_W = compute_per_phase_heating(result, data)

        return Re_k2, Im_k2, Mtot_kg, CMR2, perPhase_W

    except Exception as exc:
        log.debug(f'TidalPy failed: {exc}')
        return np.nan, np.nan, Mtot_kg, CMR2, {}


# ============================================================
# Step 3: Log-likelihood
# ============================================================

def yao_heat_flux_mWm2(Tb_K, D_iceIh_km, eta_Ih_Pas, R_body_m):
    D_m = D_iceIh_km * 1e3
    if D_m < 1e3: return 0.0
    rBot_m = R_body_m - D_m
    f = rBot_m / R_body_m
    DeltaT = Tb_K - T_SURF_TITAN_K
    if DeltaT <= 0: return 0.0
    gamma = EACT_IH_JMOL / (R_GAS * Tb_K)
    theta_m = 1.0 - 1.23 / (gamma * f**1.5)
    if theta_m <= 0 or theta_m >= 1: return 0.0
    T_m = T_SURF_TITAN_K + theta_m * DeltaT
    DeltaT_v = R_GAS * T_m**2 / EACT_IH_JMOL
    eta_m = eta_Ih_Pas * np.exp(EACT_IH_JMOL / (R_GAS * T_m) - EACT_IH_JMOL / (R_GAS * Tb_K))
    Ra_fullDT = (ALPHA_ICE_K * RHO_ICE_KGM3 * G_TITAN_MS2 * DeltaT * D_m**3 / (eta_Ih_Pas * KAPPA_ICE_M2S))
    Ra_crit = 20.9 * gamma**4
    if Ra_fullDT < Ra_crit: return 0.0
    Ra_m = (ALPHA_ICE_K * RHO_ICE_KGM3 * G_TITAN_MS2 * DeltaT_v * D_m**3 / (eta_m * KAPPA_ICE_M2S))
    Phi_c = K_ICE_WMK * DeltaT / D_m
    q_bot = (1.46 * Ra_m**0.27 / f**1.78 * (DeltaT_v / DeltaT)**1.21 * Phi_c)
    return q_bot * (rBot_m / R_body_m)**2 * 1e3

def log_likelihood(theta, grid_cache, tb_vals, d_vals):
    alpha, log10_zeta, log10_eta_Ih, log10_eta_HP, log10_eta_sil, Tb_K, D_hydro_km, rho_sil, rho_core, f_core = theta
    
    if rho_core < rho_sil: return -1e30

    re, im, mtot, cmr2, _ = forward_model(theta, grid_cache, tb_vals, d_vals)
    if np.isnan(re): return -1e30

    idx_tb = int(np.argmin(np.abs(tb_vals - Tb_K)))
    idx_d  = int(np.argmin(np.abs(d_vals  - D_hydro_km)))
    data = grid_cache.get((float(tb_vals[idx_tb]), float(d_vals[idx_d])))
    D_ih = data['D_iceIh_km']
    q_yao = yao_heat_flux_mWm2(Tb_K, D_ih, 10**log10_eta_Ih, data['R_body_m'])

    chi2 = 0.0
    chi2 += ((re - OBS['Re_k2']) / OBS['Re_k2_err'])**2
    chi2 += ((im - OBS['Im_k2']) / OBS['Im_k2_err'])**2
    chi2 += ((cmr2 - OBS['CMR2']) / OBS['CMR2_err'])**2
    chi2 += ((mtot - MTOT_OBS) / MTOT_ERR)**2
    chi2 += ((q_yao - OBS['q_surface_mWm2']) / OBS['q_surface_err_mWm2'])**2
    
    return -0.5 * chi2

# ============================================================
# Step 4: Run MCMC (pocoMC)
# ============================================================

def run_mcmc(grid_cache, tb_vals, d_vals):
    import pocomc as pc
    from scipy.stats import uniform
    from PlanetProfile.Inference.mcmc_common import run_pocomc_sampler

    prior = pc.Prior([
        uniform(loc=0.15, scale=0.30),  # alpha: [0.15, 0.45]
        uniform(loc=-3.0, scale=5.0),   # log10_zeta: [-3, 2]
        uniform(loc=10.0, scale=6.0),   # log10_eta_Ih: [10, 16]
        uniform(loc=10.0, scale=6.0),   # log10_eta_HP: [10, 16]
        uniform(loc=18.0, scale=4.0),   # log10_eta_sil: [18, 22]
        uniform(loc=TB_MIN, scale=TB_MAX - TB_MIN),
        uniform(loc=D_MIN, scale=D_MAX - D_MIN),
        uniform(loc=1800.0, scale=1700.0), # rho_sil: [1800, 3500]
        uniform(loc=1800.0, scale=6700.0), # rho_core: [1800, 8500]
        uniform(loc=0.0, scale=0.8),    # f_core: [0, 0.8]
    ])

    def _log_like(theta):
        return log_likelihood(theta, grid_cache, tb_vals, d_vals)

    samples, log_likes, sampler = run_pocomc_sampler(
        prior, _log_like, n_effective=N_EFF, random_state=RANDOM_STATE
    )
    return samples, log_likes, sampler


# ============================================================
# Step 5: Re-evaluate heating
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
# Main
# ============================================================

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PPTest49 Andrade + Yao2014 4km-Clathrate MCMC')
    parser.add_argument('--rebuild-grid', action='store_true', help='Force rebuild of structural grid')
    parser.add_argument('--replot', type=str, help='Re-plot existing results from pkl')
    parser.add_argument('--no-plots', action='store_true', help='Skip plot generation')
    args = parser.parse_args()

    grid_cache, tb_vals, d_vals = build_or_load_grid(force_rebuild=args.rebuild_grid)

    # Derive output names
    pkl_name = 'test49_clathrate4km_mcmc_results.pkl'
    out_pkl = os.path.join(OUTPUT_DIR, pkl_name)

    if args.replot:
        with open(args.replot, 'rb') as f:
            res = pickle.load(f)
        samples = res['samples']
        log_likes = res['log_likes']
    else:
        samples, log_likes, sampler = run_mcmc(grid_cache, tb_vals, d_vals)
        with open(out_pkl, 'wb') as f:
            pickle.dump({'samples': samples, 'log_likes': log_likes, 'params': PARAM_NAMES}, f)
        log.info(f'Results saved to {out_pkl}')

    if not args.no_plots:
        log.info('Evaluating posterior for plotting...')
        eval_idx, k2_res, mtot_res, cmr2_res, heat_res = evaluate_heating(
            samples, grid_cache, tb_vals, d_vals
        )
        
        from PlanetProfile.Inference import mcmc_plots as mp
        out = lambda name: os.path.join(OUTPUT_DIR, f'titan_test49_clathrate4km_{name}.png')
        title = 'Titan (Yao 2014, 4km Clathrate)'
        
        mp.plot_corner(samples, PARAM_LABELS, title + ' — Posterior', out('corner'))
        
        totals = np.array([sum(h.values()) if h else 1e-30 for h in heat_res])
        f_sil = np.array([h.get('Sil', 0.0) for h in heat_res]) / np.where(totals > 0, totals, 1)
        mp.plot_k2_scatter_by(
            k2_res, color_values=f_sil, colorbar_label='Silicate heating fraction',
            obs_re=OBS['Re_k2'], obs_im=OBS['Im_k2'],
            obs_re_err=OBS['Re_k2_err'], obs_im_err=OBS['Im_k2_err'],
            title=title + r' — $k_2$ Posterior', output_path=out('k2_scatter')
        )
        
        # Additional plots can be added here following Test48 pattern
        log.info('Plotting complete.')

    log.info('Done.')
