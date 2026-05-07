"""
MCMC exploration: Andrade rheology, all-ice Titan with Yao 2014 spherical Ih convection (PPTest50).

Explores Titan's interior structure under the hypothesis of NO liquid ocean —
the entire hydrosphere is ice (Ih + III + V + VI), as in Petricca et al. (2025).

With NO_OCEAN_EXCEPT_INNER_ICES=True and Fe_CORE=False, the ice structure is
self-consistently determined by PP and barely varies with Tb_K (< 1 km over the
entire valid range 240–251 K). The all-ice hypothesis produces a UNIQUE
structure: D_hsphere ≈ 494 km, rho_sil ≈ 2555 kg/m³, CMR2 ≈ 0.343.

The MCMC samples rheology only (7D) at this fixed structure. Per-phase HP ice
viscosities are sampled independently since each phase has distinct shear moduli
and dissipation behaviour (Petricca et al. 2025 find ~1e12 Pa s for HP ices).

Parameter space (7D):
  alpha:           Andrade exponent                [0.15, 0.45]
  log10(zeta):     Andrade zeta (Pa s^alpha)       [-3, 2]
  log10(eta_Ih):   Ice Ih viscosity (Pa s)         [12, 16]
  log10(eta_III):  Ice III viscosity (Pa s)        [10, 16]
  log10(eta_V):    Ice V viscosity (Pa s)          [10, 16]
  log10(eta_VI):   Ice VI viscosity (Pa s)         [10, 16]
  log10(eta_sil):  Silicate viscosity (Pa s)       [18, 22]

Fixed structure (Tb_K=250 K, computed by PP):
  D_iceIh  ≈ 156 km
  D_iceIII ≈ 83 km
  D_iceV   ≈ 155 km
  D_iceVI  ≈ 69 km
  D_total  ≈ 494 km
  rho_sil  ≈ 2555 kg/m³
  CMR2     ≈ 0.343

Observational constraints (Petricca et al. 2025, Cassini radio science):
  Re(k2)   = 0.608 +/- 0.048
  |Im(k2)| = 0.135 +/- 0.035

Usage:
  mamba activate PPcl
  python PlanetProfile/Test/Test46_mcmc_allice.py
  python PlanetProfile/Test/Test46_mcmc_allice.py --build-structure
  python PlanetProfile/Test/Test46_mcmc_allice.py --replot
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

# TidalPy imports
from TidalPy.RadialSolver import radial_solver, build_rs_input_from_data
from TidalPy.tides.multilayer.heating import calc_radial_volumetric_tidal_heating_from_rs_solution
from TidalPy.rheology import Andrade, Elastic

# ============================================================
# Configuration
# ============================================================
OUTPUT_DIR = os.path.join(_pp_root, 'PlanetProfile', 'Test', 'mcmc_results')
os.makedirs(OUTPUT_DIR, exist_ok=True)

STRUCTURE_CACHE_PATH = os.path.join(OUTPUT_DIR, 'titan_allice_yao2014_structure_cache.pkl')

# Observational constraints (Petricca et al. 2025, Cassini radio science)
OBS = {
    'Re_k2': 0.608,  'Re_k2_err': 0.048,
    'Im_k2': 0.135,  'Im_k2_err': 0.035,
    'CMR2':  0.343,  'CMR2_err':  0.001,
}

# Titan reference
TITAN_R_M = 2574.73e3
TITAN_M_KG = 1.3452e23

# Fixed structural Tb_K (all-ice structure is insensitive to Tb_K within valid range)
STRUCTURE_TB_K = 250.0

# MCMC settings
N_EFF = 500
RANDOM_STATE = 42
N_REEVAL = 500

PARAM_NAMES = ['alpha', 'log10_zeta', 'log10_eta_Ih', 'log10_eta_III',
               'log10_eta_V', 'log10_eta_VI', 'log10_eta_sil']
PARAM_LABELS = [
    r'$\alpha$',
    r'$\log_{10}\zeta$',
    r'$\log_{10}\eta_\mathrm{Ih}$',
    r'$\log_{10}\eta_\mathrm{III}$',
    r'$\log_{10}\eta_\mathrm{V}$',
    r'$\log_{10}\eta_\mathrm{VI}$',
    r'$\log_{10}\eta_\mathrm{sil}$',
]
N_DIM = 7


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
    """Run PlanetProfile for one Tb_K point and extract arrays.

    Tb_K controls the thermal boundary at the base of ice Ih, which determines
    the entire HP ice structure below. With Fe_CORE=False, PP determines rho_sil
    self-consistently from mass balance.
    """
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

    # Set ice Ih base temperature (controls entire ice structure)
    Planet.Bulk.Tb_K = Tb_K

    configParams.Gravity.backend = 'tidalpy'
    configParams.CALC_NEW = True
    configParams.CALC_NEW_GRAVITY = True
    configParams.NO_SAVEFILE = True
    configParams.SKIP_PLOTS = True

    Planet, Params = RunPP(Planet, configParams)
    Params.CALC_NEW_GRAVITY = True
    Planet, Params = SetupGravity(Planet, Params)

    # Extract structural arrays (same pattern as Test41)
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

    # T profile (needed for Arrhenius Ih viscosity). Reduced.T_K is inside-out.
    try:
        T_K_reduced = np.asarray(Planet.Reduced.T_K, dtype=np.float64)
        T_K = np.flipud(T_K_reduced) if T_K_reduced.size == r_m.size else np.full(r_m.size, np.nan)
    except AttributeError:
        T_K = np.full(r_m.size, np.nan)
        log.warning('Planet.Reduced.T_K not available — Arrhenius Ih viscosity will be skipped.')

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

    # Layer boundaries
    changeIndices = np.max(Planet.Reduced.changeIndices) - np.flipud(Planet.Reduced.changeIndices)
    n_layers = len(changeIndices) - 1

    # Thin-layer padding
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
        new_r, new_rho, new_K, new_mu, new_eta, new_phases, new_bv = \
            [], [], [], [], [], [], []
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
                new_ci.append(new_ci[-1] + MIN_POINTS)
            else:
                new_r.append(r_m[s:e])
                new_rho.append(rho[s:e])
                new_K.append(K_Pa[s:e])
                new_mu.append(mu_Pa[s:e])
                new_eta.append(eta_Pa_base[s:e])
                new_phases.append(phases[s:e])
                new_bv.append(bulk_visc[s:e])
                new_ci.append(new_ci[-1] + (e - s))

        r_m = np.concatenate(new_r)
        rho = np.concatenate(new_rho)
        K_Pa = np.concatenate(new_K)
        mu_Pa = np.concatenate(new_mu)
        eta_Pa_base = np.concatenate(new_eta)
        phases = np.concatenate(new_phases)
        bulk_visc = np.concatenate(new_bv)
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

    # Layer thicknesses for diagnostics
    D_iceIh_km = 0.0
    D_iceIII_km = 0.0
    D_iceV_km = 0.0
    D_iceVI_km = 0.0
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

    # CMR2 from PP
    CMR2_pp = Planet.Bulk.Cmeasured if hasattr(Planet.Bulk, 'Cmeasured') else np.nan
    try:
        CMR2_pp = float(Planet.CMR2mean)
    except (AttributeError, TypeError):
        pass

    # Total hydrosphere thickness (derived, not input)
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


def build_or_load_structure(force_rebuild=False):
    """Build or load the single all-ice structure from PP."""
    if not force_rebuild and os.path.exists(STRUCTURE_CACHE_PATH):
        try:
            with open(STRUCTURE_CACHE_PATH, 'rb') as f:
                data = pickle.load(f)
            if data is not None and 'r_m' in data:
                log.info(f'Loaded cached structure from {STRUCTURE_CACHE_PATH}')
                log.info(f'  D_hsphere={data["D_hsphere_km"]:.1f} km, '
                         f'rho_sil={data["rhoSil_kgm3"]:.0f} kg/m³, '
                         f'CMR2={data["CMR2"]:.4f}')
                return data
        except (EOFError, pickle.UnpicklingError):
            log.warning('Cache corrupted, rebuilding.')

    log.info(f'Building all-ice structure at Tb_K={STRUCTURE_TB_K:.1f} K ...')
    data = _build_single_structure(STRUCTURE_TB_K)
    _save_cache(data, STRUCTURE_CACHE_PATH)
    log.info(f'  D_hsphere={data["D_hsphere_km"]:.1f} km, '
             f'D_Ih={data["D_iceIh_km"]:.1f}, D_III={data["D_iceIII_km"]:.1f}, '
             f'D_V={data["D_iceV_km"]:.1f}, D_VI={data["D_iceVI_km"]:.1f} km')
    log.info(f'  rho_sil={data["rhoSil_kgm3"]:.0f} kg/m³, CMR2={data["CMR2"]:.4f}')
    log.info(f'  Saved to {STRUCTURE_CACHE_PATH}')
    return data


# ============================================================
# Step 2: Forward model
# ============================================================

def forward_model(theta, structure, return_heating=False):
    """Compute k2 and per-phase heating for given rheology parameters.

    Args:
        theta: [alpha, log10_zeta, log10_eta_Ih, log10_eta_III,
                log10_eta_V, log10_eta_VI, log10_eta_sil]
        structure: dict from build_or_load_structure()
    Returns:
        (Re_k2, Im_k2, perPhase_W)
    """
    alpha, log10_zeta, log10_eta_Ih, log10_eta_III, log10_eta_V, \
        log10_eta_VI, log10_eta_sil = theta
    eta_Ih = 10 ** log10_eta_Ih
    eta_III = 10 ** log10_eta_III
    eta_V = 10 ** log10_eta_V
    eta_VI = 10 ** log10_eta_VI
    eta_sil = 10 ** log10_eta_sil

    data = structure

    # Apply per-phase viscosity overrides
    eta_mod = data['eta_Pa_base'].copy()
    phases = data['phases']
    ci = data['changeIndices']
    n_layers = data['n_layers']

    for i in range(n_layers):
        s, e = ci[i], ci[i + 1]
        ph = int(phases[s])
        if ph == 1:
            eta_mod[s:e] = eta_Ih
        elif ph == 3:
            eta_mod[s:e] = eta_III
        elif ph == 5:
            eta_mod[s:e] = eta_V
        elif ph == 6:
            eta_mod[s:e] = eta_VI
        elif 50 <= ph < 100:
            eta_mod[s:e] = eta_sil

    # Apply Arrhenius T-dependence to Ice Ih. Sampled eta_Ih is the BASAL
    # viscosity at T=Tb; cold stagnant lid is much stiffer under Yao scaling.
    # η(T) = η_Ih * exp(E/R * (1/T - 1/Tb))  with E = 60 kJ/mol (Yao Table 5).
    Tb_K_val = data.get('Tb_K', STRUCTURE_TB_K)
    T_K_prof = data.get('T_K')
    EACT_IH_JMOL = 60e3
    R_GAS = 8.314462
    if T_K_prof is not None:
        T_K_arr = np.asarray(T_K_prof)
        if T_K_arr.shape == eta_mod.shape:
            for i in range(n_layers):
                s, e = ci[i], ci[i + 1]
                if int(phases[min(s, len(phases) - 1)]) == 1:
                    T_layer = T_K_arr[s:e]
                    if np.all(np.isfinite(T_layer)) and np.all(T_layer > 0):
                        exponent = (EACT_IH_JMOL / R_GAS) * (1.0 / T_layer - 1.0 / Tb_K_val)
                        eta_mod[s:e] *= np.exp(exponent)

    # Build Andrade rheology
    zeta_pa = 10 ** log10_zeta
    zeta_tp = zeta_pa ** (1.0 / alpha)
    shear = []
    for rp in data['region_phases']:
        base = rp.replace('_conv', '')
        shear.append(Elastic() if base in ('0', 'Clath') else Andrade(args=(alpha, zeta_tp)))
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
            return np.nan, np.nan, {}

        k2 = complex(result.k)
        Re_k2 = k2.real
        Im_k2 = k2.imag

        perPhase_W = {}
        if return_heating and data['eccentricity'] > 0:
            hp = calc_radial_volumetric_tidal_heating_from_rs_solution(
                data['eccentricity'], data['omega'], data['a_m'],
                data['host_mass'], result, perform_checks=False
            )
            rr = np.asarray(result.radius_array)
            r_m_local = data['r_m']
            for i_layer in range(n_layers):
                s_idx, e_idx = ci[i_layer], ci[i_layer + 1]
                ph = int(phases[s_idx])
                phase_str = {0: '0', 1: 'Ih', 3: 'III', 5: 'V', 6: 'VI'}.get(
                    ph, PhaseConv(ph, liq='0'))
                r_lo, r_hi = r_m_local[s_idx], r_m_local[e_idx - 1]
                mask = (rr >= r_lo - 1.0) & (rr <= r_hi + 1.0)
                if np.any(mask):
                    lr, lh = rr[mask], hp[mask]
                    if len(lr) > 1:
                        power = np.trapz(lh * 4.0 * np.pi * lr**2, lr)
                    else:
                        power = lh[0] * (4.0/3.0) * np.pi * (r_hi**3 - r_lo**3)
                    perPhase_W[phase_str] = perPhase_W.get(phase_str, 0) + power

        return Re_k2, Im_k2, perPhase_W

    except Exception as exc:
        log.debug(f'TidalPy failed: {exc}')
        return np.nan, np.nan, {}


# ============================================================
# Step 3: Log-likelihood
# ============================================================

def log_likelihood(theta, structure):
    """Gaussian log-likelihood on k2 (CMR2 is fixed at ~0.343, matching obs)."""
    Re_k2, Im_k2, _ = forward_model(theta, structure)
    if np.isnan(Re_k2):
        return -1e30
    chi2 = (
        ((Re_k2 - OBS['Re_k2']) / OBS['Re_k2_err']) ** 2 +
        ((abs(Im_k2) - OBS['Im_k2']) / OBS['Im_k2_err']) ** 2
    )
    return -0.5 * chi2


# ============================================================
# Step 4: MCMC
# ============================================================

def run_mcmc(structure):
    import pocomc as pc
    from scipy.stats import uniform

    log.info(f'Starting pocoMC MCMC ({N_DIM}D, n_eff={N_EFF})')

    prior = pc.Prior([
        uniform(loc=0.15,  scale=0.30),                   # alpha: [0.15, 0.45]
        uniform(loc=-3.0,  scale=5.0),                    # log10_zeta: [-3, 2]
        uniform(loc=12.0,  scale=4.0),                    # log10_eta_Ih: [12, 16]
        uniform(loc=12.0,  scale=4.0),                    # log10_eta_III: [12, 16] (brackets HP Maxwell peak)
        uniform(loc=12.0,  scale=4.0),                    # log10_eta_V: [12, 16]
        uniform(loc=12.0,  scale=4.0),                    # log10_eta_VI: [12, 16]
        uniform(loc=18.0,  scale=4.0),                    # log10_eta_sil: [18, 22]
    ])

    def _log_like(theta):
        return log_likelihood(theta, structure)

    t0 = time.time()
    sampler = pc.Sampler(
        prior=prior,
        likelihood=_log_like,
        n_effective=N_EFF,
        random_state=RANDOM_STATE,
    )
    sampler.run()
    elapsed = time.time() - t0
    log.info(f'MCMC completed in {elapsed/60:.1f} min')

    raw_x = sampler.results['x']
    raw_logl = sampler.results['logl']
    # pocoMC returns (n_iters, n_walkers, n_dim); flatten to 2D
    if raw_x.ndim == 3:
        samples = raw_x.reshape(-1, raw_x.shape[-1])
        log_prob = raw_logl.reshape(-1)
    else:
        samples = raw_x
        log_prob = raw_logl
    return samples, log_prob


# ============================================================
# Step 5: Plotting
# ============================================================

def make_plots(samples, log_prob, structure):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    PHASE_COLORS = {
        'Clath': '#D4F1F9', 'Ih': '#AEE1F8',
        'III': '#C97BAE', 'V': '#9B59B6', 'VI': '#6C3483',
        'Sil': '#C8A96E',
    }
    PHASE_LABELS = {
        'Clath': 'Clathrate', 'Ih': 'Ice Ih',
        'III': 'Ice III', 'V': 'Ice V', 'VI': 'Ice VI',
        'Sil': 'Silicate',
    }

    # --- Corner plot ---
    try:
        import corner
        fig = corner.corner(samples, labels=PARAM_LABELS,
                            quantiles=[0.16, 0.5, 0.84], show_titles=True)
        out_corner = os.path.join(OUTPUT_DIR, 'allice_yao2014_andrade_corner.png')
        fig.savefig(out_corner, dpi=150, bbox_inches='tight')
        plt.close(fig)
        log.info(f'Corner plot saved: {out_corner}')
    except ImportError:
        log.warning('corner not installed, skipping corner plot')

    # --- Re-evaluate subset for heating ---
    n_eval = min(N_REEVAL, len(samples))
    rng = np.random.default_rng(42)
    eval_idx = rng.choice(len(samples), size=n_eval, replace=False)

    heating_results = []
    k2_results = []
    for i in eval_idx:
        Re_k2_i, Im_k2_i, heat_i = forward_model(
            samples[i], structure, return_heating=True)
        heating_results.append(heat_i)
        k2_results.append((Re_k2_i, Im_k2_i))

    # Posterior arrays used by multiple plots
    re_arr = np.array([r[0] for r in k2_results])
    im_arr = np.array([abs(r[1]) for r in k2_results])
    valid = np.isfinite(re_arr)

    totals = np.array([sum(h.values()) if h else 1e-30 for h in heating_results])
    safe_tot = np.where(totals > 1e-30, totals, 1e-30)
    heating_fracs = {ph: np.array([h.get(ph, 0.0) for h in heating_results]) / safe_tot
                     for ph in ['Ih', 'III', 'V', 'VI', 'Sil', 'Clath']}
    f_sil = heating_fracs['Sil']

    # --- k2 scatter colored by f_sil with Petricca 1σ/2σ ellipses ---
    from matplotlib.patches import Ellipse
    fig, ax = plt.subplots(figsize=(8, 6))
    sc = ax.scatter(re_arr[valid], im_arr[valid], c=f_sil[valid], cmap='RdYlBu_r',
                    s=10, alpha=0.6, vmin=0, vmax=1)
    plt.colorbar(sc, ax=ax, label='Silicate heating fraction')
    for nσ, ls, lw, lbl in [(1, '--', 2, r'obs 1$\sigma$'), (2, ':', 1, r'obs 2$\sigma$')]:
        ax.add_patch(Ellipse((OBS['Re_k2'], OBS['Im_k2']),
                             2*nσ*OBS['Re_k2_err'], 2*nσ*OBS['Im_k2_err'],
                             fill=False, color='red', linewidth=lw, linestyle=ls, label=lbl))
    ax.set_xlabel(r'$\mathrm{Re}(k_2)$', fontsize=12)
    ax.set_ylabel(r'$|\mathrm{Im}(k_2)|$', fontsize=12)
    ax.set_title(f'All-Ice Titan: $k_2$ Posterior (D_hsphere={structure["D_hsphere_km"]:.0f} km, '
                 f'CMR2={structure["CMR2"]:.4f})', fontsize=11)
    ax.legend(loc='best', fontsize=9)
    out_k2 = os.path.join(OUTPUT_DIR, 'allice_yao2014_andrade_k2_scatter.png')
    fig.savefig(out_k2, dpi=200, bbox_inches='tight')
    plt.close(fig)
    log.info(f'k2 scatter saved: {out_k2}')

    # --- Combined heating figure: per-parameter scatter + cumulative fraction bar ---
    # Layout: 2×4 grid of per-parameter heating-fraction scatters (7 params fit in 8 slots)
    # + bottom full-width cumulative stacked bar (sorted by f_sil ASC).
    import matplotlib.gridspec as gridspec
    eval_samples = samples[eval_idx]
    fig = plt.figure(figsize=(20, 13))
    gs = gridspec.GridSpec(3, 4, height_ratios=[1, 1, 1.2], hspace=0.35, wspace=0.28)
    scatter_axes = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(4)]
    ax_cum = fig.add_subplot(gs[2, :])

    phase_colors_pp = {'Ih': 'C0', 'III': 'C1', 'V': 'C2', 'VI': 'C3', 'Sil': 'C4'}
    phases_scatter = ['Ih', 'III', 'V', 'VI', 'Sil']

    for ip, (pname, plabel) in enumerate(zip(PARAM_NAMES, PARAM_LABELS)):
        ax = scatter_axes[ip]
        x = eval_samples[:, ip]
        for ph in phases_scatter:
            ax.scatter(x, heating_fracs[ph], s=4, alpha=0.35,
                       color=phase_colors_pp[ph], label=PHASE_LABELS.get(ph, ph))
        ax.set_xlabel(plabel)
        ax.set_ylabel('Heating fraction')
        ax.set_ylim(-0.05, 1.05)
        if ip == 0:
            ax.legend(fontsize=8, loc='center right')

    # 8th slot: legend/text panel with structure info
    ax_info = scatter_axes[7]
    ax_info.axis('off')
    ax_info.text(0.02, 0.95,
                 'Fixed structure (PPTest50):\n'
                 f'  D_hsphere = {structure["D_hsphere_km"]:.1f} km\n'
                 f'  D_iceIh   = {structure["D_iceIh_km"]:.1f} km\n'
                 f'  D_iceIII  = {structure["D_iceIII_km"]:.1f} km\n'
                 f'  D_iceV    = {structure["D_iceV_km"]:.1f} km\n'
                 f'  D_iceVI   = {structure["D_iceVI_km"]:.1f} km\n'
                 f'  CMR2      = {structure["CMR2"]:.4f}\n'
                 f'  ρ_sil     = {structure["rhoSil_kgm3"]:.0f} kg/m³\n'
                 f'  Tb_K      = {structure.get("Tb_K", STRUCTURE_TB_K):.1f} K\n'
                 '\nAndrade + Yao 2014 Ih conv +\nKalousova HP + Arrhenius Ih',
                 family='monospace', fontsize=9, va='top', transform=ax_info.transAxes)

    # --- Cumulative stacked heating fractions, sorted by (f_sil ASC) ---
    PHASE_COLORS_STACK = {
        'Ih':  '#AEE1F8',
        'III': '#C97BAE', 'V':   '#9B59B6', 'VI': '#6C3483',
        'Sil': '#C8A96E',
    }
    stack_phases = ['Ih', 'III', 'V', 'VI', 'Sil']

    order = np.argsort(f_sil)
    n_models = len(order)
    x = np.arange(n_models)
    bottom = np.zeros(n_models)
    for ph in stack_phases:
        h_vals = heating_fracs[ph][order]
        ax_cum.bar(x, h_vals, bottom=bottom, width=1.0,
                   color=PHASE_COLORS_STACK[ph], edgecolor='none',
                   label=PHASE_LABELS.get(ph, ph))
        bottom += h_vals

    ax_cum.set_ylim(0, 1.0)
    ax_cum.set_xlim(0, n_models - 1)
    ax_cum.set_ylabel('Cumulative heating fraction', fontsize=11)
    ax_cum.set_title('Per-model heating fractions  (models sorted ascending by $f_\\mathrm{sil}$)',
                     fontsize=11)
    n_ticks = 8
    tick_positions = np.linspace(0, n_models - 1, n_ticks, dtype=int)
    tick_labels = [f'{f_sil[order][i]:.2f}' for i in tick_positions]
    ax_cum.set_xticks(tick_positions)
    ax_cum.set_xticklabels(tick_labels, fontsize=8)
    ax_cum.set_xlabel(r'$f_\mathrm{sil}$   [models sorted ascending]', fontsize=11)
    handles, lbls = ax_cum.get_legend_handles_labels()
    ax_cum.legend(handles[::-1], lbls[::-1],
                  loc='lower left', fontsize=9, ncol=5)

    fig.suptitle(f'All-Ice Titan Yao 2014 MCMC — Heating vs Parameters (top) '
                 f'and Per-Model Fractions (bottom)   [CMR2={structure["CMR2"]:.4f}]',
                 fontsize=13)
    out_heat = os.path.join(OUTPUT_DIR, 'allice_yao2014_andrade_heating.png')
    fig.savefig(out_heat, dpi=150, bbox_inches='tight')
    plt.close(fig)
    log.info(f'Heating plot saved: {out_heat}')

    # --- Summary statistics ---
    log.info('=== Posterior Summary (All-Ice, 5D Rheology) ===')
    log.info(f'  Fixed structure: D_hsphere={structure["D_hsphere_km"]:.1f} km, '
             f'CMR2={structure["CMR2"]:.4f}, rho_sil={structure["rhoSil_kgm3"]:.0f} kg/m³')
    log.info(f'  Layers: Ih={structure["D_iceIh_km"]:.1f}, '
             f'III={structure["D_iceIII_km"]:.1f}, '
             f'V={structure["D_iceV_km"]:.1f}, '
             f'VI={structure["D_iceVI_km"]:.1f} km')
    for i, name in enumerate(PARAM_NAMES):
        med = np.median(samples[:, i])
        lo, hi = np.percentile(samples[:, i], [16, 84])
        log.info(f'  {name}: {med:.3f} [{lo:.3f}, {hi:.3f}]')
    # k2 summary
    re_valid = re_arr[valid]
    im_valid = im_arr[valid]
    log.info(f'  Re(k2): {np.median(re_valid):.4f} '
             f'[{np.percentile(re_valid, 16):.4f}, {np.percentile(re_valid, 84):.4f}]')
    log.info(f'  |Im(k2)|: {np.median(im_valid):.4f} '
             f'[{np.percentile(im_valid, 16):.4f}, {np.percentile(im_valid, 84):.4f}]')


# ============================================================
# Main
# ============================================================
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PPTest46 Andrade all-ice MCMC')
    parser.add_argument('--rebuild-structure', action='store_true',
                        help='Force rebuild of structural cache')
    parser.add_argument('--build-structure', action='store_true',
                        help='Build structure then exit without running MCMC')
    parser.add_argument('--no-plots', action='store_true',
                        help='Skip plot generation')
    parser.add_argument('--replot', action='store_true',
                        help='Load existing pkl and regenerate plots')
    args = parser.parse_args()

    # Build / load single all-ice structure
    structure = build_or_load_structure(
        force_rebuild=args.rebuild_structure)

    if args.build_structure:
        log.info('--build-structure: exiting after structure build.')
        sys.exit(0)

    pkl_path = os.path.join(OUTPUT_DIR, 'allice_yao2014_andrade_mcmc_results.pkl')

    if args.replot:
        log.info(f'Loading results from {pkl_path}')
        with open(pkl_path, 'rb') as f:
            results = pickle.load(f)
        samples, log_prob = results['samples'], results['log_prob']
    else:
        samples, log_prob = run_mcmc(structure)
        results = {'samples': samples, 'log_prob': log_prob,
                   'param_names': PARAM_NAMES, 'obs': OBS,
                   'structure_info': {
                       'Tb_K': STRUCTURE_TB_K,
                       'D_hsphere_km': structure['D_hsphere_km'],
                       'CMR2': structure['CMR2'],
                       'rhoSil_kgm3': structure['rhoSil_kgm3'],
                   }}
        _save_cache(results, pkl_path)
        log.info(f'MCMC results saved to {pkl_path}')

    if not args.no_plots:
        make_plots(samples, log_prob, structure)

    log.info('Done.')
