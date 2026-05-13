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

STRUCTURE_CACHE_PATH = os.path.join(OUTPUT_DIR, 'titan_allice_yao2014_structure_grid.pkl')

# Observational constraints (Petricca et al. 2025, Cassini radio science)
OBS = {
    'Re_k2': 0.608,  'Re_k2_err': 0.048,
    'Im_k2': 0.135,  'Im_k2_err': 0.035,
    'CMR2':  0.343,  'CMR2_err':  0.001,
}

# Titan reference
TITAN_R_M = 2574.73e3
TITAN_M_KG = 1.3452e23

# Tb_K grid for Option-2 structural interpolation.  Upper bound = Ih-III-L triple
# point minus ε=0.2 K (grid-resolution-safe for nIceI=200).  Lower bound = ~2 K
# depression simulating solute-driven triple-point lowering (e.g., NaCl ~15 ppt).
# 9 points at 0.246 K spacing gives <0.1% linear-interp error.
# Probe at upper edge: PbI=207.68 MPa, CMR2=0.343, no Ih cell overshoots Tm(P).
STRUCTURE_TB_GRID = np.linspace(249.0, 250.965, 9)
STRUCTURE_TB_K = float(STRUCTURE_TB_GRID[-1])  # back-compat alias for logging

# MCMC settings
N_EFF = 500
RANDOM_STATE = 42
N_REEVAL = 500

PARAM_NAMES = ['alpha', 'log10_zeta', 'log10_eta_Ih', 'log10_eta_III',
               'log10_eta_V', 'log10_eta_VI', 'log10_eta_sil', 'Tb_K']
PARAM_LABELS = [
    r'$\alpha$',
    r'$\log_{10}\zeta$',
    r'$\log_{10}\eta_\mathrm{Ih}$',
    r'$\log_{10}\eta_\mathrm{III}$',
    r'$\log_{10}\eta_\mathrm{V}$',
    r'$\log_{10}\eta_\mathrm{VI}$',
    r'$\log_{10}\eta_\mathrm{sil}$',
    r'$T_b$ (K)',
]
N_DIM = 8


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
    # Disable gravity backend during RunPP: Yao2014 spherical convection produces
    # a thick stagnant lid discretized into very few slices, which TidalPy rejects
    # ("Layer N has 2 slices when at least 5 are required").  We skip the gravity
    # call here and re-run SetupGravity manually after the thin-layer padding step
    # below — identical to the strategy used in hybrid_structure_cache.py line 220.
    configParams.CALC_NEW_GRAVITY = False
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

    # T profile (needed for Arrhenius Ih viscosity).  Planet.T_K is the primary
    # (outside-in) thermal profile populated by RunPP.  Planet.Reduced.T_K does
    # NOT exist on this PlanetStruct, so we interpolate Planet.T_K onto the ALMA
    # r_m grid via ascending-r sort — orientation-agnostic.
    try:
        T_K_primary = np.asarray(Planet.T_K, dtype=np.float64)
        r_m_primary = np.asarray(Planet.r_m[:T_K_primary.size], dtype=np.float64)
        sort_idx = np.argsort(r_m_primary)
        T_K = np.interp(r_m, r_m_primary[sort_idx], T_K_primary[sort_idx])
    except (AttributeError, ValueError) as _exc:
        T_K = np.full(r_m.size, np.nan)
        log.warning(f'Planet.T_K extraction failed ({_exc}) — Arrhenius Ih viscosity will be skipped.')

    # P(r) profile (needed for no-ocean safeguard: compare Ih-cell T against Tm_Ih(P))
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
                # Single-slice layer (Yao 2014 stagnant lid): span from the
                # previous layer's top radius up to this slice's r and fill
                # with constant layer-property values.  T_K is interpolated
                # linearly between the previous slice and this one so the
                # Arrhenius profile respects the real surface→basal gradient.
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

    Returns:
        dict with keys:
            'Tb_K_grid': 1D array of Tb values.
            'structures': list of per-Tb structure dicts (from _build_single_structure).

    All grid structures must share identical `changeIndices`, `layer_types`,
    `n_layers`, `region_phases` — otherwise linear interpolation between
    bracketing Tb grid points in `_interp_structure` is ill-defined.  PP's
    fixed nClath/nIceI/nSilMax discretization should guarantee this; the
    check is included for safety and produces a helpful error if it fails.
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
                    f'than the first grid point. Option-2 interpolation requires '
                    f'identical cell layouts across all Tb samples.  Likely cause: '
                    f'thin-layer padding logic produced different minimum-point '
                    f'expansion at this Tb.  Narrow the Tb grid span or investigate.'
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


# Keys whose values are 1D numpy arrays and are linearly interpolated between Tb grid points.
_INTERP_ARRAY_KEYS = ('r_m', 'rho', 'K_Pa', 'mu_Pa', 'eta_Pa_base',
                      'T_K', 'P_MPa', 'bulk_visc')
# Keys whose values are scalars (float) and are linearly interpolated.
_INTERP_SCALAR_KEYS = ('Tb_K', 'CMR2', 'rhoSil_kgm3', 'D_hsphere_km',
                       'D_iceIh_km', 'D_iceIII_km', 'D_iceV_km', 'D_iceVI_km')


def _interp_structure(Tb_sampled, grid):
    """Linearly interpolate structure arrays between bracketing Tb grid points.

    Non-interpolated keys (phases, changeIndices, layer_types, region_phases,
    n_layers, omega, eccentricity, host_mass, a_m, R_body_m, Mtot_kg) are taken
    unchanged from the lower bracket — they are assumed grid-invariant by design.
    """
    Tb_grid = grid['Tb_K_grid']
    structs = grid['structures']
    Tb_clamped = float(np.clip(Tb_sampled, Tb_grid[0], Tb_grid[-1]))
    # Find bracketing indices
    j = int(np.searchsorted(Tb_grid, Tb_clamped))
    if j <= 0:
        return {k: (v.copy() if isinstance(v, np.ndarray) else v)
                for k, v in structs[0].items()}
    if j >= len(Tb_grid):
        return {k: (v.copy() if isinstance(v, np.ndarray) else v)
                for k, v in structs[-1].items()}
    t0, t1 = Tb_grid[j - 1], Tb_grid[j]
    if t1 == t0:
        w = 0.0
    else:
        w = (Tb_clamped - t0) / (t1 - t0)
    s0, s1 = structs[j - 1], structs[j]
    out = {}
    for k in _INTERP_ARRAY_KEYS:
        out[k] = (1.0 - w) * s0[k] + w * s1[k]
    for k in _INTERP_SCALAR_KEYS:
        out[k] = (1.0 - w) * s0[k] + w * s1[k]
    # layer_upper_radii is a tuple of floats — interp elementwise
    out['layer_upper_radii'] = tuple(
        (1.0 - w) * float(a) + w * float(b)
        for a, b in zip(s0['layer_upper_radii'], s1['layer_upper_radii'])
    )
    # Carry-through keys (identical across grid)
    for k in ('phases', 'changeIndices', 'n_layers', 'layer_types', 'region_phases',
              'omega', 'eccentricity', 'host_mass', 'a_m', 'R_body_m', 'Mtot_kg'):
        out[k] = s0[k]
    return out


# ============================================================
# Step 2: Forward model
# ============================================================

def forward_model(theta, structure_grid, return_heating=False):
    """Compute k2 and per-phase heating for given rheology + Tb parameters.

    Args:
        theta: [alpha, log10_zeta, log10_eta_Ih, log10_eta_III,
                log10_eta_V, log10_eta_VI, log10_eta_sil, Tb_K]
        structure_grid: dict from build_or_load_structure_grid().
                        Single-Tb dicts (legacy, no 'Tb_K_grid' key) still accepted
                        — treated as frozen structure with Tb = data['Tb_K'].
    Returns:
        (Re_k2, Im_k2, perPhase_W)
    """
    alpha, log10_zeta, log10_eta_Ih, log10_eta_III, log10_eta_V, \
        log10_eta_VI, log10_eta_sil, Tb_sampled = theta
    eta_Ih = 10 ** log10_eta_Ih
    eta_III = 10 ** log10_eta_III
    eta_V = 10 ** log10_eta_V
    eta_VI = 10 ** log10_eta_VI
    eta_sil = 10 ** log10_eta_sil

    # Interpolate structure at sampled Tb (Option-2 grid interp)
    if isinstance(structure_grid, dict) and 'Tb_K_grid' in structure_grid:
        data = _interp_structure(Tb_sampled, structure_grid)
    else:
        # Legacy single-structure input — no Tb variation captured
        data = structure_grid

    # --- No-ocean safeguard ---
    # Reject sample if any Ih cell's T would exceed the Ih-L melting curve.
    # Linearized Tm_Ih(P) = 273.16 - 0.068*P_MPa (accurate to <1 K over 0-210 MPa).
    # 0.1 K safety margin under the melt line keeps us in the solid-Ih stability
    # field at every grid cell, respecting the no-ocean model assumption.
    phases = data['phases']
    P_MPa_prof = data.get('P_MPa')
    T_K_prof = data.get('T_K')
    if P_MPa_prof is not None and T_K_prof is not None:
        P_arr = np.asarray(P_MPa_prof)
        T_arr = np.asarray(T_K_prof)
        if P_arr.shape == T_arr.shape == phases.shape:
            Ih_mask = (phases == 1)
            if np.any(Ih_mask) and np.all(np.isfinite(P_arr[Ih_mask])):
                Tm_Ih_lin = 273.16 - 0.068 * P_arr[Ih_mask]
                if np.any(T_arr[Ih_mask] >= Tm_Ih_lin - 0.1):
                    # Ocean would form — reject
                    return np.nan, np.nan, {}

    # Apply per-phase viscosity overrides
    eta_mod = data['eta_Pa_base'].copy()
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
    # Uses SAMPLED Tb (not cached) so the Arrhenius reference moves consistently
    # with the interpolated T(r) profile.
    Tb_K_val = float(Tb_sampled)
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
        uniform(loc=10.0,  scale=6.0),                    # log10_eta_Ih:  [10, 16]  (Petricca low-η admissible)
        uniform(loc=10.0,  scale=6.0),                    # log10_eta_III: [10, 16]
        uniform(loc=10.0,  scale=6.0),                    # log10_eta_V:   [10, 16]
        uniform(loc=10.0,  scale=6.0),                    # log10_eta_VI:  [10, 16]
        uniform(loc=18.0,  scale=4.0),                    # log10_eta_sil: [18, 22]
        uniform(loc=float(STRUCTURE_TB_GRID[0]),
                scale=float(STRUCTURE_TB_GRID[-1] - STRUCTURE_TB_GRID[0])),
                                                          # Tb_K: [249, 250.965]  (triple-pt depression band)
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

def make_plots(samples, log_prob, structure_grid):
    """Generate diagnostic plots via the shared mcmc_plots toolkit.

    Data preparation (posterior subsample selection, per-sample forward-model
    re-evaluation, heating arrays) stays here because it depends on Test50's
    body-specific structure dict schema.  Body-agnostic plotting is delegated
    to PlanetProfile.Inference.mcmc_plots.

    N_DIM=8 note: plot_heating_vs_parameters builds its 2x5 scatter grid from
    eval_samples cols 0..4 (5 panels), extra_xvals (0 panels — we pass []),
    then cols 6+ of eval_samples (2 panels: col 6 = log10_eta_sil, col 7 = Tb_K).
    That gives 7 active panels out of 10 slots; the remaining 3 are set invisible
    by the function.  No padding is needed.
    """
    import matplotlib
    matplotlib.use('Agg')
    from PlanetProfile.Inference import mcmc_plots as mp

    # --- Re-evaluate subset for heating and k2 ---
    n_eval = min(N_REEVAL, len(samples))
    rng = np.random.default_rng(RANDOM_STATE)
    eval_idx = rng.choice(len(samples), size=n_eval, replace=False)
    eval_idx.sort()

    heating_results = []
    k2_results_list = []
    for i in eval_idx:
        Re_k2_i, Im_k2_i, heat_i = forward_model(
            samples[i], structure_grid, return_heating=True)
        heating_results.append(heat_i)
        k2_results_list.append((Re_k2_i, abs(Im_k2_i)))

    k2_results = np.array(k2_results_list)  # shape (n_eval, 2)

    totals = np.array([sum(h.values()) if h else 1e-30 for h in heating_results])
    safe_tot = np.where(totals > 1e-30, totals, 1e-30)
    f_sil = np.array([h.get('Sil', 0.0) for h in heating_results]) / safe_tot

    eval_samples = samples[eval_idx]

    # Representative mid-grid structure for headline diagnostic numbers
    if isinstance(structure_grid, dict) and 'structures' in structure_grid:
        rep = structure_grid['structures'][len(structure_grid['structures']) // 2]
    else:
        rep = structure_grid

    out = lambda name: os.path.join(OUTPUT_DIR, f'allice_yao2014_andrade_{name}.png')
    title_prefix = (f'All-Ice Titan: D_hsphere={rep["D_hsphere_km"]:.0f} km, '
                    f'CMR2={rep["CMR2"]:.4f}')

    # --- 1. Corner plot ---
    mp.plot_corner(
        samples, PARAM_LABELS,
        title=title_prefix + ' — Posterior',
        output_path=out('corner'),
    )

    # --- 2. k2 scatter coloured by silicate heating fraction ---
    mp.plot_k2_scatter_by(
        k2_results, color_values=f_sil,
        colorbar_label='Silicate heating fraction',
        obs_re=OBS['Re_k2'], obs_im=OBS['Im_k2'],
        obs_re_err=OBS['Re_k2_err'], obs_im_err=OBS['Im_k2_err'],
        title=title_prefix + r' — $k_2$ Posterior',
        output_path=out('k2_scatter'),
        cmap='RdYlBu_r', vmin=0, vmax=1,
    )

    # --- 3. Heating vs parameters + cumulative bar ---
    # extra_xvals=[] because Test50 has no D_ocean or D_iceIh derived from a
    # grid lookup — the structure is fixed.  plot_heating_vs_parameters will
    # use cols 0..4 and col 6 of eval_samples (6 of 10 scatter slots), hiding
    # the remaining 4 slots.  eval_d_ocean is all-zeros (no ocean layer) which
    # makes the secondary sort in the cumulative bar trivial; f_sil is the
    # primary sort key.
    mp.plot_heating_vs_parameters(
        eval_samples, heating_results, PARAM_LABELS,
        extra_xvals=[], extra_xlabels=[],
        output_path=out('heating'),
        cumulative_bar=True,
        eval_d_ocean=np.zeros(n_eval),
        title=(title_prefix + ' — Tidal Heating vs Parameters (top) '
               'and Per-Model Fractions (bottom)'),
    )

    # --- Summary statistics ---
    log.info('=== Posterior Summary (All-Ice, 8D Rheology+Tb) ===')
    log.info(f'  Representative mid-grid structure (Tb={rep["Tb_K"]:.3f} K): '
             f'D_hsphere={rep["D_hsphere_km"]:.1f} km, '
             f'CMR2={rep["CMR2"]:.4f}, rho_sil={rep["rhoSil_kgm3"]:.0f} kg/m3')
    log.info(f'  Layers: Ih={rep["D_iceIh_km"]:.1f}, '
             f'III={rep["D_iceIII_km"]:.1f}, '
             f'V={rep["D_iceV_km"]:.1f}, '
             f'VI={rep["D_iceVI_km"]:.1f} km')
    for i, name in enumerate(PARAM_NAMES):
        med = np.median(samples[:, i])
        lo, hi = np.percentile(samples[:, i], [16, 84])
        log.info(f'  {name}: {med:.3f} [{lo:.3f}, {hi:.3f}]')
    re_arr = k2_results[:, 0]
    im_arr = k2_results[:, 1]
    valid = np.isfinite(re_arr)
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

    # Build / load Tb grid of all-ice structures
    structure_grid = build_or_load_structure_grid(
        force_rebuild=args.rebuild_structure)

    if args.build_structure:
        log.info('--build-structure: exiting after structure build.')
        sys.exit(0)

    pkl_path = os.path.join(OUTPUT_DIR, 'allice_yao2014_andrade_mcmc_results.pkl')

    # Representative structure for diagnostic headline numbers (mid-grid Tb)
    rep_structure = structure_grid['structures'][len(structure_grid['structures']) // 2]

    if args.replot:
        log.info(f'Loading results from {pkl_path}')
        with open(pkl_path, 'rb') as f:
            results = pickle.load(f)
        samples, log_prob = results['samples'], results['log_prob']
    else:
        samples, log_prob = run_mcmc(structure_grid)
        results = {'samples': samples, 'log_prob': log_prob,
                   'param_names': PARAM_NAMES, 'obs': OBS,
                   'structure_info': {
                       'Tb_K_grid': structure_grid['Tb_K_grid'],
                       'Tb_K_mid': float(rep_structure['Tb_K']),
                       'D_hsphere_km': rep_structure['D_hsphere_km'],
                       'CMR2': rep_structure['CMR2'],
                       'rhoSil_kgm3': rep_structure['rhoSil_kgm3'],
                   }}
        _save_cache(results, pkl_path)
        log.info(f'MCMC results saved to {pkl_path}')

    if not args.no_plots:
        make_plots(samples, log_prob, structure_grid)

    log.info('Done.')
