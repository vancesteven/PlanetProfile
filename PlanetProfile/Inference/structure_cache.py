"""
Structure caching utilities for inference.

Handles loading and validation of pre-computed Planet structure data.
Supports both single-structure and multi-structure (grid) caching for
parameter spaces that include structural variation (e.g., Tb_K, ocean salinity).

Author: PlanetProfile Team
Date: 2026-04-29
"""
import os
import sys
import pickle
import time
import logging
import importlib
import numpy as np
from typing import Dict, Optional, Any, List, Tuple, Sequence
from pathlib import Path

log = logging.getLogger('PlanetProfile')

# Phase constant (from PlanetProfile.Utilities.defineStructs.Constants)
PHASE_CLATH = 30  # Clathrate phase ID

# TidalPy requires at least this many radial points per layer
_TIDALPY_MIN_PTS = 5


def _expand_thin_reduced_layers(Planet):
    """Expand any layers in Planet.Reduced that have fewer than _TIDALPY_MIN_PTS points.

    Must be called after PlanetProfile() and before SetupGravity().  TidalPy
    requires ≥5 radial slices per layer; layers with fewer points cause
    "N slices when at least 5 are required" errors.  This can happen when a
    porous-silicate phase transition creates a 1-point boundary segment in the
    reduced model.

    Operates on the surface-down (decreasing-r) convention used by Planet.Reduced.
    """
    if not hasattr(Planet, 'Reduced') or Planet.Reduced.changeIndices is None:
        return Planet

    ci = list(Planet.Reduced.changeIndices)
    n_layers = len(ci) - 1
    if n_layers == 0:
        return Planet

    min_pts = _TIDALPY_MIN_PTS

    # Quick check: any thin layers?
    if all(ci[k + 1] - ci[k] >= min_pts for k in range(n_layers)):
        return Planet  # nothing to do

    r_in = np.array(Planet.Reduced.r_m, dtype=float)
    rho_in = np.array(Planet.Reduced.rho_kgm3, dtype=float)
    phase_in = list(Planet.Reduced.phase)
    iConv_in = list(Planet.Reduced.iConv)

    has_seismic = (hasattr(Planet.Reduced, 'Seismic')
                   and hasattr(Planet.Reduced.Seismic, 'VP_kms')
                   and Planet.Reduced.Seismic.VP_kms is not None)
    has_visc = (hasattr(Planet.Reduced, 'eta_Pas')
                and Planet.Reduced.eta_Pas is not None)

    if has_seismic:
        vp_in = np.array(Planet.Reduced.Seismic.VP_kms, dtype=float)
        vs_in = np.array(Planet.Reduced.Seismic.VS_kms, dtype=float)
        gs_in = np.array(Planet.Reduced.Seismic.GS_GPa, dtype=float)
    if has_visc:
        eta_in = np.array(Planet.Reduced.eta_Pas, dtype=float)

    new_r, new_rho, new_phase, new_iConv = [], [], [], []
    new_ci = [0]
    if has_seismic:
        new_vp, new_vs, new_gs = [], [], []
    if has_visc:
        new_eta = []

    for k in range(n_layers):
        s, e = ci[k], ci[k + 1]
        n_pts = e - s

        if n_pts < min_pts:
            # Surface-down: r decreases within [s, e).
            r_outer = r_in[s]      # largest r of this layer (index s)
            # Determine the inner radial boundary to expand toward.
            if k < n_layers - 1:
                # Non-innermost thin layer: expand down to just above the next
                # layer's outer edge, leaving a small margin to avoid overlap.
                r_lower = r_in[ci[k + 1]]
                # Place min_pts points from r_outer to (min_pts steps above r_lower)
                r_exp = np.linspace(r_outer, r_lower, min_pts + 1)[:-1]
            else:
                # Innermost layer: expand toward the planetary center (r → 0).
                # Use r_outer/min_pts as the inner endpoint.
                r_inner = max(r_outer / min_pts, 1.0)
                r_exp = np.linspace(r_outer, r_inner, min_pts)

            # Properties are constant (replicate the single point's values).
            new_r.extend(r_exp.tolist())
            new_rho.extend([rho_in[s]] * min_pts)
            new_phase.extend([phase_in[s]] * min_pts)
            new_iConv.extend([iConv_in[s]] * min_pts)
            if has_seismic:
                new_vp.extend([vp_in[s]] * min_pts)
                new_vs.extend([vs_in[s]] * min_pts)
                new_gs.extend([gs_in[s]] * min_pts)
            if has_visc:
                new_eta.extend([eta_in[s]] * min_pts)
            new_ci.append(new_ci[-1] + min_pts)
            log.info(f'TidalPy prep: expanded layer {k} from {n_pts} to {min_pts} points '
                     f'(r = {r_exp[0]/1e3:.1f}-{r_exp[-1]/1e3:.1f} km)')
        else:
            new_r.extend(r_in[s:e].tolist())
            new_rho.extend(rho_in[s:e].tolist())
            new_phase.extend(phase_in[s:e])
            new_iConv.extend(iConv_in[s:e])
            if has_seismic:
                new_vp.extend(vp_in[s:e].tolist())
                new_vs.extend(vs_in[s:e].tolist())
                new_gs.extend(gs_in[s:e].tolist())
            if has_visc:
                new_eta.extend(eta_in[s:e].tolist())
            new_ci.append(new_ci[-1] + n_pts)

    Planet.Reduced.r_m = new_r
    Planet.Reduced.rho_kgm3 = new_rho
    Planet.Reduced.phase = new_phase
    Planet.Reduced.iConv = new_iConv
    Planet.Reduced.changeIndices = new_ci
    if has_seismic:
        Planet.Reduced.Seismic.VP_kms = new_vp
        Planet.Reduced.Seismic.VS_kms = new_vs
        Planet.Reduced.Seismic.GS_GPa = new_gs
    if has_visc:
        Planet.Reduced.eta_Pas = new_eta

    return Planet


def resolve_cache_path(filepath) -> Path:
    """Resolve a structure-cache path that may be repo-relative.

    Config JSONs store repo-relative paths (e.g.
    'PlanetProfile/Test/mcmc_results/...') so result pickles stay portable
    across machines. Those resolve fine when the CWD is the repo root (CLI
    convention) but not when a caller runs from elsewhere (e.g. the Streamlit
    app runs from PlanetProfileApp/). If the path as given doesn't exist and
    isn't absolute, fall back to resolving against the repo root inferred
    from this package's location.
    """
    p = Path(filepath)
    if p.exists() or p.is_absolute():
        return p
    repo_root = Path(__file__).resolve().parents[2]
    candidate = repo_root / p
    if candidate.exists():
        return candidate
    return p


def load_structure_cache(
    filepath: str,
    validate_bodyname: Optional[str] = None
) -> Dict[str, Any]:
    """
    Load cached Planet structure from file.

    Args:
        filepath: Path to cache file (.pkl)
        validate_bodyname: If provided, check that cached structure matches this body

    Returns:
        Dict with structural arrays and metadata

    Raises:
        FileNotFoundError: If cache file doesn't exist
        ValueError: If validation fails
    """
    filepath = resolve_cache_path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Structure cache not found: {filepath}")

    log.info(f"Loading structure cache from {filepath}")

    with open(filepath, 'rb') as f:
        data = pickle.load(f)

    # Validate structure
    required_keys = ['r_m', 'rho', 'K_Pa', 'mu_Pa', 'eta_Pa_base', 'phases',
                     'changeIndices', 'n_layers', 'omega', 'eccentricity']
    missing = [k for k in required_keys if k not in data]
    if missing:
        raise ValueError(f"Cache missing required keys: {missing}")

    # Validate bodyname if requested
    if validate_bodyname is not None:
        cached_bodyname = data.get('bodyname', 'Unknown')
        if cached_bodyname.lower() != validate_bodyname.lower():
            log.warning(f"Cache bodyname '{cached_bodyname}' doesn't match "
                       f"requested '{validate_bodyname}'")

    log.info(f"  Loaded: {data['n_layers']} layers, {len(data['r_m'])} radial points")

    return data


def build_structure_from_pptest(
    test_module_name: str,
    rheology: str = 'andrade',
    force_rebuild: bool = False
) -> Dict[str, Any]:
    """
    Build structural model from a PPTest configuration.

    This runs PlanetProfile once to compute the self-consistent interior structure,
    then extracts and caches the arrays needed for fast forward model evaluation.

    Args:
        test_module_name: Full module path (e.g., 'PlanetProfile.Test.PPTest41')
        rheology: 'andrade' or 'maxwell' (affects TidalPy setup)
        force_rebuild: If True, force recomputation even if cached

    Returns:
        Dict with structural arrays and metadata
    """
    from PlanetProfile.Main import PlanetProfile
    from PlanetProfile.Gravity.Gravity import SetupGravity
    from PlanetProfile.Utilities.defineStructs import EOSlist, Constants

    log.info(f"Building structure from {test_module_name}...")
    t0 = time.time()

    # Clear EOS cache and import test module
    EOSlist.loaded.clear()
    from PlanetProfile.GetConfig import Params as configParams

    if test_module_name in sys.modules:
        importlib.reload(sys.modules[test_module_name])
    else:
        importlib.import_module(test_module_name)

    mod = sys.modules[test_module_name]
    Planet = mod.Planet

    # Inference framework always treats viscosity as a free MCMC parameter;
    # Arrhenius temperature-dependence (if any) is applied separately via
    # arrhenius_params.  Force it off here so eta_Pa_base in the cache
    # reflects plain reference viscosities, not temperature-corrected ones.
    Planet.Do.ARRHENIUS_VISCOSITY = False

    # Configure for structure calculation
    configParams.Gravity.backend = 'tidalpy'
    if rheology == 'maxwell':
        configParams.Gravity.rheology_models = {
            '0': 'newton', 'Ih': 'maxwell', 'Ih_conv': 'maxwell',
            'II': 'maxwell', 'III': 'maxwell', 'III_conv': 'maxwell',
            'IV': 'maxwell', 'V': 'maxwell', 'V_conv': 'maxwell', 'VI': 'maxwell',
            'Sil': 'maxwell', 'Fe': 'elastic', 'Clath': 'elastic', 'Clath_conv': 'maxwell'
        }
    configParams.CALC_NEW = True
    configParams.CALC_NEW_GRAVITY = True
    configParams.NO_SAVEFILE = True
    configParams.SKIP_PLOTS = True
    configParams.CALC_CONDUCT = True  # Required for induction layer extraction

    # Slightly loosen MoI validation so the cache includes structures at all
    # Tb_K even if CMR² drifts a little from the target.  Do NOT use ±1.0:
    # CalcMoIWithEOS zero-initialises C_kgm2 and a wide window captures those
    # zero entries in CMR2inds, causing an IndexError at the iValid lookup.
    Planet.Bulk.CuncertaintyLower = 0.05
    Planet.Bulk.CuncertaintyUpper = 0.05

    # Run PlanetProfile
    Planet, Params = PlanetProfile(Planet, configParams)
    Planet = _expand_thin_reduced_layers(Planet)
    Params.CALC_NEW_GRAVITY = True
    Planet, Params = SetupGravity(Planet, Params)

    # Extract structure
    data = extract_structure_from_planet(Planet, Params)
    data['source_module'] = test_module_name
    data['rheology'] = rheology

    log.info(f"Structure built in {time.time()-t0:.1f}s: {data['n_layers']} layers, "
             f"{len(data['r_m'])} points")

    return data


def extract_structure_from_planet(
    Planet,
    Params,
    min_layer_points: int = 5
) -> Dict[str, Any]:
    """
    Extract structural arrays from Planet object for caching.

    This function pulls out the radial profiles, layer boundaries, and
    orbital parameters needed for fast forward model evaluation.

    Args:
        Planet: PlanetProfile Planet object (after SetupGravity)
        Params: PlanetProfile Params object
        min_layer_points: Minimum points per layer (pad thin layers)

    Returns:
        Dict with all arrays needed for forward_model_k2()
    """
    from PlanetProfile.Thermodynamics.Electrical import PhaseConv
    from PlanetProfile.Utilities.defineStructs import Constants

    # Extract arrays from ALMA model
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

    # Compute bulk modulus from VP
    K_Pa = rho * VP_ms**2 - (4.0 / 3.0) * mu_Pa

    # Handle NaN/inf in K_Pa (estimate from Poisson's ratio)
    nan_mask = ~np.isfinite(K_Pa) | (K_Pa <= 0)
    if np.any(nan_mask):
        for i in np.where(nan_mask)[0]:
            ph = int(phases[i])
            if ph >= 50 and ph < 100:  # Silicate
                nu = 0.25
            elif ph >= 100:  # Fe core
                nu = 0.29
            else:  # Ice/liquid
                nu = 0.33
            K_Pa[i] = 2.0 * mu_Pa[i] * (1.0 + nu) / (3.0 * (1.0 - 2.0 * nu))

    K_Pa = np.maximum(K_Pa, 1e6)  # Minimum bulk modulus

    # Layer boundaries
    changeIndices = np.max(Planet.Reduced.changeIndices) - np.flipud(Planet.Reduced.changeIndices)
    n_layers = len(changeIndices) - 1

    # Region phase labels (for rheology assignment)
    _orig_iConv = np.flipud(Planet.Reduced.iConv)
    _orig_changeIndices = changeIndices.copy()
    region_phases = []

    for i_layer in range(n_layers):
        start = _orig_changeIndices[i_layer]
        phase = phases[start]

        # Clathrates: map to single phase
        if phase >= PHASE_CLATH and phase < PHASE_CLATH + 10:
            phase = PHASE_CLATH

        convection = _orig_iConv[start]
        phase_str = PhaseConv(phase, liq='0')
        if convection:
            phase_str += '_conv'

        region_phases.append(phase_str)

    # Bulk viscosity (zeros for now, could be extended)
    bulk_visc = np.zeros_like(eta_Pa_base)

    # Pad thin layers for TidalPy stability
    needs_padding = any(
        changeIndices[i+1] - changeIndices[i] < min_layer_points
        for i in range(n_layers)
    )

    if needs_padding:
        new_r, new_rho, new_K, new_mu, new_eta, new_phases, new_bv = \
            [], [], [], [], [], [], []
        new_ci = [0]

        for i_layer in range(n_layers):
            s, e = changeIndices[i_layer], changeIndices[i_layer + 1]
            n_pts = e - s

            if n_pts < min_layer_points and n_pts >= 2:
                # Interpolate to min_layer_points
                r_layer = r_m[s:e]
                r_interp = np.linspace(r_layer[0], r_layer[-1], min_layer_points)

                new_r.append(r_interp)
                new_rho.append(np.interp(r_interp, r_layer, rho[s:e]))
                new_K.append(np.interp(r_interp, r_layer, K_Pa[s:e]))
                new_mu.append(np.interp(r_interp, r_layer, mu_Pa[s:e]))
                new_eta.append(np.interp(r_interp, r_layer, eta_Pa_base[s:e]))
                new_phases.append(np.full(min_layer_points, phases[s]))
                new_bv.append(np.zeros(min_layer_points))
                new_ci.append(new_ci[-1] + min_layer_points)
            else:
                # Keep as-is
                new_r.append(r_m[s:e])
                new_rho.append(rho[s:e])
                new_K.append(K_Pa[s:e])
                new_mu.append(mu_Pa[s:e])
                new_eta.append(eta_Pa_base[s:e])
                new_phases.append(phases[s:e])
                new_bv.append(bulk_visc[s:e])
                new_ci.append(new_ci[-1] + (e - s))

        # Concatenate padded arrays
        r_m = np.concatenate(new_r)
        rho = np.concatenate(new_rho)
        K_Pa = np.concatenate(new_K)
        mu_Pa = np.concatenate(new_mu)
        eta_Pa_base = np.concatenate(new_eta)
        phases = np.concatenate(new_phases)
        bulk_visc = np.concatenate(new_bv)
        changeIndices = np.array(new_ci)

    # Layer metadata
    layer_upper_radii = []
    layer_types = []
    for i_layer in range(n_layers):
        end = changeIndices[i_layer + 1]
        layer_upper_radii.append(r_m[end - 1])
        layer_types.append('liquid' if phases[changeIndices[i_layer]] == 0 else 'solid')

    # Orbital parameters
    omega = Planet.Bulk.meanMotion_radps
    ecc = Planet.Bulk.eccentricity
    host_mass = Constants.parentMass_kg.get(Planet.parent, 1.898e27)  # Default to Jupiter
    a_m = (Constants.G * host_mass / omega**2) ** (1.0 / 3.0)

    # Phase map for heating aggregation
    phase_map = {0: '0', 1: 'Ih', 2: 'II', 3: 'III', 5: 'V', 6: 'VI'}

    # Temperature profile (if available, for Arrhenius)
    T_K = None
    if hasattr(Planet, 'T_K') and Planet.T_K is not None:
        # Interpolate onto structure grid
        # (This requires matching Planet's original radial grid to extracted r_m)
        # For now, set to None if not directly available
        pass

    # Build output dict
    data = {
        'r_m': np.ascontiguousarray(r_m),
        'rho': np.ascontiguousarray(rho),
        'K_Pa': np.ascontiguousarray(K_Pa),
        'mu_Pa': np.ascontiguousarray(mu_Pa),
        'eta_Pa_base': eta_Pa_base,
        'phases': phases,
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
        'phase_map': phase_map,
        'T_K': T_K,
        'bodyname': Planet.name,
    }

    _cmr2 = getattr(Planet, 'CMR2mean', None)
    _cmr2 = float(_cmr2) if (_cmr2 is not None and np.isfinite(_cmr2)) else np.nan
    data['CMR2mean'] = _cmr2
    data['CMR2'] = _cmr2  # Backward-compatible alias used by current inference likelihoods.
    data['Mtot_kg'] = float(getattr(Planet, 'Mtot_kg', getattr(Planet.Bulk, 'M_kg', np.nan)))
    data['rhoSilWithCore_kgm3'] = float(getattr(Planet.Sil, 'rhoSilWithCore_kgm3', np.nan))
    data['rhoSil_kgm3'] = data['rhoSilWithCore_kgm3']

    # Ocean and Ice Ih thicknesses derived from layer profile
    D_ocean_m = 0.0
    D_iceIh_m = 0.0
    for i_layer in range(n_layers):
        s, e = changeIndices[i_layer], changeIndices[i_layer + 1]
        thick = r_m[e - 1] - r_m[s]
        ph = region_phases[i_layer]
        if ph == '0':
            D_ocean_m += thick
        elif ph.startswith('Ih'):
            D_iceIh_m += thick
    data['D_ocean_km'] = D_ocean_m / 1e3
    data['D_iceIh_km'] = D_iceIh_m / 1e3

    # -------------------------------------------------------------------------
    # Induction layer data (composition-agnostic)
    # PP's conductivity pipeline (SwConduct / RktConduct / by-ion) has already
    # populated Planet.sigma_Sm.  SetupInduction converts that to the compact
    # layer representation used by MoonMag's AeList.
    # -------------------------------------------------------------------------
    R_body_m = float(Planet.Bulk.R_m)
    _rSigChange_m = None
    _sigmaLayers_Sm = None
    _sigma_ocean_mean_Sm = np.nan

    try:
        _mag = getattr(Planet, 'Magnetic', None)
        _rsc = getattr(_mag, 'rSigChange_m', None)
        _sig = getattr(_mag, 'sigmaLayers_Sm', None)

        if _rsc is None or not hasattr(_rsc, '__len__') or len(_rsc) == 0:
            from PlanetProfile.MagneticInduction.MagneticInduction import SetupInduction
            _old_cc = getattr(Params, 'CALC_CONDUCT', False)
            Params.CALC_CONDUCT = True
            Planet, Params = SetupInduction(Planet, Params)
            Params.CALC_CONDUCT = _old_cc
            _rsc = getattr(getattr(Planet, 'Magnetic', None), 'rSigChange_m', None)
            _sig = getattr(getattr(Planet, 'Magnetic', None), 'sigmaLayers_Sm', None)

        if _rsc is not None and len(_rsc) > 0:
            _rSigChange_m = np.ascontiguousarray(_rsc, dtype=np.float64)
            _sigmaLayers_Sm = np.ascontiguousarray(_sig, dtype=np.float64)

            # Mean ocean conductivity (diagnostic): ocean occupies the band
            # from the ice-shelf base down to the seafloor.
            if D_ocean_m > 0:
                r_ocean_top = R_body_m - D_iceIh_m
                r_ocean_bot = r_ocean_top - D_ocean_m
                ocean_mask = (
                    (_rSigChange_m >= r_ocean_bot - 1e3) &
                    (_rSigChange_m <= r_ocean_top + 1e3)
                )
                if np.any(ocean_mask):
                    _sigma_ocean_mean_Sm = float(np.mean(_sigmaLayers_Sm[ocean_mask]))

    except Exception as e:
        log.warning(f"Could not extract induction layers from Planet: {e}")

    data['rSigChange_m'] = _rSigChange_m
    data['sigmaLayers_Sm'] = _sigmaLayers_Sm
    data['wOcean_ppt'] = float(getattr(Planet.Ocean, 'wOcean_ppt', np.nan))
    data['sigma_ocean_mean_Sm'] = _sigma_ocean_mean_Sm
    data['R_body_m'] = R_body_m

    return data


def _apply_structure_grid_parameter(Planet, param_name: str, param_value: float) -> None:
    """Apply a structure-grid control parameter to a Planet object before PP run."""
    if param_name == 'Tb_K':
        Planet.Bulk.Tb_K = float(param_value)
    elif param_name == 'rhoSilInput_kgm3':
        # Scott Chang Monte Carlo semantics: direct silicate-density control.
        # This intentionally trades EOS self-consistency for a controlled mass/MoI axis.
        Planet.Sil.rhoSilWithCore_kgm3 = float(param_value)
        Planet.Do.CONSTANT_INNER_DENSITY = True
    else:
        raise ValueError(f"Unsupported structure-grid parameter: {param_name}")


def _grid_key(values: Sequence[float]) -> Tuple[float, ...]:
    """Normalize floating grid coordinates for stable pickle keys."""
    return tuple(float(v) for v in values)


def save_structure_cache(
    data: Dict[str, Any],
    filepath: str
) -> None:
    """
    Save structure cache to file.

    Args:
        data: Structure dict from build_structure_from_pptest() or extract_structure_from_planet()
        filepath: Output path (.pkl extension)
    """
    import tempfile
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    # Atomic write: dump to temp file, then rename — prevents corruption on kill
    fd, tmp_path = tempfile.mkstemp(dir=filepath.parent, suffix='.pkl.tmp')
    try:
        with os.fdopen(fd, 'wb') as f:
            pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(tmp_path, filepath)
    except BaseException:
        os.unlink(tmp_path)
        raise

    file_size_mb = filepath.stat().st_size / (1024 * 1024)
    log.info(f"Structure cache saved to {filepath} ({file_size_mb:.1f} MB)")


def build_structure_grid(
    test_module_name: str,
    param_name: str,
    param_values: List[float],
    cache_path: str,
    rheology: str = 'andrade',
    force_rebuild: bool = False
) -> Dict[float, Dict[str, Any]]:
    """
    Build grid of structures varying a single parameter (e.g., Tb_K, ocean salinity).

    This is used for parameter spaces where structure itself varies, requiring
    nearest-neighbor lookup during MCMC/SBI.

    Args:
        test_module_name: Full module path (e.g., 'PlanetProfile.Test.PPTest42')
        param_name: Parameter to vary (e.g., 'Tb_K', 'oceanComp')
        param_values: Values to compute (e.g., [250, 255, 260, ..., 275])
        cache_path: Path to cache file (.pkl)
        rheology: 'andrade' or 'maxwell'
        force_rebuild: If True, recompute all grid points

    Returns:
        Dict mapping param_value -> structure_data
    """
    cache_path = Path(cache_path)

    # Load existing cache
    if cache_path.exists() and not force_rebuild:
        log.info(f"Loading structure grid from {cache_path}")
        with open(cache_path, 'rb') as f:
            cache = pickle.load(f)
    else:
        cache = {}

    # Determine which grid points need computation
    missing = [val for val in param_values if val not in cache]
    if not missing:
        log.info(f"All {len(param_values)} grid points already cached")
        return cache

    log.info(f"Computing {len(missing)}/{len(param_values)} new grid points for {param_name}...")

    from PlanetProfile.Main import PlanetProfile
    from PlanetProfile.Gravity.Gravity import SetupGravity
    from PlanetProfile.Utilities.defineStructs import EOSlist, Constants
    from PlanetProfile.GetConfig import Params as configParams

    t0 = time.time()

    for i_val, param_val in enumerate(missing):
        log.info(f"  [{i_val+1}/{len(missing)}] {param_name} = {param_val}")

        try:
            # Clear and reload module
            EOSlist.loaded.clear()

            if test_module_name in sys.modules:
                importlib.reload(sys.modules[test_module_name])
            else:
                importlib.import_module(test_module_name)

            mod = sys.modules[test_module_name]
            Planet = mod.Planet

            # Inference framework treats viscosity as a free MCMC parameter;
            # disable Arrhenius so eta_Pa_base in the cache is not pre-corrected.
            Planet.Do.ARRHENIUS_VISCOSITY = False

            # Set parameter value
            if param_name == 'Tb_K':
                Planet.Bulk.Tb_K = param_val
            else:
                # Add other parameter types as needed
                log.warning(f"Parameter '{param_name}' not yet supported for grid variation")
                continue

            # Configure and run
            configParams.Gravity.backend = 'tidalpy'
            if rheology == 'maxwell':
                configParams.Gravity.rheology_models = {
                    '0': 'newton', 'Ih': 'maxwell', 'Ih_conv': 'maxwell',
                    'II': 'maxwell', 'III': 'maxwell', 'III_conv': 'maxwell',
                    'IV': 'maxwell', 'V': 'maxwell', 'V_conv': 'maxwell', 'VI': 'maxwell',
                    'Sil': 'maxwell', 'Fe': 'elastic', 'Clath': 'elastic', 'Clath_conv': 'maxwell'
                }
            configParams.CALC_NEW = True
            configParams.CALC_NEW_GRAVITY = True
            configParams.NO_SAVEFILE = True
            configParams.SKIP_PLOTS = True
            configParams.CALC_CONDUCT = True  # Required for induction layer extraction

            # ±0.05 keeps the window away from zero (zero-init'd C_kgm2 entries
            # inside CalcMoIWithEOS would otherwise land in CMR2inds and break the
            # iValid lookup at the Fe_CORE=False branch with a wider window).
            Planet.Bulk.CuncertaintyLower = 0.05
            Planet.Bulk.CuncertaintyUpper = 0.05

            try:
                Planet, Params = PlanetProfile(Planet, configParams)
            except Exception as e_pp:
                # PlanetProfile's internal SetupGravity can fail with a TidalPy
                # thin-layer error when a porous-silicate phase transition occupies
                # only 1 radial point in Reduced.  Planet is mutated in-place so
                # Planet.Reduced is fully populated at the point of failure.
                _is_thin_layer = ('slices when at least' in str(e_pp)
                                  or 'has 2 slices' in str(e_pp))
                if not _is_thin_layer:
                    raise
                if not (hasattr(Planet, 'Reduced')
                        and Planet.Reduced.changeIndices is not None):
                    raise RuntimeError(
                        f"PlanetProfile failed before Reduced model was built: {e_pp}"
                    ) from e_pp
                Params = configParams
                log.info(f"    Internal thin-layer error; expanding layers and retrying gravity")
            Planet = _expand_thin_reduced_layers(Planet)
            Params.CALC_NEW_GRAVITY = True
            Planet, Params = SetupGravity(Planet, Params)

            # Extract and cache
            structure_data = extract_structure_from_planet(Planet, Params)
            structure_data['param_name'] = param_name
            structure_data['param_value'] = param_val
            cache[param_val] = structure_data

            log.info(f"    {structure_data['n_layers']} layers, {len(structure_data['r_m'])} points")
            # Save incrementally so a killed process doesn't lose all progress
            save_structure_cache(cache, cache_path)

        except Exception as e:
            log.warning(f"    FAILED: {e}")
            continue

    elapsed = time.time() - t0
    n_in_grid = sum(1 for val in param_values if val in cache)
    log.info(f"Structure grid complete: {n_in_grid}/{len(param_values)} succeeded "
             f"(+{len(missing)} attempted in {elapsed/60:.1f} min)")

    # Final save (already saved incrementally, but update the summary log entry)
    save_structure_cache(cache, cache_path)

    return cache


def build_structure_nd_grid(
    test_module_name: str,
    grid_params: Dict[str, Sequence[float]],
    cache_path: str,
    rheology: str = 'andrade',
    force_rebuild: bool = False
) -> Dict[str, Any]:
    """
    Build an N-dimensional structure cache.

    Cache schema:
      {
        'grid_metadata': {
            'schema_version': 2,
            'grid_parameters': ['Tb_K', 'rhoSilInput_kgm3'],
            'grid_values': {'Tb_K': [...], 'rhoSilInput_kgm3': [...]},
            'assumptions': [...]
        },
        'grid_cache': {
            (Tb_K, rhoSilInput_kgm3): structure_data,
            ...
        }
      }
    """
    from itertools import product
    from PlanetProfile.Main import PlanetProfile
    from PlanetProfile.Gravity.Gravity import SetupGravity
    from PlanetProfile.Utilities.defineStructs import EOSlist
    from PlanetProfile.GetConfig import Params as configParams

    cache_path = Path(cache_path)
    param_names = list(grid_params.keys())
    param_values = {
        name: [float(v) for v in values]
        for name, values in grid_params.items()
    }

    if cache_path.exists() and not force_rebuild:
        log.info(f"Loading ND structure grid from {cache_path}")
        with open(cache_path, 'rb') as f:
            cache_data = pickle.load(f)
    else:
        cache_data = {
            'grid_metadata': {
                'schema_version': 2,
                'grid_parameters': param_names,
                'grid_values': param_values,
                'rheology': rheology,
                'source_module': test_module_name,
                'assumptions': [
                    'rhoSilInput_kgm3 maps to Planet.Sil.rhoSilWithCore_kgm3.',
                    'rhoSilInput_kgm3 forces Planet.Do.CONSTANT_INNER_DENSITY=True.',
                    'This gives direct silicate-density control but is not fully EOS self-consistent.',
                ],
            },
            'grid_cache': {},
        }

    grid_cache = cache_data.setdefault('grid_cache', {})
    all_points = list(product(*(param_values[name] for name in param_names)))
    missing = [point for point in all_points if _grid_key(point) not in grid_cache]

    if not missing:
        log.info(f"All {len(all_points)} ND grid points already cached")
        return cache_data

    log.info(f"Computing {len(missing)}/{len(all_points)} ND grid points for {param_names}")
    t0 = time.time()

    for i_point, point in enumerate(missing):
        theta = dict(zip(param_names, point))
        log.info(f"  [{i_point + 1}/{len(missing)}] {theta}")

        try:
            EOSlist.loaded.clear()

            if test_module_name in sys.modules:
                importlib.reload(sys.modules[test_module_name])
            else:
                importlib.import_module(test_module_name)

            Planet = sys.modules[test_module_name].Planet
            Planet.Do.ARRHENIUS_VISCOSITY = False
            for name, value in theta.items():
                _apply_structure_grid_parameter(Planet, name, value)

            configParams.Gravity.backend = 'tidalpy'
            if rheology == 'maxwell':
                configParams.Gravity.rheology_models = {
                    '0': 'newton', 'Ih': 'maxwell', 'Ih_conv': 'maxwell',
                    'II': 'maxwell', 'III': 'maxwell', 'III_conv': 'maxwell',
                    'IV': 'maxwell', 'V': 'maxwell', 'V_conv': 'maxwell', 'VI': 'maxwell',
                    'Sil': 'maxwell', 'Fe': 'elastic', 'Clath': 'elastic', 'Clath_conv': 'maxwell'
                }
            configParams.CALC_NEW = True
            configParams.CALC_NEW_GRAVITY = True
            configParams.NO_SAVEFILE = True
            configParams.SKIP_PLOTS = True
            configParams.CALC_CONDUCT = True
            Planet.Bulk.CuncertaintyLower = 0.05
            Planet.Bulk.CuncertaintyUpper = 0.05

            try:
                Planet, Params = PlanetProfile(Planet, configParams)
            except Exception as e_pp:
                _is_thin_layer = ('slices when at least' in str(e_pp)
                                  or 'has 2 slices' in str(e_pp))
                if not _is_thin_layer:
                    raise
                if not (hasattr(Planet, 'Reduced')
                        and Planet.Reduced.changeIndices is not None):
                    raise RuntimeError(
                        f"PlanetProfile failed before Reduced model was built: {e_pp}"
                    ) from e_pp
                Params = configParams
                log.info("    Internal thin-layer error; expanding layers and retrying gravity")
            Planet = _expand_thin_reduced_layers(Planet)
            Params.CALC_NEW_GRAVITY = True
            Planet, Params = SetupGravity(Planet, Params)

            structure_data = extract_structure_from_planet(Planet, Params)
            structure_data['grid_coordinates'] = theta
            structure_data['param_name'] = ','.join(param_names)
            structure_data['param_value'] = _grid_key(point)
            grid_cache[_grid_key(point)] = structure_data

            log.info(f"    {structure_data['n_layers']} layers, {len(structure_data['r_m'])} points")
            save_structure_cache(cache_data, cache_path)

        except Exception as e:
            log.warning(f"    FAILED: {e}")
            continue

    elapsed = time.time() - t0
    log.info(f"ND structure grid complete: {len(grid_cache)}/{len(all_points)} succeeded "
             f"in {elapsed/60:.1f} min")
    save_structure_cache(cache_data, cache_path)

    return cache_data


def validate_structure_cache(
    data: Dict[str, Any],
    expected_bodyname: Optional[str] = None,
    expected_n_layers: Optional[int] = None
) -> Tuple[bool, List[str]]:
    """
    Validate structure cache for common issues.

    Args:
        data: Structure dict from load_structure_cache()
        expected_bodyname: Expected planet name (optional)
        expected_n_layers: Expected number of layers (optional)

    Returns:
        (is_valid, warnings) tuple
    """
    warnings = []

    # Check bodyname
    if expected_bodyname is not None:
        cached_name = data.get('bodyname', 'Unknown')
        if cached_name.lower() != expected_bodyname.lower():
            warnings.append(f"Bodyname mismatch: cache has '{cached_name}', "
                          f"expected '{expected_bodyname}'")

    # Check number of layers
    if expected_n_layers is not None:
        actual_n_layers = data.get('n_layers', 0)
        if actual_n_layers != expected_n_layers:
            warnings.append(f"Layer count mismatch: cache has {actual_n_layers}, "
                          f"expected {expected_n_layers}")

    # Check for NaN/inf in critical arrays
    for key in ['r_m', 'rho', 'K_Pa', 'mu_Pa', 'eta_Pa_base']:
        arr = data.get(key)
        if arr is not None and not np.all(np.isfinite(arr)):
            warnings.append(f"Array '{key}' contains NaN or inf values")

    # Check orbital parameters
    if data.get('omega', 0) <= 0:
        warnings.append("Invalid orbital frequency (omega <= 0)")

    if data.get('eccentricity', 0) < 0 or data.get('eccentricity', 0) >= 1:
        warnings.append(f"Invalid eccentricity: {data.get('eccentricity')}")

    is_valid = len(warnings) == 0
    return is_valid, warnings
