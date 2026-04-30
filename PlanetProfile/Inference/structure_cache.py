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
from typing import Dict, Optional, Any, List, Tuple
from pathlib import Path

log = logging.getLogger('PlanetProfile')

# Phase constant (from PlanetProfile.Utilities.defineStructs.Constants)
PHASE_CLATH = 30  # Clathrate phase ID


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
    filepath = Path(filepath)
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

    # Configure for structure calculation
    configParams.Gravity.backend = 'tidalpy'
    if rheology == 'maxwell':
        configParams.Gravity.rheology_models = 'maxwell'
    configParams.CALC_NEW = True
    configParams.CALC_NEW_GRAVITY = True
    configParams.NO_SAVEFILE = True
    configParams.SKIP_PLOTS = True

    # Run PlanetProfile
    Planet, Params = PlanetProfile(Planet, configParams)
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

    return data


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
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    with open(filepath, 'wb') as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

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
                configParams.Gravity.rheology_models = 'maxwell'
            configParams.CALC_NEW = True
            configParams.CALC_NEW_GRAVITY = True
            configParams.NO_SAVEFILE = True
            configParams.SKIP_PLOTS = True

            Planet, Params = PlanetProfile(Planet, configParams)
            Params.CALC_NEW_GRAVITY = True
            Planet, Params = SetupGravity(Planet, Params)

            # Extract and cache
            structure_data = extract_structure_from_planet(Planet, Params)
            structure_data['param_name'] = param_name
            structure_data['param_value'] = param_val
            cache[param_val] = structure_data

            log.info(f"    {structure_data['n_layers']} layers, {len(structure_data['r_m'])} points")

        except Exception as e:
            log.warning(f"    FAILED: {e}")
            continue

    elapsed = time.time() - t0
    n_in_grid = sum(1 for val in param_values if val in cache)
    log.info(f"Structure grid complete: {n_in_grid}/{len(param_values)} succeeded "
             f"(+{len(missing)} attempted in {elapsed/60:.1f} min)")

    # Save updated cache
    save_structure_cache(cache, cache_path)

    return cache


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
