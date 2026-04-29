"""
Forward models for Bayesian inference.

Provides fast k2 and heating calculations given parameter samples and cached
Planet structure. Supports both Andrade and Maxwell rheologies with optional
Arrhenius temperature-dependent viscosity.

Designed for HPC-scale MCMC/SBI: ~0.1s per evaluation.

Author: PlanetProfile Team
Date: 2026-04-29
"""
import numpy as np
from typing import Dict, Tuple, Optional, Any, List
import logging

log = logging.getLogger('PlanetProfile')


# ============================================================================
# Parameter Hook System for Extensibility
# ============================================================================

_PARAMETER_HOOKS = {}


def parameter_hook(*param_ids: str):
    """
    Decorator to register parameter application function.

    Multiple param_ids can be specified for grouped parameters that should be
    applied together (e.g., 'alpha' and 'log10_zeta' for Andrade).

    The decorated function receives theta_dict and structure_data, and must
    return modified structure_data.

    Example:
        @parameter_hook('alpha', 'log10_zeta')
        def apply_andrade_params(theta_dict, structure_data):
            # Reads both alpha and log10_zeta from theta_dict
            return modified_structure_data
    """
    def decorator(func):
        for param_id in param_ids:
            _PARAMETER_HOOKS[param_id] = func
        return func
    return decorator


def apply_parameters(theta_dict: Dict[str, float], structure_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Apply all parameters to structure data via registered hooks.

    Deduplicates hook calls: each function is only called once even if it's
    registered for multiple parameters (e.g., Andrade hook for alpha + log10_zeta).

    Args:
        theta_dict: Parameter values keyed by parameter ID
        structure_data: Cached structural arrays

    Returns:
        Modified structure_data with parameter overrides applied
    """
    modified_structure = structure_data.copy()
    called_functions = set()

    for param_id in theta_dict.keys():
        if param_id in _PARAMETER_HOOKS:
            hook_func = _PARAMETER_HOOKS[param_id]

            # Deduplicate: only call each function once
            if hook_func not in called_functions:
                modified_structure = hook_func(theta_dict, modified_structure)
                called_functions.add(hook_func)

    return modified_structure


# ============================================================================
# Parameter Hook Implementations
# ============================================================================

@parameter_hook('alpha', 'log10_zeta')
def apply_andrade_params(theta_dict: Dict[str, float], structure_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Apply Andrade rheology parameters (alpha, zeta).

    Stores Andrade parameters in structure_data for later use by forward model.
    Both alpha and log10_zeta must be present in theta_dict.
    """
    if 'alpha' not in theta_dict or 'log10_zeta' not in theta_dict:
        raise ValueError("Andrade rheology requires both 'alpha' and 'log10_zeta' parameters")

    alpha = theta_dict['alpha']
    log10_zeta = theta_dict['log10_zeta']
    zeta_pa = 10 ** log10_zeta

    # Store for forward model
    structure_data['rheology_type'] = 'andrade'
    structure_data['andrade_alpha'] = alpha
    structure_data['andrade_zeta_pa'] = zeta_pa
    structure_data['andrade_zeta_tp'] = zeta_pa ** (1.0 / alpha) if zeta_pa != 1.0 else 1.0

    return structure_data


@parameter_hook('log10_eta_Ih', 'log10_eta_HP', 'log10_eta_sil')
def apply_viscosity_params(theta_dict: Dict[str, float], structure_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Apply viscosity overrides for ice and silicate layers.

    Modifies eta_Pa_base array based on phase IDs.
    """
    eta_mod = structure_data['eta_Pa_base'].copy()
    phases = structure_data['phases']
    changeIndices = structure_data['changeIndices']
    n_layers = structure_data['n_layers']

    # Convert log viscosities
    eta_values = {}
    if 'log10_eta_Ih' in theta_dict:
        eta_values['Ih'] = 10 ** theta_dict['log10_eta_Ih']
    if 'log10_eta_HP' in theta_dict:
        eta_values['HP'] = 10 ** theta_dict['log10_eta_HP']
    if 'log10_eta_sil' in theta_dict:
        eta_values['sil'] = 10 ** theta_dict['log10_eta_sil']

    # Apply per-phase
    for i_layer in range(n_layers):
        s, e = changeIndices[i_layer], changeIndices[i_layer + 1]
        ph = int(phases[s])

        if ph == 1 and 'Ih' in eta_values:  # Ice Ih
            eta_mod[s:e] = eta_values['Ih']
        elif ph in (3, 5, 6) and 'HP' in eta_values:  # HP ices III, V, VI
            eta_mod[s:e] = eta_values['HP']
        elif ph >= 50 and ph < 100 and 'sil' in eta_values:  # Silicate
            eta_mod[s:e] = eta_values['sil']
        # Fe core, clathrates, liquid: keep baseline

    structure_data['eta_Pa'] = eta_mod
    return structure_data


@parameter_hook('log10_mu_Ih')
def apply_shear_modulus(theta_dict: Dict[str, float], structure_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Apply shear modulus override for Maxwell rheology.

    Typically used for Ice Ih layers in Maxwell viscoelastic model.
    """
    if 'log10_mu_Ih' in theta_dict:
        mu_Ih = 10 ** theta_dict['log10_mu_Ih']

        # Store for Maxwell rheology setup
        structure_data['rheology_type'] = 'maxwell'
        structure_data['mu_Ih_override'] = mu_Ih

        # Apply to Ice Ih layers
        mu_mod = structure_data['mu_Pa'].copy()
        phases = structure_data['phases']
        changeIndices = structure_data['changeIndices']
        n_layers = structure_data['n_layers']

        for i_layer in range(n_layers):
            s, e = changeIndices[i_layer], changeIndices[i_layer + 1]
            ph = int(phases[s])
            if ph == 1:  # Ice Ih
                mu_mod[s:e] = mu_Ih

        structure_data['mu_Pa'] = mu_mod

    return structure_data


@parameter_hook('Tb_K')
def apply_bottom_temperature(theta_dict: Dict[str, float], structure_data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Apply bottom temperature parameter.

    NOTE: Varying Tb_K requires pre-computed structure grid, not yet implemented.
    Future implementation will load from grid cache based on Tb_K value.

    For now, this raises an error to prevent misuse.
    """
    raise ValueError(
        "Parameter 'Tb_K' requires structure grid caching (not yet implemented). "
        "To vary Tb_K, generate a grid of structure caches with prepare_structure_variants.py "
        "using different Tb_K values, then implement grid interpolation in this hook."
    )


# ============================================================================
# Flexible Forward Model (Dict-based Parameters)
# ============================================================================

def forward_model_k2_flexible(
    theta_dict: Dict[str, float],
    structure_data: Dict[str, Any],
    return_heating: bool = False,
    arrhenius_params: Optional[Dict[str, Any]] = None
) -> Tuple[float, float, Optional[Dict[str, float]]]:
    """
    Flexible forward model using dict-based parameters.

    This is the extensible entry point that uses the parameter hook system.
    New parameters can be added by registering hooks - no changes needed here.

    Args:
        theta_dict: Parameters as dict (e.g., {'alpha': 0.3, 'log10_eta_Ih': 14.0})
        structure_data: Cached structural arrays from load_structure_cache()
        return_heating: If True, compute per-phase tidal heating (slower)
        arrhenius_params: Optional Arrhenius temperature-dependent viscosity params

    Returns:
        (Re_k2, Im_k2, perPhase_W) where perPhase_W is None if return_heating=False
        Returns (np.nan, np.nan, None) if TidalPy solver fails
    """
    # Import TidalPy components (lazy to avoid import overhead)
    try:
        from TidalPy.rheology import Andrade, Maxwell, Elastic
        from TidalPy.RadialSolver import radial_solver
        from TidalPy.RadialSolver import build_rs_input_from_data
        from TidalPy.tides.multilayer.heating import calc_radial_volumetric_tidal_heating_from_rs_solution
    except ImportError as e:
        log.error(f"TidalPy not available: {e}")
        return np.nan, np.nan, None

    # Apply parameter hooks
    modified_structure = apply_parameters(theta_dict, structure_data)

    # Extract modified arrays
    eta_mod = modified_structure.get('eta_Pa', modified_structure['eta_Pa_base'].copy())
    phases = modified_structure['phases']
    changeIndices = modified_structure['changeIndices']
    n_layers = modified_structure['n_layers']

    # Apply Arrhenius temperature-dependence if requested
    if arrhenius_params is not None:
        eta_mod = apply_arrhenius_viscosity(
            eta_mod,
            phases,
            changeIndices,
            n_layers,
            modified_structure.get('T_K'),
            arrhenius_params
        )

    # Build rheology models per layer based on applied parameters
    rheology_type = modified_structure.get('rheology_type', 'maxwell')

    if rheology_type == 'andrade':
        alpha = modified_structure['andrade_alpha']
        zeta_tp = modified_structure['andrade_zeta_tp']
        shear_models = []
        for rp in modified_structure['region_phases']:
            base_phase = rp.replace('_conv', '')
            if base_phase in ('0', 'Clath'):  # Liquid and clathrates: elastic
                shear_models.append(Elastic())
            else:
                shear_models.append(Andrade(args=(alpha, zeta_tp)))

    elif rheology_type == 'maxwell':
        shear_models = []
        for rp in modified_structure['region_phases']:
            base_phase = rp.replace('_conv', '')
            if base_phase in ('0', 'Clath'):
                shear_models.append(Elastic())
            else:
                shear_models.append(Maxwell())
    else:
        raise ValueError(f"Unknown rheology type '{rheology_type}'. Must be 'andrade' or 'maxwell'")

    # Bulk rheology: always elastic
    bulk_models = [Elastic() for _ in shear_models]

    # Run TidalPy radial solver
    try:
        build_data = build_rs_input_from_data(
            modified_structure['omega'],
            modified_structure['r_m'],
            modified_structure['rho'],
            modified_structure['K_Pa'],
            modified_structure.get('mu_Pa', structure_data['mu_Pa']),
            modified_structure['bulk_visc'],
            np.ascontiguousarray(eta_mod),
            modified_structure['layer_upper_radii'],
            modified_structure['layer_types'],
            tuple([False] * n_layers),  # is_static_bylayer
            tuple([False] * n_layers),  # is_incompressible_bylayer
            tuple(shear_models),
            tuple(bulk_models),
            perform_checks=False,
            warnings=False
        )

        result = radial_solver(
            build_data.radius_array,
            build_data.density_array,
            build_data.complex_bulk_modulus_array,
            build_data.complex_shear_modulus_array,
            build_data.frequency,
            build_data.planet_bulk_density,
            build_data.layer_types,
            build_data.is_static_bylayer,
            build_data.is_incompressible_bylayer,
            build_data.upper_radius_bylayer_array,
            degree_l=2,
            solve_for=('tidal',),
            verbose=False,
            raise_on_fail=False
        )

        if not result.success:
            return np.nan, np.nan, None

        # Extract k2
        k2 = complex(result.k)
        Re_k2 = k2.real
        Im_k2 = k2.imag

        # Compute per-phase heating if requested
        perPhase_W = None
        if return_heating and modified_structure['eccentricity'] > 0:
            perPhase_W = compute_heating(result, modified_structure)

        return Re_k2, Im_k2, perPhase_W

    except Exception as e:
        log.debug(f"Forward model failed: {e}")
        return np.nan, np.nan, None


# ============================================================================
# Legacy Array-based Forward Model (Backward Compatible)
# ============================================================================

def forward_model_k2(
    theta: np.ndarray,
    structure_data: Dict[str, Any],
    rheology: str = 'andrade',
    return_heating: bool = False,
    arrhenius_params: Optional[Dict[str, Any]] = None
) -> Tuple[float, float, Optional[Dict[str, float]]]:
    """
    Compute tidal Love number k2 (and optionally heating) for given parameters.

    LEGACY INTERFACE: Backward compatible array-based parameter format.
    Internally wraps forward_model_k2_flexible() with dict conversion.

    Args:
        theta: Parameter vector. Format depends on rheology:
            - 'andrade': [alpha, log10_zeta, log10_eta_Ih, log10_eta_HP, log10_eta_sil]
            - 'maxwell': [log10_eta_Ih, log10_eta_HP, log10_eta_sil]
            - With Tb_K variation: append Tb_K as final parameter
        structure_data: Cached structural arrays from load_structure_cache()
        rheology: 'andrade' or 'maxwell'
        return_heating: If True, compute per-phase tidal heating (slower)
        arrhenius_params: Optional dict with {'activation_energy_kJ_mol': Dict[phase: E],
                                             'reference_temp_K': float, ...}

    Returns:
        (Re_k2, Im_k2, perPhase_W) where perPhase_W is None if return_heating=False
        Returns (np.nan, np.nan, None) if TidalPy solver fails
    """
    # Convert theta array to dict based on rheology
    if rheology == 'andrade':
        if len(theta) == 5:
            theta_dict = {
                'alpha': theta[0],
                'log10_zeta': theta[1],
                'log10_eta_Ih': theta[2],
                'log10_eta_HP': theta[3],
                'log10_eta_sil': theta[4]
            }
        else:
            raise ValueError(f"Andrade rheology expects 5 parameters, got {len(theta)}")

    elif rheology == 'maxwell':
        if len(theta) == 3:
            theta_dict = {
                'log10_eta_Ih': theta[0],
                'log10_eta_HP': theta[1],
                'log10_eta_sil': theta[2]
            }
        else:
            raise ValueError(f"Maxwell rheology expects 3 parameters, got {len(theta)}")
    else:
        raise ValueError(f"Unknown rheology '{rheology}'. Must be 'andrade' or 'maxwell'")

    # Call flexible version
    return forward_model_k2_flexible(
        theta_dict,
        structure_data,
        return_heating=return_heating,
        arrhenius_params=arrhenius_params
    )


def compute_heating(
    rs_solution,
    structure_data: Dict[str, Any]
) -> Dict[str, float]:
    """
    Compute per-phase tidal heating from TidalPy radial solver solution.

    Args:
        rs_solution: TidalPy RadialSolver result object
        structure_data: Cached structural arrays

    Returns:
        Dict mapping phase names (e.g., 'Ih', 'III', 'V', 'VI', 'Sil') to total
        heating in that phase (Watts)
    """
    from TidalPy.tides.multilayer.heating import calc_radial_volumetric_tidal_heating_from_rs_solution
    from PlanetProfile.Thermodynamics.Electrical import PhaseConv

    # Compute volumetric heating profile
    heating_profile = calc_radial_volumetric_tidal_heating_from_rs_solution(
        structure_data['eccentricity'],
        structure_data['omega'],
        structure_data['a_m'],
        structure_data['host_mass'],
        rs_solution,
        perform_checks=False
    )

    result_radii = np.asarray(rs_solution.radius_array)
    r_m_local = structure_data['r_m']
    phases = structure_data['phases']
    changeIndices = structure_data['changeIndices']
    n_layers = structure_data['n_layers']
    phase_map = structure_data['phase_map']

    # Integrate heating over each layer
    perPhase_W = {}
    for i_layer in range(n_layers):
        s_idx = changeIndices[i_layer]
        e_idx = changeIndices[i_layer + 1]
        layer_phase = int(phases[s_idx])

        # Map numeric phase to string
        phase_str = phase_map.get(layer_phase, PhaseConv(layer_phase, liq='0'))

        # Find radii in this layer
        r_lo = r_m_local[s_idx]
        r_hi = r_m_local[e_idx - 1]
        mask = (result_radii >= r_lo - 1.0) & (result_radii <= r_hi + 1.0)

        if np.any(mask):
            lr = result_radii[mask]
            lh = heating_profile[mask]

            # Integrate heating: ∫ heating * 4πr² dr
            if len(lr) > 1:
                total_power = np.trapz(lh * 4.0 * np.pi * lr**2, lr)
            else:
                # Single point: use layer volume
                V = (4.0/3.0) * np.pi * (r_hi**3 - r_lo**3)
                total_power = lh[0] * V

            # Accumulate (same phase may appear in multiple layers)
            perPhase_W[phase_str] = perPhase_W.get(phase_str, 0) + total_power

    return perPhase_W


def apply_arrhenius_viscosity(
    eta_base: np.ndarray,
    phases: np.ndarray,
    changeIndices: np.ndarray,
    n_layers: int,
    T_K: Optional[np.ndarray],
    arrhenius_params: Dict[str, Any]
) -> np.ndarray:
    """
    Apply Arrhenius temperature-dependent viscosity scaling.

    η(T) = η_ref * exp(E/R * (1/T - 1/T_ref))

    Args:
        eta_base: Base viscosity array (Pa·s)
        phases: Phase ID array
        changeIndices: Layer boundary indices
        n_layers: Number of layers
        T_K: Temperature profile (K). If None, returns eta_base unchanged
        arrhenius_params: Dict with:
            - 'activation_energy_kJ_mol': Dict[phase_name: E_kJ_mol]
            - 'reference_temp_K': float (default 273.15 K)
            - 'R_J_mol_K': float (default 8.314 J/mol/K)

    Returns:
        Modified viscosity array with Arrhenius scaling applied
    """
    if T_K is None:
        log.warning("Temperature profile not available, Arrhenius viscosity not applied")
        return eta_base

    eta_mod = eta_base.copy()
    activation_energies = arrhenius_params.get('activation_energy_kJ_mol', {})
    T_ref = arrhenius_params.get('reference_temp_K', 273.15)
    R = arrhenius_params.get('R_J_mol_K', 8.314)

    # Phase name mapping (numeric -> string)
    phase_names = {1: 'Ih', 3: 'III', 5: 'V', 6: 'VI'}

    for i_layer in range(n_layers):
        s, e = changeIndices[i_layer], changeIndices[i_layer + 1]
        ph = int(phases[s])

        # Check if this phase has Arrhenius parameters
        phase_name = phase_names.get(ph)
        if phase_name is None:
            continue  # Not an ice phase with Arrhenius

        E_kJ_mol = activation_energies.get(phase_name)
        if E_kJ_mol is None:
            continue  # No activation energy specified for this phase

        # Apply Arrhenius scaling
        E_J_mol = E_kJ_mol * 1000.0
        T_layer = T_K[s:e]

        # Avoid division by zero or negative temperatures
        T_layer = np.maximum(T_layer, 1.0)

        # η(T) = η_ref * exp(E/R * (1/T - 1/T_ref))
        exponent = (E_J_mol / R) * (1.0 / T_layer - 1.0 / T_ref)
        eta_mod[s:e] *= np.exp(exponent)

    return eta_mod


def evaluate_heating_on_posterior(
    samples: np.ndarray,
    structure_data: Dict[str, Any],
    rheology: str = 'andrade',
    n_eval: Optional[int] = None,
    random_state: int = 42,
    arrhenius_params: Optional[Dict[str, Any]] = None
) -> Tuple[np.ndarray, np.ndarray, List[Dict[str, float]]]:
    """
    Re-evaluate forward model on posterior samples to get k2 and heating.

    This is called after MCMC/SBI completes to compute heating for a subset
    of posterior samples (since heating calculation is slow).

    Args:
        samples: Posterior samples (n_samples, n_params)
        structure_data: Cached structural arrays
        rheology: 'andrade' or 'maxwell'
        n_eval: Number of samples to evaluate (default: all)
        random_state: Random seed for subsampling
        arrhenius_params: Optional Arrhenius parameters

    Returns:
        (sample_indices, k2_results, heating_results)
        - sample_indices: Indices of evaluated samples (n_eval,)
        - k2_results: Array of (Re_k2, Im_k2) for each sample (n_eval, 2)
        - heating_results: List of per-phase heating dicts (n_eval,)
    """
    n_samples = len(samples)
    if n_eval is None or n_eval > n_samples:
        n_eval = n_samples

    # Subsample posterior
    rng = np.random.RandomState(random_state)
    if n_eval < n_samples:
        idx = rng.choice(n_samples, n_eval, replace=False)
        idx.sort()
    else:
        idx = np.arange(n_samples)

    log.info(f"Evaluating {n_eval} posterior samples for heating...")

    k2_results = []
    heating_results = []

    for i, si in enumerate(idx):
        theta = samples[si]
        Re_k2, Im_k2, perPhase_W = forward_model_k2(
            theta,
            structure_data,
            rheology=rheology,
            return_heating=True,
            arrhenius_params=arrhenius_params
        )
        k2_results.append((Re_k2, Im_k2))
        heating_results.append(perPhase_W if perPhase_W is not None else {})

        if (i + 1) % 100 == 0:
            log.info(f"  {i+1}/{n_eval} samples evaluated")

    log.info(f"Heating evaluation complete: {n_eval} samples")

    return idx, np.array(k2_results), heating_results


def create_log_likelihood(
    observables: Dict[str, Tuple[float, float]],
    structure_data: Dict[str, Any],
    rheology: str = 'andrade',
    arrhenius_params: Optional[Dict[str, Any]] = None
) -> callable:
    """
    Create log-likelihood function for MCMC/SBI sampling.

    Returns a function that accepts theta and returns log-likelihood value.

    Args:
        observables: Dict of {observable_name: (value, uncertainty)}
            Supported observables: 'Re_k2', 'Im_k2', 'k2', 'abs_Im_k2'
        structure_data: Cached structural arrays
        rheology: 'andrade' or 'maxwell'
        arrhenius_params: Optional Arrhenius parameters

    Returns:
        Callable log_likelihood(theta) -> float
    """
    def log_likelihood(theta):
        """Gaussian log-likelihood on tidal observables."""
        Re_k2, Im_k2, _ = forward_model_k2(
            theta,
            structure_data,
            rheology=rheology,
            return_heating=False,
            arrhenius_params=arrhenius_params
        )

        if np.isnan(Re_k2):
            return -1e30  # TidalPy solver failed

        # Compute chi-squared
        chi2 = 0.0

        if 'Re_k2' in observables:
            obs_val, obs_err = observables['Re_k2']
            chi2 += ((Re_k2 - obs_val) / obs_err)**2

        if 'Im_k2' in observables:
            obs_val, obs_err = observables['Im_k2']
            chi2 += ((Im_k2 - obs_val) / obs_err)**2

        if 'abs_Im_k2' in observables:
            obs_val, obs_err = observables['abs_Im_k2']
            chi2 += ((abs(Im_k2) - obs_val) / obs_err)**2

        if 'k2' in observables:
            # Magnitude constraint: |k2| = sqrt(Re^2 + Im^2)
            obs_val, obs_err = observables['k2']
            k2_mag = np.sqrt(Re_k2**2 + Im_k2**2)
            chi2 += ((k2_mag - obs_val) / obs_err)**2

        return -0.5 * chi2

    return log_likelihood
