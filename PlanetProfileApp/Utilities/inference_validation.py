"""
Validation functions for the Inference tab.

Returns validation tuples: (is_valid: bool, message: str, severity: str|None)

Severity levels:
- 'error': Invalid, blocks execution
- 'warning': Valid but potentially problematic
- 'info': Valid with helpful context
- None: Valid, no message

Pattern adapted from app_helpers.py validation functions
"""
from typing import Tuple, Optional, Dict, Any


def validate_prior_bounds(
    param_name: str,
    prior_type: str,
    bounds: Optional[Tuple[float, float]] = None,
    mean: Optional[float] = None,
    std: Optional[float] = None
) -> Tuple[bool, str, Optional[str]]:
    """
    Validate prior specification for a parameter.

    Args:
        param_name: Parameter name (e.g., 'alpha', 'log10_eta_Ih')
        prior_type: 'uniform', 'normal', 'log-uniform'
        bounds: (min, max) for uniform/log-uniform priors
        mean: Mean for normal prior
        std: Std deviation for normal prior

    Returns:
        (is_valid, message, severity)
    """
    # Check prior type
    valid_types = ['uniform', 'normal', 'log-uniform']
    if prior_type not in valid_types:
        return (False, f"Invalid prior type '{prior_type}'. Must be one of {valid_types}.",
                'error')

    # Validate uniform/log-uniform priors
    if prior_type in ['uniform', 'log-uniform']:
        if bounds is None or len(bounds) != 2:
            return (False, f"Uniform prior requires bounds (min, max).", 'error')

        min_val, max_val = bounds
        if min_val >= max_val:
            return (False, f"Prior bounds invalid: min ({min_val}) >= max ({max_val}).",
                    'error')

        # Physical bounds checks
        if param_name == 'alpha':
            if not (0 < min_val < 1 and 0 < max_val < 1):
                return (False, "Andrade exponent alpha must be in (0, 1).", 'error')
            if not (0.2 <= min_val <= 0.4 and 0.2 <= max_val <= 0.4):
                return (True, "Typical Andrade alpha range is [0.2, 0.4] (Petricca 2025).",
                        'warning')

        elif 'eta' in param_name.lower() or 'log10_eta' in param_name.lower():
            # Viscosity bounds
            if 'log10' in param_name:
                if not (5 <= min_val <= 25 and 5 <= max_val <= 25):
                    return (False, "Log10(viscosity) must be in [5, 25] (Pa·s).", 'error')
            else:
                if min_val < 0 or max_val < 0:
                    return (False, "Viscosity must be positive.", 'error')

        elif param_name == 'log10_zeta':
            if not (-5 <= min_val <= 5 and -5 <= max_val <= 5):
                return (True, "Typical Andrade zeta range is [-2, 2] (Petricca 2025).",
                        'warning')

        elif 'Tb' in param_name or 'T_b' in param_name:
            # Basal temperature
            if not (200 <= min_val <= 300 and 200 <= max_val <= 300):
                return (False, "Basal temperature must be in [200, 300] K for Titan.",
                        'error')

    # Validate normal priors
    elif prior_type == 'normal':
        if mean is None or std is None:
            return (False, "Normal prior requires mean and std.", 'error')

        if std <= 0:
            return (False, f"Standard deviation must be positive (got {std}).", 'error')

        if std > abs(mean):
            return (True, f"Large std ({std}) relative to mean ({mean}). Prior may be too wide.",
                    'warning')

    return (True, "Prior bounds valid.", None)


def validate_observables(
    observables: Dict[str, Tuple[float, float]]
) -> Tuple[bool, str, Optional[str]]:
    """
    Validate observable constraints.

    Args:
        observables: Dict of {name: (value, uncertainty)}

    Returns:
        (is_valid, message, severity)
    """
    if not observables:
        return (False, "At least one observable is required.", 'error')

    for obs_name, (value, uncertainty) in observables.items():
        # Check for valid observable types
        valid_obs = ['Re_k2', 'Im_k2', 'k2', 'h2', 'l2', 'Bx', 'By', 'Bz']
        if obs_name not in valid_obs:
            return (True, f"Observable '{obs_name}' not in standard list: {valid_obs}.",
                    'warning')

        # Check value is numeric
        if not isinstance(value, (int, float)):
            return (False, f"Observable '{obs_name}' value must be numeric (got {type(value)}).",
                    'error')

        # Check uncertainty is positive
        if uncertainty <= 0:
            return (False, f"Observable '{obs_name}' uncertainty must be positive (got {uncertainty}).",
                    'error')

        # Physical bounds for k2
        if 'k2' in obs_name.lower():
            if not (0 < value < 2):
                return (False, f"Tidal Love number k2 must be in (0, 2). Got {value}.",
                        'error')

            if 'Re' in obs_name:
                if not (0.5 < value < 0.7):
                    return (True, f"Typical Titan Re(k2) is ~0.608 (Petricca 2025). Got {value}.",
                            'info')

            if 'Im' in obs_name:
                if value < 0:
                    return (False, f"Im(k2) constraint should use absolute value. Got {value}.",
                            'error')
                if not (0.1 < value < 0.2):
                    return (True, f"Typical Titan |Im(k2)| is ~0.135 (Petricca 2025). Got {value}.",
                            'info')

    return (True, "Observables valid.", None)


def validate_sampler_settings(
    sampler_type: str,
    settings: Dict[str, Any]
) -> Tuple[bool, str, Optional[str]]:
    """
    Validate sampler hyperparameters.

    Args:
        sampler_type: 'mcmc' or 'sbi'
        settings: Dict of sampler-specific settings

    Returns:
        (is_valid, message, severity)
    """
    if sampler_type not in ['mcmc', 'sbi']:
        return (False, f"Invalid sampler type '{sampler_type}'. Must be 'mcmc' or 'sbi'.",
                'error')

    if sampler_type == 'mcmc':
        # Validate MCMC settings
        n_eff = settings.get('n_effective', 500)
        if not isinstance(n_eff, int) or n_eff < 100:
            return (False, f"n_effective must be integer >= 100. Got {n_eff}.",
                    'error')

        if n_eff < 500:
            return (True, f"n_effective={n_eff} may give under-sampled posterior. Recommend >= 500.",
                    'warning')

        random_state = settings.get('random_state', 42)
        if not isinstance(random_state, int) or random_state < 0:
            return (False, f"random_state must be non-negative integer. Got {random_state}.",
                    'error')

        n_reeval = settings.get('n_reeval', 500)
        if not isinstance(n_reeval, int) or n_reeval < 10:
            return (False, f"n_reeval must be integer >= 10. Got {n_reeval}.",
                    'error')

        if n_reeval > n_eff:
            return (True, f"n_reeval ({n_reeval}) > n_effective ({n_eff}). Will use all samples.",
                    'info')

    elif sampler_type == 'sbi':
        # Validate SBI settings
        num_sims = settings.get('num_simulations', 10000)
        if not isinstance(num_sims, int) or num_sims < 1000:
            return (False, f"num_simulations must be integer >= 1000. Got {num_sims}.",
                    'error')

        if num_sims > 100000:
            return (True, f"num_simulations={num_sims} will take significant time (~hours).",
                    'warning')

        num_rounds = settings.get('num_rounds', 3)
        if not isinstance(num_rounds, int) or num_rounds < 1:
            return (False, f"num_rounds must be positive integer. Got {num_rounds}.",
                    'error')

    return (True, "Sampler settings valid.", None)


def validate_structure_path(
    structure_path: Optional[str],
    bodyname: str
) -> Tuple[bool, str, Optional[str]]:
    """
    Validate structure cache path exists and matches body.

    Args:
        structure_path: Path to structure cache file (.pkl or .txt)
        bodyname: Planet name (e.g., 'Titan')

    Returns:
        (is_valid, message, severity)
    """
    if structure_path is None:
        return (False, f"Structure file required. Run PlanetProfile for {bodyname} first.",
                'error')

    import os
    if not os.path.exists(structure_path):
        return (False, f"Structure file not found: {structure_path}",
                'error')

    # Check file extension
    if not (structure_path.endswith('.pkl') or structure_path.endswith('.txt')):
        return (False, f"Structure file must be .pkl or .txt. Got: {structure_path}",
                'error')

    # Check bodyname in path (basic validation)
    if bodyname.lower() not in os.path.basename(structure_path).lower():
        return (True, f"Structure file name doesn't contain '{bodyname}'. Verify it's correct.",
                'warning')

    return (True, "Structure file valid.", None)


def validate_full_config(
    param_space: Dict[str, Dict],
    observables: Dict[str, Tuple[float, float]],
    sampler_type: str,
    sampler_settings: Dict[str, Any],
    structure_path: Optional[str],
    bodyname: str
) -> Tuple[bool, str, Optional[str]]:
    """
    Validate entire inference configuration.

    Args:
        param_space: Dict of {param: {prior_type, bounds/mean/std}}
        observables: Dict of {obs: (value, uncertainty)}
        sampler_type: 'mcmc' or 'sbi'
        sampler_settings: Sampler hyperparameters
        structure_path: Path to structure cache
        bodyname: Planet name

    Returns:
        (is_valid, message, severity)
    """
    # Check parameter space is not empty
    if not param_space:
        return (False, "Parameter space is empty. Select at least one parameter to vary.",
                'error')

    # Validate each prior
    for param_name, prior_spec in param_space.items():
        prior_type = prior_spec.get('prior_type', 'uniform')
        bounds = prior_spec.get('bounds')
        mean = prior_spec.get('mean')
        std = prior_spec.get('std')

        is_valid, message, severity = validate_prior_bounds(
            param_name, prior_type, bounds, mean, std
        )

        if not is_valid:
            return (False, f"Parameter '{param_name}': {message}", severity)

    # Validate observables
    is_valid, message, severity = validate_observables(observables)
    if not is_valid:
        return (False, message, severity)

    # Validate sampler settings
    is_valid, message, severity = validate_sampler_settings(sampler_type, sampler_settings)
    if not is_valid:
        return (False, message, severity)

    # Validate structure
    is_valid, message, severity = validate_structure_path(structure_path, bodyname)
    if not is_valid:
        return (False, message, severity)

    return (True, "Configuration valid. Ready to run inference.", None)
