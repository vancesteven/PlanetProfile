"""
PlanetProfile Inference Module

Bayesian inference for constraining planetary interior parameters using
observed data (tidal Love numbers, magnetic induction, etc.).

Supports:
- MCMC (Markov Chain Monte Carlo) via pocoMC
- SBI (Simulation-Based Inference) via neural posterior estimation [future]

Main entry point: run_inference() from inference_core.py
CLI entry point: run_inference_cli.py for HPC/subprocess execution

Example usage:
    >>> from PlanetProfile.Inference import InferenceConfig, run_inference
    >>> config = InferenceConfig(
    ...     mode='mcmc',
    ...     bodyname='Titan',
    ...     param_space={'alpha': {'prior_type': 'uniform', 'bounds': [0.2, 0.4]}},
    ...     observables={'Re_k2': (0.608, 0.048)},
    ...     sampler_settings={'n_effective': 500}
    ... )
    >>> result = run_inference(config)
    >>> print(result.get_summary_stats())
"""

from .inference_core import (
    InferenceConfig,
    InferenceResult,
    run_inference,
    validate_config
)

from .forward_models import (
    forward_model_k2,
    compute_heating,
    evaluate_heating_on_posterior,
    create_log_likelihood,
    apply_arrhenius_viscosity
)

from .structure_cache import (
    load_structure_cache,
    build_structure_from_pptest,
    extract_structure_from_planet,
    save_structure_cache,
    build_structure_grid,
    validate_structure_cache
)

from .mcmc_analysis import (
    reanalyze_k2_from_pickle,
    plot_k2_posteriors,
)

# Lazy import runners to avoid dependency chains
def _get_mcmc_runner():
    from .mcmc_runner import MCMCRunner
    return MCMCRunner

def _get_sbi_runner():
    from .sbi_runner import SBIRunner
    return SBIRunner

__all__ = [
    # Core API
    'InferenceConfig',
    'InferenceResult',
    'run_inference',
    'validate_config',

    # Forward models
    'forward_model_k2',
    'compute_heating',
    'evaluate_heating_on_posterior',
    'create_log_likelihood',
    'apply_arrhenius_viscosity',

    # Structure caching
    'load_structure_cache',
    'build_structure_from_pptest',
    'extract_structure_from_planet',
    'save_structure_cache',
    'build_structure_grid',
    'validate_structure_cache',

    # Post-hoc analysis
    'reanalyze_k2_from_pickle',
    'plot_k2_posteriors',

    # Runners (lazy imports)
    '_get_mcmc_runner',
    '_get_sbi_runner',
]
