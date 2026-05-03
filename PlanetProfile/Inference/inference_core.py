"""
Core inference infrastructure for PlanetProfile.

Defines configuration dataclasses, result structures, and main API entry point.
Designed for both interactive GUI use and headless HPC execution.

Author: PlanetProfile Team
Date: 2026-04-29
"""
import json
import pickle
import hashlib
import numpy as np
from dataclasses import dataclass, field, asdict
from typing import Dict, Tuple, List, Optional, Any, Callable
from pathlib import Path


@dataclass
class InferenceConfig:
    """
    Configuration for Bayesian inference run.

    Attributes:
        mode: 'mcmc' or 'sbi'
        bodyname: Planet name (e.g., 'Titan')
        param_space: Dict of {param_name: {prior_type, bounds/mean/std}}
        observables: Dict of {observable_name: (value, uncertainty)}
        sampler_settings: Sampler-specific hyperparameters
        structure_cache_path: Path to pre-computed Planet structure (optional)
        random_state: Random seed for reproducibility
        metadata: Additional run metadata (tags, description, etc.)
    """
    mode: str  # 'mcmc' or 'sbi'
    bodyname: str
    param_space: Dict[str, Dict[str, Any]]
    observables: Dict[str, Tuple[float, float]]
    sampler_settings: Dict[str, Any]
    structure_cache_path: Optional[str] = None
    random_state: int = 42
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate configuration after initialization."""
        # Mode validation
        if self.mode not in ['mcmc', 'sbi']:
            raise ValueError(f"mode must be 'mcmc' or 'sbi', got '{self.mode}'")

        # Parameter space validation
        if not self.param_space:
            raise ValueError("param_space cannot be empty")

        # Observables validation
        if not self.observables:
            raise ValueError("observables cannot be empty")

        for obs_name, (value, uncertainty) in self.observables.items():
            if uncertainty <= 0:
                raise ValueError(f"Observable '{obs_name}' uncertainty must be positive, got {uncertainty}")

    def to_dict(self) -> Dict:
        """Convert to dictionary for serialization."""
        return asdict(self)

    def to_json(self, filepath: str) -> None:
        """Save configuration to JSON file."""
        with open(filepath, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def from_dict(cls, data: Dict) -> 'InferenceConfig':
        """Load configuration from dictionary."""
        return cls(**data)

    @classmethod
    def from_json(cls, filepath: str) -> 'InferenceConfig':
        """Load configuration from JSON file."""
        with open(filepath, 'r') as f:
            data = json.load(f)
        return cls.from_dict(data)

    def generate_hash(self) -> str:
        """
        Generate unique hash for this configuration.
        Used for cache keys and result identification.
        """
        # Include all computational parameters
        hash_data = {
            'mode': self.mode,
            'bodyname': self.bodyname,
            'param_space': sorted(self.param_space.items()),
            'observables': sorted(self.observables.items()),
            'sampler_settings': sorted(self.sampler_settings.items()),
            'random_state': self.random_state
        }
        hash_str = json.dumps(hash_data, sort_keys=True)
        return hashlib.md5(hash_str.encode()).hexdigest()[:16]


@dataclass
class InferenceResult:
    """
    Results from Bayesian inference run.

    Attributes:
        config: Original configuration used for this run
        samples: Posterior samples array (n_samples, n_params)
        log_likelihoods: Log-likelihood values for each sample
        param_names: Parameter names in order (matches samples columns)
        param_labels: LaTeX-formatted labels for plotting
        k2_results: Re(k2) and Im(k2) for evaluated samples (n_eval, 2)
        cmr2_results: C/MR² values for all posterior samples (n_samples,)
        D_ocean_results: Ocean thickness (km) for all posterior samples (n_samples,)
        D_iceIh_results: Ice Ih shell thickness (km) for all posterior samples (n_samples,)
        heating_results: Per-phase heating for evaluated samples
        convergence_metrics: R-hat, ESS, acceptance rate, etc.
        metadata: Runtime info, versions, timestamps, etc.
    """
    config: InferenceConfig
    samples: np.ndarray
    log_likelihoods: np.ndarray
    param_names: List[str]
    param_labels: List[str]
    k2_results: Optional[np.ndarray] = None
    cmr2_results: Optional[np.ndarray] = None
    D_ocean_results: Optional[np.ndarray] = None
    D_iceIh_results: Optional[np.ndarray] = None
    heating_results: Optional[List[Dict[str, float]]] = None
    convergence_metrics: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self):
        """Validate result structure after initialization."""
        n_samples, n_params = self.samples.shape

        # Validate dimensions
        if len(self.log_likelihoods) != n_samples:
            raise ValueError(f"log_likelihoods length ({len(self.log_likelihoods)}) "
                           f"must match samples ({n_samples})")

        if len(self.param_names) != n_params:
            raise ValueError(f"param_names length ({len(self.param_names)}) "
                           f"must match samples dimensions ({n_params})")

        if len(self.param_labels) != n_params:
            raise ValueError(f"param_labels length ({len(self.param_labels)}) "
                           f"must match samples dimensions ({n_params})")

        # Validate k2_results if present
        if self.k2_results is not None:
            if self.k2_results.shape[1] != 2:
                raise ValueError(f"k2_results must have shape (n_eval, 2), "
                               f"got {self.k2_results.shape}")

        # Validate heating_results if present
        if self.heating_results is not None:
            if not isinstance(self.heating_results, list):
                raise ValueError("heating_results must be a list of dicts")

    def save(self, filepath: str) -> None:
        """
        Save result to pickle file.

        Args:
            filepath: Output path (.pkl extension)
        """
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)

        with open(filepath, 'wb') as f:
            pickle.dump(self, f, protocol=pickle.HIGHEST_PROTOCOL)

    @classmethod
    def load(cls, filepath: str) -> 'InferenceResult':
        """
        Load result from pickle file.

        Args:
            filepath: Input path (.pkl file)

        Returns:
            InferenceResult object
        """
        with open(filepath, 'rb') as f:
            result = pickle.load(f)

        if not isinstance(result, cls):
            raise TypeError(f"Loaded object is not InferenceResult, got {type(result)}")

        return result

    def get_summary_stats(self) -> Dict[str, np.ndarray]:
        """
        Compute summary statistics for each parameter.

        Returns:
            Dict with keys: median, mean, std, q16, q84 (1σ), q025, q975 (2σ)
            Each value is array of length n_params
        """
        return {
            'median': np.median(self.samples, axis=0),
            'mean': np.mean(self.samples, axis=0),
            'std': np.std(self.samples, axis=0),
            'q16': np.percentile(self.samples, 16, axis=0),
            'q84': np.percentile(self.samples, 84, axis=0),
            'q025': np.percentile(self.samples, 2.5, axis=0),
            'q975': np.percentile(self.samples, 97.5, axis=0)
        }

    def get_best_fit(self) -> np.ndarray:
        """
        Get best-fit parameters (maximum likelihood).

        Returns:
            1D array of parameter values
        """
        best_idx = np.argmax(self.log_likelihoods)
        return self.samples[best_idx]


def run_inference(
    config: InferenceConfig,
    progress_callback: Optional[Callable[[Dict], None]] = None
) -> InferenceResult:
    """
    Main entry point for Bayesian inference.

    Dispatches to appropriate runner (MCMC or SBI) based on config.mode.

    Args:
        config: Inference configuration
        progress_callback: Optional callback for progress updates.
                          Called with dict containing:
                          - iteration: current iteration
                          - n_total: total iterations
                          - acceptance_rate: current acceptance rate (MCMC)
                          - n_samples: current number of posterior samples

    Returns:
        InferenceResult object with posterior samples and diagnostics

    Raises:
        ValueError: If config.mode is invalid
        ImportError: If required sampler library not installed
    """
    if config.mode == 'mcmc':
        # Import here to avoid dependency if not using MCMC
        try:
            from .mcmc_runner import MCMCRunner
        except ImportError as e:
            raise ImportError(
                f"MCMC mode requires pocoMC library. Install with: pip install pocomc\n"
                f"Original error: {e}"
            )

        runner = MCMCRunner(config)
        return runner.run(progress_callback)

    elif config.mode == 'sbi':
        # Import here to avoid dependency if not using SBI
        try:
            from .sbi_runner import SBIRunner
        except ImportError as e:
            raise ImportError(
                f"SBI mode requires sbi library. Install with: pip install sbi torch\n"
                f"Original error: {e}"
            )

        runner = SBIRunner(config)
        return runner.run(progress_callback)

    else:
        raise ValueError(f"Invalid mode '{config.mode}'. Must be 'mcmc' or 'sbi'.")


def validate_config(config: InferenceConfig) -> Tuple[bool, List[str]]:
    """
    Validate inference configuration for common issues.

    Args:
        config: Configuration to validate

    Returns:
        (is_valid, error_messages) tuple
    """
    errors = []

    # Check mode
    if config.mode not in ['mcmc', 'sbi']:
        errors.append(f"Invalid mode '{config.mode}'. Must be 'mcmc' or 'sbi'.")

    # Check parameter space
    if not config.param_space:
        errors.append("Parameter space is empty.")

    for param, spec in config.param_space.items():
        prior_type = spec.get('prior_type', 'uniform')

        if prior_type not in ['uniform', 'normal', 'log-uniform']:
            errors.append(f"Parameter '{param}': invalid prior_type '{prior_type}'")

        if prior_type in ['uniform', 'log-uniform']:
            bounds = spec.get('bounds')
            if bounds is None or len(bounds) != 2:
                errors.append(f"Parameter '{param}': missing or invalid bounds")
            elif bounds[0] >= bounds[1]:
                errors.append(f"Parameter '{param}': bounds min >= max")

        elif prior_type == 'normal':
            if spec.get('mean') is None or spec.get('std') is None:
                errors.append(f"Parameter '{param}': missing mean or std for normal prior")
            elif spec.get('std', 1) <= 0:
                errors.append(f"Parameter '{param}': std must be positive")

    # Check observables
    if not config.observables:
        errors.append("Observables dictionary is empty.")

    for obs, (value, uncertainty) in config.observables.items():
        if uncertainty <= 0:
            errors.append(f"Observable '{obs}': uncertainty must be positive")

    # Check sampler settings
    if config.mode == 'mcmc':
        n_eff = config.sampler_settings.get('n_effective', 500)
        if n_eff < 100:
            errors.append(f"MCMC n_effective ({n_eff}) should be >= 100")

    elif config.mode == 'sbi':
        n_sims = config.sampler_settings.get('num_simulations', 10000)
        if n_sims < 1000:
            errors.append(f"SBI num_simulations ({n_sims}) should be >= 1000")

    # Check structure cache if provided
    if config.structure_cache_path:
        cache_path = Path(config.structure_cache_path)
        if not cache_path.exists():
            errors.append(f"Structure cache not found: {config.structure_cache_path}")

    return (len(errors) == 0, errors)
