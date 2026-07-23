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
    param_groups: Dict[str, List[str]] = field(default_factory=dict)
    """
    param_groups: Maps a *group key* (sampled scalar) to a list of *member*
    parameter names that all receive the same sampled value.

    Example::

        param_groups = {'log10_eta_HP': ['log10_eta_III', 'log10_eta_V', 'log10_eta_VI']}

    The group key must appear in ``param_space`` and must NOT appear in
    ``param_space`` for any of its members.  The forward model receives the
    expanded dict (group key → sampled value AND each member → same value).
    """
    fixed_params: Dict[str, float] = field(default_factory=dict)
    """
    fixed_params: Parameters passed to the forward model as constants (not
    sampled).  Must NOT appear in ``param_space`` or as group keys in
    ``param_groups``.

    Example::

        fixed_params = {'Tb_K': 250.965}
    """
    derived_params: Dict[str, Any] = field(default_factory=dict)
    """
    derived_params: Parameters whose value is computed from other sampled
    parameters via a physical constraint, then optionally bounds-checked.
    Stored at parse time; the derivation itself is performed by the runner.

    Schema (v2)::

        derived_params = {
            'rho_sil_kgm3': {
                'derivation': 'mass_conservation',
                'bounds': [2200.0, 3500.0],
                'reject_if_outside_bounds': True,
            },
        }

    Currently the only registered derivation is ``mass_conservation`` for
    ``rho_sil_kgm3``, wired in ``MCMCRunner._make_flexible_log_likelihood``
    (Phase C1 Stage 2).  Unknown derivations are tolerated at parse time
    and ignored by the runner.
    """

    ocean_overrides: Dict[str, Any] = field(default_factory=dict)
    """
    ocean_overrides: Cache-builder-only knob carried through the unified
    config schema. Consumed by ``cache_builder.build_phase_c1_cache`` to
    swap ``Planet.Ocean.comp`` / ``Planet.Ocean.wOcean_ppt`` (e.g. NaCl
    100 ppt for Callisto). The MCMC runner ignores this field — it lives
    here only so a single body config file can drive both the cache build
    and the inference run without schema-skew between the two stages.

    Example::

        ocean_overrides = {'comp': 'NaCl', 'wOcean_ppt': 100.0}
    """
    # Arrhenius viscosity params (e.g. {'E_Ih_J_per_mol': 60000}); empty = no Arrhenius correction.
    arrhenius_params: Dict[str, Any] = field(default_factory=dict)
    # Importable module path for the PP body template (e.g. 'PlanetProfile.Test.PPTest48'); used by plot_structure_wedge_pp.
    planet_template_module: Optional[str] = None
    induction_bounds: Dict[str, Any] = field(default_factory=dict)
    """
    induction_bounds: One-sided constraints on the complex induction response
    Ae per excitation label — support cuts, NOT Gaussian observables (they
    never become SBI x-channels; the flow's training support and the MCMC
    likelihood support both exclude the violating region).

    Schema (ratified 2026-07-12, Europa campaign)::

        induction_bounds = {
            'synodic': {
                'amp_min': 0.7,      # reject when |Ae| < amp_min
                'im_abs_max': 0.4,   # reject when |Im(Ae)| > im_abs_max
            },
        }

    Both keys optional per label. Enforced as a hard rejection (the -1e30
    sentinel) in MCMCRunner's likelihood and as training-time support
    rejection in generate_sbi_dataset (counted in n_rejected_support).
    Because Ae depends only on Tb through the cached conductivity profile,
    each bound is exactly a Tb-support cut.
    """
    gravity_forward_model: Optional[str] = None
    """
    gravity_forward_model: When 'clairaut_hydrostatic', the observables
    'C20'/'C22' (and 'J2') are COMPUTED per sample — fluid Love number k_f
    by Clairaut integration over the sample's composite density profile
    (the same assembled core + mass-conserved silicate + cached hydrosphere
    the CMR2 recompute uses), then the hydrostatic pair C22_h = k_f q_r/4,
    C20_h = -3.324 C22_h (Tricarico rapid-rotation ratio), unnormalized at
    the GC21 1565 km reference radius — plus the sampled non-hydrostatic
    offsets dC20_nh/dC22_nh when present in param_space. When None
    (default; every pre-v4 config), 'J2'/'C22' keep their legacy meaning:
    read from the structure cache as-is. See
    plans/europa-clipper-v4-geodesy-plan.md and gravity_obs.py.
    """

    def __post_init__(self):
        """Validate configuration after initialization."""
        # Mode validation
        if self.mode not in ['mcmc', 'sbi']:
            raise ValueError(f"mode must be 'mcmc' or 'sbi', got '{self.mode}'")

        # Parameter space validation
        if not self.param_space:
            raise ValueError("param_space cannot be empty")

        # param_groups validation: members must not also appear in param_space
        all_members = []
        for group_key, members in self.param_groups.items():
            if group_key not in self.param_space:
                raise ValueError(
                    f"param_groups key '{group_key}' must appear in param_space"
                )
            for m in members:
                if m in self.param_space:
                    raise ValueError(
                        f"param_groups member '{m}' (group '{group_key}') must NOT "
                        f"also appear in param_space — it is expanded at runtime"
                    )
                all_members.append(m)

        # fixed_params must not overlap with param_space or group keys
        for fp_key in self.fixed_params:
            if fp_key in self.param_space:
                raise ValueError(
                    f"fixed_params key '{fp_key}' must NOT appear in param_space"
                )
            if fp_key in self.param_groups:
                raise ValueError(
                    f"fixed_params key '{fp_key}' must NOT be a param_groups group key"
                )

        # Observables validation
        if not self.observables:
            raise ValueError("observables cannot be empty")

        for obs_name, obs_spec in self.observables.items():
            value, uncertainty = obs_spec[0], obs_spec[1]
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
        # Legacy configs nest arrhenius_params inside sampler_settings.
        # Hoist it to the top-level field so the hash and dataclass field agree.
        data = dict(data)
        ss = data.get('sampler_settings')
        if isinstance(ss, dict) and 'arrhenius_params' in ss and not data.get('arrhenius_params'):
            data['sampler_settings'] = dict(ss)
            data['arrhenius_params'] = data['sampler_settings'].pop('arrhenius_params')
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
            'random_state': self.random_state,
            'arrhenius_params': sorted((self.arrhenius_params or {}).items()),
        }
        # Included ONLY when set so every pre-existing config hash (artifact
        # provenance, cache reuse) is unchanged by the field's introduction.
        if self.induction_bounds:
            hash_data['induction_bounds'] = sorted(self.induction_bounds.items())
        # Same include-only-when-set discipline (v4 geodesy, 2026-07-19).
        if self.gravity_forward_model:
            hash_data['gravity_forward_model'] = self.gravity_forward_model
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
    # Complex h2 per evaluated sample (n_eval, 2) = [Re_h2, Im_h2 signed];
    # None on results predating 2026-07-18.
    h2_results: Optional[np.ndarray] = None
    cmr2_results: Optional[np.ndarray] = None
    D_ocean_results: Optional[np.ndarray] = None
    D_iceIh_results: Optional[np.ndarray] = None
    # Total hydrosphere thickness (km) per sample: ice Ih + ocean + any HP
    # ices, from the cache's D_hsphere_km. Lets the GUI derive the rock
    # mantle thickness exactly (R_body - D_hsphere - R_core) even for
    # bodies with HP-ice layers. None on results predating 2026-07.
    D_hsphere_results: Optional[np.ndarray] = None
    # Aggregate high-pressure-ice (III + V + VI) shell thickness (km) per
    # sample, from the cache's D_iceIII/V/VI_km. 0.0 for bodies with no HP
    # ices (e.g. Europa seawater, basal-ocean P < 200 MPa). Lets the globe
    # draw the HP-ice shell between the ocean and the silicate floor. None
    # on results predating 2026-07-20.
    D_iceHP_results: Optional[np.ndarray] = None
    # Per-phase high-pressure-ice thicknesses (km), per sample, as a dict
    # {'III': (n,), 'V': (n,), 'VI': (n,)} from the cache's D_iceIII/V/VI_km.
    # Lets the globe draw each ice phase as its own labeled, proportional
    # shell (Ih is the textured surface; ocean uses D_ocean_results). None
    # on results predating 2026-07-20.
    D_icePhase_results: Optional[Dict[str, np.ndarray]] = None
    # Model-predicted unnormalized C20/C22 per sample (v4 geodesy configs
    # with gravity_forward_model='clairaut_hydrostatic'; includes the
    # sampled non-hydrostatic offsets). None otherwise / on older pkls.
    c20_results: Optional[np.ndarray] = None
    c22_results: Optional[np.ndarray] = None
    # Hydrostatic-REFERENCE C/MR² per sample: the C/MR² the model's
    # (offset-included) C22 would imply IF the body were hydrostatic, via
    # the Radau–Darwin inverse (gravity_obs.cmr2_from_c22_rd), computed with
    # the SAME omega/R_body/M/ref-radius the forward gravity map used. Pairs
    # with cmr2_results (the ACTUAL C/MR² from the structure moment integral):
    # the gap between the two is the inferred non-hydrostaticity. RD carries a
    # ~1% k_f systematic (≈0.0035 in C/MR²) — a floor the display must annotate
    # rather than read as physical. None unless the gravity forward model is
    # active / on pkls predating 2026-07-22.
    cmr2_hydro_results: Optional[np.ndarray] = None
    heating_results: Optional[List[Dict[str, float]]] = None
    convergence_metrics: Dict[str, Any] = field(default_factory=dict)
    metadata: Dict[str, Any] = field(default_factory=dict)
    # Importance weights from pocoMC for each posterior sample.
    # Shape: (n_samples,) with values normalised so sum ≈ 1.
    # None for older pkls loaded without this field.
    weights: Optional[np.ndarray] = None

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

        Uses ``nanargmax`` because amortized SBI results may NaN-pad the
        log-likelihood outside the recomputed subset (see
        ``SBIRunner._condition_and_package`` ``n_derived``); a plain
        ``argmax`` would return the first NaN-padded row. Falls back to the
        posterior median when every likelihood is NaN.

        Returns:
            1D array of parameter values
        """
        ll = np.asarray(self.log_likelihoods, dtype=np.float64)
        if not np.any(np.isfinite(ll)):
            return np.median(self.samples, axis=0)
        best_idx = int(np.nanargmax(ll))
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

        if prior_type not in ['uniform', 'normal', 'log-uniform',
                               'truncated_gaussian']:
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

        elif prior_type == 'truncated_gaussian':
            # v5 ice-thickness prior: needs mean, std, and [low, high] bounds.
            if spec.get('mean') is None or spec.get('std') is None:
                errors.append(f"Parameter '{param}': missing mean or std for "
                              "truncated_gaussian prior")
            elif spec.get('std', 1) <= 0:
                errors.append(f"Parameter '{param}': std must be positive")
            bounds = spec.get('bounds')
            if bounds is None or len(bounds) != 2:
                errors.append(f"Parameter '{param}': truncated_gaussian requires "
                              "[low, high] bounds")
            elif bounds[0] >= bounds[1]:
                errors.append(f"Parameter '{param}': bounds min >= max")

    # Check observables
    if not config.observables:
        errors.append("Observables dictionary is empty.")

    for obs, obs_spec in config.observables.items():
        value, uncertainty = obs_spec[0], obs_spec[1]
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
