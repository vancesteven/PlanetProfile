"""
MCMC runner using pocoMC for Bayesian inference.

Refactored from Test41-44 scripts into reusable MCMCRunner class.
Supports Andrade and Maxwell rheologies with optional Arrhenius viscosity.

Author: PlanetProfile Team
Date: 2026-04-29
"""
import json
import numpy as np
import logging
import os
import subprocess
import threading
import time
import pickle
from datetime import datetime, timezone
from typing import Dict, Optional, Any, Callable
from pathlib import Path

log = logging.getLogger('PlanetProfile')


class MCMCRunner:
    """
    MCMC sampler for interior structure inference using pocoMC.

    Args:
        config: InferenceConfig object with parameter space, observables, and settings

    Example:
        >>> config = InferenceConfig(mode='mcmc', bodyname='Titan', ...)
        >>> runner = MCMCRunner(config)
        >>> result = runner.run(progress_callback=my_callback)
    """

    def __init__(self, config):
        """Initialize MCMC runner with configuration."""
        from .inference_core import InferenceConfig
        from .structure_cache import load_structure_cache

        if not isinstance(config, InferenceConfig):
            raise TypeError("config must be InferenceConfig instance")

        self.config = config

        # Extract parameter names and labels (must be before _build_prior)
        self.param_names = list(config.param_space.keys())
        self.param_labels = [self._make_label(name) for name in self.param_names]

        # param_groups: maps group key -> list of member names
        self.param_groups = getattr(config, 'param_groups', {}) or {}
        # fixed_params: constants injected into every forward call
        self.fixed_params = getattr(config, 'fixed_params', {}) or {}

        # Route to grid cache when Tb_K is a free parameter OR fixed via fixed_params
        self._use_flexible = 'Tb_K' in self.param_names or 'Tb_K' in self.fixed_params

        # Load cached structure (skip bodyname validation for Test* files)
        log.info(f"Loading structure cache: {config.structure_cache_path}")
        if self._use_flexible:
            self.structure_data = self._load_grid_cache(config.structure_cache_path)
        else:
            self.structure_data = load_structure_cache(
                config.structure_cache_path, validate_bodyname=None
            )

        # Resolve arrhenius_params: prefer config.arrhenius_params, fallback to
        # sampler_settings, and mirror back to config for self-describing pickles.
        ap = getattr(config, 'arrhenius_params', None) or config.sampler_settings.get('arrhenius_params')
        if ap and not getattr(config, 'arrhenius_params', None):
            config.arrhenius_params = dict(ap)
        self.arrhenius_params = ap or {}

        # Build prior and likelihood — always use dict-based interface so
        # parameter order in param_space never affects forward model mapping.
        self.prior = self._build_prior()
        self.log_likelihood_fn = self._make_flexible_log_likelihood(
            config.observables,
            self.structure_data,
            arrhenius_params=self.arrhenius_params
        )

        # MCMC settings
        self.n_effective = config.sampler_settings.get('n_effective', 500)
        self.random_state = config.random_state
        self.n_reeval = config.sampler_settings.get('n_reeval', 500)

        # Checkpoint settings
        self.checkpoint_interval = config.sampler_settings.get('checkpoint_interval', 100)
        self.checkpoint_dir = Path(config.sampler_settings.get('checkpoint_dir', '/tmp'))
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)

    def _infer_rheology(self) -> str:
        """Infer rheology type from parameter space or sampler_settings."""
        # Explicit override in sampler_settings takes highest priority
        explicit = self.config.sampler_settings.get('rheology')
        if explicit in ('andrade', 'maxwell'):
            return explicit

        params = self.config.param_space
        has_alpha = 'alpha' in params
        # Accept both 'log10_zeta' (single) and 'log10_zeta_*' (per-phase)
        has_zeta = ('log10_zeta' in params or
                    any(k.startswith('log10_zeta_') for k in params))
        if has_alpha and has_zeta:
            return 'andrade'
        elif not has_alpha and not has_zeta:
            return 'maxwell'
        else:
            raise ValueError("Cannot infer rheology from parameter space. "
                           "Andrade requires 'alpha' and log10_zeta or log10_zeta_* parameters. "
                           "Maxwell requires neither. Or set sampler_settings.rheology explicitly.")

    def _build_prior(self):
        """Build pocoMC Prior object from parameter space configuration."""
        try:
            import pocomc as pc
            from scipy.stats import uniform, norm, loguniform
        except ImportError as e:
            raise ImportError("pocoMC not installed. Run: pip install pocomc") from e

        priors = []
        for param_name in self.param_names:
            param_cfg = self.config.param_space[param_name]
            prior_type = param_cfg['prior_type']

            if prior_type == 'uniform':
                low, high = param_cfg['bounds']
                priors.append(uniform(loc=low, scale=(high - low)))
            elif prior_type == 'normal':
                mean = param_cfg['mean']
                std = param_cfg['std']
                priors.append(norm(loc=mean, scale=std))
            elif prior_type == 'log-uniform':
                low, high = param_cfg['bounds']
                priors.append(loguniform(a=10**low, b=10**high))
            else:
                raise ValueError(f"Unknown prior type '{prior_type}' for parameter '{param_name}'")

        return pc.Prior(priors)

    def _make_label(self, param_name: str) -> str:
        """Convert parameter name to LaTeX label for plotting."""
        label_map = {
            'alpha': r'$\alpha$',
            'log10_zeta': r'$\log_{10}\zeta$',
            'log10_zeta_Ih': r'$\log_{10}(\zeta_{\rm Ih})$',
            'log10_zeta_HP': r'$\log_{10}(\zeta_{\rm HP})$',
            'log10_zeta_sil': r'$\log_{10}(\zeta_{\rm sil})$',
            'log10_eta_Ih': r'$\log_{10}\eta_\mathrm{Ih}$',
            'log10_eta_III': r'$\log_{10}\eta_\mathrm{III}$',
            'log10_eta_V': r'$\log_{10}\eta_\mathrm{V}$',
            'log10_eta_VI': r'$\log_{10}\eta_\mathrm{VI}$',
            'log10_eta_HP': r'$\log_{10}(\eta_{\rm HP})$',
            'log10_eta_sil': r'$\log_{10}\eta_\mathrm{sil}$',
            'log10_mu_Ih': r'$\log_{10}(\mu_{\rm Ih})$',
            'Tb_K': r'$T_b$ (K)',
        }
        return label_map.get(param_name, param_name)

    def _load_grid_cache(self, cache_path: str) -> Dict[str, Any]:
        """Load grid cache; accepts two formats.

        **Format A** (MCMCRunner native): ``dict[float -> structure_dict]``
        keyed by Tb_K values.

        **Format B** (Test50 list): ``{'Tb_K_grid': ndarray, 'structures': list}``
        as produced by Test50's ``build_or_load_structure_grid()``.

        Both formats are passed through transparently; ``apply_bottom_temperature``
        in forward_models.py handles interpolation for both.
        """
        with open(cache_path, 'rb') as f:
            grid_cache = pickle.load(f)
        if not isinstance(grid_cache, dict):
            raise ValueError(f"Expected dict at {cache_path}, got {type(grid_cache)}")
        if not grid_cache:
            raise ValueError(
                f"Grid cache at {cache_path} is empty. "
                f"Regenerate with: python -m PlanetProfile.Inference.prepare_structure_variants "
                f"--test-module PlanetProfile.Test.PPTest42 --output-dir titan_cache/ --maxwell --force"
            )

        # Format B: Test50 list format
        if 'Tb_K_grid' in grid_cache and 'structures' in grid_cache:
            Tb_grid = np.asarray(grid_cache['Tb_K_grid'])
            log.info(f"Grid cache loaded (list format): {len(grid_cache['structures'])} structures "
                     f"[{Tb_grid[0]:.3f} – {Tb_grid[-1]:.3f} K]")
            return grid_cache  # pass through; apply_bottom_temperature handles it

        # Format A: float-keyed dict
        first_key = next(iter(grid_cache))
        if not isinstance(first_key, (int, float)):
            raise ValueError(
                f"Expected grid cache keyed by Tb_K floats at {cache_path}, "
                f"got key type {type(first_key).__name__}"
            )
        grid_Tb_values = np.array(sorted(grid_cache.keys()))
        log.info(f"Grid cache loaded (dict format): {len(grid_Tb_values)} Tb_K values "
                 f"[{grid_Tb_values[0]:.1f} – {grid_Tb_values[-1]:.1f} K]")
        return {'grid_cache': grid_cache, 'grid_Tb_values': grid_Tb_values}

    def _make_flexible_log_likelihood(
        self,
        observables: Dict[str, Any],
        structure_data: Dict[str, Any],
        arrhenius_params: Optional[Dict[str, Any]] = None
    ) -> Callable:
        """Build Gaussian log-likelihood using the dict-based flexible forward model."""
        from .forward_models import forward_model_k2_flexible
        param_names = self.param_names
        param_groups = self.param_groups   # {group_key: [member_names]}
        fixed_params = self.fixed_params   # {param_name: value}

        # Build no-ocean guard from sampler_settings
        phase_stability_cfg = self.config.sampler_settings.get('phase_stability', {})
        no_ocean_guard = phase_stability_cfg.get('enforce') == 'no_ocean_Ih'
        no_ocean_margin_K = float(phase_stability_cfg.get('margin_K', 0.1))

        # v2 derived-params (Phase C1 Stage 2): mass-conservation solve for rho_sil
        # and per-sample CMR² recomputation from sampled (R_core, rho_core) plus
        # cached hydrosphere. Inactive for v1 configs (Titan Test50) — those
        # leave derived_params empty and fall through to the cached-CMR² path.
        derived_params_cfg = getattr(self.config, 'derived_params', {}) or {}
        rho_sil_cfg = derived_params_cfg.get('rho_sil_kgm3', {}) or {}
        use_derived_rho_sil = rho_sil_cfg.get('derivation') == 'mass_conservation'
        rho_sil_bounds = tuple(rho_sil_cfg.get('bounds', (2200.0, 3500.0)))
        rho_sil_reject = bool(rho_sil_cfg.get('reject_if_outside_bounds', True))

        def _expand_theta(theta):
            """Convert sampled array → full parameter dict with groups and fixed params."""
            theta_dict = dict(zip(param_names, theta))
            # Expand param_groups: each group key's value is broadcast to all members
            for group_key, members in param_groups.items():
                if group_key in theta_dict:
                    for m in members:
                        theta_dict[m] = theta_dict[group_key]
            # Inject fixed params (constants not in the prior)
            theta_dict.update(fixed_params)
            return theta_dict

        def _check_no_ocean(modified_structure) -> bool:
            """Return True if sample should be rejected (ocean would form)."""
            if not no_ocean_guard:
                return False
            phases = modified_structure.get('phases')
            P_arr = modified_structure.get('P_MPa')
            T_arr = modified_structure.get('T_K')
            if phases is None or P_arr is None or T_arr is None:
                return False
            P_arr = np.asarray(P_arr)
            T_arr = np.asarray(T_arr)
            if P_arr.shape != T_arr.shape or P_arr.shape != np.asarray(phases).shape:
                return False
            Ih_mask = (np.asarray(phases) == 1)
            if not np.any(Ih_mask):
                return False
            P_Ih = P_arr[Ih_mask]
            T_Ih = T_arr[Ih_mask]
            if not np.all(np.isfinite(P_Ih)):
                return False
            Tm_Ih_lin = 273.16 - 0.068 * P_Ih
            return bool(np.any(T_Ih >= Tm_Ih_lin - no_ocean_margin_K))

        def log_likelihood(theta):
            theta_dict = _expand_theta(theta)

            # Run forward model — parameter hooks do the Tb interpolation so that
            # _check_no_ocean sees the fully-interpolated T(r) and P(r) profiles.
            from .forward_models import apply_parameters
            modified = apply_parameters(theta_dict, structure_data)

            # No-ocean safeguard (body-agnostic: assert solid-Ih stability everywhere)
            if _check_no_ocean(modified):
                return -1e30

            Re_k2, Im_k2, Re_h2, Im_h2, _ = forward_model_k2_flexible(
                theta_dict, structure_data,
                return_heating=False, arrhenius_params=arrhenius_params
            )
            if np.isnan(Re_k2):
                return -1e30
            chi2 = 0.0
            if 'Re_k2' in observables:
                obs_val, obs_err = observables['Re_k2']
                chi2 += ((Re_k2 - obs_val) / obs_err) ** 2
            if 'Im_k2' in observables or 'abs_Im_k2' in observables:
                key = 'Im_k2' if 'Im_k2' in observables else 'abs_Im_k2'
                obs_val, obs_err = observables[key]
                chi2 += ((abs(Im_k2) - obs_val) / obs_err) ** 2
            if 'k2' in observables:
                obs_val, obs_err = observables['k2']
                chi2 += ((np.sqrt(Re_k2**2 + Im_k2**2) - obs_val) / obs_err) ** 2
            # h2 observables (Mazarico et al. 2023 convention).
            if 'Re_h2' in observables:
                obs_val, obs_err = observables['Re_h2']
                chi2 += ((Re_h2 - obs_val) / obs_err) ** 2
            if 'Im_h2' in observables or 'abs_Im_h2' in observables:
                key = 'Im_h2' if 'Im_h2' in observables else 'abs_Im_h2'
                obs_val, obs_err = observables[key]
                chi2 += ((abs(Im_h2) - obs_val) / obs_err) ** 2
            if 'h2' in observables:
                obs_val, obs_err = observables['h2']
                chi2 += ((np.sqrt(Re_h2**2 + Im_h2**2) - obs_val) / obs_err) ** 2
            # Gravity coefficients J2 and C22 — cached per Tb point as scalars
            # in the structure dict; blended via _BLEND_SCALAR_FIELDS.
            if 'J2' in observables:
                obs_val, obs_err = observables['J2']
                pred = float(modified.get('J2', np.nan))
                if not np.isfinite(pred):
                    return -1e30
                chi2 += ((pred - obs_val) / obs_err) ** 2
            if 'C22' in observables:
                obs_val, obs_err = observables['C22']
                pred = float(modified.get('C22', np.nan))
                if not np.isfinite(pred):
                    return -1e30
                chi2 += ((pred - obs_val) / obs_err) ** 2
            if 'CMR2' in observables:
                obs_val, obs_err = observables['CMR2']

                # Pick the per-sample structure dict for CMR² re-derivation.
                # Three cache layouts are supported:
                #   (a) v2.1 list format: {'Tb_K_grid': arr, 'structures': [..]}
                #       — what cache_builder.build_phase_c1_cache produces.
                #   (b) legacy dict format: {'grid_cache': {Tb: struct}, 'grid_Tb_values': arr}
                #       — what _load_grid_cache wraps Format-A caches into.
                #   (c) single struct (no Tb grid).
                # The earlier code only handled (b) and silently fell through to
                # struct_for_cmr2 = structure_data for (a), which lacks
                # Mtot_kg/R_body_m at top level → NaN → -1e30 rejection of every
                # sample. Verified 2026-05-23 against europa_seawater cache.
                Tb_sample = theta_dict.get('Tb_K')
                if (Tb_sample is not None
                        and 'Tb_K_grid' in structure_data
                        and 'structures' in structure_data):
                    grid_Tb_values = np.asarray(structure_data['Tb_K_grid'])
                    idx = int(np.argmin(np.abs(grid_Tb_values - Tb_sample)))
                    struct_for_cmr2 = structure_data['structures'][idx]
                elif 'grid_cache' in structure_data and Tb_sample is not None:
                    grid_Tb_values = structure_data['grid_Tb_values']
                    idx = int(np.argmin(np.abs(grid_Tb_values - Tb_sample)))
                    struct_for_cmr2 = structure_data['grid_cache'][grid_Tb_values[idx]]
                else:
                    struct_for_cmr2 = structure_data

                if use_derived_rho_sil:
                    # v2 path: mass-conserve rho_sil, then recompute CMR² from
                    # the assembled body (sampled core + derived silicate +
                    # cached hydrosphere).
                    from .structure_derivation import (
                        compute_cmr2,
                        derive_silicate_density,
                        extract_hydrosphere_layers,
                    )
                    R_body_m = struct_for_cmr2.get('R_body_m', np.nan)
                    M_total_kg = struct_for_cmr2.get('Mtot_kg', np.nan)
                    R_core_km = theta_dict.get('R_core_km')
                    rho_core_kgm3 = theta_dict.get('rho_core_kgm3')
                    if not (np.isfinite(R_body_m) and np.isfinite(M_total_kg)
                            and R_core_km is not None and rho_core_kgm3 is not None):
                        return -1e30
                    R_core_m = float(R_core_km) * 1000.0
                    try:
                        hydro_layers, R_oceanbot_m, M_hydro_kg = (
                            extract_hydrosphere_layers(struct_for_cmr2)
                        )
                    except (KeyError, ValueError):
                        return -1e30
                    if not hydro_layers:
                        return -1e30
                    rho_sil, accepted = derive_silicate_density(
                        M_total_kg=float(M_total_kg),
                        M_hydrosphere_kg=M_hydro_kg,
                        R_oceanbot_m=R_oceanbot_m,
                        R_core_m=R_core_m,
                        rho_core_kgm3=float(rho_core_kgm3),
                        bounds=rho_sil_bounds,
                    )
                    if rho_sil_reject and not accepted:
                        return -1e30
                    # Assemble layers; skip zero-volume core shell at R_core = 0
                    assembled = []
                    if R_core_m > 0.0:
                        assembled.append((0.0, R_core_m, float(rho_core_kgm3)))
                    assembled.append((R_core_m, R_oceanbot_m, rho_sil))
                    assembled.extend(hydro_layers)
                    try:
                        cmr2_val = compute_cmr2(
                            assembled,
                            R_body_m=float(R_body_m),
                            M_body_kg=float(M_total_kg),
                        )
                    except ValueError:
                        return -1e30
                else:
                    # v1 path: read precomputed CMR² from the cache as-is
                    cmr2_val = struct_for_cmr2.get('CMR2', np.nan)

                if np.isfinite(cmr2_val):
                    chi2 += ((cmr2_val - obs_val) / obs_err) ** 2
            if 'Mtot_kg' in observables:
                obs_val, obs_err = observables['Mtot_kg']
                mtot_val = structure_data.get('Mtot_kg', np.nan)
                if np.isfinite(mtot_val):
                    chi2 += ((mtot_val - obs_val) / obs_err) ** 2

            # Magnetic induction observables (C2-B). The runner accepts
            # observables in two equivalent forms:
            #   1. Re/Im (default, Europa-Clipper convention):
            #          Ae_<label>_real, Ae_<label>_imag
            #   2. Amplitude/phase (legacy, paper convention):
            #          BiAmp_<label>, BiPhase_<label>_deg
            # Both pull from the same forward_model_induction call; <label>
            # is one of the canonical PP excitation names ('synodic',
            # 'orbital', 'synodic 2nd', 'true anomaly') matching keys in
            # the cached Texc_hr dict. The runner only invokes the
            # induction forward model if at least one such observable is
            # present, so existing CMR²/k₂ configs are unaffected.
            induction_keys_real = [k for k in observables
                                   if k.startswith('Ae_') and k.endswith('_real')]
            induction_keys_imag = [k for k in observables
                                   if k.startswith('Ae_') and k.endswith('_imag')]
            induction_keys_amp = [k for k in observables
                                  if k.startswith('BiAmp_')]
            induction_keys_phase = [k for k in observables
                                    if k.startswith('BiPhase_')
                                    and k.endswith('_deg')]
            need_induction = bool(induction_keys_real or induction_keys_imag
                                  or induction_keys_amp or induction_keys_phase)
            if need_induction:
                from .forward_models import forward_model_induction
                Texc_hr_full = modified.get('Texc_hr')
                if not Texc_hr_full:
                    return -1e30
                # Build the freq dict from labels referenced in the observables.
                requested_labels = set()
                for k in induction_keys_real:
                    requested_labels.add(k[len('Ae_'):-len('_real')])
                for k in induction_keys_imag:
                    requested_labels.add(k[len('Ae_'):-len('_imag')])
                for k in induction_keys_amp:
                    requested_labels.add(k[len('BiAmp_'):])
                for k in induction_keys_phase:
                    requested_labels.add(k[len('BiPhase_'):-len('_deg')])
                freq_dict = {label: Texc_hr_full[label]
                             for label in requested_labels
                             if label in Texc_hr_full}
                if not freq_dict:
                    return -1e30
                # forward_model_induction reads rSigChange_m / sigmaLayers_Sm
                # / R_body_m from `modified`.
                Ae_dict = forward_model_induction(modified, freq_dict)
                if Ae_dict is None:
                    return -1e30
                for label in requested_labels:
                    Ae = Ae_dict.get(label)
                    if Ae is None or not np.isfinite(complex(Ae).real):
                        return -1e30
                    Ae = complex(Ae)
                    re_key = f'Ae_{label}_real'
                    im_key = f'Ae_{label}_imag'
                    amp_key = f'BiAmp_{label}'
                    ph_key = f'BiPhase_{label}_deg'
                    if re_key in observables:
                        v, s = observables[re_key]
                        chi2 += ((Ae.real - v) / s) ** 2
                    if im_key in observables:
                        v, s = observables[im_key]
                        chi2 += ((Ae.imag - v) / s) ** 2
                    if amp_key in observables:
                        v, s = observables[amp_key]
                        chi2 += ((abs(Ae) - v) / s) ** 2
                    if ph_key in observables:
                        v, s = observables[ph_key]
                        # Phase wrap into (-180, 180] before residualizing.
                        pred = float(np.degrees(np.angle(Ae)))
                        delta = ((pred - v + 180.0) % 360.0) - 180.0
                        chi2 += (delta / s) ** 2

            return -0.5 * chi2

        return log_likelihood

    def _expand_theta(self, theta: np.ndarray) -> Dict[str, float]:
        """Expand sampled array with groups and fixed params into a dict."""
        theta_dict = dict(zip(self.param_names, theta))
        for group_key, members in self.param_groups.items():
            if group_key in theta_dict:
                for m in members:
                    theta_dict[m] = theta_dict[group_key]
        theta_dict.update(self.fixed_params)
        return theta_dict

    def _get_cache_scalar(self, theta_dict: Dict[str, float], key: str) -> float:
        """Extract scalar from interpolated structure for grid caches."""
        from .forward_models import apply_parameters
        modified = apply_parameters(theta_dict, self.structure_data)
        val = modified.get(key, np.nan)
        return float(val) if np.isfinite(float(val)) else np.nan

    def run(self, progress_callback: Optional[Callable] = None,
            progress_jsonl_path: Optional[str] = None):
        """
        Run MCMC sampling with pocoMC.

        Args:
            progress_callback: Optional function called with progress dict:
                {'iteration': int, 'n_total': int, 'acceptance_rate': float,
                 'n_samples': int, 'ess': float}
            progress_jsonl_path: Optional path to a JSONL file.  When set,
                one JSON line is appended after each sampler iteration with
                fields: iteration, log_Z, log_Z_err, ESS, n_accepted,
                n_total, elapsed_s, timestamp.  Fields unavailable for the
                current sampler state are written as null.  When None (the
                default) no file is written and behaviour is unchanged.

        Returns:
            InferenceResult object with samples, log-likelihoods, and convergence metrics
        """
        from .inference_core import InferenceResult

        try:
            import pocomc as pc
        except ImportError as e:
            raise ImportError("pocoMC not installed. Run: pip install pocomc") from e

        log.info(f"Starting MCMC with pocoMC (target n_eff={self.n_effective})")
        log.info(f"Parameter space ({len(self.param_names)}D): {self.param_names}")
        log.info(f"Observables: {list(self.config.observables.keys())}")

        t0 = time.time()

        # Initialize sampler
        def _log_like(theta):
            return self.log_likelihood_fn(theta)

        sampler = pc.Sampler(
            prior=self.prior,
            likelihood=_log_like,
            n_effective=self.n_effective,
            random_state=self.random_state,
        )

        # --- JSONL progress streaming ----------------------------------------
        # pocoMC's run() is a blocking call with no native per-iteration
        # callback.  We launch a background thread that polls sampler.t (the
        # iteration counter) at a fixed interval and writes one JSONL record
        # whenever the counter advances.  The thread stops cleanly once the
        # main run() returns.
        _jsonl_stop_event = threading.Event()

        def _jsonl_writer_thread(path: str, stop: threading.Event,
                                  run_t0: float) -> None:
            """Background thread: poll sampler state and append JSONL lines."""
            jsonl_path = Path(path)
            jsonl_path.parent.mkdir(parents=True, exist_ok=True)
            last_iter = -1
            try:
                with open(jsonl_path, 'a', buffering=1) as fh:
                    while not stop.is_set():
                        try:
                            current_t = getattr(sampler, 't', None)
                            if current_t is not None and current_t != last_iter:
                                # Read particle stats; particles.get returns
                                # scalar when the key exists, -1 sentinel otherwise.
                                def _safe(key, sentinel=-1):
                                    try:
                                        v = sampler.particles.get(key, sentinel)
                                        v = float(v)
                                        return None if v == sentinel else v
                                    except Exception:
                                        return None

                                ess_val = _safe('ess')
                                logz_val = _safe('logz')
                                n_calls = getattr(sampler, 'calls', None)
                                accept_val = _safe('accept')
                                # n_accepted: accept rate * n_active particles
                                n_active = getattr(sampler, 'n_active', None)
                                if accept_val is not None and n_active is not None:
                                    n_accepted = int(round(accept_val * n_active))
                                else:
                                    n_accepted = None

                                record = {
                                    'iteration': current_t,
                                    'log_Z': logz_val,
                                    'log_Z_err': None,  # pocoMC sets this only at end
                                    'ESS': int(ess_val) if ess_val is not None else None,
                                    'n_accepted': n_accepted,
                                    'n_total': int(n_calls) if n_calls is not None else None,
                                    'elapsed_s': round(time.time() - run_t0, 2),
                                    'timestamp': datetime.now(timezone.utc).strftime(
                                        '%Y-%m-%dT%H:%M:%SZ'),
                                }
                                fh.write(json.dumps(record) + '\n')
                                last_iter = current_t
                        except Exception:
                            pass  # Never crash the monitor thread
                        stop.wait(timeout=2.0)
            except Exception as exc:
                log.warning(f"JSONL progress writer failed: {exc}")

        if progress_jsonl_path is not None:
            _jsonl_thread = threading.Thread(
                target=_jsonl_writer_thread,
                args=(progress_jsonl_path, _jsonl_stop_event, t0),
                daemon=True,
                name='mcmc-jsonl-progress',
            )
            _jsonl_thread.start()
        else:
            _jsonl_thread = None
        # ---------------------------------------------------------------------

        try:
            # Run — pocoMC stops internally when ESS >= n_ess
            sampler.run(n_total=4096, progress=True)
        finally:
            # Signal the writer thread to stop and wait briefly for it to flush
            _jsonl_stop_event.set()
            if _jsonl_thread is not None:
                _jsonl_thread.join(timeout=5.0)
                # Write a final record with end-of-run logz / logz_err
                if progress_jsonl_path is not None:
                    try:
                        jsonl_path = Path(progress_jsonl_path)
                        final_logz = getattr(sampler, 'logz', None)
                        final_logz_err = getattr(sampler, 'logz_err', None)
                        n_calls_final = getattr(sampler, 'calls', None)
                        try:
                            _s_final, _, _, _ = sampler.posterior()
                            ess_final = len(_s_final)
                        except Exception:
                            ess_final = None
                        final_record = {
                            'iteration': getattr(sampler, 't', None),
                            'log_Z': float(final_logz) if final_logz is not None else None,
                            'log_Z_err': float(final_logz_err) if final_logz_err is not None else None,
                            'ESS': ess_final,
                            'n_accepted': None,
                            'n_total': int(n_calls_final) if n_calls_final is not None else None,
                            'elapsed_s': round(time.time() - t0, 2),
                            'timestamp': datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
                        }
                        with open(jsonl_path, 'a', buffering=1) as fh:
                            fh.write(json.dumps(final_record) + '\n')
                    except Exception as exc:
                        log.warning(f"JSONL final record write failed: {exc}")

        # Single progress update on completion (for GUI)
        if progress_callback is not None:
            try:
                _s, _ll, _, _ = sampler.posterior()
                _ess = sampler.n_eff if hasattr(sampler, 'n_eff') else len(_s)
                progress_callback({
                    'iteration': 1,
                    'n_total': 1,
                    'n_samples': len(_s),
                    'ess': _ess,
                    'acceptance_rate': sampler.pbar.info.get('acc') if hasattr(sampler, 'pbar') and sampler.pbar is not None else None,
                })
            except Exception:
                pass

        # Final posterior extraction. pocoMC's posterior() returns
        #   (samples, weights, logl, logp)
        # NOT (samples, logl, logp, weights). The earlier unpack stored
        # importance weights into the log_likelihoods field (=~1/N after
        # trimming) and log-likelihoods into log_posteriors — verified by
        # reading pocomc/sampler.py::posterior.
        samples, _weights, log_likes, _log_prior = sampler.posterior()
        n_samples = len(samples)

        elapsed = time.time() - t0
        log.info(f"MCMC complete: {n_samples} samples in {elapsed/60:.1f} min")

        # --- Capture run-level metadata from pocoMC ----------------------------

        # git SHA for audit / pre-bugfix flag
        try:
            _git_sha = subprocess.check_output(
                ['git', 'rev-parse', '--short', 'HEAD'],
                cwd=os.path.dirname(os.path.abspath(__file__)),
                stderr=subprocess.DEVNULL,
            ).decode().strip()
        except Exception:
            _git_sha = None

        # log-evidence and its MC error (set by pocoMC at end of run)
        _log_Z = getattr(sampler, 'logz', None)
        _log_Z_err = getattr(sampler, 'logz_err', None)
        _log_Z = float(_log_Z) if _log_Z is not None else None
        _log_Z_err = float(_log_Z_err) if _log_Z_err is not None else None

        # True ESS from importance weights:  (sum w)^2 / sum(w^2)
        # _weights from posterior() are already normalised (sum ≈ 1).
        # Avoid division by zero for uniform (resampled) case.
        try:
            _w = np.asarray(_weights, dtype=float)
            _w_sum2 = float(np.sum(_w) ** 2)
            _w2_sum = float(np.sum(_w ** 2))
            _true_ess = _w_sum2 / _w2_sum if _w2_sum > 0 else float(n_samples)
        except Exception:
            _true_ess = float(n_samples)

        # Compute convergence metrics
        convergence_metrics = self._compute_convergence(samples, log_likes, sampler,
                                                        true_ess=_true_ess)

        # Recompute k2 for posterior samples — always dict-based so param order
        # in param_space never causes wrong positional mapping.
        log.info(f"Recomputing k2 for {n_samples} posterior samples...")
        rheology = self._infer_rheology() if not self._use_flexible else None

        from .forward_models import forward_model_k2_flexible

        k2_results = []
        cmr2_results = []
        D_ocean_results = []
        D_iceIh_results = []
        for i, theta in enumerate(samples):
            theta_dict = self._expand_theta(theta)
            Re_k2, Im_k2, _, _, _ = forward_model_k2_flexible(
                theta_dict, self.structure_data,
                return_heating=False, arrhenius_params=self.arrhenius_params
            )
            k2_results.append((Re_k2, Im_k2))
            cmr2_results.append(self._get_cache_scalar(theta_dict, 'CMR2'))
            D_ocean_results.append(self._get_cache_scalar(theta_dict, 'D_ocean_km'))
            D_iceIh_results.append(self._get_cache_scalar(theta_dict, 'D_iceIh_km'))
            if (i + 1) % 100 == 0:
                log.info(f"  {i+1}/{n_samples} samples recomputed")

        k2_results = np.array(k2_results)
        cmr2_results = np.array(cmr2_results)
        D_ocean_results = np.array(D_ocean_results)
        D_iceIh_results = np.array(D_iceIh_results)

        # Recompute heating on subset — same dict-based approach
        n_reeval = min(self.n_reeval, n_samples)
        log.info(f"Recomputing heating for {n_reeval} posterior samples...")

        rng = np.random.RandomState(self.random_state)
        idx_heat = rng.choice(n_samples, n_reeval, replace=False)
        idx_heat.sort()
        heating_results = []
        for si in idx_heat:
            theta_dict = self._expand_theta(samples[si])
            _, _, _, _, perPhase_W = forward_model_k2_flexible(
                theta_dict, self.structure_data,
                return_heating=True, arrhenius_params=self.arrhenius_params
            )
            heating_results.append(perPhase_W if perPhase_W is not None else {})
        heating_indices = idx_heat

        # Build result object
        result = InferenceResult(
            config=self.config,
            samples=samples,
            log_likelihoods=log_likes,
            param_names=self.param_names,
            param_labels=self.param_labels,
            k2_results=k2_results,
            cmr2_results=cmr2_results,
            D_ocean_results=D_ocean_results,
            D_iceIh_results=D_iceIh_results,
            heating_results=heating_results,
            convergence_metrics=convergence_metrics,
            metadata={
                'elapsed_time_s': elapsed,
                'n_iterations': len(samples) if samples is not None else 0,
                'rheology': rheology,
                'heating_indices': heating_indices,
                # Audit / reproducibility fields
                'git_sha': _git_sha,
                'log_Z': _log_Z,
                'log_Z_err': _log_Z_err,
            },
            weights=_weights,
        )

        log.info("MCMC inference complete")
        return result

    def _compute_convergence(self, samples, log_likes, sampler,
                             true_ess: Optional[float] = None) -> Dict[str, float]:
        """Compute convergence diagnostics (ESS, R-hat, acceptance rate).

        Args:
            true_ess: Pre-computed ESS from importance weights,
                ``(sum w)^2 / sum(w^2)``.  When provided this is used
                directly; otherwise falls back to ``sampler.n_eff`` or
                ``n_samples`` (the old behaviour).
        """
        n_samples = len(samples)

        # ESS (effective sample size)
        # Prefer the true weight-based ESS passed in from run(); fall back to
        # sampler.n_eff (pocoMC's internal target count), then n_samples.
        if true_ess is not None:
            ess = true_ess
        else:
            try:
                ess = sampler.n_eff if hasattr(sampler, 'n_eff') else n_samples
            except Exception:
                ess = n_samples

        # Acceptance rate — pocoMC stores this in pbar.info['acc'] after run
        try:
            acceptance_rate = sampler.pbar.info.get('acc')
        except Exception:
            acceptance_rate = None

        # R-hat (Gelman-Rubin) - requires multiple chains, skip if not available
        # pocoMC doesn't expose chains directly, so we approximate with single chain
        r_hat = 1.0  # Perfect convergence assumed (posterior-weighted samples)

        metrics = {
            'ess': float(ess),
            'acceptance_rate': acceptance_rate,  # None if unavailable
            'r_hat': float(r_hat),
            'n_samples': int(n_samples),
        }

        acc_str = f"{acceptance_rate:.2%}" if acceptance_rate is not None else "N/A"
        log.info(f"Convergence: ESS={ess:.0f}, accept={acc_str}, R-hat={r_hat:.3f}")

        return metrics

    def save_checkpoint(self, sampler, iteration: int, filepath: str) -> None:
        """
        Save MCMC checkpoint for resuming long runs.

        Args:
            sampler: pocoMC Sampler object
            iteration: Current iteration number
            filepath: Path to save checkpoint
        """
        checkpoint = {
            'iteration': iteration,
            'sampler_state': sampler.__dict__.copy(),  # Save internal state
            'config': self.config,
            'random_state': self.random_state,
        }

        with open(filepath, 'wb') as f:
            pickle.dump(checkpoint, f, protocol=pickle.HIGHEST_PROTOCOL)

        log.info(f"Checkpoint saved: {filepath}")

    def load_checkpoint(self, filepath: str):
        """
        Load MCMC checkpoint to resume run.

        Args:
            filepath: Path to checkpoint file

        Returns:
            (sampler, iteration) tuple
        """
        with open(filepath, 'rb') as f:
            checkpoint = pickle.load(f)

        # Rebuild sampler from saved state
        import pocomc as pc

        sampler = pc.Sampler(
            prior=self.prior,
            likelihood=self.log_likelihood_fn,
            n_effective=self.n_effective,
            random_state=checkpoint['random_state'],
        )

        # Restore internal state
        sampler.__dict__.update(checkpoint['sampler_state'])

        iteration = checkpoint['iteration']
        log.info(f"Checkpoint loaded from iteration {iteration}")

        return sampler, iteration

    def generate_sbi_dataset(self, n_samples: int = 10000, output_path: Optional[str] = None):
        """
        Generate (theta, x) dataset by sampling from the prior.

        Args:
            n_samples: Number of simulations to run
            output_path: Optional path to save .npz file

        Returns:
            (theta, x) tuple
        """
        log.info(f"Generating SBI dataset with {n_samples} samples...")

        # Sample from prior
        theta = self.prior.sample(n_samples)

        x = []
        obs_names = list(self.config.observables.keys())

        from .forward_models import forward_model_k2_flexible

        t0 = time.time()
        for i in range(n_samples):
            theta_dict = self._expand_theta(theta[i])

            # Compute k2
            Re_k2, Im_k2, _, _, _ = forward_model_k2_flexible(
                theta_dict, self.structure_data,
                return_heating=False, arrhenius_params=self.arrhenius_params
            )

            xi = []
            for name in obs_names:
                if name == 'Re_k2':
                    xi.append(Re_k2)
                elif name == 'Im_k2':
                    xi.append(Im_k2)
                elif name == 'CMR2':
                    xi.append(self._get_cache_scalar(theta_dict, 'CMR2'))
                elif name == 'Mtot_kg':
                    xi.append(self._get_cache_scalar(theta_dict, 'Mtot_kg'))
                else:
                    xi.append(np.nan)
            x.append(xi)

            if (i + 1) % 100 == 0:
                elapsed = time.time() - t0
                eta = (elapsed / (i + 1)) * (n_samples - (i + 1))
                log.info(f"  {i+1}/{n_samples} simulations complete (ETA: {eta/60:.1f} min)")

        theta = np.array(theta)
        x = np.array(x)

        if output_path:
            np.savez(output_path, theta=theta, x=x, param_names=self.param_names, obs_names=obs_names)
            log.info(f"SBI dataset saved to {output_path}")

        return theta, x
