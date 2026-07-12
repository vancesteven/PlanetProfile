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

        # Lazily-loaded CMR2 discretization-offset anchor sidecar (Test52
        # Phase 2, plan decision 1). None until _load_cmr2_offset_sidecar()
        # is first called; _cmr2_offset_sidecar_checked distinguishes "not
        # yet checked" from "checked, no sidecar file exists" so the lookup
        # (and its one-time log line) only ever happens once per runner
        # instance.
        self._cmr2_offset_sidecar = None
        self._cmr2_offset_sidecar_checked = False
        # Diagnostic side-channel set by _derive_cmr2_via_mass_conservation
        # on every call: None on success, else a short reason string
        # ('missing_inputs', 'hydrosphere_error', 'hydrosphere_empty',
        # 'rho_sil_bounds', 'density_inversion', 'cmr2_valueerror'). Lets
        # generate_sbi_dataset's support-guard counting distinguish a
        # density-inversion rejection from other reject-causes without
        # duplicating the derivation logic or changing this method's
        # public None/float return contract.
        self._last_cmr2_reject_reason: Optional[str] = None

        # Route to grid cache when Tb_K is a free parameter OR fixed via fixed_params
        self._use_flexible = 'Tb_K' in self.param_names or 'Tb_K' in self.fixed_params

        # Load cached structure (skip bodyname validation for Test* files).
        # Paths in configs are repo-relative; resolve against the repo root
        # when the CWD is elsewhere (e.g. the Streamlit app).
        from .structure_cache import resolve_cache_path
        self._resolved_cache_path = resolve_cache_path(config.structure_cache_path)
        log.info(f"Loading structure cache: {self._resolved_cache_path}")
        if self._use_flexible:
            self.structure_data = self._load_grid_cache(str(self._resolved_cache_path))
        else:
            self.structure_data = load_structure_cache(
                str(self._resolved_cache_path), validate_bodyname=None
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

        # Pre-compute induction Ae for every grid point at init time so the
        # likelihood can look up by index instead of calling AeResponse
        # (mpmath, ~1.5 s/call) on every sample.  Only done when induction
        # observables are present and a grid cache is loaded.
        self._ae_grid_cache: Dict[int, Optional[Dict[str, complex]]] = {}
        self._precompute_ae_grid(config.observables)

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

    def _precompute_ae_grid(self, observables: Dict[str, Any]) -> None:
        """Pre-compute InducedAeList for every grid point at init time.

        AeResponse uses mpmath arbitrary-precision arithmetic (~1.5 s/call).
        Since the conductivity profile depends only on the discrete Tb grid
        (not on rheological parameters), we compute once and cache by index.
        The per-likelihood call then does a dict lookup instead of recomputing.
        """
        induction_obs = [k for k in observables
                         if k.startswith('Ae_') or k.startswith('BiAmp_')
                         or k.startswith('BiPhase_')]
        # induction_bounds labels (support cuts, e.g. Europa's
        # |Ae_synodic| > 0.7) need the Ae grid too, even with no Gaussian
        # induction observable configured.
        bound_labels = set(getattr(self.config, 'induction_bounds', {}) or {})
        if not induction_obs and not bound_labels:
            return

        # Only applies to v2.1 grid format {'Tb_K_grid': ..., 'structures': [...]}
        sd = self.structure_data
        if not (isinstance(sd, dict) and 'Tb_K_grid' in sd and 'structures' in sd):
            return

        from .forward_models import forward_model_induction

        structures = sd['structures']
        n = len(structures)
        log.info(f"Pre-computing induction Ae for {n} grid points ...")

        # Collect all requested frequency labels from observables + bounds.
        requested_labels: set = set(bound_labels)
        for k in observables:
            if k.startswith('Ae_') and k.endswith('_real'):
                requested_labels.add(k[len('Ae_'):-len('_real')])
            elif k.startswith('Ae_') and k.endswith('_imag'):
                requested_labels.add(k[len('Ae_'):-len('_imag')])
            elif k.startswith('BiAmp_'):
                requested_labels.add(k[len('BiAmp_'):])
            elif k.startswith('BiPhase_') and k.endswith('_deg'):
                requested_labels.add(k[len('BiPhase_'):-len('_deg')])

        for i, struct in enumerate(structures):
            Texc_hr_full = struct.get('Texc_hr') if isinstance(struct, dict) else None
            if not Texc_hr_full:
                self._ae_grid_cache[i] = None
                continue
            freq_dict = {lbl: Texc_hr_full[lbl]
                         for lbl in requested_labels if lbl in Texc_hr_full}
            if not freq_dict:
                self._ae_grid_cache[i] = None
                continue
            result = forward_model_induction(struct, freq_dict, nn=1, do_parallel=False)
            self._ae_grid_cache[i] = result
            log.debug(f"  Grid point {i+1}/{n}: Ae computed")

        log.info("Induction Ae pre-computation complete.")

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

                if use_derived_rho_sil:
                    # v2 path: mass-conserve rho_sil from the sampled core,
                    # then recompute CMR² from the assembled body (sampled
                    # core + derived silicate + cached hydrosphere). Factored
                    # into _derive_cmr2_via_mass_conservation (restructure
                    # only, 2026-07 CMR2-reporting-path fix — logic below is
                    # unchanged, just moved into a method so run()'s
                    # cmr2_results and generate_sbi_dataset's CMR2 column can
                    # share the identical derivation instead of silently
                    # falling back to the no-core cache placeholder).
                    derived = self._derive_cmr2_via_mass_conservation(theta_dict)
                    if derived is None:
                        return -1e30
                    cmr2_val = derived
                else:
                    # v1 path: read precomputed CMR² from the cache as-is.
                    # Pick the per-sample structure dict for CMR² lookup.
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
            # induction_bounds (ratified 2026-07-12): one-sided support cuts
            # on Ae per label (amp_min: reject |Ae| < amp_min; im_abs_max:
            # reject |Im Ae| > im_abs_max). Hard rejection, not chi^2 terms.
            induction_bounds = getattr(self.config, 'induction_bounds', {}) or {}
            need_induction = bool(induction_keys_real or induction_keys_imag
                                  or induction_keys_amp or induction_keys_phase
                                  or induction_bounds)
            if need_induction:
                # Collect all requested frequency labels.
                requested_labels = set(induction_bounds)
                for k in induction_keys_real:
                    requested_labels.add(k[len('Ae_'):-len('_real')])
                for k in induction_keys_imag:
                    requested_labels.add(k[len('Ae_'):-len('_imag')])
                for k in induction_keys_amp:
                    requested_labels.add(k[len('BiAmp_'):])
                for k in induction_keys_phase:
                    requested_labels.add(k[len('BiPhase_'):-len('_deg')])

                # Fast path: look up pre-computed Ae from the grid cache
                # (avoids ~1.5 s/call mpmath cost on every likelihood eval).
                Tb_sample = theta_dict.get('Tb_K')
                if (self._ae_grid_cache
                        and Tb_sample is not None
                        and 'Tb_K_grid' in structure_data
                        and 'structures' in structure_data):
                    grid_Tb = np.asarray(structure_data['Tb_K_grid'])
                    grid_idx = int(np.argmin(np.abs(grid_Tb - Tb_sample)))
                    Ae_dict = self._ae_grid_cache.get(grid_idx)
                else:
                    # Fallback: compute on-the-fly (single-struct or no-grid configs).
                    from .forward_models import forward_model_induction
                    Texc_hr_full = modified.get('Texc_hr')
                    if not Texc_hr_full:
                        return -1e30
                    freq_dict = {label: Texc_hr_full[label]
                                 for label in requested_labels
                                 if label in Texc_hr_full}
                    if not freq_dict:
                        return -1e30
                    Ae_dict = forward_model_induction(modified, freq_dict)

                if Ae_dict is None:
                    return -1e30
                for label in requested_labels:
                    Ae = Ae_dict.get(label)
                    if Ae is None or not np.isfinite(complex(Ae).real):
                        return -1e30
                    Ae = complex(Ae)
                    bounds_spec = induction_bounds.get(label)
                    if bounds_spec:
                        amp_min = bounds_spec.get('amp_min')
                        if amp_min is not None and abs(Ae) < float(amp_min):
                            return -1e30
                        im_abs_max = bounds_spec.get('im_abs_max')
                        if im_abs_max is not None and abs(Ae.imag) > float(im_abs_max):
                            return -1e30
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

    def _load_cmr2_offset_sidecar(self) -> Optional[Dict[str, np.ndarray]]:
        """Lazily load the per-Tb CMR2 discretization-offset sidecar, once.

        Test52 Phase 2 (plans/test52-differentiated-titan-plan.md, decision
        1). ``structure_derivation.compute_cmr2`` evaluated on the reduced
        (core + derived-silicate + cached-hydrosphere) shell stack
        under-reports CMR2 by a small, stable, core-independent offset
        relative to the PP-native CMR2mean for the cached template (pure
        hydrosphere radial-discretization error). ``cache_builder`` records
        this offset per Tb grid point in a sidecar JSON next to the
        structure cache pickle, named ``<cache_stem>_offsets.json`` (e.g.
        ``titan_diff_noocean_structure_grid.pkl`` ->
        ``titan_diff_noocean_structure_grid_offsets.json``).

        Looked up exactly once per ``MCMCRunner`` instance and cached on
        ``self._cmr2_offset_sidecar``/``self._cmr2_offset_sidecar_checked``
        (loop callers hit this every likelihood evaluation, so re-reading
        the file per-call would be wasteful and is unnecessary — the
        sidecar is immutable for the lifetime of a run).

        Returns
        -------
        None
            if ``config.structure_cache_path`` is unset or no sidecar file
            exists next to it. This is the case for every config that
            predates Test52 (including all Callisto configs) — by design,
            those configs get NO offset applied and their CMR2 numerics are
            therefore byte-identical to before this change. Only configs
            whose cache has been (re)built with a sidecar (currently just
            Test52) are affected. This is an explicit, documented decision,
            not an oversight: Callisto's own +0.00095-scale offset is not
            yet applied because its cache has not been regenerated with a
            sidecar (tracked separately; see plan decision 1 disclosure
            note re: "+0.23sigma shift" upon Callisto regeneration).
        dict
            ``{'Tb_K_grid': ndarray, 'offsets': ndarray}`` otherwise, sorted
            to match the sidecar's own grid ordering (verified ascending at
            build time; not re-sorted here).
        """
        if self._cmr2_offset_sidecar_checked:
            return self._cmr2_offset_sidecar

        self._cmr2_offset_sidecar_checked = True
        cache_path_str = getattr(self.config, 'structure_cache_path', None)
        if not cache_path_str:
            return None

        cache_path = getattr(self, '_resolved_cache_path', None) or Path(cache_path_str)
        sidecar_path = cache_path.with_name(cache_path.stem + '_offsets.json')
        if not sidecar_path.exists():
            log.debug(
                f"No CMR2 offset sidecar at {sidecar_path} — proceeding "
                f"without an anchor correction (unchanged pre-Test52 "
                f"behavior for this config)."
            )
            return None

        with open(sidecar_path, 'r') as f:
            data = json.load(f)
        Tb_grid = np.asarray(data['Tb_K_grid'], dtype=float)
        offsets = np.asarray(data['CMR2_offset_per_point'], dtype=float)
        sidecar = {'Tb_K_grid': Tb_grid, 'offsets': offsets}
        self._cmr2_offset_sidecar = sidecar
        log.info(
            f"CMR2 offset anchor sidecar found and will be applied: "
            f"{sidecar_path} (n={len(offsets)} Tb points, "
            f"mean={float(np.mean(offsets)):.6f}, "
            f"std={float(np.std(offsets)):.2e})"
        )
        return sidecar

    def _derive_cmr2_via_mass_conservation(
        self, theta_dict: Dict[str, float]
    ) -> Optional[float]:
        """Re-derive CMR² for v2 (mass-conservation) configs.

        Verbatim factoring (2026-07 CMR2-reporting-path fix; pure
        restructure, no numerical change) of the CMR² branch that used to
        live inline in ``_make_flexible_log_likelihood``'s ``'CMR2'``
        observable handling, reached only when
        ``config.derived_params.rho_sil_kgm3.derivation ==
        'mass_conservation'``. Given the sampled core (``R_core_km``,
        ``rho_core_kgm3``), re-derives silicate density via mass
        conservation (``structure_derivation.derive_silicate_density``)
        against the cached hydrosphere
        (``structure_derivation.extract_hydrosphere_layers``), then
        assembles core + silicate + hydrosphere shells and recomputes CMR²
        (``structure_derivation.compute_cmr2``).

        This is now the ONLY place this derivation sequence is
        implemented: both the likelihood (``_make_flexible_log_
        likelihood``, via this same method) and the reporting paths
        (``_compute_model_cmr2``, used by ``run()``'s ``cmr2_results`` and
        ``generate_sbi_dataset``'s ``CMR2`` column) call it, so they can
        never again diverge. (Bug: reporting previously called
        ``_get_cache_scalar('CMR2')`` unconditionally, which ignores
        sampled core parameters and reports the cache-builder's no-core
        placeholder for v2 configs.)

        Test52 Phase 2 additions (plans/test52-differentiated-titan-plan.md):

        1. **CMR2 discretization-offset anchor.** If a sidecar offset file
           exists next to ``config.structure_cache_path`` (see
           ``_load_cmr2_offset_sidecar``), the per-Tb offset — linearly
           interpolated in ``Tb_K`` — is ADDED to the reconstructed CMR2
           before it is returned, correcting the reduced-stack
           radial-discretization under-report identified for Test52. When
           no sidecar exists (every pre-Test52 config, including all
           current Callisto configs), NO offset is added and this method's
           numerics are byte-identical to before this change — an explicit
           decision, not an oversight (see plan decision 1).
        2. **Density-inversion guard.** If
           ``derived_params.rho_sil_kgm3.density_inversion_guard`` is
           truthy in the config, a sample whose derived ``rho_sil`` exceeds
           the sampled ``rho_core_kgm3`` (gravitationally unstable: mantle
           denser than core) is hard-rejected (``None``) exactly like the
           existing ``reject_if_outside_bounds`` check. Absent or false for
           every config that predates this key (including Callisto), so
           this is a no-op there.

        Both additions set ``self._last_cmr2_reject_reason`` (see
        ``__init__``) to a short string identifying which check produced a
        ``None`` return, so ``generate_sbi_dataset`` can distinguish a
        density-inversion rejection from other reject-causes for its
        support-guard bookkeeping without duplicating this derivation.

        Returns
        -------
        None
            if the sample must be hard-rejected: missing/non-finite
            ``R_body_m``/``Mtot_kg`` on the selected structure or missing
            core params in ``theta_dict``; non-contiguous or empty
            hydrosphere; derived ``rho_sil`` outside ``bounds`` when
            ``reject_if_outside_bounds`` is True; derived ``rho_sil >
            rho_core_kgm3`` when ``density_inversion_guard`` is True; or a
            degenerate mass/radius in ``compute_cmr2``. Likelihood callers
            must treat ``None`` exactly like the original inline
            early-returns (return ``-1e30`` immediately, short-circuiting
            the rest of the likelihood). Reporting callers should map
            ``None`` to ``NaN``.
        float
            the derived CMR² otherwise (offset-corrected when a sidecar is
            present). May be non-finite in the (pathological, not hit by
            any current config) case where ``reject_if_outside_bounds`` is
            False and an upstream input is non-finite despite passing the
            earlier explicit checks — the original code guarded this with
            ``if np.isfinite(cmr2_val)`` before adding the chi² term rather
            than rejecting; callers must preserve that same guard rather
            than treating a non-finite return here as a hard reject.
        """
        from .structure_derivation import (
            compute_cmr2,
            derive_silicate_density,
            extract_hydrosphere_layers,
        )

        self._last_cmr2_reject_reason = None

        derived_params_cfg = getattr(self.config, 'derived_params', {}) or {}
        rho_sil_cfg = derived_params_cfg.get('rho_sil_kgm3', {}) or {}
        rho_sil_bounds = tuple(rho_sil_cfg.get('bounds', (2200.0, 3500.0)))
        rho_sil_reject = bool(rho_sil_cfg.get('reject_if_outside_bounds', True))
        density_inversion_guard = bool(rho_sil_cfg.get('density_inversion_guard', False))

        structure_data = self.structure_data
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

        R_body_m = struct_for_cmr2.get('R_body_m', np.nan)
        M_total_kg = struct_for_cmr2.get('Mtot_kg', np.nan)
        R_core_km = theta_dict.get('R_core_km')
        rho_core_kgm3 = theta_dict.get('rho_core_kgm3')
        if not (np.isfinite(R_body_m) and np.isfinite(M_total_kg)
                and R_core_km is not None and rho_core_kgm3 is not None):
            self._last_cmr2_reject_reason = 'missing_inputs'
            return None
        R_core_m = float(R_core_km) * 1000.0
        try:
            hydro_layers, R_oceanbot_m, M_hydro_kg = (
                extract_hydrosphere_layers(struct_for_cmr2)
            )
        except (KeyError, ValueError):
            self._last_cmr2_reject_reason = 'hydrosphere_error'
            return None
        if not hydro_layers:
            self._last_cmr2_reject_reason = 'hydrosphere_empty'
            return None
        rho_sil, accepted = derive_silicate_density(
            M_total_kg=float(M_total_kg),
            M_hydrosphere_kg=M_hydro_kg,
            R_oceanbot_m=R_oceanbot_m,
            R_core_m=R_core_m,
            rho_core_kgm3=float(rho_core_kgm3),
            bounds=rho_sil_bounds,
        )
        if rho_sil_reject and not accepted:
            self._last_cmr2_reject_reason = 'rho_sil_bounds'
            return None
        # Density-inversion guard (Test52 Phase 2, plan decision 3): reject
        # gravitationally unstable configurations where the derived
        # silicate mantle is denser than the sampled core. Only active when
        # the config explicitly sets density_inversion_guard=True (absent
        # for every pre-Test52 config, including Callisto, so this is a
        # no-op there).
        if (density_inversion_guard and np.isfinite(rho_sil)
                and rho_sil > float(rho_core_kgm3)):
            self._last_cmr2_reject_reason = 'density_inversion'
            return None
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
            self._last_cmr2_reject_reason = 'cmr2_valueerror'
            return None

        # CMR2 discretization-offset anchor (Test52 Phase 2, plan decision
        # 1): add the per-Tb offset recorded in the cache's sidecar file,
        # if one exists. No sidecar (every pre-Test52 config) -> no offset,
        # numerics unchanged.
        offset_sidecar = self._load_cmr2_offset_sidecar()
        if offset_sidecar is not None and Tb_sample is not None:
            offset = float(np.interp(
                Tb_sample, offset_sidecar['Tb_K_grid'], offset_sidecar['offsets']
            ))
            cmr2_val = cmr2_val + offset

        return cmr2_val

    def _compute_model_cmr2(self, theta_dict: Dict[str, float]) -> float:
        """Reporting-path CMR²: the single call site for ``run()``'s
        ``cmr2_results`` and ``generate_sbi_dataset``'s ``CMR2`` column.

        Dispatches on the same ``derived_params.rho_sil_kgm3.derivation ==
        'mass_conservation'`` flag the likelihood uses:

        - v2 (mass-conservation) configs (e.g. ``callisto_nacl_andrade_8D``):
          calls ``_derive_cmr2_via_mass_conservation`` — the identical
          core-sensitive derivation the likelihood evaluates for its
          ``'CMR2'`` observable term. A hard-reject (``None``) is mapped to
          ``NaN`` since there is no log-likelihood to reject on a reporting
          path.
        - v1 configs with no ``derived_params`` (e.g. Titan Test50): falls
          back to ``_get_cache_scalar(theta_dict, 'CMR2')`` — the
          cache-builder's precomputed, core-independent CMR², unchanged
          from before this fix.

        Fixes a reporting-path bug (2026-07): ``InferenceResult.cmr2_results``
        and the SBI dataset's ``CMR2`` column used to call
        ``_get_cache_scalar`` unconditionally for every config, silently
        ignoring sampled ``R_core_km``/``rho_core_kgm3`` for v2 configs and
        reporting a CMR² that could disagree with the likelihood's own CMR²
        by several sigma (verified: Callisto reported ~0.3412 vs. the true
        likelihood value ~0.3558). The likelihood itself was never affected
        by this bug — only these two reporting paths.
        """
        derived_params_cfg = getattr(self.config, 'derived_params', {}) or {}
        rho_sil_cfg = derived_params_cfg.get('rho_sil_kgm3', {}) or {}
        use_derived_rho_sil = rho_sil_cfg.get('derivation') == 'mass_conservation'
        if use_derived_rho_sil:
            derived = self._derive_cmr2_via_mass_conservation(theta_dict)
            return derived if derived is not None else np.nan
        return self._get_cache_scalar(theta_dict, 'CMR2')

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
            cmr2_results.append(self._compute_model_cmr2(theta_dict))
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

    def generate_sbi_dataset(self, n_samples: int = 10000, output_path: Optional[str] = None,
                             apply_support_guard: bool = False, imag_convention: str = 'signed',
                             drop_nonfinite: bool = False, seed: Optional[int] = None,
                             provenance: Optional[dict] = None,
                             obs_noise: bool = False, noise_seed: Optional[int] = None):
        """
        Generate (theta, x) dataset by sampling from the prior.

        Args:
            n_samples: Number of simulations to run
            output_path: Optional path to save .npz file
            apply_support_guard: If True, apply the same no-ocean phase-stability
                guard used by ``_make_flexible_log_likelihood`` (see
                ``sampler_settings['phase_stability']``) and drop any sampled
                row for which it would trigger, so the SBI training-set support
                matches the MCMC posterior support. Rejected rows are excluded
                from the returned/saved arrays and counted in
                ``n_rejected_support``. Default False preserves current
                behavior byte-identically (no rows dropped).
            imag_convention: 'signed' (default, unchanged behavior) stores the
                signed Im(k2) for observable name 'Im_k2' and silently NaN-fills
                any unrecognized observable name (including 'abs_Im_k2'), exactly
                as before. 'abs' stores |Im(k2)| for both the 'Im_k2' and
                'abs_Im_k2' observable names (matching the MCMC likelihood's and
                GUI's |Im k2| convention). Any other value raises ValueError.
                Induction channels (``Ae_<label>_real/_imag``,
                ``BiAmp_<label>``, ``BiPhase_<label>_deg``) are computed from
                the same precomputed Tb-grid Ae cache the likelihood uses
                (2026-07-10; previously these were silently NaN-filled).
                NaN when the cache/label is unavailable — pair with
                ``drop_nonfinite=True`` for SBI datasets. Note the BiPhase
                channel is in degrees; ``obs_noise`` adds unwrapped Gaussian
                noise there (the likelihood wraps residuals, which only
                matters near +/-180 deg).
            drop_nonfinite: If True, drop any row whose x vector contains a
                non-finite value (e.g. TidalPy forward-model failures), counted
                in ``n_rejected_nonfinite``. Default False preserves current
                behavior byte-identically (NaN rows kept). If the combined
                rejection rate (support guard + non-finite) exceeds 20%, a loud
                warning is logged (not raised).
            seed: Optional int to seed prior sampling reproducibly. pocoMC's
                installed ``Prior`` class (checked at
                site-packages/pocomc/prior.py, pocomc 1.2.6) exposes only
                ``.rvs(size)`` — there is no ``.sample()`` method and ``.rvs``
                accepts no ``random_state``/seed argument of its own (it calls
                each underlying scipy-frozen distribution's ``.rvs(size=size)``
                with no random_state, so draws depend on the global numpy RNG
                state). Reproducible seeding is therefore done via
                ``np.random.seed(seed)`` immediately before sampling. Default
                None leaves the global RNG state untouched, preserving current
                behavior.
            provenance: Optional dict merged (additively) into the saved .npz
                metadata when any of the above opt-in options are used. Ignored
                when all new kwargs are left at their defaults.
            obs_noise: If True (requires imag_convention='abs'), add Gaussian
                observation noise to every x column after the forward model:
                x_k += N(0, sigma_k) with sigma_k taken from
                ``config.observables[name][1]``. The noise is added AFTER the
                |Im k2| fold and the result is NOT re-folded, matching the
                committed diagonal-Gaussian likelihood's implied generative
                model (data = |Im_model| + noise, may be negative near zero).
                Required for NPE training to target the same posterior the
                MCMC likelihood defines (ratified 2026-07-09; without it the
                flow targets the singular noiseless conditional and fails
                SBC/crosscheck). Default False preserves existing behavior
                byte-identically.
            noise_seed: Seed for a dedicated numpy Generator used only for
                observation noise (independent of the prior-draw ``seed`` so
                the same theta set can be re-noised reproducibly). None draws
                from a fresh nondeterministic Generator.

        Returns:
            (theta, x) tuple. When ``apply_support_guard`` or
            ``drop_nonfinite`` is used, ``theta``/``x`` contain fewer than
            ``n_samples`` rows (rejected rows excluded). Rejection statistics
            are never returned inline (to keep this return signature stable);
            they are attached to ``self.last_sbi_dataset_stats`` (dict) and,
            when ``output_path`` is given, saved into the .npz as additive
            keys. The default call (no new kwargs) saves exactly the original
            4 .npz keys (theta, x, param_names, obs_names) unchanged.
        """
        if imag_convention not in ('signed', 'abs'):
            raise ValueError(
                f"imag_convention must be 'signed' or 'abs', got '{imag_convention}'"
            )
        if obs_noise and imag_convention != 'abs':
            raise ValueError(
                "obs_noise=True requires imag_convention='abs': the noise model "
                "is data = |Im_model| + noise, matching the likelihood's |Im k2| "
                "convention; with signed Im the generative model is ambiguous."
            )

        log.info(f"Generating SBI dataset with {n_samples} samples...")

        # Reproducible seeding: pocoMC's Prior.rvs() has no seed/random_state
        # argument of its own (verified against installed pocomc/prior.py),
        # so we seed the global numpy RNG that the underlying scipy frozen
        # distributions draw from.
        if seed is not None:
            np.random.seed(seed)

        # Sample from prior
        theta = self.prior.rvs(n_samples)

        x = []
        theta_kept = []
        obs_names = list(self.config.observables.keys())

        from .forward_models import apply_parameters, forward_model_k2_flexible

        # Support guard: identical logic to
        # MCMCRunner._make_flexible_log_likelihood's nested _check_no_ocean
        # (kept in sync manually — that guard is a closure local to the
        # likelihood builder and is not otherwise importable without
        # refactoring the likelihood, which is out of scope here).
        phase_stability_cfg = self.config.sampler_settings.get('phase_stability', {})
        no_ocean_guard = apply_support_guard and phase_stability_cfg.get('enforce') == 'no_ocean_Ih'
        no_ocean_margin_K = float(phase_stability_cfg.get('margin_K', 0.1))

        # Density-inversion guard (Test52 Phase 2, plan decision 3): wired
        # into the support-guard counting exactly like the no-ocean guard
        # above, but only takes effect when BOTH apply_support_guard=True
        # AND the config's derived_params.rho_sil_kgm3.density_inversion_
        # guard is truthy (absent for every pre-Test52 config, including
        # Callisto). With apply_support_guard's default of False, or for
        # any config without this key, density_inversion_guard_active is
        # False and the block below never executes — behavior for every
        # existing config/call site is therefore unchanged.
        _derived_params_cfg = getattr(self.config, 'derived_params', {}) or {}
        _rho_sil_cfg = _derived_params_cfg.get('rho_sil_kgm3', {}) or {}
        density_inversion_guard_active = (
            apply_support_guard
            and _rho_sil_cfg.get('derivation') == 'mass_conservation'
            and bool(_rho_sil_cfg.get('density_inversion_guard', False))
        )

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

        def _induction_channel_value(name, theta_dict):
            """Model value for an induction observable channel, mirroring the
            likelihood's grid-cache lookup (nearest Tb index into the Ae grid
            precomputed at __init__). Returns NaN when the cache/label is
            unavailable — with drop_nonfinite=True such rows are rejected and
            counted, closing the old landmine where these channels were
            silently NaN-filled in SBI datasets.
            """
            Tb_sample = theta_dict.get('Tb_K')
            if (not self._ae_grid_cache or Tb_sample is None
                    or 'Tb_K_grid' not in self.structure_data):
                return np.nan
            grid_Tb = np.asarray(self.structure_data['Tb_K_grid'])
            Ae_dict = self._ae_grid_cache.get(
                int(np.argmin(np.abs(grid_Tb - Tb_sample))))
            if Ae_dict is None:
                return np.nan
            if name.startswith('Ae_') and name.endswith('_real'):
                label, part = name[len('Ae_'):-len('_real')], 'real'
            elif name.startswith('Ae_') and name.endswith('_imag'):
                label, part = name[len('Ae_'):-len('_imag')], 'imag'
            elif name.startswith('BiAmp_'):
                label, part = name[len('BiAmp_'):], 'amp'
            elif name.startswith('BiPhase_') and name.endswith('_deg'):
                label, part = name[len('BiPhase_'):-len('_deg')], 'phase'
            else:
                return np.nan
            Ae = Ae_dict.get(label)
            if Ae is None:
                return np.nan
            Ae = complex(Ae)
            if part == 'real':
                return Ae.real
            if part == 'imag':
                return Ae.imag
            if part == 'amp':
                return abs(Ae)
            return float(np.degrees(np.angle(Ae)))

        # induction_bounds support cuts (ratified 2026-07-12): reject rows
        # violating the one-sided Ae constraints so the SBI training support
        # matches the MCMC likelihood support exactly (same precedent as the
        # no-ocean guard). Uses the same Tb-grid Ae lookup as the likelihood.
        induction_bounds_cfg = (getattr(self.config, 'induction_bounds', {}) or {}
                                if apply_support_guard else {})

        def _check_induction_bounds(theta_dict) -> bool:
            """True when the sample violates an induction bound (reject)."""
            for label, spec in induction_bounds_cfg.items():
                Tb_sample = theta_dict.get('Tb_K')
                if (not self._ae_grid_cache or Tb_sample is None
                        or 'Tb_K_grid' not in self.structure_data):
                    return True  # bound configured but unevaluable -> reject
                grid_Tb = np.asarray(self.structure_data['Tb_K_grid'])
                Ae_dict = self._ae_grid_cache.get(
                    int(np.argmin(np.abs(grid_Tb - Tb_sample))))
                Ae = None if Ae_dict is None else Ae_dict.get(label)
                if Ae is None:
                    return True
                Ae = complex(Ae)
                amp_min = spec.get('amp_min')
                if amp_min is not None and abs(Ae) < float(amp_min):
                    return True
                im_abs_max = spec.get('im_abs_max')
                if im_abs_max is not None and abs(Ae.imag) > float(im_abs_max):
                    return True
            return False

        n_rejected_support = 0
        n_rejected_nonfinite = 0

        t0 = time.time()
        for i in range(n_samples):
            theta_dict = self._expand_theta(theta[i])

            if apply_support_guard:
                modified = apply_parameters(theta_dict, self.structure_data)
                if _check_no_ocean(modified):
                    n_rejected_support += 1
                    continue
                if induction_bounds_cfg and _check_induction_bounds(theta_dict):
                    n_rejected_support += 1
                    continue

            cmr2_precomputed = None
            if density_inversion_guard_active:
                # Evaluate the CMR2 derivation now (it internally applies
                # the density-inversion guard and records the cause in
                # self._last_cmr2_reject_reason). Only an inversion-caused
                # rejection is counted/dropped as a support rejection here;
                # any other None-cause (e.g. rho_sil bounds) keeps its
                # existing behavior of flowing through as NaN in the CMR2
                # column below, unchanged from before this fix.
                cmr2_precomputed = self._compute_model_cmr2(theta_dict)
                if (not np.isfinite(cmr2_precomputed)
                        and self._last_cmr2_reject_reason == 'density_inversion'):
                    n_rejected_support += 1
                    continue

            # Compute k2 (and h2 — same forward call returns both)
            Re_k2, Im_k2, Re_h2, Im_h2, _ = forward_model_k2_flexible(
                theta_dict, self.structure_data,
                return_heating=False, arrhenius_params=self.arrhenius_params
            )

            xi = []
            for name in obs_names:
                if name == 'Re_k2':
                    xi.append(Re_k2)
                elif name == 'Im_k2':
                    xi.append(abs(Im_k2) if imag_convention == 'abs' else Im_k2)
                elif name == 'abs_Im_k2' and imag_convention == 'abs':
                    xi.append(abs(Im_k2))
                elif name == 'k2':
                    xi.append(np.sqrt(Re_k2 ** 2 + Im_k2 ** 2))
                # h2 channels (2026-07-10; previously silently NaN-filled).
                # Mirrors the likelihood: Re_h2 signed; the Im channel is
                # |Im h2| under 'abs' (both spellings), signed 'Im_h2' under
                # 'signed'; 'h2' is the modulus.
                elif name == 'Re_h2':
                    xi.append(Re_h2)
                elif name == 'Im_h2':
                    xi.append(abs(Im_h2) if imag_convention == 'abs' else Im_h2)
                elif name == 'abs_Im_h2' and imag_convention == 'abs':
                    xi.append(abs(Im_h2))
                elif name == 'h2':
                    xi.append(np.sqrt(Re_h2 ** 2 + Im_h2 ** 2))
                elif name == 'CMR2':
                    xi.append(cmr2_precomputed if cmr2_precomputed is not None
                               else self._compute_model_cmr2(theta_dict))
                elif name == 'Mtot_kg':
                    xi.append(self._get_cache_scalar(theta_dict, 'Mtot_kg'))
                elif (name.startswith('Ae_') or name.startswith('BiAmp_')
                        or name.startswith('BiPhase_')):
                    xi.append(_induction_channel_value(name, theta_dict))
                else:
                    xi.append(np.nan)

            if drop_nonfinite and not np.all(np.isfinite(xi)):
                n_rejected_nonfinite += 1
                continue

            x.append(xi)
            theta_kept.append(theta[i])

            if (i + 1) % 100 == 0:
                elapsed = time.time() - t0
                eta = (elapsed / (i + 1)) * (n_samples - (i + 1))
                log.info(f"  {i+1}/{n_samples} simulations complete (ETA: {eta/60:.1f} min)")

        theta = np.array(theta_kept) if (apply_support_guard or drop_nonfinite) else np.array(theta)
        x = np.array(x)

        obs_noise_meta = None
        if obs_noise and len(x):
            # Diagonal Gaussian observation noise, sigma per observable from the
            # config. Added after the |Im| fold, NOT re-folded (see docstring).
            sigmas = np.array([float(self.config.observables[name][1])
                               for name in obs_names])
            noise_rng = np.random.default_rng(noise_seed)
            x = x + noise_rng.normal(0.0, 1.0, size=x.shape) * sigmas
            obs_noise_meta = {
                'type': 'gaussian_diagonal',
                'sigma': {name: float(s) for name, s in zip(obs_names, sigmas)},
                'noise_seed': noise_seed,
                'refold_im': False,
            }

        n_rejected_total = n_rejected_support + n_rejected_nonfinite
        rejection_rate = (n_rejected_total / n_samples) if n_samples else 0.0
        if (apply_support_guard or drop_nonfinite) and rejection_rate > 0.20:
            log.warning(
                f"generate_sbi_dataset: rejection_rate={rejection_rate:.1%} exceeds 20% "
                f"(n_rejected_support={n_rejected_support}, "
                f"n_rejected_nonfinite={n_rejected_nonfinite}, n_requested={n_samples}). "
                f"Training-set support may be significantly reduced — inspect before proceeding."
            )

        stats = {
            'n_requested': n_samples,
            'n_rejected_support': n_rejected_support,
            'n_rejected_nonfinite': n_rejected_nonfinite,
            'rejection_rate': rejection_rate,
            'imag_convention': imag_convention,
            'seed': seed,
        }
        if obs_noise_meta is not None:
            stats['obs_noise'] = obs_noise_meta
        self.last_sbi_dataset_stats = stats

        using_new_options = (apply_support_guard or imag_convention != 'signed'
                             or drop_nonfinite or seed is not None or provenance is not None
                             or obs_noise)

        if output_path:
            if not using_new_options:
                # Exact original save path — untouched.
                np.savez(output_path, theta=theta, x=x, param_names=self.param_names, obs_names=obs_names)
            else:
                # Additive keys only; old readers of the 4 original keys are unaffected.
                from .parameter_registry import PARAMETER_REGISTRY

                param_bounds = []
                for name in self.param_names:
                    cfg = self.config.param_space[name]
                    bounds = cfg.get('bounds')
                    if bounds is not None:
                        param_bounds.append([float(bounds[0]), float(bounds[1])])
                    else:
                        # No natural hard bound for this prior_type (e.g. 'normal').
                        param_bounds.append([np.nan, np.nan])
                param_bounds = np.array(param_bounds)

                param_units = [
                    (PARAMETER_REGISTRY[name].units or '') if name in PARAMETER_REGISTRY else ''
                    for name in self.param_names
                ]

                try:
                    config_hash = self.config.generate_hash()
                except Exception:
                    config_hash = ''

                try:
                    git_sha = subprocess.check_output(
                        ['git', 'rev-parse', '--short', 'HEAD'],
                        cwd=os.path.dirname(os.path.abspath(__file__)),
                        stderr=subprocess.DEVNULL,
                    ).decode().strip()
                except Exception:
                    git_sha = ''

                save_kwargs = dict(
                    theta=theta, x=x, param_names=self.param_names, obs_names=obs_names,
                    param_bounds=param_bounds, param_units=param_units,
                    config_hash=config_hash, git_sha=git_sha,
                    seed=(seed if seed is not None else -1),
                    n_requested=n_samples,
                    n_rejected_nonfinite=n_rejected_nonfinite,
                    n_rejected_support=n_rejected_support,
                    imag_convention=imag_convention,
                )
                if provenance:
                    save_kwargs['provenance'] = json.dumps(provenance)
                if obs_noise_meta is not None:
                    save_kwargs['obs_noise'] = json.dumps(obs_noise_meta)
                np.savez(output_path, **save_kwargs)
            log.info(f"SBI dataset saved to {output_path}")

        return theta, x
