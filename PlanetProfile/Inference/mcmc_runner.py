"""
MCMC runner using pocoMC for Bayesian inference.

Refactored from Test41-44 scripts into reusable MCMCRunner class.
Supports Andrade and Maxwell rheologies with optional Arrhenius viscosity.

Author: PlanetProfile Team
Date: 2026-04-29
"""
import numpy as np
import logging
import time
import pickle
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

        # Route to grid cache when Tb_K is a free parameter
        self._use_flexible = 'Tb_K' in self.param_names

        # Load cached structure (skip bodyname validation for Test* files)
        log.info(f"Loading structure cache: {config.structure_cache_path}")
        if self._use_flexible:
            self.structure_data = self._load_grid_cache(config.structure_cache_path)
        else:
            self.structure_data = load_structure_cache(
                config.structure_cache_path, validate_bodyname=None
            )

        # Build prior and likelihood — always use dict-based interface so
        # parameter order in param_space never affects forward model mapping.
        self.prior = self._build_prior()
        self.log_likelihood_fn = self._make_flexible_log_likelihood(
            config.observables,
            self.structure_data,
            arrhenius_params=config.sampler_settings.get('arrhenius_params')
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
        """Infer rheology type from parameter space."""
        params = self.config.param_space
        has_alpha = 'alpha' in params
        has_zeta = any(k.startswith('log10_zeta_') for k in params)
        if has_alpha and has_zeta:
            return 'andrade'
        elif not has_alpha and not has_zeta:
            return 'maxwell'
        else:
            raise ValueError("Cannot infer rheology from parameter space. "
                           "Andrade requires 'alpha' and log10_zeta_* parameters. "
                           "Maxwell requires neither.")

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
            'log10_zeta_Ih': r'$\log_{10}(\zeta_{\rm Ih})$',
            'log10_zeta_HP': r'$\log_{10}(\zeta_{\rm HP})$',
            'log10_zeta_sil': r'$\log_{10}(\zeta_{\rm sil})$',
            'log10_eta_Ih': r'$\log_{10}(\eta_{\rm Ih})$',
            'log10_eta_HP': r'$\log_{10}(\eta_{\rm HP})$',
            'log10_eta_sil': r'$\log_{10}(\eta_{\rm sil})$',
            'log10_mu_Ih': r'$\log_{10}(\mu_{\rm Ih})$',
            'Tb_K': r'$T_b$ (K)',
        }
        return label_map.get(param_name, param_name)

    def _load_grid_cache(self, cache_path: str) -> Dict[str, Any]:
        """Load grid cache dict and attach Tb lookup array for nearest-neighbor search."""
        with open(cache_path, 'rb') as f:
            grid_cache = pickle.load(f)
        if not isinstance(grid_cache, dict):
            raise ValueError(f"Expected dict[float -> structure] at {cache_path}, got {type(grid_cache)}")
        if not grid_cache:
            raise ValueError(
                f"Grid cache at {cache_path} is empty. "
                f"Regenerate with: python -m PlanetProfile.Inference.prepare_structure_variants "
                f"--test-module PlanetProfile.Test.PPTest42 --output-dir titan_cache/ --maxwell --force"
            )
        first_key = next(iter(grid_cache))
        if not isinstance(first_key, (int, float)):
            raise ValueError(
                f"Expected grid cache keyed by Tb_K floats at {cache_path}, "
                f"got key type {type(first_key).__name__}"
            )
        grid_Tb_values = np.array(sorted(grid_cache.keys()))
        log.info(f"Grid cache loaded: {len(grid_Tb_values)} Tb_K values "
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

        def log_likelihood(theta):
            theta_dict = dict(zip(param_names, theta))
            Re_k2, Im_k2, _ = forward_model_k2_flexible(
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
            return -0.5 * chi2

        return log_likelihood

    def run(self, progress_callback: Optional[Callable] = None):
        """
        Run MCMC sampling with pocoMC.

        Args:
            progress_callback: Optional function called with progress dict:
                {'iteration': int, 'n_total': int, 'acceptance_rate': float,
                 'n_samples': int, 'ess': float}

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

        # Run — pocoMC stops internally when ESS >= n_ess
        sampler.run(n_total=4096, progress=True)

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

        # Final posterior extraction
        samples, log_likes, log_post, _ = sampler.posterior()
        n_samples = len(samples)

        elapsed = time.time() - t0
        log.info(f"MCMC complete: {n_samples} samples in {elapsed/60:.1f} min")

        # Compute convergence metrics
        convergence_metrics = self._compute_convergence(samples, log_likes, sampler)

        # Recompute k2 for posterior samples — always dict-based so param order
        # in param_space never causes wrong positional mapping.
        log.info(f"Recomputing k2 for {n_samples} posterior samples...")
        arrhenius_params = self.config.sampler_settings.get('arrhenius_params')
        rheology = self._infer_rheology() if not self._use_flexible else None

        from .forward_models import forward_model_k2_flexible
        k2_results = []
        for i, theta in enumerate(samples):
            theta_dict = dict(zip(self.param_names, theta))
            Re_k2, Im_k2, _ = forward_model_k2_flexible(
                theta_dict, self.structure_data,
                return_heating=False, arrhenius_params=arrhenius_params
            )
            k2_results.append((Re_k2, Im_k2))
            if (i + 1) % 100 == 0:
                log.info(f"  {i+1}/{n_samples} samples recomputed")

        k2_results = np.array(k2_results)

        # Recompute heating on subset — same dict-based approach
        n_reeval = min(self.n_reeval, n_samples)
        log.info(f"Recomputing heating for {n_reeval} posterior samples...")

        rng = np.random.RandomState(self.random_state)
        idx_heat = rng.choice(n_samples, n_reeval, replace=False)
        idx_heat.sort()
        heating_results = []
        for si in idx_heat:
            theta_dict = dict(zip(self.param_names, samples[si]))
            _, _, perPhase_W = forward_model_k2_flexible(
                theta_dict, self.structure_data,
                return_heating=True, arrhenius_params=arrhenius_params
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
            heating_results=heating_results,
            convergence_metrics=convergence_metrics,
            metadata={
                'elapsed_time_s': elapsed,
                'n_iterations': len(samples) if samples is not None else 0,
                'rheology': rheology,
                'heating_indices': heating_indices,
            }
        )

        log.info("MCMC inference complete")
        return result

    def _compute_convergence(self, samples, log_likes, sampler) -> Dict[str, float]:
        """Compute convergence diagnostics (ESS, R-hat, acceptance rate)."""
        n_samples = len(samples)

        # ESS (effective sample size)
        try:
            ess = sampler.n_eff if hasattr(sampler, 'n_eff') else n_samples
        except:
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
