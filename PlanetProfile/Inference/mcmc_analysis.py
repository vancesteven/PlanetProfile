"""
Post-hoc TidalPy re-analysis of MCMC posterior samples.

Loads a saved ``InferenceResult`` (or legacy dict) pickle, reloads its structure
cache, and re-runs the TidalPy radial-solver forward model
(``forward_model_k2_flexible``) on the posterior samples to produce posteriors
of (Re_k2, Im_k2, Re_h2, Im_h2) without re-running MCMC.

This is useful when:
  - The pickle was produced before a forward-model bug fix (e.g. the per-layer
    rheology fix that switched Fe→Elastic, ocean→Newton on 2026-05-22) and the
    stored ``k2_results`` are stale.
  - You want h₂ posteriors but only k₂ was originally cached.
  - You want a quick re-evaluation on a thinned subset of samples.

The function is intentionally a *pure post-processing* step — it does not
modify the pickle, does not rebuild the structure cache, and does not change
sampler state.

Author: PlanetProfile Team
Date: 2026-06-10
"""
from __future__ import annotations

import logging
import pickle
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pickle adapters
# ---------------------------------------------------------------------------

def _extract_samples_and_meta(loaded: Any) -> Dict[str, Any]:
    """Extract a uniform view of a posterior pickle.

    Accepts either an ``InferenceResult`` instance (current format) or a raw
    ``dict`` (legacy format used by Test48 / Test49 / etc.).

    Returns a dict with keys:
        samples         (np.ndarray, n_samples × n_params)
        param_names     (List[str])
        param_groups    (Dict[str, List[str]])    — empty if absent
        fixed_params    (Dict[str, float])        — empty if absent
        structure_cache_path (Optional[str])      — None if absent
        bodyname        (Optional[str])

    Raises ValueError if samples or param_names cannot be located.
    """
    out: Dict[str, Any] = {
        'param_groups': {},
        'fixed_params': {},
        'structure_cache_path': None,
        'bodyname': None,
    }

    # InferenceResult instance (duck-typed: has .config and .samples)
    if hasattr(loaded, 'samples') and hasattr(loaded, 'config'):
        cfg = loaded.config
        out['samples'] = np.asarray(loaded.samples)
        out['param_names'] = list(loaded.param_names)
        out['param_groups'] = dict(getattr(cfg, 'param_groups', {}) or {})
        out['fixed_params'] = dict(getattr(cfg, 'fixed_params', {}) or {})
        out['structure_cache_path'] = getattr(cfg, 'structure_cache_path', None)
        out['bodyname'] = getattr(cfg, 'bodyname', None)
        # Auto-load arrhenius_params: prefer top-level field, fall back to sampler_settings.
        ap = getattr(cfg, 'arrhenius_params', None)
        if not ap:
            ss = getattr(cfg, 'sampler_settings', {}) or {}
            ap = ss.get('arrhenius_params')
        out['arrhenius_params'] = dict(ap) if ap else None
        return out

    # Legacy dict format
    if isinstance(loaded, dict):
        if 'samples' not in loaded:
            raise ValueError(
                "Pickle is a dict but has no 'samples' key. Cannot re-analyse."
            )
        out['samples'] = np.asarray(loaded['samples'])
        # param_names may live at top level or inside a config-like sub-dict
        if 'param_names' in loaded:
            out['param_names'] = list(loaded['param_names'])
        else:
            raise ValueError(
                "Pickle dict is missing 'param_names'. Cannot reconstruct "
                "theta_dict for forward-model re-evaluation."
            )
        out['param_groups'] = dict(loaded.get('param_groups', {}) or {})
        out['fixed_params'] = dict(loaded.get('fixed_params', {}) or {})
        out['structure_cache_path'] = loaded.get('structure_cache_path', None)
        out['bodyname'] = loaded.get('bodyname', None)
        out['arrhenius_params'] = loaded.get('arrhenius_params')
        return out

    raise TypeError(
        f"Loaded object is neither InferenceResult nor dict: {type(loaded).__name__}"
    )


def _load_structure_cache(
    cache_path: str,
    bodyname: Optional[str] = None,
) -> Dict[str, Any]:
    """Load a structure cache, accepting either the v1 single-structure format
    or the v2.1 grid (Tb, D) format. Mirrors MCMCRunner._load_grid_cache logic.

    The forward model accepts both formats — apply_parameters dispatches on the
    'grid_cache' / 'Tb_K_grid' keys when sampling Tb_K. We pass through the
    raw loaded dict whenever it looks grid-shaped.
    """
    from .structure_cache import load_structure_cache

    p = Path(cache_path)
    if not p.exists():
        raise FileNotFoundError(f"Structure cache not found: {cache_path}")

    # Peek at the pickle to decide which loader to call.
    with open(p, 'rb') as f:
        head = pickle.load(f)

    if not isinstance(head, dict):
        raise ValueError(
            f"Structure cache at {cache_path} is not a dict (got {type(head).__name__})"
        )

    # v2.1 list format: top-level Tb_K_grid + structures
    if 'Tb_K_grid' in head and 'structures' in head:
        log.info(
            f"Structure cache: list-format grid with "
            f"{len(head['structures'])} structures (Tb range "
            f"{min(head['Tb_K_grid']):.2f}-{max(head['Tb_K_grid']):.2f} K)"
        )
        return head

    # v2 dict-of-(Tb, D) grid: keys are tuples of floats
    first_key = next(iter(head)) if head else None
    if isinstance(first_key, tuple):
        log.info(f"Structure cache: 2-D grid with {len(head)} (Tb, D) points")
        # Wrap so apply_bottom_temperature finds it under 'grid_cache'
        Tb_vals = np.array(sorted({k[0] for k in head}))
        return {'grid_cache': head, 'grid_Tb_values': Tb_vals}

    # 'grid_cache' envelope wrapper — already in the runner shape
    if 'grid_cache' in head and isinstance(head['grid_cache'], dict):
        log.info("Structure cache: pre-wrapped grid_cache envelope")
        return head

    # Otherwise treat as v1 single-structure cache and use the validated loader
    return load_structure_cache(cache_path, validate_bodyname=bodyname)


def _expand_theta(
    theta: np.ndarray,
    param_names: Sequence[str],
    param_groups: Dict[str, List[str]],
    fixed_params: Dict[str, float],
) -> Dict[str, float]:
    """Mirror MCMCRunner._expand_theta — group expansion + fixed merge."""
    theta_dict = dict(zip(param_names, theta))
    for group_key, members in (param_groups or {}).items():
        if group_key in theta_dict:
            for m in members:
                theta_dict[m] = theta_dict[group_key]
    theta_dict.update(fixed_params or {})
    return theta_dict


def _infer_rheology(param_names: Sequence[str]) -> str:
    """Infer rheology type from parameter names. Andrade requires alpha;
    if absent, fall back to maxwell. The forward model itself raises if the
    inferred type doesn't match the structure_data flags set by hooks.
    """
    pn = set(param_names)
    if 'alpha' in pn or any(n.startswith('log10_zeta') for n in pn):
        return 'andrade'
    return 'maxwell'


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def reanalyze_k2_from_pickle(
    pickle_path: str,
    structure_cache_path: Optional[str] = None,
    n_eval: Optional[int] = None,
    rheology_type: Optional[str] = None,
    arrhenius_params: Optional[Dict[str, Any]] = None,
    return_heating: bool = False,
    random_state: int = 42,
    progress_every: int = 100,
) -> Dict[str, Any]:
    """Re-run TidalPy on the posterior of a saved MCMC pickle.

    Parameters
    ----------
    pickle_path
        Path to an ``InferenceResult`` pickle produced by ``MCMCRunner.run()``,
        or a legacy dict pickle that stores ``samples`` and ``param_names`` at
        the top level.
    structure_cache_path
        Path to the structure cache (.pkl). If ``None``, the function reads
        ``config.structure_cache_path`` from the pickle. If both are missing it
        raises ``ValueError``.
    n_eval
        Number of posterior samples to evaluate. ``None`` means evaluate all.
        Subsampling is reproducible via ``random_state``.
    rheology_type
        ``'andrade'`` or ``'maxwell'``. ``None`` lets the function infer from
        the parameter names. The TidalPy hook system tags
        ``structure_data['rheology_type']`` based on the parameter dict, so
        this kwarg only affects ambiguous cases (e.g. samples with neither
        alpha nor mu).
    arrhenius_params
        Optional Arrhenius temperature-dependent viscosity dict, identical to
        what the original MCMC run used. If ``None``, auto-loaded from the
        pickle's InferenceConfig (preferred top-level field or
        sampler_settings fallback). Pass explicitly to override the pickled
        value. If still not found, ``None`` evaluates the rheology at the
        cached reference viscosity.
    return_heating
        If True, also compute per-phase tidal heating (slower).
    random_state
        Seed for the subsampling RNG.
    progress_every
        Log a progress line every N samples. Set to 0 to silence.

    Returns
    -------
    dict with keys:
        sample_indices   (n_eval,) int — indices into the original samples array
        Re_k2, Im_k2     (n_eval,) float — tidal Love number k₂ components
        Re_h2, Im_h2     (n_eval,) float — radial-displacement Love number h₂
        heating          List[Dict[str, float]] | None — per-phase heating (W)
        n_failed         int — count of samples for which TidalPy returned NaN
        metadata         dict — inputs + run info for provenance

    Raises
    ------
    FileNotFoundError
        If ``pickle_path`` or the resolved ``structure_cache_path`` does not
        exist.
    ValueError
        If the pickle is missing required fields.
    """
    pickle_path = str(pickle_path)
    if not Path(pickle_path).exists():
        raise FileNotFoundError(f"Pickle not found: {pickle_path}")

    log.info(f"Re-analysing posterior from {pickle_path}")
    with open(pickle_path, 'rb') as f:
        loaded = pickle.load(f)

    meta = _extract_samples_and_meta(loaded)
    samples = meta['samples']
    param_names = meta['param_names']
    param_groups = meta['param_groups']
    fixed_params = meta['fixed_params']
    n_samples, _ = samples.shape

    # Auto-load arrhenius_params from pickle if not explicitly passed
    if arrhenius_params is None:
        arrhenius_params = meta.get('arrhenius_params')

    # Resolve structure cache path
    cache_path = structure_cache_path or meta['structure_cache_path']
    if cache_path is None:
        raise ValueError(
            "structure_cache_path is None. Pass it explicitly to "
            "reanalyze_k2_from_pickle() — the pickle does not record it."
        )
    cache_path = str(cache_path)

    structure_data = _load_structure_cache(cache_path, bodyname=meta['bodyname'])

    # Subsample
    if n_eval is None or n_eval > n_samples:
        n_eval_eff = n_samples
        idx = np.arange(n_samples)
    else:
        rng = np.random.RandomState(random_state)
        idx = rng.choice(n_samples, n_eval, replace=False)
        idx.sort()
        n_eval_eff = n_eval

    rheo = rheology_type or _infer_rheology(param_names)

    # Lazy import — TidalPy is imported eagerly inside forward_models.
    from .forward_models import forward_model_k2_flexible

    log.info(
        f"Evaluating {n_eval_eff}/{n_samples} samples "
        f"(rheology={rheo}, heating={return_heating}, arrhenius="
        f"{'on' if arrhenius_params else 'off'})"
    )

    Re_k2 = np.full(n_eval_eff, np.nan)
    Im_k2 = np.full(n_eval_eff, np.nan)
    Re_h2 = np.full(n_eval_eff, np.nan)
    Im_h2 = np.full(n_eval_eff, np.nan)
    heating: Optional[List[Dict[str, float]]] = [] if return_heating else None
    n_failed = 0

    for k, si in enumerate(idx):
        theta_dict = _expand_theta(samples[si], param_names, param_groups, fixed_params)
        rk, ik, rh, ih, perPhase = forward_model_k2_flexible(
            theta_dict,
            structure_data,
            return_heating=return_heating,
            arrhenius_params=arrhenius_params,
        )
        Re_k2[k] = rk
        Im_k2[k] = ik
        Re_h2[k] = rh
        Im_h2[k] = ih
        if heating is not None:
            heating.append(perPhase if perPhase is not None else {})
        if not np.isfinite(rk):
            n_failed += 1
        if progress_every and ((k + 1) % progress_every == 0):
            log.info(
                f"  {k+1}/{n_eval_eff} samples evaluated "
                f"(failed so far: {n_failed})"
            )

    if n_failed:
        log.warning(
            f"{n_failed}/{n_eval_eff} samples returned NaN from TidalPy "
            f"radial_solver — check rheology/arrhenius inputs."
        )

    return {
        'sample_indices': idx,
        'Re_k2': Re_k2,
        'Im_k2': Im_k2,
        'Re_h2': Re_h2,
        'Im_h2': Im_h2,
        'heating': heating,
        'n_failed': n_failed,
        'metadata': {
            'pickle_path': pickle_path,
            'structure_cache_path': cache_path,
            'rheology_type': rheo,
            'arrhenius_params': arrhenius_params,
            'n_samples_total': n_samples,
            'n_eval': n_eval_eff,
            'random_state': random_state,
            'bodyname': meta['bodyname'],
            'param_names': list(param_names),
        },
    }


def plot_k2_posteriors(
    reanalysis: Dict[str, Any],
    output_path: str,
    observation: Optional[Tuple[float, float]] = None,
    title_prefix: str = '',
) -> None:
    """Quick 4-panel histogram of (Re_k2, Im_k2, Re_h2, Im_h2) posteriors.

    Parameters
    ----------
    reanalysis
        Dict returned by ``reanalyze_k2_from_pickle``.
    output_path
        Absolute path to save the PNG (dpi=150).
    observation
        Optional ``(value, sigma)`` for Re(k₂) — drawn as a vertical band on
        the Re(k₂) panel only.
    title_prefix
        Prepended to the figure suptitle (e.g. ``'Titan Test48: '``).
    """
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    components = [
        ('Re_k2', r'$\mathrm{Re}(k_2)$', axes[0, 0]),
        ('Im_k2', r'$\mathrm{Im}(k_2)$', axes[0, 1]),
        ('Re_h2', r'$\mathrm{Re}(h_2)$', axes[1, 0]),
        ('Im_h2', r'$\mathrm{Im}(h_2)$', axes[1, 1]),
    ]
    for key, label, ax in components:
        arr = reanalysis[key]
        finite = arr[np.isfinite(arr)]
        if finite.size == 0:
            ax.text(0.5, 0.5, 'all NaN', ha='center', va='center',
                    transform=ax.transAxes)
        else:
            ax.hist(finite, bins=40, color='#1E90FF', alpha=0.75, edgecolor='#333')
            med = np.median(finite)
            ax.axvline(med, color='k', lw=1.2, ls='--',
                       label=f'median={med:.3g}')
            ax.legend(fontsize=8, loc='best')
        ax.set_xlabel(label)
        ax.set_ylabel('count')

    if observation is not None:
        v, s = observation
        ax = axes[0, 0]
        ax.axvspan(v - s, v + s, color='red', alpha=0.18, label='1σ obs')
        ax.axvline(v, color='red', lw=1.0)

    md = reanalysis['metadata']
    fig.suptitle(
        f"{title_prefix}k₂/h₂ post-hoc re-analysis "
        f"(n_eval={md['n_eval']}, rheo={md['rheology_type']}, "
        f"failed={reanalysis['n_failed']})",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    log.info(f"Saved {output_path}")
    plt.close(fig)
