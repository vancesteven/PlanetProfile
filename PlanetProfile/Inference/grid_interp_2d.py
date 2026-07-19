"""Shared 2D (Tb × log₁₀ w) bilinear interpolant for the Europa Clipper v3
salinity-sampled structure cache.

The v3 cache is a regular grid in ``Tb_K`` × ``wOcean_ppt`` (schema ``v3.0``,
built by :func:`cache_builder.build_tbw_grid_cache`), row-major
``structures[i_Tb * n_w + i_w]``. Salinity is sampled as ``log10_wOcean_ppt``,
so interpolation is bilinear in (Tb, **log₁₀ w**) — the natural axis for a
Jeffreys-prior scale parameter and the axis on which the log-spaced grid is
uniform.

**Why one shared helper (scientific-reviewer, binding).** The MCMC likelihood
and the SBI support guard MUST evaluate the identical interpolant. If they
diverge, the training-set support ≠ the reference-MCMC support and the 2D
Tb↔w degeneracy gate fails for a spurious reason. Every consumer therefore
computes its blend from the SAME ``bilinear_weights`` output:

- structural scalars (CMR2, D_ocean_km, …) and complex Ae are blended by the
  same 4 corner indices + weights;
- ``None`` corners (failed/ice-extinct/no-ocean builds — present in TWO
  opposite corners of the tilted valid band, low-Tb/low-w and
  high-Tb/high-w) trigger a nearest-valid fallback; if no corner is valid the
  sample is rejected.

This module is deliberately free of PlanetProfile imports so it is unit-tested
standalone.
"""
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import numpy as np


def wOcean_ppt_from_theta(theta_dict) -> float:
    """Salinity in ppt from a theta dict, exponentiating ``log10_wOcean_ppt``
    EXACTLY once.

    v3 samples salinity as ``log10_wOcean_ppt`` (Jeffreys prior). Every 2D
    lookup consumer (structural blend in forward_models, Ae blend in
    mcmc_runner) routes through here so the ``10**θ`` conversion lives in one
    place — the scientific-reviewer flagged double-log as a correctness risk.
    A plain ``wOcean_ppt`` key, if present, is returned as-is.
    """
    if 'log10_wOcean_ppt' in theta_dict:
        return float(10.0 ** float(theta_dict['log10_wOcean_ppt']))
    if 'wOcean_ppt' in theta_dict:
        return float(theta_dict['wOcean_ppt'])
    raise KeyError(
        "2D (Tb × w) cache requires 'log10_wOcean_ppt' (or 'wOcean_ppt') in "
        "theta_dict; neither found."
    )


def is_2d_cache(structure_data) -> bool:
    """True if ``structure_data`` is a v3.0 2D (Tb × w) cache.

    The presence of ``wOcean_ppt_grid`` is the format discriminant: legacy
    1D v1/v2 caches (``Tb_K_grid`` + ``structures`` only) return False so they
    stay servable through the existing nearest/blend paths.
    """
    return (
        isinstance(structure_data, dict)
        and "wOcean_ppt_grid" in structure_data
        and "Tb_K_grid" in structure_data
        and "structures" in structure_data
    )


def bilinear_weights(
    Tb_grid: np.ndarray,
    w_grid: np.ndarray,
    Tb_sample: float,
    w_sample_ppt: float,
) -> Tuple[List[int], List[float]]:
    """Bilinear corner flat-indices and weights for a (Tb, w) query.

    Interpolates in (Tb, **log₁₀ w**). The sample is clamped to the grid box
    (no extrapolation). Returns four row-major flat indices
    ``i_Tb * n_w + i_w`` for the bracketing cell corners
    ``(lo,lo), (lo,hi), (hi,lo), (hi,hi)`` and their bilinear weights (summing
    to 1). Callers map indices through ``structures`` / ``_ae_grid_cache`` and
    apply :func:`resolve_none_corners` before blending.

    Parameters
    ----------
    Tb_grid, w_grid
        1-D ascending grid axes (K, ppt). ``w_grid`` is log-spaced; weighting
        is done in log₁₀ w.
    Tb_sample, w_sample_ppt
        Query point. ``w_sample_ppt`` is a salinity in ppt (i.e. the caller
        has already exponentiated ``10**log10_wOcean_ppt`` — do NOT pass the
        log value here; that would double-log).
    """
    Tb_grid = np.asarray(Tb_grid, dtype=float)
    w_grid = np.asarray(w_grid, dtype=float)
    n_Tb, n_w = Tb_grid.size, w_grid.size

    # --- Tb axis (linear) ---
    Tb_c = float(np.clip(Tb_sample, Tb_grid[0], Tb_grid[-1]))
    j = int(np.searchsorted(Tb_grid, Tb_c))
    if j <= 0:
        i_lo_T = i_hi_T = 0
        f_T = 0.0
    elif j >= n_Tb:
        i_lo_T = i_hi_T = n_Tb - 1
        f_T = 0.0
    else:
        i_lo_T, i_hi_T = j - 1, j
        t0, t1 = Tb_grid[i_lo_T], Tb_grid[i_hi_T]
        f_T = 0.0 if t1 == t0 else (Tb_c - t0) / (t1 - t0)

    # --- w axis (log₁₀) ---
    logw_grid = np.log10(w_grid)
    logw_c = float(np.clip(np.log10(w_sample_ppt), logw_grid[0], logw_grid[-1]))
    k = int(np.searchsorted(logw_grid, logw_c))
    if k <= 0:
        i_lo_w = i_hi_w = 0
        f_w = 0.0
    elif k >= n_w:
        i_lo_w = i_hi_w = n_w - 1
        f_w = 0.0
    else:
        i_lo_w, i_hi_w = k - 1, k
        w0, w1 = logw_grid[i_lo_w], logw_grid[i_hi_w]
        f_w = 0.0 if w1 == w0 else (logw_c - w0) / (w1 - w0)

    def flat(i_T, i_w):
        return i_T * n_w + i_w

    corners = [
        flat(i_lo_T, i_lo_w),
        flat(i_lo_T, i_hi_w),
        flat(i_hi_T, i_lo_w),
        flat(i_hi_T, i_hi_w),
    ]
    weights = [
        (1.0 - f_T) * (1.0 - f_w),
        (1.0 - f_T) * f_w,
        f_T * (1.0 - f_w),
        f_T * f_w,
    ]
    return corners, weights


def resolve_none_corners(
    corners: Sequence[int],
    weights: Sequence[float],
    valid: Sequence[bool],
) -> Optional[Tuple[List[int], List[float]]]:
    """Drop ``None`` corners and renormalize; return None if all invalid.

    ``valid[i]`` is True when ``structures[corners[i]]`` (or the corresponding
    ``_ae_grid_cache`` entry) is not None. Behavior:

    - all four valid → the corners/weights are returned unchanged (true
      bilinear);
    - some valid → those are kept and their weights renormalized to sum to 1
      (a nearest-valid partial blend — degrades gracefully toward the nearest
      built node as the query approaches a ``None`` corner of the tilted band);
    - none valid → returns ``None``, signalling the caller to REJECT the
      sample (the (Tb, w) query sits entirely in an unbuilt corner, e.g.
      ice-extinct at high w / thick-shell at low w).

    Keeping this policy in one place is what guarantees the likelihood and the
    support guard treat ``None`` corners identically.
    """
    kept_idx = [i for i in range(len(corners)) if valid[i]]
    if not kept_idx:
        return None
    kept_corners = [corners[i] for i in kept_idx]
    kept_w = np.array([weights[i] for i in kept_idx], dtype=float)
    total = kept_w.sum()
    if total <= 0:
        # Degenerate weights (query exactly on kept nodes with zero weight):
        # fall back to uniform over the kept corners.
        kept_w = np.ones(len(kept_idx), dtype=float)
        total = kept_w.sum()
    kept_w = kept_w / total
    return kept_corners, kept_w.tolist()


def blend_complex(
    values: Sequence[Optional[complex]],
    weights: Sequence[float],
) -> Optional[complex]:
    """Weighted blend of complex Ae values (Re and Im blended separately).

    ``values``/``weights`` are already ``None``-resolved and renormalized by
    :func:`resolve_none_corners` (i.e. no None entries remain and weights sum
    to 1). Returns None if any surviving value is None (defensive) — callers
    treat that as reject.
    """
    acc = 0.0 + 0.0j
    for v, wt in zip(values, weights):
        if v is None:
            return None
        acc += complex(v) * wt
    return acc


def blend_scalar(
    values: Sequence[float],
    weights: Sequence[float],
) -> float:
    """Weighted blend of a real structural scalar (e.g. CMR2, D_ocean_km)."""
    acc = 0.0
    for v, wt in zip(values, weights):
        acc += float(v) * wt
    return acc
