"""
Probe Tb_K prior bounds against the cached structure grid.

For each Tb_K grid point in a body's structure cache, evaluate the k2/h2
forward model with a small set of fiducial theta values drawn from the
configured prior. Reports per-Tb success rates and recommends a tightened
``[Tb_min, Tb_max]`` bound that excludes Tb grid points where the
TidalPy ``radial_solver`` integrator fails for most of the prior.

This is a generalizable diagnostic: when adding a new body or changing
the cache cell layout, run this probe before launching production MCMC
to surface integrator-failing prior regions instead of letting them
silently truncate the posterior at -1e30.

Usage
-----
    # Set NUMBA_CACHE_DIR before invoking — see run_inference_cli.py.
    NUMBA_CACHE_DIR=$TMPDIR/pp_numba_cache python -m PlanetProfile.Inference.probe_tb_prior \
        --config PlanetProfile/Inference/configs/europa_seawater_andrade_7D.json

    # Optional: more fiducials for tighter sensitivity check
    python -m PlanetProfile.Inference.probe_tb_prior \
        --config <cfg> --n-fiducials 9 --threshold 0.8

Output
------
A table of (Tb_K, n_finite/n_fiducials, success_rate) plus a recommended
``Tb_K`` prior bound that contains only grid points whose success rate
meets ``--threshold``.

Author: PlanetProfile Team
Date:   2026-05-22
"""
from __future__ import annotations

# Set NUMBA_CACHE_DIR before any TidalPy import in case this module is
# imported directly (i.e. not through run_inference_cli.py).
import os
import tempfile
if not os.environ.get('NUMBA_CACHE_DIR'):
    _default = os.path.join(tempfile.gettempdir(), 'pp_numba_cache')
    try:
        os.makedirs(_default, exist_ok=True)
        os.environ['NUMBA_CACHE_DIR'] = _default
    except OSError:
        pass

import argparse
import json
import logging
import pickle
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional

import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.forward_models import forward_model_k2_flexible

log = logging.getLogger('PlanetProfile.Inference.probe_tb_prior')


# ---------------------------------------------------------------------------
# Fiducial-theta sampling
# ---------------------------------------------------------------------------
def _prior_median(prior_cfg: Dict[str, Any]) -> float:
    """Return the median (or mean) of a prior spec — supports uniform/gaussian."""
    ptype = prior_cfg.get('prior_type', 'uniform').lower()
    if ptype == 'uniform':
        lo, hi = prior_cfg['bounds']
        return 0.5 * (float(lo) + float(hi))
    if ptype in ('normal', 'gaussian'):
        return float(prior_cfg.get('mean', 0.0))
    raise ValueError(f"Unsupported prior_type for probe: {ptype!r}")


def _prior_sample(prior_cfg: Dict[str, Any], rng: np.random.Generator) -> float:
    ptype = prior_cfg.get('prior_type', 'uniform').lower()
    if ptype == 'uniform':
        lo, hi = prior_cfg['bounds']
        return float(rng.uniform(float(lo), float(hi)))
    if ptype in ('normal', 'gaussian'):
        return float(rng.normal(float(prior_cfg.get('mean', 0.0)),
                                float(prior_cfg.get('std', 1.0))))
    raise ValueError(f"Unsupported prior_type for probe: {ptype!r}")


def build_fiducial_thetas(
    param_space: Dict[str, Dict[str, Any]],
    n_fiducials: int = 5,
    seed: int = 0,
    skip_param: str = 'Tb_K',
) -> List[Dict[str, float]]:
    """
    Build N fiducial theta dicts from the prior, fixing one parameter
    (default ``Tb_K``) to be filled in by the probe loop.

    First fiducial is the prior median (deterministic baseline). Remaining
    fiducials are drawn pseudorandomly with the supplied seed.
    """
    fiducials: List[Dict[str, float]] = []
    # First: prior medians
    median_theta = {p: _prior_median(spec) for p, spec in param_space.items()
                    if p != skip_param}
    fiducials.append(median_theta)

    # Remaining: random draws
    rng = np.random.default_rng(seed)
    for _ in range(max(0, n_fiducials - 1)):
        theta = {p: _prior_sample(spec, rng) for p, spec in param_space.items()
                 if p != skip_param}
        fiducials.append(theta)
    return fiducials


# ---------------------------------------------------------------------------
# Probe
# ---------------------------------------------------------------------------
def probe_cache(
    config: InferenceConfig,
    n_fiducials: int = 5,
    seed: int = 0,
    threshold: float = 0.8,
) -> Dict[str, Any]:
    """
    Probe each Tb_K grid point in the cache referenced by ``config`` with
    ``n_fiducials`` theta draws. Return a dict of per-Tb results and the
    recommended tightened Tb_K bounds.
    """
    cache_path = config.structure_cache_path
    if cache_path is None:
        raise ValueError("config.structure_cache_path is required for probing")
    cache_path = str(cache_path)
    if not Path(cache_path).exists():
        raise FileNotFoundError(f"Structure cache not found: {cache_path}")

    with open(cache_path, 'rb') as f:
        cache = pickle.load(f)
    Tb_grid = np.asarray(cache['Tb_K_grid'], dtype=np.float64)

    # Configured prior bounds (for context)
    Tb_prior_lo, Tb_prior_hi = config.param_space['Tb_K']['bounds']

    # Build fiducials, expand param_groups, attach fixed_params, ensure
    # any required core parameters are populated
    fiducials = build_fiducial_thetas(
        config.param_space, n_fiducials=n_fiducials, seed=seed, skip_param='Tb_K'
    )

    # Each fiducial gets group-key expansion + fixed_params. The inner
    # loop just appends Tb_K and calls forward_model_k2_flexible.
    def _expanded(theta: Dict[str, float]) -> Dict[str, float]:
        out = dict(theta)
        for group_key, members in config.param_groups.items():
            if group_key in out:
                for m in members:
                    out[m] = out[group_key]
        for k, v in config.fixed_params.items():
            out[k] = v
        return out

    arrhenius_params = config.sampler_settings.get('arrhenius_params')

    rows: List[Dict[str, Any]] = []
    for Tb in Tb_grid:
        finite_count = 0
        per_fid_results = []
        for theta in fiducials:
            full = _expanded(theta)
            full['Tb_K'] = float(Tb)
            try:
                Re_k2, Im_k2, Re_h2, Im_h2, _ = forward_model_k2_flexible(
                    full, cache, return_heating=False,
                    arrhenius_params=arrhenius_params,
                )
            except Exception as exc:  # noqa: BLE001
                log.debug(f"Tb={Tb}: forward call raised: {exc}")
                Re_k2 = np.nan
            ok = bool(np.isfinite(Re_k2))
            per_fid_results.append(ok)
            finite_count += int(ok)
        rows.append({
            'Tb_K': float(Tb),
            'n_finite': finite_count,
            'n_total': len(fiducials),
            'success_rate': finite_count / len(fiducials),
            'per_fiducial_finite': per_fid_results,
        })

    # Recommend a contiguous bound containing grid points with
    # success_rate >= threshold. If none qualify, return None for both.
    qualifying_idx = [i for i, r in enumerate(rows) if r['success_rate'] >= threshold]
    if qualifying_idx:
        # Pick the largest contiguous run of qualifying indices.
        runs: List[List[int]] = []
        cur: List[int] = [qualifying_idx[0]]
        for idx in qualifying_idx[1:]:
            if idx == cur[-1] + 1:
                cur.append(idx)
            else:
                runs.append(cur)
                cur = [idx]
        runs.append(cur)
        best = max(runs, key=len)
        Tb_min = rows[best[0]]['Tb_K']
        Tb_max = rows[best[-1]]['Tb_K']
        recommendation = {
            'Tb_K_min': Tb_min, 'Tb_K_max': Tb_max,
            'n_qualifying_grid_points': len(best),
            'n_total_grid_points': len(rows),
            'threshold': threshold,
        }
    else:
        recommendation = {
            'Tb_K_min': None, 'Tb_K_max': None,
            'n_qualifying_grid_points': 0,
            'n_total_grid_points': len(rows),
            'threshold': threshold,
            'note': ('No Tb grid point reached the success threshold. '
                     'The integrator failure spans the full sampled '
                     'range — Tb tightening alone will not unblock '
                     'this body.'),
        }

    return {
        'bodyname': config.bodyname,
        'cache_path': cache_path,
        'configured_Tb_bounds': [float(Tb_prior_lo), float(Tb_prior_hi)],
        'n_fiducials': len(fiducials),
        'rows': rows,
        'recommendation': recommendation,
    }


def format_report(result: Dict[str, Any]) -> str:
    lines = []
    lines.append(f"Body: {result['bodyname']}")
    lines.append(f"Cache: {result['cache_path']}")
    lines.append(f"Configured Tb_K bounds: {result['configured_Tb_bounds']}")
    lines.append(f"Fiducials per Tb: {result['n_fiducials']}")
    lines.append('')
    lines.append(f"{'Tb_K':>10} {'n_finite':>9} {'rate':>6}")
    lines.append('-' * 30)
    for r in result['rows']:
        lines.append(
            f"{r['Tb_K']:>10.4f} {r['n_finite']:>4}/{r['n_total']:<3} "
            f"{r['success_rate']:>6.2%}"
        )
    lines.append('')
    rec = result['recommendation']
    if rec['Tb_K_min'] is not None:
        lines.append(
            f"RECOMMENDED Tb_K bound: [{rec['Tb_K_min']:.4f}, {rec['Tb_K_max']:.4f}]  "
            f"({rec['n_qualifying_grid_points']}/{rec['n_total_grid_points']} grid "
            f"points at ≥{rec['threshold']:.0%} success)"
        )
    else:
        lines.append("RECOMMENDED Tb_K bound: NONE — see note:")
        lines.append(f"  {rec.get('note','')}")
    return '\n'.join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--config', required=True,
                        help='Path to InferenceConfig JSON file')
    parser.add_argument('--n-fiducials', type=int, default=5,
                        help='Number of theta draws per Tb grid point '
                             '(default 5; first is the prior median)')
    parser.add_argument('--seed', type=int, default=0,
                        help='RNG seed for fiducial draws (default 0)')
    parser.add_argument('--threshold', type=float, default=0.8,
                        help='Success-rate threshold for recommended bound '
                             '(default 0.80)')
    parser.add_argument('--json-out',
                        help='Write the full result dict as JSON to this path')
    args = parser.parse_args()

    logging.basicConfig(level=logging.WARNING,
                        format='%(asctime)s [%(levelname)s] %(message)s')

    config = InferenceConfig.from_json(args.config)
    result = probe_cache(config,
                         n_fiducials=args.n_fiducials,
                         seed=args.seed,
                         threshold=args.threshold)
    print(format_report(result))

    if args.json_out:
        with open(args.json_out, 'w') as f:
            json.dump(result, f, indent=2, default=float)
        print(f"\nFull JSON written to {args.json_out}")

    return 0


if __name__ == '__main__':
    raise SystemExit(_main())
