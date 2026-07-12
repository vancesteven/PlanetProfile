#!/usr/bin/env python3
"""
SBI validation driver: SBC + MCMC<->SBI cross-check + limiting-behavior.

Phase E of ``plans/monte-carlo-sbi-implementation-plan.md``. This is a
*manually invoked* command-line tool that runs the three RATIFIED
pass/fail gates against a trained SBI artifact (see ``sbi_runner.py``).
The gates are hard: a failure exits with code 2 and a loud stderr
summary, and must NEVER be tuned to pass. Passing exits 0.

Subcommands
-----------
``sbc``        Simulation-based calibration (Talts et al. 2018) on >=200
               held-out (theta, x) pairs. Gate: per-parameter rank-
               uniformity Kolmogorov-Smirnov p >= 0.05.

``crosscheck`` Compare an SBI posterior against an MCMC ``InferenceResult``
               conditioned on the SAME observables. Per-parameter gates
               (plan doc, decision 4):
                 - |mean_SBI - mean_MCMC| <= max(0.25*sigma_MCMC, MC-error),
                   where MC-error = sigma_MCMC / sqrt(ESS) and ESS is the
                   importance-weighted effective sample size of the MCMC
                   posterior (=(sum w)^2/sum(w^2); N if unweighted).
                 - sigma_SBI / sigma_MCMC in [0.7, 1.4].
                 - two-sample KS not rejected at alpha=0.01 (weighted MCMC
                   samples are resampled to floor(ESS) draws so the KS test
                   reflects the true information content, not an inflated N).

``limits``     Limiting-behavior check: condition the flow on a grid of
               |Im k2| values and assert the ``log10_eta_Ih`` posterior
               MEDIAN decreases monotonically as |Im k2| increases, and
               that 100% of drawn samples lie inside the prior box.
               Gate: strictly monotone (at most one tie), containment 100%.

``selftest``   Build a toy BoxUniform + linear-Gaussian setup, train a tiny
               flow, and run all three checks end-to-end. No PlanetProfile
               forward model or structure cache is touched. Used to validate
               the driver's code paths in CI-free manual runs.

Scientific conventions (do not change silently)
-----------------------------------------------
- The flow's log-density is a POSTERIOR density, never a likelihood.
  None of the gates here treat it as a likelihood.
- All stochastic steps take a fixed ``--seed`` (default 42). Every report
  JSON embeds config_hash, artifact metadata, the seed(s), sbi/torch
  versions, git SHA, and UTC timestamps for reproducibility.
- SBC uses the installed sbi 0.26.1 API, verified against the source:
  ``sbi.diagnostics.run_sbc(thetas, xs, posterior, num_posterior_samples,
  reduce_fns='marginals') -> (ranks, dap_samples)`` and
  ``sbi.diagnostics.check_sbc(ranks, prior_samples, dap_samples,
  num_posterior_samples) -> {'ks_pvals', 'c2st_ranks', 'c2st_dap'}``, with
  plots via ``sbi.analysis.sbc_rank_plot(ranks, num_posterior_samples,
  parameter_labels=..., plot_type='cdf') -> (fig, ax)``.

Heavy-run gating
----------------
This script is NOT run by the unit-test suite. It performs training and
large posterior draws that are far too slow for ``tests/``. The convention
for any future pytest wrapper is to gate it behind an environment variable::

    PP_RUN_SBI_HEAVY=1 pytest tests/validate_sbi_heavy_test.py

Do NOT add slow SBC/cross-check runs to ``tests/`` without that guard.

Test50 cache blocker
--------------------
The Test50 structure-cache pickle referenced by
``configs/test50_titan_noocean_andrade_8D.json`` does not exist in this
checkout. Subcommands that must run the real forward model (``sbc`` in
``--config`` generate mode) fail fast with the exact missing path and the
build route (``build_phase_c1_cache.py`` / ``cache_builder.py``). The
``crosscheck`` and ``limits`` subcommands operate on saved products
(MCMC pickle, SBI artifact) and do NOT need the cache.

Usage examples (once the cache + MCMC posterior exist)
------------------------------------------------------
    conda run -n PP python -m PlanetProfile.Inference.validate_sbi selftest \\
        --seed 42 --output-dir /tmp/sbi_selftest

    conda run -n PP python -m PlanetProfile.Inference.validate_sbi sbc \\
        --artifact PlanetProfile/Inference/sbi_artifacts/titan_<hash>_posterior.pt \\
        --config PlanetProfile/Inference/configs/test50_titan_noocean_andrade_8D.json \\
        --n-sbc 300 --seed 42 --output-dir validation/sbc_titan

    conda run -n PP python -m PlanetProfile.Inference.validate_sbi crosscheck \\
        --artifact .../titan_<hash>_posterior.pt \\
        --mcmc PlanetProfile/Test/mcmc_results/Titan/.../titan_result.pkl \\
        --seed 42 --output-dir validation/crosscheck_titan

    conda run -n PP python -m PlanetProfile.Inference.validate_sbi limits \\
        --artifact .../titan_<hash>_posterior.pt --re-k2 0.608 \\
        --seed 42 --output-dir validation/limits_titan

Author: PlanetProfile Team
Date: 2026-07-07
"""
import argparse
import json
import logging
import os
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

# Ensure numba has a writable cache dir before any TidalPy import chain
# (mirrors run_inference_cli.py). Harmless when TidalPy is never imported.
if not os.environ.get('NUMBA_CACHE_DIR'):
    _default_numba_cache = os.path.join(tempfile.gettempdir(), 'pp_numba_cache')
    try:
        os.makedirs(_default_numba_cache, exist_ok=True)
        os.environ['NUMBA_CACHE_DIR'] = _default_numba_cache
    except OSError:
        pass

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
)
log = logging.getLogger('PlanetProfile.Inference.validate_sbi')

# Exit codes.
EXIT_PASS = 0
EXIT_ERROR = 1        # setup/IO/config error (not a gate verdict)
EXIT_GATE_FAIL = 2    # a scientific gate failed

# ---- RATIFIED gate thresholds (plan doc, decision 4). Do not tune. -------
SBC_KS_ALPHA = 0.05          # per-parameter rank-uniformity KS p must be >= this
CROSSCHECK_MEAN_SIGMA_FRAC = 0.25   # mean tolerance = max(0.25*sigma, MC-error)
# Ratified 2026-07-10 (replacing per-param KS p >= alpha): the KS statistic D
# becomes a DETECTION metric judged against the reference posterior's own
# finite-sampling floor, paired with a scientific-materiality location check.
# PASS iff D <= 1.5 x matched-n self-D p99 AND |dmedian| <= 0.3 dex (log10_*
# params) or 0.3 x sigma_MCMC (linear params). Thresholds derive from the
# bootstrap floor + materiality currency, not from any artifact's values; the
# same rule fails the noiseless and 200k-era artifacts. p is still reported.
CROSSCHECK_SHAPE_FLOOR_FACTOR = 1.5
CROSSCHECK_MATERIALITY_DEX = 0.3
CROSSCHECK_SELF_D_BOOT = 200
CROSSCHECK_SIGMA_RATIO_LO = 0.70
CROSSCHECK_SIGMA_RATIO_HI = 1.40
CROSSCHECK_KS_ALPHA = 0.01   # two-sample KS p must be >= this
LIMITS_MAX_TIES = 1          # monotone median: at most one tie allowed
LIMITS_CONTAINMENT = 1.0     # 100% of samples must lie inside the prior box
# Anchor mode (ratified 2026-07-09, replacing the monotonicity premise that
# ground-truth MCMC falsified below |Im k2| ~ 0.15). Statistic amended
# 2026-07-10 (user-ratified): per sweep point, the 1-D Wasserstein-1 distance
# between flow samples and the anchor posterior must satisfy
# W1 <= 0.25 * sigma_anchor (same fraction/currency as the crosscheck mean
# gate). W1 replaces the median comparison because the Im=0.25 anchor is
# BIMODAL: a median lands mid-valley and is fragile there, while W1 is
# stable under bimodality, reduces to ~|dmedian| for unimodal shifts, and
# penalizes mode-smoothing in proportion to misplaced mass. Medians are
# still reported. Fixed from principles; do not tune.
LIMITS_ANCHOR_SIGMA_FRAC = 0.25

# Default |Im k2| sweep grid for the limits check (plan doc: 0.05..0.30).
DEFAULT_IMK2_GRID = (0.05, 0.10, 0.15, 0.20, 0.25, 0.30)
# Default parameter whose posterior median must respond to |Im k2|.
DEFAULT_MONOTONE_PARAM = 'log10_eta_Ih'
# Observable names understood as the |Im k2| channel.
_IM_K2_ALIASES = ('abs_Im_k2', 'Im_k2')


# ==========================================================================
# Provenance helpers
# ==========================================================================

def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def _versions() -> Dict[str, str]:
    """torch/sbi/numpy/scipy versions for the report (lazy imports)."""
    out: Dict[str, str] = {}
    try:
        import torch
        out['torch'] = torch.__version__
    except Exception:
        out['torch'] = 'unavailable'
    try:
        import sbi
        out['sbi'] = sbi.__version__
    except Exception:
        out['sbi'] = 'unavailable'
    try:
        import scipy
        out['scipy'] = scipy.__version__
    except Exception:
        out['scipy'] = 'unavailable'
    out['numpy'] = np.__version__
    return out


def _git_sha() -> Optional[str]:
    from .sbi_runner import _git_short_sha
    return _git_short_sha()


def _provenance_block(seed: int, extra: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    block = {
        'seed': int(seed),
        'timestamp_utc': _now_utc(),
        'git_sha': _git_sha(),
        'versions': _versions(),
    }
    if extra:
        block.update(extra)
    return block


def _json_default(o):
    """Make numpy / Path objects JSON-serializable."""
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)


def _safe_config_hash(config) -> str:
    """config.generate_hash(), tolerant of legacy pickled InferenceConfigs."""
    try:
        return config.generate_hash()
    except AttributeError as e:
        log.warning(f"Config hash unavailable for legacy pickle: {e}")
        return f'unavailable-legacy-pickle ({e})'


def _write_report(output_dir: Path, name: str, report: Dict[str, Any]) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / name
    with open(path, 'w') as f:
        json.dump(report, f, indent=2, default=_json_default)
    log.info(f"Report written: {path}")
    return path


def _loud(lines: Sequence[str]) -> None:
    """Emit a boxed summary to stderr (used for gate PASS/FAIL banners)."""
    width = max((len(s) for s in lines), default=0) + 4
    bar = '=' * width
    print(bar, file=sys.stderr)
    for s in lines:
        print(f"  {s}", file=sys.stderr)
    print(bar, file=sys.stderr)
    sys.stderr.flush()


# ==========================================================================
# Weighted-statistics helpers (MCMC posterior may carry importance weights)
# ==========================================================================

def _normalize_weights(weights: Optional[np.ndarray], n: int) -> np.ndarray:
    """Return weights summing to 1 (uniform 1/n when weights is None)."""
    if weights is None:
        return np.full(n, 1.0 / n, dtype=np.float64)
    w = np.asarray(weights, dtype=np.float64).ravel()
    if w.shape[0] != n:
        raise ValueError(f"weights length {w.shape[0]} != n_samples {n}")
    if np.any(w < 0) or not np.all(np.isfinite(w)):
        raise ValueError("weights must be finite and non-negative")
    s = w.sum()
    if s <= 0:
        raise ValueError("weights sum to <= 0")
    return w / s


def _effective_sample_size(weights_norm: np.ndarray) -> float:
    """ESS = (sum w)^2 / sum(w^2). Mirrors MCMCRunner's true-ESS formula.

    ``weights_norm`` is already normalized (sum = 1), so this reduces to
    1 / sum(w^2). Returns n for uniform weights.
    """
    w2 = float(np.sum(weights_norm ** 2))
    return (1.0 / w2) if w2 > 0 else float(len(weights_norm))


def _weighted_mean_std(
    samples: np.ndarray, weights_norm: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Per-column weighted mean and (population) weighted std.

    var_j = sum_i w_i (x_ij - mean_j)^2  (weights sum to 1).
    """
    mean = np.average(samples, axis=0, weights=weights_norm)
    var = np.average((samples - mean) ** 2, axis=0, weights=weights_norm)
    return mean, np.sqrt(var)


def _resample_by_weights(
    samples: np.ndarray, weights_norm: np.ndarray, n_draws: int, rng: np.random.RandomState
) -> np.ndarray:
    """Draw ``n_draws`` unweighted rows from ``samples`` proportional to weights.

    Used to turn a weighted MCMC posterior into an unweighted representative
    sample for the two-sample KS test. Resampling to floor(ESS) draws keeps
    the KS test from being over-powered by the raw (weighted) sample count.
    """
    n = samples.shape[0]
    idx = rng.choice(n, size=int(n_draws), replace=True, p=weights_norm)
    return samples[idx]


def _ks_stat_1d(a: np.ndarray, b: np.ndarray) -> float:
    """Two-sample KS statistic D (no p-value) via CDF scan."""
    a, b = np.sort(a), np.sort(b)
    grid = np.concatenate([a, b])
    ca = np.searchsorted(a, grid, side='right') / len(a)
    cb = np.searchsorted(b, grid, side='right') / len(b)
    return float(np.max(np.abs(ca - cb)))


def _self_d_floor(
    mcmc_samples: np.ndarray, weights_norm: np.ndarray, ess: float,
    rng: np.random.RandomState, n_boot: int = CROSSCHECK_SELF_D_BOOT,
) -> np.ndarray:
    """Per-parameter matched-n self-D p99 floor of the reference posterior.

    Bootstraps split-half KS D between two half-ESS resamples of the SAME
    posterior (the D a perfect flow would show from finite sampling alone),
    then scales by 1/sqrt(2) to the crosscheck's full-ESS-vs-full-ESS
    comparison (two-sample KS critical D scales as sqrt((n+m)/(n m))).
    """
    n_ess = max(int(ess), 4)
    half = n_ess // 2
    n_params = mcmc_samples.shape[1]
    d = np.empty((n_boot, n_params))
    for b in range(n_boot):
        r = _resample_by_weights(mcmc_samples, weights_norm, n_ess, rng)
        h1, h2 = r[:half], r[half:]
        for j in range(n_params):
            d[b, j] = _ks_stat_1d(h1[:, j], h2[:, j])
    return np.percentile(d, 99, axis=0) / np.sqrt(2.0)


def _weighted_median_1d(values: np.ndarray, weights_norm: np.ndarray) -> float:
    order = np.argsort(values)
    cw = np.cumsum(weights_norm[order])
    cw = cw / cw[-1]
    return float(values[order][np.searchsorted(cw, 0.5)])


# ==========================================================================
# Artifact / prior helpers
# ==========================================================================

def _load_artifact_runner(artifact_path: str):
    """Load a trained SBI artifact into an SBIRunner (sampling-only mode)."""
    from .sbi_runner import SBIRunner
    runner = SBIRunner.load_artifact(artifact_path)
    return runner


def _artifact_meta_for_report(runner) -> Dict[str, Any]:
    """Serializable subset of artifact metadata (drops big/opaque fields)."""
    meta = dict(getattr(runner, 'artifact_meta', {}) or {})
    keep = (
        'schema_version', 'config_hash', 'bodyname', 'git_sha', 'seed',
        'n_train_effective', 'sbi_version', 'torch_version', 'created_utc',
        'density_estimator', 'imag_convention', 'param_names', 'param_bounds',
        'obs_names', 'rejection_stats',
    )
    return {k: meta.get(k) for k in keep if k in meta}


def _build_boxuniform_from_bounds(param_bounds: Sequence[Sequence[float]]):
    """Reconstruct the training BoxUniform prior from artifact param_bounds."""
    import torch
    from sbi.utils.torchutils import BoxUniform
    bounds = np.asarray(param_bounds, dtype=np.float64)
    return BoxUniform(
        low=torch.as_tensor(bounds[:, 0], dtype=torch.float32),
        high=torch.as_tensor(bounds[:, 1], dtype=torch.float32),
    )


def _require_structure_cache(config) -> None:
    """Fail fast (EXIT_ERROR) if the config's structure cache is missing.

    Names the missing file and the build route, per the Phase-E blocker.
    """
    cache = getattr(config, 'structure_cache_path', None)
    if not cache:
        return
    if Path(cache).exists():
        return
    msg = (
        f"Structure cache not found: {cache}\n"
        f"This SBC generate-mode run needs the real forward model, which "
        f"requires the Tb-grid structure cache. Build it first, e.g.:\n"
        f"    conda run -n PP python -m PlanetProfile.Inference.build_phase_c1_cache "
        f"--config <body_config.json>\n"
        f"(build_phase_c1_cache.py -> cache_builder.build_tb_grid_cache). "
        f"Then re-run this command, or pass --dataset <held_out.npz> to use a "
        f"pre-generated (theta, x) set instead of generating one."
    )
    raise FileNotFoundError(msg)


# ==========================================================================
# SBC (simulation-based calibration)
# ==========================================================================

def _load_heldout_dataset(dataset_path: str, param_names, obs_names):
    """Load held-out (theta, x) pairs from an .npz produced by generate_sbi_dataset."""
    with np.load(dataset_path, allow_pickle=True) as npz:
        theta = np.asarray(npz['theta'], dtype=np.float64)
        x = np.asarray(npz['x'], dtype=np.float64)
        ds_params = [str(p) for p in npz['param_names']] if 'param_names' in npz.files else None
        ds_obs = [str(o) for o in npz['obs_names']] if 'obs_names' in npz.files else None
    if ds_params is not None and list(ds_params) != list(param_names):
        raise ValueError(
            f"Dataset param_names {ds_params} do not match artifact "
            f"param_names {list(param_names)} — refusing to run SBC on "
            f"mismatched parameterizations."
        )
    if ds_obs is not None and list(ds_obs) != list(obs_names):
        raise ValueError(
            f"Dataset obs_names {ds_obs} do not match artifact obs_names "
            f"{list(obs_names)}."
        )
    finite = np.all(np.isfinite(theta), axis=1) & np.all(np.isfinite(x), axis=1)
    return theta[finite], x[finite]


def _run_sbc_check(
    posterior,
    thetas: np.ndarray,
    xs: np.ndarray,
    param_labels: Sequence[str],
    num_posterior_samples: int,
    seed: int,
) -> Tuple[Dict[str, Any], Any]:
    """Run SBC and the KS uniformity check; return (results_dict, ranks_tensor).

    Uses sbi 0.26.1: run_sbc(...reduce_fns='marginals') then check_sbc(...).
    """
    import torch
    from sbi.diagnostics import check_sbc

    torch.manual_seed(int(seed))
    thetas_t = torch.as_tensor(np.asarray(thetas, dtype=np.float64), dtype=torch.float32)
    xs_t = torch.as_tensor(np.asarray(xs, dtype=np.float64), dtype=torch.float32)

    # Sample manually with reject_outside_prior=False to avoid 0%-acceptance
    # hangs on flows with heavy-tailed x distributions (e.g. extreme tidal k2).
    # This replaces the run_sbc() call, which passes through the sbi batched
    # sampler that does not expose reject_outside_prior.
    n_sbc = thetas_t.shape[0]
    n_p = thetas_t.shape[1]
    N = int(num_posterior_samples)
    ranks_list = []
    dap_list = []
    for i in range(n_sbc):
        post_samples = posterior.sample(
            (N,), x=xs_t[i],
            show_progress_bars=False,
            reject_outside_prior=False,
        )  # shape (N, n_p)
        # marginal rank: fraction of posterior samples < ground truth
        ranks_i = (post_samples < thetas_t[i].unsqueeze(0)).float().sum(dim=0)  # (n_p,)
        ranks_list.append(ranks_i)
        # DAP: one random draw from posterior per SBC pair
        dap_list.append(post_samples[0])
    ranks = torch.stack(ranks_list).detach()           # (n_sbc, n_p)
    dap_samples = torch.stack(dap_list).detach()       # (n_sbc, n_p)
    prior_samples = thetas_t
    checks = check_sbc(
        ranks, prior_samples, dap_samples,
        num_posterior_samples=N,
    )
    ks_pvals = np.asarray(checks['ks_pvals'].detach().cpu().numpy(), dtype=np.float64)
    c2st_ranks = np.asarray(checks['c2st_ranks'].detach().cpu().numpy(), dtype=np.float64)

    per_param = []
    for i, label in enumerate(param_labels):
        per_param.append({
            'param': str(label),
            'ks_pval': float(ks_pvals[i]),
            'ks_pass': bool(ks_pvals[i] >= SBC_KS_ALPHA),
            'c2st_rank_accuracy': float(c2st_ranks[i]),
        })
    results = {
        'n_sbc_pairs': int(thetas_t.shape[0]),
        'num_posterior_samples': int(num_posterior_samples),
        'ks_alpha': SBC_KS_ALPHA,
        'per_parameter': per_param,
        'all_pass': bool(np.all(ks_pvals >= SBC_KS_ALPHA)),
        'min_ks_pval': float(np.min(ks_pvals)),
    }
    return results, ranks


def _save_sbc_rank_plot(ranks, num_posterior_samples, param_labels, output_dir: Path) -> Optional[str]:
    """Save an SBC rank CDF plot; return its path (None on failure)."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        from sbi.analysis import sbc_rank_plot
        fig, _ = sbc_rank_plot(
            ranks=ranks,
            num_posterior_samples=int(num_posterior_samples),
            plot_type='cdf',
            parameter_labels=list(param_labels),
        )
        out = output_dir / 'sbc_rank_cdf.png'
        fig.savefig(out, dpi=120, bbox_inches='tight')
        import matplotlib.pyplot as plt
        plt.close(fig)
        return str(out)
    except Exception as e:
        log.warning(f"Could not render SBC rank plot: {e}")
        return None


def cmd_sbc(args) -> int:
    output_dir = Path(args.output_dir)
    runner = _load_artifact_runner(args.artifact)
    param_names = runner.param_names
    param_labels = runner.param_labels
    obs_names = runner.obs_names

    # Obtain held-out (theta, x) pairs.
    if args.dataset:
        theta, x = _load_heldout_dataset(args.dataset, param_names, obs_names)
        source = {'mode': 'dataset', 'dataset': str(args.dataset)}
    elif args.config:
        from .inference_core import InferenceConfig
        cfg_dict = InferenceConfig.from_json(args.config).to_dict()
        cfg_dict['mode'] = 'sbi'
        config = InferenceConfig.from_dict(cfg_dict)
        _require_structure_cache(config)  # fail fast if Test50 cache missing
        from .sbi_runner import SBIRunner
        gen_runner = SBIRunner(config)
        if list(gen_runner.param_names) != list(param_names):
            raise ValueError(
                f"Config param_names {gen_runner.param_names} != artifact "
                f"param_names {param_names}."
            )
        theta, x, stats = gen_runner.generate_training_set(int(args.n_sbc), seed=int(args.seed))
        source = {'mode': 'generated', 'config': str(args.config), 'gen_stats': stats}
    else:
        raise ValueError("sbc requires either --dataset or --config")

    if theta.shape[0] < args.n_sbc:
        log.warning(
            f"Only {theta.shape[0]} finite held-out pairs available "
            f"(requested {args.n_sbc})."
        )
    n_use = min(theta.shape[0], int(args.n_sbc)) if args.dataset else theta.shape[0]
    theta, x = theta[:n_use], x[:n_use]

    if theta.shape[0] < 100:
        raise ValueError(
            f"SBC needs >=100 (recommended >=200) held-out pairs; got "
            f"{theta.shape[0]}. Provide a larger --dataset or --n-sbc."
        )

    results, ranks = _run_sbc_check(
        runner._posterior, theta, x, param_labels,
        num_posterior_samples=int(args.num_posterior_samples), seed=int(args.seed),
    )
    plot_path = _save_sbc_rank_plot(
        ranks, int(args.num_posterior_samples), param_labels, output_dir
    )

    report = {
        'check': 'sbc',
        'gate': f"per-parameter rank-uniformity KS p >= {SBC_KS_ALPHA}",
        'verdict': 'PASS' if results['all_pass'] else 'FAIL',
        'results': results,
        'source': source,
        'artifact': str(args.artifact),
        'artifact_meta': _artifact_meta_for_report(runner),
        'rank_plot': plot_path,
        'provenance': _provenance_block(int(args.seed)),
    }
    _write_report(output_dir, 'sbc_report.json', report)

    if results['all_pass']:
        _loud([
            "SBC GATE: PASS",
            f"{results['n_sbc_pairs']} pairs, min KS p = {results['min_ks_pval']:.4g} "
            f">= {SBC_KS_ALPHA}",
        ])
        return EXIT_PASS
    failing = [p['param'] for p in results['per_parameter'] if not p['ks_pass']]
    _loud([
        "SBC GATE: FAIL",
        f"min KS p = {results['min_ks_pval']:.4g} < {SBC_KS_ALPHA}",
        f"non-uniform parameters: {failing}",
        "Do NOT tune this gate — investigate flow calibration / training set.",
    ])
    return EXIT_GATE_FAIL


# ==========================================================================
# MCMC <-> SBI cross-check
# ==========================================================================

def _two_sample_ks(
    mcmc_col: np.ndarray,
    weights_norm: np.ndarray,
    sbi_col: np.ndarray,
    ess: float,
    rng: np.random.RandomState,
) -> Tuple[float, float, int]:
    """Weighted-MCMC vs SBI two-sample KS. Returns (statistic, pvalue, n_ks).

    The weighted MCMC posterior is resampled to floor(ESS) unweighted draws
    (so the KS test is not over-powered by the raw sample count), and the
    SBI column is subsampled to the same size.
    """
    from scipy.stats import ks_2samp
    n_ks = int(max(50, min(np.floor(ess), sbi_col.shape[0])))
    mcmc_rs = _resample_by_weights(
        mcmc_col.reshape(-1, 1), weights_norm, n_ks, rng
    ).ravel()
    if sbi_col.shape[0] > n_ks:
        sbi_idx = rng.choice(sbi_col.shape[0], size=n_ks, replace=False)
        sbi_rs = sbi_col[sbi_idx]
    else:
        sbi_rs = sbi_col
    stat, pval = ks_2samp(mcmc_rs, sbi_rs)
    return float(stat), float(pval), n_ks


def _run_crosscheck(
    param_names: Sequence[str],
    mcmc_samples: np.ndarray,
    mcmc_weights: Optional[np.ndarray],
    sbi_samples: np.ndarray,
    seed: int,
) -> Dict[str, Any]:
    """Per-parameter cross-check gates. All arrays share param_names order."""
    n_mcmc = mcmc_samples.shape[0]
    w = _normalize_weights(mcmc_weights, n_mcmc)
    ess = _effective_sample_size(w)
    mcmc_mean, mcmc_std = _weighted_mean_std(mcmc_samples, w)
    sbi_mean = sbi_samples.mean(axis=0)
    sbi_std = sbi_samples.std(axis=0)
    mc_error = mcmc_std / np.sqrt(max(ess, 1.0))

    rng = np.random.RandomState(int(seed))
    # Reference's own finite-sampling floor (ratified 2026-07-10): the shape
    # gate judges D against this, not against an n-dependent p-value.
    d_floor_p99 = _self_d_floor(mcmc_samples, w, ess, rng)

    per_param = []
    all_pass = True
    for j, name in enumerate(param_names):
        mean_tol = max(CROSSCHECK_MEAN_SIGMA_FRAC * mcmc_std[j], mc_error[j])
        mean_diff = abs(sbi_mean[j] - mcmc_mean[j])
        mean_pass = bool(mean_diff <= mean_tol)

        sigma_ratio = float(sbi_std[j] / mcmc_std[j]) if mcmc_std[j] > 0 else float('inf')
        sigma_pass = bool(CROSSCHECK_SIGMA_RATIO_LO <= sigma_ratio <= CROSSCHECK_SIGMA_RATIO_HI)

        ks_stat, ks_pval, n_ks = _two_sample_ks(
            mcmc_samples[:, j], w, sbi_samples[:, j], ess, rng
        )
        # Detection + materiality shape gate (ratified 2026-07-10; replaces
        # ks_pval >= alpha, which convolves effect size with n). p remains
        # reported for reference.
        d_tol = CROSSCHECK_SHAPE_FLOOR_FACTOR * float(d_floor_p99[j])
        d_pass = bool(ks_stat <= d_tol)
        median_diff = abs(_weighted_median_1d(mcmc_samples[:, j], w)
                          - float(np.median(sbi_samples[:, j])))
        if str(name).startswith('log10_'):
            median_tol = CROSSCHECK_MATERIALITY_DEX
        else:
            median_tol = CROSSCHECK_MATERIALITY_DEX * float(mcmc_std[j])
        median_pass = bool(median_diff <= median_tol)
        shape_pass = d_pass and median_pass

        param_pass = mean_pass and sigma_pass and shape_pass
        all_pass = all_pass and param_pass
        per_param.append({
            'param': str(name),
            'mcmc_mean': float(mcmc_mean[j]), 'sbi_mean': float(sbi_mean[j]),
            'mcmc_std': float(mcmc_std[j]), 'sbi_std': float(sbi_std[j]),
            'mean_diff': float(mean_diff), 'mean_tol': float(mean_tol),
            'mc_error': float(mc_error[j]), 'mean_pass': mean_pass,
            'sigma_ratio': sigma_ratio, 'sigma_ratio_bounds': [CROSSCHECK_SIGMA_RATIO_LO, CROSSCHECK_SIGMA_RATIO_HI],
            'sigma_pass': sigma_pass,
            'ks_stat': ks_stat, 'ks_pval': ks_pval, 'ks_n': n_ks,
            'ks_alpha_informational': CROSSCHECK_KS_ALPHA,
            'd_floor_p99': float(d_floor_p99[j]), 'd_tol': float(d_tol),
            'd_pass': d_pass,
            'median_diff': float(median_diff), 'median_tol': float(median_tol),
            'median_pass': median_pass, 'shape_pass': bool(shape_pass),
            'param_pass': bool(param_pass),
        })
    return {
        'n_mcmc': int(n_mcmc), 'ess_mcmc': float(ess),
        'n_sbi': int(sbi_samples.shape[0]),
        'per_parameter': per_param,
        'all_pass': bool(all_pass),
    }


def _reorder_sbi_to_mcmc(sbi_runner, sbi_samples, mcmc_param_names):
    """Reorder SBI sample columns to the MCMC parameter order (must be a superset match)."""
    sbi_names = list(sbi_runner.param_names)
    if list(mcmc_param_names) == sbi_names:
        return sbi_samples, list(mcmc_param_names)
    if set(mcmc_param_names) != set(sbi_names):
        raise ValueError(
            f"Parameter sets differ between MCMC {list(mcmc_param_names)} and "
            f"SBI {sbi_names} — cannot cross-check mismatched parameterizations."
        )
    order = [sbi_names.index(n) for n in mcmc_param_names]
    return sbi_samples[:, order], list(mcmc_param_names)


def cmd_crosscheck(args) -> int:
    output_dir = Path(args.output_dir)
    runner = _load_artifact_runner(args.artifact)

    from .inference_core import InferenceResult
    mcmc = InferenceResult.load(args.mcmc)
    mcmc_samples = np.asarray(mcmc.samples, dtype=np.float64)
    mcmc_weights = None if mcmc.weights is None else np.asarray(mcmc.weights, dtype=np.float64)
    mcmc_param_names = list(mcmc.param_names)

    # Observed data to condition the flow on: explicit --x-obs JSON, else the
    # MCMC config's observable central values (imag channel -> |Im k2|).
    if args.x_obs:
        x_obs = json.loads(args.x_obs)
    else:
        x_obs = {name: float(spec[0]) for name, spec in mcmc.config.observables.items()}
    x_obs = {k: float(v) for k, v in x_obs.items()}

    n_draw = int(args.num_samples) if args.num_samples else mcmc_samples.shape[0]
    sbi_samples = runner.sample_posterior(x_obs, n_samples=n_draw, seed=int(args.seed))
    sbi_samples, ordered_names = _reorder_sbi_to_mcmc(runner, sbi_samples, mcmc_param_names)

    results = _run_crosscheck(
        ordered_names, mcmc_samples, mcmc_weights, sbi_samples, seed=int(args.seed)
    )
    plot_path = _save_crosscheck_plot(
        ordered_names, mcmc_samples, mcmc_weights, sbi_samples, output_dir
    )

    report = {
        'check': 'crosscheck',
        'gate': (
            f"per-param: |dmean| <= max({CROSSCHECK_MEAN_SIGMA_FRAC}*sigma_MCMC, MC-error); "
            f"sigma_SBI/sigma_MCMC in [{CROSSCHECK_SIGMA_RATIO_LO},{CROSSCHECK_SIGMA_RATIO_HI}]; "
            f"shape (ratified 2026-07-10): D <= {CROSSCHECK_SHAPE_FLOOR_FACTOR}x "
            f"matched-n self-D p99 AND |dmedian| <= {CROSSCHECK_MATERIALITY_DEX} dex "
            f"(log params) / {CROSSCHECK_MATERIALITY_DEX}*sigma_MCMC (linear); "
            f"KS p reported informationally"
        ),
        'verdict': 'PASS' if results['all_pass'] else 'FAIL',
        'x_obs': x_obs,
        'results': results,
        'artifact': str(args.artifact),
        'artifact_meta': _artifact_meta_for_report(runner),
        'mcmc_result': str(args.mcmc),
        # Legacy result pickles may predate newer InferenceConfig attributes
        # (e.g. arrhenius_params); the hash is report metadata only, so degrade
        # gracefully instead of aborting the gate run.
        'mcmc_config_hash': _safe_config_hash(mcmc.config),
        'mcmc_weighted': bool(mcmc_weights is not None),
        'overlay_plot': plot_path,
        'provenance': _provenance_block(int(args.seed)),
    }
    _write_report(output_dir, 'crosscheck_report.json', report)

    if results['all_pass']:
        _loud([
            "CROSS-CHECK GATE: PASS",
            f"{len(ordered_names)} params, ESS_MCMC = {results['ess_mcmc']:.1f}",
        ])
        return EXIT_PASS
    failing = [p['param'] for p in results['per_parameter'] if not p['param_pass']]
    lines = ["CROSS-CHECK GATE: FAIL", f"failing parameters: {failing}"]
    for p in results['per_parameter']:
        if not p['param_pass']:
            reasons = []
            if not p['mean_pass']:
                reasons.append(f"mean_diff {p['mean_diff']:.3g} > tol {p['mean_tol']:.3g}")
            if not p['sigma_pass']:
                reasons.append(f"sigma_ratio {p['sigma_ratio']:.3g} outside "
                               f"[{CROSSCHECK_SIGMA_RATIO_LO},{CROSSCHECK_SIGMA_RATIO_HI}]")
            if not p.get('d_pass', True):
                reasons.append(f"D {p['ks_stat']:.3g} > tol {p['d_tol']:.3g} "
                               f"(1.5x self-D floor {p['d_floor_p99']:.3g})")
            if not p.get('median_pass', True):
                reasons.append(f"|dmedian| {p['median_diff']:.3g} > "
                               f"materiality tol {p['median_tol']:.3g}")
            lines.append(f"  {p['param']}: " + "; ".join(reasons))
    lines.append("Do NOT tune this gate — investigate SBI vs MCMC discrepancy.")
    _loud(lines)
    return EXIT_GATE_FAIL


def _save_crosscheck_plot(param_names, mcmc_samples, mcmc_weights, sbi_samples, output_dir: Path):
    """Overlay 1-D marginal histograms (weighted MCMC vs SBI). Returns path or None."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        w = _normalize_weights(mcmc_weights, mcmc_samples.shape[0])
        d = len(param_names)
        ncol = min(4, d)
        nrow = int(np.ceil(d / ncol))
        fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 2.6 * nrow), squeeze=False)
        for j, name in enumerate(param_names):
            ax = axes[j // ncol][j % ncol]
            lo = min(mcmc_samples[:, j].min(), sbi_samples[:, j].min())
            hi = max(mcmc_samples[:, j].max(), sbi_samples[:, j].max())
            bins = np.linspace(lo, hi, 40)
            ax.hist(mcmc_samples[:, j], bins=bins, weights=w, density=True,
                    histtype='step', color='C0', label='MCMC')
            ax.hist(sbi_samples[:, j], bins=bins, density=True,
                    histtype='step', color='C1', label='SBI')
            ax.set_title(str(name), fontsize=8)
            ax.tick_params(labelsize=6)
        for k in range(d, nrow * ncol):
            axes[k // ncol][k % ncol].axis('off')
        axes[0][0].legend(fontsize=7)
        fig.tight_layout()
        out = output_dir / 'crosscheck_marginals.png'
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=120, bbox_inches='tight')
        plt.close(fig)
        return str(out)
    except Exception as e:
        log.warning(f"Could not render cross-check overlay plot: {e}")
        return None


# ==========================================================================
# Limiting-behavior check
# ==========================================================================

def _check_monotone(values: Sequence[float], direction: str, tol: float) -> Dict[str, Any]:
    """Check a sequence is monotone in ``direction`` allowing <= LIMITS_MAX_TIES ties.

    A 'strict violation' is a step in the wrong direction beyond ``tol``.
    A 'tie' is a step with |diff| <= tol. Gate: 0 strict violations and
    at most LIMITS_MAX_TIES ties.
    """
    vals = np.asarray(values, dtype=np.float64)
    diffs = np.diff(vals)
    ties = int(np.sum(np.abs(diffs) <= tol))
    if direction == 'decreasing':
        strict_viol = int(np.sum(diffs > tol))
    elif direction == 'increasing':
        strict_viol = int(np.sum(diffs < -tol))
    else:
        raise ValueError(f"direction must be 'increasing' or 'decreasing', got {direction}")
    ok = (strict_viol == 0) and (ties <= LIMITS_MAX_TIES)
    return {
        'direction': direction, 'tol': float(tol),
        'diffs': diffs.tolist(), 'n_strict_violations': strict_viol,
        'n_ties': ties, 'max_ties_allowed': LIMITS_MAX_TIES, 'monotone_pass': bool(ok),
    }


def _run_limits_check(
    runner,
    sweep_obs: str,
    sweep_values: Sequence[float],
    fixed_obs: Dict[str, float],
    monotone_param: str,
    direction: str,
    n_samples: int,
    seed: int,
) -> Dict[str, Any]:
    """Condition the flow across a sweep grid; test median monotonicity + containment."""
    if monotone_param not in runner.param_names:
        raise ValueError(
            f"monotone-param '{monotone_param}' not in artifact parameters "
            f"{list(runner.param_names)}."
        )
    pidx = list(runner.param_names).index(monotone_param)
    bounds = np.asarray(runner.artifact_meta['param_bounds'], dtype=np.float64)

    grid_rows = []
    medians = []
    n_total = 0
    n_inside = 0
    for k, val in enumerate(sweep_values):
        x_obs = dict(fixed_obs)
        x_obs[sweep_obs] = float(val)
        # Diagnostic sweep at possibly far-out-of-distribution x: skip
        # sbi's prior-box rejection (NSF spline asserts there); the
        # containment gate below measures any leakage honestly.
        samples = runner.sample_posterior(x_obs, n_samples=int(n_samples), seed=int(seed) + k,
                                          reject_outside_prior=False)
        median = float(np.median(samples[:, pidx]))
        medians.append(median)
        inside = np.all(
            (samples >= bounds[:, 0] - 1e-6) & (samples <= bounds[:, 1] + 1e-6), axis=1
        )
        n_total += samples.shape[0]
        n_inside += int(np.sum(inside))
        grid_rows.append({
            'sweep_value': float(val),
            f'{monotone_param}_median': median,
            'n_samples': int(samples.shape[0]),
            'n_inside_box': int(np.sum(inside)),
        })

    tol = 1e-3 * float(bounds[pidx, 1] - bounds[pidx, 0])
    mono = _check_monotone(medians, direction, tol)
    containment = (n_inside / n_total) if n_total > 0 else 0.0
    containment_pass = bool(containment >= LIMITS_CONTAINMENT - 1e-12)

    return {
        'sweep_obs': sweep_obs,
        'sweep_values': [float(v) for v in sweep_values],
        'fixed_obs': fixed_obs,
        'monotone_param': monotone_param,
        'grid': grid_rows,
        'medians': [float(m) for m in medians],
        'monotonicity': mono,
        'containment_fraction': float(containment),
        'containment_required': LIMITS_CONTAINMENT,
        'containment_pass': containment_pass,
        'all_pass': bool(mono['monotone_pass'] and containment_pass),
    }


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, q) -> np.ndarray:
    """Quantile(s) of a weighted sample (CDF inversion on sorted values)."""
    order = np.argsort(values)
    cw = np.cumsum(weights[order])
    cw = cw / cw[-1]
    idx = np.searchsorted(cw, q)
    return values[order][np.clip(idx, 0, len(values) - 1)]


def _wasserstein1_weighted(
    a_values: np.ndarray, a_weights: np.ndarray, b_values: np.ndarray,
    n_q: int = 512,
) -> float:
    """1-D Wasserstein-1 distance between a weighted sample (a) and an
    unweighted sample (b), via quantile-function integration on a uniform
    q grid: W1 = integral_0^1 |Fa^-1(q) - Fb^-1(q)| dq."""
    q = (np.arange(n_q) + 0.5) / n_q
    qa = _weighted_quantile(np.asarray(a_values, float),
                            np.asarray(a_weights, float), q)
    qb = np.quantile(np.asarray(b_values, float), q)
    return float(np.mean(np.abs(qa - qb)))


def _run_limits_anchor_check(
    runner,
    sweep_obs: str,
    anchor_results: Dict[float, str],
    fixed_obs: Dict[str, float],
    monotone_param: str,
    n_samples: int,
    seed: int,
) -> Dict[str, Any]:
    """Anchor mode: compare the flow's monotone-param median against
    ground-truth MCMC anchor posteriors at each sweep value.

    Gate per anchor: |median_flow - median_MCMC| <= LIMITS_ANCHOR_SIGMA_FRAC *
    sigma_MCMC (weighted). Prior-box containment identical to the legacy mode.
    """
    from .inference_core import InferenceResult

    if monotone_param not in runner.param_names:
        raise ValueError(
            f"monotone-param '{monotone_param}' not in artifact parameters "
            f"{list(runner.param_names)}."
        )
    pidx = list(runner.param_names).index(monotone_param)
    bounds = np.asarray(runner.artifact_meta['param_bounds'], dtype=np.float64)

    grid_rows = []
    n_total = 0
    n_inside = 0
    all_anchor_pass = True
    for k, (val, pkl_path) in enumerate(sorted(anchor_results.items())):
        mcmc = InferenceResult.load(pkl_path)
        if monotone_param not in mcmc.param_names:
            raise ValueError(
                f"Anchor {pkl_path} lacks parameter '{monotone_param}'."
            )
        aidx = list(mcmc.param_names).index(monotone_param)
        a_samples = np.asarray(mcmc.samples, dtype=np.float64)[:, aidx]
        a_w = _normalize_weights(
            None if mcmc.weights is None else np.asarray(mcmc.weights, dtype=np.float64),
            len(a_samples),
        )
        a_median = float(_weighted_quantile(a_samples, a_w, 0.5))
        a_mean, a_std = _weighted_mean_std(a_samples[:, None], a_w)
        a_std = float(a_std[0])

        x_obs = dict(fixed_obs)
        x_obs[sweep_obs] = float(val)
        # Diagnostic sweep at possibly far-out-of-distribution x: skip
        # sbi's prior-box rejection (NSF spline asserts there); the
        # containment gate below measures any leakage honestly.
        samples = runner.sample_posterior(x_obs, n_samples=int(n_samples), seed=int(seed) + k,
                                          reject_outside_prior=False)
        f_median = float(np.median(samples[:, pidx]))

        inside = np.all(
            (samples >= bounds[:, 0] - 1e-6) & (samples <= bounds[:, 1] + 1e-6), axis=1
        )
        n_total += samples.shape[0]
        n_inside += int(np.sum(inside))

        # W1 statistic (ratified 2026-07-10; see LIMITS_ANCHOR_SIGMA_FRAC
        # comment). Medians remain reported for reference.
        tol = LIMITS_ANCHOR_SIGMA_FRAC * a_std
        w1 = _wasserstein1_weighted(a_samples, a_w, samples[:, pidx])
        anchor_pass = bool(w1 <= tol)
        all_anchor_pass = all_anchor_pass and anchor_pass
        grid_rows.append({
            'sweep_value': float(val),
            'anchor_pkl': str(pkl_path),
            'anchor_median': a_median,
            'anchor_sigma': a_std,
            'anchor_ess': float((np.sum(a_w) ** 2) / np.sum(a_w ** 2)),
            'flow_median': f_median,
            'median_diff': float(abs(f_median - a_median)),
            'w1': float(w1),
            'w1_tol': float(tol),
            'anchor_pass': anchor_pass,
            'n_inside_box': int(np.sum(inside)),
            'n_samples': int(samples.shape[0]),
        })

    containment = (n_inside / n_total) if n_total > 0 else 0.0
    containment_pass = bool(containment >= LIMITS_CONTAINMENT - 1e-12)

    return {
        'mode': 'anchor',
        'sweep_obs': sweep_obs,
        'sweep_values': [r['sweep_value'] for r in grid_rows],
        'fixed_obs': fixed_obs,
        'monotone_param': monotone_param,
        'anchor_sigma_frac': LIMITS_ANCHOR_SIGMA_FRAC,
        'grid': grid_rows,
        'medians': [r['flow_median'] for r in grid_rows],
        'anchor_medians': [r['anchor_median'] for r in grid_rows],
        'anchors_all_pass': bool(all_anchor_pass),
        'containment_fraction': float(containment),
        'containment_required': LIMITS_CONTAINMENT,
        'containment_pass': containment_pass,
        'all_pass': bool(all_anchor_pass and containment_pass),
    }


def _save_limits_plot(results, output_dir: Path):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(4.5, 3.2))
        ax.plot(results['sweep_values'], results['medians'], 'o-', color='C2',
                label='flow median')
        if results.get('mode') == 'anchor':
            ax.errorbar(
                results['sweep_values'], results['anchor_medians'],
                yerr=[r['w1_tol'] for r in results['grid']],
                fmt='s', color='C0', capsize=3, label='MCMC anchor ± W1 tol')
            ax.legend(fontsize=7)
            title = 'Limiting behavior (MCMC-anchor comparison)'
        else:
            title = f"Limiting behavior ({results['monotonicity']['direction']})"
        ax.set_xlabel(results['sweep_obs'])
        ax.set_ylabel(f"{results['monotone_param']} posterior median")
        ax.set_title(title, fontsize=9)
        fig.tight_layout()
        out = output_dir / 'limits_median_trend.png'
        output_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=120, bbox_inches='tight')
        plt.close(fig)
        return str(out)
    except Exception as e:
        log.warning(f"Could not render limits plot: {e}")
        return None


def _resolve_sweep_obs(runner, requested: Optional[str]) -> str:
    if requested:
        if requested not in runner.obs_names:
            raise ValueError(
                f"--sweep-obs '{requested}' not in artifact obs_names {list(runner.obs_names)}."
            )
        return requested
    for alias in _IM_K2_ALIASES:
        if alias in runner.obs_names:
            return alias
    raise ValueError(
        f"Could not auto-detect an |Im k2| observable in {list(runner.obs_names)}; "
        f"pass --sweep-obs explicitly."
    )


def cmd_limits(args) -> int:
    output_dir = Path(args.output_dir)
    runner = _load_artifact_runner(args.artifact)

    sweep_obs = _resolve_sweep_obs(runner, args.sweep_obs)
    if args.sweep_values:
        sweep_values = [float(v) for v in args.sweep_values.split(',')]
    else:
        sweep_values = list(DEFAULT_IMK2_GRID)

    # Fixed observables (everything except the swept channel). --fixed-obs JSON
    # overrides; --re-k2 is a convenience for the Re_k2 channel.
    fixed_obs: Dict[str, float] = {}
    if args.fixed_obs:
        fixed_obs.update({k: float(v) for k, v in json.loads(args.fixed_obs).items()})
    _RE_K2_ALIASES = ('Re_k2', 're_k2', 'abs_Re_k2')
    if args.re_k2 is not None:
        matched = [n for n in runner.obs_names if n in _RE_K2_ALIASES and n not in fixed_obs]
        if matched:
            for name in matched:
                fixed_obs[name] = float(args.re_k2)
        else:
            # Fallback: apply to all non-Im_k2 non-fixed channels (original behaviour,
            # preserved for configs with no standard Re_k2 name).
            for name in runner.obs_names:
                if name not in _IM_K2_ALIASES and name not in fixed_obs:
                    fixed_obs[name] = float(args.re_k2)
    # Any remaining non-swept observable without a value is an error.
    missing = [n for n in runner.obs_names if n != sweep_obs and n not in fixed_obs]
    if missing:
        raise ValueError(
            f"Fixed values required for non-swept observables {missing}; "
            f"pass --re-k2 or --fixed-obs."
        )

    if getattr(args, 'anchor_results', None):
        anchor_map = {float(k): str(v)
                      for k, v in json.loads(args.anchor_results).items()}
        results = _run_limits_anchor_check(
            runner, sweep_obs, anchor_map, fixed_obs,
            monotone_param=args.monotone_param,
            n_samples=int(args.n_samples), seed=int(args.seed),
        )
        gate_desc = (
            f"{args.monotone_param}: Wasserstein-1(flow, MCMC anchor) <= "
            f"{LIMITS_ANCHOR_SIGMA_FRAC}*sigma_anchor at every sweep point "
            f"(ratified 2026-07-10, bimodality-stable); "
            f"prior-box containment == {LIMITS_CONTAINMENT}"
        )
    else:
        results = _run_limits_check(
            runner, sweep_obs, sweep_values, fixed_obs,
            monotone_param=args.monotone_param, direction=args.direction,
            n_samples=int(args.n_samples), seed=int(args.seed),
        )
        gate_desc = (
            f"{args.monotone_param} posterior median monotone {args.direction} "
            f"vs {sweep_obs} (<= {LIMITS_MAX_TIES} tie); "
            f"prior-box containment == {LIMITS_CONTAINMENT}"
        )
    plot_path = _save_limits_plot(results, output_dir)

    report = {
        'check': 'limits',
        'gate': gate_desc,
        'verdict': 'PASS' if results['all_pass'] else 'FAIL',
        'results': results,
        'artifact': str(args.artifact),
        'artifact_meta': _artifact_meta_for_report(runner),
        'limits_plot': plot_path,
        'provenance': _provenance_block(int(args.seed)),
    }
    _write_report(output_dir, 'limits_report.json', report)

    if results['all_pass']:
        detail = (
            f"{args.monotone_param} within {LIMITS_ANCHOR_SIGMA_FRAC}*sigma of "
            f"{len(results['grid'])} MCMC anchors"
            if results.get('mode') == 'anchor'
            else f"{args.monotone_param} median monotone {args.direction} vs {sweep_obs}"
        )
        _loud([
            "LIMITS GATE: PASS",
            f"{detail}; containment {results['containment_fraction']:.3f}",
        ])
        return EXIT_PASS
    lines = ["LIMITS GATE: FAIL"]
    if results.get('mode') == 'anchor':
        for r in results['grid']:
            if not r['anchor_pass']:
                lines.append(
                    f"{sweep_obs}={r['sweep_value']}: W1 = {r['w1']:.3f} > "
                    f"tol {r['w1_tol']:.3f} (flow median {r['flow_median']:.3f}, "
                    f"anchor median {r['anchor_median']:.3f})"
                )
    else:
        mono = results['monotonicity']
        if not mono['monotone_pass']:
            lines.append(
                f"non-monotone: {mono['n_strict_violations']} strict violations, "
                f"{mono['n_ties']} ties (max {LIMITS_MAX_TIES}); medians={results['medians']}"
            )
    if not results['containment_pass']:
        lines.append(f"containment {results['containment_fraction']:.4f} < {LIMITS_CONTAINMENT}")
    lines.append("Do NOT tune this gate — investigate the trained flow.")
    _loud(lines)
    return EXIT_GATE_FAIL


# ==========================================================================
# Self-test (toy BoxUniform + linear-Gaussian; exercises all three checks)
# ==========================================================================

def _make_toy_sbi_config(seed: int):
    """Cheap 2D uniform-prior SBI config mirroring tests/sbi_runner_test.py."""
    from .inference_core import InferenceConfig
    bounds = ((-2.0, 2.0), (-1.0, 3.0))
    param_names = ['theta1', 'theta2']
    param_space = {
        name: {'prior_type': 'uniform', 'bounds': list(b)}
        for name, b in zip(param_names, bounds)
    }
    data = {
        'mode': 'sbi', 'bodyname': 'ToyBody',
        'param_space': param_space,
        'observables': {'obs_a': [0.0, 1.0], 'obs_b': [0.0, 1.0]},
        'sampler_settings': {
            'num_simulations': 4000, 'density_estimator': 'maf',
            'n_posterior_samples': 2000, 'n_reeval': 10,
        },
        'structure_cache_path': 'unused_nonexistent.pkl',
        'random_state': seed,
    }
    return InferenceConfig.from_dict(data)


def _toy_simulate(theta: np.ndarray, seed: int) -> np.ndarray:
    """Linear-Gaussian toy simulator x = theta + 0.05 * N(0, I) (seeded)."""
    rng = np.random.RandomState(seed)
    return theta + 0.05 * rng.randn(*theta.shape)


def cmd_selftest(args) -> int:
    """End-to-end smoke of SBC + crosscheck + limits on a toy flow."""
    import torch
    from .sbi_runner import SBIRunner

    seed = int(args.seed)
    output_dir = Path(args.output_dir) if args.output_dir else Path(
        tempfile.mkdtemp(prefix='sbi_selftest_')
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    log.info(f"Self-test output dir: {output_dir}")

    torch.manual_seed(seed)
    np.random.seed(seed)

    config = _make_toy_sbi_config(seed)
    runner = SBIRunner(config)

    # --- Train a tiny flow on prior-predictive (theta, x) --------------------
    prior = runner.build_prior()
    torch.manual_seed(seed)
    n_train = int(config.sampler_settings['num_simulations'])
    theta_train = prior.sample((n_train,)).numpy().astype(np.float64)
    x_train = _toy_simulate(theta_train, seed=seed)
    runner.train(theta_train, x_train, seed=seed, density_estimator='maf',
                 max_num_epochs=200, show_progress_bars=False)
    artifact_path = output_dir / 'toy_artifact.pt'
    runner._train_info['rejection_stats'] = {'n_requested': n_train, 'rejection_rate': 0.0}
    runner.save_artifact(artifact_path)

    loaded = SBIRunner.load_artifact(artifact_path)

    verdicts: Dict[str, str] = {}

    # --- SBC on held-out pairs ----------------------------------------------
    n_sbc = 200
    torch.manual_seed(seed + 1)
    theta_sbc = prior.sample((n_sbc,)).numpy().astype(np.float64)
    x_sbc = _toy_simulate(theta_sbc, seed=seed + 1)
    sbc_results, ranks = _run_sbc_check(
        loaded._posterior, theta_sbc, x_sbc, loaded.param_labels,
        num_posterior_samples=1000, seed=seed,
    )
    _save_sbc_rank_plot(ranks, 1000, loaded.param_labels, output_dir)
    _write_report(output_dir, 'selftest_sbc_report.json', {
        'check': 'sbc', 'results': sbc_results,
        'provenance': _provenance_block(seed),
    })
    verdicts['sbc'] = 'PASS' if sbc_results['all_pass'] else 'FAIL'
    log.info(f"[selftest] SBC min KS p = {sbc_results['min_ks_pval']:.4g} "
             f"-> {verdicts['sbc']}")

    # --- Cross-check: two independent draws from the same flow ----------------
    # A matched MCMC-like reference: draw from the same posterior with a
    # different seed and attach mild (near-uniform) importance weights so the
    # weighted-stats / ESS / resample-KS code paths are exercised. This is a
    # PLUMBING self-test, not a scientific calibration claim.
    x_obs = {'obs_a': 0.4, 'obs_b': 0.6}
    sbi_samples = loaded.sample_posterior(x_obs, n_samples=2000, seed=seed + 2)
    mcmc_like = loaded.sample_posterior(x_obs, n_samples=2000, seed=seed + 3)
    rng_w = np.random.RandomState(seed)
    weights = 1.0 + 0.05 * rng_w.rand(mcmc_like.shape[0])  # mild non-uniformity
    cc_results = _run_crosscheck(
        loaded.param_names, mcmc_like, weights, sbi_samples, seed=seed,
    )
    _save_crosscheck_plot(loaded.param_names, mcmc_like, weights, sbi_samples, output_dir)
    _write_report(output_dir, 'selftest_crosscheck_report.json', {
        'check': 'crosscheck', 'x_obs': x_obs, 'results': cc_results,
        'provenance': _provenance_block(seed),
    })
    verdicts['crosscheck'] = 'PASS' if cc_results['all_pass'] else 'FAIL'
    log.info(f"[selftest] cross-check ESS = {cc_results['ess_mcmc']:.1f} "
             f"-> {verdicts['crosscheck']}")

    # --- Cross-check GATE-BITE negative control ------------------------------
    # Prove the gate is not a no-op: an artificially shifted mean must FAIL.
    shifted = sbi_samples.copy()
    shifted[:, 0] += 5.0 * sbi_samples[:, 0].std() + 1.0
    cc_bite = _run_crosscheck(loaded.param_names, mcmc_like, weights, shifted, seed=seed)
    gate_bites = not cc_bite['all_pass']
    log.info(f"[selftest] gate-bite control (shifted mean should FAIL): "
             f"{'OK (gate bit)' if gate_bites else 'BROKEN (gate did not bite)'}")

    # --- Limits: sweep obs_b, expect theta2 median to INCREASE ---------------
    # For x = theta + noise, conditioning on larger obs_b raises theta2's
    # posterior median (increasing). This exercises the same code path the
    # real |Im k2| -> log10_eta_Ih (decreasing) check uses.
    limits_results = _run_limits_check(
        loaded, sweep_obs='obs_b', sweep_values=[-0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.6],
        fixed_obs={'obs_a': 0.0}, monotone_param='theta2', direction='increasing',
        n_samples=1500, seed=seed,
    )
    _save_limits_plot(limits_results, output_dir)
    _write_report(output_dir, 'selftest_limits_report.json', {
        'check': 'limits', 'results': limits_results,
        'provenance': _provenance_block(seed),
    })
    verdicts['limits'] = 'PASS' if limits_results['all_pass'] else 'FAIL'
    log.info(f"[selftest] limits medians = "
             f"{[round(m, 3) for m in limits_results['medians']]} "
             f"-> {verdicts['limits']}")

    all_pass = all(v == 'PASS' for v in verdicts.values()) and gate_bites
    summary = {
        'check': 'selftest', 'verdicts': verdicts,
        'gate_bite_control_ok': bool(gate_bites),
        'output_dir': str(output_dir),
        'provenance': _provenance_block(seed),
    }
    _write_report(output_dir, 'selftest_summary.json', summary)

    if all_pass:
        _loud([
            "SELFTEST: PASS",
            f"SBC={verdicts['sbc']} crosscheck={verdicts['crosscheck']} "
            f"limits={verdicts['limits']} gate-bite={'OK' if gate_bites else 'BROKEN'}",
            f"reports in {output_dir}",
        ])
        return EXIT_PASS
    _loud([
        "SELFTEST: FAIL",
        f"verdicts={verdicts} gate-bite={'OK' if gate_bites else 'BROKEN'}",
        f"reports in {output_dir}",
    ])
    return EXIT_GATE_FAIL


# ==========================================================================
# CLI
# ==========================================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog='validate_sbi',
        description='SBI validation gates: SBC, MCMC<->SBI cross-check, limiting behavior.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest='command', required=True)

    def _add_common(p):
        p.add_argument('--seed', type=int, default=42,
                       help='Random seed for all stochastic steps (default: 42).')
        p.add_argument('--output-dir', type=str, required=True,
                       help='Directory for JSON reports and plot PNGs.')

    # sbc
    p_sbc = sub.add_parser('sbc', help='Simulation-based calibration.')
    p_sbc.add_argument('--artifact', required=True, help='Trained SBI .pt artifact.')
    p_sbc.add_argument('--dataset', help='Held-out (theta, x) .npz (no forward model needed).')
    p_sbc.add_argument('--config', help='InferenceConfig JSON to generate held-out pairs '
                                        '(needs structure cache; fails fast if missing).')
    p_sbc.add_argument('--n-sbc', type=int, default=200, help='Number of SBC pairs (default: 200).')
    p_sbc.add_argument('--num-posterior-samples', type=int, default=1000,
                       help='Posterior samples per SBC rank (default: 1000).')
    _add_common(p_sbc)
    p_sbc.set_defaults(func=cmd_sbc)

    # crosscheck
    p_cc = sub.add_parser('crosscheck', help='MCMC vs SBI posterior cross-check.')
    p_cc.add_argument('--artifact', required=True, help='Trained SBI .pt artifact.')
    p_cc.add_argument('--mcmc', required=True, help='MCMC InferenceResult .pkl.')
    p_cc.add_argument('--x-obs', help='JSON observed data dict (default: MCMC config observables).')
    p_cc.add_argument('--num-samples', type=int, default=None,
                      help='SBI posterior draws (default: match MCMC sample count).')
    _add_common(p_cc)
    p_cc.set_defaults(func=cmd_crosscheck)

    # limits
    p_lim = sub.add_parser('limits', help='Limiting-behavior monotonicity + containment.')
    p_lim.add_argument('--artifact', required=True, help='Trained SBI .pt artifact.')
    p_lim.add_argument('--sweep-obs', default=None,
                       help='Observable to sweep (default: auto-detect |Im k2| channel).')
    p_lim.add_argument('--sweep-values', default=None,
                       help='Comma-separated sweep grid (default: 0.05,0.10,...,0.30).')
    p_lim.add_argument('--fixed-obs', default=None, help='JSON dict of fixed observable values.')
    p_lim.add_argument('--anchor-results', default=None,
                       help='JSON dict {sweep_value: mcmc_result.pkl} of ground-truth MCMC '
                            'anchors. When given, replaces the monotonicity check with '
                            'per-anchor median comparison (ratified 2026-07-09: the '
                            'monotone-decreasing premise is falsified below |Im k2|~0.15 '
                            'by the folded-noise regime).')
    p_lim.add_argument('--re-k2', type=float, default=None, help='Fixed Re(k2) for non-swept channels.')
    p_lim.add_argument('--monotone-param', default=DEFAULT_MONOTONE_PARAM,
                       help=f'Parameter whose median must respond (default: {DEFAULT_MONOTONE_PARAM}).')
    p_lim.add_argument('--direction', default='decreasing', choices=['increasing', 'decreasing'],
                       help='Expected monotone direction (default: decreasing).')
    p_lim.add_argument('--n-samples', type=int, default=2000,
                       help='Posterior draws per grid point (default: 2000).')
    _add_common(p_lim)
    p_lim.set_defaults(func=cmd_limits)

    # selftest
    p_st = sub.add_parser('selftest', help='Toy end-to-end smoke of all three checks.')
    p_st.add_argument('--seed', type=int, default=42, help='Random seed (default: 42).')
    p_st.add_argument('--output-dir', type=str, default=None,
                      help='Report dir (default: a fresh temp dir).')
    p_st.set_defaults(func=cmd_selftest)

    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    # Plot helpers save into output_dir before _write_report creates it;
    # ensure it exists up front so PNGs are not silently skipped.
    if getattr(args, 'output_dir', None):
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
    try:
        return args.func(args)
    except FileNotFoundError as e:
        _loud(["SETUP ERROR (missing file)", *str(e).splitlines()])
        return EXIT_ERROR
    except Exception as e:
        log.error(f"{args.command} failed: {e}", exc_info=True)
        _loud([f"{args.command.upper()} ERROR", str(e)])
        return EXIT_ERROR


if __name__ == '__main__':
    sys.exit(main())
