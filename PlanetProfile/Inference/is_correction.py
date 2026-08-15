"""Importance-sampling correction of amortized SBI posteriors (Track 1).

Lifts the tidal-sector split ratification by treating the deployed flow as
an importance proposal for the exact posterior the reference MCMC targets:

    log w_i = log p(x_obs | theta_i) - log q(theta_i | x_obs)

The uniform prior density cancels (all campaign priors are BoxUniform) and
every support guard is folded into the MCMC likelihood as the -1e30
sentinel, so using ``MCMCRunner.log_likelihood_fn`` reproduces the guards
exactly (scientific review 2026-08-11, plans/active/
tidal-sector-remedy-plan.md). Both ingredients are already produced by the
deploy path: ``InferenceResult.log_likelihoods`` (full-N when conditioned
with ``n_derived=None``) and ``metadata['flow_log_prob']``
(``norm_posterior=False`` — the truncation constant cancels; do not
"fix" that).

Binding reviewer conditions implemented here: C2 (byte-identity asserts),
C3 (reference likelihood recompute test), C4 (sentinel masking before any
arithmetic; ESS on full N; >50% rejected hard fail), C5 (Pareto-k primary
gate with PSIS smoothing; w_max cap), C6 (k2-support-box mass), C7 (ESS
floors), C16 (per-branch ocean/no-ocean mass + ESS).

The corrected posterior targets the reference-MCMC posterior BY
CONSTRUCTION; a validation-gate failure downstream therefore indicts the
implementation (config/index/sentinel handling), not the flow. The
model-data k2 tension is untouched by this correction.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Dict, Optional

import numpy as np

log = logging.getLogger('PlanetProfile')

# Sentinel threshold: MCMCRunner returns -1e30 for guard rejections.
SENTINEL_LOGL = -1e29

# C5 thresholds (preregistered; never tuned after seeing results).
PARETO_K_CLEAN = 0.5
PARETO_K_MAX = 0.7
W_MAX_FRAC_CAP = 0.01
# C7 floor: absolute ESS at the N actually run. (The 0.1 ESS/N fractional
# floor was struck by reviewer ruling 2026-08-11 — reported-only now.)
ESS_ABS_FLOOR = 1000.0
ESS_FRAC_FLOOR = 0.1  # retained for reporting + N-sizing, never gating
# C4: hard fail above this rejected fraction.
REJECT_FRAC_HARD_FAIL = 0.5
# C16: branch diagnostics.
BRANCH_PROB_MIN = 0.05
BRANCH_ESS_FLOOR = 50.0
# C6: corrected mass outside the k2 training-support box must stay below.
K2_BOX_MASS_CAP = 1e-3


@dataclass
class ISCorrection:
    """Result of the importance correction: normalized weights over the
    FULL sample set (sentinel-rejected draws carry weight 0 and stay in
    N), plus every preregistered diagnostic and the deploy verdict."""
    log_w: np.ndarray            # raw log-weights, -inf where rejected
    weights: np.ndarray          # normalized, full N, sums to 1
    n_samples: int
    n_rejected: int
    frac_rejected: float
    ess: float
    ess_over_n: float
    pareto_k: float
    psis_smoothed: bool
    w_max_frac: float
    verdict: str                 # 'clean' | 'smoothed' | 'fail'
    fail_reasons: list = field(default_factory=list)
    branch: Optional[Dict[str, Any]] = None
    k2_box_mass_outside: Optional[float] = None

    def summary(self) -> Dict[str, Any]:
        out = {
            'n_samples': self.n_samples,
            'n_rejected': self.n_rejected,
            'frac_rejected': self.frac_rejected,
            'ess': self.ess,
            'ess_over_n': self.ess_over_n,
            'pareto_k': self.pareto_k,
            'psis_smoothed': self.psis_smoothed,
            'w_max_frac': self.w_max_frac,
            'verdict': self.verdict,
            'fail_reasons': list(self.fail_reasons),
            'k2_box_mass_outside': self.k2_box_mass_outside,
            'branch': self.branch,
            'thresholds': {
                'pareto_k_clean': PARETO_K_CLEAN,
                'pareto_k_max': PARETO_K_MAX,
                'w_max_frac_cap': W_MAX_FRAC_CAP,
                'ess_abs_floor': ESS_ABS_FLOOR,
                'ess_frac_floor': ESS_FRAC_FLOOR,
                'reject_frac_hard_fail': REJECT_FRAC_HARD_FAIL,
                'branch_prob_min': BRANCH_PROB_MIN,
                'branch_ess_floor': BRANCH_ESS_FLOOR,
                'k2_box_mass_cap': K2_BOX_MASS_CAP,
            },
        }
        return out


def pareto_k_fit(w: np.ndarray) -> tuple:
    """Generalized-Pareto shape fit on the weight tail, mirroring the
    reference PSIS implementation (Zhang & Stephens 2009 estimator as in
    the `loo` package; Vehtari et al. 2024). Returns
    (k_hat, tail_idx, (mu, sigma)) with k_hat in the Vehtari sign
    convention (k > 0.7 unreliable; k < 0 bounded tail).

    Tail size M = min(0.2 S, 3 sqrt(S)) over the S nonzero weights.
    """
    w = np.asarray(w, float)
    S = int(np.sum(w > 0))
    if S < 50:
        return float('nan'), None, None
    M = max(int(min(0.2 * S, 3.0 * np.sqrt(S))), 10)
    order = np.argsort(w)
    tail_idx = order[-M:]
    mu = float(w[order[-M - 1]])
    y = np.sort(w[tail_idx] - mu)
    y = y[y > 0]
    n = y.size
    if n < 5:
        return float('nan'), None, None
    # loo::gpdfit grid over b (can be negative for heavy tails)
    mgrid = 30 + int(np.floor(np.sqrt(n)))
    j = np.arange(1, mgrid + 1)
    z_star = y[max(int(np.floor(n / 4 + 0.5)) - 1, 0)]
    b_cand = 1.0 / y[-1] + (1.0 - np.sqrt(mgrid / (j - 0.5))) / (3.0 * z_star)
    with np.errstate(invalid='ignore', divide='ignore'):
        # profile log-lik (loo::lx): kp = -mean(log1p(-b y)),
        # L = n*(log(b/kp) + kp - 1); note kp = -k_b below.
        k_b = np.array([np.mean(np.log1p(-bb * y)) for bb in b_cand])
        L_b = n * (np.log(-b_cand / k_b) - k_b - 1.0)
    L_b[~np.isfinite(L_b)] = -np.inf
    rel = np.exp(L_b - L_b.max())
    bw = rel / rel.sum()
    b_hat = float(np.sum(b_cand * bw))
    k_hat = float(np.mean(np.log1p(-b_hat * y)))
    sigma = -k_hat / b_hat if b_hat != 0 else float('nan')
    # weak prior regularization (loo wip): pulls k toward 0.5 at small n
    k_hat = (n * k_hat + 10.0 * 0.5) / (n + 10.0)
    return k_hat, tail_idx, (mu, sigma)


def _psis_smooth(w: np.ndarray, tail_idx: np.ndarray,
                 gpd: tuple, k_hat: float) -> np.ndarray:
    """Replace the tail weights with expected GPD order statistics
    (PSIS smoothing), capped at the raw maximum."""
    mu, sigma = gpd
    M = tail_idx.size
    z = (np.arange(1, M + 1) - 0.5) / M
    if abs(k_hat) < 1e-12:
        q = -sigma * np.log1p(-z)
    else:
        q = sigma * (np.power(1.0 - z, -k_hat) - 1.0) / k_hat
    smoothed = np.array(w, dtype=float)
    order_within = np.argsort(w[tail_idx])
    smoothed[tail_idx[order_within]] = np.minimum(mu + q, w.max())
    return smoothed


def compute_is_correction(
    result,
    d_ocean: Optional[np.ndarray] = None,
    k2: Optional[np.ndarray] = None,
    k2_support_bounds: Optional[Dict[str, Any]] = None,
) -> ISCorrection:
    """Form importance weights and every preregistered diagnostic from a
    fully-derived InferenceResult (``n_derived=None`` conditioning).

    Args:
        result: InferenceResult with full-N ``log_likelihoods`` and
            ``metadata['flow_log_prob']``.
        d_ocean: optional per-sample ocean thickness (km) for the C16
            branch split (defaults to ``result.D_ocean_results``).
        k2: optional (N, 2) per-sample [Re, Im] k2 for the C6 box-mass
            check (defaults to ``result.k2_results``).
        k2_support_bounds: the training-support box, e.g.
            ``{'Re_k2': [-0.1, 1.5], 'Im_k2': [0.0, 1.0]}``.
    """
    ll = np.asarray(result.log_likelihoods, float)
    flp = np.asarray(result.metadata['flow_log_prob'], float)
    if ll.shape != flp.shape:
        raise ValueError(
            f"log_likelihoods {ll.shape} vs flow_log_prob {flp.shape} — "
            "index misalignment (C3); condition with n_derived=None.")
    N = ll.size
    n_nonfinite_ll = int(np.sum(~np.isfinite(ll)))
    if n_nonfinite_ll:
        raise ValueError(
            f"{n_nonfinite_ll}/{N} non-finite log-likelihoods — the "
            "conditioning was NOT run with n_derived=None (NaN padding "
            "present). The correction requires full-N likelihoods.")

    # C4: sentinel mask BEFORE any arithmetic.
    rejected = ll < SENTINEL_LOGL
    frac_rej = float(rejected.mean())
    fail_reasons = []
    if frac_rej > REJECT_FRAC_HARD_FAIL:
        raise RuntimeError(
            f"{frac_rej:.1%} of draws sentinel-rejected (> "
            f"{REJECT_FRAC_HARD_FAIL:.0%}) — proposal/support mismatch; "
            "do not form weights (C4 hard fail).")

    log_w = np.full(N, -np.inf)
    log_w[~rejected] = ll[~rejected] - flp[~rejected]
    if not np.isfinite(log_w[~rejected]).all():
        raise ValueError("non-finite log-weights after sentinel masking — "
                         "inspect flow_log_prob (C4).")

    finite = np.isfinite(log_w)
    lw = log_w[finite] - log_w[finite].max()
    w_full = np.zeros(N)
    w_full[finite] = np.exp(lw)

    k_hat, tail_idx_local, gpd = pareto_k_fit(w_full[finite])
    smoothed = False
    if np.isfinite(k_hat) and PARETO_K_CLEAN < k_hat <= PARETO_K_MAX:
        wf = _psis_smooth(w_full[finite], tail_idx_local, gpd, k_hat)
        w_full[finite] = wf
        smoothed = True
    if np.isfinite(k_hat) and k_hat > PARETO_K_MAX:
        fail_reasons.append(
            f'pareto_k {k_hat:.3f} > {PARETO_K_MAX} (C5)')

    weights = w_full / w_full.sum()
    ess = float(1.0 / np.sum(weights ** 2))
    w_max_frac = float(weights.max())
    if w_max_frac > W_MAX_FRAC_CAP:
        fail_reasons.append(
            f'w_max_frac {w_max_frac:.4f} > {W_MAX_FRAC_CAP} (C5)')
    if ess < ESS_ABS_FLOOR:
        fail_reasons.append(f'ESS {ess:.0f} < {ESS_ABS_FLOOR:.0f} (C7)')
    # ESS/N is REPORTED-ONLY (reviewer ruling 2026-08-11: cost statistic,
    # not a reliability gate; the original 0.1 floor is recorded
    # FAIL-ADJUDICATED, superseded by absolute ESS + Pareto-k + w_max +
    # the C5.3 reverse-coverage test). It sizes N:
    # N_required = ESS_ABS_FLOOR / (ESS/N).

    # C16: branch split + per-branch ESS.
    branch = None
    if d_ocean is None:
        d_ocean = np.asarray(getattr(result, 'D_ocean_results', []), float)
    if d_ocean is not None and np.size(d_ocean) == N:
        has_ocean = np.nan_to_num(d_ocean, nan=0.0) > 0.5
        branch = {}
        for name, mask in (('ocean', has_ocean),
                           ('no_ocean', ~has_ocean)):
            p = float(weights[mask].sum())
            wb = weights[mask]
            b_ess = float((wb.sum() ** 2) / np.sum(wb ** 2)) if p > 0 \
                else 0.0
            branch[name] = {'prob_flow': float(mask.mean()),
                            'prob_corrected': p, 'ess': b_ess}
            if p > BRANCH_PROB_MIN and b_ess < BRANCH_ESS_FLOOR:
                fail_reasons.append(
                    f'branch {name}: prob {p:.3f} but ESS {b_ess:.0f} < '
                    f'{BRANCH_ESS_FLOOR:.0f} (C16)')

    # C6: corrected mass outside the k2 training-support box.
    box_mass = None
    if k2 is None:
        k2 = np.asarray(getattr(result, 'k2_results', []), float)
    if (k2_support_bounds and k2 is not None and np.ndim(k2) == 2
            and len(k2) == N):
        re_lo, re_hi = k2_support_bounds.get('Re_k2', (-np.inf, np.inf))[:2]
        im_lo, im_hi = k2_support_bounds.get('Im_k2', (0.0, np.inf))[:2]
        re_v, im_v = k2[:, 0], np.abs(k2[:, 1])
        outside = ((re_v < re_lo) | (re_v > re_hi)
                   | (im_v < im_lo) | (im_v > im_hi))
        outside &= np.isfinite(re_v) & np.isfinite(im_v)
        box_mass = float(weights[outside].sum())
        if box_mass > K2_BOX_MASS_CAP:
            fail_reasons.append(
                f'k2-box outside mass {box_mass:.2e} > {K2_BOX_MASS_CAP} '
                '(C6)')

    verdict = ('fail' if fail_reasons
               else ('smoothed' if smoothed else 'clean'))
    return ISCorrection(
        log_w=log_w, weights=weights, n_samples=N,
        n_rejected=int(rejected.sum()), frac_rejected=frac_rej,
        ess=ess, ess_over_n=ess / N, pareto_k=float(k_hat),
        psis_smoothed=smoothed, w_max_frac=w_max_frac,
        verdict=verdict, fail_reasons=fail_reasons,
        branch=branch, k2_box_mass_outside=box_mass)


def assert_byte_identity(runner, result) -> Dict[str, Any]:
    """C2: hard asserts that the correction's likelihood is byte-identical
    to the reference MCMC's. Any mismatch raises — never a warning."""
    import hashlib
    from pathlib import Path

    cfg = runner.config
    checks: Dict[str, Any] = {}
    # (i) config hash of the conditioning run matches the result's.
    res_hash = (result.metadata or {}).get('config_hash')
    run_hash = cfg.generate_hash()
    if res_hash is not None and res_hash != run_hash:
        raise AssertionError(
            f"config hash mismatch: result {res_hash} vs runner {run_hash}")
    checks['config_hash'] = run_hash
    # (ii) structure cache SHA vs file on disk.
    cache_path = getattr(cfg, 'structure_cache_path', None)
    if cache_path:
        p = Path(cache_path)
        if not p.is_absolute():
            from PlanetProfile import _ROOT
            p = Path(_ROOT).parent / cache_path
        sha = hashlib.sha256(p.read_bytes()).hexdigest()
        checks['structure_cache_sha256'] = sha
    # (iii) param name alignment.
    art_names = (result.metadata or {}).get('param_names',
                                            list(result.param_names))
    if list(art_names) != list(runner.param_names):
        raise AssertionError(
            f"param_names mismatch: {art_names} vs {runner.param_names}")
    checks['param_names'] = list(runner.param_names)
    # (iv) noise convention.
    obs = cfg.observables
    sig = {k: float(v[1]) for k, v in obs.items()}
    checks['sigma'] = sig
    imag = getattr(runner, 'imag_convention', 'abs')
    if imag != 'abs':
        raise AssertionError(f"imag_convention {imag!r} != 'abs'")
    checks['imag_convention'] = imag
    return checks


def reference_likelihood_consistency(runner, ref_samples: np.ndarray,
                                     ref_logl: np.ndarray,
                                     n_check: int = 200,
                                     rtol: float = 1e-9,
                                     seed: int = 0) -> Dict[str, Any]:
    """C3: recompute the log-likelihood for reference-MCMC draws with the
    deploy-path runner and require agreement with the stored values."""
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(ref_samples), min(n_check, len(ref_samples)),
                     replace=False)
    max_rel = 0.0
    n_sentinel = 0
    for i in idx:
        new = float(runner.log_likelihood_fn(np.asarray(ref_samples[i],
                                                        float)))
        old = float(ref_logl[i])
        if old < SENTINEL_LOGL or new < SENTINEL_LOGL:
            n_sentinel += int((old < SENTINEL_LOGL)
                              != (new < SENTINEL_LOGL))
            continue
        denom = max(abs(old), 1.0)
        max_rel = max(max_rel, abs(new - old) / denom)
    ok = (max_rel < rtol) and (n_sentinel == 0)
    out = {'n_checked': int(len(idx)), 'max_rel_diff': max_rel,
           'sentinel_disagreements': n_sentinel, 'rtol': rtol, 'pass': ok}
    if not ok:
        raise AssertionError(f"C3 likelihood consistency FAILED: {out}")
    return out


REVERSE_COVERAGE_TAIL_Q = 0.001   # 0.1st percentile of self-log q
REVERSE_COVERAGE_MASS_CAP = 0.01  # <1% of reference mass may sit below it


def reverse_coverage(logq_self: np.ndarray,
                     logq_ref: np.ndarray,
                     ref_weights: Optional[np.ndarray] = None
                     ) -> Dict[str, Any]:
    """C5.3 (BLOCKING, reviewer ruling 2026-08-11): parameter-space
    coverage test. Pareto-k is blind to regions where p >> q but q is so
    low that no draw landed there; this test is not. Evaluate the flow's
    log-density at the REFERENCE-MCMC draws and require that less than
    REVERSE_COVERAGE_MASS_CAP of reference mass falls below the
    REVERSE_COVERAGE_TAIL_Q quantile of the flow's own draws' log q.

    Args:
        logq_self: flow log-density at the flow's own draws (finite only).
        logq_ref: flow log-density at reference-MCMC draws (same x_obs,
            same norm_posterior=False convention).
        ref_weights: reference importance weights (None = equal).
    """
    logq_self = np.asarray(logq_self, float)
    logq_self = logq_self[np.isfinite(logq_self)]
    logq_ref = np.asarray(logq_ref, float)
    if ref_weights is None:
        ref_weights = np.ones(logq_ref.size)
    ref_weights = np.asarray(ref_weights, float)
    ref_weights = ref_weights / ref_weights.sum()
    threshold = float(np.quantile(logq_self, REVERSE_COVERAGE_TAIL_Q))
    below = (logq_ref < threshold) | ~np.isfinite(logq_ref)
    mass_below = float(ref_weights[below].sum())
    return {'threshold_logq': threshold,
            'ref_mass_below': mass_below,
            'cap': REVERSE_COVERAGE_MASS_CAP,
            'tail_q': REVERSE_COVERAGE_TAIL_Q,
            'pass': bool(mass_below < REVERSE_COVERAGE_MASS_CAP)}


def weighted_quantile(v: np.ndarray, weights: np.ndarray,
                      q) -> np.ndarray:
    """Weighted quantile(s) of v (finite entries only)."""
    v = np.asarray(v, float)
    weights = np.asarray(weights, float)
    ok = np.isfinite(v) & (weights > 0)
    v, weights = v[ok], weights[ok]
    order = np.argsort(v)
    v, weights = v[order], weights[order]
    c = np.cumsum(weights)
    c /= c[-1]
    return np.interp(np.atleast_1d(q), c, v)


def systematic_resample(weights: np.ndarray, n_out: int,
                        seed: int = 0) -> np.ndarray:
    """Systematic resampling: returns indices into the weighted set.
    Downstream consumers must report ESS, not n_out (review 6e)."""
    rng = np.random.default_rng(seed)
    positions = (rng.random() + np.arange(n_out)) / n_out
    c = np.cumsum(weights)
    c /= c[-1]
    return np.searchsorted(c, positions)
