"""Unit tests for the Track 1 importance-sampling correction
(plans/active/tidal-sector-remedy-plan.md, reviewer conditions C2-C16)."""
import sys
import types
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from PlanetProfile.Inference.is_correction import (  # noqa: E402
    ISCorrection, compute_is_correction, pareto_k_fit,
    systematic_resample, weighted_quantile, SENTINEL_LOGL)


def _mk_result(ll, flp, d_ocean=None, k2=None):
    return types.SimpleNamespace(
        log_likelihoods=np.asarray(ll, float),
        metadata={'flow_log_prob': np.asarray(flp, float)},
        D_ocean_results=(np.asarray(d_ocean, float)
                         if d_ocean is not None else np.array([])),
        k2_results=(np.asarray(k2, float) if k2 is not None
                    else np.array([])),
    )


def test_pareto_k_sign_convention():
    rng = np.random.default_rng(1)
    # Bounded/light tail (exponential weights): k well below 0.7
    w_light = rng.exponential(size=5000)
    k_light, _, _ = pareto_k_fit(w_light)
    assert k_light < 0.5, k_light
    # Heavy tail: weights ~ Pareto(alpha=1) has k ~ 1 -> must exceed 0.7
    w_heavy = (1.0 / rng.random(5000))
    k_heavy, _, _ = pareto_k_fit(w_heavy)
    assert k_heavy > 0.7, k_heavy


def test_ideal_proposal_passes_clean():
    # q == p: all weights equal -> ESS == N, k tiny, verdict clean
    rng = np.random.default_rng(2)
    N = 5000
    ll = rng.normal(-10, 1e-9, N)
    flp = ll - 3.0  # constant offset cancels in normalization
    res = _mk_result(ll, flp)
    c = compute_is_correction(res)
    assert c.verdict == 'clean', (c.verdict, c.fail_reasons)
    assert c.ess == pytest.approx(N, rel=1e-6)
    assert c.frac_rejected == 0.0


def test_sentinel_masked_before_arithmetic():
    rng = np.random.default_rng(3)
    N = 4000
    ll = rng.normal(-5, 0.3, N)
    ll[:100] = -1e30  # guard rejections
    flp = rng.normal(-4, 0.3, N)
    res = _mk_result(ll, flp)
    c = compute_is_correction(res)
    assert c.n_rejected == 100
    assert np.all(c.weights[:100] == 0.0)
    # ESS defined on full N (rejects included in the count)
    assert c.n_samples == N


def test_hard_fail_on_majority_rejection():
    N = 1000
    ll = np.full(N, -1e30)
    ll[:400] = -5.0
    res = _mk_result(ll, np.zeros(N))
    with pytest.raises(RuntimeError, match='sentinel-rejected'):
        compute_is_correction(res)


def test_nan_padding_rejected():
    # n_derived subset conditioning (NaN padding) must be refused
    ll = np.array([-5.0, np.nan, -6.0, np.nan])
    res = _mk_result(ll, np.zeros(4))
    with pytest.raises(ValueError, match='n_derived=None'):
        compute_is_correction(res)


def test_branch_split_and_gate():
    rng = np.random.default_rng(4)
    N = 6000
    ll = rng.normal(-5, 0.2, N)
    flp = rng.normal(-5, 0.2, N)
    d_oc = np.where(rng.random(N) < 0.9, 100.0, 0.0)
    res = _mk_result(ll, flp, d_ocean=d_oc)
    c = compute_is_correction(res)
    assert c.branch is not None
    assert c.branch['ocean']['prob_corrected'] == pytest.approx(0.9,
                                                               abs=0.03)
    assert (c.branch['ocean']['prob_corrected']
            + c.branch['no_ocean']['prob_corrected']) == pytest.approx(1.0)
    # concentrate ALL corrected mass on 20 no-ocean draws -> per-branch
    # ESS floor must trip (C16)
    ll2 = np.full(N, -80.0)
    ll2[:20] = 0.0
    d2 = np.zeros(N); d2[:20] = 0.0  # the favored draws are no-ocean
    d2[20:] = 100.0
    res2 = _mk_result(ll2, np.zeros(N), d_ocean=d2)
    c2 = compute_is_correction(res2)
    assert any('branch no_ocean' in r for r in c2.fail_reasons), \
        c2.fail_reasons


def test_k2_box_mass_gate():
    rng = np.random.default_rng(5)
    N = 3000
    ll = rng.normal(-5, 0.2, N)
    flp = rng.normal(-5, 0.2, N)
    k2 = np.column_stack([np.full(N, 0.5), np.full(N, 0.1)])
    k2[:30, 0] = 2.0  # outside Re box [-0.1, 1.5] with ~1% mass
    res = _mk_result(ll, flp, k2=k2)
    c = compute_is_correction(
        res, k2_support_bounds={'Re_k2': [-0.1, 1.5],
                                'Im_k2': [0.0, 1.0]})
    assert c.k2_box_mass_outside == pytest.approx(30 / N, rel=0.5)
    assert any('k2-box' in r for r in c.fail_reasons)


def test_weighted_quantile_and_resample():
    v = np.array([1.0, 2.0, 3.0, 4.0])
    w = np.array([0.0, 0.0, 1.0, 1.0])
    med = weighted_quantile(v, w / w.sum(), 0.5)[0]
    assert med in (3.0, 4.0) or 3.0 <= med <= 4.0
    idx = systematic_resample(np.array([0.5, 0.5, 0.0]), 1000, seed=0)
    assert set(np.unique(idx)) <= {0, 1}
    counts = np.bincount(idx, minlength=3)
    assert abs(counts[0] - 500) < 60


def test_shape_mismatch_raises():
    res = _mk_result(np.zeros(10), np.zeros(9))
    with pytest.raises(ValueError, match='misalignment'):
        compute_is_correction(res)


def test_reverse_coverage_gate():
    from PlanetProfile.Inference.is_correction import reverse_coverage
    rng = np.random.default_rng(6)
    logq_self = rng.normal(-10, 1, 20000)
    # covered reference: same distribution -> ~0.1% below threshold
    ref_ok = rng.normal(-10, 1, 5000)
    r1 = reverse_coverage(logq_self, ref_ok)
    assert r1['pass'], r1
    # uncovered reference: 10% of mass sits far below the flow's support
    ref_bad = np.concatenate([rng.normal(-10, 1, 4500),
                              np.full(500, -50.0)])
    r2 = reverse_coverage(logq_self, ref_bad)
    assert not r2['pass'] and r2['ref_mass_below'] >= 0.09, r2


def test_ess_over_n_reported_not_gating():
    # low ESS/N with healthy absolute ESS must NOT fail (reviewer ruling
    # 2026-08-11: fractional floor struck)
    rng = np.random.default_rng(7)
    N = 40000
    # displaced proposal: ESS/N ~ exp(-2.7) ~ 0.067 -> ESS ~ 2680
    ll = rng.normal(0, 1, N) * 1.6
    flp = np.zeros(N)
    res = types.SimpleNamespace(
        log_likelihoods=ll, metadata={'flow_log_prob': flp},
        D_ocean_results=np.array([]), k2_results=np.array([]))
    c = compute_is_correction(res)
    assert c.ess >= 1000, c.ess
    assert c.ess_over_n < 0.1
    assert not any('ESS/N' in r for r in c.fail_reasons), c.fail_reasons
