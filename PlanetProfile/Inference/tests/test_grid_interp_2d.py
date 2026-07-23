"""Unit tests for the shared 2D (Tb × log w) bilinear interpolant (v3 salinity).

Pure-numpy; no PP/EOS import, so these run fast and standalone:

    mamba run -n PPcl pytest PlanetProfile/Inference/tests/test_grid_interp_2d.py -v
"""
import numpy as np
import pytest

from PlanetProfile.Inference.grid_interp_2d import (
    is_2d_cache,
    bilinear_weights,
    resolve_none_corners,
    blend_complex,
    blend_scalar,
)

# Small regular Tb grid × log-spaced w grid mirroring the v3 layout.
TB = np.array([260.0, 261.0, 262.0, 263.0])          # n_Tb = 4
W = np.array([0.1, 1.0, 10.0, 100.0])                # n_w = 4 (log-spaced)


def test_is_2d_cache_discriminates_format():
    assert is_2d_cache({"Tb_K_grid": TB, "wOcean_ppt_grid": W, "structures": []})
    # legacy 1D v2 cache -> False (stays servable on old path)
    assert not is_2d_cache({"Tb_K_grid": TB, "structures": []})
    assert not is_2d_cache({"grid_cache": {}, "grid_Tb_values": TB})
    assert not is_2d_cache(None)


def test_weights_sum_to_one_and_rowmajor():
    corners, weights = bilinear_weights(TB, W, 261.5, 3.0)
    assert abs(sum(weights) - 1.0) < 1e-12
    # Tb=261.5 brackets rows 1,2; w=3.0 (log10=0.477) brackets cols 1,2.
    n_w = W.size
    expected = [1 * n_w + 1, 1 * n_w + 2, 2 * n_w + 1, 2 * n_w + 2]
    assert corners == expected


def test_exact_node_recovers_that_node():
    # Query exactly on grid node (Tb=262, w=10) -> all weight on that corner.
    corners, weights = bilinear_weights(TB, W, 262.0, 10.0)
    n_w = W.size
    target = 2 * n_w + 2  # i_Tb=2, i_w=2
    total_on_target = sum(wt for c, wt in zip(corners, weights) if c == target)
    assert abs(total_on_target - 1.0) < 1e-12


def test_interp_in_log_w_not_linear_w():
    # Midpoint between w=1 and w=10 in LOG space is w=10**0.5≈3.162.
    # A query at 3.162 should split ~50/50 across the two w-corners.
    corners, weights = bilinear_weights(TB, W, 260.0, 10 ** 0.5)
    # Tb=260 is the first node -> all Tb weight on row 0; w splits 0.5/0.5.
    n_w = W.size
    wcol_lo = 0 * n_w + 1  # w=1.0
    wcol_hi = 0 * n_w + 2  # w=10.0
    w_lo = sum(wt for c, wt in zip(corners, weights) if c == wcol_lo)
    w_hi = sum(wt for c, wt in zip(corners, weights) if c == wcol_hi)
    assert abs(w_lo - 0.5) < 1e-9
    assert abs(w_hi - 0.5) < 1e-9


def test_clamp_outside_box_no_extrapolation():
    # Below the Tb box and below the w box -> clamps to node (0,0).
    corners, weights = bilinear_weights(TB, W, 100.0, 1e-6)
    n_w = W.size
    total_on_00 = sum(wt for c, wt in zip(corners, weights) if c == 0)
    assert abs(total_on_00 - 1.0) < 1e-12
    # Above the box -> clamps to the top corner.
    corners, weights = bilinear_weights(TB, W, 999.0, 1e9)
    top = (TB.size - 1) * n_w + (W.size - 1)
    total_on_top = sum(wt for c, wt in zip(corners, weights) if c == top)
    assert abs(total_on_top - 1.0) < 1e-12


def test_none_corner_low_side_renormalizes():
    # Low-Tb/low-w corner None (the tilted-band low corner). Query interior.
    corners, weights = bilinear_weights(TB, W, 261.5, 3.0)
    valid = [False, True, True, True]  # first corner (lo,lo) is None
    res = resolve_none_corners(corners, weights, valid)
    assert res is not None
    kept_corners, kept_w = res
    assert len(kept_corners) == 3
    assert abs(sum(kept_w) - 1.0) < 1e-12
    assert corners[0] not in kept_corners


def test_none_corner_high_side_renormalizes():
    # High-Tb/high-w corner None (the OTHER tilted-band corner) — reviewer
    # flagged BOTH corners must be handled, not just low-w.
    corners, weights = bilinear_weights(TB, W, 261.5, 3.0)
    valid = [True, True, True, False]  # last corner (hi,hi) is None
    res = resolve_none_corners(corners, weights, valid)
    assert res is not None
    kept_corners, kept_w = res
    assert len(kept_corners) == 3
    assert abs(sum(kept_w) - 1.0) < 1e-12
    assert corners[3] not in kept_corners


def test_all_none_rejects():
    corners, weights = bilinear_weights(TB, W, 261.5, 3.0)
    res = resolve_none_corners(corners, weights, [False, False, False, False])
    assert res is None


def test_blend_complex_matches_manual_bilinear():
    # Full 4-corner bilinear on a known complex field == analytic value.
    corners, weights = bilinear_weights(TB, W, 261.5, 10 ** 0.5)
    # define Ae = f(i_Tb, i_w) linear in indices so bilinear is exact
    n_w = W.size
    def ae(flat):
        i_T, i_w = divmod(flat, n_w)
        return complex(i_T, i_w)
    valid = [True] * 4
    res = resolve_none_corners(corners, weights, valid)
    kept_corners, kept_w = res
    got = blend_complex([ae(c) for c in kept_corners], kept_w)
    # Tb=261.5 -> i_T midpoint 1.5; w=10**0.5 -> i_w midpoint 1.5
    assert abs(got.real - 1.5) < 1e-9
    assert abs(got.imag - 1.5) < 1e-9


def test_blend_scalar_partial_corner():
    # With one corner dropped, blend is over the renormalized survivors.
    vals = [10.0, 20.0, 30.0, 40.0]
    corners = [0, 1, 2, 3]
    weights = [0.25, 0.25, 0.25, 0.25]
    res = resolve_none_corners(corners, weights, [True, True, True, False])
    kept_corners, kept_w = res
    kept_vals = [vals[c] for c in kept_corners]
    got = blend_scalar(kept_vals, kept_w)
    assert abs(got - (10 + 20 + 30) / 3) < 1e-9


# --------------------------------------------------------------------------
# Integration: 2D cache through a real MCMCRunner. Skipped unless the small
# 3×4 fixture cache exists (built by /tmp/v3_build_full_2d_cache.py's sibling
# 3×4 build). These assert the shared-interpolant wiring, not just the helper.
# --------------------------------------------------------------------------
import json
import math
import os

_FIXTURE_3x4 = "/tmp/v3_test_3x4.pkl"
_V3_CFG = "PlanetProfile/Inference/configs/europa_seawater_andrade_clipper_v3_8D.json"
_needs_fixture = pytest.mark.skipif(
    not (os.path.exists(_FIXTURE_3x4) and os.path.exists(_V3_CFG)),
    reason="requires the 3×4 v3.0 fixture cache and the v3 config",
)


def _make_runner_on_3x4():
    from PlanetProfile.Inference.inference_core import InferenceConfig
    from PlanetProfile.Inference.mcmc_runner import MCMCRunner
    with open(_V3_CFG) as f:
        cfg_dict = json.load(f)
    cfg_dict["structure_cache_path"] = _FIXTURE_3x4
    return MCMCRunner(InferenceConfig.from_dict(cfg_dict))


def _theta(runner, Tb, log10_w):
    d = {'alpha': 0.3, 'log10_zeta': 0.0, 'log10_eta_Ih': 14.0,
         'log10_eta_sil': 20.0, 'Tb_K': Tb, 'log10_wOcean_ppt': log10_w,
         'R_core_km': 400.0, 'rho_core_kgm3': 6000.0}
    return np.array([d[n] for n in runner.param_names], dtype=float)


@_needs_fixture
def test_runner_loads_2d_cache_8d():
    r = _make_runner_on_3x4()
    assert len(r.param_names) == 8
    assert 'log10_wOcean_ppt' in r.param_names
    assert is_2d_cache(r.structure_data)


@_needs_fixture
def test_support_cut_is_2d_region_in_w():
    """|Ae_synodic| rises with salinity at fixed Tb -> support cut carves a
    (Tb, w) region: low-w rejected, high-w passes."""
    r = _make_runner_on_3x4()
    # low salinity at mid-Tb: |Ae_syn| < 0.7 -> reject sentinel
    ll_low = r.log_likelihood_fn(_theta(r, 266.0, math.log10(5.0)))
    assert ll_low < -1e29
    # high salinity: |Ae_syn| > 0.7 -> finite
    ll_high = r.log_likelihood_fn(_theta(r, 269.0, 2.0))  # w=100 ppt
    assert np.isfinite(ll_high) and ll_high > -1e29


def test_unservable_sample_error_on_all_none_corner():
    """Regression (2026-07-18): a draw in an unbuilt tilted-band corner (all 4
    bilinear corners None) must raise the DEDICATED UnservableSampleError from
    _apply_bottom_temperature_2d — a ValueError subclass the likelihood/dataset
    catch to hard-reject, distinct from the rheology hooks' genuine ValueErrors.
    The pocoMC prior samples the full rectangular box including the unbuilt
    corners, so an uncaught raise here crashed the reference MCMC on sample ~0.
    """
    from PlanetProfile.Inference.forward_models import (
        _apply_bottom_temperature_2d, UnservableSampleError)
    # 2×2 cache with every structure None (fully unbuilt).
    sd = {
        'Tb_K_grid': np.array([260.0, 261.0]),
        'wOcean_ppt_grid': np.array([0.1, 100.0]),
        'structures': [None, None, None, None],
    }
    theta = {'Tb_K': 260.5, 'log10_wOcean_ppt': 0.0}  # w=1 ppt, interior
    with pytest.raises(UnservableSampleError):
        _apply_bottom_temperature_2d(theta, sd)
    # It IS a ValueError subclass (so legacy `except ValueError` still catches),
    # but the subclass lets rheology-config ValueErrors stay distinguishable.
    assert issubclass(UnservableSampleError, ValueError)


@_needs_fixture
def test_likelihood_rejects_unbuilt_corner_not_crash():
    """Regression (2026-07-18): the MCMC likelihood must return the -1e30 reject
    sentinel for an unbuilt-corner sample, NOT propagate UnservableSampleError.
    Uses the low-w corner of the fixture, which is unbuilt in the tilted band."""
    r = _make_runner_on_3x4()
    # log10_w = -1 -> w=0.1 ppt, the low-w edge that is None across the band.
    ll = r.log_likelihood_fn(_theta(r, 260.0, -1.0))
    assert ll <= -1e29  # graceful reject, no exception escaped


@_needs_fixture
def test_ae_monotone_in_salinity():
    """Shared interpolant: |Ae_synodic| increases with w at fixed Tb."""
    r = _make_runner_on_3x4()
    amps = []
    for lw in [-0.5, 0.5, 1.5]:
        ae = r._blended_ae_dict(r._expand_theta(_theta(r, 266.0, lw)))
        amps.append(abs(ae['synodic']))
    assert amps[0] < amps[1] < amps[2]


@_needs_fixture
def test_salinity_exponentiated_once():
    from PlanetProfile.Inference.grid_interp_2d import wOcean_ppt_from_theta
    r = _make_runner_on_3x4()
    td = r._expand_theta(_theta(r, 266.0, math.log10(5.0)))
    assert abs(wOcean_ppt_from_theta(td) - 5.0) < 1e-9


@_needs_fixture
def test_cmr2_sidecar_on_2d_cache_no_nameerror():
    """Regression (scientific-review condition, 2026-07-18): the CMR2-offset
    sidecar path on a 2D cache must not NameError on Tb_sample. The offset is
    interpolated in Tb only (sidecar assumed w-independent). Short-circuit
    ordering hid this because no test previously exercised sidecar + 2D.
    """
    r = _make_runner_on_3x4()
    # Inject a fake Tb-only sidecar so the offset branch executes.
    fake_sidecar = {
        'Tb_K_grid': np.array([259.5, 265.0, 271.0]),
        'offsets': np.array([0.001, 0.002, 0.003]),
    }
    r._load_cmr2_offset_sidecar = lambda: fake_sidecar
    # A valid high-w sample (passes the support cut / mass conservation).
    theta = r._expand_theta(_theta(r, 266.0, 2.0))
    cmr2 = r._derive_cmr2_via_mass_conservation(theta)
    # Must not raise; offset applied -> finite CMR2 (or a clean None reject,
    # never a NameError).
    assert cmr2 is None or np.isfinite(cmr2)


@_needs_fixture
def test_support_guard_matches_likelihood_reject_randomized():
    """Randomized equivalence sweep (scientific-review condition, 2026-07-18):
    the likelihood's reject decision must agree with an INDEPENDENT support
    check computed from the same _blended_ae_dict, sample-for-sample — proving
    the training support (SBI guard) == likelihood support. A divergence here
    is exactly the failure mode the shared interpolant exists to prevent.
    """
    r = _make_runner_on_3x4()
    bounds = r.config.induction_bounds.get('synodic', {})
    amp_min = bounds.get('amp_min')
    im_abs_max = bounds.get('im_abs_max')
    assert amp_min is not None  # fixture config carries the Galileo cut

    rng = np.random.default_rng(4848)
    n = 200
    n_checked = 0
    for _ in range(n):
        Tb = float(rng.uniform(259.5, 271.0))
        lw = float(rng.uniform(-1.0, 2.0))
        theta_arr = _theta(r, Tb, lw)
        theta_d = r._expand_theta(theta_arr)
        ae = r._blended_ae_dict(theta_d)
        # Independent support decision from the SAME blended Ae.
        if ae is None or 'synodic' not in ae:
            indep_reject = True
        else:
            A = complex(ae['synodic'])
            indep_reject = (abs(A) < float(amp_min)) or (
                im_abs_max is not None and abs(A.imag) > float(im_abs_max))
        # Likelihood reject sentinel.
        ll = r.log_likelihood_fn(theta_arr)
        ll_reject = ll < -1e29
        # The likelihood may also reject for non-induction reasons (e.g. derived
        # rho_sil out of bounds), so a likelihood-reject need not imply a
        # support-reject; but a support-reject MUST force a likelihood-reject.
        if indep_reject:
            assert ll_reject, (
                f"support says reject but likelihood accepted at "
                f"Tb={Tb:.3f} log10w={lw:.3f} (|Ae_syn| check)")
            n_checked += 1
    assert n_checked > 0  # the sweep actually exercised the reject branch


# --------------------------------------------------------------------------
# v5 ice-thickness reparameterization: D_iceIh(Tb, w) -> Tb inversion.
# Pure-numpy on a synthetic monotone-decreasing D field; no cache/EOS import.
# --------------------------------------------------------------------------
from PlanetProfile.Inference.grid_interp_2d import (
    forward_d_iceIh,
    invert_d_iceIh_to_Tb,
    d_iceIh_column_at_w,
)

# 5-node Tb grid × 4-node log-w grid. D decreases with Tb (warmer=thinner) and
# decreases with w at fixed Tb (saltier=thinner) — the real physics.
_TB5 = np.array([260.0, 262.0, 264.0, 266.0, 268.0])
_W4 = np.array([0.1, 1.0, 10.0, 100.0])


def _synthetic_D_flat(none_corners=()):
    """Row-major D_iceIh, D = 40 - 3*(Tb-260) - 2*log10(w), some None corners."""
    n_w = _W4.size
    flat = []
    for i_T, Tb in enumerate(_TB5):
        for i_w, w in enumerate(_W4):
            if (i_T, i_w) in none_corners:
                flat.append(None)
            else:
                flat.append(40.0 - 3.0 * (Tb - 260.0) - 2.0 * np.log10(w))
    return flat


def test_v5_forward_matches_bilinear_on_D():
    """forward_d_iceIh is the same bilinear operator used for other scalars."""
    D_flat = _synthetic_D_flat()
    # interior query; D linear in Tb and log10 w so bilinear is exact
    val = forward_d_iceIh(_TB5, _W4, D_flat, 263.0, 10 ** 0.5)
    # analytic: 40 - 3*3 - 2*0.5 = 40 - 9 - 1 = 30
    assert abs(val - 30.0) < 1e-9


def test_v5_inversion_round_trip_exact():
    """forward(invert(D, w), w) == D on the valid band, to bisection tol."""
    D_flat = _synthetic_D_flat()
    rng = np.random.default_rng(0)
    for _ in range(200):
        Tb = rng.uniform(_TB5[0], _TB5[-1])
        w = 10 ** rng.uniform(-1, 2)
        D = forward_d_iceIh(_TB5, _W4, D_flat, Tb, w)
        Tb_inv = invert_d_iceIh_to_Tb(_TB5, _W4, D_flat, D, w)
        assert Tb_inv is not None
        D_re = forward_d_iceIh(_TB5, _W4, D_flat, Tb_inv, w)
        assert abs(D_re - D) < 1e-5


def test_v5_inversion_edge_rejects_not_clips():
    """Targets outside [D_min(w), D_max(w)] return None (reject), not the edge."""
    D_flat = _synthetic_D_flat()
    for w in (0.1, 1.0, 10.0, 100.0):
        Tbs, Ds = d_iceIh_column_at_w(_TB5, _W4, D_flat, w)
        dmin, dmax = Ds.min(), Ds.max()
        assert invert_d_iceIh_to_Tb(_TB5, _W4, D_flat, dmin - 5.0, w) is None
        assert invert_d_iceIh_to_Tb(_TB5, _W4, D_flat, dmax + 5.0, w) is None
        mid = invert_d_iceIh_to_Tb(_TB5, _W4, D_flat, 0.5 * (dmin + dmax), w)
        assert mid is not None and _TB5[0] <= mid <= _TB5[-1]


def test_v5_inversion_monotone_thicker_is_colder():
    """At fixed w, a thicker ice target inverts to a colder Tb."""
    D_flat = _synthetic_D_flat()
    Tb_thin = invert_d_iceIh_to_Tb(_TB5, _W4, D_flat, 25.0, 10.0)
    Tb_thick = invert_d_iceIh_to_Tb(_TB5, _W4, D_flat, 34.0, 10.0)
    assert Tb_thick < Tb_thin


def test_v5_inversion_survives_none_corner():
    """A None at one tilted-band corner still inverts (nearest-valid in w),
    and the round-trip stays consistent with the same forward operator."""
    D_flat = _synthetic_D_flat(none_corners=((0, 0),))  # low-Tb/low-w None
    w = 0.1
    Tbs, Ds = d_iceIh_column_at_w(_TB5, _W4, D_flat, w)
    assert Tbs.size >= 2  # column still usable via nearest-valid in w
    Dt = 0.5 * (Ds.min() + Ds.max())
    Tb_inv = invert_d_iceIh_to_Tb(_TB5, _W4, D_flat, Dt, w)
    assert Tb_inv is not None
    D_re = forward_d_iceIh(_TB5, _W4, D_flat, Tb_inv, w)
    assert abs(D_re - Dt) < 1e-4


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
