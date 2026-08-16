"""E1 (r6 build blocker): the zb-placement invariant must be ACHIEVABLE.

r6 measured 8/8 real ocean builds rejected by the zb-placement invariant, by
5-27x, with residuals to -3.4 km at the posterior heart. Three mechanisms, all
pinned here:

1. ``GetIceShellTFreeze`` returned the LAST root-finder iterate, not the
   minimum-|residual| one. The terminal iterate is whatever the root-finder
   probed last; at Enceladus |d(zb)/d(Tb)| ~ 120 km/K it was routinely several
   km worse than an iterate the solver had already computed and thrown away.
2. That root-find's tolerance was ``Planet.TfreezeRes_K`` = 0.05 K -- a
   phase-diagram STEP SIZE shared with GetTfreeze/SetupInit, worth ~6 km of zb
   at Enceladus. It is now DERIVED per node from a km tolerance the caller
   declares (``Bulk.zbTol_km``), leaving TfreezeRes_K -- and therefore every
   Titan/Europa path -- untouched.
3. The builder scored the achieved zb with ``D_iceIh_km``, which undercounts
   the shell by exactly one radial cell (measured 0.438-0.495 km at Enceladus,
   3.5-4x the finest segment's whole 0.125 km budget). The invariant now reads
   the ice/ocean interface depth from the node's own arrays.

Every test here is fast and structural. The production-tolerance evidence
(real builds through ``build_zbw_grid_cache``, one node per config grid
segment) is recorded in the r7 report.
"""
import copy
import types

import numpy as np
import pytest

from PlanetProfile.Inference import cache_builder
from PlanetProfile.Inference.cache_builder import (
    build_zbw_grid_cache, _ocean_node_zb_km)
from PlanetProfile.Thermodynamics import LayerPropagators
from PlanetProfile.Utilities.defineStructs import BulkSubstruct


R_BODY_M = 2.521e5

# The iterate sequence measured on the real (zb = 25 km, w = 20 ppt) node --
# the "posterior heart" node r6 recorded failing by -3.4 km.
R6_T_SEQ = [270.6359375, 272.0203125, 271.84010128680694,
            271.88685729407996, 271.9118572940801]
R6_ZB_SEQ = [135.81958461365758, 8.415103487097731, 30.810515495697107,
             25.115671410365508, 22.035091610841235]


def _struct(zb_km, d_iceIh_km=None):
    """Ocean-branch node whose ARRAYS say zb_km and whose D_iceIh_km may not."""
    return {
        'r_m': np.array([1.0e5, R_BODY_M - float(zb_km) * 1e3, R_BODY_M]),
        'rho': np.array([2350.0, 1005.0, 925.0]),
        'phases': np.array([50, 0, 1]),
        'D_iceIh_km': float(zb_km if d_iceIh_km is None else d_iceIh_km),
        'D_ocean_km': 30.0,
        'CMR2': 0.3345,
        'Tb_K': 272.3,
        'Mtot_kg': 1.08022e20,
        'R_body_m': R_BODY_M,
    }


# ---------------------------------------------------------------------------
# Mechanism 3 -- what the invariant reads
# ---------------------------------------------------------------------------

def test_achieved_zb_is_the_ice_ocean_interface_not_D_iceIh():
    s = _struct(25.0, d_iceIh_km=24.505)  # the measured Enceladus undercount
    assert _ocean_node_zb_km(s) == pytest.approx(25.0, abs=1e-12)
    assert abs(s['D_iceIh_km'] - 25.0) > 0.125  # and it would have failed


def test_achieved_zb_is_nan_without_a_liquid_cell():
    s = _struct(25.0)
    s['phases'] = np.array([50, 1, 1])
    assert np.isnan(_ocean_node_zb_km(s))


def test_a_one_cell_undercount_no_longer_rejects_a_good_node(monkeypatch,
                                                             tmp_path):
    """THE E1 PIN. A node placed EXACTLY on target, whose D_iceIh_km is short
    by the measured one-cell offset, survives at PRODUCTION tolerance. Under
    the pre-E1 reading this same node was rejected by 3.5x."""
    def _fake(module, Tb_K, **kwargs):
        zb = kwargs['bulk_overrides']['zb_approximate_km']
        return _struct(zb, d_iceIh_km=zb - 0.44)

    monkeypatch.setattr(cache_builder, 'build_single_structure', _fake)
    monkeypatch.setattr(cache_builder, '_save_cache_atomic',
                        lambda cache, path: None)
    cache = build_zbw_grid_cache(
        'PlanetProfile.Default.Enceladus.PPEnceladus',
        zb_km_grid=[25.0, 25.25], wOcean_ppt_grid=[5.0, 20.0],
        output_path=str(tmp_path / 'x.pkl'),
        bulk_overrides={'Cuncertainty': 0.08},
        progress=False,
    )
    assert cache['zb_tol_km'] == pytest.approx(0.125)  # grid_step / 2
    assert cache['n_zb_placement_rejected'] == 0
    assert all(s is not None for s in cache['structures'])
    for s in cache['structures']:
        assert abs(s['zb_residual_km']) < 0.125


# ---------------------------------------------------------------------------
# Mechanism 2 -- the tolerance the root-find is given
# ---------------------------------------------------------------------------

def test_builder_declares_a_km_tolerance_to_the_root_find(monkeypatch,
                                                          tmp_path):
    """Every node is built with Bulk.zbTol_km = 0.4 x the invariant -- the
    same 0.4x convention ``_build_frozen_node`` already uses."""
    seen = []

    def _fake(module, Tb_K, **kwargs):
        seen.append(kwargs['bulk_overrides'])
        return _struct(kwargs['bulk_overrides']['zb_approximate_km'])

    monkeypatch.setattr(cache_builder, 'build_single_structure', _fake)
    monkeypatch.setattr(cache_builder, '_save_cache_atomic',
                        lambda cache, path: None)
    build_zbw_grid_cache(
        'PlanetProfile.Default.Enceladus.PPEnceladus',
        zb_km_grid=[25.0, 25.25], wOcean_ppt_grid=[5.0, 20.0],
        output_path=str(tmp_path / 'x.pkl'),
        bulk_overrides={'Cuncertainty': 0.08},
        progress=False,
    )
    assert seen
    for bo in seen:
        assert bo['zbTol_km'] == pytest.approx(0.4 * 0.125)


def test_bulk_default_leaves_every_other_body_on_the_legacy_path():
    """zbTol_km must default to None, or tightening Enceladus would silently
    re-tune every body through the shared TfreezeRes_K."""
    assert BulkSubstruct().zbTol_km is None


def test_derived_xtol_converts_km_to_K_through_the_measured_slope():
    calls = []

    def solver(T):
        calls.append(T)
        return -100.0 * (T - 271.0)   # 100 km of zb per K

    xtol = LayerPropagators._zbDrivenXtol(solver, 271.0, 272.0, 0.05)
    assert calls == [271.0, 272.0]      # measured on the established bracket
    assert xtol == pytest.approx(
        0.05 / (LayerPropagators._ZB_XTOL_SAFETY * 100.0))
    # 5e-5 K against TfreezeRes_K's 0.05 K: three orders of magnitude, which
    # is the size of the defect.
    assert xtol < 0.05


def test_derived_xtol_floors_on_a_degenerate_bracket():
    assert LayerPropagators._zbDrivenXtol(
        lambda T: 1.0, 271.0, 272.0, 0.05) == LayerPropagators._ZB_XTOL_FLOOR_K


# ---------------------------------------------------------------------------
# Mechanism 1 -- which iterate GetIceShellTFreeze returns
#
# These drive the REAL GetIceShellTFreeze. Only its two external dependencies
# are faked: GetTfreeze (the phase-diagram bracket) and IceLayers (the column
# build). The IceShellResidual class, the retention rule and the return path
# are the production ones.
# ---------------------------------------------------------------------------

def _fake_planet(zb_target_km, zbTol_km=None):
    return types.SimpleNamespace(
        Bulk=types.SimpleNamespace(Tb_K=None,
                                   zb_approximate_km=float(zb_target_km),
                                   zbTol_km=zbTol_km),
        Do=types.SimpleNamespace(ICEIh_THICKNESS=True),
        Ocean=types.SimpleNamespace(meltEOS=object()),
        TfreezeRes_K=0.05,
        TfreezeLower_K=260.0,
        TfreezeUpper_K=273.0,
        PfreezeLower_MPa=0.1,
        PfreezeUpper_MPa=20.0,
    )


@pytest.fixture
def bracket(monkeypatch):
    """GetTfreeze returns the two ends of the measured Enceladus bracket."""
    def _fake_GetTfreeze(meltEOS, P_MPa, TLower_K, TUpper_K=None, TRes_K=None):
        return R6_T_SEQ[1] if (TRes_K or 0) < 0 else R6_T_SEQ[0]

    monkeypatch.setattr(LayerPropagators, 'GetTfreeze', _fake_GetTfreeze)


class _IceResult:
    PbI_MPa = 1.0

    def __init__(self, zb_km):
        self.zb_km = zb_km


def test_solver_returns_the_minimum_residual_iterate(monkeypatch, bracket):
    """Replays r6's measured (zb = 25 km, w = 20 ppt) iterate sequence through
    the production function with a scripted root-finder.

    Last iterate misses by 2.965 km (23x the 0.125 km budget); the best misses
    by 0.116 km. The returned Planet must carry the best one's Tb_K.
    """
    seq = {'i': 0}

    def _fake_IceLayers(Planet, Params):
        i = seq['i']
        seq['i'] += 1
        return _IceResult(R6_ZB_SEQ[i])

    def _fake_GetZero(f, bracket=None, xtol=None):
        for T in R6_T_SEQ:
            f(T)
        return types.SimpleNamespace(root=R6_T_SEQ[-1], function_calls=len(R6_T_SEQ))

    monkeypatch.setattr(LayerPropagators, 'IceLayers', _fake_IceLayers)
    monkeypatch.setattr(LayerPropagators, 'GetZero', _fake_GetZero)

    out = LayerPropagators.GetIceShellTFreeze(_fake_planet(25.0), object())

    assert out.Bulk.Tb_K == pytest.approx(R6_T_SEQ[3])   # the BEST iterate
    assert out.Bulk.Tb_K != pytest.approx(R6_T_SEQ[-1])  # not the LAST
    assert out.Do.ICEIh_THICKNESS is True
    # The regression this test exists for, in the units the invariant uses:
    assert abs(25.0 - R6_ZB_SEQ[-1]) > 23 * 0.125
    assert abs(25.0 - R6_ZB_SEQ[3]) < 0.125


def test_legacy_path_keeps_TfreezeRes_K_as_the_root_find_tolerance(
        monkeypatch, bracket):
    """Byte-identity guard for Titan/Europa: with zbTol_km unset, the xtol
    handed to the root-finder is still the shared phase-diagram step, and the
    zb early-exit is still the historical 0.01 km."""
    seen = {}

    def _fake_IceLayers(Planet, Params):
        return _IceResult(25.0 + 120.0 * (R6_T_SEQ[3] - Planet.Bulk.Tb_K))

    def _fake_GetZero(f, bracket=None, xtol=None):
        seen['xtol'] = xtol
        seen['bracket'] = bracket
        f(bracket[0])
        return types.SimpleNamespace(root=bracket[0], function_calls=1)

    monkeypatch.setattr(LayerPropagators, 'IceLayers', _fake_IceLayers)
    monkeypatch.setattr(LayerPropagators, 'GetZero', _fake_GetZero)

    LayerPropagators.GetIceShellTFreeze(_fake_planet(25.0), object())
    assert seen['xtol'] == 0.05          # == Planet.TfreezeRes_K, unchanged
    assert seen['bracket'] == [R6_T_SEQ[0], R6_T_SEQ[1]]


def test_zb_driven_path_tightens_the_tolerance_by_three_orders(monkeypatch,
                                                               bracket):
    """With zbTol_km declared, the xtol is derived from the bracket slope."""
    seen = {}

    def _fake_IceLayers(Planet, Params):
        return _IceResult(25.0 + 120.0 * (R6_T_SEQ[3] - Planet.Bulk.Tb_K))

    def _fake_GetZero(f, bracket=None, xtol=None):
        seen['xtol'] = xtol
        f(bracket[0])
        return types.SimpleNamespace(root=bracket[0], function_calls=1)

    monkeypatch.setattr(LayerPropagators, 'IceLayers', _fake_IceLayers)
    monkeypatch.setattr(LayerPropagators, 'GetZero', _fake_GetZero)

    LayerPropagators.GetIceShellTFreeze(_fake_planet(25.0, zbTol_km=0.05),
                                        object())
    assert seen['xtol'] == pytest.approx(
        0.05 / (LayerPropagators._ZB_XTOL_SAFETY * 120.0))
    assert seen['xtol'] < 0.05 / 1000.0


def test_zb_driven_path_reaches_production_tolerance_with_real_brentq(
        monkeypatch, bracket):
    """End-to-end with scipy's actual root-finder on a smooth, convex zb(Tb)
    of the measured Enceladus scale: legacy xtol misses the 0.125 km budget,
    the derived xtol meets it."""
    def _fake_IceLayers(Planet, Params):
        # zb(T) = 25 km at the root, |d zb/d T| ~ 120 km/K, convex like the
        # real melting curve.
        dT = R6_T_SEQ[3] - Planet.Bulk.Tb_K
        return _IceResult(25.0 + 120.0 * dT + 40.0 * dT ** 2)

    monkeypatch.setattr(LayerPropagators, 'IceLayers', _fake_IceLayers)

    out = LayerPropagators.GetIceShellTFreeze(
        _fake_planet(25.0, zbTol_km=0.05), object())
    achieved = 25.0 + 120.0 * (R6_T_SEQ[3] - out.Bulk.Tb_K) \
        + 40.0 * (R6_T_SEQ[3] - out.Bulk.Tb_K) ** 2
    assert abs(achieved - 25.0) < 0.125
