"""Opt-in repair of the rigid Van Hoolst et al. (2008) ocean libration branch.

``librations(..., ocean=True, rigid=True, formulation='vh2008_repaired')``
repairs, TOGETHER (they must not land piecemeal; freeze-doc B2 /
k2-campaign plan Physics step 2, gate C7/P1b):

* K_int: VH08 Eq. (16) carries (4 pi G / 5)(8 pi / 15); the shipped code
  drops the 8 pi / 15. Repaired as VH08 Eq. (9) / VH13 Eq. (37) with the
  exact region-1 moment difference (Bc - Ac) - (Bi - Ai)_pressure.
* the ocean-pressure figure moments: linearized (8 pi / 15) rho beta r^5
  (VH08 Eq. 14/54) replaced by the exact ellipsoid moment rho_o * f,
  f = (4 pi / 15) a b c (a^2 - b^2), consistent with ``moi_ellps``.

The DEFAULT (``'vh2008_shipped'``) is unchanged and pinned here, because the
frozen Cassini campaign applies a +0.8 % correction calibrated against it.
"""
import pickle
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tests"))

from PlanetProfile.Gravity import pyALMA3Updated  # noqa: E402
from PlanetProfile.Gravity.Gravity import (  # noqa: E402
    ALMA_TIME_UNIT_S,
    LibrationModelInputs,
)
from PlanetProfile.Gravity.Librations import (  # noqa: E402
    compute_eccentricities,
    librations,
    radial2abc,
)
from PlanetProfile.Utilities.defineStructs import Constants  # noqa: E402
from librations_test import (  # noqa: E402
    CACHE_PATH as SMOKE_CACHE,
    ECCENTRICITY,
    OMEGA_VH16_RADPS,
    _enceladus_b2prime_structure,
    _vh16_interior,
)

PROD_CACHE = (
    REPO
    / "PlanetProfile/Test/mcmc_results/Enceladus/Cassini_isostasy"
    / "enceladus_seawater_zbw_grid_2d.pkl"
)
REPAIRED = dict(formulation="vh2008_repaired")


def _cases():
    r, rho = _vh16_interior(20.0)
    yield "vh16", r, rho, OMEGA_VH16_RADPS, ECCENTRICITY
    yield ("b2prime", *_enceladus_b2prime_structure())
    for zb in (10.0, 25.0, 40.0):
        r, rho = _vh16_interior(zb)
        yield f"vh16_zb{zb:g}", r, rho, OMEGA_VH16_RADPS, ECCENTRICITY


# --- default path is unchanged -------------------------------------------

# repr() values recorded on the pre-change code (genai @ 132ead8b).
SHIPPED_PINS_M = {
    "vh16": 524.3804491164386,
    "b2prime_hydrostatic": 399.586798187253,
    "b2prime_surface_only": 442.5725181368896,
    "b2prime_C2_0.5": 432.01560840065963,
    "b2prime_C2_1.0": 421.5376281578341,
}


def test_default_is_the_shipped_formulation_exactly():
    """Omitting ``formulation`` and passing 'vh2008_shipped' are identical."""
    for _, r, rho, om, e in _cases():
        for kw in ({}, dict(H22_obs_m=857.0),
                   dict(H22_obs_m=857.0, compensation_C2=1.0)):
            assert librations(r, rho, om, e, **kw) == librations(
                r, rho, om, e, formulation="vh2008_shipped", **kw)


def test_shipped_values_are_pinned_to_the_pre_change_code():
    r, rho = _vh16_interior(20.0)
    got = {"vh16": librations(r, rho, OMEGA_VH16_RADPS, ECCENTRICITY)}
    r, rho, om, e = _enceladus_b2prime_structure()
    got["b2prime_hydrostatic"] = librations(r, rho, om, e)
    got["b2prime_surface_only"] = librations(r, rho, om, e, H22_obs_m=857.0)
    for c2 in (0.5, 1.0):
        got[f"b2prime_C2_{c2}"] = librations(
            r, rho, om, e, H22_obs_m=857.0, compensation_C2=c2)
    for key, val in SHIPPED_PINS_M.items():
        assert got[key] == pytest.approx(val, rel=1e-12), key


def test_unknown_formulation_raises():
    r, rho, om, e = _enceladus_b2prime_structure()
    with pytest.raises(ValueError):
        librations(r, rho, om, e, formulation="vh2008")


def test_formulation_does_not_touch_other_branches():
    r, rho, om, e = _enceladus_b2prime_structure()
    for kw in (dict(ocean=False), dict(rigid=False, y1=np.full(3, 1e-3))):
        assert librations(r, rho, om, e, **kw) == librations(
            r, rho, om, e, **REPAIRED, **kw)


# --- the repair itself ----------------------------------------------------

def test_repaired_rigid_is_the_rigid_limit_of_the_deformable_branch():
    """Gate P1b, analytic form: VH13 with y1 = 0 (no tidal deformation)
    must reduce to the repaired VH08 rigid branch.

    The deformable branch builds K_int from (Bi - Ai) + (Bob - Aob) with the
    exact moments (VH13 Eq. 37), so with every k2_j = 0 it is the rigid
    problem with consistent coefficients. The shipped rigid branch misses
    this by 2.8-3.6 %; the repaired one matches to rounding.
    """
    for name, r, rho, om, e in _cases():
        rigid_rep = librations(r, rho, om, e, **REPAIRED)
        rigid_ship = librations(r, rho, om, e)
        deform_y0 = librations(r, rho, om, e, rigid=False, y1=np.zeros(3))
        assert deform_y0 == pytest.approx(rigid_rep, rel=1e-12), name
        assert abs(deform_y0 / rigid_ship - 1) > 0.01, name


def test_repair_size_at_the_vh16_marker():
    """Total repaired/shipped at the VH16 marker is +3.57 %. (Measured
    separately: the 8pi/15 K_int piece alone is +1.88 %, matching the B2
    adjudication's ~+1.9 %; the exact pressure moments supply the rest.)"""
    r, rho = _vh16_interior(20.0)
    ship = librations(r, rho, OMEGA_VH16_RADPS, ECCENTRICITY)
    rep = librations(r, rho, OMEGA_VH16_RADPS, ECCENTRICITY, **REPAIRED)
    assert rep / ship - 1 == pytest.approx(0.03572, abs=2e-4)


def test_repaired_vh16_yellow_marker_lands_in_published_band():
    """VH16 Figs. 2-3 yellow marker; observed 2-sigma band 466-590 m.

    Repaired: 543.11 m (shipped 524.38 m). The band does not discriminate.
    """
    r, rho = _vh16_interior(20.0)
    lib_m = librations(r, rho, OMEGA_VH16_RADPS, ECCENTRICITY, **REPAIRED)
    assert 466.0 <= lib_m <= 590.0
    assert lib_m == pytest.approx(543.110, abs=0.01)


def test_repaired_libration_decreases_with_shell_thickness():
    amps = [librations(*_vh16_interior(z), OMEGA_VH16_RADPS, ECCENTRICITY,
                       **REPAIRED) for z in np.arange(5.0, 46.0, 5.0)]
    assert np.all(np.diff(amps) < 0.0)


# --- elastic branch vs repaired rigid on cache nodes ----------------------

def _alma_node(structure, mu_scale=1.0):
    r = np.asarray(structure["r_m"], float)
    rho = np.asarray(structure["rho"], float)
    mu = np.asarray(structure["mu_Pa"], float) * mu_scale
    eta = np.asarray(structure["eta_Pa_base"], float)
    ph = np.asarray(structure["phases"])
    rheo = np.where(ph == 0, "newton", "elastic").tolist()
    model = pyALMA3Updated.build_model(r, rho, mu, eta, rheo,
                                       np.zeros((r.size, 2)), norm=False)
    om = float(structure["omega"])
    _, _, _, rf = pyALMA3Updated.love_numbers(
        model, [2], [(2 * np.pi / om) / ALMA_TIME_UNIT_S], "tidal",
        "periodic", 4, verbose=False, internal=True)
    M, R = float(structure["Mtot_kg"]), float(structure["R_body_m"])
    y1 = np.real(rf[:, 0]) * M / R / (Constants.G * M / R ** 2)
    rr, rrho, oi, ry = LibrationModelInputs(r, ph, rho, y1)
    kw = dict(omega=om, ecc=float(structure["eccentricity"]), ocean=True,
              ocean_idx=oi)
    return {
        "shipped": librations(rr, rrho, **kw),
        "repaired": librations(rr, rrho, **REPAIRED, **kw),
        "elastic": librations(rr, rrho, rigid=False, y1=ry, **kw),
    }


@lru_cache(maxsize=None)
def _smoke_node(mu_scale=1.0):
    if not SMOKE_CACHE.exists():
        pytest.skip("smoke cache absent")
    with SMOKE_CACHE.open("rb") as stream:
        cache = pickle.load(stream)
    i = int(np.argmin(np.abs(np.asarray(cache["Tb_K_grid"]) - 272.46)))
    return _alma_node(cache["structures"][i], mu_scale)


@lru_cache(maxsize=None)
def _prod_node(index, mu_scale=1.0):
    if not PROD_CACHE.exists():
        pytest.skip("production cache absent")
    with PROD_CACHE.open("rb") as stream:
        cache = pickle.load(stream)
    return _alma_node(cache["structures"][index], mu_scale)


def test_repaired_branch_has_the_vh16_physical_direction():
    """The strict-xfail expectation of librations_test.py, met by the
    repaired branch: elasticity LOWERS libration (VH16; VH13 sec. 5).

    Measured on the adjudicated Tb = 272.46 K smoke node: elastic/repaired
    = 0.97380 (-2.62 %, the adjudicated real elastic reduction), while the
    shipped rigid value is unchanged at 0.110905164 deg.
    """
    v = _smoke_node()
    assert v["elastic"] < v["repaired"]
    assert v["elastic"] / v["repaired"] == pytest.approx(0.973798, rel=1e-4)
    # The shipped pin of librations_test.py is untouched by this change.
    assert np.degrees(v["shipped"] / 252.1e3) == pytest.approx(
        0.110905164, rel=1e-6)


@pytest.mark.skipif(not PROD_CACHE.exists(), reason="production cache absent")
@pytest.mark.parametrize("index", [1220, 1540])
def test_repaired_physical_direction_on_production_nodes(index):
    """zb 24.75 / 26.75 km at w = 3.455 ppt: elastic is 1.9-2.1 % BELOW
    repaired rigid (shipped rigid sits 0.9 % below elastic)."""
    v = _prod_node(index)
    assert v["elastic"] < v["repaired"]
    assert -0.025 < v["elastic"] / v["repaired"] - 1 < -0.015


def test_gate_p1b_stiff_elastic_limit_reproduces_repaired_rigid():
    """Gate P1b: the deformable branch with ALMA y1 converges to the
    repaired rigid result as the solid layers stiffen (mu x 100 -> -0.027 %,
    mu x 1e4 -> -3e-6), and NOT to the shipped rigid result (+3.5 %)."""
    stiff = _smoke_node(1e4)
    assert stiff["elastic"] == pytest.approx(stiff["repaired"], rel=1e-5)
    assert abs(stiff["elastic"] / stiff["shipped"] - 1) > 0.03
    mid = _smoke_node(1e2)
    assert abs(mid["elastic"] / mid["repaired"] - 1) < 5e-4


# --- general layer stacks (ocean_idx > 1, multi-layer shell) -------------

def _vh16_mass_kg(zb_km=25.0):
    r, rho = _vh16_interior(zb_km)
    shells = np.diff(np.concatenate([[0.0], r]) ** 3)
    return 4 * np.pi / 3 * float(np.sum(rho * shells))


def _mass_conserving(r, rho, fill_idx, mass_kg):
    """Set rho[fill_idx] so the stack has the given total mass."""
    rho = np.array(rho, float)
    shells = 4 * np.pi / 3 * np.diff(np.concatenate([[0.0], r]) ** 3)
    rho[fill_idx] = 0.0
    rho[fill_idx] = (mass_kg - np.sum(rho * shells)) / shells[fill_idx]
    return rho


def _split_core_stack():
    """Two-layer interior below the ocean, so ocean_idx = 2."""
    r, rho = _vh16_interior(25.0)
    r4 = np.array([0.6 * r[0], r[0], r[1], r[2]])
    rho4 = _mass_conserving(r4, [1.3 * rho[0], 0.0, rho[1], rho[2]], 1,
                            _vh16_mass_kg())
    return r4, rho4, 2


def _two_layer_shell_stack():
    """Two-layer ice shell above the ocean (ocean_idx = 1)."""
    r, rho = _vh16_interior(25.0)
    r4 = np.array([r[0], r[1], r[2] - 10e3, r[2]])
    rho4 = _mass_conserving(r4, [0.0, rho[1], 940.0, 920.0], 0,
                            _vh16_mass_kg())
    return r4, rho4, 1


def test_repaired_moments_close_the_whole_body_identity():
    """(Bs - As) + Bsp + (Bc - Ac) - Bip = sum(B - A) exactly when the
    pressure moments are the exact ellipsoid moments (the shipped linearized
    ones miss it by ~0.5 %), for every stack shape."""
    from PlanetProfile.Gravity.Librations import moi_ellps
    stacks = [(*_vh16_interior(25.0), 1), _split_core_stack(),
              _two_layer_shell_stack()]
    for r, rho, oi in stacks:
        ep, eq = compute_eccentricities(r, rho, OMEGA_VH16_RADPS)
        a, b, c = radial2abc(r, ep, eq)
        A, B, _ = moi_ellps(a, b, c, rho)
        f = 4 / 15 * np.pi * a * b * c * (a ** 2 - b ** 2)
        lhs = ((B[oi + 1:] - A[oi + 1:]).sum() + rho[oi] * f[oi]
               + (B[:oi] - A[:oi]).sum() - rho[oi] * f[oi - 1])
        assert lhs == pytest.approx((B - A).sum(), rel=1e-13)


def test_repaired_rigid_limit_with_a_layered_interior():
    """ocean_idx = 2: the repaired branch indexes the ocean-bottom moment at
    ocean_idx - 1 and still equals the y1 = 0 deformable branch exactly."""
    r, rho, oi = _split_core_stack()
    kw = dict(omega=OMEGA_VH16_RADPS, ecc=ECCENTRICITY, ocean_idx=oi)
    rep = librations(r, rho, formulation="vh2008_repaired", **kw)
    d0 = librations(r, rho, rigid=False, y1=np.zeros(r.size), **kw)
    assert d0 == pytest.approx(rep, rel=1e-12)
    assert abs(rep / librations(r, rho, **kw) - 1) > 0.01


def test_repaired_rigid_with_a_two_layer_shell():
    """The repaired outer factor sums rho_j (beta_top - beta_bot) over every
    shell layer (VH08 Eq. 12). The deformable branch uses the shell-base
    density for the whole shell, so it differs here, at 1.6e-5 for a
    940/920 kg/m^3 shell (measured 2026-11-03). With a uniform shell the
    two agree exactly (test above)."""
    r, rho, oi = _two_layer_shell_stack()
    kw = dict(omega=OMEGA_VH16_RADPS, ecc=ECCENTRICITY, ocean_idx=oi)
    rep = librations(r, rho, formulation="vh2008_repaired", **kw)
    d0 = librations(r, rho, rigid=False, y1=np.zeros(r.size), **kw)
    assert 0 < abs(d0 / rep - 1) < 1e-4
    # Collapsing the shell to one density recovers the exact agreement.
    rho_u = _mass_conserving(r, [0.0, rho[1], 930.0, 930.0], 0,
                             _vh16_mass_kg())
    rep_u = librations(r, rho_u, formulation="vh2008_repaired", **kw)
    d0_u = librations(r, rho_u, rigid=False, y1=np.zeros(r.size), **kw)
    assert d0_u == pytest.approx(rep_u, rel=1e-12)


# --- H&M observed-figure paths under the repaired formulation -------------

def test_h22_paths_run_under_repaired_formulation():
    r, rho, om, e = _enceladus_b2prime_structure()
    base = librations(r, rho, om, e, **REPAIRED)
    surf = librations(r, rho, om, e, H22_obs_m=857.0, **REPAIRED)
    c2 = {C2: librations(r, rho, om, e, H22_obs_m=857.0,
                         compensation_C2=C2, **REPAIRED)
          for C2 in (0.0, 0.5, 1.0)}
    assert all(np.isfinite(v) and v > 0 for v in (base, surf, *c2.values()))
    # C2 -> 0 is surface-only exactly, as in the shipped formulation.
    assert c2[0.0] == surf
    # Same ordering as the adjudicated bracket: hyd < C2=1 < C2=0.5 < C2=0.
    assert base < c2[1.0] < c2[0.5] < c2[0.0]
    # Passing the hydrostatic H22 is a no-op.
    ep, eq = compute_eccentricities(r, rho, om)
    a, b, _ = radial2abc(r, ep, eq)
    noop = librations(r, rho, om, e, H22_obs_m=(a[-1] - b[-1]) / 6, **REPAIRED)
    assert noop == pytest.approx(base, rel=1e-12)


def test_repaired_shell_bracket_is_exactly_the_drho_sum():
    """Regression pin only: this re-assembles the implementation's own
    algebra. The independent check is the y1 = 0 deformable comparison.

    Under the repair Ks/(3 n^2) = rho_ice f_top s + (rho_ocean - rho_ice)
    f_base EXACTLY (H&M Eq. 19; the shipped branch meets it only to ~1 %
    because Bsp_Asp is linearized), Kc uses (rho_core - rho_ocean) f_core,
    and K_int is VH13 Eq. (37). A hand-assembled solve with those
    coefficients must reproduce the repaired surface-only result.
    """
    from PlanetProfile.Gravity.Librations import G, moi_ellps
    r, rho, om, e = _enceladus_b2prime_structure()
    ep, eq = compute_eccentricities(r, rho, om)
    a, b, c = radial2abc(r, ep, eq)
    Ai, Bi, Ci = moi_ellps(a, b, c, rho)
    f = 4 / 15 * np.pi * a * b * c * (a ** 2 - b ** 2)
    beta = (a - b) / a
    s = 857.0 / ((a[-1] - b[-1]) / 6)
    Ks = 3 * om ** 2 * (rho[2] * f[2] * s + (rho[1] - rho[2]) * f[1])
    BmA_int = (rho[0] - rho[1]) * f[0]
    assert Bi[0] - Ai[0] == pytest.approx(rho[0] * f[0], rel=1e-13)
    Kc = 3 * om ** 2 * BmA_int
    K = 4 * np.pi * G / 5 * (rho[1] * beta[1] + rho[2] * (beta[2] - beta[1])) * BmA_int
    Cs, Cc = Ci[2], Ci[0]
    A = np.array([[-om ** 2 * Cs + Ks + 2 * K, -2 * K],
                  [-2 * K, -om ** 2 * Cc + Kc + 2 * K]])
    expected = abs(np.linalg.solve(A, [2 * e * Ks, 2 * e * Kc])[0] * r[-1])
    got = librations(r, rho, om, e, H22_obs_m=857.0, **REPAIRED)
    assert got == pytest.approx(expected, rel=1e-11)
