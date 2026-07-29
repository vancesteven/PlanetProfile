"""Validation of the CoolProp NH3-H2O ocean EOS
(Thermodynamics/NH3/NH3Props.py + the HydroEOS 'NH3' branch).

Anchors:
- density/Cp vs the Melinder (2010) aqueous-ammonia polynomials
  (CoolProp INCOMP::MAM) at low pressure;
- pure-water-limit density vs SeaFreeze water1 across ocean pressures;
- LIQUIDUS vs the Melinder (2010) experimental freezing curve at 1 bar
  (the Melinder-anchored activity model, scientific review 2026-07-28)
  and the shifted phase diagram (freezing point drops with w; ices are
  pure H2O); negative-excess sign guard against the legacy CoolProp
  mixture activity;
- end-to-end GetOceanEOS('NH3', ...) construction with physical
  property profiles.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

pytest.importorskip('CoolProp')

from PlanetProfile.Thermodynamics.NH3.NH3Props import (
    NH3Props, _molefrac_NH3, _ln_aw_anchored, muLiquidusCurve_K,
    pureMeltingCurve_K)


def test_anchored_negative_excess():
    # Regression guard against ever reverting to the CoolProp mixture
    # activity: NH3-H2O requires gamma_w < 1 (negative excess), i.e.
    # ln a_w < ln x_w, at every composition, pressure, and temperature.
    for X in (0.01, 0.05, 0.10, 0.15, 0.20):
        u = _molefrac_NH3(X * 1000.0)
        for P in (0.1, 100.0, 500.0, 1000.0):
            for T in (250.0, 270.0, 300.0):
                assert _ln_aw_anchored(P, T, u, 270.0) \
                    < np.log(1.0 - u), (X, P, T)


def test_melinder_anchor_1bar():
    # The binding accuracy test (would have caught the 3.5 K CoolProp
    # bias): 1-bar liquidus depression within 0.05 K of the Melinder
    # (2010) experimental freezing curve, SELF-referenced to Melinder's
    # own zero (273.1218 K != SeaFreeze's 273.1527 K).
    from CoolProp.CoolProp import PropsSI
    Tf0 = PropsSI('T_FREEZE', 'P', 101325, 'Q', 0, 'INCOMP::MAM[1e-7]')
    Tm, _ = pureMeltingCurve_K(np.array([0.1]))
    for X in (0.02, 0.03, 0.05, 0.075, 0.10, 0.125, 0.15):
        Tf = PropsSI('T_FREEZE', 'P', 101325, 'Q', 0,
                     f'INCOMP::MAM[{X}]')
        Tliq, TmP = muLiquidusCurve_K(np.array([0.1]), X * 1000.0)
        dep = float(TmP[0] - Tliq[0])
        assert abs(dep - (Tf0 - Tf)) < 0.05, (X, dep, Tf0 - Tf)


def test_pure_and_colligative_limits():
    # w = 0 reproduces the pure melting curve; at trace w the anchored
    # model must match the exact ideal-colligative solve (gamma -> 1 as
    # u^2, so the excess contribution is ~1e-6 K there).
    import seafreeze.seafreeze as sf
    P = np.array([0.1, 100.0, 900.0])
    Tliq, Tm = muLiquidusCurve_K(P, 0.0)
    assert np.allclose(Tliq, Tm, atol=1e-9)
    w = 0.1  # ppt
    u = _molefrac_NH3(w)
    Tliq, Tm = muLiquidusCurve_K(np.array([0.1]), w)
    dep = float(Tm[0] - Tliq[0])

    def _dG(T):
        pts = np.array([(0.1, T)], dtype='f,f').astype(object)
        Gi = float(np.ravel(sf.getProp(pts, 'Ih').G)[0])
        Gw = float(np.ravel(sf.getProp(pts, 'water1').G)[0])
        return (Gi - Gw) * 0.018015

    lo, hi = float(Tm[0]) - 1.0, float(Tm[0])
    for _ in range(40):  # ideal-colligative bisection: ln a_w = ln x_w
        mid = 0.5 * (lo + hi)
        if _dG(mid) - 8.31446 * mid * np.log(1.0 - u) > 0:
            hi = mid
        else:
            lo = mid
    dep_ideal = float(Tm[0]) - 0.5 * (lo + hi)
    assert abs(dep - dep_ideal) < 1e-3, (dep, dep_ideal)


def test_molefrac():
    # 10 wt% NH3 -> x ~ 0.105
    assert _molefrac_NH3(100.0) == pytest.approx(0.1052, abs=0.001)


def test_props_vs_melinder_low_pressure():
    from CoolProp.CoolProp import PropsSI
    P = np.array([1.0, 2.0, 3.0])
    T = np.array([255.0, 265.0, 275.0, 285.0])
    for w in (50.0, 100.0, 200.0):
        rho, Cp, alpha, VP, KS = NH3Props(P, T, w)
        X = w / 1000.0
        errs = []
        for j, Tj in enumerate(T):
            try:
                r_m = PropsSI('D', 'T', Tj, 'P', 2e6, f'INCOMP::MAM[{X}]')
            except Exception:
                continue  # outside Melinder validity
            errs.append(abs(rho[1, j] / r_m - 1))
        assert errs, 'no Melinder anchor points evaluated'
        assert max(errs) < 0.02, (w, errs)
        assert np.all(np.isfinite(rho)) and np.all(rho > 800)
        assert np.all(Cp > 2000) and np.all(Cp < 6000)
        assert np.all(VP > 0.8) and np.all(VP < 4.0)   # km/s
        assert np.all(KS > 0.5) and np.all(KS < 20.0)  # GPa


def test_pure_limit_vs_seafreeze():
    import seafreeze.seafreeze as sf
    P = np.array([10.0, 100.0, 500.0, 1000.0])
    T = np.array([260.0, 280.0, 300.0])
    rho, *_ = NH3Props(P, T, 1.0)  # 1 ppt ~ pure water
    for i, Pi in enumerate(P):
        for j, Tj in enumerate(T):
            if Tj < 273.0 and Pi > 200.0:
                # pure water is SOLID here (ice V/VI): CoolProp's liquid
                # branch is invalid and NH3Props fills a smooth
                # metastable placeholder — never sampled by profiles
                # (the phase function marks it ice). Skip the anchor.
                continue
            pts = np.array([(Pi, Tj)], dtype='f,f').astype(object)
            r_sf = float(np.ravel(sf.getProp(pts, 'water1').rho)[0])
            assert abs(rho[i, j] / r_sf - 1) < 0.02, (Pi, Tj)
    # the filled grid must stay monotone in P on every isotherm (no
    # spline-poisoning spurious roots)
    assert np.all(np.diff(rho, axis=0) > -1e-6)


def test_density_increases_with_ammonia_dilution_inverse():
    # NH3 is lighter than water: at fixed (P, T), more NH3 -> lower rho
    P = np.array([50.0, 100.0])
    T = np.array([260.0, 270.0])
    rho_lo, *_ = NH3Props(P, T, 20.0)
    rho_hi, *_ = NH3Props(P, T, 150.0)
    assert np.all(rho_hi < rho_lo)


def test_get_ocean_eos_end_to_end():
    from PlanetProfile.Thermodynamics.HydroEOS import GetOceanEOS
    P = np.linspace(0.5, 400.0, 60)
    T = np.linspace(240.0, 300.0, 50)
    eos = GetOceanEOS('NH3', 150.0, P, T, None)

    def _s(v):
        return float(np.ravel(v)[0])

    rho = _s(eos.fn_rho_kgm3(np.array([100.0]), np.array([260.0])))
    assert 900 < rho < 1150
    Cp = _s(eos.fn_Cp_JkgK(np.array([100.0]), np.array([260.0])))
    assert 3000 < Cp < 6000
    VP, VS, KS, GS = eos.fn_Seismic(np.array([100.0]), np.array([260.0]))
    assert 1.0 < _s(VP) < 3.0 and _s(VS) == 0.0

    # Freezing-point depression: straddle just below the pure-water
    # melting point — pure water frozen, NH3 ocean still liquid.
    # Anchored-model depression at 150 ppt is ~23 K at low P; probe at
    # a conservative 10 K below the pure freezing point.
    dT = 20.0
    eos0 = GetOceanEOS('PureH2O', 0.0, P, T, None)
    Ptest = np.array([10.0])

    def _ph(e, t_):
        return int(np.ravel(e.fn_phase(Ptest, np.array([t_])))[0])

    # find pure-water freezing T at 10 MPa on the grid
    Tscan = np.linspace(255., 275., 201)
    ph0 = np.array([_ph(eos0, t_) for t_ in Tscan])
    Tfreeze0 = Tscan[np.argmax(ph0 == 0)]  # first liquid
    Tmid = Tfreeze0 - 0.5 * dT
    assert _ph(eos0, Tmid) != 0   # pure: frozen
    assert _ph(eos, Tmid) == 0    # NH3: liquid


def test_mu_liquidus_dilute_colligative():
    # 5 wt% at 1 bar: Melinder gives 6.05 K (self-referenced); the
    # anchored liquidus must land on it (ideal colligative ~5.9; the
    # removed L-K polynomial gave a provably-wrong 2.8 K).
    Tliq, Tm = muLiquidusCurve_K(np.array([0.5]), 50.0)
    dT = float(Tm[0] - Tliq[0])
    assert 5.7 < dT < 6.4, dT


def test_mu_liquidus_high_pressure_ices():
    # Depression must apply under the HP ices too (V at ~620 MPa, VI at
    # ~900 MPa) and grow with pressure at fixed w (CG2010 shape factor,
    # saturating at ~1.4x by 1 GPa).
    P = np.array([100.0, 620.0, 900.0])
    Tliq, Tm = muLiquidusCurve_K(P, 150.0)
    dT = Tm - Tliq
    assert np.all(np.isfinite(dT)), dT
    assert np.all(dT > 20.0) and np.all(dT < 36.0), dT
    assert dT[2] > dT[0], dT  # deeper -> larger depression
    # P-shape sanity (tanh saturation bound): dep(1 GPa)/dep(1 bar)
    Tl2, Tm2 = muLiquidusCurve_K(np.array([0.1, 1000.0]), 100.0)
    ratio = float((Tm2[1] - Tl2[1]) / (Tm2[0] - Tl2[0]))
    assert 1.2 <= ratio <= 1.5, ratio


def test_liquidus_solver_gates_machine_b():
    """Machine B defect-report gate 1, tightened for the analytic
    activity (2026-07-28): dTliq(P) smooth within the ice-Ih segment
    (< 0.05 K off its local median — no repairs exist to hide behind),
    monotone increasing in w, and the corrected target table
    ~1.1/3.6/6.2/24 K at w=10/30/50/150, 100 MPa (Melinder-anchored;
    the old ~1.1/3/5/15 targets encoded the wrong-sign CoolProp
    activity)."""
    P = np.arange(10.0, 200.0, 8.0)  # ice-Ih segment only
    Tliq, Tm = muLiquidusCurve_K(P, 30.0)
    dep = Tm - Tliq
    assert np.all(np.isfinite(dep))
    med = np.array([np.median(dep[max(0, k-2):k+3])
                    for k in range(dep.size)])
    assert np.max(np.abs(dep - med)) <= 0.05
    assert dep.max() < 5.0
    deps = []
    for w, target, tol in ((10., 1.13, 0.15), (50., 6.21, 0.3),
                           (150., 24.25, 0.8)):
        Tl, T0 = muLiquidusCurve_K(np.array([100.0]), w)
        d = float(T0[0] - Tl[0])
        assert abs(d - target) < tol, (w, d)
        deps.append(d)
    assert np.all(np.diff(deps) > 0)  # strictly increasing in w


def test_campaign_rectangle_buildable():
    # Domain guard for the Titan NH3 campaign: every corner of the
    # provisional Phase-1 rectangle Tb in [248, 257] K x w in
    # [30, 100] ppt must lie on the ice-Ih liquidus branch (below the
    # 1-bar liquidus, above the liquidus at the Ih/III corner
    # ~210 MPa), so a shell + ocean both exist.
    for w in (30.0, 100.0):
        Tl, Tm = muLiquidusCurve_K(np.array([0.1, 209.0]), w)
        Tb_max, Tb_min = float(Tl[0]), float(Tl[1])
        for Tb in (248.0, 257.0):
            assert Tb_min < Tb < Tb_max, (w, Tb, Tb_min, Tb_max)


def test_legacy_coolprop_mode_retained():
    # The legacy CoolProp-activity solver stays importable/selectable
    # for regression comparison; it under-depresses (~3.2 K at w=30,
    # 10 MPa vs 3.5 K anchored).
    from PlanetProfile.Thermodynamics.NH3.NH3Props import (
        _muLiquidusCoolProp_K)
    Tl, Tm = muLiquidusCurve_K(np.array([10.0]), 30.0)
    Tl_cp, Tm_cp = _muLiquidusCoolProp_K(np.array([10.0]), 30.0)
    d_anch = float(Tm[0] - Tl[0])
    d_cp = float(Tm_cp[0] - Tl_cp[0])
    assert np.isfinite(d_cp) and 2.8 < d_cp < d_anch, (d_cp, d_anch)


def test_single_liquid_band_machine_b():
    """Machine B gates 2/4: along fixed T, the phase function must show
    exactly ONE contiguous liquid band in P (no phantom liquid cells
    inside the ice-Ih field, no interleaved liquid/ice-V stacks)."""
    from PlanetProfile.Thermodynamics.HydroEOS import GetOceanEOS
    P = np.linspace(0.5, 900.0, 160)
    T = np.linspace(240.0, 305.0, 60)
    for w in (30.0, 100.0):
        eos = GetOceanEOS('NH3', w, P, T, None)
        for Tv in (245.0, 250.0, 258.0, 266.0):
            ph = np.ravel(eos.fn_phase(P, np.full_like(P, Tv))).astype(int)
            liq = (ph == 0).astype(int)
            n_bands = int(np.sum(np.diff(liq) == 1) + (liq[0] == 1))
            assert n_bands <= 1, (w, Tv, ph.tolist())


def test_hp_ice_straddle_end_to_end():
    # At 900 MPa (ice VI regime): between the pure and NH3 liquidi the
    # NH3 ocean is liquid while pure water is ice VI; below both, the
    # NH3 case freezes to ice VI (not Ih).
    from PlanetProfile.Thermodynamics.HydroEOS import GetOceanEOS
    P = np.linspace(0.5, 1100.0, 80)
    T = np.linspace(240.0, 320.0, 60)
    eos0 = GetOceanEOS('PureH2O', 0.0, P, T, None)
    eos = GetOceanEOS('NH3', 150.0, P, T, None)

    def _ph(e, t_):
        return int(np.ravel(e.fn_phase(np.array([900.0]),
                                       np.array([t_])))[0])

    Tscan = np.linspace(280., 310., 301)
    ph0 = np.array([_ph(eos0, t_) for t_ in Tscan])
    Tf0 = Tscan[np.argmax(ph0 == 0)]
    assert _ph(eos0, Tf0 - 4.0) == 6      # pure: ice VI
    assert _ph(eos, Tf0 - 4.0) == 0       # NH3: still liquid
    # anchored depression at 150 ppt / 900 MPa is ~31.8 K, so probe
    # well below it for the eventual freeze
    assert _ph(eos, Tf0 - 40.0) == 6      # NH3: eventually ice VI
