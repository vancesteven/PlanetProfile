"""Validation of the CoolProp NH3-H2O ocean EOS
(Thermodynamics/NH3/NH3Props.py + the HydroEOS 'NH3' branch).

Anchors:
- density/Cp vs the Melinder (2010) aqueous-ammonia polynomials
  (CoolProp INCOMP::MAM) at low pressure;
- pure-water-limit density vs SeaFreeze water1 across ocean pressures;
- Leliwa-Kopystynski (2002) liquidus depression values and the shifted
  phase diagram (freezing point drops with w; ices are pure H2O);
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
    NH3Props, NH3liquidusShift_K, _molefrac_NH3)


def test_liquidus_shift_values():
    # dT = 53.8 X + 650 X^3
    assert NH3liquidusShift_K(0.0) == 0.0
    assert NH3liquidusShift_K(100.0) == pytest.approx(6.03, abs=0.01)
    assert NH3liquidusShift_K(229.0) == pytest.approx(
        53.8 * 0.229 + 650 * 0.229 ** 3, rel=1e-12)


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

    # Freezing-point depression (mu mode >= L-K): straddle just below
    # the pure-water melting point — pure water frozen, NH3 ocean
    # still liquid.
    dT = NH3liquidusShift_K(150.0)
    assert dT == pytest.approx(10.26, abs=0.01)  # L-K value (lower bound)
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
    # 5 wt% at 1 bar: ideal colligative limit ~5.7 K; the mu-based
    # liquidus must land near it (the L-K polynomial gives 2.8 K).
    from PlanetProfile.Thermodynamics.NH3.NH3Props import muLiquidusCurve_K
    Tliq, Tm = muLiquidusCurve_K(np.array([0.5]), 50.0)
    dT = float(Tm[0] - Tliq[0])
    assert 4.5 < dT < 6.5, dT


def test_mu_liquidus_high_pressure_ices():
    # Depression must apply under the HP ices too (V at ~620 MPa, VI at
    # ~900 MPa) and grow with pressure at fixed w (stronger non-ideality)
    from PlanetProfile.Thermodynamics.NH3.NH3Props import muLiquidusCurve_K
    P = np.array([100.0, 620.0, 900.0])
    Tliq, Tm = muLiquidusCurve_K(P, 150.0)
    dT = Tm - Tliq
    assert np.all(np.isfinite(dT)), dT
    assert np.all(dT > 8.0) and np.all(dT < 30.0), dT
    assert dT[2] > dT[0], dT  # deeper -> larger depression here


def test_liquidus_solver_gates_machine_b():
    """Machine B defect-report gates (plans/HANDOFF-2026-07-26-nh3-
    liquidus-defect.md): gate 1 — dTliq(P) smooth (no isolated point
    > 0.3 K off its local median), no silent L-K pins, no cold spikes,
    reviewer target table ~1.1/3/5/15 K at w=10/30/50/150."""
    from PlanetProfile.Thermodynamics.NH3.NH3Props import (
        muLiquidusCurve_K)
    P = np.arange(10.0, 220.0, 8.0)
    Tliq, Tm = muLiquidusCurve_K(P, 30.0)
    dep = Tm - Tliq
    assert np.all(np.isfinite(dep))
    med = np.array([np.median(dep[max(0, k-2):k+3])
                    for k in range(dep.size)])
    assert np.max(np.abs(dep - med)) <= 0.3
    assert not np.any(np.abs(dep - NH3liquidusShift_K(30.0)) < 0.02)
    assert dep.max() < 8.0
    for w, target, tol in ((10., 1.1, 0.6), (50., 5.4, 1.2),
                           (150., 15., 3.0)):
        Tl, T0 = muLiquidusCurve_K(np.array([100.0]), w)
        assert abs(float(T0[0] - Tl[0]) - target) < tol


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
    assert _ph(eos, Tf0 - 28.0) == 6      # NH3: eventually ice VI
