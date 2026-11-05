"""Tests for the angular-mean density correction (PlanetProfile/Thermodynamics/AngularMeanDensity.py).

Synthetic EOS objects only -- no Perple_X table or full run required.
"""
import numpy as np
import pytest

from PlanetProfile.Thermodynamics.AngularMeanDensity import (
    AngularMeanDensityEOS, AngularMeanDensityOffset, WrapSilicateDensityCorrection)


class PolyEOS:
    """rho = 3000 + a*P + b*T + c*T**2; carries a domain key like the PlanetProfile wrappers."""
    key = 'synthetic'
    Pmin, Pmax = 0.0, 5000.0

    def __init__(self, a=0.01, b=-0.1, c=0.0):
        self.a, self.b, self.c = a, b, c

    def fn_rho_kgm3(self, P_MPa, T_K, grid=False):
        P, T = np.asarray(P_MPa, float), np.asarray(T_K, float)
        if grid:
            P, T = P[:, None], T[None, :]
        return 3000 + self.a * P + self.b * T + self.c * T**2

    def fn_Seismic(self, P_MPa, T_K, grid=False):
        return (np.ones_like(np.asarray(P_MPa, float)),) * 4


P = np.array([100.0, 1000.0, 3000.0]); D = np.array([-2.0, -6.0, -1.0])


def test_correction_added_and_clamped_above_moho():
    base = PolyEOS(); e = AngularMeanDensityEOS(base, P, D, P_moho_MPa=500.0)
    Pq = np.array([200.0, 1000.0, 2000.0]); Tq = np.full(3, 1200.0)
    expected = base.fn_rho_kgm3(Pq, Tq) + np.array([0.0, -6.0, -3.5])
    assert np.allclose(e.fn_rho_kgm3(Pq, Tq), expected)


def test_fill_values_outside_levels():
    e = AngularMeanDensityEOS(PolyEOS(), P, D, P_moho_MPa=0.0, fill_low=-9.0, fill_high=-4.0)
    assert np.allclose(e.correction([50.0, 4000.0]), [-9.0, -4.0])
    z = AngularMeanDensityEOS(PolyEOS(), P, D, P_moho_MPa=0.0)
    assert np.allclose(z.correction([50.0, 4000.0]), [0.0, 0.0])


def test_grid_mode_corrects_along_pressure_axis():
    base = PolyEOS(); e = AngularMeanDensityEOS(base, P, D, P_moho_MPa=0.0)
    Pq = np.array([100.0, 3000.0]); Tq = np.array([800.0, 1200.0, 1600.0])
    out = e.fn_rho_kgm3(Pq, Tq, grid=True)
    assert np.allclose(out - base.fn_rho_kgm3(Pq, Tq, grid=True), np.array([[-2.0] * 3, [-1.0] * 3]))


def test_other_properties_and_domain_key_are_delegated():
    base = PolyEOS(); e = AngularMeanDensityEOS(base, P, D)
    assert e.key == 'synthetic' and e.Pmax == 5000.0
    assert e.fn_Seismic(P, P)[0].shape == P.shape


def test_wrap_returns_base_when_correction_is_absent_or_zero():
    base = PolyEOS()
    assert WrapSilicateDensityCorrection(base, None) is base
    assert WrapSilicateDensityCorrection(base, dict(P_MPa=P, drho_kgm3=np.zeros(3))) is base
    assert isinstance(WrapSilicateDensityCorrection(base, dict(P_MPa=P, drho_kgm3=D)), AngularMeanDensityEOS)


def test_invalid_correction_raises():
    with pytest.raises(ValueError):
        AngularMeanDensityEOS(PolyEOS(), P[::-1], D)
    with pytest.raises(ValueError):
        AngularMeanDensityEOS(PolyEOS(), P, D[:2])


def _field(nLev=4, shape=(6, 12), seed=1):
    rng = np.random.default_rng(seed)
    T3 = 1300 + 150 * rng.standard_normal((nLev,) + shape)
    w = np.cos(np.linspace(-1.4, 1.4, shape[0]))[:, None] * np.ones(shape)
    return T3, w


def _structure():
    z = np.linspace(0.0, 1400.0, 50)
    return z, 3.3 * z, np.full(z.size, 1300.0)


def test_offset_is_zero_for_density_linear_in_temperature():
    T3, w = _field(); z, Pn, _ = _structure()
    zl = np.array([100.0, 400.0, 800.0, 1200.0])
    Tbar = np.array([(T3[k] * w).sum() / w.sum() for k in range(4)])
    # radial temperature equal to the weighted mean of each level
    zn = np.concatenate([zl, [0.0, 1400.0]]); Tn = np.concatenate([Tbar, [Tbar[0], Tbar[-1]]])
    Pn = 3.3 * zn
    Pl, d = AngularMeanDensityOffset(PolyEOS(c=0.0), zl, T3, w, zn, Pn, Tn)
    assert np.allclose(Pl, 3.3 * zl) and np.allclose(d, 0.0, atol=1e-9)


def test_offset_equals_weighted_variance_for_quadratic_density():
    T3, w = _field(); c = 2e-4
    zl = np.array([100.0, 400.0, 800.0, 1200.0])
    Tbar = np.array([(T3[k] * w).sum() / w.sum() for k in range(4)])
    var = np.array([((T3[k] - Tbar[k])**2 * w).sum() / w.sum() for k in range(4)])
    zn = np.concatenate([zl, [0.0, 1400.0]]); Tn = np.concatenate([Tbar, [Tbar[0], Tbar[-1]]])
    _, d = AngularMeanDensityOffset(PolyEOS(c=c), zl, T3, w, zn, 3.3 * zn, Tn)
    assert np.allclose(d, c * var)
