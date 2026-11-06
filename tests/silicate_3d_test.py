"""Tests for EvaluateSilicate3D (PlanetProfile/Thermodynamics/Silicate3D.py).

Synthetic EOS only -- no Perple_X table or full run required.
"""
import numpy as np
import pytest

from PlanetProfile.Thermodynamics.Silicate3D import EvaluateSilicate3D, LevelMean


class LinEOS:
    Pmin, Pmax, Tmin, Tmax = 0.0, 5000.0, 300.0, 2200.0

    def fn_rho_kgm3(self, P, T, grid=False):
        return 3300 + 0.02 * np.asarray(P) - 0.1 * np.asarray(T)

    def fn_Seismic(self, P, T, grid=False):
        P, T = np.asarray(P, float), np.asarray(T, float)
        return 8 + 1e-4 * P, 4.5 - 1e-4 * T, 100 + 0.01 * P, 60 - 0.01 * T


def _case(shape=(3, 4, 8)):
    rng = np.random.default_rng(0)
    T3 = 1300 + 100 * rng.standard_normal(shape)
    w = np.ones(shape[1:]); w[0] = 2.0
    zl = np.array([100.0, 500.0, 900.0])
    zn = np.linspace(0.0, 1400.0, 30)
    return T3, w, zl, zn, 3.3 * zn


def test_shapes_and_level_pressures():
    T3, w, zl, zn, Pn = _case()
    o = EvaluateSilicate3D(LinEOS(), zl, T3, w, zn, Pn)
    assert np.allclose(o['P_MPa'], 3.3 * zl)
    for k in ('rho_kgm3', 'VP_kms', 'VS_kms', 'KS_GPa', 'GS_GPa'):
        assert o[k].shape == T3.shape and o[k + '_ref'].shape == (3,) and o[k + '_mean'].shape == (3,)


def test_values_match_direct_evaluation():
    T3, w, zl, zn, Pn = _case(); e = LinEOS()
    o = EvaluateSilicate3D(e, zl, T3, w, zn, Pn)
    P3 = np.broadcast_to((3.3 * zl)[:, None, None], T3.shape)
    assert np.allclose(o['VS_kms'], e.fn_Seismic(P3, T3)[1])
    assert np.allclose(o['rho_kgm3'], e.fn_rho_kgm3(P3, T3))


def test_reference_uses_weighted_mean_temperature():
    T3, w, zl, zn, Pn = _case(); e = LinEOS()
    o = EvaluateSilicate3D(e, zl, T3, w, zn, Pn)
    assert np.allclose(o['Tbar_K'], LevelMean(T3, w))
    # linear EOS: the level mean of the 3-D field equals the reference at the mean temperature
    assert np.allclose(o['rho_kgm3_mean'], o['rho_kgm3_ref'])


def test_given_reference_temperature_is_used():
    T3, w, zl, zn, Pn = _case(); e = LinEOS()
    o = EvaluateSilicate3D(e, zl, T3, w, zn, Pn, Tbar_K=np.array([1000.0, 1100.0, 1200.0]))
    assert np.allclose(o['GS_GPa_ref'], 60 - 0.01 * np.array([1000.0, 1100.0, 1200.0]))


def test_points_outside_domain_are_counted():
    T3, w, zl, zn, Pn = _case(); T3[0, 0, 0] = 2500.0
    o = EvaluateSilicate3D(LinEOS(), zl, T3, w, zn, Pn)
    assert o['nOutsideDomain'] == 1
