"""Tests for the seismic properties of an explicit structural crust (Seismic.AssignCrustSeismic).

Default: structure-only crust, seismic properties NaN. With Crust.VP_kms and Crust.VS_kms set:
prescribed velocities, GS = rho VS^2 and KS = rho (VP^2 - 4/3 VS^2) at the crust density.
Uses a minimal stand-in for Planet -- no table or full run required.
"""
from types import SimpleNamespace

import numpy as np
import pytest

from PlanetProfile.Thermodynamics.Seismic import AssignCrustSeismic
from PlanetProfile.Utilities.defineStructs import Constants


def _planet(VP=None, VS=None, QS=None):
    mat = np.array([Constants.materialCrust] * 3 + [Constants.materialMantle] * 2)
    n = mat.size
    return SimpleNamespace(
        material=mat, rho_kgm3=np.array([2600.0, 2700.0, 2750.0, 3300.0, 3350.0]),
        Cp_JkgK=np.ones(n), alpha_pK=np.ones(n),
        Seismic=SimpleNamespace(VP_kms=np.zeros(n), VS_kms=np.zeros(n), QS=np.zeros(n),
                                KS_GPa=np.zeros(n), GS_GPa=np.zeros(n)),
        Crust=SimpleNamespace(VP_kms=VP, VS_kms=VS, QS=QS,
                              seismic_status='CRUST_SEISMIC_PROPERTIES_UNAVAILABLE'))


def test_default_crust_is_structure_only():
    P = _planet(); AssignCrustSeismic(P)
    for a in (P.Seismic.VP_kms, P.Seismic.VS_kms, P.Seismic.QS, P.Seismic.KS_GPa, P.Seismic.GS_GPa):
        assert np.all(np.isnan(a[:3])) and np.all(a[3:] == 0.0)
    assert np.all(np.isnan(P.Cp_JkgK[:3])) and np.all(P.Cp_JkgK[3:] == 1.0)
    assert P.Crust.seismic_status == 'CRUST_SEISMIC_PROPERTIES_UNAVAILABLE'


def test_prescribed_velocities_give_moduli_at_crust_density():
    P = _planet(VP=6.5, VS=3.7, QS=500.0); AssignCrustSeismic(P)
    rho = P.rho_kgm3[:3]
    assert np.allclose(P.Seismic.GS_GPa[:3], rho * 3700.0**2 * 1e-9)
    assert np.allclose(P.Seismic.KS_GPa[:3], rho * (6500.0**2 - 4 / 3 * 3700.0**2) * 1e-9)
    assert np.all(P.Seismic.QS[:3] == 500.0) and np.all(P.Seismic.VP_kms[3:] == 0.0)
    assert np.all(np.isnan(P.Cp_JkgK[:3]))
    assert P.Crust.seismic_status == 'CRUST_SEISMIC_PROPERTIES_PRESCRIBED'


def test_only_one_velocity_leaves_nan():
    P = _planet(VP=6.5); AssignCrustSeismic(P)
    assert np.all(np.isnan(P.Seismic.GS_GPa[:3]))


def test_unphysical_velocities_raise():
    with pytest.raises(ValueError):
        AssignCrustSeismic(_planet(VP=4.0, VS=3.7))
