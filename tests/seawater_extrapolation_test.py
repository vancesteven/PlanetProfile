"""Guard the documented high-salinity GSW Seawater extrapolation.

Europa v3 Gate 2 bounded the PSS-78 extrapolation systematic at 18.5%
above 42 ppt. Enceladus retains 70 and 100 ppt nodes to display the
induction-excluded region, so rho, Cp, and sigma must remain finite and
monotone in salinity at its representative 3 MPa, 271-274 K conditions.
"""

import numpy as np

from PlanetProfile.Thermodynamics.HydroEOS import GetOceanEOS


SALINITIES_PPT = (42.0, 70.0, 100.0)
TEMPERATURES_K = np.arange(271.0, 275.0)
PRESSURES_MPA = np.full(TEMPERATURES_K.shape, 3.0)


def _properties_at(salinity_ppt):
    """Evaluate within a spline grid that brackets the target points."""
    eos = GetOceanEOS(
        "Seawater",
        salinity_ppt,
        np.arange(2.0, 6.0),
        np.arange(270.0, 276.0),
        None,
        FORCE_NEW=True,
    )
    return (
        np.asarray(eos.fn_rho_kgm3(PRESSURES_MPA, TEMPERATURES_K)),
        np.asarray(eos.fn_Cp_JkgK(PRESSURES_MPA, TEMPERATURES_K)),
        np.asarray(eos.fn_sigma_Sm(PRESSURES_MPA, TEMPERATURES_K)),
    )


def test_high_salinity_seawater_properties_are_finite_and_monotone():
    """rho/sigma increase and Cp decreases from 42 through 100 ppt."""
    rho, heat_capacity, conductivity = (
        np.stack(values)
        for values in zip(*(_properties_at(w) for w in SALINITIES_PPT))
    )

    for values in (rho, heat_capacity, conductivity):
        assert np.isfinite(values).all()

    assert np.all(np.diff(rho, axis=0) > 0.0)
    assert np.all(np.diff(heat_capacity, axis=0) < 0.0)
    assert np.all(np.diff(conductivity, axis=0) > 0.0)
