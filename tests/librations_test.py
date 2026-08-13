"""Regression tests for the Enceladus libration forward model.

The published-value case uses the yellow-marker interior from Van Hoolst,
Baland & Trinh (2016), Figs. 2-3 (doi:10.1016/j.icarus.2016.05.025): a
20 km shell, 42.1 km ocean, 190 km core, rho_shell=925 kg/m3, and
rho_ocean=1050 kg/m3. Their observed two-sigma equatorial displacement band
is 466-590 m.

The elastic checks use the committed Tb=272.46 K Enceladus cache node. They
encode the 2026-08-12 scientific adjudication recorded in
plans/CODEX-QUEUE.md: the current rigid branch has a known K_int
normalization defect, while the elastic y1 convention is validated by the
partial-Love-number sum.
"""
import pickle
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from PlanetProfile.Gravity import pyALMA3Updated  # noqa: E402
from PlanetProfile.Gravity.Gravity import (  # noqa: E402
    ALMA_TIME_UNIT_S,
    LibrationModelInputs,
)
from PlanetProfile.Gravity.Librations import librations  # noqa: E402
from PlanetProfile.Utilities.defineStructs import Constants  # noqa: E402

R_ENCELADUS_M = 252.1e3
M_VH16_KG = 1.0794e20
OMEGA_VH16_RADPS = 2.0 * np.pi / (1.37 * 86400.0)
ECCENTRICITY = 0.0047

CACHE_PATH = (
    REPO
    / "PlanetProfile/Test/mcmc_results/Enceladus/Cassini_smoke"
    / "enceladus_seawater_smoke_grid.pkl"
)


def _vh16_interior(shell_thickness_km):
    """Return the mass-conserving yellow-marker radial/density model."""
    core_radius_m = 190.0e3
    shell_base_m = R_ENCELADUS_M - shell_thickness_km * 1e3
    rho_shell = 925.0
    rho_ocean = 1050.0
    volume_factor = 4.0 * np.pi / 3.0
    rho_core = (
        M_VH16_KG
        - rho_ocean
        * volume_factor
        * (shell_base_m**3 - core_radius_m**3)
        - rho_shell
        * volume_factor
        * (R_ENCELADUS_M**3 - shell_base_m**3)
    ) / (volume_factor * core_radius_m**3)
    return (
        np.array([core_radius_m, shell_base_m, R_ENCELADUS_M]),
        np.array([rho_core, rho_ocean, rho_shell]),
    )


def _libration_deg(radial_m, rho_kgm3, *, rigid, y1=None):
    displacement_m = librations(
        radial_m,
        rho_kgm3,
        OMEGA_VH16_RADPS,
        ECCENTRICITY,
        rigid=rigid,
        ocean=True,
        ocean_idx=1,
        y1=y1,
    )
    return float(np.degrees(displacement_m / radial_m[-1]))


@lru_cache(maxsize=1)
def _cache_node_diagnostics():
    """Evaluate the adjudicated cache node once for all elastic checks."""
    with CACHE_PATH.open("rb") as stream:
        cache = pickle.load(stream)

    index = int(np.argmin(np.abs(np.asarray(cache["Tb_K_grid"]) - 272.46)))
    assert cache["Tb_K_grid"][index] == pytest.approx(272.46)
    structure = cache["structures"][index]

    radius_m = np.asarray(structure["r_m"], dtype=float)
    density_kgm3 = np.asarray(structure["rho"], dtype=float)
    shear_modulus_pa = np.asarray(structure["mu_Pa"], dtype=float)
    viscosity_pas = np.asarray(structure["eta_Pa_base"], dtype=float)
    phases = np.asarray(structure["phases"])

    # This cache node contains an elastic silicate interior and ice shell
    # separated by a Newtonian liquid ocean, matching the default PP model.
    rheology = np.where(phases == 0, "newton", "elastic").tolist()
    alma_model = pyALMA3Updated.build_model(
        radius_m,
        density_kgm3,
        shear_modulus_pa,
        viscosity_pas,
        rheology,
        np.zeros((radius_m.size, 2)),
        norm=False,
    )

    omega_radps = float(structure["omega"])
    period_kyr = (2.0 * np.pi / omega_radps) / ALMA_TIME_UNIT_S
    _, _, love_k, radial_functions = pyALMA3Updated.love_numbers(
        alma_model,
        [2],
        [period_kyr],
        "tidal",
        "periodic",
        4,
        verbose=False,
        internal=True,
    )
    alma_k2 = np.complex128(love_k[0, 0, -1])

    body_mass_kg = float(structure["Mtot_kg"])
    body_radius_m = float(structure["R_body_m"])
    surface_gravity_ms2 = Constants.G * body_mass_kg / body_radius_m**2
    y1 = (
        np.real(radial_functions[:, 0])
        * body_mass_kg
        / body_radius_m
        / surface_gravity_ms2
    )
    reduced_r, reduced_rho, ocean_idx, reduced_y1 = LibrationModelInputs(
        radius_m, phases, density_kgm3, y1
    )

    kwargs = dict(
        omega=omega_radps,
        ecc=float(structure["eccentricity"]),
        ocean=True,
        ocean_idx=ocean_idx,
    )
    rigid_m = librations(reduced_r, reduced_rho, rigid=True, **kwargs)
    elastic_m = librations(
        reduced_r, reduced_rho, rigid=False, y1=reduced_y1, **kwargs
    )

    # Van Hoolst et al. (2013), Eq. 24, as implemented in the elastic
    # libration branch. This independent sum guards the y1 convention.
    r_core, r_ocean_top, r_body = reduced_r
    rho_core, rho_ocean, rho_shell = reduced_rho
    y_core, y_ocean_top, y_surface = reduced_y1
    factor = 4.0 / 5.0 * np.pi * Constants.G / r_body**3
    partial_k2 = (
        factor * rho_shell * (r_body**4 * y_surface - r_ocean_top**4 * y_ocean_top)
        + factor * rho_ocean * r_ocean_top**4 * y_ocean_top
        - factor * rho_ocean * r_core**4 * y_core
        + factor * rho_core * r_core**4 * y_core
    )

    return {
        "rigid_deg": float(np.degrees(rigid_m / body_radius_m)),
        "elastic_deg": float(np.degrees(elastic_m / body_radius_m)),
        "alma_k2": alma_k2,
        "partial_k2": float(partial_k2),
    }


def test_vh16_yellow_marker_lands_in_published_band():
    """The VH16 yellow-marker interior reproduces the 466-590 m band.

    The current result is 524.38 m. This band intentionally does not
    discriminate the known rigid K_int normalization: the adjudicated
    corrected result, 532.70 m, also lies inside the published band.
    """
    radial_m, rho_kgm3 = _vh16_interior(shell_thickness_km=20.0)
    libration_m = np.radians(
        _libration_deg(radial_m, rho_kgm3, rigid=True)
    ) * radial_m[-1]
    assert 466.0 <= libration_m <= 590.0


def test_libration_decreases_with_shell_thickness():
    """Rigid libration decreases over the fixed-core 5-45 km sweep."""
    amplitudes_deg = []
    for shell_thickness_km in np.arange(5.0, 46.0, 5.0):
        radial_m, rho_kgm3 = _vh16_interior(shell_thickness_km)
        amplitudes_deg.append(_libration_deg(radial_m, rho_kgm3, rigid=True))
    assert np.all(np.diff(amplitudes_deg) < 0.0)


def test_adjudicated_rigid_elastic_discrepancy_is_pinned():
    """Pin the known defect until the reviewed K_int repair lands.

    The 2026-08-12 adjudication found the rigid branch's K_int is missing
    8*pi/15. This pin must be updated in the same commit as any K_int fix.
    """
    values = _cache_node_diagnostics()
    assert values["rigid_deg"] == pytest.approx(0.110905164, rel=1e-6)
    assert values["elastic_deg"] == pytest.approx(0.111814157, rel=1e-6)
    assert values["elastic_deg"] / values["rigid_deg"] == pytest.approx(
        1.008196, rel=1e-4
    )


@pytest.mark.xfail(
    strict=True,
    reason="rigid K_int normalization bug; update the discrepancy pin with its fix",
)
def test_elasticity_has_the_vh16_physical_direction():
    """Elasticity should reduce libration once the rigid branch is repaired."""
    values = _cache_node_diagnostics()
    assert values["elastic_deg"] < values["rigid_deg"]


def test_partial_love_numbers_reproduce_alma_k2():
    """The y1-derived partial Love numbers sum to ALMA Re(k2) within 0.1%."""
    values = _cache_node_diagnostics()
    relative_error = abs(
        values["partial_k2"] - values["alma_k2"].real
    ) / abs(values["alma_k2"].real)
    assert relative_error < 1e-3
