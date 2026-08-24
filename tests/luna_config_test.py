"""Configuration guards for the lunar PlanetProfile body definition.

PPLuna.py was previously a verbatim copy of the Io configuration on every
branch. These tests pin the lunar bulk parameters and assert no Io value or
token survives, so that regression cannot recur silently.

Configuration-level only: nothing here runs a structural model, and no
structural output is asserted (a scientific 1-D lunar run is not valid until
the prescribed-temperature mode exists).
"""
import importlib.util
from pathlib import Path

import pytest

PPLUNA = Path(__file__).resolve().parents[1] / 'PlanetProfile' / 'Default' / 'Luna' / 'PPLuna.py'


@pytest.fixture(scope='module')
def planet():
    spec = importlib.util.spec_from_file_location('PPLuna_under_test', PPLUNA)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.Planet


@pytest.fixture(scope='module')
def source():
    return PPLUNA.read_text()


# --- body identity and bulk parameters ---------------------------------

def test_body_is_luna(planet):
    assert planet.name == 'Luna'


def test_radius_is_lunar(planet):
    # Neumann (2013) LOLA mean radius, as adopted in Williams et al. (2014) Table 2.
    assert planet.Bulk.R_m == pytest.approx(1737.15e3)


def test_mass_is_lunar(planet):
    # GM/G with GM from Williams, Boggs & Folkner (2013) DE430 Table 5.
    assert planet.Bulk.M_kg == pytest.approx(7.3458e22, rel=1e-4)


def test_moi_is_whole_moon_value(planet):
    # Konopliv et al. (1998) whole-Moon average moment. Must NOT be the
    # solid-Moon value (0.392728), which is not comparable to PlanetProfile's
    # whole-body shell sum.
    assert planet.Bulk.Cmeasured == pytest.approx(0.3931)
    assert planet.Bulk.Cuncertainty == pytest.approx(0.0002)
    assert planet.Bulk.Cmeasured != pytest.approx(0.392728, abs=1e-6)


def test_dry_body_with_iron_core(planet):
    assert planet.Do.NO_H2O is True
    assert planet.Do.Fe_CORE is True


# --- crust --------------------------------------------------------------

def test_crust_enabled_and_structure_only(planet):
    assert planet.Crust.enabled is True
    assert planet.Crust.structure_only is True
    assert planet.Crust.bulk_density_includes_porosity is True


def test_crust_thickness_is_nominal_40km(planet):
    assert planet.Crust.thickness_m == pytest.approx(40.0e3)


def test_crust_density_reproduces_garcia_m3_shell_mass(planet):
    """The single bulk density must reproduce the mass of the layered M3
    profile it summarizes, or the structural mass budget is wrong."""
    import numpy as np
    R = planet.Bulk.R_m
    rM = R - planet.Crust.thickness_m
    z = np.linspace(0.0, 40.0, 4001)
    rho = np.interp(z, [0.0, 20.0, 40.0], [2600.0, 2760.0, 2890.0])  # Garcia 2019 M3
    r = R - z * 1e3
    mExact = abs(np.trapz(rho * 4 * np.pi * r ** 2, z * 1e3))
    mUsed = planet.Crust.bulk_density_kgm3 * 4 / 3 * np.pi * (R ** 3 - rM ** 3)
    # Tolerance reflects rounding the derived mean to 2751.4 kg/m^3 (1 decimal).
    # Carrying more digits would assert precision Garcia M3's own inputs
    # (2600 / 2760 / 2890 kg/m^3) do not support. 1e-4 relative is ~4e17 kg on
    # a 4.08e21 kg crust -- far below any structurally meaningful difference.
    assert mUsed == pytest.approx(mExact, rel=1e-4)


def test_crust_elastic_properties_not_yet_connected(planet):
    """Expected at this stage: the crust is structure-only. This test should
    be updated, not deleted, when Garcia Vp/Vs/Ks/G are wired in."""
    assert planet.Crust.seismic_status == 'CRUST_SEISMIC_PROPERTIES_UNAVAILABLE'


# --- EOS wiring ---------------------------------------------------------

def test_mantle_eos_is_warren_unified(planet):
    assert planet.Sil.mantleEOS == 'moon_w05_preferred_stx24_unified.tab'


def test_core_eos_is_traceable_product(planet):
    assert planet.Core.coreEOS == 'Fe-S_3D_EOS_5slice_traceable_v1.mat'
    assert planet.Core.wFe_ppt == 950


def test_porous_rock_disabled(planet):
    # Required by STRUCTURE_ONLY_CRUST_MODE and appropriate for a solid-matrix
    # mantle EOS plus an effective-density crust.
    assert planet.Do.POROUS_ROCK is False


# --- no Io residue ------------------------------------------------------

@pytest.mark.parametrize('token', [
    "PlanetStruct('Io')",   # body constructor
    'PPIo',                 # docstring header
    'CV3hy1wt',             # hydrated CV3 chondrite mantle EOS
    '1845.9e-6',            # Io J2
    '553.7e-6',             # Io C22
    '1821.49e3',            # Io radius
    '8.932e22',             # Io mass
    '0.37685',              # Io C/MR^2
])
def test_no_io_token_survives(source, token):
    assert token not in source, f'Io configuration residue found: {token}'


@pytest.mark.parametrize('attr,bad', [
    ('Bulk.R_m', 1821.49e3),
    ('Bulk.M_kg', 8.932e22),
    ('Bulk.Tsurf_K', 110),
    ('Bulk.Cmeasured', 0.37685),
    ('Bulk.qSurf_Wm2', 140e-3),
])
def test_no_io_value_survives(planet, attr, bad):
    obj = planet
    for part in attr.split('.'):
        obj = getattr(obj, part)
    assert obj != bad, f'Io value survives in {attr}'
