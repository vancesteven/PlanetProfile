#!/usr/bin/env python3
"""Regression tests for isolated failed-pixel repair in 3D lateral models."""
import importlib.util
import numpy as np
from pathlib import Path


ROOT = Path(__file__).parent


def load_lateral_structure_helpers():
    spec = importlib.util.spec_from_file_location(
        'lateral_structure_direct',
        ROOT / 'PlanetProfile' / 'Lateral' / 'LateralStructure.py'
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.RepairLateralField


def load_compute_tidal_heating():
    if importlib.util.find_spec('scipy') is None:
        return None, 'scipy is not installed'

    try:
        from PlanetProfile.Lateral.TidalHeating3D import ComputeTidalHeating3D
        return ComputeTidalHeating3D, None
    except Exception as exc:
        return None, exc


RepairLateralField = load_lateral_structure_helpers()
ComputeTidalHeating3D, TIDAL_IMPORT_ERROR = load_compute_tidal_heating()


class Struct:
    pass


def make_struct(**kwargs):
    obj = Struct()
    for key, value in kwargs.items():
        setattr(obj, key, value)
    return obj


def make_lateral(n_pix=4):
    theta_rad = np.deg2rad([20.0, 35.0, 50.0, 65.0])[:n_pix]
    phi_rad = np.deg2rad([0.0, 90.0, 180.0, 270.0])[:n_pix]
    # Equal-area pixels for test (4π steradians / nPix)
    pixel_area = 4.0 * np.pi / n_pix
    # Uniform ice thickness for test (25 km)
    dIce_m = np.full(n_pix, 25.0e3)
    return make_struct(
        nPix=n_pix,
        theta_rad=theta_rad,
        phi_rad=phi_rad,
        pixArea_sr=pixel_area,
        dIce_m=dIce_m,
        repairedColumnMask=None,
        failedColumnMask=None,
        repairLog=None,
        maxRepairFrac=0.5,
    )


def test_repair_lateral_field():
    lateral = make_lateral()
    field = np.array([np.nan, 10.0, 20.0, 30.0])
    repaired = RepairLateralField(
        lateral, field, 'testField', invalidMask=np.array([True, False, False, False])
    )

    assert np.all(np.isfinite(repaired))
    assert repaired[0] > 0.0
    assert lateral.repairedColumnMask[0]
    assert lateral.repairLog['testField']['indices'] == [0]


def test_repair_lateral_field_rejects_clustered_failures():
    lateral = make_lateral()
    lateral.maxRepairFrac = 0.25
    field = np.array([np.nan, np.nan, 20.0, 30.0])

    try:
        RepairLateralField(
            lateral, field, 'clusteredField',
            invalidMask=np.array([True, True, False, False])
        )
    except ValueError as exc:
        assert 'exceeding maxRepairFrac' in str(exc)
    else:
        raise AssertionError('clustered failures should not be automatically repaired')


def make_column(valid=True):
    if not valid:
        return make_struct(invalidReason='ZeroDivisionError: division by zero')

    return make_struct(
        invalidReason='Valid',
        Steps=make_struct(nSurfIce=3, nHydro=3),
        Seismic=make_struct(GS_GPa=np.array([3.6, 3.5, 3.4])),
        eta_Pas=np.array([1.0e15, 3.0e14, 1.0e14]),
        T_K=np.array([110.0, 180.0, 260.0]),
        phase=np.array([1, 1, 1]),
    )


def test_compute_tidal_heating_repairs_failed_column():
    if ComputeTidalHeating3D is None:
        print(f'skipping tidal heating integration test: {TIDAL_IMPORT_ERROR}')
        return

    lateral = make_lateral()
    planet = make_struct(
        Lateral=lateral,
        Bulk=make_struct(
            eccentricity=0.0094,
            meanMotion_radps=2.05e-5,
            M_kg=4.8e22,
            R_m=1560.8e3,
        ),
        Magnetic=make_struct(Texc_hr=None),
    )
    columns = np.array([
        make_column(valid=False),
        make_column(valid=True),
        make_column(valid=True),
        make_column(valid=True),
    ], dtype=object)

    planet = ComputeTidalHeating3D(planet, make_struct(), columns, rheology='maxwell')

    assert np.all(np.isfinite(planet.Lateral.HtidalIce_Wm3))
    assert np.all(planet.Lateral.HtidalIce_Wm3 > 0.0)
    assert planet.Lateral.failedColumnMask[0]
    assert planet.Lateral.repairedColumnMask[0]
    assert planet.Lateral.repairLog['HtidalIce_Wm3']['indices'] == [0]


def main():
    test_repair_lateral_field()
    test_repair_lateral_field_rejects_clustered_failures()
    test_compute_tidal_heating_repairs_failed_column()
    print('lateral failure repair tests passed')


if __name__ == '__main__':
    main()
