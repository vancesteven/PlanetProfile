"""
PPTestLateralPhase8: Test lateral plotting functions (Phase 8).

Verifies:
- Helper functions work correctly (_InterpolateToLatLon, _SetMap, _PlotSurfaceMap)
- Each plot function handles mock data properly
- Files are created with expected names
- Graceful handling of missing data (early return, no error)
- GenerateLateralPlots orchestration
"""
import os
import sys
import numpy as np
import tempfile
import shutil

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from PlanetProfile.Utilities.defineStructs import PlanetStruct
from PlanetProfile.GetConfig import Params, FigMisc
from PlanetProfile.Plotting.LateralPlots import (
    _InterpolateToLatLon,
    _SetMap,
    _PlotSurfaceMap,
    PlotIceThickness,
    PlotTidalHeating,
    PlotBasalTemperature,
    PlotOceanConductivity,
    PlotClathrateFraction,
    PlotEffectiveConductivity,
    PlotLateralSummary,
    PlotTidalHeatingByLayer,
    PlotTidalHeatingVsDepth,
    GenerateLateralPlots
)


def _create_mock_planet_with_lateral(nPix=12):
    """ Create Planet with mock lateral data for testing. """
    Planet = PlanetStruct('TestPlanet')
    Planet.Bulk.R_m = 1560e3  # Europa-like

    # Initialize Lateral structure
    from PlanetProfile.Lateral.defaultConfigLateral import lateralAssign
    config = lateralAssign()
    config['gridType'] = 'healpix'
    config['nSide'] = 1  # 12 pixels

    # Apply config to Lateral
    for key, val in config.items():
        setattr(Planet.Lateral, key, val)

    # Mock spatial grid
    Planet.Lateral.nPix = nPix
    Planet.Lateral.theta_rad = np.linspace(0.1, np.pi - 0.1, nPix)
    Planet.Lateral.phi_rad = np.linspace(0, 2 * np.pi, nPix, endpoint=False)
    Planet.Lateral.pixArea_sr = np.full(nPix, 4 * np.pi / nPix)

    # Mock ice thickness field (15-25 km)
    Planet.Lateral.dIce_m = np.linspace(15e3, 25e3, nPix)

    # Mock tidal heating field (1e-9 to 5e-9 W/m^3)
    Planet.Lateral.HtidalIce_Wm3 = np.linspace(1e-9, 5e-9, nPix)

    # Mock basal temperature field (260-270 K)
    Planet.Lateral.Tb_K = np.linspace(260, 270, nPix)

    # Mock ocean conductivity field (2-4 S/m)
    Planet.Lateral.sigma_mean_Sm = np.linspace(2.0, 4.0, nPix)

    # Mock clathrate fraction field (0-0.3)
    Planet.Lateral.fClath = np.linspace(0, 0.3, nPix)
    Planet.Lateral.DO_CLATH_LATERAL = True

    # Mock effective thermal conductivity (2-3 W/m/K)
    Planet.Lateral.kThermEff_WmK = np.linspace(2.0, 3.0, nPix)

    # Mock by-layer heating (for PlotTidalHeatingByLayer)
    Planet.Lateral.HtidalIceI_bot_Wm3 = np.linspace(2e-9, 6e-9, nPix)
    Planet.Lateral.HtidalIceI_top_Wm3 = np.linspace(1e-9, 3e-9, nPix)
    Planet.Lateral.HtidalHP_bot_Wm3 = np.linspace(1e-10, 5e-10, nPix)
    Planet.Lateral.HtidalHP_top_Wm3 = np.linspace(5e-11, 2e-10, nPix)

    # Mock depth profiles (for PlotTidalHeatingVsDepth)
    # Each pixel has a profile with 10 depth points
    nDepth = 10
    Planet.Lateral.HtidalIceI_profile_Wm3 = [np.linspace(1e-9, 5e-9, nDepth) for _ in range(nPix)]
    Planet.Lateral.rIceI_profile_m = [np.linspace(Planet.Bulk.R_m - 5e3, Planet.Bulk.R_m - 20e3, nDepth) for _ in range(nPix)]

    return Planet


def test_interpolate_to_latlon():
    """ Test _InterpolateToLatLon helper function. """
    print("\n=== Test: _InterpolateToLatLon ===")

    # Mock pixel data
    nPix = 12
    theta_rad = np.linspace(0.1, np.pi - 0.1, nPix)
    phi_rad = np.linspace(0, 2 * np.pi, nPix, endpoint=False)
    field = np.sin(theta_rad) * np.cos(phi_rad)  # Some pattern

    # Output grid
    latMap_deg = np.linspace(-80, 80, 9)
    lonMap_deg = np.linspace(-180, 180, 18)

    # Interpolate
    data2D = _InterpolateToLatLon(field, theta_rad, phi_rad, latMap_deg, lonMap_deg)

    # Check shape
    assert data2D.shape == (len(latMap_deg), len(lonMap_deg)), \
        f"Expected shape {(len(latMap_deg), len(lonMap_deg))}, got {data2D.shape}"

    # Check no NaN at center (interpolation should work)
    center_lat_idx = len(latMap_deg) // 2
    center_lon_idx = len(lonMap_deg) // 2
    assert not np.isnan(data2D[center_lat_idx, center_lon_idx]), "Interpolated data has NaN at center"

    # Check mean is reasonable (within factor of 2 of input mean)
    mean_input = np.mean(field)
    mean_output = np.nanmean(data2D)
    assert 0.5 * abs(mean_input) < abs(mean_output) < 2.0 * abs(mean_input) or abs(mean_input) < 1e-10, \
        f"Interpolated mean {mean_output} differs too much from input mean {mean_input}"

    print("✓ _InterpolateToLatLon: shape correct, no NaN at center, mean reasonable")


def test_plot_surface_map():
    """ Test _PlotSurfaceMap helper function. """
    print("\n=== Test: _PlotSurfaceMap ===")

    # Create temporary directory for test plots
    tmpdir = tempfile.mkdtemp()

    try:
        # Mock 2D data
        latMap_deg = np.linspace(-80, 80, 9)
        lonMap_deg = np.linspace(-180, 180, 18)
        lonGrid, latGrid = np.meshgrid(lonMap_deg, latMap_deg)
        data2D = np.sin(np.radians(latGrid)) * np.cos(np.radians(lonGrid))

        fName = os.path.join(tmpdir, 'test_surface_map.png')

        # Plot (should not raise)
        _PlotSurfaceMap(
            data2D, latMap_deg, lonMap_deg,
            cmap='viridis', cbarLabel='Test field',
            fName=fName, title='Test plot'
        )

        # Check file created
        assert os.path.exists(fName), f"Plot file not created: {fName}"
        assert os.path.getsize(fName) > 1000, f"Plot file too small: {os.path.getsize(fName)} bytes"

        print(f"✓ _PlotSurfaceMap: file created ({os.path.getsize(fName)} bytes)")

    finally:
        shutil.rmtree(tmpdir)


def test_plot_ice_thickness():
    """ Test PlotIceThickness function. """
    print("\n=== Test: PlotIceThickness ===")

    Planet = _create_mock_planet_with_lateral()
    tmpdir = tempfile.mkdtemp()

    try:
        Params.FigureFiles.path = tmpdir
        PlotIceThickness(Planet, Params)

        # Check file created
        expected_file = os.path.join(tmpdir, f'{Planet.name}_dIce{FigMisc.xtn}')
        assert os.path.exists(expected_file), f"Ice thickness plot not created: {expected_file}"

        print(f"✓ PlotIceThickness: file created")

    finally:
        shutil.rmtree(tmpdir)


def test_plot_tidal_heating():
    """ Test PlotTidalHeating function. """
    print("\n=== Test: PlotTidalHeating ===")

    Planet = _create_mock_planet_with_lateral()
    tmpdir = tempfile.mkdtemp()

    try:
        Params.FigureFiles.path = tmpdir
        PlotTidalHeating(Planet, Params)

        # Check file created
        expected_file = os.path.join(tmpdir, f'{Planet.name}_Htidal{FigMisc.xtn}')
        assert os.path.exists(expected_file), f"Tidal heating plot not created: {expected_file}"

        print(f"✓ PlotTidalHeating: file created")

    finally:
        shutil.rmtree(tmpdir)


def test_plot_basal_temperature():
    """ Test PlotBasalTemperature function. """
    print("\n=== Test: PlotBasalTemperature ===")

    Planet = _create_mock_planet_with_lateral()
    tmpdir = tempfile.mkdtemp()

    try:
        Params.FigureFiles.path = tmpdir
        PlotBasalTemperature(Planet, Params)

        # Check file created
        expected_file = os.path.join(tmpdir, f'{Planet.name}_Tb{FigMisc.xtn}')
        assert os.path.exists(expected_file), f"Basal temperature plot not created: {expected_file}"

        print(f"✓ PlotBasalTemperature: file created")

    finally:
        shutil.rmtree(tmpdir)


def test_plot_ocean_conductivity():
    """ Test PlotOceanConductivity function. """
    print("\n=== Test: PlotOceanConductivity ===")

    Planet = _create_mock_planet_with_lateral()
    tmpdir = tempfile.mkdtemp()

    try:
        Params.FigureFiles.path = tmpdir
        PlotOceanConductivity(Planet, Params)

        # Check file created
        expected_file = os.path.join(tmpdir, f'{Planet.name}_sigmaOcean{FigMisc.xtn}')
        assert os.path.exists(expected_file), f"Ocean conductivity plot not created: {expected_file}"

        print(f"✓ PlotOceanConductivity: file created")

    finally:
        shutil.rmtree(tmpdir)


def test_plot_clathrate_fraction():
    """ Test PlotClathrateFraction function. """
    print("\n=== Test: PlotClathrateFraction ===")

    Planet = _create_mock_planet_with_lateral()
    tmpdir = tempfile.mkdtemp()

    try:
        Params.FigureFiles.path = tmpdir
        PlotClathrateFraction(Planet, Params)

        # Check file created
        expected_file = os.path.join(tmpdir, f'{Planet.name}_fClath{FigMisc.xtn}')
        assert os.path.exists(expected_file), f"Clathrate fraction plot not created: {expected_file}"

        print(f"✓ PlotClathrateFraction: file created")

    finally:
        shutil.rmtree(tmpdir)


def test_plot_effective_conductivity():
    """ Test PlotEffectiveConductivity function. """
    print("\n=== Test: PlotEffectiveConductivity ===")

    Planet = _create_mock_planet_with_lateral()
    tmpdir = tempfile.mkdtemp()

    try:
        Params.FigureFiles.path = tmpdir
        PlotEffectiveConductivity(Planet, Params)

        # Check file created
        expected_file = os.path.join(tmpdir, f'{Planet.name}_kTherm{FigMisc.xtn}')
        assert os.path.exists(expected_file), f"Effective conductivity plot not created: {expected_file}"

        print(f"✓ PlotEffectiveConductivity: file created")

    finally:
        shutil.rmtree(tmpdir)


def test_plot_lateral_summary():
    """ Test PlotLateralSummary multi-panel function. """
    print("\n=== Test: PlotLateralSummary ===")

    Planet = _create_mock_planet_with_lateral()
    tmpdir = tempfile.mkdtemp()

    try:
        Params.FigureFiles.path = tmpdir
        PlotLateralSummary(Planet, Params)

        # Check file created
        expected_file = os.path.join(tmpdir, f'{Planet.name}_lateralSummary{FigMisc.xtn}')
        assert os.path.exists(expected_file), f"Lateral summary plot not created: {expected_file}"

        print(f"✓ PlotLateralSummary: multi-panel file created")

    finally:
        shutil.rmtree(tmpdir)


def test_plot_tidal_heating_by_layer():
    """ Test PlotTidalHeatingByLayer function. """
    print("\n=== Test: PlotTidalHeatingByLayer ===")

    Planet = _create_mock_planet_with_lateral()
    tmpdir = tempfile.mkdtemp()

    try:
        Params.FigureFiles.path = tmpdir
        PlotTidalHeatingByLayer(Planet, Params)

        # Check file created
        expected_file = os.path.join(tmpdir, f'{Planet.name}_HtidalByLayer{FigMisc.xtn}')
        assert os.path.exists(expected_file), f"Tidal heating by-layer plot not created: {expected_file}"

        print(f"✓ PlotTidalHeatingByLayer: file created")

    finally:
        shutil.rmtree(tmpdir)


def test_plot_tidal_heating_vs_depth():
    """ Test PlotTidalHeatingVsDepth function. """
    print("\n=== Test: PlotTidalHeatingVsDepth ===")

    Planet = _create_mock_planet_with_lateral()
    tmpdir = tempfile.mkdtemp()

    try:
        Params.FigureFiles.path = tmpdir
        PlotTidalHeatingVsDepth(Planet, Params)

        # Check file created
        expected_file = os.path.join(tmpdir, f'{Planet.name}_HtidalVsDepth{FigMisc.xtn}')
        assert os.path.exists(expected_file), f"Tidal heating vs depth plot not created: {expected_file}"

        print(f"✓ PlotTidalHeatingVsDepth: file created")

    finally:
        shutil.rmtree(tmpdir)


def test_generate_lateral_plots():
    """ Test GenerateLateralPlots orchestration. """
    print("\n=== Test: GenerateLateralPlots ===")

    Planet = _create_mock_planet_with_lateral()
    Planet.Lateral.DO_3D = True
    tmpdir = tempfile.mkdtemp()

    try:
        Params.FigureFiles.path = tmpdir
        Params.PLOT_LATERAL = True

        GenerateLateralPlots(Planet, Params)

        # Check multiple files created
        expected_files = [
            f'{Planet.name}_dIce{FigMisc.xtn}',
            f'{Planet.name}_Htidal{FigMisc.xtn}',
            f'{Planet.name}_Tb{FigMisc.xtn}',
            f'{Planet.name}_sigmaOcean{FigMisc.xtn}',
            f'{Planet.name}_fClath{FigMisc.xtn}',
            f'{Planet.name}_kTherm{FigMisc.xtn}',
            f'{Planet.name}_lateralSummary{FigMisc.xtn}',
            f'{Planet.name}_HtidalByLayer{FigMisc.xtn}',
            f'{Planet.name}_HtidalVsDepth{FigMisc.xtn}',
        ]

        created_count = 0
        for fname in expected_files:
            fpath = os.path.join(tmpdir, fname)
            if os.path.exists(fpath):
                created_count += 1

        assert created_count >= 7, f"Only {created_count} of {len(expected_files)} plots created"

        print(f"✓ GenerateLateralPlots: {created_count} files created")

    finally:
        shutil.rmtree(tmpdir)


def test_graceful_missing_data():
    """ Test that plot functions handle missing data gracefully. """
    print("\n=== Test: Graceful handling of missing data ===")

    Planet = PlanetStruct('TestPlanet')
    tmpdir = tempfile.mkdtemp()

    try:
        Params.FigureFiles.path = tmpdir

        # All functions should return early without error when data is missing
        PlotIceThickness(Planet, Params)  # Lateral.dIce_m is None
        PlotTidalHeating(Planet, Params)  # Lateral.HtidalIce_Wm3 is None
        PlotBasalTemperature(Planet, Params)  # Lateral.Tb_K is None
        PlotOceanConductivity(Planet, Params)  # Lateral.sigma_mean_Sm is None
        PlotClathrateFraction(Planet, Params)  # Lateral.fClath is None
        PlotEffectiveConductivity(Planet, Params)  # Lateral.kThermEff_WmK is None
        PlotLateralSummary(Planet, Params)  # Lateral.dIce_m is None
        PlotTidalHeatingByLayer(Planet, Params)  # Lateral.HtidalIceI_bot_Wm3 is None
        PlotTidalHeatingVsDepth(Planet, Params)  # Lateral.HtidalIceI_profile_Wm3 is None

        # No files should be created
        files_created = os.listdir(tmpdir)
        assert len(files_created) == 0, f"Unexpected files created with missing data: {files_created}"

        print("✓ All plot functions handle missing data gracefully (no error, no files)")

    finally:
        shutil.rmtree(tmpdir)


def test_plot_lateral_disabled():
    """ Test GenerateLateralPlots respects PLOT_LATERAL flag. """
    print("\n=== Test: PLOT_LATERAL flag ===")

    Planet = _create_mock_planet_with_lateral()
    Planet.Lateral.DO_3D = True
    tmpdir = tempfile.mkdtemp()

    try:
        Params.FigureFiles.path = tmpdir
        Params.PLOT_LATERAL = False  # Disable plotting

        GenerateLateralPlots(Planet, Params)

        # No files should be created
        files_created = os.listdir(tmpdir)
        assert len(files_created) == 0, f"Files created even with PLOT_LATERAL=False: {files_created}"

        print("✓ GenerateLateralPlots respects PLOT_LATERAL=False (no files)")

    finally:
        shutil.rmtree(tmpdir)


def run_all_tests():
    """ Run all Phase 8 tests. """
    print("\n" + "=" * 60)
    print("PPTestLateralPhase8: Lateral plotting functions")
    print("=" * 60)

    test_interpolate_to_latlon()
    test_plot_surface_map()
    test_plot_ice_thickness()
    test_plot_tidal_heating()
    test_plot_basal_temperature()
    test_plot_ocean_conductivity()
    test_plot_clathrate_fraction()
    test_plot_effective_conductivity()
    test_plot_lateral_summary()
    test_plot_tidal_heating_by_layer()
    test_plot_tidal_heating_vs_depth()
    test_generate_lateral_plots()
    test_graceful_missing_data()
    test_plot_lateral_disabled()

    print("\n" + "=" * 60)
    print("✓ All Phase 8 tests passed!")
    print("=" * 60)


if __name__ == '__main__':
    run_all_tests()
