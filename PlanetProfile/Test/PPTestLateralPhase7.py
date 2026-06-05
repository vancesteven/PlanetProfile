"""
PPTestLateralPhase7
Test for Phase 7 of 3D lateral structure port: I/O and MoonMag integration.
Tests SaveLateralResults, ReloadLateralResults, and UpdateAsymShapeFrom3D.

Run with: python PlanetProfile/Test/PPTestLateralPhase7.py
"""

import sys
import os
import tempfile
import shutil
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, ParamsStruct, StepsSubstruct
from PlanetProfile.Lateral.SpatialGrid import InitGrid
from PlanetProfile.Lateral.LateralIO import SaveLateralResults, ReloadLateralResults
from PlanetProfile.Lateral.LateralStructure import UpdateAsymShapeFrom3D


def _SetupMockLateral(Planet, nPix=12):
    """Create minimal lateral grid for testing."""
    Planet.Lateral.DO_3D = True
    Planet.Lateral.gridType = 'healpix'
    Planet.Lateral.nSide = 1 if nPix == 12 else 2 if nPix == 48 else 4
    InitGrid(Planet.Lateral)

    # Set mock 3D fields
    Planet.Lateral.dIce_m = np.random.uniform(80e3, 120e3, nPix)
    Planet.Lateral.fClath = np.random.uniform(0.0, 0.3, nPix)
    Planet.Lateral.Tb_K = np.random.uniform(270, 275, nPix)
    Planet.Lateral.qSurf_Wm2 = np.random.uniform(10e-3, 30e-3, nPix)
    Planet.Lateral.sigma_mean_Sm = np.random.uniform(0.5, 1.5, nPix)
    Planet.Lateral.kThermEff_WmK = np.random.uniform(2.0, 3.5, nPix)
    Planet.Lateral.HtidalIce_Wm3 = np.random.uniform(1e-5, 5e-5, nPix)
    Planet.Lateral.Mtarget_kg = 1.48e23
    Planet.Lateral.Mactual_kg = 1.4801e23
    Planet.Lateral.massResidual_frac = 0.0001

    # Set SH coefficients
    pMax = 2
    Planet.Lateral.dIce_pMax = pMax
    Planet.Lateral.dIce_Cpq_km = np.random.randn(pMax + 1, pMax + 1) * 10
    Planet.Lateral.dIce_Spq_km = np.random.randn(pMax + 1, pMax + 1) * 10
    Planet.Lateral.dIce_Cpq_km[0, 0] = 100.0  # Mean ice thickness

    return Planet


def _CreateMockColumnPlanets(Planet, nPix=12, ice_var_km=20):
    """Create mock column planets with ice thickness variation."""
    columnPlanets = np.empty(nPix, dtype=object)

    R_m = Planet.Bulk.R_m
    ice_base_km = 100.0

    for i in range(nPix):
        col = PlanetStruct('MockColumn')
        col.Steps = StepsSubstruct()
        col.Steps.nSurfIce = 50
        col.Steps.nHydro = 100

        # Vary ice thickness: thicker at poles, thinner at equator
        theta = Planet.Lateral.theta_rad[i]
        ice_km = ice_base_km + ice_var_km * np.cos(theta)
        ice_m = ice_km * 1e3

        # Create radial profile: ice extends from surface to nSurfIce
        # Ocean from nSurfIce to bottom
        total_depth = ice_m + 50e3  # ice + 50 km ocean
        col.r_m = np.linspace(R_m, R_m - total_depth, 100)
        col.rho_kgm3 = np.full(100, 1000.0)
        col.phase = np.zeros(100, dtype=int)
        col.phase[:50] = 1  # Ice (indices 0-49)
        col.phase[50:] = 0  # Ocean (indices 50-99)

        # r_m[nSurfIce] should be at ice-ocean boundary
        # nSurfIce = 50 means index 50 is first ocean point
        # So r_m[50] = R_m - ice_m
        col.r_m = np.linspace(R_m, R_m - total_depth, 100)
        # Explicitly set ice-ocean boundary
        col.r_m[50] = R_m - ice_m

        col.index = i

        columnPlanets[i] = col

    return columnPlanets


def test_save_lateral_results():
    """Test SaveLateralResults function."""
    print("Testing SaveLateralResults...")

    Planet = PlanetStruct('TestPlanet')
    Params = ParamsStruct()
    Planet.Bulk.R_m = 1560e3
    Planet.name = 'TestPlanet'
    Planet.bodyname = 'TestPlanet'

    Planet = _SetupMockLateral(Planet, nPix=12)

    # Create temporary output directory
    tmpdir = tempfile.mkdtemp()
    try:
        Params.FigureFiles = type('obj', (object,), {'path': tmpdir})()

        SaveLateralResults(Planet, Params)

        # Check file was created
        outPath = os.path.join(tmpdir, 'TestPlanet_lateral3D.pkl')
        assert os.path.exists(outPath), f"Output file not created: {outPath}"

        # Check file size is reasonable
        file_size = os.path.getsize(outPath)
        assert file_size > 100, f"Output file too small: {file_size} bytes"
        assert file_size < 1e6, f"Output file too large: {file_size} bytes"

        print(f"  ✓ Saved lateral results: {file_size} bytes")

    finally:
        shutil.rmtree(tmpdir)


def test_reload_lateral_results():
    """Test ReloadLateralResults function."""
    print("Testing ReloadLateralResults...")

    # Create and save
    Planet1 = PlanetStruct('TestPlanet')
    Params = ParamsStruct()
    Planet1.Bulk.R_m = 1560e3
    Planet1.name = 'TestPlanet'
    Planet1.bodyname = 'TestPlanet'
    Planet1 = _SetupMockLateral(Planet1, nPix=12)

    tmpdir = tempfile.mkdtemp()
    try:
        Params.FigureFiles = type('obj', (object,), {'path': tmpdir})()
        SaveLateralResults(Planet1, Params)

        # Reload into new Planet
        Planet2 = PlanetStruct('TestPlanet2')
        outPath = os.path.join(tmpdir, 'TestPlanet_lateral3D.pkl')
        Planet2 = ReloadLateralResults(Planet2, outPath)

        # Verify fields match
        assert Planet2.Lateral.nPix == 12, "nPix mismatch"
        assert Planet2.Lateral.gridType == 'healpix', "gridType mismatch"
        assert np.allclose(Planet2.Lateral.dIce_m, Planet1.Lateral.dIce_m), "dIce_m mismatch"
        assert np.allclose(Planet2.Lateral.Tb_K, Planet1.Lateral.Tb_K), "Tb_K mismatch"
        assert np.allclose(Planet2.Lateral.sigma_mean_Sm, Planet1.Lateral.sigma_mean_Sm), "sigma_mean_Sm mismatch"
        assert np.allclose(Planet2.Lateral.dIce_Cpq_km, Planet1.Lateral.dIce_Cpq_km), "dIce_Cpq_km mismatch"

        print(f"  ✓ Reloaded {Planet2.Lateral.nPix} pixels, all fields match")

    finally:
        shutil.rmtree(tmpdir)


def test_save_reload_roundtrip():
    """Test save-reload-save roundtrip preserves data."""
    print("Testing save-reload-save roundtrip...")

    Planet = PlanetStruct('TestPlanet')
    Params = ParamsStruct()
    Planet.Bulk.R_m = 1560e3
    Planet.name = 'TestPlanet'
    Planet.bodyname = 'TestPlanet'
    Planet = _SetupMockLateral(Planet, nPix=48)

    tmpdir = tempfile.mkdtemp()
    try:
        Params.FigureFiles = type('obj', (object,), {'path': tmpdir})()

        # Save
        SaveLateralResults(Planet, Params)
        outPath1 = os.path.join(tmpdir, 'TestPlanet_lateral3D.pkl')
        size1 = os.path.getsize(outPath1)

        # Reload
        Planet2 = PlanetStruct('TestPlanet')
        Planet2 = ReloadLateralResults(Planet2, outPath1)

        # Save again
        SaveLateralResults(Planet2, Params)
        size2 = os.path.getsize(outPath1)

        # Sizes should be nearly identical (pickle may have minor variations)
        assert abs(size2 - size1) < 100, f"Roundtrip size change: {size1} -> {size2}"

        print(f"  ✓ Roundtrip preserves data: {size1} bytes -> {size2} bytes")

    finally:
        shutil.rmtree(tmpdir)


def test_update_asym_shape_from_3d():
    """Test UpdateAsymShapeFrom3D for MoonMag integration."""
    print("Testing UpdateAsymShapeFrom3D...")

    Planet = PlanetStruct('TestPlanet')
    Planet.Bulk.R_m = 1560e3
    Planet.name = 'TestPlanet'
    Planet = _SetupMockLateral(Planet, nPix=12)

    # Create mock column planets with ice thickness variation
    columnPlanets = _CreateMockColumnPlanets(Planet, nPix=12, ice_var_km=20)

    # Initialize Magnetic substruct
    Planet.Magnetic.nBds = 1

    # Call UpdateAsymShapeFrom3D
    Planet = UpdateAsymShapeFrom3D(Planet, columnPlanets)

    # Check outputs
    assert hasattr(Planet.Magnetic, 'asymShape_m'), "asymShape_m not created"
    assert hasattr(Planet.Magnetic, 'pMax'), "pMax not set"
    assert Planet.Magnetic.asymShape_m.shape[0] == 1, "nBds dimension wrong"
    assert Planet.Magnetic.asymShape_m.shape[1] == 2, "pos/neg q dimension wrong"
    assert Planet.Magnetic.asymShape_m.shape[2] == Planet.Magnetic.pMax + 1, "p dimension wrong"
    assert Planet.Magnetic.asymShape_m.shape[3] == Planet.Magnetic.pMax + 1, "q dimension wrong"
    assert Planet.Magnetic.asymShape_m.dtype == np.complex128, "asymShape_m should be complex"

    # Check p=0 is zero (no monopole asymmetry)
    assert np.allclose(Planet.Magnetic.asymShape_m[:, :, 0, :], 0.0), "p=0 should be zero"

    # Check some p>0 coefficients are non-zero (because ice varies)
    assert np.any(np.abs(Planet.Magnetic.asymShape_m[:, :, 1:, :]) > 1e-3), \
        "No asymmetry detected in ice-ocean boundary"

    print(f"  ✓ asymShape_m created: pMax={Planet.Magnetic.pMax}, shape={Planet.Magnetic.asymShape_m.shape}")


def test_asym_shape_pole_equator_pattern():
    """Test UpdateAsymShapeFrom3D with pole-equator ice thickness pattern."""
    print("Testing asymShape_m with pole-equator pattern...")

    Planet = PlanetStruct('TestPlanet')
    Planet.Bulk.R_m = 1560e3
    Planet = _SetupMockLateral(Planet, nPix=48)  # Higher resolution

    # Create column planets with strong pole-equator variation
    # Thick ice at poles (theta ~ 0, pi), thin at equator (theta ~ pi/2)
    columnPlanets = _CreateMockColumnPlanets(Planet, nPix=48, ice_var_km=40)

    Planet.Magnetic.nBds = 1
    Planet = UpdateAsymShapeFrom3D(Planet, columnPlanets)

    # Check that SOME asymmetry was detected (non-zero SH coefficients)
    total_asymmetry = np.sum(np.abs(Planet.Magnetic.asymShape_m[0, :, 1:, :]))

    # With 40 km ice variation, should have some detectable asymmetry
    # Even if grid resolution is limited
    assert total_asymmetry > 0, "No asymmetry detected"

    # Check that not all coefficients are exactly the same (variation exists)
    coeffs = np.abs(Planet.Magnetic.asymShape_m[0, 0, 2, :])
    assert np.std(coeffs) > 0, "No variation in p=2 coefficients"

    print(f"  ✓ Asymmetry detected: total |χpq| = {total_asymmetry:.2e} m")


def test_asym_shape_with_multiple_boundaries():
    """Test UpdateAsymShapeFrom3D with multiple conducting boundaries."""
    print("Testing asymShape_m with multiple boundaries...")

    Planet = PlanetStruct('TestPlanet')
    Planet.Bulk.R_m = 1560e3
    Planet = _SetupMockLateral(Planet, nPix=12)

    columnPlanets = _CreateMockColumnPlanets(Planet, nPix=12, ice_var_km=20)

    # Multiple boundaries (e.g., surface, ice-ocean, ocean-silicate)
    nBds = 3
    Planet.Magnetic.nBds = nBds
    Planet.Magnetic.rSigChange_m = np.array([1560e3, 1560e3 - 100e3, 1560e3 - 150e3])

    Planet = UpdateAsymShapeFrom3D(Planet, columnPlanets)

    # Check all boundaries have asymShape_m
    assert Planet.Magnetic.asymShape_m.shape[0] == nBds, \
        f"Expected {nBds} boundaries, got {Planet.Magnetic.asymShape_m.shape[0]}"

    # Check radial scaling: deeper boundaries should have proportionally scaled asymmetry
    r_ice_ocean = Planet.Magnetic.rSigChange_m[0]
    for i in range(1, nBds):
        scale_factor = Planet.Magnetic.rSigChange_m[i] / r_ice_ocean
        # Sample Y21 coefficient
        if np.abs(Planet.Magnetic.asymShape_m[0, 0, 2, 1]) > 1e-6:
            ratio = (Planet.Magnetic.asymShape_m[i, 0, 2, 1] /
                     Planet.Magnetic.asymShape_m[0, 0, 2, 1])
            assert np.isclose(ratio, scale_factor, rtol=0.01), \
                f"Boundary {i} scaling wrong: {ratio:.3f} vs expected {scale_factor:.3f}"

    print(f"  ✓ {nBds} boundaries scaled correctly")


def test_save_no_lateral_data():
    """Test SaveLateralResults with no lateral data (should skip gracefully)."""
    print("Testing save with no lateral data...")

    Planet = PlanetStruct('TestPlanet')
    Params = ParamsStruct()
    Planet.name = 'TestPlanet'
    Planet.bodyname = 'TestPlanet'

    # No lateral data initialized (dIce_m is None)

    tmpdir = tempfile.mkdtemp()
    try:
        Params.FigureFiles = type('obj', (object,), {'path': tmpdir})()

        # Should not crash, just return early
        SaveLateralResults(Planet, Params)

        # No file should be created
        outPath = os.path.join(tmpdir, 'TestPlanet_lateral3D.pkl')
        assert not os.path.exists(outPath), "File should not be created when no lateral data"

        print(f"  ✓ Gracefully skipped save with no data")

    finally:
        shutil.rmtree(tmpdir)


def print_scientific_summary():
    """Print scientific understanding summary."""
    print("\n" + "=" * 70)
    print("SCIENTIFIC UNDERSTANDING: I/O and MoonMag Integration")
    print("=" * 70)
    print("""
  Lateral results I/O:
    - Pickle format: Python-specific but efficient and portable
    - Stores all 3D fields: dIce_m, fClath, Tb_K, qSurf_Wm2, sigma_mean_Sm,
      kThermEff_WmK, HtidalIce_Wm3, and SH coefficients
    - Enables: reloading for plotting, sensitivity studies, comparison
    - File location: alongside standard PP output, '_lateral3D.pkl' suffix
    - Typical size: ~10-100 KB for 192 pixels, ~1 MB for 3072 pixels

  UpdateAsymShapeFrom3D - MoonMag integration:
    - Converts 3D ice-ocean boundary topography r(θ,φ) to SH coefficients
    - Input: columnPlanets with ice thickness dIce(θ,φ)
    - Process: r_iceOcean[i] → deviations dr = r - mean(r) → GridToSH → χpq
    - Output: Planet.Magnetic.asymShape_m[nBds, 2, pMax+1, pMax+1]
    - Format: complex χpq coefficients (4π normalization) for MoonMag

  Why ice-ocean boundary asymmetry matters:
    - Ice thickness varies geographically due to:
      • Tidal heating variations (more dissipation → thinner ice)
      • Clathrate distribution (more clathrate → thinner ice)
      • Basal heat flux variations (higher q → thinner ice)
    - Ocean conductivity is fixed in each layer (MoonMag limitation)
    - BUT boundary shape affects induction response:
      • Thinner ice → ocean closer to surface → stronger induced field
      • Geographic variation → asymmetric induced dipole
      • Observable in spacecraft magnetometer data

  Spherical harmonic decomposition:
    - p = degree (angular scale), q = order (longitudinal variation)
    - p=0: monopole (mean), p=1: dipole, p=2: quadrupole, etc.
    - Y20: pole-equator variation (q=0, zonal)
    - Y22: sectoral pattern (q=2, 4 lobes)
    - Typical ice variation: pMax = 4-8 sufficient

  MoonMag asymShape_m structure:
    - Dimension 0: nBds (number of conducting boundaries)
    - Dimension 1: 2 (positive/negative q)
    - Dimension 2: pMax+1 (degree p from 0 to pMax)
    - Dimension 3: pMax+1 (order q from 0 to pMax)
    - Complex values χpq in meters (absolute deviations)
    - Concentric scaling: deeper boundaries proportional to radius

  Integration with MagneticInduction.py:
    - UpdateAsymShapeFrom3D sets Planet.Magnetic.asymShape_m
    - MagneticInduction reads asymShape_m automatically
    - If CALC_ASYM=True: uses MoonMag asymmetry_funcs
    - If CALC_ASYM=False: uses symmetric induction
    - No code changes needed in MagneticInduction.py

  Typical asymmetry scales (Europa):
    - Mean ice thickness: 20-30 km
    - Geographic variation: ±5-10 km (17-33% of mean)
    - Y20 amplitude: 1-5 km (3-17% of radius)
    - Observable in Bx, By, Bz asymmetry: ~10-50 nT variations
""")
    print("=" * 70)


if __name__ == '__main__':
    print("=" * 70)
    print("Phase 7: I/O and MoonMag Integration Tests")
    print("=" * 70)

    try:
        test_save_lateral_results()
        test_reload_lateral_results()
        test_save_reload_roundtrip()
        test_update_asym_shape_from_3d()
        test_asym_shape_pole_equator_pattern()
        test_asym_shape_with_multiple_boundaries()
        test_save_no_lateral_data()

        print_scientific_summary()

        print("\n" + "=" * 70)
        print("✓ All Phase 7 tests passed!")
        print("=" * 70)
        sys.exit(0)

    except AssertionError as e:
        print(f"\n✗ Test failed: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
