"""
PPTestLateralPhase2
Test for Phase 2 of 3D lateral structure port: Spatial grid and spherical harmonics.
Tests grid initialization, area normalization, and SH transforms.

Run with: python PlanetProfile/Test/PPTestLateralPhase2.py
"""

import sys
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, LateralSubstruct
from PlanetProfile.Lateral.SpatialGrid import (
    InitGrid, SHtoGrid, GridToSH, IntegrateOverSphere, _assoc_legendre_4pi
)

# Check if healpy is available
try:
    import healpy as hp
    HEALPY_AVAILABLE = True
except ImportError:
    HEALPY_AVAILABLE = False
    print("Warning: healpy not available. Skipping HEALPix tests.")


def test_healpix_grid():
    """Test HEALPix grid initialization."""
    if not HEALPY_AVAILABLE:
        print("  Skipping HEALPix tests (healpy not installed)")
        return

    print("Testing HEALPix grid initialization...")

    lateral = LateralSubstruct()
    lateral.gridType = 'healpix'
    lateral.nSide = 8

    InitGrid(lateral)

    # Verify nPix = 12 * nSide^2
    expected_nPix = 12 * lateral.nSide**2
    assert lateral.nPix == expected_nPix, f"nPix should be {expected_nPix}, got {lateral.nPix}"
    assert lateral.nPix == 768, "nSide=8 should give 768 pixels"

    # Verify arrays are set
    assert lateral.theta_rad is not None, "theta_rad should be set"
    assert lateral.phi_rad is not None, "phi_rad should be set"
    assert len(lateral.theta_rad) == lateral.nPix, "theta_rad length mismatch"
    assert len(lateral.phi_rad) == lateral.nPix, "phi_rad length mismatch"

    # Verify theta range [0, pi]
    assert np.min(lateral.theta_rad) >= 0, "theta should be >= 0"
    assert np.max(lateral.theta_rad) <= np.pi, "theta should be <= pi"

    # Verify phi range [0, 2pi)
    assert np.min(lateral.phi_rad) >= 0, "phi should be >= 0"
    assert np.max(lateral.phi_rad) < 2 * np.pi, "phi should be < 2pi"

    # Verify pixel area (equal-area property)
    assert isinstance(lateral.pixArea_sr, (float, np.floating)), "HEALPix pixArea_sr should be scalar"
    expected_area = 4 * np.pi / lateral.nPix
    assert np.abs(lateral.pixArea_sr - expected_area) < 1e-10, "Pixel area should be 4pi/nPix"

    # Verify total area = 4pi
    total_area = lateral.nPix * lateral.pixArea_sr
    assert np.abs(total_area - 4 * np.pi) < 1e-8, f"Total area should be 4pi, got {total_area}"

    print("  ✓ HEALPix grid initialized correctly")
    print(f"  ✓ nPix = {lateral.nPix}, area per pixel = {lateral.pixArea_sr:.6e} sr")
    print(f"  ✓ Total area = {total_area:.10f} (4π = {4*np.pi:.10f})")


def test_latlon_grid():
    """Test lat-lon grid initialization."""
    print("Testing lat-lon grid initialization...")

    lateral = LateralSubstruct()
    lateral.gridType = 'latlon'
    lateral.nLat = 37
    lateral.nLon = 72

    InitGrid(lateral)

    # Verify nPix = nLat * nLon
    expected_nPix = lateral.nLat * lateral.nLon
    assert lateral.nPix == expected_nPix, f"nPix should be {expected_nPix}, got {lateral.nPix}"
    assert lateral.nPix == 2664, "37×72 should give 2664 pixels"

    # Verify arrays are set
    assert lateral.theta_rad is not None, "theta_rad should be set"
    assert lateral.phi_rad is not None, "phi_rad should be set"
    assert len(lateral.theta_rad) == lateral.nPix, "theta_rad length mismatch"
    assert len(lateral.phi_rad) == lateral.nPix, "phi_rad length mismatch"

    # Verify pixel area is array (non-equal-area)
    assert isinstance(lateral.pixArea_sr, np.ndarray), "latlon pixArea_sr should be array"
    assert len(lateral.pixArea_sr) == lateral.nPix, "pixArea_sr length mismatch"

    # Verify areas vary with latitude (sin(theta) weighting)
    area_min = np.min(lateral.pixArea_sr)
    area_max = np.max(lateral.pixArea_sr)
    assert area_max > area_min * 10, "lat-lon areas should vary significantly (poles vs equator)"

    # Verify total area = 4pi (within 0.1% for numerical integration)
    total_area = np.sum(lateral.pixArea_sr)
    rel_error = np.abs(total_area - 4 * np.pi) / (4 * np.pi)
    assert rel_error < 1e-3, f"Total area should be 4pi, got {total_area} (error {rel_error:.2%})"

    print("  ✓ Lat-lon grid initialized correctly")
    print(f"  ✓ nPix = {lateral.nPix}, area range = [{area_min:.6e}, {area_max:.6e}] sr")
    print(f"  ✓ Area ratio (equator/pole) = {area_max/area_min:.1f}×")
    print(f"  ✓ Total area = {total_area:.10f} (4π = {4*np.pi:.10f})")


def test_spherical_harmonic_roundtrip():
    """Test that GridToSH → SHtoGrid round-trip recovers original field."""
    print("Testing spherical harmonic round-trip...")

    # Use lat-lon grid with higher resolution for better SH accuracy
    lateral = LateralSubstruct()
    lateral.gridType = 'latlon'
    lateral.nLat = 37  # Higher resolution for better numerical integration
    lateral.nLon = 72
    InitGrid(lateral)

    # Create a simple test field: Y_2,1 (degree 2, order 1)
    # This is a well-defined spherical harmonic
    pMax = 4
    Cpq_orig = np.zeros((pMax + 1, pMax + 1))
    Spq_orig = np.zeros((pMax + 1, pMax + 1))

    # Set a few coefficients
    Cpq_orig[0, 0] = 100.0  # Mean value
    Cpq_orig[2, 1] = 10.0   # Degree 2, order 1 cosine
    Spq_orig[2, 1] = 5.0    # Degree 2, order 1 sine
    Cpq_orig[4, 2] = 2.0    # Degree 4, order 2 cosine

    # Forward: SH → Grid
    field = SHtoGrid(Cpq_orig, Spq_orig, pMax, lateral.theta_rad, lateral.phi_rad)

    # Backward: Grid → SH
    Cpq_recovered, Spq_recovered = GridToSH(field, lateral.theta_rad, lateral.phi_rad,
                                             lateral.pixArea_sr, pMax)

    # Check recovery (use reasonable tolerances for numerical SH integration on lat-lon grids)
    # Compute relative error only where original coefficients are nonzero
    Cpq_abs_error = np.abs(Cpq_orig - Cpq_recovered)
    Spq_abs_error = np.abs(Spq_orig - Spq_recovered)

    # Find nonzero elements
    Cpq_nonzero_mask = np.abs(Cpq_orig) > 1e-10
    Spq_nonzero_mask = np.abs(Spq_orig) > 1e-10

    if np.any(Cpq_nonzero_mask):
        Cpq_rel_errors = Cpq_abs_error[Cpq_nonzero_mask] / np.abs(Cpq_orig[Cpq_nonzero_mask])
        Cpq_max_rel_error = np.max(Cpq_rel_errors)
    else:
        Cpq_max_rel_error = 0.0

    if np.any(Spq_nonzero_mask):
        Spq_rel_errors = Spq_abs_error[Spq_nonzero_mask] / np.abs(Spq_orig[Spq_nonzero_mask])
        Spq_max_rel_error = np.max(Spq_rel_errors)
    else:
        Spq_max_rel_error = 0.0

    # Lat-lon grids with moderate resolution have ~5% SH integration error
    assert Cpq_max_rel_error < 0.10, f"Cpq recovery error too large: {Cpq_max_rel_error:.2%}"
    assert Spq_max_rel_error < 0.10, f"Spq recovery error too large: {Spq_max_rel_error:.2%}"

    # Forward again to verify field recovery
    field_recovered = SHtoGrid(Cpq_recovered, Spq_recovered, pMax,
                               lateral.theta_rad, lateral.phi_rad)
    field_rel_error = np.max(np.abs(field - field_recovered)) / np.max(np.abs(field))
    assert field_rel_error < 0.01, f"Field recovery error too large: {field_rel_error:.2%}"

    print("  ✓ Spherical harmonic round-trip successful")
    print(f"  ✓ Cpq max relative error: {Cpq_max_rel_error:.2%}")
    print(f"  ✓ Spq max relative error: {Spq_max_rel_error:.2%}")
    print(f"  ✓ Field relative error: {field_rel_error:.2%}")


def test_legendre_normalization():
    """Test 4π normalization of full spherical harmonics (including angular part)."""
    print("Testing 4π-normalized spherical harmonics...")

    # Create a lat-lon grid for integration
    nLat = 91
    nLon = 180
    lateral = LateralSubstruct()
    lateral.gridType = 'latlon'
    lateral.nLat = nLat
    lateral.nLon = nLon
    InitGrid(lateral)

    cos_theta = np.cos(lateral.theta_rad)

    # Test orthogonality: integral(Y_nm * Y_n'm' * dOmega) = 4pi * delta
    # For real SH: Y_nm = P_nm(cos θ) × cos(mφ) or sin(mφ)
    n1, m1 = 2, 1
    n2, m2 = 2, 1  # Same as n1, m1 → should give 4π

    P1 = _assoc_legendre_4pi(n1, m1, cos_theta)
    P2 = _assoc_legendre_4pi(n2, m2, cos_theta)

    # Full spherical harmonic includes angular part
    Y1_cos = P1 * np.cos(m1 * lateral.phi_rad)
    Y2_cos = P2 * np.cos(m2 * lateral.phi_rad)

    # Integrate Y1 * Y2 over sphere
    integral = np.sum(Y1_cos * Y2_cos * lateral.pixArea_sr)

    expected = 4 * np.pi
    rel_error = np.abs(integral - expected) / expected

    # Lat-lon grids have ~1% integration error
    assert rel_error < 0.05, f"Normalization error too large: {rel_error:.2%}"

    # Test different degrees should give ~0 (orthogonality)
    n3, m3 = 3, 1
    P3 = _assoc_legendre_4pi(n3, m3, cos_theta)
    Y3_cos = P3 * np.cos(m3 * lateral.phi_rad)
    integral_orthog = np.sum(Y1_cos * Y3_cos * lateral.pixArea_sr)

    # Should be ~0 but allow for numerical integration error
    assert np.abs(integral_orthog) < 0.5, f"Orthogonality violated: {integral_orthog}"

    print("  ✓ 4π-normalization verified")
    print(f"  ✓ ∫ Y_21² dΩ = {integral:.6f} (expected {expected:.6f}, error {rel_error:.2%})")
    print(f"  ✓ ∫ Y_21 × Y_31 dΩ = {integral_orthog:.6f} (should be ~0)")


def test_integration():
    """Test integration over sphere."""
    print("Testing sphere integration...")

    lateral = LateralSubstruct()
    lateral.gridType = 'latlon'
    lateral.nLat = 37
    lateral.nLon = 72
    InitGrid(lateral)

    # Test 1: Constant field should integrate to 4πR²
    R_m = 1e6  # 1000 km
    const_field = np.ones(lateral.nPix) * 100.0  # 100 units

    integral = IntegrateOverSphere(const_field, lateral.pixArea_sr, R_m)
    expected = 4 * np.pi * R_m**2 * 100.0

    error = np.abs(integral - expected) / expected
    # Lat-lon grids have ~0.1% integration error
    assert error < 1e-3, f"Integration error too large: {error:.2%}"

    print("  ✓ Integration of constant field correct")
    print(f"  ✓ Integral = {integral:.6e}, expected = {expected:.6e}, error = {error:.2%}")

    # Test 2: Field varying with latitude
    # Field = sin²(θ) should integrate to (8π/3) R²
    field_sin2 = np.sin(lateral.theta_rad)**2
    integral_sin2 = IntegrateOverSphere(field_sin2, lateral.pixArea_sr, R_m)
    expected_sin2 = (8 * np.pi / 3) * R_m**2

    error_sin2 = np.abs(integral_sin2 - expected_sin2) / expected_sin2
    assert error_sin2 < 0.01, f"sin²(θ) integration error: {error_sin2:.2%}"

    print("  ✓ Integration of sin²(θ) correct")
    print(f"  ✓ Integral = {integral_sin2:.6e}, expected = {expected_sin2:.6e}, error = {error_sin2:.2%}")


def test_scientific_understanding():
    """Document scientific understanding from Phase 2."""
    print("\nScientific understanding (Phase 2):")
    print("  HEALPix vs lat-lon trade-offs:")
    print("    - HEALPix: equal-area (4π/nPix per pixel), better for SH, no pole issues")
    print("    - Lat-lon: simpler, non-equal-area (equator ~100× larger than poles)")
    print("  4π normalization:")
    print("    - ∫|Y_pq|² dΩ = 4π (full sphere solid angle)")
    print("    - C_00 = mean value of field (degree 0 = constant)")
    print("    - Consistent with MoonMag, NASA gravity fields")
    print("  Pixel areas:")
    print("    - HEALPix: scalar, all equal")
    print("    - Lat-lon: array, A(θ) ∝ sin(θ) × ΔθΔφ")
    print("    - Integration: Σ f(θ,φ) × A_i → ∫∫ f dΩ")
    print("  Spherical harmonics:")
    print("    - Y_pq = N_pq × P_pq(cos θ) × [cos(qφ) or sin(qφ)]")
    print("    - Degree p: spatial wavelength (p=0 constant, p=2 quadrupole)")
    print("    - Order q: longitudinal variation (q=0 zonal, q>0 sectoral/tesseral)")


if __name__ == '__main__':
    print("="*70)
    print("PPTestLateralPhase2: Testing spatial grid and spherical harmonics")
    print("="*70)

    try:
        if HEALPY_AVAILABLE:
            test_healpix_grid()
        test_latlon_grid()
        test_spherical_harmonic_roundtrip()
        test_legendre_normalization()
        test_integration()
        test_scientific_understanding()

        print("\n" + "="*70)
        print("✓ All Phase 2 tests passed!")
        print("="*70)
        sys.exit(0)

    except AssertionError as e:
        print(f"\n✗ Test failed: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
