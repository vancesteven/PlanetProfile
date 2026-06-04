"""
PPTestLateralPhase6
Test for Phase 6 of 3D lateral structure port: Clathrate lateral variation.
Tests InitClathrateLateral, ComputeEffectiveConductivity, and EstimateIceThicknessFromClathrate.

Run with: python PlanetProfile/Test/PPTestLateralPhase6.py
"""

import sys
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct
from PlanetProfile.Lateral.ClathrateLateral import (
    InitClathrateLateral,
    ComputeEffectiveConductivity,
    EstimateIceThicknessFromClathrate,
    K_CLATH_WmK
)


def _SetupMockLateralGrid(Planet, nPix=12):
    """Create a minimal lateral grid for testing without full InitLateralStructure."""
    from PlanetProfile.Lateral.SpatialGrid import InitGrid

    Planet.Lateral.DO_3D = True
    Planet.Lateral.gridType = 'healpix'

    # HEALPix requires nSide parameter (nPix = 12 * nSide^2)
    # nPix=12 -> nSide=1, nPix=48 -> nSide=2, nPix=192 -> nSide=4
    if nPix == 12:
        Planet.Lateral.nSide = 1
    elif nPix == 48:
        Planet.Lateral.nSide = 2
    elif nPix == 192:
        Planet.Lateral.nSide = 4
    else:
        raise ValueError(f"nPix={nPix} not supported for healpix (use 12, 48, or 192)")

    # Initialize grid
    InitGrid(Planet.Lateral)

    return Planet


def test_uniform_clathrate():
    """Test uniform clathrate initialization."""
    print("Testing uniform clathrate initialization...")

    Planet = PlanetStruct('TestPlanet')
    Planet.Bulk.R_m = 1560e3
    Planet.zb_km = 100.0
    Planet.Bulk.Tb_K = 273.0
    Planet.Bulk.Tsurf_K = 100.0
    Planet.Bulk.clathMaxDepth_m = 10e3  # 10 km clathrate layer

    # Set up minimal lateral grid
    Planet = _SetupMockLateralGrid(Planet, nPix=12)

    # Initialize uniform clathrate
    Planet = InitClathrateLateral(Planet)

    # Expected clathrate fraction
    expected_fClath = 10e3 / (100e3)  # 0.1

    assert Planet.Lateral.fClath is not None, "fClath not initialized"
    assert len(Planet.Lateral.fClath) == 12, "fClath wrong size"
    assert np.allclose(Planet.Lateral.fClath, expected_fClath, atol=1e-6), \
        f"Uniform fClath incorrect: got {Planet.Lateral.fClath[0]:.3f}, expected {expected_fClath:.3f}"
    assert np.min(Planet.Lateral.fClath) >= 0.0, "fClath below 0"
    assert np.max(Planet.Lateral.fClath) <= 1.0, "fClath above 1"

    print(f"  ✓ Uniform clathrate: fClath = {np.mean(Planet.Lateral.fClath):.3f}")


def test_clathrate_from_sh():
    """Test clathrate initialization from spherical harmonics."""
    print("Testing clathrate from SH coefficients...")

    Planet = PlanetStruct('TestPlanet')
    Planet.Bulk.R_m = 1560e3
    Planet.zb_km = 100.0
    Planet.Bulk.Tb_K = 273.0
    Planet.Bulk.Tsurf_K = 100.0

    # Set up lateral grid
    Planet = _SetupMockLateralGrid(Planet, nPix=48)

    # Create SH coefficients for clathrate fraction
    # Y20 pattern (pole-equator variation): 0.3 at poles, 0.1 at equator
    pMax = 2
    Cpq = np.zeros((pMax + 1, pMax + 1))
    Spq = np.zeros((pMax + 1, pMax + 1))
    Cpq[0, 0] = 0.2  # Uniform baseline
    Cpq[2, 0] = 0.1  # Y20 adds variation

    Planet.Lateral.fClath_Cpq = Cpq
    Planet.Lateral.fClath_Spq = Spq
    Planet.Lateral.fClath_pMax = pMax

    # Initialize clathrate from SH
    Planet = InitClathrateLateral(Planet)

    assert Planet.Lateral.fClath is not None, "fClath not initialized from SH"
    assert len(Planet.Lateral.fClath) == 48, "fClath wrong size"

    # Check range
    fClath_min = np.min(Planet.Lateral.fClath)
    fClath_max = np.max(Planet.Lateral.fClath)
    fClath_mean = np.mean(Planet.Lateral.fClath)

    assert fClath_min >= 0.0, f"fClath below 0: {fClath_min}"
    assert fClath_max <= 1.0, f"fClath above 1: {fClath_max}"
    assert 0.05 < fClath_mean < 0.35, f"Mean fClath unexpected: {fClath_mean:.3f}"
    assert fClath_max > fClath_min, "No variation in SH clathrate field"

    print(f"  ✓ SH clathrate: mean={fClath_mean:.3f}, range=[{fClath_min:.3f}, {fClath_max:.3f}]")


def test_clathrate_clamping():
    """Test clamping of clathrate fraction to [0, 1]."""
    print("Testing clathrate clamping...")

    Planet = PlanetStruct('TestPlanet')
    Planet.Bulk.R_m = 1560e3
    Planet.zb_km = 50.0  # 50 km ice
    Planet.Bulk.Tb_K = 273.0
    Planet.Bulk.Tsurf_K = 100.0
    Planet.Bulk.clathMaxDepth_m = 60e3  # 60 km > ice thickness, should clamp to 1.0

    # Set up lateral grid
    Planet = _SetupMockLateralGrid(Planet, nPix=12)

    # Initialize clathrate (should clamp to 1.0)
    Planet = InitClathrateLateral(Planet)

    assert np.max(Planet.Lateral.fClath) <= 1.0, "fClath not clamped at upper bound"
    assert np.allclose(Planet.Lateral.fClath, 1.0), "fClath not clamped to 1.0 when exceeding ice thickness"

    print(f"  ✓ Clamping works: fClath = {np.mean(Planet.Lateral.fClath):.3f} (clamped to 1.0)")


def test_effective_conductivity():
    """Test effective thermal conductivity calculation."""
    print("Testing effective thermal conductivity...")

    T_mean = 200.0  # K
    k_ice = 651.0 / T_mean  # ~3.255 W/m/K

    # Test pure ice (fClath = 0)
    k_eff_ice = ComputeEffectiveConductivity(0.0, T_mean)
    assert np.isclose(k_eff_ice, k_ice, rtol=1e-6), \
        f"Pure ice conductivity wrong: got {k_eff_ice:.3f}, expected {k_ice:.3f}"

    # Test pure clathrate (fClath = 1)
    k_eff_clath = ComputeEffectiveConductivity(1.0, T_mean)
    assert np.isclose(k_eff_clath, K_CLATH_WmK, rtol=1e-6), \
        f"Pure clathrate conductivity wrong: got {k_eff_clath:.3f}, expected {K_CLATH_WmK:.3f}"

    # Test 50-50 mixture
    k_eff_mix = ComputeEffectiveConductivity(0.5, T_mean)
    k_expected = 0.5 * K_CLATH_WmK + 0.5 * k_ice
    assert np.isclose(k_eff_mix, k_expected, rtol=1e-6), \
        f"Mixed conductivity wrong: got {k_eff_mix:.3f}, expected {k_expected:.3f}"

    # Test temperature dependence
    T_cold = 100.0  # K
    k_eff_cold = ComputeEffectiveConductivity(0.5, T_cold)
    assert k_eff_cold > k_eff_mix, \
        "Colder ice should have higher conductivity"

    # Test array input
    fClath_array = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    k_eff_array = ComputeEffectiveConductivity(fClath_array, T_mean)
    assert len(k_eff_array) == len(fClath_array), "Array output wrong size"
    assert k_eff_array[0] > k_eff_array[-1], \
        "Conductivity should decrease with increasing clathrate fraction"

    print(f"  ✓ Effective conductivity: k_ice={k_ice:.2f} W/m/K, k_clath={K_CLATH_WmK:.2f} W/m/K")
    print(f"    k_eff(50% clath) = {k_eff_mix:.2f} W/m/K")


def test_ice_thickness_from_clathrate():
    """Test self-consistent ice thickness estimation."""
    print("Testing ice thickness estimation from clathrate...")

    Planet = PlanetStruct('TestPlanet')
    Planet.Bulk.R_m = 1560e3
    Planet.zb_km = 100.0  # Reference ice thickness
    Planet.Bulk.Tb_K = 273.0
    Planet.Bulk.Tsurf_K = 100.0

    # Set up lateral grid
    Planet = _SetupMockLateralGrid(Planet, nPix=48)

    # Create clathrate pattern: more at poles, less at equator
    pMax = 2
    Cpq = np.zeros((pMax + 1, pMax + 1))
    Spq = np.zeros((pMax + 1, pMax + 1))
    Cpq[0, 0] = 0.2  # Baseline
    Cpq[2, 0] = 0.15  # Pole-equator variation

    Planet.Lateral.fClath_Cpq = Cpq
    Planet.Lateral.fClath_Spq = Spq
    Planet.Lateral.fClath_pMax = pMax

    Planet = InitClathrateLateral(Planet)

    # Estimate ice thickness
    dIce_m = EstimateIceThicknessFromClathrate(Planet)

    # Check outputs
    assert dIce_m is not None, "Ice thickness not computed"
    assert len(dIce_m) == 48, "Ice thickness wrong size"
    assert np.all(dIce_m > 0), "Negative ice thickness"
    assert Planet.Lateral.kThermEff_WmK is not None, "Effective conductivity not stored"

    # Physical consistency: higher clathrate → lower k_eff → thinner ice
    idx_high_clath = np.argmax(Planet.Lateral.fClath)
    idx_low_clath = np.argmin(Planet.Lateral.fClath)

    assert Planet.Lateral.kThermEff_WmK[idx_high_clath] < Planet.Lateral.kThermEff_WmK[idx_low_clath], \
        "Higher clathrate should give lower conductivity"
    assert dIce_m[idx_high_clath] < dIce_m[idx_low_clath], \
        "Higher clathrate should give thinner ice"

    dIce_mean_km = np.mean(dIce_m) / 1e3
    dIce_min_km = np.min(dIce_m) / 1e3
    dIce_max_km = np.max(dIce_m) / 1e3

    # Sanity check: should be order of reference thickness
    assert 50 < dIce_mean_km < 150, f"Mean ice thickness unrealistic: {dIce_mean_km:.1f} km"

    print(f"  ✓ Ice thickness: mean={dIce_mean_km:.1f} km, range=[{dIce_min_km:.1f}, {dIce_max_km:.1f}] km")
    print(f"    Higher clathrate → thinner ice confirmed")


def test_ice_thickness_with_specified_heat_flux():
    """Test ice thickness estimation with specified basal heat flux."""
    print("Testing ice thickness with specified heat flux...")

    Planet = PlanetStruct('TestPlanet')
    Planet.Bulk.R_m = 1560e3
    Planet.zb_km = 100.0
    Planet.Bulk.Tb_K = 273.0
    Planet.Bulk.Tsurf_K = 100.0

    # Set up lateral grid with uniform clathrate
    Planet = _SetupMockLateralGrid(Planet, nPix=12)
    Planet.Lateral.fClath = np.full(12, 0.3)  # 30% clathrate

    # Specify heat flux
    q_base = 20e-3  # 20 mW/m^2

    dIce_m = EstimateIceThicknessFromClathrate(Planet, q_base_Wm2=q_base)

    # Analytical check: d = k * dT / q
    T_mean = (Planet.Bulk.Tb_K + Planet.Bulk.Tsurf_K) / 2
    k_eff = ComputeEffectiveConductivity(0.3, T_mean)
    dT = Planet.Bulk.Tb_K - Planet.Bulk.Tsurf_K
    d_expected = k_eff * dT / q_base

    assert np.allclose(dIce_m, d_expected, rtol=1e-6), \
        f"Ice thickness calculation wrong: got {dIce_m[0]/1e3:.1f} km, expected {d_expected/1e3:.1f} km"

    print(f"  ✓ Specified heat flux: q={q_base*1e3:.1f} mW/m², d_ice={dIce_m[0]/1e3:.1f} km")


def test_clathrate_no_default():
    """Test clathrate initialization when no default is set."""
    print("Testing clathrate with no default...")

    Planet = PlanetStruct('TestPlanet')
    Planet.Bulk.R_m = 1560e3
    Planet.zb_km = 100.0
    Planet.Bulk.Tb_K = 273.0
    Planet.Bulk.Tsurf_K = 100.0
    # No clathMaxDepth_m set

    Planet = _SetupMockLateralGrid(Planet, nPix=12)
    Planet = InitClathrateLateral(Planet)

    assert Planet.Lateral.fClath is not None, "fClath not initialized"
    assert np.allclose(Planet.Lateral.fClath, 0.0), \
        "Default fClath should be 0 when no clathMaxDepth_m"

    print(f"  ✓ No default: fClath = {np.mean(Planet.Lateral.fClath):.3f}")


def print_scientific_summary():
    """Print scientific understanding summary."""
    print("\n" + "=" * 70)
    print("SCIENTIFIC UNDERSTANDING: Clathrate Lateral Variation")
    print("=" * 70)
    print("""
  Physical mechanism:
    - Clathrate hydrates (gas cages in H₂O) have lower thermal conductivity
      than ice Ih due to trapped gas molecules scattering phonons
    - k_clath ≈ 0.5 W/m/K (constant, weak T-dependence)
    - k_ice ≈ 651/T K W/m/K (strong T-dependence, ~3.3 W/m/K at 200K)
    - Volume-weighted mixing: k_eff = f_clath × k_clath + (1 - f_clath) × k_ice

  Geographic variation:
    - Clathrate fraction f_clath(θ, φ) from spherical harmonics
    - Typical patterns: more clathrate at poles (colder surface)
    - Range: 0 (pure ice) to 1 (pure clathrate), typically 0-0.5

  Self-consistent ice thickness:
    - Conductive heat balance: q = k_eff × (T_melt - T_surf) / d_ice
    - Solve for d_ice: d_ice = k_eff × ΔT / q
    - Higher f_clath → lower k_eff → thinner ice (for fixed q)
    - Lower k_eff → steeper temperature gradient

  Why self-consistency matters:
    - Clathrate formation alters thermal structure
    - Must account for conductivity change in thickness estimate
    - Affects ocean-ice energy balance
    - Important for habitability (ice thickness controls ocean access)

  Typical values (Europa-like):
    - f_clath: 0-0.3 (0-30% clathrate in upper ice shell)
    - Ice thickness: 10-30 km (varies with f_clath and q)
    - Basal heat flux: 10-50 mW/m²
    - k_eff: 1.5-3.0 W/m/K (decreases with f_clath)

  Integration with 3D model:
    - Clathrate pattern set via SH coefficients (like ice thickness, heating)
    - Each column gets unique f_clath value
    - Affects thermal profile calculation in each column
    - Coupled to tidal heating (clathrate has different Q than ice)
""")
    print("=" * 70)


if __name__ == '__main__':
    print("=" * 70)
    print("Phase 6: Clathrate Lateral Variation Tests")
    print("=" * 70)

    try:
        test_uniform_clathrate()
        test_clathrate_from_sh()
        test_clathrate_clamping()
        test_effective_conductivity()
        test_ice_thickness_from_clathrate()
        test_ice_thickness_with_specified_heat_flux()
        test_clathrate_no_default()

        print_scientific_summary()

        print("\n" + "=" * 70)
        print("✓ All Phase 6 tests passed!")
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
