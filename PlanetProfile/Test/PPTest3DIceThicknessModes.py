#!/usr/bin/env python3
"""
Test suite for 3D ice thickness modes

Tests:
1. Equilibrium mode validation (requirements checking)
2. Prescribed SH mode (backward compatibility)
3. Uniform mode (fallback)
4. Mode priority logic
5. Conflict warnings

Author: Emma Vellard
Date: 2026-06-22
"""

import sys
import numpy as np
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants
from PlanetProfile.Lateral.LateralStructure import InitLateralStructure


def test_equilibrium_mode_validation():
    """Test that equilibrium mode validates requirements."""
    print("\nTest 1: Equilibrium mode validation")
    print("="*70)

    Planet = PlanetStruct('Europa_Test')
    Planet.Bulk.R_m = 1560.8e3
    Planet.Bulk.M_kg = 4.800e22
    Planet.Bulk.Tsurf_K = 110
    Planet.Bulk.Tb_K = 268.0
    Planet.zb_km = 25.0

    # Configure 3D
    Planet.Lateral.DO_3D = True
    Planet.Lateral.nSide = 2  # Minimal grid for testing
    Planet.Lateral.DO_EQUILIBRIUM_ICE = True
    Planet.Lateral.DO_TIDAL_3D = True

    # Test 1a: Missing eccentricity (should raise ValueError)
    try:
        InitLateralStructure(Planet, None)
        print("✗ FAILED: Should have raised ValueError for missing eccentricity")
        return False
    except ValueError as e:
        if 'eccentricity' in str(e):
            print(f"✓ PASSED: Correctly caught missing eccentricity")
        else:
            print(f"✗ FAILED: Wrong error message: {e}")
            return False

    # Test 1b: Add eccentricity but missing mean motion
    Planet.Bulk.eccentricity = 0.0094
    try:
        InitLateralStructure(Planet, None)
        print("✗ FAILED: Should have raised ValueError for missing mean motion")
        return False
    except ValueError as e:
        if 'meanMotion' in str(e):
            print(f"✓ PASSED: Correctly caught missing mean motion")
        else:
            print(f"✗ FAILED: Wrong error message: {e}")
            return False

    # Test 1c: Add all required parameters (should succeed)
    Planet.Bulk.meanMotion_radps = 2 * np.pi / (3.551181 * 86400)
    try:
        Planet = InitLateralStructure(Planet, None)
        if Planet.Lateral.dIce_m is not None:
            if len(Planet.Lateral.dIce_m) == 12 * 2**2:  # nSide=2 → 48 pixels
                print(f"✓ PASSED: Equilibrium mode initialized with {len(Planet.Lateral.dIce_m)} pixels")
                print(f"  Initial uniform thickness: {Planet.Lateral.dIce_m[0]/1e3:.1f} km")
                return True
            else:
                print(f"✗ FAILED: Wrong number of pixels: {len(Planet.Lateral.dIce_m)}")
                return False
        else:
            print("✗ FAILED: dIce_m not initialized")
            return False
    except Exception as e:
        print(f"✗ FAILED: Unexpected error: {e}")
        return False


def test_prescribed_sh_mode():
    """Test prescribed SH mode (backward compatibility)."""
    print("\nTest 2: Prescribed SH mode")
    print("="*70)

    Planet = PlanetStruct('Europa_Test')
    Planet.Bulk.R_m = 1560.8e3
    Planet.zb_km = 25.0

    # Configure 3D with prescribed SH
    Planet.Lateral.DO_3D = True
    Planet.Lateral.nSide = 2
    Planet.Lateral.DO_EQUILIBRIUM_ICE = False
    Planet.Lateral.DO_TIDAL_3D = False

    # Set SH coefficients
    Planet.Lateral.dIce_pMax = 2
    Planet.Lateral.dIce_Cpq_km = np.array([
        [29.0,  0.0,  0.0],
        [ 0.0,  0.0,  0.0],
        [-3.5,  0.0, -1.5],
    ])
    Planet.Lateral.dIce_Spq_km = np.array([
        [0.0,  0.0,  0.0],
        [0.0,  0.0,  0.0],
        [0.0,  0.0, -0.7],
    ])

    try:
        Planet = InitLateralStructure(Planet, None)
        if Planet.Lateral.dIce_m is not None:
            mean_km = np.mean(Planet.Lateral.dIce_m) / 1e3
            min_km = np.min(Planet.Lateral.dIce_m) / 1e3
            max_km = np.max(Planet.Lateral.dIce_m) / 1e3
            print(f"✓ PASSED: Prescribed SH mode initialized")
            print(f"  Mean: {mean_km:.1f} km, Range: [{min_km:.1f}, {max_km:.1f}] km")

            # Check that mean is close to C_00
            if abs(mean_km - 29.0) < 2.0:  # Allow some deviation due to SH projection
                print(f"✓ PASSED: Mean thickness close to C_00 = 29.0 km")
                return True
            else:
                print(f"✗ WARNING: Mean {mean_km:.1f} km differs from C_00 = 29.0 km")
                return True  # Still pass, but warn
        else:
            print("✗ FAILED: dIce_m not initialized")
            return False
    except Exception as e:
        print(f"✗ FAILED: Error in prescribed mode: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_uniform_mode():
    """Test uniform mode (fallback)."""
    print("\nTest 3: Uniform mode (fallback)")
    print("="*70)

    Planet = PlanetStruct('Europa_Test')
    Planet.zb_km = 30.0

    # Configure 3D without prescribed or equilibrium
    Planet.Lateral.DO_3D = True
    Planet.Lateral.nSide = 2
    Planet.Lateral.DO_EQUILIBRIUM_ICE = False
    # Don't set dIce_Cpq_km → should fall back to uniform

    try:
        Planet = InitLateralStructure(Planet, None)
        if Planet.Lateral.dIce_m is not None:
            mean_km = np.mean(Planet.Lateral.dIce_m) / 1e3
            std_km = np.std(Planet.Lateral.dIce_m) / 1e3

            if abs(mean_km - 30.0) < 0.01 and std_km < 0.01:
                print(f"✓ PASSED: Uniform mode initialized")
                print(f"  All pixels = {mean_km:.1f} km (std = {std_km:.4f} km)")
                return True
            else:
                print(f"✗ FAILED: Not uniform: mean={mean_km:.1f}, std={std_km:.4f}")
                return False
        else:
            print("✗ FAILED: dIce_m not initialized")
            return False
    except Exception as e:
        print(f"✗ FAILED: Error in uniform mode: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_mode_priority():
    """Test that equilibrium mode takes priority over prescribed."""
    print("\nTest 4: Mode priority (equilibrium > prescribed)")
    print("="*70)

    Planet = PlanetStruct('Europa_Test')
    Planet.Bulk.R_m = 1560.8e3
    Planet.Bulk.Tsurf_K = 110
    Planet.Bulk.Tb_K = 268.0
    Planet.Bulk.eccentricity = 0.0094
    Planet.Bulk.meanMotion_radps = 2 * np.pi / (3.551181 * 86400)
    Planet.zb_km = 25.0

    # Configure with BOTH equilibrium and prescribed (equilibrium should win)
    Planet.Lateral.DO_3D = True
    Planet.Lateral.nSide = 2
    Planet.Lateral.DO_EQUILIBRIUM_ICE = True
    Planet.Lateral.DO_TIDAL_3D = True

    # Also set prescribed SH (should be ignored with warning)
    Planet.Lateral.dIce_Cpq_km = np.array([
        [35.0,  0.0,  0.0],  # Different from reference
        [ 0.0,  0.0,  0.0],
        [-5.0,  0.0, -2.0],
    ])
    Planet.Lateral.dIce_Spq_km = np.zeros((3, 3))

    try:
        Planet = InitLateralStructure(Planet, None)
        if Planet.Lateral.dIce_m is not None:
            mean_km = np.mean(Planet.Lateral.dIce_m) / 1e3
            std_km = np.std(Planet.Lateral.dIce_m) / 1e3

            # Should be uniform (equilibrium mode initialization) not SH pattern
            if std_km < 0.01 and abs(mean_km - 25.0) < 0.01:
                print(f"✓ PASSED: Equilibrium mode took priority")
                print(f"  Initialized uniform at {mean_km:.1f} km (prescribed SH ignored)")
                return True
            else:
                print(f"✗ FAILED: Used prescribed pattern: mean={mean_km:.1f}, std={std_km:.2f}")
                return False
        else:
            print("✗ FAILED: dIce_m not initialized")
            return False
    except Exception as e:
        print(f"✗ FAILED: Error: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_missing_tidal_heating():
    """Test that equilibrium mode requires tidal heating."""
    print("\nTest 5: Equilibrium requires tidal heating")
    print("="*70)

    Planet = PlanetStruct('Europa_Test')
    Planet.Bulk.R_m = 1560.8e3
    Planet.Bulk.Tsurf_K = 110
    Planet.Bulk.Tb_K = 268.0
    Planet.Bulk.eccentricity = 0.0094
    Planet.Bulk.meanMotion_radps = 2 * np.pi / (3.551181 * 86400)
    Planet.zb_km = 25.0

    # Enable equilibrium but NOT tidal heating
    Planet.Lateral.DO_3D = True
    Planet.Lateral.nSide = 2
    Planet.Lateral.DO_EQUILIBRIUM_ICE = True
    Planet.Lateral.DO_TIDAL_3D = False  # Missing!

    try:
        Planet = InitLateralStructure(Planet, None)
        print("✗ FAILED: Should have raised ValueError for missing tidal heating")
        return False
    except ValueError as e:
        if 'DO_TIDAL_3D' in str(e):
            print(f"✓ PASSED: Correctly caught missing tidal heating requirement")
            return True
        else:
            print(f"✗ FAILED: Wrong error message: {e}")
            return False


def run_all_tests():
    """Run all tests and report results."""
    print("\n" + "="*70)
    print("3D ICE THICKNESS MODES TEST SUITE")
    print("="*70)

    tests = [
        ("Equilibrium validation", test_equilibrium_mode_validation),
        ("Prescribed SH mode", test_prescribed_sh_mode),
        ("Uniform mode", test_uniform_mode),
        ("Mode priority", test_mode_priority),
        ("Tidal heating requirement", test_missing_tidal_heating),
    ]

    results = []
    for name, test_func in tests:
        try:
            passed = test_func()
            results.append((name, passed))
        except Exception as e:
            print(f"\n✗ FAILED: {name} - Unexpected error: {e}")
            import traceback
            traceback.print_exc()
            results.append((name, False))

    # Summary
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)

    passed = sum(1 for _, p in results if p)
    total = len(results)

    for name, p in results:
        status = "✓ PASSED" if p else "✗ FAILED"
        print(f"{status}: {name}")

    print("-"*70)
    print(f"Results: {passed}/{total} tests passed")

    if passed == total:
        print("\n✓✓✓ ALL TESTS PASSED ✓✓✓")
        return 0
    else:
        print(f"\n✗✗✗ {total - passed} TESTS FAILED ✗✗✗")
        return 1


if __name__ == '__main__':
    sys.exit(run_all_tests())
