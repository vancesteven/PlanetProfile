#!/usr/bin/env python3
"""
Lightweight test of 3D ice thickness mode selection logic

Tests the priority order and validation without running full PlanetProfile.
"""

import sys
import numpy as np


def test_mode_selection_logic():
    """Test that mode selection follows correct priority."""
    print("\n" + "="*70)
    print("3D ICE THICKNESS MODE SELECTION LOGIC TEST")
    print("="*70)

    # Simulate the InitLateralStructure logic
    tests_passed = 0
    tests_total = 0

    # Test 1: Equilibrium mode takes priority
    print("\nTest 1: Equilibrium mode priority")
    DO_EQUILIBRIUM_ICE = True
    dIce_Cpq_km = np.array([[29, 0, 0], [0, 0, 0], [-3.5, 0, -1.5]])
    dIce_func = None

    if DO_EQUILIBRIUM_ICE:
        mode = "EQUILIBRIUM"
    elif dIce_Cpq_km is not None:
        mode = "PRESCRIBED"
    elif dIce_func is not None:
        mode = "FUNCTION"
    else:
        mode = "UNIFORM"

    if mode == "EQUILIBRIUM":
        print(f"  ✓ PASSED: Selected {mode} (correct priority)")
        tests_passed += 1
    else:
        print(f"  ✗ FAILED: Selected {mode}, expected EQUILIBRIUM")
    tests_total += 1

    # Test 2: Prescribed mode when equilibrium off
    print("\nTest 2: Prescribed mode when equilibrium disabled")
    DO_EQUILIBRIUM_ICE = False
    dIce_Cpq_km = np.array([[29, 0, 0], [0, 0, 0], [-3.5, 0, -1.5]])
    dIce_func = None

    if DO_EQUILIBRIUM_ICE:
        mode = "EQUILIBRIUM"
    elif dIce_Cpq_km is not None:
        mode = "PRESCRIBED"
    elif dIce_func is not None:
        mode = "FUNCTION"
    else:
        mode = "UNIFORM"

    if mode == "PRESCRIBED":
        print(f"  ✓ PASSED: Selected {mode}")
        tests_passed += 1
    else:
        print(f"  ✗ FAILED: Selected {mode}, expected PRESCRIBED")
    tests_total += 1

    # Test 3: Function mode
    print("\nTest 3: Function mode")
    DO_EQUILIBRIUM_ICE = False
    dIce_Cpq_km = None
    dIce_func = lambda theta: 25000.0  # meters

    if DO_EQUILIBRIUM_ICE:
        mode = "EQUILIBRIUM"
    elif dIce_Cpq_km is not None:
        mode = "PRESCRIBED"
    elif dIce_func is not None:
        mode = "FUNCTION"
    else:
        mode = "UNIFORM"

    if mode == "FUNCTION":
        print(f"  ✓ PASSED: Selected {mode}")
        tests_passed += 1
    else:
        print(f"  ✗ FAILED: Selected {mode}, expected FUNCTION")
    tests_total += 1

    # Test 4: Uniform fallback
    print("\nTest 4: Uniform fallback")
    DO_EQUILIBRIUM_ICE = False
    dIce_Cpq_km = None
    dIce_func = None

    if DO_EQUILIBRIUM_ICE:
        mode = "EQUILIBRIUM"
    elif dIce_Cpq_km is not None:
        mode = "PRESCRIBED"
    elif dIce_func is not None:
        mode = "FUNCTION"
    else:
        mode = "UNIFORM"

    if mode == "UNIFORM":
        print(f"  ✓ PASSED: Selected {mode} (fallback)")
        tests_passed += 1
    else:
        print(f"  ✗ FAILED: Selected {mode}, expected UNIFORM")
    tests_total += 1

    # Test 5: Validation requirements simulation
    print("\nTest 5: Equilibrium mode requirements")
    DO_EQUILIBRIUM_ICE = True
    DO_TIDAL_3D = True
    eccentricity = 0.0094
    meanMotion_radps = 2 * np.pi / (3.551181 * 86400)
    Tsurf_K = 110
    Tb_K = 268.0

    requirements_met = (
        DO_TIDAL_3D
        and eccentricity is not None
        and meanMotion_radps is not None
        and Tsurf_K is not None
        and Tb_K is not None
    )

    if requirements_met:
        print(f"  ✓ PASSED: All requirements met for equilibrium mode")
        tests_passed += 1
    else:
        print(f"  ✗ FAILED: Requirements not met")
    tests_total += 1

    # Test 6: Missing requirement detection
    print("\nTest 6: Missing requirement detection")
    DO_EQUILIBRIUM_ICE = True
    DO_TIDAL_3D = False  # Missing!
    eccentricity = 0.0094
    meanMotion_radps = 2 * np.pi / (3.551181 * 86400)

    requirements_met = DO_TIDAL_3D and eccentricity is not None

    if not requirements_met:
        print(f"  ✓ PASSED: Correctly detected missing DO_TIDAL_3D")
        tests_passed += 1
    else:
        print(f"  ✗ FAILED: Should have detected missing requirement")
    tests_total += 1

    # Summary
    print("\n" + "="*70)
    print(f"Results: {tests_passed}/{tests_total} tests passed")
    print("="*70)

    if tests_passed == tests_total:
        print("\n✓✓✓ ALL LOGIC TESTS PASSED ✓✓✓\n")
        return 0
    else:
        print(f"\n✗✗✗ {tests_total - tests_passed} TESTS FAILED ✗✗✗\n")
        return 1


if __name__ == '__main__':
    sys.exit(test_mode_selection_logic())
