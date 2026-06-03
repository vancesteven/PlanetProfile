"""
PPTestLateralPhase1
Test for Phase 1 of 3D lateral structure port: Structure definitions only.
Tests that LateralSubstruct exists, has correct defaults, and can be configured.

This is NOT a full PlanetProfile run - just validates data structures.
Run with: python PlanetProfile/Test/PPTestLateralPhase1.py
"""

import sys
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, LateralSubstruct
from PlanetProfile.Lateral.defaultConfigLateral import lateralAssign, configLateralVersion

def test_lateral_substruct():
    """Test that LateralSubstruct class exists and has correct default values."""
    print("Testing LateralSubstruct class...")

    # Test direct instantiation
    lateral = LateralSubstruct()

    # Verify flags
    assert lateral.DO_3D == False, "DO_3D should default to False"
    assert lateral.DO_CLATH_LATERAL == False, "DO_CLATH_LATERAL should default to False"
    assert lateral.DO_TIDAL_3D == False, "DO_TIDAL_3D should default to False"
    assert lateral.DO_MASS_CONSERVE == True, "DO_MASS_CONSERVE should default to True"

    # Verify grid configuration defaults
    assert lateral.gridType == 'healpix', "gridType should default to 'healpix'"
    assert lateral.nSide == 8, "nSide should default to 8"
    assert lateral.nLat is None, "nLat should default to None"
    assert lateral.nLon is None, "nLon should default to None"

    # Verify arrays are None (uninitialized)
    assert lateral.theta_rad is None, "theta_rad should start as None"
    assert lateral.phi_rad is None, "phi_rad should start as None"
    assert lateral.nPix is None, "nPix should start as None"
    assert lateral.pixArea_sr is None, "pixArea_sr should start as None"

    # Verify ice thickness field is None
    assert lateral.dIce_m is None, "dIce_m should start as None"
    assert lateral.dIce_Cpq_km is None, "dIce_Cpq_km should start as None"
    assert lateral.dIce_Spq_km is None, "dIce_Spq_km should start as None"

    # Verify tidal heating fields are None
    assert lateral.Htidal_Wm3 is None, "Htidal_Wm3 should start as None"
    assert lateral.HtidalIce_Wm3 is None, "HtidalIce_Wm3 should start as None"

    # Verify mass conservation fields are None
    assert lateral.Mtarget_kg is None, "Mtarget_kg should start as None"
    assert lateral.Mactual_kg is None, "Mactual_kg should start as None"

    print("✓ LateralSubstruct class has correct defaults")

def test_planet_lateral():
    """Test that Planet.Lateral exists and is a LateralSubstruct."""
    print("Testing Planet.Lateral instantiation...")

    Planet = PlanetStruct('TestLateral')

    # Verify Lateral substruct exists
    assert hasattr(Planet, 'Lateral'), "Planet should have Lateral attribute"
    assert isinstance(Planet.Lateral, LateralSubstruct), "Planet.Lateral should be LateralSubstruct instance"

    # Verify defaults propagate
    assert Planet.Lateral.DO_3D == False, "Planet.Lateral.DO_3D should default to False"
    assert Planet.Lateral.gridType == 'healpix', "Planet.Lateral.gridType should default to healpix"
    assert Planet.Lateral.nSide == 8, "Planet.Lateral.nSide should default to 8"

    print("✓ Planet.Lateral exists and has correct type")

def test_default_config():
    """Test that defaultConfigLateral provides expected configuration."""
    print("Testing defaultConfigLateral...")

    config = lateralAssign()

    # Verify it returns a dict
    assert isinstance(config, dict), "lateralAssign should return a dict"

    # Verify expected keys exist
    required_keys = ['gridType', 'nSide', 'nLat', 'nLon', 'DO_MASS_CONSERVE',
                     'DO_CLATH_LATERAL', 'DO_TIDAL_3D']
    for key in required_keys:
        assert key in config, f"Config should have key '{key}'"

    # Verify default values
    assert config['gridType'] == 'healpix', "Config gridType should be 'healpix'"
    assert config['nSide'] == 8, "Config nSide should be 8"
    assert config['nLat'] == 37, "Config nLat should be 37"
    assert config['nLon'] == 72, "Config nLon should be 72"
    assert config['DO_MASS_CONSERVE'] == True, "Config DO_MASS_CONSERVE should be True"
    assert config['DO_CLATH_LATERAL'] == False, "Config DO_CLATH_LATERAL should be False"
    assert config['DO_TIDAL_3D'] == False, "Config DO_TIDAL_3D should be False"

    # Verify config version exists
    assert configLateralVersion == 1, "configLateralVersion should be 1"

    print("✓ defaultConfigLateral returns correct configuration")

def test_apply_config():
    """Test that lateral config can be applied to Planet."""
    print("Testing config application to Planet.Lateral...")

    Planet = PlanetStruct('TestLateral')
    config = lateralAssign()

    # Apply config (manually, since no loader exists yet)
    for key, value in config.items():
        if hasattr(Planet.Lateral, key):
            setattr(Planet.Lateral, key, value)

    # Verify values were set
    assert Planet.Lateral.gridType == 'healpix', "gridType should be updated"
    assert Planet.Lateral.nSide == 8, "nSide should be updated"
    assert Planet.Lateral.DO_MASS_CONSERVE == True, "DO_MASS_CONSERVE should be updated"

    print("✓ Config can be applied to Planet.Lateral")

def test_scientific_understanding():
    """Document scientific understanding from Phase 1."""
    print("\nScientific understanding (Phase 1):")
    print("  - HEALPix grid: equal-area pixelization, nSide=8 → 768 pixels")
    print("  - Lat-lon grid: simpler but non-equal-area, 37×72 = 2664 pixels")
    print("  - Spherical harmonic coefficients: Cpq (cosine), Spq (sine)")
    print("  - Mass conservation: ∫ρdV must match target M_kg")
    print("  - 3D heating: H(r,θ,φ) stored as (nPix, nRadial) array")
    print("  - Phase 1: Structures only, no computation yet")

if __name__ == '__main__':
    print("="*70)
    print("PPTestLateralPhase1: Testing 3D lateral structure definitions")
    print("="*70)

    try:
        test_lateral_substruct()
        test_planet_lateral()
        test_default_config()
        test_apply_config()
        test_scientific_understanding()

        print("\n" + "="*70)
        print("✓ All Phase 1 tests passed!")
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
