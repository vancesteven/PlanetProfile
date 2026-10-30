#!/usr/bin/env python3
"""
Test script to verify all utility functions have correct signatures
and can be imported without errors.
"""

import sys
import os

# Add the app directory to path
app_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, app_dir)

def test_imports():
    """Test that all utility modules can be imported"""
    print("Testing imports...")

    try:
        from Utilities import app_helpers
        print("✓ app_helpers imported")
    except Exception as e:
        print(f"✗ app_helpers import failed: {e}")
        return False

    try:
        from Utilities import presets
        print("✓ presets imported")
    except Exception as e:
        print(f"✗ presets import failed: {e}")
        return False

    try:
        from Utilities import help_system
        print("✓ help_system imported")
    except Exception as e:
        print(f"✗ help_system import failed: {e}")
        return False

    try:
        from Utilities import session_manager
        print("✓ session_manager imported")
    except Exception as e:
        print(f"✗ session_manager import failed: {e}")
        return False

    return True


def test_function_signatures():
    """Test that key functions have the expected signatures"""
    print("\nTesting function signatures...")

    from Utilities.app_helpers import (
        show_tips_carousel,
        create_progress_indicator,
        validate_radius,
        validate_mass,
        format_scientific,
        create_metric_card,
        estimate_runtime
    )

    from Utilities.presets import (
        get_presets_for_body,
        get_preset,
        apply_preset,
        get_preset_summary
    )

    from Utilities.help_system import (
        get_help,
        create_help_button
    )

    from Utilities.session_manager import (
        show_session_manager_ui,
        show_recent_sessions
    )

    # Test show_tips_carousel with optional parameter
    import inspect
    sig = inspect.signature(show_tips_carousel)
    params = list(sig.parameters.keys())
    assert 'tips' in params or len(params) == 0, f"show_tips_carousel signature wrong: {params}"
    print("✓ show_tips_carousel signature OK")

    # Test create_progress_indicator with multiple parameter options
    sig = inspect.signature(create_progress_indicator)
    params = list(sig.parameters.keys())
    assert 'current_step' in params or 'steps_complete' in params, f"create_progress_indicator missing parameters: {params}"
    print("✓ create_progress_indicator signature OK")

    # Test get_preset_summary with flexible parameters
    sig = inspect.signature(get_preset_summary)
    params = list(sig.parameters.keys())
    assert 'preset_config' in params or ('body_name' in params and 'preset_name' in params), \
        f"get_preset_summary signature wrong: {params}"
    print("✓ get_preset_summary signature OK")

    # Test validate functions
    sig = inspect.signature(validate_radius)
    assert len(sig.parameters) >= 1, "validate_radius needs at least 1 parameter"
    print("✓ validate_radius signature OK")

    sig = inspect.signature(validate_mass)
    assert len(sig.parameters) >= 1, "validate_mass needs at least 1 parameter"
    print("✓ validate_mass signature OK")

    print("\nAll function signature tests passed!")
    return True


def test_help_system():
    """Test that help system has required data"""
    print("\nTesting help system...")

    from Utilities.help_system import HELP_TEXT, FAQ

    # Check that HELP_TEXT is a dict
    assert isinstance(HELP_TEXT, dict), "HELP_TEXT should be a dict"
    print(f"✓ HELP_TEXT has {len(HELP_TEXT)} entries")

    # Check that FAQ is a dict
    assert isinstance(FAQ, dict), "FAQ should be a dict"
    print(f"✓ FAQ has {len(FAQ)} entries")

    return True


def test_presets():
    """Test that presets are properly structured"""
    print("\nTesting presets...")

    from Utilities.presets import PRESETS

    assert isinstance(PRESETS, dict), "PRESETS should be a dict"
    print(f"✓ PRESETS has {len(PRESETS)} bodies")

    # Check structure of first preset
    for body_name, body_presets in PRESETS.items():
        assert isinstance(body_presets, dict), f"{body_name} presets should be a dict"
        for preset_name, preset_config in body_presets.items():
            assert 'description' in preset_config, f"{body_name}/{preset_name} missing description"
            assert 'settings' in preset_config, f"{body_name}/{preset_name} missing settings"
            print(f"✓ {body_name}/{preset_name} structure OK")
        break  # Just check first body

    return True


def main():
    """Run all tests"""
    print("="*60)
    print("Testing PlanetProfileApp Utilities")
    print("="*60)

    all_passed = True

    # Test imports
    if not test_imports():
        print("\n✗ Import tests FAILED")
        all_passed = False
    else:
        print("\n✓ All imports successful")

    # Test function signatures
    if not test_function_signatures():
        print("\n✗ Function signature tests FAILED")
        all_passed = False

    # Test help system
    if not test_help_system():
        print("\n✗ Help system tests FAILED")
        all_passed = False

    # Test presets
    if not test_presets():
        print("\n✗ Presets tests FAILED")
        all_passed = False

    print("\n" + "="*60)
    if all_passed:
        print("✓ ALL TESTS PASSED")
        print("="*60)
        return 0
    else:
        print("✗ SOME TESTS FAILED")
        print("="*60)
        return 1


if __name__ == "__main__":
    sys.exit(main())
