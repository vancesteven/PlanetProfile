#!/usr/bin/env python3
"""
Test script for PlanetProfileApp figure generation functionality.
Run this to verify that the interactive figure system works correctly.
"""

import os
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent))

def test_imports():
    """Test that all required modules can be imported."""
    print("Testing imports...")

    try:
        import streamlit as st
        print("✓ streamlit")
    except ImportError as e:
        print(f"✗ streamlit: {e}")
        return False

    try:
        import pandas as pd
        print("✓ pandas")
    except ImportError as e:
        print(f"✗ pandas: {e}")
        return False

    try:
        import plotly.graph_objects as go
        print("✓ plotly")
        plotly_available = True
    except ImportError as e:
        print(f"✗ plotly (optional): {e}")
        plotly_available = False

    try:
        from Utilities.figure_generator import (
            load_profile_data,
            create_hydrosphere_plot,
            create_seismic_plot,
            PLOTLY_AVAILABLE
        )
        print("✓ figure_generator module")
        print(f"  Plotly available in module: {PLOTLY_AVAILABLE}")
    except ImportError as e:
        print(f"✗ figure_generator module: {e}")
        return False

    return True


def test_data_loading():
    """Test loading sample data if available."""
    print("\nTesting data loading...")

    from Utilities.figure_generator import load_profile_data

    # Look for any profile data file in parent directory
    parent_dir = Path(__file__).parent.parent

    # Check common body directories
    for body in ['Europa', 'Ganymede', 'Callisto', 'Titan', 'Enceladus']:
        body_dir = parent_dir / body
        if body_dir.exists():
            txt_files = list(body_dir.glob('*.txt'))
            txt_files = [f for f in txt_files if 'liquidOcean' not in f.name and 'Core' not in f.name]

            if txt_files:
                print(f"  Found {len(txt_files)} profile files in {body}/")
                test_file = txt_files[0]

                try:
                    df = load_profile_data(str(test_file))
                    print(f"  ✓ Loaded {test_file.name}")
                    print(f"    Shape: {df.shape}")
                    print(f"    Columns: {len(df.columns)}")
                    return True, df, body
                except Exception as e:
                    print(f"  ✗ Failed to load {test_file.name}: {e}")
                    return False, None, None

    print("  No profile data files found in common body directories")
    return None, None, None


def test_figure_generation(df, body_name):
    """Test generating figures from loaded data."""
    print("\nTesting figure generation...")

    from Utilities.figure_generator import (
        create_hydrosphere_plot,
        create_seismic_plot,
        create_custom_plot,
        PLOTLY_AVAILABLE
    )

    if not PLOTLY_AVAILABLE:
        print("  ⚠ Plotly not available, skipping figure generation tests")
        return False

    try:
        # Test hydrosphere plot
        fig = create_hydrosphere_plot(df, body_name)
        if fig:
            print(f"  ✓ Generated hydrosphere plot for {body_name}")
        else:
            print(f"  ✗ Hydrosphere plot returned None")
            return False

        # Test seismic plot
        fig = create_seismic_plot(df, body_name)
        if fig:
            print(f"  ✓ Generated seismic plot for {body_name}")
        else:
            print(f"  ✗ Seismic plot returned None")
            return False

        # Test custom plot
        if len(df.columns) >= 2:
            fig = create_custom_plot(
                df,
                df.columns[0],
                [df.columns[1]],
                title="Test Custom Plot"
            )
            if fig:
                print(f"  ✓ Generated custom plot")
            else:
                print(f"  ✗ Custom plot returned None")
                return False

        return True

    except Exception as e:
        print(f"  ✗ Figure generation failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_figure_saving(df, body_name):
    """Test saving figures to file."""
    print("\nTesting figure saving...")

    from Utilities.figure_generator import (
        create_hydrosphere_plot,
        save_plotly_figure,
        PLOTLY_AVAILABLE
    )

    if not PLOTLY_AVAILABLE:
        print("  ⚠ Plotly not available, skipping figure saving tests")
        return False

    import tempfile

    try:
        fig = create_hydrosphere_plot(df, body_name)

        # Test HTML save
        with tempfile.NamedTemporaryFile(suffix='.html', delete=False) as tmp:
            success = save_plotly_figure(fig, tmp.name, 'html')
            if success and os.path.exists(tmp.name):
                size = os.path.getsize(tmp.name)
                print(f"  ✓ Saved HTML ({size/1024:.1f} KB)")
                os.unlink(tmp.name)
            else:
                print(f"  ✗ Failed to save HTML")
                return False

        print("  ✓ Figure saving works")
        return True

    except Exception as e:
        print(f"  ✗ Figure saving failed: {e}")
        return False


def main():
    """Run all tests."""
    print("="*60)
    print("PlanetProfileApp Figure Generation Test Suite")
    print("="*60)

    # Test imports
    if not test_imports():
        print("\n❌ Import tests failed")
        return 1

    print("\n✅ All imports successful")

    # Test data loading
    success, df, body = test_data_loading()

    if success is None:
        print("\n⚠️  No test data available (run a simulation first)")
        print("✅ Module tests passed")
        return 0
    elif not success:
        print("\n❌ Data loading failed")
        return 1

    print("\n✅ Data loading successful")

    # Test figure generation
    if test_figure_generation(df, body):
        print("\n✅ Figure generation successful")
    else:
        print("\n❌ Figure generation failed")
        return 1

    # Test figure saving
    if test_figure_saving(df, body):
        print("\n✅ Figure saving successful")
    else:
        print("\n❌ Figure saving failed")
        return 1

    print("\n" + "="*60)
    print("✅ ALL TESTS PASSED")
    print("="*60)
    return 0


if __name__ == '__main__':
    sys.exit(main())
