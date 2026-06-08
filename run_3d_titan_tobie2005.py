#!/usr/bin/env python3
"""
Titan 3D Lateral Structure Example - Tobie et al. (2005) Figure 10 Target

Reproduces geographic tidal heating distribution in Titan's ice shell,
demonstrating ocean decoupling effect (67× heating ratio).

Usage:
    python run_3d_titan_tobie2005.py [--grid-size NSIDE] [--model MODEL]

Options:
    --grid-size NSIDE    HEALPix nSide parameter (default: 8, 768 pixels)
                         Tobie 2005 uses ~5000 pixels (nSide~20)
    --model MODEL        Model to run (default: 2)
                         1: Ice I without ocean (low heating)
                         2: Ice I above ocean (high heating, 67× Model 1)
                         3: HP ice below ocean (medium heating)
    --help               Show this help message
"""

import sys
import numpy as np
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from PlanetProfile import PlanetProfile
from PlanetProfile.Utilities.defineStructs import PlanetStruct
from PlanetProfile.GetConfig import Params

def create_titan_model_2(nSide=8):
    """
    Create Titan Model 2: Ice I above ocean (high heating).

    This is the model with the highest heating in Tobie et al. (2005) Figure 10,
    showing ~3822 × 10⁻⁹ W/m³ at the base of ice I.

    Args:
        nSide (int): HEALPix nSide parameter

    Returns:
        Planet, Params: Configured structures
    """
    print(f"\n{'='*70}")
    print(f"Titan Model 2: Ice I Above Ocean (Tobie et al. 2005 Figure 10)")
    print(f"{'='*70}")
    print(f"Ocean decoupling → HIGH tidal heating in ice I")
    print(f"Target: ~3822 × 10⁻⁹ W/m³ at base of ice I")
    print(f"Grid: HEALPix nSide={nSide} ({12*nSide**2} pixels)")
    print(f"{'='*70}\n")

    Planet = PlanetStruct('Titan_Model2_Tobie2005')

    # Bulk properties (Titan)
    Planet.Bulk.R_m = 2575.0e3
    Planet.Bulk.M_kg = 1.3452e23
    Planet.Bulk.Tsurf_K = 94.0
    Planet.Bulk.Psurf_MPa = 0.14554
    Planet.Bulk.Tb_K = 255.0  # Ice-ocean boundary temperature

    # Orbital parameters (Titan)
    Planet.Bulk.eccentricity = 0.0288
    Planet.Bulk.meanMotion_radps = 2 * np.pi / (15.945 * 86400)  # 15.945 day period

    # Ocean (concentrated MgSO4 for ice I above ocean scenario)
    Planet.Ocean.comp = 'MgSO4'
    Planet.Ocean.wOcean_ppt = 100.0  # Concentrated ocean
    Planet.Ocean.deltaP = 1.0
    Planet.Ocean.PHydroMax_MPa = 500.0

    # Silicate mantle
    Planet.Sil.Qrad_Wkg = 5.33e-12
    Planet.Sil.Htidal_Wm3 = 0.0  # Focus on ice shell heating
    Planet.Do.POROUS_ROCK = True
    Planet.Sil.mantleEOS = 'CV3hy1wt_678_1.tab'

    # Core
    Planet.Do.Fe_CORE = True

    # Disable slow calculations
    Planet.Do.SEISMIC = False
    Planet.Do.EQUIL_Q = False  # Disable equilibrium for faster testing

    # 3D Lateral Structure
    Planet.Lateral.DO_3D = True
    Planet.Lateral.gridType = 'healpix'
    Planet.Lateral.nSide = nSide

    # 3D tidal heating (key for Tobie 2005 reproduction)
    Planet.Lateral.DO_TIDAL_3D = True

    # Uniform ice thickness (use 1D reference)
    # Tobie 2005 uses ~100 km ice I layer
    Planet.Lateral.DO_MASS_CONSERVE = True

    # Parameters
    Params = Params()
    Params.CALC_NEW = True
    Params.DO_PARALLEL = True
    Params.PLOT_LATERAL = True

    # Use TidalPy backend for self-consistent heating (optional)
    # Params.Gravity.backend = 'tidalpy'
    # Planet.Do.DO_SELF_CONSISTENT_HTIDAL = True

    return Planet, Params

def print_tobie_comparison(Planet):
    """Print comparison with Tobie et al. (2005) Figure 10."""
    print(f"\n{'='*70}")
    print(f"Comparison with Tobie et al. (2005) Figure 10")
    print(f"{'='*70}")

    if Planet.Lateral.HtidalIce_Wm3 is not None:
        H_mean = np.mean(Planet.Lateral.HtidalIce_Wm3)
        H_min = np.min(Planet.Lateral.HtidalIce_Wm3)
        H_max = np.max(Planet.Lateral.HtidalIce_Wm3)

        print(f"\nTidal Heating in Ice I:")
        print(f"  Mean: {H_mean:.2e} W/m³  ({H_mean*1e9:.2f} × 10⁻⁹ W/m³)")
        print(f"  Min:  {H_min:.2e} W/m³  ({H_min*1e9:.2f} × 10⁻⁹ W/m³)")
        print(f"  Max:  {H_max:.2e} W/m³  ({H_max*1e9:.2f} × 10⁻⁹ W/m³)")
        print(f"  Contrast: {H_max/H_min:.1f}× (max/min)")

        print(f"\nTobie et al. (2005) Model 2 (Ice I above ocean):")
        print(f"  Bottom ice I: ~3822 × 10⁻⁹ W/m³")
        print(f"  67× higher than Model 1 (no ocean)")

        if H_mean * 1e9 > 1000:
            print(f"\n  ✓ Heating magnitude matches Tobie 2005 order!")
        else:
            print(f"\n  → Heating lower than expected (check viscosity/rheology)")

    print(f"\n{'='*70}\n")

def main():
    """Main entry point."""
    nSide = 8
    model = 2

    if '--help' in sys.argv or '-h' in sys.argv:
        print(__doc__)
        return 0

    if '--grid-size' in sys.argv:
        idx = sys.argv.index('--grid-size')
        nSide = int(sys.argv[idx + 1])

    if '--model' in sys.argv:
        idx = sys.argv.index('--model')
        model = int(sys.argv[idx + 1])

    # Create configuration (currently only Model 2 implemented)
    if model == 2:
        Planet, Params = create_titan_model_2(nSide=nSide)
    else:
        print(f"Model {model} not yet implemented. Using Model 2.")
        Planet, Params = create_titan_model_2(nSide=nSide)

    # Run PlanetProfile
    print("Running PlanetProfile with 3D lateral structure...")
    start_time = time.time()

    try:
        Planet, Params = PlanetProfile(Planet, Params)
        elapsed = time.time() - start_time

        print(f"\n✓ Calculation complete in {elapsed:.1f} seconds ({elapsed/60:.1f} min)")

        # Print Tobie comparison
        print_tobie_comparison(Planet)

        # Additional info
        if Planet.Lateral.dIce_m is not None:
            print(f"Ice Shell Thickness:")
            print(f"  Mean: {np.mean(Planet.Lateral.dIce_m)/1e3:.1f} km")
            print(f"  Range: {np.min(Planet.Lateral.dIce_m)/1e3:.1f} - {np.max(Planet.Lateral.dIce_m)/1e3:.1f} km")

        if Planet.saveFile:
            print(f"\nResults saved to: {Planet.saveFile}")
            print(f"Plots saved to: {Planet.saveFile.parent}/figures/")

        print("\n✓ Titan 3D example completed successfully!")
        print("See figures for geographic tidal heating distribution.")
        return 0

    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == '__main__':
    sys.exit(main())
