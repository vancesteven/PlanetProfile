#!/usr/bin/env python3
"""
3D Lateral Structure Example Runner

Demonstrates PlanetProfile 3D lateral structure capability with a simple Europa-like example.

Usage:
    python run_3d_example.py [--grid-size NSIDE] [--no-tidal] [--no-plot]

Options:
    --grid-size NSIDE    HEALPix nSide parameter (default: 2, 48 pixels)
                         Options: 1 (12 pix), 2 (48 pix), 4 (192 pix), 8 (768 pix)
    --no-tidal           Disable 3D tidal heating (faster)
    --no-plot            Skip lateral plotting (data still saved)
    --help               Show this help message
"""

import sys
import numpy as np
import time
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from PlanetProfile.Main import PlanetProfile
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants
from PlanetProfile import GetConfig

def create_3d_example(nSide=2, do_tidal=True, do_plot=True):
    """
    Create a 3D lateral structure example configuration.

    Args:
        nSide (int): HEALPix nSide parameter (nPix = 12*nSide^2)
        do_tidal (bool): Whether to compute 3D tidal heating
        do_plot (bool): Whether to generate lateral plots

    Returns:
        Planet: Configured Planet object
        Params: Configured Params object
    """
    print(f"\n{'='*60}")
    print(f"3D Lateral Structure Example")
    print(f"{'='*60}")
    print(f"Grid: HEALPix nSide={nSide} ({12*nSide**2} pixels)")
    print(f"Tidal heating: {'Enabled' if do_tidal else 'Disabled'}")
    print(f"Plotting: {'Enabled' if do_plot else 'Disabled'}")
    print(f"{'='*60}\n")

    # Create Planet structure
    Planet = PlanetStruct('Europa_3D_Example')

    Planet.PfreezeUpper_MPa = 150

    # Bulk properties (Europa-like)
    Planet.Bulk.R_m = 1561.0e3
    Planet.Bulk.M_kg = 4.7991e22
    Planet.Bulk.Tsurf_K = 110
    Planet.Bulk.Psurf_MPa = 0.0
    Planet.Bulk.Cmeasured = 0.346
    Planet.Bulk.Cuncertainty = 0.005
    Planet.Bulk.Tb_K = 268.4

    # Orbital parameters (Europa)
    Planet.Bulk.eccentricity = 0.0094
    Planet.Bulk.meanMotion_radps = 2 * np.pi / (3.55 * 86400)

    # Layer step settings
    Planet.Steps.nIceI = 50
    Planet.Steps.nSilMax = 50
    Planet.Steps.nCore = 10
    Planet.Steps.iSilStart = Planet.Steps.nIceI

    # Ocean (Seawater)
    Planet.Ocean.comp = 'Seawater'
    Planet.Ocean.wOcean_ppt = Constants.stdSeawater_ppt
    Planet.Ocean.deltaP = 1.0
    Planet.Ocean.PHydroMax_MPa = 250.0

    # Silicate mantle (simple EOS to avoid PyALMA issues)
    Planet.Sil.Qrad_Wkg = 5.33e-12
    Planet.Sil.Htidal_Wm3 = 1e-10
    Planet.Do.POROUS_ROCK = False
    Planet.Sil.mantleEOS = 'CV3hy1wt_678_1.tab'
    Planet.Sil.rhoSilWithCore_kgm3 = 3539.0

    # Seismic properties
    Planet.Seismic.lowQDiv = 1.0

    # Core
    Planet.Do.Fe_CORE = True
    Planet.Core.rhoFe_kgm3 = 8000.0
    Planet.Core.rhoFeS_kgm3 = 5150.0
    Planet.Core.rhoPoFeFCC = 5455.0
    Planet.Core.QScore = 1e4
    Planet.Core.coreEOS = 'Fe-S_3D_EOS.mat'
    Planet.Core.wFe_ppt = 850
    Planet.Core.xFeSmeteoritic = 0.0405
    Planet.Core.xFeS = 0.55  # Moderate value to avoid PyALMA negative k2
    Planet.Core.xFeCore = 0.0279
    Planet.Core.xH2O = 0.0035

    # Disable slow calculations for faster demo
    Planet.Do.SEISMIC = False
    Planet.Do.GRAVITY = True  # Need this for 3D lateral
    Planet.Do.MAGNETIC = False

    # 3D Lateral Structure Configuration
    Planet.Lateral.DO_3D = True
    Planet.Lateral.gridType = 'healpix'
    Planet.Lateral.nSide = nSide

    # 3D tidal heating
    Planet.Lateral.DO_TIDAL_3D = do_tidal

    # Add ice thickness variation (Y_20 pole-equator pattern)
    Planet.Lateral.dIce_pMax = 2
    Planet.Lateral.dIce_Cpq_km = np.array([
        [25.0, 0.0, 0.0],  # C_00 = mean 25 km
        [0.0, 0.0, 0.0],   # p=1 (no degree-1)
        [5.0, 0.0, 0.0]    # C_20 = 5 km variation (thicker at poles)
    ])
    Planet.Lateral.dIce_Spq_km = np.zeros((3, 3))

    # Mass conservation
    Planet.Lateral.DO_MASS_CONSERVE = True

    # Configure parameters
    Params = GetConfig.Params
    Params.CALC_NEW = True
    Params.DO_PARALLEL = True  # Enable parallel processing
    Params.PLOT_LATERAL = do_plot

    return Planet, Params

def print_results(Planet):
    """Print summary of 3D results."""
    print(f"\n{'='*60}")
    print(f"3D Results Summary")
    print(f"{'='*60}")

    if Planet.Lateral.dIce_m is not None:
        print(f"\nGrid:")
        print(f"  Pixels: {Planet.Lateral.nPix}")
        print(f"  Grid type: {Planet.Lateral.gridType}")

        print(f"\nIce Thickness:")
        print(f"  Mean: {np.mean(Planet.Lateral.dIce_m)/1e3:.2f} km")
        print(f"  Min: {np.min(Planet.Lateral.dIce_m)/1e3:.2f} km")
        print(f"  Max: {np.max(Planet.Lateral.dIce_m)/1e3:.2f} km")
        print(f"  Range: {(np.max(Planet.Lateral.dIce_m) - np.min(Planet.Lateral.dIce_m))/1e3:.2f} km")

        print(f"\nBasal Temperature:")
        print(f"  Mean: {np.mean(Planet.Lateral.Tb_K):.2f} K")
        print(f"  Min: {np.min(Planet.Lateral.Tb_K):.2f} K")
        print(f"  Max: {np.max(Planet.Lateral.Tb_K):.2f} K")
        print(f"  Range: {np.max(Planet.Lateral.Tb_K) - np.min(Planet.Lateral.Tb_K):.2f} K")

        print(f"\nSurface Heat Flux:")
        print(f"  Mean: {np.mean(Planet.Lateral.qSurf_Wm2)*1e3:.2f} mW/m²")
        print(f"  Min: {np.min(Planet.Lateral.qSurf_Wm2)*1e3:.2f} mW/m²")
        print(f"  Max: {np.max(Planet.Lateral.qSurf_Wm2)*1e3:.2f} mW/m²")

        if Planet.Lateral.HtidalIce_Wm3 is not None:
            print(f"\nTidal Heating (Ice):")
            print(f"  Mean: {np.mean(Planet.Lateral.HtidalIce_Wm3):.2e} W/m³")
            print(f"  Min: {np.min(Planet.Lateral.HtidalIce_Wm3):.2e} W/m³")
            print(f"  Max: {np.max(Planet.Lateral.HtidalIce_Wm3):.2e} W/m³")

        if Planet.Lateral.DO_MASS_CONSERVE:
            print(f"\nMass Conservation:")
            print(f"  Target mass: {Planet.Lateral.Mtarget_kg:.4e} kg")
            print(f"  Actual mass: {Planet.Lateral.Mactual_kg:.4e} kg")
            print(f"  Residual: {Planet.Lateral.massResidual_frac*100:.4f}%")

    print(f"\n{'='*60}\n")

def main():
    """Main entry point."""
    # Parse command line arguments
    nSide = 2
    do_tidal = True
    do_plot = True

    if '--help' in sys.argv or '-h' in sys.argv:
        print(__doc__)
        return 0

    if '--grid-size' in sys.argv:
        idx = sys.argv.index('--grid-size')
        nSide = int(sys.argv[idx + 1])

    if '--no-tidal' in sys.argv:
        do_tidal = False

    if '--no-plot' in sys.argv:
        do_plot = False

    # Create configuration
    Planet, Params = create_3d_example(nSide=nSide, do_tidal=do_tidal, do_plot=do_plot)

    # Run PlanetProfile
    print("Running PlanetProfile with 3D lateral structure...")
    start_time = time.time()

    try:
        Planet, Params = PlanetProfile(Planet, Params)
        elapsed = time.time() - start_time

        print(f"\n✓ Calculation complete in {elapsed:.1f} seconds")

        # Print results
        print_results(Planet)

        # Save location
        if Planet.saveFile:
            print(f"Results saved to: {Planet.saveFile}")

        if do_plot:
            from pathlib import Path
            saveDir = Path(Planet.saveFile).parent if isinstance(Planet.saveFile, str) else Planet.saveFile.parent
            print(f"Plots saved to: {saveDir}/figures/")

        print("\n✓ 3D lateral structure example completed successfully!")
        return 0

    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == '__main__':
    sys.exit(main())
