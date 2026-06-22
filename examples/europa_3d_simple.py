#!/usr/bin/env python3
"""
Simple example: Europa 3D lateral structure with equilibrium ice thickness

Demonstrates the recommended workflow for 3D calculations:
1. Set up baseline Europa parameters
2. Enable equilibrium ice thickness mode (DO_EQUILIBRIUM_ICE=True)
3. Run PlanetProfile
4. Inspect results

This is the cleanest way to use 3D capabilities for scientific studies.

Author: Emma Vellard
Date: 2026-06-18
"""

import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from PlanetProfile.Main import PlanetProfile
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants
from PlanetProfile import GetConfig


def main():
    """Run Europa 3D with equilibrium ice thickness."""

    print("="*70)
    print("Europa 3D Lateral Structure - Equilibrium Ice Thickness Mode")
    print("="*70)

    # 1. Create baseline Europa configuration
    Planet = PlanetStruct('Europa_3D_Example')

    # Bulk properties (Anderson et al. 1998)
    Planet.Bulk.R_m = 1560.8e3
    Planet.Bulk.M_kg = 4.800e22
    Planet.Bulk.Tsurf_K = 110
    Planet.Bulk.Tb_K = 268.0
    Planet.Bulk.Cmeasured = 0.346
    Planet.Bulk.Cuncertainty = 0.005

    # Orbital parameters (REQUIRED for equilibrium mode)
    Planet.Bulk.eccentricity = 0.0094
    Planet.Bulk.meanMotion_radps = 2 * np.pi / (3.551181 * 86400)

    # Ocean composition
    Planet.Ocean.comp = 'Seawater'
    Planet.Ocean.wOcean_ppt = Constants.stdSeawater_ppt

    # Silicate properties
    Planet.Sil.Qrad_Wkg = 5.33e-12
    Planet.Sil.Htidal_Wm3 = 1e-10

    # Enable necessary calculations
    Planet.Do.GRAVITY = True
    Planet.Do.SEISMIC = False  # Disable for speed
    Planet.Do.MAGNETIC = False

    # 2. Configure 3D with equilibrium ice thickness mode
    print("\nConfiguring 3D lateral structure...")
    Planet.Lateral.DO_3D = True
    Planet.Lateral.gridType = 'healpix'
    Planet.Lateral.nSide = 4  # 192 pixels (moderate resolution)

    # Enable 3D tidal heating (required for equilibrium mode)
    Planet.Lateral.DO_TIDAL_3D = True

    # Enable equilibrium ice thickness solver (RECOMMENDED for science)
    Planet.Lateral.DO_EQUILIBRIUM_ICE = True
    Planet.Lateral.equilibriumTol_m = 100.0     # 100 m convergence tolerance
    Planet.Lateral.equilibriumMaxIter = 10      # Max iterations
    Planet.Lateral.kThermIce_WmK = 2.3          # Ice thermal conductivity

    # Optional: Override basal heat flux (otherwise computed from Sil properties)
    # Planet.Lateral.qBasal_Wm2 = 0.007  # 7 mW/m² for Europa

    # Enable mass conservation
    Planet.Lateral.DO_MASS_CONSERVE = True

    print(f"  Grid: {Planet.Lateral.gridType} with nSide={Planet.Lateral.nSide} "
          f"({12*Planet.Lateral.nSide**2} pixels)")
    print(f"  Ice thickness mode: EQUILIBRIUM (physics-based)")
    print(f"  Tidal heating: 3D (Maxwell rheology)")
    print(f"  Convergence tolerance: {Planet.Lateral.equilibriumTol_m:.0f} m")

    # 3. Run PlanetProfile
    print("\nRunning PlanetProfile...")
    Params = GetConfig.Params
    Params.CALC_NEW = True
    Params.DO_PARALLEL = True
    Params.PLOT_LATERAL = True  # Generate lateral plots

    Planet, Params = PlanetProfile(Planet, Params)

    # 4. Inspect results
    print("\n" + "="*70)
    print("RESULTS")
    print("="*70)

    print(f"\nReference 1D Model:")
    print(f"  Ice thickness: {Planet.zb_km:.2f} km")
    print(f"  Ocean thickness: {Planet.D_km:.2f} km")

    print(f"\n3D Equilibrium Ice Thickness:")
    dIce_km = Planet.Lateral.dIce_m / 1e3
    print(f"  Mean: {np.mean(dIce_km):.2f} km")
    print(f"  Range: [{np.min(dIce_km):.2f}, {np.max(dIce_km):.2f}] km")
    print(f"  Peak-to-peak variation: {np.max(dIce_km) - np.min(dIce_km):.2f} km")
    print(f"  Std deviation: {np.std(dIce_km):.2f} km")

    if Planet.Lateral.qSurf_Wm2 is not None:
        qSurf_mWm2 = Planet.Lateral.qSurf_Wm2 * 1e3
        print(f"\nSurface Heat Flux:")
        print(f"  Mean: {np.mean(qSurf_mWm2):.2f} mW/m²")
        print(f"  Range: [{np.min(qSurf_mWm2):.2f}, {np.max(qSurf_mWm2):.2f}] mW/m²")
        print(f"  Variation: {(np.max(qSurf_mWm2)-np.min(qSurf_mWm2))/np.mean(qSurf_mWm2)*100:.1f}%")

    if Planet.Lateral.HtidalIce_Wm3 is not None:
        H_tidal = Planet.Lateral.HtidalIce_Wm3
        print(f"\nTidal Heating (Ice Shell):")
        print(f"  Mean: {np.mean(H_tidal):.2e} W/m³")
        print(f"  Range: [{np.min(H_tidal):.2e}, {np.max(H_tidal):.2e}] W/m³")
        print(f"  Spatial variation: {(np.max(H_tidal)/np.min(H_tidal)):.2f}×")

    print(f"\nConvergence:")
    print(f"  Iterations: {Planet.Lateral.equilibriumIterations}")
    print(f"  Final residual: {Planet.Lateral.equilibriumResidual_m:.1f} m")

    if Planet.Lateral.massResidual_frac is not None:
        print(f"  Mass conservation residual: {Planet.Lateral.massResidual_frac*100:.4f}%")

    print(f"\nOutputs:")
    output_dir = Planet.saveFile.replace('.txt', '')
    print(f"  Figures: {output_dir}/figures/")
    print(f"  Data: {output_dir}/figures/*_lateral3D.npz")

    print("\n" + "="*70)
    print("SUCCESS")
    print("="*70)
    print("\nExpected equilibrium pattern (Tobie et al. 2003):")
    print("  - Thicker ice at sub/anti-Jovian points (0°, 180° longitude)")
    print("  - Thinner ice at mid-latitudes and poles")
    print("  - ~5 km peak-to-peak variation for Europa")
    print("  - Nearly uniform surface heat flux")

    return 0


if __name__ == '__main__':
    sys.exit(main())
