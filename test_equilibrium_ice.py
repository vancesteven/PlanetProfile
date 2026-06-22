#!/usr/bin/env python3
"""
Test script for self-consistent equilibrium ice shell thickness solver.

Reproduces the physics of Tobie et al. (2003, JGR doi:10.1029/2003JE002099)
Figure 12a: equilibrium ice thickness distribution on Europa resulting from
steady-state heat balance between conductive cooling and tidal heating.

Expected Result:
- Equatorial thickening at 0° and 180° longitude (sub/anti-Jovian)
- Mid-latitude thinning
- ~5 km peak-to-peak variation for a ~20-25 km mean shell
- Nearly uniform surface heat flux (~35-40 mW/m² for Europa)

Author: Emma Vellard
Date: 2026-06-17
"""

import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from PlanetProfile.Main import PlanetProfile
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants
from PlanetProfile import GetConfig


def create_europa_equilibrium_test():
    """
    Create Europa configuration for equilibrium ice thickness test.

    Based on Tobie et al. (2003) parameters and Europa observations.
    """
    Planet = PlanetStruct('Europa_Equilibrium_Test')

    # Bulk planetary properties (Anderson et al. 1998)
    Planet.Bulk.R_m = 1560.8e3  # Mean radius
    Planet.Bulk.M_kg = 4.800e22  # Mass
    Planet.Bulk.Tsurf_K = 110    # Surface temperature
    Planet.Bulk.Psurf_MPa = 0.0  # Negligible atmosphere
    Planet.Bulk.Cmeasured = 0.346  # Normalized MoI
    Planet.Bulk.Cuncertainty = 0.005
    Planet.Bulk.Tb_K = 268.0  # Basal temperature (adjust for desired mean thickness)

    # Orbital parameters (Europa in Jupiter system)
    Planet.Bulk.eccentricity = 0.0094  # Orbital eccentricity
    Planet.Bulk.meanMotion_radps = 2 * np.pi / (3.551181 * 86400)  # Period = 3.551181 days

    # Layer step settings
    Planet.Steps.nIceI = 100  # Reduced for faster iteration
    Planet.Steps.nSilMax = 200
    Planet.Steps.nCore = 10
    Planet.Steps.iSilStart = Planet.Steps.nIceI

    # Ocean composition (Seawater)
    Planet.Ocean.comp = 'Seawater'
    Planet.Ocean.wOcean_ppt = Constants.stdSeawater_ppt
    Planet.Ocean.deltaP = 1.0
    Planet.Ocean.deltaT = 0.1
    Planet.Ocean.PHydroMax_MPa = 350.0

    # Silicate mantle (Hussmann & Spohn 2004)
    Planet.Sil.Qrad_Wkg = 5.33e-12  # Radiogenic heating
    Planet.Sil.Htidal_Wm3 = 1e-10    # Tidal heating in silicates
    Planet.Do.POROUS_ROCK = False

    # Core properties (Fe-FeS system)
    Planet.Do.Fe_CORE = True
    Planet.Sil.mantleEOS = 'CM_hydrous_differentiated_Ganymede_Core080Fe020S_excluding_fluid_properties.tab'
    Planet.Sil.rhoSilWithCore_kgm3 = 3539.0

    Planet.Core.rhoFe_kgm3 = 8000.0
    Planet.Core.rhoFeS_kgm3 = 5150.0
    Planet.Core.rhoPoFeFCC = 5455.0
    Planet.Core.QScore = 1e4
    Planet.Core.coreEOS = 'Fe-S_3D_EOS.mat'
    Planet.Core.wFe_ppt = 800
    Planet.Core.xFeS = 0.882
    Planet.Core.xFeCore = 0.0279
    Planet.Core.xH2O = 0.0035

    # Seismic properties
    Planet.Seismic.lowQDiv = 1.0

    # Magnetic induction
    Planet.Bulk.J2 = 435.5e-6
    Planet.Bulk.C22 = 131.0e-6
    Planet.Magnetic.ionosBounds_m = 100e3
    Planet.Magnetic.sigmaIonosPedersen_Sm = 30/100e3

    # Enable necessary calculations
    Planet.Do.SEISMIC = False  # Disable for speed
    Planet.Do.GRAVITY = True   # Need for tidal calculations
    Planet.Do.MAGNETIC = False  # Disable for speed

    # 3D lateral structure configuration
    Planet.Lateral.DO_3D = True
    Planet.Lateral.gridType = 'healpix'
    Planet.Lateral.nSide = 4  # 192 pixels - moderate resolution for testing

    # Enable 3D tidal heating
    Planet.Lateral.DO_TIDAL_3D = True

    # Enable equilibrium ice thickness solver (RECOMMENDED MODE for science)
    # This computes self-consistent thickness from steady-state heat balance:
    #   k * (Tb - Tsurf) / d_ice = q_basal + H_tidal(pixel) * d_ice
    # where tidal heating H_tidal depends on ice viscosity profile, which
    # depends on d_ice. The solver iterates until thickness converges.
    Planet.Lateral.DO_EQUILIBRIUM_ICE = True
    Planet.Lateral.equilibriumTol_m = 100.0  # 100 m convergence tolerance
    Planet.Lateral.equilibriumMaxIter = 10  # Max iterations
    Planet.Lateral.kThermIce_WmK = 2.3  # Ice thermal conductivity (Tobie et al. 2003, Table 2)

    # Optional: Override basal heat flux (W/m²) from silicate mantle
    # If not set, computed from Planet.Sil.Qrad_Wkg and Planet.Sil.Htidal_Wm3
    # Planet.Lateral.qBasal_Wm2 = 0.007  # 7 mW/m² for Europa

    # Enable mass conservation after equilibrium is found
    Planet.Lateral.DO_MASS_CONSERVE = True

    # Do NOT set dIce_Cpq_km when using equilibrium mode
    # The equilibrium solver will compute thickness from physics
    # Initial guess is uniform from reference 1D model (Planet.zb_km)

    return Planet


def main():
    """Main entry point."""
    print("="*70)
    print("Europa Equilibrium Ice Thickness Test")
    print("Tobie et al. (2003) Physics")
    print("="*70)

    # Create Europa configuration
    Planet = create_europa_equilibrium_test()

    # Configure parameters
    Params = GetConfig.Params
    Params.CALC_NEW = True
    Params.DO_PARALLEL = True
    Params.PLOT_LATERAL = True  # Generate lateral plots

    print("\nRunning equilibrium ice thickness calculation...")
    print(f"Grid: nSide={Planet.Lateral.nSide} ({12*Planet.Lateral.nSide**2} pixels)")
    print(f"Tolerance: {Planet.Lateral.equilibriumTol_m:.0f} m")
    print(f"Max iterations: {Planet.Lateral.equilibriumMaxIter}")
    print()

    try:
        Planet, Params = PlanetProfile(Planet, Params)

        print("\n" + "="*70)
        print("EQUILIBRIUM ICE THICKNESS RESULTS")
        print("="*70)

        if Planet.Lateral.dIce_m is not None:
            dIce_km = Planet.Lateral.dIce_m / 1e3
            print(f"\nIce Thickness:")
            print(f"  Mean: {np.mean(dIce_km):.2f} km")
            print(f"  Range: [{np.min(dIce_km):.2f}, {np.max(dIce_km):.2f}] km")
            print(f"  Peak-to-peak: {np.max(dIce_km) - np.min(dIce_km):.2f} km")
            print(f"  Std dev: {np.std(dIce_km):.2f} km")

        if Planet.Lateral.qSurf_Wm2 is not None:
            qSurf_mWm2 = Planet.Lateral.qSurf_Wm2 * 1e3
            print(f"\nSurface Heat Flux:")
            print(f"  Mean: {np.mean(qSurf_mWm2):.2f} mW/m²")
            print(f"  Range: [{np.min(qSurf_mWm2):.2f}, {np.max(qSurf_mWm2):.2f}] mW/m²")
            print(f"  Std dev: {np.std(qSurf_mWm2):.2f} mW/m²")
            print(f"  Variation: {(np.max(qSurf_mWm2)-np.min(qSurf_mWm2))/np.mean(qSurf_mWm2)*100:.1f}%")

        if Planet.Lateral.HtidalIce_Wm3 is not None:
            H_tidal = Planet.Lateral.HtidalIce_Wm3
            print(f"\nTidal Heating:")
            print(f"  Mean: {np.mean(H_tidal):.2e} W/m³")
            print(f"  Range: [{np.min(H_tidal):.2e}, {np.max(H_tidal):.2e}] W/m³")

        if Planet.Lateral.equilibriumIterations is not None:
            print(f"\nConvergence:")
            print(f"  Iterations: {Planet.Lateral.equilibriumIterations}")
            print(f"  Final residual: {Planet.Lateral.equilibriumResidual_m:.1f} m")

        if Planet.Lateral.massResidual_frac is not None:
            print(f"  Mass conservation residual: {Planet.Lateral.massResidual_frac*100:.4f}%")

        print(f"\nOutputs:")
        print(f"  Figures: {Planet.saveFile.replace('.txt', '')}/figures/")
        print(f"  Data: {Planet.saveFile.replace('.txt', '')}/figures/*_lateral3D.npz")

        print("\n" + "="*70)
        print("TEST COMPLETE")
        print("="*70)
        print("\nExpected pattern (Tobie et al. 2003, Fig 12a):")
        print("  - Thicker ice at sub/anti-Jovian points (0°, 180° longitude)")
        print("  - Thinner ice at mid-latitudes")
        print("  - ~5 km peak-to-peak variation")
        print("  - Nearly uniform surface heat flux (~35-40 mW/m²)")

        return 0

    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())
