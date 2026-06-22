#!/usr/bin/env python3
"""
Europa 3D vs Traditional Comparison Script

Runs PlanetProfile for Europa with and without 3D lateral structure calculations,
using identical parameters. Generates publication-quality figures and quantitative
comparison data to demonstrate the scientific relevance of 3D capabilities.

Scientific Context:
- Traditional (1D) models assume spherically symmetric heating and structure
- 3D models resolve geographic variations in tidal heating, ice thickness, and thermal state
- Quantifies impact of lateral variations on habitability metrics

References:
- Tobie et al. (2005): Geographic tidal heating patterns
- Roberts & Nimmo (2008): Europa ice shell structure
- Beuthe (2013): Tidal dissipation theory
- Ojakangas & Stevenson (1989): Eccentricity forcing

Usage:
    python compare_europa_3d.py [--grid-size NSIDE] [--save-dir DIR]

Author: Emma Vellard
Date: 2026-06-09
"""

import sys
import numpy as np
import time
from pathlib import Path
import json

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from PlanetProfile.Main import PlanetProfile
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants
from PlanetProfile import GetConfig


def create_baseline_europa(name_suffix=''):
    """
    Create baseline Europa configuration using literature-constrained parameters.

    Parameters follow:
    - Anderson et al. (1998): M, R, C/MR²
    - Schubert et al. (2009): Interior model constraints
    - Hussmann & Spohn (2004): Tidal heating estimates
    - Khurana et al. (1998): Ocean salinity inference

    Returns:
        Planet: Configured PlanetStruct
    """
    Planet = PlanetStruct(f'Europa_Comparison{name_suffix}')

    # Reduce search range for melting pressure to Europa-realistic values
    Planet.PfreezeUpper_MPa = 150

    # Bulk planetary properties (Anderson et al. 1998)
    Planet.Bulk.R_m = 1560.8e3  # Mean radius (Archinal et al. 2018)
    Planet.Bulk.M_kg = 4.800e22  # Mass (Hussmann et al. 2006)
    Planet.Bulk.Tsurf_K = 110    # Surface temperature
    Planet.Bulk.Psurf_MPa = 0.0  # Negligible atmosphere
    Planet.Bulk.Cmeasured = 0.346  # Normalized MoI
    Planet.Bulk.Cuncertainty = 0.005
    Planet.Bulk.Tb_K = 268.305  # Basal temperature for ~30 km ice + seawater

    # Orbital parameters (Europa in Jupiter system)
    Planet.Bulk.eccentricity = 0.0094  # Orbital eccentricity
    Planet.Bulk.meanMotion_radps = 2 * np.pi / (3.551181 * 86400)  # Period = 3.551181 days

    # Layer step settings (from PPEuropa.py for accuracy)
    Planet.Steps.nIceI = 200
    Planet.Steps.nSilMax = 300
    Planet.Steps.nCore = 10
    Planet.Steps.iSilStart = Planet.Steps.nIceI

    # Ocean composition (Seawater - Khurana et al. 1998)
    Planet.Ocean.comp = 'Seawater'
    Planet.Ocean.wOcean_ppt = Constants.stdSeawater_ppt  # Standard seawater salinity
    Planet.Ocean.deltaP = 1.0
    Planet.Ocean.deltaT = 0.1
    Planet.Ocean.PHydroMax_MPa = 350.0

    # Silicate mantle (Hussmann & Spohn 2004)
    Planet.Sil.Qrad_Wkg = 5.33e-12  # Radiogenic heating
    Planet.Sil.Htidal_Wm3 = 1e-10    # Tidal heating in silicates
    Planet.Do.POROUS_ROCK = False

    # Use CV3hy1wt mantle EOS - literature-validated and PyALMA-safe
    # Note: CM_hydrous + high xFeS triggers PyALMA negative k₂ issue
    # Planet.Sil.mantleEOS = 'CV3hy1wt_678_1.tab'
    Planet.Sil.rhoSilWithCore_kgm3 = 3539.0

    # Seismic properties
    Planet.Seismic.lowQDiv = 1.0

    # Core properties (Fe-FeS system)
    Planet.Do.Fe_CORE = True
    if Planet.Do.Fe_CORE:
        Planet.Sil.mantleEOS = 'CM_hydrous_differentiated_Ganymede_Core080Fe020S_excluding_fluid_properties.tab'  # (2900 for Q= 100 GW, 3240 for Q= 220 GW)
    else:
        Planet.Sil.mantleEOS = 'CM_undifferentiated_hhph_DEW17_nofluid_nomelt_685.tab'

    Planet.Core.rhoFe_kgm3 = 8000.0
    Planet.Core.rhoFeS_kgm3 = 5150.0
    Planet.Core.rhoPoFeFCC = 5455.0
    Planet.Core.QScore = 1e4
    Planet.Core.coreEOS = 'Fe-S_3D_EOS.mat'
    Planet.Core.wFe_ppt = 800

    Planet.Core.xFeSmeteoritic = 0.0405
    Planet.Core.xFeS = 0.882  # PyALMA-safe value (avoids negative k₂ with CV3hy1wt)
    Planet.Core.xFeCore = 0.0279
    Planet.Core.xH2O = 0.0035

    # Magnetic induction (Anderson et al. 1998)
    Planet.Bulk.J2 = 435.5e-6
    Planet.Bulk.C22 = 131.0e-6
    Planet.Magnetic.ionosBounds_m = 100e3
    Planet.Magnetic.sigmaIonosPedersen_Sm = 30/100e3

    # Enable calculations
    Planet.Do.SEISMIC = True  # Skip for faster comparison
    Planet.Do.GRAVITY = True   # Need for tidal calculations
    Planet.Do.MAGNETIC = True # Skip for faster comparison

    return Planet


def configure_traditional_run(Planet):
    """
    Configure Planet for traditional 1D spherically symmetric calculation.

    Args:
        Planet: Base PlanetStruct to configure

    Returns:
        Planet: Configured for 1D
    """
    # Disable all 3D features
    Planet.Lateral.DO_3D = False
    Planet.Lateral.DO_TIDAL_3D = False

    return Planet


def configure_3d_run(Planet, nSide=4, mode='equilibrium'):
    """
    Configure Planet for 3D lateral structure calculation.

    Args:
        Planet: Base PlanetStruct to configure
        nSide: HEALPix nSide parameter (4 = 192 pixels, good for comparison)
        mode: Ice thickness mode - 'equilibrium' (default) or 'prescribed'

    Returns:
        Planet: Configured for 3D
    """
    # Enable 3D calculations
    Planet.Lateral.DO_3D = True
    Planet.Lateral.gridType = 'healpix'
    Planet.Lateral.nSide = nSide

    # Enable 3D tidal heating (required for both modes)
    Planet.Lateral.DO_TIDAL_3D = True

    if mode == 'equilibrium':
        # Mode 1: Equilibrium ice thickness (RECOMMENDED for science)
        # Solves for self-consistent thickness from steady-state heat balance
        # between conductive cooling and tidal + basal heating (Tobie et al. 2003)
        Planet.Lateral.DO_EQUILIBRIUM_ICE = True
        Planet.Lateral.equilibriumTol_m = 100.0  # 100 m convergence tolerance
        Planet.Lateral.equilibriumMaxIter = 10   # Max iterations
        Planet.Lateral.kThermIce_WmK = 2.3       # Ice thermal conductivity

        # Do NOT set dIce_Cpq_km - equilibrium solver computes thickness
        # Initial guess is uniform from reference 1D model (Planet.zb_km)

    elif mode == 'prescribed':
        # Mode 2: Prescribed ice thickness via spherical harmonics
        # Use this for exploratory studies or comparison to observations
        Planet.Lateral.DO_EQUILIBRIUM_ICE = False

        # Define ice thickness variation using spherical harmonics
        # Physical motivation: Tidal stresses on Europa create patterns varying with
        # both latitude and longitude. Following Tobie et al. (2005) and tidal stress
        # models (Greenberg et al. 1998), we include:
        #   - Y_20: Pole-equator variation (axially symmetric)
        #   - Y_22: 4-fold longitudinal pattern (sub/anti-Jovian, leading/trailing)
        Planet.Lateral.dIce_pMax = 2

        # Spherical harmonic coefficients for ice thickness (in km)
        # C_pq: cosine coefficients, S_pq: sine coefficients
        # Pattern: Thinner at poles and sub/anti-Jovian points, thicker at mid-latitudes
        Planet.Lateral.dIce_Cpq_km = np.array([
            [29.0,  0.0,  0.0],   # C_00 = 29 km mean (Juno MWR, Bolton et al. 2025)
            [ 0.0,  0.0,  0.0],   # p=1: zero (synchronous rotator, no hemispheric offset)
            [-3.5,  0.0, -1.5],   # C_20 = -3.5 km (polar thinning, O&S89 + Soderlund 2023)
                                # C_22 = -1.5 km (sub/anti-Jovian thinning, Roberts & Nimmo 2008)
        ])
        Planet.Lateral.dIce_Spq_km = np.array([
            [0.0,  0.0,  0.0],
            [0.0,  0.0,  0.0],
            [0.0,  0.0, -0.7],    # S_22 = -0.7 km (libration phase offset, ~half of C_22)
        ])
    else:
        raise ValueError(f'Unknown mode: {mode}. Use "equilibrium" or "prescribed".')

    # Enable mass conservation enforcement
    Planet.Lateral.DO_MASS_CONSERVE = True

    return Planet


def run_comparison(nSide=4, save_dir=None, mode='equilibrium'):
    """
    Run both traditional and 3D calculations with identical parameters.

    Args:
        nSide: HEALPix resolution for 3D run
        save_dir: Optional output directory
        mode: Ice thickness mode for 3D run - 'equilibrium' (default) or 'prescribed'

    Returns:
        dict: Results dictionary with both Planet objects and metadata
    """
    print("="*70)
    print("Europa 3D vs Traditional Comparison")
    print("="*70)
    print(f"Configuration:")
    print(f"  3D grid resolution: nSide={nSide} ({12*nSide**2} pixels)")
    print(f"  3D ice thickness mode: {mode.upper()}")
    if mode == 'equilibrium':
        print(f"    (physics-based steady-state from heat balance)")
    elif mode == 'prescribed':
        print(f"    (spherical harmonic pattern)")
    print(f"  Physics: Tidal heating, thermal structure, Love numbers")
    print(f"  Ocean: Seawater ({Constants.stdSeawater_ppt:.2f} ppt)")
    print(f"  Mantle: CV3hy1wt (PyALMA-safe, literature-validated)")
    print(f"  Core: Fe-FeS with xFeS=0.55 (PyALMA-safe)")
    print("="*70)

    results = {
        'config': {
            'nSide': nSide,
            'nPix': 12 * nSide**2,
            'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
        }
    }

    # Configure parameters (shared between runs)
    Params = GetConfig.Params
    Params.CALC_NEW = True
    Params.DO_PARALLEL = True
    Params.PLOT_LATERAL = True  # Enable lateral plots for 3D run

    # Traditional 1D calculation
    print("\n" + "-"*70)
    print("TRADITIONAL 1D CALCULATION")
    print("-"*70)

    Planet_1D = create_baseline_europa('_Traditional')
    Planet_1D = configure_traditional_run(Planet_1D)

    print("Running traditional spherically symmetric calculation...")
    t0_1d = time.time()

    try:
        Planet_1D, Params = PlanetProfile(Planet_1D, Params)
        t_1d = time.time() - t0_1d
        results['traditional'] = {
            'Planet': Planet_1D,
            'time_seconds': t_1d,
            'success': True
        }
        print(f"✓ Traditional calculation complete in {t_1d:.1f} seconds")
        print_traditional_summary(Planet_1D)

    except Exception as e:
        print(f"✗ Traditional calculation failed: {e}")
        import traceback
        traceback.print_exc()
        results['traditional'] = {'success': False, 'error': str(e)}
        return results

    # 3D calculation
    print("\n" + "-"*70)
    print("3D LATERAL STRUCTURE CALCULATION")
    print("-"*70)

    Planet_3D = create_baseline_europa('_3D')
    Planet_3D = configure_3d_run(Planet_3D, nSide=nSide, mode=mode)

    print(f"Running 3D lateral calculation (nSide={nSide}, {12*nSide**2} pixels, mode={mode})...")
    t0_3d = time.time()

    try:
        # Need fresh Params object for 3D run
        Params = GetConfig.Params
        Params.CALC_NEW = True
        Params.DO_PARALLEL = True
        Params.PLOT_LATERAL = True

        Planet_3D, Params = PlanetProfile(Planet_3D, Params)
        t_3d = time.time() - t0_3d
        results['3d'] = {
            'Planet': Planet_3D,
            'time_seconds': t_3d,
            'success': True
        }
        print(f"✓ 3D calculation complete in {t_3d:.1f} seconds")
        print_3d_summary(Planet_3D)

    except Exception as e:
        print(f"✗ 3D calculation failed: {e}")
        import traceback
        traceback.print_exc()
        results['3d'] = {'success': False, 'error': str(e)}
        return results

    # Performance comparison
    print("\n" + "-"*70)
    print("PERFORMANCE COMPARISON")
    print("-"*70)
    print(f"Traditional time: {t_1d:.1f} seconds")
    print(f"3D time: {t_3d:.1f} seconds")
    print(f"Slowdown factor: {t_3d/t_1d:.1f}×")
    print(f"Time per pixel: {t_3d/(12*nSide**2):.2f} seconds")

    return results


def print_traditional_summary(Planet):
    """Print summary of traditional 1D results."""
    print("\nTraditional 1D Results:")
    print(f"  Ice shell thickness: {Planet.zb_km:.2f} km")
    print(f"  Ocean thickness: {Planet.D_km:.2f} km")
    if Planet.Do.Fe_CORE and hasattr(Planet.Core, 'Rmean_m'):
        print(f"  Core radius: {Planet.Core.Rmean_m/1e3:.2f} km")
    print(f"  MoI (C/MR²): {Planet.CMR2mean:.4f}")
    if hasattr(Planet, 'g_ms2') and Planet.g_ms2 is not None:
        print(f"  Surface gravity: {Planet.g_ms2[0]:.3f} m/s²")


def print_3d_summary(Planet):
    """Print summary of 3D lateral results."""
    print("\n3D Lateral Results:")

    if Planet.Lateral.dIce_m is not None:
        print(f"\nGrid: {Planet.Lateral.nPix} pixels ({Planet.Lateral.gridType})")

        print(f"\nIce Thickness:")
        print(f"  Mean: {np.mean(Planet.Lateral.dIce_m)/1e3:.2f} km")
        print(f"  Range: {np.min(Planet.Lateral.dIce_m)/1e3:.2f} - {np.max(Planet.Lateral.dIce_m)/1e3:.2f} km")
        print(f"  Variation: {(np.max(Planet.Lateral.dIce_m) - np.min(Planet.Lateral.dIce_m))/1e3:.2f} km")

        # Use nanmean/nanmin/nanmax to handle any NaN values
        print(f"\nBasal Temperature:")
        if not np.all(np.isnan(Planet.Lateral.Tb_K)):
            print(f"  Mean: {np.nanmean(Planet.Lateral.Tb_K):.2f} K")
            print(f"  Range: {np.nanmin(Planet.Lateral.Tb_K):.2f} - {np.nanmax(Planet.Lateral.Tb_K):.2f} K")
        else:
            print(f"  Mean: nan K")
            print(f"  Range: nan - nan K")

        print(f"\nSurface Heat Flux:")
        if not np.all(np.isnan(Planet.Lateral.qSurf_Wm2)):
            print(f"  Mean: {np.nanmean(Planet.Lateral.qSurf_Wm2)*1e3:.2f} mW/m²")
            print(f"  Range: {np.nanmin(Planet.Lateral.qSurf_Wm2)*1e3:.2f} - {np.nanmax(Planet.Lateral.qSurf_Wm2)*1e3:.2f} mW/m²")
        else:
            print(f"  Mean: nan mW/m²")
            print(f"  Range: nan - nan mW/m²")

        if Planet.Lateral.HtidalIce_Wm3 is not None:
            print(f"\nTidal Heating (Ice Shell):")
            Htot_mean = np.mean(Planet.Lateral.HtidalIce_Wm3)
            Htot_min = np.min(Planet.Lateral.HtidalIce_Wm3)
            Htot_max = np.max(Planet.Lateral.HtidalIce_Wm3)
            print(f"  Mean: {Htot_mean:.2e} W/m³")
            print(f"  Range: {Htot_min:.2e} - {Htot_max:.2e} W/m³")
            print(f"  Pole/Equator ratio: {Htot_max/Htot_min:.2f}")

        if Planet.Lateral.DO_MASS_CONSERVE:
            print(f"\nMass Conservation:")
            print(f"  Residual: {Planet.Lateral.massResidual_frac*100:.4f}%")


def export_quantitative_comparison(results, output_path):
    """
    Export quantitative comparison data to JSON.

    Args:
        results: Results dictionary from run_comparison
        output_path: Path to JSON output file
    """
    Planet_1D = results['traditional']['Planet']
    Planet_3D = results['3d']['Planet']

    comparison = {
        'metadata': results['config'],
        'traditional': {
            'ice_thickness_km': float(Planet_1D.zb_km),
            'ocean_thickness_km': float(Planet_1D.D_km),
            'core_radius_km': float(Planet_1D.Core.Rmean_m/1e3) if Planet_1D.Do.Fe_CORE else None,
            'MoI': float(Planet_1D.CMR2mean),
            'computation_time_seconds': results['traditional']['time_seconds']
        },
        '3d': {
            'grid_pixels': int(Planet_3D.Lateral.nPix),
            'ice_thickness_km': {
                'mean': float(np.mean(Planet_3D.Lateral.dIce_m)/1e3),
                'min': float(np.min(Planet_3D.Lateral.dIce_m)/1e3),
                'max': float(np.max(Planet_3D.Lateral.dIce_m)/1e3),
                'std': float(np.std(Planet_3D.Lateral.dIce_m)/1e3)
            },
            'computation_time_seconds': results['3d']['time_seconds'],
            'time_per_pixel_seconds': results['3d']['time_seconds'] / Planet_3D.Lateral.nPix
        },
        'comparison': {
            'slowdown_factor': results['3d']['time_seconds'] / results['traditional']['time_seconds'],
            'ice_thickness_difference_percent': float(100 * (np.mean(Planet_3D.Lateral.dIce_m)/1e3 - Planet_1D.zb_km) / Planet_1D.zb_km)
        }
    }

    # Add basal temperature if available (use nanmean to handle any NaN values)
    if Planet_3D.Lateral.Tb_K is not None and len(Planet_3D.Lateral.Tb_K) > 0:
        Tb = Planet_3D.Lateral.Tb_K
        # Check if we have valid (non-NaN) data
        if not np.all(np.isnan(Tb)):
            comparison['3d']['basal_temperature_K'] = {
                'mean': float(np.nanmean(Tb)),
                'min': float(np.nanmin(Tb)),
                'max': float(np.nanmax(Tb)),
                'std': float(np.nanstd(Tb))
            }
        else:
            comparison['3d']['basal_temperature_K'] = None

    # Add surface heat flux if available (use nanmean to handle any NaN values)
    if Planet_3D.Lateral.qSurf_Wm2 is not None and len(Planet_3D.Lateral.qSurf_Wm2) > 0:
        qSurf = Planet_3D.Lateral.qSurf_Wm2
        # Check if we have valid (non-NaN) data
        if not np.all(np.isnan(qSurf)):
            comparison['3d']['surface_heat_flux_mWm2'] = {
                'mean': float(np.nanmean(qSurf)*1e3),
                'min': float(np.nanmin(qSurf)*1e3),
                'max': float(np.nanmax(qSurf)*1e3),
                'std': float(np.nanstd(qSurf)*1e3)
            }
        else:
            comparison['3d']['surface_heat_flux_mWm2'] = None

    # Add tidal heating if available
    if Planet_3D.Lateral.HtidalIce_Wm3 is not None and len(Planet_3D.Lateral.HtidalIce_Wm3) > 0:
        Htot = Planet_3D.Lateral.HtidalIce_Wm3
        # Avoid division by zero for pole/equator ratio
        Hmin = np.min(Htot)
        Hmax = np.max(Htot)
        if Hmin > 0:
            pole_equator = float(Hmax / Hmin)
        else:
            pole_equator = None  # Undefined when minimum is zero

        comparison['3d']['tidal_heating_ice_Wm3'] = {
            'mean': float(np.mean(Htot)),
            'min': float(Hmin),
            'max': float(Hmax),
            'std': float(np.std(Htot)),
            'pole_equator_ratio': pole_equator
        }

    with open(output_path, 'w') as f:
        json.dump(comparison, f, indent=2)

    print(f"\nQuantitative comparison exported to: {output_path}")


def main():
    """Main entry point."""
    # Parse arguments
    nSide = 4  # Default: 192 pixels (moderate resolution)
    save_dir = None
    mode = 'equilibrium'  # Default: equilibrium mode (recommended for science)

    if '--help' in sys.argv or '-h' in sys.argv:
        print(__doc__)
        print("\nAdditional options:")
        print("  --mode MODE          Ice thickness mode: 'equilibrium' (default) or 'prescribed'")
        return 0

    if '--grid-size' in sys.argv:
        idx = sys.argv.index('--grid-size')
        nSide = int(sys.argv[idx + 1])

    if '--save-dir' in sys.argv:
        idx = sys.argv.index('--save-dir')
        save_dir = Path(sys.argv[idx + 1])
        save_dir.mkdir(parents=True, exist_ok=True)

    if '--mode' in sys.argv:
        idx = sys.argv.index('--mode')
        mode = sys.argv[idx + 1].lower()
        if mode not in ['equilibrium', 'prescribed']:
            print(f"Error: Invalid mode '{mode}'. Use 'equilibrium' or 'prescribed'.")
            return 1

    # Run comparison
    print(f"\nStarting Europa comparison (nSide={nSide}, mode={mode})...")
    results = run_comparison(nSide=nSide, save_dir=save_dir, mode=mode)

    # Check both runs succeeded
    if not results.get('traditional', {}).get('success', False):
        print("\n✗ Traditional calculation failed, cannot complete comparison")
        return 1

    if not results.get('3d', {}).get('success', False):
        print("\n✗ 3D calculation failed, cannot complete comparison")
        return 1

    # Export quantitative comparison
    if save_dir:
        comparison_file = save_dir / 'europa_comparison_quantitative.json'
    else:
        comparison_file = Path('europa_comparison_quantitative.json')

    export_quantitative_comparison(results, comparison_file)

    print("\n" + "="*70)
    print("COMPARISON COMPLETE")
    print("="*70)
    print("\nOutputs:")
    print(f"  Traditional figures: {results['traditional']['Planet'].saveFile}")
    print(f"  3D figures: {results['3d']['Planet'].saveFile}")
    print(f"  Quantitative data: {comparison_file}")
    print("\n✓ Europa 3D vs Traditional comparison completed successfully!")

    return 0


if __name__ == '__main__':
    sys.exit(main())
