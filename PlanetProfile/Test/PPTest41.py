"""
PPTest41
Andrade no-ocean Titan: MCMC exploration of tidal dissipation parameters.

Configuration for the Andrade rheology MCMC exploration on a no-ocean
Titan model. The interior structure (density, shear modulus, layer
thicknesses) is computed self-consistently by PlanetProfile and held
fixed during the MCMC; only viscosities and Andrade parameters are
varied per sample.

MCMC parameter space (from Petricca et al. 2025, Extended Data Table 2):
  alpha (Andrade exponent):       [0.2, 0.4]
  log10(zeta) (Andrade compliance): [-2, 2]
  log10(eta_Ih) (Ice Ih viscosity): [12, 16]  Pa s
  log10(eta_HP) (HP ice viscosity): [10, 18]  Pa s
  log10(eta_sil) (Silicate visc.):  [18, 22]  Pa s

Observational constraints (Petricca et al. 2025):
  Re(k2) = 0.608 +/- 0.048
  |Im(k2)| = 0.135 +/- 0.035

Uses TidalPy backend for tidal calculations. The MCMC runner script
(Test41_mcmc_andrade_no_ocean.py) calls TidalPy directly with cached
structural arrays for ~100x speedup over full PlanetProfile reruns.

Future refinement:
  - Durham et al. (2001) temperature/pressure-dependent HP ice viscosities
    with separate eta for phases III, V, VI
  - Arakawa et al. (1994) partial-melt viscosity reduction

For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test41')

# Titan analog
Planet.parent = 'Saturn'

""" Bulk planetary settings """
Planet.Bulk.R_m = 2574.73e3
Planet.Bulk.M_kg = 1.3452e23
Planet.Bulk.Tsurf_K = 94
Planet.Bulk.Psurf_MPa = 0.15
Planet.Bulk.Cmeasured = 0.343
Planet.Bulk.Cuncertainty = 0.001
Planet.Do.NONHYDROSTATIC = False
Planet.Bulk.Tb_K = Constants.TtripleIh_III_L_K - 5
Planet.Bulk.eccentricity = 0.0288
Planet.Bulk.meanMotion_radps = 4.56e-6

""" No ocean, but compute inner HP ices """
Planet.Do.HYDROSPHERE_THICKNESS = True
Planet.Bulk.Dhsphere_m = 450e3
Planet.Do.BOTTOM_ICEIII = False
Planet.Do.BOTTOM_ICEV = False
Planet.Do.NO_ICE_CONVECTION = False
Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES = True

""" Clathrate cap """
Planet.Do.CLATHRATE = True
Planet.Bulk.clathType = 'top'
Planet.Bulk.clathMaxThick_m = 10e3
Planet.Steps.nClath = 30

""" Enable Kalousova convection and Arrhenius viscosity """
Planet.Do.KALOUSOVA_CONVECTION = True
Planet.Do.ARRHENIUS_VISCOSITY = True
Planet.Do.DO_SELF_CONSISTENT_HTIDAL = True

""" Melt-bearing ice viscosity (defaults; MCMC overrides per-phase) """
Planet.Ocean.etaMeltKalousova_Pas = {1: 1e12, 3: 1e13, 5: 1e13, 6: 1e13}

""" Ice tidal heating — converged value from PPTest38 """
Planet.Ocean.HtidalIce_Wm3 = 1.15e-7

""" Andrade rheology parameters (defaults; MCMC overrides alpha and zeta) """
Planet.Gravity.andradExponent = 0.2
Planet.Gravity.andrade_zeta = {
    'Ih': 0.005, 'III': 0.05, 'V': 0.05, 'VI': 0.05,
    'Sil': 0.005, 'Fe': 1.0, 'Clath': 1.0, 'default': 1.0
}

""" Layer step settings """
Planet.Steps.nIceI = 200
Planet.Steps.nSilMax = 200
Planet.Steps.nCore = 10
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'PureH2O'
Planet.Ocean.wOcean_ppt = 0
Planet.Ocean.deltaP = 8.0
Planet.Ocean.deltaT = 0.5
Planet.Ocean.phaseType = 'lookup'
Planet.Ocean.PHydroMax_MPa = 1800.0
Planet.Ocean.THydroMax_K = 350.0

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 1.5e-12
Planet.Sil.Htidal_Wm3 = 1e-10
Planet.Sil.etaRock_Pas = [1e19, 1e19]
Planet.Do.POROUS_ROCK = False
Planet.Do.CONSTANT_INNER_DENSITY = True
Planet.Sil.porosType = 'Han2014'
Planet.Sil.HtidalMin_Wm3 = 1e-10
Planet.Sil.phiRockMax_frac = 0.90
Planet.Sil.Pclosure_MPa = 4000
Planet.Sil.mantleEOS = 'Comet_67P_CG_v7_excluding_fluid_properties.tab'
Planet.Sil.rhoSilWithCore_kgm3 = 3539.0

""" Core assumptions """
Planet.Do.Fe_CORE = False
Planet.Core.rhoFe_kgm3 = 8000.0
Planet.Core.rhoFeS_kgm3 = 5150.0
Planet.Core.rhoPoFeFCC = 5455.0
Planet.Core.QScore = 1e4
Planet.Core.coreEOS = 'Fe-S_3D_EOS.mat'
Planet.Core.wFe_ppt = 700
Planet.Core.xFeSmeteoritic = 0.0405
Planet.Core.xFeS = 0.55
Planet.Core.xFeCore = 0.0279
Planet.Core.xH2O = 0.0035

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0

""" Magnetic induction """
Planet.Bulk.J2 = 33.089e-6
Planet.Bulk.C22 = 10.385e-6
Planet.Bulk.C21 = 0.513e-6
Planet.Bulk.S22 = -0.064e-6
Planet.Bulk.S21 = 0.612e-6
Planet.Magnetic.ionosBounds_m = [100e3, 250e3]
Planet.Magnetic.sigmaIonosPedersen_Sm = [1e-16, 80/150e3]
