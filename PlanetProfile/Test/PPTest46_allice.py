"""
PPTest46_allice
All-ice Titan template for the Andrade MCMC (Test46_mcmc_allice.py).

Derived from PPTest41 (no-ocean Titan, Kalousova convection, Arrhenius
viscosity). Key difference: D_hsphere_m and rhoSilWithCore_kgm3 are set at
runtime by the MCMC script, rather than fixed.

Interior structure:
  - Ice Ih shell (surface conductive + convective interior)
  - HP ice (III, V, VI) via NO_OCEAN_EXCEPT_INNER_ICES
  - Single silicate mantle (CONSTANT_INNER_DENSITY, Fe_CORE=False)
  - No liquid ocean

The MCMC script (Test46_mcmc_allice.py) calls PlanetProfile with varied
(D_hsphere_km, rho_sil) to build the structural grid, then samples Andrade
rheology and viscosity parameters against Petricca et al. (2025) constraints.

For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test46_allice')

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
Planet.Bulk.Tb_K = Constants.TtripleIh_III_L_K - 5  # default; overridden by MCMC
Planet.Bulk.eccentricity = 0.0288
Planet.Bulk.meanMotion_radps = 4.56e-6

""" No ocean, but compute inner HP ices """
Planet.Do.HYDROSPHERE_THICKNESS = False
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

""" Ice tidal heating """
Planet.Ocean.HtidalIce_Wm3 = 1.15e-7

""" Andrade rheology parameters (defaults; MCMC overrides alpha and zeta) """
Planet.Gravity.andradExponent = 0.3
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
Planet.Sil.rhoSilWithCore_kgm3 = 3300.0  # default; overridden by MCMC

""" Core assumptions """
Planet.Do.Fe_CORE = False

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0

""" Magnetic induction """
Planet.Bulk.J2 = 33.089e-6
Planet.Bulk.C22 = 10.385e-6
Planet.Magnetic.peaks_Hz = np.array([4.56e-6])
Planet.Magnetic.fOrb_radps = Planet.Bulk.meanMotion_radps
