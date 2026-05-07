"""
PPTest50
No-ocean Titan with Yao 2014 spherical-shell Ih convection (for Test50 MCMC).

Derived from PPTest46_allice by adding SPHERICAL_CONVECTION=True so that the
Ice Ih shell uses Yao et al. (2014) scaling instead of Deschamps & Sotin
(2001). Everything else — clathrate cap, Kalousova HP convection, Arrhenius
viscosity, no-ocean structure — is unchanged.

Used by Test50_mcmc_andrade_noocean_yao2014.py to test Petricca's hypothesis
that no-ocean + coupled HP ices can supply Im(k2) via HP dissipation.
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test50')

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
Planet.Bulk.Tb_K = Constants.triplePointT_K - 5  # default; overridden by MCMC
Planet.Bulk.eccentricity = 0.0288
Planet.Bulk.meanMotion_radps = 4.56e-6

""" No ocean, but compute inner HP ices """
Planet.Do.HYDROSPHERE_THICKNESS = False
Planet.Do.BOTTOM_ICEIII = False
Planet.Do.BOTTOM_ICEV = False
Planet.Do.NO_ICE_CONVECTION = False
Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES = True

""" Yao 2014 spherical convection for Ice Ih """
Planet.Do.SPHERICAL_CONVECTION = True

""" Clathrate cap """
Planet.Do.CLATHRATE = True
Planet.Bulk.clathType = 'top'
Planet.Bulk.clathMaxThick_m = 10e3
Planet.Steps.nClath = 30

""" Kalousova HP convection + Arrhenius viscosity """
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
