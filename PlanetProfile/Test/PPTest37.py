"""
PPTest37
Titan no-ocean with Kalousova HP ice convection and clathrate underplate
Same as PPTest35 but with clathType='bottom' (clathrate at base of ice Ih)
For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test37')

""" Bulk planetary settings """
Planet.Bulk.R_m = 2574.73e3  # Archinal et al. (2018)
Planet.Bulk.M_kg = 1.3452e23  # Jacobson et al. (2006)
Planet.Bulk.Tsurf_K = 94
Planet.Bulk.Psurf_MPa = 0.15
Planet.Bulk.Cmeasured = 0.343  # Petricca et al. (2025)
Planet.Bulk.Cuncertainty = 0.001
Planet.Do.NONHYDROSTATIC = False
Planet.Bulk.Tb_K = Constants.triplePointT_K - 5  # Below triple point for no-ocean
Planet.Bulk.eccentricity = 0.0288  # Titan orbital eccentricity
Planet.Bulk.meanMotion_radps = 4.56e-6  # Titan mean motion

""" No ocean, but compute inner HP ices """
Planet.Do.HYDROSPHERE_THICKNESS = True
Planet.Bulk.Dhsphere_m = 450e3
Planet.Do.BOTTOM_ICEIII = False
Planet.Do.BOTTOM_ICEV = False
Planet.Do.NO_ICE_CONVECTION = False
Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES = True

""" Clathrate underplate at base of ice Ih per Kalousova & Sotin (2020) """
Planet.Do.CLATHRATE = True
Planet.Bulk.clathType = 'bottom'
Planet.Bulk.clathMaxThick_m = 10e3
Planet.Bulk.qSurf_Wm2 = 6e-3  # 6 mW/m^2 surface heat flux
Planet.Steps.nClath = 30

""" Enable Kalousova convection and Arrhenius viscosity for HP ices """
Planet.Do.KALOUSOVA_CONVECTION = True
Planet.Do.ARRHENIUS_VISCOSITY = True
Planet.Do.DO_SELF_CONSISTENT_HTIDAL = True

""" Melt-bearing ice viscosity overrides per Petricca et al. (2025).
    Two-phase (partial melt) convection from Kalousova dramatically lowers
    the effective viscosity of all ice phases to ~1e12 Pa·s average. """
Planet.Ocean.etaMeltKalousova_Pas = {1: 1e12, 3: 5e12, 5: 5e12, 6: 5e12}

""" Ice tidal heating from Petricca et al. (2025) mid-range estimate """
Planet.Ocean.HtidalIce_Wm3 = 5e-9

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
Planet.Sil.etaRock_Pas = [1e22, 1e22]  # Fixed viscosity per Petricca et al. (2025)
# Rock porosity
Planet.Do.POROUS_ROCK = False
Planet.Do.CONSTANT_INNER_DENSITY = True
Planet.Sil.porosType = 'Han2014'
Planet.Sil.HtidalMin_Wm3 = 1e-10
Planet.Sil.phiRockMax_frac = 0.90
Planet.Sil.Pclosure_MPa = 4000
# Mantle equation of state model
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
Planet.Bulk.J2 = 33.089e-6  # Durante et al. (2019)
Planet.Bulk.C22 = 10.385e-6
Planet.Bulk.C21 = 0.513e-6
Planet.Bulk.S22 = -0.064e-6
Planet.Bulk.S21 = 0.612e-6
Planet.Magnetic.ionosBounds_m = [100e3, 250e3]
Planet.Magnetic.sigmaIonosPedersen_Sm = [1e-16, 80/150e3]
