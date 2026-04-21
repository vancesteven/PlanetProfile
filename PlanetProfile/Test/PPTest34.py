"""
PPTest34
Titan with Kalousova et al. (2018) HP ice convection and subsurface ocean
Tests HP ice convection with Titan parameters and ice tidal heating
Orbital parameters from Petricca et al. (2025)
For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test34')

""" Bulk planetary settings """
Planet.Bulk.R_m = 2574.73e3  # Archinal et al. (2018)
Planet.Bulk.M_kg = 1.3452e23  # Jacobson et al. (2006)
Planet.Bulk.Tsurf_K = 94
Planet.Bulk.Psurf_MPa = 0.15
Planet.Bulk.Cmeasured = 0.33  # Accounting for Gao and Stevenson (2013)
Planet.Bulk.Cuncertainty = 0.01
Planet.Bulk.Tb_K = 255.0
Planet.Bulk.eccentricity = 0.0288  # Titan orbital eccentricity
Planet.Bulk.meanMotion_radps = 4.56e-6  # Titan mean motion

""" Enable Kalousova convection """
Planet.Do.KALOUSOVA_CONVECTION = True

""" Ice tidal heating from Petricca et al. (2025) estimates """
Planet.Ocean.HtidalIce_Wm3 = 1e-9

""" Layer step settings """
Planet.Steps.nIceI = 100
Planet.Steps.nIceIIILitho = 20
Planet.Steps.nIceVLitho = 20
Planet.Steps.nSilMax = 100
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
# Rock porosity
Planet.Do.POROUS_ROCK = False
Planet.Do.CONSTANT_INNER_DENSITY = True
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
