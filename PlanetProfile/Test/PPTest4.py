"""
PPTest4
Europa-like, pure water model with whole-shell clathrate ice
For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test4')

Planet.PfreezeUpper_MPa = 150

""" Bulk planetary settings """
Planet.Bulk.R_m = 1561.0e3
Planet.Bulk.M_kg = 4.7991e22
Planet.Bulk.Tsurf_K = 110
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.346
Planet.Bulk.Cuncertainty = 0.005
Planet.Bulk.Tb_K = 297.8
Planet.Do.CLATHRATE = True
Planet.Steps.nClath = 30  # Gets ignored
Planet.Bulk.clathType = 'whole'
Planet.Bulk.clathMaxThick_m = 10e3  # Gets ignored

""" Layer step settings """
Planet.Steps.nIceI = 200
Planet.Steps.nRefRho = 30
Planet.Steps.nSilMax = 300
Planet.Steps.nCore = 10
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'PureH2O'
Planet.Ocean.wOcean_ppt = 0.0
Planet.Ocean.deltaP = 1.0
Planet.Ocean.PHydroMax_MPa = 250.0

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 5.33e-13
Planet.Sil.Htidal_Wm3 = 1e-11
# Rock porosity
Planet.Do.POROUS_ROCK = False
# Mantle equation of state model
Planet.Sil.mantleEOS = 'CV3hy1wt_678_1.tab'
Planet.Sil.rhoSilWithCore_kgm3 = 3539.0

""" Core assumptions """
Planet.Do.Fe_CORE = True
Planet.Core.rhoFe_kgm3 = 8000.0
Planet.Core.rhoFeS_kgm3 = 5150.0
Planet.Core.rhoPoFeFCC = 5455.0
Planet.Core.QScore = 1e4
Planet.Core.coreEOS = 'Fe-S_3D_EOS.mat'
Planet.Core.wFe_ppt = 850

Planet.Core.xFeSmeteoritic = 0.0405
Planet.Core.xFeS = 0.55
Planet.Core.xFeCore = 0.0279
Planet.Core.xH2O = 0.0035

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0

""" Magnetic induction """
Planet.Bulk.J2 = 435.5e-6
Planet.Bulk.C22 = 131.0e-6
Planet.Magnetic.ionosBounds_m = 100e3
Planet.Magnetic.sigmaIonosPedersen_Sm = 30/100e3
