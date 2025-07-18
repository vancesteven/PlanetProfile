"""
PPTest28
Europa-like, test using constant thermal conductivty for ice layers and ocean layers, specifying activation energy for diffusion of ice phases Ih-VI, and using a different viscosity for rock and corelayers
For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test28')

Planet.PfreezeUpper_MPa = 150

""" Bulk planetary settings """
Planet.Bulk.R_m = 1561.0e3
Planet.Bulk.M_kg = 4.7991e22
Planet.Bulk.Tsurf_K = 110
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.346
Planet.Bulk.Cuncertainty = 0.005
Planet.Bulk.Tb_K = 268.4

""" Layer step settings """
Planet.Steps.nIceI = 50
Planet.Steps.nSilMax = 50
Planet.Steps.nCore = 10
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = "CustomSolutionSeawater = NH4+: 0.5, Cl-: 0.5657647, Na+: 0.4860597, Mg+2: 0.0547421, Ca+2: 0.0106568, K+: 0.0105797, SO4-2: 0.0292643"
Planet.Ocean.wOcean_ppt = 0
Planet.Ocean.deltaP = 1.0
Planet.Ocean.PHydroMax_MPa = 250.0
Planet.Ocean.kThermIce_WmK = {phase: 2 for phase in ['Ih', 'II', 'III', 'V', 'VI', 'Clath']} # New setting
Planet.Ocean.kThermWater_WmK = 0.6 # New setting
Planet.Ocean.Eact_kJmol = [np.nan, 50, 50, 50, 50, 50, 50] # New setting

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 5.33e-12
Planet.Sil.Htidal_Wm3 = 1e-10
# Rock porosity
Planet.Do.POROUS_ROCK = False
# Mantle equation of state model
Planet.Sil.mantleEOS = 'CV3hy1wt_678_1.tab'
Planet.Sil.etaRock_Pas = [1e10, 1e5] # New setting

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
Planet.Core.etaFeSolid_Pas = 1e20 # New setting
Planet.Core.etaFeLiquid_Pas = 1e15 # New setting

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0
