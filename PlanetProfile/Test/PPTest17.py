"""
PPTest17
Ganymede-like, MgSO4 variable ppt ocean comp
For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test17')

""" Bulk planetary settings """
Planet.Bulk.R_m = 2634.1e3
Planet.Bulk.M_kg = 1.4819e23
Planet.Bulk.Tsurf_K = 110
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.3115
Planet.Bulk.Cuncertainty = 0.0028
Planet.Bulk.Tb_K = 258.0

""" Layer step settings """
Planet.Steps.nIceI = 50
Planet.Steps.nSilMax = 50
Planet.Steps.nCore = 10
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'MgSO4'
Planet.Ocean.wOcean_ppt = 0
Planet.Ocean.deltaP = 5.0
Planet.Ocean.PHydroMax_MPa = 1600.0
Planet.Ocean.THydroMax_K = 380.0
Planet.Do.VARIABLE_COMP_OCEAN = True
Planet.Ocean.Pstratified_MPa = np.arange(70, 280, 30)
Planet.Ocean.wStratified_ppt = np.linspace(10, 100, np.size(Planet.Ocean.Pstratified_MPa)+1)

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 5.33e-14
Planet.Sil.Htidal_Wm3 = 1e-18
# Rock porosity
Planet.Do.POROUS_ROCK = False
# Mantle equation of state model
Planet.Sil.mantleEOS = 'CM_hydrous_differentiated_Ganymede_Core080Fe020S_excluding_fluid_properties.tab'
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
