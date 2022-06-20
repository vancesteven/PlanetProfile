"""
PPTest5
Io-like, waterless world
For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test5')

""" Bulk planetary settings """
Planet.Do.NO_H2O = True
Planet.Bulk.qSurf_Wm2 = 140e-3
Planet.Bulk.R_m = 1821.3e3
Planet.Bulk.M_kg = 8.9319e22
Planet.Bulk.Tsurf_K = 110
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.37824
Planet.Bulk.Cuncertainty = 0.00022

""" Layer step settings """
Planet.Steps.nSilMax = 300
Planet.Steps.nCore = 10

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 5.33e-12
Planet.Sil.Htidal_Wm3 = 1e-9
# Rock porosity
Planet.Do.POROUS_ROCK = False
Planet.Sil.poreH2Orho_kgm3 = 0.0
# Mantle equation of state model
Planet.Sil.mantleEOS = 'CV3hy1wt_678_1.tab'
Planet.Do.CONSTANT_INNER_DENSITY = True  # Forcing this on because we don't have a good Io-specific composition dialed in yet
Planet.Sil.rhoSilWithCore_kgm3 = 3290.0

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
