"""
PPTest21
Titan-like, explicit specification of approximate hydrosphere thickness, no explicit ice III or V underplate but rather 
checked for in OceanLayers with Planet.Do.NO_OCEAN Planet.Bulk.Tb_K = Constants.triplePointT_K setting
For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test21')

""" Bulk planetary settings """
Planet.Bulk.R_m = 2574.73e3  # Value from mean radius in Archinal et al. (2018): https://doi.org/10.1007/s10569-017-9805-5
Planet.Bulk.M_kg = 1.3452e23  # Value from Jacobson et al. (2006): https://doi.org/10.1086/508812
Planet.Bulk.Tsurf_K = 94
Planet.Bulk.Psurf_MPa = 0.15
# Planet.Bulk.Cmeasured = 0.341  # Value from Durante et al. (2019): https://doi.org/10.1016/j.icarus.2019.03.003
# Planet.Bulk.Cmeasured = 0.3434  # Value from Flavio
Planet.Bulk.Cmeasured = 0.33  # Value from Flavio but accounting for Gao and Stevenson (2013)
Planet.Bulk.Cuncertainty = 0.01  # No uncertainty is reported by Durante et al.
Planet.Do.NONHYDROSTATIC = False
Planet.Bulk.Tb_K = Constants.triplePointT_K - 20 # Set the Tb_K to be slightly below the triple point to ensure we are in the ice Ih to high pressure ice phase transition
# Planet.Do.ICEIh_THICKNESS = True
# Planet.Bulk.zb_approximate_km = 300 # The approximate ice shell thickness desired
Planet.Do.HYDROSPHERE_THICKNESS = True
Planet.Bulk.Dhsphere_m = 450e3
Planet.Do.BOTTOM_ICEIII = False
Planet.Do.BOTTOM_ICEV = False
Planet.Do.NO_ICE_CONVECTION = False
Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES = True

""" Layer step settings """
Planet.Steps.nIceI = 200
Planet.Steps.nSilMax = 200
Planet.Steps.nCore = 10
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
# Planet.Ocean.comp = 'MgSO4'
Planet.Ocean.comp = 'NaCl'
Planet.Ocean.wOcean_ppt = 100
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
Planet.Sil.porosType = 'Han2014'
Planet.Sil.HtidalMin_Wm3 = 1e-10  # Only needed for non-Han 2014 porosTypes
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
Planet.Bulk.J2 = 33.089e-6  # Values from Durante et al. (2019): https://doi.org/10.1016/j.icarus.2019.03.003
Planet.Bulk.C22 = 10.385e-6
Planet.Bulk.C21 = 0.513e-6
Planet.Bulk.S22 = -0.064e-6
Planet.Bulk.S21 = 0.612e-6
Planet.Magnetic.ionosBounds_m = [100e3, 250e3]
Planet.Magnetic.sigmaIonosPedersen_Sm = [1e-16, 80/150e3]