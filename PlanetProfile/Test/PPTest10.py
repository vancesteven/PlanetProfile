"""
PPTest10
Enceladus-like, 1/3 Seawater model including porosity, with no iron core
For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test10')

Planet.PfreezeLower_MPa = 0.1
Planet.PfreezeUpper_MPa = 20

""" Bulk planetary settings """
Planet.Bulk.R_m = 252.1e3
Planet.Bulk.M_kg = 1.08022e20
Planet.Bulk.Tsurf_K = 75
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.335
Planet.Bulk.Cuncertainty = 0.001
Planet.Bulk.Tb_K = 272.353

""" Layer step settings """
Planet.Steps.nIceI = 50
Planet.Steps.nRefRho = 30
Planet.Steps.nSilMax = 120
Planet.Steps.nPoros = 15
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'Seawater'
Planet.Ocean.wOcean_ppt = 1/3 * Constants.stdSeawater_ppt
Planet.Ocean.deltaP = 0.1
Planet.Ocean.PHydroMax_MPa = 25.0

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 5.33e-12
Planet.Sil.Htidal_Wm3 = 3.38e-8

# Rock porosity
Planet.Do.POROUS_ROCK = True
Planet.Sil.porosType = 'Han2014'
Planet.Sil.HtidalMin_Wm3 = 1e-9  # Only needed for non-Han 2014 porosTypes
Planet.Sil.phiRockMax_frac = 0.32
# Mantle equation of state model
Planet.Sil.mantleEOS = 'Comet_67P_CG_v7_excluding_fluid_properties.tab'
Planet.Sil.rhoSilWithCore_kgm3 = 2700.0

""" Core assumptions """
Planet.Do.Fe_CORE = False

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0

""" Magnetic induction """
Planet.Magnetic.peaks_Hz = np.array([4.946e-5, 2.473e-5, 3.259e-6])
Planet.Magnetic.fOrb_radps = 2*np.pi/3.55/86400
Planet.Magnetic.ionosBounds_m = 100e3
Planet.Magnetic.sigmaIonosPedersen_Sm = 30/100e3
