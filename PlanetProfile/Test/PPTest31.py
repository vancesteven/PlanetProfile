"""
PPTest31
Europa-like, Seawater model with specified hydrosphere pressure at seafloor with constant ice, ocean, and inner properties
For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test31')

Planet.PfreezeUpper_MPa = 150

""" Bulk planetary settings """
Planet.Bulk.R_m = 1561.0e3
Planet.Bulk.M_kg = 4.7991e22
Planet.Bulk.Tsurf_K = 110
Planet.Bulk.Psurf_MPa = 0.0
#Planet.Bulk.Cmeasured = 0.346 This is to be inferred from bulk mass, radius, and hydrosphere pressure at seafloor, so it is explicitly not set
#Planet.Bulk.Cuncertainty = 0.005 This is to be inferred from bulk mass, radius, and hydrosphere pressure at seafloor, so it is explicitly not set
Planet.Bulk.Tb_K = 268.4

""" Layer step settings """
Planet.Steps.nIceI = 50
Planet.Steps.nSilMax = 50
Planet.Steps.nCore = 10
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Do.ConstantProps['Ocean'] = True # (New) Setting constant properties for the ocean layer
Planet.Ocean.ConstantProps.rho_kgm3 = 1200.0 # (New) Setting the density of the ocean layer
Planet.Bulk.PbSet_MPa = 50.0 # (New) Setting the ice I / ocean boundary pressure to be at 50 MPa
Planet.Do.SPECIFY_HYDROSPHERE_SEAFLOOR_PRESSURE = True # Setting the hydrosphere pressure at the seafloor to be explicitly set
Planet.Ocean.PHydroSeafloorSet_MPa = 250.0 # Setting the hydrosphere pressure at the seafloor to be explicitly set

""" Ice assumptions/settings """
Planet.Do.ConstantProps['Ice'] = True # (New) Setting constant properties for the ice layer
Planet.Ocean.IceConstantProps['Ih'].rho_kgm3 = 900.0 # (New) Setting the density of the ice layer

""" Inner assumptions/settings """
Planet.Do.ConstantProps['Inner'] = True # (New) Setting constant properties for the inner layer
Planet.Sil.rhoSilWithCore_kgm3 = 3539.0 # (New) Setting the density of the inner layer

Planet.Do.Fe_CORE = True
Planet.Core.rhoFe_kgm3 = 8000.0
Planet.Core.rhoFeS_kgm3 = 5150.0
Planet.Core.rhoPoFeFCC = 5455.0
Planet.Core.QScore = 1e4
Planet.Core.wFe_ppt = 850

Planet.Core.xFeSmeteoritic = 0.0405
Planet.Core.xFeS = 0.55
Planet.Core.xFeCore = 0.0279
Planet.Core.xH2O = 0.0035

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0
