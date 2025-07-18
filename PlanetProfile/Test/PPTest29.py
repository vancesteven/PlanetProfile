"""
PPTest29
Europa-like, test non self consistent modeling by specifying the following features:
Planet.Do.NON_SELF_CONSISTENT = True
Planet.Bulk.Tsurf_K: Surface temperature
Planet.Gravity.andradExponent: Andrade exponent (for andrade rheology to calculate tidal love numbers)

Planet.dzIceI_km: Ice Ih shell thickness
Planet.Ocean.rhoCondMean_kgm3['Ih']: Ice Conductive Mean Density
Planet.Ocean.rhoConvMean_kgm3['Ih']: Ice Convective Mean Density
Planet.Ocean.GScondMean_GPa: Ice conducting shear modulus
Planet.Ocean.GSconvMean_GPa: Ice convecting shear modulus
Planet.Ocean.kThermIce_WmK['Ih']: Ice Thermal Conductivity
Planet.etaMelt_Pas: Melting point ice viscosity
Planet.Ocean.Eact_kJmol: Creep activation energy

Planet.D_km: Ocean layer thickness
Planet.Ocean.rhoMean_kgm3: Mean density of ocean
Planet.Ocean.kThermWater_WmK: Ocean Thermal Conductivity

Planet.Sil.Rmean_m: Silicate Radius
Planet.Sil.rhoMean_kgm3: Silicate Mean Density
Planet.Sil.kTherm_WmK: Silicate Thermal Conductivity
Planet.Sil.etaRock_Pas: Silicate viscosity
Planet.Sil.GSmean_GPa: Silicate shear modulus

Planet.Core.Rmean_m: Core radius
Planet.Core.rhoMean_kgm3: Core mean density
Planet.Core.kTherm_WmK: Core thermal conductivity
"""

import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test29')

""" General settings """
Planet.Do.NON_SELF_CONSISTENT = True

""" Bulk planetary settings """
Planet.Bulk.R_m = 1561.0e3
Planet.Bulk.M_kg = 4.7991e22
Planet.Bulk.Tsurf_K = 110
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.346
Planet.Bulk.Cuncertainty = 0.005

""" Layer step settings """
Planet.Steps.nIceI = 3
Planet.Steps.nHydro = 1 + Planet.Steps.nIceI
Planet.Steps.nSil = 1
Planet.Steps.nCore = 1
Planet.Ocean.deltaP = 1.0 # Sets the EOS lookup step size
Planet.Ocean.deltaT = 1.0
Planet.Bulk.TfreezeUpper_K = 273

""" Ice assumptions/settings"""
Planet.Do.NO_ICE_CONVECTION = False
Planet.dzIceI_km = 30.0
Planet.Ocean.rhoCondMean_kgm3['Ih'] = 910.0
Planet.Ocean.rhoConvMean_kgm3['Ih'] = 1000.0
Planet.Ocean.kThermIce_WmK['Ih'] = 2.0
Planet.Ocean.Eact_kJmol['Ih'] = 50
Planet.Ocean.GScondMean_GPa['Ih'] = 100.0
Planet.Ocean.GSconvMean_GPa['Ih'] = 100.0
Planet.etaMelt_Pas = 1e14

""" Ocean assumptions/settings """
Planet.D_km = 70
Planet.Ocean.rhoMean_kgm3 = 1000.0
Planet.Ocean.kThermWater_WmK = 0.6

""" Core assumptions """
Planet.Do.Fe_CORE = True
Planet.Core.Rmean_m = 5.98e5
Planet.Core.rhoMean_kgm3 = 5647
Planet.Core.kTherm_WmK = 100.0

""" Silicate Mantle """
Planet.Sil.Rmean_m = Planet.Bulk.R_m - Planet.dzIceI_km - Planet.D_km - Planet.Core.Rmean_m
Planet.Sil.rhoMean_kgm3 = 3290
Planet.Sil.kTherm_WmK = 100.0
Planet.Sil.etaRock_Pas = 1e10
Planet.Sil.GSmean_GPa = 50.0