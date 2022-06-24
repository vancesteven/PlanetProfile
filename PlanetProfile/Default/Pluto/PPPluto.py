"""
PPPluto
Contains all body-specific parameters and information for PlanetProfile models of this body.
Import as a module and access information assigned to the attributes of the Planet struct.
Note that this file expects to be imported from the directory above.
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Pluto')

# Reduce search range for melting pressure to values realistic for Europa
Planet.PfreezeUpper_MPa = 150

""" Bulk planetary settings """
Planet.Bulk.R_m = 1195.0e3
Planet.Bulk.M_kg = 1.314e22
Planet.Bulk.Tsurf_K = 44
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.31  # Based on rough forward models with fixed layer densities in Hussmann et al. (2006): https://doi.org/10.1016/j.icarus.2006.06.005
Planet.Bulk.Cuncertainty = 0.03
Planet.Bulk.Tb_K = 345.0  # Thickness <100 km is unlikely re: Kamata
Planet.Do.CLATHRATE = True
Planet.Steps.nClath = 30
Planet.Bulk.clathType = 'bottom'
Planet.Bulk.clathMaxThick_m = 10.0e3  # Estimate based on Kamata et al. (2019): https://doi.org/10.1038/s41561-019-0369-8
Planet.Bulk.qSurf_Wm2 = 14.0e-3  # Maximum of ~13 mW/m^2 based on Conrad et al. (2021): https://doi.org/10.1029/2020JE006641
# Highly porous surface ice for top ~5 km re: Kamata

""" Layer step settings """
Planet.Steps.nIceI = 200
Planet.Steps.nSilMax = 120
Planet.Steps.nCore = 10
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'MgSO4'
Planet.Ocean.wOcean_ppt = 50
Planet.Ocean.deltaP = 0.25
Planet.Ocean.deltaT = 0.1
Planet.Ocean.PHydroMax_MPa = 250.0

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 4.5e-12
Planet.Sil.Htidal_Wm3 = 1e-10
# Rock porosity
Planet.Do.POROUS_ROCK = True
Planet.Sil.phiRockMax_frac = 0.35
Planet.Sil.mantleEOS = 'Comet_67P_CG_v7_excluding_fluid_properties.tab'
Planet.Sil.rhoSilWithCore_kgm3 = 2700.0

""" Core assumptions """
Planet.Do.Fe_CORE = False

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0

""" Magnetic induction """
Planet.Bulk.J2 = 435.5e-6  # Copied from Europa model
Planet.Bulk.C22 = 131.0e-6
