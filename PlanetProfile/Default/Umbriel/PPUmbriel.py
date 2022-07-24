"""
PPUmbriel
Contains all body-specific parameters and information for PlanetProfile models of this body.
Import as a module and access information assigned to the attributes of the Planet struct.
Note that this file expects to be imported from the directory above.
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Umbriel')

Planet.PfreezeLower_MPa = 5
Planet.PfreezeUpper_MPa = 45

""" Bulk planetary settings """
Planet.Bulk.R_m = 584.7e3
Planet.Bulk.M_kg = 1.172e21
Planet.Bulk.Tsurf_K = 60
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.319  # Based on rough forward models with fixed layer densities in Hussmann et al. (2006): https://doi.org/10.1016/j.icarus.2006.06.005
Planet.Bulk.Cuncertainty = 0.030
Planet.Bulk.Tb_K = 270.93

""" Layer step settings """
Planet.Steps.nIceI = 50
Planet.Steps.nSilMax = 80
Planet.Steps.nPoros = 8
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'Seawater'
Planet.Ocean.wOcean_ppt = 10  # Copied from Enceladus model
Planet.Ocean.deltaP = 0.3
Planet.Ocean.PHydroMax_MPa = 75.0

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 4.5e-12
Planet.Sil.Htidal_Wm3 = 3.38e-10

# Rock porosity
Planet.Do.POROUS_ROCK = True
Planet.Sil.porosType = 'Han2014'
Planet.Sil.HtidalMin_Wm3 = 1e-9  # Only needed for non-Han 2014 porosTypes
Planet.Sil.phiRockMax_frac = 0.10
Planet.Sil.Pclosure_MPa = 450
# Mantle equation of state model
Planet.Sil.mantleEOS = 'CI_undifferentiated_hhph_DEW17_nofluid_nomelt_685.tab'
Planet.Sil.rhoSilWithCore_kgm3 = 2700.0

""" Core assumptions """
Planet.Do.Fe_CORE = False

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0

""" Magnetic induction """
Planet.Bulk.J2 = 1454e-6  # Copied from Dione model
Planet.Bulk.C22 = 363e-6
Planet.Magnetic.ionosBounds_m = None
Planet.Magnetic.sigmaIonosPedersen_Sm = None
