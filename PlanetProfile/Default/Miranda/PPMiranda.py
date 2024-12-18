"""
PPMiranda
Contains all body-specific parameters and information for PlanetProfile models of this body.
Import as a module and access information assigned to the attributes of the Planet struct.
Note that this file expects to be imported from the directory above.
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Miranda')

Planet.PfreezeLower_MPa = 0.1
Planet.PfreezeUpper_MPa = 20

""" Bulk planetary settings """
Planet.Bulk.R_m = 235.8e3  # Value from mean radius in Archinal et al. (2018): https://doi.org/10.1007/s10569-017-9805-5
Planet.Bulk.M_kg = 0.659e20  # From Jacobson et al. (1992): https://adsabs.harvard.edu/pdf/1992AJ....103.2068J
Planet.Bulk.Tsurf_K = 60
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.346  # Based on rough forward models with fixed layer densities in Hussmann et al. (2006): https://doi.org/10.1016/j.icarus.2006.06.005
Planet.Bulk.Cuncertainty = 0.03
Planet.Bulk.Tb_K = 272.356

""" Layer step settings """
Planet.Steps.nIceI = 50
Planet.Steps.nSilMax = 120
Planet.Steps.nPoros = 8
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'Seawater'
Planet.Ocean.wOcean_ppt = 10  # Copied from Enceladus model due to similarities in orbital period, size, and density
Planet.Ocean.deltaP = 0.1
Planet.Ocean.PHydroMax_MPa = 25.0

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 4.5e-12
Planet.Sil.Htidal_Wm3 = 3.38e-10

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
Planet.Bulk.J2 = 5435.2e-6  # Copied from Enceladus model due to similarities in orbital period, size, and density
Planet.Bulk.C22 = 1549.8e-6
Planet.Magnetic.ionosBounds_m = None
Planet.Magnetic.sigmaIonosPedersen_Sm = None
