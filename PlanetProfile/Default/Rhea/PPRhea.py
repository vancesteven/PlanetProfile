"""
PPRhea
Contains all body-specific parameters and information for PlanetProfile models of this body.
Import as a module and access information assigned to the attributes of the Planet struct.
Note that this file expects to be imported from the directory above.
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Rhea')

Planet.PfreezeLower_MPa = 5.0
Planet.PfreezeUpper_MPa = 150

""" Bulk planetary settings """
Planet.Bulk.R_m = 764.5e3
Planet.Bulk.M_kg = 2.310e21
Planet.Bulk.Tsurf_K = 75  # Estimate based on eyeballing diurnal variation graphs from Howett et al. (2010): https://doi.org/10.1016/j.icarus.2009.07.016
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.3721  # From Iess et al. (2007): https://doi.org/10.1016/j.icarus.2007.03.027
Planet.Bulk.Cuncertainty = 0.0036  # Estimating large uncertainty based on Hussmann et al. (2004) modeling approach
Planet.Bulk.Tb_K = 265.0

""" Layer step settings """
Planet.Steps.nIceI = 50
Planet.Steps.nSilMax = 120
Planet.Steps.nPoros = 8
Planet.Steps.iSilStart = 1  # Start from near surface to allow for fully frozen body

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'Seawater'
Planet.Ocean.wOcean_ppt = 10  # Copying Enceladus model value, as we lack information about this body
Planet.Ocean.deltaP = 0.5
Planet.Ocean.deltaT = 0.1
Planet.Ocean.PHydroMax_MPa = 175.0
Planet.Ocean.phaseType = 'recalc'

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 5.33e-12
Planet.Sil.Htidal_Wm3 = 3.38e-10

# Rock porosity
Planet.Do.POROUS_ROCK = True
Planet.Sil.porosType = 'Han2014'
Planet.Sil.HtidalMin_Wm3 = 1e-9  # Only needed for non-Han 2014 porosTypes
Planet.Sil.phiRockMax_frac = 0.92
Planet.Sil.phiRangeMult = 15
Planet.Sil.Pclosure_MPa = 550
# Mantle equation of state model
Planet.Sil.mantleEOS = 'Comet_67P_CG_v7_excluding_fluid_properties.tab'

""" Core assumptions """
Planet.Do.Fe_CORE = False

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0

""" Magnetic induction """
Planet.Bulk.J2 = 794.7e-6  # Values from Iess et al. (2007): https://doi.org/10.1016/j.icarus.2007.03.027
Planet.Bulk.C22 = 235.26e-4
Planet.Magnetic.ionosBounds_m = None
Planet.Magnetic.sigmaIonosPedersen_Sm = None
