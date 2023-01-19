"""
PPTethys
Contains all body-specific parameters and information for PlanetProfile models of this body.
Import as a module and access information assigned to the attributes of the Planet struct.
Note that this file expects to be imported from the directory above.
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Tethys')

Planet.PfreezeLower_MPa = 2
Planet.PfreezeUpper_MPa = 75

""" Bulk planetary settings """
Planet.Bulk.R_m = 531.0e3  # Value from mean radius in Archinal et al. (2018): https://doi.org/10.1007/s10569-017-9805-5
Planet.Bulk.M_kg = 6.17449e20  # From Jacobson et al. (2006): https://adsabs.harvard.edu/pdf/1992AJ....103.2068J
Planet.Bulk.Tsurf_K = 68.25  # From Gyalay and Nimmo (2023): https://doi.org/10.1029/2022JE007550
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.34  # From Gyalay and Nimmo (2023): https://doi.org/10.1029/2022JE007550
Planet.Bulk.Cuncertainty = 0.01
Planet.Do.NONHYDROSTATIC = True
Planet.Bulk.Tb_K = 270.0

""" Layer step settings """
Planet.Steps.nIceI = 80
Planet.Steps.nSilMax = 60
Planet.Steps.nPoros = 8
Planet.Steps.iSilStart = round(Planet.Steps.nIceI*0.6)  # Start from 1 (near surface) to allow for fully frozen body

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'Seawater'
Planet.Ocean.wOcean_ppt = 10  # Copying Enceladus model value, as we lack information about this body
Planet.Ocean.deltaP = 1.0
Planet.Ocean.deltaT = 0.1
Planet.Ocean.PHydroMax_MPa = 80.0
Planet.Ocean.phaseType = 'recalc'
# Ice porosity --- values from or consistent with Gyalay and Nimmo (2023): https://doi.org/10.1029/2022JE007550
Planet.Do.POROUS_ICE = True
Planet.Ocean.phiMax_frac['Ih'] = 0.4
Planet.Ocean.Pclosure_MPa['Ih'] = 40
Planet.Ocean.porosType['Ih'] = 'Han2014'

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 4.5e-12
Planet.Sil.Htidal_Wm3 = 3.38e-10

# Rock porosity
Planet.Do.POROUS_ROCK = True
Planet.Sil.porosType = 'Han2014'
Planet.Sil.HtidalMin_Wm3 = 1e-9  # Only needed for non-Han 2014 porosTypes
Planet.Sil.phiRockMax_frac = 0.6
Planet.Sil.phiRangeMult = 5
Planet.Sil.Pclosure_MPa = 850
# Mantle equation of state model
Planet.Sil.mantleEOS = 'Comet_67P_CG_v7_excluding_fluid_properties.tab'

""" Core assumptions """
Planet.Do.Fe_CORE = False

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0

""" Magnetic induction """
Planet.Bulk.J2 = 794.7e-6  # Copied from Rhea model
Planet.Bulk.C22 = 235.26e-4
Planet.Magnetic.ionosBounds_m = None
Planet.Magnetic.sigmaIonosPedersen_Sm = None
