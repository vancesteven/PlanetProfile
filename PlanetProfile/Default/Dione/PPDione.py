"""
PPDione
Contains all body-specific parameters and information for PlanetProfile models of this body.
Import as a module and access information assigned to the attributes of the Planet struct.
Note that this file expects to be imported from the directory above.
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Dione')

Planet.PfreezeLower_MPa = 0.1
Planet.PfreezeUpper_MPa = 25

""" Bulk planetary settings """
Planet.Bulk.R_m = 561.4e3
Planet.Bulk.M_kg = 1.0955e21
Planet.Bulk.Tsurf_K = 75
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.33  # From Hemingway et al. (2016) LPSC abstract: https://www.hou.usra.edu/meetings/lpsc2016/pdf/1314.pdf
Planet.Bulk.Cuncertainty = 0.01
Planet.Bulk.Tb_K = 270.93

""" Layer step settings """
Planet.Steps.nIceI = 50
Planet.Steps.nRefRho = 30
Planet.Steps.nSilMax = 120
Planet.Steps.nPoros = 8
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'Seawater'
Planet.Ocean.wOcean_ppt = 10  # Copying Enceladus model value, as we lack information about this body
Planet.Ocean.deltaP = 0.1
Planet.Ocean.deltaT = 0.1
Planet.Ocean.PHydroMax_MPa = 75.0
Planet.Ocean.phaseType = 'recalc'

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 5.33e-12
Planet.Sil.Htidal_Wm3 = 3.38e-10

# Rock porosity
Planet.Do.POROUS_ROCK = True
Planet.Sil.porosType = 'Han2014'
Planet.Sil.HtidalMin_Wm3 = 1e-9  # Only needed for non-Han 2014 porosTypes
Planet.Sil.phiRockMax_frac = 0.50
# Mantle equation of state model
Planet.Sil.mantleEOS = 'Comet_67P_CG_v7_excluding_fluid_properties.tab'

""" Core assumptions """
Planet.Do.Fe_CORE = False

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0

""" Magnetic induction """
Planet.Bulk.J2 = 1454e-6  # From Hemingway et al. (2016) LPSC abstract: https://www.hou.usra.edu/meetings/lpsc2016/pdf/1314.pdf
Planet.Bulk.C22 = 363e-6
Planet.Magnetic.ionosBounds_m = None
Planet.Magnetic.sigmaIonosPedersen_Sm = None
