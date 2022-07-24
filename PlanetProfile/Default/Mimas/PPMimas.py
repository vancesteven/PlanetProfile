"""
PPMimas
Contains all body-specific parameters and information for PlanetProfile models of this body.
Import as a module and access information assigned to the attributes of the Planet struct.
Note that this file expects to be imported from the directory above.
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Mimas')

Planet.PfreezeLower_MPa = 0.1
Planet.PfreezeUpper_MPa = 25

""" Bulk planetary settings """
Planet.Bulk.R_m = 198.8e3
Planet.Bulk.M_kg = 3.79e19
Planet.Bulk.Tsurf_K = 80  # "Canonical" value used by Rhoden and Walker (2022): https://doi.org/10.1016/j.icarus.2021.114872
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.355  # From Hussmann et al. (2006): http://dx.doi.org/10.1016/j.icarus.2006.06.005
Planet.Bulk.Cuncertainty = 0.010  # Estimating large uncertainty, as Mimas is regarded to not be in hydrostatic equilibrium (e.g. Tajeddine et al. (2016): https://doi.org/10.1126/science.1255299)
Planet.Bulk.Tb_K = 272.5  # Ice shell should be 24-31 km to match libration constraints, see Tajeddine et al. (2014): https://doi.org/10.1126/science.1255299

""" Layer step settings """
Planet.Steps.nIceI = 50
Planet.Steps.nSilMax = 120
Planet.Steps.nPoros = 8
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'Seawater'
Planet.Ocean.wOcean_ppt = 10  # Copying Enceladus model value, as we lack information about this body
Planet.Ocean.deltaP = 0.1
Planet.Ocean.deltaT = 0.1
Planet.Ocean.PHydroMax_MPa = 25.0
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
Planet.Bulk.J2 = None
Planet.Bulk.C22 = None
Planet.Magnetic.ionosBounds_m = None
Planet.Magnetic.sigmaIonosPedersen_Sm = None
