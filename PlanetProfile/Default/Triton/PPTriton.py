"""
PPTriton
Contains all body-specific parameters and information for PlanetProfile models of this body.
Import as a module and access information assigned to the attributes of the Planet struct.
Note that this file expects to be imported from the directory above.
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Triton')

# Reduce search range for melting pressure for smaller body
Planet.PfreezeUpper_MPa = 150

""" Bulk planetary settings """
Planet.Bulk.R_m = 1353.4e3
Planet.Bulk.M_kg = 2.140e22
Planet.Bulk.Tsurf_K = 38
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.31  # Copied from Pluto model, based on rough forward models with fixed layer densities in Hussmann et al. (2006): https://doi.org/10.1016/j.icarus.2006.06.005 . These authors could not find a fit for Triton.
Planet.Bulk.Cuncertainty = 0.03
Planet.Bulk.Tb_K = 266.0

""" Layer step settings """
Planet.Steps.nIceI = 200
Planet.Steps.nSilMax = 120
Planet.Steps.nCore = 10
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'MgSO4'
Planet.Ocean.wOcean_ppt = 10
Planet.Ocean.deltaP = 0.25
Planet.Ocean.deltaT = 0.1
Planet.Ocean.PHydroMax_MPa = 250.0

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 4.5e-12
Planet.Sil.Htidal_Wm3 = 0
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
Planet.Magnetic.ionosBounds_m = [250e3, 450e3]  # Approximation based on Tyler et al. (1989): https://doi.org/10.1126/science.246.4936.1466
Planet.Magnetic.sigmaIonosPedersen_Sm = [1e-16, 0.05]
