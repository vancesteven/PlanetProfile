"""
PPTest17
Miranda-like, undifferentiated body
For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test17')

""" Bulk planetary settings """
Planet.Do.PARTIAL_DIFFERENTIATION = True
Planet.Bulk.qSurf_Wm2 = 10e-3
Planet.Bulk.R_m = 235.8e3
Planet.Bulk.M_kg = 0.659e20
Planet.Bulk.Tsurf_K = 60
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.346
Planet.Bulk.Cuncertainty = 0.03

""" Layer step settings """
Planet.Steps.nSilMax = 120
Planet.Steps.nPoros = 16

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 5.33e-12
Planet.Sil.Htidal_Wm3 = 1e-14
# Rock porosity
Planet.Sil.phiRockMax_frac = 0.97
Planet.Sil.Pclosure_MPa = 50
Planet.Sil.poreComp = 'Seawater'
Planet.Sil.wPore_ppt = 10
# Mantle equation of state model
Planet.Sil.mantleEOS = 'CV_undifferentiated_v4_687_DEW17_nofluid_nomelt_v2.tab'

""" Core assumptions """
Planet.Do.Fe_CORE = False

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0

""" Magnetic induction """
Planet.Bulk.J2 = 5435.2e-6
Planet.Bulk.C22 = 1549.8e-6
Planet.Magnetic.ionosBounds_m = None
Planet.Magnetic.sigmaIonosPedersen_Sm = None
