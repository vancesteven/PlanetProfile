"""
PPIo
Contains all body-specific parameters and information for PlanetProfile models of this body.
Import as a module and access information assigned to the attributes of the Planet struct.
Note that this file expects to be imported from the directory above.
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Io')

""" Bulk planetary settings """
Planet.Do.NO_H2O = True
Planet.Bulk.qSurf_Wm2 = 140e-3
# Planet.Bulk.qSurf_Wm2 = 2240e-3  # From Lainey et al. (2009): DOI: 10.1038/nature08108
Planet.Bulk.R_m = 1821.3e3
Planet.Bulk.M_kg = 8.9319e22
Planet.Bulk.Tsurf_K = 110
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.37685  # Value from Anderson et al. (2001): https://doi.org/10.1029/2000JE001367
Planet.Bulk.Cuncertainty = 0.00035

""" Layer step settings """
Planet.Steps.nSilMax = 300
Planet.Steps.nCore = 10

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 5.33e-12
Planet.Sil.Htidal_Wm3 = 5e-9
# Planet.Sil.Htidal_Wm3 = 3.77e-6  # Inferred from Lainey et al. (2009), with the global rate averaged over the silicate volume.
# Rock porosity
Planet.Do.POROUS_ROCK = True
Planet.Sil.phiRockMax_frac = 0.70
Planet.Sil.Pclosure_MPa = 750
# Mantle equation of state model
Planet.Sil.mantleEOS = 'CV3hy1wt_678_1.tab'

""" Core assumptions """
Planet.Do.Fe_CORE = True
Planet.Core.rhoFe_kgm3 = 8000.0
Planet.Core.rhoFeS_kgm3 = 5150.0
Planet.Core.rhoPoFeFCC = 5455.0
Planet.Core.QScore = 1e4
Planet.Core.coreEOS = 'Fe-S_3D_EOS.mat'
Planet.Core.wFe_ppt = 875

Planet.Core.xFeSmeteoritic = 0.0405
Planet.Core.xFeS = 0.55
Planet.Core.xFeCore = 0.0279
Planet.Core.xH2O = 0.0035

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0

""" Magnetic induction """
Planet.Magnetic.J2 = 1845.9e-6  # Values from Anderson et al. (2001): https://doi.org/10.1029/2000JE001367
Planet.Magnetic.C22 = 553.7e-6
Planet.Magnetic.ionosBounds_m = 100e3
Planet.Magnetic.sigmaIonosPedersen_Sm = 30/100e3
