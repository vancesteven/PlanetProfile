"""
PPCallisto
Contains all body-specific parameters and information for PlanetProfile models of this body.
Import as a module and access information assigned to the attributes of the Planet struct.
Note that this file expects to be imported from the directory above.
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Callisto')

""" Bulk planetary settings """
Planet.Bulk.R_m = 2410.3e3
Planet.Bulk.M_kg = 1.0759e23
Planet.Bulk.Tsurf_K = 110
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.3549  # Value from Anderson et al. (2001): https://doi.org/10.1006/icar.2001.6664
Planet.Bulk.Cuncertainty = 0.0042
Planet.Bulk.Tb_K = 262.0

""" Layer step settings """
Planet.Steps.nIceI = 200
Planet.Steps.nRefRho = 30
Planet.Steps.nSilMax = 300
Planet.Steps.nCore = 10
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'MgSO4'
Planet.Ocean.wOcean_ppt = 100
Planet.Ocean.deltaP = 4.0
Planet.Ocean.deltaT = 0.5
Planet.Ocean.phaseType = 'lookup'
Planet.Ocean.PHydroMax_MPa = 2000.0
Planet.Ocean.THydroMax_K = 350.0

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 5.33e-14
Planet.Sil.Htidal_Wm3 = 1e-18
# Rock porosity
Planet.Do.POROUS_ROCK = False
# Mantle equation of state model
Planet.Sil.mantleEOS = 'CI_undifferentiated_hhph_DEW17_nofluid_nomelt_685.tab'
Planet.Sil.rhoSilWithCore_kgm3 = 3539.0

""" Core assumptions """
Planet.Do.Fe_CORE = True
Planet.Core.rhoFe_kgm3 = 8000.0
Planet.Core.rhoFeS_kgm3 = 5150.0
Planet.Core.rhoPoFeFCC = 5455.0
Planet.Core.QScore = 1e4
Planet.Core.coreEOS = 'Fe-S_3D_EOS.mat'
Planet.Core.wFe_ppt = 750

Planet.Core.xFeSmeteoritic = 0.0405
Planet.Core.xFeS = 0.55
Planet.Core.xFeCore = 0.0279
Planet.Core.xH2O = 0.0035

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0

""" Magnetic induction """
Planet.Magnetic.J2 = 32.7e-6  # Cnm values from Anderson et al. (2001): https://doi.org/10.1006/icar.2001.6664
Planet.Magnetic.C22 = 10.2e-6
Planet.Magnetic.S22 = -1.1e-6
Planet.Magnetic.C21 = 0.0
Planet.Magnetic.S21 = 0.0
Planet.Magnetic.ionosBounds_m = 100e3
Planet.Magnetic.sigmaIonosPedersen_Sm = 800/100e3
