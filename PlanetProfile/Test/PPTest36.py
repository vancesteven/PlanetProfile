"""
PPTest36
Europa-like with 3D lateral tidal dissipation (Tobie et al. 2005 Figure 10 comparison)
Tests the full 3D pipeline: RunLateral3D -> ComputeTidalHeating3D -> PlotTidalHeatingByLayer
Exercises DO_3D + DO_TIDAL_3D with degree-2 ice thickness variation
For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test36')

Planet.PfreezeUpper_MPa = 150

""" Bulk planetary settings """
Planet.Bulk.R_m = 1561.0e3
Planet.Bulk.M_kg = 4.7991e22
Planet.Bulk.Tsurf_K = 110
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.346
Planet.Bulk.Cuncertainty = 0.005
Planet.Bulk.Tb_K = 267.0  # Gives ~42 km ice, Ra/RaCrit ~ 2.5 (convecting)
Planet.Bulk.eccentricity = 0.0094  # Europa orbital eccentricity
Planet.Bulk.meanMotion_radps = 2.048e-5  # Europa: 2*pi / (3.551 days * 86400)

""" Layer step settings """
Planet.Steps.nIceI = 50
Planet.Steps.nSilMax = 50
Planet.Steps.nCore = 10
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'Seawater'
Planet.Ocean.wOcean_ppt = Constants.stdSeawater_ppt
Planet.Ocean.deltaP = 1.0
Planet.Ocean.deltaT = 0.1
Planet.Ocean.PHydroMax_MPa = 250.0

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 5.33e-12
Planet.Sil.Htidal_Wm3 = 1e-10
Planet.Do.POROUS_ROCK = False
Planet.Sil.mantleEOS = 'CV3hy1wt_678_1.tab'
Planet.Sil.rhoSilWithCore_kgm3 = 3539.0

""" Core assumptions """
Planet.Do.Fe_CORE = True
Planet.Core.rhoFe_kgm3 = 8000.0
Planet.Core.rhoFeS_kgm3 = 5150.0
Planet.Core.rhoPoFeFCC = 5455.0
Planet.Core.QScore = 1e4
Planet.Core.coreEOS = 'Fe-S_3D_EOS.mat'
Planet.Core.wFe_ppt = 850
Planet.Core.xFeSmeteoritic = 0.0405
Planet.Core.xFeS = 0.55
Planet.Core.xFeCore = 0.0279
Planet.Core.xH2O = 0.0035

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0

""" Magnetic induction """
Planet.Bulk.J2 = 435.5e-6
Planet.Bulk.C22 = 131.0e-6
Planet.Magnetic.ionosBounds_m = 100e3
Planet.Magnetic.sigmaIonosPedersen_Sm = 30/100e3

""" 3D Lateral structure settings """
Planet.Lateral.DO_3D = True
Planet.Lateral.DO_TIDAL_3D = True
Planet.Lateral.gridType = 'latlon'
Planet.Lateral.nLat = 18
Planet.Lateral.nLon = 36
Planet.Lateral.DO_MASS_CONSERVE = True

# Prescribe degree-2 ice thickness variation via SH coefficients
# C00 = mean ice thickness, matching the ~42 km from the reference model (Tb=267 K)
# C20 = 3 km gives ~10 km peak-to-peak pole-to-equator variation (4pi-normalized)
Planet.Lateral.dIce_pMax = 2
pMax = Planet.Lateral.dIce_pMax
Planet.Lateral.dIce_Cpq_km = np.zeros((pMax + 1, pMax + 1))
Planet.Lateral.dIce_Spq_km = np.zeros((pMax + 1, pMax + 1))
Planet.Lateral.dIce_Cpq_km[0, 0] = 42.0  # Mean ice thickness ~42 km
Planet.Lateral.dIce_Cpq_km[2, 0] = 3.0   # Degree-2, order-0: ~10 km peak-to-peak (4pi-normalized)
