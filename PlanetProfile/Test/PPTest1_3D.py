"""
PPTest1_3D
Europa-like 3D lateral structure example
Demonstrates minimal 3D configuration with small grid (12 pixels, fast)
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test1_3D')

Planet.PfreezeUpper_MPa = 150

""" Bulk planetary settings """
Planet.Bulk.R_m = 1561.0e3
Planet.Bulk.M_kg = 4.7991e22
Planet.Bulk.Tsurf_K = 110
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.346
Planet.Bulk.Cuncertainty = 0.005
Planet.Bulk.Tb_K = 268.4

# Orbital parameters for tidal heating (Europa)
Planet.Bulk.eccentricity = 0.0094
Planet.Bulk.meanMotion_radps = 2 * np.pi / (3.55 * 86400)  # 3.55 day period

""" Layer step settings """
Planet.Steps.nIceI = 50
Planet.Steps.nSilMax = 50
Planet.Steps.nCore = 10
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'Seawater'
Planet.Ocean.wOcean_ppt = Constants.stdSeawater_ppt
Planet.Ocean.deltaP = 1.0
Planet.Ocean.PHydroMax_MPa = 250.0

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 5.33e-12
Planet.Sil.Htidal_Wm3 = 1e-10
# Rock porosity
Planet.Do.POROUS_ROCK = False
# Mantle equation of state model (use simpler EOS to avoid PyALMA issues)
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
Planet.Core.xFeS = 0.55  # Keep moderate to avoid PyALMA negative k2 issue
Planet.Core.xFeCore = 0.0279
Planet.Core.xH2O = 0.0035

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0

""" 3D Lateral Structure Settings """
# Enable 3D lateral structure
Planet.Lateral.DO_3D = True

# Grid configuration (small grid for fast testing)
Planet.Lateral.gridType = 'healpix'  # Equal-area grid
Planet.Lateral.nSide = 1  # 12 pixels (minimal for testing)

# Enable 3D tidal heating
Planet.Lateral.DO_TIDAL_3D = True

# Add ice thickness variation via spherical harmonics
# Example: Y_20 pattern (pole-equator variation)
Planet.Lateral.dIce_pMax = 2
Planet.Lateral.dIce_Cpq_km = np.array([
    [25.0, 0.0, 0.0],  # C_00 = mean thickness 25 km
    [0.0, 0.0, 0.0],   # p=1 (no degree-1)
    [3.0, 0.0, 0.0]    # C_20 = 3 km pole-equator variation (thicker at poles)
])
Planet.Lateral.dIce_Spq_km = np.zeros((3, 3))  # All sine coefficients zero

# Mass conservation
Planet.Lateral.DO_MASS_CONSERVE = True
