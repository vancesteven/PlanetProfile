"""
PPEuropaInductOgram
Contains all body-specific parameters and information for PlanetProfile models of this body.
Import as a module and access information assigned to the attributes of the Planet struct.
Note that this file expects to be imported from the directory above.
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Europa')

# Reduce search range for melting pressure to values realistic for Europa
Planet.PfreezeUpper_MPa = 150

""" Bulk planetary settings """
Planet.Bulk.R_m = 1560.8e3
Planet.Bulk.M_kg = 4.800e22
Planet.Bulk.Tsurf_K = 110
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.346
Planet.Bulk.Cuncertainty = 0.005
Planet.Bulk.Tb_K = 269.2

""" Layer step settings """
Planet.Steps.nIceI = 100
Planet.Steps.nSilMax = 150
Planet.Steps.nCore = 10
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'Seawater'
#Planet.Ocean.comp = 'MgSO4'
Planet.Ocean.wOcean_ppt = Constants.stdSeawater_ppt
#Planet.Ocean.wOcean_ppt = 100
Planet.Ocean.phaseType = 'calculate'
Planet.Ocean.deltaP = 2.0
Planet.Ocean.PHydroMax_MPa = 350.0

""" Silicate Mantle """
Planet.Do.CONSTANT_INNER_DENSITY = True
Planet.Sil.Qrad_Wkg = 5.33e-12  # Estimate from Hussmann and Spohn (2004): https://doi.org/10.1016/j.icarus.2004.05.020, calculated from total heating rate, silicate density, and layer radii
Planet.Sil.Htidal_Wm3 = 1e-10  # Approximate max. tidal heating in silicates as modeled by Tobie et al. (2003): https://doi.org/10.1029/2003JE002099
# Rock porosity
Planet.Do.POROUS_ROCK = False
# Mantle equation of state model
Planet.Sil.mantleEOS = 'CV3hy1wt_678_1.tab'  # (2900 for Q= 100 GW, 3240 for Q= 220 GW)
Planet.Sil.rhoSilWithCore_kgm3 = 3539.0  # This is the 1-bar, 275 K value from CV3hy1wt_678_1.tab

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
Planet.Magnetic.SCera = 'Galileo'
Planet.Magnetic.extModel = 'JRM33C2020'
