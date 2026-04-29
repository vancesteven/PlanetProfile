"""
PPTest41NoClath
Andrade no-ocean Titan WITHOUT clathrate cap: MCMC exploration variant.

Identical to PPTest41 except Planet.Do.CLATHRATE = False.
Used for inference structure cache generation to compare clathrate vs. no-clathrate models.

For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test41NoClath')

""" Body """
Planet.name = 'Titan'
Planet.Bulk.R_m = 2574.73e3
Planet.Bulk.M_kg = 1.3452e23
Planet.Bulk.Tsurf_K = 94.0
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.3438
Planet.Bulk.Cuncertainty = 0.0005
Planet.parent = 'Saturn'

""" No ocean, but compute inner HP ices """
Planet.Do.HYDROSPHERE_THICKNESS = True
Planet.Bulk.Dhsphere_m = 450e3
Planet.Do.BOTTOM_ICEIII = False
Planet.Do.BOTTOM_ICEV = False
Planet.Do.NO_ICE_CONVECTION = False
Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES = True

""" NO clathrate cap (difference from PPTest41) """
Planet.Do.CLATHRATE = False
Planet.Bulk.clathType = 'top'
Planet.Bulk.clathMaxThick_m = 0.0
Planet.Steps.nClath = 0

""" Enable Kalousova convection and Arrhenius viscosity """
Planet.Do.KALOUSOVA_CONVECTION = True
Planet.Do.ARRHENIUS_VISCOSITY = True

""" Layer step settings """
Planet.Steps.nIceI = 100
Planet.Steps.nSilMax = 60
Planet.Steps.nCore = 10

""" Surface heat flux (from Cassini) """
Planet.Bulk.qSurf_Wm2 = 4.8e-3

""" Silicate and core properties """
Planet.Sil.EOS = 'Comet_67P_CG_v7_excluding_fluid_properties.tab'
Planet.Sil.rhoSilWithCore_kgm3 = 2500.0
Planet.Core.rhoFe_kgm3 = 8000.0
Planet.Core.rhoFeS_kgm3 = 5150.0
Planet.Core.xFeS = 0.25

""" Ocean composition """
Planet.Ocean.comp = 'PureH2O'
Planet.Ocean.w_ocean_pct = 0.0
Planet.Ocean.deltaP = 0.0

""" Gravity (TidalPy backend for MCMC) """
Planet.Do.GRAVITY = True
