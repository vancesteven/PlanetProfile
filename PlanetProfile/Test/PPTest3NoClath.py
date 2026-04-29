"""
PPTest3NoClath
Europa-like, pure water model WITHOUT clathrate underplate layer.

Identical to PPTest3 except Planet.Do.CLATHRATE = False.
Used for inference structure cache generation to compare clathrate vs. no-clathrate models.

For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test3NoClath')

Planet.PfreezeLower_MPa = 0
Planet.PfreezeUpper_MPa = 150

""" Bulk planetary settings """
Planet.Bulk.R_m = 1561.0e3
Planet.Bulk.M_kg = 4.7991e22
Planet.Bulk.Tsurf_K = 110
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.346
Planet.Bulk.Cuncertainty = 0.005
Planet.Bulk.Tb_K = 295.0

""" NO clathrate (difference from PPTest3) """
Planet.Do.CLATHRATE = False
Planet.Bulk.clathType = 'bottom'
Planet.Bulk.clathMaxThick_m = 0.0
Planet.Steps.nClath = 0

Planet.Bulk.qSurf_Wm2 = 15.0e-3

""" Layer step settings """
Planet.Steps.nIceI = 50
Planet.Steps.nSilMax = 50
Planet.Steps.nCore = 10

""" Ocean composition """
Planet.Ocean.comp = 'PureH2O'
Planet.Ocean.w_ocean_pct = 0.0
Planet.Ocean.deltaP = 0.0

""" Silicate and core properties """
Planet.Sil.EOS = 'CV_hhph_DEW17_fluid_properties_new.tab'
Planet.Sil.rhoSilWithCore_kgm3 = 3539.0
Planet.Core.rhoFe_kgm3 = 8000.0
Planet.Core.rhoFeS_kgm3 = 5150.0
Planet.Core.xFeS = 0.0

""" Gravity """
Planet.Do.GRAVITY = True
