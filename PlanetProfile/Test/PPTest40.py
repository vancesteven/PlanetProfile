"""
PPTest40
Ocean-bearing Titan: Andrade vs Maxwell rheology comparison.

Combines PPTest34 (ocean-bearing Titan structure) with PPTest38
(tidal heating features: Kalousova convection, Arrhenius viscosity,
self-consistent tidal heating, per-phase Andrade zeta).

This test creates a Titan model WITH a subsurface ocean between ice Ih
and HP ices (III/V/VI). The ocean decouples the ice shell from the
silicate mantle, changing the tidal response relative to the no-ocean
models in PPTest35/38/39.

Usage:
  Default (Andrade):
    Params.Gravity.backend = 'tidalpy'
    Planet, Params = PlanetProfile(Planet, Params)

  Maxwell comparison:
    Params.Gravity.backend = 'tidalpy'
    Params.Gravity.rheology_models = {
        '0': 'newton', 'Ih': 'maxwell', 'Ih_conv': 'maxwell',
        'II': 'maxwell', 'III': 'maxwell', 'III_conv': 'maxwell',
        'V': 'maxwell', 'V_conv': 'maxwell', 'VI': 'maxwell',
        'Sil': 'maxwell', 'Fe': 'elastic', 'Clath': 'elastic',
        'Clath_conv': 'maxwell'
    }

  Maxwell with tuned silicate viscosity (matches Petricca et al. k2):
    Same as above, plus:
    Planet.Sil.etaRock_Pas = [5.6e14, 5.6e14]

When DO_SELF_CONSISTENT_HTIDAL is True and TidalPy per-phase data is
available, CheckIceTidalHeatingConsistency updates both:
  - Ocean.HtidalIce_Wm3 (from TidalPy ice-phase heating)
  - Sil.Htidal_Wm3 (from TidalPy silicate heating)
Re-run for full convergence with updated heating distribution.

Physical context:
  With an ocean, the ice Ih shell responds nearly independently of
  the HP ice + silicate interior. The ocean layer (liquid, zero shear
  modulus) transmits only normal stress, so tidal deformation of the
  outer shell is controlled by its own thickness and rigidity.

  Under Andrade rheology, dissipation occurs primarily in HP ices
  (transient creep amplified by low zeta). Under Maxwell, HP ices
  with eta ~ 1e12-1e13 are in the low-omega*tau regime (nearly
  Newtonian), so dissipation requires either:
    (a) Much lower ice viscosity (approaching fluid), or
    (b) Silicate dissipation near the Maxwell peak (eta ~ mu/omega)

  The sweep results (PPTest39 analysis) show that matching Re(k2)~0.6
  and Im(k2)~-0.12 under Maxwell requires eta_sil ~ 5e14 Pa-s (4
  orders below typical mantle viscosity), shifting all dissipation
  from HP ices to silicates.

For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test40')

# Titan analog with ocean
Planet.parent = 'Saturn'  # Needed for TidalPy per-layer heating (parent body mass)

""" Bulk planetary settings """
Planet.Bulk.R_m = 2574.73e3
Planet.Bulk.M_kg = 1.3452e23
Planet.Bulk.Tsurf_K = 94
Planet.Bulk.Psurf_MPa = 0.15
Planet.Bulk.Cmeasured = 0.33
Planet.Bulk.Cuncertainty = 0.01
Planet.Do.NONHYDROSTATIC = False
Planet.Bulk.Tb_K = 255.0  # Above ice Ih melting at base → creates ocean
Planet.Bulk.eccentricity = 0.0288
Planet.Bulk.meanMotion_radps = 4.56e-6

""" Ocean — subsurface liquid water layer """
Planet.Do.HYDROSPHERE_THICKNESS = True
Planet.Bulk.Dhsphere_m = 450e3
Planet.Do.BOTTOM_ICEIII = False
Planet.Do.BOTTOM_ICEV = False
Planet.Do.NO_ICE_CONVECTION = False

""" Clathrate cap """
Planet.Do.CLATHRATE = True
Planet.Bulk.clathType = 'top'
Planet.Bulk.clathMaxThick_m = 10e3
Planet.Steps.nClath = 30

""" Enable Kalousova convection and self-consistent tidal heating """
Planet.Do.KALOUSOVA_CONVECTION = True
Planet.Do.ARRHENIUS_VISCOSITY = True
Planet.Do.DO_SELF_CONSISTENT_HTIDAL = True

""" Melt-bearing ice viscosity """
Planet.Ocean.etaMeltKalousova_Pas = {1: 1e12, 3: 1e13, 5: 1e13, 6: 1e13}

""" Ice tidal heating — initial guess; self-consistent loop will override """
Planet.Ocean.HtidalIce_Wm3 = 1e-8

""" Andrade rheology parameters (active under Andrade, inert under Maxwell) """
Planet.Gravity.andradExponent = 0.2
Planet.Gravity.andrade_zeta = {
    'Ih': 0.005, 'III': 0.05, 'V': 0.05, 'VI': 0.05,
    'Sil': 0.005, 'Fe': 1.0, 'Clath': 1.0, 'default': 1.0
}

""" Layer step settings — ensure all layers have >= 5 points for TidalPy """
Planet.Steps.nIceI = 200
Planet.Steps.nIceIIILitho = 20
Planet.Steps.nIceVLitho = 20
Planet.Steps.nSilMax = 200
Planet.Steps.nCore = 10
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'PureH2O'
Planet.Ocean.wOcean_ppt = 0
Planet.Ocean.deltaP = 8.0
Planet.Ocean.deltaT = 0.5
Planet.Ocean.phaseType = 'lookup'
Planet.Ocean.PHydroMax_MPa = 1800.0
Planet.Ocean.THydroMax_K = 350.0

""" Silicate Mantle """
Planet.Sil.Qrad_Wkg = 1.5e-12
Planet.Sil.Htidal_Wm3 = 1e-10
Planet.Sil.etaRock_Pas = [1e19, 1e19]
Planet.Do.POROUS_ROCK = False
Planet.Do.CONSTANT_INNER_DENSITY = True
Planet.Sil.porosType = 'Han2014'
Planet.Sil.HtidalMin_Wm3 = 1e-10
Planet.Sil.phiRockMax_frac = 0.90
Planet.Sil.Pclosure_MPa = 4000
Planet.Sil.mantleEOS = 'Comet_67P_CG_v7_excluding_fluid_properties.tab'
Planet.Sil.rhoSilWithCore_kgm3 = 3539.0

""" Core assumptions """
Planet.Do.Fe_CORE = False
Planet.Core.rhoFe_kgm3 = 8000.0
Planet.Core.rhoFeS_kgm3 = 5150.0
Planet.Core.rhoPoFeFCC = 5455.0
Planet.Core.QScore = 1e4
Planet.Core.coreEOS = 'Fe-S_3D_EOS.mat'
Planet.Core.wFe_ppt = 700
Planet.Core.xFeSmeteoritic = 0.0405
Planet.Core.xFeS = 0.55
Planet.Core.xFeCore = 0.0279
Planet.Core.xH2O = 0.0035

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0

""" Magnetic induction """
Planet.Bulk.J2 = 33.089e-6
Planet.Bulk.C22 = 10.385e-6
Planet.Bulk.C21 = 0.513e-6
Planet.Bulk.S22 = -0.064e-6
Planet.Bulk.S21 = 0.612e-6
Planet.Magnetic.ionosBounds_m = [100e3, 250e3]
Planet.Magnetic.sigmaIonosPedersen_Sm = [1e-16, 80/150e3]
