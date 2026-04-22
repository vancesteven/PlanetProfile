"""
PPTest39
Maxwell rheology test on Titan no-ocean model.
Identical to PPTest38 except uses Maxwell rheology for all solid layers
(instead of Andrade). Tests the claim by Petricca et al. (2025, Nature)
that tidal dissipation results are robust and not strongly dependent on
the chosen rheological model.

Usage:
  This test defines the Planet structure. To activate Maxwell rheology,
  the test runner must also set:
    Params.Gravity.rheology_models = {
        '0': 'newton', 'Ih': 'maxwell', 'Ih_conv': 'maxwell',
        'II': 'maxwell', 'III': 'maxwell', 'III_conv': 'maxwell',
        'V': 'maxwell', 'V_conv': 'maxwell', 'VI': 'maxwell',
        'Sil': 'maxwell', 'Fe': 'elastic', 'Clath': 'elastic',
        'Clath_conv': 'maxwell'
    }
  Andrade parameters (andradExponent, andrade_zeta) are retained but
  inert under Maxwell rheology — they only affect the Andrade creep term.

Petricca et al. (2025) key results for no-ocean Titan:
  Observed (their reanalysis): Re(k2) = 0.341 +/- 0.349, Im(k2) = 0.048 +/- 0.035
  Model (no ocean, Andrade):   Re(k2) ~ 0.6, Im(k2) ~ 0.12, Q ~ 5
  They claim Maxwell gives similar posterior distributions.

  Their model structure (Extended Data Table 2):
    - Ice Ih conductive lid: ~3 km clathrate + ~160 km ice Ih
    - HP ice layer: III + V + VI, total ~320 km, viscosity ~10^12 Pa·s near melting
    - Silicate mantle: ~480 km, rho ~ 2500 kg/m^3
    - No iron core in best-fit model

  Tidal dissipation:
    - Total E_dot from Im(k2): ~3500 GW
    - Concentrated in HP ice near melting point
    - Dissipation drives vigorous convection in HP ices, transporting
      melt pockets and organic-rich material upward

  Their claim that results are not strongly dependent on rheology
  implies that the dominant dissipation mechanism is Maxwell relaxation
  in the low-viscosity HP ice, not Andrade transient creep. This test
  checks whether that is the case.

Physical expectation:
  For Maxwell rheology, dissipation peaks when omega*tau_M ~ 1,
  i.e., eta ~ mu/omega. For Titan (omega = 4.56e-6 rad/s):
    - eta_peak = mu/omega ~ 4e9/4.56e-6 ~ 8.8e14 Pa·s
  The convecting HP ices have eta ~ 10^12-10^13 Pa·s, which is
  well below the Maxwell peak. This means Maxwell dissipation in
  the HP ices is dominated by the viscous (Newtonian) limit:
    Im(mu*) ~ omega*eta*mu^2/(mu^2 + omega^2*eta^2) ~ omega*eta
  Andrade creep adds dissipation beyond this Newtonian baseline.
  Whether Andrade vs Maxwell matters depends on how much the
  Andrade transient term contributes relative to the Maxwell term.

For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test39')

# Titan analog: same as PPTest38
Planet.parent = 'Saturn'  # Needed for TidalPy per-layer heating (parent body mass)

""" Bulk planetary settings """
Planet.Bulk.R_m = 2574.73e3
Planet.Bulk.M_kg = 1.3452e23
Planet.Bulk.Tsurf_K = 94
Planet.Bulk.Psurf_MPa = 0.15
Planet.Bulk.Cmeasured = 0.343
Planet.Bulk.Cuncertainty = 0.001
Planet.Do.NONHYDROSTATIC = False
Planet.Bulk.Tb_K = Constants.triplePointT_K - 5
Planet.Bulk.eccentricity = 0.0288
Planet.Bulk.meanMotion_radps = 4.56e-6

""" No ocean, but compute inner HP ices """
Planet.Do.HYDROSPHERE_THICKNESS = True
Planet.Bulk.Dhsphere_m = 450e3
Planet.Do.BOTTOM_ICEIII = False
Planet.Do.BOTTOM_ICEV = False
Planet.Do.NO_ICE_CONVECTION = False
Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES = True

""" Clathrate cap """
Planet.Do.CLATHRATE = True
Planet.Bulk.clathType = 'top'
Planet.Bulk.clathMaxThick_m = 10e3
Planet.Steps.nClath = 30

""" Enable Kalousova convection """
Planet.Do.KALOUSOVA_CONVECTION = True
Planet.Do.ARRHENIUS_VISCOSITY = True
Planet.Do.DO_SELF_CONSISTENT_HTIDAL = True

""" Melt-bearing ice viscosity """
Planet.Ocean.etaMeltKalousova_Pas = {1: 1e12, 3: 1e13, 5: 1e13, 6: 1e13}

""" Ice tidal heating — converged value from Andrade model """
Planet.Ocean.HtidalIce_Wm3 = 1.15e-7

""" Andrade parameters — retained for reference but inert under Maxwell.
    To activate, set Params.Gravity.rheology_models back to 'andrade'. """
Planet.Gravity.andradExponent = 0.2
Planet.Gravity.andrade_zeta = {
    'Ih': 0.005, 'III': 0.05, 'V': 0.05, 'VI': 0.05,
    'Sil': 0.005, 'Fe': 1.0, 'Clath': 1.0, 'default': 1.0
}

""" Layer step settings """
Planet.Steps.nIceI = 200
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
