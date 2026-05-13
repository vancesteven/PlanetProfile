"""
PPTest35
Titan no-ocean with Kalousova HP ice convection and ice tidal heating
Based on Petricca et al. (2025): tidal dissipation in HP ice layer, no subsurface ocean
Template: Titan/PPTitanPureWaterNoOcean.py with Kalousova convection enabled
For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test35')

""" Bulk planetary settings """
Planet.Bulk.R_m = 2574.73e3  # Archinal et al. (2018)
Planet.Bulk.M_kg = 1.3452e23  # Jacobson et al. (2006)
Planet.Bulk.Tsurf_K = 94
Planet.Bulk.Psurf_MPa = 0.15
Planet.Bulk.Cmeasured = 0.343  # Petricca et al. (2025)
Planet.Bulk.Cuncertainty = 0.001
Planet.Do.NONHYDROSTATIC = False
Planet.Bulk.Tb_K = Constants.TtripleIh_III_L_K - 5  # Below triple point for no-ocean
Planet.Bulk.eccentricity = 0.0288  # Titan orbital eccentricity
Planet.Bulk.meanMotion_radps = 4.56e-6  # Titan mean motion

""" No ocean, but compute inner HP ices """
Planet.Do.HYDROSPHERE_THICKNESS = True
Planet.Bulk.Dhsphere_m = 450e3
Planet.Do.BOTTOM_ICEIII = False
Planet.Do.BOTTOM_ICEV = False
Planet.Do.NO_ICE_CONVECTION = False
Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES = True

""" Clathrate cap per Kalousova & Sotin (2020) — insulating lid reduces
    stagnant lid thickness from ~42 km to ~15 km """
Planet.Do.CLATHRATE = True
Planet.Bulk.clathType = 'top'
Planet.Bulk.clathMaxThick_m = 10e3
Planet.Steps.nClath = 30

""" Enable Kalousova convection and Arrhenius viscosity for HP ices """
Planet.Do.KALOUSOVA_CONVECTION = True
Planet.Do.ARRHENIUS_VISCOSITY = True
Planet.Do.DO_SELF_CONSISTENT_HTIDAL = True

""" Melt-bearing ice viscosity overrides per Petricca et al. (2025).
    Two-phase (partial melt) convection from Kalousova dramatically lowers
    the effective viscosity of all ice phases.  HP ice etaMelt controls the
    viscosity seen by the gravity model (Arrhenius scaling from this reference
    value).  eHP=1e13 gives peak dissipation in HP ices at Titan's orbital
    frequency (omega*tau_M ~ 1 when combined with Andrade transient creep). """
Planet.Ocean.etaMeltKalousova_Pas = {1: 1e12, 3: 1e13, 5: 1e13, 6: 1e13}

""" Ice tidal heating — converged value from self-consistent k2 iteration.
    Starting from 5e-9, iteration converges to ~1.15e-7 in 3 steps. """
Planet.Ocean.HtidalIce_Wm3 = 1.15e-7

""" Andrade rheology parameters for tidal dissipation (Petricca et al. 2025).
    Per-phase zeta: HP ices (III/V/VI) get moderate zeta to optimize dissipation,
    while Ih/Sil get small zeta for overall body compliance.
    alpha=0.2 (Andrade exponent), zeta<1 amplifies transient creep. """
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
Planet.Sil.etaRock_Pas = [1e19, 1e19]  # Moderate viscosity; Petricca range 1e18-1e22
# Rock porosity
Planet.Do.POROUS_ROCK = False
Planet.Do.CONSTANT_INNER_DENSITY = True
Planet.Sil.porosType = 'Han2014'
Planet.Sil.HtidalMin_Wm3 = 1e-10
Planet.Sil.phiRockMax_frac = 0.90
Planet.Sil.Pclosure_MPa = 4000
# Mantle equation of state model
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
Planet.Bulk.J2 = 33.089e-6  # Durante et al. (2019)
Planet.Bulk.C22 = 10.385e-6
Planet.Bulk.C21 = 0.513e-6
Planet.Bulk.S22 = -0.064e-6
Planet.Bulk.S21 = 0.612e-6
Planet.Magnetic.ionosBounds_m = [100e3, 250e3]
Planet.Magnetic.sigmaIonosPedersen_Sm = [1e-16, 80/150e3]
