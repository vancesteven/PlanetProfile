"""
PPTest46
Andrade Europa MCMC base structure: Bayesian inference of tidal dissipation
parameters for Europa, incorporating updated MOI from Gomez Casajus et al. (2021).

Interior structure:
  - Thin Ice Ih shell (~25 km, no HP ice phases)
  - Liquid ocean (Seawater, ~100 km)
  - Silicate mantle (CONSTANT_INNER_DENSITY)
  - Fe-FeS metallic core (Fe_CORE=True, xFeS=0.55 → rhoCore ≈ 5865 kg/m³)

This is the 'base' reference structure (Variant B: intermediate core).
Four reference structures spanning Petricca et al. (2025) credible core range:

  A. No core      (Fe_CORE=False, rhoSil=3000): partially undifferentiated
  B. This file    (Fe_CORE=True,  rhoSil=3300, xFeS=0.55, rhoCore≈5865): intermediate
  C. Dense core   (Fe_CORE=True,  rhoSil=3500, xFeS=0.10, rhoCore≈7391): small Fe-rich
  D. Sulfide core (Fe_CORE=True,  rhoSil=3200, xFeS=0.88, rhoCore≈5301): large FeS-rich

rhoCore mixing rule (molar-weighted harmonic mean, LayerPropagators.py:1828-1830):
  rhoCore = rhoFeS * rhoFe * (xFeS*M_FeS + (1-xFeS)*M_Fe)
            / (xFeS*M_FeS*rhoFe + (1-xFeS)*M_Fe*rhoFeS)

Observational constraints:
  CMR2    = 0.3547 +/- 0.0024   (Gomez Casajus et al. 2021)
  Re(k2)  = 0.328  +/- 0.15     (theoretical estimate; Europa Clipper will constrain)
  Im(k2)  = 0.015  +/- 0.010    (theoretical estimate; placeholder)

MCMC parameter space (Andrade, Ice Ih only — no HP ice phases for thin shell):
  alpha (Andrade exponent):          [0.2, 0.4]
  log10(zeta_Ih) (compliance, Ih):   [-2, 2]
  log10(zeta_sil) (compliance, Sil): [-2, 2]
  log10(eta_Ih) (Ice Ih viscosity):  [12, 16]  Pa s
  log10(eta_sil) (Silicate visc.):   [18, 22]  Pa s

References:
  Gomez Casajus et al. (2021): MOI = 0.3547 +/- 0.0024
  Petricca et al. (2025): core radius 49-269 km (1σ), rhoSil 2906-3630 kg/m³
  Howell (2021): CBE ice thickness 24.3+22.8/-1.5 km; CBE η_Ih ≈ 10^14.7 Pa s

For testing purposes
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Test46')

# Europa around Jupiter
Planet.parent = 'Jupiter'

""" Bulk planetary settings """
Planet.Bulk.R_m = 1560.8e3           # Archinal et al. (2018)
Planet.Bulk.M_kg = 4.800e22           # Hussmann et al. (2006)
Planet.Bulk.Tsurf_K = 110
Planet.Bulk.Psurf_MPa = 0.0
Planet.Bulk.Cmeasured = 0.3547        # Gomez Casajus et al. (2021) reanalysis
Planet.Bulk.Cuncertainty = 0.0024
Planet.Do.NONHYDROSTATIC = False

# Tb_K ≈ 271.0 K → ~25 km ice shell (CBE from Howell 2021).
# Upper bound: 272.7 K (melting-point depression at ~15 MPa base pressure).
# Lower bound: ~269 K (thicker shell limit for MCMC prior).
Planet.Bulk.Tb_K = 271.0

# Europa orbital parameters
Planet.Bulk.eccentricity = 0.0094
Planet.Bulk.meanMotion_radps = 2.048e-5  # 2π / (3.551181 d × 86400 s/d)

""" No HP ice phases — thin shell is purely Ice Ih """
Planet.Do.BOTTOM_ICEIII = False
Planet.Do.BOTTOM_ICEV = False
Planet.Do.NO_ICE_CONVECTION = False

""" No clathrate cap (thin ice, uncertain for Europa) """
Planet.Do.CLATHRATE = False

""" Andrade rheology parameters (defaults; MCMC overrides alpha, zeta, eta) """
Planet.Gravity.andradExponent = 0.3
Planet.Gravity.andrade_zeta = {
    'Ih': 0.005, 'Sil': 0.005, 'Fe': 1.0, 'default': 1.0
}

""" Ice tidal heating — representative value; MCMC structure is fixed """
Planet.Ocean.HtidalIce_Wm3 = 1e-7

""" Layer step settings """
Planet.Steps.nIceI = 200
Planet.Steps.nSilMax = 300
Planet.Steps.nCore = 10
Planet.Steps.iSilStart = Planet.Steps.nIceI

""" Hydrosphere assumptions/settings """
Planet.Ocean.comp = 'Seawater'
Planet.Ocean.wOcean_ppt = Constants.stdSeawater_ppt
Planet.Ocean.deltaP = 1.0
Planet.Ocean.deltaT = 0.1
# Max pressure at seafloor: ~30 MPa (ice) + ~140 MPa (ocean) = ~170 MPa
Planet.Ocean.PHydroMax_MPa = 300.0
Planet.Ocean.THydroMax_K = 320.0

""" Silicate Mantle — Variant B: CM-scale composition (Petricca 2025) """
Planet.Sil.Qrad_Wkg = 5.33e-12       # Hussmann & Spohn (2004)
Planet.Sil.Htidal_Wm3 = 1e-10
Planet.Sil.etaRock_Pas = [1e19, 1e19]
Planet.Do.POROUS_ROCK = False
Planet.Do.CONSTANT_INNER_DENSITY = True
Planet.Sil.mantleEOS = 'CM_hydrous_differentiated_Ganymede_Core080Fe020S_excluding_fluid_properties.tab'
# Silicate density for Variant B (CM-scale mantle, Petricca 2025 range: 2906-3630 kg/m³)
Planet.Sil.rhoSilWithCore_kgm3 = 3300.0

""" Core assumptions — Variant B: intermediate Fe-FeS core """
Planet.Do.Fe_CORE = True
Planet.Core.rhoFe_kgm3 = 8000.0
Planet.Core.rhoFeS_kgm3 = 5150.0
Planet.Core.rhoPoFeFCC = 5455.0
Planet.Core.QScore = 1e4
Planet.Core.coreEOS = 'Fe-S_3D_EOS.mat'
Planet.Core.xFeSmeteoritic = 0.0405
# xFeS = 0.55 → rhoCore ≈ 5865 kg/m³ (molar-weighted harmonic mixing)
# Core radius self-consistently determined from mass budget given rhoSil and rhoCore.
# Expected ~150 km based on Petricca 2025 credible range (49-269 km 1σ).
Planet.Core.xFeS = 0.55
Planet.Core.xFeCore = 0.0279
Planet.Core.xH2O = 0.0035

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0

""" Magnetic induction (Europa Galileo gravity harmonics) """
Planet.Bulk.J2 = 435.5e-6            # Anderson et al. (1998)
Planet.Bulk.C22 = 131.0e-6
Planet.Magnetic.ionosBounds_m = 100e3
Planet.Magnetic.sigmaIonosPedersen_Sm = 30 / 100e3

# ============================================================================
# Reference structure variants (create separate PPTest files for each)
# ============================================================================
# Variant A — No core (partially undifferentiated, upper limit per Petricca 2025)
#   Planet.Do.Fe_CORE = False
#   Planet.Sil.rhoSilWithCore_kgm3 = 3000.0
#   Planet.Bulk.Tb_K = 271.0
#
# Variant B — This file: intermediate Fe-FeS core
#   Planet.Sil.rhoSilWithCore_kgm3 = 3300.0, xFeS = 0.55 → rhoCore ≈ 5865 kg/m³
#
# Variant C — Small, dense Fe core (mostly iron, minimal sulfur)
#   Planet.Sil.rhoSilWithCore_kgm3 = 3500.0, xFeS = 0.10 → rhoCore ≈ 7391 kg/m³
#
# Variant D — Large, sulfide-rich core (maximum sulfur content)
#   Planet.Sil.rhoSilWithCore_kgm3 = 3200.0, xFeS = 0.88 → rhoCore ≈ 5301 kg/m³
