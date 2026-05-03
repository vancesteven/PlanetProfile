"""
PPTest47
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
  log10(zeta):                       [-2, 2]
  log10(eta_Ih):                     [12, 16]
  log10(eta_sil):                    [18, 22]
"""
from copy import deepcopy

from PlanetProfile.Utilities.defineStructs import Constants, PlanetStruct

Planet = PlanetStruct('Test47')

# ── Bulk properties (Europa) ──────────────────────────────────────────────────
Planet.Bulk.R_m         = 1560.8e3
Planet.Bulk.M_kg        = 4.7998e22
Planet.Bulk.Tsurf_K     = 110.0
Planet.Bulk.Psurf_MPa   = 0.0
Planet.Bulk.Cmeasured   = 0.3547
Planet.Bulk.Cuncertainty = 0.0024
Planet.Bulk.Tb_K        = 271.0

# ── Ocean / hydrosphere ───────────────────────────────────────────────────────
Planet.Ocean.comp       = 'Seawater'
Planet.Ocean.wOcean_ppt = Constants.stdSeawater_ppt

Planet.Do.NO_ICE_CONVECTION = True

# ── No HP ice for thin shell ───────────────────────────────────────────────────
Planet.Do.BOTTOM_ICEIII = False
Planet.Do.BOTTOM_ICEV   = False

# ── Silicate mantle ────────────────────────────────────────────────────────────
Planet.Do.CONSTANT_INNER_DENSITY = True
Planet.Sil.rhoSilWithCore_kgm3   = 3300.0
Planet.Sil.Qrad_Wkg              = 0.0
Planet.Sil.Htidal_Wm3            = 0.0

# ── Fe-FeS core (Variant B: intermediate) ─────────────────────────────────────
Planet.Do.Fe_CORE = True
Planet.Core.xFeS  = 0.55

# ── Tidal / orbital parameters ────────────────────────────────────────────────
Planet.Bulk.Tidal_rad_s  = 2 * 3.14159265358979 / (3.551181 * 86400)   # Europa's orbital period
Planet.Bulk.Tidal_body   = 'Jupiter'
Planet.Bulk.Tidal_ecc    = 0.0094

# ── Andrade rheology parameters ───────────────────────────────────────────────
Planet.Gravity.andradExponent = 0.3
Planet.Gravity.andrade_zeta   = {
    'Ih':  0.005,
    'Sil': 0.005,
    'Fe':  1.0,
    'default': 1.0,
}

# ── Layer resolution ───────────────────────────────────────────────────────────
Planet.Steps.nOceanMax   = 100
Planet.Steps.nSurfIce    = 50
Planet.Steps.nSil        = 200
Planet.Steps.nCore       = 100
