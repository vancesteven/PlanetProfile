"""
PPLuna
Lunar dry-body reference configuration for PlanetProfile.

Scope and status
----------------
This is the Stage-1 lunar structural configuration: a genuine Moon, replacing
the Io template that previously occupied this file. It defines the spherically
symmetric reference PlanetProfile is responsible for -- P(r), rho(r), g(r),
m(r), M, C/MR^2 and core radius -- for the architecture

    surface -> observational crust -> Moho -> stx24 mantle -> CMB -> Fe-S core

Key points, and what is deliberately NOT settled here:

* Mantle EOS is the unified Warren (2005) Preferred Ti-free + stx24 table
  (quartz enabled, 1000-50000 bar, 300-2200 K, 300x300, resolution-converged).
  Swapping to Hauri 2015 or Warren Alt 5 requires changing only
  Sil.mantleEOS -- no other structural parameter.
* The observational crust (Garcia 2019 M3) is separate from, and must not be
  confused with, the mantle EOS. It is a structure-only layer: it contributes
  to mass, pressure, gravity and moment of inertia, but its seismic properties
  are deliberately unset (NaN) until a later stage connects Garcia Vp/Vs/Ks/G.
  The Perple_X thermodynamic highlands crust is a validation diagnostic only
  and is NOT used as the crust property source.
* Thermal structure will be supplied externally in Stage 2
  (PRESCRIBED_SILICATE_T). Until then the legacy conductive shell path runs,
  and NO temperature-dependent structural output from this configuration
  should be interpreted scientifically. Tsurf/qSurf/Qrad below exist to let
  the code initialize and are labeled PROVISIONAL.
* The core is PROVISIONAL_CORE_CONFIGURATION_PENDING_FINAL_SELECTION. The
  Fe-S EOS product itself is repository-traceable and adopted, but the wider
  core configuration (S content beyond wFe, radius) is not scientifically
  frozen. Core radius must be solved structurally, never tuned by hand.
* LOV3D is not part of this configuration.

Parameter provenance is recorded per line below and in
lunar-bsm-interior-models/docs/LUNAR_PLANETPROFILE_STAGE0_STAGE1.md.
Labels: [OBS] observational constraint, [ASSUM] model assumption,
[PROV] provisional input, [DERIV] derived value.
"""
import numpy as np
from PlanetProfile.Utilities.defineStructs import PlanetStruct, Constants

Planet = PlanetStruct('Luna')

""" Bulk planetary settings """
Planet.Do.NO_H2O = True   # [ASSUM] Dry body: no hydrosphere.
Planet.Do.Fe_CORE = True  # [ASSUM] Fe-bearing core present (Weber et al. 2011; Williams et al. 2014).

# [ASSUM] Strict EOS-domain policy, set explicitly here rather than inherited
# from Params: a scientific lunar run must FAIL if any mantle or core EOS
# evaluation falls outside the tabulated Perple_X domain, instead of silently
# saturating at the table edge. The lunar central pressure sits close to the
# unified Warren table's 5 GPa ceiling, so this is not a hypothetical risk.
# Strict EOS-domain policy: a scientific lunar run must FAIL if any
# PHYSICALLY ACTIVE mantle or core node falls outside the tabulated Perple_X
# domain, rather than silently saturating at the table edge. The lunar central
# pressure sits close to the unified Warren table's 5 GPa ceiling, so this is
# not hypothetical.
#
# Strictness applies to physical nodes only. PlanetProfile propagates every
# candidate silicate profile down toward r = 0 as scratch storage for the
# core-radius search; those unused nodes routinely exceed the ceiling. They are
# clamped during propagation (they cannot contaminate active nodes, since the
# recursion is strictly downward and active nodes are always shallower) and the
# policy is then enforced on each candidate's active extent once nSilFinal is
# known. See docs/LUNAR_PLANETPROFILE_WARREN_NUMERICAL_CLOSURE.md.
Planet.Do.EOS_DOMAIN_POLICY = 'error'

# [OBS] Mean radius 1737.15 km. Neumann, G. A. (2013), LRO LOLA topographic
# data products, as adopted in Williams et al. (2014), JGR Planets 119,
# 1546-1578, Table 2. (Smith et al. 2010, GRL 37, L18204 give a closely
# similar but distinct 1737.153 +/- 0.010 km.)
Planet.Bulk.R_m = 1737.15e3

# [DERIV] M = GM/G = 7.3458e22 kg, from GM = 4902.80007 +/- 0.00014 km^3/s^2
# (Williams, Boggs & Folkner 2013, JPL IOM 335-JW,DB,WF-20130722-016, Table 5;
# the DE430 ephemeris value adopted in Williams et al. 2014 Table 1) and
# G = 6.6743e-11 m^3/kg/s^2 (CODATA 2018). Quoted to 5 significant figures --
# the precision the division supports given G's uncertainty, not GM's.
Planet.Bulk.M_kg = 7.3458e22

# [OBS] Whole-Moon MEAN moment-of-inertia factor I/MR^2 = 0.3931 +/- 0.0002.
# Konopliv, A. S., Binder, A. B., Hood, L. L., et al. (1998), Science 281,
# 1476-1480 ("the average moment I/MR^2 = 0.3931 +/- 0.0002").
#
# NAMING WARNING (Stage E2a/E2a.1, 2026-09-01). The attribute is called
# `Cmeasured` for historical/API reasons, but the value stored here is the
# MEAN moment I = (A+B+C)/3, NOT the polar moment C. Konopliv et al. give both
# in the same sentence: polar C/MR^2 = 0.3932 and mean I/MR^2 = 0.3931. The mean
# moment is the physically correct comparator, because a spherically symmetric
# model has A = B = C and cannot predict C separately from I; using the polar
# value would introduce a (2/3)J2 ~ 1.35e-4 bias. The value is right; only the
# attribute name is misleading.
# TODO(deprecation): rename to something like Bulk.meanMomentMeasured. NOT done
# here -- `Cmeasured` appears ~135 times across PlanetProfile, including the app,
# the test suite and many body configurations, so the rename is a separate
# authorized refactor, not a cosmetic edit.
#
# UNCERTAINTY WARNING. The quoted +/-0.0002 is, verbatim in the same sentence,
# "an uncertainty of five times the formal error". It is NOT a Gaussian 1-sigma.
# Downstream code that forms (predicted - Cmeasured)/Cuncertainty is computing a
# LEGACY-BAND UNIT, not a number of standard deviations.
#
# This is the whole-body value, which is the physically correct target for
# PlanetProfile's shell sum over every layer, fluid core included. Do NOT
# substitute Williams et al. (2014)'s solid-Moon mean moment (0.392728 at
# R = 1738.0 km, or 0.393112 at R = 1737.151 km), which excludes the fluid
# core's own contribution and is not comparable to a whole-body shell sum. The
# modern whole-Moon comparator is a MODEL-DEPENDENT RANGE, I_m/MR^2 in
# [0.39311, 0.39338] at R = 1737.151 km (Williams et al. 2014, Table 8).
Planet.Bulk.Cmeasured = 0.3931
Planet.Bulk.Cuncertainty = 0.0002

# [PROV] Surface temperature 250 K, the value adopted for the crustal
# conductive geotherm in Laneuville et al. (2013), DOI 10.1002/jgre.20103,
# and used throughout this project's crust thermal work. Provisional: it is
# a modelling convention for the geotherm, not a measured global mean, and
# becomes an input to the prescribed T(r) in Stage 2.
Planet.Bulk.Tsurf_K = 250
Planet.Bulk.Psurf_MPa = 0.0  # [ASSUM] No atmosphere.

# [PROV] Global mean surface heat flow ~12 mW/m^2, Warren & Rasmussen (1987),
# DOI 10.1029/JB092iB05p03453 (best global-mean estimate after correcting the
# Apollo 17 boundary setting). Provisional and largely superseded once
# Stage 2 supplies T(r) externally: with a prescribed thermal profile this
# value no longer sets the mantle temperature structure.
Planet.Bulk.qSurf_Wm2 = 12e-3

""" Layer step settings """
# [DERIV] Production structural resolution, set by a convergence study, not by
# convenience. PlanetProfile's silicate recursion is first-order in the radial
# step h, verified: successive refinements 300 -> 600 -> 1200 -> 2400 change
# R_core by 2.895, 1.448, 0.724 km -- ratios of exactly 2.00. The pre-registered
# convergence criteria are met only at 1200 (checked against 2400); 300 and 600
# both fail on the pressure profile and 300 also on mass residual.
# See docs/LUNAR_PLANETPROFILE_WARREN_NUMERICAL_CLOSURE.md.
#
# NOTE for composition comparisons: |R_core - R_core_converged| equals exactly
# one radial interval at every resolution, so R_core is resolved only to h.
# At 1200 that is 1.448 km. A Hauri-vs-Warren core-radius difference smaller
# than that is NOT resolved and must not be interpreted.
Planet.Steps.nSilMax = 1200
Planet.Steps.nCore = 10     # [ASSUM] Core radial resolution.

""" Observational crust (Garcia 2019 M3 nominal; structure-only) """
# Nominal crust for this project is Garcia et al. (2019), Space Sci. Rev.,
# Table 11, model M3 -- selected for global representativeness and its native
# 40 km Moho. Garcia 2019 M2 is retained as a sensitivity case (not wired
# here). See lunar-bsm-interior-models/docs/MOHO_PHASE_SPACE_AND_CRUST_REFERENCE_CLOSURE.md.
Planet.Crust.enabled = True
Planet.Crust.structure_only = True
Planet.Crust.bulk_density_includes_porosity = True
Planet.Crust.radiogenic_mode = 'fixed_specific_rate'
Planet.Crust.material_id = Constants.materialCrust

# [OBS] 40 km mean crustal thickness; GRAIL range 34-43 km, Wieczorek et al.
# (2013), Science 339, 671-675. 40 km is M3's own native Moho depth and this
# project's adopted nominal value.
Planet.Crust.thickness_m = 40.0e3

# [DERIV] 2751.4 kg/m^3: spherical-shell-volume-weighted mean of Garcia 2019
# M3's own sampled densities (2600 at 0 km, 2760 at 20 km, 2890 at 40 km),
# linearly interpolated between those points over 0-40 km. Reproduces the
# exact integrated M3 crust mass (4.078079e21 kg) to 1e-15 relative, so using
# a single bulk density here loses no mass relative to the layered profile.
# This is an effective (porosity-inclusive) seismic density, NOT a
# Perple_X solid-matrix grain density.
Planet.Crust.bulk_density_kgm3 = 2751.4

# [PROV] 1.5 W/m/K effective crustal conductivity, Laneuville et al. (2013)
# (~5 km of 0.3 W/m/K megaregolith over ~35 km of 3 W/m/K crust). A model
# assumption, not a measurement; irrelevant to structure once Stage 2
# prescribes T(r).
Planet.Crust.thermal_conductivity_WmK = 1.5

""" Prescribed silicate temperature (Stage 2) """
# PlanetProfile's own conductive silicate propagation is not validated for the
# Moon (the legacy shell-conduction branch is retained for Fe_CORE bodies), so
# the silicate temperature is supplied externally instead. Only the temperature
# assignment is replaced: P, rho, g, layer mass, MoI and the core-radius search
# are still solved by PlanetProfile, with the EOS evaluated at this T.
#
# The resulting model is hydrostatically self-consistent CONDITIONAL ON this
# prescribed thermal field. It is NOT thermally self-consistent.
#
# Source: Santangelo et al. (2025) preferred model BU16_C_th04_K_027_r15_th6,
# area-weighted global mean over the native 180x360 lateral grid -- the same
# reduction used for every radial validation in this study. Extraction and
# provenance: lunar-bsm-interior-models/inputs/thermal/
# santangelo_preferred_radial_temperature.csv
#
# Coverage: 0-1350 km depth. It spans the crust (surface value 250.00 K,
# independently consistent with Bulk.Tsurf_K above), so no separate crustal
# thermal bridge is required. It does NOT reach a plausible CMB; any mantle
# node deeper than 1350 km has no prescribed temperature and will raise rather
# than extrapolate.
Planet.Do.PRESCRIBED_SILICATE_T = True

from pathlib import Path as _Path
from PlanetProfile.Thermodynamics.PrescribedT import LoadPrescribedTProfileCSV as _LoadT
_TPROFILE = _Path('/Users/emma/Documents/Code/lunar-bsm-interior-models/inputs/thermal/'
                  'santangelo_preferred_radial_temperature.csv')
# deepContinuation='isothermal' is REQUIRED here, and is not a scientific
# choice about the deep Moon. PlanetProfile propagates every candidate
# silicate profile from its start radius down to r = 0 as scratch space for
# the core-radius search (Silicates.py sets rSilEnd_m = 0 when Fe_CORE is
# set), so it queries T at depths no finite profile can cover. Holding the
# deepest tabulated temperature constant below 1350 km lets that scratch
# space evaluate; every such use is COUNTED on the profile object, and the
# run is only interpretable if the FINAL selected model is afterwards shown
# not to depend on it (i.e. its deepest mantle node lies above 1350 km).
# The tail is in any case nearly isothermal: 1709.6 -> 1713.9 K over the
# last 24 km, a gradient of ~0.06 K/km.
Planet.Sil.T_of_z_K = _LoadT(_TPROFILE,
                             source='Santangelo et al. (2025) BU16_C_th04_K_027_r15_th6, area-weighted global mean',
                             deepContinuation='isothermal')
Planet.Sil.TprescribedSource = str(_TPROFILE)

""" Silicate mantle """
# [ASSUM] Unified Warren (2005) Preferred Ti-free + stx24 table: quartz
# enabled, P 1000-50000 bar, T 300-2200 K, 300x300, verified
# resolution-converged (lunar-bsm-interior-models/docs/UNIFIED_STX24_RESOLUTION_CONVERGENCE.md).
# Change ONLY this line to swap composition to Hauri 2015 or Warren Alt 5.
Planet.Sil.mantleEOS = 'moon_w05_preferred_stx24_unified.tab'

# [PROV] Bulk-silicate radiogenic heating. Placeholder pending the Stage-2
# prescribed thermal profile, which supersedes it as the control on mantle
# temperature. Not a sourced lunar value and must not be interpreted as one.
Planet.Sil.Qrad_Wkg = 0.0
Planet.Sil.Htidal_Wm3 = 0.0  # [ASSUM] Present-day lunar tidal heating negligible for structure.

# Porosity is not modeled: the crust's effective density already includes it,
# and the mantle EOS is a solid-matrix calculation.
Planet.Do.POROUS_ROCK = False

""" Core -- PROVISIONAL_CORE_CONFIGURATION_PENDING_FINAL_SELECTION """
# The EOS product is repository-traceable and adopted (Fe-5 wt% S, wFe_ppt=950)
# -- see lunar-bsm-interior-models/docs/FE_S_CORE_EOS_PRODUCT.md and the
# decision to reject the shipped Fe-S EOS in favour of a self-built traceable
# product. The wider core configuration is NOT frozen: sulfur content beyond
# wFe, and the core radius, remain open. Core radius must be SOLVED against
# mass and C/MR^2, never tuned by hand.
Planet.Core.coreEOS = 'Fe-S_3D_EOS_5slice_traceable_v1.mat'
Planet.Core.wFe_ppt = 950  # [PROV] Fe-5 wt% S, the adopted lunar core composition slice.
Planet.Core.wS_ppt = 50    # [PROV] Complement of wFe_ppt.

# [PROV] Endmember densities and minimum core density, inherited PlanetProfile
# defaults retained so the mass/MoI search has bounds. Not lunar-specific
# determinations; they set the maximum admissible core size only.
Planet.Core.rhoFe_kgm3 = 8000.0
Planet.Core.rhoFeS_kgm3 = 5150.0
Planet.Core.rhoMin_kgm3 = 4500.0

# [PROV] Core seismic quality factor, PlanetProfile convention. Not used
# structurally; retained for downstream seismic bookkeeping.
Planet.Core.QScore = 1e4

""" Seismic properties of solids """
Planet.Seismic.lowQDiv = 1.0  # [ASSUM] PlanetProfile default.
