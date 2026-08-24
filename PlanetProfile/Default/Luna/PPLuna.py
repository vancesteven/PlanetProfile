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

# [OBS] Whole-Moon mean moment-of-inertia factor I/MR^2 = 0.3931 +/- 0.0002.
# Konopliv, A. S., Binder, A. B., Hood, L. L., et al. (1998), Science 281,
# 1476-1480 ("the average moment I/MR^2 = 0.3931 +/- 0.0002").
# This is the whole-body value, which is the physically correct target for
# PlanetProfile's C/MR^2 -- computed as a plain shell sum over every layer,
# fluid core included. Do NOT substitute Williams et al. (2014)'s
# solid-Moon mean moment (0.392728), which excludes the fluid core's own
# contribution and is not comparable to a whole-body shell sum.
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
Planet.Steps.nSilMax = 300  # [ASSUM] Silicate radial resolution.
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
