# Titan NH3-ocean campaign spec (Machine B)

User-directed 2026-07-25. The NH3 blocker is resolved: CoolProp
NH3-H2O ocean EOS is live on genai (Ocean.comp='NH3'), with the
mu-based liquidus (chemical-potential equality vs SeaFreeze ice G) as
the default melt mode, working under all ice phases Ih–VI
(tests/coolprop_nh3_test.py). User ratification: the
LOW MoI of the achievable Titan interiors is acceptable — this
campaign uses the free-gravity C20/C22 configuration (CMR2 dropped as
observable, mass-conservation rho_sil), so classic MoI matching is
NOT required and its current failure mode for Titan does not block.

## SUPERSEDED IN PART 2026-08-02 — JOINT no-ocean+ocean posterior governs

User decision (firm 2026-07-30, re-affirmed 2026-08-02; Machine A reconciled
2026-08-02): the NH3 artifact is a JOINT no-ocean + ocean posterior over the
FULL Tb in [249, 263] K x w in [1, 70] ppt range. Frozen nodes build as real
no-ocean interiors (`retry_frozen_as_no_ocean=True`, `has_ocean=False`) and
remain in support. The "provisional replacement rectangle" below
([248,257] x [30,100]) is RETIRED — it conditioned on "an ocean exists" and
exceeded the 70 ppt NH3 ceiling. `test_campaign_rectangle_buildable` remains
in the test suite as a liquidus-solvability guard only, not campaign policy.
Authoritative execution state: `plans/STATUS-2026-08-01-machineB-joint-nh3.md`
and `plans/MACHINE-B-HANDOFF.md` §1. The activity-model correction below
still applies in full.

## IMPORTANT — activity model corrected 2026-07-28 (scientific review)

The water-activity model inside the mu-equality liquidus changed from
the CoolProp estimated mixture to a Melinder-anchored Redlich-Kister
excess on the Choukroun & Grasset (2010) P,T shape
(NH3_ACTIVITY_MODE='melinder-CG' in Thermodynamics/NH3/NH3Props.py).
Reasons: the CoolProp mixture's excess term has the WRONG SIGN
(gamma_w > 1; NH3-H2O requires < 1), under-depressing the liquidus by
9-37% over 30-150 ppt (3.5 K too warm at w=100), and its flash fails
across much of the 200-400 MPa HP-ice band. The anchored model is
analytic (no flash) and reproduces the Melinder (2010) experimental
freezing curve to < 0.05 K over X in [0.02, 0.175].

Consequences for this campaign:
- Corrected 1-bar depressions: ~1.1 K at 10 ppt, 3.5 K at 30, 6.1 K
  at 50, 13.6 K at 100, 23.2 K at 150 ppt (x ~1.34 by 1 GPa). The old
  "expected physics" numbers below CoolProp-era values are superseded.
- The previously proposed Phase-1 rectangle Tb in [241,259] K x
  w in [30,100] ppt is NOT buildable: the ice-Ih branch bottoms out
  at ~247.8 K for w=30 (pure-water Ih/III triple limit), and
  (Tb=259, w=100) leaves a ~6 km shell (degenerate).
- PROVISIONAL replacement: Tb in [248, 257] K x w in [30, 100] ppt
  (corner D_iceIh ~162/106/95/25 km). Guarded by
  tests/coolprop_nh3_test.py::test_campaign_rectangle_buildable.
- Phase 0 MUST BE RE-RUN under the corrected model before committing
  the Phase-1 cache — the earlier scan was computed on the defective
  liquidus and would bake a 24-32 km systematic D_iceIh bias into the
  SBI training set. Check region_phases homogeneity at the new
  corners: (248, 30) sits near the Ih/III triple (PbI ~203 MPa; ice
  III appears immediately) while (257, 100) has no HP ice.

## Machine B prerequisites

1. `pip install CoolProp` in PPcl (CheckCompat pins >= 8.0.0).
2. Pull genai (>= commit with Thermodynamics/NH3/). SeaFreeze must
   already be 1.1.3.
3. Sanity: `python -m pytest tests/coolprop_nh3_test.py` (9 pass,
   ~1 min; the mu-liquidus solves make EOS construction ~1 min per
   (comp, w) label, cached in EOSlist thereafter).

## Phase 0 — feasibility scan (cheap, do first)

Map the buildable (Tb, w_NH3) domain for Titan WITH an ocean:
- Body: Titan defaults + Ocean.comp='NH3'; interior per the test52
  differentiated-Titan recipe (the same reduced-stack/no-MoI-matching
  path that built titan_diff_noocean_structure_grid.pkl — do NOT use
  the classic FindInnerWithMoIAndEOS route, which currently fails for
  Titan at any composition with achievable C/MR2 ~0.317 vs 0.341).
- Scan Tb in ~[244, 263] K x w in {0, 10, 30, 50, 100, 150} ppt.
  Expected physics (Melinder-anchored model, 2026-07-28): depression
  ~1.1 K at 10 ppt, ~6.1 K at 50 ppt, ~23 K at 150 ppt at low P,
  growing ~1.2-1.35x under the HP ices — the valid ocean-bearing Tb
  band shifts DOWN with w accordingly, and the ice-Ih branch bottoms
  out at the pure-water Ih/III triple (Tb_min ~247.8 K at w=30).
  Record which nodes produce ocean + ice-I shell + (any) HP ices, and
  D_iceIh / D_ocean per node.
- Deliverable: the (Tb, w) rectangle + node spacing for Phase 1 (aim
  ~12 Tb x 6 w built nodes, transitions pinned as in cache_builder).

## Phase 1 — 2D structure cache (Tb x w_NH3)

- Same v3.0 2D cache format as the Europa (Tb x wOcean_ppt) caches so
  ALL existing interpolation machinery (forward_models 2D blend,
  grid_interp_2d, wOcean_ppt_from_theta with log10_wOcean_ppt or
  wOcean_ppt) works unchanged — w axis is NH3 ppt instead of seawater
  ppt; note that in cache metadata.
- No Ae sidecar needed: NH3(aq) is a nonelectrolyte (sigma placeholder
  ~1e-5 S/m) — no induction channel for this campaign (document as an
  explicit assumption; a conductive-impurity extension is future work).
- Path suggestion:
  PlanetProfile/Test/mcmc_results/Titan/Test54_nh3_ocean/
  titan_nh3_structure_grid_2d.pkl (+ offsets sidecar if the reduced
  stack needs the CMR2 discretization anchor like Test52 did).

## Phase 2 — config + reference MCMC

- Config: copy titan_freegrav_noocean.json →
  test54_titan_nh3_freegrav.json. Changes:
  - structure_cache_path → the Phase 1 cache.
  - param_space: keep the freegrav 12D pattern (Andrade rheology
    params, per-phase log10_eta, Tb_K OR thick-ice reparam, R_core_km,
    rho_core_kgm3, dC20_nh, dC22_nh) PLUS the salinity axis:
    log10_wOcean_ppt over the built w range (or linear wOcean_ppt if
    the grid is linear — match the cache axis).
  - Observables: unchanged from freegrav (k2 + free C20/C22, CMR2
    dropped). No induction_bounds (no channel). NO amp_min-style
    support cuts anywhere (per the open-support policy ratified for
    Europa v7 — Cochrane2025 posture).
- Reference MCMC on the new cache for the crosscheck gate.

## Phase 3 — SBI artifact + gates

- Standard recipe (1m posterior), SBC + crosscheck + limits; grid-walk
  anchors must span the full built (Tb, w) rectangle including the
  near-frozen edge (thin/no ocean) and the saltiest column.
- Ship: cache + reference result + artifact + gate reports + INDEX.md
  entry. Machine A then wires the GUI slot (appears under the
  Cassini–Titan analysis as the newest version).

## Known caveats to carry in metadata

- BULK properties (rho/Cp/alpha/sound speed) use the CoolProp
  ESTIMATED mixture (no fitted NH3-H2O binary pair in CoolProp 8.x):
  rho ~1%, Cp 2-5% vs Melinder. The LIQUIDUS uses the separate
  Melinder-anchored activity model (see above) — an intentional
  split, same as 2018 PP (error budget: a 5% Cp error moves the
  adiabat ~1.2 K over 800 MPa, smaller than the 3.5 K liquidus bias
  the split removes).
- The "L-K (2002)" polynomial (53.8X + 650X^3) has been DELETED: its
  attribution could not be verified anywhere and its dilute-limit
  slope is provably ~2x too shallow (thermodynamic requirement
  109 K/mass-fraction; it gave 53.8). Do not resurrect it.
- Ices are pure H2O (NH3 fully rejected to liquid); no NH3 hydrates
  (ammonia dihydrate etc.) in the phase diagram — acceptable above
  ~176 K and w below the eutectic (X ~ 0.33), state it.
- Liquidus temperatures below SeaFreeze water1 validity (239.6 K) are
  clamped to the floor (affects w >= ~100 ppt near the Ih/III corner
  only; no EOS grid samples below the floor).
