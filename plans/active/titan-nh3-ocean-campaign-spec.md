# Titan NH3-ocean campaign spec (Machine B)

User-directed 2026-07-25. The NH3 blocker is resolved: CoolProp
NH3-H2O ocean EOS is live on genai (Ocean.comp='NH3'), with the
mu-based liquidus (chemical-potential equality, SeaFreeze ice G vs
CoolProp water activity) as the default melt mode — validated against
the ideal colligative dilute limit and working under all ice phases
Ih–VI (tests/coolprop_nh3_test.py, 9 tests). User ratification: the
LOW MoI of the achievable Titan interiors is acceptable — this
campaign uses the free-gravity C20/C22 configuration (CMR2 dropped as
observable, mass-conservation rho_sil), so classic MoI matching is
NOT required and its current failure mode for Titan does not block.

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
- Scan Tb in ~[240, 263] K x w in {0, 10, 30, 50, 100, 150} ppt.
  Expected physics: mu-mode depression is ~1.1 K at 10 ppt, ~5.4 K at
  50 ppt, ~15 K at 150 ppt (larger under HP ices) — the valid
  ocean-bearing Tb band shifts DOWN with w accordingly. Record which
  nodes produce ocean + ice-I shell + (any) HP ices, and D_iceIh /
  D_ocean per node.
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

- CoolProp mixture uses an ESTIMATED reducing rule (no fitted NH3-H2O
  binary pair exists in CoolProp 8.x): rho ~1%, Cp 2-5% vs Melinder;
  documented in Thermodynamics/NH3/NH3Props.py.
- L-K (2002) polynomial retained only as fallback; it underestimates
  the depression ~2x at low-mid X (flagged in-code; verify original
  coefficients if it is ever promoted back).
- Ices are pure H2O (NH3 fully rejected to liquid); no NH3 hydrates
  (ammonia dihydrate etc.) in the phase diagram — acceptable above
  ~176 K and w below the eutectic, state it.
