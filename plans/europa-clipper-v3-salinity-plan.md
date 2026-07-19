# Europa Clipper v3: ocean salinity as an 8th sampled parameter (Machine B)

Author: Machine A, 2026-07-14 (user-directed: "train for a range of
salinities from 0.1 to 100 ppt (range widened from 1 ppt, user 2026-07-18)", seawater composition retained). Machine B
runs everything compute-intensive.

Status 2026-07-19: Phases 1–3 `verified` (Machine B, commits bfc7f3f1 →
f228810f; gates in validation_reports/europa_clipper_v3_1m/, artifact
RATIFIABLE). Phase 3.5 cloud plot `verified` (Machine A, AppTest + PNG
against the committed reference pickle). Phase 4 GUI slot `verified`
end-to-end: 2D cache committed (Machine B, f7a03572), then Machine A ran
the full GUI sampling path via AppTest — 10k draws/4 min, posterior
matches the reference MCMC (Tb 264.46 vs 264.3 K, w 36.6 vs 38.7 ppt,
corr(Tb, log10 w) −0.988 vs −0.986), per-sample cloud plot rendered.
ONLY user ratification remains; HF deploy after that.

## Objective

The v2 Clipper artifact fixes ocean salinity at the Test51 default — a
known oversight: salinity drives ocean conductivity (hence all 14 Bind
channels), density, and thicknesses, so v2 posteriors implicitly condition
on one assumed ocean. v3 samples it:

- New parameter `log10_wOcean_ppt`, prior uniform [-1, 2] (= 0.1-100 ppt,
  log-uniform per scientific review: a 2-decade scale parameter with a
  linear-uniform prior would put ~90% of mass in 10-100 ppt and drag the
  posterior high wherever the data are weak; Jeffreys-style log-uniform is
  the honest weak prior). Composition `Seawater` (GSW) throughout.
  Pre-register dual-space posterior reporting (w and log10 w) plus a
  prior-sensitivity check regardless.
- Everything else identical to `europa_seawater_andrade_clipper_v2.json`
  (7 params -> 8; same 19 observables incl. 14 Bind channels at
  sigma = 1.5 nT; same Galileo synodic support cut; same
  derived_params / priors).

**Scientific caveat (config metadata, quantified):** GSW/PSS-78 is
calibrated to SP <= 42; 42-100 ppt is an extrapolation (verified monotone
and plausible in magnitude, but ~O(10-20%) systematic vs hypersaline brine
data — a sigma systematic NOT inside the 1.5 nT noise budget; there is no
wmax guard on the GSW path). Phase 1 includes a HARD numeric gate: compare
sigma(w), rho(w), Tfreeze(w) at w in {50, 70, 100} against an independent
hypersaline reference (published artificial-seawater / NaCl brine data),
record divergences; do NOT substitute GSW self-monotonicity (trivially
true) for this check.

## Phase 1 — 2D structure grid (cache rebuild)

The current cache is 1D in Tb (`Tb_K_grid` + `structures[i]`). Salinity
changes the hydrosphere structure, so the cache becomes 2D:

1. Builder reality check (review 2026-07-14): the existing API is
   `build_tb_grid_cache` / `build_single_structure` (NOT the plan's earlier
   `build_phase_c1_cache` name), and it writes a 1D schema with
   TRANSITION-REFINED (irregular) Tb nodes — per-w refinement would yield
   different Tb nodes per salinity, which CANNOT be stacked row-major.
   **Decision (binding): fixed REGULAR dense Tb grid with
   `detect_transitions=False`** — e.g. 0.125 K spacing over the Tb box
   (~92 nodes) x ~16 log-spaced w — and a verification that CMR2/D_ocean
   near ice transitions do not degrade vs the refined v2 cache at w=35
   (that refinement existed for a reason; quantify what a regular grid
   loses before committing). A new 2D builder function is a real rewrite,
   not a loop — budget for it.
2. Schema: `Tb_K_grid` + `wOcean_ppt_grid`; `structures` row-major
   [i_Tb * n_w + i_w] with BOTH `Tb_K` and `wOcean_ppt` stored per entry;
   bump cache `schema_version`; absent `wOcean_ppt_grid` -> legacy 1D path
   so v1/v2 artifacts stay servable.
3. Freezing point depends on salinity: at high w the valid-Tb window
   SHIFTS (~6-7 K lower freezing at 100 ppt; higher at 1 ppt), so a fixed
   Tb box is mostly-invalid in one corner (range now reaches 0.1 ppt, nearly-fresh water: freezing point ~pure-ice; conductivity ~0.01 S/m -> the |Ae_synodic|>0.7 support cut will carve away much of the low-w plane; expect large rejected regions, report them). Store failures as explicit
   `None` entries (rejected at sampling, never silently skipped) AND
   report D_ocean/D_iceIh maps over (Tb, w) confirming the box still
   resolves the valid region at w = 1 and w = 100; widen the Tb box if
   not.
4. w-axis resolution near the physics: Ae moves fastest at LOW w (synodic
   de-saturation) and across the |Ae_synodic| = 0.7 support contour;
   nearest-neighbor over 16 log-w points staircases that boundary. Either
   bilinear interpolation in (Tb, log w) for Ae/structure scalars, or
   w-refinement near the 0.7 and no-ocean contours. Pre-registered
   convergence check: rebuild a sub-box at 2x w-resolution and confirm the
   support boundary and a fiducial posterior move by less than the gate
   tolerances.
5. Expected cost: ~92 x 16 = ~1500 PP runs with detect_transitions=False
   (cheaper per run than refined builds). Machine B.

## Phase 2 — code: 2D cache lookup (small)

`mcmc_runner` touchpoints (all currently nearest-Tb):
- `_get_cache_scalar` / `_expand_theta`-consumers: nearest node in
  (Tb, log w) — same nearest-neighbor convention as today, now 2D.
- `_precompute_ae_grid` + `_induction_channel_value`: Ae grid keyed by
  (i_Tb, i_w); conductivity now varies with BOTH.
- `induction_bounds` support cut: unchanged logic, evaluated on the 2D
  grid — the |Ae_synodic| > 0.7 cut now carves a (Tb, w) region (low
  salinity + thin ocean fails: physical, intended).
- `_compute_model_cmr2` mass-conservation path: hydrosphere columns come
  from the (Tb, w) node — no logic change beyond the lookup.
- Unit tests: 2D nearest lookup; legacy 1D cache still loads; support cut
  rejects a known low-w/low-Tb corner.

## Phase 3 — config + reference + training (Machine B)

1. `europa_seawater_andrade_clipper_v3_8D.json`: v2 copy + `wOcean_ppt`
   {uniform, [1, 100]} + new cache path + metadata (GSW extrapolation
   caveat; this plan). Fresh seeds: train 42 / data 48 / noise 4848.
2. Fiducial Bind central values: recompute at the fiducial model
   (Tb ~ 266 K, w = 35.165 ppt) from the NEW cache — do not reuse v2
   values blindly (the new grid's nearest node may differ slightly).
3. Reference MCMC on the v3 config (nautilus, as v2). COMMIT the reference
   pickle — v2's lived in /tmp on Machine B and is unavailable to Machine
   A. Target location Test/mcmc_results/Europa/ — USER PERMISSION GRANTED
   (2026-07-18); commit the reference pickle there.
4. 1M dataset (support guard ON, drop_nonfinite) -> nsf train -> artifact.
5. Gates (ratified rules): SBC (now 8 params), crosscheck, limits. Anchor
   design: (Tb, w) grid-walk, 4 Tb x 3 w on-manifold anchors PLUS at least
   one anchor PAIR on a single iso-Ae_synodic contour (same synodic Ae,
   different Tb/w) — a rectangular grid never tests the Tb<->w degeneracy
   ridge. Choose w-anchors wide enough that the ORBITAL Ae actually varies
   (at saturated synodic it is the only channel that moves with w). The v2
   lesson stands: anchors move all frequencies jointly from one physical
   (Tb, w).
6. **2D degeneracy gate (pre-registered, the central v3 deliverable):**
   the flow must reproduce the reference MCMC's JOINT Tb-w posterior at
   the fiducial — correlation coefficient + 2D credible-region match (2D
   W1 or equivalent), not just marginals: marginal SBC/W1 can all pass
   while the degeneracy tilt is wrong. Physics expectation: synodic Ae is
   saturated (0.75-0.94), so w is carried mostly by the orbital channels
   at |Bind| ~ 2-9 nT vs 1.5 nT noise — w will be weakly identified and
   Tb<->w degenerate (both raise ocean conductance). Also run the
   prior-sensitivity check from the Objective section; if the w posterior
   tracks the prior, say so plainly in the INDEX row.
7. Expect the Tb shape defect (right-skew vs near-Gaussian flow) to
   persist; report it the same way. Watch wOcean's SBC rank histogram for
   prior-dominance artifacts.

## MCMC side (user 2026-07-18)

The SAME sampled `log10_wOcean_ppt` (U[-1,2]) applies to Europa MCMC
configs once the 2D (Tb, w) cache + lookup land (Phases 1-2) — author a
`europa_seawater_andrade_mcmc_8D_salinity.json` alongside the SBI config
and register the parameter in the GUI registry (done Machine A-side:
parameter_registry log10_wOcean_ppt). Titan no-ocean configs are exempt
(no ocean to salinize).

## Phase 3.5 — GUI: k2/h2/Ae complex-plane CLOUD plot (absorbed from
Machine B's 2026-07-18 design spec; `not implemented`, Machine A)

Replaces the interim per-Tb-node "k2 complex plane by model" panel once a
salinity-sampled artifact exists. One connected path per posterior SAMPLE
(full interior model: rheology + Tb + composition) linking that sample's
dimensionless complex signals on a single Re-Im plane: k2, h2, Ae(synodic),
Ae(synodic 2nd), Ae(orbital) — 5 signals for 3-frequency runs, 3 for
Galileo. Color = salinity, alpha low (cloud form, user-confirmed). The
current plot's "Ae depends only on Tb" statement is a v2 fixed-seawater
limitation, NOT physics — must not survive into v3.

## Phase 4 — Machine A (after artifacts land)

GUI slot (config_path REQUIRED), derived-|Ae| guard recomputed for the
validated envelope from the v3 grid-walk, Tb default truncation if the
support edge persists, salinity posterior appears automatically in the
parameter panels. AppTest + user ratification before any deploy. User directive
(2026-07-16): this work returns to the MAIN repo (artifact + gates +
INDEX + slot on origin/genai) BEFORE any Hugging Face deployment — the
public app only ever consumes the repo (DEPLOYING.md release discipline).

## Non-goals

- Non-seawater compositions (MgSO4/NaCl salinity axes) — future work; the
  2D cache schema should not hardcode 'Seawater' anywhere new.
- Retrofit of v1/v2 artifacts: they stay as-is (fixed-salinity scope
  documented in their INDEX rows).
