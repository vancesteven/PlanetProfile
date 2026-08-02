# STATUS 2026-08-01 — Machine B → Machine A handoff (Titan joint NH3 SBI)

Author: Machine B (genai worktree). For: Machine A planning + code changes that
**do not collide** with Machine B's active work.

## Headline decision (user, firm 2026-07-30)

The Titan **NH3** free-gravity SBI artifact must produce a **JOINT no-ocean +
ocean posterior**. Frozen (Tb, w) grid nodes must build as REAL no-ocean
interiors (own k2 / C20 / C22 via the same per-sample forward model) and
**appear** in the posterior — NOT be rejected as out-of-support. Keep the FULL
Tb ∈ [249,263] K × w ∈ [1,70] NH3 ppt range. Do NOT truncate Tb to force an
all-ocean cache.

Rationale: the user wants to see what a distribution spanning both no-ocean and
ocean Titan interiors looks like. The current (ocean-only) build path silently
drops frozen models → posterior is implicitly conditioned on "an ocean exists."

## What Machine B is actively changing (DO NOT TOUCH these — merge conflicts)

Task #68 (blocks Phase B #51). Machine B is scoping/implementing the joint
build. Files Machine B owns for this work:

- `PlanetProfile/Inference/cache_builder.py` — adding a frozen-node retry via a
  new `do_overrides` channel (sets `Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES=True`
  on the deep-copied Planet) so a genuinely-frozen NH3 node builds as a real
  no-ocean structure instead of being stored as `None`. Node tagged
  `has_ocean=True/False` in the structure dict.
- The NH3 SBI config (`configs/test54_titan_nh3_freegrav.json`, to be created
  from `titan_freegrav_noocean.json`) — must NOT set
  `sampler_settings.phase_stability.enforce = "no_ocean_Ih"` (that guard would
  reject ocean OR no-ocean samples; the joint config keeps neither-direction
  enforcement). density_inversion_guard + k2_support_bounds still apply.
- The NH3 2D structure cache under
  `PlanetProfile/Test/mcmc_results/Titan/Test54_nh3_ocean/`.
- `plans/scripts/titanG_*` NH3 driver variants and any `/tmp` build scripts.

Reviewer (scientific-reviewer) sign-off: **PASS WITH CONCERNS** (2026-08-01),
guardrails folded into the implementation. No full 1M build launched yet.

**Code implemented in Machine B's working tree 2026-08-01 (dispatch verified,
real-PP smoke run) — NOT yet committed to genai (still uncommitted as of
2026-08-02; will be committed/pushed once the 3×4 validation confirms it):**
- `defineStructs.py`: new `NoIceLiquidTransitionError(ValueError)`.
- `LayerPropagators.py:116`: raises that typed subclass at the no-Ih→liquid
  branch (message unchanged; `except ValueError` callers unaffected).
- `cache_builder.py`: new `do_overrides` channel (`setattr(Planet.Do, …)`)
  threaded through all 4 builder fns, and `build_tbw_grid_cache` gains a
  `retry_frozen_as_no_ocean` flag (default **False** → every existing Europa
  v3/v5/v6/v7 2D cache behaves exactly as before). When True, a node raising the
  typed frozen exception is retried once as `NO_OCEAN_EXCEPT_INNER_ICES`; other
  failures still store `None`. Cache dict adds `retry_frozen_as_no_ocean` +
  `n_no_ocean_nodes`; each structure carries a `has_ocean` bool.

Machine A impact: **none forced** — the builder-API additions are backward
compatible (new kwargs default off). If Machine A reads 2D caches, be aware
structures may now carry `has_ocean` and the cache may carry the two new
metadata keys.

NOTE (verified this session): the four observables {C20, C22, Re_k2, Im_k2} are
**recomputed per sample** from the interpolated structure (k2 via TidalPy;
C20/C22 via the Clairaut integral in `mcmc_runner._derive_gravity_pair`), NOT
read from cached values. `grid_interp_2d` + `forward_models` already treat any
stored structure as valid support and already fall back to nearest-corner across
the ocean/no-ocean `region_phases` discontinuity. So **no change is needed to
`grid_interp_2d.py`, `forward_models.py`, or `mcmc_runner.py`** for the joint
posterior — Machine A may read them but Machine B does not expect to edit them.

## What Machine A can safely plan + change now (no collision)

**Phase C — GUI slot wiring (task #52).** Owned by Machine A per the campaign
spec. Target: `PlanetProfileApp/pages/Inference.py`.

1. Add four Titan slots to `_SBI_ARTIFACT_SLOTS` (~line 1049): no-ocean (Phase A,
   already shipped) + NH3 (joint) + MgSO4 + NaCl. Each: label /
   `bodyname='Titan'` / `config_path` / `cache_path` / `default_obs` (config
   centrals, **no CMR2**) / `x_obs_limits` (bound Im_k2) / `scope_note`.
2. **Shared-widget-key hazard (must fix):** the observable-input widgets use a
   global key `amort_obs_Re_k2` (~line 1740) shared with Clipper. With multiple
   Titan slots + Clipper all exposing Re_k2, keys MUST be namespaced per slot id
   (same fix applied for Clipper slots in commit ec711672). This is
   self-contained to `Inference.py` and does not touch Machine B's files —
   **safe to implement now.**
3. For the NH3 slot `scope_note`: state it is a JOINT no-ocean+ocean posterior
   (frozen Titan interiors included), gravity provenance Petricca 2025, CMR2
   dropped (C22 double-count), induction+h2 dropped (no clean Cassini signal),
   NH3 ocean w∈[1,70] ppt. Machine B will supply the final cache_path,
   config_path, gate status, and confirmed default_obs centrals when the NH3
   artifact is built — leave those as clearly-marked TODO placeholders for now.

Recommend Machine A: implement the key-namespacing fix (#2) and stub the slot
scaffolding (#1, #3) with placeholders, verify in-browser that existing slots
still load and widget state no longer collides. Do NOT hardcode NH3 artifact
paths yet — Machine B will hand those over.

## Status vocabulary reminder

`verified` / `implemented, unverified` / `not implemented` only. UI changes are
`implemented, unverified` until observed in-browser (screenshot / `/verify`).

---

## ADDENDUM 2026-08-02 — Machine B → A: queue reconciliation needed (Titan NH3)

Machine B (Claude Opus 4.8, PPcl env) read the coordination refresh in
`9cfc3c14` (new roster, `STATUS.md`, `MACHINE-B-HANDOFF.md`). Thanks — the
lane structure is clear and the backward-compatible builder API note above
still holds. One reconciliation item before Machine B proceeds:

**`MACHINE-B-HANDOFF.md` §1 and this joint-NH3 decision disagree, and the
handoff does not reference the joint decision at all.** Side by side:

| | This file / Task #68 (user firm 2026-07-30, **re-affirmed 2026-08-02**) | `MACHINE-B-HANDOFF.md` §1 (2026-08-01/02) |
|---|---|---|
| Goal | JOINT no-ocean + ocean posterior; frozen nodes build as real no-ocean interiors and appear in the posterior | "Titan NH3 Phase 0 RE-RUN" — feasibility scan → 2D cache → SBI |
| Tb × w | FULL Tb ∈ [249,263] K × w ∈ **[1,70]** ppt (spans the frozen wedge) | provisional rectangle Tb ∈ [248,257] K × w ∈ **[30,100]** ppt |
| Frozen interiors | kept in posterior (the whole point) | not mentioned; the rectangle is chosen to sit in ocean-favorable T,w and avoid the frozen wedge |

Two concrete concerns with the §1 rectangle as written:
1. It would **re-condition the posterior on "an ocean exists"** — the exact
   outcome the user asked us to eliminate. Frozen Titan interiors would again
   be dropped.
2. **w = 100 ppt exceeds the NH3 physical ceiling.** NH3 7 wt% = 70 ppt
   (mass fraction of solution); 100 ppt is off the composition range the user
   set for NH3 (see the campaign spec's per-composition salinity table).

**User adjudication (2026-08-02):** Task #68 (this joint posterior, full range,
frozen nodes kept) **governs**. NOT the §1 rectangle.

**Requested of Machine A (queue owner):** please reconcile
`MACHINE-B-HANDOFF.md` §1 to the joint-posterior design — i.e. replace the
"Phase 0 rectangle re-run" framing with: joint no-ocean+ocean cache over the
full Tb ∈ [249,263] K × w ∈ [1,70] ppt range via
`build_tbw_grid_cache(..., retry_frozen_as_no_ocean=True)`, on the corrected
Melinder-anchored NH3 model (`6c5ee2af`, which local HEAD includes). If Machine
A believes the corrected-model re-run should instead supersede the joint
decision, please say so here and we'll take it back to the user — but Machine B
is proceeding on #68 until then.

**Machine B state right now (no full 1M launched):**
- Code for #68 is implemented in Machine B's **working tree only — NOT yet
  committed to genai** (see the section above): typed
  `NoIceLiquidTransitionError`, `do_overrides` channel, `retry_frozen_as_no_ocean`
  flag (default off → all Europa caches unchanged). It will be committed and
  pushed once the 3×4 validation confirms it, so Machine A will not see these
  three edited files (`cache_builder.py`, `LayerPropagators.py`,
  `defineStructs.py`) on origin until then.
- Per-node behavior **verified** on real PP: frozen (249,1) NH3 node raises the
  typed exception → retries → builds a real no-ocean interior (D_ocean=0, ice
  Ih/III/V/VI) → forward models give finite Re_k2=0.113, |Im_k2|=0.069,
  C22=8.42e-6, C20=−28.06e-6, and the support guard **admits** it
  (`imag_convention='abs'` globally). Cited: `/tmp/nh3_frozen_reachability.json`.
- **Running now (~5 h, serial):** a 3×4 joint validation cache
  (Tb {249,251,253} × w {1,10,30,60}) exercising reviewer validations 1
  (frozen-band confinement), 2 (w-invariance of frozen x), 3 (boundary-placement
  interval), 5 (region_phases disjointness). Result → `/tmp/nh3_joint_3x4_validate.json`.
- After the 3×4 passes + scientific-reviewer interprets: author
  `configs/test54_titan_nh3_freegrav.json`, then get **user go-ahead** before the
  expensive full-range production cache + 1M NH3 build.

Nothing here changes the "Machine A can safely plan Phase-C GUI wiring now"
guidance above — that remains collision-free.
