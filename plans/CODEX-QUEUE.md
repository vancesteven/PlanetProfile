# Codex 5.6 queue (Machine A delegate lane)

Updated: 2026-08-05 at genai `0ccb3524`. Curated by the Claude model manager
(Fable 5). Codex 5.6 works ONLY tasks listed here unless the user directs
otherwise. Instructions: `AGENTS.md` (binding), which itself binds `CLAUDE.md`.

## Protocol

1. Read `plans/STATUS.md`, then this file. Claim a task by setting its status
   to `in progress (Codex)` and commit that edit.
2. Work on the `genai` branch. Commit locally with clear messages. Do NOT
   push — the manager reviews and pushes. Note commit hashes in your report.
3. Report inline under the task using the CLAUDE.md status vocabulary:
   `verified` (cite artifact: test name, PDF path, AppTest), `implemented,
   unverified`, or `not implemented` + blockers. Never `done`/`fixed`.
4. Escalate instead of proceeding when a task turns out to touch scientific
   assumptions, `PlanetProfile/Test/`, gate thresholds, or anything outside
   its written scope. Write the blocker in the report and stop.
5. Environment: `mamba activate PP` on this machine (PPcl is Machine B's).

## Tasks

### C7 — Warn when the hydrosphere exceeds the MgSO4 EOS table range [status: in progress (Codex)]

From the classic-MoI recon adjudication (`plans/active/titan-classic-moi-recon.md`
item 3): the MgSO4(aq) property table ends at 800 MPa, and thick Titan
hydrospheres reach ~1.8 GPa — above 800 MPa properties are silently
nearest-clamped. Add a clear, once-per-run log warning when an MgSO4 ocean
EOS is evaluated beyond its table's pressure ceiling, naming the ceiling,
the requested maximum pressure, and the clamping behavior. Warning ONLY —
do not change any EOS value, bound, or scientific assumption (extending the
table is a separate, reviewer-gated roadmap item). Pattern reference: the
NH3 branch's Pmax reset warning in `Thermodynamics/HydroEOS.py` (~line 296)
and the `_warn_once` mechanics used for NaCl. Scope: the MgSO4 path in
`HydroEOS.py`/its property loader; find the actual clamp site first and
name it in your report. Verification: a targeted unit-style test or REPL
run showing (a) the warning fires exactly once for a P grid crossing the
ceiling, (b) no warning for an in-range grid, (c) returned property values
are bit-identical before/after the change. `verified` requires all three
cited.

## Completed

### C5 — Deploy script: derive serve-time cache list from the slot registry [status: verified]

`plans/scripts/build_deploy_branch.sh` hand-maintains the list of structure
cache pkls to ship; it has now missed a newly wired slot twice (v5 cache in
July, NH3 cache `d53385f1`). Replace the hardcoded `for f in ...` list with
a step that extracts every `cache_path` from
`PlanetProfileApp/pages/Inference.py::_SBI_ARTIFACT_SLOTS` at build time
(small python helper invoked by the script is fine — e.g. import the module
with streamlit stubbed, or a conservative regex over the registry block;
prefer the import if clean). Keep the existing sidecar handling
(`.ae_sidecar.pkl`, `_offsets.json`) for every derived path. Placeholder
slots with `cache_path: None` must be skipped. Verification: run the script
build-only (comment/skip the git push or run in a scratch clone — do NOT
force-push app-deploy during testing), then diff the staged pkl set against
the registry cache_paths: exact match = `verified`; cite the diff output.

**Report (Codex, 2026-08-04):** `verified`. Added a literal-safe AST helper
that extracts unique non-None `cache_path` values without importing or
executing the Streamlit page. The deploy script now copies that derived set,
retains `.ae_sidecar.pkl` and `_offsets.json` handling per path, fails on a
missing registry cache, and blocks when staged primary PKLs differ from the
registry. Added `--build-only`, which commits the temporary snapshot but
skips both local and remote `app-deploy` pushes. Evidence from
`plans/scripts/build_deploy_branch.sh --build-only`: “Verified staged SBI
cache set: 6 registry paths, exact match”; no `diff -u` output; 299 MB, 528
files; “Build-only verification passed; skipped all app-deploy pushes”.
`bash -n`, helper `py_compile`, and `git diff --check` passed. Claim commit
`d8503fdb`; implementation commit `a28e11f0`.

### C6 — Wedge cold-edge legend polish [status: verified]

From C3 manager verification (cosmetic, non-blocking):
`PlanetProfileApp/Utilities/radial_profiles.py::build_perplex_geotherm_figure`
— the cold-edge `axvspan` legend swatch is white-on-white (invisible in the
legend box) and the 25 K band is faint against wide native T domains. Give
the span a visible edge (e.g. hatch pattern or light-blue fill with a
colored edgecolor) so both the on-plot band and its legend entry read
clearly, without obscuring the heatmap underneath. Bump `_GLOBE_FIG_VER`
(the figure changes shape). Verification: render CV3 + CI figures from the
Europa v4 reference posterior (the C3 verification recipe), save PNGs, and
confirm band + legend visibility; AppTest mineralogy test must stay green.
Top solo status: `implemented, unverified` + PNG paths for manager visual
check.

**Report (Codex, 2026-08-04):** `implemented, unverified`. Replaced the
white-on-white span with a pale-blue, blue-edged hatched patch at the same
exact `Tmin` to `Tmin + 25 K` interval and raised the red silicate path above
it. Bumped `_GLOBE_FIG_VER` 7→8. Structural assertions now require a hatch,
nonzero outline, and non-white fill. AppTest evidence:
`tests/app_globe_panel_test.py::test_apptest_mineralogy_tab` passed (1 passed,
one pre-existing TidalPy deprecation warning). Rendered Europa-v4
posterior-median artifacts:
`/private/tmp/ppgenai-c6-legend/europa_v4_CV3_cold_edge.png` and
`/private/tmp/ppgenai-c6-legend/europa_v4_CI_cold_edge.png`; both contain the
visible on-plot band and outlined/hatched legend swatch. `py_compile` and
`git diff --check` passed. Manager visual confirmation remains required.
Claim commit `d8503fdb`; implementation commit `5284a318`.

**Manager verification (Claude, 2026-08-04):** `verified`. Inspected the
rendered CV3 PNG: hatched pale-blue band with blue edge clearly visible at
the table cold edge, legend swatch readable, silicate path above the band;
suite 8/8. C5 also accepted: extractor code reviewed (literal-eval only,
path-traversal guards, blocking staged-vs-registry invariant) + Codex's
cited build-only evidence. Both commits pushed (rebased over B's PPC batch;
one STATUS.md merge resolved by the manager).

### C4 — Amortized widget-key namespacing + Titan slot stubs [status: verified]

Source: Machine B addendum `plans/STATUS-2026-08-01-machineB-joint-nh3.md`,
section "What Machine A can safely plan + change now". Target:
`PlanetProfileApp/pages/Inference.py` ONLY. Two parts:

1. **Widget-key namespacing (the must-fix):** the amortized observable-input
   widgets use global keys (e.g. `amort_obs_Re_k2`, ~line 1740) shared
   across slots. With multiple Titan slots + Clipper all exposing Re_k2,
   namespace the keys per slot id. Follow the precedent of commit
   `ec711672` (per-slot namespaced widget keys for Clipper).
2. **Titan slot stubs:** add four Titan entries to `_SBI_ARTIFACT_SLOTS`
   (~line 1049): no-ocean (Phase A, already shipped — real paths) + NH3
   (joint) + MgSO4 + NaCl. For the three not-yet-built slots use
   clearly-marked TODO placeholders for `cache_path` / `config_path` /
   `default_obs` centrals and gate status — do NOT hardcode NH3 artifact
   paths; Machine B hands those over. NH3 `scope_note` must state: JOINT
   no-ocean+ocean posterior (frozen Titan interiors included), gravity
   provenance Petricca 2025, CMR2 dropped (C22 double-count), induction+h2
   dropped (no clean Cassini signal), NH3 ocean w in [1,70] ppt. Gate the
   placeholder slots so they render as "awaiting artifact" rather than
   erroring.

FORBIDDEN (Machine B owns, collision risk): `cache_builder.py`,
`LayerPropagators.py`, `defineStructs.py`, `configs/test54_*`,
`PlanetProfile/Test/mcmc_results/Titan/Test54_nh3_ocean/`,
`plans/scripts/titanG_*`.

Verification: AppTest — existing slots still load; per-slot keys present and
distinct; placeholder Titan slots render their awaiting-artifact state
without exceptions. In-browser confirmation stays with the manager, so top
solo status is `implemented, unverified` + AppTest evidence.

**Report (Codex, 2026-08-02):** `implemented, unverified`. The earlier
`ec711672` repair already namespaced observable/sigma/truncation/Ae keys by
artifact; C4 preserves that behavior and gives the shipped Phase-A slot an
explicit stable `slot_id` namespace. Registered the NH₃ joint, MgSO₄, and
NaCl Titan placeholders with `artifact_filename`, `config_path`,
`cache_path`, and `default_obs` unset and gate fields marked TODO. Placeholder
discovery is deliberate and selectable, but returns before any file stat,
config read, or runner load; the Generate button is disabled. The shipped
Phase-A slot remains the default. The NH₃ scope note records the required
joint frozen/ocean support, Petricca-2025 gravity, CMR² and induction/h₂
exclusions, and `[1,70] ppt` range. AppTest structural evidence using
`AppTest.from_file(PlanetProfileApp/pages/Inference.py)`: all six shipped
top-level slots rendered without exceptions and with enabled Generate;
Phase-A, Test50, and Clipper `Re_k2` keys were distinct; Clipper retained
`Re_k2 = 0.23`; all three placeholders rendered `Awaiting artifact`, TODO
gate status, and disabled Generate without exceptions. `py_compile` and
`git diff --check` passed. Browser confirmation was not performed and remains
with the manager. Claim commit `220ec593`; implementation commit `f4090c33`.

**Manager verification (Claude, 2026-08-02):** `verified`. Independent
scripted AppTest: Cassini–Titan group lists all five versions with the
shipped Phase-A slot as default; each of the three placeholders renders its
own label, the awaiting-artifact info, TODO gate status, and zero observable
widgets; the NH3 scope note carries the joint/Petricca/[1,70] ppt content;
Titan vs Clipper Re_k2 keys are disjoint with correct per-slot values
(0.608 / 0.23). One earlier apparent NH3-scope failure was a manager test
harness bug (stale selectbox handle), not a code defect. Full
`tests/app_globe_panel_test.py` suite 7/7 after the follow-on gating fix.

**Side finding during verification (manager, fixed in `cb597490`):** the
v5/v6 slot registrations (`7211f002`) predate the adjudication hold and
made unratified v6 the silent Clipper default. Fixed by gating those six
slots with `artifact_status='not_ratified'` via the C4 status mechanism;
v4 is the Clipper default again.

### C3 — Perple_X geotherm-overlay figure in Mineralogy tab [status: verified]

`PlanetProfileApp/pages/Inference.py` Mineralogy tab + helpers in
`PlanetProfileApp/Utilities/radial_profiles.py`. Add an expander "Geotherm vs
Perple_X table" showing, for a user-selected composition from
`PERPLEX_SILICATE_TABLES`: the table's grain rho(P,T) as a heatmap over its
native P,T domain, with the selected draw's silicate (P,T) path overlaid as a
line, and the table's cold-edge band (T < Tmin + 25 K) shaded. Reuse the
existing table-loading path (`GetInnerEOS` RectBivariateSpline machinery via
`mineralogy_for_sample`'s loaders — do not re-parse .tab files by hand).
Read the `developing-with-streamlit` conventions already used in the file
(figure token caching, widget keys namespaced per tab). Verification:
AppTest asserting the expander, selectbox, and figure object exist and the
overlay path P,T ranges match the sampled profile = structural verification;
in-browser/visual confirmation stays with the manager, so top status you can
claim alone is `implemented, unverified` + AppTest evidence.

**Report (Codex, 2026-08-02):** `implemented, unverified`. Added the exact
expander and a Mineralogy-namespaced composition selectbox backed by
`PERPLEX_SILICATE_TABLES`. The heatmap uses `GetInnerEOS` and the loaded
RectBivariateSpline knots over the native P,T endpoints, with bounded cached
sampling; it overlays the selected draw's existing silicate P,T path and
shades `T < Tmin + 25 K`. The crisp figure cache token includes the result,
sample, composition, and figure version. Structural evidence:
`tests/app_globe_panel_test.py::test_apptest_mineralogy_tab` passed (1 passed,
one pre-existing TidalPy deprecation warning), asserting the expander,
selectbox, figure exports/object, grid shape/native endpoints, and exact path
range equality with `mineralogy_for_sample`. `py_compile` and
`git diff --check` passed. Browser/visual confirmation was not performed and
remains with the manager. Claim commit `e41d725a`; implementation commit
`2a31c345`.

**Manager verification (Claude, 2026-08-02):** `verified`. Re-ran
`tests/app_globe_panel_test.py::test_apptest_mineralogy_tab` (1 passed,
19 s). Rendered the figure from the Europa v4 reference posterior for CV3
and CI tables and inspected the PNGs: heatmap spans the native domain,
posterior-median silicate path (131–3472 MPa, 265–1592 K) overlays in the
correct position with inverted pressure axis, cold-edge band present at the
table's cold edge. Cosmetic only (non-blocking, future polish): the white
cold-edge legend swatch is invisible on the white legend background, and the
25 K band is faint against wide native T domains.

### C2 — Titan classic-MoI reconnaissance [status: verified]

Report-only reconnaissance of the silicate EOS, core, porosity, and
hydrosphere parameter space. Deliverable:
`plans/active/titan-classic-moi-recon.md`.

**Report (Codex, 2026-08-02):** `verified`. Targeted full-profile runs
exercised `FindInnerWithMoIAndEOS`; a direct 1-km core-radius integration used
the production `structure_derivation` helpers. Total runtime was under one
minute. The report renders with Python-Markdown (7 tables, 2 code blocks) and
`git diff --check` passed. Finding for manager adjudication: current
self-consistent families top out near 0.3395--0.3398; numerical 0.341 crossings
require either the reduced uniform-density stack or an MgSO4/CV3 boundary case
outside current EOS/density safeguards. Claim commit `e8c73142`;
implementation commit `400f84e5`.

### C1 — Artifact ledger refresh [status: verified]

Updated `PlanetProfile/Inference/sbi_artifacts/INDEX.md` so every shipped
artifact has a current row: v5 trio, v6 trio, v7 open-|Ae|, titan_freegrav —
state `delivered, not ratified (gate adjudication open; see plans/STATUS.md)`.
Did not change deployed/retired/vetoed states; corrected stale GUI references
and malformed table cells.

**Report (Codex, 2026-08-02):** `verified`. Python-Markdown rendered both
tables. A structural audit cross-checked `Path.glob('*.pt')` against artifact
row names: 16 files, 16 rows, zero missing, extra, duplicate, or malformed
rows. `git diff --check` passed. Claim commit `7bac3207`; implementation commit
`49d39eba`.
