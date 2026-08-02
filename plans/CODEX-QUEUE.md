# Codex 5.6 queue (Machine A delegate lane)

Updated: 2026-08-02 at genai `e9df7e86`. Curated by the Claude model manager
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

### C4 — Amortized widget-key namespacing + Titan slot stubs [status: in progress (Codex)]

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

## Completed

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
