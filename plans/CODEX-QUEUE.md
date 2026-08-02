# Codex 5.6 queue (Machine A delegate lane)

Updated: 2026-08-02 at genai `400f84e5`. Curated by the Claude model manager
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

### C3 — Perple_X geotherm-overlay figure in Mineralogy tab [status: in progress (Codex)]

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

## Completed

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
