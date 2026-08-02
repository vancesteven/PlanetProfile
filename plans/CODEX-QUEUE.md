# Codex 5.6 queue (Machine A delegate lane)

Updated: 2026-08-02 at genai `9cfc3c14`. Curated by the Claude model manager
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

### C1 — Artifact ledger refresh [status: in progress (Codex)]

Update `PlanetProfile/Inference/sbi_artifacts/INDEX.md` so every shipped
artifact has a current row: v5 trio, v6 trio, v7 open-|Ae|, titan_freegrav —
state `delivered, not ratified (gate adjudication open; see plans/STATUS.md)`.
Do not change rows for deployed/retired/vetoed artifacts except to fix stale
cross-references. Doc-only; no code. Verification: `verified` = the file
renders correctly and every `.pt` in `sbi_artifacts/` has exactly one row
(cross-check with `ls`).

### C2 — Titan classic-MoI reconnaissance [status: queued]

Report-only, NO code changes. Classic MoI-matching Titan runs historically
fail with achievable C/MR^2 ~0.317 vs Cmeasured 0.341. Survey the Sil/core
parameter space to answer: which (silicate EOS table, core size/density,
porosity, hydrosphere) combinations can reach 0.341? Use short targeted runs
or direct MoI integration over candidate structures — NOT a heavy sweep
(Machine B owns heavy compute; keep total runtime under ~30 min). Deliver:
`plans/active/titan-classic-moi-recon.md` with the achievable C/MR^2 range
per configuration family, what limits it, and whether any physically
defensible configuration closes the gap. Cite exact configs/commands so runs
are reproducible.

### C3 — Perple_X geotherm-overlay figure in Mineralogy tab [status: queued]

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

(none yet)
