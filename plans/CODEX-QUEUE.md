# Codex 5.6 queue (Machine A delegate lane)

Updated: 2026-08-12 (C16t verified; C17 next, then C18→C19 in order, one implementation commit per task). Curated by the Claude model manager
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

_Second doc audit (2026-08-11, full report in the audit agent record;
findings cited per task) + program-review items. All are doc/test/code
tasks with NO scientific adjudication — copy authoritative wording, never
re-adjudicate; escalate per protocol 4 on anything ambiguous._

### C18 — Enceladus libration regression test (freeze blocker B3) [status: queued]

`PlanetProfile/Gravity/Librations.py` has ZERO test coverage and is about
to carry the Enceladus campaign's dominant observable (0.092 ± 0.003 deg,
3.3% relative). Write `tests/librations_test.py`: (1) reproduce a
PUBLISHED Enceladus libration for a stated interior — target Van Hoolst
et al. 2016 or Thomas et al. 2016 parameters (find the cleanest published
(interior, amplitude) pair; document the source and tolerances); (2) the
rigid-vs-elastic branch consistency check (rigid=True vs False within
0.1%; recorded values 0.11091/0.11094 deg); (3) monotonicity: libration
decreases with increasing shell thickness over 5-45 km for a fixed
Enceladus-like interior. Escalate per protocol 4 if reproducing the
published value requires ANY physics change — report the discrepancy,
do not fix. Report `verified` with pytest output.

### C19 — Enceladus induction plumbing (freeze blockers B8-B10) [status: queued]

1. B8: add the synodic/PPO excitation — `Texc_hr['Enceladus']['synodic']
   = 15.559` in `PlanetProfile/MagneticInduction/Moments.py:55` (cite
   Saur et al. 2024 §2.2, 1-2 nT driver) + the corresponding Be table
   row (derive from the existing Be1xyz_Enceladus.txt pattern; if the
   synodic Be amplitude cannot be sourced from the committed table or
   the paper, ESCALATE rather than invent).
2. B9: fix Enceladus Be-file resolution — `GetBexc`
   (MagneticInduction.py:676-687) looks for
   `Be1xyz_Enceladus_Cassini_Cassini11noMP.txt`, which does not exist;
   `inductionData/Be1xyz_Enceladus.txt` does (B0 z = -336.5 nT, sane vs
   Saur ~330 nT). Fix by file rename/copy or loader-path correction —
   whichever matches the conventions of the other bodies; NO value
   changes.
3. B10: formalize the Saur Fig. 20 conductivity acceptance test —
   `tests/enceladus_conductivity_test.py`: Seawater sigma at P=3 MPa,
   w in {5,10,20} ppt, T in {271.15, 272.15, 274.15} K must land in
   0.45-1.79 S/m (band verified by the reviewer 2026-08-12) and stay
   under 2 S/m at the plume band; plus sigma(w=70 ppt) in 4.5-6 S/m
   bracketing Saur's ~5 S/m ceiling.
Verification: pytest green + one InducedAeList smoke at both Enceladus
periods returning finite Ae. Report `verified`.

### C13 — INDEX.md + DEPLOYING.md + plans-index currency (audit items 1.1-1.3, 3.1-3.2, 2.4, archive) [status: verified]

1. `PlanetProfile/Inference/sbi_artifacts/INDEX.md`: bump the audit
   stamp; REWRITE the "Deployment gate" section to describe the actual
   adjudicated/split-ratification regime (gate FAIL + recorded
   adjudication + scope note may deploy; cite the FAIL-ADJUDICATED
   vocabulary rule) instead of "all three PASS"; move the three wired
   Titan joint rows (NH3/MgSO4/NaCl) into the deployed-artifacts table
   (or retitle tables wired/not-wired); add an SBI config-hash column
   entry for the two salts alongside the existing "(mcmc cfg hash)"
   values, labeled.
2. `DEPLOYING.md`: "six" -> "eight" caches, add the two Titan salt cache
   paths, size ~300 -> ~320 MB (also fix the same size figure in
   `plans/scripts/build_deploy_branch.sh` header comment).
3. `plans/active/README.md`: refresh (stamp; C-queue state; v5/v6
   finalized; NH3 Phase B complete; ADD tidal-sector-remedy-plan.md as
   the live plan). Add to `plans/archive/README.md` index:
   europa-clipper-v4-geodesy-plan.md, HANDOFF-2026-07-26-nh3-liquidus-
   defect.md, STATUS-2026-07-20.md, STATUS-2026-08-01-machineB-joint-
   nh3.md, plans/active/titan-classic-moi-recon.md, plans/active/
   radiogenic-inventory-and-mineralogy-plan.md (index-based archiving,
   no file moves).
Verification: link check + `git diff --check`; report `verified`.

**Report (Codex, 2026-08-12):** `verified`. Rewrote the artifact deployment
gate to preserve raw FAIL results plus recorded scientific adjudication and
the `FAIL-ADJUDICATED-ACCEPTABLE` vocabulary; promoted the wired NH3, MgSO4,
and NaCl joint Titan rows into the deployed table; and labeled both SBI and
MCMC config hashes for the salts from their committed manifests. Updated the
deploy snapshot documentation to the exact eight-cache registry set and
~320 MB size, including the build-script header. Refreshed the active plan
index for the finalized v5/v6 rulings, completed Titan composition campaign,
live tidal-sector remedy and Enceladus plans, and added all six requested
provenance entries to the archive index without moving files. Verification
output: “Changed Markdown local-reference check: 22/22 references resolve”;
“Deploy cache documentation check: 8/8 registry paths exact”; “Salt
config-hash provenance check: 4/4 hashes match manifests”; “Titan joint row
promotion check: 3/3 appear exactly once”; `git diff --check: PASS`. Claim
commit `a1246a18`; implementation commit `411167fb`.

**Manager review (2026-08-12): C13 ACCEPTED.** Independently verified:
INDEX stamp/gate-policy rewrite/salt-row promotion/labeled dual hashes;
DEPLOYING.md 8-path list matches the registry extractor exactly (both
salt caches present, ~320 MB in doc AND script header); plans/active/
README refreshed incl. the two new Enceladus docs (good initiative);
6/6 archive-index additions present. Queue continues at C14.

### C14 — Salt RATIFICATION.md consolidation + gate-manifest schema (audit 5.1-5.2) [status: verified]

1. Write `validation_reports/titan_freegrav_{mgso4,nacl}_1m/
   RATIFICATION.md` mirroring the NH3 template: verdict, sector split,
   per-gate outcomes, adjudications (NaCl log10_eta_V
   FAIL-ADJUDICATED-ACCEPTABLE; limits monotone gated-FAIL override),
   reviewer agents + dates, countersign date 2026-08-10. STRICTLY
   copy-only from committed sources (INDEX rows, STATUS 2026-08-10
   entries, MACHINE-B-HANDOFF §0.15/0.16, gate reports). The manager
   countersigns the files after your commit — flag any statement you
   cannot source to a committed document instead of writing it.
2. Extend the gate driver's manifest writing (plans/scripts/
   titanG_run_gates.py or equivalent) to record: all three gate exit
   codes, repo-relative paths, git SHA, torch/sbi versions, artifact
   sha256. Future runs only; do not rewrite committed manifests.
Verification: RATIFICATION files render; every claim carries a source
pointer; manifest schema change unit-smoked. Report `verified`.

**Report (Codex, 2026-08-12):** `verified`. Added copy-only MgSO4 and NaCl
`RATIFICATION.md` files with the 2026-08-10 split verdict, sector scope,
raw per-gate outcomes, limits gated-FAIL override, NaCl
`log10_eta_V = FAIL-ADJUDICATED-ACCEPTABLE`, reviewer agent/date records,
and the underlying campaign countersign date. Every factual paragraph and
gate-table row cites the committed INDEX, STATUS, MACHINE-B-HANDOFF §0.16,
or a local JSON report; each file explicitly awaits the manager's post-commit
file-level countersign. Extended only the future-run salt gate-manifest path:
schema v2 records all three gate outcomes (merging partial reruns),
repo-relative artifact/config/reference/report/log paths, full Git SHA,
torch/sbi versions, UTC generation time, and artifact SHA-256. Existing
committed manifests were not rewritten. Verification: focused pytest
`tests/titanG_gate_manifest_test.py` — 1 passed; both Markdown documents
rendered with gate tables and 9/9 source links resolved; claim-to-report
assertions PASS; committed-manifest unchanged check PASS; script `--help`,
`py_compile`, and `git diff --check` PASS. Claim commit `cc53dc5d`;
implementation commit `f282b669`.

**Manager review (2026-08-12): C14 ACCEPTED + RATIFICATION files
COUNTERSIGNED.** File-level verification done: every gate-table number
re-checked against the per-gate JSONs (MgSO4 SBC 1024 kept / min BH-adj
p 0.1559 recomputed; NaCl 0.352-vs-0.30 eta_V, 0.00722-sigma containment)
and the 2026-08-10 adjudication records; source links 9/9; manifest
schema-v2 test passes; committed manifests untouched. The two
RATIFICATION.md files are now the canonical citation targets for the salt
campaigns. Queue continues at C15.

### C15 — Methodology + GUI docs refresh (audit 2.1-2.2) [status: verified]

1. `docs_ai/AMORTIZED_SBI_METHODOLOGY.md`: add the PPC/pushforward gate
   as standard equipment (four-way table, 0.5 sigma_obs flag); add a
   Track 1 IS-correction section (weights, binding reliability set:
   Pareto-k/absolute-ESS/w_max/reverse-coverage; pointer to
   PlanetProfile/Inference/is_correction.py and plans/active/
   tidal-sector-remedy-plan.md); refresh the status stamp. Copy-only
   from the remedy plan + module docstring.
2. `PlanetProfileApp/README.md`: document the newer GUI capabilities
   (interactive globe incl. static no-WebGL render + per-phase shells,
   k2 measurement-ellipse zoom, figure-provenance exports, sector
   warnings / split-ratification display, slot scope notes).
Verification: link check; strings match module constants. Report
`verified`.

**Report (Codex, 2026-08-12):** `verified`. Refreshed the amortized-SBI
methodology status and made posterior-predictive/pushforward validation
standard equipment with the observed/prior-PP/reference-MCMC-PP/SBI-PP
four-way table, the 0.5-sigma_obs median flag, and the corrected pipeline's
weighted KS/W1 requirement. Added the Track 1 exact-likelihood importance
correction, weight equation, source pointers, binding Pareto-k/absolute-ESS/
w_max/reverse-coverage set, supporting branch/support guards, and the rule
that successful correction validates the MCMC target rather than erasing
shared model-data tension. Updated the GUI README with the interactive globe,
exact “Static render (no WebGL)” fallback, proportional Ice III/V/VI shells,
exact k2 ellipse-zoom control, provenance-bearing figure exports, visible
scope/gate text, and result-level split-sector warnings; removed its stale
nonexistent `FIGURE_CHANGES.md` link. Verification output: both Markdown
documents rendered; methodology local links 4/4 and README local links 0/0
resolve; 10/10 IS constants match `is_correction.py`; GUI text matches 3/3
control labels, 2/2 slot keys, 4/4 export-provenance keys, and 3/3 phase names;
`git diff --check: PASS`. Claim commit `12c9ed53`; implementation commit
`87f4ca31`.

**Manager review (2026-08-12): C15 ACCEPTED.** Spot-verified: methodology
doc now carries the PPC four-way table + 0.5-sigma flag and the Track 1
IS-correction section with the binding reliability set matching
is_correction.py constants; GUI README documents the globe (static
no-WebGL fallback, per-phase shells), k2 ellipse zoom, provenance
exports, and sector warnings with exact control labels; stale
FIGURE_CHANGES.md link removed. Queue continues at C16t.

### C16t — Slot-level end-to-end reproduction regression test (program review item) [status: verified]

New test `tests/slot_reproduction_test.py`: for every registry slot with
an existing artifact+config+cache, assert (1) the slot config's hash
matches an artifact (or the slot's documented hash), (2) cache file SHA
matches its committed manifest where a manifest exists, (3) for ONE
designated fast slot (Titan Test50), a fixed-seed 1k-draw conditioning
reproduces committed summary statistics within tolerance (commit the
reference stats file on first run). Marked slow-safe: full-suite runtime
under ~5 min; heavier per-slot draws are out of scope. This turns the
traceability principle into a testable claim. Report `verified` with
pytest output.

**Report (Codex, 2026-08-12):** `verified`. Added
`tests/slot_reproduction_test.py` and the first-run Test50 reference record
`tests/data/slot_reproduction_test50.json`. The regression discovers the
literal GUI registry without executing the Streamlit page; all 13 locally
complete slots match their artifact config hash after the documented
MCMC-to-SBI mode switch, and each config cache path matches its slot. Seven
committed manifest records covering five ready-slot caches reproduce their
installed cache SHA-256. The designated Test50 slot reproduces mean, sample
standard deviation, median, and 5th/95th percentiles for all eight parameters
from 1,000 posterior draws at seed 1601 (`rtol=1e-4`, `atol=1e-5`), with the
artifact and cache SHA-256 pinned in the reference record. Verification:
`python -m pytest -q tests/slot_reproduction_test.py` = **15 passed, 2
warnings in 2.20s** (pre-existing TidalPy deprecation and nflows/torch API
warnings); `py_compile`, JSON parse, and `git diff --check` PASS. Claim commit
`cd1c007c`; implementation commit `798dc937`.

### C17 — Expose degree-2 C20/C22 as a standard PlanetProfile run output (core-parity, STRATEGY exception (b)) [status: in progress (Codex)]

The direct-Clairaut calculation exists in
`PlanetProfile/Inference/gravity_obs.py`; only its invocation path is
inference-only. Wire it into the standard run output path (Planet
attribute + printout parity with other derived quantities) behind a Do
flag defaulting OFF (no behavior change unless enabled), so a single
PlanetProfile run can reproduce the inference gravity observables.
DO NOT alter the calculation itself (protocol 4: escalate if any
physics change seems needed). Verification: targeted test comparing the
new standard-run output against gravity_obs directly on one Europa +
one Titan case; suite green. Report `verified`.

_Prior triage (C9-C12, closed): triage of the C8 doc-doctor report (2 PASS / 10 FINDING,
`validation_reports/doc_doctor/2026-08-06_first_pass.md`). Already resolved
elsewhere: v5/v6 scope-note superseded language + INDEX rows (manager
ratification session 2026-08-07); v4 direct-Clairaut exception and the
Geotherm-tab deviation are now RECORDED in `plans/STRATEGY.md` (items 5, 8 —
no code change wanted). Work C9→C12 in order; one commit per task._

### C9 — Ledger + routing corrections (doc-doctor items 2, 4, 9, 10) [status: verified]

Report-driven doc fixes; no code, no science:
1. `PlanetProfile/Inference/sbi_artifacts/INDEX.md`: (a) Test50 row — add a
   direct gate pointer into `validation_reports/` (find the Test50 gate
   report path); (b) Titan freegrav no-ocean row — replace "not wired for
   deployment" with the actual wired/deployed state (slot at
   `Inference.py:1097`); (c) NH3 row — replace the unqualified
   "RATIFIED (verified)" and the benign-eta-medians crosscheck reduction
   with the split-ratification wording from
   `validation_reports/titan_freegrav_nh3_1m/RATIFICATION.md` (structure
   sector ratified; tidal k2/zeta/eta sector QUARANTINED, MCMC
   authoritative; crosscheck FAIL not to be read as benign-nuisances-only).
   Copy the authoritative wording — do NOT re-adjudicate.
2. `git add` the untracked routing trees so entry points resolve in a clean
   checkout: `plans/README.md`, `plans/archive/`,
   `plans/scripts/{active,archive,reproducibility}/` (doc-doctor item 9).
   Then extend `plans/README.md` routing with `plans/CODEX-QUEUE.md`,
   `plans/STRATEGY.md`, `plans/DOC-DOCTOR.md`.
3. `DEPLOYING.md`: update ~205 MB → ~300 MB; replace the "two manually
   listed caches / copy list" text with the registry-derived cache set
   (six paths, exact-equality enforcement,
   `plans/scripts/build_deploy_branch.sh:55-96`); document `--no-push` and
   `--build-only`.
Verification: markdown-link check on the touched files + `git status`
clean-checkout resolution of every routed path. Report `verified` with the
check output.

**Report (Codex, 2026-08-07):** `verified`. Corrected the Test50 ledger with
the committed Machine A cross-version gate path; recorded the Titan
free-gravity no-ocean slot's actual 2026-07-24 deployment; and replaced the
NH3 row's superseded benign-nuisance framing with the authoritative split
ratification (gravity/structure verified, tidal/dissipation not verified,
MCMC authoritative). Added the five previously untracked routing index files
and linked `CODEX-QUEUE.md`, `STRATEGY.md`, and `DOC-DOCTOR.md` from
`plans/README.md`. Updated `DEPLOYING.md` to ~300 MB, documented
`--build-only`/`--no-push`, and listed the six registry-derived caches plus
the exact-equality invariant. Verification output: “Markdown local-reference
check: 6/6 changed references resolve”; “Deploy cache documentation check:
6/6 paths exactly match registry”; “Clean-checkout routing check: 9/9 routed
entry points tracked”; `git diff --cached --check: PASS`. Claim commit
`9cee9c90`; implementation commit `8f969076`.

### C10 — GUI assumption-text contradictions (doc-doctor items 6, 7) [status: verified]

Three factual corrections + one caption, text-only (no behavior change):
1. `Inference.py:2517-2523` shared build-up expander: "nearest Tb grid
   structure" → between-node bilinear blending (match
   `run_assumptions.py:77-86`).
2. `run_assumptions.py:115-119`: gravity-pair configs claim Radau–Darwin
   C20/C22; v4+ uses direct Clairaut (`gravity_obs.py`, v4 config `:198`).
   Make the text config-conditional (say what the loaded config declares),
   not a blanket claim.
3. `run_assumptions.py:48-51`: "emcee ensemble sampler" → pocoMC
   (preconditioned Monte Carlo), matching `Inference.py:2566-2570`.
4. Add the bilinear-blending deviation to the captions of the affected
   results panels (radial profiles `Inference.py:3906-3910`, plus
   wedge/heating/mineralogy captions where the rho_sil disclosure already
   sits) — one short clause each, mirroring the existing disclosure style.
Verification: AppTest render of the expander + one results page asserting
the new strings; app suite still green. Report `verified`.

**Report (Codex, 2026-08-07):** `verified`. The shared build-up expander now
describes transition-aware between-node structure blending (linear in Tb for
1D list caches and bilinear in Tb/log-salinity for 2D caches). Per-run text
now names pocoMC and reports the loaded config's gravity forward model;
`clairaut_hydrostatic` is described as direct Clairaut integration rather
than Radau–Darwin. Added the required short bilinear-blending clause to the
radial-profile, wedge, heating, and mineralogy captions. AppTest asserts the
new expander text, all four result-caption disclosures, pocoMC, and
config-conditional Clairaut wording: `tests/app_globe_panel_test.py` — 8
passed (10 pre-existing warnings) with a task-specific writable Numba cache.
`py_compile` and `git diff --check` passed. Claim commit `4776a320`;
implementation commit `7f5f3836`.

### C11 — Scope-note numerical priors + config-exact centrals (doc-doctor item 3) [status: verified]

1. Add numerical headline-prior ranges to the scope notes missing them:
   Test50 (`Inference.py:1088`), Galileo v1.1 (`:1274`), v4 (`:1320`),
   v5/v6 baseline + ablation notes (at least salinity and the
   non-hydrostatic offset boxes). Copy ranges from each slot's config —
   configs are the source of truth; do not invent or round beyond config
   precision. Test50: also enumerate dropped channels with reasons (from
   its config/campaign doc).
2. Replace the three v5 registry C20/C22 centrals with config-exact values
   from `europa_clipper_v5_geodesy_11D.json:101-108`
   (-4.57888786233524e-4 / 1.3775234242885802e-4) at
   `Inference.py:1356-1374,1410-1415,1435-1451`.
CAUTION: do not alter any other number; if a config value seems wrong,
escalate (protocol 4) — no silent corrections.
Verification: AppTest asserting a prior range renders in each touched scope
note; diff shows only scope-note strings + the three centrals. Report
`verified`.

**Report (Codex, 2026-08-07):** `verified`. Added config-exact headline
prior ranges to all nine requested Test50, Galileo v1.1, v4, v5, and v6
scope notes. Test50 now states its two-channel conditioning and the campaign's
documented CMR2 exclusion (core-blind 2.5e-5 span = 0.025 observational
sigma; production MCMC used the same two channels). Replaced C20/C22 in all
three v5 registry entries with the exact config literals
`-4.57888786233524e-4` / `1.3775234242885802e-4`. A direct AppTest rendered
9/9 touched scope notes with their expected ranges; an AST/config comparison
reported “3/3 slots exact for C20 and C22”; `py_compile` and `git diff
--check` passed. The selectively staged C11 patch contained only scope-note
strings and the three central pairs; the shared manager process incorporated
it into its concurrent commit `607c6d34` alongside unrelated WebGL work.
Claim commit after the manager's concurrent history rewrite: `6b173c24`.

### C12 — Exported-figure provenance (doc-doctor item 12) [status: verified]

`PlanetProfileApp/Utilities/crisp_figs.py`: embed provenance in exports —
(a) `savefig` metadata (PDF/SVG/PNG support differs; set title/subject or
XMP-ish keys where the backend allows) carrying model version (slot label +
artifact filename), conditioning summary, and app git SHA if cheaply
available; (b) filenames `<figure>_<slot-short>_<date>.<ext>` instead of
bare `k2.pdf`/`corner.png`. Keep figure TITLES unchanged (visual layout is
frozen; provenance goes in metadata/filename only). Thread the slot/
conditioning info from the results context with minimal plumbing.
Verification: download one PDF + one PNG via AppTest or direct call,
inspect embedded metadata (pikepdf/PIL), cite the inspection output.
Report `verified`.

**Report (Codex, 2026-08-07):** `verified`. Every crisp figure export now
receives the active result's model-slot label, artifact filename, and exact
conditioning values; old result pickles recover their slot from the artifact
registry, while custom MCMC runs receive an explicit fallback label. PDF,
SVG, and PNG metadata add that context plus the cheaply resolved app git SHA,
without changing the matplotlib figure title. Download names now follow
`<figure>_<slot-short>_<UTC-date>.<ext>`. Direct export inspection produced
`k2_clipper-v5-geodesy_2026-08-07.pdf`; pikepdf reported `/Creator =
PlanetProfileApp` and `/Subject = Model slot: 1D · Clipper–Europa — v5
geodesy 11D; artifact: europa_clipper_v5_geodesy_11D_posterior_1m.pt;
conditioning: C20=-0.0004355, C22=0.0001309, Re_k2=0.23; app git SHA:
c12-test-sha`. PIL reported the same provenance in PNG `Description`, with
`Software = PlanetProfileApp`; SVG inspection found the same subject and the
unchanged visual title. Automatic SHA resolution matched current HEAD
(`9d76375e` at inspection time). Representative results-page AppTest
`tests/app_globe_panel_test.py::test_apptest_europa_v4_globe_panel` passed
(1 passed, one pre-existing TidalPy deprecation warning). `py_compile`,
`git diff --check`, and staged diff checks passed. Claim commit `93a93c4e`;
implementation commit `da59b261`.

**Manager review (2026-08-07): C9–C12 ACCEPTED and pushed.** C12
independently re-verified on this machine: SVG dc:description, PDF
/Subject, and PNG Description all carry slot label + artifact filename +
conditioning + live git SHA (auto-resolved 7b4c9f9c); filename stem
`k2_<slot>_<date>` confirmed; app suite 9/9 with C12 in tree. Queue empty.
Process note: C11's staged Inference.py edits were swept into the
manager's concurrent globe commit `607c6d34` (shared worktree) — content
verified, but future concurrent lanes must claim files or stash before
manager commits touch the same file.

## Completed

### C8 — First doc-doctor pass [status: verified]

Execute the full checklist in `plans/DOC-DOCTOR.md` against the current
tree (all wired GUI slots incl. gated ones, README/DEPLOYING, the NH3 and
v4 campaigns for methods-extractability). Deliverable: a report at
`validation_reports/doc_doctor/2026-08-06_first_pass.md` with per-item
PASS / FINDING (file:line) / N/A and a short findings summary at top.
Report-only — do NOT fix findings in this task; the manager triages them
into follow-up tasks. Verification: `verified` = every checklist item has
an entry with evidence cited (grep/paths), and the report renders.

**Report (Codex, 2026-08-06):** `verified`. Audited all 13 GUI registry
entries (11 concrete artifacts plus two intentional placeholders), root and
plans README routing, DEPLOYING versus the deploy script, all checklist
reproducibility/organization/caption requirements, and methods provenance for
Titan NH3 joint and Europa v4. The report at
`validation_reports/doc_doctor/2026-08-06_first_pass.md` records every item
1-12 as PASS/FINDING/N/A with file:line evidence and a findings summary; no
findings were fixed. Result: 2 PASS, 10 FINDING, with two placeholder
subchecks N/A. Machine checks found 11/11 concrete config/cache paths exact,
all concrete files present, three v5 observable-central rounding differences,
all root README local targets present, and six deploy caches. Python-Markdown
rendered the report successfully (16,976 Markdown chars to 20,696 HTML chars;
12/12 checklist headings and one table); `git diff --check` passed. Claim
commit `16be13e1`; report commit `94ded63f`.

### C7 — Warn when the hydrosphere exceeds the MgSO4 EOS table range [status: verified]

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

**Report (Codex, 2026-08-05):** `verified`. The actual clamp site is
`PlanetProfile/Thermodynamics/MgSO4/MgSO4Props.py::MgSO4Props`, where
`ResetNearestExtrap` nearest-clamps pressure to
`MgSO4propsLookup.Pmax` before interpolation. Added a module-level
once-per-run diagnostic immediately before that existing path; it names the
800 MPa ceiling, requested maximum, and nearest-table-value behavior without
changing inputs, outputs, bounds, or scientific assumptions. Focused pytest
artifact `tests/mgso4_props_warning_test.py`: 2 passed, covering exactly one
warning across repeated crossing grids, zero warnings for an in-range grid,
and bit-identical results versus the explicit legacy clamp. A separate
pre-edit/post-edit REPL comparison for P `[0, 400, 900, 1800]` MPa and T
`[260, 290]` K asserted equal dtype, shape, `array_equal(equal_nan=True)`,
and SHA-256 bytes for all six returns; property hashes matched at rho
`ca4bbd6319488374265a9db34ed518c873f8908eb3e248e5ae3919ab53f40dde`, Cp
`cbb6e443ae0085ad363f89be0459abe8d0a12264ac66c2d27d23db8cb5def456`,
alpha `609ee380f382810e648ce531ad813fd440c0e9becae78ed58dda547aa595b0d8`,
and kTherm
`3c7c1c06d5b4bf4315f995036fca0e82eaab468da3246e94d1bd3e1daddd53b1`.
`py_compile` and `git diff --check` passed. Claim commit `379e0bf2`;
implementation commit `904dd704`.

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
