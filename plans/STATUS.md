# PlanetProfile genai — current status

Updated: 2026-08-02 at genai `f4090c33`. Single current-state source of truth.
Freshness rule: any session that pushes commits, integrates Machine B
artifacts, or changes a queue must refresh this file's `Updated:` line and the
affected sections in the same session. If this file is >7 days older than
`git log -1`, refresh it before relying on it.

## Roster and lanes (2026-08-01)

| Lane | Models | Role | Queue |
|---|---|---|---|
| Machine A — Claude Code | Fable 5 (manager) + Opus 5 / Sonnet 5 subagents | planning, scientific review/adjudication, integration, GUI, smokes, cross-machine coordination | this file |
| Machine A — Codex 5.6 | Codex 5.6 | delegated implementation + reconnaissance tasks | `plans/CODEX-QUEUE.md` (instructions: `AGENTS.md`) |
| Machine B — Claude Code | Opus 4.8 / Sonnet 4.6 / Haiku 4.5 | heavy compute: cache builds, reference MCMC, SBI training, gates | `plans/MACHINE-B-HANDOFF.md` |

Claude Fable 5 on Machine A is the model manager: it owns adjudication and
decides what enters the Codex and Machine B queues. Scientific verdicts are
never delegated to Codex.

## SBI slot states (GUI `_SBI_ARTIFACT_SLOTS`)

| Slot | State | Notes |
|---|---|---|
| Titan (Andrade, no-ocean) Test50 8D | DEPLOYED | unchanged |
| Titan freegrav no-ocean | delivered | artifact + gate reports in repo |
| Titan NH₃ joint no-ocean+ocean | **awaiting artifact** | GUI placeholder; Machine B Task #68 in flight |
| Titan MgSO₄ ocean | **awaiting artifact** | GUI placeholder; paths/centrals/gates unset |
| Titan NaCl ocean | **awaiting artifact** | GUI placeholder; paths/centrals/gates unset |
| Galileo–Europa v1.1 8D | DEPLOYED 2026-07-20 | gates 8/8 |
| Clipper–Europa v4 geodesy 11D | DEPLOYED 2026-07-21 | user-ratified |
| Clipper–Europa v5 thick-ice ablation trio | **delivered, NOT ratified** | baseline gates: sbc 0, limits 2, crosscheck 2; nok2 arm sbc 2 |
| Clipper–Europa v6 freegrav trio | **delivered, NOT ratified** | baseline gates: sbc 0, limits 2, crosscheck 2; noinduction arm sbc 2 |
| Clipper–Europa v7 open-\|Ae\| 11D | **delivered, NOT ratified** | gates: sbc 0, limits 2, crosscheck 2 |
| v1 / v2 | RETIRED | replaced by v1.1 / v4 |
| v3 | VETOED | wrong k2 bounds; 2D cache reused by v4+ |

## OPEN ADJUDICATION — v5/v6/v7 gate failures are systematic (2026-08-01)

All three delivered baselines share the identical failure signature
(`validation_reports/v{5,6,7}_gate_summary.json` + per-campaign reports):

- **SBC passes** (rank-uniformity, BH-FDR) in every baseline arm.
- **Limits FAIL — prior-box containment**, not monotonicity: v7 containment
  0.57 vs required 1.0; leakage concentrated at sweep Im_k2 >= 0.2 with
  Re_k2 pinned 0.23 (2000-sample draws fall to 846/307/124 inside box).
  Likely off-manifold conditioning (huge loss angle) — gate design question,
  not necessarily flow defect.
- **Crosscheck FAIL — same 2 params each time**: `D_iceIh_km` (v7: SBI mean
  62.9 vs MCMC 61.2 km, diff 1.69 vs tol 1.19; shape fail) and
  `log10_wOcean_ppt` (v7: mean diff 0.040 dex vs tol 0.026). Small absolute
  biases (~1.7 km, ~10% salinity) but beyond preregistered tolerance.

Machine B correctly stopped without tuning. NO GUI wiring of v5/v6/v7 until
Machine A scientific review adjudicates: retrain vs gate-recalibration vs
reference-MCMC scrutiny. v5 structure cache + Ae sidecar and the v7 reference
MCMC now exist on Machine A (`Test52_seawater_v5/`, `Test54_seawater_v7/`),
so A-side re-analysis is unblocked.

## NH3 activity model corrected (2026-07-28, commit 6c5ee2af)

CoolProp mixture excess had the WRONG SIGN (under-depressed liquidus 9–37%).
Replaced by Melinder-anchored Redlich-Kister on the Choukroun & Grasset P,T
shape (`NH3_ACTIVITY_MODE='melinder-CG'`, default). "L-K 2002" polynomial
DELETED (unattributable; dilute slope 2x too shallow). Consequence: Titan NH3
campaign Phase 0 MUST re-run under the corrected model; provisional Phase-1
rectangle Tb [248,257] K x w [30,100] ppt. Spec:
`plans/active/titan-nh3-ocean-campaign-spec.md`.

## Landed on genai since the 8e91d440 HF snapshot (selection)

- NH3: CoolProp NH3-H2O ocean EOS; Melinder-anchored liquidus correction +
  rewritten `tests/coolprop_nh3_test.py` (15 green).
- GUI: heating-tab radiogenic inventory selector (BSE/CI + age slider);
  Mineralogy tab (post-hoc Perple_X grain-density consistency, tolerance band
  + porosity headroom + cold-edge flags) — both AppTest-verified. The
  Perple_X native-domain density heatmap + selected-draw geotherm overlay
  (Codex C3) is `verified` 2026-08-02: AppTest + manager visual inspection
  of rendered CV3/CI figures. Titan amortized Phase-A/NH₃/MgSO₄/NaCl slot
  scaffolding (Codex C4) is `implemented, unverified`: six shipped slots and
  three awaiting-artifact states pass structural AppTest; manager browser
  confirmation remains.
- Fixes: Titan "C/MR^2 = None" message (SetCMR2strings guard + SetupInit
  call); PREM porosity-table path; NaCl `_warn_once` broadcast crash;
  `bulk_overrides` threading in `build_tbw_grid_cache`.
- Repro: frozen PPcl environment.yml; scipy 1.17.1 / numpy 2.4.6 bumps;
  v5/v6/v7 artifacts, caches, reference results, gate reports.

## Public deployment boundary

Hugging Face Space (vsteven-planetprofile.hf.space) still at `8e91d440`
(2026-07-21) — now MANY commits stale. Redeploy is a USER action (deploy
one-liner + HF_TOKEN). Nothing since the snapshot is deploy-blocked except
v5/v6/v7 slots (not wired pending adjudication).

## Awaiting

- Machine A (Claude): v5/v6/v7 gate adjudication (scientific-reviewer pass);
  then wire ratified slots; maintain queues.
- Machine A (Codex 5.6): queue empty — C1/C2/C3 are `verified`; C4 is
  `implemented, unverified` and awaits manager browser confirmation. Next
  tasks to be curated by the manager; candidate: C2's classic-MoI finding
  needs manager adjudication first.
- Machine B: Task #68 Titan JOINT no-ocean+ocean NH3 posterior in flight —
  3x4 validation cache running (~5 h) on the corrected activity model;
  builder edits (`retry_frozen_as_no_ocean` etc.) live in B's working tree,
  to be pushed after validation. Machine A reconciled the queue to the joint
  design 2026-08-02 (see `plans/STATUS-2026-08-01-machineB-joint-nh3.md` +
  `MACHINE-B-HANDOFF.md` §1): full Tb [249,263] K x w [1,70] ppt, frozen
  nodes kept in support. HOLD on any v5/v6/v7 retraining until adjudication.
- USER: HF redeploy; browser click-throughs (globe sample-picker, amortized
  body banner, selector dropdown wrap CSS — all `implemented, unverified`).

## Reference docs

- Machine B queue: `plans/MACHINE-B-HANDOFF.md`
- Codex lane: `plans/CODEX-QUEUE.md`, `AGENTS.md`
- Active/archive routing: `plans/active/README.md`, `plans/archive/README.md`
- Constraints per mission-body: `plans/mission-body-constraints.md`
- Gate reports: `validation_reports/` (+ `v5/v6/v7_gate_summary.json`)
- Artifact ledger: `PlanetProfile/Inference/sbi_artifacts/INDEX.md`
- NH3 defect trail: `plans/HANDOFF-2026-07-26-nh3-liquidus-defect.md`
