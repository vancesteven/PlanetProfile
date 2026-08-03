# PlanetProfile genai — current status

Updated: 2026-08-04 (Machine A: split ratification COUNTERSIGNED — PPC
independently reproduced end-to-end on A, SBI pushforward 0.542/0.042 ==
B's 0.541/0.042; results-panel sector warning added + slot verified
functional post-cache; PPC requirement extended to all Europa artifacts,
deployed v4/v1.1 first — see MACHINE-B-HANDOFF §0.5). Single current-state
source of truth.
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
| Titan NH₃ joint no-ocean+ocean | **SPLIT ratification 2026-08-03: gravity/structure VERIFIED, tidal (k₂) sector NOT verified**; GUI slot wired (in-browser render pending) | `titan_freegrav_nh3_posterior_1m.pt`; SBC pass; PPC shows flow under-updates k₂ (SBI pushforward at prior-pred median vs MCMC/data) — MCMC reference authoritative for k₂/dissipation, do NOT quote SBI Re_k₂/Im_k₂/ζ/η. See amended `validation_reports/titan_freegrav_nh3_1m/RATIFICATION.md`. Real slot replaced placeholder; headless load+sample verified |
| Titan MgSO₄ ocean | **awaiting artifact** | GUI placeholder; paths/centrals/gates unset |
| Titan NaCl ocean | **awaiting artifact** | GUI placeholder; paths/centrals/gates unset |
| Galileo–Europa v1.1 8D | DEPLOYED 2026-07-20 | gates 8/8 |
| Clipper–Europa v4 geodesy 11D | DEPLOYED 2026-07-21 | user-ratified |
| Clipper–Europa v5 thick-ice ablation trio | **delivered, NOT ratified — GUI-gated** | baseline gates: sbc 0, limits 2, crosscheck 2; nok2 arm sbc 2; slots render hold warning, conditioning disabled (`cb597490`) |
| Clipper–Europa v6 freegrav trio | **delivered, NOT ratified — GUI-gated** | baseline gates: sbc 0, limits 2, crosscheck 2; noinduction arm sbc 2; was silently the Clipper default — now gated, v4 default restored (`cb597490`) |
| Clipper–Europa v7 open-\|Ae\| 11D | **delivered, NOT ratified** | gates: sbc 0, limits 2, crosscheck 2; never wired into GUI |
| v1 / v2 | RETIRED | replaced by v1.1 / v4 |
| v3 | VETOED | wrong k2 bounds; 2D cache reused by v4+ |

## ADJUDICATED 2026-08-02 — v5/v6/v7 gates (record: plans/active/europa-v5v6v7-gate-adjudication.md)

Manager + independent Opus-5 adversarial review. Outcomes:

- **Limits FAILs**: superseded-gate artifact (pre-`3fa5a3fc` code demanded
  containment==1.0 pooled over off-manifold sweep points; unreachable for
  NSF flows, deploy path rejects out-of-box). Provisionally overturned —
  official regeneration at HEAD with the constructed sweep grid + anchor
  mode required (B1/B6) before any PASS is recorded.
- **v5: NOT ratifiable.** D_iceIh_km shape-EXCESS failure survives HEAD's
  own decomposition (1.12x tol, headline science param) + zeta_Ih mean
  1.074x + SBC underpowered (n=108 < 200 spec minimum).
- **v6: conditionally ratifiable.** Clean pass at HEAD (0 shape-excess,
  0 mean fails); blocked only on adequately powered SBC (n=102) + report
  regeneration (B1/B2/B6).
- **v7: crosscheck FAIL adjudicated a REFERENCE-side artifact.** True v7
  posterior provably identical to v5's at the fiducial (support strictly
  enlarged; newly admitted region carries 0.000000 posterior mass in both
  references); the 1.06 km reference disagreement has the wrong sign for a
  support effect and the rigid-translation shape of nested-sampling
  log-volume error. v7 flow passes 22/22 mean+sigma gates against the v5
  reference. Which reference wandered is OPEN (v5's log_Z_err 0.281 makes
  it the more suspect run). Blocked on B3/B4/B5.
- **Top open scientific risk:** every campaign with statistical power flags
  the same D–w degeneracy pair (v4 SBC FAIL on Tb+log10_w at n=1000; v5
  shape-excess on D_iceIh; v7 means on D_iceIh+log10_w). B2/B5 test it.
- **Governance:** deployed v4's ledger row corrected — its own reports
  contain undisclosed crosscheck mean breaches (zeta_Ih 1.23x, rho_core
  1.18x) and an SBC FAIL (Tb p=4e-4, log10_w p=0.034, n=1000). Whether v4
  needs re-adjudication is a USER decision.

GUI gating (`cb597490`) stays for all v5/v6 slots; v7 remains unwired.
Machine B follow-up queue B1–B7 in `plans/MACHINE-B-HANDOFF.md` §0.

## NH3 activity model corrected (2026-07-28, commit 6c5ee2af)

CoolProp mixture excess had the WRONG SIGN (under-depressed liquidus 9–37%).
Replaced by Melinder-anchored Redlich-Kister on the Choukroun & Grasset P,T
shape (`NH3_ACTIVITY_MODE='melinder-CG'`, default). "L-K 2002" polynomial
DELETED (unattributable; dilute slope 2x too shallow). Consequence: Titan NH3
campaign Phase 0 MUST re-run under the corrected model; provisional Phase-1
rectangle Tb [248,257] K x w [30,100] ppt. Spec:
`plans/active/titan-nh3-ocean-campaign-spec.md`.

## Titan NH3 JOINT no-ocean+ocean artifact — SPLIT ratification (2026-08-03)

Task #68 delivered: frozen (Tb,w) grid nodes build as REAL no-ocean interiors
and appear in the posterior (not rejected), per the user's 2026-07-30 decision.
Full campaign complete; reviewer status is SPLIT (gravity/structure verified,
tidal sector not — see the PPC subsection above):

- Cache `Test54_nh3_ocean/titan_nh3_joint_structure_grid_2d.pkl` (sha256
  `3d837cf8…`, `retry_frozen_as_no_ocean=True`); 1M gen 689,845 kept; reference
  MCMC 4461 samples; 1M nsf flow `titan_freegrav_nh3_posterior_1m.pt`.
- Gates: **SBC PASS** (353 pairs, min BH-adj KS p 0.772). limits FAIL is N/A for
  a joint mixture (HP ices carry tidal response in the no-ocean branch;
  documented, NOT tuned). crosscheck FAIL is **NOT benign-nuisance-only** — the
  PPC proved the per-param gate is blind to a joint tidal-sector under-update.
- **Split status:** gravity/structure (C20, C22, R_core, rho_core, salinity, Tb,
  dC20/dC22 null) VERIFIED; tidal (k2/zeta/eta) NOT verified — MCMC reference
  authoritative. Amended verdict + PPC finding in
  `validation_reports/titan_freegrav_nh3_1m/RATIFICATION.md`.
- **MgSO4/NaCl inherit the split-status discipline + the new PPC/pushforward-gate
  requirements** (NOT the old benign-nuisance framing). Committed 2026-08-03.

GUI: NH3 joint slot wired into `_SBI_ARTIFACT_SLOTS`
(`PlanetProfileApp/pages/Inference.py`) 2026-08-02 as the newest Cassini–Titan
version; widget keys already namespaced by artifact filename (no shared-key
collision). Headless load+sample verified (artifact loads at deployed path,
2000 finite draws at the slot fiducial, salinity+Tb span the full joint plane).
In-browser render (`implemented, unverified`) still pending.

**BLOCKER found 2026-08-03 (Machine A, user-reported):** the structure cache
`Test54_nh3_ocean/titan_nh3_joint_structure_grid_2d.pkl` was NOT committed —
the slot errors ("Structure cache not found") on every machine but B, and as
the newest Titan version it had become the DEFAULT, breaking the Titan model
selector for users. Repairs on `255585c0` (verified, AppTest 8/8): version
default skips slots whose artifact/cache is missing locally (NH3 stays
selectable with awaiting-cache guidance, Generate disabled); NH3 tag added to
the dropdown-visible label tail; wedge now labels the ocean with the run's
asserted composition (NH3) instead of Titan's MgSO4 body default, threaded
via the slot -> result config (training JSON/config hash untouched).
**Machine B: commit + push the Test54 cache pkl** (v5 precedent, 21 MB) —
until then the NH3 slot is non-functional on Machine A and any HF deploy.
Manager countersign of the B-side ratification + in-browser verify follow
once the cache lands.

Next in Phase B: MgSO4 (seeds 73/7373/73, salinity log10[0,2.288]), then NaCl
(74/7474/74, log10[0,2.477]) — each own config + 2D cache + artifact.

**Posterior-predictive + interior diagnostic (2026-08-03, in gate runner).**
`plans/scripts/titanG_ppc_interior_check.py` (wired into `titanG_nh3_run_gates.py`
→ inherited by MgSO4/NaCl gate scripts) pushes SBI posterior draws back through
the identical theta→x forward loop (new `theta_override` kwarg on
`generate_sbi_dataset`; forward path validated exact, max|fwd−stored|=0) and
reports posterior-predictive coverage vs x_obs + the prior-predictive percentile
of x_obs. NH3 result: x_obs interior on all 4 channels (C20 42/C22 60/Re_k2 80/
Im_k2 86 pctile) — the flow interpolates, confirming the prior verdict. BUT a
NEW finding surfaced: SBI posterior-predictive |Im_k2| median 0.042 vs the MCMC
reference 0.093 (obs 0.135) — the SBI flow under-updates the tidal sector,
traced to log10_zeta/eta_Ih/eta_III shifts that pass the per-param crosscheck by
the thinnest margins (90–92% of tol) and compound in the nonlinear k2
pushforward. The per-param crosscheck gate provably cannot see this (it never
pushes forward to k2).

**Reviewer verdict (2026-08-03): ratification SPLIT.** Gravity/structure sector
VERIFIED; tidal (k₂/dissipation) sector NOT verified — MCMC reference
authoritative, do not quote SBI Re_k2/Im_k2/zeta/eta at the Titan datum. SBC
(prior-averaged) PASS coexists with local miscalibration at the informative
datum (no paradox). RATIFICATION.md, gate_summary.json, GUI scope_note, and the
`project_joint_mixture_gate_scope` memory all amended to the split status (the
old "benign eta nuisance" framing is superseded). **MgSO4/NaCl requirements
before ratification:** run the PPC + prior-pred interior recheck; promote a
pushforward-observable crosscheck to a FIRST-CLASS gate (must flag the known NH3
k2 miss); diagnose the flow under-update before more compute (don't assume "more
epochs" fixes it). **Do not start MgSO4 compute until the flow-diagnosis
approach is chosen.** Report: `validation_reports/titan_freegrav_nh3_1m/ppc/`.
Deployed sampling path confirmed `reject_outside_prior=True` (reviewer MINOR
closed).

## Landed on genai since the 8e91d440 HF snapshot (selection)

- NH3: CoolProp NH3-H2O ocean EOS; Melinder-anchored liquidus correction +
  rewritten `tests/coolprop_nh3_test.py` (15 green).
- GUI: heating-tab radiogenic inventory selector (BSE/CI + age slider);
  Mineralogy tab (post-hoc Perple_X grain-density consistency, tolerance band
  + porosity headroom + cold-edge flags) — both AppTest-verified. The
  Perple_X native-domain density heatmap + selected-draw geotherm overlay
  (Codex C3) is `verified` 2026-08-02: AppTest + manager visual inspection
  of rendered CV3/CI figures. Titan amortized Phase-A/NH₃/MgSO₄/NaCl slot
  scaffolding (Codex C4) is `verified` 2026-08-02 (manager AppTest: Phase-A
  default, placeholder states, disjoint per-slot widget keys). Follow-on
  manager fix `cb597490` gates the unratified v5/v6 slots (adjudication
  hold enforcement; suite 7/7).
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
- Machine A (Codex 5.6): queue empty — C1/C2/C3/C4 all `verified`. Next
  tasks to be curated by the manager; candidate: C2's classic-MoI finding
  needs manager adjudication first.
- Machine B: Task #68 Titan JOINT no-ocean+ocean NH3 posterior in flight —
  3x4 validation cache running (~5 h) on the corrected activity model;
  builder edits (`retry_frozen_as_no_ocean` etc.) live in B's working tree,
  to be pushed after validation. Machine A reconciled the queue to the joint
  design 2026-08-02 (see `plans/STATUS-2026-08-01-machineB-joint-nh3.md` +
  `MACHINE-B-HANDOFF.md` §1): full Tb [249,263] K x w [1,70] ppt, frozen
  nodes kept in support. HOLD on any v5/v6/v7 retraining until adjudication.
- USER: HF redeploy ON HOLD (user 2026-08-04) until the v4/v1.1 PPC verdict
  from Machine B — it decides whether the public v4 slot needs a sector
  warning first. Browser click-throughs still pending (globe sample-picker,
  amortized body banner, selector dropdown wrap CSS, NH3 slot render — all
  `implemented, unverified`).

## Reference docs

- Machine B queue: `plans/MACHINE-B-HANDOFF.md`
- Codex lane: `plans/CODEX-QUEUE.md`, `AGENTS.md`
- Active/archive routing: `plans/active/README.md`, `plans/archive/README.md`
- Constraints per mission-body: `plans/mission-body-constraints.md`
- Gate reports: `validation_reports/` (+ `v5/v6/v7_gate_summary.json`)
- Artifact ledger: `PlanetProfile/Inference/sbi_artifacts/INDEX.md`
- NH3 defect trail: `plans/HANDOFF-2026-07-26-nh3-liquidus-defect.md`
