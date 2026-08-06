# PlanetProfile genai — current status

Big-picture program map: `plans/STRATEGY.md` (production line states,
what current work buys). This file is the tactical log.

Updated: 2026-08-06 (Machine B: §0.10 priority-(a) v5/v6 B1/B2/B6 gates EXECUTED
vs the FRESH pooled n_eff~2000 v5 reference, reviewer-adjudicated.
**v5 = PASS WITH CONCERNS** — ships for its DEPLOYED D_iceIh/ocean/salinity
deliverable with a mandatory scope-note (Machine A finalizes + countersigns):
SBC FAILs on dC22_nh ONLY (BH-adj p=0.0185; 621 kept pairs — the B2 n→1500 fix
EXPOSED a miscalibration the old n=108 masked), a real but LOCALIZED
v5-flow-specific undertraining on the ~100:1-compressed nuisance (rank-CDF mild
monotone; NOT free-gravity-family — v6's identical nuisance PASSES p=0.61);
deployed marginals doubly validated (SBC-clean D_iceIh p=0.996 + crosscheck PASS
even without the 0.36 km floor). **v6 = ALL GATES PASS.** B6 DEFERRED (empty
anchor set, below the 0.15 falsified boundary). **#4 architecture-pilot EXECUTED — capacity/embedding
ELIMINATED** (reviewer-adjudicated): uniform negative across D0/A/B/C, more
capacity made concentration WORSE. The Im_k2 gap PARTITIONS — 0.042→~0.093 is a
real flow-defect (training-signal/identifiability, NOT representational) and
~0.093→0.135 is a data-limit (obs sits in the upper tail of the EXACT posterior;
reachable frac≥obs ~0.19, not 0.5). Reviewer: MgSO4/NaCl proceed on the DEPLOYED
architecture under split-status (re-verify quarantine per-composition); next NH3
tidal diagnostic is upstream identifiability, NOT a bigger flow. FOR MACHINE A:
finalize v5 ratification + countersign; deploy v6. See MACHINE-B-HANDOFF §0.11
(v5/v6), §0.12 (#4), §0.13 (MgSO4/NaCl build prereqs). MgSO4/NaCl cache builds:
reviewer SIGN OFF on NaCl extrap_ocean (no-op≡clamp) + SIGN OFF WITH CONDITIONS
on the per-composition Tb grid; empirical onset probes run — NaCl is a ~39 K
DISJOINT ocean/frozen diagonal (→30-40 fine Tb nodes, not NH3's 12), MgSO4
~15 K [240,255] K with a REAL (fine-Pfreeze-confirmed) non-monotone frozen
island at the w=194 (2-molal) cap. HELD on reviewer adjudication of whether that
low-Tb ocean island is physical (keep) or pathological (exclude) before sizing
the MgSO4 grid + running gate-3 (half-cell Tb-shift). #4 corrected 4-arm×3-seed
re-run in flight (recovers best_val/train_val_gap; A@72/172 ~0.043, no pass —
corroborates the eliminated-capacity verdict).
(§0.10 priority-b, architecture-independent). MgSO4/NaCl build gates now CLEARED
(reviewer a3c9ed8ae24664527, 2026-08-06): MgSO4 w=194 island = PATHOLOGICAL
(Margules melting-monotonicity violation → excluded by construction + post-build
has_ocean invariant); NaCl w=290 monotone + retry-corner discriminator PASS;
MgSO4 extrap_ocean=True adopted (physics verified monotone/stiffening to 1500
MPa), ceiling raised 1200→1400 MPa to preserve the w=194 salinity cap, with
stratified eos_extrapolated flags + dρ/dP band invariant + T-spot-check wired
into the NEW committed build driver plans/scripts/titanG_build_ocean_cache.py.
Both configs authored + BOTH PRODUCTION CACHES BUILT, VALIDATED, and PUSHED for
Machine A review: NaCl (40 Tb × 15 w = 600 nodes; 314 ocean/219 no-ocean/67 None;
tilted-diagonal has_ocean map matches onset probe) at
Test54_nacl_ocean/titan_nacl_joint_structure_grid_2d.pkl; MgSO4 (17 Tb × 16 w =
272 nodes; 191 ocean/81 no-ocean/0 None) at
Test54_mgso4_ocean/titan_mgso4_joint_structure_grid_2d.pkl. All MgSO4 reviewer
guards report CLEAN in-manifest: island invariant OK (0 has_ocean=True at
w≥180 & Tb<248), extrap-ceiling OK (deepest w=194 node 1371<1400 MPa; 0 rejects),
dρ/dP∈[0.14,0.17] in all extrap columns, ocean-base T≤300 K everywhere; 45 deep +
7 mild eos_extrapolated flags on the high-w tail (carry the extrapolation caveat
in any high-salinity claim). NEXT (Machine B, no user input): per-comp gate-3
(half-cell Tb-shift → ocean-fraction, MgSO4 in the w≥150 corner), then 1M dataset
gens (seeds NaCl 74/7474/74, MgSO4 73/7373/73). Flow TRAINING still HELD on
per-composition quarantine re-verification (§0.10). Prior NH3 #1 → PERSIST
(salinity ELIMINATED) below.)
Prior: 2026-08-06 at genai `7bef8ba2` (Machine A Codex: C8 first
doc-doctor pass `verified`; report `94ded63f` covers all 12 checklist items,
2 PASS / 10 FINDING, with no findings fixed. Codex queue empty; manager triage
next).
Prior: 2026-08-05 (Machine B: NH3 follow-ups #3, #2 AND #1 EXECUTED, reviewer
PASS on all. #1 salinity-fixed pilot + matched-N control → PERSIST: fixing
salinity at matched N did NOT recover the Im_k2 update (banded 0.0313 ≈ control
0.0321, +0.02σ_obs, wrong sign for salinity); size effect small/monotone/toward-obs
→ salinity axis ELIMINATED as the driver of gap (a); remaining candidate is
capacity/embedding (#4). STOPPED per protocol — remedy selection is manager +
reviewer + user. MgSO4/NaCl stay HELD. See below).
Earlier: 2026-08-05 at genai `5ae7d17c`. Manager judgement §0.9 issued:
B3 accepted — the 1.06 km v5–v7 reference disagreement was an n_eff=500
resolution artifact; both old references wandered (fresh truth ~61.4–61.8 km);
§0.7 crosschecks retarget the pooled fresh n_eff=2000 references with a
preregistered 0.36 km shape-excess floor, and the earlier "v7 flow 22/22 vs
v5 reference" claim is corrected (void — wandered target). NH3 separators
accepted: joint mixture and noise-swamping ruled out for the Im_k2 miss;
follow-ups reordered #3 (pushforward plot) → #2 (matched-resolution NH3
reference — the 0.100 MCMC-pp ceiling inherits the discredited resolution
class) → #1 only if a material gap survives; MgSO4/NaCl release criteria
preregistered on the #2 outcome (details MACHINE-B-HANDOFF §0.9).
**#3 EXECUTED 2026-08-05 (reviewer PASS):** pushforward-shape diagnostic on
BOTH the capped anchor AND the deployed 1M flow confirms concentration-failure
and EXCLUDES wrong-mode — deployed concentration ratio 1.215 (>1: posterior
not narrower than prior), high-k2 tail retained (frac≥ceiling 0.233), but
frac≥obs moved only 0.128→0.161 (ceiling-independent under-update evidence).
Nonfinite-drop bias closed (<0.25% of draws, conservative direction).
**#2 EXECUTED 2026-08-05 (reviewer PASS):** matched-resolution NH3 reference MCMC
(n_eff=2000/n_active=1024, seeds 72+172, config unmutated) re-measured the
weighted |Im_k2| ceiling = **0.1037** (n_eff=500 was 0.0999; +0.0038 = +0.11σ_obs
and *away* from SBI — the ceiling is NOT an n_eff=500 artifact, the decisive
contrast with B3). SBI-pp 0.0423 vs matched MCMC-pp 0.1037 → **gap +1.76σ_obs ≫
0.5σ → SURVIVES**. With #3 (concentration-failure, wrong-mode excluded) the
under-update is a genuine flow-training deficiency. **Decision: run #1
(salinity-fixed retrain) with reviewer+user sign-off; MgSO4/NaCl PROCEED
falsified, stay HELD.** Reviewer note: matched MCMC-pp itself sits 0.89σ below obs
(ordinary model/data tension, NOT a remediation target); #1 targets only the
1.76σ SBI-vs-MCMC gap. See `validation_reports/FLOW_UNDERUPDATE_DIAGNOSIS.md`
(Follow-ups #3 + #2) and `validation_reports/nh3_diagnosis/matched_reference/`.
**#1 EXECUTED 2026-08-05 (reviewer PASS): PERSIST → salinity ELIMINATED as the
driver.** Design-review reviewer (PASS WITH CONCERNS) required a sample-size-matched,
salinity-VARYING control (the fixed-salinity band keeps only ~9%/60,039 rows vs
~690k — an ~11× cut biases toward under-concentration on its own); folded in +
band re-documented as a symmetric ±0.084-dex slice about 12.6 ppt. Results (Im_k2
pp-median): banded(N=60k,fixed) 0.0313 ≈ control(N=60k,varying) 0.0321 → salinity
axis IN vs OUT = +0.02σ_obs, wrong sign for the salinity hypothesis; control(60k)
vs S1(643k varying) 0.0392 = −0.20σ (pure size effect, toward obs, discharges the
confound); banded−ceiling(0.1037) = +2.07σ still-open flow gap. Reviewer PASS:
salinity cleared for gap (a) [SBI-pp vs MCMC-pp ceiling], remaining candidate
capacity/embedding (#4); size gain strongly sublinear (11× data → +0.007; ceiling
needs +0.065) so more ocean-only data won't close it; banded pp 5–95 [0.0023,0.340]
spans past obs → expressiveness/concentration signature, not support failure.
STOPPED per protocol; #4 remedy selection is manager + reviewer + user.
`validation_reports/nh3_diagnosis/f1_salinity_fixed/`. Also
2026-08-05: e-notation display fix for C20/C22-scale parameters in the
corner plot + summary table (`12bdbf56`, verified). Earlier-entry history:
`git log -p -- plans/STATUS.md`.
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
| Galileo–Europa v1.1 8D | DEPLOYED 2026-07-20 | gates 8/8; PPC 2026-08-04: 0/5 flagged, no tidal warning (tidal channels non-informative by construction) |
| Clipper–Europa v4 geodesy 11D | DEPLOYED 2026-07-21 | user-ratified; PPC 2026-08-04: 0/21 flagged, no tidal warning (Re_k2 informative extreme-tail datum, SBI tracks MCMC 1.43σ update to 0.04σ — decisive clean evidence) |
| Clipper–Europa v5 thick-ice ablation trio | **delivered, NOT ratified — GUI-gated** | baseline gates: sbc 0, limits 2, crosscheck 2; nok2 arm sbc 2; slots render hold warning, conditioning disabled (`cb597490`). PPC 2026-08-04 (baseline): 0/21 flagged, non-trivially clean (C20 22.6σ update tracked 0.12σ); flow-fidelity only, does NOT clear ratification blocks |
| Clipper–Europa v6 freegrav trio | **delivered, NOT ratified — GUI-gated** | baseline gates: sbc 0, limits 2, crosscheck 2; noinduction arm sbc 2; was silently the Clipper default — now gated, v4 default restored (`cb597490`). PPC 2026-08-04 (baseline): 0/20 flagged (gravity loosened by design 0.4–0.5σ, Re_k2 1.66σ + synodic 2.8–3.1σ tracked ≤0.07σ); flow-fidelity only |
| Clipper–Europa v7 open-\|Ae\| 11D | **delivered, NOT ratified** | gates: sbc 0, limits 2, crosscheck 2; never wired into GUI. PPC 2026-08-04: 0/21 flagged (Bind_synodic_x_real 52.9σ update tracked 0.18σ); gaps batch-largest 0.10–0.27σ carried to B5; PPC orthogonal to the B3 parameter-space reference-wander block |
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
  hold enforcement; suite 7/7). Perple_X cold-edge band/legend polish (C6)
  is `implemented, unverified`; CV3/CI PNGs and AppTest evidence are recorded
  in `plans/CODEX-QUEUE.md`.
- Fixes: Titan "C/MR^2 = None" message (SetCMR2strings guard + SetupInit
  call); PREM porosity-table path; NaCl `_warn_once` broadcast crash;
  `bulk_overrides` threading in `build_tbw_grid_cache`.
- Repro: frozen PPcl environment.yml; scipy 1.17.1 / numpy 2.4.6 bumps;
  v5/v6/v7 artifacts, caches, reference results, gate reports.

## Public deployment boundary

Hugging Face Space (vsteven-planetprofile.hf.space) is CURRENT at genai
`d53385f1`, deployed 2026-08-04 (user ran the upload one-liner; USER
confirmed the Space running). The batch replaces the 2026-07-21 `8e91d440`
snapshot and adds: Titan NH3 joint slot (+ its cache, sector warning, wedge
NH3 composition), v5/v6 not-ratified GUI gating with v4 as Clipper default,
mineralogy tab + Perple_X geotherm overlay, heating radiogenic-inventory
selector, ready-slot default fallback, NaCl `_warn_once` fix, and the NH3
activity-model correction underneath. Deploy-script cache list now includes
Test54 (`d53385f1`). C5 subsequently made future deploy snapshots derive and
exact-diff the full cache set from `_SBI_ARTIFACT_SLOTS`; build-only evidence
is `verified`, but no new deploy was pushed.

## Awaiting

- Machine A (Claude): v5/v6/v7 gate adjudication (scientific-reviewer pass);
  then wire ratified slots; maintain queues.
- Machine A (Codex 5.6): queue empty — C8 first doc-doctor pass is `verified`
  at `validation_reports/doc_doctor/2026-08-06_first_pass.md` (12/12 items,
  2 PASS / 10 FINDING; no fixes). Findings await manager triage into scoped
  follow-up tasks.
- **Machine A (Claude) — DECISION REQUESTED: NH3 flow under-update diagnosis
  is complete (`711c2fd2`); MgSO4/NaCl compute is HELD pending your call.**
  Diagnosis (reviewer PASS WITH CONCERNS, corrections folded): #3 x-norm scale
  RULED OUT; obs-vector width RULED OUT as causal (no-ocean control assimilates
  informative Re_k2 3.92σ→0.05σ, 0/4 flagged); mechanism localized to the
  ocean-admitting apparatus (joint mixture AND the co-varying salinity axis —
  NOT the mixture alone, since the control also adds a 13th param + removes the
  phase guard); bimodal mode-assignment is a CANDIDATE, not established; #1
  noise swamping cleared for Re_k2 but UNTESTED for the dominant Im_k2 miss.
  Three paths for you/the user to pick, cheapest-first: (1) run two cheap
  SEPARATORS before any build — an ocean-only-with-salinity PPC (splits salinity
  axis from bimodality) + the #1 reduced-noise Im_k2 pilot; (2) proceed with
  MgSO4/NaCl under split-status as-is (tidal sector MCMC-authoritative,
  pushforward gate fires by construction); (3) remediate the conditional first
  (mixture-aware embedding / has_ocean-labelled conditional / sequential round).
  Full report `validation_reports/FLOW_UNDERUPDATE_DIAGNOSIS.md`; rationale in
  `MACHINE-B-HANDOFF.md` §0.6 RESOLVED block. Machine B awaits the decision;
  will not start #1 pilot or any MgSO4/NaCl work unilaterally.
- Machine B: idle after the diagnosis push. HOLD on any v5/v6/v7 retraining
  until adjudication; MgSO4/NaCl HELD pending the decision above.
- USER: live-app click-throughs DONE 2026-08-05 — one change requested and
  shipped (`01649f3f`, verified): Europa induction conditioning defaults
  Re Ae = 0.9 for the synodic band + second harmonic (imaginary parts and
  the orbital band keep config-derived values; inputs stay editable). The
  Space runs the d53385f1 snapshot; the Ae-default change reaches it at the
  next redeploy. v5/v6 ratification (user priority) is executing via
  MACHINE-B-HANDOFF §0.7 — B3 done, B1/B2 against the fresh references
  next.

## Reference docs

- Machine B queue: `plans/MACHINE-B-HANDOFF.md`
- Codex lane: `plans/CODEX-QUEUE.md`, `AGENTS.md`
- Active/archive routing: `plans/active/README.md`, `plans/archive/README.md`
- Constraints per mission-body: `plans/mission-body-constraints.md`
- Gate reports: `validation_reports/` (+ `v5/v6/v7_gate_summary.json`)
- Artifact ledger: `PlanetProfile/Inference/sbi_artifacts/INDEX.md`
- NH3 defect trail: `plans/HANDOFF-2026-07-26-nh3-liquidus-defect.md`
