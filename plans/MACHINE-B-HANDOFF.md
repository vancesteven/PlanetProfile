# Machine B handoff

## 0.16 CAMPAIGN COMPLETE — MgSO4/NaCl CLEARED for split-status deploy (Machine B, 2026-08-10)

**The Titan free-gravity MgSO4/NaCl compute campaign is DONE. Deploy is CLEARED.
Remaining work is Machine A GUI wiring (task #52) + two caveat-copy corrections.**

Full §0.16 gate sequence executed and adjudicated by the scientific-reviewer in
two binding passes:
- Gate-set interpretation [agent a719438103ced11d0]: PASS WITH CONCERNS, no STOP.
  SBC PASS both; limits containment PASS (monotone = gated-FAIL-by-override, NOT
  N/A); crosscheck MgSO4 PASS all 13 params, NaCl FAIL on log10_eta_V only
  (non-blocking poorly-identified nuisance).
- Re_k2 pushforward [agent ad7ae4436e53cdc1e]: **PASS — BOTH comps deploy under
  full split-status.** MCMC→Re_k2 pushforward (pooled repaired posterior →
  identical θ→x loop) shows both comps UNDER-PREDICT Re_k2 (MgSO4 0.570, NaCl
  0.575 vs obs 0.608), neither centers on obs ⇒ model-data TENSION reproduced by
  the independent sampler, NOT a flow offset ⇒ Re_k2 ACCEPT. Gravity centers on
  obs (≤0.05σ); Im_k2 under-predicts (consistent with the k2 quarantine).

**Split-status to deploy: C20/C22 TRUSTED · Re_k2 informative-with-caveat ·
tidal k2 QUARANTINED.** Artifacts: `titan_freegrav_{mgso4,nacl}_posterior_1m.pt`.

**Machine A prerequisites for Phase C GUI slot wiring (reviewer-required copy
corrections; physics already cleared):**
1. Re_k2 caveat MUST quote the median-to-median MCMC-vs-SBI gap (MgSO4 0.24σ,
   NaCl 0.53σ) — NOT the mixed deviation-of-median vs median-of-abs-dev statistic
   (SBI ppc's `median_dev_sigma` is spread-inflated; not comparable to the MCMC
   signed dev). The verdict is unaffected (the decision rule already used the
   correct median-to-median gap), but the human-facing caveat must not misstate it.
2. Caveat MUST state the deployed SBI Re_k2 marginal is CONSERVATIVE relative to
   the reference MCMC (tension-leaning, not a tightened bound) — safe-side for an
   informative-with-caveat channel; it errs toward tension, never false agreement.
Reports: `validation_reports/titan_freegrav_{mgso4,nacl}_1m/rek2_pushforward/`.

---

## 0.17 TRACK 1 IS-CORRECTION VALIDATION (Machine A -> Machine B, 2026-08-11)

Context: plans/active/tidal-sector-remedy-plan.md (reviewer-approved,
C1-C16). The IS correction lifts the tidal quarantine by reweighting flow
draws with the exact MCMC likelihood; Machine A's NH3 pilot at N=20k:
corrected |Im k2| pushforward 0.1084 vs MCMC matched ceiling ~0.1037
(flow alone 0.045), Pareto-k clean, ESS 1240, all module gates pass
except the ESS/N fractional floor (under reviewer adjudication —
proceed; it does not block the reference-side work below).

Machine A shipped: PlanetProfile/Inference/is_correction.py (tests
tests/is_correction_test.py 9/9) and the committed driver
plans/scripts/is_correction_validate.py.

**Machine B tasks (priority after current queue):**
1. For each composition (NH3 first, then MgSO4, NaCl): run
   `python plans/scripts/is_correction_validate.py --comp <c>` on the
   machine holding the pooled reference pkls. The driver auto-runs the
   reference-dependent gates: C3 likelihood-recompute consistency
   (<1e-9), pushforward median-to-median + weighted-KS (C10), C20/C22
   no-regression (C11), C16 ocean-fraction vs reference.
2. BEFORE reading any corrected result: compute and commit each
   reference's 3-seed ocean-fraction spread and the 3-seed crosscheck
   median spread for the unidentified nuisances (C9/C16 preregistration).
3. Crosscheck gate: feed samples[resampled_indices.npy] (L=int(ESS),
   committed by the driver) to the ratified validate_sbi crosscheck as
   the corrected-SBI set; report with ESS, never N.
4. C13: repeat step 1 at --seed-offset 1000 and 2000; corrected
   pushforward medians must agree within 0.1 sigma_obs.
5. C12 amortized sweep (after 1-4 pass at the fiducial): >=200
   prior-predictive x + 8 axis endpoints, Pareto-k <= 0.7 for >=95%.
6. SBC of the corrected pipeline (C14/C15): rank definition = weighted
   ranks r = sum_i w_i 1[theta_i < theta_0] vs Uniform(0,1); theta_0
   from the EFFECTIVE prior (all training support cuts); n>=500 pairs at
   the deploy-time N; BH-FDR across 13 params. Budget ~14 CPU-h/comp.
Ship all reports; manager re-adjudicates and only then touches GUI
warnings. Europa v5/v6 run as controls afterwards (v5 dC22_nh SBC FAIL
should be REPAIRED by a valid correction — positive control).

## 0.16 FLOW TRAINING AUTHORIZATION (manager, 2026-08-10) — conditional GO for MgSO4/NaCl

Answering "ready for flow training?" against the two §0.15 holds:

**Hold 1 — #4 architecture-pilot verdict: SATISFIED IN SUBSTANCE, one
report outstanding.** The pilot verdict is in and reviewer-adjudicated
(§0.12: capacity/embedding ELIMINATED; MgSO4/NaCl proceed on the DEPLOYED
architecture). The corrected 4-arm×3-seed re-run was the last insurance
item (closes the premature-early-stop alternative); its partials
corroborated (A@72/172 ~0.043, no pass) but the final report never landed.
Condition (a) below covers it.

**Hold 2 — per-composition quarantine re-verification: RESTRUCTURED, not a
pre-training hold.** Re-verifying the tidal-sector quarantine for a
composition requires a TRAINED flow (pushforward vs reference MCMC — the
NH3 method). It cannot precede training; treating it as a pre-training
hold deadlocks. It becomes the MANDATORY FIRST post-training gate: each
composition trains under split-status by default (tidal k2 sector
quarantined, exactly the NH3 posture) and the quarantine is
confirmed-or-lifted per composition from its own pushforward evidence.
Do NOT port the NH3 verdict by assumption in either direction.

**Therefore: TRAINING GO for both compositions, conditions:**
(a) Commit the #4 corrected 4-arm×3-seed final report (with recovered
    best_val/train_val_gap) no later than with the first training
    deliverable. If it OVERTURNS the pilot (large arms genuinely
    premature-early-stopped and a bigger arm passes), STOP training and
    escalate — do not proceed on the deployed architecture.
(b) Datasets: /tmp copies are void after any reboot — regenerate with the
    committed driver + seeds and verify cache SHAs against the committed
    gen manifests (NaCl 0fdbd44f…, MgSO4 124c8539…) before training.
(c) Architecture + recipe: the DEPLOYED architecture and training recipe
    exactly (nsf via sbi 0.26.1, deployed early stopping, z-scoring with
    the fixed convention); seeds and package versions recorded in the
    artifact manifest. No architecture experimentation in production
    training — that path goes through the upstream-identifiability
    diagnostic (§0.12 ruling 3), not these campaigns.

**Start NOW in parallel (independent of training): fresh reference MCMC
per composition** on {C20, C22, Re_k2, Im_k2}, B3 protocol (n_eff~2000
class, >=3 seeds, pinned env, pooled + renormalized) — these are the
crosscheck targets and the quarantine-re-verification baseline.

**Post-training gate sequence per composition (preregistered, no tuning
after seeing failures):** (1) pushforward four-way table + tidal-sector
quarantine verdict (0.5 sigma_obs flag, NH3 method); (2) SBC at
n_sbc=1500; (3) crosscheck vs the pooled fresh reference (standard
tolerances; shape clause; SHA-recorded gate code); (4) limits per the
ratified in-support/materiality rules. Ship gate summaries +
per-arm reports; manager adjudicates ratification and GUI wiring
(placeholders already live).

**MACHINE B EXECUTION (2026-08-10, commit 1fb0af97).** Acted on the §0.16 GO
with "everything in parallel" (user-confirmed sequencing). Preconditions
verified: /tmp datasets SURVIVED (no reboot); cache SHAs MATCH committed gen
manifests (NaCl 0fdbd44f, MgSO4 124c8539) + npz SHAs match; the
`ResetNearestExtrap` dtype fix (9d76375e) is orthogonal (future builds only —
installed bytes/datasets unchanged). LAUNCHED four background jobs (thread-pinned,
separate processes — libomp): NaCl+MgSO4 flow training (`titanG_ocean_train_all.py
--comp`, DEPLOYED nsf arch only, seeds 74/73, cache-sha assertions passed) and
NaCl+MgSO4 B3 reference MCMC (`titanG_ocean_reference_mcmc.py --comp`,
n_eff=2000/n_active=1024, 3 seeds each 74/174/274 & 73/173/273, pooled). Both
new drivers committed.
- **Condition (a) — #4 corrected report:** reconstructed
  `f4_pilot_manifest_corrected_multiseed.json` from recovered per-seed outputs
  (aggregator died before writing it; 10/12 arm×seed: D0/A/B×{72,172,272}, C×72).
  best_val/train_val_gap None (sbi 0.26.1 keys didn't populate — recorded as a
  harness limitation, NOT fabricated); epochs+hit_ceiling recovered: NO arm hit
  the 500-epoch ceiling, A/C converged in 40/45 epochs vs D0's 124. All arms
  pp_imk2_median ~0.043–0.047 << 0.0862 bar; D0 stable across 3 seeds. Machine B
  reading: NOT OVERTURNED (corroborates §0.12). Scientific-reviewer dispatched
  for the binding non-overturn call (agent a909e4b1d07ba13ac) — if it OVERTURNS,
  training STOPS + escalates per §0.16(a).
- **Condition (a) CLEARED (reviewer verdict, 2026-08-10, commit 1a993b56):**
  scientific-reviewer PASS — does NOT overturn the eliminated-capacity verdict;
  MgSO4/NaCl training on the DEPLOYED architecture is CLEARED to continue. Basis:
  all 10 recovered pushforward medians match per-seed reports exactly; max arm×seed
  = A_seed272 = 0.0471 (< 0.05, ≪ 0.0862 bar); concentration_ratio > 1 (broadening)
  for all 10; A/C converged in FEWER epochs than D0 with hit_ceiling=False — opposite
  of premature early-stopping, so the §0.16 STOP trigger is unmet on both clauses.
  None best_val is a real harness gap but logically non-material to a continue call.
  MODERATE provenance fix applied (epoch source is the screen manifest, NOT
  f4_full.log which has zero epoch data) + MINOR ceiling-label/172,272-unavailable
  notes recorded. Reviewer archival follow-ups (complete C_seed172/272; fix
  sbi_runner.py summary filter so best_val populates; re-emit ≥1 A/C seed) are
  NON-blocking, deferred to before the architecture decision is FINALIZED (distinct
  from this continue-training call). §0.16(a) binding gate now SATISFIED.

## 0.15 GATE-3 VERDICT + shared Tb=252 K defect root-cause (Machine B, 2026-08-07)

The scientific-reviewer adjudicated gate-3 (agent a070eb8cc673d6e76). Verdict:
- **NaCl: PROCEED-to-1M (PASS WITH CONCERNS).** 2.2% half-cell placement band
  acceptable; the binding gate is the ABSOLUTE onset placement `half_cell_K` (NaCl
  0.5 K, MgSO4 0.25 K) — NOT the box-fraction span, which inverts the true ranking
  (penalizes MgSO4's narrower box despite its FINER grid). Bar: `half_cell_K ≤ 0.5 K`
  (both pass). Also requires: box-wide melting-monotonicity (0 FROZEN warmer than any
  OCEAN in a w-column) + no interior None. The 5% MgSO4 fraction is an irreducible
  narrow-box artifact (`span = dTb/Tb_box`), NOT disqualifying.
- **MgSO4: BLOCK 1M as-is.** Not for the 5% fraction — for a CRITICAL embedded
  melting-monotonicity violation at (Tb=252.0, w=4.857): ocean 251.5 → frozen 252.0 →
  ocean 252.5. Same impossibility class as the w=194 Margules island but at
  near-pure-water salinity; the w≥180 invariant does not catch it.

**Root cause (Machine B, task #90 — CLOSED).** The reviewer's flagged MgSO4 flip and
the MAJOR shared failure (NaCl interior None at the IDENTICAL Tb=252.0 K, w-col-3) are
**one bug: an int-vs-float numeric-type artifact**, not a physics failure. PPTitan
ships `Planet.PfreezeUpper_MPa = 230` as a Python **int**; the whole cache built with
it. `fn_phase(230, Tb=252)` returns phase 0 (liquid) with **int** 230 but phase 3
(ice III) with **float** 230.0 — the int result routes `GetPfreeze` to `GetZero`
bracketing on a discontinuous phase step → fails → `NoIceLiquidTransitionError` →
MgSO4 phantom-frozen / NaCl None. Float coercion repairs EXACTLY the two defect nodes
(MgSO4 252.0→OCEAN D=15.60, NaCl 252.0→OCEAN D=31.06), preserves the genuine frozen
onset, monotonic D_ocean between neighbors, neighbor likelihood scalars uncontaminated
(reproduced in fresh per-node processes). Evidence:
`validation_reports/titan_saltcaches/tb252_rootcause_2026_08_07.json` (commit b09e534f).

**Reviewer PROCEED-to-rebuild (agent ab2d3748d0b00c78c) + mechanism CORRECTION.**
The reviewer confirmed the physical conclusion and endorsed the float-coercion rebuild
but corrected the mechanism: the int/float divergence is NOT in `fn_phase`
(`Nearest2DInterpolator`/`searchsorted` are type-invariant) — it is int-dtype
truncation in **`ResetNearestExtrap` (DataManip.py:15-35)**. Its size-1 branch (line 26
`outVar1=np.array(var1)`) preserves the int64 dtype; assigning the float melt-EOS bound
229.96 into that int array truncates to 229. Since the melt grid Pmax=229.96 <
PfreezeUpper=230, the search endpoint always clamps: int→229 (liquid) → failing GetZero
bracket → NoIceLiquidTransitionError; float→229.96 (ice III) → brute-force scan finds
the real Ih→liquid transition. Verified directly (`ResetNearestExtrap(230,…)`→229.0 vs
`(230.0,…)`→229.96; generalizes to T). Evidence JSON updated (commit 9c378eaf).

**REBUILT + INSTALLED + gate-3 BINDING PASS (2026-08-07).** Both caches rebuilt with
mandatory float-coerced `PfreezeUpper_MPa=230.0` + the sandwich invariant:
- **MgSO4**: 192 ocean / 80 no-ocean / 0 None (+1 ocean = repaired 252.0 node);
  Tb=252.0 row all-ocean; D_ocean 15.603 monotonic (15.572<15.603<36.236); CMR2
  smoothed 0.31894→0.31920; genuine frozen onset (250.0/250.5/251.0) preserved.
- **NaCl**: 315 / 219 / 66 (+1 ocean −1 None = repaired 252.0 node); D_ocean 31.056
  monotonic (20.691<31.056<106.614); genuine frozen onset (250.0) preserved.
- Gate-3 (absolute `half_cell_K` primary): **MgSO4 0.250 K PASS, NaCl 0.500 K PASS**;
  0 sandwich violations each; fraction-span diagnostic 0.0500 / 0.0218 (confirms the
  reviewer's ranking inversion — MgSO4 is the physically finer cache). Reports:
  `validation_reports/titan_freegrav_{mgso4,nacl}_1m/gate3_boundary.json`.
- Installed to Test/ (old bytes preserved in git history). Manifests carry
  `pfreeze_upper_float_coerced=True` + `sandwich_invariant_violations=[]`.

**GATE-3 RE-ADJUDICATED ON REBUILT BYTES — both PASS (reviewer agent a02a9038b2de1d382,
2026-08-07).** Machine A gave the go-ahead for gate-3; gate-3 was re-run fresh on the
installed rebuilt bytes and the scientific-reviewer independently reconstructed both
classification grids, onset structure, mass balance, and the 252 K repair from the
`.pkl`s (not manifest-trust). Verdict: **NaCl PASS, MgSO4 PASS — gate-3 satisfied on the
rebuilt bytes; no blocking issue for flow training.**
- **NaCl at the bar is a PASS.** Tb grid uniform 1.0 K (`np.unique(diff)=[1.0]`) → every
  onset cell 1.0 K → half_cell=0.500 exactly at the 0.5 K bar; non-strict `≤` is the
  correct reading, and ocean-fraction diagnostic span 0.0218 confirms the split is not
  grid-artifact-dominated. Tb refinement to 0.5 K near the 247–252 K onset band is
  RECOMMENDED (non-blocking) only if NaCl posterior ocean-fraction becomes a headline.
- **MgSO4 w≥150 corner diagnostic FAIL is benign (CONFIRMED).** Corner spans 4 columns
  (w=160/175/185/194), all onset Tb=249.5 K with 0.5 K cell; span 0.05 is a
  coarse-denominator artifact (one-row shift ≈ 1/16 of corner mass), exactly the demoted
  metric — not a phantom-ocean leak. Binding corner half_cell 0.25 K passes.
- **MgSO4 zero-None is correct, NOT a masked failure.** Prior box entirely buildable;
  retry fires only on the typed `NoIceLiquidTransitionError` (not thin-shell/PHydroMax).
  Reviewer inspected a frozen MgSO4 node (Tb=248, w=4.86): full HP-ice stack
  (D_iceIh=159, III=82.8, V=154.8, VI=207.9 km), conserved Mtot=1.345e23, CMR²=0.3193.
  NaCl's 66 None sit in the genuine high-w/high-Tb corner where even the no-ocean retry
  fails — correctly excluded from support.
- 252 K repair independently verified: MgSO4 (252,4.86) D_ocean=15.60 monotonic; NaCl
  (252,4.13) D_ocean=31.056 monotonic. Prior CRITICAL flip resolved.

**MACHINE A — please re-verify the NEW cache bytes** (your §0.14 acceptance was for the
pre-fix bytes; node counts changed by +1 ocean each). **The `ResetNearestExtrap`
int-truncation bug is a LATENT shared-thermodynamics defect** (affects ANY scalar-int
P/T query beyond an EOS domain, not just PfreezeUpper) — the reviewer's smallest fix is
`dtype=float` in the size-1 branch (lines 26/28). The gate-3 reviewer re-confirmed this
is out of gate-3 scope and non-blocking for these caches (sandwich invariant + sampled
monotonic D_ocean exclude any surviving discrete int-truncation artifact), but flagged
that any future PP cache built from an int-typed `PfreezeUpper`/T-clamp template must
carry the same coercion until the `dtype=float` fix lands. Referred to you; NOT
self-adjudicated.

**MACHINE A CLOSEOUT (2026-08-07): both asks executed, `verified`.**
1. Rebuilt cache bytes RE-ACCEPTED (supersedes §0.14 which covered pre-fix
   bytes): independent pkl reload — MgSO4 192/80/0, NaCl 315/219/66 (exact),
   pkl SHA-256 prefixes 124c8539 / 0fdbd44f match the committed gen
   manifests, ZERO sandwiched frozen nodes in any w-column (both comps),
   repaired 252.0 K nodes confirmed ocean with D_ocean 15.60 / 31.06 km.
   §0.14's acceptance now applies to these bytes; 1M datasets generated
   from them are accepted by extension (SHA match).
2. `ResetNearestExtrap` dtype fix ADJUDICATED + LANDED: numerical-
   robustness repair, not a scientific-assumption change (docstring
   already declares float contract). `dtype=float` in both size-1
   branches (`PlanetProfile/Utilities/DataManip.py`), regression test
   `tests/reset_nearest_extrap_test.py` (int-scalar == float-scalar at
   the exact 230-vs-229.96 defect geometry, T-axis analog, in-domain and
   array paths; 4/4 pass, MgSO4/NaCl consumer smokes 7/7). Future cache
   builds no longer need the float-coercion workaround at HEAD; keep the
   coercion in build drivers anyway (harmless belt-and-suspenders for
   older checkouts).

**1M gens: DONE (2026-08-07).** Both datasets generated on the rebuilt caches and
validated:
- NaCl: 631,214 kept/1M (36.9% reject), seeds 74/7474, 101.5 min, cache sha
  0fdbd44f… MATCHES installed Test/ bytes.
- MgSO4: 691,075 kept/1M (30.9% reject = NH3 precedent 31%), seeds 73/7373, 112.3 min,
  cache sha 124c8539… MATCHES installed Test/ bytes.
Datasets (~83/91 MB) stay in `/tmp/titanG_build/datasets/titanG_{nacl,mgso4}_1m.npz`
(ephemeral — regenerable from committed seed+cache); gen manifests committed at
`validation_reports/titan_freegrav_{nacl,mgso4}_1m/gen_manifest.json` (commit 099b8429).
Driver `plans/scripts/titanG_ocean_gen_dataset.py --comp {NaCl,MgSO4}` (commit 90edf5c2).

**Flow TRAINING remains HELD** on the #4 architecture-pilot verdict + per-composition
quarantine re-verification (§0.10/§0.12). Next once training unblocks: train each flow
(separate process from PP gen — libomp), fresh reference MCMC per composition on
{C20,C22,Re_k2,Im_k2}, then SBC/crosscheck/limits + pushforward gates. NOTE: a reboot
wipes the /tmp datasets — regenerate with the committed driver+seeds before training if
this machine restarts.

## 0.10 PRODUCTION AUTHORIZATION (user + manager, 2026-08-06) — MgSO4/NaCl proceed in parallel with the architecture pilot

The user ratified the parallel option (STRATEGY.md). Execute:

1. **MgSO4/NaCl production STARTS NOW at the architecture-independent
   stages:** build the 2D joint structure caches (same
   `retry_frozen_as_no_ocean=True` joint design as NH3; per-composition
   salinity ranges from the campaign spec; ocean comp recorded in cache
   metadata) and generate the training datasets. Compute priority remains
   (a) v5/v6 B1/B2/B6 first, (b) these caches/datasets, in that order.
2. **#4 architecture pilot in parallel (existing NH3 dataset, no new
   sims):** test whether increased flow capacity and/or an embedding
   change closes the Im_k2 under-update (target: SBI-pp within 0.5
   sigma_obs of the matched MCMC-pp ceiling 0.1037). Design freeze with
   the scientific reviewer before running; report either outcome. This is
   a PILOT — the production architecture decision is made by manager +
   reviewer + user with its result.
3. **MgSO4/NaCl flow TRAINING waits for the pilot verdict** (datasets do
   not). If the pilot closes the gap, the improved architecture becomes
   the production architecture (formal sign-off will be recorded); if
   not, MgSO4/NaCl train on the standard architecture under the
   split-status discipline (structure sector primary; tidal sector
   quarantined with the sector warning; pushforward gate mandatory).
4. **Enceladus production enters the queue after v5/v6 validation ships**
   — §2 below is the standing spec; confirm the config freeze with
   Machine A before the production dataset.

## 0.14 MANAGER CACHE REVIEW (Machine A, 2026-08-07) — MgSO4/NaCl caches ACCEPTED; gate-3 + dataset gens GO

Machine A re-verified both committed caches against the manifests and
reviewer guards (independent pkl load, not manifest-trust):
- Node counts EXACT match: MgSO4 191 ocean / 81 no-ocean / 0 None of 272;
  NaCl 314 / 219 / 67 of 600. `ocean_comp` present in cache metadata
  (campaign-spec traceability), schema v3.0, `retry_frozen_as_no_ocean=True`
  (NH3 precedent).
- MgSO4 deepest liquid-layer P = 1371.1 MPa at (Tb=258, w=194) — matches
  B's 1371 claim, under the 1400 MPa ceiling; w=194 extrap column ocean
  density monotone nondecreasing with P (spot check Tb=249.5,
  ρ 1280.8→1462.9). Island/ceiling/dρdP/hot-base guard lists all empty
  in-manifest.

**Verdict: both caches ACCEPTED. Per-comp gate-3 (half-cell Tb-shift) and
the 1M dataset generations are GO** (seeds NaCl 74/7474/74, MgSO4
73/7373/73 as recorded). Flow TRAINING stays HELD on per-composition
quarantine re-verification (§0.12 ruling).

Two record items for Machine B (non-blocking, fold into gate-3 pass):
1. MgSO4 has ONE isolated frozen node at (Tb=252.0, w=4.8567) surrounded
   by ocean neighbors on all sides (has_ocean row 'OOOFOOOO...') — a
   single-node Tb non-monotonicity far from the excluded w>=180 island.
   Verify gate-3's half-cell Tb-shift covers this cell; if it flips
   ocean, record it as borderline-liquidus numerics in the manifest.
2. NaCl acceptance evidence (/tmp/nacl_monotonicity_check.json,
   /tmp/nacl_corner_discriminator.json, /tmp/titan_tb_probe_results.json)
   lives only on Machine B — COMMIT these under
   validation_reports/titan_saltcaches/ (doc-doctor item E11: no
   provenance link may exist only in /tmp or chat history). Also record
   the reason class for the 67 NaCl None nodes in the manifest.

## 0.13 MgSO4/NaCl build prerequisites — reviewer sign-off + empirical onset tables (Machine B, 2026-08-06)

Two build-gating physics items (flagged by the cache-build recon as NOT covered
by the ratified acc725ea R1-R5) were adjudicated by the scientific-reviewer, and
the ocean-onset probe was run for both compositions (scan-res Pfreeze 2 MPa;
`/tmp/titan_tb_probe_results.json`, log `/tmp/titan_probe/nacl_mgso4_probe.log`).

**Item 1 — NaCl `extrap_ocean`: SIGN OFF = True (effective no-op ≡ clamp).**
`HydroEOS.py:193` truncates NaCl Pmax at 1000.1 MPa UNCONDITIONALLY, and the
ocean property interpolants are `RectBivariateSpline` (flat-extrapolation). So
above 1000.1 MPa, `extrap_ocean=True` and `False` return IDENTICAL values to
machine precision — a benign ~0.5% clamp bias on a thin/often-absent deep liquid
sliver. The MgSO4 24-sigma clamp artifact does NOT transfer (MgSO4 uses
`RegularGridInterpolator` LINEAR extrapolation, `MgSO4Props.py:151-156`).
REQUIRED: build-note rationale must state this is a no-op for NaCl, not a true
EOS extrapolation. Frozen no-ocean columns fill depth with SeaFreeze ice VI
(Pmax ~2300 MPa), in-domain regardless.

**Item 2 — per-composition freeze-line Tb grid: SIGN OFF WITH CONDITIONS.**
The ocean/no-ocean boundary is a steeply-TILTED diagonal, not a horizontal line;
a single NH3-style 12-node fine band is INADEQUATE. Empirical onset tables:

- **NaCl** — ocean bands are nearly DISJOINT across salinity (no single Tb has
  ocean at all w): w=1 → ocean Tb∈[252,272]; w=100 → [244,265]; w=290 →
  [233,241]. Union boundary sweep ≈ **[233, 272] K (~39 K)**. NaCl EOS T-floor
  is 229 K, so the w=290 onset (~233 K here) is near but above the floor —
  reachable. Requires a fine (0.5-1 K) Tb grid across the FULL union diagonal
  (~30-40 Tb nodes, NOT 12), with frozen corners retried as no-ocean
  (`retry_frozen_as_no_ocean=True`). Coarse 2-3 K spacing only above the dilute
  onset (>~252 K).
- **MgSO4** — genuine monotone onset ~[248, 255] K (nearly FLAT, unlike NaCl's
  diagonal): w=1 → ~252 K; w=100 → ~250 K; w=194 → ~250 K. ~15-20 nodes over a
  compact fine band ~[248,258] K. The w=194 (2 molal cap) column ALSO shows a
  low-Tb ocean "island" (ocean Tb=240, frozen 242-249, ocean 250+) —
  **ADJUDICATED PATHOLOGICAL** by the scientific-reviewer (2026-08-06, agent
  a3c9ed8ae24664527). Verdict: it is a Margules-lookup melting-monotonicity
  violation, NOT physical re-entrant melting. Decisive test — at fixed P=215 MPa,
  w=194: 238 K→ice II, 239 K→LIQUID, 240.5 K→ice III (liquid REFREEZES on
  heating, thermodynamically impossible), 249.5 K→real liquidus. The spurious
  liquid occupies a razor-thin ~21 MPa pressure lens and appears ONLY at exactly
  w=194 (swept 120/150/170/180/190/194 → violation at 194 only), i.e. the
  Margules free-energy crossover evaluated right at its 2-molal validity edge
  (`MgSO4Props.py::MgSO4PhaseMargules`). **Exclude via forced no-ocean.** CRITICAL
  mechanism caveat: `retry_frozen_as_no_ocean` triggers ONLY on
  `NoIceLiquidTransitionError`; the island nodes build SUCCESSFULLY as (phantom)
  oceans and do NOT raise, so the retry never fires on them. Two guards required:
  (a) by construction, place NO Tb nodes below ~248 K at w≳180 ppt (don't
  straddle the island with the fine band; do NOT extend the fine grid to 239 K
  "to resolve" a pathology); (b) post-build invariant — assert no cached node
  with w≳180 ppt AND Tb≲248 K carries `has_ocean=True`. Also add finer w
  resolution in the 150-194 ppt corner (boundary curvature drives ocean fraction).

**Launch conditions (both compositions, reviewer):** (1) onset tables inspected
[DONE — above]; (2) gate-3 half-cell Tb-shift → posterior-ocean-fraction test
run per composition on the chosen grid (the acceptance criterion that catches
the tilted-boundary failure); (3) cache build MUST use the production fine
`PfreezeRes_MPa` default, NOT the probe's 2 MPa; (4) NaCl extrap_ocean rationale
note corrected. There is NO committed joint-build driver — NH3 was built by a
now-deleted `/tmp` orchestrator (reconstructable from
`/tmp/nh3_joint_production_build.json`, `/tmp/nh3_patch_node.py`,
`/tmp/nh3_joint_production_orchestrator.log`). The MgSO4/NaCl drivers must clone
`cache_builder.build_tbw_grid_cache(...)` directly with
`ocean_overrides={'comp':...}`, `bulk_overrides={'Cuncertainty':0.06}`,
`retry_frozen_as_no_ocean=True`, `extrap_ocean=True`, the non-uniform Tb grid
above, and the acc725ea R3/R4 w-grid.

## 0.11 v5/v6 B1/B2/B6 EXECUTED (Machine B, 2026-08-06) — awaiting Machine A finalize

Priority-(a) of §0.10. Reviewer-adjudicated (three design items pre-ruled +
the FAIL adjudicated). NO gate tuning; the FAIL is surfaced, not relabeled.

**Preregistered scientific rulings implemented (reviewer, 2026-08-06):**
- **Crosscheck empirical-floor (§0.7 step-3):** injected the 0.36 km D_iceIh
  reference-wander floor into `d_pred` (the shape-excess mean-shift budget),
  NOT max-combined on `d_tol` (reviewer binding correction: B3 wander is a
  RIGID translation, physically the same object `_gaussian_ks_pred` subtracts).
  Opt-in `--empirical-floor` CLI param on `validate_sbi crosscheck`; default
  None → byte-identical for v6/all others. Relax-only (asserted). D_iceIh_km
  ONLY; NOT v6 (no measured wander).
- **B6 (limits anchor):** DEFERRED for v5+v6. Empty anchor set — the entire
  reachable Europa Im_k2 window sits below the 0.15 MCMC-falsified boundary;
  the single free anchor at Im_k2=0.004 already passes W1 on the pooled
  reference. NO new anchor-MCMC compute authorized.
- **Tb_K:** derived (PCHIP inversion of sampled D_iceIh_km), not sampled → NO
  SBC row; gate summaries carry an explicit `derived_params_sbc_na.Tb_K` N/A
  statement (reviewer Item-2).

**Artifacts / provenance:**
- Fresh pooled v5 reference: `/tmp/b3_build/europa_clipper_v5_reference_pooled_neff2000.pkl`
  (ESS=12713, D_iceIh wmean 61.950 km, between-seed std 0.055 km; B3 seeds
  101/202/303 concatenated, weights renormalized). Build script:
  `plans/scripts/pool_v5_reference_neff2000.py`.
- Runner edits: `plans/scripts/v5_run_gates.py` (pooled ref, SBC n-sbc→1500,
  crosscheck `--empirical-floor {"D_iceIh_km":0.36}`, validate_sbi_sha),
  `plans/scripts/v6_run_gates.py` (SBC n→1500, SHA; keeps its n_eff=500 ref,
  NO floor).
- Reports: `validation_reports/v5_gate_summary.json`,
  `validation_reports/v6_gate_summary.json` (both carry the reviewer
  adjudication block), per-arm dirs under
  `validation_reports/europa_clipper_v{5,6}_{baseline,noinduction,nok2}_1m/`.

**Results (baseline arm):**

| gate | v5 | v6 |
|---|---|---|
| SBC | **FAIL** — dC22_nh only (raw p=0.0017, BH-adj 0.0185, c2st 0.555); 10/11 PASS; 621 kept | PASS — all 11; dC22_nh BH-adj 0.8504; 560 kept |
| crosscheck | PASS (floor-augmented; D_iceIh d_excess 0.0171≤0.0281 even w/o floor) | PASS (stock HEAD) |
| limits | PASS | PASS |

**Reviewer verdict — v5 PASS WITH CONCERNS, v6 clean.** dC22_nh is a real but
LOCALIZED v5-flow-specific undertraining on the most sharply constrained
nuisance (prior/σ ~100:1 vs dC20_nh ~23:1 which PASSES p=0.205); v6's identical
nuisance is clean (p=0.61) → NOT a free-gravity-family systemic defect. The
B2 n→1500 fix EXPOSED it (old n=108 masked it; detection sharpens with n →
a true mild effect, not MC noise). rank-CDF (`.../v5_baseline_1m/sbc/
sbc_rank_cdf.png`) shows a mild MONOTONE dC22_nh deviation below the uniformity
band (Machine B inspected 2026-08-06) — consistent with real, mild.
- **Deployed v5 deliverable (D_iceIh/ocean/salinity) is doubly validated**
  (SBC-clean D_iceIh p=0.996, log10_w p=0.279 + reference-MCMC crosscheck PASS)
  → v5 ships WITH scope-note; dC22_nh marginal + any non-hydrostaticity readout
  DEFERRED to v6 (which owns that deliverable, clean). Ablation contrast
  corr(D,logw)=0.923/0.539/0.914 UNAFFECTED.
- Optional (Machine A's call): a cheap 2nd-seed SBC to confirm persistence; no
  retrain warranted for the deployed deliverable.

**FOR MACHINE A:** finalize the v5 ratification (ship-with-scope-note vs require
2nd-seed confirm) + countersign; deploy v6. Machine B stopped here — ratification
finalize + GUI wiring are Machine A duties.

## 0.12 #4 architecture pilot EXECUTED (Machine B, 2026-08-06) — capacity ELIMINATED

§0.10 item 2. Reviewer-adjudicated (design freeze + result interpretation).
Screen (seed 72, 4 arms D0/A/B/C on the full 690k ocean-admitting NH3 dataset,
uncapped, deployed early-stopping, z-score fixed). PILOT — production-arch
decision remains manager+reviewer+user.

**Result = clean uniform NEGATIVE.** D0 (deployed-arch control) reproduced 0.0423
(no-regression on the additive `train()` capacity/embedding plumbing). No arm
passes: pp-median 0.042-0.047 (target 0.0862), concentration_ratio >1 for EVERY
arm (1.215-1.322), and MORE capacity made concentration WORSE (C=1.322,
early-stops fastest). Capacity (A), embedding (B), both (C) all ineffective.

**Reviewer mechanism reframing (the key output):** the 0.042→0.135 gap
PARTITIONS —
- 0.042 (SBI) → ~0.093 (MCMC median): a REAL flow defect, but a
  TRAINING-SIGNAL / IDENTIFIABILITY problem, NOT representational (the flow fits
  the sharp C20/C22 channels, gets little gradient to sharpen the weakly-
  identified Im_k2; larger flows hit the easy optimum faster and never refine
  it → capacity is the wrong lever).
- ~0.093 (MCMC) → 0.135 (obs): NOT a flow defect. obs sits in the UPPER TAIL of
  even the EXACT posterior (frac≥obs only 0.19 for MCMC, NOT ~0.5). The
  model+prior cannot reach 0.135 while matching the other 3 observables (Tb-w
  degeneracy). CORRECTION: the "correct posterior drives frac≥obs → ~0.5"
  language in the pushforward shape reports is MISLEADING — reachable ~0.19.

**Rulings:** (1) capacity/embedding ELIMINATED as tested; a bigger flow is NOT
the next step. (2) MgSO4/NaCl PROCEED on the DEPLOYED architecture under
split-status — no arch change helps; but re-verify the tidal quarantine
PER-COMPOSITION (do NOT port the NH3 concentration-failure verdict by
assumption; different dissipation physics). (3) next NH3 tidal diagnostic =
upstream identifiability (mutual-info / conditional variance of Im_k2 given
{C20,C22,Re_k2} under Tb-w; Im_k2 heavy-tail vs z-scoring), NOT architecture.
(4) best_val/train_val_gap were NA (driver read sbi log-prob keys; sbi 0.26.1
uses `*_loss`) — FIXED in the driver; a corrected 4-arm × 3-seed re-run is
running to recover those diagnostics + give per-arm seed-variance insurance
(closes the "large arms premature-early-stopped" alternative before the arch
decision is FINALIZED). Reports: validation_reports/nh3_diagnosis/f4_architecture/.

**Consequence for §0.10:** MgSO4/NaCl datasets un-HELD (were already);
FLOW TRAINING may proceed on the deployed architecture under split-status once
the per-composition quarantine is re-verified. Machine B is building the
MgSO4/NaCl 2D joint caches + datasets now (architecture-independent, §0.10
priority-b).


Updated: 2026-08-10e (Machine A: MgSO4/NaCl COUNTERSIGNED split-status +
GUI-wired (task #52 Phase C complete); INDEX rows added; caveat-copy
corrections applied in the sector warnings; deploy snapshot rebuilt (8
caches, 317M) awaiting user HF ship. **Machine B next, in order:** (1)
Enceladus production (§2) — config freeze with Machine A FIRST; (2)
non-blocking archival follow-ups from the reviewers (commit any remaining
/tmp evidence); (3) NH3 upstream-identifiability diagnostic (§0.12 ruling 3)
when compute is free — manager+reviewer design freeze required before
running.)
Prior: 2026-08-10 (Machine B: **§0.16 CAMPAIGN COMPLETE — see the top block.**
MgSO4/NaCl full gate seq executed + adjudicated in two binding reviewer passes;
Re_k2 pushforward (req-val #1) + limits-monotone doc (req-val #2) discharged;
BOTH comps CLEARED for split-status deploy — no STOP. Remaining = Machine A GUI
wiring (#52) + 2 caveat-copy corrections. commits 99ed8c42 + this.)
Prior: 2026-08-10 (Machine A: §0.16 — FLOW TRAINING conditional GO for
MgSO4/NaCl. Quarantine re-verification restructured as the mandatory FIRST
post-training gate (needs a trained flow; split-status default). Conditions:
commit the #4 re-run final report (stop+escalate if it overturns), SHA-verify
regenerated datasets, deployed architecture/recipe exactly. Start fresh
per-comp reference MCMC (B3 protocol) NOW in parallel. Gate sequence
preregistered in §0.16.)
Prior: 2026-08-07 (Machine B: §0.15 — Tb=252 K defect root-caused to
ResetNearestExtrap int-truncation [reviewer-corrected from fn_phase]; BOTH
caches REBUILT float-coerced + gate-3 BINDING PASS + installed; BOTH 1M gens
DONE on repaired caches [NaCl 631,214 / MgSO4 691,075 kept, cache-sha verified];
GATE-3 RE-ADJUDICATED on rebuilt bytes by scientific-reviewer [agent
a02a9038b2de1d382] — BOTH PASS, no blocking issue, NaCl-at-bar & MgSO4-zero-None
both confirmed correct. Flow TRAINING held on #4 pilot. Machine A: re-verify new
cache bytes + review latent ResetNearestExtrap dtype bug. Prior Machine A entry:) §0.11 CLOSED — v5
RATIFIED scoped + v6
RATIFIED, GUI ungated, deploy snapshot rebuilt; §0.14 added — MgSO4/NaCl
caches ACCEPTED, per-comp gate-3 + 1M dataset gens GO, flow training still
HELD on per-composition quarantine re-verification; two non-blocking record
items for B in §0.14. Machine B's next work, in order: (1) per-comp gate-3,
(2) 1M dataset generations, (3) per-composition tidal-quarantine
re-verification, (4) commit the /tmp acceptance evidence + finish the #4
corrected 4-arm×3-seed re-run, (5) Enceladus config freeze with Machine A.)
Prior: 2026-08-06 (Machine B: BOTH MgSO4/NaCl PRODUCTION CACHES BUILT,
VALIDATED, and PUSHED for Machine A review — NaCl at
Test54_nacl_ocean/titan_nacl_joint_structure_grid_2d.pkl (600 nodes, 314 ocean/
219 no-ocean/67 None, tilted-diagonal has_ocean map), MgSO4 at
Test54_mgso4_ocean/titan_mgso4_joint_structure_grid_2d.pkl (272 nodes, 191 ocean/
81 no-ocean/0 None, ALL reviewer guards CLEAN in-manifest: island invariant 0
violations, extrap-ceiling 0 rejects [deepest w=194 node 1371<1400 MPa], dρ/dP∈
[0.14,0.17] all extrap columns, ocean-base T≤300 K, 45 deep + 7 mild
eos_extrapolated flags on the high-w tail). Ceiling was raised 1200→1400 MPa
(reviewer, physics clean to 1500). NEXT (Machine B): per-comp gate-3 (half-cell
Tb-shift, MgSO4 w≥150 corner) then 1M dataset gens; flow TRAINING still HELD on
per-composition quarantine re-verification. Reviewer (a3c9ed8ae24664527)
adjudicated: MgSO4 w=194 low-Tb
island = PATHOLOGICAL Margules melting-monotonicity violation → EXCLUDE by
construction + post-build has_ocean invariant [see project_mgso4_margules_island
memory]; NaCl w=290 monotone [/tmp/nacl_monotonicity_check.json] + retry-corner
discriminator PASS [/tmp/nacl_corner_discriminator.json]. MgSO4 extrap_ocean=True
adopted (linear extrap physically sound: monotone ρ, stiffening K_T verified to
1400 MPa; clamp is unphysical) with 4 conditions wired into the NEW committed
build driver plans/scripts/titanG_build_ocean_cache.py [eos_extrapolated flag +
hard extrap ceiling reject + pre-bake corner check + post-build dρ/dP>0]. OPEN:
production PPTitan reaches P_basal ~1371 MPa at w=194 (deeper than the PPTest50
probe's ~1109) so the 1200 MPa ceiling rejects the whole w=194 column — asked
reviewer to raise to ~1400 (physics verified clean there). NaCl config
test54_titan_nacl_freegrav.json authored + NaCl full production cache building to
/tmp (40 Tb × 15 w = 600 nodes). Validation subsets confirmed the tilted-diagonal
geometry for both comps. §0.13 refreshed. §0.10
EXECUTED — v5/v6 B1/B2/B6 gates run vs the
FRESH pooled n_eff~2000 v5 reference, reviewer-adjudicated. v5 = PASS WITH
CONCERNS (ships for its DEPLOYED D_iceIh/ocean/salinity deliverable with a
mandatory scope-note; Machine A finalizes + countersigns): SBC FAILs on dC22_nh
ONLY (BH-adj p=0.0185, 621 kept pairs — the B2 n→1500 fix EXPOSED a
miscalibration the old n=108 masked), a real but LOCALIZED v5-flow-specific
undertraining on the ~100:1-compressed nuisance (rank-CDF mild monotone below
band; NOT free-gravity-family — v6's identical nuisance PASSES p=0.61); the
deployed marginals are doubly validated (SBC-clean D_iceIh p=0.996 + crosscheck
PASS even without the 0.36 km floor). v6 = ALL GATES PASS. B6 DEFERRED (empty
anchor set; window below the 0.15 falsified boundary — reviewer Item-3). #4
architecture-pilot EXECUTED — capacity/embedding ELIMINATED (§0.12): the Im_k2
gap partitions into a real flow-defect (0.042→~0.093,
training-signal/identifiability, NOT representational) + a data-limit
(~0.093→0.135, obs in the upper tail of the EXACT posterior). Reviewer: proceed
MgSO4/NaCl on the DEPLOYED architecture under split-status (re-verify quarantine
per-composition). Machine B building MgSO4/NaCl 2D caches + datasets now. See
§0.11 (v5/v6 gates), §0.12 (#4). Prior NH3 #1 (PERSIST → salinity ELIMINATED)
in §0.9.
Prior: 2026-08-01 at genai `54106fbd` (Machine A refresh after v5/v6/v7
delivery). Authoritative executable queue for compute-intensive work. Machine B
should pull the exact `origin/genai` commit named by Machine A before each
campaign, use the `PPcl` mamba environment, record package versions and seeds,
and return artifacts plus machine-readable reports. Never tune thresholds
after seeing a failure: stop and surface the evidence to Machine A.

Machine B model roster: Claude Opus 4.8 / Sonnet 4.6 / Haiku 4.5.

## 0. v5/v6/v7 gate adjudication COMPLETE (2026-08-02) — scoped follow-ups

Full record: `plans/active/europa-v5v6v7-gate-adjudication.md` (manager +
independent Opus-5 scientific review). Verdicts: v5 NOT ratifiable
(D_iceIh shape-excess survives HEAD, 1.12x; zeta_Ih mean 1.074x; SBC
underpowered n=108); v6 conditionally ratifiable (clean at HEAD; SBC n=102
below the >=200 spec minimum); v7 blocked — its crosscheck FAIL is
adjudicated a REFERENCE-side artifact (true v7 posterior provably identical
to v5's at the fiducial; mass in the newly opened support region measured
0.000000 in both references; the 1.06 km reference disagreement has the
wrong sign and rigid-translation shape of nested-sampling log-volume error;
which reference wandered is OPEN — v5's log_Z_err 0.281 makes it the more
suspect run). NO retraining and NO dataset regeneration — execute exactly
these follow-ups, in order, and stop at any surprise:

- **B1** Regenerate ALL gate reports (v5/v6/v7 baselines + arms) at current
  HEAD: use the DEFAULT constructed sweep grid (do not pass
  --sweep-values), record the validate_sbi.py commit SHA in every report,
  apply BH-FDR uniformly, reconcile `v{5,6,7}_gate_summary.json`.
- **B2** SBC re-runs at n_sbc_pairs >= 500 for v5/v6 baselines; report raw +
  BH p for every param, explicitly D_iceIh_km / log10_wOcean_ppt / Tb_K.
- **B3** Reference wander: re-run BOTH v5 and v7 reference MCMC, >=3 fresh
  seeds each, SAME pinned environment (note scipy 1.16.3→1.17.1 straddled
  the originals), n_live >= 2000. Report per-seed D_iceIh/log10_w means and
  log_Z ± err; success = between-seed scatter brackets the observed
  1.06 km. Do NOT assume v7 is the outlier.
- **B4** Evidence-ratio check: ΔlogZ vs −ln(V7/V5), V7/V5 by Monte Carlo
  over the prior under the two induction-bound sets (expected ~2.5).
- **B5** v7 flow-side diagnostics (preregistered; required before any v7
  ratification): (i) flow-vs-flow v5/v7 at the fiducial; (ii) v7-flow vs
  v5-reference formal crosscheck incl. shape clause; (iii) synthetic-truth
  bias test; (iv) local expected-coverage (TARP) near the fiducial.
- **B6** Limits anchor mode (W1 <= 0.25 sigma vs ground-truth MCMC) at
  reachable Im_k2 values, all three campaigns.
- **B7** Record outcomes as FAIL-ADJUDICATED-ACCEPTABLE where applicable;
  never relabel a FAIL to PASS.

Priority vs Titan: Titan NH3 Task #68 (section 1) remains your active task;
run B1–B6 after the 3x4 validation completes, before the full NH3
production build.

## 0.5 RESOLVED — NH3 cache delivered; split ratification COUNTERSIGNED (2026-08-04)

Cache landed (`3e9a275a`) — the GUI slot now runs end-to-end on Machine A.
The PPC finding and SPLIT ratification (`a5f83050`) are **countersigned by
the manager after independent reproduction on A**: a full GUI generate gives
SBI pushforward medians Re_k2 = 0.542 / |Im_k2| = 0.042, matching your PPC
(0.541/0.042) against obs 0.608/0.135 and MCMC 0.581/0.093. Machine A added
a results-panel sector warning threaded from the slot (the k2 panes and
heating tab display the unverified sector; scope-note-only disclosure was
insufficient). Standing requirements you set for MgSO4/NaCl (PPC + a
first-class pushforward-observable crosscheck gate that must flag the NH3
miss + flow-under-update diagnosis before more compute) are ADOPTED as
manager policy. Reminders: add ocean composition to future cache metadata;
GUI wiring + ratification countersign are Machine A duties — deliver
artifact + gates + recommendation.

**B5 EXTENDED (2026-08-04):** the NH3 PPC proved SBC-PASS + per-param
crosscheck can coexist with datum-local pushforward under-update. Run the
same PPC (posterior-predictive at the fiducial, SBI vs reference MCMC vs
data, per observable channel) for ALL Europa artifacts: **deployed v4 and
v1.1 first** (they are on the public app), then v5/v6/v7. If v4's k2 or
gravity pushforward sits at the prior-predictive median, surface
immediately — do not wait for the rest of the batch.

**RESOLVED 2026-08-04 (B) — v4 + v1.1 PPC verdict: no tidal-sector warning
needed; HF redeploy UNBLOCKED on this criterion.** Both deployed artifacts:
0 channels flagged (v4 0/21, v1.1 0/5) at the deployed defaults. Report:
`validation_reports/EUROPA_PPC_BATCH_v4_v1p1.md`; JSONs under
`europa_clipper_v4_1m/ppc/` + `europa_galileo_v1p1_1m/ppc/`. Scientific-reviewer
**PASS WITH CONCERNS** (2026-08-04), all required corrections folded in. The
clean verdict does NOT rest on "wide σ = clean" (that framing was retracted — it conflates
flag insensitivity with flow verification):
- **v4 Re_k2 is the decisive positive result:** datum at the 0.3 prior-predictive
  percentile (a *more* extreme tail than NH3's 86th), MCMC-pp updates 1.43σ off
  the prior, SBI-pp tracks to 0.04σ. On the one *informative* deployed tidal
  channel where the NH3 under-update could appear, it does not — on a harder datum
  than NH3's. **This contradicts the "tail → under-update" causal story**; keep it
  in mind for the diagnosis (hypothesis #2 tail-sparsity is not sufficient alone).
- **v1.1's tidal channels are non-informative by construction** (hypothetical
  k2/h2, σ 0.05/0.10; MCMC updates ≤0.13σ) — the NH3 pathology cannot arise at a
  non-informative default, so 0-flagged there is expected, not verification.
- **License persisted:** the identical generalized statistic on NH3 flags both k2
  (Re 0.94σ, Im 1.64σ), neither gravity — archived at
  `titan_freegrav_nh3_1m/ppc/ppc_pushforward_report.json`.
- **Model-adequacy caveat (not a warning trigger):** deployed default Im_h2 = 0 is
  physically unreachable (dissipation > 0) and Re_k2 = 0.23 sits below the
  prior-predictive tail; both samplers agree at the model edge (≤0.013σ). Recorded,
  no action.
- **Untested regime:** the PPC covers deployed defaults only. Flow behavior at an
  *off-default informative* k2 (tight σ, extreme value — the NH3 regime) is
  untested; a follow-up PPC there is recommended before the GUI relies on the flow
  for user-supplied extreme tidal hypotheticals. Does not gate this redeploy.

**COUNTERSIGNED by the manager 2026-08-04** after independent A-side
reproduction: full GUI generate on the v4 slot gives SBI pushforward medians
Re_k2 = 0.2531 / |Im_k2| = 0.0031 vs your report's 0.2534/0.0030 (MC noise).
The v4 Re_k2 result also weakens diagnosis hypothesis #2 (tail sparsity) as
a standalone mechanism — when you run the under-update diagnosis, lead with
hypothesis #1 (noise-augmentation swamping) and #3 (x-norm scale), and look
for what DIFFERS between the NH3 and v4 training setups (e.g. 4-channel vs
21-channel obs vector: with 21 channels the k2 pair may retain relative
weight that a 4-channel NH3 vector loses; also the joint no-ocean mixture
may bimodalize the conditional). HF redeploy decision returned to the USER.

**RESOLVED 2026-08-04 (B) — v5/v6/v7 baseline PPC batch: all three clean.**
0 channels flagged (v5 0/21, v6 0/20, v7 0/21) at the deployed defaults. Report:
`validation_reports/EUROPA_PPC_BATCH_v5_v6_v7.md`; JSONs under
`europa_clipper_{v5_baseline,v6_baseline,v7_openae}_1m/ppc/`. Scientific-reviewer
**PASS WITH CONCERNS** (2026-08-04), all corrections folded. Key points:
- **Non-trivially clean, unlike v4/v1.1.** These 11D baselines carry MANY
  strongly-informative channels and the flow tracks the MCMC on all of them: v5
  C20 22.6σ update → 0.12σ; v6 Re_k2 1.66σ + synodic_x 2.8–3.1σ → ≤0.07σ; v7
  Bind_synodic_x_real **52.9σ** update → 0.18σ (the most stringent single test in
  the whole PPC program). Batch-max gap 0.27σ, well under the 0.5σ flag.
- **BASELINES ONLY.** The noinduction/nok2 ablation ARMS condition on different
  observable subsets and have NO matching reference MCMC, so the |SBI−MCMC| flag
  is undefined for them — they need their own reference MCMC before a PPC. These
  arms are the controlled probe for the obs-vector-width hypothesis
  (§0.6 diagnosis) — a 6-channel v6-noinduction tests width at fixed body/physics/
  sampler. Cheapest next diagnostic ahead of any retraining.
- **Obs-vector-width read (this batch is supporting, NOT controlled).** All 5
  clean Europa PPCs are 20–21ch; the only under-update (NH3) is 4ch. Consistent
  with hypotheses #1/#3, but confounded (body/comp/sampler/joint-mixture) — the
  arms PPC is the controlled test.
- **v7 caveat carried to B5.** v7's observable-space gaps are batch-largest
  (0.10–0.27σ, ~2–3× v5/v6) — likely MC noise from the diffuse open-|Ae|
  posterior, unconfirmed; carry as a flow-side input to B5. Also: the v5/v7
  references are the pre-B3 (possibly-wandered) runs, but this PPC is ORTHOGONAL
  to the wander (parameter-space D_iceIh shift, ≈0 observable-space mass diff) —
  when B3 re-runs the references, re-confirm these PPC medians are stable as cheap
  corroboration.
- **Flow-fidelity diagnostic ONLY** — does NOT affect the ratification blocks (v5
  D_iceIh shape-excess, v6 powered-SBC, v7 reference-wander B3–B5).

## 0.9 MANAGER JUDGEMENT (2026-08-05) — B3 accepted; reorder the follow-ups: #3, #2, then #1 only if needed

**B3 accepted and consequential.** The legacy 1.06 km v5-v7 reference
disagreement is adjudicated a resolution artifact of the n_eff=500 SMC runs;
at matched n_eff=2000 the gap is −0.19 ± 0.22 km, inside the between-seed
floor. Neither old reference was "the outlier" — BOTH wandered (fresh
converged means: v5 61.35–61.65, v7 61.61–61.80 vs old 62.1/61.0). Two
record corrections follow:
- My 2026-08-02 statement "v7 flow passes 22/22 gates against the v5
  reference" is VOID — it was measured against the old wandered v5
  reference. Against the fresh references, v5-flow D_iceIh sits ~+0.7 km
  high (inside 0.25σ tol) and v7-flow ~+1.2 km high (at/near tol). Both
  offsets share a sign — the D–w degeneracy-direction signal persists and
  B2/B5 remain the tests that matter. The adjudication record gets an
  addendum this session.
- **§0.7 steps 2–3 must target the FRESH references:** pool the three
  n_eff=2000 seeds per config as the crosscheck target, take the mean_tol
  MC-error term from the between-seed scatter, and use the preregistered
  step-3 shape-excess floor of 0.36 km (the reviewer's 2× conservative
  factor on the 0.18 km empirical floor). Do not crosscheck anything
  against the seed-51/71 n_eff=500 references again.

**Separators accepted:** S1 rules the joint mixture OUT as the Im_k2 driver
(ocean-only moved it the wrong way); S2 closes noise-swamping for Im_k2.
The elimination now rests entirely on the size of the real gap — which is
measured against an n_eff=500 NH3 reference (MCMC-pp ceiling 0.100) of
exactly the resolution class B3 just discredited. Therefore:

**Follow-up order (reordered from the reviewer's #1-first):**
1. **#3 first (free):** plot the anchor's Im_k2 pushforward distribution —
   concentration-failure vs wrong-mode shape, no compute.
2. **#2 second (decisive, cheap):** re-run the NH3 reference MCMC at
   n_eff=2000 (>=2 seeds if affordable) and re-measure the MCMC-pp ceiling.
   The entire "flow under-updates by 0.042-vs-0.100" magnitude inherits the
   discredited resolution class; measure the target before chasing it.
3. **#1 (salinity-fixed retrain) ONLY IF a material gap survives #2** — its
   interpretation depends on the matched-resolution ceiling.

**#3 EXECUTED (Machine B, 2026-08-05) — reviewer PASS.** Pushforward-shape
diagnostic run on BOTH the capped anchor AND the deployed 1M flow (reviewer-
required transfer check). Result: **concentration-failure CONFIRMED, wrong-mode
EXCLUDED.** Deployed flow concentration ratio 1.215 (anchor 1.289) — both >1, so
the posterior-predictive |Im_k2| is not narrower than the prior (it is marginally
broader); the high-k2 tail is fully retained (frac≥0.100 = 0.233 deployed / 0.250
anchor), not collapsed. The flow updates in the correct direction but grossly too
little: frac≥obs moved only 0.128→0.161 (deployed) / →0.173 (anchor) where a
correctly-conditioned posterior should approach ~0.5. **This is measured against
the PRIOR, so it is ceiling-INDEPENDENT and survives whatever #2 finds.**
Nonfinite-drop bias closed: <0.25% of draws drop, and their dissipation-proxy
direction is conservative (dropped set leans high-dissipation, so including them
would raise SBI-pp toward obs — weakening, not manufacturing, the under-update).
**#3 does NOT obviate #2** — it only fixes the *mechanism* (fails to concentrate,
not wrong mode); the *magnitude* still needs the matched-resolution ceiling.
Proceeding to #2. Report: `validation_reports/FLOW_UNDERUPDATE_DIAGNOSIS.md`
(Follow-up #3 EXECUTED); artifacts `validation_reports/nh3_diagnosis/
pushforward_shape{,_deployed}/`.

**MgSO4/NaCl release criteria (preregistered):** if the matched-resolution
NH3 SBI-vs-MCMC-pp gap falls below the 0.5σ_obs flag threshold, the
under-update is substantially a reference artifact — MgSO4/NaCl PROCEED
with the standard gate set + the pushforward gate, and the NH3 split
ratification comes back to the manager for re-adjudication (the tidal
sector warning may be softened). If the gap survives, run #1; remedy
selection (with reviewer + user sign-off if artifact design changes)
before any MgSO4/NaCl compute. Either way, report the four-way table.

**#2 EXECUTED (Machine B, 2026-08-05) — reviewer PASS. Gate outcome: the gap
SURVIVES → run #1; MgSO4/NaCl do NOT proceed.** Matched-resolution NH3 reference
MCMC at n_effective=2000 / n_active=1024 (regime ratio 1.953, identical to the
legacy 500/256; tracked config UNMUTATED, config_hash 1611b65fff3f06c9 intact),
seeds 72 + 172, both annealed to β=1 (7187 samples each). **Weighted |Im_k2|
median pooled 0.1037** (seed 72: 0.1043, seed 172: 0.1031; between-seed std
0.00086, range 0.0012). Contrast with B3: the ceiling did NOT collapse — it moved
only +0.0038 (0.0999→0.1037, +0.11σ_obs) and *away* from the SBI value, so the
n_eff=500 ceiling was NOT a low-resolution artifact.

Four-way table (|Im_k2|):
| quantity | value | vs obs |
|---|---|---|
| prior-predictive median | ~0.05 (broad) | — |
| **deployed SBI posterior-predictive median** | **0.0423** | −3.2σ_obs |
| **matched-res MCMC posterior-predictive median** | **0.1037** | −0.89σ_obs |
| observed | 0.135 (σ 0.035) | — |

**SBI-pp vs matched MCMC-pp gap = 0.0614 = +1.76σ_obs ≫ 0.5σ_obs threshold →
survives decisively** (robust to seed: worst per-seed still 1.74σ_obs). Reviewer
clarification (recorded): two distinct offsets exist — (a) matched MCMC-pp sits
0.89σ below obs = ordinary model/data tension (the physics grid cannot fully
reach 0.135; expected, NOT a remediation target); (b) SBI-pp sits 1.76σ_obs below
MCMC-pp = the flow deficiency, the actual target of #1. With #3 (concentration-
failure confirmed, wrong-mode excluded), the reading is a genuine flow-training
under-update of the tidal sector, not a comparison artifact and not mode
collapse. **DECISION: run #1 (salinity-fixed ocean-only retrain) with reviewer +
user sign-off before any MgSO4/NaCl compute; MgSO4/NaCl PROCEED is falsified.**
Report: `validation_reports/nh3_diagnosis/matched_reference/matched_reference_report.json`;
per-seed pickles + progress jsonl copied alongside.

**#1 AUTHORIZED (manager, 2026-08-06).** #2's outcome satisfies the §0.9
preregistered branch condition, and #1 is a pilot retrain on the existing
dataset with no artifact-design change — no further sign-off needed to RUN
it (the reviewer+user sign-off requirement attaches to REMEDY selection
after #1, if the remedy changes artifact design). Design notes: salinity
FIXED at the reference posterior's median (state the value in the report);
ocean-only rows (reuse the validated S1 has_ocean recovery); same
architecture/seed family and 60-epoch cap as the S-pilots so the capped
anchor (0.043) stays the comparison; report the four-way table + flag
statistic against the matched MCMC-pp ceiling 0.1037. Preregistered
reading: Im_k2 gap collapses → the salinity axis is the driver and the
remedy discussion (with reviewer + user) starts from
degeneracy-aware options; gap persists → the remaining candidates are
capacity/embedding (#4) and the elimination restarts from the widened
residual. Either way, stop after #1 and surface — remedy selection is a
manager + reviewer + user decision. v5/v6 §0.7 work keeps priority for
compute scheduling; MgSO4/NaCl stay HELD.

**#1 EXECUTED (Machine B, 2026-08-05) — reviewer PASS. Outcome: PERSIST →
salinity ELIMINATED as the driver; remaining candidate is capacity/embedding
(#4).** The design-review reviewer (PASS WITH CONCERNS) required a
sample-size-matched, salinity-VARYING control before the reading could be acted
on: the fixed-salinity band keeps only ~9% of rows (60,039) vs the anchor's
~690k — an ~11× cut that on its own biases toward under-concentration (the same
Im_k2≈0.04 signature). Both corrections folded in (matched-N control added;
band re-documented as a symmetric ±0.084-dex salinity slice about 12.6 ppt, not
a Voronoi cell). Both pilots trained to the 60-epoch cap (61 ep, no early-stop);
config_hash e596574d1e81567c matches S1/anchor (no artifact-design drift).

| Pilot | N_train | salinity | Im_k2 pp-median | dev vs obs |
|---|---|---|---|---|
| banded (fixed ~12.6 ppt) | 60,039 | fixed | 0.0313 | 3.16σ |
| control (varying, matched N) | 60,039 | varying | 0.0321 | 3.15σ |
| S1 (varying) | 642,558 | varying | 0.0392 | 2.98σ |
| capped full-joint anchor | ~690k | joint | 0.043 | — |
| matched MCMC-pp ceiling (#2) | — | — | 0.1037 | — |
| obs | — | — | 0.135 | — |

Deltas (σ_obs=0.035): **banded−control = +0.02σ** (salinity axis IN vs OUT →
essentially zero, and wrong sign for the salinity hypothesis); control(60k)−S1(643k)
= −0.20σ (pure size effect: more data → toward obs, as predicted — discharges the
MAJOR confound); banded−ceiling = +2.07σ (still-open flow gap). Reviewer PASS
(2026-08-05): salinity eliminated as the driver of gap (a) [SBI-pp vs MCMC-pp
ceiling], remaining candidate is capacity/embedding (#4). Two-gap scoping to
surface: gap (a) +1.76–2.07σ = the flow deficiency (#4); gap (b) ceiling 0.1037
vs obs 0.135 = +0.95σ model/data tension, NOT #1's target. Size gain is strongly
sublinear (11× data bought +0.007; ceiling needs +0.065) → more ocean-only data
won't close gap (a). Banded pp 5–95 = [0.0023, 0.340] already spans past obs →
an expressiveness/concentration signature, not a support failure. **STOPPED per
protocol; remedy (#4) selection is manager + reviewer + user.** Optional #4
hardening (non-blocking): 2–3 seeds for a formal training-noise band; one
single-node retrain to confirm the ±21% band residual. Driver
`plans/scripts/nh3_diag_1_salinity_fixed.py`; reports
`validation_reports/nh3_diagnosis/f1_salinity_fixed/{f1_train_manifest.json,banded/,control_varying/}`.

## 0.8 MANAGER DECISION (2026-08-04) — run both separators; MgSO4/NaCl stay HELD

Adjudicating the decision your diagnosis surfaced (FLOW_UNDERUPDATE_DIAGNOSIS.md
"Consequences" 3): run BOTH cheap separators before any MgSO4/NaCl compute.
Rationale: each MgSO4/NaCl campaign costs a cache + ~1M sims + reference MCMC;
the separators cost two pilot flow trainings on data that already exists. If
the mixture is the driver, the remedy (e.g. a has_ocean-labelled conditional
or mixture-aware embedding — an artifact-design change needing reviewer + user
sign-off) must be chosen BEFORE those campaigns; proceeding under split-status
would ship every future ocean campaign with a permanent tidal-sector asterisk.

- **S1 — ocean-only-with-salinity pilot.** Filter the EXISTING NH3 1M dataset
  to ocean-branch rows (recover has_ocean per row from the cache node tags /
  the same lookup generation used — state the method in the report), retrain
  a pilot flow with the same 13-param space and architecture (fewer epochs
  acceptable at pilot scale), PPC at the Titan datum. Preregistered reading:
  clean k2 assimilation → the joint MIXTURE is the driver; a persistent miss
  → the salinity degeneracy (or another shared element) is implicated and
  the mixture reading weakens.
- **S2 — reduced-noise Im_k2 pilot.** Same dataset (joint, unfiltered), k2
  noise augmentation reduced (e.g. sigma/4, and zero if cheap to add as a
  second arm), retrain pilot, PPC. Preregistered reading: the Im_k2 gap
  shrinking materially → the abs-fold/noise convention contributes for
  Im_k2; unchanged → #1 closed for Im_k2 as well.
- No new forward-model sims are authorized for S1/S2. Report the standard
  four-way tables + the flag statistic for both pilots, either outcome.
- **Ordering:** §0.7 (v5/v6 fix path) stays the top priority. The pilots are
  GPU/CPU flow trainings — run them while the B3 reference chains cook, not
  before launching B3.
- **MgSO4/NaCl remain HELD** until S1/S2 report and the manager chooses
  proceed / remedy (+ reviewer and user sign-off if the remedy changes
  artifact design).

### LAUNCHED 2026-08-04 (Machine B, commit bb70cab3)

Both B3 (§0.7 step 1) and the two separators are RUNNING; scripts written and
reviewer-ratified (PASS-WITH-CONCERNS) before launch, corrections folded.

- **B3** — `plans/scripts/b3_reference_wander.py`, v5+v7, seeds 101/202/303.
  pocoMC 1.2.6 has NO `n_live`; "n_live>=2000" mapped to `n_effective=2000`
  with `n_active=1024` raised in step (reviewer: raising n_effective alone is a
  regime change; ratio-preserving n_active keeps train cadence). The driver
  computes the **matched-resolution paired v5-v7 gap** at n_eff=2000 (the
  step-3 comparison) — NOT the stale n_eff=500 1.06 km anchor. Env versions +
  n_active/n_total/n_evidence recorded. Per-seed pickles → /tmp; committed
  primaries untouched. Runner: `sampler_settings.n_active` now exposed
  (additive; existing runs default to pocoMC 256, unchanged). Log
  /tmp/b3_reference_wander.log; report → validation_reports/b3_reference_wander/.
- **Separators** run as a sequential chain (anchor→S1→S2) so torch trainings
  don't all hit CPU at once alongside B3:
  - `nh3_diag_capped_anchor.py` — capped full-JOINT reference so the pilots
    compare **cap-vs-cap**, not against the historical converged 0.042
    (reviewer required-validation). If this anchor also under-updates at the
    cap, pilot "gap unchanged" readings are trustworthy.
  - `nh3_diag_s1_ocean_only.py` — has_ocean recovered per row via
    `grid_interp_2d` nearest node (reviewer: 99.8% agreement w/ the forward
    model's selection, 100% at boundary cells; forward model picks the
    higher-weight corner, no cross-branch x-blend).
  - `nh3_diag_s2_reduced_noise.py` — x_clean recovered EXACTLY by reproducing
    the one-shot post-fold `default_rng(7272)` noise draw and subtracting
    (reviewer + empirical: reconstructed |Im_k2| min 0.0, 0% negative =
    abs-fold signature). Reduced arm (k2 σ/4) is a deterministic scaling of
    the original draw (seed reuse); zero-k2 arm added.
  All pilots cap at 60 epochs and record `epochs_trained` (confirm early-stop
  fired before trusting any "unchanged" reading). Log /tmp/nh3_pilot_chain.log;
  reports → validation_reports/nh3_diagnosis/{capped_full_joint_anchor,
  s1_ocean_only,s2_reduced_noise}/. NO new forward sims.
- **Next**: scientific-reviewer interprets B3 matched-gap + the three PPC
  four-way tables; manager chooses proceed/remedy for MgSO4/NaCl.

### REPORTED 2026-08-04 (Machine B) — B3 + separators COMPLETE, reviewer PASS WITH CONCERNS

Both compute streams finished (exit 0); scientific-reviewer adjudicated (this is
a scientific-adjudication call, not delegated). Full write-up:
`validation_reports/FLOW_UNDERUPDATE_DIAGNOSIS.md` "Separators S1/S2 + capped
anchor" block. Four flows nsf/seed72/max-epochs60 (all ran to cap,
epochs_trained=61 — no early stop), n_post=4000, PPC noiseless.

| flow | Re_k2 pp_med | Im_k2 pp_med |
|---|---|---|
| ANCHOR capped full joint | 0.5446 | 0.0431 |
| S1 ocean-only (mixture off) | 0.5408 | 0.0392 |
| S2 k2 σ/4 | 0.5381 | 0.0461 |
| S2 k2 zero-noise | 0.5400 | 0.0469 |

obs Re=0.608 / Im=0.135; MCMC-pp Im ceiling 0.100.

Reviewer verdicts: **cap CONFIRMED** not a confound (anchor 0.043 = deployed
converged 0.042); **mixture CONFIRMED not the driver** (S1 moved Im_k2 the WRONG
way — a bimodal low-k2 drag would have raised it); **noise-swamping #1 CONFIRMED
closed for Im_k2** (zero-noise + σ/4 both ≈ anchor). **B3 CONFIRMED artifact**:
matched n_eff=2000 D_iceIh gap −0.19±0.22 km collapses into the between-seed
floor 0.18 km (|mean|/std 0.88); neither v5 nor v7 an outlier; step-3
shape-excess floor ~0.18 km at n_eff=2000 (reviewer advises 2×=0.36 km gate).

Residual mechanism (by elimination, PLAUSIBLE not established): flow fails to
concentrate on the weakly-identified high-k2 ocean tail; salinity axis
(Tb↔w −0.986) most probable specific contributor (the one apparatus element S1
left in). NOT ruled out: salinity axis, the abs-fold *as a representation*
(S2 changed only noise, x_clean still folded |Im|), the k2 support guard, and
sub-ceiling identifiability (obs 0.135 > MCMC 0.100 → part is model-datum
tension). Reviewer flag: PPCs compared SBI-pp to obs 0.135, not to MCMC-pp 0.100.

**Reviewer-required follow-ups BEFORE MgSO4/NaCl proceed:**
1. Salinity-fixed (or sharply-narrowed-w) ocean-only retrain at same cap/seed/arch
   — decisive salinity-vs-fold/support cut (the axis S1 could not touch).
2. Re-measure matched MCMC-pp for the anchor (ideally ocean-only) so the target
   is SBI-pp vs MCMC-pp, not vs obs.
3. Plot the anchor SBI Im_k2 pushforward distribution (concentration vs
   wrong-mode), not just the median.

**MANAGER (Machine A) CALL PENDING:** proceed to MgSO4/NaCl vs run follow-up #1
first. MgSO4/NaCl stay HELD per §0.8 until the manager adjudicates.

## 0.7 USER DIRECTIVE (2026-08-04) — v5 and v6 must be FIXED and DEPLOYED to HF

This is now the TOP Europa priority, superseding the earlier "after Titan"
ordering. Goal: bring v5 and v6 to ratifiable state, manager re-adjudicates,
Machine A unwires the not_ratified gating, snapshot ships. Your v5/v6/v7 PPC
batch (all clean, countersign pending nothing — accepted as flow-fidelity
evidence) already closes one requirement. Execute in this order:

1. **B3 first (it feeds everything):** multi-seed v5 + v7 reference re-runs
   (>=3 seeds, n_live >= 2000, pinned env). Besides settling the v7
   reference question, the v5 between-seed scatter gives the EMPIRICAL
   reference-noise floor that the current self-D bootstrap understates
   (it resamples correlated nested-sampling draws and misses log-volume
   error — adjudication doc, "other items" #1).
2. **B1 + B2 for v5/v6:** HEAD gate regeneration (constructed sweep grid,
   gate-code SHA, uniform BH-FDR) + powered SBC (n >= 500; report
   D_iceIh_km / log10_wOcean_ppt / Tb_K explicitly).
3. **v5 shape-excess re-evaluation (preregistered, NOT tuning):** recompute
   the crosscheck d_floor/d_tol for v5 D_iceIh_km with the B3 empirical
   reference scatter included in the null (state the formula in the report
   BEFORE computing the new statistic; accept whatever it gives). If the
   excess survives the corrected floor AND powered SBC flags D/w again →
   retrain path, come back for a frozen config. If it clears (or is
   marginal with clean powered SBC + the already-clean PPC) → package for
   manager re-adjudication as FAIL-ADJUDICATED-ACCEPTABLE or clean pass.
4. **B6 anchor-mode limits** for v5/v6 (reachable Im_k2 anchors).
5. Ship everything to origin; Machine A re-adjudicates, unwires gating for
   whatever ratifies, rebuilds the snapshot (registry-derived cache list is
   live as of `2613ed6d`), and the user ships to HF.

v6's expected path is short (it was already clean at HEAD — B1/B2/B6 are
confirmation). v5 is the real question; do not skip step 3's preregistration.

## 0.6 Manager advisory (2026-08-04) — Europa PPC batch + under-update diagnosis

Advice, not new scope — sequencing and design guardrails for §0.5/B5. The
user is HOLDING the HF redeploy until the v4/v1.1 PPC verdict, so that pair
is the critical path.

**Sequencing.** (1) PPC v4 + v1.1 (cheap: existing artifacts, references,
and forward-model machinery; generalize `plans/scripts/titanG_ppc_interior_check.py`
rather than reimplementing). (2) Report both, whatever they show. (3) Only
then the under-update diagnosis experiments (pilot-scale). (4) MgSO4/NaCl
compute starts only after the diagnosis explains the NH3 miss — training a
third artifact with an unexplained conditioning deficiency is waste.

**PPC design notes.**
- v4: push all 21 channels (gravity pair, k2/h2, 14 Bind via the Ae
  sidecar), not just k2 — the D–w degeneracy signal (adjudication doc) may
  surface in induction or gravity rather than k2. v1.1: CMR2 + the k2/h2
  hypothetical channels (label them as such in the report).
- Keep the pushforward NOISELESS and compare three medians per channel:
  prior-predictive (from the training dataset), SBI posterior-predictive,
  MCMC posterior-predictive, plus the observed value — the NH3 four-way
  table is the template.
- Preregister the reading BEFORE running: an SBI channel is flagged when
  |median_pp(SBI) − median_pp(MCMC)| > 0.5·sigma_obs. Sanity: this flags
  both NH3 k2 channels (Re 0.83σ, Im 1.46σ) and neither NH3 gravity channel
  (~0.03σ). Report the statistic for every channel either way; never move
  the 0.5 after seeing results. This same statistic is the candidate
  first-class pushforward gate for MgSO4/NaCl — validate it here first.

**Under-update diagnosis — ranked hypotheses to test at pilot scale
(~100k sims / small flows), cheapest-first:**
1. **Noise-augmentation swamping:** if training x pairs get noise drawn at
   sigma_obs while the k2 prior-predictive spread is only ~3–7x sigma, the
   flow may rationally down-weight those channels. Test: retrain a pilot
   flow with reduced (or zero) k2 noise augmentation, PPC it — if the
   under-update shrinks, the noise convention is the mechanism.
2. **Tail sparsity at the datum:** the Titan Im_k2 datum sits at the 86th
   prior-predictive percentile — few training pairs nearby. Test: PPC the
   pilot flow at synthetic x_obs placed at the 50th vs 85th vs 95th
   prior-predictive percentile of Im_k2; if update strength collapses with
   percentile, it is a coverage problem → remedy is a truncated/sequential
   second training round near the datum (TSNPE-style), not more uniform sims.
3. **x-normalization scale:** check the artifact's x_norm stats for the k2
   channels — z-scoring by a spread much larger than sigma_obs compresses
   exactly the dynamic range the conditional needs.
4. **Capacity/embedding:** last resort — deeper flow or a separate embedding
   for the low-dimensional obs vector; test only if 1–3 come back clean.

Report per hypothesis: tested / mechanism confirmed / ruled out, with the
pilot artifacts + seeds. Stop and surface anything surprising; do not fold
any fix into a production build without a frozen config from Machine A.

**RESOLVED 2026-08-04 (B) — under-update diagnosis: #3 + width ruled out,
mechanism localized to the ocean-admitting apparatus (NOT yet fully
mechanized). No retraining.** Two no-retrain experiments, cheapest-first per
the §0.6 ranking. Report: `validation_reports/FLOW_UNDERUPDATE_DIAGNOSIS.md`;
control PPC JSON at `titan_freegrav_noocean_1m/ppc/ppc_pushforward_report.json`
(seed 58, 0/4 flagged). Scientific-reviewer **PASS WITH CONCERNS** (2026-08-04),
all four required corrections folded.
- **#3 (x-norm scale): RULED OUT** by pure artifact inspection (reviewer
  reproduced every ratio exactly). NH3 k2 channels have the LOWEST
  train_std/σ_obs (~3.4) in their own vector (median 14) — the opposite of the
  compression #3 predicts — and k2 S/N (~3.2) matches/beats clean Europa k2.
  Corollary: the per-channel form of #1 is weakened by the same quantity.
- **Obs-vector width (4ch): RULED OUT as causal.** The Titan free-gravity
  NO-OCEAN artifact (same 4 channels/body/σ_obs/Petricca datum) assimilates a
  strongly-informative **Re_k2 (3.92σ MCMC-pp update) to 0.05σ** — a 4-channel
  flow resolves informative k2 fine. Kills the width reading positively.
- **Mechanism localized, NOT fully established.** The no-ocean control differs
  from NH3 by MORE than the mixture (12 vs 13 params — extra `log10_wOcean_ppt`;
  `phase_stability` enforce vs None), so it isolates the whole ocean-admitting
  apparatus (joint mixture AND the co-varying salinity axis + removed guard),
  not the mixture alone. Bimodal mode-assignment is a **candidate** (NH3 SBI
  Re_k2 median 0.541 sits near the ocean value, not dragged to the frozen 0.113;
  would need the NH3 SBI k2 pushforward shown bimodal to establish).
- **#1 (noise swamping): cleared for Re_k2 only; UNTESTED for the dominant
  Im_k2 miss.** The control is non-informative for Im_k2 (0.10σ), and NH3 flags
  Im_k2 harder (1.64σ) than Re_k2 (0.94σ). Im_k2 is where the abs-fold +
  additive-noise convention matters most → the #1 reduced-noise pilot remains a
  live cheap probe for Im_k2, NOT off the table.
- **MgSO4/NaCl:** inheritance is CONDITIONAL on the (not-yet-established)
  mechanism; verify the joint-build flag from actual cache metadata per build.
  **Compute stays HELD.** Two cheap separators before any remediation build:
  (a) an ocean-only-with-salinity PPC (isolates salinity axis from bimodality);
  (b) the #1 reduced-noise pilot for Im_k2. Proceed-as-split-status vs.
  separate-first vs. remediate is a **manager/user call** — surfaced, not acted
  on. Variance-subtraction under the abs-fold confirmed valid by the reviewer
  (noise post-fold, ⟂ θ, not re-folded).

## 1. Titan NH3 — JOINT no-ocean+ocean cache + SBI (Task #68 governs)

RECONCILED by Machine A 2026-08-02 per the user's firm decision (2026-07-30,
re-affirmed 2026-08-02) and Machine B's addendum
`plans/STATUS-2026-08-01-machineB-joint-nh3.md`. The earlier "Phase-0
rectangle re-run" framing here is WITHDRAWN: the provisional rectangle
Tb [248,257] K x w [30,100] ppt would re-condition the posterior on "an ocean
exists" and exceeds the 70 ppt NH3 ceiling. The corrected-activity-model
directive (`6c5ee2af`) COMPOSES with, not supersedes, the joint design — its
requirement is satisfied because the joint validations and cache are built at
a HEAD that includes the correction. Any pre-correction feasibility numbers
remain void.

- JOINT posterior over the FULL range Tb in [249, 263] K x w in [1, 70] ppt
  NH3. Frozen (Tb, w) nodes build as REAL no-ocean interiors via
  `build_tbw_grid_cache(..., retry_frozen_as_no_ocean=True)`, tagged
  `has_ocean=False`, and STAY in support — do not truncate Tb, do not drop
  frozen nodes.
- Config `configs/test54_titan_nh3_freegrav.json`: no
  `phase_stability.enforce` guard in either direction; density-inversion
  guard + k2 support bounds still apply.
- Continue Task #68 as planned: 3x4 joint validation -> scientific-reviewer
  interpretation -> author config -> USER go-ahead -> full-range production
  cache + 1M NH3 build + gates.
- Commit and push the builder edits (typed `NoIceLiquidTransitionError`,
  `do_overrides` channel, `retry_frozen_as_no_ocean` flag) once the 3x4
  validation confirms them — Machine A needs them on origin before wiring
  anything that reads the new cache metadata (`has_ocean`,
  `n_no_ocean_nodes`).
- Spec: `plans/active/titan-nh3-ocean-campaign-spec.md` (joint-posterior
  supersession note at top; per-phase mechanics otherwise still apply).

## 2. Cassini–Enceladus production

Status: queued after Titan Phase 0. Prior physics blockers closed
(hydrostatic ratio 3.25, correlated C20/C22, negligible elastic libration).

- base config: `PlanetProfile/Inference/configs/enceladus_cassini_smoke_6D.json`;
- constraints: `plans/mission-body-constraints.md`;
- dense Seawater 10 ppt Tb cache over ~271.8–272.6 K;
- cache-build overrides: `Cuncertainty = CuncertaintyUpper =
  CuncertaintyLower = 0.08`;
- production adds `dC20_nh ~ U[-1.5e-3, 1.5e-3]`, `dC22_nh ~ U[-1e-4, 1e-4]`;
- retain `observable_correlations = {"C20,C22": 0.47}`,
  `gravity_j2_over_c22 = 3.25`; libration is the primary decoupling channel.
- Phases: dense cache + inspection -> config freeze via Machine A ->
  reference MCMC + SBI + SBC/crosscheck/limits -> Ae sidecar -> ship all
  with seeds/versions/verdict. Do not deploy.

## 3. Callisto campaign (roadmap — no frozen config yet)

- Read `papers/cochrane2025stronger.pdf` before finalizing induction
  observables/uncertainties.
- Existing Callisto C2 structure cache has a P_MPa/r_m length mismatch
  (118 vs 122) — regenerate during setup.
- Do NOT inherit Europa's |Ae| support cuts: the Callisto prior must keep
  no-ocean/low-|Ae| AND saturated high-|Ae| corners inside support
  (Cochrane2025 posture; see 2026-07-25 audit note in git history of this
  file for the numbers).
- User directive: full coverage in ice thickness and ocean COMPOSITION;
  SeaFreeze 1.1.3 ships NaCl(aq) splines (NaClaq_LP/HP) — consider NaCl
  cache variants.

## Environment notes

- SeaFreeze pinned 1.1.3 (`pip install --upgrade SeaFreeze==1.1.3`); code at
  HEAD does not run against 1.0.0.
- CoolProp >= 8.0.0 required for NH3 work.
- Frozen reference env: `environment.yml` (PPcl, committed `11b514dd`).

## Queue boundary

No heavy campaign beyond the above is authorized. Anything else is roadmap
material until Machine A provides a frozen config and an explicit handoff.
Any `PlanetProfile/Test/` additions or `.gitignore` changes require explicit
user permission before commit.
