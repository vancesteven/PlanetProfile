# SBI Artifact Index

Last audited: 2026-08-02

## Naming convention

GUI slot filenames are registered in
`PlanetProfileApp/pages/Inference.py::_SBI_ARTIFACT_SLOTS`. The inventory
tables below are one-to-one: every shipped `.pt` file has its own row, including
ablation and diagnostic artifacts.

`SBIRunner.save_artifact` defaults to `<bodyname>_<config_hash>_posterior.pt`; a
validated artifact is deployed by registering its filename in a GUI slot. Artifacts
embed no path — renames are safe. Artifacts are self-describing `.pt` files (schema_version
1): posterior, prior_spec, param_names/units/bounds, obs_names, imag_convention,
normalization constants, config_hash, git_sha, seed, n_train_effective, rejection_stats,
sbi/torch versions, created_utc.

## Deployment gate

An artifact may occupy a GUI slot ONLY after passing all three ratified validation gates
(`PlanetProfile/Inference/validate_sbi.py`; thresholds in
`plans/monte-carlo-sbi-implementation-plan.md` §Validation):

1. `sbc` — per-parameter rank-uniformity KS p ≥ 0.05
2. `crosscheck` — vs the matching production MCMC pickle: |Δmean| ≤ max(0.25σ, σ/√ESS),
   σ-ratio ∈ [0.7, 1.4], marginal KS α = 0.01
3. `limits` — monotone log10_eta_Ih vs |Im k2|, samples inside prior box

Pre-deploy assertions: `obs_names` matches the training config's observables,
`imag_convention == 'abs'` (the only convention the GUI accepts), `param_names` matches
the config's sampled parameters.

## Deployed artifacts

| Slot file | Source config | config_hash | git SHA | seed | n_train | Gates | Deployed |
|---|---|---|---|---|---|---|---|
| titan_andrade_noocean_posterior.pt | test50_titan_noocean_andrade_8D.json | 629afbd55a4f0ce5 | bf7c938e | train 42 / data 44 / noise 4444 | 999,816 (nsf) | ALL GREEN within domain (amended rules, 2026-07-11) | 2026-07-12 (user-ratified; artifact pushed git 4d32809b, pre-deploy assertions green) |
| europa_seawater_andrade_posterior_1m.pt | europa_seawater_andrade_7D.json | a09396bcb0d0eff5 | 3d865dc1 | train 42 / data 44 / noise 4444 | 831,750 (nsf, synodic-only support cut) | 3 gates FAIL but ALL are flow-calibration slack, NOT physics; physics-discrimination (limits W1 anchors) GREEN all 6. Deployed under domain guard + default Tb truncation (see conditions below) | 2026-07-13 (user-directed deploy w/o Machine A; scientific-reviewer DEPLOYABLE-after-remedy; GUI-verified via AppTest)**RETIRING (user 2026-07-19):** trained with Titan-derived k2 sigmas; Europa has NO measured k2 — the artifact conditions on a fake measurement. Replacement: Galileo v1.1 (honest CMR2+support-cut data; k2/h2 as labeled hypothetical channels) per plans/europa-clipper-v4-geodesy-plan.md. **RETIRED 2026-07-20**: slot replaced in place by europa_galileo_v1p1_8D_posterior_1m.pt (below); artifact retained for provenance. |
| europa_seawater_andrade_clipper_v2.pt | europa_seawater_andrade_clipper_v2.json | 46be64069a40090f | f59f43b0 | train 42 / data 47 / noise 4747 | 831,566 (nsf, synodic \|Ae\|>0.7 support cut) | SBC PASS all 7 (min KS p .19); crosscheck PASS mean/median/sigma all 7, shape FAIL Tb_K only (D .119 vs .054); limits grid-walk W1 7/8 PASS + containment .989-.997 gate-artifact | **USER-RATIFIED 2026-07-18.** GUI slot live (config_path + derived \|Ae_synodic\| guard [0.75,0.94] + Im_k2<=0.15 + Tb>=261.5 default truncation + Ae-per-frequency conditioning inputs); AppTest-verified; same-version sampling reproduction PASS (v2_sampling_reproduction_report.json); spline resample guard + Im_h2 fold fix in. Tb crosscheck shape caveat dispositioned (sub-resolution, conservative-wide; see scope note).**RETIRING (user 2026-07-19):** same Titan-derived k2-sigma defect; superseded by v4 (Mazarico-projected geodesy). **RETIRED 2026-07-20**: slot removed in the commit that integrates the v4 slot (retirement schedule, user-ratified 2026-07-19); artifact retained for provenance. |
| europa_galileo_v1p1_8D_posterior_1m.pt | europa_galileo_v1p1_8D.json | a9d70df25dabdbda | 6e995165 | train 44 / data 50 / noise 5050 | 831,903 (nsf, synodic \|Ae\|>0.7 support cut) | **ALL GATES PASS**: SBC 8/8 (min KS p .067), crosscheck 8/8 (no soft fails); zeta split decouples ice/silicate heating (corr -0.02, MCMC & SBI identical); runner-path heating reproduces stored reference exactly | **DEPLOYED 2026-07-20** (slot swap for v1 per the retirement schedule). Honest-data framing: only CMR2 (GC21) + synodic support cut are Galileo-era; k2/h2 are labeled hypothetical channels (Re_k2 [0.23, 0.05], Im_k2 [0.004, 0.05]). Machine A end-to-end GUI sampling `verified`: 10k draws; all 8 marginal means match the stored crosscheck MCMC stats within 0.5 sigma_MCMC. Guards: Im_k2 <= 0.15 (conservative carryover; no v1.1 anchor walk), Tb >= 261.5 K default truncation. Gates: validation_reports/europa_galileo_v1p1_1m/PHASE_GATES.md. |
| europa_clipper_v4_geodesy_11D_posterior_1m.pt | europa_clipper_v4_geodesy_11D.json | 0b46ae3ce56f92d5 | 6e995165 | train 43 / data 49 / noise 4949 | 326,849 kept (nsf, 67.3% support rejection, 21-D x / 11-D theta) | **DELIVERABLE VERIFIED**: Gate 4 non-hydrostatic u = dC22_nh + dC20_nh/3.324 recovered 94-95%, likelihood-dominated (36x prior shrink), u-SBC calibrated (KS p .093, conservative); Tb-w joint PASS (corr -.987 vs ref -.983, containment .640/.920); flow-fidelity PASS (PPD KS p .121). SBC/crosscheck/degeneracy FAIL only on prior-dominated interior scalars (zeta_Ih, rho_core) + diffuse Tb - diagnosed intrinsic (robustk2 ablation ruled out preprocessing; Re_k2 0.23 sits a genuine, MCMC-reproduced 1.3 sigma below the joint PPD). Reviewer: PASS, close v4, no blocker | **USER-RATIFIED + DEPLOYED 2026-07-20** (user authorized HF ship). sampling `verified`. GUI slot live (config_path, inherited Im_k2 <= 0.15 + synodic \|Ae\| [0.75,0.94] guards, no Tb truncation - 2D support); Machine A end-to-end GUI sampling: 10k draws, corr(Tb, log10w) -0.986 vs committed reference -0.982, u sigma 3.35e-7 (= gate report), per-sample c20/c22 packaged, new u-panel renders (reportable upper limit; per-component marginals marked do-not-cite). Reference MCMC: Test/mcmc_results/Europa/Test53_geodesy_v4/. Gates: validation_reports/europa_clipper_v4_1m/PHASE_GATES.md. |

**Europa "Galileo run" deployment conditions (scientific-reviewer 2026-07-13, GUI-verified):**
- **Gate verdicts (raw):** SBC FAIL (alpha p=.048, Tb p=.028; other 5 pass, e.g. eta_Ih .66);
  crosscheck vs Test51 FAIL on Tb *shape* only (D=.090>tol=.057; Tb mean/median/sigma + alpha
  all PASS — Tb is UNBIASED); limits containment 0.9938<1.0 but **W1 anchors PASS all 6**
  (Im_k2 0.00-0.15, W1 .024-.077 vs tol ~.42, no rail pileup). Reports in
  `validation_reports/europa_galileo_1m/`.
- **Why deployable despite 3 FAILs:** (1) containment is a gate-measurement artifact — it
  samples `reject_outside_prior=False`; the GUI runtime samples `=True`, so the 0.62% box leak
  never reaches a user. (2) The Tb failure is edge-smear, not structural: a matched-truncation
  crosscheck (both flow+MCMC cut at Tb>=261.5 K) drops Tb KS D from .093 to **.019 (PASS, p=.51)**.
  The NSF flow cannot represent the hard one-sided synodic support edge (Tb<~261.5 K has no
  conductive ocean, removed at training by row-rejection) and smears ~2.5-3.5% of Tb mass into
  the excluded band. (3) alpha p=.048 is near-threshold multiple-comparison noise (passes the
  full independent MCMC crosscheck).
- **Deployed guards (Inference.py `_SBI_ARTIFACT_SLOTS`):**
  (a) `x_obs_limits={'Im_k2':(0.0,0.15)}` — hard refusal beyond the W1-validated domain
  (NARROWER than Titan's 0.20; no Europa anchor ran above 0.15).
  (b) `default_truncate={'Tb_K':(261.5,None)}` — pre-applied ON, restores the induction support
  edge at sample time (reject_outside_prior does NOT re-cut it: [259.5,261.5] is inside the prior
  box). Verified: keeps 97.5% of draws, zero post-truncation mass below 261.5 K.
  (c) `scope_note` documenting the synodic-only Galileo run + Tb caveat + v2 pointer.
- **GUI verification (2026-07-13, AppTest streamlit.testing.v1):** Europa slot lists in the
  selector; Tb truncation slider defaults to (261.5, 271.0); default-truncation banner + scope
  note render; the Im_k2<=0.15 guard fires and refuses conditioning at Im_k2=0.18.
- **v2 follow-on:** 3-frequency (add Ae_synodic 2nd + Ae_orbital) is a deliberate future artifact
  requiring a fresh 1M; NOT retrofit into the Galileo run.

**Cross-version sampling validation (2026-07-13, Machine A):** both deployed
artifacts were trained on torch 2.11 (Machine B) and load on torch 2.8 (Machine A)
with the loud version-mismatch warning by design. VALIDATED for the pair
torch 2.11.0 -> 2.8.0:
- Europa: re-running the crosscheck gate on Machine A reproduces Machine B's
  committed report to four decimal places on every per-parameter statistic
  (D, dmean, verdicts; Tb-only edge-smear fail included) — flow deserialization
  is bit-consistent across these torch versions at the gate seed.
- Titan: crosscheck gate vs the Test50 production MCMC PASSES all 8 parameters
  on Machine A (max D 0.048, all location tests green).
Reports: `validation_reports/cross_version_machineA/`. The GUI slot registry marks
this pair `validated_version_pairs`, demoting the load warning to INFO for these
artifacts only; any OTHER version pair still raises the loud RuntimeWarning
(correct default caution — re-run the crosscheck gate before trusting draws).

## Candidate artifacts (NOT deployed — evidence committed for cross-machine verification)

| Artifact file | Source config | config_hash | seed | n_train | Gates | Status |
|---|---|---|---|---|---|---|
| titan_diff_noocean_andrade_test52_10D_v2.pt | test52_titan_noocean_andrade_10D.json | 2bf1f7b2d1708e28 | train 42 / data 45 / noise 4545 | 877,883 (nsf, 10xIQR x-filter) | SBC FAIL (eta_V p=.016, 9/10 pass); crosscheck shape FAIL 3 params (zeta/eta_III/Tb, all location tests pass); limits W1 PASS Im<=0.25, FAIL Im=0.30 | **implemented, unverified** — Test52 10D (Titan no-ocean + CMR2). NOT deployed. Gate reports in `validation_reports/test52_v2/`. Reviewer (opus 2026-07-11): deployable with \|Im k2\|<=0.25 restriction + lower-tail shape caveats on zeta/eta_III/Tb. Awaiting Machine A GUI guard + deployment decision. |
| europa_seawater_v3_clipper_8D_posterior_1m.pt | europa_seawater_andrade_clipper_v3_8D.json | 41097cc6d4ef64c9 | train 42 / data 48 / noise 4848 | 336,537 (nsf, \|Ae_synodic\|>0.7 support cut, 66.3% rejected) | **2D DEGENERACY GATE PASS** (corr(Tb,log10w) MCMC −.986 vs SBI −.988, \|Δ\|=.002; joint 68/95% containment .650/.943; all 8 marginals in tol); SBC PASS (min KS p .118); crosscheck 7/8 PASS, rho_core shape-only fail (NON-BLOCKING: least-identified, SBI conservatively wider, SBC clean); limits gate premise inapplicable (ocean-bearing) | **RATIFIABLE** (scientific-reviewer 2026-07-19: PASS WITH CONCERNS). v3 = v2 + sampled log10_wOcean_ppt U[−1,2] (0.1–100 ppt) over a 2D (Tb × log10 w) bilinear structure grid. GUI slot live on Machine A (render `verified` via AppTest: slot lists, scope note + Tb–w degeneracy text, Im_k2≤0.15 + inherited v2 Ae guard [0.75,0.94]); 2D cache committed (f7a03572, 21 MB, 1488 nodes); **sampling `verified` 2026-07-19** — end-to-end GUI AppTest (Generate Posterior → packaged result): 10k draws/4 min, 8 params incl. log10_wOcean_ppt, per-sample k2/h2/Ae×3 packaged, per-sample cloud plot rendered; posterior matches reference MCMC (Tb 264.46 vs 264.3 K; w 36.6 vs 38.7 ppt; corr −0.988 vs −0.986). Reference MCMC committed: Test/mcmc_results/Europa/Test52_seawater_v3/. Gates: `validation_reports/europa_clipper_v3_1m/PHASE3_GATES.md`. **VETOED by user 2026-07-19** — wrong training bounds: the k2 conditioning sigmas (0.06) are requirement-level, not the Mazarico et al. (2023) Clipper projections ((1.4–1.8)e-2), and Re_k2 central 0.25 vs the ocean-consistent 0.23. NOT deployed; GUI slot removed. Superseded by the v4 geodesy campaign (plans/europa-clipper-v4-geodesy-plan.md), which reuses this artifact's 2D (Tb × w) cache, reference-MCMC infrastructure, and gate pre-registrations unchanged. Artifact + gates retained for provenance. |
| europa_clipper_v4_geodesy_11D_posterior_1m_robustk2.pt | europa_clipper_v4_geodesy_11D.json (robust-k2 preprocessing ablation) | 0b46ae3ce56f92d5 | train 43 / data 49 / noise 4949 | 326,849 | Same-data diagnostic left gate verdicts unchanged; see `validation_reports/europa_clipper_v4_1m/PHASE_GATES.md` and `reviewer_diagnostics_123_robustk2.json`. | **Diagnostic ablation, not deployed.** Retained as evidence that k2 z-score compression was not the driver of v4 residuals. |
| europa_clipper_v5_geodesy_11D_posterior_1m.pt | europa_clipper_v5_geodesy_11D.json | 0781749945edc44b | train 51 / data 50 / noise 5050 | 381,734 | SBC PASS; limits FAIL; crosscheck FAIL. See `validation_reports/v5_gate_summary.json` and `validation_reports/europa_clipper_v5_baseline_1m/`. | **delivered, not ratified (gate adjudication open; see `plans/STATUS.md`).** Baseline arm; not wired for deployment. |
| europa_clipper_v5_noinduction_7obs_posterior_1m.pt | europa_clipper_v5_noinduction_7obs.json | bbf023a694aa95b0 | train 51 / data 50 / noise 5050 | 381,734 | SBC PASS; limits FAIL; crosscheck not applicable (no ablation reference MCMC). See `validation_reports/v5_gate_summary.json`. | **delivered, not ratified (gate adjudication open; see `plans/STATUS.md`).** No-induction ablation; not wired for deployment. |
| europa_clipper_v5_nok2_17obs_posterior_1m.pt | europa_clipper_v5_nok2_17obs.json | 077875a481f337b0 | train 51 / data 50 / noise 5050 | 381,734 | SBC FAIL; limits skipped (no k2 channel); crosscheck not applicable (no ablation reference MCMC). See `validation_reports/v5_gate_summary.json`. | **delivered, not ratified (gate adjudication open; see `plans/STATUS.md`).** No-k2/h2 ablation; not wired for deployment. |
| europa_clipper_v6_freegrav_11D_posterior_1m.pt | europa_clipper_v6_freegrav_11D.json | a8cd890fff6a0023 | train 61 / data 67 / noise 6767 | 381,093 | SBC PASS; limits FAIL; crosscheck FAIL. See `validation_reports/v6_gate_summary.json` and `validation_reports/europa_clipper_v6_baseline_1m/`. | **delivered, not ratified (gate adjudication open; see `plans/STATUS.md`).** Baseline free-gravity arm; not wired for deployment. |
| europa_clipper_v6_freegrav_noinduction_6obs_posterior_1m.pt | europa_clipper_v6_freegrav_noinduction_6obs.json | 7075e2ae5dca104b | train 61 / data 67 / noise 6767 | 381,093 | SBC FAIL; limits FAIL; crosscheck not applicable (no ablation reference MCMC). See `validation_reports/v6_gate_summary.json`. | **delivered, not ratified (gate adjudication open; see `plans/STATUS.md`).** No-induction ablation; not wired for deployment. |
| europa_clipper_v6_freegrav_nok2_16obs_posterior_1m.pt | europa_clipper_v6_freegrav_nok2_16obs.json | 70b8e622dd28c807 | train 61 / data 67 / noise 6767 | 381,093 | SBC PASS; limits skipped (no k2 channel); crosscheck not applicable (no ablation reference MCMC). See `validation_reports/v6_gate_summary.json`. | **delivered, not ratified (gate adjudication open; see `plans/STATUS.md`).** No-k2/h2 ablation; not wired for deployment. |
| europa_clipper_v7_openae_11D_posterior_1m.pt | europa_clipper_v7_openae_11D.json | 6b77483322378820 | train/data/gates 71 | 955,816 | SBC PASS; limits FAIL; crosscheck FAIL. See `validation_reports/v7_gate_summary.json` and `validation_reports/europa_clipper_v7_openae_1m/`. | **delivered, not ratified (gate adjudication open; see `plans/STATUS.md`).** Open-\|Ae\| baseline; not wired for deployment. |
| titan_freegrav_noocean_posterior_1m.pt | titan_freegrav_noocean.json | 9600b13e519baea6 | train/data/gates 71 | 877,227 | SBC PASS; limits FAIL; crosscheck FAIL. See `validation_reports/titan_freegrav_noocean_1m/titanG_gate_summary.json`. | **delivered, not ratified (gate adjudication open; see `plans/STATUS.md`).** Free-gravity Titan candidate; not wired for deployment. |

**Europa Clipper v2 scope note (scientific-reviewer 2026-07-14, COMMIT-AS-CANDIDATE):**
- **What v2 is:** the 3-frequency Clipper-era follow-on to the deployed synodic-only "Galileo run"
  (`europa_seawater_andrade_posterior_1m.pt`). Adds 14 induction observables
  `Bind_<label>_<comp>_<part>` for label in {synodic, synodic 2nd, orbital}, comp in {x,y,z}
  (kept per pruning: synodic xyz, synodic-2nd xy, orbital xy), part in {real, imag}, each
  sigma=1.5 nT (Kivelson et al. 2023). Bind = Ae*Be (unconjugated, complex per-component Be,
  FT surface-field path); signed Im (not abs-folded). See config metadata for the full
  convention + the surface-field-vs-g1 factor-of-2 assumption.
- **Gate verdicts (raw):** SBC (1000 pairs) PASS all 7 params (min KS p .19; Tb_K p .32,
  c2st .573). Crosscheck vs a fresh v2 reference MCMC (nautilus, r_hat 1.0, ESS 4259):
  FAIL on Tb_K *shape* only (D .119 > tol .054); Tb_K mean/median/sigma PASS (dmean .020,
  dmedian .010, sigma_ratio 1.184 -> SBI 18% WIDER = conservative), and all 6 other params
  PASS fully. Limits grid-walk W1 (8 on-manifold anchors, Tb 261.5-271 K): 7/8 W1 PASS,
  medians agree <=.04 K; containment .989-.997 < 1.0.
- **Why committable despite 2 shape/containment FAILs:**
  (1) The Tb_K crosscheck fail is NOT the v1 cold-edge smear: matched-truncation at Tb>=261.5 K
  is a no-op (0.0000 sub-support mass on BOTH flow and MCMC). It is a genuine within-support
  higher-moment defect — the reference Tb marginal is right-skewed/leptokurtic (skew +.70) and
  the flow renders it near-Gaussian (skew +.20). Confirmed via a dual-reference MCMC:
  MCMC(seed42)-vs-MCMC(seed7) Tb_K D=.035 ~ self-D floor .035 (both chains agree on the skew),
  so the reference is stable and the D=.119 is a real FLOW shape defect. It is sub-resolution
  (medians agree .01 K ~ .09 km ocean; tail mismatch ~1.3 km ~3.5% of ocean) and in the
  conservative (wider) direction. SBC PASS confirms the flow is calibrated for Tb across the
  prior; the fail is at a single pinned-Tb x_obs.
  (2) The single grid-walk W1 miss (Tb=263.5, W1 .0442 vs tol .0344) is the SAME flow
  shape-fidelity floor, not a discrimination failure: flow_median tracks anchor_median to
  .003 K there. It trips only at 263.5 because that anchor has the smallest sigma_anchor (.138 K)
  hence tightest .25*sigma tol. Disambiguated with W1: MCMC-vs-MCMC W1 at 263.5 = .0066 K
  (reference noise floor) vs flow-vs-MCMC .0442 K -> a ~.038 K flow shape increment, constant
  across mid-range anchors, sub-resolution. Report `limits_263p50_w1_disambiguation.json`.
  (3) Containment .989-.997 is the SAME gate-measurement artifact as the deployed v1 (v1: .9938):
  it samples `reject_outside_prior=False`; the GUI runtime samples `=True`, so the diffuse 0.72%
  NSF-spline tail leak (spread across all params, Tb_K leak 0.00%, no single param >.25%) never
  reaches a user.
- **Anchor design correction (scientific review, superseded plan):** the plan's pre-registered
  independent-per-frequency |Ae| sweep was REJECTED — Ae is a function of Tb ALONE, so the three
  excitation frequencies cannot be moved independently, and off-manifold synodic targets rail the
  reference MCMC against the hard |Ae_synodic|>0.7 cut. Replaced by single-Tb GRID-WALK anchors
  (all frequencies co-vary from one physical Ae(Tb)) as the primary W1 gate. The plan's set was
  still run as labeled EXTRAPOLATION PROBES and empirically confirms the rejection: off-manifold
  synodic probes rail (synodic_0.70 W1 .61, 0.75 .24, 0.90 .24, corner .27), and the orbital/
  synodic-2nd probe passes are INERT (synodic pins Tb; the reference barely moves) — NOT
  independent flow-quality evidence.
- **Carry-forward deploy conditions (Machine A, before any deploy — NOT done here):** runtime
  `reject_outside_prior=True`; an x_obs domain guard set to the W1-validated |Ae| envelope
  actually exercised by the grid-walk (synodic ~0.75-0.94; do not condition beyond it); Tb
  default-truncation at the induction support edge (261.5 K), as in v1; an AppTest asserting
  Bind inputs are NOT abs-folded (signed Im preserved end-to-end). Deploy gated on GUI slot +
  AppTest + user ratification. Cross-version load: trained torch 2.8.0 / sbi 0.26.1; re-run the
  crosscheck gate on Machine A before trusting draws (per the version-pair discipline above).

**Scope note (deployment condition, user-ratified 2026-07-11):**
- **Validated conditioning domain: |Im k2| <= 0.20.** SBC (1000 pairs) PASS; crosscheck
  vs the Test50 production MCMC PASS on all 8 parameters under the ratified shape gate
  (max D 0.048 vs tol 0.055, max |dmedian| 0.085 dex); W1 anchor gate PASS at
  Im = 0.05/0.10/0.15/0.16/0.17/0.18/0.20/0.30.
- **Known limitation:** directional low-viscosity (eta_Ih) bias in the bimodal posterior
  regime above Im k2 ~ 0.18; the Im = 0.25 anchor fails (W1 = 0.473 vs tol 0.349) —
  the flow overweights the low-eta mode there. The GUI enforces a HARD guard refusing
  conditioning at Im k2 > 0.20 (`x_obs_limits` in the Inference-page slot registry);
  use MCMC mode beyond the domain.
- Previous provisional 500k artifact (git 278c3bea, data seed 43/4343) superseded;
  retrievable from git history.

## Pipeline pointers

- Training runner: `PlanetProfile/Inference/sbi_runner.py` (`SBIRunner`; single-round NPE,
  MAF). Dataset generation delegates to `MCMCRunner.generate_sbi_dataset`
  (`mcmc_runner.py`) with opt-in kwargs `apply_support_guard`, `imag_convention='abs'`,
  `drop_nonfinite`, `seed`, provenance.
- libomp hazard: do NOT run dataset generation (PlanetProfile) and torch training in the
  same process. Generate `(theta, x)` to `.npz` in a PlanetProfile-only process, then
  train/save in a torch-only process. Avoid `SBIRunner.run()` end-to-end on macOS conda
  environments.
- Validation driver: `PlanetProfile/Inference/validate_sbi.py`
  (sbc / crosscheck / limits / selftest; exit 0 pass, 2 gate failure).
- GUI: `PlanetProfileApp/pages/Inference.py::_SBI_ARTIFACT_SLOTS` — slot
  availability follows the registry and artifact presence; non-`'abs'`
  artifacts fail loudly.
