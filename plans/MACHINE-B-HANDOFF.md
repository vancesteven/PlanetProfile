# Machine B handoff

Updated: 2026-08-01 at genai `54106fbd` (Machine A refresh after v5/v6/v7
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
honest basis is NOT "wide σ = clean" (that framing was retracted — it conflates
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
  arms are the genuinely controlled probe for the obs-vector-width hypothesis
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
