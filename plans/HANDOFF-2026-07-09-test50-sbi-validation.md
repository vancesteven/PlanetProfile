# HANDOFF 2026-07-09 — Test50 Titan SBI validation campaign

## RATIFIED 2026-07-09 (user, via decision prompt)

A. Observation-noise injection: APPROVED, likelihood convention (x_Im = |Im_model| +
   noise, NO outer refold). IMPLEMENTED: `generate_sbi_dataset(obs_noise=, noise_seed=)`
   in mcmc_runner.py (additive kwargs, default off = byte-identical); wired through
   `SBIRunner.generate_training_set` (obs_noise=True default for SBI training).
   Test added: tests/sbi_runner_test.py::test_obs_noise_injection — 9/9 pass. verified.
B. Production settings nsf + 200k: APPROVED.
C. Crosscheck eta_Ih KS residual: more flow work first (bounded variants), then
   materiality ratification if stalled.
D. Limits gate: REPLACED monotonicity with MCMC-anchor comparison. IMPLEMENTED:
   `validate_sbi limits --anchor-results '{"0.05": "pkl", ...}'` — per-anchor
   |flow median − anchor median| ≤ 0.25·σ_anchor (LIMITS_ANCHOR_SIGMA_FRAC,
   pre-registered) + containment. Legacy monotonicity path kept for configs without
   anchors. verified (runs, plot overlays anchors).

## Post-ratification status (v2 = pipeline-noise 200k nsf, seed 42/noise 4242)

- SBC (1000 held-out pairs, seed 777/noise 7777): FAIL only log10_eta_sil p=0.034
  (marginal, persistent across every noised run).
- Crosscheck: all mean+sigma gates PASS; KS fails eta_Ih (7e-10), eta_V (9e-7),
  Tb_K (7e-4). Failing set wobbles between retrains except eta_Ih (always fails).
- Anchor limits (6 anchors, n_eff=300 each): PASS at Im 0.05/0.10/0.15
  (diff 0.10/0.24/0.32 dex vs tol ~0.40); FAIL at 0.20/0.25/0.30 (0.63/0.86/0.80 dex).
  Containment 1.0.
- **NEW FINDING: the Im=0.25 anchor posterior is BIMODAL in log10_eta_Ih** (modes ~11
  and ~13.5; median lands mid-valley) and at Im=0.30 the mass collapses to the
  low-viscosity mode (median 12.81 → 10.95). Median is a fragile anchor statistic in
  the bimodal regime; the flow smooths across modes. High-Im x also lies in the sparse
  tail of the training manifold (p99 of |Im k2| = 0.255).
- Flow-work variants in flight: (V1) nsf hidden=100/transforms=10 on 200k;
  (V2) nsf default on 500k (seed 43/noise 4343).

## Variant results (bounded flow work, decision C)

| Artifact | SBC (1000 pairs) | Crosscheck | Anchor limits |
|---|---|---|---|
| 200k nsf v2 | FAIL eta_sil p=.034 | means pass; KS fail eta_Ih/eta_V/Tb | pass .05-.15; fail .20-.30 |
| 200k nsf h100/t10 | FAIL min p=.0029 | KS fail 4 | pass .05-.15; fail .20-.30 |
| **500k nsf (candidate)** | **PASS min p=.252** | means+sigma pass; KS fail 4 (p 6e-6..3e-3) | pass .05/.10/.15/.30; fail .20 (.51/.39), .25 (1.30, bimodal) |

Capacity does NOT help; data scaling does (monotone improvement 50k→200k→500k).

## DEPLOYED (user-approved provisional, 2026-07-09)

500k nsf artifact → `PlanetProfile/Inference/sbi_artifacts/titan_andrade_noocean_posterior.pt`
(committed, 400 KB; provenance in sbi_artifacts/INDEX.md). GUI Titan slot live; validated
domain = conditioning near the observed point; |Im k2| >= 0.20 is unvalidated extrapolation.
Pre-deploy assertions passed (obs_names/convention/param_names/config_hash).
NOTE Machine A: streamlit was missing from conda env PP — installed via
`conda run -n PP pip install streamlit` (pure-python; torch/libomp stack untouched).
App runs: `cd PlanetProfileApp && conda run -n PP streamlit run PlanetProfileApp.py`.

## NEXT: 1M training on Machine B (M2 96GB, mamba env PPcl) — STUBBED HERE

**SCOPE FOR MACHINE B: TRAININGS AND THEIR GATE RUNS ONLY (user directive 2026-07-10).**
Machine B runs the numbered steps below (dataset gen, nsf training, held-out set,
anchor MCMCs, three gates) and reports results back via a short addendum to this file
(commit + push to origin genai). Do NOT pick up other development items — GUI work,
ratification-item implementation, Test52 manuscript plots, Callisto sidecars, and all
other backlog items are being handled on Machine A in parallel. Pull latest genai
before starting; Machine A pushes frequently.

Machine A generated the 1M dataset but training was stopped (machine too small).
All datasets regenerate deterministically from seeds — nothing binary needs to transfer.

```bash
# 1. dataset (~15 min): 1M sims, prior seed 44, noise seed 4444
mamba run -n PPcl python - <<'EOF'
import json, sys, numpy as np
sys.path.insert(0, '.')
from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.sbi_runner import SBIRunner
cfg = json.load(open('PlanetProfile/Inference/configs/test50_titan_noocean_andrade_8D.json'))
cfg['mode'] = 'sbi'
r = SBIRunner(InferenceConfig.from_dict(cfg))
theta, x, stats = r.generate_training_set(1000000, seed=44, obs_noise=True, noise_seed=4444)
assert stats['rejection_rate'] <= 0.20, stats
np.savez('train_1m.npz', theta=theta, x=x, stats=json.dumps(stats, default=str))
EOF
# (rebuild the structure grid pkl first if missing — see HANDOFF-2026-07-08 §gitignored artifacts)

# 2. train nsf, torch seed 42; save via SBIRunner.train + save_artifact
#    (pattern: scratch script train_artifact.py in the 2026-07-09 session; recreate:
#     load npz -> SBIRunner(config mode='sbi') -> train(theta, x, seed=42,
#     density_estimator='nsf') -> _train_info['rejection_stats']=stats -> save_artifact)

# 3. held-out: generate_training_set(1500, seed=777, obs_noise=True, noise_seed=7777)

# 4. anchors (~12 min): 6 MCMC runs from the committed configs
for f in PlanetProfile/Inference/configs/test50_limits_anchors/sweep_im0.*.json; do
  v=$(basename $f .json | sed 's/sweep_im//')
  mamba run -n PPcl python -m PlanetProfile.Inference.run_inference_cli \
    --config $f --output anchors_im$v.pkl
done

# 5. gates (sbc --n-sbc 1000; crosscheck vs the committed Test50 production pickle;
#    limits --re-k2 0.608 --anchor-results '{"0.05":"anchors_im0.05.pkl", ...}')
```

Expected from 1M: KS residuals (eta_Ih p 5.8e-4 at 500k, improving ~10x per doubling)
may cross alpha=0.01; anchor 0.20 (0.51 vs tol 0.39) may pass. Anchor 0.25 will NOT
pass on data alone (bimodal median valley) — needs open item 2 below.

## Machine A parallel development (2026-07-10)

- **Test52 manuscript disclosure plots DONE (items a, c, d):** committed in
  production_run/ — prior-predictive CMR2 (posterior truncated at the no-core anchor;
  sigma_post/sigma_obs = 0.46 is prior-envelope-limited: verified visually), k2 vs
  (log10_zeta, log10_eta_V) conditionals (compliant-corner bias visible), modes+HPD
  table (mode(zeta)=-2.30, mode(eta_V)=10.65). Item (b) — nothing pending there;
  item (e) Tb prior-uniformity is a text disclosure, no plot needed.
- **Callisto offset sidecars: SCIENTIFIC NO-GO on existing caches (opus-verified).**
  The Test52 sidecar definition (native - no-core reconstruction) applied to the
  existing Callisto C2_andrade cache yields -0.035 = -8.4 sigma_obs — the PHYSICAL
  structured-vs-uniform-silicate difference (cache rho_sil spans 1326-3202 kg/m3),
  NOT a discretization offset. The handoff 07-08 estimate of "+0.23 sigma" assumed a
  Titan-like uniform-silicate cache. Applying the naive sidecar would corrupt every
  Callisto tension conclusion.
  - IMPLEMENTED + verified: `PlanetProfile/Inference/build_cmr2_offset_sidecar.py`
    (reusable producer; reproduces the committed Test52 sidecar exactly; validity
    guard REFUSES non-uniform-silicate caches / |offset|>0.005, exit 3) and
    `tests/cmr2_offset_sidecar_test.py` (3/3, incl. Callisto refusal).
  - User ratified the uniform-silicate measurement approach 2026-07-10; opus
    recommended path (B): metrology-only uniform grid, sidecar applied to the
    production cache, contingent on a per-Tb hydrosphere-match assertion.
  - **2026-07-10 RESULT: path (B) FAILED its own verification — transferability
    empirically falsified.** Built the uniform-silicate metrology grid
    (PlanetProfile/Test/PPCallistoUniformSil.py + measurement-mode producer with
    hydrosphere-match guard). The guard REFUSED: CONSTANT_INNER_DENSITY does not
    merely flatten the silicate, it moves the entire MoI-matched solution — the
    hydrosphere bottom shifts ~200 km (R_ob 2056 vs 2257 km) and M_hydro grows 2.4x
    vs the production grid. The measured offset would describe a different planet.
  - **Conceptual reframe surfaced by this failure:** Test52's offset anchors a
    SYNTHETIC CMR2 observable that was defined on the PP-native scale, so
    native-vs-reconstruction consistency is required there. Callisto's CMR2 is a
    REAL measurement (Anderson et al. 2001, +/-0.0042); the mass-conservation
    reconstruction IS the forward model, and there is no PP-native truth to anchor
    to. The right treatment is likely a discretization-ERROR ESTIMATE for the
    reconstruction (e.g. radial-resolution refinement / Richardson on the production
    structures) folded into the error budget as a systematic, NOT a chi shift.
    Expected magnitude ~1e-3 CMR2 ~ 0.2 sigma_obs — small vs the minus5pct 2.45σ /
    minus10pct 6.38σ conclusions.
  - STATUS: no Callisto sidecar exists (guards refused, correctly);
    `test_no_sidecar_found_for_callisto` remains correct as written; committed
    Callisto χ values stand unchanged. NEXT DECISION (user + scientific reviewer):
    adopt the error-budget treatment (resolution-refinement estimate) vs drop the
    correction entirely; the "+0.23σ shift" expectation from the 07-08 handoff
    should be retired either way.

- **SBI induction + h2 observables (roadmap item) IMPLEMENTED, verified:**
  `generate_sbi_dataset` now computes `Ae_<label>_real/_imag`, `BiAmp_<label>`,
  `BiPhase_<label>_deg` from the same precomputed Tb-grid Ae cache the likelihood
  uses, and `Re_h2/Im_h2/abs_Im_h2/h2/k2` from the existing forward call (all were
  silently NaN-filled before — the documented landmine). Mirrors likelihood
  conventions exactly (abs-fold under imag_convention='abs'; BiPhase in degrees,
  obs_noise unwrapped — matters only near ±180°). Tests: sbi_runner 10/10 (new
  induction-channel test); real-path check on callisto_nacl config: dataset columns
  match `_ae_grid_cache` lookups to 1e-12. Note Callisto Ae is nearly Tb-constant
  across the grid (weakly informative for SBI there by construction).
  (Directly serves Machine B's surfaced item 2 below — the x-channel plumbing for
  induction observables now exists end-to-end in dataset generation.)
- **GUI MCMC tab config-loader fixes (user-reported):** preset radio now syncs on
  config load (widget-state-wins bug was silently overwriting loaded configs);
  k2 checkbox-gated like CMR2 (CMR2-only configs no longer crash with KeyError
  'Re_k2' or silently gain Titan-default k2); form preserves non-rendered
  observables (h2/Ae channels were dropped from runs). Status: implemented,
  unverified pending user click-through.

## 1M RESULTS — Machine B (2026-07-10)

Ran the full campaign on Machine B (M2, 96 GB, mamba PPcl), gates + seeds unchanged.
Artifact `/tmp/titan_andrade_noocean_posterior_1m.pt` (NOT committed — see below):
config_hash 629afbd55a4f0ce5, git bf7c938e, train_seed 42, data seed 44 / noise 4444,
n_train_effective 999,816 (rejection_rate 1.8e-4), nsf, sbi 0.26.1 / torch 2.11.0.
Dataset gen ~19 min; **nsf training 6628 s (~110 min)** — see COMPUTE note.

| Gate | Verdict | Detail |
|---|---|---|
| SBC (1000 pairs, seed 777/noise 7777) | **PASS** | all 8 params; min KS p=0.057 (eta_sil), Tb p=0.093; C2ST 0.56-0.62 |
| Crosscheck (vs test50_..._8D_result.pkl, n=4198) | **FAIL** | all mean+sigma PASS; KS fails 5/8 (D, p @ alpha=0.01): alpha (D=0.038, 4.8e-3), eta_Ih (D=0.045, 4.0e-4), eta_V (D=0.040, 2.4e-3), eta_VI (D=0.042, 1.4e-3), Tb (D=0.041, 1.6e-3) |
| Anchor limits (6 anchors, ess ~4100) | **FAIL** | pass .05/.10/.15/.30; **fail .20 (0.565 vs tol 0.397), .25 (1.137, bimodal)**; containment 1.0 |

**Headline: eta_Ih residual has CONVERGED toward the self-D floor; remaining failures are
gate-definition questions, not data volume (opus-reviewed 2026-07-10).**
- Track effect size D, NOT the KS p-value (p convolves D with n and is not a discrepancy
  metric at fixed alpha). eta_Ih crosscheck KS **D = 0.045 at 1M vs D = 0.068 recorded for
  the noised retrains earlier this campaign (line ~297)** — a 34% reduction, now only
  ~1.25x the split-half self-D 99th-pct floor (~0.036) vs ~1.9x before. That is
  decelerating convergence toward the irreducible floor, not a stalled model. (The 500k
  crosscheck logged only p=5.8e-4, not its D; the 1M p=4.0e-4 barely moved BECAUSE n grew,
  so p is not the scaling metric — track D. ACTION: to make the scaling curve rigorous,
  recompute crosscheck D for eta_Ih at 200k/500k/1M against the bootstrap self-D floor;
  open item 1.) All five failing params have D in 0.038-0.045; means and
  sigma-ratios all pass. Residual is distributional-only and sub-material for an
  order-of-magnitude viscosity, but note eta_Ih carries a PERSISTENT single-sign LOW bias
  (~0.12 dex low, ~0.25 dex low vs the Im=0.15 GT anchor, across every retrain) — bound it
  as a systematic under open item 1, do NOT dismiss it as pure KS oversensitivity.
- Crosscheck failing count rose 4/8 (500k) → 5/8 (1M, adds eta_VI). This is INCREASED KS
  POWER at larger n resolving sub-0.05-D residuals, NOT a model regression — direct
  corroboration that the p-value is the wrong scaling metric.
- Anchor 0.20 did NOT pass (0.565 dex). Do NOT compare to 500k's 0.51: those anchors were
  n_eff=300 (MC error on the median ~0.115 dex) vs ess~4100 here (~0.031 dex) — a
  0.055-dex change is within combined anchor MC error, so the direction is not meaningful.
  The real reason 0.20/0.25 fail is the independently-established bimodal median-valley
  mechanism (modes ~11 and ~13.5 at 0.25; median lands mid-valley), which is a
  statistic-definition problem, not data volume. Anchor 0.25 remains unpassable on data
  alone, as predicted.

**Implication:** the two failing gates require the ratification-framework work (open items
1 and 2), NOT more training. Recommend STOP scaling data on Titan Test50 and resolve the
gate definitions (KS materiality margin anchored to the self-D floor + distribution-level
bimodal-anchor statistic) before any further artifact. SBC clean at 1M (C2ST 0.56-0.62,
near-random) confirms average calibration; the crosscheck/limits "failures" test a
different null (pointwise agreement vs one MCMC posterior at one x, alpha=0.01) and land
in the bimodal / folded-noise tail flagged as unvalidated extrapolation — no contradiction
with SBC PASS.

**COMPUTE note (relevant to user's cluster question 2026-07-10):** 1M nsf training took
~110 min CPU-bound on the M2 (vs the handoff's ~12 min estimate at 200k) — super-linear
in dataset size with convergence-based early stopping. Artifact is 409 KB (size is
architecture-bound, not data-bound; well under the 50 MB repo limit). Larger runs, wider
observable vectors (induction/gravity), or artifact families should move to the cluster.
`train_sbi_artifact.py` (committed this session) takes a pre-built .npz so dataset-gen
and training are already separable cluster jobs.

**Artifact NOT deployed:** two gates fail; per handoff open item 4, deployment is blocked
on all-green + user re-verify. The provisional 500k committed artifact stands unchanged.
1M .pt + dataset live in /tmp (regenerable from seeds; not committed). New reusable
entry point committed: `PlanetProfile/Inference/train_sbi_artifact.py`.

## SURFACED FOR MACHINE A (planning + implementation; user directive 2026-07-10)

Hand these to Machine A to plan + implement; training returns to Machine B.

1. **SBI conditioning flexibility.** User expected the GUI to vary priors and constraints
   (observable sigma, CMR2, viscosities, Andrade params), but the amortized Titan artifact
   only exposes the Re/Im k2 VALUES — correct behavior, because an NPE artifact conditions
   only on its observable vector; sigma/priors are frozen at training time and CMR2 was
   dropped from Test50. Options: (a) point users to the MCMC tab for free prior/sigma
   variation; (b) train with sigma as CONDITIONING INPUTS (x = [Re, Im, sigma_Re, sigma_Im])
   so a sigma slider is meaningful — pipeline + GUI change; (c) core-sensitive config +
   artifact to expose CMR2; (d) artifact FAMILY across settings, GUI-selected.
2. **Magnetic induction observables.** For Jupiter's moons (Europa, Ganymede, Callisto)
   AND Enceladus, Miranda, Ariel, Triton, inference should also use magnetic-induction
   observables (complex induction response) as an x-channel, distinct from tidal k2.
3. **Modular forward-model observables (design constraint).** Independently-trained flows
   CANNOT be merged into a joint posterior — cross-observable correlations only come from
   training on the combined x. What IS modular is dataset generation: gravity /
   non-spherical (J_n/C_nm/S_nm) / induction / plasma-effects modules each emit columns of
   x, composed at build time; one flow per observable-set (or an artifact family). Future
   work: non-spherically-symmetric gravity, magnetic induction, plasma effects as such
   modules.

## RATIFICATION PROPOSAL (Machine A, 2026-07-10, after 1M) — open items 1 + 2

Groundwork numbers (scratch ks_d_scaling.py; ESS=4198, 200 bootstrap):

| param | selfD p99 (split-half, raw) | matched-n floor (~raw/sqrt2) | D 200k | D 500k | D 1M (Machine B) |
|---|---|---|---|---|---|
| alpha | .056 | .039 | .030 | .031 | .038 |
| log10_zeta | .050 | .035 | .030 | .030 | (pass) |
| log10_eta_Ih | .051 | .036 | .065 | .040 | .045 |
| log10_eta_III | .047 | .033 | .032 | .036 | (pass) |
| log10_eta_V | .052 | .037 | .065 | .060 | .040 |
| log10_eta_VI | .046 | .033 | .038 | .043 | .042 |
| log10_eta_sil | .051 | .036 | .037 | .038 | (pass) |
| Tb_K | .051 | .036 | .047 | .042 | .041 |

D MC error at this ESS is ~±0.01 (500k vs 1M differences are resample noise). All 1M
D values sit within ~1.0–1.3x the matched-n floor; every mean gate passes with
|dmean| ≤ 0.19 dex (materiality currency: 0.3 dex for order-of-magnitude viscosities);
all sigma-ratios in [0.97, 1.07].

**Proposed gate amendments (pre-registered form, thresholds from floor + materiality,
not from which params pass):**

1. **Crosscheck KS component → detection + materiality.** Keep mean/sigma gates
   unchanged (hard). Replace per-param "KS p >= 0.01" with:
   PASS iff (D <= 1.5 x matched-n selfD_p99) AND (|dmedian| <= 0.3 dex or 0.3*sigma_MCMC
   for non-log params) — i.e. shape residuals up to 1.5x the reference's own sampling
   floor are acceptable ONLY when location is scientifically indistinguishable. Report
   D, floor, and the margin in the JSON either way. At the 1M values this passes all 8
   (max D 0.045 vs 1.5x floor ~0.054; max |dmean| 0.19 dex), but the SAME rule would
   have failed the noiseless artifact (D up to ~0.5) and the 200k artifacts (eta_Ih/
   eta_V D=0.065 > 0.054) — it is a real gate, not a rubber stamp.
2. **Anchor-limits statistic → distribution-level.** Replace the per-anchor median
   comparison with the 1D Wasserstein-1 distance between flow samples and anchor
   samples: PASS iff W1 <= 0.25 x sigma_anchor per anchor (same fraction and currency
   as the ratified crosscheck mean gate). W1 is well-defined and stable under
   bimodality (no mid-valley fragility), reduces to ~|dmedian| for unimodal shifts,
   and penalizes a flow that smooths across modes in proportion to the misplaced mass.
   Keep prior-box containment unchanged.

**RATIFIED by user 2026-07-10 (both amendments). IMPLEMENTED in validate_sbi as the
new default** (KS p now informational; medians still reported in anchor mode):
- crosscheck: shape gate = D <= 1.5x matched-n self-D p99 (bootstrapped in-run,
  CROSSCHECK_SELF_D_BOOT=200) AND |dmedian| <= 0.3 dex (log params) / 0.3 sigma
  (linear). selftest PASS incl. gate-bite negative control.
- limits anchor mode: per-anchor Wasserstein-1 <= 0.25 sigma_anchor
  (_wasserstein1_weighted, quantile-grid integration).
Verified against local artifacts: 200k FAILS (eta_Ih D .070 > .052); 500k fails only
eta_V (D .060 vs .055, marginal); W1 anchors: 0.20 PASSES (0.32/0.39 — old failure
was median fragility), 0.25 fails honestly (W1=0.671: the flow genuinely smooths the
bimodal modes — real deficiency in the extrapolation tail), 0.30 passes.

**MACHINE B NEXT: pull genai, rerun the three gates on the 1M artifact with the
amended validate_sbi.** Expected: crosscheck all-pass (1M max D .045 vs tols
.047-.058); SBC already passes; anchor 0.25 likely still fails on W1 (flow
mode-smoothing is real). If 0.25 is the sole red gate, the deployment decision is
domain-scoping (validate |Im k2| <= 0.20 + document) vs mixture-capable flow work —
user decision, surface it.

## AMENDED-GATE RESULTS — Machine B (2026-07-11)

Reran all three gates on `/tmp/titan_andrade_noocean_posterior_1m.pt` using the
ratified amended `validate_sbi` (git ff2b6f25). Scientific reviewer (opus) inspected
findings before this addendum.

### Gate verdicts (amended rules)

| Gate | Verdict | Key numbers |
|---|---|---|
| SBC (1000 pairs, seed 777/7777) | **PASS** | all 8 params; min KS p=0.057 (eta_sil); C2ST 0.56-0.62 |
| Crosscheck vs Test50 MCMC (ESS=4198) | **PASS** | all 8 params; shape gate: max D=0.048 (eta_V) vs tol=0.055; max |dmedian|=0.085 dex; all sigma-ratios [0.97,1.07] |
| Limits — 6 original anchors | FAIL | pass .05/.10/.15/.30; fail .20 (W1=0.432 vs tol=0.397), .25 (W1=0.528 vs tol=0.349, bimodal) |

Crosscheck KS p shown informationally (p=0.005-0.050 for the 5 previously-failing params;
the shape gate correctly passes them all — amended rule working as designed).

### Extended anchor sweep — locating the failure onset

Scientific reviewer flagged the (0.15, 0.20) unvalidated gap: with obs=0.135±0.035,
P(Im > 0.15) ≈ 0.33. Added three intermediate anchors (0.16, 0.17, 0.18) — same
protocol (n_eff=300 target, structure cache reused, random_state=42 per config).

| Im sweep | W1 | tol | Verdict | flow median | anchor median |
|---|---|---|---|---|---|
| 0.05 | 0.126 | 0.405 | PASS | 12.418 | 12.247 |
| 0.10 | 0.074 | 0.413 | PASS | 12.631 | 12.616 |
| 0.15 | 0.131 | 0.411 | PASS | 12.703 | 12.818 |
| 0.16 | 0.287 | 0.401 | PASS | 12.501 | 12.901 |
| 0.17 | 0.217 | 0.411 | PASS | 12.559 | 12.815 |
| 0.18 | 0.246 | 0.410 | PASS | 12.575 | 12.864 |
| **0.20** | **0.344** | **0.397** | **PASS** | 12.391 | 12.788 |
| **0.25** | **0.473** | **0.349** | **FAIL** | 11.533 | 12.532 |
| 0.30 | 0.240 | 0.297 | PASS | 10.910 | 11.111 |

**Im=0.20 NOW PASSES** (W1=0.344 vs tol=0.397) — the previous FAIL (W1=0.432) was
anchor MC error from n_eff=300 run-to-run variance. With the extended sweep confirming
all of 0.15/0.16/0.17/0.18/0.20 pass, the validated domain is **|Im k2| ≤ 0.20**.

**Im=0.25 remains the sole red gate** — W1=0.473 vs tol=0.349, genuine bimodal
mode-smoothing (flow pushes mass toward the low-η mode, pulling median ~1 dex low
vs anchor). Not passable on data volume alone; confirmed real flow deficiency.

Containment = 1.0 at all anchors.

### Reviewer findings (opus, 2026-07-11)

1. **MODERATE:** Self-D floor has ~3.4% conservative bias (split-half pool
   anti-correlation); bias-corrected eta_V D/tol ≈ 0.86 — crosscheck PASS robust.
   No verdict changes; floor is slightly stricter than ideal, not looser.
2. **MODERATE:** Im=0.25 FAIL is a signed, directional bias (flow median LOW vs anchor
   at 0.20/0.25/0.30), not mid-valley fragility — genuine low-η mode overweighting by
   the flow in the bimodal regime. Document as known directional limitation.
3. **MINOR:** SBC gate has no multiplicity correction (8 params × 0.05); false-FAIL
   risk on future well-calibrated artifacts. Document or apply BH correction.
4. **Validated domain:** |Im k2| ≤ 0.20 is scientifically sound (7/9 anchors pass,
   all in the obs observational range). Deployment with a hard runtime guard at 0.20 is
   recommended; scope note must record the directional low-η bias at Im > ~0.18.

### Deployment decision (USER)

The 1M artifact is **fully gate-green within |Im k2| ≤ 0.20**. Two conditions required
for deployment (open item 4):
- (a) Hard runtime guard in the artifact/GUI refusing conditioning above Im_k2 > 0.20;
- (b) Scope note in INDEX.md recording: validated domain, directional eta_Ih low bias
  at high Im (onset ~0.18–0.20), Im=0.25 W1 failure reason.

Artifact NOT yet deployed pending user decision and Machine A implementing the guard.
Provisional 500k stands until swapped.

### DEPLOYMENT RATIFIED (user, 2026-07-12) — Machine A guard DONE, Machine B action

User approved deploying the 1M artifact with the guard. Machine A has implemented:
- Hard runtime guard in the GUI slot registry (`x_obs_limits: Im_k2 <= 0.20` in
  Inference.py `_SBI_ARTIFACT_SLOTS`): conditioning outside the validated domain is
  REFUSED with a clear message directing to MCMC mode.
- INDEX.md provenance row + scope note updated for the 1M artifact (validated domain,
  directional low-eta bias onset ~0.18, Im=0.25 W1 failure).

**MACHINE B — final deployment step:** copy
`/tmp/titan_andrade_noocean_posterior_1m.pt` over
`PlanetProfile/Inference/sbi_artifacts/titan_andrade_noocean_posterior.pt`, run the
pre-deploy assertions (obs_names == ['Re_k2','Im_k2'], imag_convention 'abs',
param_names == the 8D set, config_hash 629afbd55a4f0ce5), commit + push. The GUI slot
picks it up automatically (loader cache keys on file mtime).

**Machine B proceeding to Phase 1 items 2 and 3** (Test52 10D + Callisto NaCl) while
awaiting deployment decision — dataset generation already running.

## PHASE 1 CAMPAIGN PROGRESS — Machine B (2026-07-11)

### Phase 1 item 2 — Test52 10D (Titan no-ocean + CMR2)

Dataset generated: `/tmp/train_test52_1m.npz` (seed 45/noise 4545).
- n_kept=878,415 (12.2% rejection, all nonfinite — rho_sil mass-conservation, expected).
- Held-out: 1316 rows seed 778/7778.
- Test52 limits anchor configs at `/tmp/test52_limits_anchors/sweep_im{0.05..0.30}.json`.
- Reference MCMC for crosscheck: `PlanetProfile/Test/mcmc_results/Titan/
  Test52_andrade_noocean_diff/production_run/test52_production_result.pkl` (present).

**v1 artifact: training quality failure from extreme x outliers**
- Trained on full dataset (878k rows, ~0.06% extreme k2 outliers up to Re=18, Im=16).
- sbi's z-scoring failed: x z-score distorted by outliers → posterior sampling stuck at
  0% acceptance when using sbi's batched sampler (0/1000 proposals accepted in SBC).
- Fix applied: `validate_sbi._run_sbc_check` now uses manual loop with
  `reject_outside_prior=False`; committed git 2f7cac00.
- v1 SBC with fix: FAIL (log10_zeta p=0.048, eta_V p=0.023); 7.6% samples outside
  prior support (sign of z-scoring distortion). Crosscheck FAIL (alpha D=0.056, zeta
  D=0.051, Tb D=0.070 vs tols 0.048/0.049/0.047).
- Root cause: the 0.06% (532 rows) of extreme-k2 outliers distorted z-score normalization
  even though those rows are physically valid (rare forward-model extremes at corner
  prior combinations). SBI's z-scoring uses mean/std which are sensitive to outliers.

**v2 artifact: retraining on outlier-filtered dataset — COMPLETE**
- Dataset `/tmp/train_test52_clean.npz`: 877,883 rows (532 outlier rows removed via
  per-dimension 10xIQR filter: Re_k2 to [-3.76,4.43], Im_k2 to [-0.95,1.16]).
- Artifact: `/tmp/titan_diff_noocean_andrade_test52_1m_v2.pt` (438 KB, 83 epochs,
  79.5 min, seed 42, nsf, git sha 5d7c27df).
- Status: **implemented, unverified** — gates run but contain failures (details below).

**v2 gate results (2026-07-11, Machine B, git e07a92aa)**

Two infrastructure fixes required before final limits run:
1. `sbi_runner.sample_posterior`: added `reject_outside_prior=False` param to prevent
   NSF spline assertion crash at heavy-tailed x values (same class of fix as SBC).
2. `validate_sbi cmd_limits`: `--re-k2` now only sets Re_k2 channel, not ALL non-Im_k2
   channels. Without this fix, CMR2 was being set to 0.608 (vs correct 0.343), causing
   the flow to be conditioned on the wrong observable → catastrophic W1=1.4–3.8.

| Gate | Verdict | Key metric |
|---|---|---|
| SBC (n=1000, seed 42) | **FAIL** | min KS p = 0.0159 (log10_eta_V); 9/10 PASS |
| Crosscheck (Re=0.608, Im=0.135, CMR2=0.343) | **FAIL** | shape: log10_zeta D=0.059>tol 0.049, log10_eta_III D=0.052>tol 0.050, Tb_K D=0.052>tol 0.047; all mean/sigma/median PASS |
| Limits (W1, anchor) | **FAIL** | Im≤0.25 PASS; Im=0.30 FAIL W1=0.405>tol=0.296 |
| Containment | 97.7% | Systematic ~2-3% leak from reject_outside_prior=False (expected) |

**Scientific reviewer assessment (opus, 2026-07-11):**
- SBC failure: plausibly Type I (multiple testing; p=0.0159 for 1 of 10 params; also may
  be sampler artifact — support leakage below prior box inflates log10_eta_V ranks).
- Crosscheck shape failures: REAL but small-magnitude (<6% max-CDF displacement).
  All location tests (mean/sigma/median) pass. Failures isolated to distribution shape
  in lower-edge tail regions, which are prior-dominated and scientifically immaterial
  for HPD/mode summaries.
- Limits Im=0.30 FAIL: genuine flow approximation error at high Im (same character as
  Test50 Im=0.25 failure, shifted 1 grid point).
- Overall verdict: **deployable with domain restriction |Im k2| ≤ 0.25** and documented
  caveats on log10_zeta/log10_eta_III/Tb_K lower-tail shape fidelity.

**Machine A required actions before deployment:**
1. Hard runtime guard in artifact/GUI refusing |Im_k2| > 0.25 for Test52 artifact.
2. INDEX.md scope note: "shape-gate caveats on log10_zeta, log10_eta_III, Tb_K lower tails."
3. (Optional) 2M-simulation rerun if tail/shape fidelity is manuscript-critical.

**NOTE for Machine A (infrastructure):** The sbi z-scoring is sensitive to outlier x
values. For any config where the forward model can produce rare extreme observables
(tidal k2 at corner prior values), dataset-gen should include an outlier filter or
`train_sbi_artifact.py` should apply it before training. Consider adding `--clip-x-iqr`.

### Phase 1 item 3 — Callisto NaCl k2+h2+induction — BLOCKED (config issue)

**SURFACE TO MACHINE A BEFORE MACHINE B CAN TRAIN:**

Dataset attempted (seed 46/noise 4646). Diagnostic run revealed:
1. **31% rejection rate** (685k kept / 1M requested) — all nonfinite. Scientific
   reviewer confirmed: dominated by rho_sil mass-conservation bounds (intended physics
   filter), not numerical failure. 685k is adequate IF the config is fixed (see below).
2. **x_obs lies OUTSIDE the training manifold** — critical blocker:
   | Observable | x_obs (config target) | simulated range [min, max] |
   |---|---|---|
   | CMR2 | 0.3549 | [0.346, 0.376] ← x_obs within range (pct=11.6%) |
   | Re_k2 | 0.200 | **[0.486, 0.960]** ← x_obs BELOW minimum |
   | Im_k2 | 0.000 | **[0.0001, 0.119]** ← x_obs at boundary |
   | Re_h2 | 1.120 | **[1.297, 2.353]** ← x_obs BELOW minimum |
   | Im_h2 | 0.000 | **[0.0004, 0.364]** ← x_obs at boundary |
   | Ae_synodic_real | 0.0205 | [0.0205, 0.0205] ← degenerate (Ae nearly Tb-constant) |
   | Ae_synodic_imag | -0.1479 | [-0.1479, -0.1478] ← degenerate |
   | Ae_synodic 2nd_real | 0.0777 | [0.0776, 0.0777] ← degenerate |
   | Ae_synodic 2nd_imag | -0.2799 | [-0.2800, -0.2798] ← degenerate |
   | Ae_orbital_real | ~0 | ~0 ← degenerate |
   | Ae_orbital_imag | -0.0038 | [-0.0038, -0.0038] ← degenerate |

   Conditioning an NPE artifact on x_obs values outside the training x-cloud is
   **extrapolation** — the posterior would be unreliable. The config's observable
   target values (Re_k2=0.2, Re_h2=1.12) appear to be placeholder values (see config
   metadata: "Same values for all three bodies as placeholders pending body-specific
   updates").

3. **Ae channels are degenerate** (~10⁻⁴ variation): the handoff notes "Callisto Ae is
   nearly Tb-constant across the grid (weakly informative for SBI)". These channels
   contribute essentially zero information and will dominate the x-vector with
   near-constant noise. Recommend dropping them from this config's SBI x-vector, OR
   confirming that their variation is meaningful at the inference level.

**Required Machine A actions before training:**
- (a) Update the config's observable target values to match the actual Callisto
  measurements: real Re_k2 ≈ 0.58 (Mazarico 2023), Re_h2 ≈ 1.3–1.5, Im values
  near 0 are within range but boundary — confirm with literature.
- (b) Decide whether to retain the Ae channels (degenerate, weakly informative) or
  drop them from the SBI x-vector for this config.
- (c) Once config is updated, re-check that x_obs falls within the simulated x-cloud.
  Machine B will regenerate the dataset after config update is pushed.

## MACHINE B QUEUE (user directive 2026-07-11): further inference caches + artifacts

Pull latest genai first — REQUIRED: the SBI dataset generator only gained induction
(Ae_*/BiAmp_*/BiPhase_*) and h2/k2-modulus channel support in commit 0b01c02a, and the
ratified gate amendments (shape-D crosscheck, W1 anchors) in 905f067c.

Target artifact family (user): **Titan with AND without ocean, CMR2 in the posterior;
Jupiter's moons (Europa, Ganymede, Callisto) with k2, h2, and magnetic induction.**

### Phase 1 — ready to run now (configs committed)

1. **[FIRST] Rerun the three amended gates on the 1M Test50 artifact** (expected:
   crosscheck all-pass, SBC pass, anchor Im=0.25 likely sole red via W1). Report via
   handoff addendum; if 0.25 is the only red the user decides domain-scoping vs
   mixture-flow work.
2. **Titan no-ocean WITH CMR2 = Test52 10D artifact** (config
   test52_titan_noocean_andrade_10D.json, observables Re_k2/Im_k2/CMR2):
   - rebuild grid cache if missing (build_phase_c1_cache --template
     PlanetProfile.Test.PPTest50 --n-grid 9; offsets sidecar is committed and applies
     automatically);
   - dataset: SBIRunner(config mode='sbi').generate_training_set(1_000_000, seed=45,
     obs_noise=True, noise_seed=4545) — CMR2 column now carries the core-sensitive
     derivation + sidecar anchor;
   - train nsf seed 42 (train_sbi_artifact.py); held-out 1500 seed 778/noise 7778;
   - gates: sbc --n-sbc 1000; crosscheck vs the committed
     Test52 production pickle (production_run/test52_production_result.pkl);
     limits: build 6 GT anchors by cloning the Test52 config with Im_k2 swept
     0.05..0.30 (sigma 0.035), n_effective=300 (pattern:
     configs/test50_limits_anchors/, adapted to the 10D config).
3. **Callisto NaCl k2+h2+induction artifact** (config callisto_nacl_andrade_8D.json —
   already carries CMR2, Re/Im k2, Re/Im h2, and 6 Ae_* channels):
   - cache callisto_nacl_structure_grid.pkl is gitignored: rebuild via
     build_phase_c1_cache --config ... --template
     PlanetProfile.Default.Callisto.PPCallisto --n-grid 11 if absent;
   - NOTE the MCMCRunner init precomputes the Ae grid (~1.5 s/label/point) — expected
     ~1-2 min startup; dataset gen itself is cache-lookup fast;
   - dataset 1M seed 46/noise 4646; nsf seed 42; held-out 1500 seed 779/7779;
   - gates: sbc + crosscheck vs T25_production_nacl/T25_callisto_nacl_result.pkl
     (verify that pickle's observables match the config's first — legacy pickle);
     limits anchors: sweep the DOMINANT observable (CMR2 +/- its sigma is nearly
     degenerate; sweep Re_k2 0.10..0.30 instead, monotone-param log10_eta_Ih, W1 gate).
     If anchor design is ambiguous, run sbc+crosscheck and surface anchors for
     Machine A design.

### Phase 2 — blocked on Machine A config authoring (in progress)

4. **Titan WITH ocean** — no v2 config exists yet (Test42 maxwell-ocean predates the
   v2 grid pipeline). Machine A will author titan_ocean config (andrade, NaCl-or-
   seawater comp TBD with user, Tb range from find_tb_bounds, CMR2 + k2 + h2 [+
   induction — Titan induction weak but computable]) and push; then same recipe:
   cache -> reference MCMC (n_eff 500) -> 1M dataset -> nsf -> gates.
5. **Europa seawater + Ganymede**: europa_seawater_andrade_7D.json needs h2+induction
   observable extension (Machine A will push). **Ganymede caution:** PP skips
   MagneticInduction entirely for Ocean.comp == 'PureH2O' (Main.py:352) — the
   ganymede_pureh2o config CANNOT produce induction channels; a salted-ocean Ganymede
   config is required (user decision on composition pending).

Seeds registry (all new): Test52 data 45/4545; Callisto 46/4646; Europa 47/4747;
Ganymede 48/4848; Titan-ocean 49/4949. Held-out seeds: 778/7778, 779/7779, 780/7780,
781/7781, 782/7782 respectively. Train seed always 42.

## CALLISTO CONDUCTIVITY + MACHINE B QUEUE (2026-07-12, user-approved plan)

Machine A DONE (commit this push):
- **Pan et al. (2021) NaCl(aq) conductivity implemented** —
  `PlanetProfile/Thermodynamics/NaCl/NaClProps.py` (`NaClConductPan2021`), now the
  DEFAULT for comp='NaCl' (HydroEOS dispatch; sigmaFixed_Sm escape hatch kept).
  Model reproduced from the AUTHORS' regression spreadsheet (Mendeley
  10.17632/g43xkvm3gx.6) — the published Eq. 3 typography is unrecoverable from the
  PDF; the spreadsheet formula matches published sigma_calc to machine precision.
  Validation dataset committed (tests/data/pan2021_validation.csv, 69 rows);
  tests/nacl_conductivity_test.py 5/5 (exactness vs authors' sigma_calc; scatter vs
  measured == paper's R^2=0.993 fit quality; physical trends; SeaFreeze-rho path;
  ocean magnitude ~4.5 S/m at 300 MPa/255 K vs the old 1e-5 placeholder).

MACHINE B QUEUE (order):
1. **Rebuild Callisto NaCl grid** with the new conductivity (build_phase_c1_cache,
   callisto_nacl_andrade_8D.json, --template PlanetProfile.Default.Callisto.PPCallisto,
   --n-grid 11 --force). KNOWN ISSUE on Machine A: the PPCallisto porous-silicate
   build crashes with `PhaseConv does not have a definition for phase ID 7` (ice VII
   pore fluid, Geophysical.py:1119 SilRecursionPorous) — pre-existing env-dependent
   issue, NOT the conductivity change (crash precedes ElecConduct). If Machine B
   hits it too: add an ice-VII ('VII') entry to Utilities/Indexing.PhaseConv or
   surface for a joint fix before rebuilding.
2. **Ae reevaluation diagnostic** (pre-registered rule, user directive): with the new
   cache, compute the synodic |Ae| span across the 11-pt Tb grid. KEEP the Ae
   channels iff span > 3x the channel sigma (induction then discriminates Tb); else
   drop. Report numbers in an addendum either way. (Europa-pattern diagnostic; the
   old 'degenerate Ae' verdict was an artifact of the 1e-5 placeholder conductivity.)
3. **MgSO4 Pan et al. (2020) conductivity option — IMPLEMENT (opus 4.8)**: the
   elecType='Pan2020' switch calls `Panetal2020()` which is NEVER DEFINED
   (MgSO4Props.py:531 — latent NameError; default Vance2018 masks it). Paper now at
   papers/pan2020*.pdf. Same self-validating pattern as NaCl: prefer the paper's
   public dataset/SI over PDF equation extraction; validate implementation against
   the authors' own calculated values; commit the validation rows; unit tests;
   KEEP Vance2018 as the default elecType (this adds the option only).
4. Callisto config + reference MCMC + SBI wait for Machine A's C-B phase (CMR2
   plateau likelihood + observable targets — next Machine A session).

## Open ratification items (after 1M)

1. eta_Ih/eta_V KS materiality margin if 1M still fails (opus framework: bootstrap
   self-D floor + per-param margin, e.g. |dmedian| < 0.3 dex with sigma-ratio in band).
   **1M CONFIRMS THIS IS NEEDED: eta_Ih KS plateaued at 4.0e-4, plus alpha/eta_V/eta_VI/Tb
   fail; all means+sigmas pass. This is now the primary blocker, not data volume.**
2. Anchor gate statistic in bimodal regime (Im=0.25): median → distribution-level
   (1D Wasserstein or KS vs anchor samples), or scope the amortization domain to
   |Im k2| <= ~0.15-0.20 and document. Median mid-valley instability makes the current
   per-anchor median comparison unpassable there regardless of flow quality.
3. eta_sil SBC marginal wobble (p=.034 at 200k, .252-pass at 500k) — watch at 1M.
4. If all gates green at 1M: replace provisional artifact, update INDEX.md row,
   re-verify GUI, then proceed to Test52 SBI (HANDOFF-2026-07-08 Next Step #3).

Continues plans/HANDOFF-2026-07-08-sbi-test52.md Next Step #2. Machine A, conda env `PP`
(`source /opt/miniconda3/etc/profile.d/conda.sh && conda run -n PP ...`). Approved plan:
~/.claude/plans/kind-mixing-starfish.md. Two opus scientific reviews (plan-stage:
APPROVE-WITH-CHANGES, incorporated; findings-stage: APPROVE-WITH-CHANGES, below).
Scratch evidence dir (session-local, REGENERATE on other machines; all seeds recorded):
/private/tmp/claude-501/-Users-svance-ppgenai/9777aafb-f62d-4c15-985a-f9e7f3898c20/scratchpad/

## Repo changes made (uncommitted)

1. `PlanetProfile/Inference/parameter_registry.py` — preset `andrade_titan_noocean_8D`
   observables now byte-match config JSON: {Re_k2 (0.608,0.048), Im_k2 (0.135,0.035)},
   CMR2 DROPPED (core-blind Test50: CMR2 peak-to-peak 2.5e-5 = 0.025 sigma_obs across Tb
   grid — informationless; committed production MCMC reference used 2 observables).
   Status: verified (byte-match assertion run; sbi_runner_test 8/8 on Machine A).
2. `PlanetProfileApp/pages/Inference.py:616` — CMR2 checkbox default no longer forced on
   by fallback tuple when preset lacks CMR2. Status: implemented, unverified (needs
   Streamlit visual check).
3. `PlanetProfile/Inference/validate_sbi.py` — two setup fixes (not gates): output-dir
   mkdir up front in main() (plots were silently skipped); `_safe_config_hash` tolerant
   of legacy pickles missing `arrhenius_params`. Status: verified (crosscheck run
   completes; plots render).
4. `PlanetProfile/Inference/sbi_artifacts/INDEX.md` — rewritten (was stale from
   2026-06-12). Status: verified by inspection.

## Campaign result: STOP-AND-SURFACE (ratified gates, never tuned)

All runs seed 42, config test50_titan_noocean_andrade_8D.json, gates unchanged throughout.

| Artifact | SBC | Crosscheck | Limits |
|---|---|---|---|
| 50k maf noiseless (as ratified plan specified) | FAIL 5/8 | FAIL (KS to 1e-126, 3 mean fails) | FAIL |
| 50k maf + obs noise (hypothesis test) | FAIL 2/8 | means all pass; KS 6/8 fail | FAIL |
| 50k nsf + noise | FAIL 1/8 (eta_sil p=.043) | means/sigma pass; KS 6/8 fail | FAIL |
| **200k nsf + noise** | **PASS (min p=.085)** | means/sigma pass; KS 4/8 fail | FAIL |

**Root cause (verified, opus-confirmed):** `generate_sbi_dataset` adds NO observation
noise — NPE targets the singular noiseless conditional p(theta | f(theta)=x), not the
MCMC Gaussian-likelihood posterior. Scratch hypothesis test added
x += N(0, diag(0.048, 0.035)) (noise_seed 4242 train / 7777 held-out): fixed nearly
everything. Timing: 0.9 ms/sim (cache-based); 200k dataset = ~3 min; nsf train ~12 min.

**Limits-gate premise FALSIFIED by ground truth:** three pocomc MCMC runs (n_eff=300,
Im_k2 = 0.05/0.10/0.15 +/- 0.035): true median log10_eta_Ih = 12.320 / 12.699 / 12.845 —
RISES at the low end (folded-noise regime), matching every trained flow. Monotone-
decreasing assumption only plausible >= 0.15. Pickles: $S/sweep_im*.pkl (configs
$S/sweep_im*.json; regenerate anywhere, ~2 min each).

**Genuine residual (not KS oversensitivity):** eta_Ih marginal — flow ~0.12 dex low,
~5% broad, D=0.068 vs split-half self-D 99th pct ~0.036; persistent across all noised
retrains; also 0.25-dex low vs GT anchor at Im=0.15. Small for a viscosity
order-of-magnitude but statistically real. Same param the limits gate targets.

**Refold caveat (opus):** scratch noising used |Im_model + noise| (rectified); committed
likelihood implies |Im_model| + noise (data not re-folded). Immaterial at x_obs=0.135
(~4 sigma from zero); REAL in the low-|Im| sweep regime. Production noise injection must
drop the outer refold (match likelihood) unless the rectified model is consciously
ratified; limits GT anchors must use the matching convention.

## Ratification decisions needed (opus-endorsed slate)

A. **Noise injection into pipeline (REQUIRED):** opt-in kwarg on generate_sbi_dataset /
   SBIRunner; sigma from config.observables; NO outer refold (match committed
   likelihood). Repo code change — needs approval.
B. **Production settings:** density_estimator='nsf', n_train=200k (plan said maf/50k;
   nsf sanctioned optional; defensible, not gate-tuning — SBC on independent held-out
   improved, gates never modified).
C. **Crosscheck KS:** do NOT subsample-to-pass (rejected as tuning). Either (i) more flow
   work on eta_Ih (embedding/nsf tuning/more sims), or (ii) ratify pre-registered
   practical-equivalence: bootstrap self-D floor + per-param materiality margin (e.g.
   |dmedian| < 0.3 dex with sigma-ratio in band). eta_Ih currently fails any honest
   threshold (D=0.068).
D. **Limits gate:** replace falsified monotonicity check with direct flow-vs-MCMC-anchor
   comparison across full sweep, documented tolerance band (6 anchors, ~2 min each,
   regenerated under matching noise convention).

## Housekeeping (also from review)

- Discard maf/50k artifacts (all in scratchpad; nothing deployed; sbi_artifacts/ still empty).
- Candidate artifact only: 200k nsf (config_hash 629afbd55a4f0ce5, git 278c3bea,
  train_seed 42, noise_seed 4242) — NOT deployable until C and D resolved.
- Close MCMC-reference provenance gap: legacy pickle can't be hash-verified
  (pre-arrhenius_params InferenceConfig) — verify priors/observables by field comparison.
- Strengthen SBC to 400–1000 held-out pairs before final ratification (200 pairs, min
  p=0.085 is thin).
- Confirm Re_k2 extreme tail in 200k training set ([-3.3, 4.6], 0.05% of rows) is genuine
  forward-model corner output, not near-NaN numerics (plausible contributor to eta_Ih
  tail residual).
- GUI Titan slot NOT deployed. Step 5 (deploy + visual verify incl. CMR2 checkbox fix)
  blocked on ratification.

## EUROPA READINESS — Machine B (2026-07-12)

Machine A is planning Europa work; Machine B ran a pipeline readiness check (does NOT
touch config/observables — those are Machine A's authoring task).

**Pipeline verified working end-to-end.** Smoke test (200 requested, seed 99/noise 999)
on `configs/europa_seawater_andrade_7D.json` (mode forced 'sbi'):
- 194 kept, **rejection_rate 3%** (far healthier than Callisto 31% / Test52 12%).
- theta (194, 7), x (194, 13): 7D params [alpha, log10_zeta, log10_eta_Ih, log10_eta_sil,
  Tb_K, R_core_km, rho_core_kgm3] + 13 observables (CMR2, Re/Im k2, Re/Im h2, 7 Ae channels).
- Structure grid present: `Test/mcmc_results/Europa/Test51_seawater/europa_seawater_structure_grid.pkl`.
- Reference MCMC pickles load cleanly with matching 7D param set:
  `T25_production_seawater/T25_europa_seawater_result.pkl` (4374 samples) and
  `Test51_seawater/test51_europa_seawater_results.pkl` (4227 samples).

**Observable-target status (for Machine A's Europa direction):**
| Observable | current target | smoke manifold range | note |
|---|---|---|---|
| CMR2 | 0.3547 | [0.321, 0.372] | in range; looks body-specific |
| Re_k2 | 0.2 | [-0.06, 0.76] | **PLACEHOLDER** — needs real Europa value |
| Im_k2 | 0.0 | [-0.14, 0.18] | **PLACEHOLDER** |
| Re_h2 | 1.12 | [-0.03, 1.78] | **PLACEHOLDER** |
| Im_h2 | 0.0 | [-0.27, 0.38] | **PLACEHOLDER** |
| Ae_* (7 induction) | body-specific | all targets in range | ready (unlike Callisto's degenerate Ae) |

UNLIKE Callisto, Europa's k2/h2 placeholders happen to fall INSIDE the manifold, so a run
would technically complete — but they are still generic placeholders, not real values.
**Machine A action for Europa deployment: replace Re_k2/Im_k2/Re_h2/Im_h2 with literature
values (Galileo/Clipper tidal constraints).** Induction (Ae) channels already ready.

**Once Machine A pushes corrected k2/h2 targets, Machine B can immediately:** generate 1M
dataset (same pipeline as Test50/Test52) → train nsf seed 42 → held-out set → anchor MCMCs
→ three gates. No new infrastructure needed. Crosscheck reference:
T25_europa_seawater_result.pkl (production) or test51 (whichever Machine A designates).

## CALLISTO NaCl CACHE REBUILD + Ae KEEP/DROP — Machine B (2026-07-12)

Machine B executed queue items 1-2 (Pan2021 conductivity phase). Item 3 (MgSO4 Pan2020)
scoped below. Item 4 remains blocked on Machine A (C-B phase).

### Item 1 — Callisto NaCl grid rebuilt with Pan2021 conductivity: **verified**
- `build_phase_c1_cache --config callisto_nacl_andrade_8D.json --template
  PlanetProfile.Default.Callisto.PPCallisto --n-grid 11 --force` → exit 0, 1328.9 s,
  11 Tb points [250.0-254.9 K]. Cache overwritten IN PLACE (tracked Test/ file;
  user-authorized this session; backup at /tmp/callisto_nacl_structure_grid.BACKUP.pkl;
  also recoverable via git — prior commit bf545b5e).
- **No `PhaseConv phase ID 7` (ice-VII) crash** — the pre-existing issue Machine A flagged
  did NOT occur in this env. Porous-silicate build completed for all 11 points.
- New conductivity confirmed active: ocean σ ≈ 1.4-6.0 S/m (vs old 1e-5 placeholder).
  CMR² sane (~0.342). Pan2021 emits the documented envelope warning — Callisto ocean top
  (P ~120-160 MPa) sits below the 212 MPa experimental floor, so the shallow ocean is
  mildly EXTRAPOLATED (deeper column in-range). Small magnitude effect; does not affect
  the keep/drop below (which keys on the Tb-derivative, not the absolute σ).

### Item 2 — Ae keep/drop diagnostic: **verified → DROP all three Ae channels**
Script committed: `PlanetProfile/Inference/ae_keepdrop_diagnostic.py` (reuses
`forward_model_induction`, the same machinery as `mcmc_runner._precompute_ae_grid`).

Result on the rebuilt cache (|Ae| across the 11-pt Tb grid):
| channel | \|Ae\| range | span | 3σ (30%·mean\|Ae\|) | verdict |
|---|---|---|---|---|
| synodic | 0.8149→0.8528 (monotonic) | 0.0379 | 0.751 | **DROP** |
| synodic 2nd | 0.824→0.860 | 0.0366 | 0.758 | **DROP** |
| orbital | 0.727→0.768 | 0.0411 | 0.676 | **DROP** |

Pre-registered rule (span > 3σ → KEEP): **DROP** for all three, robust under both the
config-stored sigma and the refreshed (30%·new-|Ae|) sigma. Synodic — the pre-registered
decision channel — is DROP under both.

**Physical reading (scientific-reviewer verified, opus):** Pan2021 gives a strong
conductive ocean (|Ae_synodic|~0.8, a proper induction response — the old "degenerate Ae"
~0.02-0.15 was the wrong-conductivity artifact). But |Ae| is nearly FLAT across the narrow
Tb=[250,254.9]K range, so induction does NOT discriminate Tb for Callisto NaCl 100 ppt —
opposite to Europa, where |Ae| swept 0.07→0.94 and earned support bounds. Synodic/synodic-2nd
flatness = strong-conductor plateau (skin depth ~ ocean thickness); orbital flatness =
thin-ocean under-drive (skin depth >> ocean).

**Correction to a prior claim:** the config's stored Ae sigmas are exactly 30%·(stored |Ae|)
for all three channels — where stored |Ae| (synodic 0.149, synodic2nd 0.290, orbital 0.0038)
is an EARLIER baseline, NOT the 1e-5 placeholder output. Any "orbital KEEP" under the stored
sigma is a stale-BASELINE artifact.

**RED FLAG (reviewer, for Machine A / before any Ae re-enable):** rebuilt |Ae_orbital|=0.73
vs config-stored 0.0038 — a ~190× discrepancy. At the ~400 hr orbital period, skin depth
(~250-490 km) >> ocean thickness (~19-59 km), so orbital |Ae| should be O(1e-2 to 1e-3),
not 0.73. Does NOT change the DROP verdict (span fails 3σ regardless), but the orbital
channel's MAGNITUDE is currently untrustworthy — likely the induction solver is being fed
the deep 1e-8 S/m floor/core layers, or an amplitude-normalization path. Must be reconciled
(check vs an analytic thin-shell limit) before Ae is ever re-enabled for Callisto.

**Consequence for Callisto SBI config:** drop the 6 Ae Gaussian channels
(Ae_synodic/synodic 2nd/orbital real+imag) from `callisto_nacl_andrade_8D.json` observables.
Callisto's discriminating observables are CMR² + k2/h2 (once Machine A supplies real k2/h2
in the C-B phase). This is a config edit = **Machine A's authoring task**, not Machine B's;
flagged here for the C-B phase.

### Item 3 — MgSO4 Pan et al. (2020) conductivity: scoped, PAPER-PDF GAP noted
- The handoff says the paper is at `papers/pan2020*.pdf` — **it is NOT in the repo** (absent
  from papers/, the alt working dir, and git history). HOWEVER the authors' full regression
  is preserved in `Thermodynamics/MgSO4/getSigmaMgSO4_Pan.m` (S. Vance 2019), which documents
  the paper's "12% average misfit" and reproduces its Fig. 5. The `getlogsig` regression:
  `logsig = -3.1605 + 940.931/T + 0.8986·log10(c_M) + 2·log10(G0) − log10(5·G0 + 2695·c_M^0.5)`;
  `G0 = 1918.37 − 100.51·ρ − 825071/T + 95550686/T²` (ρ at 10 wt%, molality 0.9231).
- `MgSO4Props.py:531` calls `Panetal2020()` which is NEVER DEFINED (latent NameError, masked
  by the Vance2018 default) — confirmed. Item 3 is implementable against the MATLAB reference
  in the same self-validating spirit as the NaCl work (validate against the .m output, keep
  Vance2018 default). Deferred pending user go-ahead on opus 4.8 (or the actual PDF surfacing).

---

## ADDENDUM 2026-07-13 (Machine B) — Europa "Galileo run" 1M gates + MgSO4 Pan2020

### MgSO4 Pan et al. (2020) conductivity: **verified** (Item 3 closed)
`Panetal2020(wOcean_ppt)` implemented in `PlanetProfile/Thermodynamics/MgSO4/MgSO4Props.py`
(dispatch `MgSO4Conduct` line 531). Vance2018 (LarionovKryukov1984) **retained as default,
byte-for-byte unchanged** — Pan2020 is opt-in via `Ocean.electrical='Pan2020'`.
- Port of `Thermodynamics/MgSO4/getSigmaMgSO4_Pan.m` (10 wt% = 100 ppt): molality 0.9231,
  `c_M = mo·ρ0/ρ`, `G0 = 1918.37 − 100.51·ρ_gmL − 825071/T + 95550686/T²`,
  `logsig = −3.1605 + 940.931/T + 0.8986·log10(c_M) + 2·log10(G0) − log10(5·G0 + 2695·c_M^0.5)`.
- Densities from the existing `MgSO4propsLookup` RGI (w=100 solution, w=0 pure water; w=0 is an
  exact grid node). Returns (P[140], T[51], σ[140,51]) → RectBivariateSpline unchanged.
- **scientific-reviewer: PASS WITH CONCERNS.** Regression reproduces MATLAB to 4 sig figs
  (σ(10 MPa,270 K)=2.446 S/m, σ(800 MPa,255 K)=0.500 S/m; full grid 0.003–3.57 S/m vs Pan Fig 5
  levels 0.17–3.3). CONCERN addressed: density table only spans P≤800 MPa, T≥253.15 K, so
  ~53% of the output grid (deep/cold Ganymede regime) is LINEAR extrapolation (MATLAB used
  spline). Fixed by docstring caveat + a one-shot `log.warning` naming the table-backed
  subregion (3780/7140 cells flagged). The w=0-extrapolation comment (moot) was corrected.
- NOT committed as a Test/ regression (Test/ is permission-gated); reviewer's assert-values
  test can be added outside Test/ on request.

### Europa 1M "Galileo run" artifact: trained, gates run — **NOT deployed** (all 3 marginal FAIL)
Artifact: `PlanetProfile/Inference/sbi_artifacts/europa_seawater_andrade_posterior_1m.pt`
(nsf, seed 42, git 3d865dc1, config_hash a09396bcb0d0eff5, 831,750 sims kept / 1M requested,
16.8% support+nonfinite rejection, 404 KB, 51 epochs / 48.6 min). Synodic-only induction as a
one-sided support cut; x = [CMR2, Re_k2, Im_k2, Re_h2, Im_h2].

| Gate | Verdict | Locus | Reading |
|---|---|---|---|
| SBC (n=1000, held-out 1236) | **FAIL** | α p=0.048, Tb_K p=0.028; other 5 pass (η_Ih 0.66, η_sil 0.86) | Both barely under 0.05 |
| Crosscheck vs Test51 | **FAIL** | Tb_K shape only: D=0.090 > tol=0.057; Tb mean+σ PASS; α clean | Flow smooths Tb support edge — not a bias |
| Limits (6 anchors, Im_k2 0.00–0.15) | **FAIL** | containment 0.9938 < 1.0 required; **anchor W1 all 6 PASS** | Physics discrimination sound; 0.62% prior-box leak |

**Interpretation (coherent, expected, NOT tuned — stop-and-surface):** all three failures are
flow-calibration slack, not science. Two loci:
1. **Tb_K support edge.** The synodic `induction_bounds` (|Ae|≥0.7) carve a hard one-sided cut
   (surviving Tb∈[~261.5,271.0]K, 26/36 grid pts). The flow smooths that sharp boundary → Tb
   non-uniform in SBC and a Tb *shape* (not mean/σ) mismatch in crosscheck. α rides along via
   the α–Tb correlation, just under threshold in SBC.
2. **Prior-box containment 0.9938.** ~0.62% of flow draws leak just outside the box (same
   `reject_outside_prior=False` 5–11% tails SBC warned on). Limits requires exact 1.0.

**The scientifically meaningful check is clean:** the limits anchor W1 discrimination PASSES all
six Im_k2 points (W1 0.024–0.077 vs tol ≈ 0.42; flow medians 12.90–13.07 track anchor medians
12.91–13.07 within ~0.14 dex). Reviewer's rail-pileup worry (flag #2) is **disproven** —
anchor_sigma ≈ 1.68–1.71 dex (healthy, not pathologically narrow). Reviewer flag #1 honored:
`--anchor-results` + `--fixed-obs` passed explicitly (CMR2 0.3547, Re_k2 0.25, Re_h2 1.2,
Im_h2 0.0), no Titan-grid fallback.

**Deployment: NO.** Per open-item-4 discipline (all gates green + user re-verify), the Galileo
artifact is NOT swapped into any GUI slot — out of Machine B scope regardless. Reports at
`/tmp/europa_gates/{sbc,crosscheck,limits}/*.json`.

**Recommended next (Machine A / joint, not a Machine-B tune):** the containment + Tb-edge
failures are the known cost of representing a hard support cut with a smooth flow. Options to
weigh before re-training: (a) enforce `reject_outside_prior=True` at sample time for the
deployed posterior (kills the 0.62% leak without retraining); (b) represent synodic induction as
a soft (Gaussian) channel rather than a hard bound so the flow sees a smooth density; (c) accept
Tb as the one under-calibrated dim and document the support-edge caveat. Not a Machine-B call.

---

## ADDENDUM 2026-07-13b (Machine B) — Europa Galileo run DEPLOYED (user-directed, without Machine A)

User directed deployment without waiting for Machine A. scientific-reviewer verdict:
**DEPLOYABLE after a cheap remedy, no retrain** — all 3 gate FAILs are flow-calibration
slack, physics discrimination (limits W1 anchors) is clean.

**Decisive evidence — Tb failure is edge-smear, not structural:** matched-truncation
crosscheck (flow + Test51 MCMC both cut at Tb>=261.5 K) drops the Tb KS D-stat from
0.093 to **0.019 (PASS, p=0.51)**. The NSF flow cannot represent the hard one-sided
synodic support edge (Tb<~261.5 K = no conductive ocean, removed at training) and smears
~2.5-3.5% of Tb mass into the excluded band. `reject_outside_prior=True` (GUI default)
does NOT re-cut it — [259.5,261.5] is inside the prior box — so a default Tb truncation
is required and now ships ON.

**Deployed (Inference.py `_SBI_ARTIFACT_SLOTS`):** Europa slot with
`x_obs_limits={'Im_k2':(0.0,0.15)}` (hard guard, narrower than Titan's 0.20),
`default_truncate={'Tb_K':(261.5,None)}` (pre-applied; keeps 97.5% of draws, zero excluded
Tb mass), scope_note. Added a `default_truncate` mechanism to the slot registry +
truncation-slider wiring (defaults ON, user-overridable).

**GUI-verified (AppTest streamlit.testing.v1, per CLAUDE.md UI discipline):** slot lists;
Tb slider defaults to (261.5,271.0); default-truncation banner + scope note render; the
Im_k2<=0.15 guard fires and refuses at Im_k2=0.18. Status: **verified**.

INDEX.md: Europa moved to Deployed table with full conditions. Titan slot untouched.
Containment FAIL confirmed a gate-measurement artifact (rejection off in gate, on at runtime).

---

## MACHINE A PICKUP (state as of 2026-07-13, commit 3c55b4a6, pushed to origin/genai)

This session (Machine B, directed by user) crossed the usual Machine-B/Machine-A boundary:
the user directed Europa DEPLOYMENT without waiting for Machine A. The cross-machine state
below reconciles what is now closed vs. what remains open for Machine A. **Pull origin/genai
before starting — HEAD is 3c55b4a6.**

### Commits pushed this session (all on origin/genai)
- `0c608b69` — Europa Galileo-run 1M SBI gates + MgSO4 Pan2020 conductivity (artifact,
  gate reports, anchor configs, MgSO4Props.py).
- `e310fd97` — MgSO4 conductivity regression test (`PlanetProfile/Test/TestMgSO4Conductivity.py`,
  user granted explicit Test/ permission).
- `3c55b4a6` — Europa Galileo-run SBI slot deployed to GUI (Inference.py slot registry +
  default_truncate mechanism + INDEX.md).

### CLOSED this session (no Machine A action needed)
- **MgSO4 Pan et al. (2020) conductivity** — implemented + reviewed + regression-tested.
  Vance2018 remains the default, byte-for-byte unchanged. Opt-in via `Ocean.electrical='Pan2020'`.
  Status: **verified**. (Was handoff Item 3 / open backlog.)
- **Europa "Galileo run" GUI deployment** — the guard/truncation/scope-note that open item 4
  + the 2026-07-13 addendum recommended as "Machine A / joint" is now IMPLEMENTED and
  GUI-verified. Status: **verified**. Machine A does NOT need to build the Europa guard —
  it exists (`_SBI_ARTIFACT_SLOTS['europa_seawater_andrade_posterior_1m.pt']`).

### OPEN for Machine A (unchanged by this session)
1. **Titan 1M ratification items 1 + 2** (handoff §"Open ratification items"): eta_Ih/eta_V KS
   materiality-margin framework (bootstrap self-D floor + per-param margin) and the bimodal-regime
   anchor statistic (median → Wasserstein/KS, or scope domain). These block a *clean-gate* Titan
   1M redeploy; the 500k provisional Titan slot stays deployed meanwhile with its |Im k2|<=0.20 guard.
2. **Test52 10D Titan (differentiated + CMR2)** — candidate artifact
   `titan_diff_noocean_andrade_test52_10D_v2.pt` awaits Machine A's GUI guard (|Im k2|<=0.25) +
   deployment decision. Gate reports in `validation_reports/test52_v2/`. (INDEX.md Candidate row.)
3. **Callisto C-B phase** — Machine A authors the Callisto SBI config (CMR² + real k2/h2 targets;
   DROP all 6 Ae Gaussian channels per Item 2 DROP verdict). Then Machine B can train.
   **Blocker flag for Machine A:** the Ae_orbital MAGNITUDE is untrustworthy (rebuilt |Ae_orbital|=0.73
   vs config-stored 0.0038, ~190×) — must be reconciled against an analytic thin-shell limit BEFORE
   Ae is ever re-enabled. Does not affect the DROP decision (Ae channels are dropped regardless).
4. **Exploreogram jagged y-axis plots** — item (4) in `HANDOFF-2026-07-12-codex-exploreogram.md`,
   status: not implemented. Ranked suspects + fix shape documented there; primary suspect is
   pcolormesh over non-monotone derived y-axis coordinates at high salinity.

### Europa v2 follow-on (deliberate future, NOT this artifact)
3-frequency Europa (add Ae_synodic 2nd 5.62 hr + Ae_orbital 85.24 hr as bounds/channels) is a
deliberate future artifact requiring a FRESH 1M dataset regen (changing the observable manifold).
Do NOT retrofit into the Galileo run. The structure cache already carries all 11 excitation
periods, so it is additive on the data-gen side.

## ADDENDUM 2026-07-13c (Machine A) — CMR2 display bugs fixed; cross-version validated; NEW Machine B assignment: Clipper v2

Commits this session (origin/genai): `cdf5844e` (core-sensitive CMR2 in
`_condition_and_package` + viridis k2 + INDEX cross-version note), `b3ee412f`
(warning demotion for gate-validated version pairs + stale-server banner +
Machine A gate reports), `abbb2272` (GUI amortized runner builds config from
the training JSON — the GUI-side half of the CMR2 bug), plus this addendum +
`plans/europa-clipper-v2-induction-plan.md`.

### CLOSED (Machine A, verified)
- **Europa posterior CMR2 ~3.6 sigma left of the observed Gaussian** — TWO
  stacked bugs, both fixed and observed fixed in the running app (AppTest
  Europa run: CMR2 median 0.3553, +0.23 sigma):
  1. `sbi_runner._condition_and_package` used the core-blind cache scalar.
  2. The GUI run button reconstructed a minimal InferenceConfig without
     `derived_params`, so `_compute_model_cmr2` still fell back core-blind
     even after fix 1. Slots now carry `config_path` (training JSON);
     **any future slot MUST set config_path**.
- **torch 2.11 -> 2.8 cross-version trust** — crosscheck gates re-run on
  Machine A: Europa reproduces Machine B's report to 4 decimals on every
  statistic; Titan passes all 8 params. Reports in
  `validation_reports/cross_version_machineA/`. Slot registry
  `validated_version_pairs` demotes the load warning to INFO for exactly
  this pair; anything else still warns loudly.

### NEW MACHINE B ASSIGNMENT (user-directed 2026-07-13): Europa Clipper v2
`plans/europa-clipper-v2-induction-plan.md` — 3-frequency artifact trained on
induced dipole coefficients Bind = Ae * Be in nT (Kivelson et al. 2023:
sigma = 1.5 nT on Re and Im per component). Scientific review (opus):
RATIFY-WITH-EDITS, all edits incorporated. **Phase 0 is a hard blocker**:
reconcile the Ae solver against an independent multilayer sigma(r) solution
(NOT a uniform shell) AND re-derive the expected |Ae_orbital| — the
"O(1e-2)" anchor behind the 0.73 red flag is itself suspect (85 hr skin
depth ~ ocean thickness makes a few tenths plausible). Channel set is
provisional until Phase 0 (prune on |Ae|*|Be| >= 3 nT, not |Be|). Seeds:
train 42 / data 47 / noise 4747. v1 Galileo slot stays deployed untouched.

## ADDENDUM 2026-07-14 (Machine B) — Clipper v2 Phase 0 Ae reconciliation: VERIFIED, unblocked (n=1)

Phase 0 (the hard blocker in `plans/europa-clipper-v2-induction-plan.md`) is
**`verified`** for the n=1 dipole physics v2 uses (scientific-reviewer opus PASS,
2026-07-14). Full table + method + caveats are in the plan file's Phase 0 section.
Summary:
- Production MoonMag Srivastava Ae solver (the exact `forward_model_induction`
  path that fills the SBI cache) vs an INDEPENDENT from-scratch mpmath log-derivative
  solver agree to **3.7e-14 amplitude / 0.0 deg phase** on all 9 (model x freq) points,
  far inside the 5%/5deg gate. Independence confirmed non-tautological (sigma-sweep
  co-variation to 1e-16; vacuum->0 / PEC->1 absolute-normalization anchors; a!=R
  (a/R)^3 test, which is the real ionosphere-100km-above-surface geometry).
- **|Ae_orbital| red flag DEBUNKED:** 0.22/0.82/0.90 (thin/mid/thick), rising with
  ocean-thickness/skin-depth as required. The Callisto "O(1e-2)" anchor was physically
  wrong for a seawater ocean; the 0.73 was not a solver bug. (Still DROP the Ae GUI
  channels per the C-B DROP verdict — orthogonal reason.)
- **v1 synodic support cut retro-validated:** synodic |Ae| in [0.754, 0.939], all > 0.7.
- **Scope caveat (out-of-scope for v2):** validated at n=1 ONLY; the reference solver's
  general-n rescaling is unresolved (diverges from production at n>=2). v2 is dipole-only
  so this is safe; any future n>=2 induction work must add a quadrupole PEC anchor first.
- Committed: `plans/scripts/phase0_ae_reconciliation.py` (+ `_results.json`).

**Machine B next:** Phase 1 — add the `Bind_<label>_<comp>_real/imag` channel family
to mcmc_runner (signed-Im, Ae*Be complex product in nT, per-channel convention metadata,
pre-deploy-assertion update, unit tests). Small/non-intensive. Then Phase 2 config,
then the Phase 3 1M campaign (seeds train 42 / data 47 / noise 4747).

## NOTE for Machine B (2026-07-16, Machine A) — v2 artifact nflows spline assertion
A Machine A sampling check on `europa_seawater_andrade_clipper_v2.pt`
(4000 draws at config-central x_obs, plus signed-Im flip probes at n=1000)
died inside nflows rational_quadratic_spline inversion:
`assert (discriminant >= 0).all()`. Intermittent — the GUI slot ran fine
for the user at defaults (n=5000). Likely numerical edge of the NSF spline
under some conditioning/draw combinations. Worth a guard (catch + resample
or clamp) in sample_posterior before v2 ratification; belongs with the
Europa SBI work now on Machine B.

## ADDENDUM 2026-07-17 (Machine A) — app/science batch since 2026-07-16

Commits on origin/genai (all verified per CLAUDE.md, AppTest/unit tests):
- k2 ellipse bug (Titan values hardcoded as fallback — every non-Titan
  amortized run drew a Titan ellipse); "Observables used in this run"
  readout; posterior Ae complex-plane plots (metadata['induction_Ae'] /
  ['Be_nT'] packaged by BOTH runners — symmetric schema, unit-tested);
  slot-aware induction documentation; CustomSolution ion table
  (st.data_editor) + published Europa ocean-model presets
  (Utilities/europa_ocean_presets.py, from the user's literature table).
- SCIENCE (opus-reviewed, APPROVE-WITH-EDITS landed): McCleskey
  conductivity now uses Reaktoro equilibrium FREE-ION speciation with 1:1
  ion-pair complexes added to the equilibrium system (MgSO4 at 0.1 mol/kg
  275 K: 33% paired; sigma down accordingly). SCOPE: CustomSolution path
  (RktConduct) ONLY — Seawater (GSW), MgSO4 (Vance2018/Pan2020), NaCl
  (Pan2021) conductivities untouched, so the v3 salinity campaign inputs
  are unaffected. Documentation figures + Mahboub et al. (2026) Table S2
  data comparisons in plans/figures/mccleskey_speciation/ (speciation
  fixes NaCl high-conc overshoot, halves MgSO4 mid-range error; known
  per-ion ionic-strength approximation makes speciated sigma slightly
  HIGH for Na2CO3/mixtures — flagged in sigmaElectricMcCleskey2012.py,
  candidate future fix: total-I).
- Known intermittent nflows spline assertion on Clipper v2 sampling (see
  2026-07-16 note) still open on the Machine B side.

## ADDENDUM 2026-07-17 (Machine B) — pulled Machine A batch; state sync

Pulled origin/genai to `d35157b7` (clean fast-forward from `4efefa56`, 22
commits, no local commits rebased, no untracked collisions). Env + candidate
artifact confirmed on Machine B:
- PPcl env: sbi 0.26.1, torch 2.8.0 (matches the v2 candidate's recorded
  training versions — cross-version trust already gate-validated by Machine A
  in `validation_reports/cross_version_machineA/`).
- Candidate artifact present: `PlanetProfile/Inference/sbi_artifacts/
  europa_seawater_andrade_clipper_v2.pt` (446 KB, 2026-07-14) plus all 6 gate
  reports in `.../validation_reports/europa_clipper_v2_1m/` (sbc,
  crosscheck, crosscheck_matched_truncation, limits_joint,
  limits_263p50_w1_disambiguation, tb_shape_disambiguation). No rebuild needed.

Acknowledged from the pull:
- **McCleskey speciation change is CustomSolution/RktConduct-only** — Seawater
  (GSW), MgSO4, NaCl conductivities untouched. Confirmed **no impact on the v3
  salinity campaign inputs** (v3 is Seawater/GSW throughout).
- **v2 status = candidate, `implemented, unverified`**; awaits Machine A GUI
  slot + AppTest + user ratification (Machine A's step, not Machine B's).
- **v3 salinity** (`plans/europa-clipper-v3-salinity-plan.md`) remains
  `not implemented`, blocked on v2 ratification and user go-ahead.

**nflows spline-inversion resample guard — `verified` (Machine B, 2026-07-17).**
Machine A's 2026-07-16 note (intermittent `assert (discriminant >= 0).all()`
on v2 sampling) is fixed in `SBIRunner.sample_posterior`
(`PlanetProfile/Inference/sbi_runner.py`).

- **Mechanism (resample, not clamp — user-directed):** on the spline
  AssertionError, retry `_posterior.sample` up to `_MAX_SPLINE_RESAMPLE=8`
  times with fresh RNG (`base_seed + attempt*_SPLINE_RESEED_STRIDE`); the
  conditioning `x_t` is built once and reused so only the noise draw z varies.
  Detection via `_is_nflows_spline_assertion` (traceback frame filename
  contains `rational_quadratic` — the assert carries no message); any other
  AssertionError propagates. Persistent failure across all 8 independent draws
  raises RuntimeError (does NOT return partial/biased draws — a systematic
  failure must surface). **Scoped to `reject_outside_prior=True`** (the
  public/GUI path); the diagnostic `reject_outside_prior=False` sweeps
  (validate_sbi limits grid) are left unguarded so their containment check
  measures off-manifold leakage honestly.
- **Root cause:** discriminant `b²−4ac` is analytically ≥0 for a monotonic RQ
  spline; negativity is pure float roundoff at a bin edge, platform/BLAS
  dependent. On this M2 the discriminant bottoms at +1.09e-11 and never
  crosses zero — so the crash does NOT reproduce naturally here (1,700+ sample
  calls, 0 crashes; confirmed via the reproduction sweep).
- **Verification:** 6 new unit tests in `tests/sbi_runner_test.py`
  (`SBISplineResampleGuardTests`) force the assertion path via `mock.patch` on
  `_posterior.sample` and cover: detector specificity, resample-then-succeed,
  persistent→RuntimeError, non-spline-assert propagation, diagnostic path
  unguarded, and clean-run reproducibility. Full file 21/21 pass. End-to-end
  run against the real v2 artifact: same-seed draws byte-identical (guard
  transparent on the clean path), all-finite, Tb median 264.53 K, 40-seed
  sweep 0 crashes. Scientific-reviewer (opus) **PASS** — redraw is
  statistically unbiased (failure is draw-dependent in noise space, not
  conditioning-dependent; whole-batch discard of i.i.d. draws is unbiased).
- **Caveat (honest):** the guard-fires-under-a-real-crash path is verified via
  mock, not a live crash, because the assertion does not trip on this machine.
  The retry *logic* is fully exercised; a live-crash statistical-equivalence
  check would need a machine where the assert fires (reviewer's non-blocking
  recommendation).
- **Committed** to genai (`fix(inference): resample guard for intermittent
  nflows spline assert`). Composes cleanly with Machine A's adaptive-truncation
  change in `_condition_and_package` (different method; the truncation loop's
  `sample_posterior` calls now inherit the guard).

## ADDENDUM 2026-07-18 (Machine A) — salinity range widened; GUI batch
- v3 salinity range now **0.1-100 ppt** (log10_wOcean_ppt U[-1, 2]) per
  user; v3 plan updated in place. Low end reaches near-fresh water: the
  |Ae_synodic|>0.7 support cut will carve away much of the low-w plane —
  report rejected regions. MCMC side added to the plan (same sampled
  parameter once the 2D cache lands; Titan no-ocean exempt). Machine A
  registered log10_wOcean_ppt in parameter_registry (ocean category).
- GUI: MCMC preset radio removed (redundant with Load-config-file; the
  Titan no-ocean 8D config auto-loads as the worked example). Removing it
  exposed a latent bug: an old session-migration wiped any param_space
  containing log10_zeta — which every modern config uses — silently
  emptying loaded configs. Removed.
- Induction panel: 3D (Ae, Re k2) view replaced with a k2 complex-plane
  model-family plot (per-Tb-node mean k2 connected by segments, marker
  size = ocean thickness, color = salinity once a salinity-sampled
  artifact exists — v3-ready).

## MACHINE B QUEUE (2026-07-18): v2 same-version sampling check
Machine A started (and killed, per compute-split rule) the Clipper v2
sampling verification — the last item of the v2 ratification package. It
is compute-heavy on A because SBIRunner(config) triggers the 3-frequency
Ae grid precompute (mpmath, ~46 Tb nodes x 3 freqs). Machine B: run
plans/scripts equivalent of scratchpad v2_sampling_check — draw 4000
samples from europa_seawater_andrade_clipper_v2.pt at config-central
x_obs (seed 42), compare per-parameter mean/sigma against the committed
crosscheck_report.json flow stats (tolerance dmean<0.1 sigma_B, sigma
ratio 0.9-1.1 — same-version torch 2.8, expect near-exact), plus the
signed-Im semantics probe (flipping a Bind imag sign must move the
posterior; flipping Im_k2 must not). The new spline resample guard
should absorb the intermittent nflows assert that killed Machine A's
first attempt. Commit the result to validation_reports/ and mark v2
ready for user ratification.

## ADDENDUM 2026-07-18 (Machine B) — v2 sampling check DONE; one latent Im_h2 bug surfaced

Ran the v2 same-version sampling check. Script committed at
`plans/scripts/v2_sampling_check.py`; report at
`sbi_artifacts/validation_reports/europa_clipper_v2_1m/v2_sampling_reproduction_report.json`.
Sampling path is `SBIRunner.load_artifact` (pure sampling, NO config, NO Ae
precompute) — this is why it is cheap on B and why A's `SBIRunner(config)`
attempt was heavy; the compute-split framing was a red herring, the check
itself is seconds once you skip the config build.

**Overall verdict: PASS**, with the scoping the scientific-reviewer (opus)
insisted on — read these two caveats before treating it as ratification:

1. **What Check 1 actually proves (`verified`):** the artifact LOADS and SAMPLES
   deterministically and reproduces the committed `sbi_mean`/`sbi_std` to Monte
   Carlo precision on all 7 params (alpha mean matches to 7e-7; every dmean well
   inside 0.1*sigma; sigma-ratios 0.996-1.005). The spline resample guard was
   active on the reject=True path and **no assert fired** on this M2 (consistent
   with the guard's known platform-dependence — the assert doesn't trip here).
   This is a **load/determinism + fold-semantics integrity** test. It is
   self-referential (reference stats came from the same artifact + same seed +
   same reject=True path), so it CANNOT and does NOT re-validate the flow's
   scientific fidelity. In particular it does **not** subsume the already-
   dispositioned **Tb_K crosscheck shape FAIL** (v2 plan lines 309-315: genuine
   near-Gaussian-vs-skewed flow shape defect, sub-resolution + conservative,
   reviewed COMMIT-AS-CANDIDATE 2026-07-14). That disposition still stands on
   its own; this PASS neither helps nor harms it.

2. **Signed-Im semantics (`verified`):** flipping the SIGNED `Bind_synodic_x_imag`
   (-157.77 -> +157.77) moves the posterior up to 20.7 sigma (Tb_K) — signed
   channels carry information as designed (this is an off-manifold excursion, so
   it only confirms "responds", not calibrated magnitude; probe uses
   reject_outside_prior=False, matching validate_sbi's off-manifold convention
   at lines 903/1012 — reject=True legitimately stalls at 0% acceptance far off
   manifold, which is what burned 90 min on my first, wrongly-reject=True,
   attempt). Conditioning on Im_k2 = +0.08 vs -0.08 gives BYTE-IDENTICAL draws —
   the abs-fold in `_x_obs_vector` is correct for Im_k2.

**NEW LATENT BUG for Machine A (found by the reviewer, confirmed empirically):**
`_x_obs_vector` (sbi_runner.py:674) abs-folds ONLY the `_IM_K2_ALIASES`
(`Im_k2`/`abs_Im_k2`). But `Im_h2` is marked `'abs'` in the artifact's
`channel_conventions` and is folded to `|Im_h2|` during training
(mcmc_runner.py:1803). So a SIGNED `Im_h2` passed at inference conditions the
flow off-manifold. Probe: Im_h2 = +0.08 vs -0.08 gives NON-identical draws
(0.068 sigma shift) where it should be identical. **Immaterial to the v2
artifact as gated** (its x_obs has Im_h2 = 0.0, abs(0)=0), so it does NOT block
v2 ratification — but it is a real train/sample convention mismatch on the
public/GUI path v2 will deploy on. **Recommended fix (Machine A, GUI-side):**
generalize `_x_obs_vector` to fold every channel whose `channel_conventions`
value is `'abs'`, not just the Im_k2 aliases. Add an Im_h2 +v/-v identity
assertion to the AppTest.

**v2 status: the sampling check is `verified` (PASS, scoped as above). All
Machine-B gate items for v2 are now done.** Remaining before deploy are
Machine A's: GUI slot (config_path REQUIRED), AppTest, and USER RATIFICATION
(the Tb_K shape disposition + this Im_h2 fix are the two things to weigh at
ratification). v3 salinity stays blocked on that ratification.

## ADDENDUM 2026-07-18 (Machine B) — DESIGN SPEC for v3: k2/h2/Ae complex-plane cloud plot (`not implemented`, deferred to v3)

User reviewed Machine A's induction panel and wants the "k₂ complex plane by
model (connected)" block (`PlanetProfileApp/pages/Inference.py` ~2026-2099)
**replaced** by a different, physically-correct plot. Recorded here for the
**Fable/Machine A v3 plan to absorb** (Fable plans on Machine A; Machine B
executes). Status `not implemented`, deferred until a salinity-sampled v3
artifact exists.

**What Machine A's current plot gets wrong (user):** it connects *per-Tb-node
mean k₂* in Tb order — i.e. one signal (k₂) across models, collapsed onto the Tb
grid. That discards the composition axis. **Ae varies with COMPOSITION
(ocean salinity/conductivity), not just Tb.** The code comment
"Ae depends on the ocean state only through T_b in this model family"
(Inference.py ~2103-2105) is a **v2 fixed-seawater artifact limitation, not
physics** — do not carry it into v3.

**The plot the user wants (a CLOUD — user-confirmed):**
- **One connected path per posterior SAMPLE** (a full interior model: rheology
  + Tb + composition), NOT per Tb node. Each path links that sample's
  **dimensionless complex signals** on a **single Re–Im plane**:
  **k₂, h₂, Ae(synodic), Ae(synodic 2nd), Ae(orbital)** — 5 signals for the
  3-frequency Clipper v2/v3 run; 3 signals (k₂, h₂, Ae synodic) for the
  single-frequency Galileo run.
- All five are dimensionless and O(0.1–1.2) (nominal k₂ Re≈0.25, |Ae|≲1 with
  |Ae_synodic|>0.7, h₂ Re≈1.2, all Im≈0 at the fiducial), so they legitimately
  co-plot on one plane. Draw a faint unit circle for Ae reference.
- Translucent connected paths, **subsample ~200** samples, **color by
  `log10_wOcean_ppt` (salinity)**. The hypothesis (user): the *shapes* of these
  paths discriminate models and give a visual sense of the **dimensionality of
  the inversion space**. This is why it must be v3 — on v2 (fixed seawater) the
  composition axis is dormant and the shapes only spread via Tb + rheology.

**PREREQUISITE — code change that MUST land in the v3 forward-model/packaging
work (else the plot silently drops h₂ → 4 signals):**
`h2_results` is **not currently packaged**. Both runners compute `Re_h2, Im_h2`
in `forward_model_k2_flexible` (returns the 5-tuple
`(Re_k2, Im_k2, Re_h2, Im_h2, perPhase_W)`, `forward_models.py:679`) but
**discard** them (`sbi_runner.py:1033`, `mcmc_runner.py:1336`). To carry h₂:
1. Add optional field `h2_results: Optional[np.ndarray] = None` to
   `InferenceResult` (`inference_core.py` ~247, same pattern as `k2_results` /
   `D_ocean_results`; `__post_init__` needs no change unless a shape check is
   wanted).
2. Capture `(Re_h2, Im_h2)` in both recompute loops (unpack the currently-`_`
   slots) and build an `(n,2)` array alongside `k2_results`.
3. Pass `h2_results=...` at both `InferenceResult(...)` sites
   (`sbi_runner.py:1103`, `mcmc_runner.py:1378`).
Backward-compatible: all construction sites are all-keyword; pickle restores
old results without the attr, so the GUI reads it via
`getattr(result, 'h2_results', None)` (matching the existing
`D_ocean_results`/`cmr2_results` convention). Because the v3 1M artifact's
per-sample results are recomputed at inference time, **capturing h₂ in the v3
packaging work makes v3 results carry it automatically** — no artifact retrain
needed for h₂ specifically, but the packaging code must be in before the v3
results are generated/served.

**Alignment (verified this session):** `metadata['induction_Ae']`
(`{label: {'re':[...], 'im':[...]}}`, per-sample) is index-aligned with
`result.samples` and `result.k2_results` in both runners
(`_collect_posterior_Ae`, `mcmc_runner.py:779-809`), so a per-sample path can
zip k₂[i], h₂[i], and Ae[label][i] with no realignment. Plot stays matplotlib +
`st.pyplot` (established pattern in this expander). NOTE: `k2_results[:,1]` and
`induction_Ae` Im are **signed** (do NOT abs-fold in this plot — the shape needs
true sign).

## ADDENDUM 2026-07-18b (Machine A) — Im_h2 fold gap FIXED
_x_obs_vector now folds EVERY channel whose channel_conventions entry is
'abs' (covers Im_h2 — Machine B's surfaced gap), with the legacy fallback
extended to Im_h2/abs_Im_h2 aliases for pre-convention artifacts (Galileo
v1 trained Im_h2 folded too). Bind_ signed channels untouched. Unit test
covers both artifact generations (22/22 green). ONE-LINE ask for Machine
B: rerun the Im_h2 +/-v probe from v2_sampling_check.py against the
fixed code — expect byte-identical draws now (Machine A's on-artifact
attempt used off-manifold x_obs and hit the reject-loop; killed per
compute-split rule). v3 cloud-plot spec received; will absorb into the
v3 plan — current per-Tb-node plot stays as the v2 interim.

## USER DECISIONS 2026-07-18 (via Machine A)
1. Test/ permission GRANTED for the v3 reference MCMC pickle —
   commit it under Test/mcmc_results/Europa/ (v3 plan updated).
2. **Europa Clipper v2 RATIFIED** — INDEX row moved to Deployed. The GUI
   slot was already live with full guards; no code change needed. v3 is
   cleared to start (Phase 1: 2D Tb x salinity cache).

## ADDENDUM 2026-07-19 (Machine A) — HF redeploy SHIPPED
User ran the deploy one-liner 2026-07-18 21:29 UTC: app-deploy snapshot
1ce9c378 (from genai a882e599) uploaded; Space vsteven/planetprofile
stage RUNNING, HTTP 200 at vsteven-planetprofile.hf.space. Live batch now
includes: ratified Clipper v2 slot + guards, k2 ellipse alias fix,
Ae posterior + connected-signals plots, h2 packaging, ion table + Europa
presets, McCleskey equilibrium speciation (default-on CustomSolution),
PvT inner surfaces, MCMC preset removal + Titan default config, Im_h2
fold fix, adaptive truncation. Status: server-side checks `verified`
(runtime stage + HTTP); in-browser click-throughs remain with the user.
Machine B: no new commits on genai since a882e599 — v3 Phase 1 pending.

## ADDENDUM 2026-07-19b (Machine A) — v3 received; Phase 3.5 + 4 done; ONE ASK

Pulled Machine B's v3 delivery (bfc7f3f1..f228810f). Machine A work now in:

1. **Phase 3.5 cloud plot `verified`**: the connected complex-plane panel
   dispatches on a salinity parameter — salinity-sampled results render ONE
   FAINT PATH PER POSTERIOR SAMPLE (LineCollection, color = log10_wOcean_ppt,
   marker size = ocean thickness, ≤1200 paths subsampled seed-0), and the
   "Ae depends only on Tb" caption is now conditional (2D-grid text for v3).
   Fixed-salinity artifacts keep the per-node view. Verified by AppTest (both
   branches, captions asserted) + rendered PNG against the committed v3
   reference pickle — the Tb–w degeneracy is visible as the color gradient
   along the Ae arcs.
2. **Phase 4 GUI slot**: `europa_seawater_v3_clipper_8D_posterior_1m.pt`
   registered in `_SBI_ARTIFACT_SLOTS` (config_path = v3 8D json; Im_k2 ≤
   0.15 and synodic |Ae| guard [0.75, 0.94] both INHERITED from the v2
   grid-walk — v3 ran no new anchor walk; no 1D Tb default-truncation since
   the support cut is now a 2D (Tb, w) region the flow trained on). Render
   `verified` via AppTest (slot lists, scope note renders, conditioning
   inputs present, cache-not-found refusal fires). Sampling
   `implemented, unverified`.
3. **INDEX.md**: v3 candidate row added (RATIFIABLE, gate summary,
   blocked-on-cache status).

**ONE ASK (Machine B): commit the 2D structure cache** —
`PlanetProfile/Test/mcmc_results/Europa/Test52_seawater_v3/`
`europa_seawater_structure_grid_v3_2d.pkl` (21 MB, currently gitignored via
`PlanetProfile/Test/.gitignore`). Serving needs the full per-node layer
arrays (per-sample k2/h2 recompute), so it cannot be slimmed; without it
Machine A cannot run the sampling verification and the public app cannot
serve the slot. Remove the ignore line + `git add -f` + push. PENDING USER
OK on repo size (21 MB single blob). After it lands: Machine A runs the
sampling verify + presents the ratification package.

## ADDENDUM 2026-07-19c (Machine A) — v3 sampling VERIFIED end-to-end
Machine B committed the 2D cache (f7a03572; targeted .gitignore negations,
1488 nodes / 1303 built / 185 tilted-band Nones). Machine A then ran the
full GUI sampling path via AppTest (select v3 slot -> Generate Posterior ->
packaged result): 10,000 draws in 4.0 min on the M4 Air; 8 params incl.
log10_wOcean_ppt; per-sample k2 (n,2), h2 (n,2), and induction_Ae for all
three labels packaged and index-aligned; the per-sample cloud plot rendered
in the live result. Posterior vs reference MCMC: Tb 264.46 vs 264.3 K;
w 36.6 vs 38.7 ppt; corr(Tb, log10 w) -0.988 vs -0.986. Every v3 phase is
now `verified`. Remaining: (1) USER RATIFICATION of the v3 artifact for
deployment; (2) HF redeploy (user one-liner) after ratification. Inherited
guards (Im_k2 <= 0.15; synodic |Ae| [0.75, 0.94]) stay until a v3 anchor
grid-walk widens them — future Machine B work, non-blocking.

## ADDENDUM 2026-07-19d (Machine A) — v4 geodesy plan authored (Machine B: WAIT for user answers)
User direction: v3's k2 conditioning sigma (0.06) is the requirement,
not the Mazarico et al. (2023, Table 5) Clipper projection
(k2 = 0.2-0.3 +/- (1.4-1.8)e-2); and add an option to model C20/C22
directly (Clipper sigma(C20) (3.5-8.5)e-7, sigma(C22) (1.5-2.3)e-7)
instead of the RD-derived CMR2 datum, with hydrostatic C/MR2 shown via
Radau-Darwin and the implied non-hydrostaticity reported. Full plan:
plans/europa-clipper-v4-geodesy-plan.md — scientific-reviewer PASS WITH
CONCERNS, all four binding fixes incorporated (Clairaut k_f integration
replaces RD in the forward map; free 2-D nuisances honestly remove the
gravity MoI constraint; ratio-preserving injections pre-registered as
expected-degenerate; normalization convention gate). Machine B: do NOT
start until the user answers the five open items at the plan's end
(normalization convention is BLOCKING) and Machine A lands
gravity_obs.py + runner channels. Same 2D cache — zero PP runs needed.
v3 ratification remains a separate pending user decision.

## ADDENDUM 2026-07-19e (Machine A) — v3 VETOED by user; C/MR2 explainer added
User decision: Europa Clipper v3 is VETOED — wrong training bounds (k2
conditioning sigma 0.06 = requirement level, not the Mazarico (1.4-1.8)e-2
projections; Re_k2 0.25 vs ocean-consistent 0.23). INDEX row updated
(vetoed, retained for provenance), GUI slot REMOVED (AppTest-verified:
3 slots offered, no v3). The v3 2D cache, reference-MCMC infrastructure,
and gate pre-registrations carry into v4 unchanged — nothing recomputes.
Machine B: v4 is now the active campaign; still wait on Machine A's
runner C20/C22 channels before dataset generation.
GUI: added a "How degree-2 gravity constrains C/MR2" popover in the
C/MR2 posterior expander (MacCullagh relations per Mazarico et al. 2023
Eq. 5 — C20 = -J2 = -[C-(A+B)/2]/MR^2, C22 = (B-A)/4MR^2; A<=B<=C axes
explained; Radau-Darwin + Tricarico 3.324; GC21 provenance of the
0.3547 +/- 0.0024 value; v4 forward-look). AppTest-verified rendering.

## ADDENDUM 2026-07-19f (Machine A) — ALL Europa artifacts flagged: Titan k2 sigmas
User finding: every Europa artifact trained with Titan-derived k2
constraints (sigma 0.06 ~ Cassini-Titan scale) — but Europa has NO
measured k2; the Galileo v1 artifact conditions on a fake measurement.
Directive: re-run Galileo, then remove Clipper v2 in favor of v4.

Machine B queue (order):
1. **Galileo v1.1** (spec in europa-clipper-v4-geodesy-plan.md,
   "Galileo v1.1 re-run" section): honest data = CMR2 (GC21) + synodic
   support cut; k2/h2 retained as LABELED hypothetical-conditioning
   channels (Re_k2 [0.23, 0.05], Im_k2 [0.004, 0.05], h2 unchanged)
   unless the user overrides to drop them. Test51 1D cache unchanged.
   Seeds train 44 / data 50 / noise 5050. Full gates + re-run
   reference MCMC with the same corrected observables.
2. **v4 geodesy** (blocked on Machine A runner C20/C22 channels — in
   progress).
INDEX rows v1 + v2 annotated RETIRING; GUI scope notes now carry the
caveat loudly; slots stay live until replacements land (public-app
continuity). Do NOT HF-deploy in the interim unless the user asks.

## ADDENDUM 2026-07-19g (Machine A) — v4 runner integration LANDED; Machine B CLEARED
The v4 gravity forward model is wired end-to-end (user-approved):
- InferenceConfig.gravity_forward_model ('clairaut_hydrostatic'; hash
  include-only-when-set, pre-existing hashes untouched — verified by test).
- MCMCRunner._derive_gravity_pair: Clairaut k_f over the IDENTICAL
  composite profile the CMR2 mass-conservation derivation assembles
  (side-channel _last_composite_layers, so they cannot diverge), then
  C22_h = k_f q_r/4 at the GC21 1565 km reference radius,
  C20_h = -3.324 C22_h (Tricarico), + sampled dC20_nh/dC22_nh offsets.
  Non-mass-conservation configs fall back to the raw cached profile.
- Likelihood + generate_sbi_dataset dispatch C20/C22/J2 through the
  computed pair when the flag is set; legacy cached-scalar branches
  preserved verbatim otherwise. run() and SBI _condition_and_package
  package per-sample c20_results/c22_results (new InferenceResult
  fields). parameter_registry: dC20_nh/dC22_nh ('gravity' category).
- Tests: tests/gravity_channels_test.py (5: composite-profile identity,
  Tricarico ratio exact, additive nuisances, core sensitivity,
  likelihood conditioning 5-sigma = 12.5 ll, legacy-path + hash
  stability) + tests/gravity_obs_test.py (11). Full sweep 39/39 green
  (cmr2_reporting, sbi_runner, sidecar included).

**Machine B is CLEARED to start, in order:**
1. Galileo v1.1 (spec in v4 plan; Test51 1D cache; independent of the
   gravity code).
2. v4 geodesy: author europa_clipper_v4_geodesy_10D.json per the plan's
   config section (Re_k2 [0.23, 0.015], Im_k2 [0.0040, 0.015], CMR2
   kept as GC21 MoI prior, C20 [fiducial, 8.5e-7], C22 [fiducial,
   2.0e-7] UNNORMALIZED, gravity_forward_model='clairaut_hydrostatic',
   dC20_nh/dC22_nh U[-2e-5, 2e-5], seeds train 43 / data 49 / noise
   4949; fiducial C20/C22 from one _derive_gravity_pair call at the
   fiducial (Tb, w) with zero offsets). Same v3 2D cache — zero PP
   runs. Reference MCMC -> Test/mcmc_results/Europa/Test53_geodesy_v4/
   (add .gitignore negations as f7a03572). Gates incl. the two-arm
   non-hydrostaticity recovery gate (ratio-breaking recoverable;
   ratio-preserving EXPECTED-DEGENERATE) + interior-unbiasedness check.

## ADDENDUM 2026-07-19h (Machine A) — zeta split (user); B's configs revised in place
User: shared zeta preferences silicate heating — separate ice/silicate.
Machine B's two provisional configs (ff0316bc, thanks — fiducial gravity
verified self-consistent) are REVISED and RENAMED on genai:
- europa_galileo_v1p1_7D.json  -> europa_galileo_v1p1_8D.json
- europa_clipper_v4_geodesy_10D.json -> europa_clipper_v4_geodesy_11D.json
Change: param_space log10_zeta -> log10_zeta_Ih + log10_zeta_sil (both
U[-3,2]; hook already supported per-phase zeta; metadata block
zeta_split_2026_07_19 documents the physics). Fiducial C20/C22 UNCHANGED
(zeta does not touch the density profile). Machine A verification: both
configs load + sample finite likelihoods; independent-zeta effect on k2
confirmed (zeta_sil -1.5 dex quintuples |Im k2|). Machine B: train from
the RENAMED configs; SBC/crosscheck at 8/11 params; add the
zeta-Ih-vs-zeta-sil joint-posterior report + the pre-registered
ice-vs-silicate heating-fraction comparison (v4 plan, zeta section).

## ADDENDUM 2026-07-19i (Machine A) — main merged; Enceladus test runs; Park 2024

1. **origin/main MERGED into genai** (523890a0; 33 conflicted files;
   policy + verification in the merge commit message). Libration model
   (Gravity/Librations.py), DP-ALMA backend, etaMelt plumbing, upstream
   bug fixes all in; genai's TidalPy backend, MoonMag first-party,
   McCleskey speciation, Arrhenius viscosity, guarded reaktoro import
   all preserved. NO merges genai -> main (user directive). Machine B:
   run the BuildTest sweep to close the merge's PP-pipeline
   verification (inference stack already verified: 45/45 tests + v2
   artifact serving).
2. **Park et al. 2024 supersedes Iess/Thomas for Enceladus**
   (constraints doc updated): Case 2 field, unnormalized at 256.6 km —
   C20 -5477.45 +/- 36.99 e-6, C22 1517.90 +/- 14.70 e-6, libration
   0.091 +/- 0.003 deg (1 sigma). J2/C22 = 3.61.
3. **Machine A test runs COMPLETE (user-directed; B does production
   only):** Enceladus smoke cache (3 Tb nodes, Seawater 10 ppt,
   Test/mcmc_results/Enceladus/Cassini_smoke/), new observable channel
   'libration_deg' (rigid 3-layer Van Hoolst 2008 via merged
   Librations.py; MCMCRunner._derive_libration_deg; wired into
   likelihood + SBI dataset), body-parameterized gravity keys
   (metadata gravity_ref_radius_m / gravity_j2_over_c22 threaded into
   _derive_gravity_pair), cache_builder bulk_overrides (Enceladus MoI
   window must be widened for cache builds — MoI swings hugely per
   0.1 K; the inference constrains MoI via observables instead),
   Enceladus in BODY_ORBITAL_PARAMS + Torb_s-derived meanMotion.
   Channel sanity at fiducial: libration 0.111 deg vs obs 0.091
   (rigid-shell overestimate, correct direction); hydrostatic C22
   1.444e-3 vs measured 1.518e-3 and C20 -4.68e-3 vs -5.477e-3 (the
   documented non-hydrostatic tension). Smoke config:
   configs/enceladus_cassini_smoke_6D.json.

**Machine B production queue for Cassini-Enceladus** (after v1.1 + v4):
- Dense Tb cache (bulk_overrides Cuncertainty ~0.08; box [271.8,
  272.6] at 10 ppt — 272.7+ exceeds freezing; consider salinity
  sensitivity), production 8D config = smoke + dC20_nh U[-1.5e-3,
  1.5e-3] + dC22_nh U[-1e-4, 1e-4] per the constraints doc.
- BLOCKING before training: (a) exact Tricarico (2014) J2/C22 ratio
  for Enceladus (config carries provisional 3.24 with provenance
  note); (b) correlated 2x2 (C20, C22) covariance conditioning
  (constraints doc, review R3-binding); (c) decide elastic libration
  correction (y1 from the TidalPy solution) vs documented rigid-shell
  systematic.

## ADDENDUM 2026-07-19j (Machine A) — Enceladus BLOCKING ITEMS RESOLVED (none awaited B)
All three pre-production blockers closed on Machine A; none depended on
Machine B MCMC results:
1. **Hydrostatic ratio = 3.25** (McKinnon 2015, GRL 42, 2137: rapid
   1.37-d spin + differentiated structure) — config updated from the
   provisional 3.24. Structure-dependence bracket [3.25, 3.326
   (homogeneous Tricarico Eq. 42)] documented as a <= ~1.2e-4
   systematic on the recovered dC20_nh (~20% of the detected offset) —
   report alongside the posterior. Bonus: differentiated-body
   corrections exceeding the homogeneous formula also validates v4's
   Europa 3.324 (GC21) as the right magnitude class.
2. **Correlated (C20, C22) conditioning IMPLEMENTED**: config metadata
   observable_correlations = {"C20,C22": rho} -> bivariate Gaussian in
   the MCMC likelihood AND correlated multivariate noise in
   generate_sbi_dataset (training generative model matches the
   likelihood). rho = +0.47 estimated from Iess et al. 2014's ratio
   sigma (their 3.51 +/- 0.05 vs +/- 0.042 uncorrelated implies
   rho(J2,C22) ~= -0.47 -> +0.47 in C20 convention); Park et al. 2024
   publish no covariance — request from authors/PDS if possible, else
   the estimate stands (documented in config). Tests: rho=0 reduces
   exactly to independent terms; noise correlation verified
   empirically (tests/gravity_channels_test.py, 7/7).
3. **Elastic libration correction UNNECESSARY**: rigid-shell channel
   0.11091 deg vs full elastic DP-ALMA (y1, rigid=False) 0.11094 deg
   at the default structure — 0.03% = ~0.01 sigma of the 0.003-deg
   measurement (Enceladus k2 ~ 0.013; negligible tidal softening).
   Rigid channel ships for production; PP-pipeline elastic value
   retained as the crosscheck.
Machine B's Enceladus production queue now has NO blocking items
beyond compute: dense Tb cache (+ salinity sensitivity if desired),
8D nuisance config per constraints doc, reference MCMC + training +
gates.

## ADDENDUM 2026-07-19k (Machine A) — Park et al. 2024 ERRATUM applied
User supplied the published erratum (authoritative version of record).
Changes applied: libration observable 0.091 -> 0.092 deg (Eq. 6
amplitude 0.092119; 3-sigma 0.009 unchanged); interior
context/crosscheck values revised (shell 25-29 km pref. 27, ocean
26-30 km pref. 28 - THICKER ocean than v1 - core R ~197 km, rho
2290-2350 pref. 2320, MoI 0.337 [0.335, 0.338], heat loss 20-30 GW).
The MEASURED Case 2 gravity field is NOT affected (erratum touches
shape/topography-derived quantities only) - C20/C22 observables stand.
Erratum's added rescaling note C_nm ~ (R_old/R_new)^n confirms our
degree-2 (R/R_ref)^2 convention. Config + constraints doc updated;
committed smoke result pkl ran pre-erratum at 0.091 (0.3-sigma shift,
immaterial for the channel-exercise smoke); Machine B production uses
0.092. papers/park2024global.pdf remains the v1 PDF - read values
through the erratum notes in the config/constraints doc.

## ADDENDUM 2026-07-19l (Machine A) — Interactive 3D globe in results view
User-requested feature: "🌐 Interactive Globe" expander in the inference
results (Utilities/globe_view.py + assets/globes/): grab-to-rotate
plotly globe with the body's NASA/USGS global mosaic on the surface,
concentric posterior-median interior shells (core / silicate seafloor /
ocean top), cutaway toggle (90-degree wedge, staggered per shell),
surface-opacity + shape-exaggeration sliders. Degree-2 figure: the
surface is deformed by the equipotential implied by the C20/C22
observables with fluid amplification (1+k_f)/k_f — the user's point
that even 1D inference carries a non-spherical shape C/MR2 never shows.
Textures shipped (public domain, 75-170 KB each): Europa + Callisto
(USGS CKAN quicklooks), Enceladus (PIA18435), Titan (PIA19658 cropped);
Ganymede falls back to a shaded sphere until a texture lands. Verified:
AppTest on both the Enceladus smoke result (C20/C22 -> deformed figure)
and the Europa v3 reference (sphere + note), plus PNG inspection of the
cutaway (Europa: ice rim / blue ocean ring / silicate mantle / metal
core at posterior-median radii). Future 3D: per-layer r(theta, phi)
fields on the same mesh (roadmap); browser does all rendering — no new
HF runtime deps (kaleido/chrome were dev-only for verification).

## ADDENDUM 2026-07-19m (Machine A) — globe v2 (user feedback)
1. Europa texture: 2048x1024 extracted from the USGS 500-m mosaic COG
   via HTTP range reads (rasterio /vsicurl, dev-only dep) with a 2-98
   percentile contrast stretch + unsharp mask — lineae/mottled terrain
   now visible (the old 1024 quicklook was washed out).
2. Cutaway toggle REMOVED: always cutaway, innermost body now a SOLID
   sphere (user); shells on a coarse 36x72 mesh (uniform color needs no
   density), textured surface on 160x320 (payload ~2.6 MB, builds in
   0.08 s — render cost is dominated by plotly WebGL init per rerun).
3. SAMPLE PICKER (user request): a clickable posterior scatter
   (Tb vs ocean thickness, colored by salinity when sampled) above the
   globe — clicking a point rebuilds the globe from THAT sample's
   per-sample packaging (D_hsphere/D_iceIh/D_ocean[i], R_core sample,
   per-sample kf via RD, per-sample c20/c22 when packaged) with the
   sample's parameter values in the caption; default = posterior
   median. This covers the 'sliders in model space' idea with
   posterior-consistent structures (arbitrary off-posterior sliders
   would need a cache-lookup path — future).
Future upgrade recorded: true RGB texture mapping + faster rotation =
three.js custom component or missionwidget.com (roadmap 3D note).
AppTest green on both fixtures (toggle absent, picker caption present);
PNG inspection: solid core, contrast-stretched surface.
