# Tidal-sector remedy plan — lifting the split ratification

Goal (user, 2026-08-11): the amortized flow results must be trustable as
REPRESENTATIVE OF THE MCMC across all sectors, ending the per-slot tidal
(Im k2 / zeta / eta) quarantine.

Scope boundary (established, §0.12): the flow-vs-MCMC gap is the ONLY part
we can fix. The MCMC-vs-data tension (pushforward Im k2 ~0.09 vs obs 0.135;
Re k2 ~0.57 vs 0.608) is a data/model limit that NO flow remedy touches —
"trust the flow" means "flow ≡ MCMC", not "flow matches the datum".

## Why this is fixable cheaply

1. Diagnosis (reviewer-adjudicated): training-signal/identifiability
   defect. Capacity/embedding ELIMINATED (#4 pilot); salinity, mixture,
   noise convention, x-norm, dataset size, reference resolution all
   ELIMINATED. The flow under-updates the weakly-identified tidal
   dimensions — it stays too close to the prior there.
2. Direction of the error: the deployed flows are CONSERVATIVE (wider than
   the reference MCMC) on the quarantined sector — pushforward
   concentration_ratio > 1, SBI pp more pessimistic than MCMC-pp in the
   same direction. A conservative approximate posterior is an IDEAL
   importance-sampling proposal: it covers the true posterior's mass.
3. The deploy path already computes per-sample forward observables
   (cache-interpolated k2, direct-Clairaut C20/C22) for the results
   panels — the expensive ingredient of an exact correction is already
   paid for.

## Track 1 (primary): amortized flow + importance correction (SIR)

At conditioning time, treat the flow as a proposal and correct it with the
exact likelihood the MCMC uses:

- Draw N samples theta_i ~ q(theta | x) from the deployed flow
  (DirectPosterior, reject_outside_prior as today).
- Forward-evaluate each draw through the SAME cached forward model the
  MCMC used (structure-cache interpolation -> k2, C20/C22 path per the
  config's gravity_forward_model) — this code path already runs per-sample
  for the k2 scatter / heating panels.
- log w_i = log p(x_obs | theta_i) + log p(theta_i) − log q(theta_i | x)
  with the config's Gaussian sigmas; stabilize; systematic resampling.
- Diagnostics gated in-line: importance ESS / N must exceed a floor
  (preregister 0.1 for deploy; report always); max-weight fraction cap.
  Below floor -> the panel falls back to the uncorrected flow WITH the
  existing quarantine warning (never silently).

Properties: amortization preserved (works at any user-conditioned x in
support); asymptotically exact wherever the proposal covers (and
conservativeness is exactly coverage); no retraining, no architecture
change; per-composition validation is direct (compare corrected posterior
vs the pooled B3 reference on the previously-failing statistics).

Cost estimate: forward evals are cache interpolations (~ms/sample); a 10k
draw correction is seconds-scale on the HF host. If needed, correct once
per conditioning click, cached by x.

### Validation gates (preregistered; no tuning after seeing results)

Per composition (NH3, MgSO4, NaCl — Titan first since that is where the
quarantine lives; Europa v5/v6 as no-regression controls):
1. Crosscheck of the CORRECTED posterior vs the pooled reference on ALL 13
   params under the standard tolerances — including the previously
   quarantined/failed ones (NH3 tidal set; NaCl log10_eta_V).
2. Pushforward four-way: corrected-SBI pp median must sit within
   0.5 sigma_obs of MCMC pp median on Re k2 AND Im k2 (median-to-median).
3. ESS/N >= 0.1 at the fiducial obs and across a 5-point in-support x
   sweep (correction must not silently degenerate off-fiducial).
4. SBC of the corrected pipeline at n>=500 (rank-uniformity of the
   corrected posterior; the correction is part of the inference now).
Pass -> per-composition quarantine LIFTED; sector warnings replaced by a
"flow+IS correction, validated against reference MCMC" note. Any FAIL ->
recorded, quarantine stays, escalate to Track 2 findings.

## Track 2 (parallel, informs future training): upstream identifiability

Reviewer ruling 3 (§0.12) — run on Machine B against the existing NH3/salt
datasets, no new simulation:
- Conditional-variance / mutual-information profile of Im k2 given
  {C20, C22, Re k2} under the Tb–w degeneracy: quantify how much signal a
  perfect posterior has (upper bound on any flow's update).
- Heavy-tail audit of Im k2 under the deployed z-scoring; test a
  variance-stabilizing observable transform (log / rank-gauss) by
  retraining ONE arm on NH3 with the transform and re-running the
  pushforward. If the transform materially recovers the update, it becomes
  the standard recipe for FUTURE campaigns (Enceladus onward) — cheaper
  than correcting at inference time forever.

## Track 3 (held): factorized two-stage flow

p(theta_struct | x) · p(theta_tidal | theta_struct, x) — concentrates
training gradient on the tidal conditional. Bigger lift (new training
code, new gates); held unless Tracks 1–2 both fall short.

## Sequencing

1. Reviewer design freeze on this plan (esp. Track 1 weight definition,
   ESS floor, gate set) — no implementation before sign-off.
2. Machine A implements Track 1 (deploy-path code + validation driver);
   NH3 is the pilot composition (quarantined sector + committed matched
   reference class).
3. Machine B: Track 1 validation runs per composition (references are on
   B); Track 2 in parallel as compute frees.
4. Manager re-adjudication per composition; GUI warning changes only
   after gates pass.

## Reviewer adjudication (2026-08-11, agent adeb46980da136fe0): APPROVE WITH CONDITIONS

Binding conditions C1-C16 (full text in the review; summary):
- C1 STOP-GATE first: one NH3 conditioning, n_post=10000, n_derived=All;
  log w = log_likelihoods - flow_log_prob (both ALREADY stored by the
  deploy path, sbi_runner.py:1305-1310/1423 — same likelihood closure the
  MCMC uses). Report ESS/N, absolute ESS, Pareto-k, w_max fraction,
  sentinel fraction, corrected |Im k2| pushforward median. No
  implementation before it returns.
- C2/C3 byte-identity asserts + 200-theta likelihood-recompute test vs the
  committed reference pkl (<1e-9 relative).
- C4 sentinel (-1e30) masked BEFORE arithmetic; ESS on full N; >50%
  rejected = hard fail.
- C5 Pareto-k primary gate (<=0.5 clean; 0.5-0.7 with PSIS smoothing;
  >0.7 FAIL->fallback); w_max/sum <= 0.01; reverse coverage test (<1% of
  reference-MCMC mass below the flow's own 0.1st-percentile log q).
  CORRECTION to this plan's rationale: concentration_ratio>1 is
  SBI-vs-PRIOR (restates the defect), not SBI-vs-MCMC; the
  conservativeness evidence is the MCMC-referenced spread ratio.
- C6 corrected mass outside the k2 training-support box < 1e-3 (box edges
  are 14-25 sigma from the datum; report, don't assume).
- C7 ESS >= 1000 absolute (forces N ~ 10k).
- C8 COST CORRECTION: k2 is a full TidalPy radial_solver solve per sample
  (8 ms/pass laptop, 0.2-0.6 s/eval HF vCPU), NOT a cache interpolation.
  Deployment architecture must be chosen explicitly: (i) offline
  validation + precomputed fiducial correction, uncorrected+quarantine
  for user-modified x; (ii) reduced-N with honest ESS/k reported and
  quarantine retained when floors missed; (iii) background job with
  x-keyed cache (key includes artifact/cache/config hashes + N + seed).
- C9 unidentified-nuisance adjudication rule preregistered with the
  reference 3-seed spread computed BEFORE results.
- C10 weighted-KS/W1 full-distribution gate on Re/Im k2 (median-only is
  the wrong shape test for a width defect).
- C11 C20/C22 no-regression gate.
- C12 amortized sweep: >=200 prior-predictive x + 8 axis endpoints; >=95%
  with k <= 0.7 (5 points is far too thin).
- C13 3-seed stability of corrected medians within 0.1 sigma_obs.
- C14 one preregistered SBC rank definition at deploy-time M.
- C15 SBC ground truths from the EFFECTIVE prior (all eight training
  support cuts), not the nominal box.
- C16 ocean-fraction gate: corrected has_ocean mass must match the
  reference within its 3-seed spread; per-branch ESS >= 50 above 5%
  probability. (Expected science-visible effect: the frozen branch is
  ~-10 sigma at the fiducial datum -> corrected ocean fraction ~1.0.)
Track 2 reordering: information-gain-in-nats comparison first (free);
Im k2 kurtosis measured before any transform retrain; >=3 seeds for any
transform arm; NEW Track 2c ranked above the transform — a diagnostic
flow conditioned on {Re k2, Im k2} only, separating gradient-competition
from no-information.
Preregistered readings: post-correction gate FAIL => implementation
defect (target is the MCMC's target by construction); gates PASS
validates the corrected PIPELINE, not the flow. Quarantine lift never
removes the model-data tension caveat (MCMC itself under-predicts k2).
Positive control: the correction should REPAIR v5's dC22_nh SBC FAIL.
Downstream panels must all consume corrected weights (heating subset,
wedge, k2 cloud, sample picker); weighted statistics preferred over
resampling; display ESS, not N.

## C1 STOP-GATE RESULT (Machine A, 2026-08-11): PREMISE CONFIRMED

NH3, deployed artifact (config hash e596574d), n_post=10000, full derived
recompute, seed as committed. log w = log_likelihoods − flow_log_prob
(both from the existing deploy path; sentinels masked first):

| statistic | value | bar | verdict |
|---|---|---|---|
| corrected Im k2 pp median | **0.1064** | MCMC matched ceiling ~0.1037 (flow alone 0.0439) | **defect closed** |
| corrected Re k2 pp median | 0.5936 | MCMC-pp ~0.575–0.58 | close; formal gate on B |
| Pareto-k (Zhang–Stephens, PSIS tail M=297) | **−0.18** | ≤0.5 clean | PASS |
| w_max/Σw | 0.0091 | ≤0.01 | PASS |
| ESS (full-N) | 680 | ≥1000 (C7) | short → N≈15–20k needed |
| ESS/N | 0.068 | 0.1 preregistered deploy floor | below; N-scaling covers validation, deploy floor re-visited with reviewer |
| sentinel-rejected fraction | 1.35% | <50% hard fail | PASS |
| ocean fraction flow → corrected | 0.933 → 0.927 | C16 vs reference (on B) | reviewer's collapse-to-1.0 prediction NOT borne out — surviving frozen mass is rheology-variant frozen nodes; per-branch check pending |

Wall: 2.3 min per 10k on Machine A (forward evals dominate). Artifacts:
scratchpad c1_result.json + c1_logw.npy (session), to be regenerated by the
committed validation driver in implementation.

Consequence: Track 1 proceeds to implementation under C2–C16. Deployment
architecture (C8): recommend (i)+(iii) — precomputed correction for the
fiducial x shipped with each slot, background x-keyed correction for
user-modified x with in-panel ESS/Pareto-k display and automatic
quarantine-warning fallback when floors are missed.
