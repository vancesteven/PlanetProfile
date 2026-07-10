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

## Open ratification items (after 1M)

1. eta_Ih/eta_V KS materiality margin if 1M still fails (opus framework: bootstrap
   self-D floor + per-param margin, e.g. |dmedian| < 0.3 dex with sigma-ratio in band).
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
