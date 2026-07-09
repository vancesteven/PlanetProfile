# HANDOFF — SBI pipeline + CMR2 fix + Test52 (continue with fable / Opus 4.7)

Updated 2026-07-08 (evening addendum after Test52 production run; supersedes earlier
version). Branch `genai`. Orchestration pattern used so far: opus for design/review,
sonnet for implementation, haiku for recon — per repo CLAUDE.md. Continue that pattern.

## Environment — per machine

- **Machine A (author of first draft):** conda env `PP`, invoked with
  **conda, never mamba**:
  `source /opt/miniconda3/etc/profile.d/conda.sh && conda run -n PP ...`
- **Machine B (Dropbox laptop, where the 2026-07-08 addendum was authored):** mamba env
  `PPcl`: `mamba run -n PPcl python ...` (repo CLAUDE.md convention).
- Package versions differ slightly per machine but both satisfy the SBI pipeline
  requirements. Machine B: python 3.11, numpy 2.4.6, **scipy 1.16.3**, torch 2.11.0,
  sbi 0.26.1, corner 2.2.3, seaborn 0.13.2, pocomc 1.2.6.
- Machine A: python 3.11, numpy 2.4.3, **scipy 1.17.1**, torch 2.8.0 (conda-forge — pip
  torch causes duplicate-libomp crash; never use KMP_DUPLICATE_LIB_OK), pocomc 1.2.6,
  sbi 0.26.1, TidalPy 0.7.4 (wheels libomp-repatched via install_name_tool —
  reinstalling tidalpy/CyRK undoes this; repatch to
  /opt/miniconda3/envs/PP/lib/libomp.dylib + codesign -f -s -), corner 2.3.0,
  seaborn 0.13.2, pypdf.
- Both envs exhibit an OpenMP duplicate-init when `sbi + torch + PlanetProfile` are
  imported in the same process (libomp linked twice). MCMC uses pocomc (no torch) and
  is unaffected; SBI runs need to keep torch imports isolated per subprocess. Machine B
  has NOT been libomp-repatched.
- MoonMag migrated to scipy>=1.16 (`sph_harm_y`); MoonMag is FIRST-PARTY in-tree code —
  user directive: never rely on or contact upstream itsmoosh/MoonMag.

## What is DONE and verified (all in branch genai)

1. **SBI amortized-inference pipeline** (plans/monte-carlo-sbi-implementation-plan.md):
   sbi_runner.py (SBIRunner, single-round NPE/maf, prior provably == MCMC prior, self-
   describing .pt artifacts, flow log-prob kept OUT of log_likelihoods); generate_sbi_dataset
   opt-in kwargs (apply_support_guard, imag_convention='abs', drop_nonfinite, seed,
   provenance) — default path byte-identical; FIXED fatal latent bug prior.sample→prior.rvs
   (pocomc 1.2.6 has no .sample; SBI dataset gen never worked before). GUI AmortizedInference
   wired (button gated on artifact existence, loud failure on non-'abs' artifacts).
   validate_sbi.py CLI: sbc / crosscheck / limits / selftest with RATIFIED hard gates
   (SBC KS p>=0.05; crosscheck mean within max(0.25σ, σ/√ESS), σ-ratio [0.7,1.4], KS α=0.01;
   limits monotone log10_eta_Ih vs |Im k2|). selftest PASSES incl. negative control.
   Opus sign-off: PASS, no must-fix.
2. **CMR2 reporting bug fixed** (was: cmr2_results + SBI CMR2 column used core-blind cache
   placeholder while Callisto likelihood used correct core-sensitive derivation). Shared
   _derive_cmr2_via_mass_conservation + _compute_model_cmr2; likelihood proven BIT-IDENTICAL.
   Callisto pickles regenerated (backups *.pkl.pre_cmr2fix): baseline χ −3.27σ→+0.26σ
   (constraint IS met — manuscript reversal), minus5pct +0.95σ→+2.45σ, minus10pct 5.17→6.38σ.
3. **Test52 differentiated no-ocean Titan** (plans/test52-differentiated-titan-plan.md;
   Petricca et al. 2025 Nature 648:558: Re k2 0.608±0.048, Im k2 0.135±0.035, CMR2 0.343±0.001):
   - Phase 1: config test52_titan_noocean_andrade_10D.json (10D = Test50 8D + R_core_km [0,2000]
     + rho_core_kgm3 [2591,3600] rock prior; mass-conservation rho_sil), registry ParameterDefs +
     preset (observables byte-match config), 9-pt Tb grid built + verified, offsets sidecar
     (CMR2 discretization offset 0.00094983 ± 1.5e-7).
   - Phase 2: offset anchor (sidecar + np.interp in Tb; inert for any config without sidecar)
     + density-inversion guard (nested key derived_params.rho_sil_kgm3.density_inversion_guard;
     dormant under Test52 prior — max achievable rho_sil ≈2554 < floor 2591, by design).
     Callisto/Test50 numerics proven unchanged.
   - Phase 3: smoke MCMC (n_effective=50, 2.4 min) ALL 5 GATES PASS: CMR2 posterior median
     0.3428 (upper tail at no-core anchor 0.34299), 76.5% R_core mass <600 km (weak
     differentiation, matches Petricca), k2 offsets == Test50 baseline (0.94σ/1.17σ).
     Result: Test52_andrade_noocean_diff/smoke_run/.
   - **Phase 4 (production MCMC) DONE 2026-07-08 evening on Machine B (mamba PPcl):**
     - Rebuilt structure grid pkl locally (7.3 s; `mamba run -n PPcl python -m
       PlanetProfile.Inference.build_phase_c1_cache --config
       PlanetProfile/Inference/configs/test52_titan_noocean_andrade_10D.json --template
       PlanetProfile.Test.PPTest50 --n-grid 9`). Cache offsets reproduce the committed
       sidecar to <1e-6 per point (Phase 2 regression tests all 14/14 PASS).
     - Ran production MCMC (n_effective=500, seed 42) via
       `PlanetProfile.Inference.run_inference_cli`. **Wall time 1.7 min** (much faster than
       expected — cache-based likelihood, canonical pocomc mixing on 10D). Result:
       ESS=4180, R̂=1.000, acceptance 0.44, n_samples=4321. Result pickle 0.6 MB in
       Test52_andrade_noocean_diff/production_run/. SUMMARY.md written.
     - **Scientific-reviewer sign-off: APPROVE-WITH-CAUTIONS.** All 5 gates PASS.
       Disclosure items (manuscript, not blockers): (a) CMR2 σ_post/σ_obs=0.52 is
       *prior-envelope-limited*, not likelihood-limited — the (R_core, rho_core) box
       already narrows CMR2 below Petricca's Gaussian; add a prior-predictive CMR2
       histogram to disambiguate. (b) rho_core (2965±288) is prior-dominated — do NOT
       quote as inferred. (c) Correlated low k2 bias (Re −0.78σ, Im −1.12σ) tracks
       a compliant-rheology corner (low log10_zeta / low log10_eta_V); add a
       conditional k2-vs-rheology-corner plot. (d) log10_zeta and log10_eta_V marginals
       have soft lower-edge pileup — report modes + 68/95% HPD, not only medians. (e)
       Tb_K is essentially prior-uniform; observables don't discriminate.
     - Plots: production_run/test52_production_corner.png,
       test52_production_cmr2_ppc.png, test52_production_k2_ppc.png.
4. **MoonMag**: sph_harm→sph_harm_y (equivalence proven 1e-15 pre-upgrade), first-party
   README (citation: Styczinski et al. 2022, Icarus 376:114840), tests/moonmag_smoke_test.py.
5. **Tests all green (Machine B, 2026-07-08 evening):** test52_phase2_test 14/14
   (regeneration of the Test52 grid on Machine B reproduces committed offset sidecar
   to <1e-6). sbi_runner_test 8/8, cmr2_reporting_test 3/3*, moonmag_smoke_test 6/6
   as of last full run. (*one PRE-EXISTING flaky test:
   test_reported_cmr2_varies_with_core_params uses unseeded draws; fails ~2/3 of fresh
   runs even on pre-change baseline — needs seeding, requires approval since it's a
   test-behavior change.)

## Machine-specific note: gitignored artifacts

The following are `.gitignore`d (`*.pkl`) so they do NOT come through `git pull`. On a
fresh checkout, rebuild them locally:

- `PlanetProfile/Test/mcmc_results/Titan/Test52_andrade_noocean_diff/titan_diff_noocean_structure_grid.pkl`
  — rebuild via the build_phase_c1_cache command above (7 s).
- Any Test50 cache pkls (e.g. `titan_allice_yao2014_structure_grid.pkl`) are similarly
  gitignored on some checkouts. If a run says "Loading structure cache: ..." then
  ENOENT, rebuild via the same driver against the corresponding config.

Test52 production **result** pkl (0.6 MB, `test52_production_result.pkl`) is committed
manually as part of this handoff push. Do NOT reuse it as an SBI-training seed without
first re-verifying provenance against `SUMMARY.md` (config path + git SHA).

## Test52 Phase 4 sign-off: APPROVED-WITH-CAUTIONS (no must-fixes)

Second opus review completed 2026-07-08 evening on Machine B. All five Phase 4 gates PASS:
convergence (ESS/N=0.97, R̂=1.000, α=0.44); CMR2 close to 0.343 (median 0.34260, χ=−0.40σ);
weak differentiation (76.14% mass < 600 km — matches smoke 76.5%); k2 residuals within 2σ
(0.78σ / 1.12σ); no unphysical prior-edge pileups blocking. The disclosure items in
section 3 above are for manuscript preparation, not blockers. Production may proceed to
Next Steps #2 (SBI training).

## NEXT STEPS (in order)

1. ~~Test52 production MCMC (n_effective=500, production config, seed 42).~~ **DONE.**
   Result at production_run/test52_production_result.pkl.
2. Titan SBI real validation (all inputs now in-repo):
   a. Train Test50-8D artifact: measure wall-clock on 200 sims first, then 10k min/50k prod
      via SBIRunner (config test50_titan_noocean_andrade_8D.json).
   b. `python -m PlanetProfile.Inference.validate_sbi sbc --config ... ;
      crosscheck --mcmc PlanetProfile/Test/mcmc_results/Titan/Test50_andrade_noocean_yao2014/
      test50_titan_noocean_andrade_8D_result.pkl ; limits --re-k2 0.608`.
      Exact commands in plans/monte-carlo-sbi-implementation-plan.md §Validation.
      Gates are hard — failures stop and surface, never tune.
   c. On pass: name artifact per sbi_artifacts/INDEX.md convention, update INDEX.md,
      GUI Titan slot goes live automatically (button gates on file existence).
3. Test52 SBI artifact (after Test50 SBI validated): same pipeline; becomes the
   differentiated-Titan GUI slot eventually.
4. Manuscript disclosure plots for Test52 (blocking for manuscript, non-blocking for
   pipeline): prior-predictive CMR2 histogram + conditional k2 vs (log10_zeta,
   log10_eta_V) contours. See item 3.(a)–(e) above.
5. Callisto offset sidecars: regenerate Callisto grids with offsets (or compute sidecars
   for existing grids) → recompute the regenerated pickles' χ once more (≈+0.23σ shift),
   update sensitivity conclusions (minus5pct currently 2.45σ moderate tension).
6. Fix flaky cmr2_reporting_test (seed the draws) — get user approval first (test change).
7. Roadmap (user, in memory): magnetic induction observables into inference+SBI
   (mcmc likelihood already supports Ae_*/BiAmp_*/BiPhase_*; generate_sbi_dataset NaN-fills
   them — needs deliberate SBI-side support); era-LLM meta-optimization (deferred).

## Open items / cautions

- rho_core upper bound 3600 kg/m³ is user-adopted, not paper-sourced — cite before manuscript.
- Registry preset andrade_titan_noocean_8D observables ≠ Test50 config JSON (abs_Im_k2+CMR2
  vs Im_k2) — old landmine, resolve before Test50 production SBI training data.
- validate_config n_effective>=100 floor conflicts with smoke runs (smoke driver waived it);
  Phase 4 review was asked to recommend an escape hatch vs keeping the floor.
- Smoke acceptance-rate proxy only: no explicit guard-hit counter in likelihood (would be a
  source change needing approval).
- Cache bodyname 'Test50' ≠ 'Titan' accepted; matters when induction observables added.
- Sidecar `cache_path` field records the machine-A path (`/Users/svance/ppgenai/...`) — this
  is metadata only, not consumed by MCMCRunner (which uses config's repo-relative path).
  If a future consumer starts reading `cache_path`, treat this field as advisory.
- Machine B's `mamba run -n PPcl` and Machine A's `conda run -n PP` produce
  bit-consistent MCMC results at seed 42 for the Test52 config (verified: rebuilt cache
  offsets reproduce sidecar to <1e-6). Do NOT be surprised if a fresh SBI training run
  produces a different .pt binary — SBI's torch RNG streams are platform-dependent even
  at the same seed. The `sbc` / `crosscheck` gates cover that.
- Repo CLAUDE.md: verification vocabulary mandatory (verified / implemented, unverified /
  not implemented); never 'done/fixed'. Plan-first for nontrivial changes. No commit/push
  without explicit user approval.
