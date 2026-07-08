# HANDOFF — SBI pipeline + CMR2 fix + Test52 (continue with Opus 4.7)

Updated 2026-07-08 (supersedes earlier version of this file). Branch `genai`.
Orchestration pattern used so far: opus for design/review, sonnet for implementation,
haiku for recon — per repo CLAUDE.md. Continue that pattern.

## Environment (this machine)

- conda env `PP`, invoked with **conda, never mamba**:
  `source /opt/miniconda3/etc/profile.d/conda.sh && conda run -n PP ...`
  (CLAUDE.md's `mamba activate PPcl` applies to the user's OTHER computer.)
- PP contents (all verified working together): python 3.11, numpy 2.4.3, **scipy 1.17.1**,
  torch 2.8.0 (conda-forge — pip torch causes duplicate-libomp crash; never use
  KMP_DUPLICATE_LIB_OK), pocomc 1.2.6, sbi 0.26.1, TidalPy 0.7.4 (wheels libomp-repatched
  via install_name_tool — reinstalling tidalpy/CyRK undoes this; repatch to
  /opt/miniconda3/envs/PP/lib/libomp.dylib + codesign -f -s -), corner 2.3.0,
  seaborn 0.13.2, pypdf.
- MoonMag migrated to scipy>=1.17 (`sph_harm_y`); MoonMag is FIRST-PARTY in-tree code —
  user directive: never rely on or contact upstream itsmoosh/MoonMag.

## What is DONE and verified (all in this push)

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
     NOT a converged posterior (smoke only). Result: Test52_andrade_noocean_diff/smoke_run/.
4. **MoonMag**: sph_harm→sph_harm_y (equivalence proven 1e-15 pre-upgrade), first-party
   README (citation: Styczinski et al. 2022, Icarus 376:114840), tests/moonmag_smoke_test.py.
5. **Tests all green**: sbi_runner_test 8/8, cmr2_reporting_test 3/3*, test52_phase2_test 14/14,
   moonmag_smoke_test 6/6. (*one PRE-EXISTING flaky test: test_reported_cmr2_varies_with_core_params
   uses unseeded draws; fails ~2/3 of fresh runs even on pre-change baseline — needs seeding,
   requires approval since it's a test-behavior change.)

## Test52 Phase 4 sign-off: APPROVED (no must-fixes)

Opus review completed 2026-07-08. All four items PASS: (1) numerics diff — offset added after
all rejection guards (cannot perturb them), np.interp clamped within prior, exactly one sidecar
in repo (no accidental pickup), guard reset each call; (2) smoke posterior — CMR2 constraint
bites (R_core posterior/prior std ratio 0.46), PPC edge at offset-corrected 0.34299 proving the
anchor is live in the likelihood, no pathologies (no R_core=0 pileup, no rho_sil saturation),
reproduces Petricca weak differentiation; (3) manuscript disclosures enumerated (see review in
session log; summarized under Open items); (4) production run may proceed as-is. Also:
keep validate_config's n_effective>=100 floor, waive at smoke-driver level only (no escape
hatch in the validator). Plots: smoke_run/test52_smoke_corner.png, test52_smoke_cmr2_ppc.png.

## NEXT STEPS (in order)

1. Test52 production MCMC (n_effective=500, production config, seed 42). Compare vs smoke.
   APPROVED by Phase 4 review — no blockers.
2. Titan SBI real validation (all inputs now in-repo):
   a. Train Test50-8D artifact: measure wall-clock on 200 sims first, then 10k min/50k prod
      via SBIRunner (config test50_titan_noocean_andrade_8D.json).
   b. conda run -n PP python -m PlanetProfile.Inference.validate_sbi sbc --config ... ;
      crosscheck --mcmc PlanetProfile/Test/mcmc_results/Titan/Test50_andrade_noocean_yao2014/
      test50_titan_noocean_andrade_8D_result.pkl ; limits --re-k2 0.608.
      Exact commands in plans/monte-carlo-sbi-implementation-plan.md §Validation.
      Gates are hard — failures stop and surface, never tune.
   c. On pass: name artifact per sbi_artifacts/INDEX.md convention, update INDEX.md,
      GUI Titan slot goes live automatically (button gates on file existence).
3. Test52 SBI artifact (after production MCMC): same pipeline; becomes the differentiated-
   Titan GUI slot eventually.
4. Callisto offset sidecars: regenerate Callisto grids with offsets (or compute sidecars for
   existing grids) → recompute the regenerated pickles' χ once more (≈+0.23σ shift), update
   sensitivity conclusions (minus5pct currently 2.45σ moderate tension).
5. Fix flaky cmr2_reporting_test (seed the draws) — get user approval first (test change).
6. Roadmap (user, in memory): magnetic induction observables into inference+SBI
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
- Repo CLAUDE.md: verification vocabulary mandatory (verified / implemented, unverified /
  not implemented); never 'done/fixed'. Plan-first for nontrivial changes. No commit/push
  without explicit user approval.
