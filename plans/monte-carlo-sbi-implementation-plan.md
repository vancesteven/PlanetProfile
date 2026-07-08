# Approved Plan: Monte-Carlo → SBI Amortized Inference Pipeline

Approved by user 2026-07-04. Supersedes the SBI portions of `plans/titan-mcmc-sbi-plan.md`
(deleted from working tree; retrievable via `git show d9f230f5:plans/titan-mcmc-sbi-plan.md`).

## Approved decisions

1. **Environment**: conda env `PP` on this machine (PPcl exists only on user's other computer).
   Install `torch`, `pocomc`, `sbi>=0.23,<0.27` into `PP`.
2. **Scientific fixes to `generate_sbi_dataset`** (all approved, opt-in kwargs, default path byte-identical):
   - `apply_support_guard=True` for SBI: apply `_check_no_ocean` so training-set support
     matches the MCMC posterior support (no-ocean configs).
   - `imag_convention='abs'`: store |Im k2| in x (matches likelihood and GUI convention).
   - `drop_nonfinite=True`: reject TidalPy-failure NaN rows; report
     `n_requested / n_rejected_nonfinite / n_rejected_support / rejection_rate`.
     Stop and surface if rejection_rate > 20%.
   - `seed=` + `provenance=`: seeded prior draws; .npz gains param_bounds, param_units,
     config_hash, git_sha, seed, rejection stats, imag_convention (additive keys only).
3. **era-LLM meta-optimization (roadmap Phase 4): DEFERRED** to a separate plan after SBI validation.
4. **Validation thresholds RATIFIED as pass/fail gates**:
   - SBC rank uniformity not rejected at α=0.05 per parameter (≥200 held-out pairs).
   - MCMC↔SBI cross-check on Test50 Titan no-ocean, same observables {Re_k2: 0.608, Im_k2: 0.135}:
     posterior means agree within max(0.25·σ_MCMC, MCMC MC-error);
     σ_SBI/σ_MCMC ∈ [0.7, 1.4]; 1-D marginal KS not rejected at α=0.01.
   - Limiting behavior: |Im k2| ↑ ⇒ log10_eta_Ih posterior median ↓ (monotone); posterior inside prior box.
   - Any gate failure: stop and surface, do not tune to pass.

## Architecture (from design review)

- **`PlanetProfile/Inference/sbi_runner.py` (new)**: `SBIRunner(config: InferenceConfig)`,
  `.run(progress_callback) -> InferenceResult` per `inference_core.py:360-371` contract.
  Single-round amortized NPE, MAF density estimator (`'maf'` default, `'nsf'` optional),
  identity embedding (x ∈ R²–R³). Lazy torch/sbi imports. Prior built directly from
  `config.param_space` (uniform-only first; `NotImplementedError` for normal/log-uniform).
  Training in raw theta space (log10_* params already log-space). Flow log-prob goes in
  `metadata['flow_log_prob']`; `InferenceResult.log_likelihoods` holds recomputed Gaussian
  χ² log-likelihood (posterior density ≠ likelihood). `weights=None`.
- **Artifact `.pt`** (self-describing, `sbi_artifacts/` naming per INDEX.md):
  posterior/estimator, prior_spec, param_names/units/bounds, obs_names, imag_convention,
  theta/x normalization constants, config_hash, git_sha, seed, n_train_effective,
  rejection_stats, sbi_version, torch_version, created_utc, schema_version=1.
  `load_artifact` warns loudly on sbi/torch version mismatch.
- **Dataset**: 10k min, 50k production target; measure wall-clock on 200 samples first.
- **GUI (`PlanetProfileApp/pages/AmortizedInference.py`)**: artifact resolver dict,
  enable button only when `.pt` exists, `@st.cache_resource` keyed on path,
  x_obs = {'Re_k2': re, 'abs_Im_k2': abs(im)}, corner/marginal plot + summary table.
  Untrained slots stay disabled.

## Phases

| Phase | Task | Files | Agent |
|---|---|---|---|
| A | Opt-in kwargs + provenance in `generate_sbi_dataset` | `mcmc_runner.py` (additive only) | sonnet |
| B | `sbi_runner.py` | new file | sonnet |
| C | `tests/sbi_runner_test.py`: prior-consistency, shapes/schema, abs-Im, NaN drop, support guard, seed repro, artifact roundtrip, tiny n=200 train smoke (skipUnless sbi) | new file | sonnet |
| D | GUI wiring | `AmortizedInference.py` | sonnet |
| E | SBC + MCMC-cross-check driver (gated `PP_RUN_SBI_HEAVY=1`), Test50 regression guard | new driver script | opus |
| — | Sign-off: prior mapping, conventions, guard, thresholds, benchmarks | — | opus/orchestrator |

## Do NOT touch

`inference_core.py`, `parameter_registry.py`, `forward_models.py`, `structure_cache.py`,
`structure_derivation.py`, `MCMCRunner._build_prior` / likelihood / `run`, PPTest*/Test50 modules
and cached grids, Perple_X dirs, `.claudeignore`/`.graphifyignore`, manuscripts.

## Known landmines (from code survey)

- `sbi_artifacts/INDEX.md` line ref stale: `generate_sbi_dataset` is at `mcmc_runner.py:959`, not 722.
- `parameter_registry.py` preset `andrade_titan_noocean_8D` observables (Re_k2, abs_Im_k2, CMR2)
  ≠ config JSON observables (Re_k2, Im_k2). Resolve before generating training data.
- `generate_sbi_dataset` silently NaN-fills unrecognized observable names (mcmc_runner.py:1001).
- theta excludes param_groups members and fixed_params; expansion logic (`_expand_theta`)
  lives only inside MCMCRunner — SBIRunner must reuse it, not duplicate.
- sbi ≥0.23 renamed SNPE→NPE; read installed API, never code imports from memory.
- No benchmark numbers exist yet for Test50 heat-partitioning regression — must be generated
  and blessed by user before becoming a gate.

## Validation (Phase E)

Driver: `PlanetProfile/Inference/validate_sbi.py` (manually invoked; NOT run by `tests/`).
Every subcommand takes `--seed` (default 42) and `--output-dir`, writes a JSON report plus
plot PNGs, exits `0` on pass and `2` on gate failure (`1` for setup/IO errors, e.g. a missing
structure cache). Every report embeds config_hash, artifact metadata, seeds, sbi/torch/scipy/
numpy versions, git SHA, and UTC timestamps. The gates below are the RATIFIED thresholds
(decision 4) — hard pass/fail, never tuned to pass.

Verified installed sbi 0.26.1 API (read from source, not memory):
`sbi.diagnostics.run_sbc(thetas, xs, posterior, num_posterior_samples, reduce_fns='marginals')
-> (ranks, dap_samples)`; `sbi.diagnostics.check_sbc(ranks, prior_samples, dap_samples,
num_posterior_samples) -> {'ks_pvals', 'c2st_ranks', 'c2st_dap'}`;
`sbi.analysis.sbc_rank_plot(ranks, num_posterior_samples, parameter_labels, plot_type='cdf')
-> (fig, ax)`.

1. **`sbc`** — simulation-based calibration on ≥200 held-out (θ,x) pairs. Pairs come from a
   pre-generated `--dataset .npz` (no forward model) or are generated from `--config` (needs
   the structure cache; fails fast naming the missing file + the `build_phase_c1_cache.py`
   route). **Gate:** per-parameter rank-uniformity KS p ≥ 0.05. Output: `sbc_report.json`,
   `sbc_rank_cdf.png`.
2. **`crosscheck`** — SBI artifact vs an MCMC `InferenceResult` pickle conditioned on the same
   observables. MCMC importance weights are handled: weighted mean/σ, ESS = (Σw)²/Σw² (mirrors
   `MCMCRunner`'s true-ESS), and the weighted posterior is resampled to ⌊ESS⌋ unweighted draws
   for the KS test (so it is not over-powered by the raw sample count). **Gates (per param):**
   |mean_SBI − mean_MCMC| ≤ max(0.25·σ_MCMC, σ_MCMC/√ESS); σ_SBI/σ_MCMC ∈ [0.7, 1.4]; two-sample
   KS p ≥ 0.01. Output: `crosscheck_report.json`, `crosscheck_marginals.png`.
3. **`limits`** — condition the flow on a grid of |Im k2| (default 0.05..0.30, other channels
   fixed via `--re-k2`/`--fixed-obs`). **Gates:** `log10_eta_Ih` posterior median strictly
   monotone decreasing (≤1 tie); 100% of drawn samples inside the prior box. Output:
   `limits_report.json`, `limits_median_trend.png`.

**`selftest`** — `conda run -n PP python -m PlanetProfile.Inference.validate_sbi selftest --seed 42`
builds a toy BoxUniform + linear-Gaussian setup, trains a tiny MAF, and runs all three checks
end-to-end (including a gate-bite negative control proving a shifted mean FAILs). Verified PASS
at seed 42.

Heavy-run gating: any future pytest wrapper must be gated behind `PP_RUN_SBI_HEAVY=1`; do NOT
add slow SBC/cross-check runs to `tests/`.
