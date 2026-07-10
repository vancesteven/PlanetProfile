# SBI Artifact Index

Last audited: 2026-07-08

## Naming convention

GUI slot filenames (hardcoded in `PlanetProfileApp/pages/AmortizedInference.py::_ARTIFACTS`):

```
sbi_artifacts/
  titan_andrade_noocean_posterior.pt     # Titan (Andrade, No-Ocean)  — Test50 8D
  titan_maxwell_ocean_posterior.pt       # Titan (Maxwell, Ocean)
  europa_andrade_thinshell_posterior.pt  # Europa (Andrade, Thin Shell)
```

`SBIRunner.save_artifact` defaults to `<bodyname>_<config_hash>_posterior.pt`; deploying a
validated artifact to a GUI slot means copying it to the slot filename above. Artifacts
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

| Slot file | Source config | config_hash | git SHA | seed | n_train | Gate reports | Deployed |
|---|---|---|---|---|---|---|---|
| titan_andrade_noocean_posterior.pt (**PROVISIONAL**) | test50_titan_noocean_andrade_8D.json | 629afbd55a4f0ce5 | 278c3bea | train 42 / data 43 / noise 4343 | 499,915 (nsf) | see below | 2026-07-09 |

**PROVISIONAL status (user-approved for GUI preview, 2026-07-09):** SBC PASSES
decisively (1000 held-out pairs, min KS p=0.252). Crosscheck mean+sigma gates all pass
(worst dmean 0.18 of tol 0.40); KS residual on 4 params (p 6e-6..3e-3) pending the 1M
retrain + materiality ratification. Anchor-limits passes at |Im k2| = 0.05/0.10/0.15/0.30,
fails 0.20/0.25 (0.25 is the bimodal-median regime — see
plans/HANDOFF-2026-07-09-test50-sbi-validation.md). Posteriors conditioned near the
observed point (Re k2 0.608, |Im k2| 0.135) are well-validated; treat conditioning at
|Im k2| >= 0.20 as unvalidated extrapolation until the anchor gate passes.

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
- GUI: `PlanetProfileApp/pages/AmortizedInference.py` — slot button enables automatically
  when the slot file exists; loud failure on non-`'abs'` artifacts.
