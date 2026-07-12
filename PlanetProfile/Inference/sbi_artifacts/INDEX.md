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

| Slot file | Source config | config_hash | git SHA | seed | n_train | Gates | Deployed |
|---|---|---|---|---|---|---|---|
| titan_andrade_noocean_posterior.pt | test50_titan_noocean_andrade_8D.json | 629afbd55a4f0ce5 | bf7c938e | train 42 / data 44 / noise 4444 | 999,816 (nsf) | ALL GREEN within domain (amended rules, 2026-07-11) | 2026-07-12 (user-ratified; pending Machine B artifact push) |

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
- GUI: `PlanetProfileApp/pages/AmortizedInference.py` — slot button enables automatically
  when the slot file exists; loud failure on non-`'abs'` artifacts.
