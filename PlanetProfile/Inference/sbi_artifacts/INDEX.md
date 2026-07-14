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
| titan_andrade_noocean_posterior.pt | test50_titan_noocean_andrade_8D.json | 629afbd55a4f0ce5 | bf7c938e | train 42 / data 44 / noise 4444 | 999,816 (nsf) | ALL GREEN within domain (amended rules, 2026-07-11) | 2026-07-12 (user-ratified; artifact pushed git 4d32809b, pre-deploy assertions green) |
| europa_seawater_andrade_posterior_1m.pt | europa_seawater_andrade_7D.json | a09396bcb0d0eff5 | 3d865dc1 | train 42 / data 44 / noise 4444 | 831,750 (nsf, synodic-only support cut) | 3 gates FAIL but ALL are flow-calibration slack, NOT physics; physics-discrimination (limits W1 anchors) GREEN all 6. Deployed under domain guard + default Tb truncation (see conditions below) | 2026-07-13 (user-directed deploy w/o Machine A; scientific-reviewer DEPLOYABLE-after-remedy; GUI-verified via AppTest) |

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

**Cross-version sampling validation (2026-07-13, Machine A):** the Europa artifact
was trained on torch 2.11 (Machine B) and loads on torch 2.8 (Machine A) with the
loud version-mismatch warning by design. VALIDATED: re-running the crosscheck gate
on Machine A reproduces Machine B's committed report to four decimal places on every
per-parameter statistic (D, dmean, verdicts; Tb-only edge-smear fail included) —
flow deserialization is bit-consistent across these torch versions at the gate seed.
The warning stays (correct default caution for future version pairs); for this
artifact/version pair the results are trustworthy.

## Candidate artifacts (NOT deployed — evidence committed for cross-machine verification)

| Artifact file | Source config | config_hash | seed | n_train | Gates | Status |
|---|---|---|---|---|---|---|
| titan_diff_noocean_andrade_test52_10D_v2.pt | test52_titan_noocean_andrade_10D.json | 2bf1f7b2d1708e28 | train 42 / data 45 / noise 4545 | 877,883 (nsf, 10xIQR x-filter) | SBC FAIL (eta_V p=.016, 9/10 pass); crosscheck shape FAIL 3 params (zeta/eta_III/Tb, all location tests pass); limits W1 PASS Im<=0.25, FAIL Im=0.30 | **implemented, unverified** — Test52 10D (Titan no-ocean + CMR2). NOT deployed. Gate reports in `validation_reports/test52_v2/`. Reviewer (opus 2026-07-11): deployable with |Im k2|<=0.25 restriction + lower-tail shape caveats on zeta/eta_III/Tb. Awaiting Machine A GUI guard + deployment decision. |

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
