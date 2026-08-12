# Amortized SBI for ocean-world interiors — methodology & reproduction

**Status:** active cross-body methodology (updated 2026-08-12; Europa v5/v6
and Titan free-gravity composition campaigns).
**Scope:** moon-agnostic recipe for building an amortized posterior over interior
parameters from geodesy + induction observables, then imposing literature
constraints at inference time. The worked example is Europa; the same recipe
transfers to other icy moons by swapping the config + structure cache.

This complements the MCMC-tier docs
([PlanetProfile/Inference/README.md](../PlanetProfile/Inference/README.md),
[MCMC_INFERENCE_GUIDE.md](../MCMC_INFERENCE_GUIDE.md)); it does not repeat them.

---

## 1. The core idea: broad training prior + inference-time reweighting

Train ONE amortized flow under a **broad, deliberately uninformative ("open
interpretation") prior** — the posterior Clipper *alone* could support,
independent of any single literature interpretation. Then impose a literature
constraint (e.g. Levin et al. 2025 ice thickness, or a tight Galileo/Vance
induction interpretation) at inference time by **importance reweighting**:

```
p_informed(θ | x) ∝ p_flow(θ | x) · [ π_lit(θ) / π_train(θ) ]
```

One artifact serves both the "what does Clipper say by itself?" question and any
number of literature-informed posteriors, with no retraining.

**Requirements for the reweight to be valid**
- The training prior must **densely cover** the literature mass. A reweight can
  only *remove* probability the flow already placed; it cannot create density
  where the flow saw no training samples.
- Cleanest when π_train is **uniform** over each parameter's physical range, so
  π_lit/π_train ∝ π_lit (no ratio blow-up in sparse tails). A narrow informative
  training prior *starves* the tail the literature later wants — the failure mode
  that motivated the v5 pivot.

**Transfer rule for a new moon:** put uniform priors on the parameters you have
genuine open belief about; reserve informative constraints for the reweight step.

---

## 2. Parameterize on the quantity you have belief about

Sample the interior parameter the literature constrains **directly**, not a proxy.

- Europa v5 samples `D_iceIh_km` (total ice-Ih shell thickness) `uniform[5, 80]`
  km, and *derives* `Tb_K` by inverting the cached `D_iceIh(Tb, w)` field per
  salinity column. A uniform box on Tb would be an implicit, strongly
  salinity-dependent thickness prior (freezing-point depression), contaminating
  the Tb↔salinity degeneracy — a headline deliverable.
- The inversion uses the **same** bilinear (Tb, log w) operator every downstream
  consumer uses, so training support == reference-MCMC support byte-for-byte
  (reviewer-binding). Draws outside the achievable band reject (never clip).

**Definition hazard — match the literature's definition to the model's.**
Levin et al. 2025 report `29±10 km` as the **conductive-lid** thickness; PP's
`D_iceIh_km` is the **total** shell (PP models solid-state convection, so total ⊇
conductive). A naive `N(29,10)` reweight on the total would forbid convection and
is physically wrong. The literature constraint enters as a **convection-aware
conductive→total construction** at reweight time, not as a prior on the sampled
total. Always reconcile the two definitions before reweighting.

---

## 3. Induction support: per-frequency, degree-based phase caps

Induction constraints differ by frequency. A single global Ae bound is wrong.

- Europa (Vance et al. 2021, Fig 6 / Table 2): the tight synodic interpretation
  (`|Ae|>0.70`, phase `<30°`) does **not** hold at the longer-period orbital
  band, where low-salinity thick-ice models reach `|Ae|≈0.37` at `~65°`. So bounds
  are **per label**: synodic & synodic-2nd `|Ae|>0.70 / phase<30°`; orbital
  `|Ae|>0.30 / phase<70°`.
- Cap **phase in degrees**, not via an `|Im Ae|` proxy: the proxy mismaps to phase
  at low `|Ae|`; a degree cap on `|angle(Ae)|` is amplitude-independent and
  correct. (Convention: `Ae^{-iφ}`, PP stores `Im<0` for a delay.)
- The identical predicate is applied at BOTH the MCMC likelihood site and the SBI
  support guard, so the two supports agree (regression:
  `PlanetProfile/Inference/tests/test_bind_channels.py`).

**Do NOT** try to encode a phase constraint by inflating induction σ — the
excitation amplitudes differ >20× across bands (Europa `|Be_syn|≈225` vs
`|Be_orb|≈10` nT), so uniform σ de-weights bands in the wrong proportion.

---

## 4. Paired-slice ablations (what each observable set buys)

To attribute constraining power to observable groups, generate ONE baseline
dataset and **column-slice by observable name** into ablation arms (Europa v5:
baseline 21 obs, noinduction 7, nok2 17). All arms then share identical θ draws
AND noise realizations, so posterior differences isolate the dropped observables'
information — a *paired* ablation.

Consequence: **all arms must declare identical support bounds** (induction, etc.).
The shared θ are drawn under the baseline guard; an arm with different bounds would
compare a guard-bounded flow to guard-free fresh SBC pairs and fail spuriously.
Arms are correlated — do not treat as independent for joint statements.

---

## 5. Gate suite (never tuned to pass)

Automated gates, per artifact, via `PlanetProfile.Inference.validate_sbi`:
- **SBC** — simulation-based calibration on fresh held-out pairs. The pairs are
  generated through the SAME support-guarded forward model as training (so the
  held-out support matches the trained support); sampling uses
  `reject_outside_prior=False` to avoid the 0%-acceptance hang at far-tail x.
- **limits** — limiting-behavior monotonicity + prior containment. Every
  non-swept observable MUST be pinned with `--fixed-obs` (the swept channel is
  the auto-detected |Im k2|); omitting it raises. Skipped for the nok2 arm,
  which has no k2 channel to sweep.
- **crosscheck** (baseline only) — SBI vs the reference MCMC on the same obs;
  ratification-blocking. Ablations have no reference MCMC.
- **posterior-predictive / pushforward check** — standard equipment for any
  artifact whose sampled parameters enter a nonlinear observable. Report the
  noiseless observable medians in a four-way table:

  | channel | observed | prior predictive | reference-MCMC posterior predictive | SBI posterior predictive |
  |---|---:|---:|---:|---:|
  | each conditioned observable | value ± `sigma_obs` | median | median | median |

  Flag a channel when the SBI and reference-MCMC posterior-predictive medians
  differ by more than `0.5 * sigma_obs`. This separates a flow deficiency
  (SBI vs MCMC) from shared model-data tension (MCMC vs observation). For the
  importance-corrected pipeline below, Re/Im k2 also receive weighted KS/W1
  full-distribution checks; the median flag alone is not a width test.

Reported diagnostics (NOT automated pass/fail gates — computed and recorded by
`v5_run_gates.py`, interpreted by hand):
- **Ablation / degeneracy comparison** — per arm, condition the flow on that
  arm's fiducial observables and report the `D_iceIh` posterior median and
  `corr(D_iceIh, log10_w)` (Europa: the D↔w ridge). Written to
  `v5_gate_summary.json → ablation_comparison`. This attributes constraining
  power to observable groups; it is not thresholded.
- **Coverage caveats** — any parameter region unreachable at some salinity
  (cache `D_max(w)` envelope) is reported in the coverage note, not silently
  truncated.

Gates are acceptance criteria. If a gate fails, fix the model or the design —
never loosen the gate. Diagnostics inform interpretation; they are not tuned.

---

## 6. Track 1: exact-likelihood importance correction

The deployed flow can be treated as an importance proposal for the exact
posterior targeted by the reference MCMC. For draws
`theta_i ~ q(theta | x_obs)`, the correction uses

```
log w_i = log p(x_obs | theta_i) - log q(theta_i | x_obs)
```

because the common uniform BoxUniform prior cancels. The two inputs already
come from the conditioning path: full-N `InferenceResult.log_likelihoods` and
`metadata['flow_log_prob']`. Support-guard rejections use the MCMC likelihood's
`-1e30` sentinel and are masked before weight arithmetic. The implementation is
[`PlanetProfile/Inference/is_correction.py`](../PlanetProfile/Inference/is_correction.py);
the preregistered validation and deployment decisions live in
[`plans/active/tidal-sector-remedy-plan.md`](../plans/active/tidal-sector-remedy-plan.md).

The binding reliability set is:

- Pareto-k is primary: `k <= 0.5` is clean, `0.5 < k <= 0.7` uses PSIS
  smoothing, and `k > 0.7` fails to the uncorrected/quarantined fallback.
- Absolute ESS must be at least `1000` at the N actually run. `ESS/N` remains a
  reported cost and sample-sizing diagnostic, with
  `N_required = 1000 / (ESS/N)`; it is not a reliability gate.
- The maximum normalized weight must be at most `0.01`.
- Reverse coverage is blocking: less than `0.01` of reference-MCMC mass may
  fall below the `0.001` quantile of the flow's own log-density.

Additional guards retain the full sample alignment and target definition:
sentinel-rejected fraction must not exceed `0.5`; corrected mass outside the k2
training-support box must stay below `1e-3`; and any ocean/no-ocean branch above
`0.05` probability must have branch ESS at least `50`. Downstream summaries
consume the corrected weights and report ESS rather than the resampled count.

Passing these checks validates the corrected pipeline against the MCMC target;
it does not remove model-data tension shared by both samplers. Per-composition
quarantine remains until the full preregistered crosscheck, pushforward,
coverage, branch-mass, stability, and corrected-SBC gates pass.

---

## 7. Reproduction — Europa v5 (exact commands)

Env for every step: `mamba run -n PPcl env PYTHONPATH=. \
NUMBA_CACHE_DIR=/tmp/pp_numba_cache KMP_DUPLICATE_LIB_OK=TRUE`.
Build large outputs in `/tmp`, copy to Dropbox. Seeds: data=57, noise=5757,
train=51.

```bash
# Configs (checked in): PlanetProfile/Inference/configs/europa_clipper_v5_*.json
#   geodesy_11D (baseline 21obs), noinduction_7obs, nok2_17obs

# (Cache prerequisite — build once; already present for v5)
#   python -m PlanetProfile.Inference.build_phase_c1_cache --config <baseline> \
#     --n-grid <N> --n-wgrid 16 --force
# -> Test/mcmc_results/Europa/Test52_seawater_v5/europa_seawater_structure_grid_v5_2d.pkl

# D. Reference MCMC (baseline, pocoMC) — reproducible driver
python plans/scripts/v5_reference_mcmc.py
# -> Test/.../europa_clipper_v5_reference_result.pkl + v5_reference_manifest.json

# E. One 1M baseline dataset + name-sliced ablation arms (paired)
python plans/scripts/v5_gen_dataset.py --n 1000000
# -> /tmp/v5_build/datasets/v5_{baseline,noinduction,nok2}_1m.npz + manifest

# F. Train three nsf artifacts (seed 51)
python plans/scripts/v5_train_all.py
# -> PlanetProfile/Inference/sbi_artifacts/europa_clipper_v5_*_posterior_1m.pt

# G. Gates + ablation comparison
python plans/scripts/v5_run_gates.py
# -> validation_reports/europa_clipper_v5_{baseline,noinduction,nok2}_1m/
```

**Provenance to record for each artifact (ratification):** config_hash, obs_names,
seeds, structure_cache sha256, Ae-sidecar sha256, git_sha. Where each field lives:
- SBI artifact (`.pt`) embeds `config_hash`, `obs_names`, `seed`, `git_sha`
  (via `save_artifact`) — but NOT the cache/sidecar hashes.
- `v5_train_all.py`'s train manifest records `config_hash`, `structure_cache_sha256`,
  and `ae_sidecar_sha256` keyed per arm — this is what pins an artifact to its cache.
- The reference-MCMC manifest records both cache and sidecar sha256, config, seed,
  log_Z, and git_sha.
So an artifact traces to its cache via the train manifest (same directory), not via
the `.pt` alone.

> Artifact hashes for the ratified v5 run are filled in at ratification in
> `PlanetProfile/Inference/sbi_artifacts/INDEX.md`.

---

## 8. Reproducibility hygiene

- Every artifact must be regenerable from the config + committed scripts alone.
  No ad-hoc runs without a committed command.
- Keep docs and code concise: recipe + exact commands, not prose.
- Sweep for stale files periodically: check an artifact's provenance
  (hash/date/config_hash) against the current config before trusting it; archive
  files whose assumptions have been superseded (e.g. the pre-pivot v5 artifacts
  archived 2026-07-21 when the prior changed from truncated-Gaussian to uniform).
- The scientific-reviewer verifies procedures are documented for reproducibility;
  a repro gap is a review finding.
