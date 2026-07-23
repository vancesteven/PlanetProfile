# Machine B handoff — Europa Clipper v5 + v6 SBI campaigns

**Author:** Machine A (Fable), 2026-07-22, documenting the campaigns codex staged.
**For:** Machine B (heavy compute; `PPcl` env).
**Scope:** two independent SBI training campaigns. Both are config + driver-script
+ reviewer-sign-off complete on Machine A; Machine B runs the compute and returns
artifacts + machine-readable gate reports. **Do not deploy. Do not tune thresholds
after seeing a failure — stop and surface evidence to Machine A.**

Pull the exact `origin/genai` commit Machine A names, record `PPcl` package
versions + all seeds, and note that the reference MCMC pickles are NOT committed
(see each campaign): either Machine A commits them first, or Machine B regenerates
via the `*_reference_mcmc.py` driver and returns them.

---

## Campaign v5 — thick-ice geodesy, three-arm paired ablation

**Design (reviewer PASS-WITH-CAVEATS 2026-07-20).** Reparameterizes v4's Tb as a
sampled ice thickness `D_iceIh_km` (prior truncated-gaussian μ29 σ10 [5,60]) plus
salinity, over the same 2D cache and the same 21 observables; adds a per-sample
hydrostatic-reference C/MR² for the dual-distribution GUI overlay. Three arms are
**column-sliced from ONE 1M dataset** (paired: identical θ draws + noise), so
posterior differences isolate the dropped channels' information:

- `europa_clipper_v5_geodesy_11D.json`      — baseline, 21 obs
- `europa_clipper_v5_noinduction_7obs.json` — drops the 14 `Bind_*`
- `europa_clipper_v5_nok2_17obs.json`       — drops Re/Im k2 + Re/Im h2

**Drivers** (`plans/scripts/`): `v5_fiducial_recompute.py` (fiducial gravity),
`v5_reference_mcmc.py`, `v5_gen_dataset.py` (1M baseline + name-based slicing with
the two reviewer safeguards: assert sliced obs-name lists match each ablation
config exactly and in order; report the differential drop rate from the shared
full-21 finiteness mask), `v5_train_all.py` (3 nsf artifacts), `v5_run_gates.py`
(gates for all three arms + the ablation comparison).

**Reference MCMC:** `europa_clipper_v5_reference_result.pkl` (4441 samples, ESS
4125, r̂ 1.0) was produced on Machine A but is NOT committed. Confirm with Machine
A whether to reuse a committed copy or regenerate.

**FOUR MANDATORY entries in every v5 gate report** (reviewer-required,
characterize-and-report; findings in `docs_ai/plans/v5_thick_ice_pull_findings.md`):

1. **k2 model–data tension.** Min reachable Re_k2 = 0.241 exceeds the 0.23 target
   by +0.72σ; Re_k2 decreases monotonically with ice thickness; THIS (not a real
   thick-ice measurement) drives the D pull. State the posterior D upper edge is
   prior-truncation-regularized, not data-localized (profile likelihood still
   rising at the 60 km wall: maxLL −12.4→−1.2 over D 37→59 km).
2. **Prior-truncation sensitivity.** The 60 km upper truncation is load-bearing.
   Report the flat-prior reweight (likelihood-only median 54.4 km, 97.5% 58.8 km)
   and that the reported D median (50.8–51.5 km) is a prior-mean-vs-likelihood-tail
   compromise that shifts up under looser truncation. Cheap — reweight existing
   samples, no new compute.
3. **k2 normalization confirmation (one line).** Cross-check the k2 forward map +
   the 0.23 target share the same normalization/definition (Re/Im, sign,
   unnormalized-Love) against a known Europa k2 baseline in the PP test suite.
   Pre-existing forward-model property, NOT a v5 regression — but flag explicitly
   that a future normalization mismatch would reinterpret the whole thick-ice pull.
4. **Degeneracy deliverable.** Report corr(D, w) = +0.863 TOGETHER WITH recovered
   corr(Tb_derived, w) = −0.935, so the v4 ridge is visibly preserved and the sign
   flip reads as a reparameterization artifact, not a new degeneracy. Do NOT
   present +0.863 alone as a new degeneracy.

Already verified sound (no action): cache monotonicity + D_max(w)≈96–101 km at
posterior salinities (support wall is the prior, not the cache); sampled-vs-derived
D consistency <1e-6 km; ridge preservation; convergence.

---

## Campaign v6 — free-gravity (open-interpretation) redesign

**Design (user-directed 2026-07-21/22; reviewer MODERATE).** A broad, open
"Europa-Clipper-alone" artifact; literature priors (Levin 2025 ice thickness;
Vance/Galileo induction) are applied as **separable inference-time reweights on the
one artifact**, not baked into training. Two decisive changes vs v5:

1. **Geodesy double-count eliminated.** v5 imposed CMR²=[0.3547,0.0024] as a tight
   independent observable *alongside* tight C20/C22 — but 0.3547 was itself derived
   FROM C22 via the hydrostatic Radau relation (GC21 Table 3), using the same
   relation twice and biasing the posterior toward hydrostaticity (the quantity the
   campaign measures). v6 **removes CMR² from the observables** entirely.
2. **Real, unconstrained gravity.** C20/C22 set to **GC21 Table 2 SOL-A
   (unconstrained)**: C20 = −4.3759e-4 ± 7.747e-5 (J2 = 437.59±77.47 ×1e-6),
   C22 = 1.3862e-4 ± 2.44e-6, **unnormalized at R_ref = 1565 km** (Table 2 is
   already unnormalized — NO √5/√(10/24) conversion; reviewer block-level fix,
   verified vs `papers/gomescasajus2021updated.pdf` Table 2). Measured ratio
   −C20/C22 = 3.157. `D_iceIh_km` prior widened to UNIFORM[5,80] km.

**Hybrid conditioning (state in the report):** v6 conditions GRAVITY on real
Galileo-era data (SOL-A) while k2/h2 + the 14 Bind channels remain the v5 SYNTHETIC
forecast (Mazarico Re_k2=0.23; Bind from the Tb 264.5 K / w 35.165 ppt / R_core
534.67 km fiducial). It is a real-gravity + forecast-tidal/induction hybrid, NOT a
self-consistent hydrostatic null — verified mutually consistent (fiducial
C22_h = 1.3775e-4 sits 0.36σ from SOL-A C22; fiducial C/MR² ≈ 0.355 ≈ SOL-A
hydrostatic reference).

**Configs** (`europa_clipper_v6_freegrav_*`): `11D` (baseline, 20 obs — C20/C22 +
k2/h2 + 14 Bind), `noinduction_6obs`, `nok2_16obs`. Built by
`plans/scripts/v6_build_configs.py`; guard test `tests/v6_config_test.py` (7/7,
re-run green on Machine A 2026-07-22). Drivers mirror v5: `v6_gen_dataset.py`,
`v6_reference_mcmc.py`, `v6_train_all.py`, `v6_run_gates.py` (docstrings still name
v5 — the code is v6-parameterized; confirm the config paths before running).

**Reference MCMC:** not committed — Machine A confirms reuse vs regenerate.

**v6 report must state:** the removed-CMR² rationale (no hydrostatic double-count);
the SOL-A unnormalized-at-1565 km provenance; the hybrid real-gravity/forecast-tidal
status; and — since gravity is now agnostic (wide dC20_nh/dC22_nh) — that C20/C22
carry ~zero interior C/MR² information, so the interior is constrained by k2/h2 +
induction (the honest v4/v6 physics trade).

---

## Execution order & returns

1. Confirm reference-MCMC provenance with Machine A (commit vs regenerate) for each
   campaign before gating.
2. v5: fiducial → (reference) → 1M dataset + slices → train 3 → gates (with the 4
   mandatory entries) → return.
3. v6: (reference) → 1M dataset + slices → train 3 → gates (with the v6 report
   items) → return.
4. Each return: exact commit, env/package versions, commands, wall time, per-arm
   artifacts + machine-readable gate JSONs, first full traceback per regression,
   rejection statistics, seeds, and a concise verdict. `PlanetProfile/Test/` or
   `.gitignore` additions need explicit user permission before commit.

Neither campaign is authorized to deploy. Integration into GUI slots + INDEX +
ratification is Machine A's step after the reports return.
