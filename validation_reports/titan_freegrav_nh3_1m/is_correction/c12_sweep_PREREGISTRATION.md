# C12 amortized sweep — preregistration (before any result is read)

Spec source: `plans/MACHINE-B-HANDOFF.md:383` (paraphrase) and
`plans/active/tidal-sector-remedy-plan.md:140-141,244-245` (canonical wording):
"C12 (amortization across x): >=200 prior-predictive x + 8 axis endpoints,
Pareto-k <= 0.7 for >=95%." Task 5 also requires reporting the ESS/N
distribution across x (`N_required = 1000/(ESS/N)`) and one sanity line on
the 0.9326-to-four-decimals branch coincidence.

No precedent driver exists for this exact construction (checked
`plans/scripts/*.py`, `PlanetProfile/Inference/*.py`,
`validation_reports/*/is_correction/`, `validation_reports/*/pilot/`).
This document fixes every design choice BEFORE the sweep runs, per this
campaign's "never tune after seeing a failure" discipline.

## 1. What draws the 200 "prior-predictive x"

theta ~ effective prior (`SBIRunner.generate_training_set`, which applies
the support guard so draws come from the trained-on prior, matching the C15
discipline used for SBC) -> push through the forward model -> x. This is
the standard prior-predictive convention already used by
`plans/scripts/ppc_pushforward_check.py` and
`plans/scripts/titanG_ppc_interior_check.py`.

**obs_noise=True** (not False). Rationale: C12 tests whether the
correction is reliable at x's the flow will actually be asked to condition
on, and the flow was trained on noised x. A noise-free x is off the training
manifold in a different, uninteresting way. This is a preregistered choice,
not something to revisit after seeing the k-distribution.

**N per point: 5,000** (not the fiducial's 20,000). Rationale: the fiducial
N=20k budget (~6 min wall/point) times 208 points is ~21h serial /
~1.75h at 12-way parallel (matches the "1h on 16 cores" estimate in
MACHINE-B-HANDOFF.md:312 to within core-count scaling). Pareto-k is the
gate (not absolute ESS, which is fiducial-only per the 2026-08-11 ruling
striking the ESS/N floor) and is well-estimated at N=5k. Absolute ESS is
still recorded per point for the ESS/N distribution deliverable, but is
NOT gated off-fiducial (report-only, matching the ruling's own treatment
of ESS/N). This keeps runtime to <30 min at 12-way parallelism.

## 2. The 8 axis endpoints

4 observables (`C20, C22, Re_k2, Im_k2` per
`PlanetProfile/Inference/configs/test54_titan_nh3_freegrav.json`) x 2
endpoints (mean +/- 2*sigma_train, using the deployed artifact's `x_norm`
mean/std — the repo's ratified in-support definition,
`LIMITS_INSUPPORT_Z = 2.0` in `validate_sbi.py:177`) = 8. One axis swept at
a time, the other three held at the config's fiducial centrals — the same
construction as `validate_sbi._run_limits_check`.

Im_k2's lower endpoint (mean - 2*sigma = -0.1696) is unphysical (the channel
is abs-folded, zero-bounded) and is clamped to 0.0 — same floor precedent
as `_build_default_sweep_grid` (`validate_sbi.py:1060`). This is a genuine
physical constraint, not a convenience clamp.

No joint corner point is added; 4x2=8 with no corner is the literal reading
and is what is gated. (A joint-corner point may be run and reported
separately, outside the 95% denominator, if time permits — not required.)

## 3. Denominator for the 95% pass rule

**Pooled: all 200+8=208 points**, i.e. at most floor(0.05*208)=10 points
may have Pareto-k > 0.7 (or NaN — see below) for the sweep to PASS.
Reported ALSO split as {prior-predictive-only, endpoints-only} so the
adjudicator can see whether failures concentrate at the domain edge (a
different finding than failures scattered through the bulk) even though
the gate itself is on the pooled set.

NaN Pareto-k (fewer than ~50 nonzero weights survive) counts as a FAILURE
for this rule, not an exclusion — excluding it would make the gate
gameable by degeneracy.

## 4. The config-hash / likelihood-mismatch trap (must not be gotten wrong)

`_x_obs_vector`/flow conditioning is steered by the `x_obs` argument to
`_condition_and_package`, but the LIKELIHOOD closure
(`MCMCRunner.log_likelihood_fn`, built once at `MCMCRunner.__init__` from
`config.observables`) is captured from the config at construction time and
is NOT re-derived from `x_obs`. Sweeping `x_obs` while leaving
`cfg.observables` at the fiducial would silently compute
`p(x_fiducial | theta) / q(theta | x_swept)` — a meaningless ratio that
would still produce plausible-looking Pareto-k numbers.

**Fix applied in the driver**: for each swept point, build a FRESH
`InferenceConfig` with `observables` centrals set to that point's x (sigmas
unchanged from the fiducial config) and construct a FRESH `SBIRunner`
against it, so `_get_mcmc_runner()` builds its likelihood closure from the
SAME x used for conditioning. The artifact is located and loaded via the
UNMUTATED fiducial config's hash first (so artifact matching is unaffected
by the sweep), then its `posterior`/`_train_info` are injected into the
per-point runner exactly as `is_correction_validate.py` does for the single
fiducial call. `assert_byte_identity`'s config-hash check is EXPECTED to
fail off-fiducial by construction (the observables differ) — this is
recorded explicitly per point as `"c2_hash_expected_mismatch": true`,
never silently skipped, and the OTHER three byte-identity checks (structure
cache SHA, param name alignment, imag_convention) still run and still gate.

## 5. Verification before the real sweep

Before running the 208-point sweep, the point-builder is run ONCE at the
exact fiducial x and its output (Pareto-k, ESS, ESS/N, pushforward
medians) is required to reproduce
`validation_reports/titan_freegrav_nh3_1m/is_correction/is_validation_nh3.json`
to the precision the RNG allows (same seed => identical draws => identical
diagnostics). If it does not reproduce, the sweep does not run until it
does.

## 6. Deliverables

- Per-point: x, Pareto-k, ESS, ESS/N, N_required=1000/(ESS/N), verdict,
  fail_reasons, ocean-branch prob_flow/prob_corrected (for the 0.9326
  sanity line).
- Pooled pass/fail per §3, plus the prior-predictive-only and
  endpoints-only sub-splits.
- ESS/N distribution summary (min/median/max) across all 208 points.
- One sanity line: flow-vs-corrected ocean-branch histogram comparison,
  citing why 0.9326 was stable to 4 decimals at the fiducial (already
  cleared per §0.18 Phase-1 item 2 — this sweep reports whether that
  stability holds off-fiducial too, not re-litigating the fiducial finding).
- This document, committed alongside the report, so the design record
  predates the result.

Not self-adjudicated: PASS/FAIL is reported per the rule above; the
scientific-reviewer interprets what a FAIL (if any) means for the NH3
quarantine-lift path. Gates are interpreted, never tuned to pass.
