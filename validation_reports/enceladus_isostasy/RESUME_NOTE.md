# Resume note — Enceladus isostasy campaign (Machine B)

_Written 2026-08-13 at a pause point. Working tree clean; all work committed._

## Where things stand

Machine A's §0.22 directive made **Enceladus the 1.0 foreground** (salt C16
parked post-1.0, NH3 quarantine-lift backgrounded). This session built the
frozen-config candidate and worked its physics blockers.

`PlanetProfile/Inference/configs/enceladus_cassini_isostasy_7D.json`
(filename says 7D; it is **5D** — two dead params dropped) is a
**frozen-config CANDIDATE, `not implemented`, NOT runnable and NOT ratified**.
Four scientific-reviewer passes; the last returned NEEDS-MORE-ANALYSIS.

Observables {C20, C22, C30, libration_deg}. Ocean branch samples
zb_km / log10_wOcean_ppt / rho_ice_kgm3 / compensation_C2 / libration_sys_frac;
frozen branch samples rho_rock_kgm3 / rho_ice_kgm3 / libration_sys_frac with
**zb derived by mass conservation**.

## THE ONE THING TO PICK UP FIRST

**B2' is ADJUDICATED (scientific-reviewer verdict BLOCK, 2026-08-14) and REMAINS
BLOCKING.** Read `b2prime_ADJUDICATED_drho_weighting.json` FIRST — it supersedes
both `b2prime_REOPENED_figure_coupling.json` and
`b2prime_hm_published_discriminant.json`.

**The shipped observed-figure treatment — the config's adopted "option A" — is
not a valid candidate.** `H22_obs_m` scales `(Bs - As)`, but the shell-base
interface carries weight `-rho_ice = -925` inside that difference while its
*physical* weight in `Ks` is `Delta_rho = rho_ocean - rho_ice` ~ 80-95 kg/m^3
(the `Bsp_Asp` term nearly cancels it). So scaling the difference applies an
implicit, structure-dependent, **sign-flipping** effective scale of **+0.33 to
-0.58** to that interface, while the docstring states interior interfaces stay
hydrostatic. Confirmed independently by Machine B: the identity
`Ks/(3 omega^2) = rho_ice*f_top + Delta_rho*f_base` holds to ~1e-14.

The defensible treatments are **hydrostatic** and **surface-only**, with
Delta_rho-consistent Eq.-12 bracketed by C2 between 25.99 and 27.34 km. Span at
the conditioned libration: 24.42-27.68 km over the mass-admissible ocean-thickness
subrange. **The fix is a defect repair, not a scaling guess** — which is what
this note previously said the next step must not be.

**Two previously-recorded headline numbers are ARTIFACTS of the same
mis-weighting — do not quote them:** the Eq.-12 "30.5 sigma, no zb matches"
result (corrected: **+1.66 sigma, zb = 25.99 km**), and the C2 sweep
"0.1005 -> 0.1824 deg" (corrected span across C2 in [0,1] is **1.6 sigma**,
~20x smaller). The C2 dependence is real but small, so MAJOR-4's "gravity =
compensation state, libration = sole thickness channel" split **largely
survives** after all.

**The published-answer test was run and did NOT discriminate** — all treatments
fit inside H&M's published 16-22 km libration-only band, whose +/-1 sigma
half-width (~2 km) exceeds the entire treatment spread. It did produce two
things worth having: PP's 2x2 solve is **algebraically identical to H&M Eq. 20**
(rel <= 2.2e-16, so the "not term-for-term identical" caveat was over-cautious
about the form; the COEFFICIENTS remain uncertified), and imposing H&M's own
published core parameters puts the answer inside their **joint** Table 2 band.

**Headline finding, RATIFIED:** the H&M-vs-Park shell-thickness separation is
essentially *entirely* the libration measurement revision (Thomas 0.120 -> Park
0.092) propagating through an unchanged forward model — **Delta = +5.6 km**
(19.5-19.7 -> 25.1-25.4 km at H&M's own published core parameters) against a
~5.5 km band-centre separation. The two published bands are NOT competing
answers. The campaign's stated deliverable needs restating accordingly; the
proposed framing is in the adjudication file. The frozen-branch ruling
(no-ocean interiors stay in the posterior) is **unaffected**.

**Next actions are enumerated** as `required_validation_before_acceptance` items
1-7 in the adjudication file. Item 1 (retire or repair the whole-difference
scaling) needs a USER DECISION because it moves
`tests/librations_test.py::test_h22_obs_m_reproduces_b2prime_measurement`, which
currently pins the now-suspect +1.40-1.43 sigma shift, and because the config
records a user decision to adopt option A that was taken on two premises now
falsified.

**Campaign-level risk flagged by the reviewer:** the recurring failure mode here
is not arithmetic — it is perturbing a physically meaningful quantity through a
DIFFERENCE rather than through the interface weight that actually carries it.
One root cause produced the escalation's 30.5 sigma, an agent's treatment-4
instability, AND the shipped kwarg's docstring/code divergence. Recommend adding
the Delta_rho identity as a STANDING INVARIANT on the libration module.

## Closed this session (with measurements, in validation_reports/enceladus_isostasy/)

- **B15** — Σ_model re-sourced from the H&M PDF Table 1 (5.3 / 1.7 / 4.4 e-6 +
  libration 0.00025 deg), verified against their Eq. 25 and Eq. 24. This
  **retracted** a reviewer claim that the values were ~2× too large.
- **Frozen-branch reparameterization** — onto a declared rho_rock ~ U[2200,2600]
  prior, after measuring a box-edge swing (frozen prior mass 22.7 → 43.3%) and a
  factor-2 prior tilt toward low-density rock under uniform-zb.
- **Mass-neutral rho_ice** — `isostasy.mass_neutral_shell_density`, now
  conserving against the **reduced stack's own mass** (not the template target;
  that bug injected a −0.19σ libration bias merely by enabling the nuisance).
- **Root cause of the frozen-node mass violation** — and **Titan is clean, no
  escalation** (frozen nodes +0.206 to +0.249%, ocean +0.140 to +0.366%).

## Two corrections I had to make to my OWN work — read these before trusting old numbers

1. **My blast-radius audit headline was wrong.** I used midpoint-edge
   quadrature; PP's convention is **outer-edge**. Re-run correctly, every
   residual across ~37k structures is negative or zero (max positive
   **0.0000%**), the ratified corpus spans −0.001 to −0.978%, and the claimed
   "+0.85% median in a ratified Titan cache" **does not exist**. Conclusion
   retracted: the mass gate **can** land globally.
   Knock-on: the r4 invariant's 5e-4 midpoint threshold also failed (it would
   have rejected the one healthy node it was designed to pass). Now v3:
   outer-edge, one-sided, step-scaled (~0.65% Enceladus / 0.44% Titan /
   0.36% Europa).
   **Lesson, written into the config:** an invariant may define its own
   convention; an *audit of existing data* must adopt the convention the data
   was built with.
2. **A regression test that could not catch its own bug.** My first
   mass-neutral test used a synthetic stack and passed under both old and new
   code, because a hand-built stack is mass-consistent by construction. Replaced
   with one that reduces a real cached node and asserts the coarse-graining
   residual exceeds 1e-4, so it self-invalidates if it ever stops discriminating.

## Blocked on reviewer sign-off (production physics — do not land unilaterally)

The **two-part** mass fix. The defect is one layer deeper than first diagnosed:
the guard at `Silicates.py:74-92` *does* fire and set `Do.VALID=False`, then
**resets it to True** expecting a retry — and that retry is wired only inside
the `Fe_CORE=True` branch (`LayerPropagators.py:1800`); the no-core calls at
1854 and 1926 have none. Enceladus is `Fe_CORE=False`, so the guard disarms
itself.
1. Hoist the `iSilStart` retry out of the `Fe_CORE` branch.
2. Then add the post-selection mass assertion as a tripwire.
A gate **without** part 1 turns frozen nodes from wrong into unbuildable.

**And deep-frozen Enceladus is feasible** — solving mass and MoI together gives
zb 51.78 / 53.93 / 58.29 km at rho_rock 2314 / 2334 / 2378, all landing on
C/MR² = 0.3350 (Iess 0.335 ± 0.001), inside the H&M/Park density band. So the
±0.08 MoI widening was never a physics necessity; it accommodated the solver
limitation. The user ruling that no-ocean cases stay is physically correct —
**and** that segment is currently unbuildable, which makes the two-part fix a
prerequisite, not a cleanup.

## Not blocking the build, but owed

- `rho_ice` bound rationale: now that it is a mass-neutral EOS nuisance,
  attainability is the wrong basis for U[915,935] (23% still outside the
  measured attainable range, itself sampled at a single salinity). Needs a
  SeaFreeze density error budget.
- The `Mtot_kg` field carries **two conventions** in the same package
  (`cache_builder.py:561` = target; `structure_cache.py:478` = achieved with
  target fallback). Worth reconciling; do **not** repurpose `Mtot_kg` — it is
  consumed as the conservation target in six places including two `Test/` files.
- B1 reachability, B4 zb cache mode, B12 finite-amplitude at Enceladus
  amplitudes, first-order H22 (−0.65% on H22 = 0.65σ on C22), support-edge
  mapping in (zb, w), and the code_gaps list (zb_km / compensation_C2 /
  rho_ice_kgm3 / libration_sys_frac unregistered; no C30 dispatch; no zb cache
  mode; `isostatic_hm2019` has no dispatch).

## Test baseline

40 passed + 1 strict-xfail across `librations_test.py`, `isostasy_test.py`,
`isostasy_hm2019_gate_test.py`, `enceladus_conductivity_test.py`. The xfail is
the intentional tripwire for the future `K_int` 8π/15 repair — **it must stay
xfail**, and that repair must not land alone.

## Code_gaps ledger update (Machine A, 2026-08-15)

CLOSED by the reviewed implementation pass (commits 77904840/d3819fa1/
244ccd5a/fdb81c47; manager-verified, 14 new tests + 40+1xfail baseline
green): parameter registry (zb_km, compensation_C2, rho_ice_kgm3,
libration_sys_frac), C30 dispatch, isostatic_hm2019 forward-model
dispatch (ocean branch; separate-then-sum finite-amplitude per the
ratified module convention, cross-term vs combined-FA documented <0.2
sigma), and B4 zb cache mode (build_zbw_grid_cache; solved Tb recorded
per node; real 6-node smoke build). STILL OPEN (reviewer-blocked or
ratification-time): frozen-branch builder path + MAJOR-1 mass invariant
(do not land unilaterally), Sigma_model likelihood inflation wiring,
config-schema reconciliation (the candidate JSON's top-level isostasy
keys raise TypeError in InferenceConfig.from_json — dispatch currently
reads config.metadata['isostasy']), and the B2' Delta_rho repair
(user-ruled 2026-08-15, §0.23 — B executes under reviewer sign-off).
