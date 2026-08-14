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

**B2' is reopened and is the campaign's largest open systematic.**
`b2prime_REOPENED_figure_coupling.json`.

The libration figure convention is not a scaling — it is a **structural
coupling to the compensation state**. Four treatments give matched shell
thickness 24.70 / 25.80 / 27.34 km, and H&M's own Eq.-12 treatment gives 30.5σ
with no zb matching at all. The spread (2.6 km) is comparable to the 6 km
H&M-vs-Park separation the campaign exists to adjudicate, and since the
libration is the **sole** thickness channel, this sets where the answer lands.

**Do not guess a fourth scaling.** The recommended next step is to reproduce
H&M's published Table 2 shell thicknesses from their inputs *with the libration
in the loop* — discriminating the treatments against a published answer, the
same strategy that made the B13 gravity gate trustworthy.

Also from that work: **the libration DOES depend on `compensation_C2`**
(0.1005 → 0.1824 deg as C2 goes 0 → 1), contradicting an earlier reviewer
finding of dlibration/dC2 = 0.00σ — that was measured on the hydrostatic figure
where the Airy root cannot reach the torque. So the "gravity = compensation
state, libration = thickness" split holds only in the hydrostatic formulation.

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
