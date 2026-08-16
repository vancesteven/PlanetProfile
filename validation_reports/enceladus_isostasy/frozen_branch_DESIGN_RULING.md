# Frozen-branch design ruling (scientific reviewer, 2026-08-15)

Verdict: BLOCK the config's frozen-branch design as written; APPROVE the
redesign below, conditional on feasibility smoke A1. Full ruling in the
session record; operative content:

- F1 (binding obstacle): PP selects the seafloor by argmin |C/MR2 - 0.335|
  on every path — an UNDECLARED MoI conditioning that makes the declared
  rho_rock U[2200,2600] prior unreachable (collapses to ~2334). Deeper
  than the iSilStart retry; the retry fix is DEMOTED to optional.
- F2 (the route): Do.ConstantProps['Inner'] + Do.HYDROSPHERE_THICKNESS
  gives specify-zb/derive-rho_rock with ANALYTICALLY EXACT mass closure
  (LayerPropagators.py:1571,1650); C/MR2 becomes derived, not gated.
- F4: retry_frozen_as_no_ocean is INOPERATIVE at Enceladus (needs HP ice;
  basal P ~6.8 MPa is 30x below PminHPices) — remove from config.
- F6: true frozen support is zb [46.74, 65.56] km (config's [47.45,65.26]
  was computed at rho_ice=925 only).
- Design: separate frozen arrays in the cache (schema v3.2-zbw-joint; NO
  ragged shared axis, NO w-crossing); sampling space (rho_rock, rho_ice,
  libration_sys_frac) <-> build space zb via the analytic bijection
  frozen_zb_from_mass (no iteration; MoI derived along it); branch drawn
  by TWO SEPARATE nested-sampling runs + post-hoc log-Bayes-factor at
  declared 0.5/0.5 odds (single-theta branch indicator REJECTED:
  non-nested supports would need a pseudo-prior = the C16 fragility
  relocated).
- Invariants I-F1..I-F6 (mass exact at 1e-6 on the constant-rho path;
  parameterization closure <=12 kg/m3 & <0.25 km; per-sample loop
  closure I-F5 at 3e-3 -> -inf on failure; I-F6 both-directions
  non-conditioning test: a frozen set entirely inside 0.335±0.001 FAILS
  the build).
- Cuncertainty ruling (PARTIALLY OVERTURNS RESUME_NOTE): the ±0.08 was
  doing one illegitimate job (masking mass violations — cured) and one
  legitimate one (letting the prior span). Preferred: NO MoI gate at all
  (Cuncertainty forced 0, CMR2 recorded as diagnostic). If a window is
  unavoidable: 0.015 (min non-gating 0.01201 + margin), NOT 0.001 —
  executing MODERATE_4's ±0.001 instruction as written would collapse
  the prior 10-fold and is COUNTERMANDED. Ocean branch keeps ±0.08 but
  must DECLARE the argmin-MoI selection (currently undeclared, in
  tension with cmr2_not_an_observable).
- Task list A1-A7 (Machine A), C1-C5 (Codex), B1-B2 (Machine B); A1
  feasibility smoke blocks all. B2' remains the campaign-gating item;
  this work is parallel to it except B1.

## Implementation record + manager rulings (2026-08-15)

A2-A7 IMPLEMENTED (commits 1c0ad391/96c5bfee/42ddbf26/6e74f85b/08dae94b/
21d11642; suite 79 passed + 1 xfail unchanged; manager-verified). A1
items all closed with root causes: (i) no-liquid via Dhsphere_m=0 (the
argmin then selects the ice base; iSilStart pins the sweep below it —
explains the earlier ocean contamination AND why Tb changed zb);
(ii) exact zb binding via ice-base PRESSURE (SPECIFY_ICEI_BOTTOM_PRESSURE
+ secant, residuals <=0.1 km vs the 2.1 km miss through ICEIh_THICKNESS);
(iii) the no-MoI-gate route is NOT expressible at Enceladus (the
Cuncertainty=0 path requires >=1 ocean layer) — the ruled 0.015 fallback
adopted and verified non-conditioning both structurally (window can only
reject, never re-select; rejected node acquires liquid and dies on the
no-liquid invariant) and numerically (max |CMR2-0.335| = 0.0108, margin
>=0.004). Key measurements: F5 regression exact; 3/3 real frozen nodes,
0 liquid layers, mass residual <=2.2e-16, rho_rock closure 0.0 kg/m3;
I-F6 fires both directions; I-F5 fires exactly at threshold both signs;
frozen C30 zb-independent to 0.0007 sigma.

MANAGER RULINGS on the implementer's five open items:
1. Branch-comparison asymmetry (frozen uniform-density interior vs ocean
   porous Comet_67P EOS at phi~0.31; plus the ocean branch's argmin-MoI
   selection): STAYS OPEN, routed to the r5 ratification pass — the r5
   reviewer rules equalize-vs-bound BEFORE any branch odds publish.
   Manager preliminary: bound and caveat (equalizing would re-open the
   ocean branch's ratified EOS treatment late in the campaign).
2. Frozen cache dispatch gap (no sampled-(rho_rock,rho_ice) -> node
   lookup; nearest-node snapping would eat 80% of I-F5's margin):
   ACCEPTED as a real gap — follow-up implementation dispatched
   (interpolating dispatch through the analytic bijection, never
   snapping; see A8 below).
3. "I-F2 = no-liquid invariant" naming: the implemented content is the
   full ruling's I-F3 (branch flag from the phase array) — content
   CONFIRMED correct. The full ruling's actual I-F2 (stored-array
   outer-edge mass quadrature audit) is near-redundant under analytic
   closure but ORDERED — added to A8.
4. isostasy.frozen_cmr2 diagnostic helper: ACCEPTED (tested,
   display-only).
5. Superseded test-pin rewrites: ACCEPTED (both assert the ruled
   behavior with recorded reasoning).

A8 (follow-up, dispatched): interpolating frozen-branch cache dispatch
in _struct_for_hydrosphere-equivalent (1-D in zb; I-F5 asserted on the
interpolated stack), + the stored-array I-F2 quadrature audit in the
builder.
