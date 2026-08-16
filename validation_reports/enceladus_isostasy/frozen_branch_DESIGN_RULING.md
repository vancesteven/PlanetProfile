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
