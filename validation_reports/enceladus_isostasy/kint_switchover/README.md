# K_int switch-over validation — 2026-11-10 (Linux container, working tree on 79c28066)

Clears the NO-GO conditions from the 2026-11-03 K_int review: (V1) separability
of "repaired rigid H22/C2 × per-node hydrostatic elastic factor" against a
deformable branch WITH the figure treatment, C2 ∈ [0,1] at nodes 1220/1540;
(V2) the dense-zb repair ratio at C2 = 0.81.

Driver: `plans/scripts/kint_switchover_validation.py` →
`kint_switchover_validation.json`, `f_table_checks.json`.
Test cache (88 zb × 40 w), ALMA elastic y1 per node, H22_obs_m = 857.0 (config).

## The port, and the binding-review correction
`PlanetProfile/Gravity/Librations.py`: the Δρ-consistent H22/C2 figure treatment
now composes with the deformable (VH2013, rigid=False, ocean=True) branch for
the reduced single-shell stack, opt-in via the existing kwargs. A shared
`_figure_scales()` helper feeds both branches (scales numerically identical; the
rigid branch is byte-identical to its pins).

**First-pass review (2026-11-10): BLOCK.** The initial port scaled the whole
f2 products in Ks and the moment in Ksz, i.e. it applied the figure to the
tidal-response terms. Those terms are BmA·k2/kf = k2·q·M·R² (hydrostatic
identity kf = BmA/(q·M·R²)): the torque on the tidal bulge, set by the forcing
and the rheology, independent of the static figure. The y1 = 0 anchors cannot
see this (both forms coincide at f2 = 1). **Corrected per the review:**
`Ks = 1.5ω²[(BmA_s_fig + BmA_ot_fig) − (BmA_s(1−f2s) + BmA_ot(1−f2ot))]`, and
Ksz keeps the hydrostatic BmA_s. K_int, Ki1, kf/k2/f2 and the interior terms
stay hydrostatic throughout (the K_int convention mirrors the rigid branch).
Caveat (review): VH13 is not in papers/; the identity rests on the code's own
kf definition and linear-response physics.

Anchors (tests/librations_repair_test.py, 67 passed + 1 strict tripwire xfail):
- y1 = 0 deformable+figure == repaired rigid+figure to rel 1e-10,
  C2 ∈ {None, 0.25, 0.5, 0.81, 1.0}, vh16 + b2prime stacks.
- **Discriminating test at y1 ≠ 0** (`test_figure_treatment_only_scales_static_moments`,
  via a `diagnostics` hook): Ksz figure-independent; the figure increment of Ks
  identical at y1 and y1/2 and equal to the static-moment increment. The
  defective form fails both.
- Corrected deformable+figure pins at nodes 1220/1540, C2 = 0.81:
  **432.6357784877 / 401.7141908788 m** — independently reproduced by the
  reviewer (432.64 / 401.71) before pinning.
- C2 → 0 == surface-only; no-figure deformable path pinned unchanged
  (smoke node 491.9794954534861); multi-layer-shell guards.

## V1 — separability (corrected form)
elastic/repaired no-figure F: 0.97936 (node 1220) / 0.98135 (node 1540).

| C2 | sep. err (1220) | sep. err (1540) |
|---|---|---|
| 0 | −0.113% | −0.102% |
| 0.5 | −0.079% | −0.072% |
| 0.81 | −0.057% | −0.053% |
| 1.0 | −0.043% | −0.041% |

Dense-zb sweep (80 nodes, zb 13–45 km, w 3.46 ppt, C2 = 0.81): separability
−0.080%..−0.037%, mean −0.053%. **|sep| ≤ 0.2% everywhere measured** —
the review's criterion for closing C7 by the composed channel. (Review note:
its own scan found separability growing to ~+7% at zb 5–6 km, where
F ≈ 0.67–0.73 — outside the posterior but to be stated as a caveat.)

## V2 — repair ratio and frozen-channel bias
Repair ratio at C2 = 0.81: 1.01953 / 1.01793 at the review nodes (matches the
review's +1.95%/+1.79% exactly); dense-zb range 1.0098–1.0424.

Frozen bias (shipped rigid × 1.008 vs corrected deformable+figure):
- nodes 1220/1540 at C2 = 0.81: **+0.89% / +0.85%**; C2 = 0: −0.21%/−0.17%;
  C2 = 1: +1.16%/+1.10%.
- dense-zb sweep at C2 = 0.81: +0.66%..+1.52%, mean +0.89%.

Review's impact scaling (d ln L/d zb = −0.0354/km ⇒ 0.28 km per 1%): the frozen
convention over-predicted libration by ≈ +0.8–0.9% at the posterior C2, so the
published-provisional zb posterior is **≈ 0.2–0.3 km too thick** (−0.13 km at
C2 = 0.5, −0.30 km at C2 = 1, ≈ +0.04 km at C2 = 0). Fixed-C2 scaling, not a
joint-posterior result; to be confirmed by the posterior re-run.

## Runtime channel (review ruling): option B
Repaired rigid+figure × a per-node elastic-factor table F(zb, w) =
elastic_nofig / repaired_nofig (C2- and H22-independent by construction).
The confirming review (PASS WITH CONCERNS → conditions met below) ruled
option (a): a linear ρ_rock correction in the lookup,
F·(1 + c·(ρ_rock − ρ_ref)), with c = dlnF/dρ_rock taken ALONG THE RUNTIME
MASS CLOSURE (mass_rho_rock: R_s re-solved; identity at the node's own ρ_s),
NOT a fixed-geometry density bump. That distinction matters: the
fixed-geometry FD gave −8.4e-6 per kg/m³, but along the closure manifold the
R_s response nearly cancels it — closure-consistent c ≈ +1.2e-6 per kg/m³
(sidecar range +4.6e-7..+4.1e-4, large values only at thin shells).
Builder: `plans/scripts/enceladus_libration_factor_sidecar.py` (schema
lib-factor-sidecar-v1: F, dlnF_drho_rock, rho_rock_node_kgm3, provenance with
the cache md5). Acceptance checks (`plans/scripts/kint_f_table_checks.py` →
`f_table_checks.json`, all PASS):
- LOO zb-linear interpolation, 2×-coarsened grid, zb 20–32: max 1.7e-4.
- ρ_ice ±10 kg/m³, closure-consistent, node y1: |ΔF/F| ≤ 1.4e-4.
- Required test (a)-2 — c against a fresh closure FD with ALMA y1 RECOMPUTED
  per perturbed stack (3 regions × 12 constant-property sublayers,
  volume-weighted μ/η): c ratios 0.990 / 0.988 at nodes 1220/1540 (gate 5%).
- Required test (a)-3 — linearized F at ρ_rock = 2411 vs the exact factor on
  the mass-closed stack: |diff| ≤ 7.1e-6 (gate 1e-4).

## sys_frac prescription (review, to confirm)
Keep σ = 0.004, truncation ±0.012; rewrite the rationale: separability ≤0.13%
(window), F interpolation ≤0.07%, K_int hydrostatic-β convention +0.27%(C2=0)
→ −0.02%(C2=1), ALMA Re k2 incompressibility acting on 1−F ≈ 0.1%; quadrature
≈ 0.3–0.35%. Model-form caveats NOT in sys_frac: viscoelastic rheology
(Maxwell ≤3%, Andrade/soft-base ≤~6%, one-sided), thin-shell (zb <14 km)
composed-form failure, backend compressibility spread, first-order H22
projection.

## Re-derived B2' C2 bracket (fiducial zb-match to 0.092°, D_ocean 36 km)
Cross-check: the shipped numbers reproduce to all quoted digits
(27.3395 / 26.6531 / 25.9944 km at C2 = 0 / 0.5 / 1).

| convention | C2 = 0 | C2 = 0.5 | C2 = 1 |
|---|---|---|---|
| shipped (old record) | 27.3395 | 26.6531 | 25.9944 |
| repaired, no F | 27.9998 | 27.1095 | 26.2614 |
| runtime channel (repaired × F, w 3.46 ppt) | 27.4963 | 26.6137 | 25.7742 |

Against the old runtime (shipped × 1.008 ⇒ ≈ +0.23 km on each match), the new
runtime bracket sits ≈ −0.07 km (C2 = 0) to −0.45 km (C2 = 1) thinner —
consistent with the review's −0.2..−0.3 km headline estimate at C2 ≈ 0.81.
The repaired formulation also uses EXACT ellipsoid pressure moments, which
repairs the recorded Bsp_Asp linearization caveat (the +0/−0.4 km C2 = 1
model-coefficient limit in `metadata.B2prime_signoff` condition 1).

## Status
- Port + Ks/Ksz correction + discriminating tests: **`verified`** (second
  binding review PASS WITH CONCERNS, both conditions — the closure-consistent
  ρ_rock FD and the committed checks generator — met; suite 69 passed, pins
  reproduced independently).
- V1/V2 + F-table acceptance checks: **`verified`** (all gates pass; every
  float of the validation JSON reproduced independently by the reviewer).
- Switch-over NO-GO: **cleared — GO** for the switch-over commit with the
  reviewed prescription (runtime channel B with the ρ_rock term; pins
  1.008196 → 0.973798 and the tripwire flip; multiplicative_factor retired
  for the F sidecar; sys_frac σ = 0.004 ± 0.012 kept with the rewritten
  rationale; C2 bracket re-derived above).
- Caveats (stated, not blockers): thin-shell (zb < 13–14 km) composed-form
  degradation (F down to 0.95 at 13 km, 0.67 at 5–6 km); viscoelastic
  rheology one-sided (Maxwell ≤3%, Andrade/soft-base ≤~6%); VH13 not in
  papers/ — the Ks/Ksz identity rests on the code's own kf definition and
  linear-response physics.
- Headline wording (review (v)): "at fixed C2 ≈ 0.81 (posterior C2), the
  frozen convention biased libration by +0.8–0.9% ⇒ zb ≈ 0.2–0.3 km too
  thick (≈0 at C2 = 0, −0.30 km at C2 = 1); fixed-C2 scaling, provisional
  until the posterior re-run."
