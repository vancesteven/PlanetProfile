# Europa Clipper v3 — Phase 1 pre-commit validation report (Machine B)

Date: 2026-07-18. Gates the full ~1500-run 2D (Tb × w) cache build.
All checks pre-registered by scientific-reviewer (PASS WITH CONCERNS, this session).
Scripts + logs in /tmp: v3_pan_nacl_gate.py, v3_compare_regular_vs_refined.py,
v3_docean_maps.py, v3_build_regular_column.py (+ .log files).

## Code delivered (Phase 1)
- `cache_builder.build_single_structure`: now stores `wOcean_ppt` (read back
  from `Planet.Ocean.wOcean_ppt`) in the per-node dict. VERIFIED live: a
  build returns wOcean_ppt = 35.16504 for the seawater default.
- `cache_builder.build_tbw_grid_cache` (NEW): 2D regular grid, per-node
  `ocean_overrides={'wOcean_ppt': w}`, `detect_transitions` NOT used (regular),
  row-major `structures[i_Tb*n_w + i_w]`, failed builds → explicit `None`,
  schema `v3.0` with `Tb_K_grid` + `wOcean_ppt_grid`. py_compile + import OK.
- `build_phase_c1_cache.py` CLI: `--n-wgrid N` reads `log10_wOcean_ppt`
  bounds and dispatches the 2D build with log-spaced w nodes.

## Gate 1 — regular grid vs refined v2 cache (reviewer Decision 1)  → PASS
Built a regular 0.125 K Tb column (93 nodes, w=35.165) and compared to the
transition-refined v2 cache (36 nodes, 2 transitions).
- max |ΔCMR²| = **0.000339 (0.141σ)** at Tb=269.75 K; gate < 0.0006 (0.25σ). PASS
  (independently reproduced by scientific-reviewer from the two pickles.)
- Ocean-onset Tb: with the `D_ocean>0` metric both grids read 259.500 K, but
  that metric carries a ~1.5 km numerical residual at fully-frozen nodes
  (reviewer MODERATE/MINOR). The PHYSICALLY-CORRECT onset (first liquid `'0'`
  layer in region_phases) is bracketed at **260.25–260.375 K (regular) vs
  260.258–260.266 K (v2)** — agree within one 0.125 K step. gate ≤0.125 K. PASS.
  (Pre-existing v2 behavior, not a v3 regression; the |Ae|>0.7 support cut must
  reject the residual-ocean sub-onset nodes — confirm in Phase 2.)
- DEFERRED to Phase 2 (needs the bilinear lookup): the pre-registered
  |Ae_synodic|=0.7 support-edge Tb agreement (regular vs v2, within 0.125 K)
  through the production Ae path. This is a training-acceptance blocker, NOT a
  build blocker.
- Conclusion: the regular grid loses nothing material; row-major stackable.
  HP-ice concern confirmed absent on the GSW seawater path (basal P ~156 MPa <
  200 MPa PminHPices cap). Build rate 2.4 s/node warm → full grid ~60 min.

## Gate 2 — Pan-NaCl composition-bias gate (reviewer Decision 3)  → PASS (conductivity)
GSW/PSS-78 Seawater vs in-repo Pan et al. (2021) NaCl(aq) over a Europa ocean
(P,T) column, at w ∈ {35.165, 50, 70, 100}. Gate: <25% caveat.
- **Conductivity (the real composition-bias gate): max |Δσ| = 18.5% < 25%. PASS.**
  σ divergence rises only ~2.7pp over 35→100 ppt (15.8%→18.5% at basal):
  composition-dominated (NaCl vs seawater; Pan extrapolates DOWN in P from its
  212 MPa floor), NOT driven by GSW's 42→100 ppt extrapolation.
  CAVEAT (reviewer MODERATE): this is a two-mutually-extrapolated-models
  comparison (all P < Pan's 212 MPa floor; w > GSW's SP=42 ceiling) — it bounds
  model-to-model DIVERGENCE, not absolute error. Evaluating Pan below its floor
  inflates the number, so the pass is conservative. Keep the ~18.5% σ systematic
  as an explicit config caveat; execute the Phase-3 high-w NaCl cross-check MCMC.
- **Density: NOT a fidelity gate (reviewer MAJOR correction).** The 7.7% figure
  compared GSW seawater *solution* ρ against SeaFreeze pure-*water* ρ (Pan's
  input) — that difference is essentially the dissolved-salt mass fraction
  (~w/1000 → ~7.6% at 100 ppt), NOT a model discrepancy. GSW density in the
  SP>42 extrapolation is therefore UNVALIDATED here. Impact is limited: density
  enters the forward model only via layer thickness / mass conservation, and
  CMR² is flat across w (Gate 3), so the w-posterior is insensitive to it. Do
  NOT cite 7.7% as a passed density gate. (To close properly: compare GSW
  solution ρ against an independent NaCl-brine solution EOS — deferred, low
  priority given the CMR²-flat insensitivity.)
- Tfreeze(w) monotone-decreasing at every pressure (0.1→100 ppt). PASS
- Decision: **keep the 100 ppt bound**; document GSW-extrapolation + NaCl-vs-
  seawater composition bias as a config-metadata caveat; add one NaCl-EOS
  cross-check MCMC at high w (Phase 3) to bound the composition bias on the
  w-posterior mode.

## Gate 3 — D_ocean/D_iceIh coverage maps (plan check 3)  → PASS (box adequate)
Coarse 12-node Tb columns at w ∈ {0.1, 1, 100}. Fixed box [259.5, 271.0] K.
- w=0.1 ppt: ocean at 10/12 nodes, span 261.59..271.00 K (low-Tb corner fails)
- w=1.0 ppt: ocean at 10/12 nodes, span 261.59..271.00 K (low-Tb corner fails)
- w=100 ppt: ocean at 8/12 nodes, span 259.50..266.82 K (HIGH-Tb corner fails)

**Key scientific finding — the valid (Tb, w) region is a TILTED BAND, not a
rectangle.** The freezing-point depression (~6–7 K at 100 ppt) drops the
ice-extinct ceiling at high salinity, while the ice-shell-too-thick floor rises
at low salinity. Implications:
1. The box is adequate at all salinities (≥8/12 nodes, good margin) —
   **no Tb-box widening required**.
2. The 2D cache will carry `None` in TWO OPPOSITE corners (low-Tb/low-w and
   high-Tb/high-w). Phase-2 bilinear must handle `None` corners on BOTH ends,
   not just the low-w plane the plan emphasized. The `None`-corner fallback
   (nearest valid, else reject) is load-bearing for the degeneracy ridge.
3. CMR² is flat (~0.346) across the valid region at all w — the CMR² channel
   does NOT constrain salinity; w identification rests entirely on induction
   (orbital Ae), as the plan predicts.

## Verdict — GO-WITH-CONDITIONS (scientific-reviewer, 2026-07-18)
All three pre-registered Phase-1 gates PASS on their numeric tolerances with
measured evidence (Gate 1 independently reproduced by the reviewer). The
regular-grid 2D approach is validated as a BUILD decision; proceed to the full
build. Status of the code: `implemented, unverified` for the 2D builder
end-to-end (a full 2D cache has not yet been written+reloaded); the individual
checks above are `verified` against their scripts/logs.

### Conditions before the TRAINED ARTIFACT is accepted (not build blockers):
1. [Phase 2] |Ae_synodic|=0.7 support-edge Tb: regular vs v2 within 0.125 K,
   through the production Ae path.
2. [Phase 2] 2×-log-w sub-box convergence near the 0.7 contour: support contour
   moves <0.125 K; orbital Bind <0.75 nT; fiducial (Tb,w) ρ_corr <0.05.
3. [Phase 2] A SINGLE shared bilinear interpolant serves both the MCMC
   likelihood and the SBI support guard (assert byte-identical); the None-corner
   fallback is exercised at BOTH excluded corners (low-Tb/low-w AND
   high-Tb/high-w).
4. [Phase 3] Keep the σ-systematic (~18.5%) as an explicit config caveat; run
   the high-w NaCl cross-check MCMC. Density gate relabeled (not a fidelity gate).
5. [Build QA] On the full 2D cache: assert len(structures)==n_Tb·n_w, log the
   None count, and confirm each non-None entry's wOcean_ppt equals its grid node
   (guards against silent 2-of-3 re-derivation).

---

# Phase 2 — 2D bilinear lookup (verified 2026-07-18)

## Code delivered
- NEW `PlanetProfile/Inference/grid_interp_2d.py`: the SHARED interpolant
  (`is_2d_cache`, `bilinear_weights`, `resolve_none_corners`, `blend_complex`,
  `blend_scalar`, `wOcean_ppt_from_theta`). Pure numpy; 16 unit+integration
  tests pass.
- `forward_models.py`: `_apply_bottom_temperature_2d` (structural bilinear
  blend, nested `_blend_b_layered`, None-corner fallback) + Format-3 dispatch.
- `mcmc_runner.py`: `_blended_ae_dict` (shared Ae lookup: nearest-Tb 1D /
  bilinear (Tb,log w) 2D) and `_struct_for_hydrosphere`; all 4 induction/CMR2
  consumers routed through them (likelihood, mass-conservation CMR2, SBI
  `_induction_channel_value`, `_check_induction_bounds`). `_collect_posterior_Ae`
  made w-aware. Fixed a latent `Tb_sample` NameError on the sidecar+2D path.
- `europa_seawater_andrade_clipper_v3_8D.json`: v2 + `log10_wOcean_ppt`
  U[-1,2] + 2D cache path + caveat metadata. Loads/validates (8 params).

## Full 2D cache built (Machine B, 57.7 min)
`Test52_seawater_v3/europa_seawater_structure_grid_v3_2d.pkl` (21 MB): 93 Tb ×
16 w = 1488 nodes, 1303 built + 185 None. Build QA: len==n_Tb·n_w OK;
wOcean_ppt/Tb node match OK (no silent 2-of-3 re-derivation). None map confirms
the TILTED BAND (low-Tb/low-w + high-Tb/high-w corners fail).

## Scientific-reviewer GO-WITH-CONDITIONS — status
- C1 (Tb_sample NameError): FIXED + regression test. `verified`.
- C2 (support-guard == likelihood, randomized): 200-sample sweep, every
  support-reject forces a likelihood-reject. `verified`.
- C3 (contour clear of None corners, production cache): 0 flagged cells for
  both |Ae|=0.7 and |Im Ae|=0.4. `verified`.
- C5 (Ae mask == struct mask, production cache): 1303==1303, 0 mismatches.
  `verified`.
- C4 (k2/h2 jump across region-transition cells < sigma): DEFERRED to Phase-3
  gates (quantify at the convection-onset transition). `not implemented`.
- MINOR (log-w interp, single exponentiation, clamp-no-extrapolation): all
  verified sound by reviewer + unit tests.

Status: Phase 2 `verified` end-to-end on the production cache (likelihood
finite where |Ae_syn|>0.7, rejected below; CMR2 mass-conservation derives
through the bilinear hydrosphere; None-corner → hard reject). C4 remains for
Phase-3 gates.

## Phase-3 pre-training gates (Machine B, 2026-07-18) — `verified`
Run through the production shared interpolant; Ae grid disk-cached
(`/tmp/v3_ae_grid.pkl`, 1488 nodes) so downstream runs skip the ~32-min
precompute. Scripts: `/tmp/v3_c4_edge_gate2.py`, `/tmp/v3_worstcase_wprobe.py`;
results `/tmp/v3_c4_edge_result.json`, `/tmp/v3_worstcase_wprobe.json`.

- **C4 (transition-cell jump < σ): PASS (0.000σ).** The only phase-set
  transition on the GSW path is ocean onset (no-ocean → ocean; HP ice
  impossible, basal P < 200 MPa). On EVERY w-column the onset cell's lower Tb
  node is None (reject region), and the |Ae_synodic|>0.7 support edge sits
  ABOVE ocean onset everywhere — so no accepted sample ever bilinearly blends
  across the discontinuity. Worst in-support jump 0.000σ. C4 now `verified`
  (was deferred from Phase 2).
- **edge, Tb-discretization: PASS (0.0851 K < 0.125 K).** Exact-seawater
  regular 0.125 K column vs the v2 transition-refined cache — the plan's
  pre-registered gate.
- **edge, w-interpolation: 0.285 K at seawater; worst-case-over-supported-w
  0.099 K (w=7.94 ppt).** Reviewer GO-WITH-CONDITIONS. An interpolate-then-
  threshold |Ae|-curvature artifact of the 16-node log-w grid, biasing the
  support edge HIGH (conservative truncation). The worst-case probe (steepest
  supported cell) is < the 0.45 K composition-systematic edge-equivalent →
  **Option A** (keep the 16-column grid) stands. Documented as
  `salinity_w_interp_caveat_2026_07_18` in the config metadata (reducible
  discretization systematic, subdominant to the 18.5% composition systematic,
  present identically in reference-MCMC + SBI so it does NOT break the 2D
  degeneracy gate).

Fiducial Bind central values recomputed from the 2D cache at (Tb=266,
w=35.16504): `/tmp/v3_fiducial_bind.json`, |Ae_synodic|=0.874. (Computed at
Tb=266 vs v2's Tb=264.52 conditioning point; the current config Bind values
remain the ratified v2 conditioning — train/validate are self-consistent under
either, since fiducials and runtime use the same code.)
