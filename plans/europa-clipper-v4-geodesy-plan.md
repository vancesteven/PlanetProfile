# Europa Clipper v4: Mazarico-projected geodesy (k2 sigmas + C20/C22 option)

Author: Machine A, 2026-07-19 (user-directed). Machine B runs all
compute; this plan is written to minimize Machine B token use — all
decisions are made here, B executes. Status: `not implemented`.

Scientific review 2026-07-19 (opus): **PASS WITH CONCERNS** — physics of
the forward map confirmed (RD inverse algebra exact; q_r rotational
parameter correct, tidal contribution already in the 1/4 and 10/3
factors; signs correct); four binding changes, ALL INCORPORATED below:
(R1) pin the normalization convention of Mazarico's sigmas before
freezing anything; (R2) replace Radau–Darwin with a direct Clairaut
k_f integration from the cached rho(r) — the ~1% RD systematic is
~1e-6 in C22 = 5×sigma(C22) and aliases into the non-hydrostaticity
deliverable with near-unit efficiency; (R3) state honestly that free
2-D nuisances remove the interior MoI constraint from gravity (the
plan's earlier manifold-component identifiability claim was FALSE);
(R4) the ratio-preserving injection arm is exactly degenerate with a
C̄ shift — pre-register it as expected-degenerate, not recoverable.

## Why

1. **v3's k2 conditioning values are NOT the Clipper projections.** The
   Clipper configs freeze `Re_k2 [0.25, 0.06]`, `Im_k2 [0.0, 0.06]` —
   sigma 0.06 is the *requirement*-level accuracy, not the projected
   performance. Mazarico et al. (2023, SSRv 219:30, Table 5) project
   k2 = 0.2–0.3 with sigma (1.4–1.8)e-2 (49-flyby precision < 0.02,
   confirmed by the Clipper geodesy overview paper). Amortized artifacts
   freeze sigma at training time → retrain required.
2. **C20/C22 instead of C/MR2 (user 2026-07-19).** Clipper will estimate
   degree-2 gravity independently (ratio uncertainty 0.003, correlation
   0.03) — projected sigma(C20) (3.5–8.5)e-7, sigma(C22) (1.5–2.3)e-7.
   The current CMR2 observable [0.3547, 0.0024] bakes the hydrostatic
   Radau–Darwin assumption into the DATA (Anderson et al. 1998 imposed
   J2/C22 = 10/3). Modeling C20/C22 directly moves the hydrostatic
   assumption into the forward model where it belongs, and lets the
   inference REPORT the degree of non-hydrostaticity instead of assuming
   zero. Pre-Clipper context: Anderson C22 = 131.5e-6 ± 2.5e-6
   (hydrostatic imposed); Gomez Casajus et al. 2021 reanalysis
   C22 = 139e-6 ± 25e-6 without that constraint — a 9e-6 offset at
   0.36 sigma, CONSISTENT with hydrostatic (review R-refinement: v4
   sets an upper limit on non-hydrostaticity; it does not chase a
   claimed detection).

**Honest physics trade (review R3, binding):** degree-2 gravity WITHOUT
the hydrostatic assumption genuinely cannot constrain the moment of
inertia — with two free nuisance offsets, (C20, C22) carry essentially
zero interior information after marginalization (C̄ bounded only to
~±0.02 by the nuisance prior box, vs CMR2's ±0.0024: a ~9x loss).
Interior constraints in v4 come from k2 (now 4x tighter), h2, and
induction; expect rho_core/R_core/C̄ to become markedly more
prior-dominated, and SBC/crosscheck on those to possibly degrade —
report the v3→v4 per-parameter constraint change as a headline result
("what does Clipper geodesy buy?" may honestly answer: for the
interior, the k2 tightening; for non-hydrostaticity, an upper limit
whose width is set by the C̄ uncertainty (~5–9e-6), NOT by
sigma(C22) = 2e-7).

## Physics: hydrostatic C20/C22 from the existing structures (NO cache rebuild)

The committed 2D cache stores per-node `omega` (spin rate = mean motion,
synchronous — state this in gravity_obs metadata since q_r is now
load-bearing), `R_body_m`, `Mtot_kg`, AND the full radial density
profile `r_m`/`rho` per node. Forward map (review R2, binding — Clairaut
integration, NOT Radau–Darwin):

- q_r = omega^2 R^3 / (G M)  (Europa: 4.98e-4; rotational parameter —
  the tidal contribution for a synchronous body is already inside the
  1/4 and 10/3 factors below; review-confirmed)
- **k_f per structure node by direct integration of Clairaut's
  equation** over the cached rho(r) (first-order figure theory; solve
  the Radau form d(eta)/dr with eta(0)=0, k_f = (3 − eta_s)/(1 + eta_s)
  ... using the standard surface relation; pure numpy ODE over ~100
  cached layers, microseconds per node, precomputed once per cache
  load). This ELIMINATES the ~1% Radau–Darwin systematic (~1e-6 in
  C22 = 5x sigma(C22)) that would otherwise alias into the
  non-hydrostaticity deliverable with near-unit efficiency. The RD
  closed form (y = (5/2)(1 − (3/2)C̄), k_f = (4 − y^2)/(1 + y^2))
  stays as a unit-test cross-check (agreement to ~1% required) and as
  the GUI's quick inverse readout.
- Hydrostatic, synchronous rotation: C22_h = k_f q_r / 4;
  J2_h = (10/3) C22_h; C20_h = −J2_h.
- Convention verified numerically (Machine A 2026-07-19): C̄ = 0.346 →
  k_f(RD) = 1.044 → C22_h = 130e-6, J2_h = 434e-6 — reproduces
  Anderson's hydrostatic pair.
- Interpolation: k_f is a per-NODE quantity (from rho(r)); blend it
  bilinearly in (Tb, log10 w) exactly like the other structural
  scalars (shared grid_interp_2d path), then apply the sampled-core
  correction consistently with how per-sample CMR2 handles the core —
  DESIGN NOTE for implementation: the cached rho(r) is the
  template-core profile; the sampled core (R_core, rho_core,
  mass-conservation rho_sil) modifies the interior density. The
  Clairaut integration must run on the SAME per-sample composite
  profile the CMR2 recompute uses, not the raw cached one. This is the
  one non-trivial piece of gravity_obs.py; unit-test it against the
  analytic two-layer k_f solution.

**Statistical structure (review-ratified R3):** with sigma(C22) ~ 2e-7
on a 1.3e-4 signal, a PURE hydrostatic forward model would let real
non-hydrostatic power contort the interior posterior — bias, not
constraint. Therefore sample TWO nuisance parameters:

- `dC20_nh`, `dC22_nh` — additive non-hydrostatic offsets,
  priors uniform [−2e-5, 2e-5] (~2x the nominal GC21 offset; a
  conservative broad prior for an upper-limit measurement;
  pre-register a prior-sensitivity check exactly as for salinity).
- Model: C20_model = −(10/3)(k_f q_r/4) + dC20_nh;
  C22_model = k_f q_r/4 + dC22_nh.
- The POSTERIOR of (dC20_nh, dC22_nh) IS "the degree of
  non-hydrostaticity implied" (user's ask) — report it in absolute
  units, in units of the Clipper sigmas, and as the implied
  −C20/C22 ratio vs the hydrostatic 10/3.
- Identifiability (review-corrected): the hydrostatic pair moves only
  along the 1-D k_f(C̄) manifold, but the free 2-D nuisances span BOTH
  directions of (C20, C22)-space, so after marginalization the gravity
  pair carries essentially NO interior information (see "Honest physics
  trade" above). Only the RATIO-BREAKING combination
  dC22_nh + (3/10)dC20_nh is identifiable as non-hydrostatic; the
  ratio-preserving combination is exactly degenerate with a C̄ shift.
  Report the dC22_nh–C̄ posterior correlation explicitly. Rejected
  alternatives (documented): retaining CMR2 alongside C20/C22
  double-counts (the CMR2 datum IS RD-derived from C22); a single
  ratio-breaking nuisance unphysically forces the other non-hydrostatic
  combination to zero.

## Config changes (Machine B authors trivially — copy v3 8D json)

`europa_clipper_v4_geodesy_10D.json` = v3 config with:

1. `Re_k2`: [0.25, 0.015]  (central kept at 0.25 = mid of Mazarico's
   0.2–0.3 and the ocean-bearing predictions 0.24–0.26; sigma = mid of
   (1.4–1.8)e-2)
2. `Im_k2`: [0.0, 0.015]  — **FLAG: placeholder.** Mazarico Table 5
   projects amplitude; the phase-uncertainty mapping
   sigma(Im k2) ≈ |k2|·sigma(phase, rad) needs the user's read of the
   table. Amplitude-sigma parity is the interim choice; user confirms
   or supplies the phase number before training.
3. REMOVE `CMR2` observable. ADD:
   - `C20`: [<fiducial C20_h>, 6.0e-7]  (sigma = mid of (3.5–8.5)e-7)
   - `C22`: [<fiducial C22_h>, 2.0e-7]  (sigma = mid of (1.5–2.3)e-7)
   **R1 GATE (blocking, before ANY freeze):** confirm whether Mazarico
   Table 5 sigmas are UNNORMALIZED or fully-normalized (4-pi) harmonic
   coefficients. The forward map above is unnormalized (C22_h = 130e-6,
   Anderson convention). If Table 5 is normalized: N20 = sqrt(5),
   N22 ≈ 0.6455 — sigma(C̄22) = 2e-7 normalized is 1.29e-7 unnormalized
   (35% likelihood-width error if ignored). Record the resolved
   convention + conversion in config metadata. USER: please check the
   table caption (Machine A could not access the full paper).
   Central values COMPUTED at the v3 fiducial node (Tb 264.5 K,
   w 35.165 ppt) from the new forward-model code with
   dC20_nh = dC22_nh = 0 — i.e., self-consistent hydrostatic fiducial,
   NOT Anderson's numbers. (Conditioning on a hydrostatic fiducial with
   zero-centered nuisances makes the null test clean: recovered
   non-hydrostaticity should be consistent with 0.)
4. `param_space` += `dC20_nh`, `dC22_nh` (uniform [−2e-5, 2e-5]) →
   10 sampled parameters.
5. Everything else (8 v3 params, induction observables at 1.5 nT,
   support cut, 2D cache path, GSW caveat metadata) UNCHANGED. Fresh
   seeds: train 43 / data 49 / noise 4949.

## Code (Machine A authors BEFORE B starts; cheap, cache-based)

- `PlanetProfile/Inference/gravity_obs.py` (new, pure numpy):
  `clairaut_kf(r_m, rho)` (Radau-form ODE integration over the
  per-sample composite density profile — the SAME profile the CMR2
  recompute uses, sampled core included), `radau_darwin_kf(cmr2)`
  (cross-check + GUI inverse only), `hydrostatic_c20_c22(kf, omega,
  R_m, M_kg)`, `cmr2_from_c22_rd(c22, omega, R_m, M_kg)` (GUI readout
  inverse). Unit tests: Anderson round-trip (C̄ 0.346 → 130e-6/434e-6);
  Clairaut-vs-RD agreement to ~1% on a realistic profile; analytic
  two-layer k_f solution; homogeneous-body limit (k_f = 3/2,
  C̄ = 0.4); inverse consistency; k_f monotonicity.
- `mcmc_runner`: observable channels `C20`, `C22` — per-sample CMR2
  (existing path incl. sampled core) → hydrostatic pair + nuisance
  offsets from theta. Nuisance params registered in
  `parameter_registry.py` (category 'gravity', no structure rebuild).
- Both runners package `c20_results`/`c22_results` (n,) per-sample +
  hydrostatic CMR2_h per sample for the GUI panels (same pattern as
  k2_results).
- MCMC mode in the GUI gains the option IMMEDIATELY (config-driven, no
  retrain); amortized waits for B's artifact.
- GUI (Machine A, after artifact lands): observable toggle CMR2 vs
  C20+C22 in the observables panel; results panel shows (a) posterior
  dC20_nh/dC22_nh ("non-hydrostaticity implied", with 0 marked), (b)
  RD-implied hydrostatic C/MR2 from the OBSERVED C22
  (`cmr2_from_c22_rd`) vs the posterior CMR2 distribution, (c) the
  posterior −C20/C22 ratio vs 10/3.

## Machine B execution list (the compute; single campaign)

1. Pull; verify Machine A's gravity_obs code + unit tests green.
2. Author the v4 config (above; fiducial C20/C22 from one
   `hydrostatic_c20_c22` call at the fiducial node — trivial).
3. Reference MCMC (nautilus/pocoMC as v3) on the v4 config; commit
   pickle to Test/mcmc_results/Europa/Test52_seawater_v3/ (same
   permission grant) or a Test53_geodesy_v4/ sibling dir (preferred;
   add the two .gitignore negation lines as in f7a03572).
4. 1M dataset (SAME 2D cache — no PP runs; support guard ON,
   drop_nonfinite) → nsf train → artifact
   `europa_clipper_v4_geodesy_10D_posterior_1m.pt`.
5. Gates: SBC (10 params), crosscheck, 2D Tb–w degeneracy gate (same
   pre-registration as v3), PLUS a **non-hydrostaticity recovery gate**
   (pre-registered per review R4, the central v4 deliverable). Two
   injection arms with DIFFERENT pre-registered expectations:
   - **Ratio-breaking arm** (identifiable): inject
     dC22_nh ∈ {0, ±5e-6, ±1e-5} with dC20_nh = 0 (and the mirrored
     dC20_nh-only set). Magnitudes stay ≥ one expected posterior width
     inside the ±2e-5 box (NOT +1.5e-5 = 75% of half-width — edge
     truncation biases medians; review refinement). Reference MCMC and
     flow must BOTH recover within |Δmedian| ≤ 0.25·sigma_post + W1.
     Negative injections included for asymmetry check.
   - **Ratio-preserving arm** (EXPECTED-DEGENERATE): inject
     dC20_nh = −(10/3)·dC22_nh — this is the hydrostatic-manifold
     tangent, exactly degenerate with a C̄ shift. Pre-register that
     the flow should NOT localize it: expected outcome is a broad
     nuisance posterior with a correlated C̄ shift. Gate PASSES if the
     flow reports the degeneracy (correlation structure matches the
     reference MCMC); it FAILS if the flow falsely localizes.
   - **Interior-unbiasedness check (both arms):** C̄/rho_core/R_core
     posteriors must be unbiased under injection, and the dC22_nh–C̄
     posterior correlation is reported explicitly — recovering the
     offset while silently biasing C̄ would otherwise pass.
6. Also re-run the Tb–w degeneracy gate: the much tighter k2 sigma
   (0.06 → 0.015) will sharpen rheology/Tb marginals; expect the Tb–w
   ridge to persist (induction-driven) — report corr as before.
7. Report per-parameter constraint improvement v3 → v4 (the point of
   the exercise: what does Clipper-grade geodesy buy?).

## Scope guards

- v3 artifact stays as-is (candidate, ratification pending user) — v4
  does not retrofit it; INDEX gets a v4 row when B delivers.
- No change to fundamental scientific assumptions: hydrostatic theory
  enters only as the C20/C22 forward map, with non-hydrostaticity
  explicitly sampled, and the RD systematic documented.
- Cache untouched (schema v3.0 serves v4 unchanged).
- Open numeric items for USER before B trains:
  (a) **R1 normalization gate**: are Mazarico Table 5 C20/C22 sigmas
      unnormalized or 4-pi-normalized? (blocking; see config section);
  (b) Im_k2 sigma — does Table 5 project a k2 PHASE uncertainty?
      (sigma(Im k2) ≈ |k2|·sigma(phase, rad); interim placeholder
      0.015 = amplitude parity);
  (c) confirm sigma midpoints (k2 0.015, C20 6e-7, C22 2e-7) vs
      picking a specific Mazarico tracking scenario;
  (d) nuisance prior half-width ±2e-5 OK? (GC21 offset 9e-6 at
      0.36 sigma — consistent with hydrostatic; v4 is an upper-limit
      measurement, and the honest sensitivity floor is set by the C̄
      uncertainty ~5–9e-6, not by sigma(C22) = 2e-7);
  (e) accept the MoI-constraint trade? Dropping CMR2 for free-nuisance
      C20/C22 means gravity no longer constrains the interior — that
      IS the honest physics, but it changes what v4 delivers.
