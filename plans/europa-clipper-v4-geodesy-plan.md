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

`europa_clipper_v4_geodesy_11D.json` = v3 config with:

1. `Re_k2`: [0.23, 0.015]  (USER 2026-07-19: ocean-consistent central
   per Mazarico section 3.1.1 / Moore & Schubert 2000; sigma = mid of
   the projected (1.4–1.8)e-2)
2. `Im_k2`: [0.0040, 0.015]  (USER: central = 0.23 x sin(1 deg), the
   1-degree phase-lag equivalent from section 3.1.1; sigma = same as
   amplitude — the paper projects no separate phase uncertainty)
3. KEEP `CMR2` [0.3547, 0.0024] — relabeled in metadata as the
   GALILEO-DERIVED MoI PRIOR (GC21 Table 3; see USER ANSWERS (e)),
   not a Clipper observable. ADD:
   - `C20`: [<fiducial C20_h>, 8.5e-7]  (unnormalized; top of the
     projected (3.5–8.5)e-7 range — absorbs the normalization-
     convention risk, see USER ANSWERS (a))
   - `C22`: [<fiducial C22_h>, 2.0e-7]  (unnormalized; conservative
     under either convention reading)
   Convention: UNNORMALIZED coefficients at reference radius 1565 km
   (GC21), J2_h = 3.324 x C22_h (Tricarico 2014 rapid-rotation
   correction; classical 10/3 in docstrings only). Resolved
   convention + conversions recorded in config metadata (R1 CLOSED).
   Central values COMPUTED at the v3 fiducial node (Tb 264.5 K,
   w 35.165 ppt) from the new forward-model code with
   dC20_nh = dC22_nh = 0 — i.e., self-consistent hydrostatic fiducial,
   NOT Anderson's numbers. (Conditioning on a hydrostatic fiducial with
   zero-centered nuisances makes the null test clean: recovered
   non-hydrostaticity should be consistent with 0.)
4. `param_space` += `dC20_nh`, `dC22_nh` (uniform [−2e-5, 2e-5]) →
   11 sampled parameters (incl. the zeta split below).
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
   `europa_clipper_v4_geodesy_11D_posterior_1m.pt`.
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
## USER ANSWERS 2026-07-19 (all five open items RESOLVED)

(a) **Normalization**: Mazarico et al. (2023) write the potential in
    4-pi-normalized Stokes coefficients (Eq. 2, Kaula formalism), BUT
    the user's quoted Table 5 C22 central range (1250–1400)e-7 matches
    the UNNORMALIZED literature span exactly (Jacobson 125.0e-6 →
    GC21 138.42e-6; normalized C̄22 would be ~2037e-7). GC21 Table 2
    is explicitly "unnormalized ... reference radius 1565 km".
    **Decision: v4 works in UNNORMALIZED coefficients throughout**
    (Anderson/GC21 convention), documented in config metadata with the
    conversion factors (C20: x sqrt(5) = 2.236 normalized→...; C22:
    C22_unnorm = 0.6455 x C̄22_norm). Sigma convention risk, per
    direction: if Table 5's sigma(C22) 2e-7 is normalized, unnorm is
    1.29e-7 → using 2.0e-7 is CONSERVATIVE (wider). If sigma(C20)
    6e-7 is normalized, unnorm is 13.4e-7 → using 6e-7 would be
    2.2x too TIGHT; adopt sigma(C20) = 13.4e-7-equivalent risk bound:
    use 8.5e-7 (top of the projected range) and note the residual
    factor ~1.6 risk in metadata. Non-hydrostaticity is an upper-limit
    deliverable whose floor is set by the MoI prior (below), so this
    residual sigma ambiguity does not drive the science claim.
(b) **Im_k2 sigma = 0.015** (same as amplitude; both accessible paper
    versions project amplitude only — no separate phase uncertainty).
(c) **Centrals: Re_k2 = 0.23** (ocean-consistent, Moore & Schubert
    2000 per the paper: "around 0.23 with an ocean and <0.015
    without"); **Im_k2 = 0.23 x sin(1 deg) = 0.0040** (paper section
    3.1.1: phase lag < 1 deg for a non-dissipative mantle; partial
    melt raises it to several degrees — 1 deg is the adopted
    conditioning value, user-directed).
(d) Nuisance-parameter framing explained to user (an offset sampled
    alongside the physical parameters so real non-hydrostatic power
    cannot masquerade as interior structure); design unchanged.
(e) **MoI bound RESTORED as a Galileo-derived prior (user-directed):**
    keep the existing CMR2 term [0.3547, 0.0024] — this IS the GC21
    MoI (their Table 3, "This Work") — but RELABELED in config
    metadata as a PRIOR derived from Galileo-era data (GC21's C22 via
    Radau + hydrostatic assumption), NOT a Clipper observable. No
    double-counting: it is independent of the future Clipper C20/C22
    measurement the new observables represent. This fully restores
    the interior MoI constraint the review flagged as lost (R3) while
    keeping the non-hydrostaticity readout honest. Config knob
    `moi_prior` = {source: 'gomescasajus2021 Table 3', value: 0.3547,
    sigma: 0.0024} so future users can swap Anderson (0.3475
    +/- 0.0026) or Jacobson (0.3405 +/- 0.0022) solutions.

**Additional physics refinement from GC21 (adopt):** Europa's
relatively rapid rotation modifies the hydrostatic ratio J2/C22 from
10/3 to **3.324** (Tricarico 2014); GC21 applied this in their MoI
derivation. v4 forward model uses J2_h = 3.324 x C22_h (the 10/3
value stays only in docstrings as the classical limit). The
ratio-preserving injection arm and the reported ratio diagnostic use
3.324 accordingly. Reference radius: coefficients are quoted at
1565 km (GC21) — the forward model must evaluate q_r and the
coefficients at the SAME reference radius, not R_body_m = 1560.8 km
(a (1565/1560.8)^2 ~ 0.5% effect on C22 = 3x sigma(C22): load-bearing,
document it).

## Galileo v1.1 re-run + Europa slot retirement (user 2026-07-19)

User finding: ALL Europa GUI artifacts trained with Titan-derived k2
constraints — sigma(k2) = 0.06 is the Cassini-Titan scale
(Titan: 0.608 +/- 0.048, 0.135 +/- 0.035), but Europa has NO measured
k2 whatsoever; the Galileo v1 artifact conditions on a fake
measurement (Re_k2 0.25 +/- 0.06). User directive: re-run Galileo,
then remove Clipper v2 in favor of v4.

**Galileo v1.1 spec (Machine B; reuses the Test51 1D cache unchanged):**
- Honest Galileo-era data: CMR2 [0.3547, 0.0024] (GC21 MoI, relabeled
  Galileo-derived as in v4) + the synodic |Ae| > 0.7 support cut.
  Galileo measured neither k2 nor h2.
- k2/h2 channels: RETAIN as HYPOTHETICAL-CONDITIONING channels with
  theory-motivated widths — Re_k2 [0.23, 0.05] (spans the 0.2-0.3
  ocean-consistent theory range), Im_k2 [0.004, 0.05], Re_h2
  [1.2, 0.1], Im_h2 [0.0, 0.1] — labeled in config metadata and the
  GUI scope note as "no Galileo-era measurement exists; these channels
  express hypothetical tidal data for exploration, defaults at the
  ocean-consistent theory point". Alternative (REJECTED, documented):
  dropping k2/h2 entirely is the strictly-honest Galileo dataset but
  leaves 7 params constrained by one observable + a support cut — the
  posterior is prior-dominated and the GUI demo loses its purpose;
  the wide-sigma retain keeps the tool useful with honest labeling.
  USER may override to the drop-channels variant before B trains.
- Everything else as v1 (7 params, Test51 cache, synodic-only).
  Fresh seeds: train 44 / data 50 / noise 5050. Full gate suite
  (SBC/crosscheck vs a re-run reference MCMC with the SAME corrected
  observables/limits grid-walk).

**Retirement schedule (user-ratified sequence):**
1. Machine B trains Galileo v1.1 -> gates -> Machine A slot swap
   (v1 slot replaced in place; v1 artifact -> INDEX retired row).
2. v4 trains -> gates -> Machine A GUI -> user ratification ->
   **Clipper v2 slot REMOVED** in the same commit that adds the v4
   slot (v2 -> INDEX retired row; artifact retained for provenance).
3. Until each replacement lands, the current v1/v2 slots stay live
   (public-app continuity) — their scope notes already carry validity
   caveats; do NOT deploy to HF between veto and replacement unless
   the user asks.

## Future work (user 2026-07-19; out of v4 scope)

Higher-degree gravity coefficients (l >= 3) within modeled error
bounds — as with directional Bind components, this requires a
simulation of the actual Europa Clipper trajectory/flyby geometry
(per-flyby sensitivity), not just projected sigmas. Same
infrastructure slot as the future Clipper-trajectory Bind inversion.

## Zeta split: independent ice / silicate transient creep (user 2026-07-19)

A SINGLE sampled `log10_zeta` couples the Andrade transient-creep
amplification of ice and rock: any zeta the sampler raises to fit ice
dissipation simultaneously amplifies silicate dissipation —
systematically preferencing silicate heating in the posterior heating
budget. Both new configs therefore replace `log10_zeta` with
`log10_zeta_Ih` + `log10_zeta_sil` (each U[-3, 2], the previous single
range; no HP-ice zeta — the Europa Tb box has no III/V/VI layers):

- `europa_galileo_v1p1_8D.json` (was 7D)
- `europa_clipper_v4_geodesy_11D.json` (was 10D)

The forward-model hook `apply_andrade_params` already supported
per-phase zeta (override precedence) — config-only change, verified on
Machine A: both configs load, prior draws yield finite likelihoods
(14/20 and 13/30), and perturbing each zeta independently moves k2
(zeta_sil -1.5 dex quintuples |Im k2| at the test draw — the
silicate-dominance mechanism, confirmed). Machine B's fiducial
C20/C22 remain VALID (gravity depends on the density profile, not
zeta); the fiducial-metadata interior median simply reads as the
single-zeta value for both phases. SBC/crosscheck now cover 8 and 11
params respectively; expect the two zetas to be individually less
identified than the old single zeta (they share the Im k2 budget) —
report their joint posterior (corr) alongside the marginals, and
pre-register a heating-fraction comparison (ice vs silicate) between
the split and single-zeta runs as the user-facing deliverable.
