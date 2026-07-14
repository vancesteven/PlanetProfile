# Europa Clipper-era v2: induction observables in measurement space (Machine B)

Author: Machine A, 2026-07-13. Scientific review (opus): RATIFY-WITH-EDITS
2026-07-13; all six required edits are incorporated below. Status of every
item below: `not implemented` unless marked otherwise. Machine B (opus 4.8, env PPcl) executes all
compute-intensive phases; Machine A authored this plan and will do GUI
integration + verification after artifacts land.

## Objective

Train and deploy `europa_seawater_andrade_clipper_v2` — a 3-frequency SBI
artifact whose induction observables are the per-component complex
**induced dipole coefficients expressed as equivalent surface field**:

```
Bind_comp(f) = Ae(f) * Be_comp(f)        [nT, complex]
```

NOT a local field measurement at a flyby position (that would carry the
r^-3 dipole geometry factor and spacecraft trajectory). Clipper's
data-reduction inverts flyby data to recover exactly these induced-dipole
coefficients; the Kivelson et al. 2023 1.5 nT requirement is the
post-inversion precision on Re and Im of each component, which is what we
condition on. Assumptions, stated explicitly:

- **Spherically symmetric ocean**: a scalar complex Ae per frequency is
  the exact response only for a radially layered conductor; lateral
  ocean/ionosphere structure would need a tensor/multipole response.
- **Normalization consistency**: Ae and Be must share the same n=1
  Gauss/Schmidt normalization and reference radius; Phase 1 must verify
  the convention (factor-of-2 pitfalls) so the product is meaningful.
- 1.5 nT is an instrument+inversion floor: plume, ionospheric, and plasma
  current contamination are NOT in the forward model. The observable is
  the *cleaned* induced-dipole coefficient; sub-nT differences must not be
  over-interpreted. Errors across components/Re-Im are treated as
  independent diagonal Gaussians — the flyby inversion actually correlates
  them; carrying a covariance is a documented future upgrade, not v2.

This replaces conditioning on dimensionless Ae: a flat 1.5 nT noise floor
maps to a *different* effective Ae precision per frequency/component
(sigma_Ae,eff = 1.5 nT / |Be_comp|), which only measurement-space training
represents correctly.

The v1 Galileo-run artifact (synodic-only support cut, no induction
channels in x) stays deployed untouched; v2 is a separate slot.

## User-ratified constraints carried over from v1

- Galileo synodic bound stays as a training-time support cut:
  `induction_bounds: {synodic: {amp_min: 0.7, im_abs_max: 0.4}}`
  (A > 0.7 one-sided, |Im A| < 0.4; do NOT tighten to literature 0.92/0.97).
- k2/h2 observables per Mazarico et al. 2023 assumptions (same channels and
  sigmas as `europa_seawater_andrade_7D.json`).
- CMR2 0.3547 +/- 0.0024 Gaussian, core-sensitive mass-conservation
  derivation (derived_params block identical to v1).
- Same 7D parameter space and priors as v1 (alpha, log10_zeta, log10_eta_Ih,
  log10_eta_sil, Tb_K, R_core_km, rho_core_kgm3). Priors frozen at training.

## Phase 0 — BLOCKER: Ae solver magnitude reconciliation — **VERIFIED (Machine B, 2026-07-14); UNBLOCKED for n=1**

**Status: `verified` for the n=1 dipole physics v2 uses** (scientific-reviewer
opus PASS, 2026-07-14). Gate PASSES: production MoonMag Srivastava solver
(`InducedAeList`/`AeResponse` — the exact path `forward_model_induction` uses to
fill the SBI cache) vs an INDEPENDENT from-scratch mpmath (50-digit) log-derivative
propagation solver agree to **worst 3.7e-14 relative amplitude and 0.0 deg phase**
across all 9 (model x frequency) points — far inside the 5% / 5 deg acceptance.
Script + results committed: `plans/scripts/phase0_ae_reconciliation.py`,
`plans/scripts/phase0_ae_reconciliation_results.json`.

3 models (valid conductive-ocean, Tb>=261.5, Test51_seawater cache idx 10/23/35):

| model | Tb (K) | D_ocean (km) | freq | T (hr) | \|Ae\| prod | \|Ae\| indep | d-amp % | d-phase deg |
|-------|--------|--------------|------|--------|-----------|------------|---------|-------------|
| THIN  | 261.5  | 11.2  | synodic     | 11.23 | 0.7536 | 0.7536 | 0 | 0 |
| THIN  | 261.5  | 11.2  | synodic 2nd | 5.62  | 0.8149 | 0.8149 | 0 | 0 |
| THIN  | 261.5  | 11.2  | orbital     | 85.24 | 0.2245 | 0.2245 | 0 | 0 |
| MID   | 268.0  | 67.9  | synodic     | 11.23 | 0.8969 | 0.8969 | 0 | 0 |
| MID   | 268.0  | 67.9  | synodic 2nd | 5.62  | 0.9016 | 0.9016 | 0 | 0 |
| MID   | 268.0  | 67.9  | orbital     | 85.24 | 0.8165 | 0.8165 | 0 | 0 |
| THICK | 271.0  | 106.3 | synodic     | 11.23 | 0.9390 | 0.9390 | 0 | 0 |
| THICK | 271.0  | 106.3 | synodic 2nd | 5.62  | 0.9535 | 0.9535 | 0 | 0 |
| THICK | 271.0  | 106.3 | orbital     | 85.24 | 0.8961 | 0.8961 | 0 | 0 |

Independence is genuine (not tautological): the two solvers share only the
spherical-Bessel basis; production uses Srivastava beta/gamma/delta/epsilon transfer
coefficients with a Lambda-recursion, the independent solver propagates the Riccati
log-derivative L=P'/P and matches to the outer vacuum potential. Confirmed by
(a) an ocean-sigma perturbation sweep (x0.5..x5) where |Ae| genuinely moves and both
solvers still track to ~1e-16 — an input echo would not co-vary; (b) analytic anchors
pinning ABSOLUTE normalization: vacuum -> Ae=0, PEC sphere filling body -> |Ae|=1;
(c) a discriminating PEC-sphere-at-a!=R test giving exactly (a/R)^3 in both — this is
the real geometry, since the cached models place the ionosphere boundary at
r_out=1,660,800 m, 100 km ABOVE R_body=1,560,800 m (rscaling=1.064 != 1 IS exercised).

**|Ae_orbital| red flag DEBUNKED.** Re-derived |Ae_orbital| = 0.22 (thin) / 0.82 (mid)
/ 0.90 (thick), rising monotonically with ocean-thickness / skin-depth ratio (85-hr
seawater skin depth ~162-183 km, comparable to the ocean thickness). This is the correct
shielding limit; the Callisto "O(1e-2) expected" anchor has NO physical basis for a
seawater ocean at these thicknesses. The 0.73 was not a solver bug — the expectation was
wrong. (Ae may be re-enabled on the physics; still DROP the Ae GUI channels for other
reasons per the C-B DROP verdict.)

**v1 synodic support cut retro-validated:** synodic |Ae| spans 0.754-0.939 across the
valid conductive Tb grid — all safely above the 0.7 cut.

**Scope caveat (reviewer, MODERATE, out-of-scope for v2):** reconciliation is validated
at **n=1 ONLY**. The two solvers use algebraically different degree-n surface rescaling
exponents that coincide only at n=1 (production `rscaling^(n+2)`; independent
`R^-(2n+1)`); they diverge 6.0% at n=2, 11.7% at n=3. Irrelevant to v2 (dipole only;
`forward_model_induction` defaults nn=1, Bind=Ae*Be is dipole), but any future n>=2
induction work MUST add an analytic n=2 (quadrupole PEC) anchor and resolve which
exponent is the correct general-n reference-radius rescaling before use. Do NOT claim
the reference solver is validated beyond the dipole.

--- original task spec below (satisfied) ---

The Callisto rebuild flagged |Ae_orbital| = 0.73 where O(1e-2) was expected.
Until the solver magnitude is reconciled, no Ae-derived quantity may enter
training data.

Task (Machine B): for 3 representative Europa models spanning the Tb grid
(thin ocean / mid / thick), compare the grid-cached Ae at synodic,
synodic 2nd, and orbital periods against an INDEPENDENT multilayer
spherical-induction solution fed the SAME layered sigma(r) profile
(recursive Eckhardt/Srivastava formulation; a uniform-sigma_mean shell is
NOT an adequate reference — induction is nonlinear in sigma and a real
profile effect could false-FAIL a 5% gate or hide a normalization bug).
Acceptance: |Ae_solver - Ae_reference| / |Ae_reference| < 5% in amplitude
AND < 5 deg in phase for ALL 9 (model x frequency) points. Record the
table in the handoff. FAIL -> stop, diagnose MoonMag/solver normalization
before any dataset generation. This gate also retro-validates the v1
support cut (which used synodic Ae only).

IMPORTANT: also RE-DERIVE the expected |Ae_orbital| with the trusted
reference before diagnosing the Callisto flag. At the 85 hr period the
seawater skin depth (~150 km at ~3 S/m) is comparable to Europa's ocean
thickness, so |Ae_orbital| of a few tenths is physically plausible — the
"expected O(1e-2)" anchor behind the |Ae_orbital| = 0.73 red flag is
itself suspect and must be checked, not assumed.

## Phase 1 — code: `Bind_` channel family — **VERIFIED (Machine B, 2026-07-14; commit 92eb61c3)**

**Status: `verified`.** Scientific-reviewer (opus) verdict PASS WITH CONCERNS:
Phase 1 code is correct and internally self-consistent (`Bind = Ae * Be_comp`,
unconjugated, matches PP's `Binm` line 207 + FT surface-field path line 982,
NOT the Zimmer-2000 `conj(Ae)` display convention). Abs-fold isolation,
`_parse_bind_channel` grammar, closest-period label mapping (with a 0.1-hr
margin guard), and rejection paths all confirmed sound.

Implemented: `_parse_bind_channel`, `_load_be_excitation` (closest-period
argmin + `BE_PERIOD_MATCH_TOL_HR=0.1` guard), likelihood chi² + SBI x-vector
`Bind_` branches, `_channel_conventions` artifact metadata (signed for
`Bind_`, global for k2/h2 Im). Verified: 13 unit tests (test_bind_channels.py,
gitignored per repo convention) + end-to-end `generate_sbi_dataset` on the
real Test51_seawater cache / mpmath Ae grid / Europa Be file — `Bind_` columns
finite, row value matches independent analytic `Ae*Be` to 1e-6, imag column
all-negative (no abs-folding).

### Phase-2 BLOCKERS surfaced by the Phase-1 review (external-data/GUI contract)

These do NOT affect Phase 1 mechanics (fiducial values are generated by the
same code, so train/validate is self-consistent regardless). They set the
Phase 2 fiducial VALUES and σ interpretation and MUST be resolved before the
config is frozen:

1. **Factor-of-2 (MODERATE).** Verify whether Kivelson et al. 2023's "1.5 nT
   on the induced-dipole coefficient" refers to (a) the Gauss/dipole
   coefficient g₁ = (n/(n+1))·Ae·Be = ½·Ae·Be for n=1, or (b) the peak/
   equatorial induced surface field = Ae·Be (what this code computes, = PP's
   `Bi1xyz`). If (a), insert the n/(n+1)=½ factor consistently in the
   likelihood, `_induction_channel_value`, AND the fiducial computation, and
   the σ=1.5 nT interpretation changes. Record the chosen definition in the
   artifact metadata. Kivelson 2023 not in papers/; check styczinski2023
   / vance2021magnetic (papers/) for PP's own convention, and the primary
   Kivelson reference for the "coefficient" definition + phase-sign.
2. **Phase-convention ingest (MODERATE) — CORRECTED 2026-07-14 (Phase-2
   review).** The earlier "negate Im to convert from Zimmer" note was WRONG
   for the complex per-component Be actually used. Because `Be` is complex,
   `Ae·Be` and `conj(Ae)·Be` differ in BOTH real and imaginary parts, not by
   a bare Im sign flip (verified at the fiducial: synodic-x `Ae·Be` =
   91.82−157.77j vs `conj(Ae)·Be` = 126.64−131.47j; real parts differ ~35 nT).
   The v2 convention is `Bind = Ae·Be` with complex Be, matching PP's FT
   surface-field path (`MagneticInduction.py:982`, `Be1xyzFT·Ae1FT`) — NOT
   the periodic display `Bi1xyz` (`:244`, `|Be|·conj(Ae)`), a third distinct
   quantity. Harmless internally (self-consistent). GUI contract: label
   `Bind = Ae·Be (unconjugated, complex Be, FT surface-field path)` and MUST
   NOT implement a naive "negate Im" ingest recipe — mapping an externally
   quoted coefficient requires that source's full complex definition
   (amplitude + phase reference). Documented in config metadata
   (`bind_convention_2026_07_14`); GUI input labels are a Machine A task.

### Phase-2 validation checks (from the review)

- Round-trip fiducial test: condition a zero-noise MCMC/SBI on the exact
  Phase-2 fiducial `Bind` values from the grid cache, confirm the posterior
  peaks at the generating Tb (proves internal self-consistency).
- (LOW, DONE in Phase 1) closest-period margin-guard unit test added.
- (Adjacent, pre-existing, out-of-scope) `sbi_runner._x_obs_vector` folds
  only `_IM_K2_ALIASES`, NOT `Im_h2`/`abs_Im_h2`. Since v2 carries h2
  channels, confirm the h2 inference-time fold path separately (predates
  this change; not a `Bind_` defect).

--- original Phase 1 task spec below (satisfied) ---

## Phase 1 — code: `Bind_` channel family (Machine B; small, non-intensive)

1. `mcmc_runner.py` `_induction_channel_value` + likelihood observable
   parsing: add the family
   `Bind_<label>_<comp>_real` / `Bind_<label>_<comp>_imag`, comp in
   {x, y, z}, label in canonical PP excitation names. Value =
   `(Ae(label) * Be_comp(label)).real / .imag` in nT, complex
   multiplication (this encodes the phase rotation of the excitation into
   the measured component — do NOT multiply amplitudes).
2. Be_comp source: `PlanetProfile/MagneticInduction/MoonMag/excitation/
   Be1xyz_Europa.txt` (complex Bex/Bey/Bez per row). Load once per runner
   init, keyed by exc name; coordinate frame is whatever the Be1xyz file
   uses (IAU components as produced by MoonMag) — document the frame string
   in the config metadata so the future GUI labels axes correctly.
3. Sign convention: these channels use SIGNED Im. They are not k2 aliases;
   the abs-fold in mcmc_runner is hardcoded to Im_k2/Im_h2 only, so
   dataset generation is already signed — but the ARTIFACT/GUI contract is
   not: the pre-deploy assertion requires imag_convention == 'abs' ("the
   only convention the GUI accepts"), and v2 is the first artifact with
   signed-Im channels in the x-vector. Required: per-channel convention
   metadata in the artifact (e.g. `channel_conventions: {Bind_*: 'signed',
   Im_k2: 'abs', ...}`), an updated pre-deploy assertion accepting
   signed-Im for the Bind family, a unit test that folding does not touch
   Bind channels, and a GUI AppTest confirming the app does not abs-fold
   Bind inputs.
4. Noise injection for SBI datasets: Gaussian sigma per channel from the
   config observables (1.5 nT), same generative model as the Gaussian
   likelihood. Recorded in rejection_stats.obs_noise as usual.
5. Unit tests (mirror the induction_bounds tests): channel value equals
   analytic complex product on a synthetic Ae grid; abs-folding untouched;
   NaN row-rejection when the label/cache is missing.

## Excitation constants (from Be1xyz_Europa.txt, for design + pruning)

| label        | period (hr) | \|Bex\| nT | \|Bey\| nT | \|Bez\| nT |
|--------------|------------|-----------|-----------|-----------|
| synodic      | 11.233     | 213.2     | 75.7      | 15.9      |
| synodic 2nd  | 5.617      | 17.5      | 11.7      | 1.7       |
| adjusted orbital | 85.213 | 10.4      | 3.6       | 0.4       |

Pre-registered channel pruning rule (AMENDED per scientific review
2026-07-13): prune on the SIGNAL, not the excitation — drop any component
whose max over the ocean-state (Tb-grid) range of |Ae(f)| * |Be_comp| is
< 3 nT (SNR < 2 at 1.5 nT). The excitation-only cut mis-ranks the orbital
channels, where |Ae| is smallest. The kept-channel set is therefore
PROVISIONAL until Phase 0 delivers trusted |Ae(f)| values; the table above
is the excitation-amplitude input to that computation, not the decision.
Provisional expectation: synodic {x,y,z} + synodic-2nd {x,y} + orbital
{x,y} = 7 components x (Re, Im) = 14 induction channels, minus whatever
the |Bind| cut removes (orbital-y at |Be| = 3.6 nT is the most exposed).
Full x dimension (if all 7 survive): 14 + CMR2 + Re/Im k2 + Re/Im h2 = 19.

Note 'true anomaly' (84.63 hr) has |Bez| = 11.8 nT — a real Clipper-era
signal. NOT included in v2 (user scoped synodic + synodic 2nd + orbital);
flag for a possible v3 if the user wants it.

## Phase 2 — config: `europa_seawater_andrade_clipper_v2.json` — **VERIFIED (Machine B, 2026-07-14)**

**Status: `verified`.** Config authored at
`PlanetProfile/Inference/configs/europa_seawater_andrade_clipper_v2.json`. v1
carry-over is byte-for-byte (params/priors/derived_params/induction_bounds + 5
v1 observables), confirmed by scientific-reviewer (opus). 14 Bind channels
added (7 kept components × Re/Im). Fiducial values computed via the PRODUCTION
`MCMCRunner` Ae grid + Be loader at the v1 Test51 posterior weighted-median
Tb = 264.52 K (nearest grid point 264.50 K). Pruning (pre-registered ≥3 nT SNR
over Tb ≥ 261.5 K) KEPT synodic {x,y,z}, synodic-2nd {x,y}, orbital {x,y};
DROPPED synodic-2nd-z (1.65 nT) and orbital-z (0.37 nT) — exactly the plan's
provisional expectation. **Round-trip self-consistency VERIFIED:** conditioning
the v2 likelihood on the exact fiducials minimizes the induction chi² at the
generating grid Tb (264.5 K, chi²=3e-9) with a clean monotonic well.
Scientific-reviewer verdict **RATIFY-WITH-EDITS** (2026-07-14); both edits
(phase-convention text MODERATE, line-244 citation MINOR) were
documentation-only in metadata and are incorporated — training was never
blocked. Factor-of-2 handled as a documented, revisitable ASSUMPTION
(surface-field `Ae·Be`); train/validate self-consistent under either
definition. Fiducial script + results: `/tmp/compute_v2_fiducials.py`,
`/tmp/v2_fiducials.json` (session-local, not committed — regenerable).

- Copy `europa_seawater_andrade_7D.json` (params, priors, derived_params,
  induction_bounds, cache path, template module).
- Add the 14 Bind channels. Central values: evaluate the forward model at
  the v1 posterior median (Tb ~ 266 K fiducial; Machine B computes the
  exact fiducial Bind values from the grid cache after Phase 0 passes) —
  these are placeholder conditioning defaults, the GUI exposes them as
  free values. Sigma = 1.5 for every Bind channel (Kivelson et al. 2023).
- Keep the 5 v1 channels unchanged.
- metadata: cite Kivelson et al. 2023 (1.5 nT), coordinate frame, pruning
  rule, plume/ionosphere + diagonal-covariance caveats, and this plan file.
- Galileo-cut vs Clipper-likelihood interaction (documented, not double
  counting): the hard |Ae_synodic| > 0.7 support cut is Galileo prior
  support; the Bind_synodic Gaussians are the Clipper likelihood — valid
  Bayesian updating when the conditioned values are fresh Clipper data.
  Conflict regime: conditioning Bind_synodic on values implying
  |Ae_synodic| < 0.7 makes the cut hard-reject what the likelihood pulls
  toward. The GUI must warn when the user's synodic Bind inputs imply
  |Ae| < 0.7 (consistency check via the excitation constants).

## Phase 3 — Machine B run recipe (intensive)

**Phase 3 status: candidate committed — `implemented, unverified` (Machine B, 2026-07-14).**
Artifact `sbi_artifacts/europa_seawater_andrade_clipper_v2.pt` (config_hash
46be64069a40090f, git_sha f59f43b0, seed 42, n_train 831,566, nsf, sbi 0.26.1,
torch 2.8.0). Reference MCMC (nautilus, r_hat 1.0, ESS 4259) at fiducial Tb 264.5 K.
Gates (reports in `sbi_artifacts/validation_reports/europa_clipper_v2_1m/`):
SBC PASS all 7 (min KS p .19); crosscheck PASS mean/median/sigma all 7, shape FAIL
Tb_K only (D .119 vs .054) — a genuine flow shape defect (near-Gaussian flow vs
skewed reference, confirmed by dual-reference MCMC D .035~floor .035; NOT the v1
edge-smear, matched-truncation is a no-op), sub-resolution + conservative; limits
grid-walk W1 7/8 PASS (medians agree <=.04 K over 261.5-271 K), 263.5 miss is the
same shape floor (MCMC-vs-MCMC W1 .0066 vs flow .0442), containment .989-.997 = the
v1 gate-measurement artifact. Reviewer (opus 2026-07-14): COMMIT-AS-CANDIDATE.
**Anchor design amended AGAIN (this session):** the pre-registered independent-per-
frequency sweep below (3.4) was REJECTED on review — Ae is a function of Tb ALONE, so
frequencies cannot be moved independently and off-manifold synodic targets rail the
reference MCMC. Replaced by single-Tb GRID-WALK anchors as the primary W1 gate; the
plan set was run as EXTRAPOLATION PROBES and empirically confirmed the rejection.
Scripts: `plans/scripts/phase3_*.py`. INDEX candidate row + scope note carry the
Machine-A deploy conditions. NOT deployed — awaits Machine A GUI slot + AppTest + user
ratification.

Same structure cache as v1 (`Test51_seawater/europa_seawater_structure_grid
.pkl`) — the Ae grid is derived per Tb from the same cache; no rebuild.
Seeds: train 42 / data 47 / noise 4747 (fresh data+noise seeds — v1 used
44/4444).

1. Reference MCMC on the v2 config (new channels -> old Test51 pickle is
   NOT a valid crosscheck reference). Production settings as Test51.
2. 1M dataset (`generate_sbi_dataset`, apply_support_guard=True,
   drop_nonfinite=True). Expect ~15-17% support rejection (v1: 16.8%).
3. nsf train, save artifact (schema v1 fields; torch/sbi versions recorded).
4. Gates (ratified amended rules):
   - SBC 1000 pairs, KS p >= 0.05 per param.
   - crosscheck vs the Phase-3.1 reference MCMC (mean/sigma + shape gate
     D <= 1.5x bootstrap self-D floor + |dmedian| <= 0.3 dex). Expect the
     same Tb edge-smear as v1; a matched-truncation crosscheck at
     Tb >= 261.5 K is the pre-registered remedy check.
   - limits/W1 anchors (AMENDED per scientific review 2026-07-13): a
     change in Ae(f) moves ALL components of that frequency jointly —
     anchors that move only one component sit off the physical manifold.
     Sweep design: (a) Ae_synodic in {0.70, 0.75, 0.80, 0.85, 0.90, 0.95},
     each anchor updating synodic x, y, AND z (Re, Im) jointly via the
     excitation constants; (b) a coarse independent sweep of
     Ae_synodic-2nd (3 anchors) and Ae_orbital (3 anchors spanning the
     Phase-0-reconciled range), each moving its own frequency's components
     jointly; (c) one joint-corner anchor with all three frequencies
     displaced simultaneously. Non-swept channels at fiducial. W1 <= 0.25
     sigma_anchor + containment, reject_outside_prior=False for the sweep
     only.
5. Commit artifact + validation reports + INDEX row (candidate section)
   and push. Machine A does GUI slot + AppTest verification + user
   ratification before deploy.

## GUI integration (Machine A, after artifacts land)

- Slot `europa_seawater_andrade_clipper_v2.pt` in `_SBI_ARTIFACT_SLOTS`
  with `config_path` (REQUIRED — see 2026-07-13 CMR2 dispatch bug, commit
  abbb2272), `x_obs_limits` from the anchor-validated domain,
  `default_truncate` Tb if the edge-smear reproduces, scope note.
- Observable inputs grouped per frequency with nT units and the 1.5 nT
  sigma read-only (frozen at training).
- Consistency warning when synodic Bind inputs imply |Ae_synodic| < 0.7
  (conflicts with the trained Galileo support cut — see Phase 2 note).
- AppTest asserting Bind inputs are NOT abs-folded (signed Im preserved
  end-to-end).
- Verification: AppTest click-through + CMR2/heating/k2 plot inspection
  per CLAUDE.md discipline before any `verified` status.

## Status vocabulary

Use `verified` / `implemented, unverified` / `not implemented` per
CLAUDE.md. No item may be layered on an unverified predecessor; Phase 0 is
a hard gate for everything downstream.
