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

## Phase 0 — BLOCKER: Ae solver magnitude reconciliation

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

## Phase 2 — config: `europa_seawater_andrade_clipper_v2.json`

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
