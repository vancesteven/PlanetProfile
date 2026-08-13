# Enceladus Cassini campaign — config freeze (Machine A, opened 2026-08-11)

The last missing 1.0 mission pair (STRATEGY, user-ratified). Freeze
process: this document collects the requirements → scientific-reviewer
design adjudication → frozen config JSON committed → Machine B §0.18
Phase 2+ executes (cache → 1M dataset → B3 references → training at the
Phase 4 hold point under the preregistered release rule).

## Requirements

1. **USER DIRECTIVE (2026-08-11): account for the Cassini
   magnetic-induction results of Saur et al. 2024** (PSJ,
   "Analysis of Enceladus's Time-variable Space Environment to
   Magnetically Sound its Interior"; doi:10.3847/PSJ/ad8130):
   - Induction responses expected at the ORBITAL and SYNODIC periods,
     amplitudes up to a few nT; north-polar flybys most diagnostic;
     ~1 nT accuracy constrains bulk ocean conductivity (hence salinity).
   - Freeze decisions needed (reviewer adjudicates): (a) condition on
     Saur-derived induction observables (Bind channels per frequency,
     the Europa-Clipper pattern) vs training-time support cut vs
     display-reference only — driven by what amplitude+uncertainty the
     paper actually reports as a MEASUREMENT vs an upper limit;
     (b) source the paper: PDF not yet in papers/ — obtain before the
     design pass (user to drop in papers/ or fetch);
     (c) conductivity-salinity mapping for the chosen ocean
     composition(s) at Enceladus P,T (existing PP conductivity
     machinery; MoonMag first-party for the response function).
   - Follow-on literature to check at design time: Wivell et al. 2026
     (stratified-ocean inductive response) for whether a 1D bulk
     conductivity is adequate at current data precision.
2. **Salinity axis: evaluate 2D (Tb × log10 w) sampling** (program
   review 2026-08-11) — Enceladus has the program's strongest
   composition constraints (CDA/INMS plume chemistry) AND the induction
   directive makes salinity load-bearing. The July §2 spec's fixed
   10 ppt Seawater 1D cache is the fallback only if the 2D build is
   shown infeasible. Composition choice (Seawater vs explicit
   NaCl/MgSO4 per the Titan pattern) is a design-pass question;
   plume chemistry favors an alkaline Na-Cl-CO3 ocean.
3. **Libration**: the July §2 spec's primary decoupling channel —
   retained (physical libration amplitude constrains shell decoupling
   and thickness independently of k2).
4. **Gravity**: degree-2 (J2/C22 Iess et al. 2014 class) + the
   free-gravity non-hydrostatic-offset pattern (v6/Titan precedent);
   CMR2 dropped per the double-count ruling.
5. **Tidal k2**: no measured Enceladus k2 — decide labeled-hypothetical
   channels (v1.1 precedent) vs omission; NO fake-measurement
   conditioning (v1/v2 retirement lesson).
6. **Architecture**: deployed nsf under split-status by default at the
   Phase 4 hold point unless Track 2a/2c + C16 verdicts land first
   (§0.18 preregistered release rule). No architecture experimentation
   in production training (§0.16(c)).
7. **Gates**: full preregistered battery (SBC n=1500, crosscheck vs
   pooled B3 3-seed reference, limits under the ratified in-support
   rules, PPC/pushforward four-way + 0.5 sigma_obs flag) + the IS
   correction validated per §0.17 machinery once trained.

## Status

- Design pass COMPLETE (adjudication below, 2026-08-12; Saur 2024 PDF in
  papers/). Config JSON `not implemented` — blocked on B1–B7 closing and
  the user rulings.
- USER RULINGS (2026-08-12, complete): point 1 RATIFIED — k2 omitted.
  Point 2 OVERRIDDEN by user: **KEEP the no-ocean cases** (joint
  support retained against the reviewer's ocean-only recommendation —
  the frozen branch stays in the posterior; design question dispatched:
  how no-ocean nodes enter a zb-parameterized cache, likely zb extended
  to the freeze-through boundary per w, giving a clean monotone
  no-ocean edge instead of the Tb ribbon). Point 3 RATIFIED: **zb_km
  axis** ("there should be a mapping from thickness to unique Tb" —
  yes: zb↔Tb is monotone per (w, structure); Tb solved and recorded per
  node). **Seawater RATIFIED.** **"Include libration as a parameter"**
  — interpreted as: libration is a conditioned observable AND a
  per-sample derived quantity surfaced in the GUI results (corner plot
  column + summary row), not a free sampled parameter (it is fully
  determined by structure; sampling it would break the forward model).
  Flag to user if a different reading was intended. NEW DIRECTIVE: "omit
  hypothetical k2 but fit shape model as per Hemingway and Mittal 2019"
  — shape becomes a conditioned/fitted element of the campaign. Paper
  fetched to papers/hemingway2019enceladus.pdf (Icarus 332: joint
  shape+gravity+libration isostasy; shell mean 19-24 km, 4-12 km south
  pole; ocean 30-39 km; core 191-198 km). Reviewer design pass on the
  shape channel DISPATCHED (how to honor H&M within the 1D model:
  conditioned degree-2 shape + nh offset nuisances, isostasy-coupled
  offsets, display reference, or escalate as new physics; supersedes
  the earlier "shape — out" ruling).

## Reviewer design adjudication (2026-08-12, after reading Saur 2024)

Full review in the session record; binding rulings:

**Observables (frozen): C20/C22 + libration ONLY.**
- Gravity: **Park et al. 2024 Table 8 Case 2 (erratum-checked)** — C20
  -5477.45e-6 ± 36.99e-6, C22 +1517.90e-6 ± 14.70e-6, corr 0.47,
  R_ref 256.6 km, j2_over_c22 = 3.25 (McKinnon 2015). SUPERSEDES
  Iess 2014 + Thomas 2016. dC20_nh/dC22_nh nuisance boxes retained.
  CMR2 dropped (hydrostatic double-count). Honesty note: gravity alone
  is weak here (hydrostatic C22 signal ~ nuisance box) — the interior
  information is in the libration.
- Libration: **0.092 ± 0.003 deg (Park 2024 erratum; revises Thomas
  0.120)**. Implemented (Librations.py, Van Hoolst 2008 ocean+rigid;
  wired through mcmc_runner + SBI x-vector). σ = 3.3% relative;
  ~1.1σ per km of shell — the campaign's dominant channel AND its
  single point of failure: ZERO test coverage, an unexplained 6.3σ
  fiducial offset (likely benign: fiducial shell ~21.5 km vs Park's
  27 km — never demonstrated), no systematic budget.
- **k2 OMITTED entirely** (rejects the constraints-doc hypothetical
  channel): three real channels exist so the v1.1 justification does
  not transfer; a [0.015 ± 0.02] box straddling zero carries ~no
  information and re-creates fake-measurement optics; and it would
  couple the last 1.0 campaign to the open Track 1 validation.
  ** USER RATIFICATION SOUGHT (point 1).**

**Composition: Seawater (TEOS-10/GSW), log10 w U[-1, 2].** NaCl REJECTED
for 1.0 on implementation grounds: the repo's NaCl conductivity (Pan
2021) is fitted for P in [212, 1713] MPa; Enceladus ocean is 0.5-7 MPa —
a 30-400x extrapolation that would invalidate the conductivity display
deliverable. GSW is IN-envelope at Enceladus P,T; the w>42 ppt
extrapolation carries the measured <=18.5% systematic (Europa v3 Gate 2).
Reviewer verified numerically: plume band (5-20 ppt) gives sigma
0.45-1.79 S/m, reproducing Saur Fig. 20 "<2 S/m" exactly (committed as
acceptance test B10); Saur's ~5 S/m ceiling lands at w ~ 70 ppt; keep
100 ppt bound so the excluded region is DISPLAYED as excluded.
Explicit NaCl / Na-Cl-CO3 (Reaktoro) campaigns are 2.0.

**Induction: display reference ONLY — no conditioning AND no support
cut** (new ruling: Saur's two scenarios jointly span the whole salinity
axis, so a sigma_ocean cut would covertly assert "the core is not
conductive" — inverting Saur's own preferred branch). Decisive fact:
the inducing-field model scatter (2.4 nT in Bz) EXCEEDS the candidate
signal (1-3 nT); no published amplitude+uncertainty exists. ORBITAL
period is primary (5.0 nT x A~0.47 ~ 2.4 nT) over synodic (1-2 nT x
A~0.55 ~ 0.8 nT) — my earlier ordering corrected. Deliverables:
per-sample sigma_ocean, A and Phi at orbital+synodic, |B_ind| nT vs
Saur's 1-3 nT band (v6 published-CMR2 display pattern). Core
conductivity: PP's self-consistent porous-rock sigma as one curve + an
IMPOSED 25 S/m Saur-scenario overlay as a second — NOT a sampled
nuisance (degenerate vs amplitude), NOT a config parameter. Caption
must state all six mandated items (incl. the amplitude-convention trap:
Saur's A through a 5-10 km polar crust vs our mean-shell A ~25% lower,
and the 256.6-vs-252.1 km radius conventions). Upgrade-to-conditioned
criteria preregistered in config metadata (published amplitude+
uncertainty at a named period; sigma(B_ind) < 0.5|B_ind|; plasma
systematic marginalized). Wivell 2026 stratification: CLOSED
non-blocking (Saur §4: single-period sounding sees bulk averages only).

**Support: OCEAN-ONLY (retry_frozen_as_no_ocean=False).** The Titan
joint pattern does not transfer: plume salt grains are independent
observational evidence of liquid water (Postberg 2009/2011/2018) —
conditioning on independent data, not circularity. Frozen-model
exclusion reported as a one-evaluation chi^2 for the paper.
**USER RATIFICATION SOUGHT (point 2).**

**Cache geometry: (Tb x w) is NUMERICALLY INFEASIBLE at Enceladus —
sample zb_km (shell thickness) instead.** g = 0.1134 m/s^2 puts the
whole 5-40 km shell range inside 0.27 K of Tb (0.0077 K/km) while
salinity swings Tb by 6.1 K — the ocean region of a (Tb,w) rectangle is
a razor-thin tilted ribbon; libration precision (0.2 km class) would
demand ~66,000 nodes, <1% useful. Sampling zb_km in [5,45] x log10 w
gives ~3,200 nodes, ~100% useful, well-posed interpolation, and the
reportable coordinate (direct comparison to Park's 25-29 km).
PlanetProfile already supports zb via Bulk.zb_approximate_km; the cache
BUILDER needs a zb_grid mode (blocking item B4). The July Tb window was
~45% no-ocean. **USER RATIFICATION SOUGHT (point 3).**

**Blocking items before the config JSON commits (owners assigned):**
B1 libration reachability scan across the (zb,w) box (A + reviewer);
B2 libration systematic budget in sigma_obs units (A + reviewer)
   [UPDATE 2026-08-12: the rigid-vs-elastic "0.03%" figure justifying
   the rigid deferral is FALSIFIED — C18 escalation measured 0.82-0.95%
   with elasticity RAISING the amplitude, opposite the VH16 statement.
   Under adjudication; provisionally budget the elastic term at ~0.3
   sigma_obs symmetric until the ruling lands];
B3 libration regression test vs published value — Librations.py has no
test at all (Codex-suitable); B4 zb-axis cache-builder support (A);
B5 d(libration)/d(zb) scan to set node spacing (A); B6 drop rheology
params that map to no observable after k2 removal (A); B7 Ocean.deltaT
<= ~0.002 K near the melting curve verified (A). Induction-deliverable
blockers (parallel): B8 add synodic 15.559 hr excitation to Moments.py
+ Be row; B9 fix Enceladus Be-file resolution (GetBexc looks for
Be1xyz_Enceladus_Cassini_Cassini11noMP.txt — absent); B10 formalize the
Saur Fig. 20 conductivity acceptance test (already passing informally).
Non-blocking N1-N5 recorded in the review. July §2 supersession table
in the review (incl. mission-body-constraints.md:223-232 stale-line
fix; Cuncertainty re-derivation under zb).

Frozen-config sketch (params/priors/observables/metadata) is in the
review record; commits after B1-B7 close and the user rules on the
three ratification points.

## Shape-channel addendum ruling (reviewer, 2026-08-12, after reading H&M 2019)

**Central finding: H&M 2019 fit gravity (C20, C22, C30) + libration ONLY;
shape enters as a forward-model INPUT** (observed surface relief sets the
top of the shell; the Airy-compensated basal relief is computed; the
resulting NON-hydrostatic gravity is predicted). Shape uncertainty enters
as model-covariance inflation (their Eq. 22), not as a residual.

**RULING: adopt the H&M isostasy coupling (option b).**
- Conditioning on shape coefficients directly (a) REJECTED: hydrostatic
  part double-counts the same Tricarico relation as gravity (the CMR2/C22
  argument); non-hydrostatic part with free nuisances is pure absorption
  (observed shape is 4.3-84 sigma super-hydrostatic; geoid only 8.5
  sigma — the 4.09-vs-3.42 ratio gap IS the measurement, and only the
  coupling reads it).
- dC20_nh/dC22_nh free boxes REMOVED — replaced by ONE physical nuisance
  `compensation_C2 ~ U[0,1]` (Airy fraction) + H&M Sigma_model added
  variances. This converts gravity from nearly-uninformative to a real
  shell-thickness channel: the Airy cancellation residual scales as
  (R_b/R_t)^(l+2) — ~4.2%/km, ~0.5 sigma/km on C20 (vs libration's
  ~1.1 sigma/km; complementary).
- **ADD C30 = [1.7782e-4, 3.342e-5] (Park) as a conditioned observable**
  — 100% non-hydrostatic, H&M's cleanest shell-thickness channel;
  reverses the constraints-doc "out of 1D scope" tag (that assumed no
  non-hydrostatic forward model). S22 stays out (needs lateral physics).
- Shape source: **Tajeddine 2017 primary** (sigma_model on C20 5.3e-6 —
  negligible vs gravity sigma 35.4e-6) — Nimmo 2011 would inflate
  sigma(C20) 7x and kill the channel; run Nimmo as a preregistered
  sensitivity ablation with H&M's own comparability caveat quoted.
- rho_ice / rho_ocean promoted to sampled params (Airy amplification
  rho_ice/delta_rho ~ 9.7: 2% density error = 20% root-amplitude error).
- Finite-amplitude terrain correction (Wieczorek & Phillips 1998)
  REQUIRED (~5% on J2 = 8 sigma if omitted).
- Reference-radius trap: Park gravity at 256.6 km vs Tajeddine shape at
  252.22 km — 3.6 sigma on C22 if mishandled; single unit-tested
  conversion (B14).
- Reference comparison: show BOTH H&M (shell 19-24 km) and Park (25-29
  km) bands, with the stated dependence caveat (we hybridize Park
  gravity+libration with Tajeddine shape).
- Libration systematic budget B2' EXTENDED: H&M compute libration from
  the observed NON-hydrostatic triaxial figure; Librations.py uses the
  hydrostatic figure — plausibly multi-sigma at sigma=3.3%; quantify
  (B-A)/C difference. Tajeddine-propagated libration model sigma
  0.00025 deg is fine; Nimmo's 0.004 deg would exceed the measurement
  sigma.

**New blocking items:** B11 isostatic-gravity module (H&M Eqs. 8/9/12,
degree 2+3, ~300-500 lines — NEW forward-model code on the 1.0 critical
path; schedule risk surfaced to user); B12 finite-amplitude correction;
**B13 H&M REPRODUCTION GATE** (with Iess+Thomas+Tajeddine inputs our
module must reproduce H&M's Airy solution: shell 19-24 km, ocean 35-39,
core 192-195, rho_core 2340-2410 — fallback on failure is display-only +
restored nuisance boxes); B14 radius reconciliation; B15 Sigma_model;
B16 per-interface hydrostatic-figure verification (1st-order theory
underestimates H22_hyd ~4% -> ~30% error on H22_nh; must be the full
Tricarico numeric).

**OPEN TENSION (reviewer follow-up pass queued):** the Airy coupling
presupposes an ice/ocean interface — undefined on the frozen branch —
while the USER has ruled the no-ocean cases STAY. Candidate resolution
for the reviewer to adjudicate: branch-dependent support model — ocean
nodes use the sampled-C2 Airy coupling; frozen nodes use rigid
(uncompensated) support, whose predicted non-hydrostatic gravity (~12
sigma off observed) and near-zero libration decoupling make the frozen
branch DOUBLY discriminated by real physics rather than assumption.
Also pending in that pass: zb-axis no-ocean edge design and the
"libration as a parameter" interpretation check.
