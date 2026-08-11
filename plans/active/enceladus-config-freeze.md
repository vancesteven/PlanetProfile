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

- `not implemented` — requirements collected; reviewer design pass not
  yet dispatched. Blockers: Saur 2024 PDF not in papers/; induction
  observable-vs-support decision needs the paper's measured values.
