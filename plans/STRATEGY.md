# Program strategy map

Updated: 2026-08-06. The big-picture companion to `plans/STATUS.md`
(tactical). STATUS answers "what changed"; this answers "where each
production line stands, why current work exists, and what it buys."
Refresh whenever a production line changes state.

## The program

A public Bayesian-inference GUI (Hugging Face Space) for icy-moon
interiors: one ratified 1D interior model per mission–body pair. A
"production run" is a campaign: structure cache → reference MCMC →
~1M-simulation trained flow (SBI artifact) → validation gates →
ratification → GUI slot → deploy. Machine B trains; Machine A (manager)
validates, adjudicates, integrates; Codex handles scoped implementation.

## Production lines, one row each

| Line | State | Next production action |
|---|---|---|
| Europa Galileo v1.1 | **DEPLOYED, revalidated** (PPC clean 2026-08-04) | none |
| Europa Clipper v4 | **DEPLOYED, revalidated** (PPC clean) | none |
| Europa Clipper v5/v6 | trained, delivered, **ratification in progress** (user priority) | Machine B: B1/B2/B6 validation against the fresh references → manager re-adjudication → deploy |
| Europa Clipper v7 open-\|Ae\| | trained, delivered, blocked behind v5/v6 machinery (shares B1–B6) | after v5/v6 |
| Titan no-ocean (Test50 + freegrav Phase A) | DEPLOYED / delivered | none |
| Titan NH3 joint | **DEPLOYED with split ratification**: ocean-structure science (ice thickness, ocean, salinity, core) verified and usable; tidal/dissipation sector quarantined (reference MCMC authoritative) | tidal-sector remedy decision (below) |
| Titan MgSO4 / NaCl | **HELD — the next production campaigns**; the hold is the remedy decision | user decision this week (options below) |
| Enceladus Cassini | **queued since July, untouched by any current issue** | can start as soon as Machine B has capacity |
| Callisto, Ganymede | roadmap | after Enceladus |

## What the last 24 hours were (the k2 detour) and what they bought

Trigger: the user's broad-k2 concern → a posterior-predictive check found
the Titan NH3 flow effectively ignores its k2 data (its tidal posterior
sits at the prior). The production question: one-off quirk, or a systemic
pipeline defect that would poison every future campaign?

Answer, at near-zero simulation cost (every test reused existing data):

- The deployed Europa artifacts are CLEAN — the public app needed no
  warning, and the concern does not generalize backward.
- The NH3 defect is real but narrowly scoped. Eliminated as causes:
  observation-vector width, x-normalization, the noise convention, the
  joint no-ocean mixture, the salinity axis, training-set size, and
  reference-MCMC error. Remaining: flow capacity/embedding for this
  problem class (13 parameters, strong Tb–w degeneracy, weakly identified
  tidal sector).
- Side products now standard in production: the pushforward (PPC) gate
  that catches this defect class, matched-resolution reference standards
  (n_eff=2000; the old n_eff=500 Europa references wandered ~1 km), and
  an empirical shape-tolerance floor for the v5/v6 ratification.

## What the NH3 finding implies for future production

The defect affects ONE output sector (tidal k2/dissipation) of ONE model
family (joint ocean campaigns), and the structure/salinity science — the
primary payload of the ocean-chemistry campaigns — is unaffected and
verified. Options for MgSO4/NaCl:

- **(A) Proceed under split-status:** train MgSO4/NaCl now, publish the
  structure sector, quarantine the tidal sector exactly as NH3 does. The
  gates to catch and label the defect all exist.
- **(B) One capacity/embedding pilot first (#4, cheap, existing data):**
  a larger/attention-embedded flow on the NH3 dataset; if it closes the
  Im_k2 gap, fold that architecture into MgSO4/NaCl before their flows
  train.
- **Recommended: A and B in parallel.** Dataset generation is the long
  pole of a campaign and is architecture-independent — start MgSO4/NaCl
  caches + datasets now, run the #4 pilot meanwhile, and pick the flow
  architecture at training time with the pilot's answer in hand.

## Architecture principles (user-ratified 2026-08-06)

The GUI is an engine for scientific exploration, which imposes four
standing constraints on all development:

1. **Traceability.** Every model assumption must be traceable to how the
   calculation is made: config → structure cache → artifact → gate
   reports → GUI scope note form an unbroken, documented chain.
2. **Reproducibility & core parity.** Any individual model inside an
   inference run must be reproducible as a single PlanetProfile run
   (GUI or CLI), and core PlanetProfile evolves WITH the inference work —
   no inference-only physics (NH3 set the precedent: it landed as a
   standard `Ocean.comp` first). The GUI should make switching physical
   features on/off easy. Known documented exceptions: (a) the freegrav
   rho_sil mass-conservation rescale (roadmap: self-consistent Perple_X
   refit); (b) the v4+ direct-Clairaut C20/C22 forward calculation
   (`PlanetProfile/Inference/gravity_obs.py`) is invoked only by the
   inference runner — no standard PlanetProfile-run path computes it
   (doc-doctor 2026-08-06 item 5; roadmap: expose degree-2 gravity
   coefficients as a standard output of a PlanetProfile run so GUI
   single-model runs reproduce the inference observable).
3. **Explorability.** Model assumptions must be easy to understand in the
   GUI (assumptions expanders, scope notes, sector warnings) and results
   easy to explore (the 3D globe panel direction).
4. **Uniform organization across mission-body models.** Same tab/panel
   structure per slot wherever possible. Known accepted deviation: the
   experimental k2 complex-plane / induction plot (not scheduled for
   repair; new slots should follow the standard organization, and any
   intentional deviation gets recorded here). Second accepted deviation
   (recorded 2026-08-07, doc-doctor item 8): the standard tab set uses a
   **Geotherm** tab where the checklist named Geodesy — geodesy figures
   (dual C/MR², u-panel) render in the results flow outside the tab
   strip. Accepted as the shipped organization; a future Geodesy tab
   consolidation is roadmap, not a defect.

## Documentation program

- **Format adjudication: Markdown stays the source of truth.** GitHub
  renders .md natively; raw .html files in a repo display as source, not
  pages (HTML is only useful via GitHub Pages or an external host). If a
  browsable site is wanted later, generate it FROM the .md (MkDocs/Sphinx
  on a Pages branch) rather than authoring HTML — same reasoning as the
  handoff-format adjudication (diffable, agent-readable, no build step).
- **Doc-doctor:** a recurring documentation audit
  (`plans/DOC-DOCTOR.md` checklist) run at every integration point —
  new slot wired, campaign ratified, physics feature added. Codex
  executes the pass; the manager reviews. First pass queued (C8).
- **Methods/results support:** per-campaign provenance (config, priors,
  observables, gate outcomes, versions, seeds) should be liftable into a
  paper's methods section with minimal rewriting. The doc-doctor
  checklist includes a "methods-extractability" item; a generated
  per-campaign methods snippet is roadmap.

## Standing discipline (unchanged)

Preregistered gates; FAILs adjudicated, never relabeled; Machine B stops
at surprises; manager countersigns ratifications; artifact-design changes
need reviewer + user sign-off; nothing deploys unratified.
