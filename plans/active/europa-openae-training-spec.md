# Europa open-|Ae| SBI training spec (Machine B)

User directive 2026-07-25: rebuild the Europa SBI training to permit
the full induction-amplitude range (|Ae_synodic| < 1) — no
Galileo-motivated support cut baked into training. Rationale: full
parameter-space coverage in ice thickness / ocean state
(Cochrane2025: high-salinity thick-ocean cases disfavored but not
excludable); Galileo/literature induction knowledge then enters as a
SEPARABLE inference-time reweight, completing the v5
"open-interpretation" pivot already ratified in
europa_clipper_v5_geodesy_11D.json metadata
(open_interpretation_redesign_2026_07_21).

## Why this is cheap on the structure side

NO cache rebuild. The v5 2D (Tb × w) cache + Ae sidecar already span
the full physical support (Machine A audit 2026-07-25, all built
nodes, per label):

| label       | \|Ae\| range     | \|phase\| range | survives current cut |
|-------------|------------------|-----------------|----------------------|
| synodic     | 0.057 – 0.960    | 2.2 – 62.3°     | 38% |
| synodic 2nd | 0.057 – 0.971    | 1.6 – 62.4°     | 47% |
| orbital     | 0.057 – 0.921    | 0.4 – 62.5°     | 31% |

The current `induction_bounds` (synodic/2nd: amp_min 0.7, phase ≤ 30°;
orbital: amp_min 0.3, phase ≤ 70°) discard 53–69% of the physically
built node space — including ALL fresh/near-frozen low-conductance
interiors and (via the validated envelope) the saturated
high-salinity thin-ice corner.

## Config changes (new config: europa_clipper_v7_openae_11D.json)

Copy europa_clipper_v5_geodesy_11D.json, then:

1. `induction_bounds`: REMOVE amp_min for all three labels (or set 0.0).
   Set `phase_deg_max: 70.0` for all three (covers the 62.5° cache
   maximum with margin; the true support edge is then the cache itself,
   which is what we want the flow to learn).
2. Metadata: record this spec, the audit table above, and that the
   Galileo synodic band (|Ae| 0.92–0.97, phase < 1°) plus Levin-2025
   ice thickness remain available as inference-time reweights — the
   science content of the removed cut is NOT lost, just moved to where
   it can be toggled.
3. Everything else (priors incl. D_iceIh U[5,80], log10_w U[-1,2],
   observables, sigmas, cache path, arrhenius) unchanged.

## Training + gates (Machine B)

- Dataset + flow: same recipe/scale as v5/v6 1m builds. Expect the
  accepted-simulation fraction to RISE (fewer support rejects), so
  n_train_effective goes up at equal simulation budget.
- Reference MCMC: rerun with the SAME opened induction_bounds —
  crosschecking the open flow against the old (cut) reference would
  fail by construction.
- Grid-walk anchors: extend to span |Ae_synodic| ≈ 0.06–0.96 (both new
  corners: near-frozen low-amp AND saturated high-amp). The GUI slot's
  validated envelope (`derived_ae_guards.ae_range`) then widens
  accordingly; `warn_support_below: 0.7` is DELETED for this artifact.
- Standard gates: SBC, crosscheck, limits — same pass criteria as
  v5/v6 (multiplicity-aware, SHAPE_EXCESS residual clause).
- Known risks to watch in gates: (a) posterior mass piling at the
  |Ae| saturation edge for tight synthetic conditioning (flow edge
  behavior — this is the "hard support edge" concern that originally
  motivated the cut; now the edge is the physical saturation, and the
  gate must probe conditioning near it); (b) the near-frozen corner is
  ~0 ocean — verify the no-ocean guard behavior (D_ocean → 0 nodes are
  built and legitimate here, unlike configs that hard-reject them).

## GUI follow-up (Machine A, after artifact ships)

- New slot in _SBI_ARTIFACT_SLOTS: envelope from the new anchors, no
  warn_support_below, scope note stating open-support + reweight
  philosophy; optional one-click "apply Galileo induction reweight"
  toggle (ties into the reweight roadmap in the v5 metadata).

## Out of scope here (separate handoff items already filed)

- NaCl-composition cache variants (SeaFreeze 1.1.3 NaClaq splines).
- Callisto campaign support design (do not inherit Europa cuts;
  Cochrane2025).
