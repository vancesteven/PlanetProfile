# Doc-doctor checklist

Recurring documentation audit (STRATEGY.md "Documentation program"). Run at
every integration point: a new GUI slot wired, a campaign ratified, a
physics feature added, or on request. Codex executes a pass as a queued
task; the manager reviews the findings. Each pass produces a short report:
per item PASS / FINDING (with file:line) / N/A — findings become queue
tasks, the report is committed under `validation_reports/doc_doctor/`.

## A. Traceability chain (per GUI slot)

1. The slot's `config_path` exists, and the config's `structure_cache_path`
   matches the slot's `cache_path`.
2. The artifact has a row in `PlanetProfile/Inference/sbi_artifacts/INDEX.md`
   with current state (deployed / delivered / gated / retired) and a gate
   pointer into `validation_reports/`.
3. The scope note states: observables used and dropped (with reasons),
   prior ranges for headline parameters, gate status, and any sector
   quarantine. Numbers in the scope note match the config.
4. Any FAIL-adjudicated gate is described with the adjudication pointer,
   never relabeled.

## B. Reproducibility & core parity

5. Every physical feature the campaign uses is reachable from a standard
   PlanetProfile run (body file / Do flags / Ocean.comp) — no
   inference-only physics. New exceptions require a documented roadmap
   entry in STRATEGY.md.
6. The GUI assumptions expander and the campaign spec agree on what is
   sampled, fixed, and derived.
7. Known reproducibility deviations (currently: freegrav rho_sil rescale;
   between-node bilinear blending; sampled per-phase viscosities) are
   stated in the captions of the affected panels.

## C. Organization uniformity

8. The slot's results panels follow the standard tab organization
   (globe / profiles / wedge / heating / geodesy / mineralogy / data
   where applicable). Intentional deviations are listed in STRATEGY.md
   "Architecture principles" item 4; anything else is a finding.

## D. README and entry points

9. README's documentation links resolve and point at live (not archived)
   docs; plans/README.md routing matches the actual queue files.
10. DEPLOYING.md matches the current deploy script behavior.

## E. Methods extractability

11. For each ratified campaign, a methods paragraph could be assembled
    from committed material alone: config (priors, observables, sigmas),
    cache provenance (grid, comp, builder flags), training (n_train,
    seeds, architecture, versions), gates (criteria + outcomes +
    adjudications). Flag any link in that chain that exists only in chat
    history or uncommitted state.
12. Result-figure provenance: the GUI's exported figures carry enough
    caption/metadata to cite (model version, conditioning values).
