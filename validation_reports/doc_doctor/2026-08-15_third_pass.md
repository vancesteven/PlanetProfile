# Doc-doctor third pass — 2026-08-15

Full pass of `plans/DOC-DOCTOR.md` against committed `genai` at
`e74144f2` (C24 claim), after C20–C23. Outcome: **11 PASS, 1 FINDING**.
The Cassini–Enceladus registry entry is an explicit awaiting-artifact
placeholder, so artifact/config/campaign checks are N/A for that slot until
ratification. Concurrent unstaged work in
`PlanetProfile/Inference/cache_builder.py`, `mcmc_runner.py`, and
`parameter_registry.py` was not staged, modified, or used as evidence for this
report-only pass.

## A. Traceability chain

### 1. Config/cache existence and agreement — PASS

All 13 executable registry slots have a committed artifact, config, and cache;
each config's `structure_cache_path` exactly equals its registry `cache_path`.
The registry spans the Titan entries at
`PlanetProfileApp/pages/Inference.py:1072-1340` and the Europa entries at
`PlanetProfileApp/pages/Inference.py:1341-1724`. The standing regression loads
each config, resolves all three paths, and compares the cache paths at
`tests/slot_reproduction_test.py:49-86`; the focused run passed 14 tests (13
per-slot config/artifact checks plus the cache-manifest check) with the slow
sampling reproduction deselected.

N/A for Cassini–Enceladus: the entry intentionally sets artifact, config, and
cache to `None` with `artifact_status='awaiting_artifact'`
(`PlanetProfileApp/pages/Inference.py:1725-1744`), guarded by AppTest assertions
at `tests/enceladus_placeholder_test.py:76-85`.

### 2. Artifact ledger state and gate pointer — PASS

Every executable slot has a current ledger row and a validation pointer. The
five Titan rows are at
`PlanetProfile/Inference/sbi_artifacts/INDEX.md:47-50,116`; the current Galileo
v1.1 and Clipper v4 rows are at `INDEX.md:53-54`; and the v5/v6 baseline and
ablation rows are at `INDEX.md:109-114`. A mechanical artifact-name scan found
all 13 current rows and a `validation_reports/` pointer in each.

N/A for Cassini–Enceladus because no artifact exists yet
(`PlanetProfileApp/pages/Inference.py:1725-1744`).

### 3. Scope-note completeness and numerical agreement — PASS

The Titan notes state their used/dropped channels and reasons, headline prior
ranges, gate posture, and applicable quarantine/caveats
(`PlanetProfileApp/pages/Inference.py:1088-1153,1167-1225,1234-1282,1287-1339`).
The Europa notes do the same for Galileo v1.1, v4, and the ratified v5/v6
families (`Inference.py:1360-1370,1409-1423,1432-1491,1499-1570,1587-1718`).
A mechanical comparison of every configured observable central value against
the registry `default_obs` found exact agreement for all 13 executable slots,
including all channel-ablation arms.

N/A for the awaiting Enceladus slot: its note is explicitly frozen-design
intent rather than a claim about an artifact (`Inference.py:1735-1744`).

### 4. FAIL-adjudication wording — PASS

The ledger defines the binding vocabulary and requires surviving FAILs to
remain FAILs with an adjudication and visible scope control
(`PlanetProfile/Inference/sbi_artifacts/INDEX.md:19-29`). The NH3 record retains
the crosscheck FAIL and tidal-sector quarantine (`INDEX.md:48`); MgSO4 and NaCl
retain the gated limits FAILs and the NaCl nuisance FAIL (`INDEX.md:49-50`);
v4 retains the later hard-breach record and user-decision pointer
(`INDEX.md:54`); v5 retains `FAIL-ADJUDICATED-ACCEPTABLE` (`INDEX.md:109`); and
Titan freegrav retains both adjudicated FAIL exits (`INDEX.md:116`). The visible
slot text carries the corresponding split/quarantine/adjudication wording at
`PlanetProfileApp/pages/Inference.py:1147-1153,1167-1175,1234-1245,1287-1301,
1432-1436`.

## B. Reproducibility and core parity

### 5. Standard PlanetProfile reachability — PASS

The active campaign physics is reachable through standard PlanetProfile code:
NH3 is a normal `Ocean.comp` option and hydrosphere EOS branch
(`PlanetProfile/Utilities/defineStructs.py:228`;
`PlanetProfile/Thermodynamics/HydroEOS.py:275-347`), and direct hydrostatic
C20/C22 is exposed by the default-off `Do.CALC_C20_C22` flag
(`PlanetProfile/Utilities/defineStructs.py:139`;
`PlanetProfile/Gravity/Gravity.py:41-83`). STRATEGY records that direct-Clairaut
exception as closed and retains only the documented freegrav density-rescale
exception (`plans/STRATEGY.md:126-138`).

N/A for Enceladus isostasy as a campaign feature: the GUI entry is
non-executable and awaiting its artifact/config (`Inference.py:1725-1744`), so
the frozen isostasy design does not yet make a shipped inference claim.

### 6. GUI assumptions versus campaign configuration — PASS

The per-run assumptions are generated from the loaded config's sampled, fixed,
derived, and conditioned fields
(`PlanetProfileApp/Utilities/run_assumptions.py:30-75`), and now describe the
actual cache interpolation, per-phase rheology,
pocoMC/SBI mode, and direct-Clairaut declaration (`run_assumptions.py:77-129`).
The shared expander agrees: transition-aware 1D/bilinear 2D structure handling,
per-phase viscosities, direct density-profile C/MR2 treatment, and pocoMC are
stated at `PlanetProfileApp/pages/Inference.py:2647-2715`. No sampled/fixed/
derived contradiction was found in the current campaign specs.

### 7. Known deviations in affected-panel captions — PASS

The affected result captions now disclose all three standing deviations:
sampled viscosities and interpolation in radial profiles
(`PlanetProfileApp/pages/Inference.py:4164-4172`), freegrav/mass-conservation
silicate-density handling plus bilinear structure in wedge/heating/geotherm
(`Inference.py:4244-4263,4284-4293,4389-4412,4451-4475`), and bilinear
between-node structure in mineralogy (`Inference.py:4567-4601`).

## C. Organization uniformity

### 8. Standard result tabs — PASS

All executable slots flow through the shared seven-tab renderer: Globe, Radial
profiles, Wedge, Heating, Geotherm, Mineralogy, and Data table
(`PlanetProfileApp/pages/Inference.py:4068-4077`). STRATEGY explicitly records
the complex-plane plot and the accepted Geotherm-for-Geodesy organization
(`plans/STRATEGY.md:142-151`). The focused globe/results suite passed 11 tests,
including Europa v4, heating, mineralogy, no-ocean gating, and the static-globe
fallback.

N/A for the awaiting Enceladus entry because it stops before results rendering
(`tests/enceladus_placeholder_test.py:43-74`).

## D. README and entry points

### 9. README links and routing — PASS

`plans/README.md` routes first to STATUS and then to Machine B/Codex, strategy,
doc-doctor, and active/archive/script indexes (`plans/README.md:1-17`). The
active index points to the live queue and current Enceladus/spec/remedy records
(`plans/active/README.md:8-28`). A local Markdown-link resolution check found no
missing relative targets in `README.md`, `plans/README.md`, or
`plans/active/README.md`; the root deployment entry resolves at
`README.md:182-183`.

### 10. DEPLOYING versus deploy script — PASS

The docs match the script's `--build-only`/`--no-push` behavior
(`DEPLOYING.md:38-50`; `plans/scripts/build_deploy_branch.sh:10-24,187-210`),
registry-derived cache extraction and exact-set guard
(`DEPLOYING.md:103-125`; `build_deploy_branch.sh:55-96`), and Hugging Face
upload path (`DEPLOYING.md:22-36`; `build_deploy_branch.sh:212-238`). A real
`--build-only` verification staged exactly 8 registry caches, produced a 317 MB
/ 537-file snapshot, and exited before any local or remote branch push.

## E. Methods extractability

### 11. Ratified-campaign methods chain — FINDING

PASS for the eight generated campaign families: NH3, MgSO4, NaCl, v5, v6,
Galileo v1.1, v4, and Test50 are registered with committed config, artifact,
cache, gate, and deployment sources at
`plans/scripts/generate_methods_snippet.py:33-133`. Running `all --verify`
confirmed all eight outputs' source values, hashes, parameter/observable names,
and Markdown rendering.

FINDING: the separately deployed Titan free-gravity no-ocean campaign is listed
as reviewer-verified and deployed at
`PlanetProfile/Inference/sbi_artifacts/INDEX.md:116`, but it has no campaign key
in the generator's complete table (`generate_methods_snippet.py:33-133`) and no
`validation_reports/titan_freegrav_noocean_1m/methods_snippet.md`. Its methods
chain remains manually recoverable from committed config/artifact/gate sources,
but item 11's per-ratified-campaign extractability is not yet satisfied by the
standing generated artifact. Add this campaign in a separately scoped follow-up;
no fix was made during this report-only audit.

N/A for Cassini–Enceladus until an artifact is trained and ratified.

### 12. Exported result-figure provenance — PASS

Every inference crisp export is wrapped with the active model label, artifact,
slot ID, and exact conditioning summary
(`PlanetProfileApp/pages/Inference.py:2986-3025`). The common exporter embeds
that provenance plus app git SHA into SVG/PDF/PNG metadata and uses a
model/date-aware filename (`PlanetProfileApp/Utilities/crisp_figs.py:58-88,
111-174`); its public API documents the provenance contract at
`crisp_figs.py:177-194`. The same conditioned values are also shown adjacent to
the results at `Inference.py:3110-3127`.

## Verification record

- `pytest -q tests/slot_reproduction_test.py -m 'not slow'` — 14 passed,
  1 deselected.
- `python plans/scripts/generate_methods_snippet.py all --verify` — all eight
  configured campaign snippets verified and rendered.
- `pytest -q tests/app_globe_panel_test.py` — 11 passed.
- `bash plans/scripts/build_deploy_branch.sh --build-only` — exact 8-cache,
  317 MB, 537-file snapshot; no push.
- `bash -n plans/scripts/build_deploy_branch.sh` — pass.
- Local relative-link check for the root/plans/active READMEs — pass.
