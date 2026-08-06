# Doc-doctor first pass — 2026-08-06

Audited tree: `genai` at `4506ad89` plus the report-only C8 claim
`16be13e1`. Scope: all 13 entries in
`PlanetProfileApp/pages/Inference.py::_SBI_ARTIFACT_SLOTS`, the root and plans
README entry points, `DEPLOYING.md`, and methods extractability for Titan NH3
joint and Europa Clipper v4. No findings were fixed in this pass.

## Findings summary

Checklist totals: **2 PASS, 10 FINDING**. Two placeholder-slot subchecks are
N/A by design.

1. All 11 concrete GUI slots have existing artifacts, configs, and caches;
   every config cache path matches its registry cache path. The two intentional
   MgSO4/NaCl placeholders are N/A.
2. The artifact ledger is stale for the wired Titan freegrav/NH3 and gated
   v5/v6 slots, and the Test50 row lacks a direct gate-report pointer.
3. Most GUI scope notes omit required headline-prior ranges; the three v5
   registry gravity centrals are also rounded rather than config-exact.
4. NH3 and v5/v6 gate descriptions retain superseded or relabeled language.
5. NH3 is core-reachable, but the v4 direct-Clairaut gravity calculation is
   inference-only and is not listed as a STRATEGY exception.
6. GUI assumption text contains three contradictions: nearest-node versus
   blended interpolation, Radau-Darwin versus direct Clairaut, and emcee versus
   pocoMC.
7. Captions disclose rho_sil rescaling and sampled viscosities, but do not
   explicitly disclose the known bilinear-blending deviation.
8. Slots share one organization, but it substitutes a Geotherm tab for the
   checklist's Geodesy tab without a recorded architecture exception.
9. Root README links resolve; `plans/README.md` is untracked and omits the
   Codex queue from its routing.
10. `DEPLOYING.md` still says two manually listed caches/~205 MB; the script
    derives six caches and describes a ~300 MB snapshot.
11. NH3 and v4 methods provenance is extractable from committed material.
12. Downloaded figures have generic filenames/titles and no embedded model or
    conditioning provenance.

## A. Traceability chain

### 1. Config/cache existence and identity — PASS

An AST inventory found 13 registry entries. For all 11 concrete entries the
artifact, config, and cache exist, and the config's `structure_cache_path`
equals the slot's `cache_path` exactly. Evidence by slot:

| GUI slot | Registry | Config cache field | Result |
|---|---:|---:|---|
| Titan Test50 | `PlanetProfileApp/pages/Inference.py:1073` | `PlanetProfile/Inference/configs/test50_titan_noocean_andrade_8D.json:28` | PASS |
| Titan freegrav no-ocean | `PlanetProfileApp/pages/Inference.py:1097` | `PlanetProfile/Inference/configs/titan_freegrav_noocean.json:48` | PASS |
| Titan NH3 joint | `PlanetProfileApp/pages/Inference.py:1152` | `PlanetProfile/Inference/configs/test54_titan_nh3_freegrav.json:49` | PASS |
| Titan MgSO4 placeholder | `PlanetProfileApp/pages/Inference.py:1221` | `config_path=None`, `cache_path=None` | N/A — explicitly awaiting artifact |
| Titan NaCl placeholder | `PlanetProfileApp/pages/Inference.py:1238` | `config_path=None`, `cache_path=None` | N/A — explicitly awaiting artifact |
| Galileo Europa v1.1 | `PlanetProfileApp/pages/Inference.py:1255` | `PlanetProfile/Inference/configs/europa_galileo_v1p1_8D.json:105` | PASS |
| Clipper Europa v4 | `PlanetProfileApp/pages/Inference.py:1284` | `PlanetProfile/Inference/configs/europa_clipper_v4_geodesy_11D.json:191` | PASS |
| Clipper v5 baseline | `PlanetProfileApp/pages/Inference.py:1334` | `PlanetProfile/Inference/configs/europa_clipper_v5_geodesy_11D.json:191` | PASS |
| Clipper v5 no-induction | `PlanetProfileApp/pages/Inference.py:1398` | `PlanetProfile/Inference/configs/europa_clipper_v5_noinduction_7obs.json:135` | PASS |
| Clipper v5 no-k2/h2 | `PlanetProfileApp/pages/Inference.py:1423` | `PlanetProfile/Inference/configs/europa_clipper_v5_nok2_17obs.json:175` | PASS |
| Clipper v6 baseline | `PlanetProfileApp/pages/Inference.py:1478` | `PlanetProfile/Inference/configs/europa_clipper_v6_freegrav_11D.json:187` | PASS |
| Clipper v6 no-induction | `PlanetProfileApp/pages/Inference.py:1540` | `PlanetProfile/Inference/configs/europa_clipper_v6_freegrav_noinduction_6obs.json:131` | PASS |
| Clipper v6 no-k2/h2 | `PlanetProfileApp/pages/Inference.py:1566` | `PlanetProfile/Inference/configs/europa_clipper_v6_freegrav_nok2_16obs.json:171` | PASS |

### 2. Artifact ledger state and gate pointer — FINDING

- PASS: Galileo v1.1 and v4 have deployed rows and direct gate pointers at
  `PlanetProfile/Inference/sbi_artifacts/INDEX.md:41-42`.
- FINDING: the Test50 row has current deployed status but no direct pointer
  into `validation_reports/` (`INDEX.md:38`).
- FINDING: Titan freegrav no-ocean is an active GUI entry
  (`Inference.py:1097-1145`), while its ledger row calls it a candidate “not
  wired for deployment” (`INDEX.md:101`).
- FINDING: the NH3 row says unqualified “RATIFIED (verified)”
  (`INDEX.md:91`), while the GUI and current status quarantine the tidal
  sector (`Inference.py:1158-1170`; `plans/STATUS.md:83`).
- FINDING: all six v5/v6 artifacts are wired but gated with
  `artifact_status='not_ratified'` (`Inference.py:1334-1606`); their ledger
  rows still say “not wired for deployment” (`INDEX.md:94-99`). The delivered,
  not-ratified state and gate paths are otherwise present.
- N/A: MgSO4/NaCl have no artifacts to ledger yet (`Inference.py:1221-1252`).

### 3. Scope-note completeness and numerical agreement — FINDING

- PASS: the NH3 note states used/dropped channels and reasons, the headline
  `Tb` and salinity ranges, split gate status, and tidal quarantine
  (`Inference.py:1195-1219`). Its `[249,263] K` and `1-70 ppt` numbers match
  `test54_titan_nh3_freegrav.json:9-11,43-45`.
- PASS: Titan freegrav no-ocean names its used/dropped channels, headline
  offset-prior boxes, and adjudicated gate status (`Inference.py:1122-1144`).
- FINDING: Test50 (`Inference.py:1088-1092`), Galileo v1.1
  (`Inference.py:1274-1280`), and v4 (`Inference.py:1320-1331`) do not state
  numerical headline-prior ranges. Test50 also does not enumerate dropped
  channels/reasons.
- FINDING: the v5/v6 baseline and ablation notes omit numerical priors for at
  least salinity and the non-hydrostatic offsets; sibling notes merely say
  “same priors” (`Inference.py:1383-1391,1417-1421,1459-1462,1520-1534,
  1558-1564,1601-1606`).
- FINDING: all registry observable centrals match their config exactly except
  the three v5 entries. Their registry C20/C22 values
  `-4.578888e-4/1.377523e-4` (`Inference.py:1356-1374,1410-1415,1435-1451`)
  round config values `-4.57888786233524e-4/1.3775234242885802e-4`
  (`europa_clipper_v5_geodesy_11D.json:101-108`). The differences are tiny
  (~`1.6e-5` and `2.1e-4` observational sigma) but are not config-exact.
- N/A: placeholder scope notes correctly describe fields as pending
  (`Inference.py:1233-1236,1250-1253`).

### 4. FAIL adjudication wording — FINDING

- FINDING: the NH3 ledger row reduces the crosscheck FAIL to three benign eta
  medians (`INDEX.md:91`), but the authoritative ratification explicitly says
  not to read that FAIL as “benign nuisances only” because the marginal gate
  misses the tidal-sector defect
  (`validation_reports/titan_freegrav_nh3_1m/RATIFICATION.md:30-32,110-125`).
- FINDING: the v5 scope note says “undertraining + gate mis-spec, no model
  defect” (`Inference.py:1389-1391`), while the adjudication records a surviving
  headline `D_iceIh_km` shape-excess and underpowered SBC
  (`plans/active/europa-v5v6v7-gate-adjudication.md:11-15`).
- FINDING: the v6 note calls the gate failures benign design consequences
  (`Inference.py:1525-1529`); current adjudication instead records the original
  failures as superseded artifacts, a clean HEAD crosscheck, and remaining
  B1/B2/B6 prerequisites, while requiring FAIL-adjudicated vocabulary rather
  than relabeling (`europa-v5v6v7-gate-adjudication.md:6-15,97-100`).
- PASS: the registry still blocks all v5/v6 conditioning and points to
  `plans/STATUS.md` (`Inference.py:1339-1342,1479-1482`).

## B. Reproducibility and core parity

### 5. Standard PlanetProfile reachability — FINDING

- PASS: NH3 is a standard core ocean composition, not inference-only physics:
  it is registered in `defineStructs.py:227`, implemented in the normal
  hydrosphere EOS path (`PlanetProfile/Thermodynamics/HydroEOS.py:275-350`),
  and the committed cache records `ocean_comp='NH3'` and
  `retry_frozen_as_no_ocean=True`. The campaign spec describes reproducing
  nodes through Titan `Ocean.comp='NH3'` plus ordinary Do overrides
  (`plans/active/titan-nh3-ocean-campaign-spec.md:67-95`).
- PASS-with-documented-exception: freegrav's mass-conservation rho_sil rescale
  is explicitly listed as the current exception and roadmap item
  (`plans/STRATEGY.md:80-87`).
- FINDING: v4's load-bearing direct-Clairaut C20/C22 forward calculation lives
  only under `PlanetProfile/Inference/gravity_obs.py:1-24` and is invoked by
  the inference runner. A repository search for `clairaut`, `dC20_nh`, and
  `gravity_forward_model` outside `PlanetProfile/Inference/` found no standard
  PlanetProfile-run path. STRATEGY lists only rho_sil as an exception
  (`plans/STRATEGY.md:80-87`), so the v4 feature is neither core-reachable nor
  registered as a roadmap exception.

### 6. GUI assumptions versus campaign/config — FINDING

- PASS: the per-run assumptions expander derives sampled, fixed, derived, and
  conditioned fields directly from the loaded config
  (`PlanetProfileApp/Utilities/run_assumptions.py:30-75`). NH3's joint range
  and regime treatment agree between spec and config
  (`titan-nh3-ocean-campaign-spec.md:13-25`;
  `test54_titan_nh3_freegrav.json:53-86`). V4's sampled parameters and direct
  Clairaut design agree between plan and config
  (`europa-clipper-v4-geodesy-plan.md:62-109,129-161`;
  `europa_clipper_v4_geodesy_11D.json:195-223`).
- FINDING: the shared build-up expander says each posterior sample selects the
  nearest Tb grid structure (`Inference.py:2517-2523`), while the per-run text
  and production paths blend structures between nodes
  (`run_assumptions.py:77-86`; v4 config `:218-219`).
- FINDING: the per-run text says gravity-pair configs derive C20/C22 via
  Radau-Darwin (`run_assumptions.py:115-119`), contradicting v4's direct
  Clairaut implementation (`gravity_obs.py:1-10`) and config
  (`europa_clipper_v4_geodesy_11D.json:198`).
- FINDING: the per-run text labels MCMC as the “emcee ensemble sampler”
  (`run_assumptions.py:48-51`), while the shared GUI assumptions identify the
  actual sampler as pocoMC (`Inference.py:2566-2570`).

### 7. Known deviations in affected-panel captions — FINDING

- PASS: sampled/derived rho_sil rescaling is disclosed in wedge/geotherm and
  heating/mineralogy captions (`Inference.py:4023-4032,4128-4150,
  4189-4207,4305-4325`).
- PASS: sampled phase viscosities and Arrhenius handling are stated in the
  radial-profile caption (`Inference.py:3906-3920`) and heating caption
  (`Inference.py:4128-4133`).
- FINDING: the radial caption says only “(Tb,w) structure-cache interpolation”
  (`Inference.py:3906-3910`), and the other affected panels do not identify
  the known between-node **bilinear blending** deviation. The shared expander
  actively says nearest-node selection (`Inference.py:2517-2523`). The
  per-run assumptions mention blending (`run_assumptions.py:79-86`), but that
  is not a caption on each affected results panel as checklist item 7 requires.

## C. Organization uniformity

### 8. Standard result-panel organization — FINDING

- PASS: all executable slots use the same `render_results` implementation and
  the same seven tabs (`Inference.py:3834-3843`); placeholders return before
  results and are N/A. The experimental complex-plane panel is an accepted
  exception (`plans/STRATEGY.md:91-95`).
- FINDING: the implemented tab set is Globe / Radial profiles / Wedge /
  Heating / **Geotherm** / Mineralogy / Data (`Inference.py:3840-3843`), while
  checklist item 8 specifies a **Geodesy** tab where applicable. Geodesy
  figures currently live outside that tab set, and this substitution is not
  listed as an intentional deviation in STRATEGY.

## D. README and entry points

### 9. README links and plans routing — FINDING

- PASS: an automated Markdown-link check found every local target in root
  `README.md` exists, including all five feature guides, the SBI methodology,
  app README, changelog, and `DEPLOYING.md` (`README.md:11-22,172-185`). No root
  link points into an archived plan.
- PASS in the current working directory: all paths named by
  `plans/README.md:3-12` exist and its active/archive/script classification
  matches those directories.
- FINDING: `plans/README.md` routes Machine B work but omits the active Codex
  queue `plans/CODEX-QUEUE.md` (`plans/README.md:3-12`; lane table at
  `plans/STATUS.md:67-71`). It also omits `plans/STRATEGY.md` and the new
  `plans/DOC-DOCTOR.md` entry point.
- FINDING: `git status --short` reports `plans/README.md`, `plans/archive/`, and
  the three `plans/scripts/{active,archive,reproducibility}/` routing trees as
  untracked, so these entry points do not resolve in a clean checkout at the
  audited HEAD.

### 10. DEPLOYING versus deploy script — FINDING

- PASS: the documented HF `upload_folder` path, explicit redeploy model, public
  mode, and optional `DEPLOY_REPO` behavior match
  `plans/scripts/build_deploy_branch.sh:187-237`.
- FINDING: `DEPLOYING.md:3-6` says the snapshot is ~205 MB while the script
  describes ~300 MB (`build_deploy_branch.sh:2-8`; the last recorded C5 build
  was 299 MB in `plans/CODEX-QUEUE.md`).
- FINDING: `DEPLOYING.md:94-101` says Test bulk is excluded except “the two”
  slot caches and tells maintainers to add new files to a copy list. The script
  now derives the cache set from the registry and enforces exact equality
  (`build_deploy_branch.sh:55-96`); the extractor returns **six** current cache
  paths.
- FINDING: the script supports `--no-push` and `--build-only`
  (`build_deploy_branch.sh:10,20-24,187-199`), but DEPLOYING does not document
  either verification/safe-build mode.

## E. Methods extractability

### 11. Ratified NH3 and v4 campaign methods — PASS

**Titan NH3 joint (split ratification):** committed sources provide the full
chain. Priors/observables/sigmas are in
`test54_titan_nh3_freegrav.json:1-49`; cache composition, joint-node mechanism,
and builder flags are in config metadata `:53-86`, cache keys, and
`titan-nh3-ocean-campaign-spec.md:13-25,86-115`; training provenance is in
`validation_reports/titan_freegrav_nh3_1m/RATIFICATION.md:16-22` and artifact
metadata echoed by `crosscheck/crosscheck_report.json:437-446,529-556`
(689,845, nsf, seed 72, torch 2.8.0, sbi 0.26.1); gate criteria/outcomes and
sector quarantine are in `RATIFICATION.md:27-72,110-130`. No required methods
link exists only in chat or uncommitted state.

**Europa Clipper v4:** priors/observables/sigmas and physical conventions are
in `europa_clipper_v4_geodesy_11D.json:1-223`; the reused cache's 93 x 16 grid,
Seawater build, builder flags, and QA are committed at
`plans/europa-clipper-v3-phase1-validation.md:12-20,132-145`; training
provenance is in `validation_reports/europa_clipper_v4_1m/PHASE_GATES.md:1-10`
and `crosscheck/crosscheck_report.json:358-366,456-496` (326,849, nsf, seeds
43/49/4949, 107 epochs, torch 2.8.0, sbi 0.26.1); gate outcomes and the later
governance correction are in `INDEX.md:42` and
`plans/active/europa-v5v6v7-gate-adjudication.md:63-76`. A methods paragraph
can therefore be assembled entirely from committed material, with the current
v4 caveat included.

### 12. Exported result-figure provenance — FINDING

- PASS: the live results page shows the artifact/training count and a table of
  the exact conditioning values adjacent to the figures
  (`Inference.py:2885-2930`).
- FINDING: the common exporter writes SVG/PDF/PNG with bare `savefig` calls and
  no metadata (`PlanetProfileApp/Utilities/crisp_figs.py:55-77`), then uses
  generic filenames such as `k2.pdf` or `corner.png` (`crisp_figs.py:80-91`).
  Figure titles are likewise generic (for example `Inference.py:3207,3282,
  3328,4663,4729,5065`) or identify only a selected sample/Tb, not the model
  version and complete conditioning vector (`Inference.py:3895-3905`). The
  downloaded artifact alone is therefore not citable to a model version and
  conditioning state.

## Verification record

- AST/config/filesystem audit: 13 registry entries; 11 concrete artifacts,
  configs, and caches present; 11/11 config-cache paths exact; two placeholders
  N/A.
- Registry/config observable audit: 8/11 concrete slots exact; only the three
  v5 C20/C22 centrals differ, by rounding as recorded above.
- README local-link audit: all local root README targets exist.
- Deploy extractor: six unique current cache paths.
- Report render verification is recorded in the C8 queue report.
