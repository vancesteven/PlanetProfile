# PlanetProfile testing architecture — plan (DRAFT for user approval)

Status: **Phase 0 `verified` 2026-11-03** (approved T1–T5 by the user); Phases 1–3 `not implemented`. Drafted 2026-11-02 from a repo survey plus an
independent design by an Opus architect agent; every factual claim below
was re-checked against the repo by the manager (marked ✔).

## Why now

None of the defects found this fortnight was caught by a test. Each was
found by ad-hoc investigation or by a reviewer agent:

| Defect | How it was found | Test class that would have caught it |
|---|---|---|
| ALMA viscosity in Pa·s with the period in kyr (Maxwell/Andrade ice effectively elastic, for months) | Reviewer asked whether k2 varies with η | Unit/dimensional; metamorphic (rescale η and T together → k2 invariant); non-zero sensitivity |
| Rheology hook order: `log10_mu_Ih` forced Maxwell (a705cd8d) | Reviewer read the code | Metamorphic (permute theta keys → identical output) |
| `SFphase` returned a 1-element array for scalar input (87f3ab08) | Running `phaseType='calc'` on Ganymede | Shape contract (scalar in → scalar out) |
| GUI k2 panel plotted NaN placeholders (5853ef25) | User report from the live app | GUI contract ("all-NaN → explicit no-data") |
| `AndradeExponent` broke existing configs (5a2bc744) | Reviewer | Config-compatibility corpus |
| Backend switch changed published-type outputs | Planned, but nothing would have flagged it | Golden regression with a backend tag |
| Reaktoro speciation tests failing at HEAD and before (`mccleskey_speciation_test.py`, `KeyError 'ReaktoroDatabases'`) | Full-suite run today | Any CI at all |

## Current state (✔ = manager-verified 2026-11-02)

- ✔ **Test locations.**
  - `tests/` holds 50 `*_test.py` files, flat; 494 tests collect in ~4 s.
  - `PlanetProfile/Inference/tests/` holds 6 `test_*.py` files, which
    `pytest tests` does not collect.
  - The repo root holds 6 ad-hoc `test_*.py` scripts, 2 of which raise at
    import, so `pytest` from the root aborts collection. Two more are in
    `PlanetProfileApp/`.
- ✔ **Configuration.** One marker (`slow`) and no `conftest.py`.
- ✔ **Skips.** 17 files `skipif` on local data. Most assets are
  force-added to git (63 pkl under `Test/mcmc_results`), so a full clone of
  genai skips little today. A shallow, `main` or deploy-snapshot checkout,
  or any future move of data out of git, would skip silently.
  `enceladus_sbi_guard_test.py` defaults `PP_SCRATCH` to
  `/Users/svance/scratch`, so it always skips off this Mac.
- ✔ **No test CI.** `.github/workflows/` contains only the Sphinx docs build.
- ✔ **`environment.yml` is a macOS export** (`libcxx`); it is unlikely to
  solve on Linux. `requirements.txt` is the only cross-platform spec.
- ✔ **Legacy harness `PlanetProfile/BuildTest.py`** runs `PPTest*` bodies
  end to end and checks only that they run and reload; no outputs are
  compared. Defects:
  - It counts 54 files but imports `PPTest2 … PPTest54`, and 51–54 do not
    exist. The import is outside the `try`, so `full()` should crash with
    `ModuleNotFoundError` after body 50. This is inferred from the code and
    has not been run.
  - The suffixed bodies (`3NoClath`, `41NoClath`, `46_allice`, `49_2km`)
    never run.
  - `Test/` files must not be altered without permission (CLAUDE.md), so
    the new harness wraps the bodies and does not edit them.
    `BuildTest.py` is outside `Test/`; fixing it is the user's call.
- ✔ **Testing tools absent from PPcl:** hypothesis, pytest-xdist,
  pytest-cov, mutmut, pint. Playwright lives in a separate venv.
- ✔ **Full-suite baseline (PPcl, this Mac, 2026-11-02).**
  - `pytest tests` gives **488 passed, 2 failed, 2 xfailed** in 9.3 min,
    excluding the two Enceladus cache-build smokes (216 s and 152 s), which
    passed in the earlier run.
  - The two failures are both in `mccleskey_speciation_test.py`
    (`KeyError 'ReaktoroDatabases'`) and are pre-existing: they also fail at
    5853ef25.
  - No skips on this Mac, because all the data is local. That is exactly
    the case that hides the silent-skip risk.
  - Slowest tests besides the cache smokes: CoolProp NH3 (114 s, 52 s) and
    the GUI AppTests (15–36 s each).
  - A PR tier under 10 min therefore needs the cache smokes, the NH3
    end-to-end tests and the AppTests moved to nightly, or run in
    parallel.

## Test taxonomy

Each layer lists what it catches, where it lives and its runtime budget.

1. **Unit** (`tests/unit/<pkg>/`): logic and error paths. Under 50 ms each;
   the tier runs in under 2 min.
2. **Contract / shape** (`tests/contracts/`): one parametrized test over
   every EOS/phase/property callable.
   - Inputs: a Python scalar, a 0-d array, a 1-d array, and `grid=True`.
   - Assert the output shape and dtype, and that there is no NaN inside the
     declared validity range.
   - Catches the `SFphase` class of defect. Under 1 min.
3. **Metamorphic** (`tests/metamorphic/`): the output must not change
   under transformations that should not matter.
   - Permuting theta keys (the hook-order bug).
   - Rescaling η and the forcing period together (the ALMA unit bug).
   - Splitting a layer into two identical layers.
   - Swapping the backend within the documented bands.
   - Under 5 min.
4. **Non-zero sensitivity** (`tests/sensitivity/`), generated from the
   parameter registry. Each registered parameter declares
   `affects=[observables]` with an expected sign. Perturbing it must move
   those observables by more than ε. A new hook without a declaration fails
   the suite. `tests/alma_viscosity_units_test.py` is the template.
5. **Property-based** (`tests/property/`, Hypothesis): monotonicity, bounds
   and round trips (for example `invert_d_iceIh_to_Tb` ∘ `forward_d_iceIh`
   is the identity). Uses a derandomized `ci` profile with ≤ 50 examples.
6. **Physics validation** (`tests/validation/`, with data in
   `tests/data/published/` carrying provenance headers). Agreement with
   published values: Van Hoolst 2013/2016 libration (already partly in
   `librations_test.py`), the homogeneous-sphere Love number, analytic
   conduction, Pan 2021 conductivity, Hemingway & Mittal 2019 isostasy, and
   the Saur 2024 Fig. 9 induction case (which would pin today's open 1.25×
   gap). A fast subset runs on PRs; the rest nightly.
7. **Golden regression** (`tests/golden/`, references in
   `tests/data/golden/`): per-body scalars as JSON and profiles as `.npz`,
   with a header recording generator, git SHA, environment-lock hash,
   platform, config hash, **backend** and approver. Five fast bodies run on
   PRs; all bodies nightly.
8. **Integration** (`tests/integration/`): config → PlanetProfile → reload;
   cache → forward model → likelihood; slot → artifact. Nightly.
9. **Legacy PPTest wrapper** (`tests/legacy/test_pptest_bodies.py`):
   - discovers bodies with `pkgutil.iter_modules`, not a numeric range, so
     the suffixed bodies are included;
   - runs each in `tmp_path` with `BuildTest.setFullSettings`;
   - checks that a reload matches the run, and compares both with goldens.
   Weekly.
10. **GUI** (`tests/gui/`): AppTest for logic on PRs that touch
    `PlanetProfileApp/`. The Playwright scenarios already in
    `.claude/skills/verifier-streamlit/scenarios/` run nightly on the
    self-hosted runner, with screenshots uploaded. AppTest verifies logic;
    only Playwright evidence counts as rendering verification (CLAUDE.md).
11. **Statistical gates** (`tests/statistical/`): SBC, pocoMC crosscheck
    and slot reproduction, read from the preregistered gate JSON. Weekly,
    or whenever an artifact changes.
12. **Mutation testing** (mutmut, weekly), scoped to the physics kernels in
    the risk list below. Target ≥ 70 % of mutants killed. Every survivor
    gets a new test or an "equivalent mutant" note.

Roughly 70 % of the test count should come from the cheap layers (1, 2 and
5). Layers 3, 4 and 7 are few in number but mandatory: they are the ones
that would have caught this fortnight's defects.

## Layout and conventions

```
tests/
  conftest.py        # markers, --require-data, thread/seed pinning, data_asset fixture
  unit/ contracts/ metamorphic/ sensitivity/ property/ validation/
  golden/ integration/ legacy/ gui/{apptest,browser}/ statistical/
  data/MANIFEST.toml # path, sha256, size, source, tier for every external/large asset
  data/golden/<body>/  data/golden/tolerances.toml  data/published/
  tools/regen_golden.py
```

- **Naming.** Keep the `*_test.py` suffix. Add to pyproject:
  - `testpaths = ["tests"]`;
  - `--strict-markers`;
  - `-ra`, so skip reasons are always printed;
  - `--import-mode=importlib`, which removes the per-file `sys.path` hacks.
  Move `PlanetProfile/Inference/tests/*` into `tests/unit/inference/`, and
  move the root `test_*.py` scripts to `scripts/adhoc/`.
- **Markers:** `slow`, `data`, `gui`, `browser`, `golden`, `legacy`,
  `statistical`, `selfhosted`, and `quarantine(owner, expires, issue)`.
  An expired quarantine fails the build.
- **Data dependencies fail loudly.** A `data_asset("<id>")` fixture resolves
  the asset through `MANIFEST.toml` and `PP_DATA_ROOT`, and verifies its
  sha256.
  - Locally, a missing asset is a skip with a reason.
  - Under `--require-data`, which CI always sets, a missing asset or a hash
    mismatch is an **error**.
  - A terminal-summary table lists every data skip; CI fails if any test in
    a tier it claims to run was skipped for data.
- **Isolation.** A `pp_params` fixture gives a deep-copied `Params` and a
  temporary `HOME`, as `HOME=~/scratch/pphome` does by hand today, so user
  configs never leak into tests.

## Golden references and tolerances

- **Regeneration.** One command:
  `python -m tests.tools.regen_golden --bodies … --reason "…"`. It writes
  the files, updates the manifest, and produces `golden_diff.md` listing
  each quantity's old value, new value, relative change and tolerance.
- **Approval.** A golden PR needs that diff, a CHANGELOG "RESULTS CHANGE"
  entry (the 2026-11-01 entry is the model) and CODEOWNERS approval. Any
  quantity moving beyond tolerance also needs a linked scientific-reviewer
  verdict.
- **Tolerances.** One table, `tests/data/golden/tolerances.toml`, keyed by
  quantity, with `same_platform` and `cross_platform` columns. Proposed
  starting values, to be ratified:
  - integrated quantities (mass, MoI): rtol 1e-10 / 1e-7;
  - root-found boundaries (zb, D_ocean): atol 1 m plus rtol 1e-6;
  - Love numbers: rtol 1e-6, with Im k2 atol 1e-8;
  - libration (a small difference of large terms): rtol 1e-4;
  - phases: exact;
  - NaN positions: exact.
  Never loosen a tolerance in the same PR that changes physics.
- **Stopping goldens from absorbing bugs.**
  - Goldens are always paired with independent validation and metamorphic
    tests, which regeneration cannot touch.
  - A frozen anchor set of published comparisons is refused by the regen
    tool.
  - The backend tag is part of every golden: a backend switch then shows up
    as a mismatch, not as numeric noise.

## Determinism

- `conftest` pins `OMP/MKL/OPENBLAS_NUM_THREADS=1` before numpy is
  imported, and calls `torch.set_num_threads(1)` and
  `torch.use_deterministic_algorithms(True)`.
- Each test gets an RNG seeded from a hash of its node id; the seed is
  printed on failure.
- Hypothesis profiles: `ci` (derandomized), `dev`, `nightly`.
- Goldens are generated on linux-64.
- **SBI.** Never compare weights or raw samples across platforms. Check the
  artifact sha256 against its manifest, and posterior summaries under a
  fixed seed within Monte Carlo error bands (|Δmean| < 3·sd/√n_eff).
  Check training determinism only on the same machine.

## CI (GitHub Actions)

- **Environment.**
  - Replace the 450-line osx export with a minimal `environment.in.yml`
    and `conda-lock` files for `osx-arm64` and `linux-64`.
  - Use `setup-micromamba` with an environment cache keyed on the lock
    hash.
  - Keep a pip-only `requirements.txt` job that mirrors the HF Space.
- **PR tier** (required, under 10 min): import/compile; the fast markers
  with `--require-data -n auto`; AppTest when `PlanetProfileApp/**`
  changes; the manifest/golden hash check; diff coverage ≥ 80 % on changed
  lines (≥ 90 % on kernels). Matrix: ubuntu-latest and macos-14 on Python
  3.11.
- **Nightly:** slow, data, validation, integration and all goldens, with
  data fetched by hash. Python 3.12 added here first.
- **Weekly, on a self-hosted Mac runner** (label `pp-mac`, own-branch
  triggers only, never fork PRs):
  - the legacy wrapper;
  - the statistical gates;
  - the Playwright scenarios;
  - mutation testing;
  - a deploy-snapshot smoke test (`build_deploy_branch.sh --build-only`,
    then AppTest against the staged tree).
- **Not on hosted runners:** the Perple_X tables (2.4 GB), caches over
  40 MB, SBI training, GPU work, and Machine B's multi-hour pocoMC runs.
  Their results enter CI as hashed artifacts that a test verifies.
- **Branch protection on `main`:** the PR tier is required, plus one
  review. CODEOWNERS covers `tests/data/golden/**`, `Gravity/**`,
  `Thermodynamics/**` and `Inference/forward_models.py`.
- **Flaky tests.** Quarantine with an owner and an expiry date. No retries
  on physics tests; retries are allowed only on browser tests.
- **Artifacts on failure:** `golden_diff.md`, test-produced PDFs and PNGs,
  Playwright screenshots, JUnit XML, and the data-skip table.

## Governance

- **Definition of done, mapped to the CLAUDE.md vocabulary.**
  - `verified` requires a named test or artifact that exercises the changed
    path **and fails without the change**. Prove it with a red-first CI job
    that runs the PR's new tests against the parent commit, or with a
    mutation check, as done for a705cd8d today.
  - `implemented, unverified` means the tests are skipped or quarantined,
    or there is no in-browser evidence for a UI change.
- **Division of labour.** The scientific reviewer rules on *what* must be
  true: invariants, published targets, tolerances, and the sign and size of
  sensitivities. Every quantitative ruling becomes a test or a gate file
  that cites it. Tests do not replace adjudication of new physics; the
  reviewer does not replace tests.
- **PR checklist.**
  - Does each new test fail on the old code?
  - Does it exercise the physics path rather than a mock?
  - Is the tolerance justified?
  - Is the data in the manifest?

## Risk-ordered modules (where tests go first)

1. `Gravity/Gravity.py`: backend dispatch, `AlmaViscosity`,
   `AndradeExponent`, `_tidalpy_libration_m`, the padding.
2. `Inference/forward_models.py`: parameter hooks, `apply_parameters`
   ordering, `forward_model_k2_flexible`.
3. `Thermodynamics/HydroEOS.py`: `SFphase` and the EOS callables.
4. `Thermodynamics/LayerPropagators.py` (2,677 lines, essentially untested)
   and the `Geophysical` silicate boundary conditions.
5. `Gravity/Librations.py` and `isostasy.py`.
6. `GetConfig.py`: config versioning and migration.
7. `PlanetProfileApp/` plot builders, against the no-data contract.

## Phased migration

**Phase 0 (2–4 days). Make the existing suite trustworthy.**
- Add pyproject pytest options (testpaths, strict markers, `-ra`,
  importlib mode) and the marker list.
- Move the root `test_*.py` scripts out of the collection path.
- Add `tests/conftest.py`: thread and seed pinning, the `data_asset`
  fixture and `--require-data`.
- Write `MANIFEST.toml` for the ~12 assets the current skipifs use, and
  remove the hard-coded `/Users/svance/scratch` default.
- Triage and fix, or quarantine with an owner, the currently failing
  Reaktoro tests and anything else the full baseline turns up.
- Add a PR workflow on ubuntu using `requirements.txt` (the fast subset).
- Report the BuildTest range/suffix bug (fixing it is a user decision).
- **Success:** a green PR tier, zero silent data skips in CI, and the full
  suite with no unexplained failures.

**Phase 1 (2–3 weeks). The defect classes we actually had.**
- Shape-contract suite (~1 day).
- Hook-order metamorphic test, and an explicit declared hook order (~2
  days).
- Registry-driven sensitivity tests, with an `affects=` field (~3 days).
- η/period scaling invariance, and a TidalPy-vs-ALMA band test (~2 days).
- Config-compatibility corpus (~2 days).
- Cross-platform `conda-lock` (~2 days).
- Diff-coverage gate.
- **Success:** each of the six historical defects, reintroduced as a
  mutation on a branch, is caught by a fast test (6/6); the PR tier stays
  under 10 min.

**Phase 2 (3–4 weeks). Golden references and published validation.**
- Golden tooling: the regen command, `golden_diff.md`, the tolerance
  table and CODEOWNERS.
- Goldens for 5 bodies on PRs and every body nightly, via the legacy
  wrapper.
- Physics-validation suite, including the Saur Fig. 9 benchmark.
- The self-hosted runner.
- Move new caches to hashed `.npz`/HDF5 with JSON headers outside git; the
  existing pickles migrate as they are rebuilt.
- **Success:** every PPTest body has a golden, and nightly stays green for
  two weeks with no expired quarantines.

**Phase 3 (4+ weeks). Depth.**
- Weekly mutation testing on the kernels (≥ 70 % killed).
- Statistical gates in CI, with slot reproduction for every deployed slot.
- Playwright scenarios run nightly.
- The deploy-snapshot smoke test.
- Python 3.12 in the matrix.
- **Success:** over a quarter, defects are first caught by tests rather
  than by review (tracked as the ratio of test-first to review-first
  catches), and a gravity or EOS regression is detected in under a day.

## Decisions for the user

- **T1.** Approve the plan and the phase order.
- **T2.** Permission to fix `BuildTest.py` (outside `Test/`), or leave it
  and let the new wrapper replace it.
- **T3.** Branch protection and CODEOWNERS on `main`. Who approves golden
  updates (the PI)?
- **T4.** A self-hosted CI runner on this Mac, and where large test data
  lives (Zenodo, an HF dataset repo, S3) as caches leave git.
- **T5.** Ratify the starting tolerance table above, or ask the scientific
  reviewer to propose per-quantity values.

## Phase 0 — done (2026-11-03)

- **Collection.** `pytest` from the repo root collects 578 tests: `testpaths`
  covers `tests` and `PlanetProfile/Inference/tests`. Added `-ra` and
  `--strict-markers`, and registered the markers `slow`, `gui`, `data` and
  `statistical`.
- **`tests/conftest.py`.**
  - Pins BLAS/OpenMP threads.
  - Adds `--require-data`. A skip whose reason says "absent"/"not
    present"/"not found" becomes a failure, and a terminal summary table
    lists data skips. Probed in both modes.
  - Adds the `data_asset` fixture, with an optional sha256 check.
  - Auto-marks `gui` (AppTest modules) and `slow` (the measured slow tests).
    Tiers: 535 fast, 36 gui, 7 slow.
- **Found and fixed.** A bare `test_*.py` in `.gitignore` had kept 5
  Inference test files out of git; `pytest tests` never collected the
  directory at all. 16 of its tests were failing silently because their
  hand-written runner stubs had fallen behind `MCMCRunner`.
  `PlanetProfile/Inference/tests/runner_stub.py` now binds every private
  method reachable from the likelihood (discovered from the source) plus
  `__init__`'s literal state defaults, so the stubs cannot go stale. Result:
  73 passed, 8 skipped. The ignore rule is now scoped to `/test_*.py`
  (root-level ad-hoc scripts only).
- **Hard-coded path.** `enceladus_sbi_guard_test.py` no longer defaults to
  `/Users/svance/scratch`: it uses the committed cache and artifact, and the
  real home for scratch.
- **`BuildTest.full()`** runs bodies by name (`ListTestBodies`), including the
  four suffixed variants, and no longer imports the nonexistent
  PPTest51–54. Covered by `tests/buildtest_discovery_test.py` (3 tests).
- **CI.** `.github/workflows/tests.yml`:
  - The fast tier runs on push/PR on ubuntu from `requirements.txt` with
    `--require-data`.
  - The full tier runs nightly or on demand on a self-hosted `pp-mac` runner,
    restricted to `schedule` and `workflow_dispatch` because the repo is
    public.
- **Local results.**
  - Full suite from root: 560 passed, 8 skipped, 0 failed (9.4 min).
  - Fast tier with the CI flags: **525 passed, 0 failed, 2.8 min**.
- **GitHub-hosted fast tier: `verified` 2026-11-03.**
  - Run 41689437340: **523 passed, 0 failed, 4.7 min**.
  - The first run (41687173631) failed only on CI-environment gaps, all now
    fixed: Perple_X tables not downloaded; the SBI guard's thin shard
    scratch-only (now committed under `tests/data/enceladus_sbi/`); and an
    ALMA Im k2 check too strict across platforms. `--require-data` caught
    the shard, as designed.
  - The 10 hosted skips are all non-data reasons: 7 grid fixtures, 1
    deferred NH3, and 2 Reaktoro (conda-only, so absent from
    requirements.txt). The Reaktoro tests are covered by the nightly
    `pp-mac` tier.
- **Previously noted as not yet verified on GitHub** (resolved by the run
  above): The fast tier has not run on a hosted
  runner. Expect missing-dependency fixes on first run (requirements.txt
  omits reaktoro, which is conda-only; those tests use `importorskip`). The
  `pp-mac` runner needs registering by the user. The 8 Inference skips have
  non-data reasons (fixture-specific), so `--require-data` does not fail
  them; give them `data_asset` in Phase 1.
