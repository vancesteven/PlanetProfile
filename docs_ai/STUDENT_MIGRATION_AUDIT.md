# Student migration audit — `genai` → student's main-based branch

**Generated:** 2026-06-30
**Scope:** 63 commits on `genai` ahead of `main`, classified for cherry-pickability into the student's branch.
**Method:** three sonnet sub-agents each audited a 21-commit slice; opus synthesized.
**Audience:** Steven (gatekeeper), then the student taking over Exploreogram polish + new science.

---

## TL;DR

- `main` currently has **neither** `PlanetProfileApp/` nor `PlanetProfile/Inference/`. The PPApp foundational drop (`0f56abcb`) must land before any Exploreogram patch.
- Of 63 commits: **~9 are cleanly portable** (docs + Exploreogram chain), **~12 are conditionally portable** (need partial cherry-pick or main has changed), **~36 are entangled with MCMC scaffolding** (do not port), **~6 are science-sensitive** (need scientific-reviewer if ever ported).
- **18 commits touch `PlanetProfile/Test/`** — forbidden territory per CLAUDE.md regardless of category.
- The student's recommended chain for Exploreogram-only work is short: **`0f56abcb` (foundation) → `80aa882a` (exploreogram enhancement) → `4af3f31d` → `30e5abc3` → `ab1cffaa` (partial) → `8a631cc4` → `b7da7807` (partial) → `86faa709` (partial)**.

---

## How to read this doc

- **portable-clean** — apply with `git cherry-pick` directly, low conflict risk.
- **portable-with-deps** — needs a prior commit applied first, or needs partial application (`git cherry-pick -n` then drop forbidden hunks before recommit).
- **entangled** — only meaningful inside the genai MCMC scaffolding; **do not port**.
- **science-sensitive** — touches forward physics / thermodynamics / inference math; **requires scientific-reviewer before porting**.

A commit can be both portable-with-deps and science-sensitive; in that case the science gate dominates.

---

## Recommended student port chain (Exploreogram-only)

These are in chronological order, which is also the cherry-pick order:

| # | SHA | What it adds | Strategy |
|---|---|---|---|
| 1 | `0f56abcb` | **Foundation:** entire `PlanetProfileApp/` tree including initial `Exploreogram.py` | Bedrock. Must land first if main does not already have PPApp. |
| 2 | `80aa882a` | Exploreogram derived params, Illustrator fonts, export scripts, Plotly layout, `explore_cache.py` | Watch `Utilities/defineStructs.py` hunk (+99 lines) for conflicts. |
| 3 | `04cdf820` | `.gitignore` hygiene for `CLAUDE.md` / `.claude/` | Trivial. |
| 4 | `c1112efa` | README link to PPApp | Trivial. |
| 5 | `4af3f31d` | Salinity↔conductivity secondary axis | Standalone, +155 lines in `exploreogram_plotly.py`. |
| 6 | `30e5abc3` | Per-axis log-scale toggles | Tiny, +20 lines. |
| 7 | `ab1cffaa` (**partial**) | Secondary axis fixes + pcolormesh + log-scale UI | Drop `Inference.py` and `mcmc_run_loader.py` hunks. |
| 8 | `8a631cc4` | Wire induction component selector | Clean. |
| 9 | `b7da7807` (**partial**) | Matplotlib component routing (Bx/By/Bz → zName) | Drop the `CLAUDE.md` hunk (genai-specific). |
| 10 | `86faa709` (**partial**) | `sigmaMean_Sm` added to `DERIVED_PARAMS` | Take ONLY the 2-line `Exploreogram.py` hunk; drop the 5 config JSON edits. |

### How to do a "partial" cherry-pick

```bash
git cherry-pick -n <sha>          # stage all hunks without committing
git restore --staged <forbidden-file>     # unstage forbidden file
git checkout HEAD -- <forbidden-file>     # revert the working-tree change
git commit -m "fix(gui): <reword>"
```

### Optional ports the student could consider separately

| SHA | What it adds | Why optional |
|---|---|---|
| `d4ba5560` | 1-line `int()` → `bool()` fix in `HydroEOS.py` | Standalone numpy-compat fix, valuable for main regardless. |
| `4427c4c5` | Three user-facing physics guides (`CONVECTION_GUIDE.md` etc.) | Docs-only, no dependencies. |
| `283dbab0` | `.gitignore` for inference cache pkls | Only useful if any inference port follows. |

---

## Per-commit classification (all 63 commits)

### Slice A — most recent 21 (HEAD..HEAD~21)

| SHA | Subject | Category | Classification | Student? |
|---|---|---|---|---|
| `9b022fa9` | docs consolidate to SESSION.md | docs | portable-clean | no |
| `b7da7807` | fix(gui): matplotlib induction component routing | GUI-exploreogram | portable-with-deps (partial) | **yes** |
| `8a631cc4` | fix(gui): induction component selector + axis routing | GUI-exploreogram | portable-with-deps | **yes** |
| `ab1cffaa` | feat(gui): secondary axis + pcolormesh + log-scale | GUI-exploreogram (mixed) | portable-with-deps (partial) | **yes** |
| `30e5abc3` | feat(gui): per-axis log-scale toggles | GUI-exploreogram | portable-clean | **yes** |
| `46c0723d` | feat(gui): mcmc_run_loader utility | GUI-inference | entangled | no |
| `4af3f31d` | feat(gui): linked salinity↔conductivity axis | GUI-exploreogram | portable-clean | **yes** |
| `86faa709` | feat(inference+gui): planet_template_module + sigmaMean_Sm | GUI-exploreogram (mixed) | portable-with-deps (partial) | **yes** (partial only) |
| `7b49f1a5` | fix(inference): hoist arrhenius_params | MCMC-scaffolding | entangled | no |
| `a9f24584` | chore(numpy2): migrate deprecated dtype aliases | physics-science | **science-sensitive** | maybe |
| `bb6dd87d` | perf(inference): pre-compute induction Ae | MCMC-scaffolding | entangled | no |
| `4d5509e6` | feat(inference): T25 + CMR2 + diagnostics + NumPy2 | MCMC-scaffolding | entangled (touches Test/) | no |
| `e2c53e5b` | feat(inference): Titan MCMC results + wedge plots | MCMC-scaffolding | entangled (touches Test/) | no |
| `8ed5cfce` | fix(inference): wedge plot rendering | MCMC-scaffolding | entangled | no |
| `5c050a9d` | feat(inference): plot naming + Titan 4km/2km + Europa Test51 | MCMC-scaffolding | entangled (adds Test/) | no |
| `f0a4dc65` | feat(inference): Phase C1 toolkit | MCMC-scaffolding | entangled (touches Test/) | no |
| `8c79fa8d` | fix(gui): zero-variance columns + era scoring | GUI-inference | entangled | no |
| `d9f230f5` | feat(inference): Titan tidal heat + SBI pre-training + Era | physics-science | **science-sensitive** | no |
| `909c6a73` | docs: split CLAUDE.md | docs | portable-clean | no |
| `be264f5d` | feat(inference): D_ocean_km cache key + Tb priors | MCMC-scaffolding | entangled | no |
| `c268efd9` | feat(inference): k₂/h₂ observables + three-body validation | MCMC-scaffolding | **science-sensitive** | no |

### Slice B — middle 21 (commits 22–42)

Headline: **zero Exploreogram content**. This is the MCMC v2.1 / Kalousova / Yao2014 / Test46–50 construction era.

| SHA | Subject | Category | Classification | Student? |
|---|---|---|---|---|
| `ca1b6006` | fix(inference): six MCMC bugs | MCMC-scaffolding | entangled | no |
| `cd351db9` | feat(inference): Phase C1 v2.1 caches | MCMC-scaffolding | **science-sensitive** | no |
| `578a94af` | refactor(test50): delegate to MCMCRunner | MCMC-scaffolding | entangled (touches Test/) | no |
| `f2c4154a` | feat(inference): Phase B toolkit | MCMC-scaffolding | entangled | no |
| `c5ae78bf` | docs(mcmc): MCMC_ARCHITECTURE.md | docs | entangled (touches Test/) | no |
| `247f89ff` | fix(kalousova): Ice VI placeholder | physics-science | **science-sensitive** | no |
| `80853f4f` | fix(silicate): T_center-anchored conduction | physics-science | **science-sensitive** | no |
| `fe610a30` | fix(clathrate depth) + Test50 8D + rename | physics-science | **science-sensitive** (touches Test/) | no |
| `0ff1f60a` | Document MCMC toolkit + CLAUDE.md | docs | entangled | no |
| `03b21d01` | Add Test50 Titan no-ocean MCMC | MCMC-scaffolding | entangled (Test/ only) | no |
| `23b110fa` | Extract mcmc_common + mcmc_plots, port Test48 | MCMC-scaffolding | entangled (touches Test/) | no |
| `4427c4c5` | Add CONVECTION/ARRHENIUS/TIDAL_LOVE guides | docs | portable-clean | maybe |
| `75667fae` | Update MCMC guide + .claudeignore | docs | entangled | no |
| `4cb41cd5` | Yao 2014 spherical convection + MCMC harness | physics-science | **science-sensitive** (touches Test/) | no |
| `80aa0d36` | Test46 Andrade hybrid + 3727-file cleanup | MCMC-scaffolding | entangled (huge churn, touches Test/) | no |
| `a1f93576` | Test46 + rename Europa→PPTest47 | MCMC-scaffolding | entangled (Test/ only) | no |
| `9a1212b5` | PPTest45 6D MCMC harness | MCMC-scaffolding | entangled (touches Test/) | no |
| `d4ba5560` | Fix TypeError in GetTfreeze (numpy compat) | physics-science | portable-clean | maybe |
| `120996cf` | Add PPTest45 hybrid hydrosphere path | MCMC-scaffolding | entangled (touches Test/) | no |
| `d1ebe50f` | Fix TidalPy thin-layer + Test42 params | GUI-inference | entangled (touches Test/) | no |
| `ceb9f4f7` | Fix param-order in MCMCRunner | MCMC-scaffolding | entangled | no |

### Slice C — oldest 21 (commits 43–63)

| SHA | Subject | Category | Classification | Student? |
|---|---|---|---|---|
| `0bcdf709` | Maxwell+Tb_K cache 70% acceptance | MCMC-scaffolding | entangled | no |
| `74523a68` | Phase 3 viz + tab rename + Europa tidal | GUI-inference | entangled (touches Test/) | no |
| `0f962ae3` | MCMC working: pocoMC API fix | MCMC-scaffolding | **science-sensitive** (touches Test/) | no |
| `22d6d028` | Per-phase zeta fix (3775 files, +182k lines) | MCMC-scaffolding + infra | entangled (poison commit) | no |
| `5454e6ee` | Auto-populate structure cache path | GUI-inference | entangled | no |
| `dfaf73fa` | Fix NameError in render_parameter_config | GUI-inference | entangled | no |
| `33f6432b` | Fix ImportError create_help_button | GUI-inference | entangled | no |
| `79a46a9d` | Phase 2: parameter config UI | GUI-inference | entangled | no |
| `283dbab0` | .gitignore for inference cache | infra | portable-clean | maybe |
| `196dc8f1` | Phase 2 prep: parameter registry | MCMC-scaffolding | entangled | no |
| `5c6487c9` | Phase 5: mcmc_runner bug fixes | MCMC-scaffolding | entangled (touches Test/) | no |
| `b821d401` | Phase 5: MCMCRunner with dual cache | MCMC-scaffolding | entangled (touches Test/) | no |
| `c2419132` | Titan rheology Tests 41–44 | MCMC-scaffolding | entangled (touches Test/) | no |
| `c1112efa` | README link to PPApp | docs | portable-clean | **yes** |
| `4b65e5c7` | Update docs: TidalPy + MoonMag | docs | portable-with-deps | maybe |
| `fc9a5c93` | TidalPy backend + per-phase tidal heating | physics-science | **science-sensitive** (touches Test/) | no |
| `272499b0` | Arrhenius HP + Kalousova + MoonMag + Lateral (135 files) | physics-science | **science-sensitive** (touches Test/) | no |
| `80aa882a` | Exploreogram enhancements + explore_cache.py | GUI-exploreogram | portable-with-deps | **yes** |
| `04cdf820` | .gitignore CLAUDE.md + .claude/ | infra | portable-clean | **yes** |
| `0f56abcb` | **Foundation: PlanetProfileApp/ entire tree** | GUI (foundational) | portable-with-deps | **yes** |
| `73ac449f` | CSV export + inductogram test configs | physics-science | entangled (touches Test/) | maybe |

---

## Key dependency notes

- The Inference module is monolithic: `196dc8f1` (registry) → `b821d401` (MCMCRunner) → cascades all the way up to `ca1b6006`. You cannot peel off late MCMC commits without the foundation.
- `22d6d028` is a **poison commit**: 3775 files, vendors TidalPy and `.claude-flow/`, mixes the per-phase zeta physics fix with massive repo noise. Treat its physics content as lost unless extracted by hand.
- `272499b0` is **also massive** (135 files, MoonMag vendoring, Kalousova HP-ice viscosity, Lateral module). Required prerequisite for anything that imports `PlanetProfile.MagneticInduction.MoonMag` or `PlanetProfile.Lateral`. If main has neither, this is a non-starter for any later port.
- `0f56abcb` is the **foundational PPApp drop**. Verified absent from main as of 2026-06-30. Everything in PlanetProfileApp/ assumes it.

---

## Science-sensitive cluster (require scientific-reviewer before any port)

Six commits touch forward physics or numerical baselines:

1. `a9f24584` — NumPy 2.x dtype migration across 25 core modules (induction, gravity, thermodynamics). Prior Gemini regression on `np.bool_` silently broke MATLAB-export `isinstance` checks. Per [[feedback-numpy2-baselines]] if any such memory exists, otherwise: vet numerical baselines.
2. `cd351db9` — v2.1 transition-aware caches; touches `LayerPropagators.py` HP-ice blending.
3. `247f89ff` — Kalousova Ice VI melting-curve placeholder + meltFraction=0.5.
4. `80853f4f` — Silicate `T_center`-anchored conduction (solid-sphere bodies, `c1=0`). Triggered by Pluto Test 15.
5. `fe610a30` — Clathrate-depth conduction (Turcotte-Schubert eq 4.40) + bulk Test/ rename.
6. `4cb41cd5` — Yao 2014 spherical convection + per-phase Arrhenius dispatch.
7. `fc9a5c93` — TidalPy backend + per-phase tidal heating redistribution.
8. `272499b0` — Arrhenius HP-ice + Kalousova convection + MoonMag vendoring + Lateral module.
9. `0f962ae3` — Per-phase zeta (the actual physics inside the poison commit `22d6d028`).
10. `c268efd9` — Body-specific k₂/h₂ priors (Mazarico+2023 placeholders).
11. `d9f230f5` — Titan tidal heat partitioning + aggregate-reservoir definitions.

Per [[project-kalousova-hp-ice]] and [[project-yao2014-spherical-convection]], these are the load-bearing physics commits. None are for the student's Exploreogram workstream, but some may eventually be wanted on main.

---

## What to do next

1. **Steven reviews** this audit and approves the recommended chain.
2. The student is given a **shortlist of 10 commits** (the chain above), plus a 2-paragraph explanation of "partial" cherry-pick mechanics.
3. **Before student starts:** confirm whether their branch already has PPApp foundation. If yes, skip `0f56abcb`. If no, that commit becomes step 1.
4. **scientific-reviewer pass** triggered separately for the 11 science-sensitive commits before any upstream merge to main.
5. Devin's workstream B (body-organized MCMC tab refactor) runs in parallel in the worktree at `~/devin-planetprofile-genai`; it will NOT collide with the student's Exploreogram work since Devin only touches Inference module + Inference.py page.

---

## Outstanding questions for Steven

1. Does the student's branch already have any PPApp commits, or is `0f56abcb` truly the bedrock for them?
2. Are the three docs commits (`4427c4c5`, `c1112efa`, `4b65e5c7`) worth offering to main upstream separately, or are they reserved for the genai narrative?
3. Should `d4ba5560` (1-line `bool()` fix) and `a9f24584` (NumPy 2.x migration) be split out as their own upstream PRs to main, independent of the student? They're the highest-value standalone non-MCMC fixes in the 63-commit set.
