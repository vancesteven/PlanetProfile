# Development pathway — `genai` → main

**As of 2026-06-30.** This document captures the collaborator workflow feeding `genai` features back to `main`.

## Branches in play

| Branch | Owner | Base | Commits ahead | Last activity | Role |
|---|---|---|---|---|---|
| `main` | upstream | — | — | 5 months ago | Authoritative trunk |
| `genai` | Steven | `main` | 63 | 1 day ago | Active feature-development branch (this worktree) |
| `genai-clean-port` | **Emma Vellard** | `main` | 46 | 7 days ago | **Emma's primary working branch** — clean cherry-picks of `genai` features back to main |
| `integrate-genai-to-main` | Emma Vellard | `main` | 25 | 6 weeks ago | Earlier integration attempt; appears stale (HP-ice / Kalousova / Arrhenius focus) |
| `devin/workstream-B` | Devin (Sonnet IDE) | `genai` | — | today | Worktree at `~/devin-planetprofile-genai`; body-organized MCMC tab refactor |
| `chang-updates` | Chang | — | — | recent | Separate workstream, not yet inventoried |
| (TBD) | **Mara Niesyt** | `main` | — | not yet started | Will branch from `main` for Exploreogram GUI polish + new science |

## Workflow

```
                  (upstream)
                     main
                      │
        ┌─────────────┼─────────────┬─────────────────────┐
        ▼             ▼             ▼                     ▼
      genai      genai-clean-     integrate-      mara/gui-upgrades
   (Steven)        port             genai-          (Mara, future)
                  (Emma)         to-main
                                  (Emma, stale)
```

- **Steven** develops new features (MCMC, SBI scaffolding, Exploreogram fixes, science physics) on `genai`. This branch accumulates broadly and is not expected to merge to main wholesale.
- **Emma Vellard** uses `genai-clean-port` as her primary branch, cherry-picking and re-shaping `genai` features into main-clean commits. Her commits are prefixed `[port]` to mark this. Current focus: 3D Lateral structure, parallelism, obliquity tides, regression tests.
- **Mara Niesyt** (incoming) will branch fresh from `main` for GUI/Exploreogram work. The student migration audit (`STUDENT_MIGRATION_AUDIT.md`) gives her the cherry-pick chain from `genai`.
- **Devin** runs in parallel inside `genai`'s ecosystem on the body-organized MCMC tab refactor. Stays inside `devin/workstream-B`; merges via PR to `genai`, never directly to main.
- All collaborators eventually feed `main` through their respective integration branches, not through `genai` directly.

## Coordination notes

- `genai-clean-port` and the planned Mara branch share the same base (`main`), so cherry-picks of GUI commits from `genai` will likely conflict if both branches independently port `0f56abcb` (the foundational `PlanetProfileApp/` drop). **Check first** whether Emma's `genai-clean-port` already includes `PlanetProfileApp/` before directing Mara to port it.
- Emma's `integrate-genai-to-main` (HP ice + Kalousova content) overlaps with the science-sensitive cluster identified in the student audit. If those features need to land on main, they likely already have a head-start on her branch — no need to re-port from `genai`.
- Steven's `genai` should NOT be merged to main directly. It is the development substrate, not a release branch.

## Suggested next steps

1. When Mara's branch name is known, add it to the table above.
2. Open a short coordination thread with Emma to learn (a) the current state of `PlanetProfileApp/` on `genai-clean-port`, (b) whether `integrate-genai-to-main` is abandoned or still load-bearing, (c) whether the HP-ice / Kalousova / Arrhenius science work there supersedes the corresponding commits on `genai` (which the audit flagged as science-sensitive).
3. Mara should read `STUDENT_MIGRATION_AUDIT.md` before cherry-picking. The recommended chain assumes her branch is bare-from-main; if Emma's port already has PPApp foundation, that step can be skipped.
