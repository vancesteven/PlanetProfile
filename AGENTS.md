# Codex project instructions

Read `CLAUDE.md` completely before planning or changing this repository. Its
scientific constraints, verification discipline, environment guidance, and
`PlanetProfile/Test/` permissions are binding.

Roles (2026-08-01): Machine A runs two agent lanes. Claude Code (Fable 5) is
the model manager — it owns planning, scientific review/adjudication, gate
verdicts, cross-machine coordination, and pushes to origin. Codex 5.6 is the
delegate lane: implementation and reconnaissance tasks explicitly queued in
`plans/CODEX-QUEUE.md`. Machine B (Claude Opus 4.8 / Sonnet 4.6 / Haiku 4.5)
owns heavy compute via `plans/MACHINE-B-HANDOFF.md`.

Codex workflow:

- Read `plans/STATUS.md` first for current state, then `plans/CODEX-QUEUE.md`
  for your tasks and the claim/report protocol. Work only queued tasks unless
  the user directs otherwise.
- Documents classified by `plans/archive/README.md` are provenance, not
  current instructions.
- Commit locally on `genai` with clear messages; do not push — the manager
  reviews and pushes. Record commit hashes in your queue report.
- Report status strictly in the CLAUDE.md vocabulary: `verified` (cite the
  artifact), `implemented, unverified`, or `not implemented`.
- Escalate rather than proceed when a task touches scientific assumptions,
  `PlanetProfile/Test/`, gate thresholds, or anything beyond its written
  scope. Never change a scientific assumption to make a gate pass.
- Do not start intensive campaigns locally; those belong to Machine B.
- Environment on this machine: `mamba activate PP` (PPcl is Machine B's env).

Preserve unrelated working-tree changes.
