# Plans and handoffs

Start with `plans/STATUS.md`. It is the only maintained current-state summary.
Then read `plans/MACHINE-B-HANDOFF.md` if the task concerns heavy computation.

- `plans/CODEX-QUEUE.md` is the Machine A Codex implementation queue and
  defines its claim/report protocol.
- `plans/STRATEGY.md` records active project strategy, accepted deviations,
  and explicit exceptions.
- `plans/DOC-DOCTOR.md` is the repeatable documentation-audit checklist.
- `plans/active/README.md` lists current scientific documents.
- `plans/archive/README.md` classifies historical handoffs and superseded
  campaigns.
- `plans/scripts/active/README.md` identifies operational scripts.
- `plans/scripts/reproducibility/README.md` identifies scientific evidence and
  figure-regeneration scripts.
- `plans/scripts/archive/README.md` identifies retired campaign scripts.

Historical files intentionally remain at stable paths. Config JSON, tests, code
comments, reports, and artifact provenance embed those paths; cosmetic moves
would either break links or change trained-artifact config hashes. Archive
classification is therefore index-based rather than a physical relocation.
