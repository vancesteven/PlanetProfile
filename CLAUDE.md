planetprofile genai, worktree for rapidly developing new features that other branches use to migrate back to main

for python env use "mamba activate PPcl"

current focus is debugging PlanetProfile/PlanetProfileApp exploreogram and mcmc features.

Do not alter Test/ files without permission

Do not change fundamental scientific assumptions of the code
Have the scientific-reviewer check work frequently while planning and implementing.

Use sonnet and haiku agents whenever possible. Use Opus agents for planning and review.

## Verification discipline for UI and behavioral changes

A change is NOT "done" until its specified behavior has been observed in the running system.

For PlanetProfileApp/ (Streamlit GUI) changes:
- `python -m py_compile` passing does NOT count as verification. Neither does code review or "syntax clean."
- A UI bug fix is unverified until the app is launched, the specific user interaction is reproduced, and the targeted behavior is visually confirmed (screenshot, PDF inspection, or the `/verify` skill).
- For matplotlib output: the produced PDF must be opened and the specific change inspected.
- For Plotly/Streamlit output: the chart must be rendered in-browser and the change confirmed.
- If the agent cannot launch the app or cannot visually verify (no display, no time, etc.), the agent must say so explicitly and label the change `implemented, unverified` — never `done` or `fixed`.

For non-UI behavioral changes (forward-model, inference, plotting library code):
- A targeted test or benchmark must exercise the changed code path and produce the expected numerical/structural result.
- Compiling, importing, or running a smoke pass that does not exercise the specific change is not verification.

In handoff documents, status entries must use this vocabulary:
- `verified` — observed working under the documented reproduction steps; cite the artifact (PDF path, screenshot, test name).
- `implemented, unverified` — code written and syntax-checked, but the targeted behavior has not been observed. The next session must verify before any further fixes are layered on top.
- `not implemented` — planned but no code written.

Do NOT use `done`, `fixed`, `complete`, `syntax clean`, or `review passed` as a status. Those describe intermediate states, not whether the bug is gone.

Layering new fixes on top of an `implemented, unverified` change is forbidden without first verifying the prior change. A chain of unverified fixes is a chain of unknowns.
