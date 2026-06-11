# CLAUDE.md

This file provides operating rules for Claude Code when working in this repository.

Primary goals:
- preserve scientific correctness,
- avoid silent changes to physical assumptions,
- keep code changes auditable,
- conserve context,
- produce manuscript-defensible evidence for numerical/modeling work.

## Hard constraints

- Use the `PPcl` mamba environment. Never use system Python.
  ```bash
  mamba activate PPcl
  ```

-Primary working directory: ~/Dropbox/planetprofile-genai/

Manuscripts directory: ~/Dropbox/manuscripts/

Treat manuscripts as read-only unless explicitly told to edit. Never commit or push manuscript changes without explicit instruction.

Do not scan, grep, summarize, index, or modify Perple_X directories unless explicitly instructed.
Do not modify .claudeignore or .graphifyignore without explicit approval.
Do not commit, push, merge, rebase, reset, tag, or create releases without explicit approval.
Do not modify numbered PPTest* files unless explicitly approved.
Do not silently change physical assumptions, constants, parameter meanings, units, priors, likelihoods, validation criteria, defaults, or saved-output formats.

If required input is missing, inaccessible, or ambiguous, say exactly what is missing. Do not infer unseen diffs, files, logs, benchmark outputs, or user intent.

-Scientific modeling rules

PlanetProfile is scientific software. Passing tests is necessary but not sufficient.

- For PlanetProfile, PMC, MCMC, and other modeling work:

Preserve reproducibility.
Document equations, units, priors, likelihoods, parameter bounds, validation data, and numerical tolerances.
Distinguish clearly between inputs, priors, likelihood terms, deterministic forward-model outputs, posterior samples, and derived quantities.
Do not confuse outputs such as Love numbers, moment of inertia, mass, gravity coefficients, or magnetic induction responses with priors unless they are explicitly used that way.
Flag anything that should be disclosed in a manuscript.
If numerical results change, identify whether the cause is a physical/modeling change, bug fix, default change, tolerance change, data/benchmark change, or unresolved issue.

- For modeling or inference changes, check where applicable:

unit consistency,
dimensional consistency,
expected monotonicity or limiting behavior,
benchmark agreement,
reproducibility with fixed seeds,
likelihood/prior behavior,
physically plausible outputs,
manuscript-reportable limitations.

If tests pass but the result is physically suspicious, stop and surface the concern.

Autonomy and checkpoints

Claude runs autonomously with guidance from a opus or better orchestration agent. 
The main sonnet run and any associated agents need to check with the more capable orchestrator

A change is nontrivial if it:

touches more than one file,
changes more than about 50 lines,
changes scientific/modeling logic,
changes numerical outputs,
changes inference behavior,
changes public APIs,
adds dependencies,
modifies tests to accommodate changed behavior,
changes defaults,
changes saved outputs,
affects manuscript-facing figures or claims.

For nontrivial changes, produce a plan before editing.

The plan must include:

files to touch,
files not to touch,
intended behavior change,
scientific or numerical assumptions,
tests or benchmarks to run,
success criteria,
known risks.

After you check and  approve the plan:

implement only the approved plan,
run targeted tests,
self-debug within the approved scope,
avoid unrelated cleanup,
re-plan and ask before expanding scope.

Stop and ask for approval from the human running claude before continuing if the work requires:

touching new modules,
adding dependencies,
changing scientific assumptions,
changing numerical outputs beyond expectations,
rewriting tests,
changing benchmark data,
changing saved-output formats.

Before yielding control, provide:

Summary.
Files changed.
Tests run and results.
Numerical/scientific impact.
Deviations from plan.
Remaining risks or follow-ups.

Then stop. Do not commit or push.

##Test and benchmark discipline

Do not modify tests merely to make failures disappear.

If changing tests:

explain why the old expectation was wrong or obsolete,
preserve or improve scientific coverage,
add regression tests for bugs,
do not loosen tolerances without numerical justification,
do not replace meaningful checks with weaker smoke tests,
do not delete failing tests without approval.

##For stochastic methods:

fix random seeds where practical,
report sampling settings,
separate stochastic variability from code-induced changes,
do not claim posterior convergence without diagnostics.
Module-specific rules

##PlanetProfile/Inference/*

Implement approved designs and bug fixes autonomously within scope.
Always plan first for changes to priors, likelihoods, samplers, parameter transforms, forward-model coupling, saved outputs, posterior interpretation, or numerical results.
Ensure individual models in the datasets can easily be referenced back to full planetprofile models. A philosophy of the Gui for Inference is that one should be able to inspect the outputs in a dedicate tab for that.

##Core PlanetProfile modules

Always plan first before editing core modules, including Thermodynamics/, Main.py, Plotting/, MagneticInduction/, Gravity/, hydrosphere/interior-structure code, EOS interfaces, and forward-model code.
When plotting new science, use existing plotting infrastructure if possible. For example, do not create new wedge plotting functionality when you can use the existing wedge plot functions.

##PlanetProfileApp/

Plan first for new pages, session-state changes, model execution changes, file I/O changes, or changes affecting scientific output.
Autonomous edits are acceptable for cosmetic changes, help text, typo fixes, and UI wording that does not alter behavior.
Repository navigation

##For substantial PlanetProfile work, orient before editing.

Use graphify for architectural questions:

graphify query "inference mcmc forward model parameter space"

Use llm-wiki when broader project context is needed:

cd ~/llm-wiki && python3 -m llmwiki query "PlanetProfile sbi PlanetProfileApp 3D viscosity convection" && cd -

Read README.md only when repository-level assumptions matter or session context is stale.

Skip startup queries for narrow bug fixes when the relevant file path and failure are already known.

Use:

graphify query "topic" for architecture, flow, and cross-module symbol location,
grep or ripgrep for exact strings or known-symbol searches,
targeted file reads instead of broad file dumps.

After modifying code structure, run:

graphify update .

Do not modify .graphifyignore without approval.

Context and model discipline

Context is limited scientific working memory.

Summarize context before passing it to subagents.
Do not pass full conversation history to subagents unless necessary.
Avoid re-reading files already in context.
Do not spawn a subagent for tasks completable in fewer than about five tool calls.
Prefer targeted searches and summaries over full file dumps.
Avoid dumping long logs into context.
Preserve task framing when delegating.

Use the strongest available model for:

scientific assumptions,
inference architecture,
likelihood/prior design,
manuscript-facing conclusions,
numerical-output changes,
unexpected benchmark changes.

Use cheaper or faster models for:

file summaries,
grep/search,
git status,
log triage,
mechanical edits,
simple refactors.

Do not rely on model memory for APIs, equations, parameter meanings, file contents, or benchmark behavior. Verify against source files, tests, documentation, or benchmark outputs.

For bug fixes, inspect the relevant source before editing. Never guess API names or behavior from memory.

Multi-agent failure modes

Watch for:

context collapse: subagents lose task framing,
specification gaming: tests pass without real correctness,
goal misgeneralization: scope expands beyond the approved task,
false context assumption: agents assume unseen diffs, files, logs, or results exist.

If any of these appear, stop and surface the issue before continuing.

Missing context rule

If required context is unavailable, say so directly.

Do not infer:

unseen diffs,
unprovided logs,
hidden files,
expected outputs,
benchmark results,
manuscript claims,
physical assumptions.

Provide the smallest safe next step, such as the exact file, command, snippet, benchmark, or inspection plan needed.

Final report format

End every coding session with:
```
## Summary

## Files changed

## Tests run

## Numerical/scientific impact

## Deviations from plan

## Remaining risks and follow-ups
```
Then stop.

Do not commit or push.

Architecture and testing

See PP_ARCHITECTURE_TESTING.md for repository architecture, test organization, and benchmark guidance.

If this file and PP_ARCHITECTURE_TESTING.md conflict, follow the stricter rule and surface the conflict.
