# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

For detailed architecture and module structure, see [.claude/ARCHITECTURE.md](.claude/ARCHITECTURE.md).

---

## Current Mission: `genai-clean-port` Branch

Port features from the `genai` branch into `main` via the `genai-clean-port` branch:

- **Surgical** — only bring over well-defined, self-contained functionality
- **Non-breaking** — tests must stay green, no API changes without discussion
- **Reviewable** — one feature = one focused commit

### Branch Strategy

| Branch | Role |
|--------|------|
| `main` | Stable, released code. Target for this port. |
| `genai` | Experimental branch with AI-assisted features. Source. |
| `genai-clean-port` | **Working branch** — all commits land here. |

**Critical rules:**
- Commit exclusively to `genai-clean-port`
- Never commit to `main` or `genai`
- **Never open, draft, or suggest a pull request** — Emma handles the PR
- Never force-push or rebase commits already pushed to remote

---

## Before Every Porting Session

```bash
# 1. Confirm branch
git branch --show-current   # must print "genai-clean-port"

# 2. Sync with main
git fetch origin
git rebase origin/main

# 3. Establish clean baseline
python Testing.py           # ALL tests must pass
```

---

## Porting Workflow

### 1. Identify the Feature
- What does it do physically and computationally?
- Which files in `genai` implement it?
- Any genai-only dependencies (LLM calls, new packages)?
- Conflicts with existing `main` code?

**If not cleanly separable, skip it. Do not port partial code.**

### 2. Diff and Read
```bash
git diff origin/main..origin/genai -- <path/to/file>
```
Understand every changed line before porting.

### 3. Port Minimally
Copy only code needed for the feature. **Do not copy:**
- Debug `print` statements
- Dead code, commented blocks, TODOs
- AI client code, API keys, prompts
- Imports not in `pyproject.toml`

### 4. Add/Update Tests
Every feature needs a test in `PlanetProfile/Test/`.

```bash
python Testing.py          # full suite
```

### 5. Commit Atomically
One feature = one commit:

```
[port] <short imperative description>

Ported from genai branch: <what it does and why it belongs in main>.

Files changed: <list key files>
Tests: <what the test covers>
```

### 6. Self-Review
```bash
git diff origin/main..HEAD   # review everything
python Testing.py            # confirm tests pass
```

---

## Scientific Accuracy (MANDATORY)

**PlanetProfile is peer-reviewed scientific software.** Every ported feature must be physically meaningful and literature-consistent.

### Verify Before Every Commit

- ✅ **Units and dimensions** — Cross-check against `VARIABLE_REFERENCES.md`. SI units unless explicitly stated (e.g. `_km`, `_bar`)
- ✅ **Physical regime** — Model applied only within validated P-T-composition range?
- ✅ **EOS consistency** — Thermodynamic code consistent with EOS hierarchy (SeaFreeze → GSW → Reaktoro → Perple_X). Don't mix EOS outputs without checking phase boundaries
- ✅ **Numerical stability** — No division by zero, log of negatives, ill-conditioned operations
- ✅ **Literature grounding** — Name the paper in commit message and inline comment. Add to `REFERENCES.md` if missing
- ✅ **Conservation laws** — Layer calculations must conserve mass and be pressure-consistent

### Flag Scientific Concerns

If something looks physically suspect:

```python
# SCIENCE-REVIEW NEEDED: <description of concern>
# Reference: <paper or equation if known>
```

Call it out in the commit message.

---

## Things That MUST NEVER Be Ported

- Network calls to LLM APIs (OpenAI, Anthropic, etc.)
- Prompt templates or system-prompt strings
- API keys, tokens, credentials
- Auto-generation of planetary input files via LLM
- Features that bypass `CALC_NEW` / caching logic

**When in doubt, stop and ask.**

---

## Code Conventions

### Variable Naming
- SI units in names: `_m`, `_kg`, `_K`, `_MPa`, `_kgm3`
- Arrays: `P_MPa`, `T_K`, `rho_kgm3`
- Indices: `indsLiq`, `indsI`, `indsSil`
- See `VARIABLE_REFERENCES.md` for canonical names

### Python Style
- **Python 3.11** required for developers
- **NumPy-style docstrings** (Sphinx reads them)
- Use `log.debug/info/warning` — **no bare `print()` calls**
- No f-strings with walrus operators (3.8-3.11 compatibility)
- All physical quantities in SI units internally

### Key Flags in `Params`
- `CALC_NEW` — Force recalculation vs. reload cache
- `DO_PARALLEL` — Enable multiprocessing (must gate all parallel code)
- `CALC_NEW_REF` / `CALC_NEW_INDUCT` / `CALC_NEW_GRAVITY`

### Memory Management
Clear EOS lists between tests to prevent memory issues:
```python
EOSlist.loaded = {}
EOSlist.loaded['CustomSolutionEOS'] = {}
EOSlist.loaded['ReaktoroDatabases'] = {}
EOSlist.ranges = {}
```

---

## Documentation Updates

For every ported feature:
- [ ] Update `CHANGELOG.md` under `[Unreleased]` section
- [ ] Add new `Planet`/`Params` fields to `VARIABLE_REFERENCES.md` with units
- [ ] Update `docs/` if user-facing
- [ ] Add paper to `REFERENCES.md` if citing new literature
- [ ] Build docs: `cd docs && rm -rf stubs/ && make clean && make html`

---

## Dependencies

**No new dependencies without maintainer approval.**

Before adding imports not in `pyproject.toml`:
1. Check if existing utilities cover the need
2. Discuss with maintainers (GitHub issue or Slack)
3. Update `pyproject.toml`, `MANIFEST.in`, README if approved

---

## Quick Checklist Before Every Commit

- [ ] `git branch --show-current` → `genai-clean-port`
- [ ] `python Testing.py` passes with no errors
- [ ] No genai-only imports or LLM calls
- [ ] **Scientific accuracy verified** (units, regime, EOS, stability, literature, conservation)
- [ ] Scientific concerns flagged with `# SCIENCE-REVIEW NEEDED:`
- [ ] Test covers the feature
- [ ] `CHANGELOG.md` updated under `[Unreleased]`
- [ ] `VARIABLE_REFERENCES.md` updated if new fields added
- [ ] `REFERENCES.md` updated if new paper cited
- [ ] Commit follows `[port] ...` format
- [ ] Docs build cleanly if touched
- [ ] **No pull request opened** — Emma handles PR

---

## Useful Commands

```bash
# Run full test suite
python Testing.py

# Run PlanetProfile (smoke test)
python PlanetProfileCLI.py Europa

# Diff genai vs main for specific file
git diff origin/main..origin/genai -- PlanetProfile/Path/File.py

# List all changed files
git diff --name-only origin/main..origin/genai

# Build docs
cd docs && rm -rf stubs/ && make clean && make html

# Build package
python -m build
```

---

## Contacts

- **Owner/Maintainer:** Dr. Steven D. Vance — steven.d.vance@jpl.nasa.gov
- **Lead Developers:** Dr. Mohit Melwani Daswani, Scott Chang (JPL)
- **Port Developer:** Emma Vellard — emma.vellard@outlook.fr
- **Main Repo:** https://github.com/vancesteven/PlanetProfile (submit PRs here)
- **Mirror:** https://github.com/NASA-Planetary-Science/PlanetProfile (do NOT submit PRs)
- **Docs:** https://vancesteven.github.io/PlanetProfile
