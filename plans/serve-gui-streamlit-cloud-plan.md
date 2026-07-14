# Serve PlanetProfileApp on Streamlit Community Cloud (public, amortized-only)

Author: Machine A, 2026-07-14. User decisions: target = **Streamlit Community
Cloud** (share.streamlit.io, free public URL); heavy modes = **amortized-only
public** (MCMC + exploreogram compute hidden from visitors). All phases are
light work (Machine A). Status: `not implemented` throughout.

## Constraints (measured 2026-07-14)

- Cloud gives ~2.7 GB RAM, shared CPU, deploys from a public GitHub
  repo/branch, no git-LFS. Repo working set is **3.0 GB tracked**, dominated
  by root-level `Thermodynamics/Perple_X` (2.5 GB legacy tables the app never
  reads at serve time) — a dedicated slim deploy branch is required.
- Serve-time payload is small: SBI artifacts 1.2 MB total, structure-grid
  pkls ~1.2 MB, `SPICE/` 360 KB, `PlanetProfile/MagneticInduction` 8.3 MB,
  in-package `PlanetProfile/Thermodynamics` 191 MB (import-time deps; keep).
- Pip on Linux resolves plain `torch` to CUDA wheels (~2.5 GB) — must pin the
  CPU wheel. CPU wheel reports `torch.__version__ == '2.8.0+cpu'`, which
  breaks the artifact `validated_version_pairs` string match — Phase A.3.
- macOS libomp patching (memory: pp-env-libomp-fixes) is NOT needed on Linux;
  TidalPy 0.7.4 / CyRK install clean from PyPI.

## Phase A — public-mode gating (code)

1. Mode flag: `_public_mode()` helper — true when `st.secrets` or env sets
   `PP_PUBLIC_MODE=1`. Local runs unaffected (flag absent).
2. When public: Inference page forces amortized exec mode (MCMC radio hidden,
   banner explains: "public demo — full MCMC available in the repo");
   n_reeval options capped to [100, 250] (default 100 — 500 forward evals is
   minutes on Cloud's shared CPU); n_posterior_samples ≤ 20,000. Exploreogram
   run button disabled with the same banner. Browse-saved-runs stays.
3. Torch version normalization: in `sbi_runner._load_artifact_dict`, compare
   versions with the local build suffix stripped (`'2.8.0+cpu' -> '2.8.0'`)
   for BOTH the warning check and the validated-pairs demotion, keeping the
   raw string in logs. Unit test: `+cpu` suffix matches validated pair; a
   genuinely different version still warns.
4. Verify (AppTest, PP_PUBLIC_MODE=1): MCMC radio absent, caps in effect,
   Europa amortized run completes, guards + widget reseed still work.

## Phase B — packaging

1. Root `requirements.txt` (Cloud reads it): `--extra-index-url
   https://download.pytorch.org/whl/cpu`, then pins matching the validated
   env: torch==2.8.0, sbi==0.26.1, TidalPy==0.7.4, CyRK, streamlit==1.59.1,
   numpy==2.4.3, scipy==1.17.1, matplotlib, corner, pandas, seafreeze, gsw,
   spiceypy. Derive the exact closure by scraping app imports against
   `pip list` in the PP env; fail the phase if anything imports outside it.
2. `.streamlit/config.toml`: theme + `server.maxMessageSize` bump (large
   matplotlib payloads).
3. Cloud advanced settings: Python 3.11.

## Phase C — slim deploy branch

1. Script `plans/scripts/build_deploy_branch.sh`: assembles an ORPHAN branch
   `app-deploy` (single commit, no history) containing only:
   `PlanetProfileApp/`, `PlanetProfile/` (with `Test/**` pruned to the two
   slot structure-grid pkls the `_SBI_ARTIFACT_SLOTS` cache_paths name),
   `SPICE/`, `requirements.txt`, `.streamlit/`, `LICENSE`, minimal README.
   Excludes root `Thermodynamics/`, `papers/`, `plans/`, `docs/`, Test bulk.
   Expected clone ~250-300 MB. Re-running the script + `push -f app-deploy`
   is the ONLY redeploy path — public app is isolated from genai dev churn.
2. Local verification gate (blocking, per CLAUDE.md): bare-clone `app-deploy`
   to a temp dir, fresh venv from requirements.txt, `streamlit run` +
   AppTest click-through (Europa amortized end-to-end, plots, guards), and
   measure peak RSS with psutil — must stay < ~2 GB or Cloud will OOM.
3. Contingency if Cloud's clone still drags the full 885 MB pack: move
   `app-deploy` to a separate deploy repo (`PlanetProfileApp-deploy`) pushed
   by the same script; decide only after the first deploy attempt.

## Phase D — deploy + cloud verification

1. share.streamlit.io → New app → repo `vancesteven/PlanetProfile`, branch
   `app-deploy`, entrypoint `PlanetProfileApp/PlanetProfileApp.py`,
   Python 3.11. Pick a stable subdomain.
2. Cloud click-through checklist (manual; AppTest can't run there): cold-start
   time; Europa amortized run end-to-end; layer-thickness + heating + k2 +
   C/MR² panels render; Im_k2=0.18 guard refuses; slot switch reseeds
   defaults; NO torch version warning (validated pair + `+cpu` normalization);
   a second browser session runs while the first is active.
3. Note sleep/wake behavior (Cloud sleeps idle apps; first visitor wakes it,
   ~1-2 min cold start) in the app's intro page.

## Phase E — ops + docs

- README: public URL badge + "public demo is amortized-only" note.
- Handoff + INDEX pointer to this plan; record deployed commit SHA in the
  deploy script output.
- Redeploys are deliberate: verify on genai first (AppTest), then rebuild
  `app-deploy`.

## Risks

- RAM: torch+TidalPy+streamlit RSS unmeasured — Phase C.2 measures before any
  deploy. If > ~2 GB: drop to Hugging Face Spaces (16 GB) with the same
  branch, one-line change of target.
- CPU: shared vCPU makes 100 forward evals ~30-90 s; caps chosen accordingly.
- Version drift: Cloud rebuilds on requirements change only; pins keep the
  torch pair inside the validated 2.11→2.8 envelope.
