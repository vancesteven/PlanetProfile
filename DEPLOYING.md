# Deploying / updating the public PlanetProfileApp

The public GUI is served as a **Hugging Face Docker Space** from an
auto-generated snapshot — never from this repo's dev branches. Snapshot =
app + PlanetProfile package + serve-time data only (~205 MB; the full repo
is ~3 GB and does not fit hosted serving).

Public mode (`PP_PUBLIC_MODE=1`, baked into the Dockerfile) exposes
amortized SBI inference, single forward-model runs (one at a time via a
global run lock), the precomputed demo library, and configuration
save/load as files. MCMC and new exploreogram grids stay hidden (see
`PlanetProfileApp/Utilities/app_mode.py`). The image is conda-based
(miniforge) so reaktoro/CustomSolution compositions work. Server-side
session storage is hidden publicly (shared container); uploaded
configuration files are strictly validated (allowlisted keys, scalar
values, size cap — see `session_manager.validate_session_state`).

## Updating the live app (the usual case)

1. Verify the change on `genai` first (AppTest / unit suite — see
   CLAUDE.md verification discipline).
2. Rebuild + push the snapshot:

   ```bash
   HF_SPACE_ID=vsteven/planetprofile HF_TOKEN=hf_XXXX \
     bash plans/scripts/build_deploy_branch.sh
   ```

   Upload goes through `huggingface_hub.upload_folder` — NOT git push: HF
   rejects git pushes containing files > 10 MiB without LFS (the EOS
   tables are up to 43 MB); the hub API handles large files server-side.
   HF_TOKEN is a **Write** token from hf.co/settings/tokens (revocable
   after the upload; the Space keeps running). The Space rebuilds
   automatically (~5-10 min). Live Space:
   https://huggingface.co/spaces/vsteven/planetprofile
   (app URL: https://vsteven-planetprofile.hf.space).

The script also refreshes the `app-deploy` branch here and, with
`DEPLOY_REPO=...`, the `PlanetProfileApp-deploy` GitHub mirror. The live
app never tracks pushes to `genai`/`main` — redeploys are always this
deliberate script run.

## First-time Space setup

1. huggingface.co → New Space → name `planetprofile`, License of choice,
   SDK **Docker** (blank template), CPU basic, Public. NOTE (policy change
   2026-07): compute Spaces (Docker/Gradio) require a PRO subscription;
   only static Spaces are free. The Space can also be created without the
   web UI: `HfApi(token=...).create_repo('<user>/planetprofile',
   repo_type='space', space_sdk='docker')`.
2. Push the snapshot (command above). The snapshot contains the
   `Dockerfile`, HF YAML frontmatter in `README.md`, and pinned
   `requirements.txt` — no other Space configuration needed.
3. Watch the first build in the Space's Logs tab; confirm the app loads,
   then run the click-through checklist in
   `plans/serve-gui-streamlit-cloud-plan.md` Phase D.2.

## Demo library (precomputed results shipped with the app)

Public visitors cannot compute, so the snapshot ships read-only demo data:

- **Exploreogram grids**: every `output/exploreograms/cache/*.pkl` with a
  `.meta.json` sidecar appears in the page's "Precomputed grid library"
  selector (public mode only). Generate on a dev machine by running the
  grid in the local app (or headlessly via AppTest), then write the
  sidecar (widgets dict + label) — see `save_to_cache(..., meta=...)` in
  `PlanetProfileApp/Utilities/explore_cache.py`. Bigger grids: Machine B.
- **Reference PlanetProfile runs**: `<Body>/` output dirs at repo root
  (profile .txt + `figures/*.pdf`), browsed by the Outputs page (defaults
  to Europa on fresh sessions). Generate with a normal local PP run.

Both are picked up automatically by the snapshot script on the next
deliberate redeploy.

## What's in / out of the snapshot

In: `PlanetProfileApp/`, `PlanetProfile/` (minus `Test/**` except the two
slot structure-grid pkls named by `_SBI_ARTIFACT_SLOTS`), `SPICE/`,
`requirements.txt`, `packages.txt`, `.streamlit/`, generated `README.md` +
`Dockerfile`. Out: `Thermodynamics/` legacy tables (2.5 GB), `papers/`,
`plans/`, `docs/`, Test bulk. Add new serve-time data files to the script's
copy list when a new artifact slot needs them.

## History

Streamlit Community Cloud was the original target (2026-07-14) but app
provisioning failed at the account level (app object never created, no
logs; platform status green; slim repo did not help) — details in
`plans/serve-gui-streamlit-cloud-plan.md`. The `app-deploy` branch and
GitHub mirror remain usable if that is ever revisited.
