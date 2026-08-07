# Deploying / updating the public PlanetProfileApp

The public GUI is served as a **Hugging Face Docker Space** from an
auto-generated snapshot — never from this repo's dev branches. Snapshot =
app + PlanetProfile package + serve-time data only (~300 MB; the full repo
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

   Two verification modes avoid an origin deployment while testing:

   - `--build-only` assembles and verifies the temporary snapshot, then exits
     before updating any local or remote deploy branch or uploading to the
     Space.
   - `--no-push` updates the local `app-deploy` branch but skips the
     `origin/app-deploy` push. Explicit `DEPLOY_REPO` or Hugging Face upload
     settings still take effect, so leave those unset for a local-only check.

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

## Release discipline for inference artifacts

New SBI/amortized artifacts (e.g. the Europa salinity campaign) deploy to
the public app ONLY after they have landed on the main repository
(origin/genai: artifact + gate reports + INDEX row + GUI slot, verified per
CLAUDE.md) — never directly from a compute machine. The public app is a
consumer of the repo, not a second source of truth.

## Local development notes

- PDF previews on the Outputs page need the poppler binaries
  (`pdftoppm`): `conda install -n PP -c conda-forge poppler`. The HF
  image installs poppler-utils via apt. Without it the app degrades to
  download buttons.
- The speciation-vs-depth figure (PlotHydrosphereSpecies) is produced
  automatically by PlanetProfile for CustomSolution ocean compositions
  (Params.PLOT_SPECIES_HYDROSPHERE, on by default) — Seawater/MgSO4/NaCl
  runs do not produce it by design.

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

In: `PlanetProfileApp/`, `PlanetProfile/` (minus `Test/**` except the six
structure-grid caches derived from non-placeholder `cache_path` entries in
`_SBI_ARTIFACT_SLOTS`), `SPICE/`, `requirements.txt`, `packages.txt`,
`.streamlit/`, generated `README.md` + `Dockerfile`. The current registry set
is:

- `PlanetProfile/Test/mcmc_results/Europa/Test51_seawater/europa_seawater_structure_grid.pkl`
- `PlanetProfile/Test/mcmc_results/Europa/Test52_seawater_v3/europa_seawater_structure_grid_v3_2d.pkl`
- `PlanetProfile/Test/mcmc_results/Europa/Test52_seawater_v5/europa_seawater_structure_grid_v5_2d.pkl`
- `PlanetProfile/Test/mcmc_results/Titan/Test50_andrade_noocean_yao2014/titan_allice_yao2014_structure_grid.pkl`
- `PlanetProfile/Test/mcmc_results/Titan/Test52_andrade_noocean_diff/titan_diff_noocean_structure_grid.pkl`
- `PlanetProfile/Test/mcmc_results/Titan/Test54_nh3_ocean/titan_nh3_joint_structure_grid_2d.pkl`

For every derived primary cache, the script also copies any adjacent
`.ae_sidecar.pkl` and `_offsets.json` sidecars. It then compares the staged
primary PKLs with the registry-derived list and aborts unless the sets are
exactly equal (`plans/scripts/build_deploy_branch.sh:55-96`). New artifact
slots therefore declare `cache_path` in the registry; there is no manual cache
copy list to maintain.

Out: `Thermodynamics/` legacy tables (2.5 GB), `papers/`, `plans/`, `docs/`,
and the remainder of `PlanetProfile/Test/`.

## History

Streamlit Community Cloud was the original target (2026-07-14) but app
provisioning failed at the account level (app object never created, no
logs; platform status green; slim repo did not help) — details in
`plans/serve-gui-streamlit-cloud-plan.md`. The `app-deploy` branch and
GitHub mirror remain usable if that is ever revisited.
