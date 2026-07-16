#!/usr/bin/env bash
# Build the slim orphan `app-deploy` branch for Streamlit Community Cloud.
#
# Assembles a serve-time-only snapshot (~300 MB vs the 3.0 GB full tree:
# excludes root Thermodynamics/Perple_X 2.5 GB, papers/, plans/, docs/,
# and PlanetProfile/Test bulk except the structure-grid pkls the GUI slot
# registry names), commits it as a single orphan commit, and pushes it to
# origin/app-deploy (force: the branch has no history by design).
#
# Usage: plans/scripts/build_deploy_branch.sh [--no-push]
# Rerun after verifying genai — this is the ONLY redeploy path; the public
# app never tracks dev pushes.
set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel)"
SRC_SHA="$(git -C "$REPO_ROOT" rev-parse --short HEAD)"
STAGE="$(mktemp -d /tmp/pp-app-deploy.XXXXXX)"
trap 'rm -rf "$STAGE"' EXIT

echo "Staging deploy snapshot of ${SRC_SHA} in ${STAGE}"

(
  cd "$REPO_ROOT"
  rsync -a --exclude='__pycache__' --exclude='*.pyc' \
    --exclude='PlanetProfile/Test' \
    PlanetProfile PlanetProfileApp SPICE \
    PlanetProfileCLI.py \
    requirements.txt packages.txt .streamlit LICENSE \
    "$STAGE"/
  # App-side config defaults (GetConfig copies these into place on first
  # touch; shipping them avoids a first-run write in a fresh container).
  if [[ -d UserConfigs ]]; then
    rsync -a --exclude='__pycache__' UserConfigs "$STAGE"/
  fi
  # Demo library (user request 2026-07-15): precomputed exploreogram grids
  # (loaded read-only by the public app's demo-library selector) and
  # reference PlanetProfile run outputs (browsed by the Outputs page).
  # Regenerate with plans/scripts (gen scripts documented in DEPLOYING.md).
  if [[ -d output/exploreograms/cache ]]; then
    rsync -a output/exploreograms/cache "$STAGE"/output/exploreograms/
  fi
  for body in Europa; do
    if [[ -d "$body" ]]; then
      rsync -a --exclude='inductionData' "$body" "$STAGE"/
    fi
  done
)

# Serve-time data from Test/: only the structure grids the
# _SBI_ARTIFACT_SLOTS cache_paths reference (plus their dirs).
for f in \
  "PlanetProfile/Test/mcmc_results/Europa/Test51_seawater/europa_seawater_structure_grid.pkl" \
  "PlanetProfile/Test/mcmc_results/Titan/Test50_andrade_noocean_yao2014/titan_allice_yao2014_structure_grid.pkl" \
  ; do
  mkdir -p "$STAGE/$(dirname "$f")"
  cp "$REPO_ROOT/$f" "$STAGE/$f"
done

# README with Hugging Face Spaces YAML frontmatter (harmless on GitHub).
# HF retired the native streamlit SDK; this deploys as a Docker Space.
cat > "$STAGE/README.md" <<EOF
---
title: PlanetProfile
emoji: 🪐
colorFrom: blue
colorTo: indigo
sdk: docker
app_port: 7860
pinned: false
short_description: Icy-moon interior structure + Bayesian inference GUI
---

# PlanetProfileApp — deployment snapshot

Auto-generated serve-only snapshot of
[PlanetProfile](https://github.com/vancesteven/PlanetProfile) genai
commit ${SRC_SHA} (plans/scripts/build_deploy_branch.sh). Do not develop
here — changes belong on the source branch. Update procedure:
DEPLOYING.md in the source repo.

Entrypoint: \`PlanetProfileApp/PlanetProfileApp.py\` · public demo mode
(\`PP_PUBLIC_MODE=1\`, baked into the Dockerfile): amortized inference
only, heavy compute hidden.
EOF

# Dockerfile for Hugging Face Docker Space (uid-1000 non-root per HF).
# Conda-based (miniforge) so reaktoro — conda-forge-only, no pip wheel —
# ships in the image and CustomSolution ocean compositions work publicly.
cat > "$STAGE/Dockerfile" <<'EOF'
FROM condaforge/miniforge3:latest

# poppler: pdf2image previews. texlive set: the PlanetProfile usetex
# preamble (UserConfigs/configPPplots.py) needs siunitx + mhchem
# (texlive-science) AND the stix font + upgreek (texlive-fonts-extra —
# no smaller Debian package carries stix.sty); cm-super + dvipng for
# matplotlib's usetex pipeline.
RUN apt-get update \
 && apt-get install -y --no-install-recommends poppler-utils \
      texlive-latex-base texlive-latex-extra texlive-science \
      texlive-fonts-recommended texlive-fonts-extra cm-super dvipng \
 && rm -rf /var/lib/apt/lists/*

# HF runs Spaces as uid 1000. Ubuntu-Noble bases (miniforge) already ship
# an 'ubuntu' user at uid 1000 — creating another fails the build
# ("UID 1000 is not unique"); ensure one exists and give it a writable HOME.
RUN (id -u 1000 >/dev/null 2>&1 || useradd -m -u 1000 user) \
 && mkdir -p /home/appuser && chown 1000 /home/appuser

# reaktoro from conda-forge (pip cannot install it). A DEDICATED env —
# never mutate miniforge's base (its pinned newer python makes a
# python=3.11 downgrade unresolvable, failing the build).
RUN mamba create -y -n pp python=3.11 reaktoro -c conda-forge \
 && mamba clean -afy
ENV PATH=/opt/conda/envs/pp/bin:$PATH

WORKDIR /app
COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

COPY --chown=1000 . /app
# COPY --chown owns the contents, but the /app directory node itself is
# root's (created by WORKDIR): the app must be able to create dirs in its
# CWD (UserConfigs/, run outputs <Body>/, output/exploreograms/).
RUN chown 1000 /app
USER 1000
ENV HOME=/home/appuser

# Public demo mode (see PlanetProfileApp/Utilities/app_mode.py):
# amortized inference + single capped forward runs; MCMC and new
# exploreogram grids stay disabled.
ENV PP_PUBLIC_MODE=1

EXPOSE 7860
CMD ["streamlit", "run", "PlanetProfileApp/PlanetProfileApp.py", \
     "--server.port=7860", "--server.address=0.0.0.0", "--server.headless=true"]
EOF

du -sh "$STAGE"

cd "$STAGE"
git init -q -b app-deploy
git add -A
git -c user.name="deploy-script" -c user.email="noreply@planetprofile" \
  commit -q -m "deploy snapshot from genai ${SRC_SHA}"

echo "Snapshot commit: $(git rev-parse --short HEAD) ($(git ls-files | wc -l | tr -d ' ') files)"

# Push into the local repo (branch is never checked out there), then to origin.
git push -q -f "$REPO_ROOT" app-deploy:app-deploy
if [[ "${1:-}" != "--no-push" ]]; then
  git -C "$REPO_ROOT" push -f origin app-deploy
  echo "Pushed origin/app-deploy (source ${SRC_SHA})"
else
  echo "Local app-deploy updated; skipped origin push (--no-push)"
fi

# Contingency C.3 (plan): Streamlit Cloud clones the ENTIRE repo pack
# (~885 MB across all branches), which can kill provisioning silently.
# Set DEPLOY_REPO to also push the snapshot as `main` of a dedicated
# lightweight repo, e.g.:
#   DEPLOY_REPO=https://github.com/vancesteven/PlanetProfileApp-deploy.git \
#     plans/scripts/build_deploy_branch.sh
if [[ -n "${DEPLOY_REPO:-}" ]]; then
  git push -f "$DEPLOY_REPO" app-deploy:main
  echo "Pushed snapshot to ${DEPLOY_REPO} main (source ${SRC_SHA})"
fi

# Hugging Face Docker Space (primary serving target — see DEPLOYING.md).
# Uploaded via huggingface_hub (NOT git push: HF's pre-receive hook rejects
# files > 10 MiB without LFS, and the EOS tables are up to 43 MB;
# upload_folder handles large files server-side with no local git-lfs).
# Usage:
#   HF_SPACE_ID=<user>/planetprofile HF_TOKEN=hf_xxx \
#     plans/scripts/build_deploy_branch.sh
if [[ -n "${HF_SPACE_ID:-}" ]]; then
  : "${HF_TOKEN:?HF_TOKEN (write token from hf.co/settings/tokens) is required}"
  PYBIN="${HF_PYTHON:-/opt/miniconda3/envs/PP/bin/python}"
  STAGE="$STAGE" HF_SPACE_ID="$HF_SPACE_ID" HF_TOKEN="$HF_TOKEN" SRC_SHA="$SRC_SHA" \
  "$PYBIN" - <<'PYEOF'
import os
from huggingface_hub import HfApi
api = HfApi(token=os.environ['HF_TOKEN'])
api.upload_folder(
    folder_path=os.environ['STAGE'],
    repo_id=os.environ['HF_SPACE_ID'],
    repo_type='space',
    commit_message=f"deploy snapshot from genai {os.environ['SRC_SHA']}",
    delete_patterns=['**'],   # clean sync: remove files dropped upstream
    ignore_patterns=['.git/*', '.git*'],
)
print(f"Uploaded snapshot to Space {os.environ['HF_SPACE_ID']} "
      f"(source {os.environ['SRC_SHA']})")
PYEOF
fi
