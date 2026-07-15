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
    requirements.txt packages.txt .streamlit LICENSE \
    "$STAGE"/
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

cat > "$STAGE/README.md" <<EOF
# PlanetProfileApp — deployment snapshot

Auto-generated serve-only snapshot of
[PlanetProfile](https://github.com/vancesteven/PlanetProfile) genai
commit ${SRC_SHA} (plans/scripts/build_deploy_branch.sh). Do not develop
here — changes belong on the source branch.

Entrypoint: \`PlanetProfileApp/PlanetProfileApp.py\` · public mode via
\`PP_PUBLIC_MODE=1\` in app secrets (amortized inference only).
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
