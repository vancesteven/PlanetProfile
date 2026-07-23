# Perple_X EOS cache-miss crash on the HF app — for Machine A to consider

**Status:** diagnosed, NOT fixed. No code changed. Machine A may have a better approach.
**Date:** 2026-07-21 (genai). **Severity:** breaks all CustomSolution / forward-model
runs on the public Hugging Face Space.

## Symptom (reported from the live HF app)

A custom solution run failed:

```
File ".../PlanetProfile/Thermodynamics/InnerEOS.py", line 160, in __init__
    = np.loadtxt(self.fpath, skiprows=nHeaders, unpack=True)
FileNotFoundError: /home/appuser/.cache/PlanetProfile/Perple_X/
  CM_hydrous_differentiated_Ganymede_Core080Fe020S_excluding_fluid_properties.tab not found.
```

Not an API error. Not a missing repo file.

## Root cause (verified this session)

The requested table **ships in the package** at
`PlanetProfile/Thermodynamics/EOStables/Perple_X/` (13 MB real text table) and is live
on GitHub `main` (raw URL returns HTTP 200). The bug is *where the loader looks*:

1. `PerplexEOSStruct.__init__` (`InnerEOS.py:86`) loads tables **only** from the user
   cache dir `_PERPLEXCACHE = ~/.cache/PlanetProfile/Perple_X/`. Cache miss → straight
   to `np.loadtxt` → `FileNotFoundError`. **No package-local fallback, no
   download-on-miss.**
2. That cache is seeded **only** by `DownloadPerplexFiles()` in `install.py`, called
   **only** by `PPinstall()` (manual `python -m PlanetProfile.install`). **Nothing runs
   it at container/app startup** → empty cache on the HF Space.
3. The HF Docker image ships the in-package `EOStables/Perple_X/` (164 MB, via the
   `PlanetProfile` rsync in `plans/scripts/build_deploy_branch.sh`; only
   `PlanetProfile/Test` is excluded) to `/app/...`, but never copies it into
   `~/.cache/PlanetProfile/Perple_X/`.

**Why it surfaced only now:** the public app's other paths (amortized SBI, precomputed
demos) never touch the silicate EOS. A **CustomSolution forward run is the first path
to call `GetPlanetInnerEOS`** → first cache miss → crash. The GUI compounds it:
`CoreSettings.py:134` builds its EOS dropdown by listing the in-package EOStables dir,
so it offers users exactly the tables the loader can't find in the cache.

**Latent secondary bug:** `DownloadPerplexFiles`'s guard is all-or-nothing
(`if any(*.tab) in cache: return`) — a partially-seeded cache never fetches a single
missing table.

## Candidate approaches (unranked — Machine A to weigh)

- **A. Runtime package-local fallback in the loader (my leaning).** On cache miss in
  `PerplexEOSStruct.__init__`, copy `EOSfname` from
  `os.path.join(_ROOT, 'Thermodynamics', 'EOStables', 'Perple_X')` into `_PERPLEXCACHE`
  (reuse `CopyOnlyIfNeeded` from `PlanetProfile.__init__`), then load. Offline, no
  network, per-file (also dodges the all-or-nothing guard). Put it before the
  `'3D_EOS' in self.fpath` split so `Fe-S_3D_EOS.mat` benefits too. Improve the
  not-found-anywhere error to name both paths + the install hint.
  - *Con:* silently duplicates 164 MB into the cache on first use; adds a copy on the
    cold path.

- **B. Deploy-time cache seeding.** In the Dockerfile (`build_deploy_branch.sh`), copy
  `EOStables/Perple_X/*` → `$HOME/.cache/PlanetProfile/Perple_X/` after `COPY . /app`.
  Warm cache, no runtime copy.
  - *Con:* HF-specific; doesn't help fresh pip installs; grows the image; if done via
    `install.py` it hits the flaky GitHub API at build time.

- **C. Load directly from the package dir, skip the cache entirely** for shipped
  tables (treat `_PERPLEXCACHE` as write-through only for downloaded/generated tables).
  - *Con:* larger change to the loader's path logic and `tableKey` reuse; needs care
    that user-generated tables still resolve.

- **D. Fix `DownloadPerplexFiles` guard + call it at app startup.** Make the guard
  per-file and invoke a seed step when the app boots.
  - *Con:* keeps a network dependency on the critical path (API/raw GitHub); worst for
    offline/air-gapped and for HF cold starts.

## Files involved (for whoever implements)

- `PlanetProfile/Thermodynamics/InnerEOS.py` — `PerplexEOSStruct.__init__`, cache path
  at line 86, loads at 123 (`.mat`) / 160 (`.tab`). `_ROOT` already imported (line 6).
- `PlanetProfile/__init__.py` — `_PERPLEXCACHE` (line 42), `CopyOnlyIfNeeded` (line 19).
- `PlanetProfile/install.py` — `DownloadPerplexFiles()` (line 12), all-or-nothing guard
  (line 16), API/raw-GitHub URLs (26/38).
- `PlanetProfileApp/pages/CoreSettings.py:134` — GUI EOS dropdown source dir.
- `plans/scripts/build_deploy_branch.sh` — HF snapshot/Dockerfile builder (rsync at 24,
  Dockerfile heredoc from 95).

## Constraints to honor

- No change to physics, EOS parsing, units, or interpolants — this is an I/O-path bug.
- Any redeploy of the HF Space is a separate, user-initiated `build_deploy_branch.sh`
  run per `DEPLOYING.md` (opt-in).
- Verify on a cold cache locally (move the table out of `~/.cache/.../Perple_X/`, run
  the Europa `PPTest` default that uses this Ganymede table at
  `PlanetProfile/Default/Europa/PPTest.py:56`) before calling any fix `verified`; HF
  end-to-end stays `implemented, unverified` until observed on the Space.
