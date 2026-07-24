"""Titan free-gravity (no-ocean) reference MCMC (pocoMC) — reproducible driver.

Runs the full {C20, C22, Re_k2, Im_k2} reference posterior with pocoMC against
titan_freegrav_noocean.json. This is the ratification-blocking reference the
no-ocean SBI flow is crosschecked against (validate_sbi crosscheck).

The existing test52 no-ocean production pickle CANNOT serve as this reference:
it conditions on {Re_k2, Im_k2, CMR2}, not the free-gravity observable set
{C20, C22, Re_k2, Im_k2}. A fresh reference on the new observables is required.

Determinism
-----------
- MCMC seed from config sampler_settings.random_state if present, else --seed
  (default 71, the Titan no-ocean train seed).
- Forward model reads the pinned test52 no-ocean 1D structure cache named in the
  config; its sha256 is recorded. (No Ae sidecar exists for this cache; _sha256
  returns None gracefully.)

Outputs
-------
- Reference pickle -> Test/mcmc_results/Titan/Test52_andrade_noocean_diff/
  freegrav_reference/ (NEW subdir; does NOT touch the existing production_run/
  pickle). Written to /tmp first, then copied to Dropbox.
- Progress JSONL + run manifest (seed, config_hash, cache sha256, log_Z, elapsed,
  git_sha, sample count).

NOTE: writing under Test/ — the repo rule forbids ALTERING existing Test/ files
without permission. This creates a NEW subdirectory + files only; it does not
modify any existing Test/ file. If even new-file creation under Test/ needs
sign-off, pass --dst to redirect outside Test/.

Run:
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/titanG_reference_mcmc.py
"""
import argparse, hashlib, json, os, shutil, subprocess, time
from pathlib import Path

import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.mcmc_runner import MCMCRunner

ROOT = Path(__file__).resolve().parents[2]
CFG = ROOT / "PlanetProfile/Inference/configs/titan_freegrav_noocean.json"
REF_DST = (ROOT / "PlanetProfile/Test/mcmc_results/Titan/"
           "Test52_andrade_noocean_diff/freegrav_reference/"
           "titan_freegrav_noocean_reference_result.pkl")
TMP_DIR = Path("/tmp/titanG_build/reference")
DEFAULT_SEED = 71


def _sha256(path):
    if not path or not os.path.exists(path):
        return None
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _git_sha():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip()
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(CFG))
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--dst", default=str(REF_DST),
                    help="final pickle path; written via /tmp first")
    args = ap.parse_args()

    TMP_DIR.mkdir(parents=True, exist_ok=True)
    cfg = InferenceConfig.from_json(args.config)
    seed = int((cfg.sampler_settings or {}).get("random_state", args.seed))

    print(f"[ref] config={args.config}")
    print(f"[ref] observables ({len(cfg.observables)}): {list(cfg.observables)}")
    print(f"[ref] param_space: {list(cfg.param_space.keys())}")
    print(f"[ref] n_effective={cfg.sampler_settings.get('n_effective')}  seed={seed}")

    runner = MCMCRunner(cfg)
    runner.random_state = seed

    cache_path = cfg.structure_cache_path
    cache_sha = _sha256(cache_path)
    print(f"[ref] structure cache sha256 = {cache_sha}")

    jsonl = TMP_DIR / "titanG_reference_progress.jsonl"
    t0 = time.time()
    print(f"[ref] starting pocoMC ... (progress -> {jsonl})")
    result = runner.run(progress_jsonl_path=str(jsonl))
    dt = time.time() - t0

    samples = np.asarray(result.samples)
    md = getattr(result, "metadata", {}) or {}
    n_samples = int(samples.shape[0])
    log_Z = md.get("log_Z")
    print(f"[ref] done in {dt/60:.1f} min: samples={samples.shape} log_Z={log_Z}")

    tmp_pkl = TMP_DIR / "titan_freegrav_noocean_reference_result.pkl"
    result.save(str(tmp_pkl))

    manifest = {
        "kind": "titanG_noocean_reference_mcmc",
        "config": os.path.relpath(args.config, ROOT),
        "config_hash": cfg.generate_hash(),
        "seed": seed,
        "n_effective": cfg.sampler_settings.get("n_effective"),
        "n_samples": n_samples,
        "log_Z": log_Z, "log_Z_err": md.get("log_Z_err"),
        "elapsed_time_s": dt,
        "obs_names": list(cfg.observables),
        "param_names": list(result.param_names),
        "structure_cache_path": str(cache_path),
        "structure_cache_sha256": cache_sha,
        "git_sha": _git_sha(),
        "note": ("Titan free-gravity no-ocean reference. Observable set "
                 "{C20,C22,Re_k2,Im_k2} (induction+h2 dropped). Distinct from "
                 "the test52 production_run pickle, which conditions on "
                 "{Re_k2,Im_k2,CMR2}."),
    }
    tmp_manifest = TMP_DIR / "titanG_reference_manifest.json"
    with open(tmp_manifest, "w") as f:
        json.dump(manifest, f, indent=2, default=str)

    dst = Path(args.dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(tmp_pkl, dst)
    shutil.copy2(tmp_manifest, dst.parent / "titanG_reference_manifest.json")
    shutil.copy2(jsonl, dst.parent / "titanG_reference_progress.jsonl")
    print(f"[ref] reference pickle -> {dst}")
    print(f"[ref] manifest         -> {dst.parent / 'titanG_reference_manifest.json'}")


if __name__ == "__main__":
    main()
