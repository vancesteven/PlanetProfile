"""v7 open-|Ae| reference MCMC (pocoMC) — reproducible driver.

Runs the full 21-observable reference posterior with pocoMC against the v7
open-support config (amp_min cut removed; phase_deg_max 70deg all labels). This
is the ratification-blocking reference the v7 SBI flow is crosschecked against.

CRITICAL (spec): the reference MUST use the SAME opened induction_bounds as the
flow. Crosschecking the open flow against the old v5 (cut) reference would fail
by construction — different support = different posterior.

Determinism
-----------
- MCMC seed: config sampler_settings.random_state if present, else --seed
  (default 71, the v7 train seed).
- Reads the pinned v5 (Tb x w) structure cache + Ae sidecar named in the v7
  config (reused unchanged); their sha256 are recorded.

Outputs
-------
- Reference pickle -> a NEW Test54_seawater_v7/ dir (new output; does not alter
  existing Test/ fixtures). Written to /tmp first, then copied to Dropbox.
- A JSONL progress log + a run manifest alongside the pickle.

Run
---
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/v7_reference_mcmc.py
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.mcmc_runner import MCMCRunner

ROOT = Path(__file__).resolve().parents[2]
CFG = ROOT / "PlanetProfile/Inference/configs/europa_clipper_v7_openae_11D.json"
REF_DST = (ROOT / "PlanetProfile/Test/mcmc_results/Europa/Test54_seawater_v7/"
           "europa_clipper_v7_reference_result.pkl")
TMP_DIR = Path("/tmp/v7_build/reference")
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
    ap.add_argument("--dst", default=str(REF_DST))
    args = ap.parse_args()

    TMP_DIR.mkdir(parents=True, exist_ok=True)
    cfg = InferenceConfig.from_json(args.config)
    seed = int((cfg.sampler_settings or {}).get("random_state", args.seed))

    print(f"[ref] config={args.config}")
    print(f"[ref] observables ({len(cfg.observables)}): {list(cfg.observables)}")
    print(f"[ref] induction_bounds: {json.dumps(cfg.induction_bounds)}")
    print(f"[ref] param_space: {list(cfg.param_space.keys())}")
    print(f"[ref] n_effective={cfg.sampler_settings.get('n_effective')}  seed={seed}")

    runner = MCMCRunner(cfg)
    runner.random_state = seed

    cache_path = cfg.structure_cache_path
    sidecar_path = str(cache_path) + ".ae_sidecar.pkl"
    cache_sha = _sha256(cache_path)
    sidecar_sha = _sha256(sidecar_path)
    print(f"[ref] structure cache sha256 = {cache_sha}")
    print(f"[ref] Ae sidecar     sha256 = {sidecar_sha}")

    jsonl = TMP_DIR / "v7_reference_progress.jsonl"
    t0 = time.time()
    print(f"[ref] starting pocoMC ... (progress -> {jsonl})")
    result = runner.run(progress_jsonl_path=str(jsonl))
    dt = time.time() - t0

    samples = np.asarray(result.samples)
    md = getattr(result, "metadata", {}) or {}
    n_samples = int(samples.shape[0])
    log_Z = md.get("log_Z")
    print(f"[ref] done in {dt/60:.1f} min: samples={samples.shape} log_Z={log_Z}")

    tmp_pkl = TMP_DIR / "europa_clipper_v7_reference_result.pkl"
    result.save(str(tmp_pkl))

    manifest = {
        "kind": "v7_openae_reference_mcmc",
        "config": os.path.relpath(args.config, ROOT),
        "config_hash": cfg.generate_hash(),
        "seed": seed,
        "n_effective": cfg.sampler_settings.get("n_effective"),
        "n_samples": n_samples,
        "log_Z": log_Z,
        "log_Z_err": md.get("log_Z_err"),
        "elapsed_time_s": dt,
        "obs_names": list(cfg.observables),
        "param_names": list(result.param_names),
        "induction_bounds": cfg.induction_bounds,
        "structure_cache_path": str(cache_path),
        "structure_cache_sha256": cache_sha,
        "ae_sidecar_sha256": sidecar_sha,
        "git_sha": _git_sha(),
        "prior_note": ("OPEN-|Ae| config: amp_min support cut REMOVED, "
                       "phase_deg_max 70deg for synodic/synodic-2nd/orbital. "
                       "D_iceIh uniform[5,80] km, log10_w uniform[-1,2]. Same "
                       "opened bounds as the v7 flow (crosscheck validity)."),
    }
    tmp_manifest = TMP_DIR / "v7_reference_manifest.json"
    with open(tmp_manifest, "w") as f:
        json.dump(manifest, f, indent=2, default=str)

    dst = Path(args.dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(tmp_pkl, dst)
    shutil.copy2(tmp_manifest, dst.parent / "v7_reference_manifest.json")
    shutil.copy2(jsonl, dst.parent / "v7_reference_progress.jsonl")
    print(f"[ref] reference pickle -> {dst}")
    print(f"[ref] manifest         -> {dst.parent / 'v7_reference_manifest.json'}")
    print(f"[ref] progress log     -> {dst.parent / 'v7_reference_progress.jsonl'}")


if __name__ == "__main__":
    main()
