"""Titan free-gravity NH3 JOINT reference MCMC (pocoMC) — reproducible driver.

Runs the full {C20, C22, Re_k2, Im_k2} reference posterior with pocoMC against
test54_titan_nh3_freegrav.json (the JOINT no-ocean + ocean NH3 config). This is
the ratification-blocking reference the NH3 SBI flow is crosschecked against
(validate_sbi crosscheck).

Determinism
-----------
- MCMC seed from config sampler_settings.random_state if present, else --seed
  (default 72, the NH3 train seed).
- Forward model reads the pinned 2D NH3 structure cache named in the config; its
  sha256 is recorded.

Outputs
-------
- Reference pickle -> validation_reports/titan_freegrav_nh3_1m/reference/
  (OUTSIDE Test/ — no Test/ write). Written to /tmp first, then copied.
- Progress JSONL + run manifest (seed, config_hash, cache sha256, log_Z, elapsed,
  git_sha, sample count).

Run:
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/titanG_nh3_reference_mcmc.py
"""
import argparse, hashlib, json, os, shutil, subprocess, time
from pathlib import Path

import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.mcmc_runner import MCMCRunner

ROOT = Path(__file__).resolve().parents[2]
CFG = ROOT / "PlanetProfile/Inference/configs/test54_titan_nh3_freegrav.json"
REF_DST = (ROOT / "validation_reports/titan_freegrav_nh3_1m/reference/"
           "titan_freegrav_nh3_reference_result.pkl")
TMP_DIR = Path("/tmp/titanG_build/nh3_reference")
DEFAULT_SEED = 72


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

    jsonl = TMP_DIR / "titanG_nh3_reference_progress.jsonl"
    t0 = time.time()
    print(f"[ref] starting pocoMC ... (progress -> {jsonl})")
    result = runner.run(progress_jsonl_path=str(jsonl))
    dt = time.time() - t0

    samples = np.asarray(result.samples)
    md = getattr(result, "metadata", {}) or {}
    n_samples = int(samples.shape[0])
    log_Z = md.get("log_Z")
    print(f"[ref] done in {dt/60:.1f} min: samples={samples.shape} log_Z={log_Z}")

    tmp_pkl = TMP_DIR / "titan_freegrav_nh3_reference_result.pkl"
    result.save(str(tmp_pkl))

    manifest = {
        "kind": "titanG_nh3_joint_reference_mcmc",
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
        "note": ("Titan free-gravity NH3 JOINT no-ocean+ocean reference. "
                 "Observable set {C20,C22,Re_k2,Im_k2}. 2D (Tb x w) NH3 cache "
                 "with retry_frozen_as_no_ocean=True; log10_wOcean_ppt axis."),
    }
    tmp_manifest = TMP_DIR / "titanG_nh3_reference_manifest.json"
    with open(tmp_manifest, "w") as f:
        json.dump(manifest, f, indent=2, default=str)

    dst = Path(args.dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(tmp_pkl, dst)
    shutil.copy2(tmp_manifest, dst.parent / "titanG_nh3_reference_manifest.json")
    shutil.copy2(jsonl, dst.parent / "titanG_nh3_reference_progress.jsonl")
    print(f"[ref] reference pickle -> {dst}")
    print(f"[ref] manifest         -> {dst.parent / 'titanG_nh3_reference_manifest.json'}")


if __name__ == "__main__":
    main()
