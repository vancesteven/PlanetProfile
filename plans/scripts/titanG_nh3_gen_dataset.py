"""Titan free-gravity NH3 JOINT (no-ocean + ocean) SBI dataset generation.

Phase B of the Titan free-gravity campaign (plans/fluffy-snacking-fountain.md).
Mirrors titanG_gen_dataset.py but points at test54_titan_nh3_freegrav.json — the
JOINT no-ocean + ocean NH3 config (Task #68). Same observable vector
{C20, C22, Re_k2, Im_k2}; the difference is the 2D (Tb x w) NH3 structure cache
(retry_frozen_as_no_ocean=True) + the log10_wOcean_ppt salinity axis.

Seeds (NH3, per config metadata sbi_seeds_2026_08_02): data=72, noise=7272.
Outputs (in /tmp, copied to Dropbox artifacts by the train step):
  /tmp/titanG_build/datasets/titanG_nh3_1m.npz  (theta, x, stats)
  /tmp/titanG_build/datasets/titanG_nh3_gen_manifest.json

Run:
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/titanG_nh3_gen_dataset.py --n 1000000
"""
import argparse, json, os, hashlib, time
import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.sbi_runner import SBIRunner

CFG = "PlanetProfile/Inference/configs/test54_titan_nh3_freegrav.json"
OUTDIR = "/tmp/titanG_build/datasets"
DATA_SEED, NOISE_SEED = 72, 7272


def _sha256(path):
    if not path or not os.path.exists(path):
        return None
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1_000_000)
    args = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)

    obs_names = list(InferenceConfig.from_json(CFG).observables.keys())
    print(f"[gen] Titan NH3 joint obs ({len(obs_names)}): {obs_names}")
    assert obs_names == ["C20", "C22", "Re_k2", "Im_k2"], \
        f"unexpected observable set: {obs_names}"

    cfg = json.load(open(CFG)); cfg["mode"] = "sbi"
    cache_path = cfg.get("structure_cache_path")
    print(f"[gen] structure cache: {cache_path} sha256={_sha256(cache_path)}")
    runner = SBIRunner(InferenceConfig.from_dict(cfg))
    assert list(runner.obs_names) == obs_names, "runner obs order != config order"

    t0 = time.time()
    print(f"[gen] generating {args.n:,} sims (seed data={DATA_SEED}, "
          f"noise={NOISE_SEED}, support_guard ON, drop_nonfinite ON)...")
    theta, x, stats = runner.generate_training_set(
        args.n, seed=DATA_SEED, obs_noise=True, noise_seed=NOISE_SEED)
    dt = time.time() - t0
    print(f"[gen] done in {dt/60:.1f} min: theta={theta.shape}, x={x.shape}")
    print(f"[gen] rejection stats: {json.dumps(stats, default=str)}")
    assert x.shape[1] == len(obs_names), \
        f"x has {x.shape[1]} cols, expected {len(obs_names)}"

    out = os.path.join(OUTDIR, "titanG_nh3_1m.npz")
    np.savez_compressed(
        out, theta=theta.astype(np.float64), x=x.astype(np.float64),
        stats=json.dumps(stats, default=str))

    manifest = {
        "kind": "titanG_nh3_joint_dataset",
        "config": CFG, "config_hash": runner.config.generate_hash(),
        "n_requested": args.n, "n_kept": int(len(theta)),
        "data_seed": DATA_SEED, "noise_seed": NOISE_SEED,
        "gen_seconds": dt, "stats": stats,
        "obs_names": obs_names, "n_obs": len(obs_names),
        "npz": out, "npz_sha256": _sha256(out),
        "structure_cache_path": cache_path,
        "structure_cache_sha256": _sha256(cache_path),
        "method": ("JOINT no-ocean + ocean NH3 posterior (Task #68). Observable "
                   "vector {C20,C22,Re_k2,Im_k2}; 2D (Tb x w) NH3 cache with "
                   "retry_frozen_as_no_ocean=True; log10_wOcean_ppt salinity axis "
                   "[0,1.845]. No phase_stability.enforce -> both regimes admitted."),
    }
    with open(os.path.join(OUTDIR, "titanG_nh3_gen_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2, default=str)
    print(f"[gen] {len(obs_names)} obs -> {out}")
    print(f"[gen] manifest -> {OUTDIR}/titanG_nh3_gen_manifest.json")


if __name__ == "__main__":
    main()
