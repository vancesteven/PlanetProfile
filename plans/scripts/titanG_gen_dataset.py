"""Titan free-gravity (no-ocean) SBI dataset generation — SINGLE-ARM.

Unlike the Europa v6 campaign, the Titan free-gravity campaign has NO ablation
siblings: the observable vector is exactly {C20, C22, Re_k2, Im_k2} (induction +
h2 dropped — Saturn's field is spin-aligned and the ionosphere screens the
inducing signal; there is no measured Titan h2). So this generates ONE 1M
training set from the no-ocean config, with no name-based column slicing.

Config: titan_freegrav_noocean.json (test52 10D no-ocean interior + free-gravity
geodesy: CMR2 dropped, C20/C22 = Petricca 2025, agnostic dC20_nh/dC22_nh offsets;
12 sampled params). Reuses the EXISTING test52 no-ocean 1D structure cache — no
ocean EOS, no salinity axis, no rebuild.

Seeds (Titan no-ocean, independent of Europa v6's 67/6767): data=71, noise=7171.
Outputs (in /tmp, copied to Dropbox by the caller if desired):
  /tmp/titanG_build/datasets/titanG_noocean_1m.npz  (theta, x, stats)
  /tmp/titanG_build/datasets/titanG_noocean_gen_manifest.json

Run:
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/titanG_gen_dataset.py --n 1000000
"""
import argparse, json, os, hashlib, time
import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.sbi_runner import SBIRunner

CFG = "PlanetProfile/Inference/configs/titan_freegrav_noocean.json"
OUTDIR = "/tmp/titanG_build/datasets"
DATA_SEED, NOISE_SEED = 71, 7171


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
    print(f"[gen] Titan no-ocean obs ({len(obs_names)}): {obs_names}")
    assert obs_names == ["C20", "C22", "Re_k2", "Im_k2"], \
        f"unexpected observable set: {obs_names}"

    cfg = json.load(open(CFG)); cfg["mode"] = "sbi"
    cache_path = cfg.get("structure_cache_path")
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

    out = os.path.join(OUTDIR, "titanG_noocean_1m.npz")
    np.savez_compressed(
        out, theta=theta.astype(np.float64), x=x.astype(np.float64),
        stats=json.dumps(stats, default=str))

    manifest = {
        "kind": "titanG_noocean_dataset",
        "config": CFG, "config_hash": runner.config.generate_hash(),
        "n_requested": args.n, "n_kept": int(len(theta)),
        "data_seed": DATA_SEED, "noise_seed": NOISE_SEED,
        "gen_seconds": dt, "stats": stats,
        "obs_names": obs_names, "n_obs": len(obs_names),
        "npz": out, "npz_sha256": _sha256(out),
        "structure_cache_path": cache_path,
        "structure_cache_sha256": _sha256(cache_path),
        "method": ("SINGLE-ARM: {C20,C22,Re_k2,Im_k2} only. No ablation slicing "
                   "(induction + h2 dropped for Titan). Reuses the test52 "
                   "no-ocean 1D structure cache unchanged."),
    }
    with open(os.path.join(OUTDIR, "titanG_noocean_gen_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2, default=str)
    print(f"[gen] {len(obs_names)} obs -> {out}")
    print(f"[gen] manifest -> {OUTDIR}/titanG_noocean_gen_manifest.json")


if __name__ == "__main__":
    main()
