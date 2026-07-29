"""v7 open-|Ae| SBI dataset generation (single 21-observable arm).

Machine B, spec plans/active/europa-openae-training-spec.md (user directive
2026-07-25). Generates ONE 1M open-support training set from the v7 config
europa_clipper_v7_openae_11D.json — the induction amp_min support cut is
REMOVED (phase_deg_max 70deg for all three labels; the cache's physical
saturation is now the only support edge). Reuses the v5 2D (Tb x w) cache +
Ae sidecar UNCHANGED (Machine A audit: the cache already spans |Ae| 0.057-0.97
and phase up to 62.5deg), so NO cache rebuild.

Unlike v5 this is a SINGLE arm: no noinduction/nok2 ablation slicing (the v7
deliverable is the one open-interpretation flow; Galileo/Levin knowledge enters
later as a separable inference-time reweight, not a training-time cut).

Seeds (v7, distinct from v5 data=57/noise=5757 so the runs stay cleanly
separable): data=77, noise=7777.
Outputs (in /tmp, copied to Dropbox by the caller):
  /tmp/v7_build/datasets/v7_openae_1m.npz  (theta, x, stats)
  /tmp/v7_build/datasets/v7_gen_manifest.json

Run:
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/v7_gen_dataset.py --n 1000000
"""
import argparse, json, os, hashlib, time
import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.sbi_runner import SBIRunner

CFG = "PlanetProfile/Inference/configs/europa_clipper_v7_openae_11D.json"
OUTDIR = "/tmp/v7_build/datasets"
DATA_SEED, NOISE_SEED = 77, 7777


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
    ap.add_argument("--outdir", default=OUTDIR)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    obs_names = list(InferenceConfig.from_json(CFG).observables.keys())
    print(f"[gen] v7 open-|Ae| obs ({len(obs_names)}): {obs_names}")
    print(f"[gen] induction_bounds: "
          f"{json.dumps(InferenceConfig.from_json(CFG).induction_bounds)}")

    cfg = json.load(open(CFG)); cfg["mode"] = "sbi"
    runner = SBIRunner(InferenceConfig.from_dict(cfg))
    assert list(runner.obs_names) == obs_names, "runner obs order != config order"

    cache_path = cfg.get("structure_cache_path")
    sidecar_path = (str(cache_path) + ".ae_sidecar.pkl") if cache_path else None

    t0 = time.time()
    print(f"[gen] generating {args.n:,} open-support sims "
          f"(seed data={DATA_SEED}, noise={NOISE_SEED}, "
          f"support_guard={runner._support_guard_active()}, drop_nonfinite ON)...")
    theta, x, stats = runner.generate_training_set(
        args.n, seed=DATA_SEED, obs_noise=True, noise_seed=NOISE_SEED)
    dt = time.time() - t0
    print(f"[gen] done in {dt/60:.1f} min: theta={theta.shape}, x={x.shape}")
    print(f"[gen] rejection stats: {json.dumps(stats, default=str)}")
    assert x.shape[1] == len(obs_names), \
        f"x has {x.shape[1]} cols, expected {len(obs_names)}"

    out = os.path.join(args.outdir, "v7_openae_1m.npz")
    np.savez_compressed(
        out, theta=theta.astype(np.float64), x=x.astype(np.float64),
        stats=json.dumps(stats, default=str))

    def sha_npz(path):
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()

    kept = int(len(theta))
    manifest = {
        "kind": "v7_openae_single_arm",
        "config": CFG,
        "config_hash": InferenceConfig.from_json(CFG).generate_hash(),
        "n_requested": args.n, "n_kept": kept,
        "kept_fraction": kept / max(args.n, 1),
        "data_seed": DATA_SEED, "noise_seed": NOISE_SEED,
        "gen_seconds": dt, "stats": stats,
        "obs_names": obs_names,
        "npz": out, "npz_sha256": sha_npz(out),
        "structure_cache_path": cache_path,
        "structure_cache_sha256": _sha256(cache_path),
        "ae_sidecar_sha256": _sha256(sidecar_path),
        "method": ("single 21-obs open-|Ae| training set; amp_min support cut "
                   "removed (phase_deg_max 70deg all labels). Reuses v5 2D cache "
                   "unchanged. No ablation slicing (v7 is one flow). "
                   "Spec: plans/active/europa-openae-training-spec.md."),
    }
    with open(os.path.join(args.outdir, "v7_gen_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2, default=str)
    print(f"[gen] kept {kept:,}/{args.n:,} ({100.0*kept/max(args.n,1):.2f}%) "
          f"-> {out}")
    print(f"[gen] manifest -> {args.outdir}/v7_gen_manifest.json")


if __name__ == "__main__":
    main()
