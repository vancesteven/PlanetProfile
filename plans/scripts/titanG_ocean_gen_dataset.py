"""Titan free-gravity MgSO4 / NaCl JOINT (no-ocean + ocean) SBI dataset generation.

Phase B of the Titan free-gravity campaign (plans/fluffy-snacking-fountain.md).
Parameterized analogue of titanG_nh3_gen_dataset.py for the MgSO4 and NaCl joint
configs. Same observable vector {C20, C22, Re_k2, Im_k2}; the difference per
composition is the 2D (Tb x w) structure cache (retry_frozen_as_no_ocean=True)
and the log10_wOcean_ppt salinity range.

These run on the REBUILT caches (float-coerced PfreezeUpper_MPa, Tb=252 K int/float
defect repaired; gate-3 binding PASS 2026-08-07; commit 4ca9a9ed). The gen driver
re-hashes the cache into the manifest so a downstream consumer can confirm it trained
on the repaired bytes, not the pre-fix ones.

Seeds (per campaign spec): MgSO4 data=73, noise=7373; NaCl data=74, noise=7474.
Outputs (in /tmp, copied to Dropbox artifacts by the train step):
  /tmp/titanG_build/datasets/titanG_<comp>_1m.npz  (theta, x, stats)
  /tmp/titanG_build/datasets/titanG_<comp>_gen_manifest.json

libomp hazard: PP generation (numba) and any torch training MUST be separate
processes. Env prefix with thread pins:
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
    python plans/scripts/titanG_ocean_gen_dataset.py --comp NaCl --n 1000000
"""
import argparse, json, os, hashlib, time
import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.sbi_runner import SBIRunner

OUTDIR = "/tmp/titanG_build/datasets"
COMPS = {
    "MgSO4": {
        "cfg": "PlanetProfile/Inference/configs/test54_titan_mgso4_freegrav.json",
        "data_seed": 73, "noise_seed": 7373, "tag": "mgso4",
    },
    "NaCl": {
        "cfg": "PlanetProfile/Inference/configs/test54_titan_nacl_freegrav.json",
        "data_seed": 74, "noise_seed": 7474, "tag": "nacl",
    },
}


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
    ap.add_argument("--comp", required=True, choices=list(COMPS))
    ap.add_argument("--n", type=int, default=1_000_000)
    args = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)

    spec = COMPS[args.comp]
    CFG = spec["cfg"]
    DATA_SEED, NOISE_SEED = spec["data_seed"], spec["noise_seed"]

    obs_names = list(InferenceConfig.from_json(CFG).observables.keys())
    print(f"[gen] Titan {args.comp} joint obs ({len(obs_names)}): {obs_names}")
    assert obs_names == ["C20", "C22", "Re_k2", "Im_k2"], \
        f"unexpected observable set: {obs_names}"

    cfg = json.load(open(CFG)); cfg["mode"] = "sbi"
    cache_path = cfg.get("structure_cache_path")
    cache_hash = _sha256(cache_path)
    print(f"[gen] structure cache: {cache_path} sha256={cache_hash}")
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

    out = os.path.join(OUTDIR, f"titanG_{spec['tag']}_1m.npz")
    np.savez_compressed(
        out, theta=theta.astype(np.float64), x=x.astype(np.float64),
        stats=json.dumps(stats, default=str))

    manifest = {
        "kind": f"titanG_{spec['tag']}_joint_dataset",
        "comp": args.comp,
        "config": CFG, "config_hash": runner.config.generate_hash(),
        "n_requested": args.n, "n_kept": int(len(theta)),
        "data_seed": DATA_SEED, "noise_seed": NOISE_SEED,
        "gen_seconds": dt, "stats": stats,
        "obs_names": obs_names, "n_obs": len(obs_names),
        "npz": out, "npz_sha256": _sha256(out),
        "structure_cache_path": cache_path,
        "structure_cache_sha256": cache_hash,
        "cache_provenance": (
            "REBUILT cache with float-coerced PfreezeUpper_MPa (Tb=252 K int/float "
            "ResetNearestExtrap defect repaired); gate-3 binding PASS 2026-08-07; "
            "commit 4ca9a9ed. The sha256 above must match the installed Test/ cache."
        ),
        "method": ("JOINT no-ocean + ocean posterior. Observable vector "
                   "{C20,C22,Re_k2,Im_k2}; 2D (Tb x w) cache with "
                   "retry_frozen_as_no_ocean=True; log10_wOcean_ppt salinity axis. "
                   "No phase_stability.enforce -> both regimes admitted."),
    }
    with open(os.path.join(OUTDIR, f"titanG_{spec['tag']}_gen_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2, default=str)
    print(f"[gen] {len(obs_names)} obs -> {out}")
    print(f"[gen] manifest -> {OUTDIR}/titanG_{spec['tag']}_gen_manifest.json")


if __name__ == "__main__":
    main()
