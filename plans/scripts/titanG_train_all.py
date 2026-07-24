"""Train the Titan free-gravity (no-ocean) nsf SBI artifact — SINGLE-ARM.

Depends on plans/scripts/titanG_gen_dataset.py having produced:
  /tmp/titanG_build/datasets/titanG_noocean_1m.npz

Trains one nsf flow on (theta, x={C20,C22,Re_k2,Im_k2}). Titan no-ocean train
seed = 71 (independent of Europa v6's 61). Artifact lands in /tmp first, then is
copied to the Dropbox sbi_artifacts dir.

libomp hazard: this MUST run in a SEPARATE process from generation (torch vs
PlanetProfile/numba). Do not import both in one interpreter.

Run (after generation completes):
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/titanG_train_all.py
"""
import hashlib, json, os, time, shutil
import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.sbi_runner import SBIRunner

DSDIR = "/tmp/titanG_build/datasets"
ARTDIR = "/tmp/titanG_build/artifacts"
DROPBOX_ART = "PlanetProfile/Inference/sbi_artifacts"
TRAIN_SEED = 71
CFG = "PlanetProfile/Inference/configs/titan_freegrav_noocean.json"
ARTNAME = "titan_freegrav_noocean_posterior_1m.pt"


def _sha256(path):
    if not path or not os.path.exists(path):
        return None
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    os.makedirs(ARTDIR, exist_ok=True)
    npz = os.path.join(DSDIR, "titanG_noocean_1m.npz")
    if not os.path.exists(npz):
        raise SystemExit(f"[train] missing dataset {npz} — run titanG_gen_dataset.py first")
    with np.load(npz, allow_pickle=True) as d:
        theta = np.asarray(d["theta"], np.float64)
        x = np.asarray(d["x"], np.float64)
        stats = json.loads(str(d["stats"].item())) if "stats" in d.files else None

    cfg = json.load(open(CFG)); cfg["mode"] = "sbi"
    cache_path = cfg.get("structure_cache_path")
    cache_sha = _sha256(cache_path)
    runner = SBIRunner(InferenceConfig.from_dict(cfg))
    assert x.shape[1] == len(runner.obs_names), \
        f"x cols {x.shape[1]} != obs {len(runner.obs_names)}"

    print(f"[train] theta={theta.shape} x={x.shape} seed={TRAIN_SEED} nsf")
    t0 = time.time()
    runner.train(theta, x, seed=TRAIN_SEED, density_estimator="nsf")
    if stats is not None:
        runner._train_info["rejection_stats"] = stats
    out = os.path.join(ARTDIR, ARTNAME)
    runner.save_artifact(out)
    dt = time.time() - t0
    dst = os.path.join(DROPBOX_ART, ARTNAME)
    shutil.copy2(out, dst)

    manifest = {
        "kind": "titanG_noocean_train", "train_seed": TRAIN_SEED,
        "artifact": dst, "train_seconds": dt, "n_train": int(len(theta)),
        "n_obs": x.shape[1], "obs_names": list(runner.obs_names),
        "config": CFG, "config_hash": runner.config.generate_hash(),
        "structure_cache_path": cache_path, "structure_cache_sha256": cache_sha,
    }
    with open(os.path.join(ARTDIR, "titanG_train_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"[train] trained in {dt/60:.1f} min -> {dst}")
    print(f"[train] manifest -> {ARTDIR}/titanG_train_manifest.json")


if __name__ == "__main__":
    main()
