"""Train the three v5 nsf SBI artifacts from the sliced 1M datasets.

Depends on plans/scripts/v6_gen_dataset.py having produced:
  /tmp/v6_build/datasets/v6_{baseline,noinduction,nok2}_1m.npz

Each arm trains its own nsf flow on its own (theta, sliced-x) via
train_sbi_artifact.py's train path. v5 train seed = 51 (independent of v4's 43).
Artifacts land in /tmp first, then are copied to the Dropbox sbi_artifacts dir.

Run (after generation completes):
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/v6_train_all.py
"""
import hashlib, json, os, time, shutil
import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.sbi_runner import SBIRunner


def _sha256(path):
    if not path or not os.path.exists(path):
        return None
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

DSDIR = "/tmp/v6_build/datasets"
ARTDIR = "/tmp/v6_build/artifacts"
DROPBOX_ART = "PlanetProfile/Inference/sbi_artifacts"
TRAIN_SEED = 61

ARMS = {
    "baseline":    ("PlanetProfile/Inference/configs/europa_clipper_v6_freegrav_11D.json",
                    "europa_clipper_v6_freegrav_11D_posterior_1m.pt"),
    "noinduction": ("PlanetProfile/Inference/configs/europa_clipper_v6_freegrav_noinduction_6obs.json",
                    "europa_clipper_v6_freegrav_noinduction_6obs_posterior_1m.pt"),
    "nok2":        ("PlanetProfile/Inference/configs/europa_clipper_v6_freegrav_nok2_16obs.json",
                    "europa_clipper_v6_freegrav_nok2_16obs_posterior_1m.pt"),
}


def main():
    os.makedirs(ARTDIR, exist_ok=True)
    results = {}
    for tag, (cfgpath, artname) in ARMS.items():
        npz = os.path.join(DSDIR, f"v6_{tag}_1m.npz")
        if not os.path.exists(npz):
            raise SystemExit(f"[train] missing dataset {npz} — run v6_gen_dataset.py first")
        with np.load(npz, allow_pickle=True) as d:
            theta = np.asarray(d["theta"], np.float64)
            x = np.asarray(d["x"], np.float64)
            stats = json.loads(str(d["stats"].item())) if "stats" in d.files else None

        cfg = json.load(open(cfgpath)); cfg["mode"] = "sbi"
        cache_path = cfg.get("structure_cache_path")
        sidecar_path = (str(cache_path) + ".ae_sidecar.pkl") if cache_path else None
        cache_sha = _sha256(cache_path)
        sidecar_sha = _sha256(sidecar_path)
        runner = SBIRunner(InferenceConfig.from_dict(cfg))
        assert x.shape[1] == len(runner.obs_names), \
            f"{tag}: x cols {x.shape[1]} != obs {len(runner.obs_names)}"

        print(f"[train] {tag}: theta={theta.shape} x={x.shape} seed={TRAIN_SEED} nsf")
        t0 = time.time()
        runner.train(theta, x, seed=TRAIN_SEED, density_estimator="nsf")
        if stats is not None:
            runner._train_info["rejection_stats"] = stats
        out = os.path.join(ARTDIR, artname)
        runner.save_artifact(out)
        dt = time.time() - t0
        # copy to Dropbox
        dst = os.path.join(DROPBOX_ART, artname)
        shutil.copy2(out, dst)
        results[tag] = {"artifact": dst, "train_seconds": dt,
                        "n_train": int(len(theta)), "n_obs": x.shape[1],
                        "config": cfgpath,
                        "config_hash": runner.config.generate_hash(),
                        "structure_cache_path": cache_path,
                        "structure_cache_sha256": cache_sha,
                        "ae_sidecar_sha256": sidecar_sha}
        print(f"[train] {tag}: trained in {dt/60:.1f} min -> {dst}")

    with open(os.path.join(ARTDIR, "v6_train_manifest.json"), "w") as f:
        json.dump({"train_seed": TRAIN_SEED, "arms": results}, f, indent=2)
    print(f"[train] all done. manifest -> {ARTDIR}/v6_train_manifest.json")


if __name__ == "__main__":
    main()
