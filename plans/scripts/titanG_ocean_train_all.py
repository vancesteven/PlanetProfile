"""Train the Titan free-gravity MgSO4 / NaCl JOINT nsf SBI artifact — parameterized.

Phase B of the Titan free-gravity campaign (plans/fluffy-snacking-fountain.md);
authorized by MACHINE-B-HANDOFF §0.16 (flow-training conditional GO, 2026-08-10).
Parameterized analogue of titanG_nh3_train_all.py for the MgSO4 and NaCl joint
configs. Trains ONE nsf flow per composition on (theta, x={C20,C22,Re_k2,Im_k2})
using the DEPLOYED architecture + recipe exactly (§0.16 condition c): nsf via
sbi, deployed early stopping, independent z-scoring. NO architecture
experimentation — the eliminated-capacity verdict (§0.12) holds; MgSO4/NaCl train
on the deployed architecture under split-status (tidal k2 sector quarantined
until its own post-training pushforward re-verifies it, §0.16 Hold 2).

Depends on plans/scripts/titanG_ocean_gen_dataset.py having produced the 1M
dataset in /tmp; asserts the structure-cache sha256 matches the committed gen
manifest so the flow is trained on the repaired (float-coerced) cache bytes.

Seeds (campaign spec): MgSO4 train=73; NaCl train=74.

libomp hazard: this MUST run in a SEPARATE process from PP generation (torch vs
PlanetProfile/numba). Do not import both in one interpreter.

Run (after generation completes; datasets already in /tmp):
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
    python plans/scripts/titanG_ocean_train_all.py --comp NaCl
"""
import argparse, hashlib, json, os, time, shutil
import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.sbi_runner import SBIRunner

DSDIR = "/tmp/titanG_build/datasets"
ARTDIR = "/tmp/titanG_build/artifacts"
DROPBOX_ART = "PlanetProfile/Inference/sbi_artifacts"

COMPS = {
    "MgSO4": {
        "cfg": "PlanetProfile/Inference/configs/test54_titan_mgso4_freegrav.json",
        "train_seed": 73, "tag": "mgso4",
        "cache_sha_expect": "124c8539",
    },
    "NaCl": {
        "cfg": "PlanetProfile/Inference/configs/test54_titan_nacl_freegrav.json",
        "train_seed": 74, "tag": "nacl",
        "cache_sha_expect": "0fdbd44f",
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
    args = ap.parse_args()
    spec = COMPS[args.comp]
    CFG = spec["cfg"]
    TRAIN_SEED = spec["train_seed"]
    tag = spec["tag"]
    ARTNAME = f"titan_freegrav_{tag}_posterior_1m.pt"

    os.makedirs(ARTDIR, exist_ok=True)
    npz = os.path.join(DSDIR, f"titanG_{tag}_1m.npz")
    if not os.path.exists(npz):
        raise SystemExit(f"[train] missing dataset {npz} — run titanG_ocean_gen_dataset.py --comp {args.comp} first")
    with np.load(npz, allow_pickle=True) as d:
        theta = np.asarray(d["theta"], np.float64)
        x = np.asarray(d["x"], np.float64)
        stats = json.loads(str(d["stats"].item())) if "stats" in d.files else None

    cfg = json.load(open(CFG)); cfg["mode"] = "sbi"
    cache_path = cfg.get("structure_cache_path")
    cache_sha = _sha256(cache_path)
    print(f"[train] {args.comp}: structure cache {cache_path} sha256={cache_sha}")
    assert cache_sha and cache_sha.startswith(spec["cache_sha_expect"]), (
        f"cache sha {cache_sha} does not match expected repaired-cache prefix "
        f"{spec['cache_sha_expect']} — refusing to train on unverified bytes")

    # cross-check against the committed gen manifest (train on the SAME bytes gen used)
    gen_manifest = f"validation_reports/titan_freegrav_{tag}_1m/gen_manifest.json"
    if os.path.exists(gen_manifest):
        gm = json.load(open(gen_manifest))
        assert gm.get("structure_cache_sha256") == cache_sha, (
            f"cache sha {cache_sha} != gen_manifest sha {gm.get('structure_cache_sha256')}")
        print(f"[train] cache sha MATCHES committed gen manifest ✓")

    runner = SBIRunner(InferenceConfig.from_dict(cfg))
    assert list(runner.obs_names) == ["C20", "C22", "Re_k2", "Im_k2"], \
        f"unexpected obs order: {runner.obs_names}"
    assert x.shape[1] == len(runner.obs_names), \
        f"x cols {x.shape[1]} != obs {len(runner.obs_names)}"

    print(f"[train] {args.comp}: theta={theta.shape} x={x.shape} seed={TRAIN_SEED} nsf (DEPLOYED arch)")
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
        "kind": f"titanG_{tag}_joint_train", "comp": args.comp,
        "train_seed": TRAIN_SEED, "artifact": dst, "train_seconds": dt,
        "n_train": int(len(theta)), "n_obs": x.shape[1],
        "obs_names": list(runner.obs_names),
        "config": CFG, "config_hash": runner.config.generate_hash(),
        "structure_cache_path": cache_path, "structure_cache_sha256": cache_sha,
        "architecture": "DEPLOYED nsf (no capacity/embedding experimentation — "
                        "§0.12 eliminated-capacity verdict; §0.16 condition c)",
        "split_status": ("tidal k2 sector QUARANTINED by default until this "
                         "composition's own post-training pushforward re-verifies "
                         "it (§0.16 Hold 2 — NH3 verdict NOT ported by assumption)"),
        "dataset_npz": npz,
    }
    with open(os.path.join(ARTDIR, f"titanG_{tag}_train_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"[train] trained in {dt/60:.1f} min -> {dst}")
    print(f"[train] manifest -> {ARTDIR}/titanG_{tag}_train_manifest.json")


if __name__ == "__main__":
    main()
