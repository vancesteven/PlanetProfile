"""Titan free-gravity MgSO4 / NaCl JOINT reference MCMC (pocoMC) — multi-seed, B3 protocol.

Phase B of the Titan free-gravity campaign; authorized by MACHINE-B-HANDOFF §0.16
(2026-08-10): "Start NOW in parallel (independent of training): fresh reference
MCMC per composition on {C20,C22,Re_k2,Im_k2}, B3 protocol (n_eff~2000 class,
>=3 seeds, pinned env, pooled + renormalized) — these are the crosscheck targets
and the quarantine-re-verification baseline."

Parameterized analogue of titanG_nh3_reference_mcmc.py + the B3 multi-seed
protocol (plans/scripts/b3_reference_wander.py). Runs the full {C20,C22,Re_k2,
Im_k2} reference posterior at n_effective=2000 / n_active=1024 (the ratified
regime-preserving mapping: ratio ~2 like the originals, 4x ensemble) for >=3
fresh seeds per composition, then POOLS the equal-resolution equal-weight seed
samples into a single reference (renormalized = equal per-seed contribution).

Determinism / provenance
------------------------
- Each seed pins runner.random_state (config random_state deliberately ignored so
  seeds sweep, exactly as B3).
- Forward model reads the pinned repaired 2D structure cache; sha256 recorded and
  cross-checked against the committed gen manifest.
- Per-seed pickles + a pooled pickle + a scatter/pool manifest go to /tmp first,
  then copied to validation_reports/titan_freegrav_<comp>_1m/reference/ (NO Test/ write).

Seeds (campaign-traceable): MgSO4 {73,173,273}; NaCl {74,174,274}.

Run (separate process; env pinned):
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
    python plans/scripts/titanG_ocean_reference_mcmc.py --comp NaCl
"""
import argparse, hashlib, json, os, shutil, subprocess, time
from pathlib import Path

import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.mcmc_runner import MCMCRunner

ROOT = Path(__file__).resolve().parents[2]
COMPS = {
    "MgSO4": {
        "cfg": ROOT / "PlanetProfile/Inference/configs/test54_titan_mgso4_freegrav.json",
        "tag": "mgso4", "seeds": [73, 173, 273], "cache_sha_expect": "124c8539",
    },
    "NaCl": {
        "cfg": ROOT / "PlanetProfile/Inference/configs/test54_titan_nacl_freegrav.json",
        "tag": "nacl", "seeds": [74, 174, 274], "cache_sha_expect": "0fdbd44f",
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


def _git_sha():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip()
    except Exception:
        return None


def _run_one(cfg, seed, n_effective, n_active, tmp_dir, tag):
    runner = MCMCRunner(cfg)
    runner.random_state = int(seed)
    runner.n_effective = int(n_effective)
    if n_active is not None:
        runner.n_active = int(n_active)
    jsonl = tmp_dir / f"titanG_{tag}_ref_seed{seed}_progress.jsonl"
    print(f"[ref]   seed={seed} n_effective={n_effective} n_active={n_active} ...")
    t0 = time.time()
    res = runner.run(progress_jsonl_path=str(jsonl))
    dt = time.time() - t0
    samples = np.asarray(res.samples)
    md = getattr(res, "metadata", {}) or {}
    print(f"[ref]   seed={seed} done in {dt/60:.1f} min: samples={samples.shape} log_Z={md.get('log_Z')}")
    return res, samples, md, dt, list(res.param_names), jsonl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--comp", required=True, choices=list(COMPS))
    ap.add_argument("--seeds", type=int, nargs="+", default=None,
                    help="override the campaign seeds")
    ap.add_argument("--n-effective", type=int, default=2000)
    ap.add_argument("--n-active", type=int, default=1024,
                    help="raise in step with n_effective to preserve regime (ratio ~2)")
    args = ap.parse_args()

    spec = COMPS[args.comp]
    tag = spec["tag"]
    seeds = args.seeds if args.seeds else spec["seeds"]
    cfg_path = spec["cfg"]
    tmp_dir = Path(f"/tmp/titanG_build/{tag}_reference")
    tmp_dir.mkdir(parents=True, exist_ok=True)
    dst_dir = ROOT / f"validation_reports/titan_freegrav_{tag}_1m/reference"
    dst_dir.mkdir(parents=True, exist_ok=True)

    cfg = InferenceConfig.from_json(str(cfg_path))
    print(f"[ref] {args.comp}: config={cfg_path}")
    print(f"[ref] observables ({len(cfg.observables)}): {list(cfg.observables)}")
    print(f"[ref] param_space: {list(cfg.param_space.keys())}")
    assert list(cfg.observables) == ["C20", "C22", "Re_k2", "Im_k2"], \
        f"unexpected obs set: {list(cfg.observables)}"

    cache_path = cfg.structure_cache_path
    cache_sha = _sha256(cache_path)
    print(f"[ref] structure cache sha256 = {cache_sha}")
    assert cache_sha and cache_sha.startswith(spec["cache_sha_expect"]), (
        f"cache sha {cache_sha} != expected repaired prefix {spec['cache_sha_expect']}")
    gen_manifest = ROOT / f"validation_reports/titan_freegrav_{tag}_1m/gen_manifest.json"
    if gen_manifest.exists():
        gm = json.load(open(gen_manifest))
        assert gm.get("structure_cache_sha256") == cache_sha, \
            f"cache sha != gen manifest sha ({gm.get('structure_cache_sha256')})"
        print("[ref] cache sha MATCHES committed gen manifest ✓")

    runs = []
    param_names = None
    for sd in seeds:
        res, samples, md, dt, pnames, jsonl = _run_one(
            cfg, sd, args.n_effective, args.n_active, tmp_dir, tag)
        param_names = pnames
        tmp_pkl = tmp_dir / f"titan_freegrav_{tag}_reference_seed{sd}.pkl"
        res.save(str(tmp_pkl))
        shutil.copy2(tmp_pkl, dst_dir / tmp_pkl.name)
        shutil.copy2(jsonl, dst_dir / jsonl.name)
        seed_w = getattr(res, "weights", None)
        seed_w = None if seed_w is None else np.asarray(seed_w, dtype=float)
        runs.append({"seed": int(sd), "n_samples": int(samples.shape[0]),
                     "log_Z": md.get("log_Z"), "log_Z_err": md.get("log_Z_err"),
                     "elapsed_s": dt, "samples": samples, "weights": seed_w,
                     "pkl": tmp_pkl.name})

    # POOL: concatenate the pocoMC importance-weighted seed samples. The per-seed
    # weights are NOT uniform (pocoMC importance weights), so both samples AND
    # weights must be concatenated. Each seed's weights are renormalized to sum
    # to 1/n_seeds so every seed contributes EQUAL total posterior mass to the
    # pool (equal per-seed contribution, renormalized). The pool is the multi-seed
    # reference / crosscheck target.
    n_seeds = len(runs)
    pooled = np.concatenate([r["samples"] for r in runs], axis=0)
    if all(r["weights"] is not None for r in runs):
        pooled_w = np.concatenate(
            [r["weights"] / r["weights"].sum() / n_seeds for r in runs], axis=0)
    else:
        # no per-seed weights => equal-weight draws; assign uniform pooled weights
        pooled_w = np.full(pooled.shape[0], 1.0 / pooled.shape[0], dtype=float)
    assert pooled_w.shape[0] == pooled.shape[0], \
        f"pooled weights {pooled_w.shape[0]} != pooled samples {pooled.shape[0]}"
    pooled_pkl = tmp_dir / f"titan_freegrav_{tag}_reference_pooled.pkl"
    # save pooled via the last result object's container for schema compatibility
    res.samples = pooled
    res.weights = pooled_w
    res.save(str(pooled_pkl))
    shutil.copy2(pooled_pkl, dst_dir / pooled_pkl.name)

    # per-parameter between-seed scatter (the reference-noise floor for crosscheck tol)
    per_param = {}
    for name in param_names:
        j = param_names.index(name)
        seed_means = np.array([float(np.mean(r["samples"][:, j])) for r in runs])
        seed_stds = np.array([float(np.std(r["samples"][:, j])) for r in runs])
        per_param[name] = {
            "per_seed_mean": {int(r["seed"]): float(np.mean(r["samples"][:, j])) for r in runs},
            "per_seed_std": {int(r["seed"]): float(np.std(r["samples"][:, j])) for r in runs},
            "between_seed_mean_range": float(seed_means.max() - seed_means.min()),
            "between_seed_mean_std": float(seed_means.std(ddof=1)) if seed_means.size > 1 else None,
            "pooled_mean": float(np.mean(pooled[:, j])),
            "pooled_std": float(np.std(pooled[:, j])),
        }

    manifest = {
        "kind": f"titanG_{tag}_joint_reference_mcmc_multiseed_B3",
        "comp": args.comp,
        "config": os.path.relpath(str(cfg_path), ROOT),
        "config_hash": cfg.generate_hash(),
        "seeds": [int(s) for s in seeds],
        "n_effective": args.n_effective, "n_active": args.n_active,
        "protocol": ("B3 (§0.16): n_eff=2000 / n_active=1024 regime-preserving "
                     "(ratio ~2), >=3 fresh seeds, pooled equal-weight samples."),
        "obs_names": list(cfg.observables),
        "param_names": param_names,
        "structure_cache_path": str(cache_path),
        "structure_cache_sha256": cache_sha,
        "git_sha": _git_sha(),
        "runs": [{k: v for k, v in r.items() if k != "samples"} for r in runs],
        "n_pooled_samples": int(pooled.shape[0]),
        "pooled_pkl": pooled_pkl.name,
        "per_param_scatter": per_param,
        "note": ("Titan free-gravity {C} JOINT no-ocean+ocean multi-seed reference. "
                 "Crosscheck target + tidal-quarantine re-verification baseline "
                 "(§0.16). Trained on the repaired float-coerced cache bytes.").format(C=args.comp),
    }
    tmp_manifest = tmp_dir / f"titanG_{tag}_reference_manifest.json"
    with open(tmp_manifest, "w") as f:
        json.dump(manifest, f, indent=2, default=str)
    shutil.copy2(tmp_manifest, dst_dir / tmp_manifest.name)
    print(f"[ref] pooled reference ({pooled.shape[0]} samples) -> {dst_dir / pooled_pkl.name}")
    print(f"[ref] manifest -> {dst_dir / tmp_manifest.name}")
    for name in ("D_ocean_km", "log10_wOcean_ppt", "Tb_K"):
        if name in per_param:
            e = per_param[name]
            print(f"[ref]   {name:18s} pooled_mean={e['pooled_mean']:.4g} "
                  f"between_seed_std={e['between_seed_mean_std']}")


if __name__ == "__main__":
    main()
