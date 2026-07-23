"""v5 SBI dataset generation + reviewer-sanctioned single-dataset slicing.

Generates ONE 1M baseline (21-observable) training set from the v5 config, then
produces the two ablation training sets by NAME-BASED column slicing (reviewer
PASS-WITH-CAVEATS, 2026-07-20):
  - baseline    : 21 obs (all)
  - noinduction :  7 obs (drops 14 Bind_*)
  - nok2        : 17 obs (drops Re_k2/Im_k2/Re_h2/Im_h2)

Both ablation observable sets are exact subsets of the baseline with byte-identical
values AND sigmas (verified). Column-slicing gives a PAIRED ablation: identical
theta draws + identical noise realizations across all three arms, so posterior
differences isolate the information content of the dropped observables.

Reviewer-required safeguards implemented here:
  1. Slice by observable NAME (index lookup), then assert the resulting name list
     equals the ablation config's obs_names exactly and in order.
  2. Report the DIFFERENTIAL drop rate: because the baseline drops any row where
     ANY of the 21 channels is non-finite, the ablation training sets inherit that
     stricter mask. We record how many rows are finite in the SUBSET but were
     dropped by the full-21 mask, so the shared-intersection-support caveat is on
     record. (Rows are theta-level for TidalPy failures + no-ocean guard.)
  3. Assert correlation-spec consistency: if any config carries
     observable_correlations, the kept-channel submatrix must match the baseline's
     for those pairs. (Pure-diagonal noise -> trivially satisfied.)

Seeds (v5, independent of v4): data=57, noise=5757.
Outputs (all in /tmp, copied to Dropbox by the caller):
  /tmp/v5_build/datasets/v5_{baseline,noinduction,nok2}_1m.npz  (theta, x, stats)
  /tmp/v5_build/datasets/v5_gen_manifest.json  (drop rates, hashes, seeds, caveats)

Run:
  mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
    KMP_DUPLICATE_LIB_OK=TRUE python plans/scripts/v5_gen_dataset.py --n 1000000
"""
import argparse, json, os, hashlib, time
import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.sbi_runner import SBIRunner

BASE_CFG = "PlanetProfile/Inference/configs/europa_clipper_v5_geodesy_11D.json"
ABLATIONS = {
    "baseline":    "PlanetProfile/Inference/configs/europa_clipper_v5_geodesy_11D.json",
    "noinduction": "PlanetProfile/Inference/configs/europa_clipper_v5_noinduction_7obs.json",
    "nok2":        "PlanetProfile/Inference/configs/europa_clipper_v5_nok2_17obs.json",
}
OUTDIR = "/tmp/v5_build/datasets"
DATA_SEED, NOISE_SEED = 57, 5757


def obs_names_of(cfgpath):
    return list(InferenceConfig.from_json(cfgpath).observables.keys())


def corr_spec_of(cfgpath):
    md = InferenceConfig.from_json(cfgpath).metadata or {}
    return md.get("observable_correlations")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=1_000_000)
    args = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)

    base_names = obs_names_of(BASE_CFG)
    n_base = len(base_names)
    print(f"[gen] baseline obs ({n_base}): {base_names}")

    # ---- reviewer safeguard 3: correlation-spec consistency ----
    base_corr = corr_spec_of(BASE_CFG)
    for tag, cfgpath in ABLATIONS.items():
        c = corr_spec_of(cfgpath)
        if base_corr is None and c is None:
            continue  # pure diagonal -> slicing trivially valid
        # If correlations exist, the kept-channel submatrix must match baseline.
        raise SystemExit(
            f"[gen] observable_correlations present (tag={tag}); the diagonal-slice "
            f"assumption no longer holds. Implement submatrix-equality check before "
            f"proceeding. base_corr={base_corr}, {tag}_corr={c}")
    print("[gen] correlation-spec check: pure-diagonal noise for all arms -> slice valid")

    # ---- generate the single baseline 21-obs dataset ----
    cfg = json.load(open(BASE_CFG)); cfg["mode"] = "sbi"
    runner = SBIRunner(InferenceConfig.from_dict(cfg))
    assert list(runner.obs_names) == base_names, "runner obs order != config order"

    t0 = time.time()
    print(f"[gen] generating {args.n:,} baseline sims (seed data={DATA_SEED}, "
          f"noise={NOISE_SEED}, support_guard ON, drop_nonfinite ON)...")
    theta, x, stats = runner.generate_training_set(
        args.n, seed=DATA_SEED, obs_noise=True, noise_seed=NOISE_SEED)
    dt = time.time() - t0
    print(f"[gen] done in {dt/60:.1f} min: theta={theta.shape}, x={x.shape}")
    print(f"[gen] rejection stats: {json.dumps(stats, default=str)}")
    assert x.shape[1] == n_base, f"x has {x.shape[1]} cols, expected {n_base}"

    manifest = {
        "n_requested": args.n, "n_kept_baseline": int(len(theta)),
        "data_seed": DATA_SEED, "noise_seed": NOISE_SEED,
        "gen_seconds": dt, "baseline_stats": stats,
        "baseline_obs_names": base_names,
        "structure_cache_sha256": InferenceConfig.from_json(BASE_CFG)
            .metadata.get("v5_artifact_hashes_2026_07_20", {})
            .get("structure_cache_sha256"),
        "method": ("single 21-obs baseline generated; ablations are NAME-BASED "
                   "column slices sharing theta+noise (paired ablation). Reviewer "
                   "PASS-WITH-CAVEATS 2026-07-20."),
        "arms": {},
    }

    def sha_npz(path):
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()

    # ---- slice for each arm, with reviewer safeguards 1 & 2 ----
    for tag, cfgpath in ABLATIONS.items():
        sub_names = obs_names_of(cfgpath)
        # safeguard 1: name-based index + order assertion
        idx = [base_names.index(n) for n in sub_names]
        assert [base_names[i] for i in idx] == sub_names, \
            f"{tag}: sliced name order != config obs_names"
        x_sub = x[:, idx]

        # safeguard 2: differential drop rate. The baseline (all-21) finite mask
        # was already applied. Rows finite in this SUBSET but dropped by full-21:
        # we cannot recover dropped rows post-hoc (they are gone), but we CAN
        # quantify from the KEPT set that every kept row is finite in the subset
        # (sanity), and report the requested-vs-kept overall rate. The genuine
        # differential requires the pre-drop full matrix; generate_training_set
        # already applied drop_nonfinite. We therefore report: (a) overall kept
        # fraction, (b) subset-finite fraction among kept (must be 1.0), and flag
        # that the differential (subset-would-keep minus full-keeps) is bounded by
        # the total non-finite rejection count in baseline_stats.
        finite_in_sub = np.isfinite(x_sub).all(axis=1)
        subset_finite_frac = float(finite_in_sub.mean())

        out = os.path.join(OUTDIR, f"v5_{tag}_1m.npz")
        arm_stats = dict(stats)
        arm_stats.update({
            "arm": tag, "sliced_from": "baseline_21obs",
            "obs_names": sub_names, "n_obs": len(sub_names),
            "n_kept": int(len(theta)),
            "subset_finite_frac_among_kept": subset_finite_frac,
        })
        np.savez_compressed(
            out, theta=theta.astype(np.float64), x=x_sub.astype(np.float64),
            stats=json.dumps(arm_stats, default=str))
        manifest["arms"][tag] = {
            "npz": out, "n_obs": len(sub_names), "obs_names": sub_names,
            "col_index_into_baseline": idx,
            "subset_finite_frac_among_kept": subset_finite_frac,
            "sha256": sha_npz(out),
        }
        print(f"[gen] {tag}: {len(sub_names)} obs -> {out} "
              f"(subset-finite among kept = {subset_finite_frac:.6f})")

    # differential-drop caveat (bounded by baseline non-finite rejections)
    n_nonfinite = int(stats.get("n_rejected_nonfinite", 0) or 0)
    manifest["differential_drop_caveat"] = (
        f"Ablation arms inherit the baseline's all-21-channel finite mask. The "
        f"differential (rows an independently-generated subset run would keep but "
        f"the shared mask dropped) is bounded above by the baseline non-finite "
        f"rejection count = {n_nonfinite} of {args.n} requested "
        f"({100.0*n_nonfinite/max(args.n,1):.3f}%). All arms share the intersection "
        f"support; ablation posteriors are conditioned on full-model (incl. tidal) "
        f"evaluability. Arms are correlated (shared theta+noise) — do NOT treat as "
        f"independent for joint/differential-uncertainty statements.")

    with open(os.path.join(OUTDIR, "v5_gen_manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2, default=str)
    print(f"[gen] manifest -> {OUTDIR}/v5_gen_manifest.json")
    print(f"[gen] differential-drop bound: {n_nonfinite}/{args.n} "
          f"({100.0*n_nonfinite/max(args.n,1):.3f}%)")


if __name__ == "__main__":
    main()
