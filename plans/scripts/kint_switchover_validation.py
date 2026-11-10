"""K_int switch-over validation (binding review 2026-11-03, STATUS 2026-11-03a).

Required before the vh2008_repaired switch-over commit:
 (V1) repaired rigid H22/C2 x the per-node hydrostatic elastic factor vs the
      deformable branch WITH the figure treatment, over C2 in [0, 1] at
      production-cache nodes 1220 and 1540 (separability error budget);
 (V2) the dense-zb repair ratio at C2 = 0.81 (one w column of the cache).

Also measures, per zb node, the implied Cassini-channel bias of the frozen
convention (shipped rigid x 1.008) against the deformable+figure reference,
and the per-node elastic factor elastic/repaired (hydrostatic), whose zb
spread feeds the libration_sys_frac re-derivation.

Uses the committed TEST cache (88 zb x 40 w); runs anywhere (ALMA y1 per
node, ~ms). No PlanetProfile cache or Mac scratch needed.

  plans/scripts/kint_switchover_validation.py \
      [--cache PATH] [--nodes 1220 1540] [--c2 0 0.25 0.5 0.81 1.0] \
      [--zb-iw 20] [--out-dir validation_reports/enceladus_isostasy/kint_switchover]
"""
import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from PlanetProfile.Gravity import pyALMA3Updated  # noqa: E402
from PlanetProfile.Gravity.Gravity import (  # noqa: E402
    ALMA_TIME_UNIT_S,
    LibrationModelInputs,
)
from PlanetProfile.Gravity.Librations import librations  # noqa: E402
from PlanetProfile.Utilities.defineStructs import Constants  # noqa: E402

CACHE = (REPO / "PlanetProfile/Test/mcmc_results/Enceladus/Cassini_isostasy"
         / "enceladus_seawater_zbw_grid_2d.pkl")
CFG = (REPO / "PlanetProfile/Inference/configs"
       / "enceladus_cassini_isostasy_ocean6D_rhorock.json")
REPAIRED = dict(formulation="vh2008_repaired")


def _node_inputs(structure):
    """Reduced libration stack + ALMA elastic y1 for one cache node.

    Mirrors tests/librations_repair_test.py::_alma_node (kept in sync by
    the cross-check test there)."""
    r = np.asarray(structure["r_m"], float)
    rho = np.asarray(structure["rho"], float)
    mu = np.asarray(structure["mu_Pa"], float)
    eta = np.asarray(structure["eta_Pa_base"], float)
    ph = np.asarray(structure["phases"])
    rheo = np.where(ph == 0, "newton", "elastic").tolist()
    model = pyALMA3Updated.build_model(r, rho, mu, eta, rheo,
                                       np.zeros((r.size, 2)), norm=False)
    om = float(structure["omega"])
    _, _, _, rf = pyALMA3Updated.love_numbers(
        model, [2], [(2 * np.pi / om) / ALMA_TIME_UNIT_S], "tidal",
        "periodic", 4, verbose=False, internal=True)
    M, R = float(structure["Mtot_kg"]), float(structure["R_body_m"])
    y1 = np.real(rf[:, 0]) * M / R / (Constants.G * M / R ** 2)
    rr, rrho, oi, ry = LibrationModelInputs(r, ph, rho, y1)
    kw = dict(omega=om, ecc=float(structure["eccentricity"]), ocean=True,
              ocean_idx=oi)
    return rr, rrho, ry, kw


def _node_set(structure, h22, c2_list):
    rr, rrho, ry, kw = _node_inputs(structure)
    fig = lambda c2: dict(H22_obs_m=h22, compensation_C2=c2)
    out = {
        "shipped_nofig": librations(rr, rrho, **kw),
        "repaired_nofig": librations(rr, rrho, **REPAIRED, **kw),
        "elastic_nofig": librations(rr, rrho, rigid=False, y1=ry, **kw),
        "surface_only": {
            "shipped": librations(rr, rrho, H22_obs_m=h22, **kw),
            "repaired": librations(rr, rrho, H22_obs_m=h22, **REPAIRED, **kw),
            "deformable": librations(rr, rrho, rigid=False, y1=ry,
                                     H22_obs_m=h22, **kw),
        },
    }
    out["elastic_factor"] = out["elastic_nofig"] / out["repaired_nofig"]
    out["elastic_factor_vs_shipped"] = (out["elastic_nofig"]
                                        / out["shipped_nofig"])
    by_c2 = {}
    for c2 in c2_list:
        rig_ship = librations(rr, rrho, **fig(c2), **kw)
        rig_rep = librations(rr, rrho, **fig(c2), **REPAIRED, **kw)
        deform = librations(rr, rrho, rigid=False, y1=ry, **fig(c2), **kw)
        composed = rig_rep * out["elastic_factor"]
        by_c2[f"{c2:g}"] = {
            "shipped_fig": rig_ship,
            "repaired_fig": rig_rep,
            "deformable_fig": deform,
            "composed": composed,
            "separability_err": composed / deform - 1.0,
            "repair_ratio": rig_rep / rig_ship,
            "frozen_bias": rig_ship * 1.008 / deform - 1.0,
        }
    out["by_C2"] = by_c2
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=CACHE)
    ap.add_argument("--nodes", type=int, nargs="+", default=[1220, 1540])
    ap.add_argument("--c2", type=float, nargs="+",
                    default=[0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,
                             0.81, 0.9, 1.0])
    ap.add_argument("--zb-iw", type=int, default=20,
                    help="w column for the dense-zb sweep (V2)")
    ap.add_argument("--zb-c2", type=float, default=0.81)
    ap.add_argument("--out-dir", type=Path,
                    default=REPO / "validation_reports/enceladus_isostasy"
                    / "kint_switchover")
    args = ap.parse_args()

    with args.cache.open("rb") as stream:
        cache = pickle.load(stream)
    zb_grid = np.asarray(cache["zb_km_grid"], float)
    w_grid = np.asarray(cache["wOcean_ppt_grid"], float)
    nw = w_grid.size
    h22 = json.load(CFG.open())["metadata"][
        "libration_model_correction"]["H22_obs_m"]

    result = {"cache": str(args.cache), "H22_obs_m": h22,
              "c2_grid": args.c2, "nodes": {}, "zb_sweep": {}}

    # --- V1: separability at the review's two nodes ----------------------
    for idx in args.nodes:
        st = cache["structures"][idx]
        if st is None:
            print(f"node {idx}: None (skipped)")
            continue
        r = _node_set(st, h22, args.c2)
        r["zb_km"] = float(zb_grid[idx // nw])
        r["w_ppt"] = float(w_grid[idx % nw])
        result["nodes"][str(idx)] = r
        print(f"node {idx} (zb {r['zb_km']:g} km, w {r['w_ppt']:.3g} ppt): "
              f"elastic/repaired {r['elastic_factor']:.6f}")
        for c2, v in r["by_C2"].items():
            print(f"  C2={c2:>4}: sep_err {v['separability_err']:+.3e}  "
                  f"repair {v['repair_ratio']:.5f}  "
                  f"frozen_bias {v['frozen_bias']:+.4%}")

    # --- V2: dense-zb sweep at one w column ------------------------------
    iw = args.zb_iw
    sweep = {}
    for i_zb in range(zb_grid.size):
        idx = i_zb * nw + iw
        st = cache["structures"][idx]
        if st is None:
            continue
        rr, rrho, ry, kw = _node_inputs(st)
        fig = dict(H22_obs_m=h22, compensation_C2=args.zb_c2)
        rig_ship = librations(rr, rrho, **fig, **kw)
        rig_rep = librations(rr, rrho, **fig, **REPAIRED, **kw)
        deform = librations(rr, rrho, rigid=False, y1=ry, **fig, **kw)
        fac = (librations(rr, rrho, rigid=False, y1=ry, **kw)
               / librations(rr, rrho, **REPAIRED, **kw))
        sweep[f"{zb_grid[i_zb]:g}"] = {
            "node": idx,
            "repair_ratio": rig_rep / rig_ship,
            "elastic_factor": fac,
            "separability_err": rig_rep * fac / deform - 1.0,
            "frozen_bias": rig_ship * 1.008 / deform - 1.0,
            "deformable_fig_m": deform,
        }
    result["zb_sweep"] = {"i_w": iw, "w_ppt": float(w_grid[iw]),
                          "C2": args.zb_c2, "by_zb": sweep}
    rows = [(float(k), v) for k, v in sweep.items()]
    if rows:
        rep = np.array([v["repair_ratio"] for _, v in rows])
        fac = np.array([v["elastic_factor"] for _, v in rows])
        sep = np.array([v["separability_err"] for _, v in rows])
        bia = np.array([v["frozen_bias"] for _, v in rows])
        print(f"\nzb sweep (w {w_grid[iw]:.3g} ppt, C2 {args.zb_c2}): "
              f"{len(rows)} nodes, zb {rows[0][0]:g}-{rows[-1][0]:g} km")
        for name, arr in (("repair_ratio", rep), ("elastic_factor", fac),
                          ("separability_err", sep), ("frozen_bias", bia)):
            print(f"  {name:>17}: min {arr.min():+.5f} max {arr.max():+.5f} "
                  f"mean {arr.mean():+.5f}")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out = args.out_dir / "kint_switchover_validation.json"
    out.write_text(json.dumps(result, indent=1))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
