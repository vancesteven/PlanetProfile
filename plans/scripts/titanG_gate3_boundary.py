"""Gate-3: freeze-line half-cell-shift -> posterior-ocean-fraction sensitivity.

Reviewer launch-condition (2) for the MgSO4/NaCl Titan joint production caches
(plans/fluffy-snacking-fountain.md validation item 3; MACHINE-B-HANDOFF §0.10):

  "Boundary-placement error - quantify how a half-cell shift of the freeze line
   changes the posterior ocean fraction; gate on it; refine Tb if needed."

The ocean/no-ocean freeze boundary in the 2D (Tb x log10 w) structure cache is
located only to the Tb grid spacing. The trained flow interpolates over
(Tb, log10 w), so a sample's ocean-vs-no-ocean classification flips as the
freeze line moves within a cell. This gate measures how much the prior-mass
ocean fraction moves under a rigid +-half-cell Tb translation of the freeze
line relative to the grid -- the boundary-placement uncertainty the fixed grid
imposes on the posterior ocean fraction.

PURE cache analysis: reads the .pkl + config only, runs NO PlanetProfile
forward model, imports no numba/torch. Safe to run standalone.

Method
------
- Reshape the flat `structures` list (len n_Tb * n_w, index = i_Tb*n_w + i_w)
  into an (n_Tb, n_w) classification grid: OCEAN (has_ocean True), FROZEN
  (has_ocean False, a real no-ocean interior), or NONE (unbuildable -> not
  support, excluded from numerator AND denominator, exactly as grid_interp_2d
  rejects all-None regions).
- Put a uniform fine evaluation grid over the config prior box
  (Tb in param_space.Tb_K.bounds, log10 w in param_space.log10_wOcean_ppt.bounds).
  A uniform fine grid in (Tb, log10 w) already encodes the uniform-Tb x
  uniform-log10w (Jeffreys) prior.
- Classify each fine point by its NEAREST cache node in (Tb, log10 w) space
  (log10 w because the interpolant and salinity prior live in log10 w). Ocean
  fraction = (# fine points nearest an OCEAN node) / (# nearest a BUILDABLE
  node), i.e. NONE-nearest points are dropped from support.
- Half-cell shift: translate ALL fine points by +-0.5 * median(dTb) along Tb
  before the nearest-node lookup. This is a rigid half-cell translation of the
  freeze line relative to the grid. Report:
    f0 = ocean_frac(0), fcold = ocean_frac(-0.5 dTb), fwarm = ocean_frac(+0.5 dTb)
    span      = |fcold - fwarm|          (full half-cell placement band)
    half_span = span / 2
    max_dev   = max(|fcold-f0|, |fwarm-f0|)
- MgSO4 corner sub-analysis (--corner-w-min): repeat restricted to fine points
  with w >= corner_w_min ppt (the w>=150 phantom-ocean-leak corner the reviewer
  flagged for MgSO4).

Acceptance
----------
PROVISIONAL preregistered bar (`--threshold`, default 0.03): the half-cell
boundary-placement band `span` should be a small fraction of prior mass so the
ocean/no-ocean split is not grid-artifact-dominated. This is NOT a tuned gate;
the scientific-reviewer adjudicates the raw numbers and sets the binding bar.
Reported PASS/FAIL is advisory pending that sign-off.

Usage
-----
  micromamba run -n PPcl python plans/scripts/titanG_gate3_boundary.py \
     --cache PlanetProfile/Test/mcmc_results/Titan/Test54_nacl_ocean/titan_nacl_joint_structure_grid_2d.pkl \
     --config PlanetProfile/Inference/configs/test54_titan_nacl_freegrav.json \
     --out validation_reports/titan_freegrav_nacl_1m/gate3_boundary.json

  micromamba run -n PPcl python plans/scripts/titanG_gate3_boundary.py \
     --cache PlanetProfile/Test/mcmc_results/Titan/Test54_mgso4_ocean/titan_mgso4_joint_structure_grid_2d.pkl \
     --config PlanetProfile/Inference/configs/test54_titan_mgso4_freegrav.json \
     --corner-w-min 150 \
     --out validation_reports/titan_freegrav_mgso4_1m/gate3_boundary.json
"""
from __future__ import annotations

import argparse
import json
import pickle
from pathlib import Path

import numpy as np

OCEAN, FROZEN, NONE = 1, 0, -1


def _load_grid(cache_path):
    cache = pickle.load(open(cache_path, "rb"))
    Tb = np.asarray(cache["Tb_K_grid"], float)
    w = np.asarray(cache["wOcean_ppt_grid"], float)
    structs = cache["structures"]
    n_Tb, n_w = len(Tb), len(w)
    assert len(structs) == n_Tb * n_w, \
        f"structures len {len(structs)} != n_Tb*n_w {n_Tb*n_w}"
    cls = np.full((n_Tb, n_w), NONE, dtype=int)
    for i in range(n_Tb):
        for j in range(n_w):
            s = structs[i * n_w + j]
            if s is None:
                cls[i, j] = NONE
            elif s.get("has_ocean"):
                cls[i, j] = OCEAN
            else:
                cls[i, j] = FROZEN
    return Tb, w, cls, cache.get("ocean_comp")


def _onset_local_dTb(Tb_nodes, cls):
    """Local Tb spacing straddling the FROZEN/NONE -> OCEAN onset, per w-column.

    Reviewer ab2d3748d0b00c78c (2026-08-07): the BINDING gate-3 quantity is the
    absolute onset-temperature placement uncertainty = half the LOCAL Tb cell
    spanning the freeze transition, NOT the global median(dTb) (which understates
    the half-cell if the onset happens to sit in a coarse band) and NOT the
    ocean-fraction span (which penalizes a narrow prior box and inverts the true
    ranking). Returns (max_local_dTb, per_column_list). For each column with at
    least one OCEAN node and one non-OCEAN node below it, the onset cell is the
    (Tb[i], Tb[i+1]) pair where cls flips non-ocean -> ocean with increasing Tb;
    local dTb = Tb[i+1]-Tb[i]. half_cell_K = max over columns / 2 (worst-case
    placement uncertainty across the salinity axis).
    """
    n_T, n_w = cls.shape
    per_col = []
    for j in range(n_w):
        col = cls[:, j]
        ocean_iTs = np.where(col == OCEAN)[0]
        if ocean_iTs.size == 0:
            continue
        first_ocean = int(ocean_iTs.min())
        if first_ocean == 0:
            continue  # onset at/below grid floor; no straddling cell in-box
        # onset cell straddles the first ocean node and the node just below it
        local = float(Tb_nodes[first_ocean] - Tb_nodes[first_ocean - 1])
        per_col.append({"iw": j, "onset_Tb_K": float(Tb_nodes[first_ocean]),
                        "local_dTb_K": local})
    max_local = max((c["local_dTb_K"] for c in per_col), default=float("nan"))
    return max_local, per_col


def _ocean_frac(Tb_nodes, logw_nodes, cls, tb_eval, logw_eval, shift_tb,
                corner_logw_min=None):
    """Nearest-node ocean fraction over the fine (tb_eval x logw_eval) grid,
    with the freeze line rigidly translated by shift_tb along Tb.

    Scale Tb and log10 w to comparable units before the nearest-node lookup so
    neither axis dominates the Voronoi assignment (Tb spans ~10-40 K, log10 w
    spans ~2.3 decades).
    """
    tb_span = Tb_nodes.max() - Tb_nodes.min()
    lw_span = logw_nodes.max() - logw_nodes.min()
    tb_scale = tb_span if tb_span > 0 else 1.0
    lw_scale = lw_span if lw_span > 0 else 1.0

    TB, LW = np.meshgrid(tb_eval, logw_eval, indexing="ij")
    TB = TB.ravel()
    LW = LW.ravel()
    if corner_logw_min is not None:
        keep = LW >= corner_logw_min
        TB, LW = TB[keep], LW[keep]
        if TB.size == 0:
            return float("nan"), 0

    # translate eval points by shift_tb (rigid freeze-line shift rel. to grid)
    tb_q = (TB + shift_tb) / tb_scale
    lw_q = LW / lw_scale
    tbn = Tb_nodes / tb_scale
    lwn = logw_nodes / lw_scale

    # nearest node index in the (n_Tb, n_w) grid, separably (grids are axis-aligned)
    i_near = np.abs(tb_q[:, None] - tbn[None, :]).argmin(axis=1)
    j_near = np.abs(lw_q[:, None] - lwn[None, :]).argmin(axis=1)
    node_cls = cls[i_near, j_near]

    buildable = node_cls != NONE
    n_build = int(buildable.sum())
    if n_build == 0:
        return float("nan"), 0
    n_ocean = int((node_cls == OCEAN).sum())
    return n_ocean / n_build, n_build


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--corner-w-min", type=float, default=None,
                    help="if set, also report ocean fraction for w >= this ppt")
    ap.add_argument("--n-tb", type=int, default=4000)
    ap.add_argument("--n-w", type=int, default=800)
    ap.add_argument("--threshold", type=float, default=0.03,
                    help="secondary/diagnostic ocean-fraction span bar (reviewer adjudicates)")
    ap.add_argument("--half-cell-bar-K", type=float, default=0.5,
                    help="BINDING absolute onset-placement bar: max local half-cell in K "
                         "(reviewer ab2d3748d0b00c78c 2026-08-07: <= 0.5 K)")
    args = ap.parse_args()

    Tb, w, cls, comp = _load_grid(args.cache)
    cfg = json.load(open(args.config))
    ps = cfg["param_space"]
    tb_lo, tb_hi = ps["Tb_K"]["bounds"]
    lw_lo, lw_hi = ps["log10_wOcean_ppt"]["bounds"]

    logw_nodes = np.log10(w)
    dTb_med = float(np.median(np.diff(Tb)))
    # BINDING metric: local onset half-cell in absolute Kelvin (reviewer).
    max_local_dTb, onset_cols = _onset_local_dTb(Tb, cls)
    half_cell_K = 0.5 * max_local_dTb
    # SECONDARY diagnostic: the ocean-fraction span uses the global median dTb
    # (unchanged, so the two reports stay comparable to the pre-fix numbers).
    half = 0.5 * dTb_med

    tb_eval = np.linspace(tb_lo, tb_hi, args.n_tb)
    logw_eval = np.linspace(lw_lo, lw_hi, args.n_w)

    def frac(shift, corner=None):
        return _ocean_frac(Tb, logw_nodes, cls, tb_eval, logw_eval, shift,
                           corner_logw_min=corner)

    f0, n_build = frac(0.0)
    fcold, _ = frac(-half)
    fwarm, _ = frac(+half)
    span = abs(fcold - fwarm)
    max_dev = max(abs(fcold - f0), abs(fwarm - f0))

    # sandwich check on the classification grid (no non-ocean node between two
    # ocean nodes in a w-column) — mirrors the build-driver invariant.
    sandwich_violations = []
    for j in range(len(w)):
        col = cls[:, j]
        ocean_iTs = np.where(col == OCEAN)[0]
        if ocean_iTs.size >= 2:
            for it in range(int(ocean_iTs.min()) + 1, int(ocean_iTs.max())):
                if col[it] != OCEAN:
                    sandwich_violations.append(
                        [float(Tb[it]), float(w[j]), int(col[it])])

    binding_pass = (half_cell_K <= args.half_cell_bar_K
                    and len(sandwich_violations) == 0)

    result = {
        "kind": "titanG_gate3_freeze_line_placement",
        "comp": comp,
        "cache": args.cache,
        "config": args.config,
        "grid": {
            "n_Tb": int(len(Tb)), "n_w": int(len(w)),
            "Tb_range_K": [float(Tb.min()), float(Tb.max())],
            "w_range_ppt": [float(w.min()), float(w.max())],
            "dTb_median_K": dTb_med,
            "n_ocean_nodes": int((cls == OCEAN).sum()),
            "n_frozen_nodes": int((cls == FROZEN).sum()),
            "n_none_nodes": int((cls == NONE).sum()),
        },
        "prior_box": {"Tb_K": [tb_lo, tb_hi], "log10_wOcean_ppt": [lw_lo, lw_hi]},
        "eval_grid": {"n_tb": args.n_tb, "n_w": args.n_w,
                      "n_buildable_support_points": n_build},
        "PRIMARY_binding_onset_placement": {
            "half_cell_K": half_cell_K,
            "max_local_dTb_K": max_local_dTb,
            "bar_K": args.half_cell_bar_K,
            "pass": bool(half_cell_K <= args.half_cell_bar_K),
            "per_column_onset": onset_cols,
            "note": ("BINDING metric (reviewer ab2d3748d0b00c78c 2026-08-07): "
                     "absolute onset-temperature placement uncertainty = half the "
                     "LOCAL Tb cell spanning the freeze onset, worst-case over w. "
                     "Replaces the box-fraction span as the launch gate."),
        },
        "sandwich_invariant": {
            "violations": sandwich_violations,
            "pass": bool(len(sandwich_violations) == 0),
            "note": ("no None/frozen node between two ocean nodes in a w-column; "
                     "catches the Tb=252 K int/float clamp artifact + regressions."),
        },
        "ocean_fraction_DIAGNOSTIC": {
            "half_cell_K_global_median": half,
            "nominal": f0,
            "boundary_cold_shift_-half": fcold,
            "boundary_warm_shift_+half": fwarm,
            "placement_band_span": span,
            "placement_half_span": span / 2.0,
            "max_deviation_from_nominal": max_dev,
            "note": ("SECONDARY/diagnostic only (reviewer): the fraction span "
                     "penalizes a narrow prior box and inverts the true placement "
                     "ranking. Reported for continuity with pre-fix numbers; NOT "
                     "the binding gate."),
        },
        "threshold_provisional_span": args.threshold,
        "binding_verdict": "PASS" if binding_pass else "FAIL",
        "verdict_note": ("Binding PASS requires half_cell_K <= bar_K AND zero "
                         "sandwich violations. Ocean-fraction span is diagnostic "
                         "only. Scientific-reviewer adjudicates the final call."),
    }

    if args.corner_w_min is not None:
        c_lo = np.log10(args.corner_w_min)
        cf0, cn = frac(0.0, corner=c_lo)
        ccold, _ = frac(-half, corner=c_lo)
        cwarm, _ = frac(+half, corner=c_lo)
        cspan = abs(ccold - cwarm)
        result["corner"] = {
            "w_min_ppt": args.corner_w_min,
            "log10_w_min": c_lo,
            "n_buildable_support_points": cn,
            "ocean_fraction": {
                "nominal": cf0,
                "boundary_cold_shift_-half": ccold,
                "boundary_warm_shift_+half": cwarm,
                "placement_band_span": cspan,
                "max_deviation_from_nominal": max(abs(ccold - cf0), abs(cwarm - cf0)),
            },
            "advisory_verdict_diagnostic": "PASS" if cspan <= args.threshold else "FAIL",
            "note": ("MgSO4 w>=150 corner: catches a phantom-ocean leak near the "
                     "onset (reviewer launch-condition 2). Diagnostic span only; the "
                     "binding gate is the absolute half_cell_K above."),
        }

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    print(f"\n[gate3] report -> {args.out}")


if __name__ == "__main__":
    main()
