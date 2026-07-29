"""Titan NH3-ocean Phase 0 feasibility scan (Machine A spec 2026-07-25).

Maps the buildable (Tb, w_NH3) domain for Titan WITH a CoolProp NH3-H2O
ocean, using the reduced-stack / no-MoI-matching build path (the same
ConstantProps['Inner'] recipe that built titan_diff_noocean_structure_grid).
Reuses titan_tb_probe.probe_one so the build recipe is identical.

Spec grid: Tb in ~[240, 263] K x w in {0, 10, 30, 50, 100, 150} ppt.
Expected depression: ~1.1 K @ 10 ppt, ~5.4 K @ 50 ppt, ~15 K @ 150 ppt.

Deliverable: the (Tb, w) rectangle + node spacing for Phase 1, recording
which nodes produce ocean + ice-I shell + (any) HP ices, and
D_iceIh / D_ocean per node.

Run from repo root:
    mamba run -n PPcl env PYTHONPATH=. NUMBA_CACHE_DIR=/tmp/pp_numba_cache \
        KMP_DUPLICATE_LIB_OK=TRUE \
        python plans/scripts/titan_nh3_phase0_scan.py
"""
from __future__ import annotations

import json
import logging
import os
import sys
import time

if not os.environ.get("NUMBA_CACHE_DIR"):
    _nd = "/tmp/pp_numba_cache"
    os.makedirs(_nd, exist_ok=True)
    os.environ["NUMBA_CACHE_DIR"] = _nd
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

logging.basicConfig(level=logging.WARNING,
                    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
for _noisy in ["PlanetProfile", "SeaFreeze", "tidalpy", "TidalPy", "lbftd"]:
    logging.getLogger(_noisy).setLevel(logging.ERROR)

from plans.scripts.titan_tb_probe import probe_one  # noqa: E402

# Machine A spec grid. NH3 hard EOS Tb floor is 241 K; 240 is below it so
# start at 241. Upper end 263 per spec.
# w=0.0 dropped: 0-ppt NH3 is degenerate (== PureH2O); the pure/no-ocean
# case is Phase A. The melt-depression physics this scan probes needs w>0.
TB_GRID = [241.0, 244.0, 247.0, 250.0, 253.0, 256.0, 259.0, 262.0, 263.0]
W_GRID = [10.0, 30.0, 50.0, 100.0, 150.0]

# meltEOS resolution: 2 MPa (300 P-pts) instead of 0.05 default (12,000).
# The NH3 mu-liquidus is ~340 ms/pressure; the default made a single node
# ~68 min. See titan_tb_probe.probe_one docstring.
PFREEZE_RES_MPA = 2.0

# v2: EOS cleared per node. v1 (/tmp/titan_nh3_phase0_results.json) reused the
# EOS across the Tb axis and is contaminated (see per-node-clear comment above);
# kept on disk as evidence, NOT used as the Phase 1 deliverable.
OUT = "/tmp/titan_nh3_phase0_results_v2.json"


def main() -> None:
    rows = []
    t_start = time.time()
    n_total = len(TB_GRID) * len(W_GRID)
    n = 0
    print(f"# Titan NH3 Phase 0 scan: {len(TB_GRID)} Tb x {len(W_GRID)} w "
          f"= {n_total} nodes (PfreezeRes={PFREEZE_RES_MPA} MPa; EOS cleared "
          f"per node — matches deployed cache_builder.py:122)", flush=True)
    # EOS cache cleared on EVERY node. An earlier version reused the melt EOS
    # across the Tb axis (clear only on first Tb of each w) on the assumption
    # that the ocean+melt EOS depends only on (comp, w). That assumption is
    # EMPIRICALLY FALSE (2026-07-25 control /tmp/nh3_cache_control.py): when a
    # w-column's leading node builds a degenerate/failed melt curve, the
    # poisoned EOS leaks into every later Tb, pinning D_iceIh to a bogus
    # constant (w=10 column). Fresh-EOS builds gave the correct Tb-varying
    # structures (Tb=259/w=10: 105.5 km, not 48.45). The deployed builder
    # clears per node (cache_builder.py:122) for exactly this reason.
    for w in W_GRID:
        for i_tb, Tb in enumerate(TB_GRID):
            n += 1
            print(f"[{n}/{n_total}] NH3 Tb={Tb:.1f} w={w:.1f} "
                  f"(EOS cleared) ... ", end="", flush=True)
            t0 = time.time()
            try:
                r = probe_one(Tb, "NH3", w, clear_eos=True,
                              PfreezeRes_MPa=PFREEZE_RES_MPA)
            except Exception as e:  # noqa: BLE001
                r = {"Tb_K": Tb, "comp": "NH3", "w_ppt": w,
                     "status": "error", "error": repr(e)}
            r["wall_s"] = round(time.time() - t0, 1)
            rows.append(r)
            status = r.get("status", "?")
            d_oc = r.get("D_ocean_km")
            d_ih = r.get("D_iceIh_km")
            extra = ""
            if d_oc is not None:
                extra = f" D_ocean={d_oc:.1f} D_Ih={d_ih:.1f}"
            print(f"{status}{extra}  ({r['wall_s']}s)", flush=True)

    result = {
        "metadata": {
            "script": "plans/scripts/titan_nh3_phase0_scan.py",
            "spec": "plans/active/titan-nh3-ocean-campaign-spec.md",
            "TB_GRID": TB_GRID, "W_GRID": W_GRID,
            "comp": "NH3", "wall_total_s": round(time.time() - t_start, 1),
        },
        "rows": rows,
    }
    with open(OUT, "w") as f:
        json.dump(result, f, indent=2, default=str)
    print(f"\nWrote {OUT}  ({result['metadata']['wall_total_s']}s total)")

    # Coverage summary
    built = [r for r in rows if r.get("status") == "built"
             and r.get("D_ocean_km", 0) and r.get("D_ocean_km", 0) > 1.0]
    print(f"\nBuilt ocean-bearing nodes (D_ocean > 1 km): "
          f"{len(built)}/{n_total}")
    for w in W_GRID:
        tbs = sorted(r["Tb_K"] for r in built if r.get("w_ppt") == w)
        print(f"  w={w:6.1f} ppt: ocean Tb = {tbs}")


if __name__ == "__main__":
    main()
