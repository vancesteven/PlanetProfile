"""Per-node elastic libration factor sidecar for a zbw cache.

For every non-None OCEAN node of a v3.2-zbw-joint cache, computes

    F(zb, w)   = librations(elastic, ALMA y1) / librations(repaired rigid)

on the node's reduced stack (both WITHOUT the H22 figure treatment, so F is
C2- and H22-independent by construction -- binding review 2026-11-10, runtime
channel B), plus the linear rho_rock sensitivity

    c_rho(zb, w) = d ln F / d rho_rock      [per kg/m3]

by a +/-100 kg/m3 central difference taken ALONG THE RUNTIME MASS CLOSURE
(``mass_rho_rock``: rho_o, rho_i, R_o, R held fixed, R_s re-solved from the
reduced stack's own mass, so the unperturbed closure is an identity no-op --
confirming-review condition (a), 2026-11-10), with the node y1 held fixed.
The runtime stack carries the SAMPLED rho_rock while the node was built at
the build density, so the runtime factor is
F * (1 + c_rho * (rho_rock_runtime - rho_rock_node)).

Writes ``<cache>.lib_factor_sidecar.pkl`` with the grids, the three arrays,
NaN at unbuilt nodes, and provenance (git SHA, cache path+md5, ALMA order).

  plans/scripts/enceladus_libration_factor_sidecar.py \
      --cache PATH [--out PATH] [--drho 100] [--limit N]
"""
import argparse
import hashlib
import json
import pickle
import subprocess
import sys
import time
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

REPAIRED = dict(formulation="vh2008_repaired")
ALMA_ORDER = 4


def _node_factor(structure, drho):
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
        "periodic", ALMA_ORDER, verbose=False, internal=True)
    M, R = float(structure["Mtot_kg"]), float(structure["R_body_m"])
    y1 = np.real(rf[:, 0]) * M / R / (Constants.G * M / R ** 2)
    rr, rrho, oi, ry = LibrationModelInputs(r, ph, rho, y1)
    kw = dict(omega=om, ecc=float(structure["eccentricity"]), ocean=True,
              ocean_idx=oi)

    def factor(stack_r, stack_rho):
        return (librations(stack_r, stack_rho, rigid=False, y1=ry, **kw)
                / librations(stack_r, stack_rho, **REPAIRED, **kw))

    F = factor(rr, rrho)

    # rho_rock FD along the runtime mass closure (mass_rho_rock): hold
    # rho_o, rho_i, R_o, R and the stack's own mass, re-solve R_s. At the
    # node's own rho_s this returns the node stack exactly (identity).
    r_in = np.concatenate(([0.0], rr[:-1]))
    m_stack = float(np.sum(rrho * (4.0 * np.pi / 3.0) * (rr ** 3 - r_in ** 3)))
    rho_o, rho_i = float(rrho[1]), float(rrho[2])
    R_o, R = float(rr[1]), float(rr[2])
    a_cl = (3.0 * m_stack / (4.0 * np.pi) - rho_o * R_o ** 3
            - rho_i * (R ** 3 - R_o ** 3))

    def closed(rho_s):
        R_s = (a_cl / (rho_s - rho_o)) ** (1.0 / 3.0)
        return (np.array([R_s, R_o, R]),
                np.array([rho_s, rho_o, rho_i]))

    rho_s0 = float(rrho[0])
    c_rho = (np.log(factor(*closed(rho_s0 + drho)))
             - np.log(factor(*closed(rho_s0 - drho)))) / (2.0 * drho)
    return F, c_rho, rho_s0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--drho", type=float, default=100.0)
    ap.add_argument("--limit", type=int, default=None,
                    help="stop after N built nodes (smoke runs)")
    args = ap.parse_args()
    out = args.out or args.cache.with_name(
        args.cache.name + ".lib_factor_sidecar.pkl")

    md5 = hashlib.md5(args.cache.read_bytes()).hexdigest()
    with args.cache.open("rb") as stream:
        cache = pickle.load(stream)
    zb = np.asarray(cache["zb_km_grid"], float)
    w = np.asarray(cache["wOcean_ppt_grid"], float)
    nz, nw = zb.size, w.size
    F = np.full((nz, nw), np.nan)
    c_rho = np.full((nz, nw), np.nan)
    rho_node = np.full((nz, nw), np.nan)

    t0 = time.time()
    done = 0
    for i in range(nz):
        for j in range(nw):
            st = cache["structures"][i * nw + j]
            if st is None:
                continue
            F[i, j], c_rho[i, j], rho_node[i, j] = _node_factor(st, args.drho)
            done += 1
            if args.limit and done >= args.limit:
                break
        else:
            continue
        break

    sha = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                         capture_output=True, text=True).stdout.strip()
    sidecar = {
        "schema": "lib-factor-sidecar-v1",
        "zb_km_grid": zb, "wOcean_ppt_grid": w,
        "F": F, "dlnF_drho_rock": c_rho, "rho_rock_node_kgm3": rho_node,
        "drho_fd_kgm3": args.drho, "alma_order": ALMA_ORDER,
        "definition": "F = elastic(ALMA y1, no figure) / repaired rigid "
                      "(no figure); runtime factor = F * (1 + dlnF_drho_rock"
                      " * (rho_rock_runtime - rho_rock_node_kgm3))",
        "provenance": {"cache": str(args.cache), "cache_md5": md5,
                       "build_commit": sha,
                       "built_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ",
                                                  time.gmtime())},
    }
    with out.open("wb") as stream:
        pickle.dump(sidecar, stream)
    finite = np.isfinite(F)
    print(json.dumps({
        "out": str(out), "built_nodes": int(finite.sum()),
        "elapsed_s": round(time.time() - t0, 1),
        "F_min": float(np.nanmin(F)), "F_max": float(np.nanmax(F)),
        "c_rho_min": float(np.nanmin(c_rho)),
        "c_rho_max": float(np.nanmax(c_rho)),
        "cache_md5": md5,
    }, indent=1))


if __name__ == "__main__":
    main()
