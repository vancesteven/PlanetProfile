"""F-sidecar acceptance checks (confirming review 2026-11-10, blockers D1-D2).

Writes validation_reports/enceladus_isostasy/kint_switchover/
f_table_checks.json with:

1. ``loo``        -- leave-one-out zb-linear interpolation error of F on a
                     2x-coarsened grid, zb 20-32 (gate: max < 1e-3).
2. ``rho_ice``    -- |dF/F| for rho_ice +/-10 kg/m3 at nodes 1220/1540,
                     closure-consistent (rho_i enters the mass closure),
                     node y1 fixed (gate: < 1e-3).
3. ``c_vs_alma``  -- required test (a)-2: the sidecar's c = dlnF/drho_rock
                     at nodes 1220/1540 against a fresh closure FD with
                     ALMA y1 RECOMPUTED per perturbed stack (the node
                     profile's regions reduced to 3 layers x ``--sub``
                     constant-property sublayers, volume-weighted mu/eta,
                     so the full-vs-reduced ALMA baseline cancels in the
                     derivative). Gate: |c_alma/c_sidecar - 1| < 0.05.
4. ``linearity``  -- required test (a)-3: F * (1 + c*(rho_s - rho_ref)) at
                     rho_s = 2411 against the exact nonlinear factor on the
                     mass-closed stack (node y1). Gate: |diff| < 1e-4.

  plans/scripts/kint_f_table_checks.py [--cache PATH] [--sidecar PATH]
      [--nodes 1220 1540] [--rho-probe 2411] [--sub 12]
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
REPAIRED = dict(formulation="vh2008_repaired")


def _full_profile(st):
    return (np.asarray(st["r_m"], float), np.asarray(st["rho"], float),
            np.asarray(st["mu_Pa"], float),
            np.asarray(st["eta_Pa_base"], float), np.asarray(st["phases"]))


def _node_y1_inputs(st):
    r, rho, mu, eta, ph = _full_profile(st)
    rheo = np.where(ph == 0, "newton", "elastic").tolist()
    model = pyALMA3Updated.build_model(r, rho, mu, eta, rheo,
                                       np.zeros((r.size, 2)), norm=False)
    om = float(st["omega"])
    _, _, _, rf = pyALMA3Updated.love_numbers(
        model, [2], [(2 * np.pi / om) / ALMA_TIME_UNIT_S], "tidal",
        "periodic", 4, verbose=False, internal=True)
    M, R = float(st["Mtot_kg"]), float(st["R_body_m"])
    y1 = np.real(rf[:, 0]) * M / R / (Constants.G * M / R ** 2)
    rr, rrho, oi, ry = LibrationModelInputs(r, ph, rho, y1)
    kw = dict(omega=om, ecc=float(st["eccentricity"]), ocean=True,
              ocean_idx=oi)
    return rr, rrho, ry, kw


def _closure(rr, rrho):
    """mass_rho_rock closure around the node's own reduced stack."""
    r_in = np.concatenate(([0.0], rr[:-1]))
    m = float(np.sum(rrho * (4.0 * np.pi / 3.0) * (rr ** 3 - r_in ** 3)))
    rho_o, rho_i = float(rrho[1]), float(rrho[2])
    R_o, R = float(rr[1]), float(rr[2])

    def closed(rho_s=None, rho_ice=None):
        rho_s = float(rrho[0]) if rho_s is None else float(rho_s)
        rho_i2 = rho_i if rho_ice is None else float(rho_ice)
        a = (3.0 * m / (4.0 * np.pi) - rho_o * R_o ** 3
             - rho_i2 * (R ** 3 - R_o ** 3))
        R_s = (a / (rho_s - rho_o)) ** (1.0 / 3.0)
        return np.array([R_s, R_o, R]), np.array([rho_s, rho_o, rho_i2])
    return closed


def _vw(x, r, lo, hi):
    """Volume-weighted mean of profile x over radius region (lo, hi]."""
    r_in = np.concatenate(([0.0], r[:-1]))
    dv = r ** 3 - r_in ** 3
    m = (r > lo) & (r <= hi * (1 + 1e-12))
    return float(np.sum(x[m] * dv[m]) / np.sum(dv[m]))


def _alma_reduced_factor(st, stack_r, stack_rho, kw, sub):
    """elastic/repaired with y1 RECOMPUTED by ALMA on the closed stack,
    each of the 3 regions as ``sub`` constant-property sublayers with
    volume-weighted mu/eta from the node profile's regions."""
    r, _, mu, eta, _ = _full_profile(st)
    rr0, _, _, _ = _node_y1_inputs.cache[id(st)] if False else (None,) * 4
    bounds = [0.0] + [float(x) for x in stack_r]
    props = []
    # mu/eta regions from the NODE profile boundaries (same ordering:
    # core, ocean, shell), mapped onto the closed stack's radii.
    node_bounds = _node_bounds_cache[id(st)]
    for k in range(3):
        props.append((_vw(mu, r, node_bounds[k], node_bounds[k + 1]),
                      _vw(eta, r, node_bounds[k], node_bounds[k + 1])))
    rs, rhos, mus, etas, rheos = [], [], [], [], []
    for k in range(3):
        edges = np.linspace(bounds[k], bounds[k + 1], sub + 1)[1:]
        rs.extend(edges)
        rhos.extend([stack_rho[k]] * sub)
        mus.extend([props[k][0]] * sub)
        etas.extend([props[k][1]] * sub)
        rheos.extend(["newton" if k == 1 else "elastic"] * sub)
    rs = np.asarray(rs)
    model = pyALMA3Updated.build_model(
        rs, np.asarray(rhos), np.asarray(mus), np.asarray(etas), rheos,
        np.zeros((rs.size, 2)), norm=False)
    om = kw["omega"]
    _, _, _, rf = pyALMA3Updated.love_numbers(
        model, [2], [(2 * np.pi / om) / ALMA_TIME_UNIT_S], "tidal",
        "periodic", 4, verbose=False, internal=True)
    M = float(np.sum(np.asarray(rhos) * (4.0 * np.pi / 3.0)
                     * (rs ** 3 - np.concatenate(([0.0], rs[:-1])) ** 3)))
    R = rs[-1]
    y1 = np.real(rf[:, 0]) * M / R / (Constants.G * M / R ** 2)
    ph = np.where(np.asarray(rheos) == "newton", 0, 1)
    rr, rrho, oi, ry = LibrationModelInputs(rs, ph, np.asarray(rhos), y1)
    kw2 = dict(kw, ocean_idx=oi)
    return (librations(rr, rrho, rigid=False, y1=ry, **kw2)
            / librations(rr, rrho, **REPAIRED, **kw2))


_node_bounds_cache = {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=CACHE)
    ap.add_argument("--sidecar", type=Path, default=None)
    ap.add_argument("--nodes", type=int, nargs="+", default=[1220, 1540])
    ap.add_argument("--rho-probe", type=float, default=2411.0)
    ap.add_argument("--drho", type=float, default=100.0)
    ap.add_argument("--sub", type=int, default=12)
    ap.add_argument("--out", type=Path,
                    default=REPO / "validation_reports/enceladus_isostasy"
                    / "kint_switchover/f_table_checks.json")
    args = ap.parse_args()
    sidecar_path = args.sidecar or args.cache.with_name(
        args.cache.name + ".lib_factor_sidecar.pkl")
    with args.cache.open("rb") as s:
        cache = pickle.load(s)
    with sidecar_path.open("rb") as s:
        sc = pickle.load(s)
    zb = np.asarray(cache["zb_km_grid"], float)
    w = np.asarray(cache["wOcean_ppt_grid"], float)
    nw = w.size
    F_tab = np.asarray(sc["F"], float)
    c_tab = np.asarray(sc["dlnF_drho_rock"], float)
    rho_ref = np.asarray(sc["rho_rock_node_kgm3"], float)
    out = {"cache": str(args.cache), "sidecar": str(sidecar_path),
           "sidecar_md5_of_cache": sc["provenance"]["cache_md5"]}

    # 1. LOO zb-linear interpolation on a 2x-coarsened grid, zb 20-32.
    izbs = [i for i in range(zb.size) if 20 <= zb[i] <= 32][::2]
    iws = list(range(0, nw, 6))
    errs = []
    for j in iws:
        col = [(zb[i], F_tab[i, j]) for i in izbs
               if np.isfinite(F_tab[i, j])]
        for k in range(1, len(col) - 1):
            z0, f0 = col[k - 1]
            z1, f1 = col[k + 1]
            zm, fm = col[k]
            errs.append(abs(f0 + (f1 - f0) * (zm - z0) / (z1 - z0) - fm))
    errs = np.asarray(errs)
    out["loo"] = {"max": float(errs.max()), "mean": float(errs.mean()),
                  "n": int(errs.size), "pass": bool(errs.max() < 1e-3)}

    out["nodes"] = {}
    for idx in args.nodes:
        st = cache["structures"][idx]
        r_full, *_ = _full_profile(st)
        rr, rrho, ry, kw = _node_y1_inputs(st)
        _node_bounds_cache[id(st)] = [0.0] + [float(x) for x in rr]
        closed = _closure(rr, rrho)

        def fac(stack):
            return (librations(stack[0], stack[1], rigid=False, y1=ry, **kw)
                    / librations(stack[0], stack[1], **REPAIRED, **kw))

        i, j = idx // nw, idx % nw
        F0, c0, rref = F_tab[i, j], c_tab[i, j], rho_ref[i, j]
        rec = {"F_sidecar": float(F0), "c_sidecar": float(c0),
               "rho_ref": float(rref)}

        # 2. rho_ice +/-10, closure-consistent, node y1.
        d_ice = max(abs(fac(closed(rho_ice=rrho[2] + s)) / F0 - 1)
                    for s in (-10.0, +10.0))
        rec["rho_ice_pm10_absdF"] = float(d_ice)
        rec["rho_ice_pass"] = bool(d_ice < 1e-3)

        # 3. c with ALMA y1 recomputed per perturbed stack (test (a)-2).
        f_hi = _alma_reduced_factor(st, *closed(rho_s=rref + args.drho),
                                    kw=kw, sub=args.sub)
        f_lo = _alma_reduced_factor(st, *closed(rho_s=rref - args.drho),
                                    kw=kw, sub=args.sub)
        c_alma = (np.log(f_hi) - np.log(f_lo)) / (2 * args.drho)
        rec["c_alma_recomputed"] = float(c_alma)
        rec["c_ratio"] = float(c_alma / c0)
        rec["c_pass"] = bool(abs(c_alma / c0 - 1) < 0.05)

        # 4. linearized vs exact at rho_probe (test (a)-3), node y1.
        lin = F0 * (1.0 + c0 * (args.rho_probe - rref))
        exact = fac(closed(rho_s=args.rho_probe))
        rec["rho_probe"] = args.rho_probe
        rec["F_linearized"] = float(lin)
        rec["F_exact"] = float(exact)
        rec["linearity_absdiff"] = float(abs(lin - exact))
        rec["linearity_pass"] = bool(abs(lin - exact) < 1e-4)

        out["nodes"][str(idx)] = rec

    out["pass"] = bool(out["loo"]["pass"] and all(
        r["rho_ice_pass"] and r["c_pass"] and r["linearity_pass"]
        for r in out["nodes"].values()))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
