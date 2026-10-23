#!/usr/bin/env python
"""rho_grain,max for the mass_moi 'rho_s_above_grain_max' loose reject.

Zero-porosity rock (grain) density of the Perple_X table the Enceladus
template ACTUALLY loads (Planet.Sil.mantleEOS, resolved through
PlanetProfile._PERPLEXCACHE exactly as InnerEOS.PerplexEOSStruct does), read
at seafloor conditions. The Perple_X table is the zero-porosity solid; PP
applies the porosity (Han2014) afterwards, so the table density IS the grain
density. Value stored in the 7D config as
metadata.isostasy.rho_grain_max_kgm3 (with provenance) -- never hard-coded in
the runner.

  python plans/scripts/enceladus_rho_grain_max.py [--p-mpa-max 10]
         [--t-k-max 290] [--json-out FILE] [--table PATH]

Reported: the table path + md5 (compare between machines: the container and
the Mac may hold different tables), the density at the reference point
(P=6.1 MPa -- the PHydroMax PP prints for zb~39.5 km -- T=273.15 K) and the
max over the box P in [table Pmin, p-mpa-max], T in [273.15, t-k-max]
(the stored value; a conservative upper bound).
"""
import argparse
import hashlib
import importlib
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
TEMPLATE = 'PlanetProfile.Default.Enceladus.PPEnceladus'
NHEADERS = 13


def load_table(path):
    cols = np.loadtxt(path, skiprows=NHEADERS, unpack=True)
    first, second, rho = cols[0], cols[1], cols[2]
    with open(path) as f:
        for _ in range(NHEADERS - 1):
            f.readline()
        hdr = f.readline().strip()
    dim2 = next(i for i in range(second.size) if second[i] != second[0])
    if hdr[0] == 'T':          # T fastest
        T = first[:dim2]
        P_bar = second[::dim2]
        rho = rho.reshape(P_bar.size, dim2)       # [P, T]
    else:                      # P fastest
        P_bar = first[:dim2]
        T = second[::dim2]
        rho = rho.reshape(T.size, dim2).T         # [P, T]
    rho = np.where(rho > 0, rho, np.nan)
    return P_bar * 0.1, T, rho                    # P in MPa


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--table')
    ap.add_argument('--p-mpa-max', type=float, default=10.0)
    ap.add_argument('--t-k-max', type=float, default=290.0)
    ap.add_argument('--p-ref-mpa', type=float, default=6.1)
    ap.add_argument('--t-ref-k', type=float, default=273.15)
    ap.add_argument('--json-out')
    a = ap.parse_args()
    tmpl = importlib.import_module(TEMPLATE).Planet
    eos = tmpl.Sil.mantleEOS
    if a.table:
        path = a.table
    else:
        from PlanetProfile import _PERPLEXCACHE
        path = str(Path(_PERPLEXCACHE) / eos)
    P, T, rho = load_table(path)
    from scipy.interpolate import RegularGridInterpolator
    f = RegularGridInterpolator((P, T), rho, method='linear')
    ref = float(f([[a.p_ref_mpa, a.t_ref_k]])[0])
    Pb = np.linspace(P.min(), a.p_mpa_max, 60)
    Tb = np.linspace(273.15, a.t_k_max, 60)
    PP_, TT_ = np.meshgrid(Pb, Tb, indexing='ij')
    box = f(np.stack([PP_.ravel(), TT_.ravel()], -1))
    out = {'template': TEMPLATE, 'mantleEOS': eos, 'table_path': path,
           'table_md5': hashlib.md5(Path(path).read_bytes()).hexdigest(),
           'table_P_MPa_range': [float(P.min()), float(P.max())],
           'table_T_K_range': [float(T.min()), float(T.max())],
           'rho_ref_kgm3': ref, 'ref_P_MPa': a.p_ref_mpa, 'ref_T_K': a.t_ref_k,
           'box_P_MPa': [float(Pb.min()), float(Pb.max())],
           'box_T_K': [float(Tb.min()), float(Tb.max())],
           'rho_grain_max_kgm3': float(np.nanmax(box)),
           'rho_box_min_kgm3': float(np.nanmin(box))}
    print(json.dumps(out, indent=1))
    if a.json_out:
        Path(a.json_out).write_text(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
