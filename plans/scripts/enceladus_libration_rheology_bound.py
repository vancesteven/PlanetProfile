#!/usr/bin/env python
"""How much does ice-shell rheology move the Enceladus forced libration, and
what zb shift does that imply? (user question 2026-10-30: "libration depends
on viscosity though?")

The campaign forward model is the RIGID Van Hoolst et al. (2008) ocean
branch with a frozen +0.8 % correction and a +/-0.4 % residual
(libration_sys_frac). Here the DEFORMABLE branch (Van Hoolst et al. 2013,
librations(rigid=False, y1=...)) is evaluated with y1 from an ALMA
viscoelastic solve at the forcing period, exactly as
tests/librations_test.py::_cache_node_diagnostics does, while scanning the
shell's shear modulus mu and Maxwell viscosity eta.

KNOWN DEFECT carried, not fixed: the rigid branch's K_int lacks 8*pi/15
(pinned by test_adjudicated_rigid_elastic_discrepancy_is_pinned; the
deformable branch currently comes out ABOVE rigid, a strict-xfail
expectation). Results are therefore reported as RELATIVE changes of the
deformable branch with mu/eta (deformable(mu, eta) / deformable(reference)),
which do not depend on the rigid-branch normalization, plus the raw
deformable/rigid ratio for the record.

UNITS: pyALMA3Updated.build_model(norm=False) with the period in kyr (as
Gravity.ComputeGravityObservations passes it) needs viscosity in Pa*kyr for
the Maxwell term mu/eta to share the frequency's units. PlanetProfile's
gravity path passes Pa*s (Gravity.py:797-806), which makes every ice layer
~3e10x too viscous (elastic). This script converts eta to Pa*kyr; the
production-path defect is reported, not fixed here.

zb shift: d(libration)/d(zb) from a least-squares line through the rigid
libration of +/-3 production-cache nodes in the same w column (neighbouring
nodes alone were too noisy), so a fractional libration change f maps to
dzb = f * L / (dL/dzb).

  python plans/scripts/enceladus_libration_rheology_bound.py --out-dir DIR
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
    ALMA_TIME_UNIT_S, LibrationModelInputs)
from PlanetProfile.Gravity.Librations import librations  # noqa: E402
from PlanetProfile.Utilities.defineStructs import Constants  # noqa: E402

CACHE = (REPO / 'PlanetProfile/Test/mcmc_results/Enceladus/Cassini_isostasy/'
         'enceladus_seawater_zbw_grid_2d.pkl')
ICE = 1.0       # phase code of ice Ih in the cache
OCEAN = 0.0


def node_libration(s, mu_ice=None, eta_ice=None, shell_rheology='elastic'):
    """(rigid_deg, deformable_deg, k2) for node s, with optional ice mu/eta."""
    r = np.asarray(s['r_m'], float)
    rho = np.asarray(s['rho'], float)
    mu = np.array(s['mu_Pa'], float)
    eta = np.array(s['eta_Pa_base'], float)
    ph = np.asarray(s['phases'])
    ice = ph == ICE
    if mu_ice is not None:
        mu[ice] = mu_ice
    if eta_ice is not None:
        eta[ice] = eta_ice / ALMA_TIME_UNIT_S      # Pa*s -> Pa*kyr (see UNITS)
    rheo = np.where(ph == OCEAN, 'newton',
                    np.where(ice, shell_rheology, 'elastic')).tolist()
    model = pyALMA3Updated.build_model(r, rho, mu, eta, rheo,
                                       np.zeros((r.size, 2)), norm=False)
    om = float(s['omega'])
    per = (2 * np.pi / om) / ALMA_TIME_UNIT_S
    _, _, k, rf = pyALMA3Updated.love_numbers(model, [2], [per], 'tidal',
                                              'periodic', 4, verbose=False,
                                              internal=True)
    M, R = float(s['Mtot_kg']), float(s['R_body_m'])
    g = Constants.G * M / R ** 2
    y1 = np.real(rf[:, 0]) * M / R / g
    rr, rrho, oi, ry1 = LibrationModelInputs(r, ph, rho, y1)
    kw = dict(omega=om, ecc=float(s['eccentricity']), ocean=True, ocean_idx=oi)
    rig = librations(rr, rrho, rigid=True, **kw)
    dfm = librations(rr, rrho, rigid=False, y1=ry1, **kw)
    return (float(np.degrees(rig / R)), float(np.degrees(dfm / R)),
            complex(k[0, 0, -1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out-dir', required=True)
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    d = pickle.load(open(CACHE, 'rb'))
    zb, w, S = (np.asarray(d['zb_km_grid']), np.asarray(d['wOcean_ppt_grid']),
                d['structures'])
    nw = len(w)
    iw = int(np.argmin(np.abs(np.log10(w) - np.log10(3.5))))   # posterior-median w
    rows = []
    mus = [1e9, 2e9, 3.3e9, 4e9]
    etas = [1e12, 1e13, 3e13, 1e14, 3e14, 1e15, 1e17]
    n = 2 * np.pi / float(S[next(i for i, s in enumerate(S) if s)]['omega'])
    for iz_target in (20.0, 26.8, 35.0):
        iz = int(np.argmin(np.abs(zb - iz_target)))
        s = S[iz * nw + iw]
        if not s:
            continue
        # dL/dzb (rigid) from the neighbouring zb nodes at the same w
        nb = [S[j * nw + iw] for j in range(iz - 3, iz + 4)
              if 0 <= j < len(zb) and S[j * nw + iw]]
        zz = np.array([x['zb_km_actual'] for x in nb])
        LL = np.array([node_libration(x)[0] for x in nb])
        dLdz = float(np.polyfit(zz, LL, 1)[0])
        rig0, ref0, k0 = node_libration(s)                 # as-built elastic
        base = {'zb_km': s['zb_km_actual'], 'w_ppt': s['wOcean_ppt'],
                'rigid_deg': rig0, 'deformable_asbuilt_deg': ref0,
                'deformable_over_rigid': ref0 / rig0,
                'asbuilt_k2': [k0.real, k0.imag],
                'dL_dzb_deg_per_km': dLdz, 'scan': []}
        for mu in mus:
            _, de, k = node_libration(s, mu_ice=mu, shell_rheology='elastic')
            f = de / ref0 - 1
            base['scan'].append({'rheology': 'elastic', 'mu_Pa': mu,
                                 'deformable_deg': de, 'frac_vs_asbuilt': f,
                                 'dzb_km': f * ref0 / dLdz,
                                 'k2': [k.real, k.imag]})
            for eta in etas:
                _, de, k = node_libration(s, mu_ice=mu, eta_ice=eta,
                                          shell_rheology='maxwell')
                f = de / ref0 - 1
                base['scan'].append({'rheology': 'maxwell', 'mu_Pa': mu,
                                     'eta_Pas': eta,
                                     'maxwell_time_over_period': eta / mu / n,
                                     'deformable_deg': de,
                                     'frac_vs_asbuilt': f,
                                     'dzb_km': f * ref0 / dLdz,
                                     'k2': [k.real, k.imag]})
        rows.append(base)
        print(f"zb {base['zb_km']:.1f}: rigid {rig0:.4f} deformable {ref0:.4f} "
              f"(x{ref0 / rig0:.4f}) dL/dzb {dLdz:.5f} deg/km k2 {k0:.4f}")
        for r_ in base['scan']:
            print('   ', r_['rheology'], f"mu {r_['mu_Pa']:.1e}",
                  f"eta {r_.get('eta_Pas', float('nan')):.0e}",
                  f"frac {100 * r_['frac_vs_asbuilt']:+.2f} %",
                  f"dzb {r_['dzb_km']:+.2f} km",
                  f"k2 {r_['k2'][0]:.4f}{r_['k2'][1]:+.4f}i")
    (out / 'libration_rheology_bound.json').write_text(json.dumps(
        {'cache': str(CACHE.relative_to(REPO)), 'w_index': iw, 'rows': rows},
        indent=1))


if __name__ == '__main__':
    main()
