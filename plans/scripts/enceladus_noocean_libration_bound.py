#!/usr/bin/env python
"""Can a NO-OCEAN Enceladus with low-viscosity ice reach the observed forced
libration (Park et al. 2024 erratum: 0.092 +/- 0.003 deg)?

User question 2026-10-28 (prompted by Downey & Nimmo 2025, Sci. Adv. 11,
eadl4741, whose shell/interior coupling parameter K/C_s is the same quantity
used here). Two bounds, both computed with PlanetProfile's own libration
machinery (PlanetProfile/Gravity/Librations.py):

(A) WHOLE-BODY VISCOELASTIC (no decoupled layer). Van Hoolst et al. (2013)
    oceanless non-rigid branch, librations(ocean=False, rigid=False,
    rek2=k2) scanned over Re(k2) in [0, k_f]. Deformation lowers the
    amplitude, so the rigid value is the maximum.

(B) DECOUPLED SHELL OVER A SOFT (viscous / Maxwell) ICE LAYER. Three-layer
    stack: uniform rock core / soft ice layer of thickness h / ice shell of
    thickness zb, mass-closed to Enceladus' mass. The Van Hoolst et al. (2008)
    rigid-shell libration system (the code's ocean=True, rigid=True,
    hydrostatic branch, reproduced exactly at zero coupling and asserted)
    plus a frequency-domain viscoelastic coupling torque between shell and
    interior transmitted by the soft layer:
        K~ = 8 pi * i n * eta_eff * r1^3 r2^3 / (r2^3 - r1^3)
        eta_eff = eta mu / (mu + i n eta)          (Maxwell)
    the exact Stokes torque between rotating concentric spheres (thin-layer
    limit 8 pi eta R^4 / 3h), generalised to a Maxwell rheology. The shell
    amplitude |Gamma_s| is solved from the 2x2 complex system. eta -> 0 is the
    inviscid (ocean-like) limit; eta -> inf is elastic (mu) coupling.
    The shell is rigid here, which OVER-predicts the amplitude (shell
    elasticity lowers it, Van Hoolst et al. 2013), so (B) is an upper bound
    on the soft-layer scenario.
    The soft layer's own inertia is not carried (as in the ocean model); the
    strong-coupling limit is reported against the whole-body rigid value.

For each (zb, h, rho_rock) the script reports the viscosity eta* below which
the amplitude first reaches the observed value minus 3 sigma (0.083 deg).

  python plans/scripts/enceladus_noocean_libration_bound.py --out-dir DIR
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from PlanetProfile.Gravity.Librations import (  # noqa: E402
    librations, compute_eccentricities, radial2abc, moi_ellps)
from PlanetProfile.Gravity.Librations import G  # noqa: E402  (same constant as the module)
M = 1.08022e20            # Enceladus mass (frozen build spec)
R = 252.1e3               # mean radius (frozen build spec)
N = 2 * np.pi / (1.370218 * 86400)   # orbital = spin rate
ECC = 0.0047
OBS, SIG = 0.092, 0.003
TARGET = OBS - 3 * SIG    # 0.083 deg


def stack(zb, h, rho_rock, rho_ice=925.0, rho_soft=925.0):
    """Mass-closed 3-layer stack (core, soft layer, shell); None if infeasible."""
    r2, r1 = R - zb, R - zb - h
    a = 3 * M / (4 * np.pi) - rho_ice * (R ** 3 - r2 ** 3) - rho_soft * (r2 ** 3 - r1 ** 3)
    if a <= 0 or rho_rock <= rho_soft:
        return None
    # core radius from mass: rho_rock Rc^3 + rho_soft (r1^3 - Rc^3) = a
    Rc3 = (a - rho_soft * r1 ** 3) / (rho_rock - rho_soft)
    if not (0 < Rc3 < r1 ** 3):
        return None
    Rc = Rc3 ** (1 / 3)
    return np.array([Rc, r1, r2, R]), np.array([rho_rock, rho_soft, rho_soft, rho_ice])


def vh2008_terms(radial, rho, omega, ocean_idx):
    """Hydrostatic Van Hoolst (2008) terms exactly as Librations.librations
    (ocean=True, rigid=True, H22_obs_m=None)."""
    ep, eq = compute_eccentricities(radial, rho, omega)
    a, b, c = radial2abc(radial, ep, eq)
    Ai, Bi, Ci = moi_ellps(a, b, c, rho)
    As, Bs, Cs = Ai[ocean_idx + 1:].sum(), Bi[ocean_idx + 1:].sum(), Ci[ocean_idx + 1:].sum()
    Ac, Bc, Cc = Ai[:ocean_idx].sum(), Bi[:ocean_idx].sum(), Ci[:ocean_idx].sum()
    beta = (a - b) / a
    Bip_Aip = 8 / 15 * np.pi * rho[ocean_idx] * beta[0] * radial[0] ** 5
    Bsp_Asp = 8 / 15 * np.pi * rho[ocean_idx] * beta[ocean_idx] * radial[ocean_idx] ** 5
    Ks = 3 * omega ** 2 * ((Bs - As) + Bsp_Asp)
    Kc = 3 * omega ** 2 * ((Bc - Ac) - Bip_Aip)
    K_int = 4 * np.pi * G / 5 * (rho[-1] * beta[-1] + (rho[1] - rho[-1]) * beta[1]) \
        * ((rho[0] - rho[1]) * beta[0] * radial[0] ** 5)
    return dict(Ks=Ks, Kc=Kc, K_int=K_int, Cs=Cs, Cc=Cc,
                A=Ai.sum(), B=Bi.sum(), C=Ci.sum())


def shell_amplitude_deg(T, omega, ecc, Kv):
    """|Gamma_s| in degrees with complex coupling Kv (0 => VH2008 rigid)."""
    A = np.array([[-omega ** 2 * T['Cs'] + T['Ks'] + 2 * T['K_int'] + Kv, -(2 * T['K_int'] + Kv)],
                  [-(2 * T['K_int'] + Kv), -omega ** 2 * T['Cc'] + T['Kc'] + 2 * T['K_int'] + Kv]],
                 dtype=complex)
    b = np.array([2 * ecc * T['Ks'], 2 * ecc * T['Kc']], dtype=complex)
    return float(np.degrees(abs(np.linalg.solve(A, b)[0])))


def coupling(eta, mu, r1, r2, omega):
    eta_eff = eta * mu / (mu + 1j * omega * eta) if np.isfinite(mu) else eta
    return 8 * np.pi * 1j * omega * eta_eff * r1 ** 3 * r2 ** 3 / (r2 ** 3 - r1 ** 3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--mu', type=float, default=3.3e9, help='ice shear modulus, Pa')
    a = ap.parse_args()
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    etas = np.logspace(3, 18, 301)
    rep = {'inputs': {'M': M, 'R': R, 'n': N, 'ecc': ECC, 'obs': OBS, 'sigma': SIG,
                      'target_obs_minus_3sigma': TARGET, 'mu_ice': a.mu,
                      'maxwell_eta_mu_over_n': a.mu / N},
           'checks': {}, 'A_whole_body': {}, 'B_soft_layer': []}

    # --- (A) whole-body: rigid vs Van Hoolst 2013 deformable, Re(k2) scan ---
    A_rows = []
    for rr in (2200., 2400., 2600.):
        for zb in (20e3, 45e3, 65e3):
            s = stack(zb, 1e-3, rr)       # 2-layer equivalent (h -> 0)
            if s is None:
                continue
            radial, rho = s
            rigid = float(np.degrees(librations(radial, rho, N, ECC, rigid=True, ocean=False) / R))
            # fluid Love number of the body (upper end of Re k2)
            kf_scan = np.linspace(0.0, 1.5, 31)
            deform = [float(np.degrees(librations(radial, rho, N, ECC, rigid=False,
                                                  ocean=False, rek2=k) / R)) for k in kf_scan]
            A_rows.append({'rho_rock': rr, 'zb_km': zb / 1e3, 'rigid_deg': rigid,
                           'max_over_rek2_deg': max(deform),
                           'monotone_decreasing_in_rek2': bool(np.all(np.diff(deform) <= 1e-15)),
                           'rek2_scan': kf_scan.tolist(), 'deg': deform})
    rep['A_whole_body'] = {
        'rows': A_rows,
        'max_rigid_deg': max(r['rigid_deg'] for r in A_rows),
        'max_any_deg': max(max(r['max_over_rek2_deg'], r['rigid_deg']) for r in A_rows)}

    # --- (B) decoupled shell over a soft Maxwell / viscous ice layer ---
    worst_check = 0.0
    for rr in (2200., 2400., 2600.):
        for zb in (5e3, 10e3, 20e3, 27e3, 45e3):
            for h in (1e3, 10e3, 30e3, 100e3):
                s = stack(zb, h, rr)
                if s is None:
                    continue
                radial, rho = s
                T = vh2008_terms(radial, rho, N, ocean_idx=1)
                ref = float(np.degrees(librations(radial, rho, N, ECC, rigid=True,
                                                  ocean=True, ocean_idx=1) / R))
                inv = shell_amplitude_deg(T, N, ECC, 0.0)
                worst_check = max(worst_check, abs(inv - ref) / ref)
                r1, r2 = radial[1], radial[2]
                rows = {}
                for mu_lab, mu in (('maxwell', a.mu), ('viscous', np.inf)):
                    amp = np.array([shell_amplitude_deg(T, N, ECC, coupling(e, mu, r1, r2, N))
                                    for e in etas])
                    reach = amp >= TARGET
                    eta_star = float(etas[reach][-1]) if reach.any() else None
                    rows[mu_lab] = {'eta_star_Pas': eta_star,
                                    'amp_at_eta_1e13': float(np.interp(13, np.log10(etas), amp)),
                                    'amp_at_eta_1e14': float(np.interp(14, np.log10(etas), amp)),
                                    'amp_strong_coupling': float(amp[-1])}
                rigid_whole = float(np.degrees(librations(radial, rho, N, ECC, rigid=True,
                                                          ocean=False) / R))
                rep['B_soft_layer'].append({
                    'rho_rock': rr, 'zb_km': zb / 1e3, 'h_km': h / 1e3,
                    'inviscid_limit_deg': inv, 'librations_ocean_ref_deg': ref,
                    'whole_body_rigid_deg': rigid_whole, **rows})
    rep['checks']['zero_coupling_reproduces_librations_max_rel_diff'] = worst_check
    assert worst_check < 1e-10, f'zero-coupling replication off by {worst_check}'

    B = rep['B_soft_layer']
    es = [r['maxwell']['eta_star_Pas'] for r in B if r['maxwell']['eta_star_Pas']]
    rep['summary'] = {
        'whole_body_max_deg': rep['A_whole_body']['max_any_deg'],
        'soft_layer_cases': len(B),
        'cases_reaching_target_at_any_eta': len(es),
        'max_eta_star_maxwell_Pas': max(es) if es else None,
        'max_amp_at_eta_1e13_deg': max(r['maxwell']['amp_at_eta_1e13'] for r in B),
        'max_amp_at_eta_1e14_deg': max(r['maxwell']['amp_at_eta_1e14'] for r in B),
    }
    (out / 'noocean_libration_bound.json').write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep['summary'], indent=1))
    print('zero-coupling replication max rel diff', worst_check)
    print('(A) whole body:', [(r['rho_rock'], r['zb_km'], round(r['rigid_deg'], 4),
                              round(r['max_over_rek2_deg'], 4), r['monotone_decreasing_in_rek2'])
                             for r in A_rows])
    for r in B:
        if r['rho_rock'] == 2400.:
            print('(B) zb %4.0f h %4.0f  inviscid %.4f  eta*(Maxwell) %s  amp@1e13 %.4f  strong %.4f  whole-rigid %.4f'
                  % (r['zb_km'], r['h_km'], r['inviscid_limit_deg'],
                     '%.2e' % r['maxwell']['eta_star_Pas'] if r['maxwell']['eta_star_Pas'] else 'never',
                     r['maxwell']['amp_at_eta_1e13'], r['maxwell']['amp_strong_coupling'],
                     r['whole_body_rigid_deg']))


if __name__ == '__main__':
    main()
