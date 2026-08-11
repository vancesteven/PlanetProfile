"""B13 — the Hemingway & Mittal (2019) REPRODUCTION GATE.

The whole Enceladus shape-channel design rests on this test
(plans/active/enceladus-config-freeze.md): with H&M's inputs — Iess et
al. 2014 gravity + Thomas et al. 2016 libration + Tajeddine et al. 2017
shape, equal-pressures Airy isostasy, finite-amplitude correction,
rho_shell = 925, rho_ocean = 1020 kg/m^3 — our forward model's misfit
minimum must reproduce their preferred three-layer solution:

    shell 19-24 km, ocean 35-39 km (core 192-195 km follows from mass).

FALLBACK preregistered on failure: shape reverts to display-only and the
free dC20/dC22_nh nuisance boxes return (reviewer 2026-08-12).
"""
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from PlanetProfile.Gravity.isostasy import (  # noqa: E402
    airy_root, eccentricities_to_axes, finite_amplitude_coeffs,
    interface_gravity_coeff, triaxial_to_H2m, gravity_g)
from PlanetProfile.Gravity.Librations import (  # noqa: E402
    compute_eccentricities, librations)

# --- H&M 2019 inputs (Tajeddine block, Table 1; unnormalized) ----------
R_REF = 252.22e3          # Tajeddine mean radius = gravity reference
R_T = 252.22e3            # surface mean radius (same convention)
M_ENC = 1.08022e20
OMEGA = 5.307e-5          # rad/s
ECC = 0.0047              # orbital eccentricity
RHO_ICE = 925.0
RHO_OCEAN = 1020.0

H_OBS = {(2, 0): -3510.0, (2, 2): 857.0, (3, 0): 420.0}   # metres
C_OBS = {(2, 0): -5521e-6, (2, 2): 1574e-6, (3, 0): 118e-6}
SIG_C = {(2, 0): 35.8e-6, (2, 2): 15.9e-6, (3, 0): 23.8e-6}
LIB_OBS_DEG = 0.120       # Thomas et al. 2016
SIG_LIB = 0.014


def predict(D_shell_km: float, D_ocean_km: float):
    """Predicted (C20, C22, C30, libration_deg) for a 3-layer Enceladus."""
    R_b = R_T - D_shell_km * 1e3
    R_c = R_b - D_ocean_km * 1e3
    if R_c <= 50e3:
        return None
    V = 4.0 * np.pi / 3.0
    m_shell = RHO_ICE * V * (R_T ** 3 - R_b ** 3)
    m_ocean = RHO_OCEAN * V * (R_b ** 3 - R_c ** 3)
    rho_core = (M_ENC - m_shell - m_ocean) / (V * R_c ** 3)
    if not (1500.0 < rho_core < 4500.0):
        return None

    radial = np.array([R_c, R_b, R_T])
    rho = np.array([rho_core, RHO_OCEAN, RHO_ICE])

    # hydrostatic figures per interface (full iterative Tricarico)
    ep, eq = compute_eccentricities(radial, rho, OMEGA)
    H_hyd = []
    for i, R_i in enumerate(radial):
        a, b, c = eccentricities_to_axes(R_i, ep[i], eq[i])
        H20, H22 = triaxial_to_H2m(a, b, c)
        H_hyd.append({(2, 0): H20, (2, 2): H22})

    # non-hydrostatic surface topography + Airy root at the shell base
    H_nh_t = {lm: H_OBS.get(lm, 0.0)
              - (H_hyd[2].get(lm, 0.0) if lm[0] == 2 else 0.0)
              for lm in H_OBS}
    M_below_shell = M_ENC - m_shell
    g_t = gravity_g(M_ENC, R_T)
    g_b = gravity_g(M_below_shell, R_b)
    H_nh_b = airy_root(H_nh_t, RHO_ICE, RHO_OCEAN, g_t, g_b, C2=1.0)

    # total per-interface topography (H&M Eq. 9), finite-amplitude
    # adjusted, then the Eq. 8 interface sum with below-minus-above
    # density contrasts.
    interfaces = (
        (R_T, RHO_ICE - 0.0,
         {lm: H_hyd[2].get(lm, 0.0) + H_nh_t.get(lm, 0.0)
          for lm in H_OBS}),
        (R_b, RHO_OCEAN - RHO_ICE,
         {lm: H_hyd[1].get(lm, 0.0) + H_nh_b.get(lm, 0.0)
          for lm in H_OBS}),
        (R_c, rho_core - RHO_OCEAN,
         {lm: H_hyd[0].get(lm, 0.0) for lm in H_OBS}),
    )
    C_pred = {lm: 0.0 for lm in H_OBS}
    for R_i, drho, H_lm in interfaces:
        H_adj = finite_amplitude_coeffs(H_lm, R_i, n_theta=48, n_phi=96)
        for (l, m) in C_pred:
            C_pred[(l, m)] += interface_gravity_coeff(
                H_adj[(l, m)], drho, R_i, l, M_ENC, R_REF)

    lib_m = librations(radial, rho, OMEGA, ECC,
                       rigid=True, ocean=True, ocean_idx=1)
    lib_deg = float(np.degrees(lib_m / R_T))
    return C_pred, lib_deg, rho_core


def misfit(D_shell_km, D_ocean_km):
    out = predict(D_shell_km, D_ocean_km)
    if out is None:
        return np.inf
    C_pred, lib_deg, _ = out
    L = sum((C_pred[lm] - C_OBS[lm]) ** 2 / SIG_C[lm] ** 2
            for lm in C_OBS)
    L += (lib_deg - LIB_OBS_DEG) ** 2 / SIG_LIB ** 2
    return L


@pytest.mark.slow
def test_hm2019_reproduction_gate():
    shells = np.arange(8.0, 41.0, 1.0)
    oceans = np.arange(10.0, 61.0, 1.0)
    best = (np.inf, None, None)
    for D_s in shells:
        for D_o in oceans:
            L = misfit(D_s, D_o)
            if L < best[0]:
                best = (L, D_s, D_o)
    L_min, D_s_best, D_o_best = best
    print(f'\nB13 gate: misfit minimum L={L_min:.2f} at '
          f'shell={D_s_best:.0f} km, ocean={D_o_best:.0f} km')
    _, _, rho_core = predict(D_s_best, D_o_best)
    R_c_km = (R_T - (D_s_best + D_o_best) * 1e3) / 1e3
    print(f'   core radius {R_c_km:.0f} km, rho_core {rho_core:.0f}')
    # H&M preferred solution box (their section 4 / abstract)
    assert 19.0 <= D_s_best <= 24.0, \
        f'shell {D_s_best} outside H&M 19-24 km'
    assert 35.0 <= D_o_best <= 39.0, \
        f'ocean {D_o_best} outside H&M 35-39 km'
    assert 192.0 <= R_c_km <= 195.0
    assert 2340.0 <= rho_core <= 2410.0
