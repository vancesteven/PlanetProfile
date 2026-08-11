"""Isostatic non-hydrostatic gravity following Hemingway & Mittal (2019),
Icarus 332, 111-131 (papers/hemingway2019enceladus.pdf) and Hemingway &
Matsuyama (2017) equal-pressures Airy isostasy.

Forward model for bodies (Enceladus-class) whose observed long-wavelength
shape is strongly super-hydrostatic and isostatically compensated: the
OBSERVED surface shape is a model INPUT (never a likelihood residual);
its non-hydrostatic part is Airy-compensated at the ice/ocean interface
(H&M Eq. 12), and the resulting degree-2 and degree-3 gravity is
predicted via the interface-topography potential sum (H&M Eqs. 4, 8)
with the Wieczorek & Phillips (1998) finite-amplitude correction
(H&M Eqs. 5-6; ~5% on J2 for Enceladus if omitted).

Core-parity (user directive 2026-08-12): this module lives in core
PlanetProfile and is importable from both the inference runners and a
standard CLI PlanetProfile run (exposure alongside the standard-run
C20/C22 output; see the C17 task). No inference-only physics.

Spherical-harmonic conventions (pinned; see tests/isostasy_test.py):
- Unnormalized real harmonics with Y20 = P20(cos th) = (3 cos^2 th - 1)/2,
  Y22 = 3 sin^2 th cos 2ph, Y30 = P30(cos th) = (5 cos^3 th - 3 cos th)/2.
- Shape: r(th, ph) = R + sum_lm H_lm Y_lm; for a triaxial ellipsoid with
  semi-axes a >= b >= c (a sub-parent, c polar) to first order:
  H20 = (2/3) [c - (a + b)/2],  H22 = (a - b)/6.
  (Checks: Enceladus a-b = 5.2 km -> H22 = 867 m ~ published 857-917 m.)
- Gravity: U(r) = GM/r sum_lm (R_ref/r)^l C_lm Y_lm (H&M Eq. 1);
  an interface at radius R_i with density jump drho_i (BELOW minus
  ABOVE) and topography coefficient H_ilm contributes
  C_lm = 4 pi drho_i H_ilm R_i^(l+2) / [M (2l+1) R_ref^l]   (H&M Eq. 8).
- Airy root (H&M Eq. 12, equal pressures), with C2 the sampled
  compensation fraction (C2 = 1 pure Airy):
  H_b_lm = -C2 * H_t_lm * (rho_ice / drho_b) * (g_t / g_b),
  drho_b = rho_ocean - rho_ice > 0. The minus sign: topography is
  upward-positive on every interface; a surface high demands a
  DEPRESSED shell base (root protruding into the ocean). With the
  gravity rule above the two contributions combine to
  C_lm^nh ~ rho_ice H_t [R_t^(l+2) - C2 (g_t/g_b) R_b^(l+2)] / (...)
  whose near-cancellation carries the shell-thickness signal
  (~4.2%/km at l=2 for Enceladus).
Frozen branch (no ice/ocean interface): rigid support — surface term
only (no root; C2 irrelevant). This is a real physical prediction that
lets a joint no-ocean+ocean posterior discriminate the branches.
"""
from __future__ import annotations

import math
from typing import Dict, Optional, Tuple

import numpy as np

G = 6.674e-11

# Supported (l, m) channels and unnormalized basis normalization
# N_lm = integral Y_lm^2 dOmega.
_SUPPORTED_LM = ((2, 0), (2, 2), (3, 0))
_NORMS = {(2, 0): 4.0 * np.pi / 5.0,
          (2, 2): 48.0 * np.pi / 5.0,
          (3, 0): 4.0 * np.pi / 7.0}


def triaxial_to_H2m(a: float, b: float, c: float) -> Tuple[float, float]:
    """Unnormalized degree-2 shape coefficients of a triaxial ellipsoid
    (first order in flattening): H20 = (2/3)[c - (a+b)/2], H22 = (a-b)/6."""
    return (2.0 / 3.0) * (c - 0.5 * (a + b)), (a - b) / 6.0


def eccentricities_to_axes(R_mean: float, ep: float, eq: float
                           ) -> Tuple[float, float, float]:
    """Semi-axes (a, b, c) from Tricarico polar/equatorial eccentricities
    at volumetric mean radius R_mean (abc = R^3):
    ep^2 = (a^2 - c^2)/a^2, eq^2 = (a^2 - b^2)/a^2."""
    fp = 1.0 - ep ** 2
    fq = 1.0 - eq ** 2
    a = R_mean / (fp * fq) ** (1.0 / 6.0)
    return a, a * np.sqrt(fq), a * np.sqrt(fp)


def _basis_fields(theta: np.ndarray, phi: np.ndarray) -> Dict[tuple,
                                                              np.ndarray]:
    ct = np.cos(theta)
    st = np.sin(theta)
    return {
        (2, 0): 0.5 * (3.0 * ct ** 2 - 1.0),
        (2, 2): 3.0 * st ** 2 * np.cos(2.0 * phi),
        (3, 0): 0.5 * (5.0 * ct ** 3 - 3.0 * ct),
    }


def _grid(n_theta: int = 64, n_phi: int = 128):
    """Gauss-Legendre in cos(theta) x uniform phi quadrature."""
    x, wx = np.polynomial.legendre.leggauss(n_theta)
    theta = np.arccos(x)
    phi = np.linspace(0.0, 2.0 * np.pi, n_phi, endpoint=False)
    wphi = 2.0 * np.pi / n_phi
    TH, PH = np.meshgrid(theta, phi, indexing='ij')
    W = wx[:, None] * wphi * np.ones((1, n_phi))
    return TH, PH, W


def project_powers(H_lm: Dict[tuple, float], n_max: int,
                   n_theta: int = 64, n_phi: int = 128
                   ) -> Dict[int, Dict[tuple, float]]:
    """Spherical-harmonic coefficients (supported channels) of H^n for
    n = 1..n_max, where H(th,ph) = sum_lm H_lm Y_lm (H&M Eq. 5)."""
    TH, PH, W = _grid(n_theta, n_phi)
    basis = _basis_fields(TH, PH)
    H = np.zeros_like(TH)
    for lm, coeff in H_lm.items():
        if lm not in _NORMS:
            raise ValueError(f'unsupported (l,m)={lm}')
        H = H + coeff * basis[lm]
    out: Dict[int, Dict[tuple, float]] = {}
    Hn = np.ones_like(H)
    for n in range(1, n_max + 1):
        Hn = Hn * H
        out[n] = {lm: float(np.sum(Hn * basis[lm] * W) / _NORMS[lm])
                  for lm in _NORMS}
    return out


def finite_amplitude_coeffs(H_lm: Dict[tuple, float], R_i: float,
                            n_theta: int = 64, n_phi: int = 128
                            ) -> Dict[tuple, float]:
    """Finite-amplitude-adjusted topography coefficients H+_ilm for
    evaluation at r >= R_i (H&M Eq. 6 / Wieczorek & Phillips 1998):

      H+_ilm = sum_{n=1}^{l+3} H(n)_ilm / R_i^(n-1)
               * [prod_{j=1}^{n} (l + 4 - j)] / [n! (l + 3)]

    The n = 1 term is the raw topography coefficient."""
    l_max_n = max(l for (l, m) in _NORMS) + 3
    powers = project_powers(H_lm, l_max_n, n_theta, n_phi)
    out = {}
    for (l, m) in _NORMS:
        total = 0.0
        for n in range(1, l + 3 + 1):
            num = 1.0
            for j in range(1, n + 1):
                num *= (l + 4 - j)
            fact = float(math.factorial(n))
            term = powers[n][(l, m)] / R_i ** (n - 1) * num / (fact
                                                               * (l + 3))
            total += term
        out[(l, m)] = total
    return out


def interface_gravity_coeff(H_ilm: float, drho: float, R_i: float,
                            l: int, M_body: float, R_ref: float) -> float:
    """C_lm contribution of one interface (H&M Eq. 8 evaluated at R_ref):
    C_lm = 4 pi drho H_ilm R_i^(l+2) / [M (2l+1) R_ref^l].
    drho = rho(below) - rho(above); H upward-positive."""
    return (4.0 * np.pi * drho * H_ilm * R_i ** (l + 2)
            / (M_body * (2 * l + 1) * R_ref ** l))


def airy_root(H_t_lm: Dict[tuple, float], rho_ice: float,
              rho_ocean: float, g_t: float, g_b: float,
              C2: float = 1.0) -> Dict[tuple, float]:
    """Equal-pressures Airy basal relief (H&M Eq. 12) scaled by the
    compensation fraction C2 in [0, 1]."""
    drho = rho_ocean - rho_ice
    if drho <= 0:
        raise ValueError(f'rho_ocean {rho_ocean} must exceed rho_ice '
                         f'{rho_ice} for a floating shell')
    fac = -C2 * (rho_ice / drho) * (g_t / g_b)
    return {lm: fac * v for lm, v in H_t_lm.items()}


def gravity_g(M_enclosed: float, R: float) -> float:
    return G * M_enclosed / R ** 2


def isostatic_gravity(
    H_obs_lm: Dict[tuple, float],
    H_hyd_t_lm: Dict[tuple, float],
    R_t: float,
    R_b: Optional[float],
    rho_ice: float,
    rho_ocean: Optional[float],
    M_body: float,
    R_ref: float,
    C2: float = 1.0,
    finite_amplitude: bool = True,
    n_theta: int = 64,
    n_phi: int = 128,
) -> Dict[tuple, float]:
    """NON-hydrostatic gravity coefficients from the observed shape under
    Airy isostasy (ocean branch) or rigid support (frozen branch).

    Args:
        H_obs_lm: observed surface shape coefficients at (their) R_ref
            convention, metres. Degree-3 entries are wholly
            non-hydrostatic (H_hyd has no l=3 power).
        H_hyd_t_lm: hydrostatic surface figure coefficients (from the
            Tricarico interface figures), metres.
        R_t: mean radius of the shell top (m).
        R_b: mean radius of the shell base (m); None => FROZEN branch
            (rigid support, surface term only).
        rho_ice, rho_ocean: shell / ocean densities (kg/m^3);
            rho_ocean may be None on the frozen branch.
        M_body, R_ref: body mass (kg) and gravity reference radius (m).
        C2: Airy compensation fraction (sampled nuisance).
        finite_amplitude: apply H&M Eq. 6 (required for production;
            switchable for tests).

    Returns {(l, m): C_lm^nh} for (2,0), (2,2), (3,0).
    """
    H_nh_t = {}
    for lm in _SUPPORTED_LM:
        obs = H_obs_lm.get(lm, 0.0)
        hyd = H_hyd_t_lm.get(lm, 0.0) if lm[0] == 2 else 0.0
        H_nh_t[lm] = obs - hyd

    if finite_amplitude:
        Ht_eff = finite_amplitude_coeffs(H_nh_t, R_t, n_theta, n_phi)
    else:
        Ht_eff = dict(H_nh_t)

    out = {}
    for (l, m) in _SUPPORTED_LM:
        # surface: below = ice, above = vacuum -> drho = +rho_ice
        c = interface_gravity_coeff(Ht_eff[(l, m)], rho_ice, R_t, l,
                                    M_body, R_ref)
        out[(l, m)] = c

    if R_b is not None:
        if rho_ocean is None:
            raise ValueError('ocean branch requires rho_ocean')
        M_below_shell = M_body - rho_ice * (4.0 * np.pi / 3.0) \
            * (R_t ** 3 - R_b ** 3)
        g_t = gravity_g(M_body, R_t)
        g_b = gravity_g(M_below_shell, R_b)
        H_b = airy_root(H_nh_t, rho_ice, rho_ocean, g_t, g_b, C2)
        if finite_amplitude:
            Hb_eff = finite_amplitude_coeffs(H_b, R_b, n_theta, n_phi)
        else:
            Hb_eff = H_b
        for (l, m) in _SUPPORTED_LM:
            # base: below = ocean, above = ice -> drho = +(rho_ocean-rho_ice)
            out[(l, m)] += interface_gravity_coeff(
                Hb_eff[(l, m)], rho_ocean - rho_ice, R_b, l,
                M_body, R_ref)
    return out


def rescale_coeff(C_lm: float, l: int, R_ref_from: float,
                  R_ref_to: float) -> float:
    """Reference-radius conversion (B14): C_lm(R_to) =
    C_lm(R_from) * (R_from / R_to)^l  [from U ~ (R_ref/r)^l C_lm]."""
    return C_lm * (R_ref_from / R_ref_to) ** l
