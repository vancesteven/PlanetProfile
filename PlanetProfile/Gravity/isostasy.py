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


_RHO_INTERIOR_CEILING_KGM3 = 9000.0


def mass_neutral_shell_density(
    radial: np.ndarray,
    rho: np.ndarray,
    rho_shell_new: float,
    M_body: float,
    shell_idx: int = -1,
    interior_idx: int = 0,
) -> Tuple[np.ndarray, float]:
    """Apply a shell-density override MASS-NEUTRALLY.

    Any EOS/nuisance override of a shell's density changes that layer's
    mass unless something else in the stack compensates. This helper
    absorbs the shell mass delta into a nominated interior layer's density
    so the body's total mass is conserved EXACTLY (by construction, to
    floating-point precision) rather than left to drift with the override.

    PlanetProfile convention: ``rho[i]`` occupies the spherical shell
    ``[radial[i-1], radial[i]]`` (with an implicit inner edge of 0 for
    ``i = 0``); ``radial`` is ascending and ``radial[-1]`` is the surface.

    Args:
        radial: ascending outer radii of each layer (m).
        rho: density of each layer (kg/m^3), same length as ``radial``.
        rho_shell_new: the overridden density (kg/m^3) for the layer at
            ``shell_idx``.
        M_body: total body mass (kg) to conserve exactly.
        shell_idx: index of the layer receiving the override (default:
            outermost layer, -1).
        interior_idx: index of the layer whose density absorbs the shell
            mass delta (default: innermost layer, 0).

    Returns:
        (rho_new, mass_residual_frac): a NEW array (input ``rho`` is not
        mutated) with ``rho_new[shell_idx] = rho_shell_new`` and
        ``rho_new[interior_idx]`` solved so that
        ``sum(rho_new * layer_volume) == M_body`` exactly; and the
        achieved fractional mass residual
        ``(sum(rho_new * layer_volume) - M_body) / M_body`` (zero to
        floating-point precision, by construction).

    Raises:
        ValueError: if ``radial``/``rho`` have fewer than 2 layers, if
            ``shell_idx`` and ``interior_idx`` resolve to the same layer,
            if any layer has non-positive volume, or if the interior
            density required to restore ``M_body`` is unphysical
            (<= 0 or > 9000 kg/m^3).
    """
    radial = np.asarray(radial, dtype=float)
    rho_in = np.asarray(rho, dtype=float)
    n = radial.size
    if n < 2 or rho_in.size != n:
        raise ValueError(
            'radial and rho must be 1D arrays of equal length >= 2, got '
            f'sizes {n} and {rho_in.size}')

    shell_idx = shell_idx % n
    interior_idx = interior_idx % n
    if shell_idx == interior_idx:
        raise ValueError(
            f'shell_idx ({shell_idx}) and interior_idx ({interior_idx}) '
            'must refer to distinct layers')

    inner_edges = np.concatenate(([0.0], radial[:-1]))
    vol = (4.0 / 3.0) * np.pi * (radial ** 3 - inner_edges ** 3)
    if np.any(vol <= 0.0):
        raise ValueError('non-positive layer volume in radial stack '
                          '(radial must be strictly ascending and > 0)')

    # Mass of every layer OTHER than the shell/interior pair is held fixed;
    # the interior density is solved so the total lands on M_body exactly.
    m_fixed = float(np.sum(rho_in * vol)
                    - rho_in[shell_idx] * vol[shell_idx]
                    - rho_in[interior_idx] * vol[interior_idx])
    m_shell_new = rho_shell_new * vol[shell_idx]
    m_interior_new = M_body - m_fixed - m_shell_new
    rho_interior_new = m_interior_new / vol[interior_idx]

    if (not np.isfinite(rho_interior_new) or rho_interior_new <= 0.0
            or rho_interior_new > _RHO_INTERIOR_CEILING_KGM3):
        raise ValueError(
            'mass-neutral shell-density override is unphysical: restoring '
            f'M_body={M_body:.6g} kg with shell density {rho_shell_new:.6g} '
            f'kg/m^3 requires interior density {rho_interior_new:.6g} '
            f'kg/m^3, outside (0, {_RHO_INTERIOR_CEILING_KGM3:.0f}] kg/m^3')

    rho_out = rho_in.copy()
    rho_out[shell_idx] = rho_shell_new
    rho_out[interior_idx] = rho_interior_new

    mass_residual_frac = float(np.sum(rho_out * vol) - M_body) / M_body
    return rho_out, mass_residual_frac


# ---------------------------------------------------------------------------
# Frozen-branch analytic mass closure (frozen-branch design ruling, 3.1)
# ---------------------------------------------------------------------------
#
# On the frozen branch there is no ocean, so at Enceladus pressures the whole
# hydrosphere is ice Ih (basal P ~6.8 MPa is ~30x below PminHPices -- ruling
# F4) and the body reduces to TWO uniform layers: a rock interior of radius
# R_rock and density rho_rock, and an ice shell of thickness zb = R - R_rock
# and density rho_ice. Total mass is then
#
#     M = (4 pi / 3) [ rho_rock R_rock^3 + rho_ice (R^3 - R_rock^3) ]
#
# which is a strictly monotone, ANALYTICALLY INVERTIBLE relation between zb
# and rho_rock at fixed (M, R, rho_ice). That bijection is the frozen branch's
# sampling<->build map (ruling: sample rho_rock, DERIVE zb; no iteration, and
# C/MR^2 is derived ALONG it rather than gating it -- ruling F1/F2).
#
# Both functions here are PURE and analytic: no PlanetProfile call, no solver,
# no global state. ``strict`` selects the failure mode -- ``False`` (default)
# rejects to ``None`` for the builder's HARD-REJECT-to-None discipline,
# ``True`` raises ``ValueError`` for callers that must not silently drop a
# node (tests, the per-sample loop-closure check I-F5).

def _frozen_reject(msg: str, strict: bool):
    """Reject-to-None (non-strict) or raise (strict). Ruling 3.1."""
    if strict:
        raise ValueError(msg)
    return None


def frozen_zb_from_mass(
    M_body_kg: float,
    R_body_m: float,
    rho_rock_kgm3: float,
    rho_ice_kgm3: float,
    strict: bool = False,
) -> Optional[float]:
    """Ice-shell thickness of a FROZEN two-layer body, from mass closure.

    Inverts the two-layer mass relation exactly (no iteration)::

        R_rock^3 = (3 M / (4 pi) - rho_ice R^3) / (rho_rock - rho_ice)
        zb       = R - R_rock

    This is the frozen branch's sampling->build map: ``rho_rock`` is the
    sampled coordinate (declared U[2200, 2600] kg/m^3) and ``zb`` is
    DERIVED, so the frozen support is fixed by the density prior rather
    than by an arbitrary zb box edge (config
    ``branch_model.why_reparameterized``).

    Args:
        M_body_kg: total body mass (kg), conserved EXACTLY.
        R_body_m: surface radius (m).
        rho_rock_kgm3: uniform interior density (kg/m^3).
        rho_ice_kgm3: uniform ice-shell density (kg/m^3).
        strict: ``False`` (default) returns ``None`` on any rejection;
            ``True`` raises ``ValueError`` with the reason.

    Returns:
        Shell thickness zb in METRES, or ``None`` (non-strict) when the
        inputs admit no physical two-layer solution.

    Rejection / raise conditions:
        - any input non-finite, or ``M_body_kg``/``R_body_m`` non-positive;
        - ``rho_ice_kgm3 <= 0``;
        - ``rho_rock_kgm3 <= rho_ice_kgm3`` (no differentiated solution;
          the relation is singular at equality);
        - mean density below ice (``3M/(4 pi) < rho_ice R^3``): the body
          cannot be built from rock + ice at all, ``R_rock^3 < 0``;
        - mean density above rock (``R_rock > R``, i.e. ``zb < 0``).

    Verified against the ruling's F5 feasibility numbers (RESUME_NOTE
    "deep-frozen Enceladus is feasible", M = 1.08022e20 kg, R = 252.1 km):
    rho_ice 900 / 925 / 970 with rho_rock 2314.3 / 2334.3 / 2377.5 gives
    zb = 51.78 / 53.93 / 58.29 km, all at C/MR^2 = 0.3350.
    """
    M = float(M_body_kg)
    R = float(R_body_m)
    rr = float(rho_rock_kgm3)
    ri = float(rho_ice_kgm3)
    if not all(math.isfinite(v) for v in (M, R, rr, ri)):
        return _frozen_reject(
            f'non-finite input: M_body_kg={M}, R_body_m={R}, '
            f'rho_rock_kgm3={rr}, rho_ice_kgm3={ri}', strict)
    if M <= 0.0 or R <= 0.0:
        return _frozen_reject(
            f'M_body_kg ({M}) and R_body_m ({R}) must be positive', strict)
    if ri <= 0.0:
        return _frozen_reject(
            f'rho_ice_kgm3 ({ri}) must be positive', strict)
    if rr <= ri:
        return _frozen_reject(
            f'frozen two-layer closure requires rho_rock_kgm3 ({rr}) > '
            f'rho_ice_kgm3 ({ri}); the mass relation is singular at '
            'equality and unphysical below it', strict)

    M_over_vol = 3.0 * M / (4.0 * math.pi)      # = rho_mean * R^3
    R3 = R ** 3
    R_rock_cubed = (M_over_vol - ri * R3) / (rr - ri)
    if R_rock_cubed < 0.0:
        return _frozen_reject(
            f'mean density {M_over_vol / R3:.6g} kg/m^3 is below the ice '
            f'density {ri:.6g}: no rock/ice two-layer body exists', strict)
    if R_rock_cubed > R3:
        return _frozen_reject(
            f'mean density {M_over_vol / R3:.6g} kg/m^3 exceeds the rock '
            f'density {rr:.6g}: the closure returns a negative shell '
            'thickness', strict)

    zb_m = R - R_rock_cubed ** (1.0 / 3.0)
    if not math.isfinite(zb_m):
        return _frozen_reject('non-finite zb from mass closure', strict)
    return float(zb_m)


def frozen_rho_rock_from_zb(
    M_body_kg: float,
    R_body_m: float,
    zb_m: float,
    rho_ice_kgm3: float,
    strict: bool = False,
) -> Optional[float]:
    """Interior density of a FROZEN two-layer body, from mass closure.

    Exact inverse of :func:`frozen_zb_from_mass` at fixed
    ``(M_body_kg, R_body_m, rho_ice_kgm3)``::

        R_rock   = R - zb
        rho_rock = [3 M / (4 pi) - rho_ice (R^3 - R_rock^3)] / R_rock^3

    Used to read a grid node's zb back into the sampled coordinate
    (cache-build parameterization closure, invariant I-F3/I-F4) and to
    close the per-sample loop I-F5.

    Args:
        M_body_kg: total body mass (kg).
        R_body_m: surface radius (m).
        zb_m: ice-shell thickness (m).
        rho_ice_kgm3: uniform ice-shell density (kg/m^3).
        strict: ``False`` (default) returns ``None`` on rejection;
            ``True`` raises ``ValueError``.

    Returns:
        Interior density in kg/m^3, or ``None`` (non-strict) when the
        inputs admit no physical two-layer solution.

    Rejection / raise conditions:
        - any input non-finite, or ``M_body_kg``/``R_body_m`` non-positive;
        - ``rho_ice_kgm3 <= 0``;
        - ``zb_m`` outside ``(0, R_body_m)`` (a zero-thickness shell leaves
          no shell, a full-radius shell leaves no interior and the
          expression is singular);
        - the implied ``rho_rock <= rho_ice`` (undifferentiated or
          inverted density profile -- the same physical exclusion the
          forward direction applies).
    """
    M = float(M_body_kg)
    R = float(R_body_m)
    zb = float(zb_m)
    ri = float(rho_ice_kgm3)
    if not all(math.isfinite(v) for v in (M, R, zb, ri)):
        return _frozen_reject(
            f'non-finite input: M_body_kg={M}, R_body_m={R}, zb_m={zb}, '
            f'rho_ice_kgm3={ri}', strict)
    if M <= 0.0 or R <= 0.0:
        return _frozen_reject(
            f'M_body_kg ({M}) and R_body_m ({R}) must be positive', strict)
    if ri <= 0.0:
        return _frozen_reject(
            f'rho_ice_kgm3 ({ri}) must be positive', strict)
    if not (0.0 < zb < R):
        return _frozen_reject(
            f'zb_m ({zb}) must lie strictly inside (0, R_body_m={R})',
            strict)

    R_rock = R - zb
    R_rock3 = R_rock ** 3
    rho_rock = (3.0 * M / (4.0 * math.pi)
                - ri * (R ** 3 - R_rock3)) / R_rock3
    if not math.isfinite(rho_rock):
        return _frozen_reject('non-finite rho_rock from mass closure',
                              strict)
    if rho_rock <= ri:
        return _frozen_reject(
            f'zb_m ({zb}) implies rho_rock {rho_rock:.6g} kg/m^3 <= '
            f'rho_ice {ri:.6g} kg/m^3: undifferentiated or inverted, '
            'outside the frozen branch\'s physical domain', strict)
    return float(rho_rock)


def frozen_cmr2(
    M_body_kg: float,
    R_body_m: float,
    zb_m: float,
    rho_rock_kgm3: float,
    rho_ice_kgm3: float,
) -> float:
    """Polar moment factor C/MR^2 of the frozen two-layer body.

    ``C = (8 pi / 15) [rho_rock R_rock^5 + rho_ice (R^5 - R_rock^5)]``
    for uniform spherical layers (the leading-order figure correction is
    below the 1e-3 scale this quantity is used at).

    DIAGNOSTIC ONLY on this branch. Per ruling 3.2 the frozen build carries
    NO MoI gate -- C/MR^2 is recorded and reported, never used to select a
    structure -- so this exists to make the derived value auditable (and to
    let the I-F6 non-conditioning check measure how far the frozen set
    spreads across the Iess window).

    Returns ``nan`` on degenerate inputs rather than raising; callers that
    need a hard reject use the two closure functions above.
    """
    M = float(M_body_kg)
    R = float(R_body_m)
    R_rock = R - float(zb_m)
    if not all(math.isfinite(v) for v in
               (M, R, R_rock, float(rho_rock_kgm3), float(rho_ice_kgm3))):
        return float('nan')
    if M <= 0.0 or R <= 0.0 or not (0.0 <= R_rock <= R):
        return float('nan')
    C = (8.0 * math.pi / 15.0) * (
        float(rho_rock_kgm3) * R_rock ** 5
        + float(rho_ice_kgm3) * (R ** 5 - R_rock ** 5))
    return float(C / (M * R ** 2))
