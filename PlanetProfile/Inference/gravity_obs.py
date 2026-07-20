"""Hydrostatic degree-2 gravity observables for inference (Clipper v4).

Forward map from a radial density profile to the hydrostatic gravity
pair (C20_h, C22_h) for a synchronously rotating satellite, via direct
integration of Clairaut's equation in Radau form — NOT the Radau–Darwin
closed form, whose ~1% error in k_f (~1e-6 in C22) exceeds the
projected Europa Clipper sigma(C22) ~ 2e-7 by ~5x and would alias into
the sampled non-hydrostaticity offsets (scientific review 2026-07-19,
plans/europa-clipper-v4-geodesy-plan.md). Radau–Darwin is retained here
only as a unit-test cross-check and as the GUI's quick inverse readout.

Conventions (binding, from the v4 plan):
- UNNORMALIZED Stokes coefficients (Anderson et al. 1998 / Gomez
  Casajus et al. 2021 convention). Fully-normalized conversions:
  C20_unnorm = sqrt(5) * C̄20, C22_unnorm = sqrt(10/24) * C̄22 ≈
  0.6455 * C̄22.
- Coefficients are quoted at the GC21 reference radius (1565 km), not
  the physical radius (1560.8 km): Clm scales as (R/R_ref)^2, a ~0.5%
  effect on C22 = 3x sigma(C22) — load-bearing.
- J2_h/C22_h = 3.324 (Tricarico 2014 rapid-rotation correction, as
  applied by GC21), not the classical synchronous 10/3.
- q_r = omega^2 R^3 / (G M) uses the SPIN rate omega (= mean motion
  for synchronous Europa) and the physical radius R.
"""
from __future__ import annotations

import numpy as np

# Gravitational constant [m^3 kg^-1 s^-2]; matches PP's Constants.G to
# the precision relevant here.
_G = 6.674e-11

# Hydrostatic J2/C22 for a synchronous satellite with Europa's rotation
# rate: classical value 10/3, modified to 3.324 by the finite-rotation
# correction (Tricarico 2014), as adopted by Gomez Casajus et al. 2021.
J2_OVER_C22 = 3.324

# GC21 quote their unnormalized coefficients at this reference radius.
R_REF_GC21_M = 1565.0e3


def rotation_parameter(omega: float, R_m: float, M_kg: float) -> float:
    """q_r = omega^2 R^3 / (G M) (dimensionless; Europa ~4.98e-4)."""
    return float(omega) ** 2 * float(R_m) ** 3 / (_G * float(M_kg))


def _downsample_profile(r: np.ndarray, d: np.ndarray, max_shells: int,
                        jump_rel: float = 0.01):
    """Merge a fine layered profile into fewer mass-preserving shells
    WITHOUT smearing density interfaces: the profile is first split at
    every relative density jump > jump_rel (core/mantle/ocean/ice
    boundaries are kept exact — binning across them corrupts k_f at the
    percent level), then each smooth run is binned to its share of
    max_shells by cumulative volume with volume-weighted (mass-exact)
    densities."""
    if r.size <= max_shells:
        return r, d
    jumps = np.where(np.abs(np.diff(d))
                     > jump_rel * np.maximum(d[:-1], 1e-30))[0]
    seg_bounds = np.concatenate(([0], jumps + 1, [r.size]))
    r_edges_all = np.concatenate(([0.0], r))
    out_r, out_d = [], []
    n_seg = len(seg_bounds) - 1
    for s in range(n_seg):
        i0, i1 = seg_bounds[s], seg_bounds[s + 1]
        seg_r = r[i0:i1]
        seg_d = d[i0:i1]
        lo_edge = r_edges_all[i0]
        vol = seg_r ** 3 - np.concatenate(([lo_edge ** 3],
                                           seg_r[:-1] ** 3))
        n_take = max(1, int(round(max_shells * (i1 - i0) / r.size)))
        if seg_r.size <= n_take:
            out_r.extend(seg_r)
            out_d.extend(seg_d)
            continue
        cumvol = np.cumsum(vol)
        targets = np.linspace(0, cumvol[-1], n_take + 1)[1:]
        idx = np.unique(np.clip(np.searchsorted(cumvol, targets - 1e-12),
                                0, seg_r.size - 1))
        start = 0
        for i in idx:
            v = vol[start:i + 1].sum()
            m = (seg_d[start:i + 1] * vol[start:i + 1]).sum()
            out_r.append(seg_r[i])
            out_d.append(m / v if v > 0 else seg_d[i])
            start = i + 1
    return np.asarray(out_r), np.asarray(out_d)


def clairaut_kf(r_m: np.ndarray, rho: np.ndarray, n_sub: int = 16,
                max_shells: int = 64) -> float:
    """Fluid Love number k_f by integrating Clairaut's equation (Radau
    form) over a layered density profile.

    Radau form: r * deta/dr = -eta^2 + eta + 6 - 6*D(r)*(1 + eta),
    D(r) = rho(r) / rhobar(r), rhobar(r) = 3 m(r) / (4 pi r^3), with
    eta(0) = 0; then k_f = (3 - eta_s) / (2 + eta_s) at the surface.
    Limits: homogeneous body eta_s = 0 -> k_f = 3/2; point-mass (Roche)
    eta_s = 3 -> k_f = 0.

    Parameters
    ----------
    r_m, rho : layer arrays (any radial order; sorted ascending here).
        rho is treated as piecewise-constant per layer (the PP cache
        representation). Non-finite/nonpositive-radius entries dropped.
    n_sub : RK4 substeps per layer (integration error << the ~1e-3
        relative level that matters at sigma(C22)/C22 ~ 1.5e-3).
    """
    r = np.asarray(r_m, dtype=float)
    d = np.asarray(rho, dtype=float)
    keep = np.isfinite(r) & np.isfinite(d) & (r > 0)
    r, d = r[keep], d[keep]
    if r.size < 2:
        raise ValueError("clairaut_kf needs >= 2 profile points")
    order = np.argsort(r)
    r, d = r[order], d[order]
    # Performance: per-sample inference calls integrate this in pure
    # Python; downsampling a ~200-layer cache profile to <= max_shells
    # mass-preserving shells cuts the cost ~10x at verified sub-1e-4
    # relative k_f error (tests).
    r, d = _downsample_profile(r, d, max_shells)

    # Piecewise-constant density: layer i spans (r_edges[i], r_edges[i+1])
    # with density d[i]; innermost layer extends to r = 0.
    r_edges = np.concatenate(([0.0], r))

    # Cumulative mass at the layer edges (analytic within each shell).
    shell_m = (4.0 / 3.0) * np.pi * d * (r_edges[1:] ** 3 - r_edges[:-1] ** 3)
    m_edges = np.concatenate(([0.0], np.cumsum(shell_m)))

    def mass_at(rr: float, i_layer: int) -> float:
        return m_edges[i_layer] + (4.0 / 3.0) * np.pi * d[i_layer] * (
            rr ** 3 - r_edges[i_layer] ** 3)

    def deta_dr(rr: float, eta: float, i_layer: int) -> float:
        m = mass_at(rr, i_layer)
        if m <= 0.0 or rr <= 0.0:
            return 0.0
        rhobar = 3.0 * m / (4.0 * np.pi * rr ** 3)
        D = d[i_layer] / rhobar
        return (-eta * eta + eta + 6.0 - 6.0 * D * (1.0 + eta)) / rr

    # eta(0) = 0; start the integration a hair inside the innermost
    # layer (the RHS -> 0 smoothly as r -> 0 for a regular center).
    eta = 0.0
    for i in range(d.size):
        a = max(r_edges[i], 1e-6 * r[-1])
        b = r_edges[i + 1]
        if b <= a:
            continue
        # Adaptive per-shell retry: at extreme density contrasts eta
        # rises steeply toward its Roche bound (3) and a coarse RK4 step
        # can overshoot into the unstable quadratic branch (NaN). Redo
        # the shell with 8x substeps (up to 3 refinements) when the
        # result is non-finite or unphysical.
        eta_in = eta
        n_i = n_sub
        for _attempt in range(4):
            eta = eta_in
            h = (b - a) / n_i
            rr = a
            for _ in range(n_i):
                k1 = deta_dr(rr, eta, i)
                k2 = deta_dr(rr + 0.5 * h, eta + 0.5 * h * k1, i)
                k3 = deta_dr(rr + 0.5 * h, eta + 0.5 * h * k2, i)
                k4 = deta_dr(rr + h, eta + h * k3, i)
                eta += (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
                rr += h
            if np.isfinite(eta) and -1.0 < eta < 3.5:
                break
            n_i *= 8

    return (3.0 - eta) / (2.0 + eta)


def radau_darwin_kf(cmr2: float) -> float:
    """Closed-form Radau–Darwin k_f from C/MR^2 (cross-check + GUI
    inverse only — NOT the v4 forward map; ~1% systematic)."""
    y = 2.5 * (1.0 - 1.5 * float(cmr2))
    if y < 0:
        raise ValueError(f"C/MR^2 = {cmr2} > 2/3 is unphysical")
    y2 = y * y
    return (4.0 - y2) / (1.0 + y2)


def radau_darwin_cmr2(kf: float) -> float:
    """Inverse of radau_darwin_kf: C/MR^2 from a fluid Love number."""
    kf = float(kf)
    y = np.sqrt((4.0 - kf) / (1.0 + kf))
    return (2.0 / 3.0) * (1.0 - 0.4 * y)


def hydrostatic_c20_c22(kf: float, omega: float, R_m: float, M_kg: float,
                        R_ref_m: float = R_REF_GC21_M,
                        j2_over_c22: float = J2_OVER_C22) -> tuple:
    """Hydrostatic (C20_h, C22_h), UNNORMALIZED, at reference radius
    R_ref_m (default: GC21's Europa 1565 km). C22_h = k_f q_r / 4
    evaluated at the physical radius, rescaled by (R/R_ref)^2;
    J2_h = j2_over_c22 * C22_h; C20_h = -J2_h.

    j2_over_c22 defaults to Europa's Tricarico-corrected 3.324. The
    correction from the classical 10/3 grows with q_r (mission-body
    constraints doc, review-binding): pass the body-specific value via
    the config keys gravity_j2_over_c22 / gravity_ref_radius_m for any
    non-Europa body (e.g. Enceladus q_r ~ 6.3e-3 -> ~3.22)."""
    q_r = rotation_parameter(omega, R_m, M_kg)
    scale = (float(R_m) / float(R_ref_m)) ** 2
    c22 = 0.25 * float(kf) * q_r * scale
    return (-float(j2_over_c22) * c22, c22)


def cmr2_from_c22_rd(c22: float, omega: float, R_m: float, M_kg: float,
                     R_ref_m: float = R_REF_GC21_M) -> float:
    """GUI readout: the hydrostatic C/MR^2 an observed (unnormalized,
    R_ref-referenced) C22 implies under Radau–Darwin — the classical
    Anderson-style data reduction, shown alongside the posterior."""
    q_r = rotation_parameter(omega, R_m, M_kg)
    scale = (float(R_m) / float(R_ref_m)) ** 2
    kf = 4.0 * float(c22) / (q_r * scale)
    return radau_darwin_cmr2(kf)
