"""Libration amplitude helpers for PlanetProfile gravity calculations."""

import logging

import numpy as np

from PlanetProfile.Gravity.isostasy import airy_root, gravity_g
from PlanetProfile.Utilities.defineStructs import Constants

log = logging.getLogger('PlanetProfile')

G = Constants.G

def radial2abc(radial, ep, eq):

    factor = ((np.sqrt(1 - ep**2)) * (np.sqrt(1 - eq**2)))**(1/3)

    a = radial/factor
    b = a * (1 - eq**2) ** 0.5
    c = a * (1 - ep**2) ** 0.5

    return a, b, c


def compute_eccentricities(radial, rho, omega, max_iter=1000):
    # Decouple from caller so it does not modify the original arrays
    radial = np.array(radial, dtype=float, copy=True)
    rho    = np.array(rho,    dtype=float, copy=True)

    # Check for layers with zero thickness, the algorithm will not converge
    # and this case must be handled
    zero_idx = np.isclose(np.diff(radial), 0.0, atol=1e-12)
    zero_check = np.any(zero_idx)

    if zero_check:
        idx = np.where(zero_idx)[0][0]
        radial[idx + 1:] = radial[idx + 1:] + 1.0e-6

    # Density and radial profiles must be inwards for this formulation
    if np.all(np.diff(radial)>0):
        rho = rho[::-1]
        radial = radial[::-1]

    # Define arrays
    n_layers = len(radial)
    ep_sq = np.zeros(n_layers)
    eq_sq = np.zeros(n_layers)

    lamd = np.sqrt(omega**2 / (np.pi * G * rho[0]))

    for iteration in range(max_iter):
        sigmas = np.hstack([1, (rho[1:] - rho[:-1]) / rho[0]])
        mus = radial / radial[0]

        ep_sq_prev = ep_sq.copy()
        eq_sq_prev = eq_sq.copy()

        for i in range(n_layers):
            # Compute polar eccentricity for each layer
            num = 60 * lamd**2 + 12 * (sigmas[:i] * ep_sq[:i]).sum() + \
                12 * ((mus[i+1:]/mus[i]) ** 5 * sigmas[i+1:] * ep_sq[i+1:]).sum()

            den = 20 * sigmas[:i].sum() + 8 * sigmas[i] \
                + 20 * ((mus[i+1:]/mus[i]) ** 3 * sigmas[i+1:]).sum()

            ep_sq[i] = num/den

            # Compute equatorial eccentricity for each layer
            num = 45 * lamd**2 + 12 * (sigmas[:i] * eq_sq[:i]).sum() + \
                (12 * (mus[i+1:]/mus[i]) ** 5 * sigmas[i+1:] * eq_sq[i+1:]).sum()

            den = 20 * sigmas[:i].sum() + 8 * sigmas[i] \
                + 20 * ((mus[i+1:]/mus[i]) ** 3 * sigmas[i+1:]).sum()

            eq_sq[i] = num/den

        if (np.abs(ep_sq[0] - ep_sq_prev[0]) < 1e-15) and \
            (np.abs(eq_sq[0] - eq_sq_prev[0]) < 1e-15):
                # print(f"Converged after {iteration+1} iterations.")
                break

    else:
        log.warning("Libration eccentricity calculation did not converge within %s iterations.", max_iter)

    return np.sqrt(ep_sq)[::-1], np.sqrt(eq_sq)[::-1]

def moi_ellps(a, b, c, rho):
    abc = a * b * c
    absq = a**2 + b**2
    bcsq = b**2 + c**2
    acsq = a**2 + c**2

    Ai = np.hstack([
        1/5 * 4/3*np.pi * rho[0] * abc[0] * bcsq[0],
        1/5 * 4/3*np.pi * rho[1:] * (abc[1:] * bcsq[1:] - abc[:-1] * bcsq[:-1]),
    ])

    Bi = np.hstack([
        1/5 * 4/3*np.pi * rho[0] * abc[0] * acsq[0],
        1/5 * 4/3*np.pi * rho[1:] * (abc[1:] * acsq[1:] - abc[:-1] * acsq[:-1]),
    ])

    Ci = np.hstack([
        1/5 * 4/3*np.pi * rho[0] * abc[0] * absq[0],
        1/5 * 4/3*np.pi * rho[1:] * (abc[1:] * absq[1:] - abc[:-1] * absq[:-1]),
    ])

    return Ai, Bi, Ci

def librations(radial, rho, omega, ecc, rigid=True, ocean=True,
                       ocean_idx=1, y1=None, rek2=None, H22_obs_m=None,
                       compensation_C2=None):
    """Compute the forced libration amplitude (metres, at the surface).

    Figure convention (ocean=True, rigid=True branch only)
    -------------------------------------------------------
    PlanetProfile's historical (default) behaviour drives every layer's
    (B - A) contribution to the libration torques -- including the outer
    shell's -- from the HYDROSTATIC Tricarico figure computed internally
    by ``compute_eccentricities`` / ``radial2abc`` / ``moi_ellps``. This is
    the PP-native convention and is exactly what runs when ``H22_obs_m`` is
    left at its default of ``None``.

    Hemingway & Mittal (2019) instead set the shell's OUTER figure to the
    OBSERVED (non-hydrostatic) surface triaxial figure and obtain the
    ice/ocean interface figure as the Airy-isostatic response to that
    surface load (their Eq. 12, section 2.3), scaled by the compensation
    fraction C2. Passing ``H22_obs_m`` -- with ``compensation_C2`` for the
    isostatic root -- opts into that treatment.

    Delta_rho weighting: the B2' defect repair (user ruling 2026-08-15)
    -------------------------------------------------------------------
    The shell torque telescopes exactly into H&M's Eq. 19 density-contrast
    sum (verified to ~1e-14; standing invariant in
    tests/libration_density_contrast_test.py)::

        Ks / (3 omega^2) = rho_ice * f_top + (rho_ocean - rho_ice) * f_base
        f_i = (4 pi / 15) a_i b_i c_i (a_i^2 - b_i^2)

    so the shell-BASE interface's physical weight is the density contrast
    Delta_rho = rho_ocean - rho_ice ~ 80-95 kg/m^3, even though in PP's
    layer form it appears as a difference of two ~10x larger cancelling
    terms: ``-rho_ice * f_base`` inside ``(Bs - As)`` and
    ``+rho_ocean * f_base`` inside ``Bsp_Asp``.

    The 2026-08-13 implementation scaled the whole of ``(Bs - As)``, which
    touched only the first of those two halves and therefore applied an
    implicit, structure-dependent, SIGN-FLIPPING effective scale of +0.33
    to -0.58 to the base interface while claiming it stayed hydrostatic.
    That treatment is RESCINDED. Every figure scale here is applied to the
    interface's Delta_rho-weighted term -- to ``f_base`` wherever it
    appears, INCLUDING ``Bsp_Asp`` -- never to a moment difference.
    See validation_reports/enceladus_isostasy/
    b2prime_ADJUDICATED_drho_weighting.json and
    plans/MACHINE-B-HANDOFF.md section 0.23.

    Scope, stated precisely: the surface and the shell base are the two
    interfaces the treatment touches. Deep-interior interfaces (below the
    ocean) stay hydrostatic, as do any interior interfaces WITHIN a
    multi-layer shell. K_int is deliberately untouched: it carries the
    separately-adjudicated 8*pi/15 normalization defect whose repair is a
    post-campaign task that must not land piecemeal (see the strict-xfail
    tripwire in tests/librations_test.py).

    At the Enceladus fiducial (zb = 25 km, D_ocean = 36 km,
    rho_ocean = 1005, H22_obs = 857 m) the treatment spans +3.25 sigma_obs
    at C2 = 0 to +1.66 sigma_obs at C2 = 1 relative to hydrostatic, and the
    shell thickness matching the Park libration 0.092 deg is bracketed by
    C2 between 27.34 km (C2 = 0) and 25.99 km (C2 = 1).

    Parameters
    ----------
    H22_obs_m : float, optional
        Observed surface H22 in metres, same convention as
        ``PlanetProfile.Gravity.isostasy.triaxial_to_H2m``, i.e.
        H22 = (a - b) / 6. When ``None`` (default), the hydrostatic figure
        is used everywhere and behaviour is identical to leaving this
        parameter unset entirely. Only supported for the
        ``ocean=True, rigid=True`` branch; passing a non-default value for
        any other branch raises ``NotImplementedError``.
    compensation_C2 : float, optional
        Airy compensation fraction (H&M Eq. 12; C2 = 1 is pure Airy).
        Requires ``H22_obs_m``. When ``None`` (or 0) the shell base keeps
        its hydrostatic figure, which is the surface-only treatment and is
        the exact C2 -> 0 limit of the Eq.-12 treatment.
    """
    figure_option_used = H22_obs_m is not None or compensation_C2 is not None
    if figure_option_used and not (ocean and rigid):
        raise NotImplementedError(
            "The observed-figure libration convention (H22_obs_m / "
            "compensation_C2) is only implemented for the ocean=True, "
            "rigid=True (Van Hoolst et al. 2008 rigid) branch. Got "
            "ocean=%r, rigid=%r." % (ocean, rigid)
        )
    if compensation_C2 is not None and H22_obs_m is None:
        raise ValueError(
            "compensation_C2 requires H22_obs_m: the Airy root (H&M "
            "Eq. 12) responds to the NON-HYDROSTATIC part of the observed "
            "surface figure, which is identically zero when the surface "
            "figure is left hydrostatic."
        )
    if compensation_C2 is not None and len(rho) - ocean_idx - 1 != 1:
        raise NotImplementedError(
            "The Eq.-12 isostatic root is defined against a single shell "
            "density; got %d layers above the ocean. Reduce the stack "
            "first (see Gravity.LibrationModelInputs)."
            % (len(rho) - ocean_idx - 1)
        )

    # Compute eccentricities of the layers, semi-major axes, and moments of
    # inertia from Tricarico (2014) formulation
    ep, eq = compute_eccentricities(radial, rho, omega)
    a, b, c = radial2abc(radial, ep, eq)
    Ai, Bi, Ci = moi_ellps(a, b, c, rho)

    # Define auxiliary variables
    r_body = radial[-1]
    fpi = 4/3 * np.pi

    m_body  =  fpi * rho[0] * radial[0]**3
    m_body += (fpi *(rho[1:] * (radial[1:]**3 - radial[0:-1]**3))).sum()
    r_body  = radial[-1]

    q = omega ** 2 * r_body ** 3 / (G * m_body)

    # Oceanless body, compute solid libration amplitude
    if not ocean:
        # Infinintely rigid body libration amplitude
        if rigid:
            A = Ai.sum()
            B = Bi.sum()
            C = Ci.sum()

            # Compute libration amplitude
            gamma = 6 * ecc * (B - A) / C

            gs = abs(gamma * r_body)
        # Oceanless body deformed by tidal deformations (Van Hoolst et al. 2013)
        else:
            if rek2 is None:
                # Real part of the second degree tidal Love number
                raise ValueError("Re(k2) must be provided for oceanless non-rigid bodies.")

            # Body moments of inertia
            A = Ai.sum()
            B = Bi.sum()
            C = Ci.sum()

            # Compute effects of zonal and centrifugal deformations
            delta_c = 4 * rek2 * omega**2 * r_body**5 / 9 / G
            kf = (B - A) / (q * m_body * r_body ** 2)
            fact_ab = (kf - 5/6 * rek2) / kf

            # Computed modified moments of inertia
            Atil = A * fact_ab
            Btil = B * fact_ab
            Ctil = C + delta_c

            gamma = 6 * ecc * (Btil - Atil) / Ctil

            gs = abs(gamma * r_body)

    # Ocean present, compute libration amplitude with ocean effects
    else:
        # Infinintely rigid body libration amplitude (Van Hoolst et al. 2008)
        if rigid:
            # Outer shell moments of inertia
            As = Ai[ocean_idx+1:].sum()
            Bs = Bi[ocean_idx+1:].sum()
            Cs = Ci[ocean_idx+1:].sum()

            # Deep interior moments of inertia
            Ac = Ai[:ocean_idx].sum()
            Bc = Bi[:ocean_idx].sum()
            Cc = Ci[:ocean_idx].sum()

            alpha = 1 - 2 * c / (b + a)
            beta = (a - b) / a

            # Effects of ocean pressure
            Bip_Aip = 8/15*np.pi * rho[ocean_idx] * beta[0] * radial[0]**5
            Bsp_Asp = 8/15*np.pi * rho[ocean_idx] * beta[ocean_idx] * radial[ocean_idx]**5

            # Compute torques separately acting on outer shell and interior
            if H22_obs_m is None:
                # PP-native hydrostatic figures. Byte-identical to the
                # historical code path.
                Ks = 3 * omega**2 * ((Bs - As) + Bsp_Asp)
            else:
                # Delta_rho-consistent figure treatment (B2' repair; see the
                # docstring). Split (Bs - As) into the two interfaces that
                # carry it, so a figure rescaling acts on each interface's
                # physical density-contrast weight rather than on the
                # near-cancelling difference.
                f_fig = 4/15*np.pi * a * b * c * (a**2 - b**2)
                BmA_surface = rho[-1] * f_fig[-1]
                BmA_base = -rho[ocean_idx+1] * f_fig[ocean_idx]
                # Interfaces inside a multi-layer shell stay hydrostatic;
                # this residual is identically zero for a reduced stack.
                BmA_shell_interior = (Bs - As) - BmA_surface - BmA_base

                H22_hyd_surface_m = (a[-1] - b[-1]) / 6
                surface_figure_scale = H22_obs_m / H22_hyd_surface_m

                if compensation_C2 is None:
                    # Surface-only treatment == the C2 -> 0 Eq.-12 limit.
                    base_figure_scale = 1.0
                    H22_root_m = 0.0
                else:
                    shell_mass = (fpi * rho[ocean_idx+1:]
                                  * (radial[ocean_idx+1:]**3
                                     - radial[ocean_idx:-1]**3)).sum()
                    g_t = gravity_g(m_body, radial[-1])
                    g_b = gravity_g(m_body - shell_mass, radial[ocean_idx])
                    H22_hyd_base_m = (a[ocean_idx] - b[ocean_idx]) / 6
                    H22_root_m = airy_root(
                        {(2, 2): H22_obs_m - H22_hyd_surface_m},
                        rho[ocean_idx+1], rho[ocean_idx], g_t, g_b,
                        compensation_C2,
                    )[(2, 2)]
                    base_figure_scale = ((H22_hyd_base_m + H22_root_m)
                                         / H22_hyd_base_m)

                log.debug(
                    "librations(): H&M Delta_rho-consistent figure "
                    "treatment active -- H22_hyd_surface=%.6g m, "
                    "H22_observed=%.6g m, surface_figure_scale=%.6g, "
                    "C2=%r, H22_airy_root=%.6g m, base_figure_scale=%.6g",
                    H22_hyd_surface_m, H22_obs_m, surface_figure_scale,
                    compensation_C2, H22_root_m, base_figure_scale,
                )

                # The base scale multiplies BOTH halves of the base
                # interface term, so the net perturbation enters at the
                # physical Delta_rho weight.
                Ks = 3 * omega**2 * (
                    BmA_shell_interior
                    + surface_figure_scale * BmA_surface
                    + base_figure_scale * (BmA_base + Bsp_Asp)
                )
            Kc = 3 * omega**2 * ((Bc - Ac) - Bip_Aip)

            # Gravitational coupling torque
            K_int = 4*np.pi*G/5 * (rho[-1] * beta[-1] + (rho[1] - rho[-1]) * beta[1]) \
                    *((rho[0] - rho[1]) * beta[0] * radial[0] ** 5)

            # Frequencies of shell and core
            sigma_s = np.sqrt(Ks / Cs)
            sigma_c = np.sqrt(Kc / Cc)

            # Solve the characterisitc equation associated to the 2nd order ODE
            # that describes the libration amplitude (Eq. 23-24 in Van Hoolst
            # et al. 2008)
            A_mat = np.zeros((2, 2))
            b_mat = np.zeros(2)

            A_mat[0, 0] = -omega**2 * Cs + sigma_s**2 * Cs + 2 * K_int
            A_mat[0, 1] = -2 * K_int
            A_mat[1, 0] = -2 * K_int
            A_mat[1, 1] = -omega**2 * Cc + sigma_c**2 * Cc + 2 * K_int

            b_mat[0] = 2 * ecc * sigma_s**2 * Cs
            b_mat[1] = 2 * ecc * sigma_c**2 * Cc

            sol = np.linalg.solve(A_mat, b_mat)

            gs = sol[0] * radial[-1]
        # Body deformed by tidal deformations (Van Hoolst et al. 2013)
        else:
            if y1 is None:
                # Radial function describing the radial deformation of the body
                raise ValueError("y1 must be provided for non-rigid bodies.")

            # Get radial functions for each layer
            ys = y1[-1]
            y_ot = y1[ocean_idx]
            y_ob = y1[ocean_idx - 1]
            y_i = y1[:ocean_idx]

            # Radial interfaces positions
            r_ot = radial[ocean_idx]
            r_ob = radial[ocean_idx - 1]
            r_i = radial[:ocean_idx]

            # Density of each layer
            rho_shell = rho[ocean_idx + 1]
            rho_ocean = rho[ocean_idx]
            rho_i = rho[:ocean_idx]

            # Now compute the Love numbers for each layer as in the definition
            # given by Van Hoolst et al. (2013) Equation 24
            fact = 4/5 * np.pi * G / r_body**3
            k2s  = fact * rho_shell * (r_body**4*ys - r_ot**4*y_ot)
            k2ob = fact * rho_ocean * (-r_ob**4*y_ob)
            k2ot = fact * rho_ocean * (r_ot**4*y_ot)
            k2i  = (fact * rho_i * r_i ** 4 * y_i).sum()

            # To apply Tricarico (2014) formulation, a sphere must be created
            # inside the ocean layers to compute the MoIs and fluid Love numbers
            # This sphere is defined as the average of the outer and inner
            # ocean radius, and the density is set to the outer ocean density
            #
            # The alternative would be to apply the Clairaut formulation and use
            # the flattenings of each layer
            r_sphere = (r_ot + r_ob) / 2

            r_upd = np.insert(radial, ocean_idx, r_sphere)
            rho_upd     = np.insert(rho, ocean_idx, rho[ocean_idx])

            ocean_idx += 1

            # Compute the eccentricities of the updated radial profile and
            # set the eccenctricities of the bottom of the ocean to zero
            # (as it is a sphere)
            ep, eq = compute_eccentricities(r_upd, rho_upd, omega)

            ep[ocean_idx - 1] = 0
            eq[ocean_idx - 1] = 0

            # Compute the semi-major axes of the updated radial profile
            a, b, c = radial2abc(r_upd, ep, eq)
            beta = (a - b) / a

            # Compute the moments of inertia of the updated radial profile
            Ai, Bi, Ci = moi_ellps(a, b, c, rho_upd)
            kfi = (Bi - Ai) / (q * m_body * r_body ** 2)

            # Now compute the fluid Love numbers for each layer
            kf_s  = kfi[ocean_idx+1:].sum()
            kf_ot = kfi[ocean_idx]
            kf_ob = kfi[ocean_idx-1]
            kf_i  = kfi[:ocean_idx-1].sum()

            f2s   = (kf_s  - k2s)  / kf_s
            f2ot  = (kf_ot - k2ot) / kf_ot
            f2ob  = (kf_ob - k2ob) / kf_ob
            f2i   = (kf_i  - k2i)  / kf_i

            # Equatorial MoIs differences
            BmA_i = (Bi - Ai)[:ocean_idx-1].sum()
            BmA_ob = (Bi - Ai)[ocean_idx-1]
            BmA_ot = (Bi - Ai)[ocean_idx]
            BmA_s = (Bi - Ai)[ocean_idx+1:].sum()

            # Polar MoIs
            C_s = Ci[ocean_idx+1:].sum()
            C_i = Ci[:ocean_idx-1].sum()

            # Flattenings
            beta_s = beta[-1]
            beta_ot = beta[ocean_idx]
            beta_ob = beta[ocean_idx-1]

            # Van hoolst et al. (2013) model based on Clairaut equation
            # BmA_ob2 = -8/15 * np.pi * rho_ocean * beta[0] * radial[0]**5
            # BmA_ot2 = 8/15 * np.pi * rho_ocean * beta[2] * radial[2]**5
            # BmA_s2 = 8/15 * np.pi * rho_shell * (beta[-1] * radial[-1]**5 - beta[2] * radial[2]**5)
            # BmA_i2 = 8/15 * np.pi * rho_core * beta[0] * radial[0]**5

            # Compute torques separately acting on outer shell and interior
            Ks    = 1.5 * omega ** 2 * (BmA_s * f2s + BmA_ot * f2ot)
            Ki    = 1.5 * omega ** 2 * (BmA_i * f2i + BmA_ob * f2ob)

            # Gravitational coupling torque
            K_int = 4/5 * np.pi * G * (BmA_i + BmA_ob) * \
                (rho_ocean*beta_ot + rho_shell*(beta_s-beta_ot))

            # Gravitational and pressure coupling between interior and
            # periodic tidal bulges
            Ki1 = omega**2 * r_body**5 / G * (k2i + k2ob)
            Ki1 *= (rho_ocean * beta_ot + rho_shell * (beta_s - beta_ot))

            Ki2 = (BmA_i + BmA_ob) * 1.5 * omega**2 * r_body**2
            Ki2 *= (rho_shell * (ys / r_body - y_ot / r_ot) + rho_ocean * y_ot /r_ot)

            KiM   = 0.8 * np.pi * G * (-Ki1 + Ki2)
            Kii   = 0.8 * np.pi * G * Ki1
            Kis   = 0.8 * np.pi * G * Ki2

            Ksz   = k2s/kf_s / 4 * BmA_s * omega**2
            Kiz   = k2i/kf_i / 4 * BmA_i * omega**2

            # Assemble all the torques from Equations 47 - 52
            K1    = 2 * (Ks + K_int - Kis)
            K2    = 2 * (Kii - K_int)
            K3    = Ks - KiM + Ksz
            K4    = 2 * (Kis - K_int)
            K5    = 2 * (Ki + K_int - Kii)
            K6    = Ki + KiM + Kiz

            # Analytical solution of the characteristic equation
            # (Eq. 45-46 in Van Hoolst et al. 2013)
            A = K1 * C_i + K5 * C_s
            B = 4 * (K2 * K4 - K1 * K5) * C_i * C_s + (K1 * C_i + K5 * C_s) ** 2
            denom = 2 * C_i * C_s

            # Eigenfrequencies
            sigma1_sq = (A + np.sqrt(B)) / denom
            sigma2_sq = (A - np.sqrt(B)) / denom

            # Libration amplitudes of shell and interior
            denom = C_i * C_s * (omega**2 - sigma1_sq) * (omega**2 - sigma2_sq)
            gs = 4 * ecc * (K3 * K5 - K2 * K6 - omega**2 * K3 * C_i) / denom
            gi = 4 * ecc * (K1 * K6 - K3 * K4 - omega**2 * K6 * C_s) / denom

            gs *= r_body
            # breakpoint()

    return abs(gs)
