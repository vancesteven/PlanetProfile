"""Libration amplitude helpers for PlanetProfile gravity calculations."""

import logging

import numpy as np

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
                       ocean_idx=1, y1=None, rek2=None):
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
            Ks = 3 * omega**2 * ((Bs - As) + Bsp_Asp)
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
