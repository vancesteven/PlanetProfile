"""
3D tidal heating for laterally-varying ice shell structure.

Computes spatially resolved tidal dissipation rate H(r, theta, phi) from
local viscoelastic properties (shear modulus, viscosity) and the
degree-2 tidal strain pattern.

Supports Maxwell and Andrade rheologies:
- Maxwell: H = omega * mu * eps0^2 * f(theta,phi) / (mu^2 + omega^2 * eta^2)
- Andrade: Uses complex compliance J*(omega) with transient creep term

Based on the thin-shell formalism of:
- Ojakangas & Stevenson (1989) for the tidal strain pattern
- Tobie et al. (2003) for the dissipation model
- Beuthe (2013, 2015) for the thin-shell treatment
"""
import numpy as np
import logging
from scipy.special import gamma as gamma_func

log = logging.getLogger('PlanetProfile.Lateral.TidalHeating3D')


def TidalStrainPattern(theta_rad, phi_rad, e=0.0, obliq_rad=0.0):
    """ Compute the degree-2 tidal strain heating pattern f(theta, phi).

        For eccentricity tides in a synchronous rotator
        (Ojakangas & Stevenson 1989, Eq. 7-8). The result is normalized
        so that the spherical average is 1.

        Args:
            theta_rad: Colatitude array in radians (nPix,).
            phi_rad: Longitude array in radians (nPix,).
            e: Orbital eccentricity.
            obliq_rad: Obliquity in radians.

        Returns:
            f: Tidal strain pattern (nPix,), normalized to mean = 1.
    """
    cos2t = np.cos(2 * theta_rad)
    cos2p = np.cos(2 * phi_rad)
    sin2t = np.sin(theta_rad)**2  # 1 - cos^2(theta)

    # Eccentricity tide pattern (Ojakangas & Stevenson 1989)
    # Two terms: radial (e0) and librational (e2)
    f_e0 = (5 + 3 * cos2t)  # Radial tide component
    f_e2_cos = (5 + cos2p)  # Azimuthal component (cos)
    f_e2_sin = (5 * cos2p - 1)  # Azimuthal component (sin)

    # Combined (OS89 Eq. 7): averaged over orbit
    f = (f_e0 * f_e2_cos + sin2t * f_e2_sin) / 64

    # Normalize to spherical mean = 1
    # The mean is computed analytically for the OS89 pattern
    f_mean = np.mean(f)
    if f_mean > 0:
        f = f / f_mean

    return f


def _MaxwellDissipation(omega, mu_Pa, eta_Pas):
    """ Compute Im(mu*) for Maxwell rheology: the imaginary part of the
        complex shear modulus.

        For a Maxwell body:
            mu* = mu * i*omega*eta / (mu + i*omega*eta)
            Im(mu*) = omega * eta * mu^2 / (mu^2 + omega^2 * eta^2)

        so that H = omega * Im(mu*) * eps0^2 * f(theta, phi)
                  = omega^2 * eta * mu^2 / (mu^2 + omega^2 * eta^2) * eps0^2 * f

        Following Tobie et al. (2003) Eq. 5.

        Args:
            omega: Tidal angular frequency (rad/s).
            mu_Pa: Shear modulus array (Pa).
            eta_Pas: Viscosity array (Pa*s).

        Returns:
            D: Dissipation factor omega*Im(mu*) in Pa/s, such that
               H [W/m^3] = D * eps0^2 * f(theta, phi).
    """
    return omega**2 * eta_Pas * mu_Pa**2 / (mu_Pa**2 + omega**2 * eta_Pas**2)


def _AndradeDissipation(omega, mu_Pa, eta_Pas, alpha=0.2):
    """ Compute volumetric dissipation factor for Andrade rheology.

        The Andrade complex compliance is:
            J*(omega) = 1/mu - i/(omega*eta)
                        + (i*omega*zeta)^(-alpha) * Gamma(1+alpha) / mu

        where zeta = eta/mu is the Andrade timescale. The result is
        normalized consistently with _MaxwellDissipation: when the
        Andrade transient term vanishes, this returns the same D as
        Maxwell.

        Specifically, we compute:
            D_Andrade = D_Maxwell * Im(mu*_Andrade) / Im(mu*_Maxwell)

        where Im(mu*) = -J_imag / |J*|^2 is the imaginary part of the
        complex shear modulus.

        Args:
            omega: Tidal angular frequency (rad/s).
            mu_Pa: Shear modulus array (Pa).
            eta_Pas: Viscosity array (Pa*s).
            alpha: Andrade exponent (dimensionless, typically 0.2-0.4).

        Returns:
            D: Dissipation factor omega*Im(mu*_Andrade) in Pa/s, such that
               H [W/m^3] = D * eps0^2 * f(theta, phi).
    """
    zeta = eta_Pas / mu_Pa  # Andrade timescale (s)
    omega_zeta = omega * zeta

    Gamma_1a = gamma_func(1 + alpha)
    ow_neg_alpha = omega_zeta**(-alpha)

    # Andrade compliance: J* = J_real + i*J_imag
    # (i*omega*zeta)^(-alpha) = (omega*zeta)^(-alpha) * [cos(alpha*pi/2) - i*sin(alpha*pi/2)]
    J_real_A = 1.0 / mu_Pa + ow_neg_alpha * np.cos(alpha * np.pi / 2) * Gamma_1a / mu_Pa
    J_imag_A = -1.0 / (omega * eta_Pas) - ow_neg_alpha * np.sin(alpha * np.pi / 2) * Gamma_1a / mu_Pa
    J_abs_sq_A = J_real_A**2 + J_imag_A**2
    Im_mu_star_A = -J_imag_A / J_abs_sq_A

    # Maxwell compliance (for normalization): J* = 1/mu - i/(omega*eta)
    J_real_M = 1.0 / mu_Pa
    J_imag_M = -1.0 / (omega * eta_Pas)
    J_abs_sq_M = J_real_M**2 + J_imag_M**2
    Im_mu_star_M = -J_imag_M / J_abs_sq_M

    # Scale Maxwell dissipation by the Andrade/Maxwell Im(mu*) ratio
    D_Maxwell = _MaxwellDissipation(omega, mu_Pa, eta_Pas)
    ratio = Im_mu_star_A / Im_mu_star_M

    return D_Maxwell * ratio


def ComputeTidalHeating3D(Planet, Params, columnPlanets=None, rheology=None):
    """ Compute 3D tidal heating from local viscoelastic properties.

        Supports Maxwell and Andrade rheologies. The volumetric dissipation
        rate at each grid point is:

            H(r, theta, phi) = D(omega, mu, eta) * eps0^2 * f(theta, phi)

        where D is the rheology-dependent dissipation factor, mu = shear
        modulus, eta = viscosity, omega = tidal frequency, eps0 = tidal
        strain amplitude, and f is the degree-2 strain pattern.

        Args:
            Planet: Reference PlanetStruct with Lateral substruct.
            Params: ParamsStruct.
            columnPlanets: Array of PlanetStruct with completed hydrosphere.
                If None, uses mean properties from reference model.
            rheology: 'maxwell', 'andrade', or None (auto-detect from
                Planet.Gravity.andradExponent).

        Returns:
            Planet: Updated with Htidal_Wm3 and HtidalIce_Wm3 on Lateral.
    """
    Lateral = Planet.Lateral

    # Determine rheology
    if rheology is None:
        if hasattr(Planet, 'Gravity') and hasattr(Planet.Gravity, 'andradExponent') \
                and Planet.Gravity.andradExponent is not None:
            rheology = 'andrade'
        else:
            rheology = 'maxwell'

    alpha = None
    if rheology == 'andrade':
        alpha = getattr(Planet.Gravity, 'andradExponent', 0.2) if hasattr(Planet, 'Gravity') else 0.2
        log.info(f'Using Andrade rheology with alpha={alpha}')
    else:
        log.info('Using Maxwell rheology')

    # Tidal frequency (from excitation periods)
    if Planet.Magnetic.Texc_hr is not None and len(Planet.Magnetic.Texc_hr) > 0:
        T_hr = Planet.Magnetic.Texc_hr[0]
    else:
        T_hr = 85.2
        log.warning(f'No excitation periods found; using default T = {T_hr} hr')

    omega = 2 * np.pi / (T_hr * 3600)  # rad/s

    # Tidal strain amplitude: eps0 = (3/2) * e * omega^2 * R / g
    # where g = GM/R^2 is surface gravity. This gives the peak tidal strain
    # for a synchronous satellite (Tobie et al. 2003).
    from PlanetProfile.Utilities.defineStructs import Constants
    e = Planet.Bulk.eccentricity if hasattr(Planet.Bulk, 'eccentricity') and Planet.Bulk.eccentricity is not None else 0.01
    n = Planet.Bulk.meanMotion_radps if hasattr(Planet.Bulk, 'meanMotion_radps') and Planet.Bulk.meanMotion_radps is not None else omega
    g_surf = Constants.G * Planet.Bulk.M_kg / Planet.Bulk.R_m**2
    eps0 = 1.5 * e * n**2 * Planet.Bulk.R_m / g_surf

    # Tidal strain pattern
    f_pattern = TidalStrainPattern(Lateral.theta_rad, Lateral.phi_rad, e=e)

    # Select dissipation function
    if rheology == 'andrade':
        def _dissipation(omega, mu, eta):
            return _AndradeDissipation(omega, mu, eta, alpha=alpha)
    else:
        _dissipation = _MaxwellDissipation

    if columnPlanets is not None:
        HtidalIce = np.zeros(Lateral.nPix)
        HtidalIceI_top = np.zeros(Lateral.nPix)
        HtidalIceI_bot = np.zeros(Lateral.nPix)
        HtidalHP_top = np.zeros(Lateral.nPix)
        HtidalHP_bot = np.zeros(Lateral.nPix)
        # Full radial profiles (variable length per column)
        iceI_profiles = [None] * Lateral.nPix
        iceI_radii = [None] * Lateral.nPix
        hp_profiles = [None] * Lateral.nPix
        hp_radii = [None] * Lateral.nPix

        for i, colPlanet in enumerate(columnPlanets):
            if (hasattr(colPlanet, 'invalidReason')
                    and colPlanet.invalidReason is not None
                    and colPlanet.invalidReason != 'Valid'):
                continue

            nSurf = colPlanet.Steps.nSurfIce

            # --- Surface ice (ice I) dissipation ---
            if nSurf > 0:
                # Get shear modulus and viscosity in ice layers
                if hasattr(colPlanet, 'Seismic') and colPlanet.Seismic.GS_GPa is not None:
                    mu_Pa = colPlanet.Seismic.GS_GPa[:nSurf] * 1e9
                else:
                    mu_Pa = np.full(nSurf, 3.5e9)

                if hasattr(colPlanet, 'eta_Pas') and colPlanet.eta_Pas is not None:
                    eta = colPlanet.eta_Pas[:nSurf]
                else:
                    T_ice = colPlanet.T_K[:nSurf]
                    eta = _ArrheniusViscosity(T_ice)

                # Local dissipation rate at each radial layer
                D_local = _dissipation(omega, mu_Pa, eta)
                H_local = D_local * eps0**2 * f_pattern[i]

                # Column-integrated heating (W/m^3 averaged over ice shell)
                HtidalIce[i] = np.mean(H_local)
                # Top and bottom of surface ice layer
                HtidalIceI_top[i] = H_local[0]
                HtidalIceI_bot[i] = H_local[-1]
                # Full radial profile
                iceI_profiles[i] = H_local.copy()
                if hasattr(colPlanet, 'r_m') and colPlanet.r_m is not None:
                    iceI_radii[i] = colPlanet.r_m[:nSurf].copy()

            # --- HP ice (below ocean) dissipation ---
            nHydro = colPlanet.Steps.nHydro if hasattr(colPlanet.Steps, 'nHydro') and colPlanet.Steps.nHydro is not None else 0
            if nHydro > nSurf and hasattr(colPlanet, 'phase') and colPlanet.phase is not None:
                hp_indices = [j for j in range(nSurf, nHydro) if colPlanet.phase[j] in [3, 5, 6]]
                if len(hp_indices) > 0:
                    if hasattr(colPlanet, 'Seismic') and colPlanet.Seismic.GS_GPa is not None:
                        mu_hp = colPlanet.Seismic.GS_GPa[hp_indices] * 1e9
                    else:
                        mu_hp = np.full(len(hp_indices), 6.0e9)
                    if hasattr(colPlanet, 'eta_Pas') and colPlanet.eta_Pas is not None:
                        eta_hp = colPlanet.eta_Pas[hp_indices]
                    else:
                        T_hp = colPlanet.T_K[hp_indices]
                        eta_hp = _ArrheniusViscosity(T_hp)
                    D_hp = _dissipation(omega, mu_hp, eta_hp)
                    H_hp = D_hp * eps0**2 * f_pattern[i]
                    HtidalHP_top[i] = H_hp[0]
                    HtidalHP_bot[i] = H_hp[-1]
                    hp_profiles[i] = H_hp.copy()
                    if hasattr(colPlanet, 'r_m') and colPlanet.r_m is not None:
                        hp_radii[i] = colPlanet.r_m[hp_indices].copy()

        Lateral.HtidalIce_Wm3 = HtidalIce
        Lateral.HtidalIceI_top_Wm3 = HtidalIceI_top
        Lateral.HtidalIceI_bot_Wm3 = HtidalIceI_bot
        Lateral.HtidalHP_top_Wm3 = HtidalHP_top
        Lateral.HtidalHP_bot_Wm3 = HtidalHP_bot
        Lateral.HtidalIceI_profile_Wm3 = iceI_profiles
        Lateral.rIceI_profile_m = iceI_radii
        Lateral.HtidalHP_profile_Wm3 = hp_profiles
        Lateral.rHP_profile_m = hp_radii

    else:
        # Simplified calculation using reference model properties
        nSurf = Planet.Steps.nSurfIce
        if nSurf > 0:
            if hasattr(Planet, 'Seismic') and Planet.Seismic.GS_GPa is not None:
                mu_top = Planet.Seismic.GS_GPa[0] * 1e9
                mu_bot = Planet.Seismic.GS_GPa[nSurf - 1] * 1e9
                mu_mean = np.mean(Planet.Seismic.GS_GPa[:nSurf]) * 1e9
            else:
                mu_top = mu_bot = mu_mean = 3.5e9

            if hasattr(Planet, 'eta_Pas') and Planet.eta_Pas is not None:
                eta_top = Planet.eta_Pas[0]
                eta_bot = Planet.eta_Pas[nSurf - 1]
                eta_mean = np.exp(np.mean(np.log(Planet.eta_Pas[:nSurf])))
            else:
                T_mean = np.mean(Planet.T_K[:nSurf])
                eta_top = _ArrheniusViscosity(np.array([Planet.T_K[0]]))[0]
                eta_bot = _ArrheniusViscosity(np.array([Planet.T_K[nSurf - 1]]))[0]
                eta_mean = _ArrheniusViscosity(np.array([T_mean]))[0]

            D_ref = _dissipation(omega, np.array([mu_mean]), np.array([eta_mean]))[0]
            Lateral.HtidalIce_Wm3 = D_ref * eps0**2 * f_pattern

            D_top = _dissipation(omega, np.array([mu_top]), np.array([eta_top]))[0]
            D_bot = _dissipation(omega, np.array([mu_bot]), np.array([eta_bot]))[0]
            Lateral.HtidalIceI_top_Wm3 = D_top * eps0**2 * f_pattern
            Lateral.HtidalIceI_bot_Wm3 = D_bot * eps0**2 * f_pattern

    log.info(f'3D tidal heating ({rheology}): mean={np.mean(Lateral.HtidalIce_Wm3):.2e} W/m^3, '
             f'max={np.max(Lateral.HtidalIce_Wm3):.2e} W/m^3')
    if Lateral.HtidalIceI_top_Wm3 is not None and np.any(Lateral.HtidalIceI_top_Wm3 > 0):
        log.info(f'  Ice I top: mean={np.mean(Lateral.HtidalIceI_top_Wm3):.2e}, '
                 f'bot: mean={np.mean(Lateral.HtidalIceI_bot_Wm3):.2e} W/m^3')
    if Lateral.HtidalHP_top_Wm3 is not None and np.any(Lateral.HtidalHP_top_Wm3 > 0):
        log.info(f'  HP ice top: mean={np.mean(Lateral.HtidalHP_top_Wm3):.2e}, '
                 f'bot: mean={np.mean(Lateral.HtidalHP_bot_Wm3):.2e} W/m^3')

    return Planet


def ConvergeTidalFeedback(Planet, Params, columnPlanets, alpha=0.4, maxIter=10, tol=0.01):
    """ Iterate tidal heating <-> thermal structure to convergence.

        Tidal heating affects ice temperature (via volumetric heating in
        thermal profile), which affects viscosity (Arrhenius), which
        affects tidal heating. This feedback loop is converged with
        damped iteration.

        Args:
            Planet: Reference PlanetStruct with Lateral substruct.
            Params: ParamsStruct.
            columnPlanets: Array of PlanetStruct (modified in place).
            alpha: Damping factor for update (0 < alpha < 1).
            maxIter: Maximum number of feedback iterations.
            tol: Convergence tolerance (fractional change in max H_tidal).

        Returns:
            Planet: Updated with converged tidal heating field.
            columnPlanets: Updated with converged column profiles.
    """
    from PlanetProfile.Lateral.LateralStructure import RunLateralColumns

    Lateral = Planet.Lateral

    for iteration in range(maxIter):
        H_old = Lateral.HtidalIce_Wm3.copy() if Lateral.HtidalIce_Wm3 is not None else np.zeros(Lateral.nPix)

        # Step 1: Run column models with current tidal heating
        columnPlanets = RunLateralColumns(Planet, Params)

        # Step 2: Compute new tidal heating from updated viscosity/modulus
        Planet = ComputeTidalHeating3D(Planet, Params, columnPlanets)

        # Step 3: Damped update
        H_new = Lateral.HtidalIce_Wm3
        Lateral.HtidalIce_Wm3 = alpha * H_new + (1 - alpha) * H_old

        # Check convergence
        if np.max(H_old) > 0:
            max_change = np.max(np.abs(H_new - H_old)) / np.max(H_old)
        else:
            max_change = np.inf

        log.info(f'Tidal feedback iteration {iteration + 1}: max change = {max_change:.4f}')

        if max_change < tol:
            log.info(f'Tidal feedback converged after {iteration + 1} iterations')
            return Planet, columnPlanets

    log.warning(f'Tidal feedback did not converge after {maxIter} iterations '
                f'(max change = {max_change:.4f})')
    return Planet, columnPlanets


def _ArrheniusViscosity(T_K, eta0=1e14, Q_J=60e3, T_ref=263.0):
    """ Arrhenius viscosity for ice.

        eta = eta0 * exp(Q/R * (1/T - 1/T_ref))

        Args:
            T_K: Temperature array in Kelvin.
            eta0: Reference viscosity in Pa*s at T_ref.
            Q_J: Activation energy in J/mol.
            T_ref: Reference temperature in K.

        Returns:
            eta: Viscosity array in Pa*s.
    """
    R_gas = 8.314  # J/(mol*K)
    return eta0 * np.exp(Q_J / R_gas * (1.0 / T_K - 1.0 / T_ref))
