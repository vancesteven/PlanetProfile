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

        Combines eccentricity and obliquity tides for a synchronous rotator
        following Ojakangas & Stevenson (1989, Eq. 7-10). The result is
        normalized so that the spherical average is 1.

        Args:
            theta_rad: Colatitude array in radians (nPix,).
            phi_rad: Longitude array in radians (nPix,).
            e: Orbital eccentricity.
            obliq_rad: Obliquity in radians.

        Returns:
            f: Tidal strain pattern (nPix,), normalized to mean = 1.

        References:
            Ojakangas & Stevenson (1989), Icarus 81, 220-241
            Beuthe (2013), Icarus 223, 308-329
    """
    cost = np.cos(theta_rad)
    sint = np.sin(theta_rad)
    cos2t = np.cos(2 * theta_rad)
    sin2t = np.sin(2 * theta_rad)
    cos2p = np.cos(2 * phi_rad)
    sin2p = np.sin(2 * phi_rad)

    # Avoid singularities at poles
    eps_pole = 1e-12
    sint_safe = np.maximum(np.abs(sint), eps_pole) * np.sign(sint + eps_pole)

    # Initialize total pattern
    f = np.zeros_like(theta_rad)

    # Eccentricity tide pattern (OS89 Eq. 7-8)
    if abs(e) > 1e-10:
        # Radial (e0) and librational (e2) terms
        f_e0 = (5 + 3 * cos2t)  # Radial tide component
        f_e2_cos = (5 + cos2p)  # Azimuthal component (cos)
        f_e2_sin = (5 * cos2p - 1)  # Azimuthal component (sin)

        # Combined, orbit-averaged, proportional to e²
        f_ecc = (f_e0 * f_e2_cos + sint**2 * f_e2_sin) / 64
        f += e**2 * f_ecc

    # Obliquity tide pattern (OS89 Eq. 9-10)
    if abs(obliq_rad) > 1e-10:
        # Diurnal libration pattern from obliquity
        # Three components: radial, librational cos(φ), librational sin(φ)
        f_o0 = (3 + cost**2)  # Radial component
        f_o1_cos = sint_safe**2  # Librational cos component
        f_o1_sin = -2 * sint_safe * cost  # Librational sin component

        # Combined, orbit-averaged, proportional to obliq²
        # The sin2p term comes from the φ-dependence of obliquity forcing
        f_obliq = (f_o0 + f_o1_cos * cos2p + f_o1_sin * sin2p) / 16
        f += obliq_rad**2 * f_obliq

    # Normalize to spherical mean = 1
    f_mean = np.mean(f)
    if f_mean > 1e-15:
        f = f / f_mean
    else:
        # If both e and obliq are tiny, return uniform pattern
        f = np.ones_like(theta_rad)

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
    # Add small epsilon to prevent division by zero
    eps = 1e-20
    denominator = mu_Pa**2 + omega**2 * eta_Pas**2 + eps
    return omega**2 * eta_Pas * mu_Pa**2 / denominator


def _AndradeDissipation(omega, mu_Pa, eta_Pas, alpha=0.2, zeta_pa=1.0):
    """ Compute volumetric dissipation factor for Andrade rheology.

        The Andrade complex compliance is (PyALMA3 convention):
            J*(omega) = 1/mu - i/(omega*eta)
                        + [Gamma(1+alpha)/zeta_pa] * (1/mu) * (omega*tau_M)^(-alpha)
                          * [cos(alpha*pi/2) - i*sin(alpha*pi/2)]

        where tau_M = eta/mu is the Maxwell time and zeta_pa is the
        dimensionless Andrade creep amplitude parameter (PyALMA3 convention).
        Smaller zeta_pa means more Andrade creep; zeta_pa=1 recovers the
        simplest Andrade model where transient and Maxwell timescales are
        equal.  This is the same convention used in Gravity.py for Love
        number calculations.

        The result is normalized consistently with _MaxwellDissipation:
        when the Andrade transient term vanishes, this returns the same D
        as Maxwell.

        Args:
            omega: Tidal angular frequency (rad/s).
            mu_Pa: Shear modulus array (Pa).
            eta_Pas: Viscosity array (Pa*s).
            alpha: Andrade exponent (dimensionless, typically 0.2-0.4).
            zeta_pa: Andrade zeta in PyALMA3 convention (scalar or array
                matching mu_Pa). Divides Gamma(1+alpha) linearly.
                Default 1.0 (no amplification).

        Returns:
            D: Dissipation factor omega*Im(mu*_Andrade) in Pa/s, such that
               H [W/m^3] = D * eps0^2 * f(theta, phi).
    """
    # Add small epsilon to prevent division by zero
    eps = 1e-20

    tau_M = eta_Pas / (mu_Pa + eps)  # Maxwell time (s)
    omega_tau = omega * tau_M

    Gamma_1a = gamma_func(1 + alpha)
    ow_neg_alpha = omega_tau**(-alpha)

    # Andrade creep coefficient: Gamma(1+alpha) / zeta_pa
    # At zeta_pa=1 this equals Gamma(1+alpha); at zeta_pa=0.005 it is 200x larger.
    beta = Gamma_1a / (zeta_pa + eps)

    # Andrade compliance: J* = J_real + i*J_imag
    # (i*omega*tau_M)^(-alpha) = (omega*tau_M)^(-alpha) * [cos(alpha*pi/2) - i*sin(alpha*pi/2)]
    J_real_A = 1.0 / (mu_Pa + eps) + ow_neg_alpha * np.cos(alpha * np.pi / 2) * beta / (mu_Pa + eps)
    J_imag_A = -1.0 / (omega * eta_Pas + eps) - ow_neg_alpha * np.sin(alpha * np.pi / 2) * beta / (mu_Pa + eps)
    J_abs_sq_A = J_real_A**2 + J_imag_A**2 + eps
    Im_mu_star_A = -J_imag_A / J_abs_sq_A

    # Maxwell compliance (for normalization): J* = 1/mu - i/(omega*eta)
    J_real_M = 1.0 / (mu_Pa + eps)
    J_imag_M = -1.0 / (omega * eta_Pas + eps)
    J_abs_sq_M = J_real_M**2 + J_imag_M**2 + eps
    Im_mu_star_M = -J_imag_M / J_abs_sq_M

    # Scale Maxwell dissipation by the Andrade/Maxwell Im(mu*) ratio
    D_Maxwell = _MaxwellDissipation(omega, mu_Pa, eta_Pas)
    ratio = Im_mu_star_A / (Im_mu_star_M + eps)

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
    andrade_zeta = None
    if rheology == 'andrade':
        alpha = getattr(Planet.Gravity, 'andradExponent', 0.2) if hasattr(Planet, 'Gravity') else 0.2
        # Retrieve per-phase Andrade zeta (PyALMA3 convention)
        if hasattr(Planet.Gravity, 'andrade_zeta') and Planet.Gravity.andrade_zeta is not None:
            andrade_zeta = Planet.Gravity.andrade_zeta
            if isinstance(andrade_zeta, dict):
                log.info(f'Using Andrade rheology with alpha={alpha}, '
                         f'per-phase zeta (PyALMA3 convention): {andrade_zeta}')
            else:
                log.info(f'Using Andrade rheology with alpha={alpha}, zeta={andrade_zeta}')
        else:
            andrade_zeta = 1.0
            log.info(f'Using Andrade rheology with alpha={alpha}, zeta=1.0 (default)')
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

    # Map phase IDs to andrade_zeta dict keys
    _phase_to_zeta_key = {1: 'Ih', 2: 'II', 3: 'III', 5: 'V', 6: 'VI'}

    def _get_zeta_for_phases(phase_arr):
        """ Look up per-point zeta_pa from phase array and andrade_zeta spec. """
        if andrade_zeta is None or not isinstance(andrade_zeta, dict):
            return andrade_zeta if andrade_zeta is not None else 1.0
        zeta_out = np.ones(len(phase_arr))
        for j, ph in enumerate(phase_arr):
            key = _phase_to_zeta_key.get(int(ph), 'default')
            zeta_out[j] = andrade_zeta.get(key, andrade_zeta.get('default', 1.0))
        return zeta_out

    # Select dissipation function
    if rheology == 'andrade':
        def _dissipation(omega, mu, eta, zeta_pa=1.0):
            return _AndradeDissipation(omega, mu, eta, alpha=alpha, zeta_pa=zeta_pa)
    else:
        def _dissipation(omega, mu, eta, zeta_pa=1.0):
            return _MaxwellDissipation(omega, mu, eta)

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

                # Per-point zeta for ice I (all same phase)
                zeta_iceI = _get_zeta_for_phases(np.full(nSurf, 1))

                # Local dissipation rate at each radial layer
                D_local = _dissipation(omega, mu_Pa, eta, zeta_pa=zeta_iceI)
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
                    zeta_hp = _get_zeta_for_phases(colPlanet.phase[hp_indices])
                    D_hp = _dissipation(omega, mu_hp, eta_hp, zeta_pa=zeta_hp)
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

            zeta_Ih = _get_zeta_for_phases(np.array([1]))[0] if isinstance(andrade_zeta, dict) else (andrade_zeta if andrade_zeta is not None else 1.0)
            D_ref = _dissipation(omega, np.array([mu_mean]), np.array([eta_mean]), zeta_pa=zeta_Ih)[0]
            Lateral.HtidalIce_Wm3 = D_ref * eps0**2 * f_pattern

            D_top = _dissipation(omega, np.array([mu_top]), np.array([eta_top]), zeta_pa=zeta_Ih)[0]
            D_bot = _dissipation(omega, np.array([mu_bot]), np.array([eta_bot]), zeta_pa=zeta_Ih)[0]
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


def CalcEquilibriumIceThickness(Planet, Params, columnPlanets):
    """ Compute self-consistent equilibrium ice shell thickness per pixel.

        Solves the steady-state heat balance equation for each HEALPix pixel:
            k * (T_bot - T_surf) / d_ice = q_basal + H_tidal(pixel) * d_ice

        where tidal heating H_tidal depends on the ice viscosity profile, which
        depends on d_ice. Iterates until ice thickness converges.

        Physical model from Tobie et al. (2003, JGR doi:10.1029/2003JE002099):
        - Conductive heat transport through ice shell
        - Volumetric tidal dissipation in ice (computed via TidalPy or thin-shell)
        - Basal heat flux from silicate mantle
        - Self-consistent coupling between thickness, temperature, and heating

        Expected result: Reproduces Tobie et al. (2003) Figure 12a behavior:
        - Thicker ice at sub/anti-Jovian points (low tidal dissipation)
        - Thinner ice at poles and mid-latitudes (high dissipation)
        - ~5 km peak-to-peak variation for ~20-25 km mean shell
        - Nearly uniform surface heat flux (~35-40 mW/m² for Europa)

        Args:
            Planet: PlanetStruct with Lateral.HtidalIce_Wm3 already computed.
            Params: ParamsStruct.
            columnPlanets: Array of PlanetStruct from RunLateralColumns.

        Returns:
            Planet: Updated with Lateral.dIce_m set to equilibrium values.
                    Also sets Lateral.qSurf_Wm2, equilibriumIterations, equilibriumResidual_m.

        References:
            Tobie et al. (2003), JGR 108(E11), 5124, doi:10.1029/2003JE002099
    """
    from PlanetProfile.Lateral.LateralStructure import RunLateralColumns

    Lateral = Planet.Lateral

    # Check prerequisites
    if Lateral.HtidalIce_Wm3 is None:
        raise ValueError('CalcEquilibriumIceThickness requires Lateral.HtidalIce_Wm3 to be computed first')

    if not hasattr(Planet.Bulk, 'Tb_K') or not hasattr(Planet.Bulk, 'Tsurf_K'):
        raise ValueError('CalcEquilibriumIceThickness requires Planet.Bulk.Tb_K and Tsurf_K')

    # Parameters
    k_ice = Lateral.kThermIce_WmK  # Ice thermal conductivity (W/m/K)
    T_surf = Planet.Bulk.Tsurf_K  # Surface temperature (K)
    T_bot = Planet.Bulk.Tb_K  # Basal temperature (K)
    max_iter = Lateral.equilibriumMaxIter
    tol_m = Lateral.equilibriumTol_m
    nPix = Lateral.nPix

    # Compute basal heat flux from silicate mantle (uniform across pixels)
    # Use override if provided, otherwise compute from silicate properties
    if Lateral.qBasal_Wm2 is not None:
        q_basal = Lateral.qBasal_Wm2
        log.info(f'Using prescribed q_basal = {q_basal*1e3:.2f} mW/m² (from Lateral.qBasal_Wm2)')
    elif hasattr(Planet.Sil, 'Qrad_Wkg') and hasattr(Planet.Sil, 'Htidal_Wm3'):
        # Estimate mantle mass and volume from reference Planet
        R_planet = Planet.Bulk.R_m
        R_mantle = getattr(Planet.Sil, 'Rmean_m', 0.9 * R_planet)  # Approximate
        rho_sil = getattr(Planet.Sil, 'rhoSilWithCore_kgm3', 3000.0)
        V_mantle = (4.0/3.0) * np.pi * R_mantle**3
        M_mantle = rho_sil * V_mantle

        Q_rad_total = Planet.Sil.Qrad_Wkg * M_mantle  # Total radiogenic power (W)
        Q_tidal_total = Planet.Sil.Htidal_Wm3 * V_mantle  # Total silicate tidal power (W)
        Q_mantle_total = Q_rad_total + Q_tidal_total

        # Divide by surface area to get heat flux
        A_surface = 4.0 * np.pi * R_planet**2
        q_basal = Q_mantle_total / A_surface  # W/m²
        log.info(f'Computed q_basal = {q_basal*1e3:.2f} mW/m² from silicate heating')
    else:
        # Default fallback for Europa: ~5-10 mW/m² from silicates
        q_basal = 0.007  # W/m² (7 mW/m²)
        log.warning(f'Using default q_basal = {q_basal*1e3:.1f} mW/m² (silicate heat flux not found)')

    log.info(f'Equilibrium ice thickness solver starting:')
    log.info(f'  T_surf = {T_surf:.1f} K, T_bot = {T_bot:.3f} K')
    log.info(f'  k_ice = {k_ice:.2f} W/m/K, q_basal = {q_basal*1e3:.2f} mW/m²')
    log.info(f'  Tolerance = {tol_m:.1f} m, Max iterations = {max_iter}')

    # Initialize ice thickness to reference 1D mean
    d_ice_m = np.full(nPix, Planet.zb_km * 1e3)  # Convert km to m
    d_ice_old = d_ice_m.copy()

    # Self-consistent iteration loop
    for iteration in range(max_iter):
        # Update Lateral.dIce_m with current guess
        Lateral.dIce_m = d_ice_m.copy()

        # Re-run lateral columns with updated ice thickness
        # This recomputes thermal profiles and viscosity
        columnPlanets = RunLateralColumns(Planet, Params,
                                          checkpoint_interval=9999,  # No checkpointing during iteration
                                          resume_from_checkpoint=False)

        # Recompute tidal heating with updated thermal structure
        Planet = ComputeTidalHeating3D(Planet, Params, columnPlanets)
        H_tidal = Lateral.HtidalIce_Wm3  # Updated tidal heating per pixel (W/m³)

        # Solve quadratic heat balance per pixel:
        # H_tidal * d² + q_basal * d - k * (T_bot - T_surf) = 0
        # Positive root: d = (-q_basal + sqrt(q_basal² + 4*H_tidal*k*ΔT)) / (2*H_tidal)
        deltaT = T_bot - T_surf
        k_deltaT = k_ice * deltaT

        d_ice_new = np.zeros(nPix)
        for i in range(nPix):
            H_i = H_tidal[i]
            if H_i > 1e-15:  # Nonzero tidal heating
                # Quadratic formula
                discriminant = q_basal**2 + 4.0 * H_i * k_deltaT
                if discriminant >= 0:
                    d_ice_new[i] = (-q_basal + np.sqrt(discriminant)) / (2.0 * H_i)
                else:
                    # Negative discriminant (should not happen physically)
                    # Fall back to pure conduction
                    d_ice_new[i] = k_deltaT / q_basal
                    log.warning(f'Pixel {i}: negative discriminant, using pure conduction')
            else:
                # Zero tidal heating: pure conduction
                # q_basal * d = k * ΔT → d = k * ΔT / q_basal
                d_ice_new[i] = k_deltaT / q_basal

        # Check convergence
        residual = np.abs(d_ice_new - d_ice_old)
        max_residual = np.max(residual)

        log.info(f'  Iteration {iteration+1}/{max_iter}: max residual = {max_residual:.1f} m')

        if max_residual < tol_m:
            log.info(f'Equilibrium ice thickness converged after {iteration+1} iterations')
            Lateral.equilibriumIterations = iteration + 1
            Lateral.equilibriumResidual_m = max_residual
            d_ice_m = d_ice_new
            break

        # Update for next iteration with damping to improve stability
        damping = 0.5  # Mix 50% old, 50% new
        d_ice_m = damping * d_ice_old + (1.0 - damping) * d_ice_new
        d_ice_old = d_ice_m.copy()

    else:
        # Max iterations reached without convergence
        log.warning(f'Equilibrium ice thickness did not converge after {max_iter} iterations '
                    f'(max residual = {max_residual:.1f} m)')
        Lateral.equilibriumIterations = max_iter
        Lateral.equilibriumResidual_m = max_residual
        d_ice_m = d_ice_new  # Use final iteration result

    # Store final equilibrium ice thickness
    Lateral.dIce_m = d_ice_m

    # Compute equilibrium surface heat flux per pixel
    Lateral.qSurf_Wm2 = k_ice * deltaT / d_ice_m

    # Print summary
    log.info(f'Equilibrium ice thickness:')
    log.info(f'  Mean: {np.mean(d_ice_m)/1e3:.2f} km')
    log.info(f'  Range: [{np.min(d_ice_m)/1e3:.2f}, {np.max(d_ice_m)/1e3:.2f}] km')
    log.info(f'  Variation: {(np.max(d_ice_m) - np.min(d_ice_m))/1e3:.2f} km peak-to-peak')
    log.info(f'  Surface heat flux: mean = {np.mean(Lateral.qSurf_Wm2)*1e3:.2f} mW/m², '
             f'range = [{np.min(Lateral.qSurf_Wm2)*1e3:.2f}, {np.max(Lateral.qSurf_Wm2)*1e3:.2f}] mW/m²')

    return Planet, columnPlanets
