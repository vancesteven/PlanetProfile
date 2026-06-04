"""
Spatial grid management and spherical harmonic transforms for 3D lateral structure.

Supports HEALPix (equal-area) and regular lat-lon grids. Spherical harmonic
coefficients use 4pi-normalization consistent with MoonMag shape files.
"""
import numpy as np
import logging

log = logging.getLogger('PlanetProfile.Lateral.SpatialGrid')

# Try importing healpy; fall back gracefully
try:
    import healpy as hp
    _HEALPY_AVAILABLE = True
except ImportError:
    _HEALPY_AVAILABLE = False


def InitGrid(Lateral):
    """ Initialize spatial grid on Lateral substruct.

        Sets theta_rad, phi_rad, nPix, pixArea_sr based on gridType.

        Args:
            Lateral: LateralSubstruct instance (modified in place)
    """
    if Lateral.gridType == 'healpix':
        if not _HEALPY_AVAILABLE:
            raise ImportError(
                'healpy is required for HEALPix grids. Install with: conda install -c conda-forge healpy')
        nSide = Lateral.nSide
        Lateral.nPix = hp.nside2npix(nSide)
        Lateral.pixArea_sr = hp.nside2pixarea(nSide)
        # Get colatitude and longitude for each pixel center
        Lateral.theta_rad, Lateral.phi_rad = hp.pix2ang(nSide, np.arange(Lateral.nPix))

    elif Lateral.gridType == 'latlon':
        if Lateral.nLat is None or Lateral.nLon is None:
            raise ValueError('nLat and nLon must be set for latlon grid type.')
        nLat = Lateral.nLat
        nLon = Lateral.nLon
        Lateral.nPix = nLat * nLon
        # Colatitude from pole to pole (excluding exact poles for numerical stability)
        colat = np.linspace(np.pi / (2 * nLat), np.pi - np.pi / (2 * nLat), nLat)
        lon = np.linspace(0, 2 * np.pi * (1 - 1 / nLon), nLon)
        colat_grid, lon_grid = np.meshgrid(colat, lon, indexing='ij')
        Lateral.theta_rad = colat_grid.ravel()
        Lateral.phi_rad = lon_grid.ravel()
        # Pixel areas (solid angle elements)
        dtheta = np.pi / nLat
        dphi = 2 * np.pi / nLon
        Lateral.pixArea_sr = np.abs(np.sin(Lateral.theta_rad)) * dtheta * dphi
    else:
        raise ValueError(f'Unknown gridType: {Lateral.gridType}. Use "healpix" or "latlon".')

    log.debug(f'Initialized {Lateral.gridType} grid with {Lateral.nPix} pixels.')


def SHtoGrid(Cpq, Spq, pMax, theta_rad, phi_rad):
    """ Evaluate spherical harmonic coefficients on a grid.

        Uses 4pi-normalized real spherical harmonics, consistent with
        MoonMag shape file convention.

        Args:
            Cpq: Cosine coefficients, shape (pMax+1, pMax+1). Cpq[p,q] is degree p, order q.
            Spq: Sine coefficients, shape (pMax+1, pMax+1). Spq[p,q] is degree p, order q.
            pMax: Maximum degree of expansion.
            theta_rad: Colatitude array in radians, shape (nPix,).
            phi_rad: Longitude array in radians, shape (nPix,).

        Returns:
            field: Array of values at each grid point, shape (nPix,).
    """
    from scipy.special import lpmv

    nPix = len(theta_rad)
    field = np.zeros(nPix)
    cos_theta = np.cos(theta_rad)

    for p in range(pMax + 1):
        for q in range(p + 1):
            # 4pi-normalized associated Legendre function
            Pnm = _assoc_legendre_4pi(p, q, cos_theta)
            if q == 0:
                field += Cpq[p, q] * Pnm
            else:
                field += Pnm * (Cpq[p, q] * np.cos(q * phi_rad)
                                + Spq[p, q] * np.sin(q * phi_rad))

    return field


def GridToSH(field, theta_rad, phi_rad, pixArea_sr, pMax):
    """ Decompose a grid field into spherical harmonic coefficients.

        Uses numerical quadrature over the grid with 4pi-normalization.

        Args:
            field: Array of values at each grid point, shape (nPix,).
            theta_rad: Colatitude array in radians, shape (nPix,).
            phi_rad: Longitude array in radians, shape (nPix,).
            pixArea_sr: Area of each pixel in steradians (scalar for HEALPix, array for latlon).
            pMax: Maximum degree to decompose to.

        Returns:
            Cpq: Cosine coefficients, shape (pMax+1, pMax+1).
            Spq: Sine coefficients, shape (pMax+1, pMax+1).
    """
    nPix = len(field)
    Cpq = np.zeros((pMax + 1, pMax + 1))
    Spq = np.zeros((pMax + 1, pMax + 1))
    cos_theta = np.cos(theta_rad)

    for p in range(pMax + 1):
        for q in range(p + 1):
            Pnm = _assoc_legendre_4pi(p, q, cos_theta)
            # Numerical quadrature: integral of field * Y_pq * dOmega
            integrand_c = field * Pnm * np.cos(q * phi_rad) * pixArea_sr
            integrand_s = field * Pnm * np.sin(q * phi_rad) * pixArea_sr
            # For 4pi normalization, basis functions satisfy:
            #   integral(Y_pq * Y_p'q' * dOmega) = 4pi * delta
            # So the analysis coefficients are:
            #   Cpq = (1/4pi) * integral(field * Ypq_cos * dOmega)
            Cpq[p, q] = np.sum(integrand_c) / (4 * np.pi)
            if q > 0:
                Spq[p, q] = np.sum(integrand_s) / (4 * np.pi)

    return Cpq, Spq


def IntegrateOverSphere(field, pixArea_sr, R_m):
    """ Area-weighted integration of a field over the sphere.

        Args:
            field: Array of values at each grid point, shape (nPix,).
            pixArea_sr: Area of each pixel in steradians (scalar or array).
            R_m: Radius of the sphere in meters.

        Returns:
            integral: Total integral in units of [field] * m^2.
    """
    return R_m**2 * np.sum(field * pixArea_sr)


def _assoc_legendre_4pi(n, m, cos_theta):
    """ Compute 4pi-normalized associated Legendre function P_n^m(cos(theta)).

        The 4pi normalization satisfies:
            integral(Y_nm * Y_n'm' * dOmega) = 4pi * delta_nn' * delta_mm'

        This matches the convention used in MoonMag shape files.

        Args:
            n: Degree (int).
            m: Order (int, 0 <= m <= n).
            cos_theta: Cosine of colatitude, shape (nPix,).

        Returns:
            Pnm: 4pi-normalized associated Legendre values, shape (nPix,).
    """
    from scipy.special import lpmv
    from math import factorial

    # scipy lpmv gives the unnormalized associated Legendre function
    # (without Condon-Shortley phase)
    Pnm_unnorm = lpmv(m, n, cos_theta)

    # 4pi normalization factor:
    # N_nm = sqrt((2n+1) * (2-delta_m0) * (n-m)! / (n+m)!)
    if m == 0:
        norm = np.sqrt(2.0 * n + 1.0)
    else:
        norm = np.sqrt((2.0 * n + 1.0) * 2.0 * factorial(n - m) / factorial(n + m))

    return norm * Pnm_unnorm
