""" Radial structure built on the angular-mean density of a three-dimensional mantle.

A spherically symmetric structure is normally built with the density at the mean
temperature of each depth, rho(P, Tbar). The mass of a real three-dimensional mantle is the
angular mean of its local densities, <rho(P, T(theta, phi))>, and the two differ wherever the
density is not linear in temperature (Jensen's inequality), for instance across a phase
change. This module provides the two pieces needed to build the structure on the angular mean:

* AngularMeanDensityOffset computes Delta_rho(P) = <rho(P, T_3D)> - rho(P, Tbar) on the levels of
  a temperature field, at the pressures of a given radial structure.
* AngularMeanDensityEOS wraps a silicate EOS and adds Delta_rho(P) to its density only; every
  other property is delegated unchanged. The correction is zero above the Moho pressure, so an
  explicit crust is never altered.

Set Planet.Sil.densityCorrection = dict(P_MPa=..., drho_kgm3=..., P_moho_MPa=...,
fill_low=0.0, fill_high=0.0) and GetPlanetInnerEOS wraps the silicate EOS. Because Delta_rho
depends on the pressure profile, a self-consistent structure is obtained by iterating: build,
compute the offset at the new pressures, rebuild, until the core radius no longer changes and
the pressure profile changes by less than a chosen tolerance. With Delta_rho = 0 the wrapped
EOS gives back the unwrapped structure exactly.

Same construction as the lunar interior study (lunar-bsm-interior-models,
scripts/stage_e9a5/e9a5_lib.py and scripts/stage_e9a6/e9a6_lib.py).
"""
import numpy as np


class AngularMeanDensityEOS:
    """ Silicate EOS with a prescribed degree-0 density correction.

        rho_wrapped(P, T) = rho_EOS(P, T) + Delta_rho(P), with Delta_rho = 0 for P < P_moho_MPa.
        Delta_rho is interpolated linearly in P between the given levels; outside them it takes
        fill_low and fill_high (0 by default). Every other attribute is delegated to the base EOS.
    """
    def __init__(self, base, P_MPa, drho_kgm3, P_moho_MPa=0.0, fill_low=0.0, fill_high=0.0):
        P = np.asarray(P_MPa, dtype=float); d = np.asarray(drho_kgm3, dtype=float)
        if P.shape != d.shape or P.ndim != 1:
            raise ValueError('P_MPa and drho_kgm3 must be 1-D arrays of the same length.')
        if np.any(np.diff(P) <= 0):
            raise ValueError('P_MPa must be strictly increasing.')
        object.__setattr__(self, '_base', base)
        object.__setattr__(self, '_P', P)
        object.__setattr__(self, '_d', d)
        object.__setattr__(self, '_Pmoho', 0.0 if P_moho_MPa is None else float(P_moho_MPa))
        object.__setattr__(self, '_lo', float(fill_low))
        object.__setattr__(self, '_hi', float(fill_high))
        object.__setattr__(self, 'nDensityCalls', 0)

    def __getattr__(self, k):
        return getattr(object.__getattribute__(self, '_base'), k)

    def correction(self, P_MPa):
        P = np.asarray(P_MPa, dtype=float)
        d = np.interp(P, self._P, self._d, left=self._lo, right=self._hi)
        return np.where(P < self._Pmoho, 0.0, d)

    def fn_rho_kgm3(self, P_MPa, T_K, grid=False):
        rho = self._base.fn_rho_kgm3(P_MPa, T_K, grid)
        object.__setattr__(self, 'nDensityCalls', self.nDensityCalls + 1)
        if grid:
            # grid mode returns an (nP, nT) array; correct along the P axis
            return rho + self.correction(P_MPa)[:, None]
        return rho + self.correction(P_MPa)


def WrapSilicateDensityCorrection(EOS, densityCorrection):
    """ Return EOS wrapped with densityCorrection (dict as in Planet.Sil.densityCorrection),
        or EOS itself if densityCorrection is None or identically zero.
    """
    if densityCorrection is None:
        return EOS
    d = np.asarray(densityCorrection['drho_kgm3'], dtype=float)
    lo = float(densityCorrection.get('fill_low', 0.0)); hi = float(densityCorrection.get('fill_high', 0.0))
    if not np.any(d != 0.0) and lo == 0.0 and hi == 0.0:
        return EOS
    return AngularMeanDensityEOS(EOS, densityCorrection['P_MPa'], d,
                                 densityCorrection.get('P_moho_MPa', 0.0), lo, hi)


def AngularMeanDensityOffset(EOS, z_lev_km, T3D_K, weights, z_nodes_km, P_nodes_MPa, T_nodes_K):
    """ Delta_rho on the levels of a 3-D temperature field, at the pressures of a structure.

        Parameters
        ----------
        EOS : object with fn_rho_kgm3(P_MPa, T_K)
        z_lev_km : (nLev,) depths of the field levels (km)
        T3D_K : (nLev, ...) temperatures on each level (any map shape)
        weights : map-shaped quadrature weights (for example cell areas), same for every level
        z_nodes_km, P_nodes_MPa, T_nodes_K : radial structure (depth, pressure, mean temperature)
        Returns
        -------
        P_lev_MPa, drho_kgm3 : (nLev,) pressure of each level and its density offset
    """
    w = np.asarray(weights, dtype=float); W = w.sum()
    z = np.asarray(z_nodes_km, dtype=float); o = np.argsort(z)
    P_lev = np.interp(z_lev_km, z[o], np.asarray(P_nodes_MPa, dtype=float)[o])
    T_lev = np.interp(z_lev_km, z[o], np.asarray(T_nodes_K, dtype=float)[o])
    drho = np.empty(len(P_lev))
    for k in range(len(P_lev)):
        Tf = np.ascontiguousarray(np.asarray(T3D_K[k], dtype=float).ravel())
        Pf = np.full(Tf.size, float(P_lev[k]))
        mean = float((EOS.fn_rho_kgm3(Pf, Tf).reshape(np.shape(T3D_K[k])) * w).sum() / W)
        ref = float(EOS.fn_rho_kgm3(np.array([P_lev[k]]), np.array([T_lev[k]]))[0])
        drho[k] = mean - ref
    return P_lev, drho
