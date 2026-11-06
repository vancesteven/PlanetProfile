""" Silicate properties on a three-dimensional temperature field.

PlanetProfile's Lateral module varies the hydrosphere only; the silicate mantle and the core stay
at the radial reference. This module evaluates the mantle EOS on a 3-D temperature field
T(z, lat, lon) around a radial structure: each level takes the pressure of the structure at its
depth, and the EOS gives rho, VP, VS, KS and GS at every point of the level. This is a
pressure-temperature decomposition about the radial structure, not a 3-D hydrostatic solution:
lateral pressure variations are neglected.

    from PlanetProfile.Thermodynamics.Silicate3D import EvaluateSilicate3D
    out = EvaluateSilicate3D(Planet.Sil.EOS, z_km, T3D_K, w, zStruct_km, PStruct_MPa)
    out['VS_kms']        # same shape as T3D_K
    out['VS_kms_ref']    # radial reference VS(P, Tbar) on the levels
    out['VS_kms_mean']   # weighted mean of the 3-D field on each level

Same construction as the lunar interior study (lunar-bsm-interior-models,
scripts/stage_e7/e7_fields.py).
"""
import numpy as np

PROPS = ['rho_kgm3', 'VP_kms', 'VS_kms', 'KS_GPa', 'GS_GPa']


def LevelMean(field, weights):
    """ Weighted mean of a (nLev, ...) field on each level, with map-shaped weights. """
    w = np.asarray(weights, dtype=float)
    return np.array([(np.asarray(field[k], dtype=float) * w).sum() / w.sum() for k in range(len(field))])


def EvaluateSilicate3D(EOS, z_lev_km, T3D_K, weights, z_nodes_km, P_nodes_MPa, Tbar_K=None):
    """ Evaluate a silicate EOS on a 3-D temperature field around a radial structure.

        Parameters
        ----------
        EOS : object with fn_rho_kgm3 and fn_Seismic (returning VP, VS, KS, GS)
        z_lev_km : (nLev,) depths of the field levels (km), inside the mantle of the structure
        T3D_K : (nLev, ...) temperatures on each level
        weights : map-shaped quadrature weights, for the level means and the reference temperature
        z_nodes_km, P_nodes_MPa : depth and pressure of the radial structure (mantle nodes)
        Tbar_K : optional (nLev,) reference temperatures; default the weighted mean of each level
        Returns
        -------
        dict with P_MPa and Tbar_K (nLev,), every property in PROPS on the 3-D field, its radial
        reference '<prop>_ref' at (P, Tbar), its level mean '<prop>_mean', and the number of
        points outside the EOS domain when the EOS exposes Pmin, Pmax, Tmin, Tmax.
    """
    T3 = np.asarray(T3D_K, dtype=float)
    z = np.asarray(z_nodes_km, dtype=float); o = np.argsort(z)
    P_lev = np.interp(z_lev_km, z[o], np.asarray(P_nodes_MPa, dtype=float)[o])
    Tbar = LevelMean(T3, weights) if Tbar_K is None else np.asarray(Tbar_K, dtype=float)
    P3 = np.broadcast_to(P_lev.reshape((-1,) + (1,) * (T3.ndim - 1)), T3.shape)
    out = dict(P_MPa=P_lev, Tbar_K=Tbar)
    if all(hasattr(EOS, k) for k in ('Pmin', 'Pmax', 'Tmin', 'Tmax')):
        bad = (P3 < EOS.Pmin) | (P3 > EOS.Pmax) | (T3 < EOS.Tmin) | (T3 > EOS.Tmax)
        out['nOutsideDomain'] = int(bad.sum())
    Pf, Tf = np.ascontiguousarray(P3.ravel()), np.ascontiguousarray(T3.ravel())
    VP, VS, KS, GS = EOS.fn_Seismic(Pf, Tf)
    rho = EOS.fn_rho_kgm3(Pf, Tf)
    for k, v in zip(PROPS, (rho, VP, VS, KS, GS)):
        out[k] = np.asarray(v).reshape(T3.shape)
    VP0, VS0, KS0, GS0 = EOS.fn_Seismic(np.ascontiguousarray(P_lev), np.ascontiguousarray(Tbar))
    rho0 = EOS.fn_rho_kgm3(np.ascontiguousarray(P_lev), np.ascontiguousarray(Tbar))
    for k, v in zip(PROPS, (rho0, VP0, VS0, KS0, GS0)):
        out[k + '_ref'] = np.asarray(v)
        out[k + '_mean'] = LevelMean(out[k], weights)
    return out
