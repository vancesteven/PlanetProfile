"""
Lateral clathrate variation for 3D interior models.

Specifies spatially varying clathrate volume fraction f_clath(theta, phi)
and computes its effects on effective thermal conductivity and
self-consistent ice thickness.

Clathrate has k_therm ~ 0.5 W/m/K (constant), while ice Ih has
k_therm ~ 2-4 W/m/K (T-dependent). Higher clathrate fraction means
lower effective conductivity, steeper temperature gradient, and
thinner ice for a given basal heat flux.
"""
import numpy as np
import logging

log = logging.getLogger('PlanetProfile.Lateral.ClathrateLateral')

# Clathrate thermal conductivity (constant, from ClathrateProps.py)
K_CLATH_WmK = 0.5


def InitClathrateLateral(Planet):
    """ Initialize lateral clathrate variation from configuration.

        Reads clathrate fraction as SH coefficients or evaluates on grid.
        Clamps to [0, 1].

        Args:
            Planet: PlanetStruct with Lateral substruct configured.
                Requires either fClath_Cpq/fClath_Spq/fClath_pMax (SH mode)
                or fClath already set on Lateral.

        Returns:
            Planet: Updated with fClath field on grid.
    """
    from PlanetProfile.Lateral.SpatialGrid import SHtoGrid

    Lateral = Planet.Lateral

    if Lateral.fClath_Cpq is not None:
        # SH coefficient mode
        Lateral.fClath = SHtoGrid(
            Lateral.fClath_Cpq,
            Lateral.fClath_Spq,
            Lateral.fClath_pMax,
            Lateral.theta_rad,
            Lateral.phi_rad
        )
        log.info(f'Clathrate fraction from SH (pMax={Lateral.fClath_pMax})')

    elif Lateral.fClath is None:
        # Default: uniform from Bulk setting
        fClath0 = Planet.Bulk.clathMaxDepth_m / (Planet.zb_km * 1e3) \
            if hasattr(Planet.Bulk, 'clathMaxDepth_m') and Planet.Bulk.clathMaxDepth_m is not None \
            else 0.0
        Lateral.fClath = np.full(Lateral.nPix, fClath0)
        log.info(f'Uniform clathrate fraction: {fClath0:.3f}')

    # Clamp to [0, 1]
    Lateral.fClath = np.clip(Lateral.fClath, 0.0, 1.0)

    log.info(f'Clathrate fraction: mean={np.mean(Lateral.fClath):.3f}, '
             f'range=[{np.min(Lateral.fClath):.3f}, {np.max(Lateral.fClath):.3f}]')

    return Planet


def ComputeEffectiveConductivity(fClath, T_K_mean=200.0):
    """ Compute effective thermal conductivity for mixed ice-clathrate.

        Uses simple volume-weighted mixing (consistent with J-rule
        mixing in PlanetProfile's GetIceEOS).

        Args:
            fClath: Clathrate volume fraction (scalar or array).
            T_K_mean: Representative temperature for ice Ih conductivity.

        Returns:
            k_eff: Effective thermal conductivity in W/m/K.
    """
    # Ice Ih thermal conductivity (approximate, T-dependent)
    # Hobbs (1974): k_ice ≈ 651 / T for T in K, giving ~3.3 W/m/K at 200 K
    k_ice = 651.0 / T_K_mean

    k_eff = fClath * K_CLATH_WmK + (1.0 - fClath) * k_ice
    return k_eff


def EstimateIceThicknessFromClathrate(Planet, q_base_Wm2=None):
    """ Estimate self-consistent ice thickness from lateral clathrate variation.

        For a conductive shell with fixed basal heat flux:
            d_ice = k_eff * (T_melt - T_surf) / q_base

        Lower k_eff (more clathrate) -> thinner ice.

        Args:
            Planet: PlanetStruct with Lateral.fClath set.
            q_base_Wm2: Basal heat flux in W/m^2. If None, estimated from
                reference model.

        Returns:
            dIce_m: Estimated ice thickness at each grid point in meters.
    """
    Lateral = Planet.Lateral

    if q_base_Wm2 is None:
        # Estimate from reference model: q = k * dT / d
        k_ref = 651.0 / 200.0  # Approximate ice conductivity
        dT_ref = Planet.Bulk.Tb_K - Planet.Bulk.Tsurf_K
        d_ref = Planet.zb_km * 1e3
        q_base_Wm2 = k_ref * dT_ref / d_ref
        log.debug(f'Estimated basal heat flux: {q_base_Wm2 * 1e3:.1f} mW/m^2')

    # Temperature difference across ice shell
    dT = Planet.Bulk.Tb_K - Planet.Bulk.Tsurf_K

    # Mean temperature for conductivity estimate
    T_mean = (Planet.Bulk.Tb_K + Planet.Bulk.Tsurf_K) / 2

    # Effective conductivity at each grid point
    k_eff = ComputeEffectiveConductivity(Lateral.fClath, T_mean)
    Lateral.kThermEff_WmK = k_eff

    # Self-consistent ice thickness
    dIce_m = k_eff * dT / q_base_Wm2

    log.info(f'Estimated ice thickness from clathrate: '
             f'mean={np.mean(dIce_m)/1e3:.1f} km, '
             f'range=[{np.min(dIce_m)/1e3:.1f}, {np.max(dIce_m)/1e3:.1f}] km')

    return dIce_m
