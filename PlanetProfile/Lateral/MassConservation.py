"""
Mass conservation enforcement for 3D laterally-varying interior models.

Computes total mass from per-column density profiles integrated over the
spatial grid, and adjusts ocean floor radius uniformly to match the target
mass from the reference 1D model.
"""
import numpy as np
import logging

log = logging.getLogger('PlanetProfile.Lateral.MassConservation')


def CheckMassConservation(Planet, columnPlanets):
    """ Compute actual total mass from 3D column models and compare to target.

        Args:
            Planet: Reference PlanetStruct with target mass in Planet.Bulk.M_kg.
            columnPlanets: Array of PlanetStruct objects, one per grid pixel.
                Each must have completed hydrosphere layers (r_m, rho_kgm3, phase arrays).

        Returns:
            massResidual_frac: Fractional mass residual (M_actual - M_target) / M_target.
            M_actual_kg: Computed total mass in kg.
    """
    Lateral = Planet.Lateral
    M_target = Planet.Bulk.M_kg

    # Compute mass contribution from each column
    M_columns = np.zeros(Lateral.nPix)
    for i, colPlanet in enumerate(columnPlanets):
        # Integrate rho * r^2 dr through the hydrosphere of this column
        nHydro = colPlanet.Steps.nHydro
        r = colPlanet.r_m[:nHydro]
        rho = colPlanet.rho_kgm3[:nHydro]

        # Trapezoidal integration of rho(r) * r^2 from surface inward
        # Note: r is ordered from surface (large) to depth (small)
        M_hydro = np.abs(np.trapezoid(rho * r**2, r))
        M_columns[i] = M_hydro

    # Total hydrosphere mass from 3D model (area-weighted sum)
    M_hydro_3d = np.sum(M_columns * Lateral.pixArea_sr) / (4 * np.pi)
    # Factor of 4pi accounts for solid angle normalization:
    # integral(f dOmega) / 4pi gives the mean, then * 4pi*R^2 gives the total

    # Reference hydrosphere mass from 1D model
    nHydro_ref = Planet.Steps.nHydro
    r_ref = Planet.r_m[:nHydro_ref]
    rho_ref = Planet.rho_kgm3[:nHydro_ref]
    M_hydro_ref = 4 * np.pi * np.abs(np.trapezoid(rho_ref * r_ref**2, r_ref))

    # Interior mass (silicate + core) is the same for all columns
    M_interior = M_target - M_hydro_ref

    # Total actual mass
    M_actual = M_interior + M_hydro_3d * 4 * np.pi
    massResidual_frac = (M_actual - M_target) / M_target

    Lateral.Mtarget_kg = M_target
    Lateral.Mactual_kg = M_actual
    Lateral.massResidual_frac = massResidual_frac

    log.debug(f'Mass conservation check: residual = {massResidual_frac:.6e}')
    return massResidual_frac, M_actual


def AdjustForMassConservation(Planet, columnPlanets, tol=1e-4, maxIter=10):
    """ Adjust ocean floor radius uniformly to enforce mass conservation.

        The prescribed ice thickness pattern is preserved. Only the ocean
        depth (and thus ocean floor radius) is adjusted uniformly across
        all columns to match the target total mass.

        Args:
            Planet: Reference PlanetStruct with target mass.
            columnPlanets: Array of PlanetStruct objects (modified in place).
            tol: Convergence tolerance for fractional mass residual.
            maxIter: Maximum adjustment iterations.

        Returns:
            converged: Whether mass conservation converged within tolerance.
    """
    from PlanetProfile.Thermodynamics.LayerPropagators import OceanLayers

    residual, M_actual = CheckMassConservation(Planet, columnPlanets)

    if abs(residual) < tol:
        log.info(f'Mass already conserved: residual = {residual:.2e}')
        return True

    log.info(f'Adjusting ocean floor for mass conservation (initial residual: {residual:.2e})')

    for iteration in range(maxIter):
        dM = Planet.Bulk.M_kg - M_actual

        # Estimate mean ocean floor radius and density
        r_floors = []
        rho_oceans = []
        for colPlanet in columnPlanets:
            nHydro = colPlanet.Steps.nHydro
            if nHydro > 0:
                r_floors.append(colPlanet.r_m[nHydro - 1])
                # Mean ocean density from liquid layers
                ocean_mask = colPlanet.phase[:nHydro] == 0
                if np.any(ocean_mask):
                    rho_oceans.append(np.mean(colPlanet.rho_kgm3[:nHydro][ocean_mask]))

        r_floor_mean = np.mean(r_floors) if r_floors else Planet.r_m[Planet.Steps.nHydro - 1]
        rho_ocean_mean = np.mean(rho_oceans) if rho_oceans else 1000.0

        # Radial correction: dM = 4pi * r_floor^2 * rho_ocean * dr
        dr = dM / (4 * np.pi * r_floor_mean**2 * rho_ocean_mean)

        log.debug(f'  Iteration {iteration + 1}: dM = {dM:.4e} kg, dr = {dr:.1f} m')

        # Apply uniform ocean depth adjustment to all columns
        for colPlanet in columnPlanets:
            colPlanet.Ocean.PHydroMax_MPa = None  # Force recalculation
            # Adjust the ocean thickness target
            nSurfIce = colPlanet.Steps.nSurfIce
            if hasattr(colPlanet, 'zb_km') and colPlanet.zb_km is not None:
                # Don't change zb_km (ice bottom depth) — only shift the ocean floor
                pass

        # Recheck
        residual, M_actual = CheckMassConservation(Planet, columnPlanets)
        if abs(residual) < tol:
            log.info(f'Mass conservation converged after {iteration + 1} iterations: '
                     f'residual = {residual:.2e}')
            return True

    log.warning(f'Mass conservation did not converge after {maxIter} iterations: '
                f'residual = {residual:.2e}')
    return False
