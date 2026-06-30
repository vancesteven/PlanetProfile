"""
Mass conservation enforcement for 3D laterally-varying interior models.

Computes total mass from per-column density profiles integrated over the
spatial grid, and adjusts a common seafloor pressure to match the target
mass from the reference 1D model.
"""
import numpy as np
import logging

log = logging.getLogger('PlanetProfile.Lateral.MassConservation')


def _HydrosphereMass_kg(planet):
    """Return the full-sphere hydrosphere mass represented by one profile."""
    nHydro = int(planet.Steps.nHydro)
    layerMasses_kg = getattr(planet, 'MLayer_kg', None)
    if layerMasses_kg is not None and np.size(layerMasses_kg) >= nHydro:
        layerMasses_kg = np.asarray(layerMasses_kg[:nHydro], dtype=float)
        if np.all(np.isfinite(layerMasses_kg)) and np.any(layerMasses_kg > 0):
            return float(np.sum(layerMasses_kg))

    # Fallback for reduced or synthetic profiles that do not carry MLayer_kg.
    # The 4*pi factor converts the radial integral to a full-sphere mass.
    r_m = np.asarray(planet.r_m[:nHydro], dtype=float)
    rho_kgm3 = np.asarray(planet.rho_kgm3[:nHydro], dtype=float)
    try:
        radialIntegral = np.trapezoid(rho_kgm3 * r_m**2, r_m)
    except AttributeError:
        radialIntegral = np.trapz(rho_kgm3 * r_m**2, r_m)
    return float(4 * np.pi * np.abs(radialIntegral))


def EstimateOceanFloorPressureStep_MPa(dM_kg, rFloor_m, gFloor_ms2,
                                       damping=0.8, maxStep_MPa=2.0):
    """Estimate a stable common seafloor-pressure correction from a mass error.

    Hydrostatic column mass obeys dM/dP approximately equal to 4*pi*r^2/g.
    A negative ``dM_kg`` therefore lowers the seafloor pressure and removes
    hydrosphere mass, while a positive value deepens the ocean.
    """
    if rFloor_m <= 0 or gFloor_ms2 <= 0:
        raise ValueError('Ocean-floor radius and gravity must be positive')
    if not 0 < damping <= 1:
        raise ValueError('Mass-conservation damping must be in (0, 1]')

    dP_MPa = damping * dM_kg * gFloor_ms2 / (4 * np.pi * rFloor_m**2) * 1e-6
    return float(np.clip(dP_MPa, -maxStep_MPa, maxStep_MPa))


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

    M_columns = np.zeros(Lateral.nPix)
    for i, colPlanet in enumerate(columnPlanets):
        M_columns[i] = _HydrosphereMass_kg(colPlanet)

    pixelAreas_sr = np.asarray(Lateral.pixArea_sr, dtype=float)
    if pixelAreas_sr.size == 1:
        pixelAreas_sr = np.full(M_columns.size, pixelAreas_sr.item())
    if M_columns.size != pixelAreas_sr.size:
        raise ValueError(
            f'Mass conservation received {M_columns.size} columns for '
            f'{pixelAreas_sr.size} lateral pixel areas'
        )
    pixelWeights = pixelAreas_sr / (4 * np.pi)
    M_hydro_3d = float(np.sum(M_columns * pixelWeights))
    M_hydro_ref = _HydrosphereMass_kg(Planet)

    M_interior = M_target - M_hydro_ref

    M_actual = M_interior + M_hydro_3d
    massResidual_frac = (M_actual - M_target) / M_target

    Lateral.Mtarget_kg = M_target
    Lateral.Mactual_kg = M_actual
    Lateral.massResidual_frac = massResidual_frac

    log.debug(f'Mass conservation check: residual = {massResidual_frac:.6e}')
    return massResidual_frac, M_actual


def AdjustForMassConservation(Planet, columnPlanets, Params=None, tol=1e-4,
                              maxIter=10):
    """Adjust the common ocean-floor pressure to enforce mass conservation.

        The prescribed ice thickness pattern is preserved. Only the ocean
        depth is adjusted through a common seafloor pressure across all
        columns to match the target total mass.

        Args:
            Planet: Reference PlanetStruct with target mass.
            columnPlanets: Array of PlanetStruct objects (modified in place).
            Params: ParamsStruct used to regenerate the hydrosphere columns.
            tol: Convergence tolerance for fractional mass residual.
            maxIter: Maximum adjustment iterations.

        Returns:
            converged: Whether mass conservation converged within tolerance.
    """
    from PlanetProfile.Lateral.LateralStructure import RunLateralColumns

    residual, M_actual = CheckMassConservation(Planet, columnPlanets)

    if abs(residual) < tol:
        log.info(f'Mass already conserved: residual = {residual:.2e}')
        return True

    if Params is None:
        raise ValueError(
            'Params is required to regenerate lateral columns during mass conservation'
        )

    log.info(f'Adjusting ocean floor for mass conservation (initial residual: {residual:.2e})')

    for iteration in range(maxIter):
        dM = Planet.Bulk.M_kg - M_actual

        r_floors = []
        g_floors = []
        for colPlanet in columnPlanets:
            nHydro = int(colPlanet.Steps.nHydro)
            if nHydro > 0:
                iFloor = min(nHydro, len(colPlanet.r_m) - 1)
                r_floors.append(colPlanet.r_m[iFloor])
                g_floors.append(colPlanet.g_ms2[iFloor])

        r_floor_mean = np.mean(r_floors) if r_floors else Planet.r_m[Planet.Steps.nHydro - 1]
        g_floor_mean = np.mean(g_floors) if g_floors else Planet.g_ms2[Planet.Steps.nHydro - 1]
        dP_MPa = EstimateOceanFloorPressureStep_MPa(
            dM, r_floor_mean, g_floor_mean,
            damping=getattr(Planet.Lateral, 'massConservationDamping', 0.8),
        )

        currentP_MPa = getattr(Planet.Lateral, 'referenceHydroBottom_MPa', None)
        if currentP_MPa is None:
            currentP_MPa = float(Planet.P_MPa[Planet.Steps.nHydro - 1])
        minP_MPa = max(float(colPlanet.Pb_MPa) for colPlanet in columnPlanets)
        minP_MPa += 2 * float(Planet.Ocean.deltaP)
        nextP_MPa = max(currentP_MPa + dP_MPa, minP_MPa)

        if np.isclose(nextP_MPa, currentP_MPa, atol=1e-12):
            log.warning(
                'Mass conservation cannot reduce the ocean further without '
                'removing liquid layers beneath the thickest ice column'
            )
            break

        Planet.Lateral.referenceHydroBottom_MPa = nextP_MPa
        log.info(
            f'  Iteration {iteration + 1}: dM={dM:.4e} kg, '
            f'ocean-floor pressure {currentP_MPa:.4f} -> {nextP_MPa:.4f} MPa'
        )

        updatedColumns = RunLateralColumns(
            Planet, Params, checkpoint_interval=0, resume_from_checkpoint=False
        )
        columnPlanets[:] = updatedColumns

        residual, M_actual = CheckMassConservation(Planet, columnPlanets)
        Planet.Lateral.massConservationIterations = iteration + 1
        if abs(residual) < tol:
            log.info(f'Mass conservation converged after {iteration + 1} iterations: '
                     f'residual = {residual:.2e}')
            return True

    log.warning(f'Mass conservation did not converge after {maxIter} iterations: '
                f'residual = {residual:.2e}')
    return False
