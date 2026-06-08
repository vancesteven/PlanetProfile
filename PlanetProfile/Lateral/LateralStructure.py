"""
3D laterally-varying interior structure driver for PlanetProfile.

Orchestrates per-column hydrosphere calculations across a spatial grid,
with support for:
- Prescribed ice thickness from SH coefficients or callable functions
- Mass conservation enforcement
- Feed of ice-ocean boundary topography to MoonMag for asymmetric induction

Literature:
- Ojakangas & Stevenson (1989): Tidal heating with eccentricity forcing
- Tobie et al. (2005): Geographic tidal dissipation in Titan
"""
import numpy as np
import logging
from copy import deepcopy

log = logging.getLogger('PlanetProfile.Lateral.LateralStructure')


def InitLateralStructure(Planet, Params):
    """ Initialize 3D lateral structure from configuration.

        Reads ice thickness as SH coefficients or evaluates a callable,
        sets up the spatial grid, and prepares for column model runs.

        Args:
            Planet: PlanetStruct with Lateral substruct configured.
                Required: Either dIce_Cpq_km/dIce_Spq_km/dIce_pMax (SH mode)
                or dIce_func (callable mode) must be set.
            Params: ParamsStruct.

        Returns:
            Planet: Updated PlanetStruct with grid and ice thickness field.
    """
    from PlanetProfile.Lateral.SpatialGrid import InitGrid, SHtoGrid

    Lateral = Planet.Lateral

    InitGrid(Lateral)
    log.info(f'Lateral grid initialized: {Lateral.gridType} with {Lateral.nPix} pixels')

    if Lateral.dIce_Cpq_km is not None:
        Lateral.dIce_m = SHtoGrid(
            Lateral.dIce_Cpq_km * 1e3,
            Lateral.dIce_Spq_km * 1e3,
            Lateral.dIce_pMax,
            Lateral.theta_rad,
            Lateral.phi_rad
        )
        log.info(f'Ice thickness from SH (pMax={Lateral.dIce_pMax}): '
                 f'mean={np.mean(Lateral.dIce_m)/1e3:.1f} km, '
                 f'range=[{np.min(Lateral.dIce_m)/1e3:.1f}, {np.max(Lateral.dIce_m)/1e3:.1f}] km')

    elif Lateral.dIce_func is not None:
        Lateral.dIce_m = np.array([Lateral.dIce_func(theta) for theta in Lateral.theta_rad])
        log.info(f'Ice thickness from function: '
                 f'mean={np.mean(Lateral.dIce_m)/1e3:.1f} km, '
                 f'range=[{np.min(Lateral.dIce_m)/1e3:.1f}, {np.max(Lateral.dIce_m)/1e3:.1f}] km')
    else:
        Lateral.dIce_m = np.full(Lateral.nPix, Planet.zb_km * 1e3)
        log.info(f'Uniform ice thickness from reference: {Planet.zb_km:.1f} km')

    return Planet


def RunLateralColumns(Planet, Params):
    """ Run hydrosphere calculations for each grid column.

        Each column gets the reference Planet's interior (silicate/core)
        but a locally varying ice thickness. Only ice and ocean layers
        are recalculated per column.

        Ice thickness variation is achieved by adjusting Tb_K for each
        column: the reference model's P(z) profile gives the pressure at
        the desired ice-ocean boundary depth, and the ocean freezing curve
        gives the corresponding Tb_K at that pressure.

        Args:
            Planet: Reference PlanetStruct (must have completed full 1D model).
            Params: ParamsStruct.

        Returns:
            columnPlanets: Array of PlanetStruct, one per pixel (nPix,).
    """
    from PlanetProfile.Main import HydroOnly
    from PlanetProfile.Thermodynamics.HydroEOS import GetTfreeze, GetOceanEOS

    Lateral = Planet.Lateral
    nPix = Lateral.nPix

    log.info(f'Running {nPix} lateral columns...')

    nSurf = Planet.Steps.nSurfIce
    z_ref_m = Planet.z_m[:nSurf+1]
    P_ref_MPa = Planet.P_MPa[:nSurf+1]

    Pb_min = np.interp(np.min(Lateral.dIce_m), z_ref_m, P_ref_MPa)
    Pb_max = np.interp(np.max(Lateral.dIce_m), z_ref_m, P_ref_MPa)
    Pmelt_wide = np.arange(max(0.01, Pb_min - 5), Pb_max + 5, Planet.PfreezeRes_MPa)
    Tmelt_wide = np.arange(Planet.TfreezeLower_K, Planet.Bulk.Tb_K + 10, Planet.TfreezeRes_K)
    lateralMeltEOS = GetOceanEOS(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt,
                                  Pmelt_wide, Tmelt_wide, None,
                                  phaseType=Planet.Ocean.phaseType,
                                  FORCE_NEW=True, MELT=True)
    log.info(f'Lateral meltEOS: P=[{Pmelt_wide[0]:.1f}, {Pmelt_wide[-1]:.1f}] MPa, '
             f'T=[{Tmelt_wide[0]:.1f}, {Tmelt_wide[-1]:.1f}] K')

    columnPlanets = np.empty(nPix, dtype=object)
    nFailed = 0
    for i in range(nPix):
        colPlanet = deepcopy(Planet)

        d_ice_km = Lateral.dIce_m[i] / 1e3
        d_ice_m = d_ice_km * 1e3
        if d_ice_m <= z_ref_m[-1] and d_ice_m >= z_ref_m[0]:
            Pb_col_MPa = np.interp(d_ice_m, z_ref_m, P_ref_MPa)
        else:
            dPdz = (P_ref_MPa[-1] - P_ref_MPa[0]) / (z_ref_m[-1] - z_ref_m[0])
            Pb_col_MPa = P_ref_MPa[0] + dPdz * d_ice_m

        try:
            Tb_col_K = GetTfreeze(lateralMeltEOS, Pb_col_MPa,
                                  Planet.TfreezeLower_K, TRes_K=Planet.TfreezeRes_K)
            colPlanet.Bulk.Tb_K = Tb_col_K
        except Exception as e:
            log.warning(f'Column {i}: could not compute Tb_K for d_ice={d_ice_km:.1f} km, '
                       f'Pb={Pb_col_MPa:.1f} MPa: {e}')
            colPlanet.invalidReason = str(e)
            nFailed += 1
            columnPlanets[i] = colPlanet
            continue

        if Lateral.DO_CLATH_LATERAL and Lateral.fClath is not None:
            colPlanet.Bulk.clathType = 'bottom'
            colPlanet.Bulk.clathMaxDepth_m = d_ice_m
            colPlanet.Bulk.volumeFractionClathrate = Lateral.fClath[i]

        if Lateral.DO_TIDAL_3D and Lateral.HtidalIce_Wm3 is not None:
            colPlanet.Bulk.Htidal_Wm3 = Lateral.HtidalIce_Wm3[i]

        colPlanet.index = i
        columnPlanets[i] = colPlanet

    if nFailed > 0:
        log.warning(f'{nFailed}/{nPix} columns failed Tb_K computation')

    import time
    Params.nModels = nPix
    Params.tStart_s = time.time()

    if Params.DO_PARALLEL:
        columnPlanets = _RunColumnsParallel(columnPlanets, Params)
    else:
        columnPlanets = _RunColumnsSerial(columnPlanets, Params)

    _ExtractColumnSummaries(Lateral, columnPlanets)

    log.info(f'Lateral columns complete. Tb range: '
             f'[{np.nanmin(Lateral.Tb_K):.1f}, {np.nanmax(Lateral.Tb_K):.1f}] K')

    return columnPlanets


def _ColumnFailed(colPlanet):
    """ Check if a column was marked as failed during Tb_K computation.
        invalidReason is None (not yet run) or 'Valid' (successful) for good columns,
        and a descriptive error string for failed columns.
    """
    return (hasattr(colPlanet, 'invalidReason')
            and colPlanet.invalidReason is not None
            and colPlanet.invalidReason != 'Valid')


def _RunColumnsSerial(columnPlanets, Params):
    """ Run column models in serial. """
    from PlanetProfile.Main import HydroOnly

    for i, colPlanet in enumerate(columnPlanets):
        if _ColumnFailed(colPlanet):
            continue
        try:
            columnPlanets[i], _ = HydroOnly(colPlanet, Params)
        except Exception as e:
            log.warning(f'Column {i} failed: {e}')
            columnPlanets[i].invalidReason = str(e)
            continue

        if (i + 1) % 100 == 0:
            log.info(f'  Completed {i + 1}/{len(columnPlanets)} columns')

    return columnPlanets


def _RunColumnsParallel(columnPlanets, Params):
    """ Run column models in parallel using multiprocessing. """
    import multiprocessing as mp
    import platform

    from PlanetProfile.Main import HydroOnly

    nPix = len(columnPlanets)
    nCores = min(Params.maxCores, nPix)

    if platform.system() == 'Windows':
        mtpContext = mp.get_context('spawn')
    else:
        mtpContext = mp.get_context('fork')

    log.info(f'Running {nPix} columns on {nCores} cores...')

    pool = mtpContext.Pool(nCores)
    results = [pool.apply_async(HydroOnly, (deepcopy(col), deepcopy(Params)))
               for col in columnPlanets]
    pool.close()
    pool.join()

    for i, result in enumerate(results):
        try:
            columnPlanets[i] = result.get()[0]
        except Exception as e:
            log.warning(f'Column {i} failed: {e}')
            columnPlanets[i].invalidReason = str(e)

    return columnPlanets


def _ExtractColumnSummaries(Lateral, columnPlanets):
    """ Extract summary fields from column models into Lateral arrays. """
    nPix = len(columnPlanets)
    Lateral.Tb_K = np.full(nPix, np.nan)
    Lateral.qSurf_Wm2 = np.full(nPix, np.nan)
    Lateral.sigma_mean_Sm = np.full(nPix, np.nan)

    for i, colPlanet in enumerate(columnPlanets):
        if _ColumnFailed(colPlanet):
            continue

        nSurf = colPlanet.Steps.nSurfIce
        if nSurf > 0 and nSurf < len(colPlanet.T_K):
            Lateral.Tb_K[i] = colPlanet.T_K[nSurf - 1]

        if hasattr(colPlanet, 'qSurf_Wm2') and colPlanet.qSurf_Wm2 is not None:
            Lateral.qSurf_Wm2[i] = colPlanet.qSurf_Wm2

        nHydro = colPlanet.Steps.nHydro
        ocean_mask = colPlanet.phase[:nHydro] == 0
        if np.any(ocean_mask) and hasattr(colPlanet.Ocean, 'EOS') and colPlanet.Ocean.EOS is not None:
            try:
                P_ocean = colPlanet.P_MPa[:nHydro][ocean_mask]
                T_ocean = colPlanet.T_K[:nHydro][ocean_mask]
                sigma_ocean = colPlanet.Ocean.EOS.fn_sigma_Sm(P_ocean, T_ocean)
                Lateral.sigma_mean_Sm[i] = np.mean(sigma_ocean)
            except Exception:
                pass


def UpdateAsymShapeFrom3D(Planet, columnPlanets):
    """ Convert 3D ice-ocean boundary topography to SH coefficients
        for MoonMag asymmetric induction calculations.

        Args:
            Planet: Reference PlanetStruct (modified in place).
            columnPlanets: Array of PlanetStruct with completed hydrosphere.

        Returns:
            Planet: Updated with asymShape_m for MoonMag.
    """
    from PlanetProfile.Lateral.SpatialGrid import GridToSH
    try:
        from PlanetProfile.MagneticInduction.MoonMag.asymmetry_funcs import get_chipq_from_CSpq_single
    except ModuleNotFoundError:
        from MoonMag.asymmetry_funcs import get_chipq_from_CSpq_single

    Lateral = Planet.Lateral

    r_iceOcean = np.zeros(Lateral.nPix)
    for i, colPlanet in enumerate(columnPlanets):
        if _ColumnFailed(colPlanet):
            r_iceOcean[i] = np.nan
            continue
        nSurf = colPlanet.Steps.nSurfIce
        if nSurf > 0 and nSurf < len(colPlanet.r_m):
            r_iceOcean[i] = colPlanet.r_m[nSurf]
        else:
            r_iceOcean[i] = Planet.Bulk.R_m - Lateral.dIce_m[i]

    r_mean = np.nanmean(r_iceOcean)
    dr = r_iceOcean - r_mean
    dr[np.isnan(dr)] = 0.0

    pMax = Lateral.dIce_pMax if Lateral.dIce_pMax is not None else 4
    Cpq_m, Spq_m = GridToSH(dr, Lateral.theta_rad, Lateral.phi_rad,
                             Lateral.pixArea_sr, pMax)

    nBds = Planet.Magnetic.nBds if hasattr(Planet.Magnetic, 'nBds') and Planet.Magnetic.nBds is not None else 1
    Planet.Magnetic.pMax = pMax
    Planet.Magnetic.asymShape_m = np.zeros((nBds, 2, pMax + 1, pMax + 1), dtype=np.complex128)

    for p in range(1, pMax + 1):
        # get_chipq_from_CSpq_single expects arrays of length p+1
        chipq_p = get_chipq_from_CSpq_single(p, Cpq_m[p, :p+1], Spq_m[p, :p+1])
        for q in range(p + 1):
            for iLayer in range(nBds):
                if hasattr(Planet.Magnetic, 'rSigChange_m') and Planet.Magnetic.rSigChange_m is not None:
                    rScale = Planet.Magnetic.rSigChange_m[iLayer] / r_mean
                else:
                    rScale = 1.0
                # chipq_p shape is (2, p+1): [pos/neg q, q index]
                Planet.Magnetic.asymShape_m[iLayer, :, p, q] = chipq_p[:, q] * rScale

    log.info(f'Asymmetric shape updated from 3D model: pMax={pMax}, nBds={nBds}')
    return Planet


def RunLateral3D(Planet, Params):
    """ Full 3D lateral structure pipeline.

        1. Initialize grid and ice thickness field
        2. (Optional) Initialize lateral clathrate -> self-consistent d_ice
        3. Run lateral column models
        4. (Optional) 3D tidal heating with feedback convergence
        5. Enforce mass conservation
        6. Update asymmetric shape for MoonMag
        7. (Optional) Run MoonMag asymmetric induction
        8. (Optional) Run PyALMA3 for multi-degree Love numbers
        9. Save 3D fields

        Args:
            Planet: Reference PlanetStruct (must have completed full 1D model).
            Params: ParamsStruct.

        Returns:
            Planet: Updated with 3D lateral structure results.
    """
    log.info(f'Starting 3D lateral structure calculation for {Planet.name}')

    Planet = InitLateralStructure(Planet, Params)

    if Planet.Lateral.DO_CLATH_LATERAL:
        from PlanetProfile.Lateral.ClathrateLateral import InitClathrateLateral, EstimateIceThicknessFromClathrate
        Planet = InitClathrateLateral(Planet)
        Planet.Lateral.dIce_m = EstimateIceThicknessFromClathrate(Planet)

    columnPlanets = RunLateralColumns(Planet, Params)

    if Planet.Lateral.DO_TIDAL_3D:
        from PlanetProfile.Lateral.TidalHeating3D import ComputeTidalHeating3D, ConvergeTidalFeedback
        Planet = ComputeTidalHeating3D(Planet, Params, columnPlanets)
        Planet, columnPlanets = ConvergeTidalFeedback(Planet, Params, columnPlanets)

    if Planet.Lateral.DO_MASS_CONSERVE:
        from PlanetProfile.Lateral.MassConservation import CheckMassConservation, AdjustForMassConservation
        residual, _ = CheckMassConservation(Planet, columnPlanets)
        if abs(residual) > 1e-4:
            AdjustForMassConservation(Planet, columnPlanets)

    if Params.CALC_NEW_INDUCT and hasattr(Planet.Magnetic, 'CALC_ASYM') and Planet.Magnetic.CALC_ASYM:
        Planet = UpdateAsymShapeFrom3D(Planet, columnPlanets)

    if Params.CALC_NEW_INDUCT and hasattr(Planet.Magnetic, 'CALC_ASYM') and Planet.Magnetic.CALC_ASYM:
        from PlanetProfile.MagneticInduction.MagneticInduction import MagneticInduction
        Planet, Params = MagneticInduction(Planet, Params)
        log.info('Asymmetric induction calculation complete')

    if Params.CALC_NEW_GRAVITY and len(Params.Gravity.harmonic_degrees) > 1:
        from PlanetProfile.Gravity.Gravity import GravityParameters
        Planet, Params = GravityParameters(Planet, Params)
        log.info(f'Love numbers computed for degrees: {Params.Gravity.harmonic_degrees}')

    from PlanetProfile.Lateral.LateralIO import SaveLateralResults
    SaveLateralResults(Planet, Params)

    # Generate lateral plots if requested
    if Params.PLOT_LATERAL:
        from PlanetProfile.Plotting.LateralPlots import GenerateLateralPlots
        GenerateLateralPlots(Planet, Params)

    log.info('3D lateral structure calculation complete')
    return Planet, Params
