import numpy as np
import logging
from scipy.signal import savgol_filter
from copy import deepcopy
from scipy.optimize import root_scalar as GetZero

from PlanetProfile.Thermodynamics.IronCore import IronCoreLayers
from PlanetProfile.Thermodynamics.HydroEOS import GetPfreeze, GetTfreeze, \
    GetIceArrheniusViscosityKwargs, GetIceEOS, GetOceanEOS
from PlanetProfile.Utilities.Indexing import PhaseConv, GetPhaseIndices
from PlanetProfile.Thermodynamics.InnerEOS import GetHtidalFunc, GetphiCalc, GetInnerEOS
from PlanetProfile.Thermodynamics.Silicates import SilicateLayers
from PlanetProfile.Thermodynamics.ThermalProfiles.Convection import IceIConvectSolid, IceIConvectPorous, \
    IceIIIConvectSolid, IceIIIConvectPorous, IceVConvectSolid, IceVConvectPorous, \
    ClathShellConvectSolid, ClathShellConvectPorous
from PlanetProfile.Thermodynamics.ThermalProfiles.IceConduction import IceIWholeConductSolid, IceIWholeConductPorous, \
    IceIConductClathLidSolid, IceIConductClathLidPorous, IceIConductClathUnderplateSolid, IceIConductClathUnderplatePorous, \
    IceIIIConductSolid, IceIIIConductPorous, IceVConductSolid, IceVConductPorous
from PlanetProfile.Thermodynamics.ThermalProfiles.ThermalProfiles import ConvectionDeschampsSotin2001, \
    ConvectionKalousova2018, GetRaCrit
from PlanetProfile.Thermodynamics.Geophysical import PropogateConductionFromDepth
from PlanetProfile.Utilities.defineStructs import Constants, EOSlist, HPIcePhaseState, Timing, ResolveHPIceConvectionModel, \
    HP_ICE_PRODUCTION_ENERGY_ABS_FLOOR_W, HP_ICE_PRODUCTION_ENERGY_FRAC_TOL, \
    HP_ICE_PRODUCTION_LAYER_CLOSURE_ABS_TOL_M, HP_ICE_PRODUCTION_LAYER_CLOSURE_FRAC_TOL, \
    HP_ICE_PRODUCTION_MIN_THICKNESS_M, HP_ICE_PRODUCTION_PHASE_BOUNDARY_TOL_K
import time

# Assign logger
log = logging.getLogger('PlanetProfile')

POSTHOC_EOS_PRESSURE_MARGIN_WARN_MPA = 5.0
POSTHOC_EOS_TEMPERATURE_MARGIN_WARN_K = 1.0
POSTHOC_PHASE_BOUNDARY_MARGIN_WARN_K = 0.1
POSTHOC_RAYLEIGH_NEAR_CRITICAL_RATIO = 1.1
ACTIVE_ICE_VI_HEAT_FLUX_FRAC_TOL = 1e-10
ACTIVE_ICE_VI_HEAT_FLUX_ABS_FLOOR_WM2 = 1e-12
ACTIVE_ICE_VI_ROLLBACK_POLICY_VERSION = "active_ice_vi_candidate_rollback_v1"
ACTIVE_ICE_VI_THERMAL_UPDATE_STRATEGY_NO_OP = "no_op_reconstruction"

def IceLayers(Planet, Params):
    """ Wrapper function for ice layer propogation. Decides between self consistent and non-self consistent modeling.
    """
    Timing.setFunctionTime(time.time())
    # Early branching for non-self-consistent modeling
    if Planet.Do.NON_SELF_CONSISTENT:
        Planet, Params =  NonSelfConsistentIceLayer(Planet, Params)
    else:
        Planet, Params = SelfConsistentIceLayer(Planet, Params)
    Timing.printFunctionTimeDifference('IceLayers()', time.time())
    return Planet

def SelfConsistentIceLayer(Planet, Params):
    """ Layer propagation from surface downward through the ice using geophysics.
        Iteratively sets up the thermal profile (the density and temperature)
        of the layer with each pressure step for all ice layers including
        ice Ih, ice III, ice V by calling the necessary subfunctions

        Assigns Planet attributes:
            Steps.nSurfIce, phase, r_m, z_m, g_ms2, T_K, P_MPa, rho_kgm3, Cp_JkgK, alpha_pK, MLayer_kg, PbI_MPa, Pb_MPa
    """

    if Planet.Do.PARTIAL_DIFFERENTIATION:
        log.debug('Skipping ice layer in partially differentiated body. A future update will ' +
                  'include an option for whether to include an ice shell atop a mixed rock-and-' +
                  'ice interior. The current implementation includes only variable mixing in an ' +
                  'entirely rock+ice body.')
    else:
        if Planet.Do.ICEIh_THICKNESS:
            Planet = GetIceShellTFreeze(Planet, Params)
        else:
            Planet.Steps.nIbottom = Planet.Steps.nClath + Planet.Steps.nIceI
            Planet.Steps.nIIIbottom = Planet.Steps.nIbottom + Planet.Steps.nIceIIILitho
            Planet.Steps.nSurfIce = Planet.Steps.nIIIbottom + Planet.Steps.nIceVLitho
            # Assign phase values for near-surface ices
            Planet.phase[:Planet.Steps.nIbottom] = 1  # Ice Ih layers (some will be reassigned if Do.CLATHRATE = True)
            Planet.phase[Planet.Steps.nIbottom:Planet.Steps.nIIIbottom] = 3  # Ice III layers
            Planet.phase[Planet.Steps.nIIIbottom:Planet.Steps.nSurfIce] = 5  # Ice V layers
            # Finally, we need to set up the logical array for ice convection layers (this will be default False and overriden if we do ice convection))
            Planet.Steps.iConv = np.zeros(Planet.Steps.nSurfIce, dtype=bool)
            

            # Get the pressure consistent with the bottom of the surface ice layer that is
            # consistent with the choice of Tb_K we suppose for this model
            Planet.PbI_MPa = GetPfreeze(Planet.Ocean.meltEOS, 1, Planet.Bulk.Tb_K,
                                            PLower_MPa=Planet.PfreezeLower_MPa, PUpper_MPa=Planet.PfreezeUpper_MPa,
                                            PRes_MPa=Planet.PfreezeRes_MPa, UNDERPLATE=(Planet.Do.BOTTOM_ICEIII or Planet.Do.BOTTOM_ICEV), HPNOOCEAN=Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES,
                                            ALLOW_BROKEN_MODELS=Params.ALLOW_BROKEN_MODELS, DO_EXPLOREOGRAM=Params.DO_EXPLOREOGRAM)
            if(Planet.Do.CLATHRATE and
                    (Planet.Bulk.clathType == 'bottom' or
                    Planet.Bulk.clathType == 'whole')):
                PbClath_MPa = Planet.Ocean.ClathDissoc.PbClath_MPa()
                if not np.isnan(PbClath_MPa):
                    log.debug(f'Clathrate dissociation pressure: {PbClath_MPa:.3f} MPa.')
                    if PbClath_MPa < Planet.PbI_MPa:
                        raise ValueError('Dissociation pressure for clathrates is lower than the ice Ih ' +
                                        'melting pressure consistent with Bulk.Tb_K. This means ice Ih layers ' +
                                        'will be found underneath the clathrate layers, inconsistent with the ' +
                                        'assumption that clathrates are in contact with the ocean. Increase ' +
                                        f'Bulk.Tb_K until PbClath_MPa ({PbClath_MPa:.3f} MPa) ' +
                                        f'exceeds PbI_MPa ({Planet.PbI_MPa:.3f} MPa).')
                Planet.PbI_MPa = PbClath_MPa
            else:
                if np.isnan(Planet.PbI_MPa):
                    msg = f'No valid phase transition was found for Tb_K = {Planet.Bulk.Tb_K:.3f} K for P in the range ' + \
                        f'[{Planet.PfreezeLower_MPa:.1f} MPa, {Planet.PfreezeUpper_MPa:.1f} MPa]. ' + \
                        'This likely means Tb_K is too high and the phase at the lower end of this range matches ' + \
                        'the phase at the upper end. Try decreasing Tb_K or increasing Planet.PfreezeUpper_MPa. ' + \
                        'For this model, the ice shell will be set to zero thickness.'
                    if (not Params.DO_EXPLOREOGRAM) and (not Params.DO_INDUCTOGRAM):
                        if Planet.Bulk.Tb_K > 271:
                            log.warning(msg)
                        else:
                            if Params.ALLOW_BROKEN_MODELS:
                                Planet.Do.VALID = False
                                Planet.invalidReason = f'No valid phase transition was found for Tb_K = {Planet.Bulk.Tb_K:.3f} K for P in the range ' + \
                                                    f'[{Planet.PfreezeLower_MPa:.1f} MPa, {Planet.PfreezeUpper_MPa:.1f} MPa]. '
                            else:
                                raise ValueError(msg)
                    Planet.PbI_MPa = 0.0
                log.debug(f'Ice Ih transition pressure: {Planet.PbI_MPa:.3f} MPa.')
            if Planet.PbI_MPa > 0:
                # Now do the same for HP ices, if present, to make sure we have a possible configuration before continuing
                if Planet.Do.BOTTOM_ICEV:
                    Planet.PbIII_MPa = GetPfreeze(Planet.Ocean.meltEOS, 3, Planet.Bulk.TbIII_K,
                                PLower_MPa=Planet.PbI_MPa, PUpper_MPa=Planet.Ocean.PHydroMax_MPa,
                                PRes_MPa=Planet.PfreezeRes_MPa, UNDERPLATE=True, HPNOOCEAN=False,
                                ALLOW_BROKEN_MODELS=Params.ALLOW_BROKEN_MODELS, DO_EXPLOREOGRAM=Params.DO_EXPLOREOGRAM)
                    if(Planet.PbIII_MPa <= Planet.PbI_MPa) or np.isnan(Planet.PbIII_MPa):
                        msg = 'Ice III bottom pressure is not greater than ice I bottom pressure. ' + \
                                'This likely indicates TbIII_K is too high for the corresponding Tb_K.' + \
                                f'\nPbI_MPa = {Planet.PbI_MPa:.3f}' + \
                                f', Tb_K = {Planet.Bulk.Tb_K:.3f}' + \
                                f'\nPbIII_MPa = {Planet.PbIII_MPa:.3f}' + \
                                f', TbIII_K = {Planet.Bulk.TbIII_K:.3f}'
                        if Params.ALLOW_BROKEN_MODELS:
                            Planet.PbIII_MPa = np.nan
                            if Params.DO_EXPLOREOGRAM:
                                log.info(msg)
                            else:
                                log.error(msg)
                            Planet.Do.VALID = False
                            Planet.invalidReason = 'TbIII_K is too high compared to Tb_K'
                        else:
                            raise ValueError(msg)
                    if not np.isnan(Planet.PbIII_MPa):
                        Planet.PbV_MPa = GetPfreeze(Planet.Ocean.meltEOS, 5, Planet.Bulk.TbV_K,
                                                        PLower_MPa=Planet.PbIII_MPa, PUpper_MPa=Planet.Ocean.PHydroMax_MPa,
                                                        PRes_MPa=Planet.PfreezeRes_MPa, UNDERPLATE=False, HPNOOCEAN=False,
                                                        ALLOW_BROKEN_MODELS=Params.ALLOW_BROKEN_MODELS,
                                                        DO_EXPLOREOGRAM=Params.DO_EXPLOREOGRAM)
                    else:
                        Planet.PbV_MPa = np.nan
                    Planet.Pb_MPa = Planet.PbV_MPa
                    if(Planet.PbV_MPa <= Planet.PbIII_MPa) or np.isnan(Planet.PbV_MPa):
                        msg = 'Ice V bottom pressure is not greater than ice III bottom pressure. ' + \
                                'This likely indicates TbV_K is too high for the corresponding TbIII_K.' + \
                                f'\nPbIII_MPa = {Planet.PbIII_MPa:.3f}' + \
                                f', TbIII_K = {Planet.Bulk.TbIII_K:.3f}' + \
                                f'\nPbV_MPa = {Planet.PbV_MPa:.3f}' + \
                                f', TbV_K = {Planet.Bulk.TbV_K:.3f}'
                        if Params.ALLOW_BROKEN_MODELS:
                            Planet.PbIII_MPa = np.nan
                            if Params.DO_EXPLOREOGRAM:
                                log.info(msg)
                            else:
                                log.error(msg)
                            Planet.Do.VALID = False
                            Planet.invalidReason = 'TbV_K is too high compared to Tb_K'
                        else:
                            raise ValueError(msg)
                elif Planet.Do.BOTTOM_ICEIII:
                    Planet.PbIII_MPa = GetPfreeze(Planet.Ocean.meltEOS, 3, Planet.Bulk.TbIII_K,
                                                    PLower_MPa=Planet.PbI_MPa, PUpper_MPa=Planet.Ocean.PHydroMax_MPa,
                                                    PRes_MPa=Planet.PfreezeRes_MPa, UNDERPLATE=False, HPNOOCEAN=False,
                                                    ALLOW_BROKEN_MODELS=Params.ALLOW_BROKEN_MODELS,
                                                    DO_EXPLOREOGRAM=Params.DO_EXPLOREOGRAM)
                    if(Planet.PbIII_MPa <= Planet.PbI_MPa) or np.isnan(Planet.PbIII_MPa):
                        msg = 'Ice III bottom pressure is not greater than ice I bottom pressure. ' + \
                                'This likely indicates TbIII_K is too high for the corresponding Tb_K.' + \
                                f'\nPbI_MPa = {Planet.PbI_MPa:.3f}' + \
                                f', Tb_K = {Planet.Bulk.Tb_K:.3f}' + \
                                f'\nPbIII_MPa = {Planet.PbIII_MPa:.3f}' + \
                                f', TbIII_K = {Planet.Bulk.TbIII_K:.3f}'
                        if Params.ALLOW_BROKEN_MODELS:
                            Planet.PbIII_MPa = np.nan
                            if Params.DO_EXPLOREOGRAM:
                                log.info(msg)
                            else:
                                log.error(msg)
                            Planet.Do.VALID = False
                            Planet.invalidReason = 'TbIII_K is too high compared to Tb_K'
                        else:
                            raise ValueError(msg)
                    Planet.Pb_MPa = Planet.PbIII_MPa
                else:
                    Planet.Pb_MPa = Planet.PbI_MPa
            elif Planet.Pb_MPa == 0 and Planet.Bulk.Tsurf_K == Planet.Bulk.Tb_K:
                # This config needs to be caught in SetupInit.
                pass

            else:
                Planet.Pb_MPa = np.nan
                Planet.Do.VALID = False
                Planet.invalidReason = 'Tb_K too high compared to underplate TbIII_K and/or TbV_K'
                if not Params.ALLOW_BROKEN_MODELS:
                    raise RuntimeError('Unable to find a valid pressure corresponding to Bulk.TbX_K values. ' +
                                    f'This is usually because Bulk.Tb_K (currently {Planet.Bulk.Tb_K:.3f}) ' +
                                    'is set too high. Try decreasing Bulk.Tb_K before running again.')

            # Now, we want to check for a convective profile. First, we need to get zb_km, so we need to suppose
            # a whole-layer conductive profile. The densities will change slightly, so we depart from self-consistency
            # here. Repeated applications of IceConvect will get closer to self-consistency.

            if Planet.Pb_MPa > 0 and Planet.Pb_MPa < Planet.P_MPa[0]:
                negDeltaPmsg = f'Calculated Pb value of {Planet.Pb_MPa:.2f} MPa is less than surface pressure of {Planet.P_MPa[0]:.2f} MPa. ' + \
                    'This likely means Tb_K is set too high. Try to decrease and run again to get a valid model.'
                if Params.ALLOW_BROKEN_MODELS:
                    log.warning(negDeltaPmsg + ' ALLOW_BROKEN_MODELS is True, so execution will continue.')
                    Planet.Do.VALID = False
                    Planet.invalidReason = 'Pb_MPa is greater than Psurf_MPa'
                else:
                    raise ValueError(negDeltaPmsg)

            elif Planet.Pb_MPa > 0:
                if Planet.Do.CLATHRATE:
                    """ For ice shells insulated by a layer of clathrate at the surface or against the bottom
                        Calculates state variables of the layer with each pressure step
                        Applies different behavior based on Bulk.clathType:
                            'top': Models a conductive lid of clathrates limited to Bulk.clathMaxThick_m or eLid_m
                                (calculated for convection), whichever is less
                            'bottom': Models a clathrate layer at the ice-ocean interface with a fixed thickness
                                equal to Bulk.clathMaxThick_m. Assumes a purely conductive lid, as justified in
                                Kamata et al. (2019) for Pluto: https://doi.org/10.1038/s41561-019-0369-8
                            'whole': Models clathrates as present throughout the outer ice shell, checking for
                                convection, and assumes no ice I is present in the shell. This option is handled in IceLayers.
                    """
                    if Planet.Do.MIXED_CLATHRATE_ICE:
                        phaseIndex = Constants.phaseClath + 1
                    else:
                        phaseIndex = Constants.phaseClath
                    if Planet.Bulk.clathType == 'top':
                        log.debug('Applying clathrate lid conduction.')

                        Planet.phase[:Planet.Steps.nClath] = phaseIndex
                        if Planet.Do.POROUS_ICE:
                            Planet = IceIConductClathLidPorous(Planet, Params)
                        else:
                            Planet = IceIConductClathLidSolid(Planet, Params)
                    elif Planet.Bulk.clathType == 'bottom':
                        log.debug('Applying clathrate underplating to ice I shell.')
                        Planet.phase[Planet.Steps.nIceI:Planet.Steps.nIbottom] = phaseIndex
                        if Planet.Do.POROUS_ICE:
                            Planet = IceIConductClathUnderplatePorous(Planet, Params)
                        else:
                            Planet = IceIConductClathUnderplateSolid(Planet, Params)

                    elif Planet.Bulk.clathType == 'whole':
                        log.debug('Applying whole-shell clathrate modeling with possible convection.')
                        Planet.phase[:Planet.Steps.nIbottom] = phaseIndex
                        if Planet.Do.POROUS_ICE:
                            Planet = IceIWholeConductPorous(Planet, Params)
                        else:
                            Planet = IceIWholeConductSolid(Planet, Params)
                    else:
                        raise ValueError(f'Bulk.clathType option "{Planet.Bulk.clathType}" is not supported. ' +
                                        'Options are "top", "bottom", and "whole".')
                else:
                    if Planet.Do.POROUS_ICE:
                        Planet = IceIWholeConductPorous(Planet, Params)
                    else:
                        Planet = IceIWholeConductSolid(Planet, Params)

                log.debug('Upper ice initial conductive profile complete.')
                if not Planet.Do.NO_ICE_CONVECTION and not Planet.Bulk.clathType == 'bottom':
                    # Record zb_m to see if it gets adjusted significantly
                    zbOld_m = Planet.z_m[Planet.Steps.nIbottom-1] + 0.0
                    # Now check for convective region and get dimensions if present
                    if Planet.Do.CLATHRATE and Planet.Bulk.clathType == 'whole':
                        if Planet.Do.POROUS_ICE:
                            Planet = ClathShellConvectPorous(Planet, Params)
                        else:
                            Planet = ClathShellConvectSolid(Planet, Params)
                    else:
                        if Planet.Do.POROUS_ICE:
                            Planet = IceIConvectPorous(Planet, Params)
                        else:
                            Planet = IceIConvectSolid(Planet, Params)
                        if Planet.Bulk.clathType == 'top':
                            # Reassign clathrate/ice I transition following convection calcs
                            Planet.zClath_m =  Planet.z_m[Planet.Steps.nClath]
                    # Run IceIConvect a second time if zbI_m changed by more than a set tolerance
                    if(np.abs(Planet.z_m[Planet.Steps.nIbottom-1] - zbOld_m)/Planet.z_m[Planet.Steps.nIbottom-1] > Planet.Bulk.zbChangeTol_frac):
                        log.debug('The bottom depth of surface ice I changed by ' +
                                f'{(Planet.z_m[Planet.Steps.nIbottom-1] - zbOld_m)/1e3:.2f} km from IceIConvect, which is greater than ' +
                                f'{Planet.Bulk.zbChangeTol_frac * 100:.0f}%. running IceIConvect a second time...')
                        if Planet.Do.CLATHRATE and Planet.Bulk.clathType == 'whole':
                            if Planet.Do.POROUS_ICE:
                                Planet = ClathShellConvectPorous(Planet, Params)
                            else:
                                Planet = ClathShellConvectSolid(Planet, Params)
                        else:
                            if Planet.Do.POROUS_ICE:
                                Planet = IceIConvectPorous(Planet, Params)
                            else:
                                Planet = IceIConvectSolid(Planet, Params)
                else:
                    if Planet.Do.NO_ICE_CONVECTION:
                        log.debug('NO_ICE_CONVECTION is True -- skipping ice I convection calculations.')
                    Planet.RaConvect = np.nan
                    Planet.RaCrit = np.nan
                    Planet.Tconv_K = np.nan
                    Planet.etaConv_Pas = np.nan

                    Planet.eLid_m = Planet.z_m[Planet.Steps.nSurfIce]
                    Planet.Dconv_m = 0.0
                    Planet.deltaTBL_m = 0.0
                    # Find the surface heat flux from the conductive profile. This assumes there is no tidal heating!
                    Planet.Ocean.QfromMantle_W = Planet.kTherm_WmK[Planet.Steps.nIbottom-2] * Planet.T_K[Planet.Steps.nIbottom-2] / \
                                                (Planet.z_m[Planet.Steps.nIbottom-1] - Planet.z_m[Planet.Steps.nIbottom-2]) \
                                                * np.log(Planet.T_K[Planet.Steps.nIbottom-1]/Planet.T_K[Planet.Steps.nIbottom-2]) \
                                                * 4*np.pi*(Planet.Bulk.R_m - Planet.z_m[Planet.Steps.nIbottom-1])**2

                # Additional adiabats + possible convection in ice III and/or V underplate layers --
                # for thick, cold ice shells and saline oceans
                if Planet.Do.BOTTOM_ICEV:
                    log.debug('Modeling ice III and V underplating...')
                    Planet = IceIIIUnderplate(Planet, Params)
                    Planet = IceVUnderplate(Planet, Params)
                elif Planet.Do.BOTTOM_ICEIII:
                    log.debug('Modeling ice III underplating...')
                    Planet = IceIIIUnderplate(Planet, Params)

                # Print and save transition pressure and upper ice thickness
                Planet.zb_km = Planet.z_m[Planet.Steps.nSurfIce] / 1e3
                log.info(f'Upper ice transition pressure: {Planet.Pb_MPa:.3f} MPa, ' +
                        f'thickness: {Planet.zb_km:.3f} km.')

                # Set surface HP ice layers to have negative phase ID to differentiate from in-ocean HP ices
                indsHP = np.where(np.logical_and(abs(Planet.phase[:Planet.Steps.nSurfIce]) > 1,
                                                abs(Planet.phase[:Planet.Steps.nSurfIce]) <= 6))[0]
                Planet.phase[:Planet.Steps.nSurfIce][indsHP] = -Planet.phase[:Planet.Steps.nSurfIce][indsHP]

                # Get heat flux out of the possibly convecting region
                Planet.qCon_Wm2 = Planet.Ocean.QfromMantle_W / (4*np.pi * (Planet.Bulk.R_m - Planet.z_m[Planet.Steps.nSurfIce])**2)
                # Get heat flux at the surface, assuming Htidal = Qrad = 0 throughout the entire hydrosphere.
                Planet.qSurf_Wm2 = Planet.Ocean.QfromMantle_W / (4*np.pi * Planet.Bulk.R_m**2)

            elif Planet.Pb_MPa == 0 and Planet.Bulk.Tsurf_K == Planet.Bulk.Tb_K:
                Planet.zb_km = 0
                # This configuration should be accounted for in SetupInit.
            else:
                # Set necessary empty variables for when we have an invalid profile
                Planet.Do.VALID = False
                Planet.invalidReason = 'Pb_MPa is negative'
                Planet.zb_km, Planet.PbClathMax_MPa, Planet.PbIII_MPa, Planet.PbV_MPa, Planet.RaConvect, \
                Planet.RaCrit, Planet.Tconv_K, Planet.TconvIII_K, Planet.TconvV_K, Planet.etaConv_Pas, \
                Planet.etaConvIII_Pas, Planet.etaConvV_Pas, Planet.eLid_m, Planet.Dconv_m, Planet.deltaTBL_m, \
                Planet.qCon_Wm2, Planet.qSurf_Wm2, Planet.TclathTrans_K, Planet.Ocean.QfromMantle_W \
                    = (np.nan for _ in range(19))
                
    return Planet, Params
                
def NonSelfConsistentIceLayer(Planet, Params):
    """ Non-self-consistent ice layer modeling using specified mean properties
        instead of detailed EOS calculations. Uses user-specified layer densities,
        thicknesses, and thermal properties.
        
        Only setup for ice Ih layers at the moment #TODO

        Assigns Planet attributes:
    """
    log.debug('Using non-self-consistent ice layer modeling.')
    if Planet.Do.PARTIAL_DIFFERENTIATION or Planet.Do.NO_DIFFERENTIATION or Planet.Do.CLATHRATE: #TODO: Revisit in future
        raise ValueError('Non-self-consistent ice layer modeling is not supported for partial or no differentiation or clathrate modeling at this moment.')
    # Check that ice shell thickness is set
    if Planet.zb_km is None:
        raise ValueError('Planet.zb_km must be set for non-self-consistent ice modeling.')

    # Set unnecessary variables used in self consistent modeling to false
    """Planet.zb_km, Planet.PbClathMax_MPa, Planet.PbIII_MPa, Planet.PbV_MPa, Planet.RaConvect, \
            Planet.RaCrit, Planet.Tconv_K, Planet.TconvIII_K, Planet.TconvV_K, Planet.etaConv_Pas, \
            Planet.etaConvIII_Pas, Planet.etaConvV_Pas, Planet.eLid_m, Planet.Dconv_m, Planet.deltaTBL_m, \
            Planet.qCon_Wm2, Planet.qSurf_Wm2, Planet.TclathTrans_K, Planet.Ocean.QfromMantle_W \
                = (np.nan for _ in range(19))"""
    
    if Planet.zb_km == 0:
        Planet.Pb_MPa = Planet.Bulk.Psurf_MPa
        Planet.PbI_MPa = Planet.Bulk.Psurf_MPa
        Planet.Bulk.Tb_K = Planet.Bulk.Tsurf_K
        Planet.zb_km = 0.0
        Planet.Steps.nIceI = 0
        Planet.Steps.nIbottom = 0
        Planet.Steps.nSurfIce = 0
        Planet.Do.NO_ICE_CONVECTION = True
        Planet.Do.CLATHRATE = False
        log.debug('No ice shell present (zb_km = 0).')
        return Planet, Params
    else:
        # Set up basic layer structure - each layer is treated as one discrete layer in array
        Planet.Steps.nIbottom = Planet.Steps.nClath + Planet.Steps.nIceI
        Planet.Steps.nIIIbottom = Planet.Steps.nIbottom + Planet.Steps.nIceIIILitho
        Planet.Steps.nSurfIce = Planet.Steps.nIIIbottom + Planet.Steps.nIceVLitho
        # Assign phase values for near-surface ices
        Planet.phase[:Planet.Steps.nIbottom] = 1 # Ice Ih layers (some will be reassigned if Do.CLATHRATE = True)
        Planet.phase[Planet.Steps.nIbottom:Planet.Steps.nIIIbottom] = 3  # Ice III layers
        Planet.phase[Planet.Steps.nIIIbottom:Planet.Steps.nSurfIce] = 5  # Ice V layers
        # Finally, we need to set up the logical array for ice convection layers (this will be default False and overriden if we do ice convection))
        Planet.Steps.iConv = np.zeros(Planet.Steps.nSurfIce, dtype=bool)
        
        # Propogate conduction layers
        Planet = IceIWholeConductSolid(Planet, Params)
        
        # Calculate tempearture and thickness of convective sub-layer using Arrhenius law #TODO ask Flavio
        if not Planet.Do.NO_ICE_CONVECTION:
            Planet = IceIConvectSolid(Planet, Params)
        
        # Get the bottom pressure
        Planet.PbI_MPa = Planet.P_MPa[Planet.Steps.nIbottom]
        Planet.Pb_MPa = Planet.PbI_MPa
        
        
         # Get heat flux out of the possibly convecting region
        Planet.qCon_Wm2 = Planet.Ocean.QfromMantle_W / (4*np.pi * (Planet.Bulk.R_m - Planet.z_m[Planet.Steps.nSurfIce])**2)
        # Get heat flux at the surface, assuming Htidal = Qrad = 0 throughout the entire hydrosphere.
        Planet.qSurf_Wm2 = Planet.Ocean.QfromMantle_W / (4*np.pi * Planet.Bulk.R_m**2)



        log.info(f'Non-self-consistent ice shell thickness: {Planet.zb_km:.3f} km, ' +
                f'transition pressure: {Planet.Pb_MPa:.3f} MPa.')

    return Planet, Params


def IceIIIUnderplate(Planet, Params):
    """ Conductive and convective profile calculations for ice III layers between
        the ocean/underplate ice V and surface ice I layer.

        Assigns Planet attributes:
            PbIII_MPa, all physical layer arrays
    """

    log.debug(f'Ice III bottom phase transition pressure: {Planet.PbIII_MPa:.3f} MPa ' +
             f'at TbIII_K = {Planet.Bulk.TbIII_K:.3f} K.')

    if Planet.Do.POROUS_ICE:
        Planet = IceIIIConductPorous(Planet, Params)
    else:
        Planet = IceIIIConductSolid(Planet, Params)

    if not Planet.Do.NO_ICE_CONVECTION:
        # Record zbIII_m to see if it gets adjusted significantly
        zbIIIold_m = Planet.z_m[Planet.Steps.nIIIbottom-1] + 0.0
        # Now check for convective region and get dimensions if present
        if Planet.Do.POROUS_ICE:
            Planet = IceIIIConvectPorous(Planet, Params)
        else:
            Planet = IceIIIConvectSolid(Planet, Params)
        # Run IceIIIConvect a second time if zbIII_m changed by more than a set tolerance
        if(np.abs(Planet.z_m[Planet.Steps.nIIIbottom-1] - zbIIIold_m)/Planet.z_m[Planet.Steps.nIIIbottom-1] > Planet.Bulk.zbChangeTol_frac):
            log.debug('The bottom depth of underplate ice III changed by ' +
                    f'{(Planet.z_m[Planet.Steps.nIIIbottom-1] - zbIIIold_m)/1e3:.2f} km from IceIIIConvect, which is greater than ' +
                    f'{Planet.Bulk.zbChangeTol_frac * 100:.0f}%. running IceIIIConvect a second time...')
            if Planet.Do.POROUS_ICE:
                Planet = IceIIIConvectPorous(Planet, Params)
            else:
                Planet = IceIIIConvectSolid(Planet, Params)
    else:
        log.debug('NO_ICE_CONVECTION is True -- skipping ice III convection calculations.')
        Planet.eLidIII_m = Planet.z_m[Planet.Steps.nIIIbottom-1]
        Planet.DconvIII_m = 0.0
        Planet.deltaTBLIII_m = 0.0

    return Planet


def IceVUnderplate(Planet, Params):
    """ Conductive and convective profile calculations for ice V layers between
        the ocean and surface ice III layer.

        Assigns Planet attributes:
            PbV_MPa, all physical layer arrays
    """
    log.debug(f'Ice V bottom phase transition pressure: {Planet.PbV_MPa:.3f} MPa ' +
                             f'at TbV_K = {Planet.Bulk.TbV_K:.3f} K.')

    if Planet.Do.POROUS_ICE:
        Planet = IceVConductPorous(Planet, Params)
    else:
        Planet = IceVConductSolid(Planet, Params)

    if not Planet.Do.NO_ICE_CONVECTION:
        # Record zbV_m to see if it gets adjusted significantly
        zbVold_m = Planet.z_m[Planet.Steps.nSurfIce-1] + 0.0
        # Now check for convective region and get dimensions if present
        if Planet.Do.POROUS_ICE:
            Planet = IceVConvectPorous(Planet, Params)
        else:
            Planet = IceVConvectSolid(Planet, Params)
        # Run IceVConvect a second time if zbV_m changed by more than a set tolerance
        if(np.abs(Planet.z_m[Planet.Steps.nSurfIce-1] - zbVold_m)/Planet.z_m[Planet.Steps.nSurfIce-1] > Planet.Bulk.zbChangeTol_frac):
            log.debug('The bottom depth of underplate ice V changed by ' +
                     f'{(Planet.z_m[Planet.Steps.nSurfIce-1] - zbVold_m)/1e3:.2f} km from IceVConvect, which is greater than ' +
                     f'{Planet.Bulk.zbChangeTol_frac * 100:.0f}%. running IceVConvect a second time...')
            if Planet.Do.POROUS_ICE:
                Planet = IceVConvectPorous(Planet, Params)
            else:
                Planet = IceVConvectSolid(Planet, Params)
    else:
        log.debug('NO_ICE_CONVECTION is True -- skipping ice V convection calculations.')
        Planet.eLidV_m = Planet.z_m[Planet.Steps.nSurfIce-1]
        Planet.DconvV_m = 0.0
        Planet.deltaTBLV_m = 0.0

    return Planet



def IceVIUnderplate(Planet, Params):
    """ Conductive and convective profile calculations for ice V layers between
        the ocean and surface ice III layer.

        Assigns Planet attributes:
            PbV_MPa, all physical layer arrays
    """
    log.debug(f'Ice V bottom phase transition pressure: {Planet.PbV_MPa:.3f} MPa ' +
                             f'at TbV_K = {Planet.Bulk.TbV_K:.3f} K.')

    if Planet.Do.POROUS_ICE:
        Planet = IceVIConductPorous(Planet, Params)
    else:
        Planet = IceVIConductSolid(Planet, Params)

    if not Planet.Do.NO_ICE_CONVECTION:
        # Record zbV_m to see if it gets adjusted significantly
        zbVold_m = Planet.z_m[Planet.Steps.nSurfIce-1] + 0.0
        # Now check for convective region and get dimensions if present
        if Planet.Do.POROUS_ICE:
            Planet = IceVIConvectPorous(Planet, Params)
        else:
            Planet = IceVIConvectSolid(Planet, Params)
        # Run IceVIConvect a second time if zbV_m changed by more than a set tolerance
        if(np.abs(Planet.z_m[Planet.Steps.nSurfIce-1] - zbVold_m)/Planet.z_m[Planet.Steps.nSurfIce-1] > Planet.Bulk.zbChangeTol_frac):
            log.debug('The bottom depth of underplate ice V changed by ' +
                     f'{(Planet.z_m[Planet.Steps.nSurfIce-1] - zbVold_m)/1e3:.2f} km from IceVIConvect, which is greater than ' +
                     f'{Planet.Bulk.zbChangeTol_frac * 100:.0f}%. running IceVIConvect a second time...')
            if Planet.Do.POROUS_ICE:
                Planet = IceVIConvectPorous(Planet, Params)
            else:
                Planet = IceVIConvectSolid(Planet, Params)
    else:
        log.debug('NO_ICE_CONVECTION is True -- skipping ice V convection calculations.')
        Planet.eLidV_m = Planet.z_m[Planet.Steps.nSurfIce-1]
        Planet.DconvV_m = 0.0
        Planet.deltaTBLV_m = 0.0

    return Planet


def iceShellTFreeze(T, Planet, Params, zb_approximate_km):
    """Computes the difference between target and actual ice shell thickness for a given temperature.
    
    Args:
        T: Temperature at ice-ocean boundary (K)
        Planet: Planet object
        Params: Parameters object
        zb_approximate_km: Target ice shell thickness (km)
        
    Returns:
        Difference between target and computed ice shell thickness (km)
        
    Raises:
        ValueError: If ice computation fails (PbI_MPa == 0.0)
    """
    Planet_copy = deepcopy(Planet)
    Params_copy = deepcopy(Params)
    Planet_copy.Do.ICEIh_THICKNESS  = False
    Planet_copy.Bulk.Tb_K = T
    result = IceLayers(Planet_copy, Params_copy)
    
    if result.PbI_MPa == 0.0:
        raise ValueError('Ice layer computation failed: PbI_MPa == 0.0')
    
    return zb_approximate_km - result.zb_km


def GetIceShellTFreeze(Planet, Params):
    """
    Determines the temperature at which ice transitions to water for a given ice shell thickness.
    
    This function operates in two main stages:
    1. It establishes a valid temperature bracket [T_lower, T_upper] where the
       IceLayers computation can succeed. It starts by calculating freezing
       temperatures at the boundaries of the pressure search range (PfreezeLower_MPa
       and PfreezeUpper_MPa).
    2. It uses a root-finding algorithm on the `iceShellTFreeze` helper function
       to find the precise temperature (Tb_K) that results in the desired ice
       shell thickness (Planet.Bulk.zb_approximate_km). It iteratively adjusts
       the bracket if the initial range does not contain the solution.

    Args:
        Planet: Planet object containing ice shell properties.
        Params: Parameters object.
        
    Returns:
        float: Computed freezing temperature (K).
        
    Raises:
        ValueError: If a valid temperature bracket for the root-finding cannot be established.
    """
    zb_approximate_km = Planet.Bulk.zb_approximate_km
    
    # Deepcopy Params once to avoid repeatedly copying it inside the solver loop.
    Params_copy = deepcopy(Params)
    # Create wrapper class that we can use for GetZero function that allows us to use zb_approximate_km tolerance
    class IceShellResidual:
        def __init__(self, Planet, Params, zb_target_km, zb_tol_km=1.0):
            self.Planet = Planet
            self.zb_target_km = zb_target_km
            self.zb_tol_km = zb_tol_km
            self.best_planet = None
            self.best_T = None
            self.last_residual = None

        def __call__(self, T):
            Planet_copy = deepcopy(self.Planet)
            Planet_copy.Do.ICEIh_THICKNESS = False
            Planet_copy.Bulk.Tb_K = T

            result = IceLayers(Planet_copy, Params_copy)

            if result.PbI_MPa == 0.0:
                raise ValueError("Ice layer computation failed: PbI_MPa == 0.0")

            residual = self.zb_target_km - result.zb_km
            abs_residual = abs(residual)

            # Save the best result
            self.last_residual = abs_residual
            self.best_planet = Planet_copy
            self.best_planet.Do.ICEIh_THICKNESS = True
            self.best_T = T

            # Early exit condition
            if abs_residual < self.zb_tol_km:
                raise StopIteration("Residual within tolerance")

            return residual

    # Part 1: Establish a valid temperature bracket
    try:
        # Find the lower and upper temperature limits given the loewr and upper P freeze we specify
        # Note we have to subtract the TfreezeRes_K from the upper limit to ensure we have a valid upper Tb_K (i.e. within the phase diagram we generate)
        TupperLimit_K = GetTfreeze(Planet.Ocean.meltEOS, Planet.PfreezeLower_MPa, Planet.TfreezeLower_K, TRes_K=-Planet.TfreezeRes_K)
        TlowerLimit_K = GetTfreeze(Planet.Ocean.meltEOS, Planet.PfreezeUpper_MPa, Planet.TfreezeLower_K, TRes_K=Planet.TfreezeRes_K)
        log.debug(f"Established temperature bounds from phase diagram: [{TlowerLimit_K:.2f}, {TupperLimit_K:.2f}] K")
    except Exception as e:
        raise ValueError(f"Could not determine temperature bounds from phase diagram. Try lowering Planet.TfreezeLower_K, which represents the lower limit of the temperature search.")
    solver = IceShellResidual(Planet, Params, zb_approximate_km, zb_tol_km=0.01)
    # Find the precise freezing temperature using root finding
    try:
        sol = GetZero(solver, 
                      bracket=[TlowerLimit_K, TupperLimit_K], 
                      xtol=Planet.TfreezeRes_K)
        
        T_freeze = sol.root
        log.debug(f"Computed ice shell freezing temperature: {T_freeze:.3f} K after {sol.function_calls} iterations.")
        return solver.best_planet

    except StopIteration:
        log.debug("Stopped early because zb_km is within tolerance.")
        return solver.best_planet
    except Exception as e:
        msg = ''
        if Planet.TfreezeLower_K < TlowerLimit_K + Planet.TfreezeRes_K:
            msg += 'Try increasing Planet.PfreezeUpper_MPa'
        else:
            msg += 'Try decreasing Planet.TfreezeLower_K'
        raise ValueError(f"Root finding for corresponding Tb_K to achieve a given zb_approximate_km failed for temperature bracket [{TlowerLimit_K:.2f}, {TupperLimit_K:.2f}] K. {msg}. Also, the code is directly outptting the error: {e}")



def OceanLayers(Planet, Params):
    """ Wrapper function for ocean layer propogation. Decides between self consistent and non-self consistent modeling.
    """
    Timing.setStartingTime(time.time())
    # Early branching for non-self-consistent modeling
    if Planet.Do.NON_SELF_CONSISTENT:
        Planet, Params = NonSelfConsistentOceanLayer(Planet, Params)
    else:
        Planet, Params = SelfConsistentOceanLayer(Planet, Params)
    Timing.printFunctionTimeDifference('OceanLayers()', time.time())
    return Planet


def NonSelfConsistentOceanLayer(Planet, Params):
    """ Non-self-consistent ocean layer modeling using specified mean properties
        instead of detailed EOS calculations. Uses user-specified layer densities,
        thicknesses, and thermal properties.

        Assigns Planet attributes:
            phase, r_m, z_m, g_ms2, T_K, P_MPa, rho_kgm3, Cp_JkgK, alpha_pK, MLayer_kg
    """
    if Planet.Do.VALID and not Planet.Do.NO_OCEAN:
        log.debug('Evaluating non-self-consistent ocean layer.')
        
        
        # Case is handled in SetupInit.py
        if Planet.D_km == 0:
            Planet.Steps.nHydro = Planet.Steps.nSurfIce
            return Planet, Params
        
        else:
            # Get number of ocean layers
            nOcean = Planet.Steps.nHydro - Planet.Steps.nSurfIce
            # Set up depth profile to bottom of hydrosphere
            zOcean_m = np.linspace(Planet.zb_km * 1e3, (Planet.D_km + Planet.zb_km) * 1e3, nOcean + 1)
            Planet.z_m[Planet.Steps.nSurfIce:Planet.Steps.nHydro + 1] = zOcean_m
            # Get the oceanEOS
            POcean_MPa = np.arange(Planet.Pb_MPa, Planet.Ocean.PHydroMax_MPa, Planet.Ocean.deltaP)
            TOcean_K = np.arange(Planet.Bulk.Tb_K, Planet.Ocean.THydroMax_K, Planet.Ocean.deltaT)
            Planet.Ocean.EOS = GetOceanEOS(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, POcean_MPa, TOcean_K,
                                        Planet.Ocean.MgSO4elecType, rhoType=Planet.Ocean.MgSO4rhoType,
                                        scalingType=Planet.Ocean.MgSO4scalingType, FORCE_NEW=Params.FORCE_EOS_RECALC,
                                        phaseType=Planet.Ocean.phaseType, EXTRAP=Params.EXTRAP_OCEAN,
                                        sigmaFixed_Sm=Planet.Ocean.sigmaFixed_Sm, LOOKUP_HIRES=Planet.Do.OCEAN_PHASE_HIRES, kThermConst_WmK=Planet.Ocean.kThermWater_WmK,
                                        doConstantProps=Planet.Do.CONSTANTPROPSEOS,
                                        constantProperties=Planet.Ocean.oceanConstantProperties,
                                        propsStepReductionFactor=Planet.Ocean.propsStepReductionFactor)
            # Log and propogate the first ocean layer separately - to prevent index issues
            log.debug(f'il: {Planet.Steps.nSurfIce:d}; P_MPa: {Planet.P_MPa[Planet.Steps.nSurfIce]:.3f}; T_K: {Planet.T_K[Planet.Steps.nSurfIce]:.3f}; phase: {Planet.phase[Planet.Steps.nSurfIce]:d}')
            Planet.rhoMatrix_kgm3[Planet.Steps.nSurfIce] = Planet.Ocean.EOS.fn_rho_kgm3(Planet.P_MPa[Planet.Steps.nSurfIce], Planet.T_K[Planet.Steps.nSurfIce])
            Planet.Cp_JkgK[Planet.Steps.nSurfIce] = Planet.Ocean.EOS.fn_Cp_JkgK(Planet.P_MPa[Planet.Steps.nSurfIce], Planet.T_K[Planet.Steps.nSurfIce])
            Planet.alpha_pK[Planet.Steps.nSurfIce] = Planet.Ocean.EOS.fn_alpha_pK(Planet.P_MPa[Planet.Steps.nSurfIce], Planet.T_K[Planet.Steps.nSurfIce])
            Planet.kTherm_WmK[Planet.Steps.nSurfIce] = Planet.Ocean.EOS.fn_kTherm_WmK(Planet.P_MPa[Planet.Steps.nSurfIce], Planet.T_K[Planet.Steps.nSurfIce])
            
            # Propagate layer-top arrays
            thisMAbove_kg = np.sum(Planet.MLayer_kg[:Planet.Steps.nSurfIce-1])
            # Get constant gravity if we will be assigning it
            if Planet.Do.CONSTANT_GRAVITY:
                Planet.g_ms2[Planet.Steps.nSurfIce+1:] = Constants.G * (Planet.Bulk.M_kg - thisMAbove_kg) / Planet.r_m[Planet.Steps.nSurfIce-1]**2
            else:
                # Ensure g values to be assigned are zero since we will be adding to them
                Planet.g_ms2[Planet.Steps.nSurfIce+1:] = 0
            # Assign 0 or 1 multiplier for constant/variable gravity calcs in loop
            VAR_GRAV = int(not Planet.Do.CONSTANT_GRAVITY)

            for i in range(Planet.Steps.nSurfIce+1, Planet.Steps.nHydro+1):
                # Increment depth based on change in pressure, combined with gravity and density
                Planet.P_MPa[i] = Planet.P_MPa[i-1] + Planet.rhoMatrix_kgm3[i-1] * Planet.g_ms2[i-1] * (Planet.z_m[i] - Planet.z_m[i-1]) * 1e-6
                # Convert depth to radius
                Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
                # Calculate layer mass
                Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rhoMatrix_kgm3[i-1] * (Planet.r_m[i-1] ** 3 - Planet.r_m[i] ** 3)
                thisMAbove_kg += Planet.MLayer_kg[i-1]
                thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
                # Use remaining mass below in Gauss's law for gravity to get g at the top of this layer
                Planet.g_ms2[i] += VAR_GRAV * Constants.G * thisMBelow_kg / Planet.r_m[i] ** 2
                # Now use the present layer's properties to calculate an adiabatic thermal profile for layers below
                Planet.T_K[i] = Planet.T_K[i-1] + Planet.T_K[i-1] * Planet.alpha_pK[i-1] / \
                                Planet.Cp_JkgK[i-1] / Planet.rhoMatrix_kgm3[i-1] * ((Planet.P_MPa[i] - Planet.P_MPa[i-1]) * 1e6)
                # Calculate this layer's properties
                if i < Planet.Steps.nHydro:
                    Planet.rhoMatrix_kgm3[i] = Planet.Ocean.EOS.fn_rho_kgm3(Planet.P_MPa[i], Planet.T_K[i])
                    Planet.Cp_JkgK[i] = Planet.Ocean.EOS.fn_Cp_JkgK(Planet.P_MPa[i], Planet.T_K[i])
                    Planet.alpha_pK[i] = Planet.Ocean.EOS.fn_alpha_pK(Planet.P_MPa[i], Planet.T_K[i])
                    Planet.kTherm_WmK[i] = Planet.Ocean.EOS.fn_kTherm_WmK(Planet.P_MPa[i], Planet.T_K[i])
                if i == Planet.Steps.nHydro:
                    log.debug(f'Propogating starting point for next layer...')
                    log.debug(f'il: {i:d}; P_MPa: {Planet.P_MPa[i]:.3f}; T_K: {Planet.T_K[i]:.3f}')
                else:
                    log.debug(f'il: {i:d}; P_MPa: {Planet.P_MPa[i]:.3f}; T_K: {Planet.T_K[i]:.3f}; phase: {Planet.phase[i]:d}')
            
            # Fill in properties
            Planet.rho_kgm3[Planet.Steps.nSurfIce:Planet.Steps.nHydro+1] = Planet.rhoMatrix_kgm3[Planet.Steps.nSurfIce:Planet.Steps.nHydro+1]
           



            log.info(f'Ocean layers complete. zMax: {Planet.z_m[Planet.Steps.nHydro]/1e3:.1f} km, ' +
                 f'upper ice thickness zb: {Planet.zb_km:.3f} km')

    
    return Planet, Params


def SelfConsistentOceanLayer(Planet, Params):
    """ Geophysical and thermodynamic calculations for ocean layer
        Calculates state variables of the layer with each pressure step

        Assigns Planet attributes:
            phase, r_m, z_m, g_ms2, T_K, P_MPa, rho_kgm3, Cp_JkgK, alpha_pK, MLayer_kg
    """
    if Planet.Do.VALID and (not Planet.Do.NO_OCEAN or Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES):
        log.debug('Evaluating ocean layers.')
        if Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES:
            log.debug(f'Planet.Do.NO_OCEAN_EXCEPT_INNER_LAYERS is set to True. In this case, we will propogate the ocean layers function from the calculated bottom pressure of {Planet.Pb_MPa} MPa' +
                      f' and input bottom temperature of {Planet.Bulk.Tb_K} to calculate high pressure ices. These input conditions should be such that no liquid layers will form, else we will raise an error.')

        # Confirm that we haven't made mistakes in phase assignment in IceLayers()
        # Note that this assignment is only temporary, as we re-check the phase of
        # this layer with the ocean EOS momentarily--it could also be ice II or III.
        if Planet.phase[Planet.Steps.nSurfIce] != 0:
            raise ValueError('Phase of first "ocean" layer is not zero.')

        # Start ocean pressure at the melting pressure and assign linear pressure profile for layers
        # First, check that PHydroMax is high enough for the ice shell solution
        if Planet.Ocean.PHydroMax_MPa < Planet.Pb_MPa:
            raise ValueError(f'Ocean.PHydroMax_MPa is {Planet.Ocean.PHydroMax_MPa:.1f}, but the ' +
                             f'ice shell bottom pressure was found to be {Planet.Pb_MPa} MPa. ' +
                             f'Increase Ocean.PHydroMax_MPa to well beyond this value to avoid ' +
                             f'this problem.')
        POcean_MPa = np.arange(Planet.Pb_MPa, Planet.Ocean.PHydroMax_MPa, Planet.Ocean.deltaP)
        Planet.Steps.nOceanMax = np.size(POcean_MPa)

        # Initialize remaining local arrays
        TOcean_K, rhoOcean_kgm3, CpOcean_JkgK, alphaOcean_pK, kThermOcean_WmK \
            = (np.zeros(Planet.Steps.nOceanMax) for _ in range(5))
        TOcean_K = np.insert(TOcean_K, 0, Planet.T_K[Planet.Steps.nSurfIce])

        # Add HP ices to iceEOS if needed
        PHydroMax_MPa = Planet.Ocean.PHydroMax_MPa
        if Planet.Do.POROUS_ROCK and Planet.Sil.PHydroMax_MPa > PHydroMax_MPa:
            PHydroMax_MPa = Planet.Sil.PHydroMax_MPa
            PHPices_MPa = np.arange(POcean_MPa[0], PHydroMax_MPa, Planet.Ocean.deltaP)
        else:
            PHPices_MPa = POcean_MPa
        if PHydroMax_MPa > Constants.PminHPices_MPa:
            GetOceanHPIceEOS(Planet, Params, PHPices_MPa, minPres_MPa=Params.minPres_MPa, minTres_K=Params.minTres_K)
        
        # If we are not allowing liquid layers in the ocean, we need to make sure the first layer is a high pressure ice
        if Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES:
            initialPOcean_MPa = POcean_MPa[0]
            thisPhase = Planet.Ocean.EOS.fn_phase(initialPOcean_MPa, TOcean_K[0]).astype(np.int_)
            # Due to numerical approximation with GetZero function, thisPhase can be an ice Ih or ocean layer
            if thisPhase < 2:
                # Depending on the root answer that fn_phase finds, thisPhase can be an ocean layer. So we must get the ice phase by incrementing P by deltaP so we are on the high pressure side of the phase diagram
                thisPhase = Planet.Ocean.EOS.fn_phase(initialPOcean_MPa+Planet.Ocean.EOS.deltaP, TOcean_K[0]).astype(np.int_)
                # Past error debugging that is too complicated. We have simplified
                """log.warning(f'The first calculated phase is not a high pressure ice layer. \n' +
                                 f'When Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES is True, the layers below the initial ice propogation should be high pressure ices.\n' +
                                 f' T will be set to be lower than the melting temp temporarily and P will be set slightly higher than the melting pressure to construct the first high pressure ice layer.')
                # Increase P by deltaP temporarily so we can move into the next 'layer' of phase diagram
                initialPOcean_MPa += Planet.Ocean.deltaP
                # Get the freezing temperature and decrease by TfreezeOffset_K to ensure we stay within the high pressure phase diagram
                TOcean_K[0] = GetTfreeze(Planet.Ocean.EOS, initialPOcean_MPa, TOcean_K[0], TRes_K=Planet.Ocean.TfreezeOffset_K) - Planet.Ocean.TfreezeOffset_K
                thisPhase = Planet.Ocean.EOS.fn_phase(initialPOcean_MPa, TOcean_K[0]).astype(np.int_)
                if thisPhase == 0:
                    raise ValueError('Even after slightly increasing P and slightly decreasing T, the first phase is still a liquid. \n' +
                                     f'Generating an ocean is not desired when Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES is True. \n' +
                                     f'Please check your Planet.Bulk.Tb_K is set correctly for high pressure ices to form. Namely, try decreasing it further.')"""
            Planet.phase[Planet.Steps.nSurfIce] = thisPhase
            log.debug(f'il: {Planet.Steps.nSurfIce:d}; P_MPa: {POcean_MPa[0]:.3f}; ' +
            f'T_K: {TOcean_K[0]:.3f}; phase: {Planet.phase[Planet.Steps.nSurfIce]:d}')
            thisPhaseName = PhaseConv(thisPhase)
            rhoOcean_kgm3[0] = Planet.Ocean.iceEOS[thisPhaseName].fn_rho_kgm3(POcean_MPa[0], TOcean_K[0])
            CpOcean_JkgK[0] = Planet.Ocean.iceEOS[thisPhaseName].fn_Cp_JkgK(POcean_MPa[0], TOcean_K[0])
            alphaOcean_pK[0] = Planet.Ocean.iceEOS[thisPhaseName].fn_alpha_pK(POcean_MPa[0], TOcean_K[0])
            kThermOcean_WmK[0] = Planet.Ocean.iceEOS[thisPhaseName].fn_kTherm_WmK(POcean_MPa[0], TOcean_K[0])
            # We use GetTfreeze here to propagate the next high pressure ice layer temperature
            TOcean_K[1] = np.mean([GetTfreeze(Planet.Ocean.EOS, initialPOcean_MPa, TOcean_K[0], TRes_K=Planet.Ocean.TfreezeOffset_K) - Planet.Ocean.TfreezeOffset_K, TOcean_K[0]])
            iStart = 1
        else:
            # Do initial ocean step separately in order to catch potential Melosh layer--
            # see Melosh et al. (2004): https://doi.org/10.1016/j.icarus.2003.11.026
            # insert no_ocean first layer as ice Ih if Planet.Do.NO_OCEAN is True
            rhoOcean_kgm3[0] = Planet.Ocean.EOS.fn_rho_kgm3(POcean_MPa[0], TOcean_K[0])
            CpOcean_JkgK[0] = Planet.Ocean.EOS.fn_Cp_JkgK(POcean_MPa[0], TOcean_K[0])
            alphaOcean_pK[0] = Planet.Ocean.EOS.fn_alpha_pK(POcean_MPa[0], TOcean_K[0])
            kThermOcean_WmK[0] = Planet.Ocean.EOS.fn_kTherm_WmK(POcean_MPa[0], TOcean_K[0])

            log.debug(f'il: {Planet.Steps.nSurfIce:d}; P_MPa: {POcean_MPa[0]:.3f}; ' +
                    f'T_K: {TOcean_K[0]:.3f}; phase: {Planet.phase[Planet.Steps.nSurfIce]:d}')

            if Planet.phase[Planet.Steps.nSurfIce] == 0 and alphaOcean_pK[0] < 0 and not Planet.Do.NO_MELOSH_LAYER:
                log.info(f'Thermal expansivity alpha at the ice-ocean interface is negative. Modeling Melosh et al. conductive layer.')
                # Layer should be thin, so we just use a fixed dT/dz value
                dTdz = Planet.Ocean.QfromMantle_W / (4*np.pi * Planet.r_m[Planet.Steps.nSurfIce]**2) / kThermOcean_WmK[0]
                i = 0
                # Use a smaller pressure step to make sure we don't overshoot by a lot
                deltaPMelosh = Planet.Ocean.deltaP / 100
                # Initialize for while loop
                alphaMelosh_pK = alphaOcean_pK[0] + 0
                rhoMelosh_kgm3 = rhoOcean_kgm3[0] + 0
                gMelosh_ms2 = Planet.g_ms2[Planet.Steps.nSurfIce] + 0
                TMelosh_K = TOcean_K[0] + 0
                deltaPtop = 0
                zMelosh = 0
                while alphaMelosh_pK < 0 and Planet.phase[Planet.Steps.nSurfIce+i] == 0:
                    dz = deltaPMelosh*1e6 / gMelosh_ms2 / rhoMelosh_kgm3
                    zMelosh += dz
                    deltaPtop += deltaPMelosh
                    thisP_MPa = deltaPtop + POcean_MPa[0]
                    # Model temperature as linear and conductive in this layer
                    TMelosh_K += dTdz * dz

                    rhoMelosh_kgm3 = Planet.Ocean.EOS.fn_rho_kgm3(thisP_MPa, TMelosh_K)
                    alphaMelosh_pK = Planet.Ocean.EOS.fn_alpha_pK(thisP_MPa, TMelosh_K)

                    if alphaMelosh_pK > 0 or deltaPtop >= (i+1)*Planet.Ocean.deltaP:
                        i += 1
                        POcean_MPa[i] = thisP_MPa + 0
                        TOcean_K[i] = TMelosh_K + 0
                        rhoOcean_kgm3[i] = Planet.Ocean.EOS.fn_rho_kgm3(POcean_MPa[i], TOcean_K[i])
                        CpOcean_JkgK[i] = Planet.Ocean.EOS.fn_Cp_JkgK(POcean_MPa[i], TOcean_K[i])
                        alphaOcean_pK[i] = Planet.Ocean.EOS.fn_alpha_pK(POcean_MPa[i], TOcean_K[i])
                        kThermOcean_WmK[i] = Planet.Ocean.EOS.fn_kTherm_WmK(POcean_MPa[i], TOcean_K[i])
                        Planet.phase[Planet.Steps.nSurfIce+i] = Planet.Ocean.EOS.fn_phase(POcean_MPa[i], TOcean_K[i]).astype(np.int_)
                        log.debug(f'il: {Planet.Steps.nSurfIce+i:d}; P_MPa: {POcean_MPa[i]:.3f}; ' +
                                f'T_K: {TOcean_K[i]:.3f}; phase: {Planet.phase[Planet.Steps.nSurfIce+i]:d}')
                iStart = i
                # Reset pressure profile to use standard pressure step below Melosh layer bottom
                POcean_MPa[i+1:] = np.linspace(POcean_MPa[i], POcean_MPa[-1], Planet.Steps.nOceanMax - i)[1:]
                log.info(f'Melosh et al. layer complete, thickness {zMelosh:.1f} m.')
            else:
                if Planet.Do.NO_MELOSH_LAYER and alphaOcean_pK[0] < 0:
                    log.debug('Melosh layer is present, but Do.NO_MELOSH_LAYER is True. alpha_pK will be set to zero here.')
                    alphaOcean_pK[0] = 0
                # Now use the present layer's properties to calculate an adiabatic thermal profile for layers below
                TOcean_K[1] = TOcean_K[0] + alphaOcean_pK[0] * TOcean_K[0] / \
                                CpOcean_JkgK[0] / rhoOcean_kgm3[0] * Planet.Ocean.deltaP*1e6
                iStart = 1
        for i in range(iStart, Planet.Steps.nOceanMax):
            Planet.phase[Planet.Steps.nSurfIce+i] = Planet.Ocean.EOS.fn_phase(POcean_MPa[i], TOcean_K[i]).astype(np.int_)
            if not Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES and i < 4 and Planet.phase[Planet.Steps.nSurfIce+i] != 0:
                log.debug(f'Top ocean layers (i={i}) are not liquid. This will cause indexing problems. ' +
                          'T will be set to exceed the melting temp temporarily to construct at least 4 ocean layers.')
                Planet.THIN_OCEAN = True
                TOcean_K[i] = GetTfreeze(Planet.Ocean.EOS, POcean_MPa[i], TOcean_K[i], TRes_K = 0.00001)
                Planet.phase[Planet.Steps.nSurfIce+i] = 0
            log.debug(f'il: {Planet.Steps.nSurfIce+i:d}; P_MPa: {POcean_MPa[i]:.3f}; ' +
                      f'T_K: {TOcean_K[i]:.3f}; phase: {Planet.phase[Planet.Steps.nSurfIce+i]:d}')
            if Planet.phase[Planet.Steps.nSurfIce+i] == 0 and Planet.phase[Planet.Steps.nSurfIce+i-1] > 2:
                log.warning(f'Ocean layer {i} is a liquid layer, but the previous layer is a high pressure ice layer. ')
            if Planet.phase[Planet.Steps.nSurfIce+i] < 2:
                # Liquid water layers -- get fluid properties for the present layer but with the
                # overlaying layer's temperature. Note that we include ice Ih in these layers because
                # ice Ih layers result only from instabilities in phase diagram calculations. There should
                # not be any ice Ih below the ice--ocean interface at Tb.
                rhoOcean_kgm3[i] = Planet.Ocean.EOS.fn_rho_kgm3(POcean_MPa[i], TOcean_K[i])
                CpOcean_JkgK[i] = Planet.Ocean.EOS.fn_Cp_JkgK(POcean_MPa[i], TOcean_K[i])
                alphaOcean_pK[i] = Planet.Ocean.EOS.fn_alpha_pK(POcean_MPa[i], TOcean_K[i])
                kThermOcean_WmK[i] = Planet.Ocean.EOS.fn_kTherm_WmK(POcean_MPa[i], TOcean_K[i])
                if Planet.Do.NO_MELOSH_LAYER:
                    alphaOcean_pK[i] = np.abs(alphaOcean_pK[i])
                # Now use the present layer's properties to calculate an adiabatic thermal profile for layers below
                TOcean_K[i+1] = TOcean_K[i] + alphaOcean_pK[i] * TOcean_K[i] / \
                                CpOcean_JkgK[i] / rhoOcean_kgm3[i] * Planet.Ocean.deltaP*1e6
            else:
                # Undersea high-pressure ices -- we use GetTfreeze here to propagate the layer temperatures.
                # This is based on an assumption that the undersea HP ices are vigorously mixed by
                # two-phase convection, such that each layer is in local equilibrium with the liquid,
                # meaning each layer's temperature is equal to the melting temperature.
                # We implement this by averaging the upper layer temp with the melting temp minus two times the deltaT of the EOS, which ensures we are on the solid side of the phase diagram
                # to step more gently and avoid overshooting that causes phase oscillations.
                thisPhase = PhaseConv(Planet.phase[Planet.Steps.nSurfIce+i])
                rhoOcean_kgm3[i] = Planet.Ocean.iceEOS[thisPhase].fn_rho_kgm3(POcean_MPa[i], TOcean_K[i])
                CpOcean_JkgK[i] = Planet.Ocean.iceEOS[thisPhase].fn_Cp_JkgK(POcean_MPa[i], TOcean_K[i])
                alphaOcean_pK[i] = Planet.Ocean.iceEOS[thisPhase].fn_alpha_pK(POcean_MPa[i], TOcean_K[i])
                kThermOcean_WmK[i] = Planet.Ocean.iceEOS[thisPhase].fn_kTherm_WmK(POcean_MPa[i], TOcean_K[i])
                TOcean_K[i+1] = np.mean([GetTfreeze(Planet.Ocean.EOS, POcean_MPa[i], TOcean_K[i], TRes_K=Planet.Ocean.EOS.deltaT) - Planet.Ocean.EOS.deltaT*2, TOcean_K[i]])

        # Assign ocean layer critical properties to Planet fields
        Planet.P_MPa[Planet.Steps.nSurfIce:Planet.Steps.nSurfIce + Planet.Steps.nOceanMax] = POcean_MPa
        Planet.T_K[Planet.Steps.nSurfIce:Planet.Steps.nSurfIce + Planet.Steps.nOceanMax] = TOcean_K[:-1]
        Planet.rho_kgm3[Planet.Steps.nSurfIce:Planet.Steps.nSurfIce + Planet.Steps.nOceanMax] = rhoOcean_kgm3
        Planet.Cp_JkgK[Planet.Steps.nSurfIce:Planet.Steps.nSurfIce + Planet.Steps.nOceanMax] = CpOcean_JkgK
        Planet.alpha_pK[Planet.Steps.nSurfIce:Planet.Steps.nSurfIce + Planet.Steps.nOceanMax] = alphaOcean_pK
        Planet.kTherm_WmK[Planet.Steps.nSurfIce:Planet.Steps.nSurfIce + Planet.Steps.nOceanMax] = kThermOcean_WmK

        # Apply smoothing filter to avoid bumpiness from discretized phase diagram
        HPphases = np.logical_and(Planet.phase[:Planet.Steps.nHydroMax] > 1,
                                  Planet.phase[:Planet.Steps.nHydroMax] < 10)
        if Planet.Ocean.phaseType == 'lookup' and Planet.Do.HP_MELT_SMOOTHING and np.any(HPphases):
            if Planet.Do.FIXED_HPSMOOTH_WINDOW:
                window = Planet.Ocean.smoothingWindowOverride
            else:
                # Get the ratio of step size used in hydrosphere relative to lookup table spacing
                Pratio = Planet.Ocean.EOS.EOSdeltaP / Planet.Ocean.deltaP
                window = int(Planet.Ocean.smoothingFactor * Pratio)

            # Ensure smoothing will happen
            window = np.maximum(window, Planet.Ocean.smoothingPolyOrder + 1)

            # Ensure window is odd
            window = window + np.mod(window + 1, 2)
            log.debug(f'Applying smoothing to in-ocean HP ice profile with a window size of {window} pts ({Planet.Ocean.deltaP*window:.1f} MPa).')
            Planet.T_K[HPphases] = savgol_filter(Planet.T_K[HPphases], window, Planet.Ocean.smoothingPolyOrder)
            Planet.rho_kgm3[HPphases] = savgol_filter(Planet.rho_kgm3[HPphases], window, Planet.Ocean.smoothingPolyOrder)
            Planet.Cp_JkgK[HPphases] = savgol_filter(Planet.Cp_JkgK[HPphases], window, Planet.Ocean.smoothingPolyOrder)
            Planet.alpha_pK[HPphases] = savgol_filter(Planet.alpha_pK[HPphases], window, Planet.Ocean.smoothingPolyOrder)
            Planet.kTherm_WmK[HPphases] = savgol_filter(Planet.kTherm_WmK[HPphases], window, Planet.Ocean.smoothingPolyOrder)

        # Evaluate remaining physical quantities for ocean layers
        MAbove_kg = np.sum(Planet.MLayer_kg[:Planet.Steps.nSurfIce])
        # Get constant gravity if we will be assigning it
        if Planet.Do.CONSTANT_GRAVITY:
            Planet.g_ms2[Planet.Steps.nSurfIce:] = Constants.G * (Planet.Bulk.M_kg - MAbove_kg) / Planet.r_m[Planet.Steps.nSurfIce-1]**2
        else:
            # Ensure g values to be assigned are zero since we will be adding to them
            Planet.g_ms2[Planet.Steps.nSurfIce:] = 0
        # Assign 0 or 1 multiplier for constant/variable gravity calcs in loop
        VAR_GRAV = int(not Planet.Do.CONSTANT_GRAVITY)

        for i in range(Planet.Steps.nSurfIce, Planet.Steps.nSurfIce + Planet.Steps.nOceanMax):
            Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1]) * 1e6 / Planet.g_ms2[i-1] / \
                            Planet.rho_kgm3[i-1]
            Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
            Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
            MAbove_kg += Planet.MLayer_kg[i-1]
            MBelow_kg = Planet.Bulk.M_kg - MAbove_kg
            Planet.g_ms2[i] += VAR_GRAV * Constants.G * MBelow_kg / Planet.r_m[i]**2

        if Planet.Do.CLATHRATE:
            if Planet.Bulk.clathType == 'whole':
                zClathInfo = f', all clathrates.'
            elif Planet.Bulk.clathType == 'top':
                zClathInfo = f', including {Planet.zClath_m/1e3:.1f} km clathrate lid atop {Planet.zb_km - Planet.zClath_m/1e3:.1f} km ice Ih.'
            elif Planet.Bulk.clathType == 'bottom':
                # For underplate clathrates, Planet.zClath_m denotes the *thickness* of the layer rather than its starting depth.
                zClathInfo = f', including {Planet.zClath_m/1e3:.1f} km clathrate layer under {Planet.zb_km - Planet.zClath_m/1e3:.1f} km ice Ih.'
            else:
                raise ValueError(f'Bulk.clathType "{Planet.Bulk.clathType}" not recognized.')
        else:
            zClathInfo = '.'

        log.info(f'Ocean layers complete. zMax: {Planet.z_m[Planet.Steps.nSurfIce + Planet.Steps.nOceanMax - 1]/1e3:.1f} km, ' +
                 f'upper ice thickness zb: {Planet.zb_km:.3f} km{zClathInfo}')

        hpIceConvectionModel = ResolveHPIceConvectionModel(Planet)
        if hpIceConvectionModel != 'none':
            Planet = HPIceConvectionDiagnostics(Planet, Params, hpIceConvectionModel)

    return Planet, Params


def HPIceConvectionDiagnostics(Planet, Params, hpIceConvectionModel=None):
    """Compute opt-in diagnostics for in-ocean HP ice convection.

    This uses either Deschamps and Sotin (2001) or the opt-in Kalousova and
    Sotin (2018) scaling as a diagnostic calculator only. It does not modify
    the thermal, phase, mass, gravity, or heat-flux profile produced by the
    layer propagators.
    """
    if not Planet.Do.VALID or Planet.Do.NO_H2O:
        return Planet
    if Planet.Do.NO_OCEAN and not Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES:
        return Planet
    if Planet.Do.BOTTOM_ICEV or Planet.Do.BOTTOM_ICEIII:
        log.info('HP ice convection diagnostics skipped for underplate HP ice configuration.')
        return Planet
    if Planet.Do.NO_ICE_CONVECTION:
        log.info('HP ice convection diagnostics skipped because NO_ICE_CONVECTION is True.')
        return Planet
    if hpIceConvectionModel is None:
        hpIceConvectionModel = ResolveHPIceConvectionModel(Planet)
    useProductionDryRun = hpIceConvectionModel == 'Kalousova2018_production_experimental'
    useKalousova = hpIceConvectionModel in ('Kalousova2018_diagnostic', 'Kalousova2018_production_experimental')
    if hpIceConvectionModel not in ('DS2001_diagnostic', 'Kalousova2018_diagnostic',
                                    'Kalousova2018_production_experimental'):
        raise ValueError(f'HP ice convection model "{hpIceConvectionModel}" is not a recognized HP ice model.')
    productionMode = hpIceConvectionModel if useProductionDryRun else None

    phaseNames = {3: 'III', 5: 'V', 6: 'VI'}
    # Match genai's bottom-to-top HP ice traversal. The DS2001 diagnostic path
    # is phase-local, so this only affects logging/order until Kalousova is enabled.
    phaseOrder = (6, 5, 3)
    Planet.HPIceDiagnostics = {}
    Planet.DO_HP_MELT = False
    for phaseID, phaseName in phaseNames.items():
        _SetHPIceDiagnosticFields(Planet, phaseName, status='absent', phaseID=phaseID,
                                  productionMode=productionMode)

    iOceanStart = Planet.Steps.nSurfIce
    iOceanEnd = Planet.Steps.nSurfIce + Planet.Steps.nOceanMax
    oceanPhases = Planet.phase[iOceanStart:iOceanEnd]
    if not np.any(np.logical_and(oceanPhases > 1, oceanPhases < 10)):
        log.info('HP ice convection diagnostics enabled, but no in-ocean HP ice phases were found.')
        return Planet

    methodFamily = 'Kalousova and Sotin (2018)' if useKalousova else 'Deschamps and Sotin (2001)'
    if useProductionDryRun:
        methodFamily = 'Kalousova and Sotin (2018) experimental dry run'
    log.info(f'Computing opt-in HP ice convection diagnostics using {methodFamily}.')

    Qthrough_W = getattr(Planet.Ocean, 'QfromMantle_W', np.nan)
    if useKalousova and (Qthrough_W is None or not np.isfinite(Qthrough_W)):
        log.warning(
            'QfromMantle_W is not finite for Kalousova diagnostics. '
            'Using the local conductive heat-flux estimate inside the scaling law.'
        )
        Qthrough_W = np.nan

    for phaseID in phaseOrder:
        phaseName = phaseNames[phaseID]
        phaseInds = np.where(oceanPhases == phaseID)[0] + iOceanStart
        if len(phaseInds) == 0:
            log.info(f'HP ice {phaseName} diagnostics: phase absent.')
            continue

        iTop = phaseInds[0]
        iBot = phaseInds[-1]
        Ttop_K = Planet.T_K[iTop]
        Tb_K = Planet.T_K[iBot]
        rTop_m = Planet.r_m[iTop]
        rBot_m = Planet.r_m[iBot]
        zTop_m = Planet.z_m[iTop]
        zBot_m = Planet.z_m[iBot]
        kTop_WmK = Planet.kTherm_WmK[iTop]
        gtop_ms2 = Planet.g_ms2[iTop]
        zb_m = zBot_m - zTop_m
        Ptop_MPa = Planet.P_MPa[iTop]
        Pbot_MPa = Planet.P_MPa[iBot]
        Pmid_MPa = (Planet.P_MPa[iTop] + Planet.P_MPa[iBot]) / 2
        rhoPhase_kgm3 = _HPIcePhaseMean(Planet, 'rho_kgm3', phaseInds)
        CpPhase_JkgK = _HPIcePhaseMean(Planet, 'Cp_JkgK', phaseInds)
        alphaPhase_pK = _HPIcePhaseMean(Planet, 'alpha_pK', phaseInds)
        kThermPhase_WmK = _HPIcePhaseMean(Planet, 'kTherm_WmK', phaseInds)

        if zb_m < 1e3 or not np.all(np.isfinite([Ttop_K, Tb_K, rTop_m, kTop_WmK, gtop_ms2, Pmid_MPa])):
            reason = 'too thin' if zb_m < 1e3 else 'invalid P/T/r/k/g'
            _SetHPIceDiagnosticFields(
                Planet, phaseName, status=reason, phaseID=phaseID,
                iTop=iTop, iBot=iBot, rTop_m=rTop_m, rBot_m=rBot_m,
                zTop_m=zTop_m, zBot_m=zBot_m, thickness_m=zb_m,
                Ttop_K=Ttop_K, Tbot_K=Tb_K, Ptop_MPa=Ptop_MPa,
                Pmid_MPa=Pmid_MPa, Pbot_MPa=Pbot_MPa,
                rho_kgm3=rhoPhase_kgm3, Cp_JkgK=CpPhase_JkgK,
                alpha_pK=alphaPhase_pK, kTherm_WmK=kThermPhase_WmK,
                productionMode=productionMode,
            )
            log.info(f'HP ice {phaseName} diagnostics skipped: {reason}.')
            continue

        iceEOS = Planet.Ocean.iceEOS.get(phaseName)
        if iceEOS is None:
            _SetHPIceDiagnosticFields(
                Planet, phaseName, status='missing EOS', phaseID=phaseID,
                iTop=iTop, iBot=iBot, rTop_m=rTop_m, rBot_m=rBot_m,
                zTop_m=zTop_m, zBot_m=zBot_m, thickness_m=zb_m,
                Ttop_K=Ttop_K, Tbot_K=Tb_K, Ptop_MPa=Ptop_MPa,
                Pmid_MPa=Pmid_MPa, Pbot_MPa=Pbot_MPa,
                rho_kgm3=rhoPhase_kgm3, Cp_JkgK=CpPhase_JkgK,
                alpha_pK=alphaPhase_pK, kTherm_WmK=kThermPhase_WmK,
                productionMode=productionMode,
            )
            log.warning(f'HP ice {phaseName} diagnostics skipped: EOS was not loaded.')
            continue

        meltFraction = np.nan
        qBot_Wm2 = np.nan
        Q_in_W = np.nan
        Q_out_W = np.nan
        if useKalousova:
            method = 'Kalousova and Sotin (2018)'
            etaMelt_Pas = _GetKalousovaEtaMelt_Pas(Planet, phaseID, phaseName)
            if np.isfinite(Qthrough_W):
                Q_in_W = Qthrough_W
                qBot_Wm2 = Qthrough_W / (4*np.pi*rBot_m**2) if rBot_m > 0 else None
            else:
                qBot_Wm2 = None
            try:
                Tconv_K, etaConv_Pas, eLid_m, Dconv_m, deltaTBL_m, Qbot_W, Ra, RaCrit = \
                    ConvectionKalousova2018(
                        Ttop_K, rTop_m, kTop_WmK, Tb_K, zb_m, gtop_ms2,
                        Pmid_MPa, Planet.Ocean.EOS, iceEOS, phaseID,
                        Planet.Do.EQUIL_Q, Planet.Ocean.Eact_kJmol,
                        qBot_Wm2=qBot_Wm2, Htidal_Wm3=Planet.Ocean.HtidalIce_Wm3,
                        etaMelt_Pas=etaMelt_Pas
                    )
            except Exception as exc:
                _SetHPIceDiagnosticFields(
                    Planet, phaseName, status=f'error: {exc}', phaseID=phaseID,
                    iTop=iTop, iBot=iBot, rTop_m=rTop_m, rBot_m=rBot_m,
                    zTop_m=zTop_m, zBot_m=zBot_m, thickness_m=zb_m,
                    Ttop_K=Ttop_K, Tbot_K=Tb_K, qBot_Wm2=qBot_Wm2,
                    Ptop_MPa=Ptop_MPa, Pmid_MPa=Pmid_MPa,
                    Pbot_MPa=Pbot_MPa, method=method,
                    rho_kgm3=rhoPhase_kgm3, Cp_JkgK=CpPhase_JkgK,
                    alpha_pK=alphaPhase_pK, kTherm_WmK=kThermPhase_WmK,
                    productionMode=productionMode,
                )
                log.warning(f'HP ice {phaseName} Kalousova diagnostics failed: {exc}')
                continue

            isTemperate = (
                np.isfinite(Ra) and np.isfinite(RaCrit) and Ra > RaCrit and
                np.isfinite(eLid_m) and eLid_m > 0 and
                np.isfinite(Dconv_m) and Dconv_m > 0
            )
            meltFraction = Planet.Ocean.phiPercolationKalousova_frac if isTemperate else 0.0
            if isTemperate:
                Planet.DO_HP_MELT = True
            if np.isfinite(Qbot_W):
                Q_out_W = Qbot_W
                Qthrough_W = Qbot_W
        else:
            phaseEOS = _FixedPhaseEOS(Planet.Ocean.EOS, phaseID)
            method = 'Deschamps and Sotin (2001)'
            if phaseID == 6:
                method = 'DS2001 phase-local diagnostic fallback'
                log.info(
                    'HP ice VI diagnostics use the phase-local bottom-temperature fallback '
                    'because the production DS2001 melt lookup assumes a low-temperature '
                    'pure-water search range that is not valid for these Ice VI blocks.'
                )
                try:
                    Tconv_K, etaConv_Pas, eLid_m, Dconv_m, deltaTBL_m, Qbot_W, Ra, RaCrit = \
                        _ConvectionDeschampsSotinHPIceDiagnostic(
                            Ttop_K, rTop_m, kTop_WmK, Tb_K, zb_m, gtop_ms2,
                            Pmid_MPa, iceEOS, phaseID, Planet.Do.EQUIL_Q,
                            Planet.Ocean.Eact_kJmol
                        )
                except Exception as fallbackExc:
                    _SetHPIceDiagnosticFields(
                        Planet, phaseName, status=f'error: {fallbackExc}', phaseID=phaseID,
                        iTop=iTop, iBot=iBot, rTop_m=rTop_m, rBot_m=rBot_m,
                        zTop_m=zTop_m, zBot_m=zBot_m, thickness_m=zb_m,
                        Ttop_K=Ttop_K, Tbot_K=Tb_K, Ptop_MPa=Ptop_MPa,
                        Pmid_MPa=Pmid_MPa, Pbot_MPa=Pbot_MPa,
                        rho_kgm3=rhoPhase_kgm3, Cp_JkgK=CpPhase_JkgK,
                        alpha_pK=alphaPhase_pK, kTherm_WmK=kThermPhase_WmK,
                        method=method, productionMode=productionMode,
                    )
                    log.warning(f'HP ice {phaseName} diagnostics failed: {fallbackExc}')
                    continue
            else:
                try:
                    Tconv_K, etaConv_Pas, eLid_m, Dconv_m, deltaTBL_m, Qbot_W, Ra, RaCrit = \
                        ConvectionDeschampsSotin2001(Ttop_K, rTop_m, kTop_WmK, Tb_K, zb_m,
                                                     gtop_ms2, Pmid_MPa, phaseEOS,
                                                     iceEOS, phaseID, Planet.Do.EQUIL_Q,
                                                     Planet.Ocean.Eact_kJmol)
                except Exception as exc:
                    log.info(
                        f'HP ice {phaseName} DS2001 diagnostic lookup failed ({exc}); '
                        'using phase-local bottom-temperature fallback for diagnostics only.'
                    )
                    method = 'DS2001 phase-local diagnostic fallback'
                    try:
                        Tconv_K, etaConv_Pas, eLid_m, Dconv_m, deltaTBL_m, Qbot_W, Ra, RaCrit = \
                            _ConvectionDeschampsSotinHPIceDiagnostic(
                                Ttop_K, rTop_m, kTop_WmK, Tb_K, zb_m, gtop_ms2,
                                Pmid_MPa, iceEOS, phaseID, Planet.Do.EQUIL_Q,
                                Planet.Ocean.Eact_kJmol
                            )
                    except Exception as fallbackExc:
                        _SetHPIceDiagnosticFields(
                            Planet, phaseName, status=f'error: {fallbackExc}', phaseID=phaseID,
                            iTop=iTop, iBot=iBot, rTop_m=rTop_m, rBot_m=rBot_m,
                            zTop_m=zTop_m, zBot_m=zBot_m, thickness_m=zb_m,
                            Ttop_K=Ttop_K, Tbot_K=Tb_K, Ptop_MPa=Ptop_MPa,
                            Pmid_MPa=Pmid_MPa, Pbot_MPa=Pbot_MPa,
                            rho_kgm3=rhoPhase_kgm3, Cp_JkgK=CpPhase_JkgK,
                            alpha_pK=alphaPhase_pK, kTherm_WmK=kThermPhase_WmK,
                            method=method, productionMode=productionMode,
                        )
                        log.warning(f'HP ice {phaseName} diagnostics failed: {fallbackExc}')
                        continue

            etaMelt_Pas = Constants.etaMelt_Pas[phaseID]
            if np.isfinite(Qbot_W):
                Q_in_W = Qbot_W
                Q_out_W = Qbot_W

        meltCurveChecks = _GetIceVIMeltCurveCandidateChecks(
            Planet, phaseID, Ptop_MPa, Pmid_MPa, Pbot_MPa,
            Ttop_K, Tconv_K, Tb_K, productionMode,
        )
        _SetHPIceDiagnosticFields(
            Planet, phaseName, status='computed', phaseID=phaseID,
            iTop=iTop, iBot=iBot, rTop_m=rTop_m, rBot_m=rBot_m,
            zTop_m=zTop_m, zBot_m=zBot_m, thickness_m=zb_m,
            Ttop_K=Ttop_K, Tbot_K=Tb_K, Ptop_MPa=Ptop_MPa,
            Pmid_MPa=Pmid_MPa, Pbot_MPa=Pbot_MPa, qBot_Wm2=qBot_Wm2,
            Q_in_W=Q_in_W, Q_out_W=Q_out_W,
            rho_kgm3=rhoPhase_kgm3, Cp_JkgK=CpPhase_JkgK,
            alpha_pK=alphaPhase_pK, kTherm_WmK=kThermPhase_WmK,
            Tconv_K=Tconv_K, etaConv_Pas=etaConv_Pas,
            etaMelt_Pas=etaMelt_Pas,
            eLid_m=eLid_m, Dconv_m=Dconv_m,
            deltaTBL_m=deltaTBL_m, Ra=Ra, RaCrit=RaCrit,
            Qbot_W=Qbot_W, method=method,
            meltFraction=meltFraction, productionMode=productionMode,
            **meltCurveChecks,
        )
        meltInfo = f', meltFraction={meltFraction:.3f}' if np.isfinite(meltFraction) else ''
        log.info(
            f'HP ice {phaseName} diagnostics ({method}): Tconv={Tconv_K:.3f} K, '
            f'etaConv={etaConv_Pas:.3e} Pa s, eLid={eLid_m/1e3:.3f} km, '
            f'Dconv={Dconv_m/1e3:.3f} km, deltaTBL={deltaTBL_m/1e3:.3f} km, '
            f'Ra={Ra:.3e}, RaCrit={RaCrit:.3e}{meltInfo}.'
        )

    return Planet


def _GetKalousovaEtaMelt_Pas(Planet, phaseID, phaseName):
    """Resolve the explicit Kalousova melt-viscosity parameter for one phase."""
    etaParams = getattr(Planet.Ocean, 'etaMeltKalousova_Pas', None)
    etaMelt_Pas = None
    if isinstance(etaParams, dict):
        for key in (phaseName, phaseID, str(phaseID)):
            if key in etaParams:
                etaMelt_Pas = etaParams[key]
                break
    elif etaParams is not None:
        etaMelt_Pas = etaParams

    if etaMelt_Pas is None:
        etaMelt_Pas = Constants.etaMelt_Pas[phaseID]

    if not np.isfinite(etaMelt_Pas) or etaMelt_Pas <= 0:
        log.warning(
            f'Invalid etaMeltKalousova_Pas for ice {phaseName}: {etaMelt_Pas}. '
            'Falling back to Constants.etaMelt_Pas.'
        )
        etaMelt_Pas = Constants.etaMelt_Pas[phaseID]

    return etaMelt_Pas


def _ConvectionDeschampsSotinHPIceDiagnostic(Ttop_K, rTop_m, kTop_WmK, Tb_K, zb_m,
                                             gtop_ms2, Pmid_MPa, iceEOS, phaseID,
                                             EQUIL_Q, Eact_kJmol):
    """Phase-local DS2001-style calculator for HP ice diagnostics only.

    This mirrors the existing DS2001 scaling but avoids the melt-curve lookup
    used by the production helper. That lookup assumes a 274 K upper bound for
    pure water melting and can be outside the valid range for Ice VI diagnostic
    blocks. The fallback uses the block-bottom temperature as the local
    reference temperature and does not alter the propagated profile.
    """
    if Tb_K <= Ttop_K:
        raise ValueError('HP ice diagnostic bottom temperature is not warmer than the top temperature.')

    phaseString = PhaseConv(phaseID)
    if not np.isnan(Eact_kJmol[phaseString]):
        Eact_kJmol_use = Eact_kJmol[phaseString]
    else:
        Eact_kJmol_use = Constants.Eact_kJmol[phaseID]

    c1 = 1.43
    c2 = -0.03
    A = Eact_kJmol_use * 1e3 / Constants.R / Tb_K
    B = Eact_kJmol_use * 1e3 / 2 / Constants.R / c1
    C = c2 * (Tb_K - Ttop_K)
    Tconv_K = B * (np.sqrt(1 + 2/B*(Tb_K - C)) - 1)
    if Tconv_K < Ttop_K:
        Tconv_K = Ttop_K

    Tmelt_K = max(Tb_K, Tconv_K + 1e-6)
    etaMelt_Pas = Constants.etaMelt_Pas[phaseID]
    etaConv_Pas = etaMelt_Pas * np.exp(A * (Tmelt_K/Tconv_K - 1))

    rhoMid_kgm3 = iceEOS.fn_rho_kgm3(Pmid_MPa, Tconv_K)
    CpMid_JkgK = iceEOS.fn_Cp_JkgK(Pmid_MPa, Tconv_K)
    alphaMid_pK = iceEOS.fn_alpha_pK(Pmid_MPa, Tconv_K)
    kMid_WmK = iceEOS.fn_kTherm_WmK(Pmid_MPa, Tconv_K)
    if iceEOS.POROUS:
        log.warning('Porosity corrections are not applied in calculating HP ice diagnostic Rayleigh numbers.')

    rayleighInputs = [rhoMid_kgm3, CpMid_JkgK, alphaMid_pK, kMid_WmK, etaConv_Pas]
    if (not np.all(np.isfinite(rayleighInputs))) or alphaMid_pK <= 0 or kMid_WmK <= 0 or etaConv_Pas <= 0:
        raise ValueError('invalid HP ice diagnostic material properties')

    Ra = alphaMid_pK * CpMid_JkgK * rhoMid_kgm3**2 * gtop_ms2 * (Tb_K - Ttop_K) * zb_m**3 / etaConv_Pas / kMid_WmK
    if not np.isfinite(Ra) or Ra <= 0:
        raise ValueError('invalid HP ice diagnostic Rayleigh number')

    Radelta = 0.28 * Ra**0.21
    if Tb_K > Tconv_K:
        deltaTBL_m = (etaConv_Pas * kMid_WmK * Radelta /
                      alphaMid_pK / CpMid_JkgK / rhoMid_kgm3**2 /
                      gtop_ms2 / (Tb_K - Tconv_K))**(1/3)
        qBot_Wm2 = kMid_WmK * (Tb_K - Tconv_K) / deltaTBL_m
        qTop_Wm2 = (rTop_m - zb_m)**2 / rTop_m**2 * qBot_Wm2
        eLid_m = kTop_WmK * (Tconv_K - Ttop_K) / qTop_Wm2
    else:
        deltaTBL_m = 0.0
        eLid_m = zb_m

    RaCrit = GetRaCrit(Eact_kJmol_use, Tb_K, Ttop_K, Tconv_K)
    if Ra < RaCrit:
        eLid_m = zb_m
        deltaTBL_m = 0.0
        Tconv_K = Ttop_K

    if not EQUIL_Q:
        qBot_Wm2 = kMid_WmK * Tconv_K / eLid_m * np.log(Tb_K/Tconv_K)

    Dconv_m = zb_m - eLid_m - deltaTBL_m
    if Dconv_m < 0 and abs(Dconv_m) < 1e-6:
        Dconv_m = 0.0
    Qbot_W = qBot_Wm2 * 4*np.pi * (rTop_m - zb_m)**2

    return Tconv_K, etaConv_Pas, eLid_m, Dconv_m, deltaTBL_m, Qbot_W, Ra, RaCrit


class _FixedPhaseEOS:
    """Minimal EOS wrapper for phase-local diagnostics.

    The diagnostics evaluate each already-present HP ice block as its identified
    phase. This avoids modifying the model profile when a diagnostic convective
    temperature would otherwise trigger the ocean phase-boundary adjustment path.
    """

    def __init__(self, oceanEOS, phaseID):
        self.Tmin = getattr(oceanEOS, 'Tmin', -np.inf)
        self.phaseID = phaseID

    def fn_phase(self, P_MPa, T_K):
        return np.zeros_like(np.asarray(T_K), dtype=int) + self.phaseID


def _HPIcePhaseMean(Planet, field, indices):
    """Read one phase-local mean profile property without mutating the profile."""
    values = getattr(Planet, field, None)
    if values is None:
        return np.nan
    try:
        phaseValues = np.asarray(values, dtype=float)[indices]
    except (IndexError, TypeError, ValueError):
        return np.nan
    if phaseValues.size == 0 or np.all(~np.isfinite(phaseValues)):
        return np.nan
    return float(np.nanmean(phaseValues))


def _GetIceVIMeltCurveCandidateChecks(Planet, phaseID, Ptop_MPa, Pmid_MPa, Pbot_MPa,
                                      Ttop_K, Tconv_K, Tbot_K, productionMode,
                                      boundaryContext='in_run_candidate',
                                      candidateBoundarySource='in_run_candidate_bounds',
                                      finalProfileCoverageStatus='final_profile_not_evaluated'):
    """Read-only Ice VI melt-curve candidate checks for production dry runs."""
    checks = {
        'Tmelt_top_K': np.nan,
        'Tmelt_mid_K': np.nan,
        'Tmelt_bot_K': np.nan,
        'TmeltSource': None,
        'meltCurveStatus': 'not_evaluated',
        'phaseBoundaryResidual_K': np.nan,
        'phaseBoundaryStatus': 'not_evaluated',
        'eosBoundaryContext': 'not_evaluated',
        'eosBoundaryStatus': 'not_evaluated',
        'eosBoundaryReason': None,
        'candidateBoundarySource': None,
        'finalProfileCoverageStatus': finalProfileCoverageStatus,
    }
    if productionMode != "Kalousova2018_production_experimental" or phaseID != 6:
        return checks

    statusPrefix = 'posthoc_candidate' if boundaryContext == 'posthoc_final_profile' else 'in_run_candidate'
    checks.update({
        'eosBoundaryContext': boundaryContext,
        'candidateBoundarySource': candidateBoundarySource,
    })
    oceanEOS = getattr(getattr(Planet, 'Ocean', None), 'EOS', None)
    if oceanEOS is None or not hasattr(oceanEOS, 'fn_phase'):
        checks['meltCurveStatus'] = 'missing_melt_curve_rejected'
        checks['phaseBoundaryStatus'] = 'phase_boundary_unavailable_rejected'
        checks['eosBoundaryStatus'] = 'eos_boundary_unavailable'
        checks['eosBoundaryReason'] = 'missing_ocean_eos'
        return checks

    pressures_MPa = (Ptop_MPa, Pmid_MPa, Pbot_MPa)
    temperatures_K = (Ttop_K, Tconv_K, Tbot_K)
    if any(not np.isfinite(value) for value in pressures_MPa + temperatures_K):
        checks['meltCurveStatus'] = 'outside_eos_domain_rejected'
        checks['phaseBoundaryStatus'] = 'phase_boundary_unavailable_rejected'
        checks['eosBoundaryStatus'] = f'{statusPrefix}_boundary_nonfinite'
        checks['eosBoundaryReason'] = 'nonfinite_candidate_boundary'
        return checks

    Pmin_MPa = getattr(oceanEOS, 'Pmin', -np.inf)
    Pmax_MPa = getattr(oceanEOS, 'Pmax', np.inf)
    Tmin_K = getattr(oceanEOS, 'Tmin', -np.inf)
    Tmax_K = getattr(oceanEOS, 'Tmax', np.inf)
    if (
        min(pressures_MPa) < Pmin_MPa or max(pressures_MPa) > Pmax_MPa or
        min(temperatures_K) < Tmin_K or max(temperatures_K) > Tmax_K
    ):
        checks['meltCurveStatus'] = 'outside_eos_domain_rejected'
        checks['phaseBoundaryStatus'] = 'phase_boundary_unavailable_rejected'
        checks['eosBoundaryStatus'] = f'{statusPrefix}_boundary_outside_eos_domain'
        checks['eosBoundaryReason'] = 'candidate_boundary_outside_declared_eos_domain'
        return checks

    TsearchMax_K = np.nanmax((
        getattr(getattr(Planet, 'Ocean', None), 'THydroMax_K', np.nan),
        getattr(Planet, 'TfreezeUpper_K', np.nan),
        Tmax_K if np.isfinite(Tmax_K) else np.nan,
        max(temperatures_K) + 50.0,
    ))
    if not np.isfinite(TsearchMax_K):
        TsearchMax_K = max(temperatures_K) + 50.0
    TRes_K = getattr(Planet, 'TfreezeRes_K', 0.05)

    meltValues_K = []
    for P_MPa, T_K in zip(pressures_MPa, temperatures_K):
        try:
            phaseAtCandidate = int(np.asarray(oceanEOS.fn_phase(P_MPa, T_K)).item())
        except Exception:
            checks['meltCurveStatus'] = 'outside_eos_domain_rejected'
            checks['phaseBoundaryStatus'] = 'phase_boundary_unavailable_rejected'
            checks['eosBoundaryStatus'] = f'{statusPrefix}_phase_lookup_failed'
            checks['eosBoundaryReason'] = 'fn_phase_exception'
            return checks
        if phaseAtCandidate != 6:
            checks['meltCurveStatus'] = 'outside_eos_domain_rejected'
            checks['phaseBoundaryStatus'] = 'outside_eos_domain_rejected'
            checks['eosBoundaryStatus'] = f'{statusPrefix}_boundary_not_ice_vi'
            checks['eosBoundaryReason'] = 'candidate_boundary_not_ice_vi'
            return checks

        searchRange_K = max(TsearchMax_K - T_K, 1.0)
        try:
            Tmelt_K = GetTfreeze(
                oceanEOS, P_MPa, T_K,
                TfreezeRange_K=searchRange_K, TRes_K=TRes_K,
            )
        except Exception:
            checks['meltCurveStatus'] = 'missing_melt_curve_rejected'
            checks['phaseBoundaryStatus'] = 'phase_boundary_unavailable_rejected'
            checks['eosBoundaryStatus'] = 'melt_curve_lookup_failed'
            checks['eosBoundaryReason'] = 'GetTfreeze_exception'
            return checks
        meltValues_K.append(Tmelt_K)

    if any(not np.isfinite(value) for value in meltValues_K):
        checks['meltCurveStatus'] = 'melt_curve_nonfinite_rejected'
        checks['phaseBoundaryStatus'] = 'phase_boundary_unavailable_rejected'
        checks['eosBoundaryStatus'] = 'melt_curve_lookup_nonfinite'
        checks['eosBoundaryReason'] = 'nonfinite_GetTfreeze_result'
        return checks

    boundaryResidual_K = max(0.0, max(T_K - Tmelt_K for T_K, Tmelt_K in zip(temperatures_K, meltValues_K)))
    checks.update({
        'Tmelt_top_K': float(meltValues_K[0]),
        'Tmelt_mid_K': float(meltValues_K[1]),
        'Tmelt_bot_K': float(meltValues_K[2]),
        'TmeltSource': 'GetTfreeze',
        'meltCurveStatus': 'melt_curve_ok',
        'phaseBoundaryResidual_K': float(boundaryResidual_K),
        'phaseBoundaryStatus': 'phase_boundary_ok' if boundaryResidual_K <= 0.1 else 'phase_boundary_rejected',
        'eosBoundaryStatus': 'posthoc_candidate_in_domain' if boundaryContext == 'posthoc_final_profile' else 'in_run_candidate_boundary_in_domain',
        'eosBoundaryReason': 'candidate_boundary_in_declared_eos_domain',
        'finalProfileCoverageStatus': 'final_profile_nodes_in_domain' if boundaryContext == 'posthoc_final_profile' else finalProfileCoverageStatus,
    })
    return checks


def _PosthocIceVIResult(status, reason=None, blockers=None, **updates):
    result = {
        'posthocCandidateEvaluated': True,
        'posthocCandidateStatus': status,
        'posthocCandidateReason': reason or status,
        'posthocBoundarySource': 'finalized_profile_nodes',
        'posthocPtop_MPa': np.nan,
        'posthocPmid_MPa': np.nan,
        'posthocPbot_MPa': np.nan,
        'posthocTtop_K': np.nan,
        'posthocTmid_K': np.nan,
        'posthocTbot_K': np.nan,
        'posthocTmelt_top_K': np.nan,
        'posthocTmelt_mid_K': np.nan,
        'posthocTmelt_bot_K': np.nan,
        'posthocPhaseBoundaryResidual_K': np.nan,
        'posthocMeltCurveStatus': 'not_evaluated',
        'posthocPhaseBoundaryStatus': 'not_evaluated',
        'posthocUpdateAccepted': False,
        'posthocAcceptanceBlockers': tuple(blockers or (status,)),
        'posthocAllNodesInEOSDomain': False,
        'posthocAllNodesIceVI': False,
        'posthocEOSPressureMargin_MPa': np.nan,
        'posthocEOSTemperatureMargin_K': np.nan,
        'posthocMinPhaseBoundaryMargin_K': np.nan,
        'posthocAllNodeMeltCurveFailures': 0,
        'posthocTemperatureContrastStatus': 'not_evaluated',
        'posthocRayleighRegimeStatus': 'not_evaluated',
        'posthocThicknessStatus': 'not_evaluated',
        'posthocLayerClosureStatus': 'not_evaluated',
        'posthocViscosityStatus': 'not_evaluated',
        'posthocSensitivityRiskStatus': 'not_evaluated',
        'posthocRiskReasons': tuple(),
    }
    result.update(updates)
    return result


def _StorePosthocIceVIProductionCandidate(Planet, result):
    if not hasattr(Planet, 'HPIceDiagnostics') or Planet.HPIceDiagnostics is None:
        Planet.HPIceDiagnostics = {}
    Planet.HPIceDiagnostics.setdefault('VI', {})['posthocProductionCandidate'] = result
    for phaseName in ('III', 'V'):
        if phaseName in Planet.HPIceDiagnostics:
            Planet.HPIceDiagnostics[phaseName].setdefault('posthocProductionCandidate', _PosthocIceVIResult(
                'diagnostic_only_extrapolative',
                'posthoc_production_not_implemented_for_phase',
                blockers=('diagnostic_only_extrapolative',),
                posthocCandidateEvaluated=False,
            ))
    return result


def _PosthocFinitePositive(value):
    return value is not None and np.isfinite(value) and value > 0


def _PosthocMarginToBounds(values, lower, upper):
    values = np.asarray(values, dtype=float)
    margins = []
    if np.isfinite(lower):
        margins.append(values - lower)
    if np.isfinite(upper):
        margins.append(upper - values)
    if not margins:
        return np.inf
    margin = np.concatenate([np.asarray(item, dtype=float).reshape(-1) for item in margins])
    if margin.size == 0 or np.all(~np.isfinite(margin)):
        return np.nan
    return float(np.nanmin(margin))


def _PosthocIceVIAllNodeDiagnostics(Planet, iceVI, P_MPa, T_K):
    eos = getattr(getattr(Planet, 'Ocean', None), 'EOS', None)
    Pnodes_MPa = np.asarray(P_MPa, dtype=float)[iceVI]
    Tnodes_K = np.asarray(T_K, dtype=float)[iceVI]
    Pmin_MPa = getattr(eos, 'Pmin', -np.inf)
    Pmax_MPa = getattr(eos, 'Pmax', np.inf)
    Tmin_K = getattr(eos, 'Tmin', -np.inf)
    Tmax_K = getattr(eos, 'Tmax', np.inf)
    pressureInDomain = (Pnodes_MPa >= Pmin_MPa) & (Pnodes_MPa <= Pmax_MPa)
    temperatureInDomain = (Tnodes_K >= Tmin_K) & (Tnodes_K <= Tmax_K)
    finiteNodes = np.isfinite(Pnodes_MPa) & np.isfinite(Tnodes_K)
    phaseValues = []
    meltMargins_K = []
    meltCurveFailures = 0
    TsearchMax_K = np.nanmax((
        getattr(getattr(Planet, 'Ocean', None), 'THydroMax_K', np.nan),
        getattr(Planet, 'TfreezeUpper_K', np.nan),
        Tmax_K if np.isfinite(Tmax_K) else np.nan,
        np.nanmax(Tnodes_K) + 50.0 if Tnodes_K.size else np.nan,
    ))
    if not np.isfinite(TsearchMax_K) and Tnodes_K.size:
        TsearchMax_K = np.nanmax(Tnodes_K) + 50.0
    TRes_K = getattr(Planet, 'TfreezeRes_K', 0.05)
    for Pval_MPa, Tval_K, inDomain in zip(Pnodes_MPa, Tnodes_K, pressureInDomain & temperatureInDomain & finiteNodes):
        try:
            phaseAtNode = int(np.asarray(eos.fn_phase(float(Pval_MPa), float(Tval_K))).item())
        except Exception:
            phaseAtNode = None
        phaseValues.append(phaseAtNode)
        if not inDomain or phaseAtNode != 6:
            continue
        try:
            Tmelt_K = GetTfreeze(
                eos, float(Pval_MPa), float(Tval_K),
                TfreezeRange_K=max(TsearchMax_K - float(Tval_K), 1.0),
                TRes_K=TRes_K,
            )
        except Exception:
            meltCurveFailures += 1
            continue
        if np.isfinite(Tmelt_K):
            meltMargins_K.append(float(Tmelt_K - Tval_K))
        else:
            meltCurveFailures += 1
    allInDomain = bool(np.all(pressureInDomain & temperatureInDomain & finiteNodes))
    allIceVI = bool(all(value == 6 for value in phaseValues))
    return {
        'posthocAllNodesInEOSDomain': allInDomain,
        'posthocAllNodesIceVI': allIceVI,
        'posthocEOSPressureMargin_MPa': _PosthocMarginToBounds(Pnodes_MPa, Pmin_MPa, Pmax_MPa),
        'posthocEOSTemperatureMargin_K': _PosthocMarginToBounds(Tnodes_K, Tmin_K, Tmax_K),
        'posthocMinPhaseBoundaryMargin_K': float(np.nanmin(meltMargins_K)) if meltMargins_K else np.nan,
        'posthocAllNodeMeltCurveFailures': int(meltCurveFailures),
    }


def _PosthocIceVIRiskDiagnostics(Planet, iceVI, Pvals_MPa, Tvals_K, phaseDiagnostics, base):
    reasons = []
    fatalReasons = []

    if Tvals_K[2] <= Tvals_K[0]:
        fatalReasons.append('invalid_contrast')
    temperatureContrastStatus = 'invalid_contrast' if 'invalid_contrast' in fatalReasons else 'ok'

    ra = phaseDiagnostics.get('RaConvect', np.nan)
    raCrit = phaseDiagnostics.get('RaCrit', np.nan)
    if _PosthocFinitePositive(ra) and _PosthocFinitePositive(raCrit):
        raRatio = ra / raCrit
        if raRatio <= 1.0:
            fatalReasons.append('subcritical')
            rayleighStatus = 'subcritical'
        elif raRatio <= POSTHOC_RAYLEIGH_NEAR_CRITICAL_RATIO:
            reasons.append('near_critical')
            rayleighStatus = 'near_critical'
        else:
            rayleighStatus = 'supercritical'
    else:
        fatalReasons.append('subcritical')
        rayleighStatus = 'subcritical'

    thickness_m = np.nan
    if hasattr(Planet, 'z_m'):
        try:
            zNodes_m = np.asarray(Planet.z_m, dtype=float)[iceVI]
            thickness_m = float(np.nanmax(zNodes_m) - np.nanmin(zNodes_m))
        except (IndexError, TypeError, ValueError):
            thickness_m = np.nan
    if not np.isfinite(thickness_m):
        thickness_m = phaseDiagnostics.get('thickness_m', np.nan)
    if not _PosthocFinitePositive(thickness_m):
        fatalReasons.append('too_thin')
        thicknessStatus = 'too_thin'
    elif thickness_m < HP_ICE_PRODUCTION_MIN_THICKNESS_M:
        fatalReasons.append('too_thin')
        thicknessStatus = 'too_thin'
    else:
        thicknessStatus = 'ok'

    layerClosureResidual_m = phaseDiagnostics.get('layerClosureResidual_m', np.nan)
    if np.isfinite(layerClosureResidual_m) and np.isfinite(thickness_m) and thickness_m > 0:
        closureTol_m = max(
            HP_ICE_PRODUCTION_LAYER_CLOSURE_FRAC_TOL * thickness_m,
            HP_ICE_PRODUCTION_LAYER_CLOSURE_ABS_TOL_M,
        )
        if abs(layerClosureResidual_m) > closureTol_m:
            fatalReasons.append('boundary_layer_exceeds_layer')
            layerClosureStatus = 'boundary_layer_exceeds_layer'
        else:
            layerClosureStatus = 'ok'
    else:
        layerClosureStatus = 'not_evaluated'

    viscosityValues = (
        phaseDiagnostics.get('etaConv_Pas', np.nan),
        phaseDiagnostics.get('etaMelt_Pas', np.nan),
    )
    if any(value is not None and (not np.isfinite(value) or value <= 0) for value in viscosityValues):
        fatalReasons.append('invalid_viscosity')
        viscosityStatus = 'invalid_viscosity'
    else:
        viscosityStatus = 'ok'

    if base.get('posthocEOSPressureMargin_MPa', np.inf) <= POSTHOC_EOS_PRESSURE_MARGIN_WARN_MPA:
        reasons.append('near_eos_boundary')
    if base.get('posthocEOSTemperatureMargin_K', np.inf) <= POSTHOC_EOS_TEMPERATURE_MARGIN_WARN_K:
        reasons.append('near_eos_boundary')
    if base.get('posthocMinPhaseBoundaryMargin_K', np.inf) <= POSTHOC_PHASE_BOUNDARY_MARGIN_WARN_K:
        reasons.append('near_phase_boundary')

    if fatalReasons:
        sensitivityRiskStatus = 'high_risk_rejected'
    elif 'near_phase_boundary' in reasons:
        sensitivityRiskStatus = 'near_phase_boundary'
    elif 'near_eos_boundary' in reasons:
        sensitivityRiskStatus = 'near_eos_boundary'
    elif 'near_critical' in reasons:
        sensitivityRiskStatus = 'near_critical'
    elif reasons:
        sensitivityRiskStatus = 'requires_science_review'
    else:
        sensitivityRiskStatus = 'nominal'

    return {
        'posthocTemperatureContrastStatus': temperatureContrastStatus,
        'posthocRayleighRegimeStatus': rayleighStatus,
        'posthocThicknessStatus': thicknessStatus,
        'posthocLayerClosureStatus': layerClosureStatus,
        'posthocViscosityStatus': viscosityStatus,
        'posthocSensitivityRiskStatus': sensitivityRiskStatus,
        'posthocRiskReasons': tuple(fatalReasons + reasons),
    }


def EvaluatePosthocIceVIProductionCandidate(Planet, productionMode=None):
    """Evaluate finalized-profile Ice VI production candidacy without mutation."""
    if productionMode is None:
        productionMode = ResolveHPIceConvectionModel(Planet) if hasattr(Planet, 'Do') else "Kalousova2018_production_experimental"
    if productionMode != "Kalousova2018_production_experimental":
        return _StorePosthocIceVIProductionCandidate(
            Planet,
            _PosthocIceVIResult('posthoc_not_enabled', 'experimental_production_selector_not_enabled'),
        )
    diagnostics = getattr(Planet, 'HPIceDiagnostics', None)
    if isinstance(diagnostics, dict):
        iceVIDiagnostics = diagnostics.get('VI')
        if isinstance(iceVIDiagnostics, dict):
            candidate = iceVIDiagnostics.get('activeProductionCandidate')
            posthoc = iceVIDiagnostics.get('posthocProductionCandidate')
            if (
                isinstance(candidate, dict) and candidate.get('candidateDiscarded', False) and
                isinstance(posthoc, dict)
            ):
                return posthoc

    try:
        phase = np.asarray(getattr(Planet, 'phase'))
        P_MPa = np.asarray(getattr(Planet, 'P_MPa'), dtype=float)
        T_K = np.asarray(getattr(Planet, 'T_K'), dtype=float)
    except (AttributeError, TypeError, ValueError):
        return _StorePosthocIceVIProductionCandidate(
            Planet,
            _PosthocIceVIResult('posthoc_nonfinite_rejected', 'missing_profile_arrays'),
        )

    iceVI = np.where(phase == 6)[0]
    if iceVI.size < 2:
        return _StorePosthocIceVIProductionCandidate(
            Planet,
            _PosthocIceVIResult('posthoc_missing_ice_vi', 'requires_at_least_two_finalized_ice_vi_nodes'),
        )
    try:
        selected = np.array([iceVI[0], iceVI[iceVI.size // 2], iceVI[-1]], dtype=int)
        Pvals_MPa = P_MPa[selected]
        Tvals_K = T_K[selected]
        phaseVals = phase[selected]
    except (IndexError, TypeError, ValueError):
        return _StorePosthocIceVIProductionCandidate(
            Planet,
            _PosthocIceVIResult('posthoc_nonfinite_rejected', 'invalid_profile_array_shape'),
        )

    base = {
        'posthocPtop_MPa': float(Pvals_MPa[0]) if np.isfinite(Pvals_MPa[0]) else np.nan,
        'posthocPmid_MPa': float(Pvals_MPa[1]) if np.isfinite(Pvals_MPa[1]) else np.nan,
        'posthocPbot_MPa': float(Pvals_MPa[2]) if np.isfinite(Pvals_MPa[2]) else np.nan,
        'posthocTtop_K': float(Tvals_K[0]) if np.isfinite(Tvals_K[0]) else np.nan,
        'posthocTmid_K': float(Tvals_K[1]) if np.isfinite(Tvals_K[1]) else np.nan,
        'posthocTbot_K': float(Tvals_K[2]) if np.isfinite(Tvals_K[2]) else np.nan,
    }
    if (
        np.any(~np.isfinite(Pvals_MPa)) or np.any(~np.isfinite(Tvals_K)) or
        np.any(phaseVals != 6)
    ):
        return _StorePosthocIceVIProductionCandidate(
            Planet,
            _PosthocIceVIResult('posthoc_nonfinite_rejected', 'nonfinite_or_wrong_phase_finalized_nodes', **base),
        )

    base.update(_PosthocIceVIAllNodeDiagnostics(Planet, iceVI, P_MPa, T_K))
    if not base['posthocAllNodesInEOSDomain']:
        base['posthocMeltCurveStatus'] = 'outside_eos_domain_rejected'
        base['posthocPhaseBoundaryStatus'] = 'phase_boundary_unavailable_rejected'
        return _StorePosthocIceVIProductionCandidate(
            Planet,
            _PosthocIceVIResult('posthoc_outside_eos_domain_rejected', 'posthoc_all_nodes_outside_eos_domain', **base),
        )
    if not base['posthocAllNodesIceVI']:
        base['posthocMeltCurveStatus'] = 'outside_eos_domain_rejected'
        base['posthocPhaseBoundaryStatus'] = 'outside_eos_domain_rejected'
        return _StorePosthocIceVIProductionCandidate(
            Planet,
            _PosthocIceVIResult('posthoc_outside_eos_domain_rejected', 'posthoc_all_nodes_not_ice_vi', **base),
        )
    if (
        base.get('posthocAllNodeMeltCurveFailures', 0) > 0 or
        not np.isfinite(base.get('posthocMinPhaseBoundaryMargin_K', np.nan))
    ):
        base['posthocMeltCurveStatus'] = 'melt_curve_nonfinite_rejected'
        base['posthocPhaseBoundaryStatus'] = 'phase_boundary_unavailable_rejected'
        return _StorePosthocIceVIProductionCandidate(
            Planet,
            _PosthocIceVIResult('posthoc_nonfinite_rejected', 'posthoc_all_node_melt_curve_nonfinite_or_unavailable', **base),
        )

    meltChecks = _GetIceVIMeltCurveCandidateChecks(
        Planet, 6,
        Pvals_MPa[0], Pvals_MPa[1], Pvals_MPa[2],
        Tvals_K[0], Tvals_K[1], Tvals_K[2],
        productionMode,
        boundaryContext='posthoc_final_profile',
        candidateBoundarySource='finalized_profile_nodes',
        finalProfileCoverageStatus='final_profile_not_evaluated',
    )
    base.update({
        'posthocTmelt_top_K': meltChecks['Tmelt_top_K'],
        'posthocTmelt_mid_K': meltChecks['Tmelt_mid_K'],
        'posthocTmelt_bot_K': meltChecks['Tmelt_bot_K'],
        'posthocPhaseBoundaryResidual_K': meltChecks['phaseBoundaryResidual_K'],
        'posthocMeltCurveStatus': meltChecks['meltCurveStatus'],
        'posthocPhaseBoundaryStatus': meltChecks['phaseBoundaryStatus'],
        'eosBoundaryContext': meltChecks['eosBoundaryContext'],
        'eosBoundaryStatus': meltChecks['eosBoundaryStatus'],
        'eosBoundaryReason': meltChecks['eosBoundaryReason'],
        'candidateBoundarySource': meltChecks['candidateBoundarySource'],
        'finalProfileCoverageStatus': meltChecks['finalProfileCoverageStatus'],
    })
    if meltChecks['meltCurveStatus'] == 'outside_eos_domain_rejected':
        return _StorePosthocIceVIProductionCandidate(
            Planet,
            _PosthocIceVIResult('posthoc_outside_eos_domain_rejected', meltChecks['eosBoundaryStatus'], **base),
        )
    if meltChecks['meltCurveStatus'] == 'missing_melt_curve_rejected':
        return _StorePosthocIceVIProductionCandidate(
            Planet,
            _PosthocIceVIResult('posthoc_missing_melt_curve_rejected', meltChecks['eosBoundaryReason'], **base),
        )
    if meltChecks['meltCurveStatus'] == 'melt_curve_nonfinite_rejected':
        return _StorePosthocIceVIProductionCandidate(
            Planet,
            _PosthocIceVIResult('posthoc_nonfinite_rejected', meltChecks['eosBoundaryReason'], **base),
        )
    if meltChecks['phaseBoundaryStatus'] != 'phase_boundary_ok':
        return _StorePosthocIceVIProductionCandidate(
            Planet,
            _PosthocIceVIResult('posthoc_phase_boundary_rejected', meltChecks['phaseBoundaryStatus'], **base),
        )

    phaseDiagnostics = getattr(Planet, 'HPIceDiagnostics', {}).get('VI', {})
    base.update(_PosthocIceVIRiskDiagnostics(Planet, iceVI, Pvals_MPa, Tvals_K, phaseDiagnostics, base))
    if base['posthocSensitivityRiskStatus'] == 'high_risk_rejected':
        return _StorePosthocIceVIProductionCandidate(
            Planet,
            _PosthocIceVIResult(
                'posthoc_high_risk_rejected',
                base['posthocRiskReasons'][0] if base['posthocRiskReasons'] else 'high_risk_rejected',
                blockers=base['posthocRiskReasons'],
                **base,
            ),
        )

    energyResidual_W = phaseDiagnostics.get('energyResidual_W', np.nan)
    energyResidual_frac = phaseDiagnostics.get('energyResidual_frac', np.nan)
    Q_in_W = phaseDiagnostics.get('Q_in_W', 1.0)
    Q_out_W = phaseDiagnostics.get('Q_out_W', 1.0)
    if Q_in_W is None or not np.isfinite(Q_in_W):
        Q_in_W = 1.0
    if Q_out_W is None or not np.isfinite(Q_out_W):
        Q_out_W = 1.0
    Qscale_W = max(abs(Q_in_W), abs(Q_out_W), 1.0)
    energyAbsTol_W = max(HP_ICE_PRODUCTION_ENERGY_FRAC_TOL * Qscale_W, HP_ICE_PRODUCTION_ENERGY_ABS_FLOOR_W)
    base.update({
        'posthocEnergyResidual_W': energyResidual_W,
        'posthocEnergyResidual_frac': energyResidual_frac,
    })
    if (
        not np.isfinite(energyResidual_W) or not np.isfinite(energyResidual_frac) or
        abs(energyResidual_W) > energyAbsTol_W or abs(energyResidual_frac) > HP_ICE_PRODUCTION_ENERGY_FRAC_TOL
    ):
        return _StorePosthocIceVIProductionCandidate(
            Planet,
            _PosthocIceVIResult('posthoc_energy_residual_rejected', 'energy_residual_exceeds_tolerance', **base),
        )

    massResidual_kg = phaseDiagnostics.get('massResidual_kg', 0.0)
    massResidual_frac = phaseDiagnostics.get('massResidual_frac', 0.0)
    if massResidual_kg is None or not np.isfinite(massResidual_kg):
        massResidual_kg = 0.0
    if massResidual_frac is None or not np.isfinite(massResidual_frac):
        massResidual_frac = 0.0
    base.update({
        'posthocMassResidual_kg': massResidual_kg,
        'posthocMassResidual_frac': massResidual_frac,
    })
    if massResidual_kg != 0.0 or massResidual_frac != 0.0:
        return _StorePosthocIceVIProductionCandidate(
            Planet,
            _PosthocIceVIResult('posthoc_mass_residual_rejected', 'mass_residual_nonzero', **base),
        )

    return _StorePosthocIceVIProductionCandidate(
        Planet,
        _PosthocIceVIResult(
            'posthoc_candidate_passed',
            'posthoc_candidate_state_only',
            blockers=(),
            posthocUpdateAccepted=True,
            posthocAcceptanceBlockers=(),
            **base,
        ),
    )


def _ActiveIceVIProtectedProfileSnapshot(Planet):
    snapshot = {}
    for field in (
        'P_MPa', 'T_K', 'rho_kgm3', 'phase', 'eta_Pas',
        'qSurf_Wm2', 'qCon_Wm2', 'Mtot_kg', 'CMR2mean',
    ):
        if not hasattr(Planet, field):
            snapshot[field] = None
            continue
        value = getattr(Planet, field)
        snapshot[field] = value.copy() if isinstance(value, np.ndarray) else value
    return snapshot


def _ActiveIceVIProtectedProfileUnchanged(Planet, snapshot):
    for field, before in snapshot.items():
        if before is None:
            if hasattr(Planet, field):
                return False
            continue
        if not hasattr(Planet, field):
            return False
        after = getattr(Planet, field)
        if isinstance(before, np.ndarray):
            if not np.array_equal(after, before):
                return False
        elif after != before:
            return False
    return True


def _ActiveIceVIRejectedResult(status, reason, protectedFieldsUnchanged=True):
    return {
        'candidateCopyCreated': False,
        'candidateCopySource': 'finalized_posthoc_ice_vi_nodes',
        'candidatePhase': 'VI',
        'candidateNodeCount': 0,
        'candidateIndexStart': None,
        'candidateIndexEnd': None,
        'candidateP_MPa': np.array([], dtype=float),
        'candidateT_K': np.array([], dtype=float),
        'candidatePhaseArray': np.array([], dtype=int),
        'candidateR_m': np.array([], dtype=float),
        'candidateZ_m': np.array([], dtype=float),
        'candidateEta_Pas': np.array([], dtype=float),
        'candidateQ_in_W': np.nan,
        'candidateQ_out_W': np.nan,
        'candidateQbot_W': np.nan,
        'candidateq_in_Wm2': np.nan,
        'candidateq_out_Wm2': np.nan,
        'candidateAppliedToProfile': False,
        'candidateAccepted': False,
        'candidateDiscarded': False,
        'candidateDiscardReason': None,
        'candidateResidualsPassed': False,
        'candidateHeatPowerResidual_W': np.nan,
        'candidateEnergyResidual_W': np.nan,
        'candidateEnergyResidual_frac': np.nan,
        'candidateHeatFluxResidual_Wm2': np.nan,
        'candidateExpected_q_in_Wm2': np.nan,
        'candidateExpected_q_out_Wm2': np.nan,
        'candidateSphericalFluxScalingStatus': 'not_evaluated',
        'candidateMassResidual_kg': np.nan,
        'candidateMassResidual_frac': np.nan,
        'candidatePhaseBoundaryResidual_K': np.nan,
        'candidateLayerClosureResidual_m': np.nan,
        'candidateEOSPressureMargin_MPa': np.nan,
        'candidateEOSTemperatureMargin_K': np.nan,
        'candidateMinPhaseBoundaryMargin_K': np.nan,
        'candidateRaOverRaCrit': np.nan,
        'candidateResidualStatus': 'candidate_not_evaluated',
        'candidateResidualReasons': tuple(),
        'candidateStatus': status,
        'candidateReason': reason,
        'rollbackRequired': False,
        'rollbackApplied': False,
        'rollbackReason': None,
        'rollbackReasons': tuple(),
        'rollbackStatus': 'rollback_not_evaluated',
        'rollbackPolicyVersion': ACTIVE_ICE_VI_ROLLBACK_POLICY_VERSION,
        'thermalUpdate': {
            'candidateThermalUpdateStrategy': ACTIVE_ICE_VI_THERMAL_UPDATE_STRATEGY_NO_OP,
            'candidateThermalUpdateStatus': 'candidate_thermal_update_not_evaluated',
            'candidateThermalUpdateReasons': tuple(),
            'candidateThermalUpdateAppliedToCopy': False,
            'candidateThermalUpdateAppliedToPlanet': False,
            'candidateThermalUpdateAccepted': False,
        },
        'protectedFieldsUnchanged': protectedFieldsUnchanged,
    }


def _StoreActiveIceVIProductionCandidateCopy(Planet, result):
    if not hasattr(Planet, 'HPIceDiagnostics') or Planet.HPIceDiagnostics is None:
        Planet.HPIceDiagnostics = {}
    Planet.HPIceDiagnostics.setdefault('VI', {})['activeProductionCandidate'] = result
    for phaseName in ('III', 'V'):
        if phaseName in Planet.HPIceDiagnostics:
            Planet.HPIceDiagnostics[phaseName].setdefault('activeProductionCandidate', {
                'candidateCopyCreated': False,
                'candidateCopySource': 'not_applicable',
                'candidatePhase': phaseName,
                'candidateNodeCount': 0,
                'candidateAppliedToProfile': False,
                'candidateAccepted': False,
                'candidateStatus': 'diagnostic_only_extrapolative',
                'candidateReason': 'active_production_not_implemented_for_phase',
                'protectedFieldsUnchanged': True,
            })
    return result


def BuildActiveIceVIProductionCandidateCopy(Planet, productionMode=None):
    """Build an isolated Ice VI active-production candidate copy without mutation."""
    protectedBefore = _ActiveIceVIProtectedProfileSnapshot(Planet)
    if productionMode is None:
        productionMode = ResolveHPIceConvectionModel(Planet) if hasattr(Planet, 'Do') else "Kalousova2018_production_experimental"
    if productionMode != "Kalousova2018_production_experimental":
        return _StoreActiveIceVIProductionCandidateCopy(
            Planet,
            _ActiveIceVIRejectedResult(
                'candidate_copy_not_enabled',
                'experimental_production_selector_not_enabled',
                _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore),
            ),
        )

    diagnostics = getattr(Planet, 'HPIceDiagnostics', None)
    if not isinstance(diagnostics, dict):
        return _StoreActiveIceVIProductionCandidateCopy(
            Planet,
            _ActiveIceVIRejectedResult(
                'missing_diagnostics_rejected',
                'missing_hp_ice_diagnostics',
                _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore),
            ),
        )
    iceVIDiagnostics = diagnostics.get('VI')
    if not isinstance(iceVIDiagnostics, dict):
        return _StoreActiveIceVIProductionCandidateCopy(
            Planet,
            _ActiveIceVIRejectedResult(
                'missing_ice_vi_rejected',
                'missing_ice_vi_diagnostics',
                _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore),
            ),
        )
    candidate = iceVIDiagnostics.get('activeProductionCandidate')
    if isinstance(candidate, dict) and candidate.get('candidateDiscarded', False):
        return _ActiveIceVIPreserveTerminalDiscard(candidate)

    posthoc = iceVIDiagnostics.get('posthocProductionCandidate')
    if not isinstance(posthoc, dict) or not posthoc.get('posthocCandidateEvaluated', False):
        return _StoreActiveIceVIProductionCandidateCopy(
            Planet,
            _ActiveIceVIRejectedResult(
                'missing_posthoc_candidate_rejected',
                'requires_finalized_posthoc_ice_vi_candidate',
                _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore),
            ),
        )
    if posthoc.get('posthocSensitivityRiskStatus') == 'high_risk_rejected':
        return _StoreActiveIceVIProductionCandidateCopy(
            Planet,
            _ActiveIceVIRejectedResult(
                'high_risk_posthoc_rejected',
                'posthoc_candidate_is_high_risk',
                _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore),
            ),
        )
    if posthoc.get('posthocCandidateStatus') != 'posthoc_candidate_passed':
        return _StoreActiveIceVIProductionCandidateCopy(
            Planet,
            _ActiveIceVIRejectedResult(
                'posthoc_candidate_not_passed_rejected',
                posthoc.get('posthocCandidateReason', 'posthoc_candidate_not_passed'),
                _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore),
            ),
        )

    try:
        phase = np.asarray(getattr(Planet, 'phase'))
        P_MPa = np.asarray(getattr(Planet, 'P_MPa'), dtype=float)
        T_K = np.asarray(getattr(Planet, 'T_K'), dtype=float)
    except (AttributeError, TypeError, ValueError):
        return _StoreActiveIceVIProductionCandidateCopy(
            Planet,
            _ActiveIceVIRejectedResult(
                'missing_profile_arrays_rejected',
                'requires_finalized_profile_arrays',
                _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore),
            ),
        )

    iceVI = np.where(phase == 6)[0]
    if iceVI.size < 2:
        return _StoreActiveIceVIProductionCandidateCopy(
            Planet,
            _ActiveIceVIRejectedResult(
                'missing_ice_vi_rejected',
                'requires_at_least_two_finalized_ice_vi_nodes',
                _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore),
            ),
        )
    if np.any(phase[iceVI] != 6):
        return _StoreActiveIceVIProductionCandidateCopy(
            Planet,
            _ActiveIceVIRejectedResult(
                'non_ice_vi_nodes_rejected',
                'candidate_copy_requires_phase_6_nodes',
                _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore),
            ),
        )

    def copy_optional_array(field, dtype=float):
        if not hasattr(Planet, field):
            return np.array([], dtype=dtype)
        try:
            return np.asarray(getattr(Planet, field), dtype=dtype)[iceVI].copy()
        except (IndexError, TypeError, ValueError):
            return np.array([], dtype=dtype)

    result = {
        'candidateCopyCreated': True,
        'candidateCopySource': 'finalized_posthoc_ice_vi_nodes',
        'candidatePhase': 'VI',
        'candidateNodeCount': int(iceVI.size),
        'candidateIndexStart': int(iceVI[0]),
        'candidateIndexEnd': int(iceVI[-1]),
        'candidateP_MPa': P_MPa[iceVI].copy(),
        'candidateT_K': T_K[iceVI].copy(),
        'candidatePhaseArray': phase[iceVI].copy(),
        'candidateR_m': copy_optional_array('r_m'),
        'candidateZ_m': copy_optional_array('z_m'),
        'candidateEta_Pas': copy_optional_array('eta_Pas'),
        'candidateQ_in_W': iceVIDiagnostics.get('Q_in_W', np.nan),
        'candidateQ_out_W': iceVIDiagnostics.get('Q_out_W', np.nan),
        'candidateQbot_W': iceVIDiagnostics.get('Qbot_W', np.nan),
        'candidateq_in_Wm2': iceVIDiagnostics.get('q_in_Wm2', np.nan),
        'candidateq_out_Wm2': iceVIDiagnostics.get('q_out_Wm2', np.nan),
        'candidateAppliedToProfile': False,
        'candidateAccepted': False,
        'candidateStatus': 'candidate_copy_created',
        'candidateReason': 'candidate_copy_state_only',
        'rollbackRequired': False,
        'rollbackApplied': False,
        'rollbackReason': None,
        'protectedFieldsUnchanged': _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore),
    }
    return _StoreActiveIceVIProductionCandidateCopy(Planet, result)


def _ActiveIceVIFiniteScalar(value):
    return value is not None and np.isfinite(value)


def _ActiveIceVIFinitePositive(value):
    return _ActiveIceVIFiniteScalar(value) and value > 0


def _ActiveIceVIResidualStatus(candidate, passed, status, reasons):
    candidate['candidateResidualsPassed'] = bool(passed)
    candidate['candidateResidualStatus'] = status
    candidate['candidateResidualReasons'] = tuple(reasons)
    candidate['candidateAccepted'] = False
    candidate['candidateAppliedToProfile'] = False
    return candidate


def _ActiveIceVIReasonTuple(*values):
    reasons = []
    for value in values:
        if value is None:
            continue
        if isinstance(value, (tuple, list)):
            items = value
        else:
            items = (value,)
        for item in items:
            if item is None or item == '':
                continue
            if item not in reasons:
                reasons.append(item)
    return tuple(reasons)


def _ActiveIceVIRollbackStatus(candidate, *, status, required, discarded, applied,
                               reason=None, reasons=()):
    candidate['candidateAccepted'] = False
    candidate['candidateAppliedToProfile'] = False
    candidate['candidateDiscarded'] = bool(discarded)
    candidate['candidateDiscardReason'] = reason if discarded else None
    candidate['rollbackStatus'] = status
    candidate['rollbackRequired'] = bool(required)
    candidate['rollbackApplied'] = bool(applied)
    candidate['rollbackReason'] = reason
    candidate['rollbackReasons'] = tuple(reasons)
    candidate['rollbackPolicyVersion'] = ACTIVE_ICE_VI_ROLLBACK_POLICY_VERSION
    return candidate


def _ActiveIceVIPreserveTerminalDiscard(candidate):
    """Keep an already-discarded candidate terminal until an explicit reset exists."""
    residualStatus = candidate.get('candidateResidualStatus', 'candidate_not_evaluated')
    residualReasons = tuple(candidate.get('candidateResidualReasons', ()))
    existingStatus = candidate.get('rollbackStatus')
    existingReason = candidate.get('rollbackReason')
    existingReasons = tuple(candidate.get('rollbackReasons', ()))
    discardReason = candidate.get('candidateDiscardReason')

    if existingStatus in (None, '', 'rollback_not_evaluated'):
        if candidate.get('candidateCopyCreated', False):
            rollbackStatus = 'rollback_candidate_discarded'
            rollbackApplied = True
            reason = discardReason or existingReason or residualStatus or 'candidate_discarded'
            reasons = _ActiveIceVIReasonTuple(existingReasons, residualReasons, reason)
        else:
            rollbackStatus = 'rollback_missing_candidate'
            rollbackApplied = False
            reason = discardReason or existingReason or 'missing_candidate_copy'
            reasons = _ActiveIceVIReasonTuple(existingReasons, 'missing_candidate_copy', residualReasons, reason)
    else:
        rollbackStatus = existingStatus
        rollbackApplied = bool(candidate.get('rollbackApplied', rollbackStatus == 'rollback_candidate_discarded'))
        reason = existingReason or discardReason or residualStatus or candidate.get('candidateReason')
        reasons = _ActiveIceVIReasonTuple(existingReasons) if existingReasons else _ActiveIceVIReasonTuple(residualReasons, reason)

    if not reasons and reason is not None:
        reasons = (reason,)
    candidate['candidateAccepted'] = False
    candidate['candidateAppliedToProfile'] = False
    candidate['candidateDiscarded'] = True
    candidate['candidateDiscardReason'] = discardReason or reason
    candidate['rollbackStatus'] = rollbackStatus
    candidate['rollbackRequired'] = True
    candidate['rollbackApplied'] = rollbackApplied
    candidate['rollbackReason'] = reason
    candidate['rollbackReasons'] = reasons
    candidate['rollbackPolicyVersion'] = ACTIVE_ICE_VI_ROLLBACK_POLICY_VERSION
    return candidate


def EvaluateActiveIceVICandidateResiduals(Planet, productionMode=None):
    """Evaluate conservation residuals on the isolated Ice VI candidate copy only."""
    protectedBefore = _ActiveIceVIProtectedProfileSnapshot(Planet)
    if productionMode is None:
        productionMode = ResolveHPIceConvectionModel(Planet) if hasattr(Planet, 'Do') else "Kalousova2018_production_experimental"

    diagnostics = getattr(Planet, 'HPIceDiagnostics', None)
    if not isinstance(diagnostics, dict):
        result = _ActiveIceVIRejectedResult(
            'missing_diagnostics_rejected',
            'missing_hp_ice_diagnostics',
            _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore),
        )
        _ActiveIceVIResidualStatus(
            result, False, 'candidate_missing_required_field_rejected',
            ('missing_hp_ice_diagnostics',),
        )
        return _StoreActiveIceVIProductionCandidateCopy(Planet, result)

    iceVIDiagnostics = diagnostics.get('VI')
    if not isinstance(iceVIDiagnostics, dict):
        result = _ActiveIceVIRejectedResult(
            'missing_ice_vi_rejected',
            'missing_ice_vi_diagnostics',
            _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore),
        )
        _ActiveIceVIResidualStatus(
            result, False, 'candidate_missing_required_field_rejected',
            ('missing_ice_vi_diagnostics',),
        )
        return _StoreActiveIceVIProductionCandidateCopy(Planet, result)

    candidate = iceVIDiagnostics.get('activeProductionCandidate')
    if not isinstance(candidate, dict):
        candidate = _ActiveIceVIRejectedResult(
            'missing_candidate_copy_rejected',
            'requires_active_candidate_copy',
            _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore),
        )
        _StoreActiveIceVIProductionCandidateCopy(Planet, candidate)
    elif candidate.get('candidateDiscarded', False):
        return _ActiveIceVIPreserveTerminalDiscard(candidate)

    candidate['candidateAccepted'] = False
    candidate['candidateAppliedToProfile'] = False
    candidate['protectedFieldsUnchanged'] = _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore)

    if productionMode != "Kalousova2018_production_experimental":
        return _ActiveIceVIResidualStatus(
            candidate, False, 'candidate_not_evaluated',
            ('experimental_production_selector_not_enabled',),
        )
    if not candidate.get('candidateCopyCreated', False):
        return _ActiveIceVIResidualStatus(
            candidate, False, 'candidate_missing_required_field_rejected',
            (candidate.get('candidateReason', 'requires_active_candidate_copy'),),
        )

    failures = []

    def add_failure(status, reason):
        failures.append((status, reason))

    def candidate_array(field, dtype=float):
        if field not in candidate:
            add_failure('candidate_missing_required_field_rejected', f'missing_{field}')
            return np.array([], dtype=dtype)
        try:
            values = np.asarray(candidate[field], dtype=dtype)
        except (TypeError, ValueError):
            add_failure('candidate_missing_required_field_rejected', f'invalid_{field}')
            return np.array([], dtype=dtype)
        if values.size == 0:
            add_failure('candidate_missing_required_field_rejected', f'empty_{field}')
        return values

    Pvals_MPa = candidate_array('candidateP_MPa')
    Tvals_K = candidate_array('candidateT_K')
    phaseVals = candidate_array('candidatePhaseArray', dtype=int)
    etaVals_Pas = candidate_array('candidateEta_Pas')

    if Pvals_MPa.size and np.any(~np.isfinite(Pvals_MPa)):
        add_failure('candidate_missing_required_field_rejected', 'nonfinite_candidateP_MPa')
    if Tvals_K.size and np.any(~np.isfinite(Tvals_K)):
        add_failure('candidate_missing_required_field_rejected', 'nonfinite_candidateT_K')
    if phaseVals.size and np.any(phaseVals != 6):
        add_failure('candidate_phase_boundary_rejected', 'candidate_contains_non_ice_vi_nodes')
    if etaVals_Pas.size and (np.any(~np.isfinite(etaVals_Pas)) or np.any(etaVals_Pas <= 0)):
        add_failure('candidate_invalid_viscosity_rejected', 'candidate_profile_viscosity_invalid')

    temperatureContrast_K = np.nan
    if Tvals_K.size >= 2 and np.all(np.isfinite(Tvals_K[[0, -1]])):
        temperatureContrast_K = float(Tvals_K[-1] - Tvals_K[0])
    if not _ActiveIceVIFinitePositive(temperatureContrast_K):
        add_failure('candidate_invalid_contrast_rejected', 'candidate_temperature_contrast_not_positive')

    energyResidual_W = iceVIDiagnostics.get('energyResidual_W', np.nan)
    energyResidual_frac = iceVIDiagnostics.get('energyResidual_frac', np.nan)
    Q_in_W = candidate.get('candidateQ_in_W', iceVIDiagnostics.get('Q_in_W', np.nan))
    Q_out_W = candidate.get('candidateQ_out_W', iceVIDiagnostics.get('Q_out_W', np.nan))
    if not (_ActiveIceVIFiniteScalar(energyResidual_W) and _ActiveIceVIFiniteScalar(energyResidual_frac)):
        add_failure('candidate_missing_required_field_rejected', 'missing_energy_residual')
    if not (_ActiveIceVIFiniteScalar(Q_in_W) and _ActiveIceVIFiniteScalar(Q_out_W)):
        add_failure('candidate_missing_required_field_rejected', 'missing_heat_power_scale')
        Qscale_W = 1.0
    else:
        Qscale_W = max(abs(Q_in_W), abs(Q_out_W), 1.0)
    energyAbsTol_W = max(
        HP_ICE_PRODUCTION_ENERGY_FRAC_TOL * Qscale_W,
        HP_ICE_PRODUCTION_ENERGY_ABS_FLOOR_W,
    )
    if (
        _ActiveIceVIFiniteScalar(energyResidual_W) and
        _ActiveIceVIFiniteScalar(energyResidual_frac) and
        (abs(energyResidual_W) > energyAbsTol_W or abs(energyResidual_frac) > HP_ICE_PRODUCTION_ENERGY_FRAC_TOL)
    ):
        add_failure('candidate_energy_residual_rejected', 'energy_residual_exceeds_tolerance')
    if _ActiveIceVIFiniteScalar(Q_in_W) and _ActiveIceVIFiniteScalar(Q_out_W):
        heatPowerResidual_W = float(Q_out_W - Q_in_W)
        if abs(heatPowerResidual_W) > energyAbsTol_W:
            add_failure('candidate_heat_flux_residual_rejected', 'heat_power_residual_exceeds_tolerance')
    else:
        heatPowerResidual_W = np.nan

    q_in_Wm2 = candidate.get('candidateq_in_Wm2', np.nan)
    q_out_Wm2 = candidate.get('candidateq_out_Wm2', np.nan)
    expected_q_in_Wm2 = np.nan
    expected_q_out_Wm2 = np.nan
    sphericalFluxScalingStatus = 'not_evaluated'
    if not (_ActiveIceVIFiniteScalar(q_in_Wm2) and _ActiveIceVIFiniteScalar(q_out_Wm2)):
        heatFluxResidual_Wm2 = np.nan
        add_failure('candidate_missing_required_field_rejected', 'missing_heat_flux_residual')
    else:
        heatFluxResidual_Wm2 = float(q_out_Wm2 - q_in_Wm2)
        rVals_m = np.asarray(candidate.get('candidateR_m', np.array([], dtype=float)), dtype=float)
        finiteRadii = rVals_m.size >= 2 and np.all(np.isfinite(rVals_m)) and np.all(rVals_m > 0)
        if finiteRadii and _ActiveIceVIFiniteScalar(Q_in_W) and _ActiveIceVIFiniteScalar(Q_out_W):
            rTop_m = float(rVals_m[0])
            rBot_m = float(rVals_m[-1])
            expected_q_in_Wm2 = float(Q_in_W / (4.0 * np.pi * rBot_m**2))
            expected_q_out_Wm2 = float(Q_out_W / (4.0 * np.pi * rTop_m**2))
            qscale_Wm2 = max(
                abs(q_in_Wm2), abs(q_out_Wm2),
                abs(expected_q_in_Wm2), abs(expected_q_out_Wm2), 1.0,
            )
            heatFluxTol_Wm2 = max(
                ACTIVE_ICE_VI_HEAT_FLUX_FRAC_TOL * qscale_Wm2,
                ACTIVE_ICE_VI_HEAT_FLUX_ABS_FLOOR_WM2,
            )
            qInMatchesArea = abs(q_in_Wm2 - expected_q_in_Wm2) <= heatFluxTol_Wm2
            qOutMatchesArea = abs(q_out_Wm2 - expected_q_out_Wm2) <= heatFluxTol_Wm2
            if qOutMatchesArea:
                if qInMatchesArea:
                    sphericalFluxScalingStatus = 'spherical_area_scaled'
                else:
                    sphericalFluxScalingStatus = 'spherical_area_scaled_input_boundary_mismatch'
            else:
                sphericalFluxScalingStatus = 'spherical_flux_scaling_mismatch'
                add_failure('candidate_heat_flux_residual_rejected', 'spherical_flux_scaling_mismatch')
        else:
            qscale_Wm2 = max(abs(q_in_Wm2), abs(q_out_Wm2), 1.0)
            sphericalFluxScalingStatus = 'spherical_flux_scaling_unavailable'
            add_failure('candidate_heat_flux_residual_rejected', 'spherical_flux_scaling_unavailable')
        heatFluxTol_Wm2 = max(
            ACTIVE_ICE_VI_HEAT_FLUX_FRAC_TOL * qscale_Wm2,
            ACTIVE_ICE_VI_HEAT_FLUX_ABS_FLOOR_WM2,
        )

    massResidual_kg = iceVIDiagnostics.get('massResidual_kg', np.nan)
    massResidual_frac = iceVIDiagnostics.get('massResidual_frac', np.nan)
    if not (_ActiveIceVIFiniteScalar(massResidual_kg) and _ActiveIceVIFiniteScalar(massResidual_frac)):
        add_failure('candidate_missing_required_field_rejected', 'missing_mass_residual')
    elif massResidual_kg != 0.0 or massResidual_frac != 0.0:
        add_failure('candidate_mass_residual_rejected', 'mass_residual_nonzero')

    posthoc = iceVIDiagnostics.get('posthocProductionCandidate', {})
    phaseBoundaryResidual_K = posthoc.get(
        'posthocPhaseBoundaryResidual_K',
        iceVIDiagnostics.get('phaseBoundaryResidual_K', np.nan),
    )
    minPhaseBoundaryMargin_K = posthoc.get('posthocMinPhaseBoundaryMargin_K', np.nan)
    eosPressureMargin_MPa = posthoc.get('posthocEOSPressureMargin_MPa', np.nan)
    eosTemperatureMargin_K = posthoc.get('posthocEOSTemperatureMargin_K', np.nan)
    if not _ActiveIceVIFiniteScalar(phaseBoundaryResidual_K):
        add_failure('candidate_missing_required_field_rejected', 'missing_phase_boundary_residual')
    elif abs(phaseBoundaryResidual_K) > HP_ICE_PRODUCTION_PHASE_BOUNDARY_TOL_K:
        add_failure('candidate_phase_boundary_rejected', 'phase_boundary_residual_exceeds_tolerance')
    if not _ActiveIceVIFinitePositive(minPhaseBoundaryMargin_K):
        add_failure('candidate_phase_boundary_rejected', 'phase_boundary_margin_nonpositive')
    if not (_ActiveIceVIFinitePositive(eosPressureMargin_MPa) and _ActiveIceVIFinitePositive(eosTemperatureMargin_K)):
        add_failure('candidate_outside_eos_domain_rejected', 'eos_margin_nonpositive')

    zVals_m = np.asarray(candidate.get('candidateZ_m', np.array([], dtype=float)), dtype=float)
    if zVals_m.size >= 2 and np.all(np.isfinite(zVals_m)):
        thickness_m = float(np.nanmax(zVals_m) - np.nanmin(zVals_m))
    else:
        thickness_m = iceVIDiagnostics.get('thickness_m', np.nan)
    layerClosureResidual_m = iceVIDiagnostics.get('layerClosureResidual_m', np.nan)
    if not _ActiveIceVIFinitePositive(thickness_m):
        add_failure('candidate_missing_required_field_rejected', 'missing_candidate_thickness')
        closureTol_m = HP_ICE_PRODUCTION_LAYER_CLOSURE_ABS_TOL_M
    else:
        closureTol_m = max(
            HP_ICE_PRODUCTION_LAYER_CLOSURE_FRAC_TOL * thickness_m,
            HP_ICE_PRODUCTION_LAYER_CLOSURE_ABS_TOL_M,
        )
    if not _ActiveIceVIFiniteScalar(layerClosureResidual_m):
        add_failure('candidate_missing_required_field_rejected', 'missing_layer_closure_residual')
    elif abs(layerClosureResidual_m) > closureTol_m:
        add_failure('candidate_layer_closure_rejected', 'layer_closure_residual_exceeds_tolerance')

    RaConvect = iceVIDiagnostics.get('RaConvect', np.nan)
    RaCrit = iceVIDiagnostics.get('RaCrit', np.nan)
    if _ActiveIceVIFinitePositive(RaConvect) and _ActiveIceVIFinitePositive(RaCrit):
        raOverRaCrit = float(RaConvect / RaCrit)
        if raOverRaCrit <= 1.0:
            add_failure('candidate_subcritical_rejected', 'rayleigh_ratio_not_supercritical')
    else:
        raOverRaCrit = np.nan
        add_failure('candidate_missing_required_field_rejected', 'missing_rayleigh_ratio')

    etaConv_Pas = iceVIDiagnostics.get('etaConv_Pas', np.nan)
    etaMelt_Pas = iceVIDiagnostics.get('etaMelt_Pas', np.nan)
    if not (_ActiveIceVIFinitePositive(etaConv_Pas) and _ActiveIceVIFinitePositive(etaMelt_Pas)):
        add_failure('candidate_invalid_viscosity_rejected', 'phase_viscosity_invalid')

    candidate.update({
        'candidateHeatPowerResidual_W': heatPowerResidual_W,
        'candidateEnergyResidual_W': energyResidual_W,
        'candidateEnergyResidual_frac': energyResidual_frac,
        'candidateHeatFluxResidual_Wm2': heatFluxResidual_Wm2,
        'candidateExpected_q_in_Wm2': expected_q_in_Wm2,
        'candidateExpected_q_out_Wm2': expected_q_out_Wm2,
        'candidateSphericalFluxScalingStatus': sphericalFluxScalingStatus,
        'candidateMassResidual_kg': massResidual_kg,
        'candidateMassResidual_frac': massResidual_frac,
        'candidatePhaseBoundaryResidual_K': phaseBoundaryResidual_K,
        'candidateLayerClosureResidual_m': layerClosureResidual_m,
        'candidateEOSPressureMargin_MPa': eosPressureMargin_MPa,
        'candidateEOSTemperatureMargin_K': eosTemperatureMargin_K,
        'candidateMinPhaseBoundaryMargin_K': minPhaseBoundaryMargin_K,
        'candidateRaOverRaCrit': raOverRaCrit,
        'candidateTemperatureContrast_K': temperatureContrast_K,
        'protectedFieldsUnchanged': _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore),
    })

    if failures:
        status = failures[0][0]
        reasons = tuple(reason for _, reason in failures)
        return _ActiveIceVIResidualStatus(candidate, False, status, reasons)
    return _ActiveIceVIResidualStatus(candidate, True, 'candidate_residuals_passed', tuple())


def ApplyActiveIceVIRollbackPolicy(Planet, productionMode=None):
    """Map active Ice VI candidate residual state into metadata-only discard status."""
    if productionMode is None:
        productionMode = ResolveHPIceConvectionModel(Planet) if hasattr(Planet, 'Do') else "Kalousova2018_production_experimental"

    diagnostics = getattr(Planet, 'HPIceDiagnostics', None)
    if not isinstance(diagnostics, dict):
        result = _ActiveIceVIRejectedResult(
            'missing_diagnostics_rejected',
            'missing_hp_ice_diagnostics',
            False,
        )
        _ActiveIceVIRollbackStatus(
            result,
            status='rollback_missing_candidate',
            required=True,
            discarded=True,
            applied=False,
            reason='missing_candidate_copy',
            reasons=('missing_candidate_copy',),
        )
        return _StoreActiveIceVIProductionCandidateCopy(Planet, result)

    iceVIDiagnostics = diagnostics.get('VI')
    if not isinstance(iceVIDiagnostics, dict):
        result = _ActiveIceVIRejectedResult(
            'missing_ice_vi_rejected',
            'missing_ice_vi_diagnostics',
            False,
        )
        _ActiveIceVIRollbackStatus(
            result,
            status='rollback_missing_candidate',
            required=True,
            discarded=True,
            applied=False,
            reason='missing_candidate_copy',
            reasons=('missing_candidate_copy',),
        )
        return _StoreActiveIceVIProductionCandidateCopy(Planet, result)

    candidate = iceVIDiagnostics.get('activeProductionCandidate')
    if not isinstance(candidate, dict):
        candidate = _ActiveIceVIRejectedResult(
            'missing_candidate_copy_rejected',
            'requires_active_candidate_copy',
            False,
        )
        _StoreActiveIceVIProductionCandidateCopy(Planet, candidate)
        return _ActiveIceVIRollbackStatus(
            candidate,
            status='rollback_missing_candidate',
            required=True,
            discarded=True,
            applied=False,
            reason='missing_candidate_copy',
            reasons=('missing_candidate_copy',),
        )

    if candidate.get('candidateDiscarded', False):
        return _ActiveIceVIPreserveTerminalDiscard(candidate)

    if productionMode != "Kalousova2018_production_experimental":
        return _ActiveIceVIRollbackStatus(
            candidate,
            status='rollback_not_evaluated',
            required=False,
            discarded=False,
            applied=False,
            reason='experimental_production_selector_not_enabled',
            reasons=('experimental_production_selector_not_enabled',),
        )

    if not candidate.get('protectedFieldsUnchanged', True):
        return _ActiveIceVIRollbackStatus(
            candidate,
            status='rollback_protected_fields_changed',
            required=True,
            discarded=True,
            applied=False,
            reason='protected_fields_changed',
            reasons=('protected_fields_changed',),
        )

    if not candidate.get('candidateCopyCreated', False):
        reason = candidate.get('candidateReason', 'missing_candidate_copy')
        if candidate.get('candidateStatus') == 'high_risk_posthoc_rejected':
            return _ActiveIceVIRollbackStatus(
                candidate,
                status='rollback_candidate_discarded',
                required=True,
                discarded=True,
                applied=True,
                reason='high_risk_posthoc_candidate',
                reasons=('high_risk_posthoc_candidate', reason),
            )
        return _ActiveIceVIRollbackStatus(
            candidate,
            status='rollback_missing_candidate',
            required=True,
            discarded=True,
            applied=False,
            reason='missing_candidate_copy',
            reasons=('missing_candidate_copy', reason),
        )

    residualStatus = candidate.get('candidateResidualStatus', 'candidate_not_evaluated')
    residualReasons = tuple(candidate.get('candidateResidualReasons', ()))
    if candidate.get('candidateResidualsPassed', False) and residualStatus == 'candidate_residuals_passed':
        return _ActiveIceVIRollbackStatus(
            candidate,
            status='rollback_not_required',
            required=False,
            discarded=False,
            applied=False,
            reason=None,
            reasons=(),
        )

    reason = residualStatus or 'candidate_residuals_not_passed'
    reasons = residualReasons or (reason,)
    return _ActiveIceVIRollbackStatus(
        candidate,
        status='rollback_candidate_discarded',
        required=True,
        discarded=True,
        applied=True,
        reason=reason,
        reasons=reasons,
    )


def _ActiveIceVIThermalUpdateStatus(candidate, *, status, reasons=(),
                                    appliedToCopy=False, extra=None):
    thermalUpdate = {
        'candidateThermalUpdateStrategy': ACTIVE_ICE_VI_THERMAL_UPDATE_STRATEGY_NO_OP,
        'candidateThermalUpdateStatus': status,
        'candidateThermalUpdateReasons': tuple(reasons),
        'candidateThermalUpdateAppliedToCopy': bool(appliedToCopy),
        'candidateThermalUpdateAppliedToPlanet': False,
        'candidateThermalUpdateAccepted': False,
        'candidateUpdatedT_K': np.array([], dtype=float),
        'candidateUpdatedQ_W': np.array([], dtype=float),
        'candidateUpdatedq_Wm2': np.array([], dtype=float),
        'candidateUpdatedPhaseArray': np.array([], dtype=int),
        'candidateThermalEnergyResidual_W': np.nan,
        'candidateThermalEnergyResidual_frac': np.nan,
        'candidateThermalHeatPowerResidual_W': np.nan,
        'candidateThermalPhaseBoundaryResidual_K': np.nan,
        'candidateThermalEOSPressureMargin_MPa': np.nan,
        'candidateThermalEOSTemperatureMargin_K': np.nan,
        'candidateThermalRiskStatus': 'not_evaluated',
        'candidateThermalRiskReasons': tuple(reasons),
    }
    if extra:
        thermalUpdate.update(extra)
    candidate['thermalUpdate'] = thermalUpdate
    candidate['candidateAccepted'] = False
    candidate['candidateAppliedToProfile'] = False
    return thermalUpdate


def BuildActiveIceVIThermalUpdateCandidate(Planet, productionMode=None):
    """Build a no-op thermal-update reconstruction on the isolated Ice VI copy."""
    protectedBefore = _ActiveIceVIProtectedProfileSnapshot(Planet)
    if productionMode is None:
        productionMode = ResolveHPIceConvectionModel(Planet) if hasattr(Planet, 'Do') else "Kalousova2018_production_experimental"

    diagnostics = getattr(Planet, 'HPIceDiagnostics', None)
    if not isinstance(diagnostics, dict):
        candidate = _ActiveIceVIRejectedResult(
            'missing_diagnostics_rejected',
            'missing_hp_ice_diagnostics',
            _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore),
        )
        _StoreActiveIceVIProductionCandidateCopy(Planet, candidate)
        return _ActiveIceVIThermalUpdateStatus(
            candidate,
            status='candidate_thermal_update_missing_candidate_rejected',
            reasons=('missing_hp_ice_diagnostics',),
        )

    iceVIDiagnostics = diagnostics.get('VI')
    if not isinstance(iceVIDiagnostics, dict):
        candidate = _ActiveIceVIRejectedResult(
            'missing_ice_vi_rejected',
            'missing_ice_vi_diagnostics',
            _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore),
        )
        _StoreActiveIceVIProductionCandidateCopy(Planet, candidate)
        return _ActiveIceVIThermalUpdateStatus(
            candidate,
            status='candidate_thermal_update_missing_candidate_rejected',
            reasons=('missing_ice_vi_diagnostics',),
        )

    candidate = iceVIDiagnostics.get('activeProductionCandidate')
    if not isinstance(candidate, dict):
        candidate = _ActiveIceVIRejectedResult(
            'missing_candidate_copy_rejected',
            'requires_active_candidate_copy',
            _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore),
        )
        _StoreActiveIceVIProductionCandidateCopy(Planet, candidate)
        return _ActiveIceVIThermalUpdateStatus(
            candidate,
            status='candidate_thermal_update_missing_candidate_rejected',
            reasons=('requires_active_candidate_copy',),
        )

    candidate['candidateAccepted'] = False
    candidate['candidateAppliedToProfile'] = False
    candidate['protectedFieldsUnchanged'] = _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore)

    if productionMode != "Kalousova2018_production_experimental":
        return _ActiveIceVIThermalUpdateStatus(
            candidate,
            status='candidate_thermal_update_not_evaluated',
            reasons=('experimental_production_selector_not_enabled',),
        )

    if candidate.get('candidateDiscarded', False):
        _ActiveIceVIPreserveTerminalDiscard(candidate)
        return _ActiveIceVIThermalUpdateStatus(
            candidate,
            status='candidate_thermal_update_discarded_candidate_rejected',
            reasons=(candidate.get('candidateDiscardReason') or 'candidate_discarded',),
        )

    if not candidate.get('candidateCopyCreated', False):
        return _ActiveIceVIThermalUpdateStatus(
            candidate,
            status='candidate_thermal_update_missing_candidate_rejected',
            reasons=(candidate.get('candidateReason', 'requires_active_candidate_copy'),),
        )

    if not (
        candidate.get('candidateResidualsPassed', False) and
        candidate.get('candidateResidualStatus') == 'candidate_residuals_passed'
    ):
        return _ActiveIceVIThermalUpdateStatus(
            candidate,
            status='candidate_thermal_update_residuals_not_passed_rejected',
            reasons=tuple(candidate.get('candidateResidualReasons', ())) or
                    (candidate.get('candidateResidualStatus', 'candidate_residuals_not_passed'),),
        )

    try:
        Tvals_K = np.asarray(candidate['candidateT_K'], dtype=float)
        phaseVals = np.asarray(candidate['candidatePhaseArray'], dtype=int)
        rVals_m = np.asarray(candidate['candidateR_m'], dtype=float)
    except (KeyError, TypeError, ValueError):
        return _ActiveIceVIThermalUpdateStatus(
            candidate,
            status='candidate_thermal_update_missing_candidate_rejected',
            reasons=('missing_candidate_thermal_arrays',),
        )

    if Tvals_K.size == 0 or phaseVals.size == 0:
        return _ActiveIceVIThermalUpdateStatus(
            candidate,
            status='candidate_thermal_update_missing_candidate_rejected',
            reasons=('empty_candidate_thermal_arrays',),
        )
    if Tvals_K.size != phaseVals.size:
        return _ActiveIceVIThermalUpdateStatus(
            candidate,
            status='candidate_thermal_update_missing_candidate_rejected',
            reasons=('candidate_thermal_array_size_mismatch',),
        )
    if np.any(~np.isfinite(Tvals_K)) or np.any(phaseVals != 6):
        return _ActiveIceVIThermalUpdateStatus(
            candidate,
            status='candidate_thermal_update_nonfinite_rejected',
            reasons=('nonfinite_candidate_temperature_or_non_ice_vi_phase',),
        )
    if rVals_m.size != Tvals_K.size or np.any(~np.isfinite(rVals_m)) or np.any(rVals_m <= 0):
        return _ActiveIceVIThermalUpdateStatus(
            candidate,
            status='candidate_thermal_update_missing_radius_rejected',
            reasons=('candidate_radius_required_for_spherical_flux_reconstruction',),
        )

    Q_in_W = candidate.get('candidateQ_in_W', iceVIDiagnostics.get('Q_in_W', np.nan))
    Q_out_W = candidate.get('candidateQ_out_W', iceVIDiagnostics.get('Q_out_W', np.nan))
    if not (_ActiveIceVIFiniteScalar(Q_in_W) and _ActiveIceVIFiniteScalar(Q_out_W)):
        return _ActiveIceVIThermalUpdateStatus(
            candidate,
            status='candidate_thermal_update_missing_candidate_rejected',
            reasons=('missing_candidate_heat_power',),
        )

    Qscale_W = max(abs(Q_in_W), abs(Q_out_W), 1.0)
    heatPowerResidual_W = float(Q_out_W - Q_in_W)
    energyAbsTol_W = max(
        HP_ICE_PRODUCTION_ENERGY_FRAC_TOL * Qscale_W,
        HP_ICE_PRODUCTION_ENERGY_ABS_FLOOR_W,
    )
    if abs(heatPowerResidual_W) > energyAbsTol_W:
        return _ActiveIceVIThermalUpdateStatus(
            candidate,
            status='candidate_thermal_update_residuals_not_passed_rejected',
            reasons=('candidate_heat_power_residual_exceeds_tolerance',),
        )

    updatedT_K = np.array(Tvals_K, dtype=float, copy=True)
    updatedPhase = np.array(phaseVals, dtype=int, copy=True)
    Qthrough_W = float(0.5 * (Q_in_W + Q_out_W))
    updatedQ_W = np.full(updatedT_K.shape, Qthrough_W, dtype=float)
    updatedq_Wm2 = updatedQ_W / (4.0 * np.pi * rVals_m**2)
    energyResidual_W = candidate.get('candidateEnergyResidual_W', np.nan)
    energyResidual_frac = candidate.get('candidateEnergyResidual_frac', np.nan)
    phaseBoundaryResidual_K = candidate.get('candidatePhaseBoundaryResidual_K', np.nan)
    eosPressureMargin_MPa = candidate.get('candidateEOSPressureMargin_MPa', np.nan)
    eosTemperatureMargin_K = candidate.get('candidateEOSTemperatureMargin_K', np.nan)
    riskReasons = tuple()
    if candidate.get('candidateSphericalFluxScalingStatus') not in (
        'spherical_area_scaled',
        'spherical_area_scaled_input_boundary_mismatch',
    ):
        riskReasons = (candidate.get('candidateSphericalFluxScalingStatus', 'spherical_flux_scaling_not_evaluated'),)
    riskStatus = 'nominal' if not riskReasons else 'requires_review'

    extra = {
        'candidateUpdatedT_K': updatedT_K,
        'candidateUpdatedQ_W': updatedQ_W,
        'candidateUpdatedq_Wm2': updatedq_Wm2,
        'candidateUpdatedPhaseArray': updatedPhase,
        'candidateThermalEnergyResidual_W': energyResidual_W,
        'candidateThermalEnergyResidual_frac': energyResidual_frac,
        'candidateThermalHeatPowerResidual_W': heatPowerResidual_W,
        'candidateThermalPhaseBoundaryResidual_K': phaseBoundaryResidual_K,
        'candidateThermalEOSPressureMargin_MPa': eosPressureMargin_MPa,
        'candidateThermalEOSTemperatureMargin_K': eosTemperatureMargin_K,
        'candidateThermalRiskStatus': riskStatus,
        'candidateThermalRiskReasons': riskReasons,
        'candidateThermalProtectedFieldsUnchanged': _ActiveIceVIProtectedProfileUnchanged(Planet, protectedBefore),
    }
    return _ActiveIceVIThermalUpdateStatus(
        candidate,
        status='candidate_thermal_update_no_op_reconstruction',
        reasons=tuple(),
        appliedToCopy=True,
        extra=extra,
    )


def _SetHPIceDiagnosticFields(Planet, phaseName, status, phaseID=None, iTop=None, iBot=None,
                              rTop_m=np.nan, rBot_m=np.nan, zTop_m=np.nan, zBot_m=np.nan,
                              thickness_m=np.nan, Ttop_K=np.nan, Tbot_K=np.nan,
                              Ptop_MPa=np.nan, Pbot_MPa=np.nan, qBot_Wm2=np.nan,
                              Q_in_W=np.nan, Q_out_W=np.nan, q_in_Wm2=np.nan, q_out_Wm2=np.nan,
                              Q_internal_W=0.0, Q_latent_W=0.0, Tconv_K=np.nan,
                              etaConv_Pas=np.nan, etaMelt_Pas=np.nan, eLid_m=np.nan,
                              Dconv_m=np.nan, deltaTBL_m=np.nan, Ra=np.nan,
                              RaCrit=np.nan, Qbot_W=np.nan, Pmid_MPa=np.nan,
                              method=None, meltFraction=np.nan,
                              rho_kgm3=np.nan, Cp_JkgK=np.nan, alpha_pK=np.nan,
                              kTherm_WmK=np.nan,
                              productionMode=None, productionCandidate=False,
                              updateAccepted=False, candidateStatus="not_evaluated",
                              candidateReason=None, massResidual_kg=np.nan,
                              massResidual_frac=np.nan,
                              phaseBoundaryResidual_K=np.nan, Tmelt_top_K=np.nan,
                              Tmelt_mid_K=np.nan, Tmelt_bot_K=np.nan,
                              TmeltSource=None, meltCurveStatus="not_evaluated",
                              phaseBoundaryStatus="not_evaluated",
                              eosBoundaryContext="not_evaluated",
                              eosBoundaryStatus="not_evaluated",
                              eosBoundaryReason=None,
                              candidateBoundarySource=None,
                              finalProfileCoverageStatus="final_profile_not_evaluated",
                              viscositySource=None):
    """Single writer for top-level HP diagnostic fields and HPIceDiagnostics."""
    setattr(Planet, f'Tconv{phaseName}_K', Tconv_K)
    setattr(Planet, f'etaConv{phaseName}_Pas', etaConv_Pas)
    setattr(Planet, f'etaMelt{phaseName}_Pas', etaMelt_Pas)
    setattr(Planet, f'eLid{phaseName}_m', eLid_m)
    setattr(Planet, f'Dconv{phaseName}_m', Dconv_m)
    setattr(Planet, f'deltaTBL{phaseName}_m', deltaTBL_m)
    setattr(Planet, f'RaConvect{phaseName}', Ra)
    setattr(Planet, f'RaCrit{phaseName}', RaCrit)
    setattr(Planet, f'meltFraction{phaseName}', meltFraction)
    phaseState = HPIcePhaseState(
        phaseName=phaseName, phaseID=phaseID, present=status != 'absent',
        status=status, iTop=iTop, iBot=iBot, rTop_m=rTop_m, rBot_m=rBot_m,
        zTop_m=zTop_m, zBot_m=zBot_m, thickness_m=thickness_m,
        Ttop_K=Ttop_K, Tbot_K=Tbot_K, Ptop_MPa=Ptop_MPa, Pbot_MPa=Pbot_MPa,
        qBot_Wm2=qBot_Wm2, Qbot_W=Qbot_W, Q_in_W=Q_in_W, Q_out_W=Q_out_W,
        q_in_Wm2=q_in_Wm2, q_out_Wm2=q_out_Wm2,
        Q_internal_W=Q_internal_W, Q_latent_W=Q_latent_W,
        rho_kgm3=rho_kgm3, Cp_JkgK=Cp_JkgK, alpha_pK=alpha_pK,
        kTherm_WmK=kTherm_WmK,
        productionMode=productionMode, productionCandidate=productionCandidate,
        updateAccepted=updateAccepted, candidateStatus=candidateStatus,
        candidateReason=candidateReason, massResidual_kg=massResidual_kg,
        massResidual_frac=massResidual_frac,
        phaseBoundaryResidual_K=phaseBoundaryResidual_K,
        Tmelt_top_K=Tmelt_top_K, Tmelt_mid_K=Tmelt_mid_K,
        Tmelt_bot_K=Tmelt_bot_K, TmeltSource=TmeltSource,
        meltCurveStatus=meltCurveStatus,
        phaseBoundaryStatus=phaseBoundaryStatus,
        eosBoundaryContext=eosBoundaryContext,
        eosBoundaryStatus=eosBoundaryStatus,
        eosBoundaryReason=eosBoundaryReason,
        candidateBoundarySource=candidateBoundarySource,
        finalProfileCoverageStatus=finalProfileCoverageStatus,
        viscositySource=viscositySource,
        Tconv_K=Tconv_K,
        etaConv_Pas=etaConv_Pas, etaMelt_Pas=etaMelt_Pas,
        eLid_m=eLid_m, Dconv_m=Dconv_m, deltaTBL_m=deltaTBL_m,
        RaConvect=Ra, RaCrit=RaCrit, Pmid_MPa=Pmid_MPa, method=method,
        meltFraction=meltFraction,
        DO_HP_MELT=bool(np.isfinite(meltFraction) and meltFraction > 0),
    )
    Planet.HPIceDiagnostics[phaseName] = phaseState.as_dict()


def GetOceanHPIceEOS(Planet, Params, POcean_MPa, minPres_MPa=None, minTres_K=None):
    """ Assign EOS functions for possible high-pressure ices expected in the ocean
        based on the min/max temperatures and pressures we plan to model.

        Args:
            POcean_MPa (float, shape N): Range of pressures expected in the ocean, from
                Pb_MPa to PHydroMax_MPa
        Assigns Planet attributes:
            Ocean.iceEOS
    """

    # Generate linear sets of outside possibilities of P,T combinations for the ocean
    POceanHPices_MPa = POcean_MPa[POcean_MPa>Constants.PminHPices_MPa]
    TOceanHPices_K = np.arange(Planet.Bulk.Tb_K, Planet.Ocean.THydroMax_K, Planet.Ocean.deltaT)
    PHPicesLin_MPa = np.array([P for P in POceanHPices_MPa for _ in TOceanHPices_K])
    THPicesLin_K = np.array([T for _ in POceanHPices_MPa for T in TOceanHPices_K])

    # Stopgap measure to avoid MgSO4 calcs taking ages with the slow Margules formulation phase calcs
    # Remove this if/else block (just do the "else") when a faster phase calculation is implemented!
    if (Planet.Ocean.comp == 'MgSO4' or Planet.Sil.poreComp == 'MgSO4') and Planet.Ocean.phaseType == 'calc':
        # Just load all HP ice phases in case we need them. This part is way faster than Margules phase calcs
        Planet.Ocean.iceEOS['II'] = GetIceEOS(POceanHPices_MPa, TOceanHPices_K, 'II',
                                              porosType=Planet.Ocean.porosType['II'],
                                              phiTop_frac=Planet.Ocean.phiMax_frac['II'],
                                              Pclosure_MPa=Planet.Ocean.Pclosure_MPa['II'],
                                              phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE['II'],
                                              minPres_MPa=minPres_MPa, minTres_K=minTres_K, kThermConst_WmK=Planet.Ocean.kThermIce_WmK)
        Planet.Ocean.iceEOS['III'] = GetIceEOS(POceanHPices_MPa, TOceanHPices_K, 'III',
                                               porosType=Planet.Ocean.porosType['III'],
                                               phiTop_frac=Planet.Ocean.phiMax_frac['III'],
                                               Pclosure_MPa=Planet.Ocean.Pclosure_MPa['III'],
                                               phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE['III'],
                                               minPres_MPa=minPres_MPa, minTres_K=minTres_K, kThermConst_WmK=Planet.Ocean.kThermIce_WmK,
                                               **GetIceArrheniusViscosityKwargs(Planet, 'III'))
        Planet.Ocean.iceEOS['V'] = GetIceEOS(POceanHPices_MPa, TOceanHPices_K, 'V',
                                             porosType=Planet.Ocean.porosType['V'],
                                             phiTop_frac=Planet.Ocean.phiMax_frac['V'],
                                             Pclosure_MPa=Planet.Ocean.Pclosure_MPa['V'],
                                             phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE['V'],
                                             minPres_MPa=minPres_MPa, minTres_K=minTres_K, kThermConst_WmK=Planet.Ocean.kThermIce_WmK,
                                             **GetIceArrheniusViscosityKwargs(Planet, 'V'))
        Planet.Ocean.iceEOS['VI'] = GetIceEOS(POceanHPices_MPa, TOceanHPices_K, 'VI',
                                              porosType=Planet.Ocean.porosType['VI'],
                                              phiTop_frac=Planet.Ocean.phiMax_frac['VI'],
                                              Pclosure_MPa=Planet.Ocean.Pclosure_MPa['VI'],
                                              phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE['VI'],
                                              minPres_MPa=minPres_MPa, minTres_K=minTres_K, kThermConst_WmK=Planet.Ocean.kThermIce_WmK,
                                              **GetIceArrheniusViscosityKwargs(Planet, 'VI'))
    else:
        # Get phase of each P,T combination
        expandPhases = Planet.Ocean.EOS.fn_phase(POceanHPices_MPa, TOceanHPices_K, grid = True).flatten()
        # Check if any of them are not liquid
        #TODO Expand this implementation for underplate clathrates mixed with high pressure ices
        if np.any(expandPhases != 0):
            _, _, _, indsIceII, _, indsIceIII, _, indsIceV, _, indsIceVI, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _ \
                = GetPhaseIndices(expandPhases)

            if(np.size(indsIceII) != 0):
                log.debug('Loading ice II EOS functions for ocean layers...')
                PiceIImin_MPa = PHPicesLin_MPa[indsIceII[0]]
                PiceIImax_MPa = PHPicesLin_MPa[indsIceII[-1]]
                TiceIImin_K = np.min(THPicesLin_K[indsIceII])
                TiceIImax_K = np.max(THPicesLin_K[indsIceII])
                Planet.Ocean.iceEOS['II'] = GetIceEOS(np.linspace(PiceIImin_MPa, PiceIImax_MPa, Planet.Steps.nPsHP),
                                                      np.linspace(TiceIImin_K, TiceIImax_K, Planet.Steps.nTsHP), 'II',
                                                      porosType=Planet.Ocean.porosType['II'],
                                                      phiTop_frac=Planet.Ocean.phiMax_frac['II'],
                                                      Pclosure_MPa=Planet.Ocean.Pclosure_MPa['II'],
                                                      phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE['II'],
                                                      minPres_MPa=minPres_MPa, minTres_K=minTres_K, kThermConst_WmK=Planet.Ocean.kThermIce_WmK)
            if(np.size(indsIceIII) != 0):
                log.debug('Loading ice III EOS functions for ocean layers...')
                PiceIIImin_MPa = PHPicesLin_MPa[indsIceIII[0]]
                PiceIIImax_MPa = PHPicesLin_MPa[indsIceIII[-1]]
                TiceIIImin_K = np.min(THPicesLin_K[indsIceIII])
                TiceIIImax_K = np.max(THPicesLin_K[indsIceIII])
                Planet.Ocean.iceEOS['III'] = GetIceEOS(np.linspace(PiceIIImin_MPa, PiceIIImax_MPa, Planet.Steps.nPsHP),
                                                       np.linspace(TiceIIImin_K, TiceIIImax_K, Planet.Steps.nTsHP), 'III',
                                                       porosType=Planet.Ocean.porosType['III'],
                                                       phiTop_frac=Planet.Ocean.phiMax_frac['III'],
                                                       Pclosure_MPa=Planet.Ocean.Pclosure_MPa['III'],
                                                       phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE['III'],
                                                       minPres_MPa=minPres_MPa, minTres_K=minTres_K, kThermConst_WmK=Planet.Ocean.kThermIce_WmK,
                                                       **GetIceArrheniusViscosityKwargs(Planet, 'III'))
            if(np.size(indsIceV) != 0):
                log.debug('Loading ice V EOS functions for ocean layers...')
                PiceVmin_MPa = PHPicesLin_MPa[indsIceV[0]]
                PiceVmax_MPa = PHPicesLin_MPa[indsIceV[-1]]
                TiceVmin_K = np.min(THPicesLin_K[indsIceV])
                TiceVmax_K = np.max(THPicesLin_K[indsIceV])
                Planet.Ocean.iceEOS['V'] = GetIceEOS(np.linspace(PiceVmin_MPa, PiceVmax_MPa, Planet.Steps.nPsHP),
                                                     np.linspace(TiceVmin_K, TiceVmax_K, Planet.Steps.nTsHP), 'V',
                                                     porosType=Planet.Ocean.porosType['V'],
                                                     phiTop_frac=Planet.Ocean.phiMax_frac['V'],
                                                     Pclosure_MPa=Planet.Ocean.Pclosure_MPa['V'],
                                                     phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE['V'],
                                                     minPres_MPa=minPres_MPa, minTres_K=minTres_K, kThermConst_WmK=Planet.Ocean.kThermIce_WmK,
                                                     **GetIceArrheniusViscosityKwargs(Planet, 'V'))
            if(np.size(indsIceVI) != 0):
                log.debug('Loading ice VI EOS functions for ocean layers...')
                PiceVImin_MPa = PHPicesLin_MPa[indsIceVI[0]]
                PiceVImax_MPa = PHPicesLin_MPa[indsIceVI[-1]]
                TiceVImin_K = np.min(THPicesLin_K[indsIceVI])
                TiceVImax_K = np.max(THPicesLin_K[indsIceVI])
                Planet.Ocean.iceEOS['VI'] = GetIceEOS(np.linspace(PiceVImin_MPa, PiceVImax_MPa, Planet.Steps.nPsHP),
                                                      np.linspace(TiceVImin_K, TiceVImax_K, Planet.Steps.nTsHP), 'VI',
                                                      porosType=Planet.Ocean.porosType['VI'],
                                                      phiTop_frac=Planet.Ocean.phiMax_frac['VI'],   
                                                      Pclosure_MPa=Planet.Ocean.Pclosure_MPa['VI'],
                                                      phiMin_frac=Planet.Ocean.phiMin_frac, EXTRAP=Params.EXTRAP_ICE['VI'],
                                                      minPres_MPa=minPres_MPa, minTres_K=minTres_K, kThermConst_WmK=Planet.Ocean.kThermIce_WmK,
                                                      **GetIceArrheniusViscosityKwargs(Planet, 'VI'))
        else:
            log.debug('No high-pressure ices found in ocean layers.')

    return Planet, Params


def InnerLayers(Planet, Params):
    """ Wrapper function for inner layer propogation. Decides between self consistent and non-self consistent modeling.
    """
    Timing.setFunctionTime(time.time())
    # Early branching for non-self-consistent modeling
    if Planet.Do.NON_SELF_CONSISTENT:
        Planet, Params = NonSelfConsistentInnerLayer(Planet, Params)
    else:
        Planet, Params = SelfConsistentInnerLayer(Planet, Params)
    Timing.printFunctionTimeDifference('InnerLayers()', time.time())
    return Planet

def SelfConsistentInnerLayer(Planet, Params):
    """ Geophysical and thermodynamic calculations for silicate and core layers
        Calculates state variables of the layer with each pressure step

        Assigns Planet attributes:
            Steps.nTotal, all layer arrays
    """
    if Planet.Do.VALID:
        # Make sure the necessary EOSs have been loaded (mainly only important in parallel ExploreOgram runs)
        if not Params.SKIP_INNER and Planet.Sil.EOS.key not in EOSlist.loaded.keys():
            Planet.Sil.EOS = GetInnerEOS(Planet.Sil.mantleEOS, EOSinterpMethod=Params.lookupInterpMethod,
                                         kThermConst_WmK=Planet.Sil.kTherm_WmK, HtidalConst_Wm3=Planet.Sil.Htidal_Wm3,
                                         porosType=Planet.Sil.porosType, phiTop_frac=Planet.Sil.phiRockMax_frac,
                                         Pclosure_MPa=Planet.Sil.Pclosure_MPa, phiMin_frac=Planet.Sil.phiMin_frac,
                                         EXTRAP=Params.EXTRAP_SIL, etaSilFixed_Pas=Planet.Sil.etaRock_Pas, etaCoreFixed_Pas=[Planet.Core.etaFeSolid_Pas, Planet.Core.etaFeLiquid_Pas],
                                         TviscTrans_K=Planet.Sil.TviscTrans_K,
                                         doConstantProps=Planet.Do.CONSTANT_INNER_DENSITY, constantProperties={'rho_kgm3': Planet.Sil.rhoSilWithCore_kgm3, 'Cp_JkgK': np.nan, 'alpha_pK': np.nan, 'kTherm_WmK': Planet.Sil.kTherm_WmK,
                                                                                   'VP_kms': Planet.Sil.VPset_kms, 'VS_kms': Planet.Sil.VSset_kms, 'KS_GPa': Planet.Sil.KSset_GPa, 'GS_GPa': Planet.Sil.GSset_GPa, 'eta_Pas': Planet.Sil.etaRock_Pas,
                                                                                   'sigma_Sm': Planet.Sil.sigmaSil_Sm})
        # Iron core if present
        if not Params.SKIP_INNER and Planet.Do.Fe_CORE and Planet.Core.EOS.key not in EOSlist.loaded.keys():
            Planet.Core.EOS = GetInnerEOS(Planet.Core.coreEOS, EOSinterpMethod=Params.lookupInterpMethod, Fe_EOS=True,
                                          kThermConst_WmK=Planet.Core.kTherm_WmK, EXTRAP=Params.EXTRAP_Fe,
                                          wFeCore_ppt=Planet.Core.wFe_ppt, wScore_ppt=Planet.Core.wS_ppt, etaSilFixed_Pas=Planet.Sil.etaRock_Pas, etaCoreFixed_Pas=[Planet.Core.etaFeSolid_Pas, Planet.Core.etaFeLiquid_Pas],
                                          TviscTrans_K=Planet.Core.TviscTrans_K,
                                          doConstantProps=True, constantProperties={'rho_kgm3': Planet.Core.rhoFe_kgm3, 'Cp_JkgK': np.nan, 'alpha_pK': np.nan, 'kTherm_WmK': Planet.Core.kTherm_WmK,
                                                                                   'VP_kms': np.nan, 'VS_kms': np.nan, 'KS_GPa': np.nan, 'GS_GPa': Planet.Core.GSset_GPa, 'eta_Pas': Planet.Core.etaFeSolid_Pas,
                                                                                   'sigma_Sm': Planet.Core.sigmaCore_Sm})

        # Pore fluids if present
        if not Params.SKIP_INNER and Planet.Do.POROUS_ROCK and Planet.Sil.poreEOS.key not in EOSlist.loaded.keys():
            if Planet.Do.NO_H2O:
                Ppore_MPa, Tpore_K = (np.linspace(0, 1, 10) for _ in range(2))
            else:
                if Planet.Sil.poreComp == 'Seawater' and Planet.Sil.PHydroMax_MPa > 300:
                    log.warning('GSW yields NaN for Cp at pressures above 300 MPa. Reducing PsilMax to this value.')
                    Planet.Sil.PHydroMax_MPa = 300
                Ppore_MPa = np.linspace(Planet.Bulk.Psurf_MPa, Planet.Sil.PHydroMax_MPa, 100)
                Tpore_K = np.linspace(Planet.Bulk.Tb_K, Planet.Sil.THydroMax_K, 140)
            # Get pore fluid EOS
            Planet.Sil.poreEOS = GetOceanEOS(Planet.Sil.poreComp, Planet.Sil.wPore_ppt, Ppore_MPa, Tpore_K,
                                             Planet.Ocean.MgSO4elecType, rhoType=Planet.Ocean.MgSO4rhoType,
                                             scalingType=Planet.Ocean.MgSO4scalingType, FORCE_NEW=Params.FORCE_EOS_RECALC,
                                             phaseType=Planet.Ocean.phaseType, EXTRAP=Params.EXTRAP_OCEAN, PORE=True,
                                             sigmaFixed_Sm=Planet.Sil.sigmaPoreFixed_Sm, kThermConst_WmK=Planet.Ocean.kThermWater_WmK,
                                             propsStepReductionFactor=Planet.Ocean.propsStepReductionFactor)


            # Make sure Sil.phiRockMax_frac is set in case we're using a porosType that doesn't require it
            if Planet.Sil.phiRockMax_frac is None or Planet.Sil.porosType != 'Han2014':
                Planet.Sil.phiRockMax_frac = Planet.Sil.EOS.fn_phi_frac(0, 0)

        if Planet.Do.CONSTANT_INNER_DENSITY or Params.SKIP_INNER:
            Planet, mantleProps, coreProps = CalcMoIConstantRho(Planet, Params)
        else:
            Planet, mantleProps, coreProps = CalcMoIWithEOS(Planet, Params)

        if Planet.Steps.nHydro <= Planet.Steps.nSurfIce and not Planet.Do.NO_H2O:
            log.warning('For these run settings, the hydrosphere is entirely frozen and contains only surface ice.')
            Planet.Do.NO_OCEAN = True
        Planet.Steps.nTotal = Planet.Steps.nHydro + Planet.Steps.nSil + Planet.Steps.nCore

        log.debug('Evaluating remaining quantities for layer arrays...')
        # Extend Planet layer arrays to make space for silicate and possible core layers
        extend = np.zeros(Planet.Steps.nSil + Planet.Steps.nCore)
        Planet.P_MPa = np.concatenate((Planet.P_MPa[:Planet.Steps.nHydro], extend))
        Planet.T_K = np.concatenate((Planet.T_K[:Planet.Steps.nHydro], extend))
        Planet.r_m = np.concatenate((Planet.r_m[:Planet.Steps.nHydro], extend, [0.0]))
        Planet.phase = np.concatenate((Planet.phase[:Planet.Steps.nHydro], extend.astype(np.int_)))
        Planet.rho_kgm3 = np.concatenate((Planet.rho_kgm3[:Planet.Steps.nHydro], extend))
        Planet.Cp_JkgK = np.concatenate((Planet.Cp_JkgK[:Planet.Steps.nHydro], extend))
        Planet.alpha_pK = np.concatenate((Planet.alpha_pK[:Planet.Steps.nHydro], extend))
        Planet.kTherm_WmK = np.concatenate((Planet.kTherm_WmK[:Planet.Steps.nHydro], extend))
        Planet.g_ms2 = np.concatenate((Planet.g_ms2[:Planet.Steps.nHydro], extend))
        Planet.phi_frac = np.concatenate((Planet.phi_frac[:Planet.Steps.nHydro], extend))
        Planet.Htidal_Wm3 = np.concatenate((Planet.Htidal_Wm3[:Planet.Steps.nHydro], extend))
        Planet.z_m = np.concatenate((Planet.z_m[:Planet.Steps.nHydro], extend, [0.0]))
        Planet.Ppore_MPa = np.concatenate((Planet.Ppore_MPa[:Planet.Steps.nHydro], extend))
        Planet.rhoMatrix_kgm3 = np.concatenate((Planet.rhoMatrix_kgm3[:Planet.Steps.nHydro], extend))
        Planet.rhoPore_kgm3 = np.concatenate((Planet.rhoPore_kgm3[:Planet.Steps.nHydro], extend))
        Planet.MLayer_kg = np.concatenate((Planet.MLayer_kg[:Planet.Steps.nHydro], extend))

        # Unpack results from MoI calculations
        iOS = Planet.Steps.nHydro
        iSC = Planet.Steps.nHydro + Planet.Steps.nSil
        Planet.P_MPa[iOS:iSC], Planet.T_K[iOS:iSC], Planet.r_m[iOS:iSC], Planet.rho_kgm3[iOS:iSC], \
            Planet.g_ms2[iOS:iSC], Planet.phi_frac[iOS:iSC], Planet.Htidal_Wm3[iOS:iSC], Planet.kTherm_WmK[iOS:iSC], \
            Planet.Ppore_MPa[iOS:iSC], Planet.rhoMatrix_kgm3[iOS:iSC], Planet.rhoPore_kgm3[iOS:iSC], \
            Planet.MLayer_kg[iOS:iSC], phasePore \
            = mantleProps

        iCC = Planet.Steps.nTotal
        if Planet.Do.Fe_CORE and iCC > iSC:
            # Unpack results from MoI calculations
            Planet.P_MPa[iSC:iCC], Planet.T_K[iSC:iCC], Planet.r_m[iSC:iCC], Planet.g_ms2[iSC:iCC], Planet.rho_kgm3[iSC:iCC], \
                Planet.Cp_JkgK[iSC:iCC], Planet.alpha_pK[iSC:iCC], Planet.kTherm_WmK[iSC:iCC], Planet.MLayer_kg[iSC:iCC] \
                = coreProps

        Planet.z_m[iOS:iCC+1] = Planet.Bulk.R_m - Planet.r_m[iOS:iCC+1]

        # Assign phase values for silicates and core
        Planet.phase[iOS:iSC] = Constants.phaseSil + phasePore
        Planet.phase[iSC:iCC] = Constants.phaseFe

        # Record ocean layer thickness and properties
        Planet.phiSeafloor_frac = Planet.phi_frac[iOS]
        if Planet.Do.NO_H2O or not np.any(Planet.phase == 0):
            Planet.D_km = 0
            Planet.Pseafloor_MPa = np.nan
            # Reset ice layer thickness in the event there is no ocean
            Planet.zb_km = Planet.z_m[iOS] / 1e3
        else:
            if Planet.Do.VALID:
                # Get first index of phase changing from 0 to something different ---
                # this is the bottom of the contiguous ocean layer.
                iOceanBot = next(i for i,phase in enumerate(Planet.phase[:Planet.Steps.nHydro]) if phase == 0 and phase != Planet.phase[i+1])
                Planet.D_km = Planet.z_m[iOceanBot + 1]/1e3 - Planet.zb_km
                Planet.Pseafloor_MPa = Planet.P_MPa[iOceanBot + 1]
            else:
                Planet.D_km = np.nan
                Planet.Pseafloor_MPa = np.nan

                # Calculate total salt and water masses
        Planet.Mcore_kg = np.sum(Planet.MLayer_kg[iSC:iCC])
        Planet.Mrock_kg = 4/3*np.pi * np.sum((Planet.r_m[iOS:iSC]**3 - Planet.r_m[iOS+1:iSC+1]**3)
                                             * Planet.rhoMatrix_kgm3[iOS:iSC] * (1 - Planet.phi_frac[iOS:iSC]))
        # Next, fetch the phase IDs of the silicate layers, which are incremented when
        # they contain non-liquid phases.
        silPhases = Planet.phase[iOS:iSC]
        # Add matrix density for ice phases; rhoMatrix is not set for phase == 0, so we
        # can safely include those layers in the sum. We also must include ices in the
        # pore space of silicates.
        Planet.Mice_kg = 4/3*np.pi * (np.sum((Planet.r_m[0:iOS]**3 - Planet.r_m[1:iOS+1]**3)
                                             * Planet.rhoMatrix_kgm3[0:iOS] * (1 - Planet.phi_frac[0:iOS]))
                                      + np.sum((Planet.r_m[iOS:iSC][silPhases != Constants.phaseSil]**3
                                                - Planet.r_m[iOS+1:iSC+1][silPhases != Constants.phaseSil]**3)
                                               * Planet.rhoPore_kgm3[iOS:iSC][silPhases != Constants.phaseSil]
                                               * Planet.phi_frac[iOS:iSC][silPhases != Constants.phaseSil]))
        # The remainder is the ocean fluids, including H2O and salts and pore spaces.
        Planet.Mfluid_kg = Planet.Mtot_kg - Planet.Mcore_kg - Planet.Mrock_kg - Planet.Mice_kg
        # Get the mass contained in clathrate layers and just the trapped gas
        clathPhases = np.logical_and(Planet.phase >= Constants.phaseClath, Planet.phase < Constants.phaseClath + 10)
        Planet.Mclath_kg = np.sum(Planet.MLayer_kg[clathPhases])
        Planet.MclathGas_kg = Planet.Mclath_kg * Constants.clathGasFrac_ppt / 1e3
        # The portion just in the ocean is simple to evaluate:
        Planet.Mocean_kg = np.sum(Planet.MLayer_kg[Planet.phase == 0])
        # The difference is the amount contained in the pore space:
        Planet.MporeFluid_kg = Planet.Mfluid_kg - Planet.Mocean_kg
        # Multiply mass concentration of solute to get total mass of salt in ocean
        Planet.MoceanSalt_kg = Planet.Mfluid_kg * Planet.Ocean.wOcean_ppt / 1e3
        # Record the mass of salt in the pore space in case we want to track it separately
        Planet.MporeSalt_kg = Planet.MporeFluid_kg * Planet.Sil.wPore_ppt / 1e3
        # Combine these to get total salt content
        Planet.Msalt_kg = Planet.MoceanSalt_kg + Planet.MporeSalt_kg
        # The remainder, plus ice mass excluding gasses in clathrates,
        # is the total mass contained in water molecules for the body
        Planet.MH2O_kg = Planet.Mfluid_kg - Planet.Msalt_kg + Planet.Mice_kg - Planet.MclathGas_kg

        # Get the mean density of ocean layers and conducting/convecting upper ice layers
        Planet.VLayer_m3 = 4/3*np.pi * (Planet.r_m[:-1]**3 - Planet.r_m[1:]**3)
        Planet.Ocean.Vtot_m3 = np.sum(Planet.VLayer_m3[Planet.phase == 0])
        if Planet.Do.NO_H2O or Planet.Ocean.Vtot_m3 == 0:
            Planet.Ocean.rhoMean_kgm3 = np.nan
            Planet.Ocean.Tmean_K = np.nan
        else:
            Planet.Ocean.rhoMean_kgm3 = Planet.Mocean_kg / Planet.Ocean.Vtot_m3
            # Get average ocean temperature by summing the total heat energy in the
            # ocean and dividing by the total heat storage capacity
            oceanHeat_pK = Planet.Cp_JkgK[Planet.phase == 0] * Planet.MLayer_kg[Planet.phase == 0]
            Planet.Ocean.Tmean_K = np.sum(Planet.T_K[Planet.phase == 0] * oceanHeat_pK) / np.sum(oceanHeat_pK)

        # Get mean tidal heating in silicate layers
        Planet.Sil.HtidalMean_Wm3 = np.mean(Planet.Htidal_Wm3[iOS:iSC])

        # Check for any negative temperature gradient (indicates non-equilibrium conditions)
        gradTneg = np.where(np.diff(Planet.T_K) < 0)
        if np.any(gradTneg) and not Params.SKIP_INNER:
            log.warning(f'Negative temperature gradient starting at index {gradTneg[0][0]:d}. This indicates that ' +
                        'internal heating parameters Qrad_Wkg and/or Htidal_Wm3 are likely set too high to be consistent ' +
                        'with the heat flux into the ocean. The current configuration represents a ' +
                        'non-equilibrium state.')
        if Planet.invalidReason is None:
            Planet.invalidReason = 'Valid'

    else:
        # Set remaining quantities that are still None if the profile is invalid:
        Planet.Ocean.Tmean_K, Planet.Sil.HtidalMean_Wm3, Planet.Ocean.rhoMean_kgm3, Planet.Ocean.Vtot_m3, \
            Planet.D_km, Planet.CMR2mean, Planet.CMR2less, Planet.CMR2more, Planet.MH2O_kg, Planet.MclathGas_kg, \
            Planet.Mclath_kg, Planet.Mcore_kg, Planet.Mfluid_kg, Planet.Mice_kg, Planet.MoceanSalt_kg, \
            Planet.Mocean_kg, Planet.MporeFluid_kg, Planet.MporeSalt_kg, Planet.Mrock_kg, Planet.Msalt_kg, \
            Planet.Mtot_kg, Planet.Sil.Rmean_m, Planet.Sil.Rrange_m, Planet.Core.Rmean_m, Planet.Core.Rrange_m, \
            Planet.Sil.rhoMean_kgm3, Planet.Core.rhoMean_kgm3, Planet.Sil.GSmean_GPa, Planet.Core.GSmean_GPa, \
            Planet.Pseafloor_MPa, Planet.Sil.phiCalc_frac, Planet.phiSeafloor_frac = (np.nan for _ in range(32))
        Planet.Sil.Rtrade_m, Planet.Core.Rtrade_m, Planet.Sil.rhoTrade_kgm3 \
            = (np.array([np.nan]) for _ in range(3))
        Planet.Steps.nHydro = 1
        Planet.Steps.nSil = 0
        Planet.Steps.nTotal = 1

    return Planet, Params

def NonSelfConsistentInnerLayer(Planet, Params):
    """ Non-self-consistent inner layer modeling using specified mean properties
        instead of detailed EOS calculations. Uses user-specified layer densities,
        radii, and thermal properties.
        
        Only setup for basic silicate and core layers at the moment

        Assigns Planet attributes:
            Steps.nTotal, Sil.Rmean_m, Core.Rmean_m, Sil.rhoMean_kgm3, Core.rhoMean_kgm3,
            phase, r_m, z_m, g_ms2, T_K, P_MPa, rho_kgm3, MLayer_kg
    """
    log.debug('Evaluating non-self-consistent inner layer.')
    
    if Planet.Do.PARTIAL_DIFFERENTIATION or Planet.Do.NO_DIFFERENTIATION: #TODO: Revisit in future
        raise ValueError('Non-self-consistent inner layer modeling is not supported for partial or no differentiation modeling at this moment.')
    
    if Planet.Do.VALID:
        # Add final end layers to r_m and z_m so we can calculate masses
        Planet.r_m = np.append(Planet.r_m, [0])
        Planet.z_m = np.append(Planet.z_m, Planet.r_m[0])
        # Calculate position of top of silicate layer
        Planet.Sil.Rmean_m = Planet.Bulk.R_m - Planet.zb_km * 1000 - Planet.D_km * 1000
        dzSil_m = Planet.Sil.Rmean_m - Planet.Core.Rmean_m
        # Check if we are modeling a silicate layer
        if Planet.Sil.Rmean_m <= 0:
            Planet.Steps.nSil = 0
            pass
        else:
            # Set up derived layer arrays
            # NOTE how we don't include the 'next' layer in the derived arrays for the inner layers in order to prevent out of bounds indexing
            Planet.phase[Planet.Steps.nHydro:Planet.Steps.nHydro + Planet.Steps.nSil] = Constants.phaseSil
            Planet.z_m[Planet.Steps.nHydro:Planet.Steps.nHydro + Planet.Steps.nSil+1] = np.linspace(Planet.z_m[Planet.Steps.nHydro], Planet.z_m[Planet.Steps.nHydro] + dzSil_m, Planet.Steps.nSil+1)
            
            Planet.Sil.EOS = GetInnerEOS(Planet.Sil.mantleEOS, EOSinterpMethod=Params.lookupInterpMethod,
                                         kThermConst_WmK=Planet.Sil.kTherm_WmK, HtidalConst_Wm3=Planet.Sil.Htidal_Wm3,
                                         porosType=Planet.Sil.porosType, phiTop_frac=Planet.Sil.phiRockMax_frac,
                                         Pclosure_MPa=Planet.Sil.Pclosure_MPa, phiMin_frac=Planet.Sil.phiMin_frac,
                                         EXTRAP=Params.EXTRAP_SIL, etaSilFixed_Pas=Planet.Sil.etaRock_Pas, etaCoreFixed_Pas=[Planet.Core.etaFeSolid_Pas, Planet.Core.etaFeLiquid_Pas],
                                         TviscTrans_K=Planet.Sil.TviscTrans_K,
                                         doConstantProps=True, constantProperties={'rho_kgm3': Planet.Sil.rhoSilWithCore_kgm3, 'Cp_JkgK': np.nan, 'alpha_pK': np.nan, 'kTherm_WmK': Planet.Sil.kTherm_WmK,
                                                                                   'VP_kms': Planet.Sil.VPset_kms, 'VS_kms': Planet.Sil.VSset_kms, 'KS_GPa': Planet.Sil.KSset_GPa, 'GS_GPa': Planet.Sil.GSset_GPa, 'eta_Pas': Planet.Sil.etaRock_Pas,
                                                                                   'sigma_Sm': Planet.Sil.sigmaSil_Sm})
            
            # Propogate conduction
            Tbot_K = Planet.T_K[Planet.Steps.nHydro] #TODO Fix this to make it self consistent
            Planet = PropogateConductionFromDepth(Planet, Params, Planet.Steps.nHydro, Planet.Steps.nHydro + Planet.Steps.nSil, Tbot_K, Planet.Sil.EOS, propogateNextLayer=Planet.Do.Fe_CORE) #TODO Shouldn't use this equation to calculate tempearture
        if not Planet.Do.Fe_CORE:
            Planet.Steps.nCore = 0
        # Set up core layer
        else:
            # Set up derived layer arrays
            # NOTE how we have to reset 
            startCore = Planet.Steps.nHydro + Planet.Steps.nSil
            Planet.phase[startCore:startCore + Planet.Steps.nCore] = Constants.phaseFe
            Planet.z_m[startCore:startCore + Planet.Steps.nCore+1] = np.linspace(Planet.z_m[startCore-1]+Planet.Sil.Rmean_m, Planet.z_m[startCore-1]+Planet.Sil.Rmean_m+Planet.Core.Rmean_m, Planet.Steps.nCore+1)
            
            if Planet.Core.rhoMean_kgm3 is None:
                Planet.Core.rhoMean_kgm3 = Planet.Core.rhoFe_kgm3
                log.warning('Planet.Core.rhoMean_kgm3 is not set, using Planet.Core.rhoFe_kgm3.')
            if Planet.Core.kTherm_WmK is None:
                Planet.Core.kTherm_WmK = Constants.kThermFe_WmK
                log.warning('Planet.Core.kTherm_WmK is not set, using Constants.kThermFe_WmK.')
            # Set up derived layer arrays
            Planet.Core.EOS = GetInnerEOS(Planet.Core.coreEOS, EOSinterpMethod=Params.lookupInterpMethod, Fe_EOS=True,
                                          kThermConst_WmK=Planet.Core.kTherm_WmK, EXTRAP=Params.EXTRAP_Fe,
                                          wFeCore_ppt=Planet.Core.wFe_ppt, wScore_ppt=Planet.Core.wS_ppt, etaSilFixed_Pas=Planet.Sil.etaRock_Pas, etaCoreFixed_Pas=[Planet.Core.etaFeSolid_Pas, Planet.Core.etaFeLiquid_Pas],
                                          TviscTrans_K=Planet.Core.TviscTrans_K,
                                          doConstantProps=True, constantProperties={'rho_kgm3': Planet.Core.rhoFe_kgm3, 'Cp_JkgK': np.nan, 'alpha_pK': np.nan, 'kTherm_WmK': Planet.Core.kTherm_WmK,
                                                                                   'VP_kms': np.nan, 'VS_kms': np.nan, 'KS_GPa': np.nan, 'GS_GPa': Planet.Core.GSset_GPa, 'eta_Pas': Planet.Core.etaFeSolid_Pas,
                                                                                   'sigma_Sm': Planet.Core.sigmaCore_Sm})
            
            # Propogate conduction
            Tbot_K = Planet.T_K[startCore-1] #TODO Fix this to make it self consistent
            Planet = PropogateConductionFromDepth(Planet, Params, startCore, startCore + Planet.Steps.nCore, Tbot_K, Planet.Core.EOS, propogateNextLayer=False)
        # Calculate the mass of the final layer
        Planet.MLayer_kg[-1] = 4/3*np.pi * (Planet.r_m[-2]**3 - Planet.r_m[-1]**3) * Planet.rho_kgm3[-1]

        # Calculate total mass and CMR2
        Mtot_kg = np.sum(Planet.MLayer_kg)
        MR2_kgm2 = Planet.Bulk.M_kg * Planet.Bulk.R_m**2
        C_kgm2 = np.sum(8*np.pi/15 * Planet.rho_kgm3 * (Planet.r_m[:-1]**5 - Planet.r_m[1:]**5))
        CMR2 = C_kgm2 / MR2_kgm2
        
        # Calculate differences in mass and CMR2
        Mdiff_frac = np.abs(1 - Mtot_kg / Planet.Bulk.M_kg)
        CMR2diff_upper = Planet.Bulk.Cmeasured + Planet.Bulk.CuncertaintyUpper
        CMR2diff_lower = Planet.Bulk.Cmeasured - Planet.Bulk.CuncertaintyLower
        MdiffThresh = 0.05
        if Mdiff_frac > MdiffThresh:
            if Mtot_kg - Planet.Bulk.M_kg > 0:
                Planet.Do.VALID = False
                Planet.invalidReason = f'Mass of planet is more than {100 * MdiffThresh:g}% more than the total body mass.'
                invalidMessage = Planet.invalidReason + ' Try lowering the density of layers or increasing the ice thickness, which should have the lowest density.'
            else:
                Planet.Do.VALID = False
                Planet.invalidReason = f'Mass of planet is more than {100 * MdiffThresh:g}% less than the total body mass.'
                invalidMessage = Planet.invalidReason + ' Try increasing the density of layers or increasing the inner layers thickness, which should have the highest density.'
            if Params.ALLOW_BROKEN_MODELS:
                if Params.DO_EXPLOREOGRAM or Params.DO_INDUCTOGRAM or Params.DO_MONTECARLO:
                    log.info(invalidMessage + ' Params.ALLOW_BROKEN_MODELS is True, so calculations will proceed ' +
                                             'with many values set to nan.')
                else:
                    log.error(invalidMessage + ' Params.ALLOW_BROKEN_MODELS is True, so calculations will proceed ' +
                                             'with many values set to nan.')
                Planet.Do.STILL_CALCULATE_BROKEN_PROPERTIES = True
            else:
                raise RuntimeError(invalidMessage)
        elif CMR2diff_upper > CMR2 or CMR2 < CMR2diff_lower:
            if CMR2diff_upper > CMR2:
                Planet.Do.VALID = False
                Planet.invalidReason = f'CMR2 of planet is more than {100 * CMR2diff_upper:.1f}% more than the measured value.'
                invalidMessage = Planet.invalidReason + ' Try lowering the density of layers or increasing the inner layers thickness, which should have the highest density.'
            else:
                Planet.Do.VALID = False
                Planet.invalidReason = f'CMR2 of planet is less than {100 * CMR2diff_lower:.1f}% less than the measured value.'
                invalidMessage = Planet.invalidReason + ' Try increasing the density of layers or increasing the inner layers thickness, which should have the highest density.'
            if Params.ALLOW_BROKEN_MODELS:
                if Params.DO_EXPLOREOGRAM or Params.DO_INDUCTOGRAM or Params.DO_MONTECARLO:
                    log.info(invalidMessage + ' Params.ALLOW_BROKEN_MODELS is True, so calculations will proceed ' +
                                             'with many values set to nan.')
                else:
                    log.error(invalidMessage + ' Params.ALLOW_BROKEN_MODELS is True, so calculations will proceed ' +
                                             'with many values set to nan.')
                Planet.Do.STILL_CALCULATE_BROKEN_PROPERTIES = True
            else:
                raise RuntimeError(invalidMessage)
        if not Planet.Do.VALID and not Params.ALLOW_BROKEN_MODELS:
            Planet.CMR2mean = np.nan
            Planet.CMR2less = Planet.CMR2mean
            Planet.CMR2more = Planet.CMR2mean

        else:
            Planet.CMR2mean = CMR2
            Planet.CMR2less = CMR2
            Planet.CMR2more = CMR2
        nans = np.array([np.nan])
        Planet.Sil.rhoTrade_kgm3 = nans
        Planet.Sil.Rtrade_m = nans
        Planet.Sil.Rrange_m = np.nan
        Planet.Core.Rtrade_m = nans
        Planet.Core.Rrange_m = np.nan
            
        
        # Simplified mass calculations for compatibility
        Planet.Mtot_kg = Mtot_kg
        Planet.Mice_kg = np.sum(Planet.MLayer_kg[:Planet.Steps.nSurfIce])
        Planet.Mcore_kg = np.sum(Planet.MLayer_kg[Planet.Steps.nHydro + Planet.Steps.nSil:])
        Planet.Mrock_kg = np.sum(Planet.MLayer_kg[Planet.Steps.nHydro:Planet.Steps.nHydro + Planet.Steps.nSil])
        Planet.Mfluid_kg = Planet.Mtot_kg - Planet.Mcore_kg - Planet.Mrock_kg - Planet.Mice_kg
        #Planet.Mclath_kg = 0.0
        #Planet.MclathGas_kg = 0.0
        #Planet.Mocean_kg = Planet.Mfluid_kg
        #Planet.MporeFluid_kg = 0.0
        #Planet.MporeSalt_kg = 0.0
        #Planet.Msalt_kg = 0.0
        Planet.Mocean_kg = np.sum(Planet.MLayer_kg[Planet.phase == 0])
        Planet.MporeFluid_kg = Planet.Mfluid_kg - Planet.Mocean_kg
        #Planet.Mclath_kg = 0.0
        #Planet.MclathGas_kg = 0.0
        #Planet.MoceanSalt_kg = Planet.Mfluid_kg - Planet.MporeFluid_kg
        #Planet.MporeSalt_kg = Planet.Mfluid_kg - Planet.MporeFluid_kg
        #Planet.Msalt_kg = Planet.Mfluid_kg - Planet.MporeFluid_kg
        #Planet.MH2O_kg = Planet.Mfluid_kg
        #Planet.Mtot_kg = Planet.Mcore_kg + Planet.Mrock_kg + Planet.Mfluid_kg
        
        # Set remaining mean properties
        Planet.Ocean.Tmean_K = Planet.Bulk.Tb_K
        # Get the mean density of ocean layers and conducting/convecting upper ice layers
        Planet.VLayer_m3 = 4/3*np.pi * (Planet.r_m[:-1]**3 - Planet.r_m[1:]**3)
        Planet.Ocean.Vtot_m3 = np.sum(Planet.VLayer_m3[Planet.phase == 0])
        # Set validity flag
        if Planet.invalidReason is None:
            Planet.invalidReason = 'Valid'
        
        log.info(f'Non-self-consistent inner layers: ' +
                 f'Sil radius: {Planet.Sil.Rmean_m/1e3:.1f} km, ' +
                 f'Core radius: {Planet.Core.Rmean_m/1e3:.1f} km.')
        
    else:
        # Set remaining quantities that are still None if the profile is invalid:
        Planet.Ocean.Tmean_K, Planet.Sil.HtidalMean_Wm3, Planet.Ocean.rhoMean_kgm3, Planet.Ocean.Vtot_m3, \
            Planet.D_km, Planet.CMR2mean, Planet.CMR2less, Planet.CMR2more, Planet.MH2O_kg, Planet.MclathGas_kg, \
            Planet.Mclath_kg, Planet.Mcore_kg, Planet.Mfluid_kg, Planet.Mice_kg, Planet.MoceanSalt_kg, \
            Planet.Mocean_kg, Planet.MporeFluid_kg, Planet.MporeSalt_kg, Planet.Mrock_kg, Planet.Msalt_kg, \
            Planet.Mtot_kg, Planet.Sil.Rmean_m, Planet.Sil.Rrange_m, Planet.Core.Rmean_m, Planet.Core.Rrange_m, \
            Planet.Sil.rhoMean_kgm3, Planet.Core.rhoMean_kgm3, Planet.Sil.GSmean_GPa, Planet.Core.GSmean_GPa, \
            Planet.Pseafloor_MPa, Planet.Sil.phiCalc_frac, Planet.phiSeafloor_frac = (np.nan for _ in range(32))
        Planet.Sil.Rtrade_m, Planet.Core.Rtrade_m, Planet.Sil.rhoTrade_kgm3 \
            = (np.array([np.nan]) for _ in range(3))
        Planet.Steps.nHydro = 1
        Planet.Steps.nSil = 0
        Planet.Steps.nTotal = 1

    return Planet, Params
def CalcMoIConstantRho(Planet, Params):
    """ Find the relative sizes of silicate, core, and hydrosphere layers that are
        consistent with the measured moment of inertia, based on calculated hydrosphere
        properties and assumptions about the silicate and possible core layers.

        Assigns Planet attributes:
            CMR2mean, Sil.RsilMean_m, Sil.RsilRange_m, Core.RFeMean_m, Core.RFeRange_m, Steps.nHydro
    """
    log.debug('Finding MoI consistent with measured value for constant-density inner layers...')
    # Get MR^2 -- we will need to divide each C by this later.
    MR2_kgm2 = Planet.Bulk.M_kg * Planet.Bulk.R_m**2

    # Get final number of layers modeled in "overshoot" hydrosphere
    if Planet.Do.NO_H2O:
        Planet.Steps.iSilStart = 0
        nHydroActual = 2
    else:
        nHydroActual = Planet.Steps.nSurfIce + Planet.Steps.nOceanMax
    # Find contribution to axial moment of inertia C from each ocean layer
    dChydro_kgm2 = 8*np.pi/15 * Planet.rho_kgm3[:-1] * (Planet.r_m[:-1]**5 - Planet.r_m[1:]**5)
    # Find total mass contained above each hydrosphere layer
    MAbove_kg = np.array([np.sum(Planet.MLayer_kg[:i]) for i in range(nHydroActual)])
    # Find volume of a full sphere of silicate corresponding to each valid layer
    VsilSphere_m3 = 4/3*np.pi * Planet.r_m[Planet.Steps.iSilStart:]**3

    if Planet.Do.Fe_CORE:
        # Find core bulk density based on assumed sulfide content
        rhoCore_kgm3 = Planet.Core.rhoFeS_kgm3 * Planet.Core.rhoFe_kgm3 \
            * (Planet.Core.xFeS * Constants.m_gmol['FeS'] + (1 - Planet.Core.xFeS) * Constants.m_gmol['Fe'] ) \
            / (Planet.Core.xFeS * Constants.m_gmol['FeS'] * Planet.Core.rhoFe_kgm3 + (1 - Planet.Core.xFeS) * Constants.m_gmol['Fe'] * Planet.Core.rhoFeS_kgm3)
            # / (Planet.Core.xFeS * (Planet.Core.rhoFe_kgm3 - Planet.Core.rhoFeS_kgm3) + Planet.Core.rhoFeS_kgm3)  # Vance et al. (2014) Eq. 10
        # Calculate core volume for a silicate layer with outer radius equal to bottom of each hydrosphere layer
        # and inner radius equal to the core radius
        VCore_m3 = np.array([(Planet.Bulk.M_kg - MAbove_kg[i] - VsilSphere_m3[i-Planet.Steps.iSilStart] *
                             Planet.Sil.rhoSilWithCore_kgm3) / (rhoCore_kgm3 - Planet.Sil.rhoSilWithCore_kgm3)
                             for i in range(Planet.Steps.iSilStart, nHydroActual-1)])
        # Find values for which the silicate radius is too large
        try:
            nTooBig = next((i[0] for i, val in np.ndenumerate(VCore_m3) if val>0))
        except StopIteration:
            msg = f'Failed to find a core size consistent with rhoSil = {Planet.Sil.rhoSilWithCore_kgm3:.1f} kg/m3 ' + \
                  f'and xFeS = {Planet.Core.xFeS:.3f} for PHydroMax_MPa = {Planet.Ocean.PHydroMax_MPa:.1f}. ' #+ \
                   #'Core size will be set to zero.'
            if Params.DO_EXPLOREOGRAM:
                log.debug(msg)
            else:
                raise ValueError(msg)
            nTooBig = 0
            rCore_m = np.zeros_like(VCore_m3)
            Planet.Steps.nCore = 0
            INVALIDCORE = True
        else:
            INVALIDCORE = False
            # Calculate corresponding core radii based on above density
            rCore_m = (VCore_m3[nTooBig:]*3/4/np.pi)**(1/3)

        # Assign fixed density to an array for dual-use code looking for compatible C/MR^2
        rhoSil_kgm3 = np.ones_like(rCore_m) * Planet.Sil.rhoSilWithCore_kgm3
    else:
        # Find silicate density consistent with observed bulk mass for each radius
        rhoSil_kgm3 = np.array([(Planet.Bulk.M_kg - MAbove_kg[i]) / VsilSphere_m3[i-Planet.Steps.iSilStart]
                              for i in range(Planet.Steps.iSilStart, nHydroActual-1)])
        # Density of silicates is scaled to fit the total mass, so there is no nTooBig in this case.
        nTooBig = 0
        # Set core radius and density to zero so calculations can proceed
        rCore_m = np.zeros(nHydroActual-1 - Planet.Steps.iSilStart)
        rhoCore_kgm3 = 0
        INVALIDCORE = False

    # Calculate C for a mantle extending up to each hydrosphere layer in turn
    C_kgm2 = np.zeros(nHydroActual - 1)
    C_kgm2[Planet.Steps.iSilStart + nTooBig:] = [np.sum(dChydro_kgm2[:i + Planet.Steps.iSilStart + nTooBig + 1]) +
            8*np.pi/15 * rhoSil_kgm3[i] * (Planet.r_m[i + Planet.Steps.iSilStart + nTooBig]**5 - rCore_m[i]**5) +
            8*np.pi/15 * rhoCore_kgm3 * rCore_m[i]**5
            for i in range(nHydroActual - Planet.Steps.iSilStart - nTooBig - 1)]
    CMR2 = C_kgm2 / MR2_kgm2

    CMR2inds = [i[0] for i, valCMR2 in np.ndenumerate(CMR2)
                 if valCMR2 >= Planet.Bulk.Cmeasured - Planet.Bulk.CuncertaintyLower
                and valCMR2 <= Planet.Bulk.Cmeasured + Planet.Bulk.CuncertaintyUpper]
    if len(CMR2inds) == 0 or INVALIDCORE:
        if Planet.Do.NO_H2O:
            suggestion = '\nTry adjusting properties of silicates and core to get C/MR^2 values in range.'
        else:
            suggestion = '\nTry adjusting properties of silicates and core to get C/MR^2 values in range. ' + \
                         'Increasing PHydroMax_MPa can also lower C/MR^2 values.'
        msg = f'No MoI found matching C/MR^2 = {Planet.Bulk.Cmeasured:.3f}+{Planet.Bulk.CuncertaintyUpper:.3f}-{Planet.Bulk.CuncertaintyLower}. ' + \
                  f'Min: {np.min(CMR2[CMR2>0]):.3f}, Max: {np.max(CMR2):.3f}.'
        if Params.ALLOW_BROKEN_MODELS:
            fullMsg = msg + suggestion + ' Params.ALLOW_BROKEN_MODELS is True, so calculations will proceed with many values set to nan.'
            if Params.DO_EXPLOREOGRAM or Params.DO_INDUCTOGRAM or Params.DO_MONTECARLO:
                log.info(msg)
            else:
                log.error(fullMsg)
            Planet.Do.VALID = False
            if Planet.invalidReason is None:
                Planet.invalidReason = ''
            else:
                Planet.invalidReason += '. '
            Planet.invalidReason += f'No valid MoI was found between {Planet.Bulk.Cmeasured - Planet.Bulk.CuncertaintyLower:.4f} and {Planet.Bulk.Cmeasured + Planet.Bulk.CuncertaintyUpper:.4f}'
        else:
            raise RuntimeError(msg + suggestion)

        nans = np.array([np.nan])
        Planet.CMR2mean = np.nan
        Planet.CMR2less = Planet.CMR2mean
        Planet.CMR2more = Planet.CMR2mean
        Planet.Sil.rhoTrade_kgm3 = nans
        Planet.Sil.Rmean_m = np.nan
        Planet.Sil.Rtrade_m = nans
        Planet.Sil.Rrange_m = np.nan
        Planet.Core.Rmean_m = np.nan
        Planet.Core.Rtrade_m = nans
        Planet.Core.Rrange_m = np.nan
        Planet.Sil.rhoNoCore_kgm3 = np.nan
        Planet.Steps.nSil = Planet.Steps.nSilMax
        # Use Rset_m to indicate that we have already determined the core size in using SilicateLayers
        Planet.Core.Rset_m = np.nan
        iCMR2 = 0
        iCMR2inner = 0
        if Planet.Do.NO_H2O:
            Planet.Steps.nHydro = 0
        else:
            Planet.Steps.nHydro = Planet.Steps.nOceanMax

    else:
        if Planet.Do.HYDROSPHERE_THICKNESS:
            # CMR2inds = np.nonzero(CMR2) # remove constraint of CMR2 in specified bounds
            hydrosphere_thickness_m = Planet.Bulk.R_m  - Planet.r_m[CMR2inds]
            CMR2diff = np.abs(hydrosphere_thickness_m - Planet.Bulk.Dhsphere_m)
            # Get index of closest match in CMR2inds
            iCMR2ind = np.argmin(CMR2diff)
        else:# Find the C/MR^2 value most closely matching the measured value
            CMR2diff = np.abs(CMR2[CMR2inds] - Planet.Bulk.Cmeasured)
            # Get index of closest match in CMR2inds
            iCMR2ind = np.argmin(CMR2diff)
            # Find Planet array index corresponding to closest matching value
        CMR2indsInner = [ind - Planet.Steps.iSilStart - nTooBig for ind in CMR2inds]
        iCMR2 = CMR2inds[iCMR2ind]
        iCMR2inner = iCMR2 - Planet.Steps.iSilStart - nTooBig
        # Record the best-match C/MR^2 value
        Planet.CMR2mean = CMR2[iCMR2]
        # We don't have neighboring values because we used the MoI to calculate properties
        Planet.CMR2less = Planet.CMR2mean
        Planet.CMR2more = Planet.CMR2mean
        # Record interior sizes
        Planet.Sil.rhoTrade_kgm3 = rhoSil_kgm3[CMR2indsInner]
        Planet.Sil.Rmean_m = Planet.r_m[iCMR2]
        Planet.Sil.Rtrade_m = Planet.r_m[CMR2inds]
        Planet.Sil.Rrange_m = np.max(Planet.Sil.Rtrade_m) - np.min(Planet.Sil.Rtrade_m)
        Planet.Core.Rmean_m = rCore_m[iCMR2inner]
        Planet.Core.Rtrade_m = rCore_m[CMR2indsInner]
        Planet.Core.Rrange_m = np.max(Planet.Core.Rtrade_m) - np.min(Planet.Core.Rtrade_m)
        Planet.Sil.rhoNoCore_kgm3 = rhoSil_kgm3[iCMR2inner]
        # Now we finally know how many layers there are in the hydrosphere
        Planet.Steps.nHydro = iCMR2
        # Use Rset_m to indicate that we have already determined the core size in using SilicateLayers
        Planet.Core.Rset_m = Planet.Core.Rmean_m + 0.0

    if not Params.SKIP_INNER:
        Planet.Sil.fn_phi_frac = GetphiCalc(Planet.Sil.phiRockMax_frac, Planet.Sil.EOS.fn_phi_frac, Planet.Sil.phiMin_frac)
        Planet.Sil.fn_Htidal_Wm3 = GetHtidalFunc(Planet.Sil.Htidal_Wm3)  # Placeholder until we implement a self-consistent calc
        # Evaluate the silicate EOS for each layer
        indsSilValid, nProfiles, Psil_MPa, Tsil_K, rSil_m, rhoSilEOS_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, \
        phiSil_frac, HtidalSil_Wm3, kThermSil_WmK, Ppore_MPa, rhoSilMatrix_kgm3, rhoPore_kgm3, phasePore \
            = SilicateLayers(Planet, Params)
        nSilTooBig = nProfiles - np.size(indsSilValid)

        if Planet.Do.Fe_CORE and Planet.Steps.nCore > 0:
            # Evaluate the core EOS for each layer
            nSilFinal, Pcore_MPa, Tcore_K, rCoreEOS_m, rhoCore_kgm3, MLayerCore_kg, gCore_ms2, CpCore_JkgK, alphaCore_pK, \
                kThermCore_WmK = IronCoreLayers(Planet, Params,
                               indsSilValid, nSilTooBig, nProfiles, Psil_MPa, Tsil_K, rSil_m, MAboveSil_kg, gSil_ms2)

            MtotCore_kg = np.sum(MLayerCore_kg)
            Planet.Core.rhoMean_kgm3 = MtotCore_kg / VCore_m3[nTooBig:][iCMR2inner]

            coreProps = (Pcore_MPa, Tcore_K, rCoreEOS_m[0,:-1], gCore_ms2, rhoCore_kgm3, CpCore_JkgK, alphaCore_pK,
                         kThermCore_WmK, MLayerCore_kg)
        else:
            nSilFinal = Planet.Steps.nSilMax
            MtotCore_kg = 0
            Planet.Core.rhoMean_kgm3 = 0
            Planet.Core.Rtrade_m = np.zeros_like(Planet.Sil.Rtrade_m)
            Planet.Core.Rrange_m = 0
            coreProps = None

        Planet.Steps.nSil = nSilFinal
        # Fill core/mantle trade arrays and set mean values consistent with MoI
        MtotSil_kg = np.sum(MLayerSil_kg[0,:Planet.Steps.nSil])
        Planet.Sil.rhoMean_kgm3 = MtotSil_kg / (4/3*np.pi * (rSil_m[0,0]**3 - rSil_m[0,-1]**3))
        Planet.Mtot_kg = np.sum(Planet.MLayer_kg[:iCMR2]) + MtotSil_kg + MtotCore_kg
        if not np.isnan(Planet.CMR2mean):
            log.info(f'Found matching MoI of {Planet.CMR2mean:.4f} ' +
                     f'(C/MR^2 = {Planet.CMR2strPrint}) for ' +
                     f'R_sil = {Planet.Sil.Rmean_m / Planet.Bulk.R_m:.2f} R_{Planet.name[0]}, ' +
                     f'R_core = {Planet.Core.Rmean_m / Planet.Bulk.R_m:.2f} R_{Planet.name[0]}, ' +
                     f'rho_sil (found) = {rhoSil_kgm3[iCMR2inner]:.0f} kg/m^3, ' +
                     f'rho_sil (actual) = {Planet.Sil.rhoMean_kgm3:.0f} kg/m^3, ' +
                     f'P_HydroMax = {Planet.P_MPa[iCMR2]:.1f} MPa, ' +
                     f'M_tot = {Planet.Mtot_kg/Planet.Bulk.M_kg:.4f} M_{Planet.name[0]}.')
            log.warning('Because silicate and core properties were determined from the EOS after finding their ' +
                        'sizes by assuming constant densities, the body mass may not match the measured value.')

    else:
        if not np.isnan(Planet.CMR2mean):
            log.info(f'Found matching MoI of {Planet.CMR2mean:.4f} ' +
                     f'(C/MR^2 = {Planet.CMR2strPrint}) for ' +
                     f'R_sil = {Planet.Sil.Rmean_m / Planet.Bulk.R_m:.2f} R_{Planet.name[0]}, ' +
                     f'R_core = {Planet.Core.Rmean_m / Planet.Bulk.R_m:.2f} R_{Planet.name[0]}, ' +
                     f'rho_sil = {rhoSil_kgm3[iCMR2inner]:.0f} kg/m^3, ' +
                     f'P_HydroMax = {Planet.P_MPa[iCMR2]:.1f} MPa, ' +
                     f'M_tot = {1.0:.5f} M_{Planet.name[0]} (fixed).')
        if not Params.DO_EXPLOREOGRAM:
            log.debug('Params.SKIP_INNER is True, assigning interior properties to 0.')
        Planet.Steps.nSil = Planet.Steps.nSilMax
        Psil_MPa, Tsil_K, rhoSilEOS_kgm3, gSil_ms2, phiSil_frac, kThermSil_WmK, Ppore_MPa, rhoSilMatrix_kgm3, \
            rhoPore_kgm3, HtidalSil_Wm3, MLayerSil_kg \
            = (np.zeros((1,Planet.Steps.nSil)) for _ in range(11))
        phasePore = np.zeros((1, Planet.Steps.nSil), dtype=np.int_)
        rSil_m = np.zeros((1, Planet.Steps.nSil+1))
        rSil_m[0,0] = Planet.Sil.Rmean_m
        coreProps = (np.zeros(Planet.Steps.nCore) for _ in range(9))
        Planet.Sil.rhoMean_kgm3 = rhoSil_kgm3[iCMR2inner]
        Planet.Core.rhoMean_kgm3 = rhoCore_kgm3
        Planet.Mtot_kg = Planet.Bulk.M_kg

    if Planet.Do.POROUS_ROCK:
        Planet.Sil.phiCalc_frac = Planet.Sil.phiRockMax_frac
    else:
        Planet.Sil.phiCalc_frac = np.nan

    mantleProps = (Psil_MPa[0,:Planet.Steps.nSil], Tsil_K[0,:Planet.Steps.nSil], rSil_m[0,:Planet.Steps.nSil],
                   rhoSilEOS_kgm3[0,:Planet.Steps.nSil], gSil_ms2[0,:Planet.Steps.nSil],
                   phiSil_frac[0,:Planet.Steps.nSil], HtidalSil_Wm3[0,:Planet.Steps.nSil],
                   kThermSil_WmK[0,:Planet.Steps.nSil], Ppore_MPa[0,:Planet.Steps.nSil],
                   rhoSilMatrix_kgm3[0,:Planet.Steps.nSil], rhoPore_kgm3[0,:Planet.Steps.nSil],
                   MLayerSil_kg[0,:Planet.Steps.nSil], phasePore[0,:Planet.Steps.nSil])

    return Planet, mantleProps, coreProps


def CalcMoIWithEOS(Planet, Params):
    """ Find the relative sizes of silicate, core, and hydrosphere layers that are
        consistent with the measured moment of inertia, based on calculated hydrosphere
        properties and EOS data for assumed mantle and core compositions output by Perple_X.

        For models with a core, consistency with the MoI and body mass is determined with
        Sil.EOS and Core.EOS, Ocean.QfromMantle_W, Sil.Qrad_Wkg, and Sil.Htidal_Wm3. The radii
        of the silicate and core layers are treated as the free variables in 2 equations to match
        the 2 "unknowns", MoI and body mass.
        For models without a core, consistency with the MoI and body mass is determined with
        Sil.EOS, Ocean.QfromMantle_W, and Sil.Qrad_Wkg. Sil.Htidal_Wm3 and the radius of the
        silicate layer are treated as the free variables in 2 equations to match MoI and body mass.
        Sil.Htidal_Wm3 is started at 0 and increased until the thermal profile extends beyond the
        T domain of the Perplex_X EOS.

        Assigns Planet attributes:
            CMR2mean, Sil.RsilMean_m, Sil.RsilRange_m, Core.RFeMean_m, Core.RFeRange_m, Steps.nHydro, Steps.nSil,
            all layer arrays
    """
    log.debug('Finding MoI consistent with measured value...')
    # Get MR^2 -- we will need to divide each C by this later.
    MR2_kgm2 = Planet.Bulk.M_kg * Planet.Bulk.R_m**2

    # Find contribution to axial moment of inertia C from each ocean layer
    dCfromH2O_kgm2 = 8*np.pi/15 * Planet.rho_kgm3[:-1] * (Planet.r_m[:-1]**5 - Planet.r_m[1:]**5)

    if Planet.Do.Fe_CORE:
        Planet.Sil.fn_Htidal_Wm3 = GetHtidalFunc(Planet.Sil.Htidal_Wm3)  # Placeholder until we implement a self-consistent calc
        Planet.Sil.fn_phi_frac = GetphiCalc(Planet.Sil.phiRockMax_frac, Planet.Sil.EOS.fn_phi_frac, Planet.Sil.phiMin_frac)
        # Propagate the silicate EOS from each hydrosphere layer to the center of the body
        log.debug(f'Propagating silicate EOS for each possible mantle size ({Planet.Steps.nHydroMax-Planet.Steps.iSilStart} options)...')
        indsSilValid, nProfiles, Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, \
        phiSil_frac, HtidalSil_Wm3, kThermSil_WmK, PsilPore_MPa, rhoSilMatrix_kgm3, rhoSilPore_kgm3, phaseSilPore \
            = SilicateLayers(Planet, Params)
        if not Planet.Do.VALID and Planet.Steps.iSilStart > 1 and np.size(indsSilValid) != 0:
            Planet.Steps.iSilStart = 1
            indsSilValid, nProfiles, Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, \
            phiSil_frac, HtidalSil_Wm3, kThermSil_WmK, PsilPore_MPa, rhoSilMatrix_kgm3, rhoSilPore_kgm3, phaseSilPore \
                = SilicateLayers(Planet, Params)
        nSilTooBig = nProfiles - np.size(indsSilValid)
        # Propagate the core EOS from each silicate layer at the max core radius to the center of the body
        nSilFinal, Pcore_MPa, Tcore_K, rCore_m, rhoCore_kgm3, MLayerCore_kg, gCore_ms2, CpCore_JkgK, alphaCore_pK, \
            kThermCore_WmK = IronCoreLayers(Planet, Params,
                           indsSilValid, nSilTooBig, nProfiles, Psil_MPa, Tsil_K, rSil_m, MAboveSil_kg, gSil_ms2)

        dCfromCore_kgm2 = 8*np.pi/15 * rhoCore_kgm3 * (rCore_m[:,:-1]**5 - rCore_m[:,1:]**5)
        Ccore_kgm2 = np.array([np.sum(dCfromCore_kgm2[i,:]) for i in range(nProfiles - nSilTooBig)])

        # Get indices of valid silicate portions of the layer profile
        iValid = np.array([Planet.Steps.iSilStart + i for i in indsSilValid]).astype(np.int_)
        C_kgm2 = np.zeros(nProfiles + Planet.Steps.iSilStart)
    else:
        # Propagate the silicate EOS from each hydrosphere layer to the center of the body

        if Planet.Do.POROUS_ROCK and not Planet.Do.FIXED_POROSITY:
            # In this case, we use the user-specified phiTop_frac value as a "middle" option,
            # and vary between 1/Planet.Sil.phiRangeMult times the difference from this middle
            # value to the endpoints, 0 and 1, to get the range of porosities to model.
            # For example, the user specifies 0.22 and phiRangeMult = 5. phiTop will then
            # take values from a min value of 0.22 - 0.22 / 5 = 0.176 to a max value of
            # 0.22 + (1 - 0.22) / 5 = 0.376.
            thisHtidal_Wm3 = Planet.Sil.Htidal_Wm3
            phiMin_frac = Planet.Sil.phiRockMax_frac - Planet.Sil.phiRockMax_frac / Planet.Sil.phiRangeMult
            phiMax_frac = Planet.Sil.phiRockMax_frac + (1 - Planet.Sil.phiRockMax_frac) / Planet.Sil.phiRangeMult
            multphi_frac = (phiMax_frac/phiMin_frac)**(1/Planet.Steps.nPoros)
            log.debug(f'Propagating silicate EOS for each possible mantle size and porosity from ' +
                      f'phiVac = {phiMin_frac:.3f} to {phiMax_frac:.3f} in {Planet.Steps.nPoros} ' +
                      f'steps...')
        else:
            # In this case, we will use Sil.HtidalMin_Wm3 and Sil.deltaHtidal_logUnits to get
            # a valid set of profiles.
            HtidalStart_Wm3 = Planet.Sil.HtidalMin_Wm3
            multHtidal_Wm3 = 10**Planet.Sil.deltaHtidal_logUnits
            thisHtidal_Wm3 = 0
            nHsteps = np.floor(np.log10(Planet.Sil.HtidalMax_Wm3/HtidalStart_Wm3)/Planet.Sil.deltaHtidal_logUnits) + 1
            log.debug(f'Propagating silicate EOS for each possible mantle size and heating from Htidal = {thisHtidal_Wm3:.2e} to {Planet.Sil.HtidalMax_Wm3:.2e} W/m^3 in {nHsteps} steps...')
            phiMin_frac = Planet.Sil.phiRockMax_frac
            phiMax_frac = phiMin_frac

        nProfiles = 0

        # Start at minimum tidal heating/porosity and initialize arrays
        thisphiTop_frac = phiMax_frac + 0.0
        Planet.Sil.fn_Htidal_Wm3 = GetHtidalFunc(thisHtidal_Wm3)  # Placeholder until we implement a self-consistent calc
        Planet.Sil.fn_phi_frac = GetphiCalc(Planet.Sil.phiRockMax_frac, Planet.Sil.EOS.fn_phi_frac, Planet.Sil.phiMin_frac)
        Planet.Sil.fn_phi_frac.update(thisphiTop_frac)
        indsSilValidTemp, _, PsilTemp_MPa, TsilTemp_K, rSilTemp_m, rhoSilTemp_kgm3, MLayerSilTemp_kg, MAboveSilTemp_kg, gSilTemp_ms2, \
        phiSilTemp_frac, HtidalSilTemp_Wm3, kThermSilTemp_WmK, PporeTemp_MPa, rhoSilMatrixTemp_kgm3, rhoSilPoreTemp_kgm3, phaseSilPoreTemp \
            = SilicateLayers(Planet, Params)
        if (np.size(indsSilValidTemp) == 0):
            masses = np.sum(MLayerSilTemp_kg, axis=1)
            msg = 'No silicate mantle size was less than the total body mass for the initialization ' + \
                 f'setting of {thisHtidal_Wm3:.2e} W/m^3 tidal heating\nand the expected maximum porosity ' + \
                 f'of {thisphiTop_frac:.3f}. The min silicate mass was {np.min(masses)/Planet.Bulk.M_kg:.3f} M_P. ' + \
                  'Try adjusting run settings that affect mantle density,\nlike porosity, silicate ' + \
                  'composition, and radiogenic heat flux.'
            if Params.ALLOW_BROKEN_MODELS:
                if Params.DO_EXPLOREOGRAM:
                    log.info(msg)
                else:
                    log.error(msg + ' Params.ALLOW_BROKEN_MODELS is True, so calculations will proceed ' +
                              'with many values set to nan.')
                Planet.Do.VALID = False
                Planet.invalidReason = f'All mantle models exceeded total body mass'
            else:
                raise RuntimeError(msg)

        # Initialize empty arrays for silicate layer properties
        Psil_MPa, Tsil_K, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, phiSil_frac, HtidalSil_Wm3, \
            kThermSil_WmK, PsilPore_MPa, rhoSilMatrix_kgm3, rhoSilPore_kgm3, phaseSilPore, phiTop_frac \
            = (np.empty((0, Planet.Steps.nSilMax)) for _ in range(14))
        rSil_m = np.empty((0, Planet.Steps.nSilMax+1))
        iValid = [indsSilValidTemp + Planet.Steps.iSilStart]

        while np.size(indsSilValidTemp) != 0 and (thisHtidal_Wm3 <= Planet.Sil.HtidalMax_Wm3 and thisphiTop_frac >= phiMin_frac):
            nProfiles += 1
            # Append the newest mass matches alongside existing ones --
            # do this in the loop so that we don't append the final,
            # invalid profile
            Psil_MPa = np.vstack([Psil_MPa, PsilTemp_MPa[indsSilValidTemp,:]])
            Tsil_K = np.vstack([Tsil_K, TsilTemp_K[indsSilValidTemp,:]])
            rSil_m = np.vstack([rSil_m, rSilTemp_m[indsSilValidTemp,:]])
            rhoSil_kgm3 = np.vstack([rhoSil_kgm3, rhoSilTemp_kgm3[indsSilValidTemp,:]])
            MLayerSil_kg = np.vstack([MLayerSil_kg, MLayerSilTemp_kg[indsSilValidTemp,:]])
            MAboveSil_kg = np.vstack([MAboveSil_kg, MAboveSilTemp_kg[indsSilValidTemp,:]])
            gSil_ms2 = np.vstack([gSil_ms2, gSilTemp_ms2[indsSilValidTemp,:]])
            phiSil_frac = np.vstack([phiSil_frac, phiSilTemp_frac[indsSilValidTemp,:]])
            HtidalSil_Wm3 = np.vstack([HtidalSil_Wm3, HtidalSilTemp_Wm3[indsSilValidTemp,:]])
            kThermSil_WmK = np.vstack([kThermSil_WmK, kThermSilTemp_WmK[indsSilValidTemp,:]])
            PsilPore_MPa = np.vstack([PsilPore_MPa, PporeTemp_MPa[indsSilValidTemp,:]])
            rhoSilMatrix_kgm3 = np.vstack([rhoSilMatrix_kgm3, rhoSilMatrixTemp_kgm3[indsSilValidTemp,:]])
            rhoSilPore_kgm3 = np.vstack([rhoSilPore_kgm3, rhoSilPoreTemp_kgm3[indsSilValidTemp,:]])
            phaseSilPore = np.vstack([phaseSilPore, phaseSilPoreTemp[indsSilValidTemp,:]])
            phiTop_frac = np.append(phiTop_frac, thisphiTop_frac)

            if Planet.Do.POROUS_ROCK and not Planet.Do.FIXED_POROSITY:
                rSilOuter_m = rSilTemp_m[indsSilValidTemp, 0]
                log.debug(f'Silicate match for phiTop = {thisphiTop_frac:.3f} with ' +
                          f'rSil = {rSilOuter_m / 1e3:.1f} km ({rSilOuter_m / Planet.Bulk.R_m:.3f} R_{Planet.name[0]}) ' +
                          f'at PHydroMax_MPa = {PsilTemp_MPa[indsSilValidTemp, 0]:.1f}.')
                thisphiTop_frac = phiMax_frac / multphi_frac**nProfiles
                Planet.Sil.fn_phi_frac.update(thisphiTop_frac)
            else:
                rSilOuter_m = rSilTemp_m[indsSilValidTemp, 0]
                log.debug(f'Silicate match for Htidal = {thisHtidal_Wm3:.3e} W/m^3 with ' +
                          f'rSil = {rSilOuter_m / 1e3:.1f} km ({rSilOuter_m / Planet.Bulk.R_m:.3f} R_{Planet.name[0]}) ' +
                          f'at PHydroMax_MPa = {PsilTemp_MPa[indsSilValidTemp, 0]:.1f}.')
                thisHtidal_Wm3 = HtidalStart_Wm3 * multHtidal_Wm3**nProfiles
                Planet.Sil.fn_Htidal_Wm3 = GetHtidalFunc(thisHtidal_Wm3)  # Placeholder until we implement a self-consistent calc

            # Perform a check to avoid the final evaluation of SilicateLayers if we already know we won't use the outputs
            if (thisHtidal_Wm3 <= Planet.Sil.HtidalMax_Wm3 and thisphiTop_frac <= phiMax_frac):
                indsSilValidTemp, _, PsilTemp_MPa, TsilTemp_K, rSilTemp_m, rhoSilTemp_kgm3, MLayerSilTemp_kg, MAboveSilTemp_kg, gSilTemp_ms2, \
                phiSilTemp_frac, HtidalSilTemp_Wm3, kThermSilTemp_WmK, PporeTemp_MPa, rhoSilMatrixTemp_kgm3, rhoPoreTemp_kgm3, phaseSilPoreTemp \
                    = SilicateLayers(Planet, Params)
                # Record hydrosphere indices to include along with each profile
            iValid = np.append(iValid, indsSilValidTemp + Planet.Steps.iSilStart)

        # Mark all mass-matching profiles as valid
        indsSilValid = range(nProfiles)
        if nProfiles == np.size(iValid) - 1 or np.size(iValid) == 0:
            # Final profile is the one that violates the loop condition (or final value repeated if SilicateLayers is skipped),
            # so it is not valid -- pop it off
            iValid = iValid[:-1]
        # Now fill in values we're missing from not doing core calculations
        # We need copies of Steps.nSilMax for each possible result from the C/MR^2 search
        # to be compatible with the same infrastructure for when we have a core.
        nSilFinal = Planet.Steps.nSilMax * np.ones_like(indsSilValid).astype(np.int_)
        Ccore_kgm2 = 0

        C_kgm2 = np.zeros(Planet.Steps.nHydroMax)

    # Find contribution to axial moment of inertia C from each silicate layer
    dCfromSil_kgm2 = 8*np.pi/15 * rhoSil_kgm3 * (rSil_m[:,:-1]**5 - rSil_m[:,1:]**5)

    # Calculate C for a mantle extending up to each hydrosphere layer in turn
    Chydro_kgm2 = np.array([np.sum(dCfromH2O_kgm2[:i+1]) for i in iValid])
    Csil_kgm2 = np.array([np.sum(dCfromSil_kgm2[i,:nSilFinal[i]]) for i in indsSilValid])
    C_kgm2[iValid] = Chydro_kgm2 + Csil_kgm2 + Ccore_kgm2
    CMR2 = C_kgm2 / MR2_kgm2

    CMR2inds = [i[0] for i, valCMR2 in np.ndenumerate(CMR2)
                 if valCMR2 > Planet.Bulk.Cmeasured - Planet.Bulk.CuncertaintyLower
                and valCMR2 < Planet.Bulk.Cmeasured + Planet.Bulk.CuncertaintyUpper]

    if len(CMR2inds) == 0:
        if np.size(iValid) == 0 or ((np.min(CMR2[iValid]) < Planet.Bulk.Cmeasured) and (np.max(CMR2[iValid]) > Planet.Bulk.Cmeasured)):
            if Planet.Do.Fe_CORE:
                tweakable = 'increasing Steps.nSil or decreasing Ocean.deltaP'
            elif Planet.Do.POROUS_ROCK and not Planet.Do.FIXED_POROSITY:
                tweakable = 'increasing Steps.nPoros'
            else:
                tweakable = 'decreasing Sil.deltaHtidal'
            suggestion = f'\nTry increasing the resolution in C/MR^2 by {tweakable}.'
        elif(np.max(CMR2) > Planet.Bulk.Cmeasured):
            suggestion = '\nTry decreasing Tb_K or adjusting silicate/core composition to get lower C/MR^2 values.'
        else:
            suggestion = '\nTry adjusting properties of silicates and core to get higher C/MR^2 values.'

        if np.size(iValid) == 0:
            msg = 'No silicate mantle size was less than the total body mass for the initialization settings.'
        else:
            msg = f'No MoI found matching C/MR^2 = {Planet.CMR2strPrint}.\n' + \
                  f'Min: {np.min(CMR2[CMR2>0]):.4f}, Max: {np.max(CMR2):.4f}. '
        if Params.ALLOW_BROKEN_MODELS:
            if Params.DO_EXPLOREOGRAM:
                log.info(msg)
            else:
                log.error(msg + suggestion + ' Params.ALLOW_BROKEN_MODELS is True, so calculations will proceed ' +
                                             'with many values set to nan.')
            Planet.Do.VALID = False
            if Planet.invalidReason is None:
                Planet.invalidReason = ''
            else:
                Planet.invalidReason += '. '
            Planet.invalidReason += f'No valid MoI was found between {Planet.Bulk.Cmeasured - Planet.Bulk.CuncertaintyLower:.4f} and {Planet.Bulk.Cmeasured + Planet.Bulk.CuncertaintyUpper:.4f}'
        else:
            raise ValueError(msg + suggestion)

        nans = np.array([np.nan])
        Planet.CMR2mean = np.nan
        Planet.CMR2less = Planet.CMR2mean
        Planet.CMR2more = Planet.CMR2mean
        Planet.Sil.rhoTrade_kgm3 = nans
        Planet.Sil.Rmean_m = np.nan
        Planet.Sil.Rtrade_m = nans
        Planet.Sil.Rrange_m = np.nan
        Planet.Core.Rmean_m = np.nan
        Planet.Core.Rtrade_m = nans
        Planet.Core.Rrange_m = np.nan
        Planet.Steps.nSil = Planet.Steps.nSilMax
        iCMR2sil = 0
        nSilFinal = [-1]
        if Planet.Do.NO_H2O:
            Planet.Steps.nHydro = 0
        else:
            Planet.Steps.nHydro = Planet.Steps.nOceanMax

        Psil_MPa, Tsil_K, rSil_m, rhoSilEOS_kgm3, gSil_ms2, phiSil_frac, kThermSil_WmK, Ppore_MPa, PsilPore_MPa, \
        rhoSilMatrix_kgm3, rhoSil_kgm3, rhoSilPore_kgm3, HtidalSil_Wm3, MLayerSil_kg \
            = (np.zeros((1, Planet.Steps.nSil+1)) for _ in range(14))
        phaseSilPore = np.zeros((1, Planet.Steps.nSil+1), dtype=np.int_)
        coreProps = (np.zeros(Planet.Steps.nCore) for _ in range(9))
        Planet.Sil.phiCalc_frac = np.nan
        Planet.Sil.rhoMean_kgm3 = np.nan
        Planet.Core.rhoMean_kgm3 = np.nan
        Planet.Mtot_kg = Planet.Bulk.M_kg
        RcoreOrHtidalLine = ''

    else:
        # Find the C/MR^2 value most closely matching the measured value
        CMR2diff = np.abs(CMR2[CMR2inds] - Planet.Bulk.Cmeasured)
        # Get index of closest match in CMR2inds
        iCMR2ind = np.argmin(CMR2diff)
        # Check that we have only one best-fit option
        if np.size(iCMR2ind) > 1:
            log.warning('Multiple mass-matching profiles had the same MoI value. Using the first one. Try increasing the ' +
                     'resolution in the hydrosphere (by decreasing Ocean.deltaP) to avoid this problem.')
            iCMR2ind = iCMR2ind[0]
        # Find Planet array index corresponding to closest matching value
        iCMR2 = CMR2inds[iCMR2ind]
        # Get indices for silicate layer arrays
        if Planet.Do.Fe_CORE:
            iCMR2sil = iCMR2 - Planet.Steps.iSilStart
            CMR2indsSil = [ind - Planet.Steps.iSilStart for ind in CMR2inds]
        else:
            iCMR2sil = np.where(iCMR2 == iValid)[0][0]
            CMR2indsSil = [np.where(ind == iValid)[0][0] for ind in CMR2inds]
        # Record the best-match C/MR^2 value and neighboring values
        Planet.CMR2mean = CMR2[iCMR2]
        if Planet.Do.NO_H2O:
            log.warning('Only one profile is calculated for waterless bodies. Another method ' +
                        'of giving a trade space between total mass and MoI should replace this ' +
                        'implementation, varying surface heat flux, tidal heating, and/or porosity.')
            Planet.CMR2less = np.nan
            Planet.CMR2more = np.nan
        else:
            CMR2validSorted = np.sort(CMR2[CMR2inds])
            iSortedCMR2mean = np.where(CMR2validSorted == Planet.CMR2mean)[0][0]
            if iSortedCMR2mean == 0:
                Planet.CMR2less = Planet.CMR2mean
            else:
                Planet.CMR2less = CMR2validSorted[iSortedCMR2mean-1]
            if iSortedCMR2mean == np.size(CMR2validSorted) - 1:
                Planet.CMR2more = Planet.CMR2mean
            else:
                Planet.CMR2more = CMR2validSorted[iSortedCMR2mean+1]
        # Now we finally know how many layers there are in the hydrosphere and silicates
        Planet.Steps.nHydro = iCMR2
        Planet.Steps.nSil = nSilFinal[iCMR2sil]

        # Fill core/mantle trade arrays and set mean values consistent with MoI
        MtotSil_kg = np.sum(MLayerSil_kg[iCMR2sil,:nSilFinal[iCMR2sil]])
        VtotSil_m3 = 4/3*np.pi * (rSil_m[iCMR2sil,0]**3 - rSil_m[iCMR2sil,nSilFinal[iCMR2sil]]**3)
        Planet.Sil.rhoMean_kgm3 = MtotSil_kg / VtotSil_m3
        Planet.Sil.rhoTrade_kgm3 = np.array([np.sum(MLayerSil_kg[i,:nSilFinal[i]]) / (4/3*np.pi * (rSil_m[i,0]**3 - rSil_m[i,nSilFinal[i]-1]**3)) for i in CMR2indsSil])
        Planet.Sil.Rmean_m = Planet.r_m[iCMR2]
        Planet.Sil.Rtrade_m = Planet.r_m[CMR2inds]
        Planet.Sil.Rrange_m = np.max(Planet.Sil.Rtrade_m) - np.min(Planet.Sil.Rtrade_m)
        HtotSil_W = np.sum(HtidalSil_Wm3[iCMR2sil,:nSilFinal[iCMR2sil]] * 4/3*np.pi *
                           (rSil_m[iCMR2sil,:nSilFinal[iCMR2sil]]**3 - rSil_m[iCMR2sil,1:nSilFinal[iCMR2sil]+1]**3))
        Planet.Sil.Htidal_Wm3 = HtotSil_W / VtotSil_m3

        if Planet.Do.POROUS_ROCK and not Planet.Do.FIXED_POROSITY:
            if Planet.Do.Fe_CORE:
                Planet.Sil.phiCalc_frac = Planet.Sil.phiRockMax_frac
            else:
                Planet.Sil.phiCalc_frac = phiTop_frac[iCMR2sil]
        else:
            if Planet.Do.FIXED_POROSITY:
                Planet.Sil.phiCalc_frac = Planet.Sil.phiRockMax_frac
            else:
                Planet.Sil.phiCalc_frac = np.nan
        if Planet.Do.Fe_CORE:
            # Get indices for iron layer arrays
            iCMR2core = np.where(indsSilValid == iCMR2sil)[0][0]
            CMR2indsCore = [np.where(indsSilValid == i)[0][0] for i in CMR2indsSil]
            MtotCore_kg = np.sum(MLayerCore_kg[iCMR2core,:])
            Planet.Core.rhoMean_kgm3 = MtotCore_kg / (4/3*np.pi * rCore_m[iCMR2core,0]**3)
            Planet.Core.Rmean_m = rCore_m[iCMR2core,0]
            Planet.Core.Rtrade_m = rCore_m[CMR2indsCore,0]
            Planet.Core.Rrange_m = np.max(Planet.Core.Rtrade_m) - np.min(Planet.Core.Rtrade_m)
            RcoreOrHtidalLine = f'R_core = {Planet.Core.Rmean_m / Planet.Bulk.R_m:.2f} R, '

            # Package up core properties for returning
            coreProps = (Pcore_MPa[iCMR2core,:], Tcore_K[iCMR2core,:], rCore_m[iCMR2core,:-1],
                         gCore_ms2[iCMR2core, :], rhoCore_kgm3[iCMR2core,:], CpCore_JkgK[iCMR2core,:],
                         alphaCore_pK[iCMR2core,:], kThermCore_WmK[iCMR2core,:],
                         MLayerCore_kg[iCMR2core,:])
        else:
            if Planet.Do.POROUS_ROCK and not Planet.Do.FIXED_POROSITY:
                RcoreOrHtidalLine = f'phi_vac = {Planet.Sil.phiCalc_frac:.3f}, '
            else:
                RcoreOrHtidalLine = f'H_tidal = {Planet.Sil.Htidal_Wm3:.2e} W/m^3, '
            MtotCore_kg = 0
            Planet.Core.rhoMean_kgm3 = 0
            Planet.Core.Rmean_m = 0
            Planet.Core.Rtrade_m = np.zeros_like(Planet.Sil.Rtrade_m)
            Planet.Core.Rrange_m = 0
            coreProps = None

        Planet.Mtot_kg = np.sum(Planet.MLayer_kg[:iCMR2]) + MtotSil_kg + MtotCore_kg

    if not np.isnan(Planet.CMR2mean):
        log.info(f'Found matching MoI of {Planet.CMR2mean:.4f} ' +
                 f'(C/MR^2 = {Planet.CMR2strPrint}) for ' +
                 f'rho_sil = {Planet.Sil.rhoMean_kgm3:.0f} kg/m^3, ' +
                 f'R_sil = {Planet.Sil.Rmean_m / Planet.Bulk.R_m:.3f} R_{Planet.name[0]}, ' +
                 f'P_HydroMax = {Psil_MPa[iCMR2sil,0]:.1f} MPa, ' +
                 RcoreOrHtidalLine +
                 f'M_tot = {Planet.Mtot_kg/Planet.Bulk.M_kg:.4f} M_{Planet.name[0]}.')

    mantleProps = (Psil_MPa[iCMR2sil,:nSilFinal[iCMR2sil]], Tsil_K[iCMR2sil,:nSilFinal[iCMR2sil]],
                   rSil_m[iCMR2sil,:nSilFinal[iCMR2sil]], rhoSil_kgm3[iCMR2sil,:nSilFinal[iCMR2sil]],
                   gSil_ms2[iCMR2sil,:nSilFinal[iCMR2sil]], phiSil_frac[iCMR2sil,:nSilFinal[iCMR2sil]],
                   HtidalSil_Wm3[iCMR2sil,:nSilFinal[iCMR2sil]], kThermSil_WmK[iCMR2sil,:nSilFinal[iCMR2sil]],
                   PsilPore_MPa[iCMR2sil,:nSilFinal[iCMR2sil]], rhoSilMatrix_kgm3[iCMR2sil,:nSilFinal[iCMR2sil]],
                   rhoSilPore_kgm3[iCMR2sil,:nSilFinal[iCMR2sil]], MLayerSil_kg[iCMR2sil,:nSilFinal[iCMR2sil]],
                   phaseSilPore[iCMR2sil,:nSilFinal[iCMR2sil]])

    return Planet, mantleProps, coreProps
