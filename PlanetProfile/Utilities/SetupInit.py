""" Docstring explaining what we're doing here """

import os
import numpy as np
import logging
from datetime import datetime
from collections.abc import Iterable
from PlanetProfile import _ROOT
from PlanetProfile.GetConfig import Color, Style, FigLbl, FigMisc
from PlanetProfile.Thermodynamics.HydroEOS import GetOceanEOS, GetIceEOS, GetTfreeze
from PlanetProfile.Thermodynamics.InnerEOS import GetInnerEOS
from PlanetProfile.Thermodynamics.Reaktoro.CustomSolution import SetupCustomSolution, strip_latex_formatting_from_CustomSolutionLabel
from PlanetProfile.Thermodynamics.Clathrates.ClathrateProps import ClathDissoc
from PlanetProfile.Utilities.PPversion import ppVerNum, CheckCompat
from PlanetProfile.Utilities.defineStructs import DataFilesSubstruct, FigureFilesSubstruct, Constants
from PlanetProfile.TrajecAnalysis import _MAGdir, _scList
from PlanetProfile.TrajecAnalysis.FlybyEvents import scTargets
from PlanetProfile.TrajecAnalysis.RefileMAGdata import RefileName, MAGtoHDF5, LoadMAG

# Assign logger
log = logging.getLogger('PlanetProfile')

def SetupInit(Planet, Params):

    # Print version number
    log.debug(f'-- PlanetProfile {ppVerNum} --')
    if ppVerNum[-3:] == 'dev': log.debug('This version is in development.')
    # Print body name
    log.info(f'Body name: {Planet.name}')
    log.debug(f'Input file: {Planet.fname}')

    # Check dependency compatibility
    CheckCompat('seafreeze')  # SeaFreeze
    if Planet.Ocean.comp is not None and Planet.Ocean.comp == 'Seawater': CheckCompat('gsw')  # Gibbs Seawater
    if Planet.Do.TAUP_SEISMIC: CheckCompat('obspy')  # TauP (accessed as obspy.taup)
    if Params.CALC_NEW_INDUCT: CheckCompat('MoonMag')  # MoonMag
    if Params.CALC_NEW_GRAVITY: CheckCompat('pyalma3')

    # Check if Custom Reaktoro Solution is being used and if so then update Params with necessary parameters to plot
    if Planet.Ocean.comp is not None and 'CustomSolution' in Planet.Ocean.comp:
        CheckCompat('reaktoro')
        Planet, Params = SetupCustomSolution(Planet, Params)

    # Afford for additional MoI lower-bound uncertainty under non-hydrostatic conditions of 3% of C/MR^2,
    # in accordance with Gao and Stevenson (2013): https://doi.org/10.1016/j.icarus.2013.07.034
    if Planet.Bulk.CuncertaintyUpper is None:
        Planet.Bulk.CuncertaintyUpper = Planet.Bulk.Cuncertainty
    if Planet.Bulk.CuncertaintyLower is None:
        if Planet.Do.NONHYDROSTATIC:
            Planet.Bulk.CuncertaintyLower = Planet.Bulk.Cuncertainty + 0.03 * Planet.Bulk.Cmeasured
        else:
            Planet.Bulk.CuncertaintyLower = Planet.Bulk.Cuncertainty
    Planet = SetCMR2strings(Planet)
    # If we are not allowing ocean layers except for inner ices, we need to set NO_OCEAN to True
    if Planet.Do.NO_OCEAN_EXCEPT_INNER_ICES:
        Planet.Do.NO_OCEAN = True
    # Waterless, partially differentiated, and undifferentiated bodies. We have to do this first,
    # before filename generation, to ensure ocean comp is set.
    if Planet.Do.NO_H2O or Planet.Do.NO_DIFFERENTIATION or Planet.Do.PARTIAL_DIFFERENTIATION:
        Planet.Do.NO_OCEAN = True

        if Planet.Bulk.qSurf_Wm2 is None:
            raise ValueError('Bulk.qSurf_Wm2 must be set in order to model waterless or ' +
                             'partial/undifferentiated bodies.')
        Planet.Ocean.QfromMantle_W = Planet.Bulk.qSurf_Wm2 * 4*np.pi * Planet.Bulk.R_m**2
        Planet.qSurf_Wm2 = Planet.Bulk.qSurf_Wm2
        Planet.qCon_Wm2 = np.nan
        Planet.etaConv_Pas = np.nan
        Planet.Pb_MPa = Planet.Bulk.Psurf_MPa
        Planet.PbI_MPa = Planet.Bulk.Psurf_MPa
        Planet.Bulk.Tb_K = Planet.Bulk.Tsurf_K
        Planet.zb_km = 0.0
        Planet.Steps.nIceI = 0
        Planet.Steps.nIceIIILitho = 0
        Planet.Steps.nIceVLitho = 0
        Planet.Steps.nIbottom = 0
        Planet.Steps.nSurfIce = 0
        Planet.Steps.nOceanMax = 1
        Planet.Steps.nHydroMax = 1
        Planet.Steps.iSilStart = 0
        Planet.Ocean.deltaP = 0.1
        Planet.Do.NO_ICE_CONVECTION = True
        Planet.Do.CLATHRATE = False

        if Planet.Do.NO_H2O:
            log.info('Modeling a waterless body.')
            Planet.Ocean.comp = 'none'
            Planet.Ocean.wOcean_ppt = 0.0
            Planet.Sil.PHydroMax_MPa = Planet.Bulk.Psurf_MPa
            if Planet.Do.POROUS_ROCK:
                # Generate zero-yielding ocean "EOS" for use in porosity calculations
                # Note that there must be enough input points for creating spline
                # interpolators, even though we will not use them.
                Planet.Ocean.EOS = GetOceanEOS('none', None, np.linspace(0, 1, 10), np.linspace(0, 1, 10), None)

        else:
            if Planet.Do.NO_DIFFERENTIATION:
                log.info('Modeling an undifferentiated body.')
                Planet.Sil.Pclosure_MPa = Constants.PclosureUniform_MPa
            else:
                log.info('Modeling a partially differentiated body.')
                if Planet.Do.DIFFERENTIATE_VOLATILES:
                    raise ValueError('Do.DIFFERENTIATE_VOLATILES is not implemented yet.')

            Planet.Sil.PHydroMax_MPa = Planet.Ocean.PHydroMax_MPa
            # Make sure everything needed is set
            Planet.Do.POROUS_ROCK = True
            if Planet.Sil.phiRockMax_frac is None:
                raise ValueError('Sil.phiRockMax_frac must be set for partial/undifferentiated bodies.')
            if Planet.Do.Fe_CORE:
                raise ValueError('Do.Fe_CORE must be False for partial/undifferentiated bodies.')

            if Planet.Ocean.comp is None and Planet.Sil.poreComp is None:
                raise ValueError(f'Either Ocean.comp or Sil.poreComp must be set.')
            elif Planet.Ocean.comp is None:
                Planet.Ocean.comp = Planet.Sil.poreComp
            elif Planet.Sil.poreComp is None:
                Planet.Sil.poreComp = Planet.Ocean.comp
            else:
                if Planet.Sil.poreComp != Planet.Ocean.comp:
                    log.warning(f'Sil.poreComp ("{Planet.Sil.poreComp}") does not match ' +
                                f'Ocean.comp ("{Planet.Ocean.comp}"), but only pore fluids are ' +
                                f'modeled for partial/undifferentiated bodies. Sil.poreComp will ' +
                                f'be ignored.')
                Planet.Sil.poreComp = Planet.Ocean.comp

            if Planet.Ocean.wOcean_ppt is None and Planet.Sil.wPore_ppt is None:
                raise ValueError(f'Either Ocean.wOcean_ppt or Sil.wPore_ppt must be set.')
            elif Planet.Ocean.wOcean_ppt is None:
                Planet.Ocean.wOcean_ppt = Planet.Sil.wPore_ppt
            elif Planet.Sil.wPore_ppt is None:
                Planet.Sil.wPore_ppt = Planet.Ocean.wOcean_ppt
            else:
                if Planet.Sil.wPore_ppt != Planet.Ocean.wOcean_ppt:
                    log.warning(f'Sil.wPore_ppt ("{Planet.Sil.wPore_ppt}") does not match ' +
                                f'Ocean.wOcean_ppt ("{Planet.Ocean.wOcean_ppt}"), but only pore ' +
                                f'fluids are modeled for partial/undifferentiated bodies. ' +
                                f'Sil.wPore_ppt will be ignored.')
                Planet.Sil.wPore_ppt = Planet.Ocean.wOcean_ppt
    if Planet.Do.POROUS_ROCK:
        if Planet.Ocean.wOcean_ppt is None and Planet.Sil.wPore_ppt is None:
            raise ValueError(f'Either Ocean.wOcean_ppt or Sil.wPore_ppt must be set.')
        elif Planet.Ocean.wOcean_ppt is None:
            Planet.Ocean.wOcean_ppt = Planet.Sil.wPore_ppt
        elif Planet.Sil.wPore_ppt is None:
            Planet.Sil.wPore_ppt = Planet.Ocean.wOcean_ppt

    # Get filenames for saving/loading
    Planet, Params.DataFiles, Params.FigureFiles = SetupFilenames(Planet, Params)
    if Planet.Do.NON_SELF_CONSISTENT:
        SetupNonSelfConsistent(Planet, Params)
    else:
        # Set steps and settings for unused options to zero, check that we have settings we need
        # Core settings
        if Planet.Do.Fe_CORE:
            # (Re)set a predefined core radius, i.e. from CONSTANT_INNER_DENSITY, so that we can
            # check if it's None to see if we used that option.
            Planet.Core.Rset_m = None
        else:
            Planet.Steps.nCore = 0

        # Clathrates
        if Planet.Do.CLATHRATE:
            if Planet.Bulk.clathType is None:
                raise ValueError('Clathrate model type must be set. Options are "top", "bottom", and "whole".')
            elif Planet.Bulk.clathType == 'whole':
                # Pick whichever number of layers is greater, if the user has set Steps.nClath.
                # This allows the user to skip setting Steps.nClath separately if they just want
                # to model whole-shell clathrates with other standard run settings.
                if Planet.Steps.nClath is None: Planet.Steps.nClath = 0
                Planet.Steps.nClath = np.maximum(Planet.Steps.nIceI, Planet.Steps.nClath)
                Planet.Steps.nIceI = 0
            else:
                if Planet.Bulk.clathMaxThick_m is None:
                    raise ValueError('Bulk.clathMaxThick_m must be set for this clathType model.')
                elif Planet.Steps.nClath is None:
                    raise ValueError('Steps.nClath must be set for this clathType model.')
                elif Planet.Bulk.clathType == 'bottom':
                    if Planet.Bulk.qSurf_Wm2 is None:
                        raise ValueError('Bulk.qSurf_Wm2 must be set for this clathType model.')
                    if not Planet.Do.NO_ICE_CONVECTION:
                        log.warning('Do.NO_ICE_CONVECTION is False, but convection is incompatible with clathrate underplating. ' +
                                    'Do.NO_ICE_CONVECTION will be forced on.')
                        Planet.Do.NO_ICE_CONVECTION = True
            Planet.Ocean.ClathDissoc = ClathDissoc(Planet.Bulk.Tb_K, NAGASHIMA_CLATH_DISSOC=Planet.Do.NAGASHIMA_CLATH_DISSOC,
                                                ALLOW_BROKEN_MODELS=Params.ALLOW_BROKEN_MODELS,
                                                DO_EXPLOREOGRAM=Params.DO_EXPLOREOGRAM)
        else:
            Planet.Steps.nClath = 0
            Planet.zClath_m = 0
            Planet.Bulk.clathType = 'none'

        if not Planet.Do.NO_OCEAN:
            # In addition, perform some checks on underplating settings to be sure they make sense
            if not Planet.Do.BOTTOM_ICEIII and not Planet.Do.BOTTOM_ICEV:
                Planet.Steps.nIceIIILitho = 0
                Planet.Steps.nIceVLitho = 0
                Planet.RaConvectIII = np.nan
                Planet.RaConvectV = np.nan
                Planet.RaCritIII = np.nan
                Planet.RaCritV = np.nan
                Planet.eLidIII_m = 0
                Planet.eLidV_m = 0
                Planet.DconvIII_m = 0
                Planet.DconvV_m = 0
                Planet.deltaTBLIII_m = 0
                Planet.deltaTBLV_m = 0
            elif not Planet.Do.BOTTOM_ICEV:
                Planet.Steps.nIceVLitho = 0
                Planet.RaConvectV = np.nan
                Planet.RaCritV = np.nan
                Planet.eLidV_m = 0
                Planet.DconvV_m = 0
                Planet.deltaTBLV_m = 0
                if Planet.Ocean.PHydroMax_MPa < 209.9:
                    raise ValueError('Hydrosphere max pressure is less than the pressure of the ice I-III phase transition, ' +
                                    'but Do.BOTTOM_ICEIII is True.')
                if Planet.Bulk.Tb_K > Planet.Bulk.TbIII_K:
                    log.warning('Bottom temperature of ice I (Tb_K) is greater than bottom temperature of underplate ' +
                                'ice III (TbIII_K). This likely represents a non-equilibrium state.')
            else:
                if Planet.Ocean.PHydroMax_MPa < 344.3:
                    raise ValueError('Hydrosphere max pressure is less than the pressure of the ice III-V phase transition, ' +
                                    'but Do.BOTTOM_ICEV is True.')
                if Planet.Bulk.Tb_K > Planet.Bulk.TbIII_K:
                    log.warning('Bottom temperature of ice I (Tb_K) is greater than bottom temperature of underplate ' +
                                'ice III (TbIII_K). This likely represents a non-equilibrium state.')
                if Planet.Bulk.TbIII_K > Planet.Bulk.TbV_K:
                    log.warning('Bottom temperature of ice III (Tb_K) is greater than bottom temperature of underplate ' +
                                'ice V (TbV_K). This likely represents a non-equilibrium state.')
                if Planet.Do.CLATHRATE:
                    log.warning('Clathrates are stable under a very large range of pressures and temperatures, and this ' +
                                'may be contradictory with having underplating ice III or V.')
      
            # Make sure ocean max temp is above melting temp
            if Planet.Ocean.THydroMax_K <= Planet.Bulk.Tb_K:
                Planet.Ocean.THydroMax_K = Planet.Bulk.Tb_K + 30

        # Get ocean EOS functions
        POcean_MPa = np.arange(Planet.PfreezeLower_MPa, Planet.Ocean.PHydroMax_MPa, Planet.Ocean.deltaP)
        # Set Ocean.deltaT to use default of 0.1 K but recommend the user to set
        if Planet.Ocean.deltaT is None:
            Planet.Ocean.deltaT = 1e-1
            log.warning('Ocean.deltaT is not set--defaulting to 0.1 K. This may not be precise enough ' +
                        'for shallow oceans or fine control over ice shell thickness calculations. ' +
                        'It is recommended to set Ocean.deltaT manually in the PPBody.py file.')
        # Check ocean parameter space and load EOS
        if Planet.Ocean.THydroMax_K < Planet.Bulk.Tb_K:
            raise ValueError(f'Ocean.THydroMax_K of {Planet.Ocean.THydroMax_K} is less than Bulk.Tb_K of {Planet.Bulk.Tb_K}.')
        elif Planet.Do.ICEIh_THICKNESS:
            TOcean_K = np.arange(Planet.TfreezeLower_K, Planet.Ocean.THydroMax_K, Planet.Ocean.deltaT)
        else:
            TOcean_K = np.arange(Planet.Bulk.Tb_K, Planet.Ocean.THydroMax_K, Planet.Ocean.deltaT)
        Planet.Ocean.EOS = GetOceanEOS(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, POcean_MPa, TOcean_K,
                                       Planet.Ocean.MgSO4elecType, rhoType=Planet.Ocean.MgSO4rhoType,
                                       scalingType=Planet.Ocean.MgSO4scalingType, FORCE_NEW=Params.FORCE_EOS_RECALC,
                                       phaseType=Planet.Ocean.phaseType, EXTRAP=Params.EXTRAP_OCEAN,
                                       sigmaFixed_Sm=Planet.Ocean.sigmaFixed_Sm, LOOKUP_HIRES=Planet.Do.OCEAN_PHASE_HIRES)
        if Planet.Ocean.EOS.deltaP != Planet.Ocean.deltaP:
            log.debug(f'Updating Ocean.deltaP to match the more refined EOS.deltaP ({Planet.Ocean.EOS.deltaP}).')
            Planet.Ocean.deltaP = Planet.Ocean.EOS.deltaP
        if Planet.Ocean.EOS.deltaT != Planet.Ocean.deltaT:
            log.debug(f'Updating Ocean.deltaT to match the more refined EOS.deltaT ({Planet.Ocean.EOS.deltaT}).')
            Planet.Ocean.deltaT = Planet.Ocean.EOS.deltaT
        
        # Get separate, simpler EOS for evaluating the melting curve
        if Planet.Do.ICEIh_THICKNESS:
            Tmelt_K = np.arange(Planet.TfreezeLower_K, Planet.TfreezeUpper_K, Planet.TfreezeRes_K)
        else:
            Planet.Bulk.zb_approximate_km = np.nan
            if Planet.Do.BOTTOM_ICEV:
                # Make sure Tb values are physically reasonable
                if Planet.Bulk.Tb_K > Planet.Bulk.TbV_K or Planet.Bulk.Tb_K > Planet.Bulk.TbIII_K:
                    raise ValueError('Bulk.Tb_K must be less than underplate layer Tb values.')
                if Planet.Bulk.TbIII_K > Planet.Bulk.TbV_K:
                    raise ValueError('Bulk.TbIII_K must be less than underplate TbV_K value.')
                TmeltI_K = np.linspace(Planet.Bulk.Tb_K - 0.01, Planet.Bulk.Tb_K + 0.01, 11)
                TmeltIII_K = np.linspace(Planet.Bulk.TbIII_K - 0.01, Planet.Bulk.TbIII_K + 0.01, 11)
                TmeltV_K = np.linspace(Planet.Bulk.TbV_K - 0.01, Planet.Bulk.TbV_K + 0.01, 11)
                Tmelt_K = np.concatenate((TmeltI_K, TmeltIII_K, TmeltV_K))
                if Planet.PfreezeUpper_MPa < 700:
                    log.info(f'PfreezeUpper_MPa is set to {Planet.PfreezeUpper_MPa}, but Do.BOTTOM_ICEV is True and ' +
                             'ice V is stable up to 700 MPa. PfreezeUpper_MPa will be raised to this value.')
                    Planet.PfreezeUpper_MPa = 700
            elif Planet.Do.BOTTOM_ICEIII:
                if Planet.Bulk.Tb_K > Planet.Bulk.TbIII_K:
                    raise ValueError('Bulk.Tb_K must be less than underplate TbIII_K value.')
                TmeltI_K = np.linspace(Planet.Bulk.Tb_K - 0.01, Planet.Bulk.Tb_K + 0.01, 11)
                TmeltIII_K = np.linspace(Planet.Bulk.TbIII_K - 0.01, Planet.Bulk.TbIII_K + 0.01, 11)
                Tmelt_K = np.concatenate((TmeltI_K, TmeltIII_K))
                if Planet.PfreezeUpper_MPa < 700:
                    log.info(f'PfreezeUpper_MPa is set to {Planet.PfreezeUpper_MPa}, but Do.BOTTOM_ICEIII is True and ' +
                             'ice V is stable up to 400 MPa. PfreezeUpper_MPa will be raised to this value.')
                    Planet.PfreezeUpper_MPa = 400
            else:
                Tmelt_K = np.linspace(Planet.Bulk.Tb_K - 0.01, Planet.Bulk.Tb_K + 0.01, 11)
        Pmelt_MPa = np.arange(Planet.PfreezeLower_MPa, Planet.PfreezeUpper_MPa, Planet.PfreezeRes_MPa)
        Planet.Ocean.meltEOS = GetOceanEOS(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, Pmelt_MPa, Tmelt_K,  None,
                                           phaseType=Planet.Ocean.phaseType, FORCE_NEW=True, MELT=True,
                                           LOOKUP_HIRES=Planet.Do.OCEAN_PHASE_HIRES)

        # Make sure convection checking outputs are set if we won't be modeling them
        if Planet.Do.NO_ICE_CONVECTION:
            Planet.RaConvect = np.nan
            Planet.RaConvectIII = np.nan
            Planet.RaConvectV = np.nan
            Planet.RaCrit = np.nan
            Planet.RaCritIII = np.nan
            Planet.RaCritV = np.nan
            Planet.Tconv_K = np.nan
            Planet.TconvIII_K = np.nan
            Planet.TconvV_K = np.nan
            Planet.eLid_m = 0.0
            Planet.eLidIII_m = 0.0
            Planet.eLidV_m = 0.0
            Planet.Dconv_m = 0.0
            Planet.DconvIII_m = 0.0
            Planet.DconvV_m = 0.0
            Planet.deltaTBL_m = 0.0
            Planet.deltaTBLIII_m = 0.0
            Planet.deltaTBLV_m = 0.0

        # Porous rock
        if Planet.Do.POROUS_ROCK:
            # Make sure pore max pressure is set, since we will need it to check if we'll need HP ices
            if not Planet.Do.PORE_EOS_DIFFERENT and Planet.Sil.PHydroMax_MPa is None:
                Planet.Sil.PHydroMax_MPa = Planet.Ocean.PHydroMax_MPa

            if Planet.Sil.porosType == 'Han2014':
                if Planet.Sil.phiRockMax_frac is None:
                    log.warning('Sil.phiRockMax_frac is not set. Using arbitrary max porosity of 0.7.')
                    Planet.Sil.phiRockMax_frac = 0.7
            else:
                Planet.Do.FIXED_POROSITY = True
            if Planet.Sil.phiRangeMult <= 1:
                raise ValueError(f'Sil.phiRangeMult = {Planet.Sil.phiRangeMult}, but it must be greater than 1.')
            # Ensure we set pore comp to ocean comp if not specified as different
            if Planet.Sil.poreComp is None:
                Planet.Sil.poreComp = Planet.Ocean.comp
                Planet.Do.PORE_EOS_DIFFERENT = False
            else:
                Planet.Do.PORE_EOS_DIFFERENT = True
            if Planet.Sil.wPore_ppt is None:
                Planet.Sil.wPore_ppt = Planet.Ocean.wOcean_ppt
                Planet.Do.PORE_EOS_DIFFERENT = False
            else:
                Planet.Do.PORE_EOS_DIFFERENT = True
        else:
            if (not Planet.Do.Fe_CORE) and (not Planet.Do.CONSTANT_INNER_DENSITY):
                raise RuntimeError('Matching the body MoI requires either a core or porosity in the rock.' +
                                'Set Planet.Do.POROUS_ROCK to True and rerun to continue modeling with no core.')
            Planet.Sil.porosType = 'none'
            Planet.Sil.poreH2Orho_kgm3 = 0
            Planet.Sil.phiRockMax_frac = 0
            Planet.Sil.poreComp = Planet.Ocean.comp
            Planet.Sil.wPore_ppt = Planet.Ocean.wOcean_ppt

        # Porous ice
        if not Planet.Do.POROUS_ICE:
            Planet.Ocean.phiMax_frac = {key:0 for key in Planet.Ocean.phiMax_frac.keys()}
            Planet.Ocean.porosType = {key:None for key in Planet.Ocean.porosType.keys()}

        # Effective pressure in pore space
        if not Planet.Do.P_EFFECTIVE:
            # Peffective is calculated from Pmatrix - alpha*Ppore, so setting alpha to zero avoids the need for repeated
            # conditional checks during layer propagation -- calculations are typically faster than conditional checks.
            if Planet.Sil.alphaPeff != 0:
                log.debug('Sil.alphaPeff was not 0, but is being set to 0 because Do.P_EFFECTIVE is False.')
            Planet.Sil.alphaPeff = 0
            for phase in Planet.Ocean.alphaPeff.keys():
                if Planet.Ocean.alphaPeff[phase] != 0:
                    log.debug(f'Ocean.alphaPeff[{phase}] was not 0, but is being set to 0 because Do.P_EFFECTIVE is False.')
                Planet.Ocean.alphaPeff[phase] = 0

        # Calculate bulk density from total mass and radius, and warn user if they specified density
        if Planet.Bulk.M_kg is None:
            Planet.Bulk.M_kg = Planet.Bulk.rho_kgm3 * (4/3*np.pi * Planet.Bulk.R_m**3)
        else:
            if Planet.Bulk.rho_kgm3 is not None and not Params.DO_EXPLOREOGRAM:
                log.warning('Both bulk mass and density were specified. Only one is required--' +
                            'density will be recalculated from bulk mass for consistency.')
            Planet.Bulk.rho_kgm3 = Planet.Bulk.M_kg / (4/3*np.pi * Planet.Bulk.R_m**3)

        # Load EOS functions for deeper interior
        if not Params.SKIP_INNER:
            # Get silicate EOS
            Planet.Sil.EOS = GetInnerEOS(Planet.Sil.mantleEOS, EOSinterpMethod=Params.lookupInterpMethod,
                                        kThermConst_WmK=Planet.Sil.kTherm_WmK, HtidalConst_Wm3=Planet.Sil.Htidal_Wm3,
                                        porosType=Planet.Sil.porosType, phiTop_frac=Planet.Sil.phiRockMax_frac,
                                        Pclosure_MPa=Planet.Sil.Pclosure_MPa, phiMin_frac=Planet.Sil.phiMin_frac,
                                        EXTRAP=Params.EXTRAP_SIL, etaSilFixed_Pas=Planet.Sil.etaRock_Pas, etaCoreFixed_Pas=[Planet.Core.etaFeSolid_Pas, Planet.Core.etaFeLiquid_Pas])

            # Pore fluids if present
            if Planet.Do.POROUS_ROCK:
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
                                                sigmaFixed_Sm=Planet.Sil.sigmaPoreFixed_Sm, LOOKUP_HIRES=Planet.Do.OCEAN_PHASE_HIRES, kThermConst_WmK=Planet.Ocean.kThermWater_WmK)

                if Planet.Do.NO_DIFFERENTIATION or Planet.Do.PARTIAL_DIFFERENTIATION:
                    Planet.Ocean.EOS = Planet.Sil.poreEOS
                    for icePhase in ['Ih', 'II', 'III', 'V', 'VI']:
                        Planet.Ocean.surfIceEOS[icePhase] = GetIceEOS(Ppore_MPa, Tpore_K, icePhase,
                                                                    EXTRAP=Params.EXTRAP_ICE[icePhase],
                                                                    ICEIh_DIFFERENT=Planet.Do.ICEIh_DIFFERENT, kThermConst_WmK=Planet.Ocean.kThermIce_WmK)

                # Make sure Sil.phiRockMax_frac is set in case we're using a porosType that doesn't require it
                if Planet.Sil.phiRockMax_frac is None or Planet.Sil.porosType != 'Han2014':
                    Planet.Sil.phiRockMax_frac = Planet.Sil.EOS.fn_phi_frac(0, 0)
                if Planet.Sil.phiRockMax_frac > 1.0:
                    raise ValueError(f'A maximum rock porosity value of {Planet.Sil.phiRockMax_frac}' +
                                    'has been set. Values greater than 1 are unphysical--check input' +
                                    'file settings.')

            # Iron core if present
            if Planet.Do.Fe_CORE:
                Planet.Core.EOS = GetInnerEOS(Planet.Core.coreEOS, EOSinterpMethod=Params.lookupInterpMethod, Fe_EOS=True,
                                            kThermConst_WmK=Planet.Core.kTherm_WmK, EXTRAP=Params.EXTRAP_Fe,
                                            wFeCore_ppt=Planet.Core.wFe_ppt, wScore_ppt=Planet.Core.wS_ppt, etaSilFixed_Pas=Planet.Sil.etaRock_Pas, etaCoreFixed_Pas=[Planet.Core.etaFeSolid_Pas, Planet.Core.etaFeLiquid_Pas])

        # Ensure ionosphere bounds and conductivity are in a format we expect
        if Planet.Magnetic.ionosBounds_m is None:
            Planet.Magnetic.ionosBounds_m = [np.nan]
        elif not isinstance(Planet.Magnetic.ionosBounds_m, Iterable):
            Planet.Magnetic.ionosBounds_m = [Planet.Magnetic.ionosBounds_m]
        if Planet.Magnetic.sigmaIonosPedersen_Sm is None:
            Planet.Magnetic.sigmaIonosPedersen_Sm = [np.nan]
        elif not isinstance(Planet.Magnetic.sigmaIonosPedersen_Sm, Iterable):
            Planet.Magnetic.sigmaIonosPedersen_Sm = [Planet.Magnetic.sigmaIonosPedersen_Sm]

    # Preallocate layer physical quantity arrays
    Planet = SetupLayers(Planet)

    return Planet, Params


def SetupInversion(Params):

    # Ensure configPP/Params toggles are set as needed
    Params.SKIP_INDUCTION = False
    Params.CALC_NEW = True
    Params.CALC_CONDUCT = True
    Params.CALC_SEISMIC = False

    # Check if relevant spacecraft data has been refiled into HDF5, and create if not
    magFnames = {}
    magData = {}
    for scName in Params.Trajec.scSelect:
        if Params.Trajec.targetBody not in scTargets[scName]:
            log.warning(f'Target {Params.Trajec.targetBody} not listed in scTargets for {scName}:' +
                        f' {scTargets[scName]}.')

        # Make sure data for this spacecraft will be found
        if scName not in _scList:
            log.warning(f'"{scName}" not found in {_MAGdir}.')

        # Make sure files are present
        magFnames[scName] = RefileName(Params.Trajec.targetBody, scName,
                                       MAGdir=Params.Trajec.MAGdir)
        if Params.Trajec.FORCE_MAG_RECALC or not os.path.isfile(magFnames[scName]):
            MAGtoHDF5(Params, scName)

        # Load magnetometer data
        magData[scName] = LoadMAG(Params, magFnames[scName], scName)

    # Fill TrajecParams from selections
    Params.Trajec.flybys = {scName: [fbID for fbID in data.t_UTC.keys()]
                            for scName, data in magData.items()}
    Params.Trajec.flybyNames = {scName: {fbID: f'{Params.Trajec.fbDescrip[scName]}{fbID}'
                                         for fbID in fbList} for scName, fbList in
                                Params.Trajec.flybys.items()}
    Params.Trajec.flybyFnames = {scName: {fbID: fbName.replace(' ', '')
                                          for fbID, fbName in fbList.items()} for scName, fbList in
                                 Params.Trajec.flybyNames.items()}
    Params.Trajec.nFitParams = Params.Trajec.nFitParamsGlobal + Params.Trajec.nFitParamsFlybys * \
                               np.sum([data.nFlybys for data in magData.values()])

    datPath = Params.Trajec.targetBody
    if datPath == 'Test':
        datPath = os.path.join(_ROOT, datPath)
    figPath = os.path.join(datPath, 'figures')
    Params.Trajec.FigureFiles = FigureFilesSubstruct(figPath, '', FigMisc.xtn,
                                                     flybys=Params.Trajec.flybyFnames)
    scStr = ''.join([f'{scName}{np.size(fbList)}'
                     for scName, fbList in Params.Trajec.flybys.items()])
    fitFbase = f'{Params.Trajec.targetBody}Fit{scStr}Trajec'
    Params.Trajec.DataFiles = DataFilesSubstruct(datPath, fitFbase, None,
                                                 inductAppend=Params.Trajec.trajecAppend)

    return Params, magData


def SetupFilenames(Planet, Params, exploreAppend=None, figExploreAppend=None, monteCarloAppend=None):
    """ Generate filenames for saving data and figures.
    """
    datPath = Planet.bodyname
    if Planet.bodyname == 'Test':
        datPath = os.path.join(_ROOT, datPath)
    figPath = os.path.join(datPath, 'figures')

    # Account for wppt being None or less than or equal to zero for CustomSolution, which in this case we need to set Do parameter so we don't have Nonetype error when generating filenames
    if Planet.Ocean.wOcean_ppt is None or Planet.Ocean.wOcean_ppt < 0:
        # Flag that we are not using wOcean_ppt as independent parameter
        Planet.Do.USE_WOCEAN_PPT = False
    # Account for differing ocean/pore composition here, since we need it for filenames
    # Use ocean composition and salinity if user has not specified different ones for pore space
    Planet.Do.PORE_EOS_DIFFERENT = False
    if Planet.Sil.poreComp is not None and Planet.Sil.poreComp != Planet.Ocean.comp:
        Planet.Do.PORE_EOS_DIFFERENT = True
    elif Planet.Sil.wPore_ppt is not None and Planet.Sil.wPore_ppt != Planet.Ocean.wOcean_ppt:
        Planet.Do.PORE_EOS_DIFFERENT = True

        
    saveBase = Planet.name + 'Profile_'
    saveLabel = ''
    label = ''
    comp = Planet.Ocean.comp
    if Planet.Do.NON_SELF_CONSISTENT:
        saveLabel = 'NonSelfConsistent_'
        setStr = 'Non-self-consistent'
    else:
        if Planet.Do.NO_H2O:
            saveLabel += f'NoH2O_qSurf{Planet.Bulk.qSurf_Wm2*1e3:.1f}mWm2'
            setStr = f'$q_\mathrm{{surf}}\,{Planet.Bulk.qSurf_Wm2*FigLbl.qMult:.1f}\,\si{{{FigLbl.fluxUnits}}}$'
            label += setStr
            Planet.compStr = r'No~\ce{H2O}'

        elif Planet.Do.NO_DIFFERENTIATION:
            saveLabel += f'NoDiff_qSurf{Planet.Bulk.qSurf_Wm2*1e3:.1f}mWm2'
            setStr = f'$q_\mathrm{{surf}}\,{Planet.Bulk.qSurf_Wm2*FigLbl.qMult:.1f}\,\si{{{FigLbl.fluxUnits}}}$'
            label += setStr
            Planet.compStr = r'Undifferentiated'

        elif Planet.Do.PARTIAL_DIFFERENTIATION:
            saveLabel += f'PartialDiff_qSurf{Planet.Bulk.qSurf_Wm2*1e3:.1f}mWm2'
            setStr = f'$q_\mathrm{{surf}}\,{Planet.Bulk.qSurf_Wm2*FigLbl.qMult:.1f}\,\si{{{FigLbl.fluxUnits}}}$'
            label += setStr

            if Planet.Sil.poreComp == 'PureH2O':
                Planet.compStr = r'Pure~\ce{H2O}'
                saveLabel += f'_{Planet.Sil.poreComp}Pores'
            else:
                Planet.compStr = f'${Planet.Sil.wPore_ppt*FigLbl.wMult:.1f}\,\si{{{FigLbl.wUnits}}}$~\ce{{{Planet.Sil.poreComp}}}'
                saveLabel += f'_{Planet.Sil.poreComp}_{Planet.Sil.wPore_ppt:.1f}pptPores'

            saveLabel += f'_PorousRock_phi{Planet.Sil.phiRockMax_frac:.2f}_Pc{Planet.Sil.Pclosure_MPa:5.2e}'
            label += ' w/$\phi_\mathrm{sil}$'
        
        else:
            if Planet.Do.ICEIh_THICKNESS:
                saveLabelAppendage = f'zb{Planet.Bulk.zb_approximate_km}km'
                labelAppendage = f'$z_b\,\SI{{{Planet.Bulk.zb_approximate_km}}}{{\kilo\meter}}$'
            else:
                saveLabelAppendage = f'Tb{Planet.Bulk.Tb_K}K'
                labelAppendage = f'$T_b\,\SI{{{Planet.Bulk.Tb_K}}}{{K}}$'
            if Planet.Ocean.comp == 'PureH2O':
                saveLabel += f'{Planet.Ocean.comp}_{saveLabelAppendage}'
                setStr = f'Pure \ce{{H2O}}'
                label += f'{setStr}, {labelAppendage}'
                Planet.compStr = r'Pure~\ce{H2O}'
            elif "CustomSolution" in Planet.Ocean.comp:
                # Get text to left of = sign
                CustomSolutionLabel = Planet.Ocean.comp.split('=')[0].strip()
                # Get label for plotting
                setStr = CustomSolutionLabel.replace("CustomSolution", "")
                # Strip latex formatting for save label
                saveCustomSolutionLabel = strip_latex_formatting_from_CustomSolutionLabel(CustomSolutionLabel)
                comp = CustomSolutionLabel
                # In this case, we are using input speciation from user with no input w_ppt, so we will generate filenames without w_ppt
                if not Planet.Do.USE_WOCEAN_PPT:
                    saveLabel += f'{saveCustomSolutionLabel}_{saveLabelAppendage}'
                    #label = f'{setStr}, {labelAppendage}'
                    label = f'{setStr}'
                    Planet.compStr = f'{CustomSolutionLabel}'
                else:
                    saveLabel += f'{saveCustomSolutionLabel}_{Planet.Ocean.wOcean_ppt:.1f}ppt' + \
                            f'_{saveLabelAppendage}'
                    label = f'{Planet.Ocean.comp}_{Planet.Ocean.wOcean_ppt:.1f}ppt' + \
                            f'_{labelAppendage}'
                    Planet.compStr = f'${Planet.Ocean.wOcean_ppt*FigLbl.wMult:.1f}\,\si{{{FigLbl.wUnits}}}${CustomSolutionLabel}'
            else:
                saveLabel += f'{Planet.Ocean.comp}_{Planet.Ocean.wOcean_ppt:.1f}ppt' + \
                            f'_{saveLabelAppendage}'
                setStr = f'${Planet.Ocean.wOcean_ppt*FigLbl.wMult:.1f}\,\si{{{FigLbl.wUnits}}}$ \ce{{{Planet.Ocean.comp}}}'
                label += setStr + \
                    f', {labelAppendage}'
                Planet.compStr = f'${Planet.Ocean.wOcean_ppt*FigLbl.wMult:.1f}\,\si{{{FigLbl.wUnits}}}$~\ce{{{Planet.Ocean.comp}}}'
            if Planet.Do.CLATHRATE:
                if Planet.Do.MIXED_CLATHRATE_ICE:
                    saveLabel += '_MixedClathrates'
                    label += ' w/mixed clathrates'
                else:
                    saveLabel += '_Clathrates'
                    label += ' w/clath'
            if Planet.Do.POROUS_ICE:
                if Planet.Do.CLATHRATE and Planet.Bulk.clathType != 'bottom':
                    saveLabel += f'_PorousIce_phi{Planet.Ocean.phiMax_frac["Clath"]:.2f}_Pc{Planet.Ocean.Pclosure_MPa["Clath"]:5.2e}'
                else:
                    saveLabel += f'_PorousIce_phi{Planet.Ocean.phiMax_frac["Ih"]:.2f}_Pc{Planet.Ocean.Pclosure_MPa["Ih"]:5.2e}'
                label += ' w/$\phi_\mathrm{ice}$'
            if Planet.Do.PORE_EOS_DIFFERENT:
                if Planet.Sil.poreComp == 'PureH2O':
                    saveLabel += f'_{Planet.Sil.poreComp}Pores'
                    label += f'Pure \ce{{H2O}} pores'
                else:
                    saveLabel += f'_{Planet.Sil.poreComp}_{Planet.Sil.wPore_ppt:.1f}pptPores'
                    label += f', {Planet.Sil.wPore_ppt*FigLbl.wMult:.1f}\,\si{{{FigLbl.wUnits}}} \ce{{{Planet.Sil.poreComp}}} pores'
            elif Planet.Do.POROUS_ROCK:
                saveLabel += f'_PorousRock_phi{Planet.Sil.phiRockMax_frac:.2f}_Pc{Planet.Sil.Pclosure_MPa:5.2e}'
                label += ' w/$\phi_\mathrm{sil}$'
            if Planet.Do.HYDROSPHERE_THICKNESS:
                saveLabel += f'_{Planet.Bulk.Dhsphere_m}'
                label += f'$hydroThickness\,\SI{{{Planet.Bulk.Dhsphere_m}}}{{\meter}}$'
        if Planet.Sil.mantleEOSName is not None: saveLabel += f'_{Planet.Sil.mantleEOSName}'

    # Add time and date label
    if Params.TIME_AND_DATE_LABEL:
        current_datetime = datetime.now()
        formatted_datetime = current_datetime.strftime("%Y-%m-%d_%H-%M")
        saveLabel += f"_{formatted_datetime}"
    Planet.saveLabel = saveLabel
    Planet.tradeLabel = f'{label}, $C/MR^2\,{Planet.CMR2str}$'
    Planet.label = label
    if Params.DO_INDUCTOGRAM:
        inductBase = f'{Planet.name}_{Params.Induct.inductOtype}'
        Params.Induct.SetFlabel(Planet.bodyname)
        inductAppend = Params.Induct.fLabel
        exploreBase = None
    else:
        inductBase = None
        inductAppend = None
        if Params.DO_EXPLOREOGRAM:
            exploreBase = f'{Planet.name}ExploreOgram_{exploreAppend}_{saveLabel}'
        else:
            exploreBase = None
            if Params.DO_MONTECARLO:
                monteCarloBase = f'{Planet.name}MonteCarlo_{monteCarloAppend}_{saveLabel}'
            else:
                monteCarloBase = None
    
    DataFiles = DataFilesSubstruct(datPath, saveBase + saveLabel, comp, inductBase=inductBase,
                                   exploreAppend=exploreAppend, EXPLORE=(Params.DO_INDUCTOGRAM or
                                       Params.DO_EXPLOREOGRAM or Params.DO_MONTECARLO or Params.INDUCTOGRAM_IN_PROGRESS),
                                   inductAppend=inductAppend, monteCarloAppend=monteCarloAppend)
    FigureFiles = FigureFilesSubstruct(figPath, saveBase + saveLabel, FigMisc.xtn,
                                       comp=comp, exploreBase=exploreBase, inductBase=inductBase, monteCarloBase=monteCarloBase,
                                       exploreAppend=figExploreAppend, inductAppend=inductAppend, monteCarloAppend=monteCarloAppend)

    # Attach profile name to PlanetStruct in addition to Params
    Planet.saveFile = DataFiles.saveFile + ''

    return Planet, DataFiles, FigureFiles


def SetupLayers(Planet):
    """ Initialize layer arrays in Planet.
    """

    if not Planet.Do.NO_H2O:
        if not Planet.Do.NON_SELF_CONSISTENT:
            nOceanMax = int(Planet.Ocean.PHydroMax_MPa / Planet.Ocean.deltaP)
            Planet.Steps.nHydroMax = Planet.Steps.nClath + Planet.Steps.nIceI + Planet.Steps.nIceIIILitho + Planet.Steps.nIceVLitho + nOceanMax
            nStepsForArrays = Planet.Steps.nHydroMax
        else:
            nStepsForArrays = Planet.Steps.nTotal
    
        Planet.phase = np.zeros(nStepsForArrays, dtype=np.int_)
        Planet.P_MPa, Planet.T_K, Planet.r_m, Planet.rho_kgm3, \
            Planet.Cp_JkgK, Planet.alpha_pK, Planet.g_ms2, Planet.phi_frac, \
            Planet.sigma_Sm, Planet.z_m, Planet.MLayer_kg, Planet.VLayer_m3, Planet.kTherm_WmK, \
            Planet.Htidal_Wm3, Planet.Ppore_MPa, Planet.rhoMatrix_kgm3, Planet.rhoPore_kgm3 = \
            (np.zeros(nStepsForArrays) for _ in range(17))

    # Layer property initialization for surface
    Planet.z_m[0] = 0.0  # Set first layer depth to zero (layer properties correspond to outer radius)
    Planet.r_m[0] = Planet.Bulk.R_m  # Set first layer to planetary surface radius
    Planet.g_ms2[0] = Constants.G * Planet.Bulk.M_kg / Planet.Bulk.R_m**2  # Set first layer gravity at surface
    Planet.T_K[0] = Planet.Bulk.Tsurf_K  # Set first layer surface temp
    Planet.P_MPa[0] = Planet.Bulk.Psurf_MPa  # Set first layer to surface pressure

    return Planet


def SetCMR2strings(Planet):
    # Create strings to describe input MoI in terminal outputs and plots/tables
    if Planet.Bulk.CuncertaintyLower == Planet.Bulk.CuncertaintyUpper:
        Planet.CMR2strPrint = f'{Planet.Bulk.Cmeasured:.4f} +/- {Planet.Bulk.CuncertaintyLower:.4f}'
        Planet.CMR2str = f'{Planet.Bulk.Cmeasured}\pm{Planet.Bulk.CuncertaintyLower}'
        Planet.CMR2str5 = f'{Planet.Bulk.Cmeasured:.5f}\pm{Planet.Bulk.CuncertaintyLower:.5f}'
    else:
        Planet.CMR2strPrint = f'{Planet.Bulk.Cmeasured:.4f} + {Planet.Bulk.CuncertaintyUpper:.4f} - {Planet.Bulk.CuncertaintyLower:.4f}'
        if round(Planet.Bulk.Cmeasured, 3) == Planet.Bulk.Cmeasured and round(Planet.Bulk.CuncertaintyUpper, 3) == Planet.Bulk.CuncertaintyUpper:
            CmidStr = f'{Planet.Bulk.Cmeasured:.3f}'
            CpStr = f'+{Planet.Bulk.CuncertaintyUpper:.3f}'
            CmStr = f'-{Planet.Bulk.CuncertaintyLower:.3f}'
        elif round(Planet.Bulk.Cmeasured, 4) == Planet.Bulk.Cmeasured and round(Planet.Bulk.CuncertaintyUpper, 4) == Planet.Bulk.CuncertaintyUpper:
            CmidStr = f'{Planet.Bulk.Cmeasured:.4f}'
            CpStr = f'+{Planet.Bulk.CuncertaintyUpper:.4f}'
            CmStr = f'-{Planet.Bulk.CuncertaintyLower:.4f}'
        else:
            CmidStr = f'{Planet.Bulk.Cmeasured:.5f}'
            CpStr = f'+{Planet.Bulk.CuncertaintyUpper:.5f}'
            CmStr = f'-{Planet.Bulk.CuncertaintyLower:.5f}'
        if FigMisc.TEX_INSTALLED:
            Planet.CMR2str = f'{CmidStr}\substack{{{CpStr} \\\\ {CmStr}}}'
            Planet.CMR2str5 = f'{Planet.Bulk.Cmeasured:.5f}\substack{{+{Planet.Bulk.CuncertaintyUpper:.5f} \\\\ -{Planet.Bulk.CuncertaintyLower:.5f}}}'
        else:
            Planet.CMR2str = f'{CmidStr}^{{{CpStr}}}_{{{CmStr}}}'
            Planet.CMR2str5 = f'{Planet.Bulk.Cmeasured:.5f}^{{+{Planet.Bulk.CuncertaintyUpper:.5f}}}_{{-{Planet.Bulk.CuncertaintyLower:.5f}}}'

    return Planet


def SetupNonSelfConsistent(Planet, Params):
    """ Initialize non-self-consistent layer arrays in Planet.
    """        # Verify that all depths have been specified
    if Planet.Do.NON_SELF_CONSISTENT:
        Planet.zb_km = 0
        if Planet.dzIceI_km == np.nan or Planet.dzIceI_km < 0:
            raise ValueError('Planet.dzIceI_km must be set to non-negative thickness for non-self-consistent ice modeling.')
        elif Planet.dzIceI_km == 0:
            Planet.Steps.nIceI = 0
        else:
            Planet.zb_km += Planet.dzIceI_km
            # Check that ice shell density is set, otherwise use default
            if Planet.Ocean.rhoCondMean_kgm3['Ih'] is None:
                Planet.Ocean.rhoCondMean_kgm3['Ih'] = Constants.STP_kgm3['Ih']
                log.warning('Planet.Sil.rhoCondMean_kgm3 is not set, using Constants.STP_kgm3["Ih"].')
            else:
                Planet.Do.CONSTANTPROPSEOS = True
        
            # Set up melting point viscosity
            if Planet.etaMelt_Pas is None:
                Planet.etaMelt_Pas = Constants.etaMelt_Pas[1]
                log.warning('Planet.Sil.etaMelt_Pas is not set, using Constants.etaMelt_Pas[1].')
            else:
                Planet.Do.CONSTANTPROPSEOS = True
        
            # Set up thermal conductivity
            if Planet.Ocean.kThermIce_WmK['Ih'] is None:
                Planet.Ocean.kThermIce_WmK['Ih'] = Constants.kThermIce_WmK['Ih']
                log.warning('Planet.Sil.kThermIce_WmK is not set, using Constants.kThermIce_WmK["Ih"].')
            else:
                Planet.Do.CONSTANTPROPSEOS = True
            
            # Set up creep parameters
            if Planet.Ocean.Eact_kJmol['Ih'] is None:
                Planet.Ocean.Eact_kJmol['Ih'] = Constants.Eact_kJmol[1]
                log.warning('Planet.Sil.Eact_kJmol is not set, using Constants.Eact_kJmol[1].')
            else:
                Planet.Do.CONSTANTPROPSEOS = True
                
            # Set up shear modulus in conducting layer
            if Planet.Ocean.GScondMean_GPa['Ih'] is None:
                Planet.Ocean.GScondMean_GPa['Ih'] = Constants.GS_GPa[1]
                log.warning('Planet.Sil.GScondMean_GPa is not set, using Constants.GS_GPa[1].')
            else:
                Planet.Do.CONSTANTPROPSEOS = True
            
            if Planet.Do.CONSTANTPROPSEOS:
                Planet.Ocean.constantProperties['Ih'] = {'rho_kgm3': Planet.Ocean.rhoCondMean_kgm3['Ih'],
                                    'Cp_JkgK': Constants.Cp_JkgK['Ih'],
                                    'alpha_pK': Constants.alphaIce_pK['Ih'],
                                    'kTherm_WmK': Planet.Ocean.kThermIce_WmK['Ih'],
                                    'VP_GPa': Constants.VP_GPa[1],
                                    'VS_GPa': Constants.VS_GPa[1],
                                    'KS_GPa': Constants.KS_GPa[1],
                                    'GS_GPa': Planet.Ocean.GScondMean_GPa['Ih'],
                                    'sigma_Sm': Planet.Ocean.sigmaIce_Sm['Ih'],
                                    'eta_Pas': Constants.etaIce_Pas[0]}
        
        if Planet.Do.BOTTOM_ICEIII:
            if (Planet.dzIceIII_km == np.nan or Planet.dzIceIII_km < 0):
                raise ValueError('Planet.Do.BOTTOM_ICEIII is True, but Planet.dzIceIII_km is not set or is negative.')
            else:
                Planet.zb_km += Planet.dzIceIII_km
                            # Check that ice shell density is set, otherwise use default
                if Planet.Ocean.rhoCondMean_kgm3['III'] is None:
                    Planet.Ocean.rhoCondMean_kgm3['III'] = Constants.STP_kgm3['III']
                    log.warning('Planet.Sil.rhoCondMean_kgm3 is not set, using Constants.STP_kgm3["III"].')
                else:
                    Planet.Do.CONSTANTPROPSEOS = True
            
                # Set up melting point viscosity
                if Planet.etaMeltIII_Pas is None:
                    Planet.etaMeltIII_Pas = Constants.etaMelt_Pas[3]
                    log.warning('Planet.Sil.etaMeltIII_Pas is not set, using Constants.etaMelt_Pas[3].')
                else:
                    Planet.Do.CONSTANTPROPSEOS = True
            
                # Set up thermal conductivity
                if Planet.Ocean.kThermIce_WmK['III'] is None:
                    Planet.Ocean.kThermIce_WmK['III'] = Constants.kThermIce_WmK['III']
                    log.warning('Planet.Sil.kThermIce_WmK is not set, using Constants.kThermIce_WmK["III"].')
                else:
                    Planet.Do.CONSTANTPROPSEOS = True
                
                # Set up creep parameters
                if Planet.Ocean.Eact_kJmol['III'] is None:
                    Planet.Ocean.Eact_kJmol['III'] = Constants.Eact_kJmol[3]
                    log.warning('Planet.Sil.Eact_kJmol is not set, using Constants.Eact_kJmol[3].')
                else:
                    Planet.Do.CONSTANTPROPSEOS = True
                    
                # Set up shear modulus in conducting layer
                if Planet.Ocean.GScondMean_GPa['III'] is None:
                    Planet.Ocean.GScondMean_GPa['III'] = Constants.GS_GPa[3]
                    log.warning('Planet.Sil.GScondMean_GPa is not set, using Constants.GS_GPa[3].')
                else:
                    Planet.Do.CONSTANTPROPSEOS = True
                
                if Planet.Do.CONSTANTPROPSEOS:
                    Planet.Ocean.constantProperties['III'] = {'rho_kgm3': Planet.Ocean.rhoCondMean_kgm3['III'],
                                        'Cp_JkgK': Constants.Cp_JkgK['III'],
                                        'alpha_pK': Constants.alphaIce_pK['III'],
                                        'kTherm_WmK': Planet.Ocean.kThermIce_WmK['III'],
                                        'VP_GPa': Constants.VP_GPa[3],
                                        'VS_GPa': Constants.VS_GPa[3],
                                        'KS_GPa': Constants.KS_GPa[3],
                                        'GS_GPa': Planet.Ocean.GScondMean_GPa['III'],
                                        'sigma_Sm': Planet.Ocean.sigmaIce_Sm['III'],
                                        }
        else:
            Planet.Steps.nIceIIILitho = 0
            Planet.Steps.nIceVLitho = 0
        if Planet.Do.BOTTOM_ICEV:
            if (Planet.dzIceV_km == np.nan or Planet.dzIceV_km < 0):
                raise ValueError('Planet.Do.BOTTOM_ICEV is True, but Planet.dzIceV_km is not set or is negative.')
            else:
                Planet.zb_km += Planet.dzIceV_km
                                        # Check that ice shell density is set, otherwise use default
                if Planet.Ocean.rhoCondMean_kgm3['V'] is None:
                    Planet.Ocean.rhoCondMean_kgm3['V'] = Constants.STP_kgm3['V']
                    log.warning('Planet.Sil.rhoCondMean_kgm3 is not set, using Constants.STP_kgm3["V"].')
                else:
                    Planet.Do.CONSTANTPROPSEOS = True
            
                # Set up melting point viscosity
                if Planet.etaMeltV_Pas is None:
                    Planet.etaMeltV_Pas = Constants.etaMelt_Pas[5]
                    log.warning('Planet.Sil.etaMeltV_Pas is not set, using Constants.etaMelt_Pas[5].')
                else:
                    Planet.Do.CONSTANTPROPSEOS = True
            
                # Set up thermal conductivity
                if Planet.Ocean.kThermIce_WmK['V'] is None:
                    Planet.Ocean.kThermIce_WmK['V'] = Constants.kThermIce_WmK['V']
                    log.warning('Planet.Sil.kThermIce_WmK is not set, using Constants.kThermIce_WmK["V"].')
                else:
                    Planet.Do.CONSTANTPROPSEOS = True
                
                # Set up creep parameters
                if Planet.Ocean.Eact_kJmol['V'] is None:
                    Planet.Ocean.Eact_kJmol['V'] = Constants.Eact_kJmol[5]
                    log.warning('Planet.Sil.Eact_kJmol is not set, using Constants.Eact_kJmol[5].')
                else:
                    Planet.Do.CONSTANTPROPSEOS = True
                    
                # Set up shear modulus in conducting layer
                if Planet.Ocean.GScondMean_GPa['V'] is None:
                    Planet.Ocean.GScondMean_GPa['V'] = Constants.GS_GPa[5]
                    log.warning('Planet.Sil.GScondMean_GPa is not set, using Constants.GS_GPa[5].')
                else:
                    Planet.Do.CONSTANTPROPSEOS = True
                
                if Planet.Do.CONSTANTPROPSEOS:
                    Planet.Ocean.constantProperties['V'] = {'rho_kgm3': Planet.Ocean.rhoCondMean_kgm3['V'],
                                        'Cp_JkgK': Constants.Cp_JkgK['V'],
                                        'alpha_pK': Constants.alphaIce_pK['V'],
                                        'kTherm_WmK': Planet.Ocean.kThermIce_WmK['V'],
                                        'VP_GPa': Constants.VP_GPa[5],
                                        'VS_GPa': Constants.VS_GPa[5],
                                        'KS_GPa': Constants.KS_GPa[5],
                                        'GS_GPa': Planet.Ocean.GScondMean_GPa['V'],
                                        'sigma_Sm': Planet.Ocean.sigmaIce_Sm['V'],
                                        }
        else:
            Planet.Steps.nIceVLitho = 0
        if Planet.Do.CLATHRATE:
            Planet.Steps.nClath = 0
        else:
            Planet.Steps.nClath = 0
        Planet.Steps.nSurfIce = Planet.Steps.nIceI+Planet.Steps.nIceIIILitho+Planet.Steps.nIceVLitho+Planet.Steps.nClath
        if Planet.D_km is None or Planet.D_km < 0:
            raise ValueError('Planet.D_km must be set to non-negative thickness for non-self-consistent ocean modeling.')
        elif Planet.D_km == 0:
            Planet.Do.NO_OCEAN = True
            Planet.Steps.nOcean = 0
        else:
            # Set up ocean density, if specified
            if Planet.Ocean.rhoMean_kgm3 is None:
                Planet.Ocean.rhoMean_kgm3 = Constants.STP_kgm3['0']
                log.warning('Planet.Ocean.rhoMean_kgm3 is not set, using Constants.STP_kgm3["0"].')
            else:
                Planet.Do.CONSTANTPROPSEOS = True

            # Set up thermal conductivity, if specified
            if Planet.Ocean.kThermWater_WmK is None:
                Planet.Ocean.kThermOcean_WmK = Constants.kThermWater_WmK
                log.warning('Planet.Ocean.kThermWater_WmK is not set, using Constants.kThermWater_WmK.')
            else:
                Planet.Do.CONSTANTPROPSEOS = True
                
            # Set up electrical condonctvity, if specified
            if Planet.Ocean.sigmaFixed_Sm is None:
                Planet.Ocean.sigmaFixed_Sm = Constants.sigmaH2O_Sm
                log.warning('Planet.Ocean.sigmaFixed_Sm is not set, using Constants.sigmaH2O_Sm.')
            else:
                Planet.Do.CONSTANTPROPSEOS = True
            
            if Planet.Do.CONSTANTPROPSEOS:
                Planet.Ocean.oceanConstantProperties = {'rho_kgm3': Planet.Ocean.rhoMean_kgm3,
                                    'Cp_JkgK': Constants.CpWater_JkgK,
                                    'alpha_pK': Constants.alphaWater_pK,
                                    'kTherm_WmK': Planet.Ocean.kThermWater_WmK,
                                    'VP_kms': Constants.VPOcean_kms,
                                    'VS_kms': Constants.VSOcean_kms,
                                    'sigma_Sm': Planet.Ocean.sigmaFixed_Sm,
                                    'eta_Pas': Constants.etaH2O_Pas,
                                    }
        Planet.Steps.nHydro = Planet.Steps.nSurfIce + Planet.Steps.nOcean
        if Planet.Core.Rmean_m is None or Planet.Core.Rmean_m < 0:
            raise ValueError('Planet.Core.Rmean_m must be set to non-negative radius for non-self-consistent inner modeling.')
        elif Planet.Core.Rmean_m == 0:
            Planet.Do.Fe_CORE = False
            Planet.Steps.nCore = 0
        Planet.Sil.mantleEOS = 'none'
        
        # Check if we have specified a valid ocean composition
        if Planet.Ocean.comp is None:
            Planet.Ocean.comp = 'none'
            Planet.Ocean.wOcean_ppt = 0.0
            Planet.Sil.PHydroMax_MPa = Planet.Bulk.Psurf_MPa
        else:
            # If so, then we will use the ocean composition to query bulk temperature andEC
            pass
        Planet.Steps.nTotal = Planet.Steps.nHydro + Planet.Steps.nSil + Planet.Steps.nCore
    return Planet