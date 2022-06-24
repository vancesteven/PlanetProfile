""" Docstring explaining what we're doing here """

import os
import numpy as np
import logging
from collections.abc import Iterable
from PlanetProfile import _ROOT
from PlanetProfile.GetConfig import FigMisc, FigLbl
from PlanetProfile.Thermodynamics.HydroEOS import GetOceanEOS
from PlanetProfile.Thermodynamics.InnerEOS import GetInnerEOS
from PlanetProfile.Thermodynamics.Clathrates.ClathrateProps import ClathDissoc
from PlanetProfile.Utilities.PPversion import ppVerNum, CheckCompat
from PlanetProfile.Utilities.defineStructs import DataFilesSubstruct, FigureFilesSubstruct, Constants

# Assign logger
log = logging.getLogger('PlanetProfile')

def SetupInit(Planet, Params):

    # Print version number
    log.debug(f'-- PlanetProfile v{ppVerNum} --')
    if ppVerNum[-3:] == 'dev': log.debug('This version is in development.')

    # Check dependency compatibility
    CheckCompat('seafreeze')  # SeaFreeze
    if Planet.Ocean.comp == 'Seawater': CheckCompat('gsw')  # Gibbs Seawater
    if Planet.Do.TAUP_SEISMIC: CheckCompat('obspy')  # TauP (accessed as obspy.taup)
    if Params.CALC_NEW_INDUCT: CheckCompat('MoonMag')  # MoonMag

    # Waterless bodies. We have to do this first, before filename
    # generation, to ensure ocean comp is set.
    if Planet.Do.NO_H2O:
        log.info('Modeling a waterless body.')
        if Planet.Bulk.qSurf_Wm2 is None:
            raise ValueError('Bulk.qSurf_Wm2 must be set in order to model waterless bodies.')
        Planet.Ocean.QfromMantle_W = Planet.Bulk.qSurf_Wm2 * 4*np.pi * Planet.Bulk.R_m**2
        Planet.qSurf_Wm2 = Planet.Bulk.qSurf_Wm2
        Planet.qCon_Wm2 = np.nan
        Planet.etaConv_Pas = np.nan
        Planet.Pb_MPa = Planet.Bulk.Psurf_MPa
        Planet.PbI_MPa = Planet.Bulk.Psurf_MPa
        Planet.Sil.PHydroMax_MPa = Planet.Bulk.Psurf_MPa
        Planet.Bulk.Tb_K = Planet.Bulk.Tsurf_K
        Planet.zb_km = 0.0
        Planet.Steps.nIceI = 0
        Planet.Steps.nIbottom = 0
        Planet.Steps.nSurfIce = 0
        Planet.Steps.nOceanMax = 1
        Planet.Steps.nHydroMax = 1
        Planet.Ocean.comp = 'none'
        Planet.Ocean.wOcean_ppt = 0.0
        Planet.Ocean.deltaP = 0.1
        Planet.Do.NO_ICE_CONVECTION = True
        if Planet.Do.POROUS_ROCK:
            # Generate zero-yielding ocean "EOS" for use in porosity calculations
            # Note that there must be enough input points for creating spline
            # interpolators, even though we will not use them.
            Planet.Ocean.EOS = GetOceanEOS('none', None, np.linspace(0, 1, 4), np.linspace(0, 1, 4), None)

    # Get filenames for saving/loading
    Params.DataFiles, Params.FigureFiles = SetupFilenames(Planet, Params)

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

    if not Planet.Do.NO_H2O:
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
            if(Planet.Ocean.PHydroMax_MPa < 209.9):
                raise ValueError('Hydrosphere max pressure is less than the pressure of the ice I-III phase transition, ' +
                                 'but Do.BOTTOM_ICEIII is True.')
            if(Planet.Bulk.Tb_K > Planet.Bulk.TbIII_K):
                log.warning('Bottom temperature of ice I (Tb_K) is greater than bottom temperature of underplate ' +
                            'ice III (TbIII_K). This likely represents a non-equilibrium state.')
        else:
            if(Planet.Ocean.PHydroMax_MPa < 344.3):
                raise ValueError('Hydrosphere max pressure is less than the pressure of the ice III-V phase transition, ' +
                                 'but Do.BOTTOM_ICEV is True.')
            if(Planet.Bulk.Tb_K > Planet.Bulk.TbIII_K):
                log.warning('Bottom temperature of ice I (Tb_K) is greater than bottom temperature of underplate ' +
                            'ice III (TbIII_K). This likely represents a non-equilibrium state.')
            if(Planet.Bulk.TbIII_K > Planet.Bulk.TbV_K):
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
        # Check ocean parameter space, and prevent setup from taking forever if we have a deep ocean:
        if Planet.Ocean.THydroMax_K < Planet.Bulk.Tb_K:
            raise ValueError(f'Ocean.THydroMax_K of {Planet.Ocean.THydroMax_K} is less than Bulk.Tb_K of {Planet.Bulk.Tb_K}.')
        elif Planet.Bulk.Tb_K + 30 < Planet.Ocean.THydroMax_K:
            TOcean_K = np.concatenate((np.linspace(Planet.Bulk.Tb_K, Planet.Bulk.Tb_K + 30, int(100/Planet.Ocean.deltaT), endpoint=False),
                                       np.arange(Planet.Bulk.Tb_K + 30, Planet.Ocean.THydroMax_K, 2)))
        else:
            TOcean_K = np.arange(Planet.Bulk.Tb_K, Planet.Ocean.THydroMax_K, Planet.Ocean.deltaT)
        Planet.Ocean.EOS = GetOceanEOS(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, POcean_MPa, TOcean_K,
                                       Planet.Ocean.MgSO4elecType, rhoType=Planet.Ocean.MgSO4rhoType,
                                       scalingType=Planet.Ocean.MgSO4scalingType, FORCE_NEW=Params.FORCE_EOS_RECALC,
                                       phaseType=Planet.Ocean.phaseType, EXTRAP=Params.EXTRAP_OCEAN)
        # Get separate, simpler EOS for evaluating the melting curve
        Pmelt_MPa = np.arange(Planet.PfreezeLower_MPa, Planet.PfreezeUpper_MPa, Planet.PfreezeRes_MPa)
        Planet.Ocean.meltEOS = GetOceanEOS(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, Pmelt_MPa,
                                           np.linspace(Planet.Bulk.Tb_K - 0.01, Planet.Bulk.Tb_K + 0.01, 11), None,
                                           phaseType=Planet.Ocean.phaseType, FORCE_NEW=True, MELT=True)

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
        if not Planet.Do.PORE_EOS_DIFFERENT:
            Planet.Sil.PHydroMax_MPa = Planet.Ocean.PHydroMax_MPa

        if Planet.Sil.porosType == 'Han2014':
            if Planet.Sil.phiRockMax_frac is None:
                log.warning('Sil.phiRockMax_frac is not set. Using arbitrary max porosity of 0.7.')
                Planet.Sil.phiRockMax_frac = 0.7
        else:
            Planet.Do.FIXED_POROSITY = True
        if Planet.Sil.phiRangeMult <= 1:
            raise ValueError(f'Sil.phiRangeMult = {Planet.Sil.phiRangeMult}, but it must be greater than 1.')
    else:
        Planet.Sil.porosType = 'none'
        Planet.Sil.poreH2Orho_kgm3 = 0
        Planet.Sil.phiRockMax_frac = 0
        Planet.Sil.poreComp = Planet.Ocean.comp
        Planet.Sil.wPore_ppt = Planet.Ocean.wOcean_ppt

    # Porous ice
    if not Planet.Do.POROUS_ICE:
        Planet.Ocean.phiMax_frac = {key:0 for key in Planet.Ocean.phiMax_frac.keys()}

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
                                     EXTRAP=Params.EXTRAP_SIL)

        # Pore fluids if present
        if Planet.Do.POROUS_ROCK:
            if Planet.Do.NO_H2O:
                Ppore_MPa, Tpore_K = (np.linspace(0, 1, 4) for _ in range(2))
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
                                             phaseType=Planet.Ocean.phaseType, EXTRAP=Params.EXTRAP_OCEAN)

            # Make sure Sil.phiRockMax_frac is set in case we're using a porosType that doesn't require it
            if Planet.Sil.phiRockMax_frac is None or Planet.Sil.porosType != 'Han2014':
                Planet.Sil.phiRockMax_frac = Planet.Sil.EOS.fn_phi_frac(0, 0)

        # Iron core if present
        if Planet.Do.Fe_CORE:
            Planet.Core.EOS = GetInnerEOS(Planet.Core.coreEOS, EOSinterpMethod=Params.lookupInterpMethod, Fe_EOS=True,
                                          kThermConst_WmK=Planet.Core.kTherm_WmK, EXTRAP=Params.EXTRAP_Fe,
                                          wFeCore_ppt=Planet.Core.wFe_ppt, wScore_ppt=Planet.Core.wS_ppt)

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


def SetupFilenames(Planet, Params, exploreAppend=None, figExploreAppend=None):
    """ Generate filenames for saving data and figures.
    """
    datPath = Planet.bodyname
    if Planet.bodyname == 'Test':
        datPath = os.path.join(_ROOT, datPath)
    figPath = os.path.join(datPath, 'figures')

    # Account for differing ocean/pore composition here, since we need it for filenames
    # Use ocean composition and salinity if user has not specified different ones for pore space
    Planet.Do.PORE_EOS_DIFFERENT = False
    if Planet.Sil.poreComp is None:
        Planet.Sil.poreComp = Planet.Ocean.comp
    elif Planet.Sil.poreComp != Planet.Ocean.comp:
        Planet.Do.PORE_EOS_DIFFERENT = True
    if Planet.Sil.wPore_ppt is None:
        Planet.Sil.wPore_ppt = Planet.Ocean.wOcean_ppt
    elif Planet.Sil.wPore_ppt != Planet.Ocean.wOcean_ppt:
        Planet.Do.PORE_EOS_DIFFERENT = True
        
    saveBase = Planet.name + 'Profile_'
    saveLabel = ''
    label = ''
    if Planet.Do.NO_H2O:
        saveLabel += f'NoH2O_qSurf{Planet.Bulk.qSurf_Wm2*1e3:.1f}mWm2'
        setStr = f'$q_\mathrm{{surf}}\,{Planet.Bulk.qSurf_Wm2*FigLbl.qMult:.1f}\,\si{{{FigLbl.fluxUnits}}}$'
        label += setStr
    else:
        if Planet.Ocean.comp == 'PureH2O':
            saveLabel += f'{Planet.Ocean.comp}_Tb{Planet.Bulk.Tb_K}K'
            setStr = f'Pure \ce{{H2O}}'
            label += f'{setStr}, $T_b\,\SI{{{Planet.Bulk.Tb_K}}}{{K}}$'
        else:
            saveLabel += f'{Planet.Ocean.comp}_{Planet.Ocean.wOcean_ppt:.1f}ppt' + \
                        f'_Tb{Planet.Bulk.Tb_K}K'
            setStr = f'${Planet.Ocean.wOcean_ppt*FigLbl.wMult:.1f}\,\si{{{FigLbl.wUnits}}}$ \ce{{{Planet.Ocean.comp}}}'
            label += setStr + \
                f', $T_b\,\SI{{{Planet.Bulk.Tb_K}}}{{K}}$'
        if Planet.Do.CLATHRATE:
            saveLabel += '_Clathrates'
            label += ' w/clath'
        if Planet.Do.POROUS_ICE:
            saveLabel += '_PorousIce'
            label += ' w/$\phi_{ice}$'
        if Planet.Do.PORE_EOS_DIFFERENT:
            if Planet.Sil.poreComp == 'PureH2O':
                saveLabel += f'_{Planet.Sil.poreComp}Pores'
                label += f'Pure \ce{{H2O}} pores'
            else:
                saveLabel += f'_{Planet.Sil.poreComp}_{Planet.Sil.wPore_ppt:.1f}pptPores'
                label += f', {Planet.Sil.wPore_ppt*FigLbl.wMult:.1f}\,\si{{{FigLbl.wUnits}}} \ce{{{Planet.Sil.poreComp}}} pores'
    if Planet.Sil.mantleEOSName is not None: saveLabel += f'_{Planet.Sil.mantleEOSname}'

    Planet.saveLabel = saveLabel
    Planet.tradeLabel = f'{label}, $C/MR^2\,{Planet.Bulk.Cmeasured}\pm{Planet.Bulk.Cuncertainty}$'
    Planet.label = label
    if Params.DO_INDUCTOGRAM:
        inductBase = f'{Planet.name}_{Params.Induct.inductOtype}'
    else:
        inductBase = None
    DataFiles = DataFilesSubstruct(datPath, saveBase + saveLabel, Planet.Ocean.comp, inductBase=inductBase,
                                   exploreAppend=exploreAppend)
    FigureFiles = FigureFilesSubstruct(figPath, saveBase + saveLabel, FigMisc.xtn,
                                       comp=Planet.Ocean.comp, inductBase=inductBase,
                                       exploreAppend=figExploreAppend)

    return DataFiles, FigureFiles


def SetupLayers(Planet):
    """ Initialize layer arrays in Planet.
    """

    if not Planet.Do.NO_H2O:
        nOceanMax = int(Planet.Ocean.PHydroMax_MPa / Planet.Ocean.deltaP)
        Planet.Steps.nHydroMax = Planet.Steps.nClath + Planet.Steps.nIceI + Planet.Steps.nIceIIILitho + Planet.Steps.nIceVLitho + nOceanMax
    
    Planet.phase = np.zeros(Planet.Steps.nHydroMax, dtype=np.int_)
    Planet.P_MPa, Planet.T_K, Planet.r_m, Planet.rho_kgm3, \
        Planet.Cp_JkgK, Planet.alpha_pK, Planet.g_ms2, Planet.phi_frac, \
        Planet.sigma_Sm, Planet.z_m, Planet.MLayer_kg, Planet.VLayer_m3, Planet.kTherm_WmK, \
        Planet.Htidal_Wm3, Planet.Ppore_MPa, Planet.rhoMatrix_kgm3, Planet.rhoPore_kgm3 = \
        (np.zeros(Planet.Steps.nHydroMax) for _ in range(17))

    # Layer property initialization for surface
    Planet.z_m[0] = 0.0  # Set first layer depth to zero (layer properties correspond to outer radius)
    Planet.r_m[0] = Planet.Bulk.R_m  # Set first layer to planetary surface radius
    Planet.g_ms2[0] = Constants.G * Planet.Bulk.M_kg / Planet.Bulk.R_m**2  # Set first layer gravity at surface
    Planet.T_K[0] = Planet.Bulk.Tsurf_K  # Set first layer surface temp
    Planet.P_MPa[0] = Planet.Bulk.Psurf_MPa  # Set first layer to surface pressure

    return Planet
