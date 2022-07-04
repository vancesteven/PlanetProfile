# Import necessary Python modules
import os, sys, time, importlib
import numpy as np
import logging
from scipy.io import savemat, loadmat
from copy import deepcopy
from distutils.util import strtobool
from collections.abc import Iterable
from os.path import isfile
from glob import glob as FilesMatchingPattern

# Import all function definitions for this file
from PlanetProfile import _Defaults, _TestImport, CopyCarefully
from PlanetProfile.GetConfig import Params as configParams, FigMisc
from PlanetProfile.MagneticInduction.MagneticInduction import MagneticInduction, ReloadInduction, Benm2absBexyz
from PlanetProfile.MagneticInduction.Moments import InductionResults, Excitations as Mag
from PlanetProfile.Plotting.ProfilePlots import GeneratePlots, PlotExploreOgram
from PlanetProfile.Plotting.MagPlots import GenerateMagPlots, PlotInductOgram, \
    PlotInductOgramPhaseSpace
from PlanetProfile.Thermodynamics.LayerPropagators import IceLayers, OceanLayers, InnerLayers
from PlanetProfile.Thermodynamics.Electrical import ElecConduct
from PlanetProfile.Thermodynamics.Seismic import SeismicCalcs
from PlanetProfile.Utilities.defineStructs import Constants, FigureFilesSubstruct, PlanetStruct, ExplorationResults
from PlanetProfile.Utilities.SetupInit import SetupInit, SetupFilenames
from PlanetProfile.Utilities.SummaryTables import GetLayerMeans, PrintGeneralSummary, PrintLayerSummaryLatex, PrintLayerTableLatex

# Parallel processing
import multiprocessing as mtp
mtpFork = mtp.get_context('fork')
# Assign logger
log = logging.getLogger('PlanetProfile')

""" MAIN RUN BLOCK """
def run(bodyname=None, opt=None, fNames=None):

    # Copy global Params settings to local variable so we can add things like filenames
    Params = configParams
    
    if fNames is None and bodyname is None:
        log.info('No body name entered. Defaulting to Europa.')
        bodyname = 'Europa'
    bodyname = bodyname.capitalize()
    if bodyname != '':
        log.info(f'Body name: {bodyname}')

    if opt is not None or fNames is not None:
        Params, fNames = ExecOpts(Params, bodyname, opt, fNames=fNames)

    """ Run PlanetProfile """
    if Params.DO_INDUCTOGRAM:
        if bodyname == '':
            raise ValueError('A single body must be specified for an InductOgram.')
        else:
            Induction, Params = InductOgram(bodyname, Params)
        if not Params.SKIP_PLOTS:
            PlotInductOgram(Induction, Params)

            if Params.COMPARE:
                inductOgramFiles = FilesMatchingPattern(os.path.join(Params.DataFiles.fNameInduct+'*.mat'))
                Params.nModels = np.size(inductOgramFiles)
                InductionList = np.empty(Params.nModels, dtype=object)
                InductionList[0] = deepcopy(Induction)
                # Move the filename for this run to the front of the list
                if Params.DataFiles.inductOgramFile in inductOgramFiles:
                    inductOgramFiles.remove(Params.DataFiles.inductOgramFile)
                    inductOgramFiles.insert(0, Params.DataFiles.inductOgramFile)
                for i, reloadInduct in enumerate(inductOgramFiles[1:]):
                    InductionList[i+1] = ReloadInductOgram(bodyname, Params, fNameOverride=reloadInduct)[0]
            else:
                InductionList = [Induction]
            PlotInductOgramPhaseSpace(InductionList, Params)
    elif Params.DO_EXPLOREOGRAM:
        if bodyname == '':
            raise ValueError('A single body must be specified for an ExploreOgram.')
        else:
            Exploration, Params = ExploreOgram(bodyname, Params)
        if not Params.SKIP_PLOTS:
            if Params.COMPARE:
                exploreOgramFiles = FilesMatchingPattern(os.path.join(Params.DataFiles.fNameExplore+'*.mat'))
                Params.nModels = np.size(exploreOgramFiles)
                ExplorationList = np.empty(Params.nModels, dtype=object)
                ExplorationList[0] = deepcopy(Exploration)
                # Move the filename for this run to the front of the list
                if Params.DataFiles.exploreOgramFile in exploreOgramFiles:
                    exploreOgramFiles.remove(Params.DataFiles.exploreOgramFile)
                    exploreOgramFiles.insert(0, Params.DataFiles.exploreOgramFile)
                for i, reloadExplore in enumerate(exploreOgramFiles[1:]):
                    ExplorationList[i+1] = ReloadExploreOgram(bodyname, Params, fNameOverride=reloadExplore)[0]
            else:
                ExplorationList = [Exploration]
            PlotExploreOgram(ExplorationList, Params)
    else:
        # Set timekeeping for recording elapsed times
        tMarks = np.empty(0)
        tMarks = np.append(tMarks, time.time())
        if opt == 'reload':
            PlanetList = np.empty(0, dtype=object)
            for fName in fNames:
                loadName = os.path.join(bodyname, fName)
                Planet, Params = ReloadProfile(None, None, fnameOverride=loadName)
                PlanetList = np.append(PlanetList, deepcopy(Planet))
        else:
            Params, loadNames = LoadPPfiles(Params, fNames, bodyname=bodyname)
            PlanetList = np.empty(Params.nModels, dtype=object)
            # Run main model first, so that it always appears as 0-index
            PlanetList[0] = importlib.import_module(loadNames[0]).Planet
            PlanetList[0].index = 1
            PlanetList[0], Params = PlanetProfile(PlanetList[0], Params)
            tMarks = np.append(tMarks, time.time())
            if Params.RUN_ALL_PROFILES:
                for i,loadName in enumerate(loadNames[1:]):
                    PlanetList[i+1] = deepcopy(importlib.import_module(loadName).Planet)
                    PlanetList[i+1].index = i+2
                    PlanetList[i+1], Params = PlanetProfile(PlanetList[i+1], Params)
                    tMarks = np.append(tMarks, time.time())

        dt = np.diff(tMarks)
        if np.size(dt) > 0:
            log.debug('Elapsed time:\n' + '\n'.join([f'    {dt[i]:.3f} s for {Planet.saveLabel}' for i,Planet in enumerate(PlanetList)]))

        """ Post-processing """
        # Loading BodyProfile...txt files to plot them together
        if Params.COMPARE and not Params.RUN_ALL_PROFILES:
            fNamesToCompare = np.array(FilesMatchingPattern(os.path.join(PlanetList[0].bodyname, f'{PlanetList[0].name}Profile*.txt')))
            isProfile = [Params.DataFiles.saveFile != fName and 'mantle' not in fName for fName in fNamesToCompare]
            fProfiles = fNamesToCompare[isProfile]
            nCompare = np.size(fProfiles) + 1
            log.info('Loading comparison profiles.')
            CompareList = np.empty(nCompare, dtype=object)
            CompareList[0] = PlanetList[0]
            for i,fName in enumerate(fProfiles):
                CompareList[i+1], _ = ReloadProfile(deepcopy(CompareList[0]), Params, fnameOverride=fName)
        else:
            CompareList = PlanetList
            
        MULTIPLOT = Params.COMPARE or Params.RUN_ALL_PROFILES
            
        # Get large-scale layer properties, which are needed for table outputs and some plots
        if Params.DISP_LAYERS or Params.DISP_TABLE or ((not Params.SKIP_PLOTS) and MULTIPLOT):
            CompareList, Params = GetLayerMeans(CompareList, Params)
            
        # Plot combined figures
        if (not Params.SKIP_PLOTS) and MULTIPLOT:
            if Params.ALL_ONE_BODY:
                comparePath = os.path.join(CompareList[0].bodyname, 'figures')
                compareBase = f'{CompareList[0].name}Comparison'
            else:
                comparePath = Params.compareDir
                compareBase = 'Comparison'

            Params.FigureFiles = FigureFilesSubstruct(comparePath, compareBase, FigMisc.xtn)
            GeneratePlots(CompareList, Params)

            if Params.CALC_CONDUCT and not Params.SKIP_INDUCTION:
                GenerateMagPlots(CompareList, Params)

        # Print table outputs
        if Params.DISP_LAYERS or Params.DISP_TABLE:
            if Params.DISP_LAYERS:
                PrintLayerSummaryLatex(CompareList, Params)
                PrintLayerTableLatex(CompareList, Params)
            if Params.DISP_TABLE:
                PrintGeneralSummary(CompareList, Params)

    return

""" END MAIN RUN BLOCK """


def PlanetProfile(Planet, Params):

    if Params.CALC_NEW:
        # Initialize
        Planet, Params = SetupInit(Planet, Params)
        if not Planet.Do.NO_H2O:
            Planet = IceLayers(Planet, Params)
            Planet = OceanLayers(Planet, Params)
        Planet = InnerLayers(Planet, Params)
        Planet = ElecConduct(Planet, Params)
        Planet = SeismicCalcs(Planet, Params)

        # Save data after modeling
        if not Params.NO_SAVEFILE:
            WriteProfile(Planet, Params)
    else:
        # Reload previous run
        Planet, Params = ReloadProfile(Planet, Params)

    # Main plotting functions
    if ((not Params.SKIP_PLOTS) and not (Params.DO_INDUCTOGRAM or Params.DO_EXPLOREOGRAM)) \
        and Planet.Do.VALID:
        # Calculate large-scale layer properties
        PlanetList, Params = GetLayerMeans(np.array([Planet]), Params)
        # Plotting functions

        GeneratePlots(PlanetList, Params)
        Planet = PlanetList[0]

    # Magnetic induction calculations and plots
    if Params.CALC_CONDUCT and not Params.SKIP_INDUCTION:
        # Calculate induced magnetic moments
        Planet, Params = MagneticInduction(Planet, Params)

        # Plot induced dipole surface strength
        if ((not Params.SKIP_PLOTS) and Planet.Do.VALID) and \
            not (Params.DO_INDUCTOGRAM or Params.DO_EXPLOREOGRAM):
            GenerateMagPlots([Planet], Params)

    PrintCompletion(Planet, Params)
    return Planet, Params


def HydroOnly(Planet, Params):
    """ Wrapper for PlanetProfile function, up through hydrosphere calculations,
        for parameter exploration that needs only to redo the interior.
    """

    Planet, Params = SetupInit(Planet, Params)
    if not Planet.Do.NO_H2O:
        Planet = IceLayers(Planet, Params)
        Planet = OceanLayers(Planet, Params)

    PrintCompletion(Planet, Params)
    return Planet, Params


def InteriorEtc(Planet, Params):
    """ Wrapper for PlanetProfile function, minus hydrosphere calculations,
        for parameter exploration that needs only to redo the interior.
    """

    Planet = InnerLayers(Planet, Params)
    Planet = ElecConduct(Planet, Params)
    Planet = SeismicCalcs(Planet, Params)
    if not Params.SKIP_INDUCTION and (Params.CALC_CONDUCT and Params.CALC_NEW_INDUCT):
        # Calculate induced magnetic moments
        Planet, Params = MagneticInduction(Planet, Params)

    PrintCompletion(Planet, Params)
    return Planet, Params


def InductionOnly(Planet, Params):
    """ Wrapper for MagneticInduction function similar to above that includes PrintCompletion.
    """

    Planet, Params = MagneticInduction(Planet, Params)

    PrintCompletion(Planet, Params)
    return Planet, Params


def PrintCompletion(Planet, Params):
    """ Print a message at the end of each PlanetProfile run estimating completion time.
    """

    if Planet.index is None:
        indicator = ''
        ending = '!'
    else:
        indicator = f' {Planet.index}/{Params.nModels}'
        if Planet.index == Params.nModels:
            ending = '!'
        else:
            tNow_s = time.time()
            tTot_s = Params.nModels / Planet.index * (tNow_s - Params.tStart_s)
            tRemain_s = (Params.tStart_s + tTot_s - tNow_s)
            remain = ''
            if tRemain_s > 3600:
                remain += f' {int(tRemain_s/3600)} hr'
                tRemain_s = tRemain_s % 3600
                if tRemain_s > 60:
                    remain += f' {int(tRemain_s/60)} min'
            else:
                if tRemain_s > 60:
                    remain += f' {int(tRemain_s/60)} min'
                remain += f' {int(tRemain_s % 60)} s'

            ending = f'. Approx.{remain} remaining.'
    log.profile(f'Profile{indicator} complete{ending}')
    return


def ExecOpts(Params, bodyname, opt, fNames=None):
    """ Actions to take if opt is passed to run(). Params may be changed,
        so we return Params.
    """

    if opt is not None:
        if opt == 'clear':
            fNamesToClear = FilesMatchingPattern(os.path.join(bodyname, '*.txt'))
            if len(fNamesToClear) > 0:
                log.info(f'Clearing previous run files for {bodyname}:')
                for fName in fNamesToClear:
                    log.info(f'    {fName}')
                    os.remove(fName)
                log.info(f'{bodyname} files cleared.')
            else:
                log.warning(f'Attempted to remove previous run files for {bodyname}, but found none.')
            if not Params.CALC_NEW:
                log.warning('CALC_NEW is set to False in configPP.py, but files are being cleared. ' +
                            'CALC_NEW has been forced on for this run.')
                Params.CALC_NEW = True
        elif opt == 'compare':
            log.info('Comparing with other profiles from this body.')
            Params.COMPARE = True
        elif opt == 'all':
            log.info(f'Running all available profiles for {bodyname}.')
            Params.RUN_ALL_PROFILES = True
            Params.COMPARE = True
        elif opt == 'inductogram':
            Params.DO_INDUCTOGRAM = True
            Params.NO_SAVEFILE = True
        elif opt == 'reload':
            Params.CALC_NEW = False
        else:
            log.warning(f'Unrecognized option: {opt}. Skipping.')

    if fNames is not None:
        # Set up fNames as a list of filenames to load and run
        Params.SPEC_FILE = True
        if np.size(fNames) > 1:
            Params.RUN_ALL_PROFILES = True
            Params.COMPARE = True
        if opt == 'reload':
            for fName in fNames:
                if os.path.split(fName)[0] == '':
                    expected = os.path.join(bodyname, fName)
                else:
                    expected = fName
                if not os.path.isfile(expected):
                    default = os.path.join(_Defaults, bodyname, fName)
                    if os.path.isfile(default):
                        CopyCarefully(default, expected)
                    else:
                        log.warning(f'{expected} does not exist and no default was found at {default}. Skipping.')
                        fNames.remove(fName)

            if np.size(fNames) == 0:
                raise ValueError('None of the specified PP files were found.')
        else:
            for fName in fNames:
                if os.path.split(fName)[0] == '':
                    expected = os.path.join(bodyname, fName)
                else:
                    expected = fName
                if not os.path.isfile(expected):
                    if os.path.split(fName)[0] == '':
                        default = os.path.join(_Defaults, bodyname, fName)
                    else:
                        default = os.path.join(_Defaults, fName)
                    if os.path.isfile(default):
                        CopyCarefully(default, expected)
                    else:
                        log.warning(f'{expected} does not exist and no default was found at {default}. Skipping.')
                        fNames.remove(fName)

            if np.size(fNames) == 0:
                raise ValueError('None of the specified PP files were found.')

    return Params, fNames


def WriteProfile(Planet, Params):
    """ Write out all profile calculations to disk """
    headerLines = [
        f'MoI label = {Planet.tradeLabel}',
        f'Iron core = {Planet.Do.Fe_CORE}',
        f'Silicate EOS file = {Planet.Sil.mantleEOS}',
        f'Iron core EOS file = {Planet.Core.coreEOS}',
        f'Ocean salt = {Planet.Ocean.comp}',
        f'Pore salt = {Planet.Sil.poreComp}',
        f'wOcean_ppt = {Planet.Ocean.wOcean_ppt:.3f}',
        f'wPore_ppt = {Planet.Sil.wPore_ppt:.3f}',
        f'R_m = {Planet.Bulk.R_m:.3f}',
        f'M_kg = {Planet.Bulk.M_kg:.5e}',
        f'Cmeasured = {Planet.Bulk.Cmeasured}',
        f'Cuncertainty = {Planet.Bulk.Cuncertainty}',
        f'Psurf_MPa = {Planet.Bulk.Psurf_MPa:.3f}',
        f'Tsurf_K = {Planet.Bulk.Tsurf_K:.3f}',
        f'qSurf_Wm2 = {Planet.qSurf_Wm2:.3e}',
        f'qCon_Wm2 = {Planet.qCon_Wm2:.3e}',
        f'Tb_K = {Planet.Bulk.Tb_K}',
        f'zb_km = {Planet.zb_km:.3f}',
        f'zClath_m = {Planet.zClath_m:.3f}',
        f'D_km = {Planet.D_km:.3f}',
        f'Pb_MPa = {Planet.Pb_MPa:.3f}',
        f'PbI_MPa = {Planet.PbI_MPa:.3f}',
        f'deltaP = {Planet.Ocean.deltaP:.3f}',
        f'Mtot_kg = {Planet.Mtot_kg:.6e}',
        f'CMR2mean = {Planet.CMR2mean:.5f}',
        f'CMR2less = {Planet.CMR2less:.5f}',
        f'CMR2more = {Planet.CMR2more:.5f}',
        f'QfromMantle_W = {Planet.Ocean.QfromMantle_W:.6e}',
        f'rhoOcean_kgm3 = {Planet.Ocean.rhoMean_kgm3:.3f}',
        f'phiRockMax = {Planet.Sil.phiRockMax_frac:.3f}',
        f'Qrad_Wkg = {Planet.Sil.Qrad_Wkg:.3f}',
        f'HtidalSil_Wm3 = {Planet.Sil.HtidalMean_Wm3:.3f}',
        f'RsilMean_m = {Planet.Sil.Rmean_m:.3f}',
        f'RsilRange_m = {Planet.Sil.Rrange_m:.3f}',
        f'rhoSil_kgm3 = {Planet.Sil.rhoMean_kgm3:.3f}',
        f'RcoreMean_m = {Planet.Core.Rmean_m:.3f}',
        f'RcoreRange_m = {Planet.Core.Rrange_m:.3f}',
        f'rhoCore_kgm3 = {Planet.Core.rhoMean_kgm3:.3f}',
        f'MH2O_kg = {Planet.MH2O_kg:.5e}',
        f'Mrock_kg = {Planet.Mrock_kg:.5e}',
        f'Mcore_kg = {Planet.Mcore_kg:.5e}',
        f'Mice_kg = {Planet.Mice_kg:.5e}',
        f'Msalt_kg = {Planet.Msalt_kg:.5e}',
        f'MporeSalt_kg = {Planet.MporeSalt_kg:.5e}',
        f'Mocean_kg = {Planet.Mocean_kg:.5e}',
        f'Mfluid_kg = {Planet.Mfluid_kg:.5e}',
        f'MporeFluid_kg = {Planet.MporeFluid_kg:.5e}',
        f'Mclath_kg = {Planet.Mclath_kg:.5e}',
        f'MclathGas_kg = {Planet.MclathGas_kg:.5e}',
        f'sigmaOceanMean_Sm = {Planet.Ocean.sigmaMean_Sm:.3f}',
        f'sigmaPoreMean_Sm = {Planet.Sil.sigmaPoreMean_Sm:.3f}',
        f'sigmaPorousLayerMean_Sm = {Planet.Sil.sigmaPorousLayerMean_Sm:.3f}',
        f'Tconv_K = {Planet.Tconv_K:.3e}',
        f'etaConv_Pas = {Planet.etaConv_Pas:.3e}',
        f'RaConvect = {Planet.RaConvect:.2e}',
        f'RaConvectIII = {Planet.RaConvectIII:.2e}',
        f'RaConvectV = {Planet.RaConvectV:.2e}',
        f'RaCrit = {Planet.RaCrit:.2e}',
        f'RaCritIII = {Planet.RaCritIII:.2e}',
        f'RaCritV = {Planet.RaCritV:.2e}',
        f'eLid_m = {Planet.eLid_m:.3f}',
        f'eLidIII_m = {Planet.eLidIII_m:.3f}',
        f'eLidV_m = {Planet.eLidV_m:.3f}',
        f'Dconv_m = {Planet.Dconv_m:.3f}',
        f'DconvIII_m = {Planet.DconvIII_m:.3f}',
        f'DconvV_m = {Planet.DconvV_m:.3f}',
        f'deltaTBL_m = {Planet.deltaTBL_m:.3f}',
        f'deltaTBLIII_m = {Planet.deltaTBLIII_m:.3f}',
        f'deltaTBLV_m = {Planet.deltaTBLV_m:.3f}',
        f'Porous ice = {Planet.Do.POROUS_ICE}',
        f'Steps.nClath = {Planet.Steps.nClath:d}',
        f'Steps.nIceI = {Planet.Steps.nIceI:d}',
        f'Steps.nIceIIILitho = {Planet.Steps.nIceIIILitho:d}',
        f'Steps.nIceVLitho = {Planet.Steps.nIceVLitho:d}',
        f'Steps.nHydro = {Planet.Steps.nHydro:d}',
        f'Steps.nSil = {Planet.Steps.nSil:d}',
        f'Steps.nCore = {Planet.Steps.nCore:d}']
    colHeaders = [' P (MPa)'.ljust(24),
                  'T (K)'.ljust(24),
                  'r (m)'.ljust(24),
                  'phase ID'.ljust(8),
                  'rho (kg/m3)'.ljust(24),
                  'Cp (J/kg/K)'.ljust(24),
                  'alpha (1/K)'.ljust(24),
                  'g (m/s2)'.ljust(24),
                  'phi (void/solid frac)'.ljust(24),
                  'sigma (S/m)'.ljust(24),
                  'k (W/m/K)'.ljust(24),
                  'VP (km/s)'.ljust(24),
                  'VS (km/s)'.ljust(24),
                  'QS'.ljust(24),
                  'KS (GPa)'.ljust(24),
                  'GS (GPa)'.ljust(24),
                  'Ppore (MPa)'.ljust(24),
                  'rhoMatrix (kg/m3)'.ljust(24),
                  'rhoPore (kg/m3)'.ljust(24),
                  'MLayer (kg)'.ljust(24),
                  'VLayer (m3)'.ljust(24),
                  'Htidal (W/m3)']
    # Print number of header lines early so we can skip the rest on read-in if we want to
    Params.nHeadLines = np.size(headerLines) + 3
    headerLines = np.insert(headerLines, 0, f'  nHeadLines = {Params.nHeadLines:d}')
    with open(Params.DataFiles.saveFile,'w') as f:
        f.write(Planet.label + '\n')
        f.write('\n  '.join(headerLines) + '\n')
        f.write(' '.join(colHeaders) + '\n')
        # Now print the columnar data
        for i in range(Planet.Steps.nTotal):
            f.write(' '.join([
                f'{Planet.P_MPa[i]:24.17e}',
                f'{Planet.T_K[i]:24.17e}',
                f'{Planet.r_m[i]:24.17e}',
                f'{Planet.phase[i]:8d}',
                f'{Planet.rho_kgm3[i]:24.17e}',
                f'{Planet.Cp_JkgK[i]:24.17e}',
                f'{Planet.alpha_pK[i]:24.17e}',
                f'{Planet.g_ms2[i]:24.17e}',
                f'{Planet.phi_frac[i]:24.17e}',
                f'{Planet.sigma_Sm[i]:24.17e}',
                f'{Planet.kTherm_WmK[i]:24.17e}',
                f'{Planet.Seismic.VP_kms[i]:24.17e}',
                f'{Planet.Seismic.VS_kms[i]:24.17e}',
                f'{Planet.Seismic.QS[i]:24.17e}',
                f'{Planet.Seismic.KS_GPa[i]:24.17e}',
                f'{Planet.Seismic.GS_GPa[i]:24.17e}',
                f'{Planet.Ppore_MPa[i]:24.17e}',
                f'{Planet.rhoMatrix_kgm3[i]:24.17e}',
                f'{Planet.rhoPore_kgm3[i]:24.17e}',
                f'{Planet.MLayer_kg[i]:24.17e}',
                f'{Planet.VLayer_m3[i]:24.17e}',
                f'{Planet.Htidal_Wm3[i]:24.17e}']) + '\n')

    # Write out data from core/mantle trade
    with open(Params.DataFiles.mantCoreFile, 'w') as f:
        f.write(' '.join(['RsilTrade (m)'.ljust(24),
                          'RcoreTrade (m)'.ljust(24),
                          'rhoSilTrade (kg/m3)']) + '\n')
        for i in range(np.size(Planet.Sil.Rtrade_m)):
            line = ' '.join([
                f'{Planet.Sil.Rtrade_m[i]:24.17e}',
                f'{Planet.Core.Rtrade_m[i]:24.17e}',
                f'{Planet.Sil.rhoTrade_kgm3[i]:24.17e}\n '
            ])
            f.write(line)

    log.info(f'Profile saved to file: {Params.DataFiles.saveFile}')
    return


def ReloadProfile(Planet, Params, fnameOverride=None):
    """ Reload previously saved PlanetProfile run from disk """
    if Planet is None:
        bodyname = fnameOverride.split('Profile')[0].split(os.sep)[-1]
        Planet = PlanetStruct(bodyname)
    if Params is None:
        Params = configParams

    if fnameOverride is not None:
        Params.DataFiles.saveFile = fnameOverride
        Params.DataFiles.mantCoreFile = f'{fnameOverride[:-4]}_mantleCore.txt'
        Params.DataFiles.mantPermFile = f'{fnameOverride[:-4]}_mantlePerm.txt'
        nSkip = len(os.path.join(Planet.bodyname, f'{Planet.name}Profile_'))
        Planet.saveLabel = fnameOverride[nSkip:-4]
    else:
        Params.DataFiles, Params.FigureFiles = SetupFilenames(Planet, Params)
    log.info(f'Reloading previously saved run from file: {Params.DataFiles.saveFile}')
    log.debug(f'Steps.n settings from PP{Planet.name}.py will be ignored.')
    if not isfile(Params.DataFiles.saveFile):
        raise ValueError(f'CALC_NEW is set to False in configPP.py but the reload file at {Params.DataFiles.saveFile} ' +
                         'was not found.\nRe-run with CALC_NEW set to True to generate the profile.')

    with open(Params.DataFiles.saveFile) as f:
        # Get legend label for differentiating runs
        Planet.label = f.readline().strip()
        # Get number of header lines to read in from (and skip for columnar data)
        Params.nHeadLines = int(f.readline().split('=')[-1])
        # Get MoI-included label for tradeoff plots
        Planet.tradeLabel = f.readline().split('=')[-1].strip()
        # Get whether iron core is modeled
        Planet.Do.Fe_CORE = bool(strtobool(f.readline().split('=')[-1].strip()))
        # Get silicate mantle Perple_X EOS file
        Planet.Sil.mantleEOS = f.readline().split('=')[-1].strip()
        # Get iron core Perple_X EOS file
        Planet.Core.coreEOS = f.readline().split('=')[-1].strip() 
        # Get dissolved salt supposed for ocean (present in filename, but this is intended for future-proofing when we move to a database lookup)
        Planet.Ocean.comp = f.readline().split('=')[-1].strip()
        # Get dissolved salt supposed for pore space
        Planet.Sil.poreComp = f.readline().split('=')[-1].strip()
        # Get float values from header
        Planet.Ocean.wOcean_ppt, Planet.Sil.wPore_ppt, Planet.Bulk.R_m, Planet.Bulk.M_kg, Planet.Bulk.Cmeasured, \
        Planet.Bulk.Cuncertainty, Planet.Bulk.Psurf_MPa, Planet.Bulk.Tsurf_K, Planet.qSurf_Wm2, Planet.qCon_Wm2, \
        Planet.Bulk.Tb_K, Planet.zb_km, Planet.zClath_m, Planet.D_km, Planet.Pb_MPa, Planet.PbI_MPa, \
        Planet.Ocean.deltaP, Planet.Mtot_kg, Planet.CMR2mean, Planet.CMR2less, Planet.CMR2more, Planet.Ocean.QfromMantle_W, \
        Planet.Ocean.rhoMean_kgm3, Planet.Sil.phiRockMax_frac, Planet.Sil.Qrad_Wkg, Planet.Sil.HtidalMean_Wm3, \
        Planet.Sil.Rmean_m, Planet.Sil.Rrange_m, Planet.Sil.rhoMean_kgm3, Planet.Core.Rmean_m, Planet.Core.Rrange_m, \
        Planet.Core.rhoMean_kgm3, Planet.MH2O_kg, Planet.Mrock_kg, Planet.Mcore_kg, Planet.Mice_kg, \
        Planet.Msalt_kg, Planet.MporeSalt_kg, Planet.Mocean_kg, Planet.Mfluid_kg, Planet.MporeFluid_kg, \
        Planet.Mclath_kg, Planet.MclathGas_kg, Planet.Ocean.sigmaMean_Sm, Planet.Sil.sigmaPoreMean_Sm, \
        Planet.Sil.sigmaPorousLayerMean_Sm, Planet.Tconv_K, Planet.etaConv_Pas, Planet.RaConvect, Planet.RaConvectIII, Planet.RaConvectV, \
        Planet.RaCrit, Planet.RaCritIII, Planet.RaCritV, Planet.eLid_m, Planet.eLidIII_m, Planet.eLidV_m, \
        Planet.Dconv_m, Planet.DconvIII_m, Planet.DconvV_m, Planet.deltaTBL_m, Planet.deltaTBLIII_m, Planet.deltaTBLV_m \
            = (float(f.readline().split('=')[-1]) for _ in range(63))
        # Note porosity flags
        Planet.Do.POROUS_ICE = bool(strtobool(f.readline().split('=')[-1].strip()))
        Planet.Do.POROUS_ROCK = Planet.Sil.phiRockMax_frac > 0
        # Get integer values from header (nSteps values)
        Planet.Steps.nClath, Planet.Steps.nIceI, \
        Planet.Steps.nIceIIILitho, Planet.Steps.nIceVLitho, \
        Planet.Steps.nHydro, Planet.Steps.nSil, Planet.Steps.nCore \
            = (int(f.readline().split('=')[-1]) for _ in range(7))

        if Planet.Ocean.comp == 'none':
            Planet.Do.NO_H2O = True
        if Planet.Do.NO_H2O:
            Planet.Bulk.qSurf_Wm2 = Planet.qSurf_Wm2 + 0.0

    Planet.Steps.nIbottom = Planet.Steps.nClath + Planet.Steps.nIceI
    Planet.Steps.nIIIbottom = Planet.Steps.nIbottom + Planet.Steps.nIceIIILitho
    Planet.Steps.nSurfIce = Planet.Steps.nIIIbottom + Planet.Steps.nIceVLitho
    Planet.Steps.nTotal = Planet.Steps.nHydro + Planet.Steps.nSil + Planet.Steps.nCore
    # Read in columnar data that follows header lines -- full-body
    Planet.P_MPa, Planet.T_K, Planet.r_m, Planet.phase, Planet.rho_kgm3, Planet.Cp_JkgK, Planet.alpha_pK, \
    Planet.g_ms2, Planet.phi_frac, Planet.sigma_Sm, Planet.kTherm_WmK, Planet.Seismic.VP_kms, Planet.Seismic.VS_kms,\
    Planet.Seismic.QS, Planet.Seismic.KS_GPa, Planet.Seismic.GS_GPa, Planet.Ppore_MPa, Planet.rhoMatrix_kgm3, \
    Planet.rhoPore_kgm3, Planet.MLayer_kg, Planet.VLayer_m3, Planet.Htidal_Wm3 \
        = np.loadtxt(Params.DataFiles.saveFile, skiprows=Params.nHeadLines, unpack=True)
    Planet.r_m = np.concatenate((Planet.r_m, [0]))
    Planet.z_m = Planet.Bulk.R_m - Planet.r_m
    Planet.phase = Planet.phase.astype(np.int_)

    # Read in data for core/mantle trade
    Planet.Sil.Rtrade_m, Planet.Core.Rtrade_m, Planet.Sil.rhoTrade_kgm3, \
        = np.loadtxt(Params.DataFiles.mantCoreFile, skiprows=1, unpack=True)

    return Planet, Params


def InitBayes(bodyname, fEnd):
    """ Load in a specific profile to use as in initial prior in Bayesian analysis. """
    Params = configParams
    # Prevent unnecessary slowdowns and disk space usage
    Params.NO_SAVEFILE = True
    Params.SKIP_PLOTS = True
    # Make sure CALC_NEW settings are as desired
    Params.CALC_NEW = True
    Params.CALC_NEW_INDUCT = True
    # Quiet messages unless we're debugging
    if bodyname != 'Test':
        log.setLevel(logging.WARN+5)

    # Fetch starting parameters
    fName = f'PP{bodyname}{fEnd}.py'
    expected = os.path.join(bodyname, fName)
    if not os.path.isfile(expected):
        default = os.path.join(_Defaults, bodyname, fName)
        if os.path.isfile(default):
            CopyCarefully(default, expected)
        else:
            log.warning(f'{expected} does not exist and no default was found at {default}.')
    Planet = importlib.import_module(expected[:-3].replace(os.sep, '.')).Planet
    return Planet, Params


def UpdateRun(Planet, Params, changes=None):
    """ Wrapper for editing the settings of Planet using a dict naming the changes.

        Args:
            changes (dict): Dict of {input: value} pairs to change in Planet. Options are listed in AssignPlanetVal.

    """
    if changes is None:
        log.warning('No changes passed to UpdateRun.')
        changes = {}

    for key, value in changes.items():
        Planet = AssignPlanetVal(Planet, key, value)

    Planet, Params = PlanetProfile(Planet, Params)
    return Planet, Params


def InductOgram(bodyname, Params):
    """ Run PlanetProfile models over a variety of settings to get induction
        responses for each input. 
    """
    if Params.CALC_NEW_INDUCT:
        log.info(f'Running calculations for induct-o-gram type {Params.Induct.inductOtype}.')
        Params.CALC_CONDUCT = True
        if bodyname[:4] == 'Test':
            loadname = bodyname + ''
            bodyname = 'Test'
            bodydir = os.path.join('PlanetProfile', 'Test')
        else:
            loadname = bodyname
            bodydir = bodyname

        fName = f'PP{loadname}InductOgram.py'
        expected = os.path.join(bodydir, fName)
        if not os.path.isfile(expected):
            default = os.path.join(_Defaults, bodydir, fName)
            if os.path.isfile(default):
                CopyCarefully(default, expected)
            else:
                log.warning(f'{expected} does not exist and no default was found at {default}.')
        Planet = importlib.import_module(expected[:-3].replace(os.sep, '.')).Planet
        DataFiles, FigureFiles = SetupFilenames(Planet, Params)
        tMarks = np.empty(0)
        tMarks = np.append(tMarks, time.time())
        settings = {}
        k = 1
        if Params.Induct.inductOtype == 'sigma':
            sigmaList = np.logspace(Params.Induct.sigmaMin[bodyname], Params.Induct.sigmaMax[bodyname], Params.Induct.nSigmaPts)
            Dlist = np.logspace(Params.Induct.Dmin[bodyname], Params.Induct.Dmax[bodyname], Params.Induct.nDpts)
            Params.nModels = Params.Induct.nSigmaPts * Params.Induct.nDpts
            settings['OceanChanges'] = {'sigmaMean_Sm': sigmaList}
            settings['PlanetChanges'] = {'D_km': Dlist}
            Planet.zb_km = Params.Induct.zbFixed_km[bodyname]
            Planet.Magnetic.rSigChange_m = np.array([0, Planet.Bulk.R_m - Planet.zb_km*1e3, Planet.Bulk.R_m])
            Planet.Magnetic.sigmaLayers_Sm = np.array([Planet.Ocean.sigmaIce_Sm['Ih'], 0, Planet.Sil.sigmaSil_Sm])
            PlanetGrid = np.empty((Params.Induct.nSigmaPts, Params.Induct.nDpts), dtype=object)
            for i, sigma_Sm in enumerate(sigmaList):
                for j, D_km in enumerate(Dlist):
                    Planet.Sil.rhoMean_kgm3 = Planet.Sil.rhoSilWithCore_kgm3
                    Planet.Sil.phiRockMax_frac = 0
                    Planet.Ocean.sigmaMean_Sm = sigma_Sm
                    Planet.Ocean.sigmaTop_Sm = sigma_Sm
                    Planet.Ocean.Tmean_K = Constants.T0
                    Planet.Magnetic.sigmaLayers_Sm[1] = sigma_Sm
                    Planet.D_km = D_km
                    Planet.Magnetic.rSigChange_m[0] = Planet.Bulk.R_m - 1e3 * (D_km + Planet.zb_km)
                    Planet.index = k
                    k += 1
                    PlanetGrid[i,j] = deepcopy(Planet)

        else:
            wList = np.logspace(Params.Induct.wMin[bodyname], Params.Induct.wMax[bodyname], Params.Induct.nwPts)
            settings['OceanChanges'] = {'wOcean_ppt': wList}
            if Params.Induct.inductOtype == 'Tb':
                Params.SKIP_INNER = True
                TbList = np.linspace(Params.Induct.TbMin[bodyname], Params.Induct.TbMax[bodyname], Params.Induct.nTbPts)
                settings['BulkChanges'] = {'Tb_K': TbList}
                Params.nModels = Params.Induct.nwPts * Params.Induct.nTbPts
                PlanetGrid = np.empty((Params.Induct.nwPts, Params.Induct.nTbPts), dtype=object)
                for i, w_ppt in enumerate(wList):
                    for j, Tb_K in enumerate(TbList):
                        Planet.Ocean.wOcean_ppt = w_ppt
                        Planet.Bulk.Tb_K = Tb_K
                        Planet.index = k
                        k += 1
                        PlanetGrid[i,j] = deepcopy(Planet)
            elif Params.Induct.inductOtype == 'phi':
                phiList = np.logspace(Params.Induct.phiMin[bodyname], Params.Induct.phiMax[bodyname], Params.Induct.nphiPts)
                settings['SilChanges'] = {'phiRockMax_frac': phiList}
                Params.nModels = Params.Induct.nwPts * Params.Induct.nphiPts
                PlanetGrid = np.empty((Params.Induct.nwPts, Params.Induct.nphiPts), dtype=object)
                Params.CONSTANT_DENSITY_INNER = False
                Params.SKIP_INNER = False
                for i, w_ppt in enumerate(wList):
                    for j, phiRockMax_frac in enumerate(phiList):
                        Planet.Ocean.wOcean_ppt = w_ppt
                        Planet.Sil.phiRockMax_frac = phiRockMax_frac
                        Planet.index = k
                        k += 1
                        PlanetGrid[i,j] = deepcopy(Planet)
            elif Params.Induct.inductOtype == 'rho':
                Params.SKIP_INNER = True
                rhoList = np.linspace(Params.Induct.rhoMin[bodyname], Params.Induct.rhoMax[bodyname], Params.Induct.nrhoPts)
                settings['SilChanges'] = {'rhoSilWithCore_kgm3': rhoList}
                Params.nModels = Params.Induct.nwPts * Params.Induct.nrhoPts
                PlanetGrid = np.empty((Params.Induct.nwPts, Params.Induct.nrhoPts), dtype=object)
                for i, w_ppt in enumerate(wList):
                    for j, rhoSilWithCore_kgm3 in enumerate(rhoList):
                        Planet.Ocean.wOcean_ppt = w_ppt
                        Planet.Sil.rhoSilWithCore_kgm3 = rhoSilWithCore_kgm3
                        Planet.index = k
                        k += 1
                        PlanetGrid[i,j] = deepcopy(Planet)
            else:
                raise ValueError(f'inductOtype {Params.Induct.inductOtype} behavior not defined.')

        tMarks = np.append(tMarks, time.time())
        log.info('PlanetGrid constructed. Calculating induction responses.')
        Params.INDUCTOGRAM_IN_PROGRESS = True
        PlanetGrid = ParPlanet(PlanetGrid, Params)
        tMarks = np.append(tMarks, time.time())
        dt = tMarks[-1] - tMarks[-2]
        log.info(f'Parallel run elapsed time: {dt:.1f} s.')

        # Organize data into a format that can be plotted/saved for plotting
        Bex_nT, Bey_nT, Bez_nT = Benm2absBexyz(PlanetGrid[0,0].Magnetic.Benm_nT)
        nPeaks = np.size(Bex_nT)
        Induction = InductionResults
        Induction.bodyname = bodyname
        Induction.yName = Params.Induct.inductOtype
        Induction.Texc_hr = Mag.Texc_hr[bodyname]
        Induction.Amp = np.array([[[Planeti.Magnetic.Amp[i] for Planeti in line] for line in PlanetGrid] for i in range(nPeaks)])
        Induction.phase = np.array([[[Planeti.Magnetic.phase[i] for Planeti in line] for line in PlanetGrid] for i in range(nPeaks)])
        Induction.Bix_nT = np.array([Induction.Amp[i, ...] * Bex_nT[i] for i in range(nPeaks)])
        Induction.Biy_nT = np.array([Induction.Amp[i, ...] * Bey_nT[i] for i in range(nPeaks)])
        Induction.Biz_nT = np.array([Induction.Amp[i, ...] * Bez_nT[i] for i in range(nPeaks)])
        Induction.wOcean_ppt = np.array([[Planeti.Ocean.wOcean_ppt for Planeti in line] for line in PlanetGrid])
        Induction.oceanComp = np.array([[Planeti.Ocean.comp for Planeti in line] for line in PlanetGrid])
        Induction.Tb_K = np.array([[Planeti.Bulk.Tb_K for Planeti in line] for line in PlanetGrid])
        Induction.Tmean_K = np.array([[Planeti.Ocean.Tmean_K for Planeti in line] for line in PlanetGrid])
        Induction.rhoSilMean_kgm3 = np.array([[Planeti.Sil.rhoMean_kgm3 for Planeti in line] for line in PlanetGrid])
        Induction.phiSilMax_frac = np.array([[Planeti.Sil.phiRockMax_frac for Planeti in line] for line in PlanetGrid])
        Induction.sigmaMean_Sm = np.array([[Planeti.Ocean.sigmaMean_Sm for Planeti in line] for line in PlanetGrid])
        Induction.sigmaTop_Sm = np.array([[Planeti.Ocean.sigmaTop_Sm for Planeti in line] for line in PlanetGrid])
        Induction.D_km = np.array([[Planeti.D_km for Planeti in line] for line in PlanetGrid])
        Induction.zb_km = np.array([[Planeti.zb_km for Planeti in line] for line in PlanetGrid])
        Induction.R_m = np.array([[Planeti.Bulk.R_m for Planeti in line] for line in PlanetGrid])
        Induction.rBds_m = np.array([[Planeti.Magnetic.rSigChange_m for Planeti in line] for line in PlanetGrid])
        Induction.sigmaLayers_Sm = np.array([[Planeti.Magnetic.sigmaLayers_Sm for Planeti in line] for line in PlanetGrid])

        Params.DataFiles = DataFiles
        Params.FigureFiles = FigureFiles
        WriteInductOgram(Induction, Params)
        Induction.SetAxes(Params.Induct.inductOtype)
        Induction.SetComps(Params.Induct.inductOtype)
    else:
        log.info(f'Reloading induct-o-gram type {Params.Induct.inductOtype}.')
        Induction, Params = ReloadInductOgram(bodyname, Params)

    return Induction, Params


def WriteInductOgram(Induction, Params):
    """ Organize Induction results from an induct-o-gram run into a dict
        and print to a .mat file.
    """

    saveDict = {
        'bodyname': Induction.bodyname,
        'yName': Induction.yName,
        'Texc_hr_keys': [key for key in Induction.Texc_hr.keys()],
        'Texc_hr_values': [value for value in Induction.Texc_hr.values()],
        'Amp': Induction.Amp,
        'phase': Induction.phase,
        'Bix_nT': Induction.Bix_nT,
        'Biy_nT': Induction.Biy_nT,
        'Biz_nT': Induction.Biz_nT,
        'w_ppt': Induction.wOcean_ppt,
        'oceanComp': Induction.oceanComp,
        'Tb_K': Induction.Tb_K,
        'Tmean_K': Induction.Tmean_K,
        'rhoSilMean_kgm3': Induction.rhoSilMean_kgm3,
        'phiSilMax_frac': Induction.phiSilMax_frac,
        'sigmaMean_Sm': Induction.sigmaMean_Sm,
        'sigmaTop_Sm': Induction.sigmaTop_Sm,
        'D_km': Induction.D_km,
        'zb_km': Induction.zb_km,
        'R_m': Induction.R_m,
        'rBds_m': Induction.rBds_m,
        'sigmaLayers_Sm': Induction.sigmaLayers_Sm
    }
    savemat(Params.DataFiles.inductOgramFile, saveDict)
    log.info(f'Saved induct-o-gram {Params.DataFiles.inductOgramFile} to disk.')

    return


def ReloadInductOgram(bodyname, Params, fNameOverride=None):
    """ Reload a previously run induct-o-gram from disk.
    """

    if fNameOverride is None:
        if bodyname[:4] == 'Test':
            loadname = bodyname + ''
            bodydir = _TestImport
        else:
            loadname = bodyname
            bodydir = bodyname
        Planet = importlib.import_module(f'{bodydir}.PP{loadname}InductOgram').Planet
        Params.DataFiles, Params.FigureFiles = SetupFilenames(Planet, Params)
        loadFile = Params.DataFiles.inductOgramFile
    else:
        loadFile = fNameOverride
    if os.path.isfile(loadFile):
        reload = loadmat(loadFile)
    else:
        raise FileNotFoundError(f'Attempted to reload inductogram, but {loadFile} was not found. ' +
                                'Re-run with Params.CALC_NEW_INDUCT = True in configPP.py.')

    Induction = InductionResults
    Induction.bodyname = reload['bodyname'][0]
    Induction.yName = reload['yName'][0]
    Induction.Texc_hr = {key.strip(): value for key, value in zip(reload['Texc_hr_keys'], reload['Texc_hr_values'][0])}
    Induction.Amp = reload['Amp']
    Induction.phase = reload['phase']
    Induction.Bix_nT = reload['Bix_nT']
    Induction.Biy_nT = reload['Biy_nT']
    Induction.Biz_nT = reload['Biz_nT']
    Induction.wOcean_ppt = reload['w_ppt']
    Induction.oceanComp = reload['oceanComp']
    Induction.Tb_K = reload['Tb_K']
    Induction.Tmean_K = reload['Tmean_K']
    Induction.rhoSilMean_kgm3 = reload['rhoSilMean_kgm3']
    Induction.phiSilMax_frac = reload['phiSilMax_frac']
    Induction.sigmaMean_Sm = reload['sigmaMean_Sm']
    Induction.sigmaTop_Sm = reload['sigmaTop_Sm']
    Induction.D_km = reload['D_km']
    Induction.zb_km = reload['zb_km']
    Induction.R_m = reload['R_m']
    Induction.rBds_m = reload['rBds_m']
    Induction.sigmaLayers_Sm = reload['sigmaLayers_Sm']

    Induction.SetAxes(Params.Induct.inductOtype)
    Induction.SetComps(Params.Induct.inductOtype)

    return Induction, Params


def ParPlanet(PlanetList, Params):
    """ Run a list of PlanetProfile models over arrays of run settings.

        Args:
            PlanetList (PlanetStruct, shape N, NxM, Nx...): List of Planet objects
                over which to run in parallel. 
    """
    if Params.logParallel > logging.INFO:
        log.info('Quieting messages to avoid spam in gridded run.')
    saveLevel = log.getEffectiveLevel() + 0
    log.setLevel(Params.logParallel)
    dims = np.shape(PlanetList)
    nParDims = np.size(dims)
    if nParDims == 1 and np.size(PlanetList) == 1:
        PlanetList[0] = PlanetProfile(PlanetList[0], Params)
    else:
        if Params.INDUCTOGRAM_IN_PROGRESS and Params.Induct.inductOtype == 'sigma':
            PlanetList = GridPlanetProfileFunc(MagneticInduction, PlanetList, Params)
        else:
            PlanetList = GridPlanetProfileFunc(PlanetProfile, PlanetList, Params)

    # Return settings to what they were before we entered here
    log.setLevel(saveLevel)
    Params.INDUCTOGRAM_IN_PROGRESS = False

    return PlanetList


def ParPlanetExplore(Planet, Params, xList, yList):
    """ Run a parameter exploration over arrays of run settings, starting from a base Planet object.

        Args:
            xList, yList (float, shape nx or ny): Lists of values to use for x,y variables
    """
    # Construct PlanetGrid to use for exploration
    PlanetGrid = np.empty((Params.Explore.nx, Params.Explore.ny), dtype=object)
    nTot = Params.Explore.nx * Params.Explore.ny
    IND_SKIP_SAVE = Params.SKIP_INDUCTION and True
    k = 0
    if Params.logParallel > logging.INFO:
        log.info('Quieting messages to avoid spam in gridded run.')
    saveLevel = log.getEffectiveLevel() + 0
    log.setLevel(Params.logParallel)

    if (Params.Explore.exploreType[Params.Explore.xName] == 'hydro' and
        Params.Explore.exploreType[Params.Explore.yName] == 'hydro'):
        # In this case, we have to run the whole gamut for each model.
        for i, xVal in enumerate(xList):
            for j, yVal in enumerate(yList):
                k += 1
                Planet = AssignPlanetVal(Planet, Params.Explore.xName, xVal)
                Planet = AssignPlanetVal(Planet, Params.Explore.yName, yVal)
                Planet.index = k
                PlanetGrid[i,j] = deepcopy(Planet)
        log.info('PlanetGrid constructed. Calculating exploration responses.')
        Params.nModels = nTot
        Params.tStart_s = time.time()
        PlanetGrid = GridPlanetProfileFunc(PlanetProfile, PlanetGrid, Params)
    else:
        if (Params.Explore.exploreType[Params.Explore.xName] == 'ionos' and
            Params.Explore.exploreType[Params.Explore.yName] == 'ionos'):
            # This is the simplest case. We can run one interior model and just do
            # the MagneticInduction part again for the rest.
            log.info('Running common interior model to iterate on for ionosphere-only exploration.')
            Planet = AssignPlanetVal(Planet, Params.Explore.xName, xList[0])
            Planet = AssignPlanetVal(Planet, Params.Explore.yName, yList[0])
            Params.SKIP_INDUCTION = True
            PlanetGrid[0,0], Params = PlanetProfile(Planet, Params)
            Params.SKIP_INDUCTION = IND_SKIP_SAVE
            Params.nModels = nTot
            Params.tStart_s = time.time()
            log.info('Copying common interior model to entire explore grid.')
            for i, xVal in enumerate(xList):
                for j, yVal in enumerate(yList):
                    k += 1
                    if not (i == 0 and j == 0):
                        PlanetGrid[i,j] = deepcopy(PlanetGrid[0,0])
                        PlanetGrid[i,j] = AssignPlanetVal(PlanetGrid[i,j], Params.Explore.xName, xVal)
                        PlanetGrid[i,j] = AssignPlanetVal(PlanetGrid[i,j], Params.Explore.yName, yVal)
                        PlanetGrid[i,j].index = k
                        PrintCompletion(PlanetGrid[i,j], Params)
            log.info('PlanetGrid constructed. Calculating exploration responses.')
            Params.tStart_s = time.time()
            PlanetGrid = GridPlanetProfileFunc(InductionOnly, PlanetGrid, Params)

        elif (Params.Explore.exploreType[Params.Explore.xName] == 'ionos' or
              Params.Explore.exploreType[Params.Explore.yName] == 'ionos'):
            # In this case, we have one ionosphere exploration. First, check which dim:
            if Params.Explore.exploreType[Params.Explore.xName] == 'ionos':
                # Now figure out whether we can reuse the hydrosphere for the other dim:
                if Params.Explore.exploreType[Params.Explore.yName] == 'inner':
                    log.info('Running common hydrosphere model to iterate on for inner+ionosphere exploration.')
                    PlanetGrid[0,0] = AssignPlanetVal(Planet, Params.Explore.xName, xList[0])
                    PlanetGrid[0,0] = AssignPlanetVal(PlanetGrid[0,0], Params.Explore.yName, yList[0])
                    PlanetGrid[0,0].index = 1
                    PlanetGrid[0,0], Params = HydroOnly(PlanetGrid[0,0], Params)
                    Params.nModels = Params.Explore.ny
                    Params.tStart_s = time.time()
                    log.info('Copying common hydrosphere model to grid row.')
                    for j, yVal in enumerate(yList):
                        k += 1
                        if j != 0:
                            PlanetGrid[0,j] = deepcopy(PlanetGrid[0,0])
                            PlanetGrid[0,j] = AssignPlanetVal(PlanetGrid[0,j], Params.Explore.yName, yVal)
                            PlanetGrid[0,j].index = k
                            PrintCompletion(PlanetGrid[0,j], Params)
                    log.info('PlanetGrid row constructed. Calculating exploration responses to propagate for ionosphere exploration.')
                    Params.tStart_s = time.time()
                    Params.SKIP_INDUCTION = True
                    PlanetGrid[0,:] = GridPlanetProfileFunc(InteriorEtc, PlanetGrid[0,:], Params)
                    Params.SKIP_INDUCTION = IND_SKIP_SAVE
                    k = 0
                    Params.nModels = nTot
                    Params.tStart_s = time.time()
                    log.info('Copying hydrosphere grid row to remaining grid.')
                    for i, xVal in enumerate(xList):
                        for j, _ in enumerate(yList):
                            k += 1
                            if i != 0:
                                PlanetGrid[i,j] = deepcopy(PlanetGrid[0,j])
                                PlanetGrid[i,j] = AssignPlanetVal(PlanetGrid[i,j], Params.Explore.xName, xVal)
                                PlanetGrid[i,j].index = k
                                PrintCompletion(PlanetGrid[i,j], Params)
                    log.info('Running interior model row to iterate on for ionosphere exploration.')
                    Params.tStart_s = time.time()
                    PlanetGrid = GridPlanetProfileFunc(InductionOnly, PlanetGrid, Params)

                else:
                    # In this case, we need to run full PlanetProfile interior calcs for the non-ionos row.
                    Planet = AssignPlanetVal(Planet, Params.Explore.xName, xList[0])
                    for j, yVal in enumerate(yList):
                        k += 1
                        Planet = AssignPlanetVal(Planet, Params.Explore.yName, yVal)
                        Planet.index = k
                        PlanetGrid[0,j] = deepcopy(Planet)
                    log.info('PlanetGrid row constructed. Calculating exploration responses.')
                    Params.tStart_s = time.time()
                    Params.SKIP_INDUCTION = True
                    PlanetGrid[0,:] = GridPlanetProfileFunc(PlanetProfile, PlanetGrid[0,:], Params)
                    Params.SKIP_INDUCTION = IND_SKIP_SAVE
                    k = 0
                    Params.nModels = nTot
                    Params.tStart_s = time.time()
                    for i, xVal in enumerate(xList):
                        for j, _ in enumerate(yList):
                            k += 1
                            if i != 0:
                                PlanetGrid[i,j] = deepcopy(PlanetGrid[0,j])
                                PlanetGrid[i,j] = AssignPlanetVal(PlanetGrid[i,j], Params.Explore.xName, xVal)
                                PlanetGrid[i,j].index = k
                                PrintCompletion(PlanetGrid[i,j], Params)
                    log.info('Running interior model row to iterate on for ionosphere exploration.')
                    Params.tStart_s = time.time()
                    PlanetGrid = GridPlanetProfileFunc(InductionOnly, PlanetGrid, Params)

            else:
                # Repeat of above case, but now we have ionos calcs on the y axis.
                # First, check if we can reuse hydrosphere calcs on x axis:
                if Params.Explore.exploreType[Params.Explore.xName] == 'inner':
                    log.info('Running common hydrosphere model to iterate on for inner+ionosphere exploration.')
                    PlanetGrid[0,0] = AssignPlanetVal(Planet, Params.Explore.xName, xList[0])
                    PlanetGrid[0,0] = AssignPlanetVal(PlanetGrid[0,0], Params.Explore.yName, yList[0])
                    PlanetGrid[0,0].index = 1
                    PlanetGrid[0,0], Params = HydroOnly(PlanetGrid[0,0], Params)
                    log.info('Copying common hydrosphere to grid row.')
                    Params.nModels = Params.Explore.nx
                    Params.tStart_s = time.time()
                    for i, xVal in enumerate(xList):
                        k += 1
                        if i != 0:
                            PlanetGrid[i,0] = deepcopy(PlanetGrid[0,0])
                            PlanetGrid[i,0] = AssignPlanetVal(PlanetGrid[i,0], Params.Explore.xName, xVal)
                            PlanetGrid[i,0].index = k
                            PrintCompletion(PlanetGrid[i,0], Params)
                    log.info('PlanetGrid row constructed. Calculating exploration responses to propagate for ionosphere exploration.')
                    Params.tStart_s = time.time()
                    Params.SKIP_INDUCTION = True
                    PlanetGrid[:,0] = GridPlanetProfileFunc(InteriorEtc, PlanetGrid[:,0], Params)
                    Params.SKIP_INDUCTION = IND_SKIP_SAVE
                    k = 0
                    Params.nModels = nTot
                    Params.tStart_s = time.time()
                    log.info('Copying interior model row to entire explore grid.')
                    for i, _ in enumerate(xList):
                        for j, yVal in enumerate(yList):
                            k += 1
                            if j != 0:
                                PlanetGrid[i,j] = deepcopy(PlanetGrid[i,0])
                                PlanetGrid[i,j] = AssignPlanetVal(PlanetGrid[i,j], Params.Explore.yName, yVal)
                                PlanetGrid[i,j].index = k
                                PrintCompletion(PlanetGrid[i,j], Params)
                    log.info('Running interior model row to iterate on for ionosphere exploration.')
                    Params.tStart_s = time.time()
                    PlanetGrid = GridPlanetProfileFunc(InductionOnly, PlanetGrid, Params)

                else:
                    # In this case, we need to run full PlanetProfile interior calcs for the non-ionos row.
                    Planet = AssignPlanetVal(Planet, Params.Explore.yName, yList[0])
                    for i, xVal in enumerate(xList):
                        k += 1
                        Planet = AssignPlanetVal(Planet, Params.Explore.xName, xVal)
                        Planet.index = k
                        PlanetGrid[i,0] = deepcopy(Planet)
                    log.info('PlanetGrid row constructed. Calculating exploration responses.')
                    Params.nModels = Params.Explore.nx
                    Params.tStart_s = time.time()
                    Params.SKIP_INDUCTION = True
                    PlanetGrid[:,0] = GridPlanetProfileFunc(PlanetProfile, PlanetGrid[:,0], Params)
                    Params.SKIP_INDUCTION = IND_SKIP_SAVE
                    k = 0
                    Params.nModels = nTot
                    Params.tStart_s = time.time()
                    log.info('Copying interior calcs to entire explore grid.')
                    for i, _ in enumerate(xList):
                        for j, yVal in enumerate(yList):
                            k += 1
                            if j != 0:
                                PlanetGrid[i,j] = deepcopy(PlanetGrid[i,0])
                                PlanetGrid[i,j] = AssignPlanetVal(PlanetGrid[i,j], Params.Explore.yName, yVal)
                                PlanetGrid[i,j].index = k
                                PrintCompletion(PlanetGrid[i,j], Params)
                    log.info('Running interior model row to iterate on for ionosphere exploration.')
                    Params.tStart_s = time.time()
                    PlanetGrid = GridPlanetProfileFunc(InductionOnly, PlanetGrid, Params)

        else:
            # Finally, we have a combination of hydro and inner, or both inner.
            if (Params.Explore.exploreType[Params.Explore.xName] == 'inner' and
                Params.Explore.exploreType[Params.Explore.yName] == 'inner'):
                # In this case, we need to run 1 hydro-only model, then propagate to the whole grid.
                log.info('Running common hydrosphere model to iterate on for interior-only exploration.')
                PlanetGrid[0,0] = AssignPlanetVal(Planet, Params.Explore.xName, xList[0])
                PlanetGrid[0,0] = AssignPlanetVal(PlanetGrid[0,0], Params.Explore.yName, yList[0])
                PlanetGrid[0,0].index = 1
                PlanetGrid[0,0], Params = HydroOnly(Planet, Params)
                log.info('Copying hydrosphere calcs to entire explore grid.')
                Params.nModels = nTot
                Params.tStart_s = time.time()
                for i, xVal in enumerate(xList):
                    for j, yVal in enumerate(yList):
                        k += 1
                        if not (i == 0 and j == 0):
                            PlanetGrid[i,j] = deepcopy(PlanetGrid[0,0])
                            PlanetGrid[i,j] = AssignPlanetVal(PlanetGrid[i,j], Params.Explore.xName, xVal)
                            PlanetGrid[i,j] = AssignPlanetVal(PlanetGrid[i,j], Params.Explore.yName, yVal)
                            PlanetGrid[i,j].index = k
                            PrintCompletion(PlanetGrid[i,j], Params)
                log.info('PlanetGrid constructed. Calculating exploration responses.')
                Params.tStart_s = time.time()
                PlanetGrid = GridPlanetProfileFunc(InteriorEtc, PlanetGrid, Params)

            else:
                # Now, we finally have the case that we have one hydro and one inner.
                # First, find out which one is on the x axis:
                if Params.Explore.exploreType[Params.Explore.xName] == 'hydro':
                    # In this case, we need to first run hydrosphere calcs for the non-inner row.
                    Planet = AssignPlanetVal(Planet, Params.Explore.yName, yList[0])
                    for i, xVal in enumerate(xList):
                        k += 1
                        Planet = AssignPlanetVal(Planet, Params.Explore.xName, xVal)
                        Planet.index = k
                        PlanetGrid[i,0] = deepcopy(Planet)
                    log.info('PlanetGrid row constructed. Calculating exploration responses.')
                    Params.nModels = Params.Explore.nx
                    Params.tStart_s = time.time()
                    PlanetGrid[:,0] = GridPlanetProfileFunc(HydroOnly, PlanetGrid[:,0], Params)
                    k = 0
                    Params.nModels = nTot
                    Params.tStart_s = time.time()
                    log.info('Copying hydrosphere row to remaining explore grid.')
                    for i, _ in enumerate(xList):
                        for j, yVal in enumerate(yList):
                            k += 1
                            if j != 0:
                                PlanetGrid[i,j] = deepcopy(PlanetGrid[i,0])
                                PlanetGrid[i,j] = AssignPlanetVal(PlanetGrid[i,j], Params.Explore.yName, yVal)
                                PlanetGrid[i,j].index = k
                                PrintCompletion(PlanetGrid[i,j], Params)
                    log.info('Running hydrosphere model row to iterate on for interior exploration.')
                    Params.tStart_s = time.time()
                    PlanetGrid = GridPlanetProfileFunc(InteriorEtc, PlanetGrid, Params)

                else:
                    # Lastly, do the same but for hydro on the y axis instead:
                    Planet = AssignPlanetVal(Planet, Params.Explore.xName, xList[0])
                    for j, yVal in enumerate(yList):
                        k += 1
                        Planet = AssignPlanetVal(Planet, Params.Explore.yName, yVal)
                        Planet.index = k
                        PlanetGrid[0,j] = deepcopy(Planet)
                    log.info('PlanetGrid row constructed. Calculating exploration responses.')
                    Params.nModels = Params.Explore.ny
                    Params.tStart_s = time.time()
                    PlanetGrid[0,:] = GridPlanetProfileFunc(HydroOnly, PlanetGrid[0,:], Params)
                    k = 0
                    Params.nModels = nTot
                    Params.tStart_s = time.time()
                    log.info('Copying hydrosphere row to remaining explore grid.')
                    for i, xVal in enumerate(xList):
                        for j, _ in enumerate(yList):
                            k += 1
                            if i != 0:
                                PlanetGrid[i,j] = deepcopy(PlanetGrid[0,j])
                                PlanetGrid[i,j] = AssignPlanetVal(PlanetGrid[i,j], Params.Explore.xName, xVal)
                                PlanetGrid[i,j].index = k
                                PrintCompletion(PlanetGrid[i,j], Params)
                    log.info('Running hydrosphere model row to iterate on for interior exploration.')
                    Params.tStart_s = time.time()
                    PlanetGrid = GridPlanetProfileFunc(InteriorEtc, PlanetGrid, Params)

    # Return log settings to what they were before we entered here
    log.setLevel(saveLevel)

    return PlanetGrid


def GridPlanetProfileFunc(FuncName, PlanetGrid, Params):
    """ Wrapper for (optionally) parallel run of multiple Planet objects through the
        funcName function.
    """
    PlanetList1D = np.reshape(PlanetGrid, -1)
    if Params.DO_PARALLEL:
        # Prevent slowdowns from competing process spawning when #cores > #jobs
        nCores = np.min([Params.maxCores, np.product(np.shape(PlanetList1D)), Params.threadLimit])
        pool = mtpFork.Pool(nCores)
        parResult = [pool.apply_async(FuncName, (deepcopy(Planet),
                                                      deepcopy(Params))) for Planet in PlanetList1D]
        pool.close()
        pool.join()

        for i, result in enumerate(parResult):
            PlanetList1D[i] = result.get()[0]
    else:
        log.profile('Running grid without parallel processing. This may take some time.')
        PlanetList1D = np.array([FuncName(deepcopy(Planet), deepcopy(Params)) for Planet in PlanetList1D])[:, 0]

    PlanetGrid = np.reshape(PlanetList1D, np.shape(PlanetGrid))

    return PlanetGrid


def ExploreOgram(bodyname, Params):
    """ Run PlanetProfile models over a variety of settings to get interior
        properties for each input.
    """
    if Params.CALC_NEW:
        log.info(f'Running {Params.Explore.xName} x {Params.Explore.yName} explore-o-gram for {bodyname}.')
        if bodyname[:4] == 'Test':
            loadname = bodyname + ''
            bodyname = 'Test'
            bodydir = os.path.join('PlanetProfile', 'Test')
        else:
            loadname = bodyname
            bodydir = bodyname

        fName = f'PP{loadname}Explore.py'
        expected = os.path.join(bodydir, fName)
        if not os.path.isfile(expected):
            default = os.path.join(_Defaults, bodydir, fName)
            if os.path.isfile(default):
                CopyCarefully(default, expected)
            else:
                log.warning(f'{expected} does not exist and no default was found at {default}.')
        Planet = importlib.import_module(expected[:-3].replace(os.sep, '.')).Planet
        tMarks = np.empty(0)
        tMarks = np.append(tMarks, time.time())

        Exploration = ExplorationResults
        Exploration.xName = Params.Explore.xName
        Exploration.yName = Params.Explore.yName
        Exploration.zName = Params.Explore.zName
        DataFiles, FigureFiles = SetupFilenames(Planet, Params, exploreAppend=f'{Exploration.xName}{Exploration.yName}',
                                                figExploreAppend=Exploration.zName)
        if bodyname == 'Test':
            Params.Explore.nx = 5
            Params.Explore.ny = 5
        xList = np.linspace(Params.Explore.xRange[0], Params.Explore.xRange[1], Params.Explore.nx)
        yList = np.linspace(Params.Explore.yRange[0], Params.Explore.yRange[1], Params.Explore.ny)
        Params.nModels = Params.Explore.nx * Params.Explore.ny
        if not Params.SKIP_INNER:
            log.warning('Running explore-o-gram with interior calculations, which will be slow.')
        Params.NO_SAVEFILE = True
        Params.ALLOW_BROKEN_MODELS = True

        tMarks = np.append(tMarks, time.time())
        PlanetGrid = ParPlanetExplore(Planet, Params, xList, yList)
        tMarks = np.append(tMarks, time.time())
        dt = tMarks[-1] - tMarks[-2]
        log.info(f'Parallel run elapsed time: {dt:.1f} s.')

        # Organize data into a format that can be plotted/saved for plotting
        Exploration.bodyname = bodyname
        Exploration.wOcean_ppt = np.array([[Planeti.Ocean.wOcean_ppt for Planeti in line] for line in PlanetGrid])
        Exploration.oceanComp = np.array([[Planeti.Ocean.comp for Planeti in line] for line in PlanetGrid])
        Exploration.R_m = np.array([[Planeti.Bulk.R_m for Planeti in line] for line in PlanetGrid])
        Exploration.Tb_K = np.array([[Planeti.Bulk.Tb_K for Planeti in line] for line in PlanetGrid])
        Exploration.xFeS = np.array([[Planeti.Core.xFeS for Planeti in line] for line in PlanetGrid])
        Exploration.rhoSilInput_kgm3 = np.array([[Planeti.Sil.rhoSilWithCore_kgm3 for Planeti in line] for line in PlanetGrid])
        Exploration.silPhi_frac = np.array([[Planeti.Sil.phiRockMax_frac for Planeti in line] for line in PlanetGrid])
        Exploration.icePhi_frac = np.array([[Planeti.Ocean.phiMax_frac['Ih'] for Planeti in line] for line in PlanetGrid])
        Exploration.silPclosure_MPa = np.array([[Planeti.Sil.Pclosure_MPa for Planeti in line] for line in PlanetGrid])
        Exploration.icePclosure_MPa = np.array([[Planeti.Ocean.Pclosure_MPa['Ih'] for Planeti in line] for line in PlanetGrid])
        Exploration.ionosTop_km = np.array([[Planeti.Magnetic.ionosBounds_m[-1]/1e3 for Planeti in line] for line in PlanetGrid])
        Exploration.sigmaIonos_Sm = np.array([[Planeti.Magnetic.sigmaIonosPedersen_Sm[-1] for Planeti in line] for line in PlanetGrid])
        Exploration.Htidal_Wm3 = np.array([[Planeti.Sil.Htidal_Wm3 for Planeti in line] for line in PlanetGrid])
        Exploration.Qrad_Wkg = np.array([[Planeti.Sil.Qrad_Wkg for Planeti in line] for line in PlanetGrid])
        Exploration.rhoSilMean_kgm3 = np.array([[Planeti.Sil.rhoMean_kgm3 for Planeti in line] for line in PlanetGrid])
        Exploration.rhoCoreMean_kgm3 = np.array([[Planeti.Core.rhoMean_kgm3 for Planeti in line] for line in PlanetGrid])
        Exploration.sigmaMean_Sm = np.array([[Planeti.Ocean.sigmaMean_Sm for Planeti in line] for line in PlanetGrid])
        Exploration.sigmaTop_Sm = np.array([[Planeti.Ocean.sigmaTop_Sm for Planeti in line] for line in PlanetGrid])
        Exploration.Tmean_K = np.array([[Planeti.Ocean.Tmean_K for Planeti in line] for line in PlanetGrid])
        Exploration.D_km = np.array([[Planeti.D_km for Planeti in line] for line in PlanetGrid])
        Exploration.zb_km = np.array([[Planeti.zb_km for Planeti in line] for line in PlanetGrid])
        Exploration.Rcore_km = np.array([[Planeti.Core.Rmean_m/1e3 for Planeti in line] for line in PlanetGrid])
        Exploration.qSurf_Wm2 = np.array([[Planeti.qSurf_Wm2 for Planeti in line] for line in PlanetGrid])

        # Ensure everything is set so things will play nicely with .mat saving and plotting functions
        nans = np.nan * Exploration.R_m
        for name, attr in Exploration.__dict__.items():
            if attr is None:
                setattr(Exploration, name, nans)

        Params.DataFiles = DataFiles
        Params.FigureFiles = FigureFiles
        WriteExploreOgram(Exploration, Params)
    else:
        log.info(f'Reloading explore-o-gram for {bodyname}.')
        Exploration, Params = ReloadExploreOgram(bodyname, Params)

    return Exploration, Params


def AssignPlanetVal(Planet, name, val):
    """ Set values in Planet object based on descriptive key. Variable descriptions:
            R_m: Body surface radius in m in Planet.Bulk.R_m
            compOcean: Ocean composition string in Planet.Ocean.comp
            compSil: Silicate composition to use from available Perplex output files
            compFe: Iron core composition to use from available Perplex output files
            wOcean_ppt: Salinity in Planet.Ocean.wOcean_ppt
            Tb_K: Ocean bottom temperature in K in Planet.Bulk.Tb_K
            xFeS: Core FeS / Fe mixing ratio in Planet.Core.xFeS
            rhoSilInput_kgm3: Fixed density in silicate layers in Planet.Sil.rhoSilWithCore_kgm3 (for use with Planet.Do.CONSTANT_INNER_DENSITY)
            silPhi_frac: Vacuum-extrapolated porosity in silicates in Planet.Sil.phiRockMax_frac
            silPclosure_MPa: Pore closure pressure in silicates in Planet.Sil.Pclosure_MPa
            icePhi_frac: Vacuum porosity in ices in Planet.Ocean.phiMax_frac
            icePclosure_MPa: Pore closure pressure in ices in Planet.Ocean.Pclosure_MPa
            Htidal_Wm3: Fixed tidal heating in silicates in Planet.Sil.Htidal_Wm3
            Qrad_Wkg: Fixed radiogenic heating in silicates in Planet.Sil.Qrad_Wkg
            qSurf_Wm2: Surface heat flux for waterless bodies in Planet.Bulk.qSurf_Wm2
            ionosTop_km: Ionosphere upper limit altitude above the surface in km, used in Planet.Magnetic.ionosBounds_m.
            sigmaIonos_Sm: Ionosphere Pedersen conductivity in S/m in Planet.Magnetic.sigmaIonosPedersen_Sm.
    """

    if name == 'R_m':
        Planet.Bulk.R_m = val
    elif name == 'xFeS':
        Planet.Core.xFeS = val
        Planet.Do.CONSTANT_INNER_DENSITY = True
    elif name == 'rhoSilInput_kgm3':
        Planet.Sil.rhoSilWithCore_kgm3 = val
        Planet.Do.CONSTANT_INNER_DENSITY = True
    elif name == 'wOcean_ppt':
        Planet.Ocean.wOcean_ppt = val
    elif name == 'Tb_K':
        Planet.Bulk.Tb_K = val
    elif name == 'ionosTop_km' or name == 'sigmaIonos_Sm':
        # Make sure ionosphere top altitude and conductivity are both set and valid
        if Planet.Magnetic.ionosBounds_m is None or np.any(np.isnan(Planet.Magnetic.ionosBounds_m)):
            Planet.Magnetic.ionosBounds_m = [Constants.ionosTopDefault_km*1e3]
        elif not isinstance(Planet.Magnetic.ionosBounds_m, Iterable):
            Planet.Magnetic.ionosBounds_m = [Planet.Magnetic.ionosBounds_m]

        if Planet.Magnetic.sigmaIonosPedersen_Sm is None or np.any(np.isnan(Planet.Magnetic.sigmaIonosPedersen_Sm)):
            Planet.Magnetic.sigmaIonosPedersen_Sm = [Constants.sigmaIonosPedersenDefault_Sm]
        elif not isinstance(Planet.Magnetic.sigmaIonosPedersen_Sm, Iterable):
            Planet.Magnetic.sigmaIonosPedersen_Sm = [Planet.Magnetic.sigmaIonosPedersen_Sm]

        if name == 'ionosTop_km':
            Planet.Magnetic.ionosBounds_m[-1] = val*1e3
        else:
            Planet.Magnetic.sigmaIonosPedersen_Sm[-1] = val
    elif name == 'silPhi_frac':
        Planet.Sil.phiRockMax_frac = val
        if val == 0:
            Planet.Do.POROUS_ROCK = False
        else:
            Planet.Do.POROUS_ROCK = True
        if Planet.Sil.porosType is None or Planet.Sil.porosType == 'none':
            Planet.Sil.porosType = 'Han2014'
        Planet.Do.CONSTANT_INNER_DENSITY = False
    elif name == 'silPclosure_MPa':
        Planet.Sil.Pclosure_MPa = val
        Planet.Do.POROUS_ROCK = True
        Planet.Do.CONSTANT_INNER_DENSITY = False
    elif name == 'icePhi_frac':
        Planet.Ocean.phiMax_frac = {key: val for key in Planet.Ocean.phiMax_frac.keys()}
        if val == 0:
            Planet.Do.POROUS_ICE = False
        else:
            Planet.Do.POROUS_ICE = True
    elif name == 'icePclosure_MPa':
        Planet.Ocean.Pclosure_MPa = {key: val for key in Planet.Ocean.Pclosure_MPa.keys()}
        Planet.Do.POROUS_ICE = True
    elif name == 'Htidal_Wm3':
        Planet.Sil.Htidal_Wm3 = val
    elif name == 'Qrad_Wkg':
        Planet.Sil.Qrad_Wkg = val
    elif name == 'qSurf_Wm2':
        Planet.Bulk.qSurf_Wm2 = val
    elif name == 'compOcean':
        Planet.Ocean.comp = val
    elif name == 'compSil':
        Planet.Sil.mantleEOS = val
        Planet.Do.CONSTANT_INNER_DENSITY = False
    elif name == 'compFe':
        Planet.Core.coreEOS = val
        Planet.Do.CONSTANT_INNER_DENSITY = False
    elif name == 'wFeCore_ppt':
        Planet.Core.wFe_ppt = val
        Planet.Core.coreEOS = 'Fe-S_3D_EOS.mat'
        Planet.Do.CONSTANT_INNER_DENSITY = False
    else:
        log.warning(f'No defined behavior for Planet setting named "{name}". Returning unchanged.')

    return Planet


def WriteExploreOgram(Exploration, Params):
    """ Organize Exploration results from an explore-o-gram run into a dict
        and print to a .mat file.
    """

    saveDict = {
        'bodyname': Exploration.bodyname,
        'xName': Exploration.xName,
        'yName': Exploration.yName,
        'zName': Exploration.zName,
        'wOcean_ppt': Exploration.wOcean_ppt,
        'oceanComp': Exploration.oceanComp,
        'R_m': Exploration.R_m,
        'Tb_K': Exploration.Tb_K,
        'xFeS': Exploration.xFeS,
        'rhoSilInput_kgm3': Exploration.rhoSilInput_kgm3,
        'silPhi_frac': Exploration.silPhi_frac,
        'icePhi_frac': Exploration.icePhi_frac,
        'silPclosure_MPa': Exploration.silPclosure_MPa,
        'icePclosure_MPa': Exploration.icePclosure_MPa,
        'ionosTop_km': Exploration.ionosTop_km,
        'sigmaIonos_Sm': Exploration.sigmaIonos_Sm,
        'Htidal_Wm3': Exploration.Htidal_Wm3,
        'Qrad_Wkg': Exploration.Qrad_Wkg,
        'rhoSilMean_kgm3': Exploration.rhoSilMean_kgm3,
        'rhoCoreMean_kgm3': Exploration.rhoCoreMean_kgm3,
        'sigmaMean_Sm': Exploration.sigmaMean_Sm,
        'sigmaTop_Sm': Exploration.sigmaTop_Sm,
        'Tmean_K': Exploration.Tmean_K,
        'D_km': Exploration.D_km,
        'zb_km': Exploration.zb_km,
        'Rcore_km': Exploration.Rcore_km
    }
    savemat(Params.DataFiles.exploreOgramFile, saveDict)
    log.info(f'Saved explore-o-gram {Params.DataFiles.exploreOgramFile} to disk.')

    return


def ReloadExploreOgram(bodyname, Params, fNameOverride=None):
    """ Reload a previously run explore-o-gram from disk.
    """

    if fNameOverride is None:
        if bodyname[:4] == 'Test':
            loadname = bodyname + ''
            bodydir = _TestImport
        else:
            loadname = bodyname
            bodydir = bodyname
        Planet = importlib.import_module(f'{bodydir}.PP{loadname}Explore').Planet
        Params.DataFiles, Params.FigureFiles = SetupFilenames(Planet, Params,
                                                              exploreAppend=f'{Params.Explore.xName}{Params.Explore.yName}',
                                                              figExploreAppend=f'{Params.Explore.zName}')
        reload = loadmat(Params.DataFiles.exploreOgramFile)
    else:
        reload = loadmat(fNameOverride)

    Exploration = ExplorationResults
    Exploration.bodyname = reload['bodyname'][0]
    Exploration.xName = reload['xName'][0]
    Exploration.yName = reload['yName'][0]
    Exploration.zName = reload['zName'][0]
    Exploration.wOcean_ppt = reload['wOcean_ppt']
    Exploration.oceanComp = reload['oceanComp']
    Exploration.R_m = reload['R_m']
    Exploration.Tb_K = reload['Tb_K']
    Exploration.xFeS = reload['xFeS']
    Exploration.rhoSilInput_kgm3 = reload['rhoSilInput_kgm3']
    Exploration.silPhi_frac = reload['silPhi_frac']
    Exploration.icePhi_frac = reload['icePhi_frac']
    Exploration.silPclosure_MPa = reload['silPclosure_MPa']
    Exploration.icePclosure_MPa = reload['icePclosure_MPa']
    Exploration.ionosTop_km = reload['ionosTop_km']
    Exploration.sigmaIonos_Sm = reload['sigmaIonos_Sm']
    Exploration.Htidal_Wm3 = reload['Htidal_Wm3']
    Exploration.Qrad_Wkg = reload['Qrad_Wkg']
    Exploration.rhoSilMean_kgm3 = reload['rhoSilMean_kgm3']
    Exploration.rhoCoreMean_kgm3 = reload['rhoCoreMean_kgm3']
    Exploration.sigmaMean_Sm = reload['sigmaMean_Sm']
    Exploration.sigmaTop_Sm = reload['sigmaTop_Sm']
    Exploration.Tmean_K = reload['Tmean_K']
    Exploration.D_km = reload['D_km']
    Exploration.zb_km = reload['zb_km']
    Exploration.Rcore_km = reload['Rcore_km']

    return Exploration, Params


def LoadPPfiles(Params, fNames, bodyname=''):
    """ Loads the settings in bodyname/fName.py to run or reload a specific model
        or models.
    """
    if Params.SPEC_FILE and fNames is not None:
        if bodyname == '':
            loadNames = [fName.replace('.py', '').replace(os.sep, '.') for fName in fNames]
        else:
            loadNames = ['.'.join([bodyname, fName.replace('.py', '').replace(os.sep, '.')]) for fName in fNames]
    else:
        if bodyname[:4] == 'Test':
            # Just get a length-1 list of the matching test profile
            loadNames = [f'{_TestImport}.PP{bodyname}']
        else:
            # Get model names from body directory
            models = FilesMatchingPattern(os.path.join(f'{bodyname}', f'PP{bodyname}*.py'))
            # If we don't find any, look in the defaults
            if np.size(models) == 0:
                models = FilesMatchingPattern(os.path.join(_Defaults, bodyname, f'PP{bodyname}*.py'))
                # Copy them over where we will find them if we have them
                if np.size(models) != 0:
                    [CopyCarefully(model, os.path.join(bodyname, os.path.basename(model))) for model in models]
                else:
                    raise ValueError('No profiles were found for input filenames with the search string ' +
                                     os.path.join(f'{bodyname}', f'PP{bodyname}*.py') + '.')

            # Splitting on PP and taking the -1 chops the path out, then taking
            # all but the last 3 characters chops the .py extension
            models = [model.split('PP')[-1][:-3] for model in models]
            # If there is a standard profile, make sure it is run first
            if bodyname in models:
                models.remove(bodyname)
                models.insert(0, bodyname)
            loadNames = ['.'.join([bodyname, f'PP{model}']) for model in models]

    if (Params.RUN_ALL_PROFILES and Params.COMPARE) or Params.SPEC_FILE:
        Params.nModels = np.size(loadNames)
    else:
        Params.nModels = 1

    return Params, loadNames


def RunPPfile(bodyname, fName, Params=None):
    # Simple wrapper to just run the given profile through the main PP function
    if fName[-3:] == '.py':
        loadName = fName[:-3]
    else:
        loadName = fName
    if Params is None:
        Params = configParams
    fName = f'{loadName}.py'
    expected = os.path.join(bodyname, fName)
    if not os.path.isfile(expected):
        default = os.path.join(_Defaults, bodyname, fName)
        if os.path.isfile(default):
            CopyCarefully(default, expected)
        else:
            log.warning(f'{expected} does not exist and no default was found at {default}.')
    Planet = importlib.import_module(f'{bodyname}.{loadName}').Planet
    Planet, Params = PlanetProfile(Planet, Params)
    
    return Planet, Params


if __name__ == '__main__':
    # Command line args
    nArgs = len(sys.argv)
    clArg = None
    fNames = None
    if nArgs > 1 and ('PP' not in sys.argv[1] and '.txt' not in sys.argv[1]):
        # Body name was passed as command line argument
        bodyname = sys.argv[1]

        # Additional command line arguments
        if nArgs > 2:
            if 'PP' in sys.argv[2]:
                log.debug('PP in CL arg 2 -- interpreting as (list of) filename(s).')
                fNames = sys.argv[2:]
            elif '.txt' in sys.argv[2]:
                log.debug('.txt in CL arg 2 -- interpreting as (list of) filename(s) to reload.')
                fNames = sys.argv[2:]
                clArg = 'reload'
            else:
                clArg = sys.argv[2]
                if nArgs > 3:
                    if 'PP' in sys.argv[3]:
                        fNames = sys.argv[3:]
                    else:
                        log.warning(f'Too many command line args passed. Ignoring command "{sys.argv[3]}" and any after it.')
    elif 'PP' in sys.argv[1]:
        log.debug('PP in first CL arg -- interpreting as (list of) filename(s).')
        fNames = sys.argv[1:]
        if np.size(fNames) == 1:
            bodyname = os.path.split(sys.argv[1])[0]
        else:
            bodyname = ''
    elif '.txt' in sys.argv[1]:
        log.debug('.txt in first CL arg -- interpreting as (list of) filename(s) to reload.')
        fNames = sys.argv[1:]
        bodyname = ''
        clArg = 'reload'
    else:
        # No command line argument, ask user which body to run
        bodyname = input('Please input body name: ')


    run(bodyname=bodyname, opt=clArg, fNames=fNames)
