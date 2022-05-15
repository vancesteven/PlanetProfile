# Import necessary Python modules
import os, sys, time, importlib
import numpy as np
import logging as log
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
from PlanetProfile.Plotting.ProfilePlots import GeneratePlots, PlotInductOgram, PlotInductOgramPhaseSpace
from PlanetProfile.Thermodynamics.LayerPropagators import IceLayers, OceanLayers, InnerLayers
from PlanetProfile.Thermodynamics.Electrical import ElecConduct
from PlanetProfile.Thermodynamics.Seismic import SeismicCalcs
from PlanetProfile.Utilities.defineStructs import Constants, FigureFilesSubstruct, PlanetStruct
from PlanetProfile.Utilities.SetupInit import SetupInit, SetupFilenames
from PlanetProfile.Utilities.SummaryTables import GetLayerMeans, PrintGeneralSummary, PrintLayerSummaryLatex, PrintLayerTableLatex

# Parallel processing
import multiprocessing as mtp
mtpFork = mtp.get_context('fork')

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
        Induction, Params = InductOgram(bodyname, Params)
        if not Params.SKIP_PLOTS:
            PlotInductOgram(Induction, Params)

            if Params.COMPARE:
                inductOgramFiles = FilesMatchingPattern(os.path.join(Params.DataFiles.fNameInductOgram+'*.mat'))
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

    if not Params.SKIP_PLOTS and not Params.DO_INDUCTOGRAM:
        # Calculate large-scale layer properties
        PlanetList, Params = GetLayerMeans(np.array([Planet]), Params)
        # Plotting functions
        GeneratePlots(PlanetList, Params)
        Planet = PlanetList[0]

    # Magnetic induction calculations
    if Params.CALC_CONDUCT and Params.CALC_NEW_INDUCT:
        # Calculate induced magnetic moments
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
                expected = os.path.join(bodyname, fName)
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
                expected = os.path.join(bodyname, fName)
                if not os.path.isfile(expected):
                    default = os.path.join(_Defaults, bodyname, fName)
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
        Planet.Sil.sigmaPorousLayerMean_Sm, Planet.etaConv_Pas, Planet.RaConvect, Planet.RaConvectIII, Planet.RaConvectV, \
        Planet.RaCrit, Planet.RaCritIII, Planet.RaCritV, Planet.eLid_m, Planet.eLidIII_m, Planet.eLidV_m, \
        Planet.Dconv_m, Planet.DconvIII_m, Planet.DconvV_m, Planet.deltaTBL_m, Planet.deltaTBLIII_m, Planet.deltaTBLV_m \
            = (float(f.readline().split('=')[-1]) for _ in range(62))
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
        log.basicConfig(level=log.WARN+5, format=Params.printFmt, force=True)

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
            changes (dict): Dict of {input: value} pairs to change in Planet. Options are:
                compOcean: Ocean composition string in Planet.Ocean.comp
                compSil: Silicate composition to use from available Perplex output files
                compFe: Iron core composition to use from available Perplex output files
                wOcean_ppt: Salinity in Planet.Ocean.wOcean_ppt
                Tb_K: Ocean bottom temperature in K in Planet.Bulk.Tb_K
                silDensity_kgm3: Fixed density in silicate layers in Planet.Sil.rhoSilWithCore_kgm3 (for use with Planet.Do.CONSTANT_INNER_DENSITY)
                silPhi_frac: Vacuum-extrapolated porosity in silicates in Planet.Sil.phiRockMax_frac
                silPclosure_MPa: Pore closure pressure in silicates in Planet.Sil.Pclosure_MPa
                icePhi_frac: Vacuum porosity in surface ice Ih in Planet.Ocean.phiRockMax_frac
                icePclosure_MPa: Pore closure pressure in ice in Planet.Sil.Pclosure_MPa
                Htidal_Wm3: Fixed tidal heating in silicates in Planet.Sil.Htidal_Wm3
                Qrad_Wkg: Fixed radiogenic heating in silicates in Planet.Sil.Qrad_Wkg
                ionosBounds_m: Ionosphere boundary altitudes above the surface in m in Planet.Magnetic.ionosBounds_m. Correspond to outer radii of each sigma below of same index.
                sigmaIonos_Sm: Ionosphere Pedersen conductivities in S/m in Planet.Magnetic.sigmaIonosPedersen_Sm. Must have same length as ionosBounds_m
    """
    if changes is None:
        log.warning('No changes passed to UpdateRun.')
        changes = {}

    for key, value in changes.items():
        if key == 'wOcean_ppt':
            Planet.Ocean.wOcean_ppt = value
        elif key == 'Tb_K':
            Planet.Bulk.Tb_K = value
        elif key == 'silDensity_kgm3':
            Planet.Sil.rhoSilWithCore_kgm3 = value
            Planet.Do.CONSTANT_INNER_DENSITY = True
        elif key == 'ionosBounds_m':
            Planet.Magnetic.ionosBounds_m = value
        elif key == 'sigmaIonos_Sm':
            Planet.Magnetic.sigmaIonosPedersen_Sm = value
        elif key == 'silPhi_frac':
            Planet.Sil.phiRockMax_frac = value
            if value == 0:
                Planet.Do.POROUS_ROCK = False
            else:
                Planet.Do.POROUS_ROCK = True
            Planet.Do.CONSTANT_INNER_DENSITY = False
        elif key == 'silPclosure_MPa':
            Planet.Sil.Pclosure_MPa = value
            Planet.Do.POROUS_ROCK = True
            Planet.Do.CONSTANT_INNER_DENSITY = False
        elif key == 'icePhi_frac':
            Planet.Ocean.phiMax_frac['Ih'] = value
            if value == 0:
                Planet.Do.POROUS_ICE = False
            else:
                Planet.Do.POROUS_ICE = True
        elif key == 'icePclosure_MPa':
            Planet.Ocean.Pclosure_MPa['Ih'] = value
            Planet.Do.POROUS_ICE = True
        elif key == 'Htidal_Wm3':
            Planet.Sil.Htidal_Wm3 = value
        elif key == 'Qrad_Wkg':
            Planet.Sil.Qrad_Wkg = value
        elif key == 'compOcean':
            Planet.Ocean.comp = value
        elif key == 'compSil':
            Planet.Sil.mantleEOS = value
            Planet.Do.CONSTANT_INNER_DENSITY = False
        elif key == 'compFe':
            Planet.Core.coreEOS = value
            Planet.Do.CONSTANT_INNER_DENSITY = False
        else:
            log.warning(f'UpdateRun key not recognized: {key}. Skipping.')

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
            bodydir = _TestImport
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
        log.info(f'Parallel run elapsed time: {dt} s.')

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
        Induction.w_ppt = np.array([[Planeti.Ocean.wOcean_ppt for Planeti in line] for line in PlanetGrid])
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
        'w_ppt': Induction.w_ppt,
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
            bodyname = 'Test'
            bodydir = _TestImport
        else:
            loadname = bodyname
            bodydir = bodyname
        Planet = importlib.import_module(f'{bodydir}.PP{loadname}InductOgram').Planet
        Params.DataFiles, Params.FigureFiles = SetupFilenames(Planet, Params)
        reload = loadmat(Params.DataFiles.inductOgramFile)
    else:
        reload = loadmat(fNameOverride)

    Induction = InductionResults
    Induction.bodyname = reload['bodyname'][0]
    Induction.yName = reload['yName'][0]
    Induction.Texc_hr = {key.strip(): value for key, value in zip(reload['Texc_hr_keys'], reload['Texc_hr_values'][0])}
    Induction.Amp = reload['Amp']
    Induction.phase = reload['phase']
    Induction.Bix_nT = reload['Bix_nT']
    Induction.Biy_nT = reload['Biy_nT']
    Induction.Biz_nT = reload['Biz_nT']
    Induction.w_ppt = reload['w_ppt']
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
    if Params.logParallel > log.INFO:
        log.info('Quieting messages to avoid spam in gridded run.')
    saveLevel = log.getLogger().level + 0
    log.getLogger().setLevel(Params.logParallel)
    dims = np.shape(PlanetList)
    nParDims = np.size(dims)
    if nParDims == 1 and np.size(PlanetList) == 1:
        PlanetList[0] = PlanetProfile(PlanetList[0], Params)
    else:
        PlanetList1D = np.reshape(PlanetList, -1)
        if Params.DO_PARALLEL:
            # Prevent slowdowns from competing process spawning when #cores > #jobs
            nCores = np.minimum(Params.maxCores, np.product(dims))
            if Params.INDUCTOGRAM_IN_PROGRESS and Params.Induct.inductOtype == 'sigma':
                pool = mtpFork.Pool(nCores)
                parResult = [pool.apply_async(MagneticInduction, (deepcopy(Planet),
                                                                  deepcopy(Params))) for Planet in PlanetList1D]
                pool.close()
                pool.join()
            else:
                pool = mtpFork.Pool(nCores)
                parResult = [pool.apply_async(PlanetProfile, (deepcopy(Planet),
                                                              deepcopy(Params))) for Planet in PlanetList1D]
                pool.close()
                pool.join()

            for i, result in enumerate(parResult):
                PlanetList1D[i] = result.get()[0]
        else:
            log.profile('Running grid without parallel processing. This may take some time.')
            if Params.Induct.inductOtype == 'sigma':
                PlanetList1D = np.array([MagneticInduction(deepcopy(Planet), deepcopy(Params)) for Planet in PlanetList1D])[:, 0]
            else:
                PlanetList1D = np.array([PlanetProfile(deepcopy(Planet), deepcopy(Params)) for Planet in PlanetList1D])[:, 0]

        PlanetList = np.reshape(PlanetList1D, dims)

        # Return settings to what they were before we entered here
        log.basicConfig(level=saveLevel, format=Params.printFmt)
        Params.INDUCTOGRAM_IN_PROGRESS = False

    return PlanetList


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


def CompareProfile(Planet, Params, fname2, tol=0.01, tiny=1e-6):
    """ Checks saved data in a named file against the current Planet run
        to say whether the contents of the files differ.

        Args:
             Planet (PlanetStruct): The current run's completed data
             Params (ParamsStruct): Params for the current run
             fname2 (string): The filename for save data to compare against
             tol (float): Tolerance in % difference to consider values the same
                in comparing layer arrays
        Returns:
            None
    """

    log.info(f'Comparing current run with {fname2}...')
    Planet2 = deepcopy(Planet)
    Params2 = deepcopy(Params)
    Planet2, Params2 = ReloadProfile(Planet2, Params2, fnameOverride=fname2)

    # Avoid divide-by-zero errors
    if Planet.Ocean.wOcean_ppt == 0: Planet.Ocean.wOcean_ppt = tiny
    if Planet.Bulk.Tb_K == 0: Planet.Bulk.Tb_K = tiny
    if Planet.zb_km == 0: Planet.zb_km = tiny
    if Planet.zClath_m == 0: Planet.zClath_m = tiny
    if Planet.Pb_MPa == 0: Planet.Pb_MPa = tiny
    if Planet.PbI_MPa == 0: Planet.PbI_MPa = tiny
    if Planet.Ocean.deltaP == 0: Planet.Ocean.deltaP = tiny
    if Planet.CMR2mean == 0: Planet.CMR2mean = tiny
    if Planet.Ocean.QfromMantle_Wm2 == 0: Planet.Ocean.QfromMantle_Wm2 = tiny
    if Planet.Sil.phiRockMax_frac == 0: Planet.Sil.phiRockMax_frac = tiny
    if Planet.Sil.Rmean_m == 0: Planet.Sil.Rmean_m = tiny
    if Planet.Sil.Rrange_m == 0: Planet.Sil.Rrange_m = tiny
    if Planet.Sil.rhoMean_kgm3 == 0: Planet.Sil.rhoMean_kgm3 = tiny
    if Planet.Core.Rmean_m == 0: Planet.Core.Rmean_m = tiny
    if Planet.Core.Rrange_m == 0: Planet.Core.Rrange_m = tiny
    if Planet.Core.rhoMean_kgm3 == 0: Planet.Core.rhoMean_kgm3 = tiny

    # Avoid false positives for zero values
    if Planet2.Ocean.wOcean_ppt == 0: Planet2.Ocean.wOcean_ppt = tiny
    if Planet2.Bulk.Tb_K == 0: Planet2.Bulk.Tb_K = tiny
    if Planet2.zb_km == 0: Planet2.zb_km = tiny
    if Planet2.zClath_m == 0: Planet2.zClath_m = tiny
    if Planet2.Pb_MPa == 0: Planet2.Pb_MPa = tiny
    if Planet2.PbI_MPa == 0: Planet2.PbI_MPa = tiny
    if Planet2.Ocean.deltaP == 0: Planet2.Ocean.deltaP = tiny
    if Planet2.CMR2mean == 0: Planet2.CMR2mean = tiny
    if Planet2.Ocean.QfromMantle_Wm2 == 0: Planet2.Ocean.QfromMantle_Wm2 = tiny
    if Planet2.Sil.phiRockMax_frac == 0: Planet2.Sil.phiRockMax_frac = tiny
    if Planet2.Sil.Rmean_m == 0: Planet2.Sil.Rmean_m = tiny
    if Planet2.Sil.Rrange_m == 0: Planet2.Sil.Rrange_m = tiny
    if Planet2.Sil.rhoMean_kgm3 == 0: Planet2.Sil.rhoMean_kgm3 = tiny
    if Planet2.Core.Rmean_m == 0: Planet2.Core.Rmean_m = tiny
    if Planet2.Core.Rrange_m == 0: Planet2.Core.Rrange_m = tiny
    if Planet2.Core.rhoMean_kgm3 == 0: Planet2.Core.rhoMean_kgm3 = tiny

    # Compare header info
    same_nHeadLines = Params.nHeadLines == Params2.nHeadLines
    same_comp = Planet.Ocean.comp == Planet2.Ocean.comp
    same_Fe_CORE = Planet.Do.Fe_CORE and Planet2.Do.Fe_CORE
    same_wOcean_ppt = (Planet.Ocean.wOcean_ppt - Planet2.Ocean.wOcean_ppt) / Planet.Ocean.wOcean_ppt < tol
    same_Tb_K = (Planet.Bulk.Tb_K - Planet2.Bulk.Tb_K) / Planet.Bulk.Tb_K < tol
    same_zb_km = (Planet.zb_km > Planet2.zb_km) / Planet.zb_km < tol
    same_zClath_m = (Planet.zClath_m - Planet2.zClath_m) / Planet.zClath_m < tol
    same_Pb_MPa = (Planet.Pb_MPa - Planet2.Pb_MPa) / Planet.Pb_MPa < tol
    same_PbI_MPa = (Planet.PbI_MPa - Planet2.PbI_MPa) / Planet.PbI_MPa < tol
    same_deltaP = (Planet.Ocean.deltaP - Planet2.Ocean.deltaP) / Planet.Ocean.deltaP < tol
    same_CMR2mean = (Planet.CMR2mean - Planet2.CMR2mean) / Planet.CMR2mean < tol
    same_QfromMantle_Wm2 = (Planet.Ocean.QfromMantle_Wm2 - Planet2.Ocean.QfromMantle_Wm2) / Planet.Ocean.QfromMantle_Wm2 < tol
    same_phiRockMax_frac = (Planet.Sil.phiRockMax_frac - Planet2.Sil.phiRockMax_frac) / Planet.Sil.phiRockMax_frac < tol
    same_silRmean_m = (Planet.Sil.Rmean_m - Planet2.Sil.Rmean_m) / Planet.Sil.Rmean_m < tol
    same_silRrange_m = (Planet.Sil.Rrange_m - Planet2.Sil.Rrange_m) / Planet.Sil.Rrange_m < tol
    same_silrhoMean_kgm3 = (Planet.Sil.rhoMean_kgm3 - Planet2.Sil.rhoMean_kgm3) / Planet.Sil.rhoMean_kgm3 < tol
    same_coreRmean_m = (Planet.Core.Rmean_m - Planet2.Core.Rmean_m) / Planet.Core.Rmean_m < tol
    same_coreRrange_m = (Planet.Core.Rrange_m - Planet2.Core.Rrange_m) / Planet.Core.Rrange_m < tol
    same_corerhoMean_kgm3 = (Planet.Core.rhoMean_kgm3 - Planet2.Core.rhoMean_kgm3) / Planet.Core.rhoMean_kgm3 < tol
    same_nClath = Planet.Steps.nClath == Planet2.Steps.nClath
    same_nIceI = Planet.Steps.nIceI == Planet2.Steps.nIceI
    same_nIceIIILitho = Planet.Steps.nIceIIILitho == Planet2.Steps.nIceIIILitho
    same_nIceVLitho = Planet.Steps.nIceVLitho == Planet2.Steps.nIceVLitho
    same_nHydro = Planet.Steps.nHydro == Planet2.Steps.nHydro
    same_nSil = Planet.Steps.nSil == Planet2.Steps.nSil
    same_nCore = Planet.Steps.nCore == Planet2.Steps.nCore

    same_steps = same_nClath and same_nIceI and same_nIceIIILitho and same_nIceVLitho and same_nHydro and same_nSil \
        and same_nCore
    headers_match = same_nHeadLines and same_comp and same_Fe_CORE and same_wOcean_ppt and same_Tb_K and same_zb_km \
        and same_zClath_m and same_Pb_MPa and same_PbI_MPa and same_deltaP and same_CMR2mean and same_QfromMantle_Wm2 \
        and same_phiRockMax_frac and same_silRmean_m and same_silRrange_m and same_silrhoMean_kgm3 and same_coreRmean_m \
        and same_coreRrange_m and same_corerhoMean_kgm3 and same_steps

    if not headers_match:
        if not same_nHeadLines: log.info(f'nHeadLines differs. {Params.nHeadLines:d} | {Params2.nHeadLines:d}')
        if not same_comp: log.info(f'Ocean.comp differs. {Planet.Ocean.comp} | {Planet2.Ocean.comp}')
        if not same_Fe_CORE: log.info(f'Do.Fe_CORE differs. {Planet.Do.Fe_CORE} | {Planet2.Do.Fe_CORE}')
        if not same_wOcean_ppt: log.info(f'Ocean.wOcean_ppt differs. {Planet.Ocean.wOcean_ppt:.3f} | {Planet2.Ocean.wOcean_ppt:.3f}')
        if not same_Tb_K: log.info(f'Bulk.Tb_K differs. {Planet.Bulk.Tb_K:.4f} | {Planet2.Bulk.Tb_K:.4f}')
        if not same_zb_km: log.info(f'zb_km differs. {Planet.zb_km:.3f} | {Planet2.zb_km:.3f}')
        if not same_zClath_m: log.info(f'zClath_m differs. {Planet.zClath_m:.3f} | {Planet2.zClath_m:.3f}')
        if not same_Pb_MPa: log.info(f'Pb_MPa differs. {Planet.Pb_MPa:.3f} | {Planet2.Pb_MPa:.3f}')
        if not same_PbI_MPa: log.info(f'PbI_MPa differs. {Planet.PbI_MPa:.3f} | {Planet2.PbI_MPa:.3f}')
        if not same_deltaP: log.info(f'Ocean.deltaP differs. {Planet.Ocean.deltaP:.3f} | {Planet2.Ocean.deltaP:.3f}')
        if not same_CMR2mean: log.info(f'CMR2mean differs. {Planet.CMR2mean:.3f} | {Planet2.CMR2mean:.3f}')
        if not same_QfromMantle_Wm2: log.info(f'Ocean.QfromMantle_Wm2 differs. {Planet.Ocean.QfromMantle_Wm2:.3f} | {Planet2.Ocean.QfromMantle_Wm2:.3f}')
        if not same_phiRockMax_frac: log.info(f'Sil.phiRockMax_frac differs. {Planet.Sil.phiRockMax_frac:.3f} | {Planet2.Sil.phiRockMax_frac:.3f}')
        if not same_silRmean_m: log.info(f'Sil.Rmean_m differs. {Planet.Sil.Rmean_m:.3f} | {Planet2.Sil.Rmean_m:.3f}')
        if not same_silRrange_m: log.info(f'Sil.Rrange_m differs. {Planet.Sil.Rrange_m:.3f} | {Planet2.Sil.Rrange_m:.3f}')
        if not same_silrhoMean_kgm3: log.info(f'Sil.rhoMean_kgm3 differs. {Planet.Sil.rhoMean_kgm3:.3f} | {Planet2.Sil.rhoMean_kgm3:.3f}')
        if not same_coreRmean_m: log.info(f'Core.Rmean_m differs. {Planet.Core.Rmean_m:.3f} | {Planet2.Core.Rmean_m:.3f}')
        if not same_coreRrange_m: log.info(f'Core.Rrange_m differs. {Planet.Core.Rrange_m:.3f} | {Planet2.Core.Rrange_m:.3f}')
        if not same_corerhoMean_kgm3: log.info(f'Core.rhoMean_kgm3 differs. {Planet.Core.rhoMean_kgm3:.3f} | {Planet2.Core.rhoMean_kgm3:.3f}')
        if not same_nClath: log.info(f'Steps.nClath differs. {Planet.Steps.nClath:d} | {Planet2.Steps.nClath:d}')
        if not same_nIceI: log.info(f'Steps.nIceI differs. {Planet.Steps.nIceI:d} | {Planet2.Steps.nIceI:d}')
        if not same_nIceIIILitho: log.info(f'Steps.nIceIIILitho differs. {Planet.Steps.nIceIIILitho:d} | {Planet2.Steps.nIceIIILitho:d}')
        if not same_nIceVLitho: log.info(f'Steps.nIceVLitho differs. {Planet.Steps.nIceVLitho:d} | {Planet2.Steps.nIceVLitho:d}')
        if not same_nHydro: log.info(f'Steps.nHydro differs. {Planet.Steps.nHydro:d} | {Planet2.Steps.nHydro:d}')
        if not same_nSil: log.info(f'Steps.nSil differs. {Planet.Steps.nSil:d} | {Planet2.Steps.nSil:d}')
        if not same_nCore: log.info(f'Steps.nCore differs. {Planet.Steps.nCore:d} | {Planet2.Steps.nCore:d}')
    else:
        log.info('All header values match!')

    if same_steps:
        log.info(f'Some arrays differ. The first index more than {tol*100:.2f}% different will be printed.')
        iIO = Planet.Steps.nSurfIce
        iOS = Planet.Steps.nHydro
        iSC = Planet.Steps.nHydro+Planet.Steps.nSil
        iCC = Planet.Steps.nTotal

        # Avoid divide-by-zero errors
        Planet.P_MPa[Planet.P_MPa==0] = tiny
        Planet.T_K[Planet.T_K==0] = tiny
        Planet.r_m[Planet.r_m==0] = tiny
        Planet.rho_kgm3[Planet.rho_kgm3==0] = tiny
        Planet.Cp_JkgK[Planet.Cp_JkgK==0] = tiny
        Planet.alpha_pK[Planet.alpha_pK==0] = tiny
        Planet.g_ms2[Planet.g_ms2==0] = tiny
        Planet.phi_frac[Planet.phi_frac==0] = tiny
        Planet.sigma_Sm[Planet.sigma_Sm==0] = tiny
        Planet.kTherm_WmK[Planet.kTherm_WmK==0] = tiny
        Planet.Seismic.VP_kms[Planet.Seismic.VP_kms==0] = tiny
        Planet.Seismic.VS_kms[Planet.Seismic.VS_kms==0] = tiny
        Planet.Seismic.QS[Planet.Seismic.QS==0] = tiny
        Planet.Seismic.KS_GPa[Planet.Seismic.KS_GPa==0] = tiny
        Planet.Seismic.GS_GPa[Planet.Seismic.GS_GPa==0] = tiny
        Planet.Ppore_MPa[Planet.Ppore_MPa==0] = tiny
        Planet.rhoMatrix_kgm3[Planet.rhoMatrix_kgm3==0] = tiny
        Planet.rhoPore_kgm3[Planet.rhoPore_kgm3==0] = tiny

        # Avoid false positives from zero values
        Planet2.P_MPa[Planet2.P_MPa==0] = tiny
        Planet2.T_K[Planet2.T_K==0] = tiny
        Planet2.r_m[Planet2.r_m==0] = tiny
        Planet2.rho_kgm3[Planet2.rho_kgm3==0] = tiny
        Planet2.Cp_JkgK[Planet2.Cp_JkgK==0] = tiny
        Planet2.alpha_pK[Planet2.alpha_pK==0] = tiny
        Planet2.g_ms2[Planet2.g_ms2==0] = tiny
        Planet2.phi_frac[Planet2.phi_frac==0] = tiny
        Planet2.sigma_Sm[Planet2.sigma_Sm==0] = tiny
        Planet2.kTherm_WmK[Planet2.kTherm_WmK==0] = tiny
        Planet2.Seismic.VP_kms[Planet2.Seismic.VP_kms==0] = tiny
        Planet2.Seismic.VS_kms[Planet2.Seismic.VS_kms==0] = tiny
        Planet2.Seismic.QS[Planet2.Seismic.QS==0] = tiny
        Planet2.Seismic.KS_GPa[Planet2.Seismic.KS_GPa==0] = tiny
        Planet2.Seismic.GS_GPa[Planet2.Seismic.GS_GPa==0] = tiny
        Planet2.Ppore_MPa[Planet2.Ppore_MPa==0] = tiny
        Planet2.rhoMatrix_kgm3[Planet2.rhoMatrix_kgm3==0] = tiny
        Planet2.rhoPore_kgm3[Planet2.rhoPore_kgm3==0] = tiny

        # Compare layer calculations
        diff_P_MPa = [i[0] for i, val in np.ndenumerate(Planet.P_MPa) if abs(val-Planet2.P_MPa[i[0]])/val>tol]
        diff_T_K = [i[0] for i, val in np.ndenumerate(Planet.T_K) if abs(val-Planet2.T_K[i[0]])/val>tol]
        diff_r_m = [i[0] for i, val in np.ndenumerate(Planet.r_m) if abs(val-Planet2.r_m[i[0]])/val>tol]
        diff_phase = [i[0] for i, val in np.ndenumerate(Planet.phase) if Planet2.phase[i[0]]!=val]
        diff_rho_kgm3 = [i[0] for i, val in np.ndenumerate(Planet.rho_kgm3) if abs(val-Planet2.rho_kgm3[i[0]])/val>tol]
        diff_Cp_JkgK = [i[0] for i, val in np.ndenumerate(Planet.Cp_JkgK) if abs(val-Planet2.Cp_JkgK[i[0]])/val>tol]
        diff_alpha_pK = [i[0] for i, val in np.ndenumerate(Planet.alpha_pK) if abs(val-Planet2.alpha_pK[i[0]])/val>tol]
        diff_g_ms2 = [i[0] for i, val in np.ndenumerate(Planet.g_ms2) if abs(val-Planet2.g_ms2[i[0]])/val>tol]
        diff_phi_frac = [i[0] for i, val in np.ndenumerate(Planet.phi_frac) if abs(val-Planet2.phi_frac[i[0]])/val>tol]
        diff_sigma_Sm = [i[0] for i, val in np.ndenumerate(Planet.sigma_Sm) if abs(val-Planet2.sigma_Sm[i[0]])/val>tol]
        diff_kTherm_WmK = [i[0] for i, val in np.ndenumerate(Planet.kTherm_WmK) if abs(val-Planet2.kTherm_WmK[i[0]])/val>tol]
        diff_VP_kms = [i[0] for i, val in np.ndenumerate(Planet.Seismic.VP_kms) if abs(val-Planet2.Seismic.VP_kms[i[0]])/val>tol]
        diff_VS_kms = [i[0] for i, val in np.ndenumerate(Planet.Seismic.VS_kms) if abs(val-Planet2.Seismic.VS_kms[i[0]])/val>tol]
        diff_QS = [i[0] for i, val in np.ndenumerate(Planet.Seismic.QS) if abs(val-Planet2.Seismic.QS[i[0]])/val>tol]
        diff_KS_GPa = [i[0] for i, val in np.ndenumerate(Planet.Seismic.KS_GPa) if abs(val-Planet2.Seismic.KS_GPa[i[0]])/val>tol]
        diff_GS_GPa = [i[0] for i, val in np.ndenumerate(Planet.Seismic.GS_GPa) if abs(val-Planet2.Seismic.GS_GPa[i[0]])/val>tol]
        diff_Ppore_MPa = [i[0] for i, val in np.ndenumerate(Planet.Ppore_MPa) if abs(val-Planet2.Ppore_MPa[i[0]])/val>tol]
        diff_rhoMatrix_kgm3 = [i[0] for i, val in np.ndenumerate(Planet.rhoMatrix_kgm3) if abs(val-Planet2.rhoMatrix_kgm3[i[0]])/val>tol]
        diff_rhoPore_kgm3 = [i[0] for i, val in np.ndenumerate(Planet.rhoPore_kgm3) if abs(val-Planet2.rhoPore_kgm3[i[0]])/val>tol]

        same_P_MPa = len(diff_P_MPa) == 0
        same_T_K = len(diff_T_K) == 0
        same_r_m = len(diff_r_m) == 0
        same_phase = len(diff_phase) == 0
        same_rho_kgm3 = len(diff_rho_kgm3) == 0
        same_Cp_JkgK = len(diff_Cp_JkgK) == 0
        same_alpha_pK = len(diff_alpha_pK) == 0
        same_g_ms2 = len(diff_g_ms2) == 0
        same_phi_frac = len(diff_phi_frac) == 0
        same_sigma_Sm = len(diff_sigma_Sm) == 0
        same_kTherm_WmK = len(diff_kTherm_WmK) == 0
        same_VP_kms = len(diff_VP_kms) == 0
        same_VS_kms = len(diff_VS_kms) == 0
        same_QS = len(diff_QS) == 0
        same_KS_GPa = len(diff_KS_GPa) == 0
        same_GS_GPa = len(diff_GS_GPa) == 0
        same_Ppore_MPa = len(diff_Ppore_MPa) == 0
        same_rhoMatrix_kgm3 = len(diff_rhoMatrix_kgm3) == 0
        same_rhoPore_kgm3 = len(diff_rhoPore_kgm3) == 0

        layers_match = same_P_MPa and same_T_K and same_r_m and same_phase and same_rho_kgm3 and same_Cp_JkgK \
            and same_alpha_pK and same_g_ms2 and same_phi_frac and same_sigma_Sm and same_kTherm_WmK and same_VP_kms \
            and same_VS_kms and same_QS and same_KS_GPa and same_GS_GPa and same_Ppore_MPa and same_rhoMatrix_kgm3 \
            and same_rhoPore_kgm3
        all_match = headers_match and layers_match

        if not layers_match:
            if not same_P_MPa: log.info(f'P_MPa differs in position {diff_P_MPa[0]:.3f}: {Planet.P_MPa[diff_P_MPa[0]]:.3f} | {Planet2.P_MPa[diff_P_MPa[0]]:.3f}')
            if not same_T_K: log.info(f'T_K differs in position {diff_T_K[0]:.3f}: {Planet.T_K[diff_T_K[0]]:.3f} | {Planet2.T_K[diff_T_K[0]]:.3f}')
            if not same_r_m: log.info(f'r_m differs in position {diff_r_m[0]:.3f}: {Planet.r_m[diff_r_m[0]]:.3f} | {Planet2.r_m[diff_r_m[0]]:.3f}')
            if not same_phase: log.info(f'phase differs in position {diff_phase[0]:d}: {Planet.phase[diff_phase[0]]:d} | {Planet2.phase[diff_phase[0]]:d}')
            if not same_rho_kgm3: log.info(f'rho_kgm3 differs in position {diff_rho_kgm3[0]:.3f}: {Planet.rho_kgm3[diff_rho_kgm3[0]]:.3f} | {Planet2.rho_kgm3[diff_rho_kgm3[0]]:.3f}')
            if not same_Cp_JkgK: log.info(f'Cp_JkgK differs in position {diff_Cp_JkgK[0]:.3f}: {Planet.Cp_JkgK[diff_Cp_JkgK[0]]:.3f} | {Planet2.Cp_JkgK[diff_Cp_JkgK[0]]:.3f}')
            if not same_alpha_pK: log.info(f'alpha_pK differs in position {diff_alpha_pK[0]:.3f}: {Planet.alpha_pK[diff_alpha_pK[0]]:.3f} | {Planet2.alpha_pK[diff_alpha_pK[0]]:.3f}')
            if not same_g_ms2: log.info(f'g_ms2 differs in position {diff_g_ms2[0]:.3f}: {Planet.g_ms2[diff_g_ms2[0]]:.3f} | {Planet2.g_ms2[diff_g_ms2[0]]:.3f}')
            if not same_phi_frac: log.info(f'phi_frac differs in position {diff_phi_frac[0]:.3f}: {Planet.phi_frac[diff_phi_frac[0]]:.3f} | {Planet2.phi_frac[diff_phi_frac[0]]:.3f}')
            if not same_sigma_Sm: log.info(f'sigma_Sm differs in position {diff_sigma_Sm[0]:.3f}: {Planet.sigma_Sm[diff_sigma_Sm[0]]:.3f} | {Planet2.sigma_Sm[diff_sigma_Sm[0]]:.3f}')
            if not same_kTherm_WmK: log.info(f'kTherm_WmK differs in position {diff_kTherm_WmK[0]:.3f}: {Planet.kTherm_WmK[diff_kTherm_WmK[0]]:.3f} | {Planet2.kTherm_WmK[diff_kTherm_WmK[0]]:.3f}')
            if not same_VP_kms: log.info(f'VP_kms differs in position {diff_VP_kms[0]:.3f}: {Planet.Seismic.VP_kms[diff_VP_kms[0]]:.3f} | {Planet2.Seismic.VP_kms[diff_VP_kms[0]]:.3f}')
            if not same_VS_kms: log.info(f'VS_kms differs in position {diff_VS_kms[0]:.3f}: {Planet.Seismic.VS_kms[diff_VS_kms[0]]:.3f} | {Planet2.Seismic.VS_kms[diff_VS_kms[0]]:.3f}')
            if not same_QS: log.info(f'QS differs in position {diff_QS[0]:.3f}: {Planet.Seismic.QS[diff_QS[0]]:.3f} | {Planet2.Seismic.QS[diff_QS[0]]:.3f}')
            if not same_KS_GPa: log.info(f'KS_GPa differs in position {diff_KS_GPa[0]:.3f}: {Planet.Seismic.KS_GPa[diff_KS_GPa[0]]:.3f} | {Planet2.Seismic.KS_GPa[diff_KS_GPa[0]]:.3f}')
            if not same_GS_GPa: log.info(f'GS_GPa differs in position {diff_GS_GPa[0]:.3f}: {Planet.Seismic.GS_GPa[diff_GS_GPa[0]]:.3f} | {Planet2.Seismic.GS_GPa[diff_GS_GPa[0]]:.3f}')
            if not same_Ppore_MPa: log.info(f'Ppore_MPa differs in position {diff_Ppore_MPa[0]:.3f}: {Planet.Ppore_MPa[diff_Ppore_MPa[0]]:.3f} | {Planet2.Ppore_MPa[diff_Ppore_MPa[0]]:.3f}')
            if not same_rhoMatrix_kgm3: log.info(f'rhoMatrix_kgm3 differs in position {diff_rhoMatrix_kgm3[0]:.3f}: {Planet.rhoMatrix_kgm3[diff_rhoMatrix_kgm3[0]]:.3f} | {Planet2.rhoMatrix_kgm3[diff_rhoMatrix_kgm3[0]]:.3f}')
            if not same_rhoPore_kgm3: log.info(f'rhoPore_kgm3 differs in position {diff_rhoPore_kgm3[0]:.3f}: {Planet.rhoPore_kgm3[diff_rhoPore_kgm3[0]]:.3f} | {Planet2.rhoPore_kgm3[diff_rhoPore_kgm3[0]]:.3f}')
        elif not all_match:
            log.info('All layer calculations match!')
        else:
            log.info('All values match! Both profiles are identical.')
    else:
        log.warning('Unable to compare layer calculations, as Steps values do not match.')

    return


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
