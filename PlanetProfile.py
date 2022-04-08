# Import necessary Python modules
import os, sys
import numpy as np
import importlib
import logging as log
from copy import deepcopy
from distutils.util import strtobool
from os.path import isfile
from glob import glob as FilesMatchingPattern
from config import Params as configParams

# Import all function definitions for this file
from Utilities.SetupInit import SetupInit, SetupFilenames
from Utilities.defineStructs import FigureFilesSubstruct
from Thermodynamics.LayerPropagators import IceLayers, OceanLayers, InnerLayers
from Thermodynamics.FromLiterature.Electrical import ElecConduct
from Seismic import SeismicCalcs
from Utilities.PrintLayerTable import PrintLayerTable
from Plotting.ProfilePlots import GeneratePlots

""" MAIN RUN BLOCK """
def main():

    # Copy global Params settings to local variable to we can add things like filenames
    Params = configParams
    # Set up message logging and apply verbosity level
    if Params.VERBOSE:
        logLevel = log.DEBUG
    elif Params.QUIET:
        logLevel = log.WARN
    else:
        logLevel = log.INFO
    log.basicConfig(level=logLevel, format=Params.printFmt)
    log.debug('Printing verbose runtime messages. Toggle with Params.VERBOSE in config.py.')

    # Command line args
    nArgs = len(sys.argv)
    if nArgs > 1:
        # Body name was passed as command line argument
        bodyname = sys.argv[1]

    else:
        # No command line argument, ask user which body to run
        bodyname = input('Please input body name: ')
        if bodyname == '':
            log.info('No body name entered. Defaulting to Europa.')
            bodyname = 'Europa'

    bodyname = bodyname.capitalize()
    log.info(f'Body name: {bodyname}')

    if Params.DEBUG:
        # Compare to test Matlab output
        Planet, Params = RunPPfile(bodyname, f'PP{bodyname}')
        CompareProfile(Planet, Params, os.path.join('Europa', 'EuropaProfile_Seawater_0WtPpt_Tb269.800K_mat.txt'))
        exit()

    # Get model names from body directory.
    models = FilesMatchingPattern(os.path.join(f'{bodyname}', f'PP{bodyname}*.py'))
    # Splitting on PP and taking the -1 chops the path out, then taking
    # all but the last 3 characters chops the .py extension
    models = [model.split('PP')[-1][:-3] for model in models]

    # Additional command line arguments
    if nArgs > 2:
        if sys.argv[2] == 'clear':
            fNamesToClear = FilesMatchingPattern(os.path.join(bodyname, '*.txt'))
            if len(fNamesToClear) > 0:
                log.info(f'Clearing previous run files for {bodyname}:')
                for fName in fNamesToClear:
                    log.info(f'    {fName}')
                    os.remove(fName)
                log.info(f'{bodyname} files cleared.')
            else:
                log.warning(f'Attempted to remove previous run files for {bodyname}, but found none.')
        elif sys.argv[2] == 'compare':
            log.info('Comparing with other profiles from this body.')
            Params.COMPARE = True
        elif sys.argv[2] == 'all':
            log.info(f'Running all available profiles for {bodyname}.')
            Params.COMPARE = True
        else:
            raise ValueError(f'Unrecognized command: {sys.argv[2]}')
        if nArgs > 3:
            log.warning(f'Too many command line args passed. Ignoring command "{sys.argv[3]}" and any after it.')

    """ Run PlanetProfile """
    PlanetList = np.empty(np.size(models), dtype=object)
    # Run main model first, so that it always appears as 0-index
    PlanetList[0] = importlib.import_module(f'{bodyname}.PP{models[0]}').Planet
    PlanetList[0], Params = PlanetProfile(PlanetList[0], Params)
    if Params.RUN_ALL_PROFILES:
        for i,model in enumerate(models[1:]):
            PlanetList[i+1] = deepcopy(importlib.import_module(f'{bodyname}.PP{model}').Planet)
            PlanetList[i+1], Params = PlanetProfile(PlanetList[i+1], Params)

        # Plot combined figures
        if not Params.SKIP_PLOTS and Params.COMPARE:
            Params.FigureFiles = FigureFilesSubstruct(
                bodyname, os.path.join('figures', f'{bodyname}Comparison'), Params.xtn)
            GeneratePlots(PlanetList, Params)
    else:
        PlanetList = PlanetList[:1]

    """ Post-processing """

    # Loading BodyProfile...txt files to plot them together
    if Params.COMPARE and not Params.RUN_ALL_PROFILES:
        fNamesToCompare = np.array(FilesMatchingPattern(os.path.join(bodyname, f'{bodyname}Profile*.txt')))
        isProfile = [Params.DataFiles.saveFile != fName and 'mantle' not in fName for fName in fNamesToCompare]
        fProfiles = fNamesToCompare[isProfile]
        nCompare = np.size(fProfiles) + 1
        log.info('Loading comparison profiles.')
        CompareList = np.empty(nCompare, dtype=object)
        CompareList[0] = PlanetList[0]
        for i,fName in enumerate(fProfiles):
            CompareList[i+1], _ = ReloadProfile(deepcopy(CompareList[0]), Params, fnameOverride=fName)

        # Plot combined figures
        if not Params.SKIP_PLOTS:
            Params.FigureFiles = FigureFilesSubstruct(
                bodyname, os.path.join('figures', f'{bodyname}Comparison'), Params.xtn)
            GeneratePlots(CompareList, Params)

    if Params.DISP_LAYERS:
        PrintLayerTable(PlanetList, Params)
    if Params.DISP_TABLE:
        PrintLayerTableLatex(PlanetList, Params)

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
        WriteProfile(Planet, Params)
    else:
        # Reload previous run
        Planet, Params = ReloadProfile(Planet, Params)

    if not Params.SKIP_PLOTS:
        # Plotting functions
        GeneratePlots(np.array([Planet]), Params)

    log.info('Run complete!')
    return Planet, Params


def WriteProfile(Planet, Params):
    """ Write out all profile calculations to disk """
    Params.nHeadLines = 41  # Increment as new header lines are added
    with open(Params.DataFiles.saveFile,'w') as f:
        f.write(Planet.label + '\n')
        # Print number of header lines early so we can skip the rest on read-in if we want to
        f.write(f'  nHeadLines = {Params.nHeadLines:d}\n')
        f.write(f'  Ocean salt = {Planet.Ocean.comp}\n')
        f.write(f'  Iron core = {Planet.Do.Fe_CORE}\n')
        f.write(f'  Salinity(ppt) = {Planet.Ocean.wOcean_ppt:.3f}\n')
        f.write(f'  Tb_K = {Planet.Bulk.Tb_K}\n')
        f.write(f'  zb_km = {Planet.zb_km:.3f}\n')
        f.write(f'  zClath_m = {Planet.zClath_m:.3f}\n')
        f.write(f'  D_km = {Planet.D_km:.3f}\n')
        f.write(f'  Pb_MPa = {Planet.Pb_MPa:.3f}\n')
        f.write(f'  PbI_MPa = {Planet.PbI_MPa:.3f}\n')
        f.write(f'  deltaP = {Planet.Ocean.deltaP:.3f}\n')
        f.write(f'  Mtot_kg = {Planet.Mtot_kg:.6e}\n')
        f.write(f'  CMR2mean = {Planet.CMR2mean:.5f}\n')
        f.write(f'  CMR2less = {Planet.CMR2less:.5f}\n')
        f.write(f'  CMR2more = {Planet.CMR2more:.5f}\n')
        f.write(f'  QfromMantle_W = {Planet.Ocean.QfromMantle_W:.6e}\n')
        f.write(f'  phiRockMax = {Planet.Sil.phiRockMax_frac:.3f}\n')
        f.write(f'  RsilMean_m = {Planet.Sil.Rmean_m:.3f}\n')
        f.write(f'  RsilRange_m = {Planet.Sil.Rrange_m:.3f}\n')
        f.write(f'  rhoSil_kgm3 = {Planet.Sil.rhoMean_kgm3:.3f}\n')
        f.write(f'  RcoreMean_m = {Planet.Core.Rmean_m:.3f}\n')
        f.write(f'  RcoreRange_m = {Planet.Core.Rrange_m:.3f}\n')
        f.write(f'  rhoCore_kgm3 = {Planet.Core.rhoMean_kgm3:.3f}\n')
        f.write(f'  MH2O_kg = {Planet.MH2O_kg:.6e}\n')
        f.write(f'  Mrock_kg = {Planet.Mrock_kg:.6e}\n')
        f.write(f'  Mcore_kg = {Planet.Mcore_kg:.6e}\n')
        f.write(f'  Mice_kg = {Planet.Mice_kg:.6e}\n')
        f.write(f'  Msalt_kg = {Planet.Msalt_kg:.6e}\n')
        f.write(f'  MporeSalt_kg = {Planet.MporeSalt_kg:.6e}\n')
        f.write(f'  Mocean_kg = {Planet.Mocean_kg:.6e}\n')
        f.write(f'  Mfluid_kg = {Planet.Mfluid_kg:.6e}\n')
        f.write(f'  MporeFluid_kg = {Planet.MporeFluid_kg:.6e}\n')
        f.write(f'  Steps.nClath = {Planet.Steps.nClath:d}\n')
        f.write(f'  Steps.nIceI = {Planet.Steps.nIceI:d}\n')
        f.write(f'  Steps.nIceIIILitho = {Planet.Steps.nIceIIILitho:d}\n')
        f.write(f'  Steps.nIceVLitho = {Planet.Steps.nIceVLitho:d}\n')
        f.write(f'  Steps.nHydro = {Planet.Steps.nHydro:d}\n')
        f.write(f'  Steps.nSil = {Planet.Steps.nSil:d}\n')
        f.write(f'  Steps.nCore = {Planet.Steps.nCore:d}\n')
        f.write(f' '.join(['P (MPa)'.ljust(24),
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
                           'MLayer (kg)']) + '\n')
        # Now print the columnar data
        for i in range(Planet.Steps.nTotal):
            line = \
                f'{Planet.P_MPa[i]:24.17e} ' + \
                f'{Planet.T_K[i]:24.17e} ' + \
                f'{Planet.r_m[i]:24.17e} ' + \
                f'{Planet.phase[i]:8d} ' + \
                f'{Planet.rho_kgm3[i]:24.17e} ' + \
                f'{Planet.Cp_JkgK[i]:24.17e} ' + \
                f'{Planet.alpha_pK[i]:24.17e} ' + \
                f'{Planet.g_ms2[i]:24.17e} ' + \
                f'{Planet.phi_frac[i]:24.17e} ' + \
                f'{Planet.sigma_Sm[i]:24.17e} ' + \
                f'{Planet.kTherm_WmK[i]:24.17e} ' + \
                f'{Planet.Seismic.VP_kms[i]:24.17e} ' + \
                f'{Planet.Seismic.VS_kms[i]:24.17e} ' + \
                f'{Planet.Seismic.QS[i]:24.17e} ' + \
                f'{Planet.Seismic.KS_GPa[i]:24.17e} ' + \
                f'{Planet.Seismic.GS_GPa[i]:24.17e} ' + \
                f'{Planet.Ppore_MPa[i]:24.17e} ' + \
                f'{Planet.rhoMatrix_kgm3[i]:24.17e} ' + \
                f'{Planet.rhoPore_kgm3[i]:24.17e} ' + \
                f'{Planet.MLayer_kg[i]:24.17e}\n '
            f.write(line)

    # Write out data from core/mantle trade
    with open(Params.DataFiles.mantCoreFile, 'w') as f:
        f.write(' '.join(['RsilTrade (m)'.ljust(24),
                          'RcoreTrade (m)'.ljust(24),
                          'rhoSilTrade (kg/m3)']) + '\n')
        for i in range(np.size(Planet.Sil.Rtrade_m)):
            line = \
                f'{Planet.Sil.Rtrade_m[i]:24.17e} ' + \
                f'{Planet.Core.Rtrade_m[i]:24.17e} ' + \
                f'{Planet.Sil.rhoTrade_kgm3[i]:24.17e}\n '
            f.write(line)

    log.info(f'Profile saved to file: {Params.DataFiles.saveFile}')
    return


def ReloadProfile(Planet, Params, fnameOverride=None):
    """ Reload previously saved PlanetProfile run from disk """

    if fnameOverride is not None:
        Params.DataFiles.saveFile = fnameOverride
        Params.DataFiles.mantCoreFile = f'{fnameOverride[:-4]}_mantleCore.txt'
        Params.DataFiles.mantPermFile = f'{fnameOverride[:-4]}_mantlePerm.txt'
        nSkip = len(os.path.join(Planet.name, f'{Planet.name}Profile_'))
        Planet.label = fnameOverride[nSkip:-4]
    else:
        Params.DataFiles, Params.FigureFiles = SetupFilenames(Planet, Params)
    log.info(f'Reloading previously saved run from file: {Params.DataFiles.saveFile}')
    log.debug(f'Steps.n settings from PP{Planet.name}.py will be ignored.')
    if not isfile(Params.DataFiles.saveFile):
        raise ValueError('Params.CALC_NEW is set to False in config.py but the reload file was not found.\n' +
                         'Re-run with CALC_NEW set to True to generate the profile.')

    with open(Params.DataFiles.saveFile) as f:
        # Get legend label for differentiating runs
        Planet.label = f.readline().strip()
        # Get number of header lines to read in from (and skip for columnar data)
        Params.nHeadLines = int(f.readline().split('=')[-1])
        # Get dissolved salt supposed for ocean (present in filename, but this is intended for future-proofing when we move to a database lookup)
        Planet.Ocean.comp = f.readline().split('=')[-1].strip()
        # Get whether iron core is modeled
        Planet.Do.Fe_CORE = bool(strtobool(f.readline().split('=')[-1].strip()))
        # Get float values from header
        Planet.Ocean.wOcean_ppt, Planet.Bulk.Tb_K, Planet.zb_km, Planet.zClath_m, Planet.D_km, \
        Planet.Pb_MPa, Planet.PbI_MPa, Planet.Ocean.deltaP, Planet.Mtot_kg, Planet.CMR2mean, Planet.CMR2less,\
        Planet.CMR2more, Planet.Ocean.QfromMantle_W, Planet.Sil.phiRockMax_frac, Planet.Sil.Rmean_m, \
        Planet.Sil.Rrange_m, Planet.Sil.rhoMean_kgm3, Planet.Core.Rmean_m, Planet.Core.Rrange_m, \
        Planet.Core.rhoMean_kgm3, Planet.MH2O_kg, Planet.Mrock_kg, Planet.Mcore_kg, Planet.Mice_kg, \
        Planet.Msalt_kg, Planet.MporeSalt_kg, Planet.Mocean_kg, Planet.Mfluid_kg, Planet.MporeFluid_kg \
            = (float(f.readline().split('=')[-1]) for _ in range(29))
        # Get integer values from header (nSteps values)
        Planet.Steps.nClath, Planet.Steps.nIceI, \
        Planet.Steps.nIceIIILitho, Planet.Steps.nIceVLitho, \
        Planet.Steps.nHydro, Planet.Steps.nSil, Planet.Steps.nCore \
            = (int(f.readline().split('=')[-1]) for _ in range(7))

    Planet.Steps.nIbottom = Planet.Steps.nClath + Planet.Steps.nIceI
    Planet.Steps.nIIIbottom = Planet.Steps.nIbottom + Planet.Steps.nIceIIILitho
    Planet.Steps.nVbottom = Planet.Steps.nIIIbottom + Planet.Steps.nIceVLitho
    Planet.Steps.nTotal = Planet.Steps.nHydro + Planet.Steps.nSil + Planet.Steps.nCore
    # Read in columnar data that follows header lines -- full-body
    Planet.P_MPa, Planet.T_K, Planet.r_m, Planet.phase, Planet.rho_kgm3, Planet.Cp_JkgK, Planet.alpha_pK, \
    Planet.g_ms2, Planet.phi_frac, Planet.sigma_Sm, Planet.kTherm_WmK, Planet.Seismic.VP_kms, Planet.Seismic.VS_kms,\
    Planet.Seismic.QS, Planet.Seismic.KS_GPa, Planet.Seismic.GS_GPa, Planet.Ppore_MPa, Planet.rhoMatrix_kgm3, \
    Planet.rhoPore_kgm3, Planet.MLayer_kg \
        = np.loadtxt(Params.DataFiles.saveFile, skiprows=Params.nHeadLines, unpack=True)
    Planet.r_m = np.concatenate((Planet.r_m, [0]))
    Planet.z_m = Planet.Bulk.R_m - Planet.r_m
    Planet.phase = Planet.phase.astype(np.int_)

    # Read in data for core/mantle trade
    Planet.Sil.Rtrade_m, Planet.Core.Rtrade_m, Planet.Sil.rhoTrade_kgm3, \
        = np.loadtxt(Params.DataFiles.mantCoreFile, skiprows=1, unpack=True)

    # Read in data for mantle permeability properties
    #Planet.Sil.perm1, Planet.Sil.perm2, Planet.Sil.perm3, Planet.Sil.perm4, \
    #Planet.Sil.perm5 \
    #    = np.loadtxt(Params.DataFiles.permFile, skiprows=1, unpack=True)

    return Planet, Params


def RunPPfile(bodyname, fName):
    """ Loads the settings in bodyname/fName.py to run or reload a specific model. """
    bodyname = bodyname.capitalize()
    Planet = deepcopy(importlib.import_module(f'{bodyname}.{fName}').Planet)
    Planet, Params = PlanetProfile(Planet, configParams)

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


if __name__ == '__main__': main()
