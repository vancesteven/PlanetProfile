""" Load in default settings, then overwrite them with the user's settings in top-level dir. """
from warnings import warn
import os, sys
import time
import matplotlib
import spiceypy as spice
import logging
import multiprocessing as mtp
from functools import partial, partialmethod
from PlanetProfile import _DefaultList, _ROOT
import importlib
import MoonMag.symmetry_funcs, MoonMag.asymmetry_funcs

""" Get default config settings """
from PlanetProfile.defaultConfig import configAssign, configVersion
from PlanetProfile.MagneticInduction.defaultConfigInduct import inductAssign, configInductVersion
from PlanetProfile.TrajecAnalysis.defaultConfigTrajec import trajecAssign, configTrajecVersion
from PlanetProfile.Plotting.defaultConfigPlots import plotAssign, configPlotsVersion
from PlanetProfile.CustomSolution.defaultConfigCustomSolution import customSolutionAssign, configCustomSolutionVersion
from PlanetProfile.Gravity.defaultConfigGravity import gravityAssign, configGravityVersion
from PlanetProfile.MonteCarlo.defaultConfigMonteCarlo import montecarloAssign, configMonteCarloVersion
from PlanetProfile.Inversion.defaultConfigInversion import inversionAssign, configInversionVersion
Params, ExploreParams = configAssign()
defLSK = os.path.join(Params.spiceDir, Params.spiceTLS)
if not os.path.isfile(defLSK):
    raise FileNotFoundError(f'Leapseconds kernel was not found at {defLSK}. This likely means PlanetProfile ' +
                            f'has not been fully installed. Run the install script with the following command:\n' +
                            f'python -m PlanetProfile.install')
spice.furnsh(defLSK)
SigParams, ExcSpecParams, InductParams, _ = inductAssign()
TrajecParams = trajecAssign()
CustomSolutionParams = customSolutionAssign()
GravityParams = gravityAssign()
MonteCarloParams = montecarloAssign()
InversionParams = inversionAssign()
Color, Style, FigLbl, FigSize, FigMisc = plotAssign()


# Create parallel printout log level
logging.PROFILE = logging.WARN + 5
Params.logParallel = logging.PROFILE + 0
logging.addLevelName(logging.PROFILE, 'PROFILE')
logging.Logger.profile = partialmethod(logging.Logger.log, logging.PROFILE)
logging.profile = partial(logging.log, logging.PROFILE)
if Params.VERBOSE:
    # Allow debug messages to be printed if VERBOSE is selected
    Params.logParallel -= 30
elif Params.QUIET:
    # Allow progress printout to be silenced if QUIET is selected
    Params.logParallel += 10
    
# === Custom Log Level: PERFORMANCE ===
logging.TIMING = logging.INFO + 5
logging.addLevelName(logging.TIMING, 'TIMING')
logging.Logger.timing = partialmethod(logging.Logger.log, logging.TIMING)
logging.timing = partial(logging.log, logging.TIMING)
    
# Set up message logging, child loggers, and apply verbosity level
log = logging.getLogger('PlanetProfile')
timeLog = logging.getLogger('PlanetProfile.Timing')
if Params.VERBOSE:
    logLevel = logging.DEBUG
elif Params.QUIET:
    logLevel = logging.WARN
else:
    logLevel = logging.INFO
if Params.QUIET_MOONMAG:
    logLevelMoonMag = logging.WARNING
else:
    logLevelMoonMag = logLevel
if Params.QUIET_LBF:
    logLevelLBF = logging.ERROR
else:
    logLevelLBF = logLevel
if Params.TIMING:
    timingLogLevel = logging.TIMING
else:
    timingLogLevel = logging.TIMING + 1

stream = logging.StreamHandler(sys.stdout)
stream.setFormatter(logging.Formatter(Params.printFmt))
log.setLevel(logLevel)
log.addHandler(stream)
timeLog.setLevel(timingLogLevel)
timeStream = logging.StreamHandler(sys.stdout)
timeStream.setFormatter(logging.Formatter(Params.printFmt))
timeLog.addHandler(timeStream)
timeLog.propagate = False # Prevent double logging to parent handlers
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('PIL').setLevel(logging.WARNING)
logging.getLogger('MoonMag').setLevel(logLevelMoonMag)
logging.getLogger('lbftd').setLevel(logLevelLBF)
log.debug('Printing verbose runtime messages. Toggle with Params.VERBOSE in configPP.py.')


""" Load user settings to allow for configuration """
def loadUserSettings(configModule: str = ''):
    if configModule != '':
        configModule = 'UserConfigs.' + configModule
    else:
        configModule = 'UserConfigs'
    configPP = importlib.import_module(configModule + '.configPP')
    userConfigAssign = configPP.configAssign
    userParams, userExploreParams = userConfigAssign()
    userConfigVersion = configPP.configVersion
    # Check config file versions and warn user if they differ
    if configVersion != userConfigVersion:
        warn(f'User configPP file is version {userConfigVersion}, but the default file is ' +
                        f'version {configVersion}. Some settings may be missing; default values will be used. ' +
                        'To align the file version, delete configPP.py and run again, or execute reset.py ' +
                        'with python -m PlanetProfile.reset')
    if hasattr(userParams, 'spiceTLS') and hasattr(userParams, 'spiceDir'):
        userLSK = os.path.join(userParams.spiceDir, userParams.spiceTLS)
        if not os.path.isfile(userLSK):
            raise FileNotFoundError(f'Leapseconds kernel was not found at {userLSK}. This likely means PlanetProfile ' +
                                    f'has not been fully installed. Run the install script with the following command:\n' +
                                    f'python -m PlanetProfile.install')
        spice.furnsh(userLSK)
    else:
        defLSK = os.path.join(Params.spiceDir, Params.spiceTLS)
        if not os.path.isfile(defLSK):
            raise FileNotFoundError(f'Leapseconds kernel was not found at {defLSK}. This likely means PlanetProfile ' +
                                    f'has not been fully installed. Run the install script with the following command:\n' +
                                    f'python -m PlanetProfile.install')
        spice.furnsh(defLSK)
    configInduct = importlib.import_module(configModule + '.configPPinduct')
    userConfigInductVersion = configInduct.configInductVersion
    configTrajec = importlib.import_module(configModule + '.configPPtrajec')
    userConfigTrajecVersion = configTrajec.configTrajecVersion
    configPlots = importlib.import_module(configModule + '.configPPplots')
    userConfigPlotsVersion = configPlots.configPlotsVersion
    configCustomSolution = importlib.import_module(configModule + '.configPPcustomsolution')
    userConfigCustomSolutionVersion = configCustomSolution.configCustomSolutionVersion
    configGravity = importlib.import_module(configModule + '.configPPgravity')
    userConfigGravityVersion = configGravity.configGravityVersion
    configMonteCarlo = importlib.import_module(configModule + '.configPPmontecarlo')
    userConfigMonteCarloVersion = configMonteCarlo.configMonteCarloVersion
    configInversion = importlib.import_module(configModule + '.configPPinversion')
    userConfigInversionVersion = configInversion.configInversionVersion

    # Check sub-config file versions and warn user if they differ
    if configInductVersion != userConfigInductVersion:
        warn(f'User configPPinduct file is version {userConfigInductVersion}, but the default file is ' +
            f'version {configInductVersion}. Some settings may be missing; default values will be used. ' +
            f'To align the file version, delete configPPinduct.py and run again, or execute reset.py ' +
            f'with python -m PlanetProfile.reset')
    if configPlotsVersion != userConfigPlotsVersion:
        warn(f'User configPPplots file is version {userConfigPlotsVersion}, but the default file is ' +
            f'version {configPlotsVersion}. Some settings may be missing; default values will be used. ' +
            f'To align the file version, delete configPPplots.py and run again, or execute reset.py ' +
            f'with python -m PlanetProfile.reset')
    if configCustomSolutionVersion != userConfigCustomSolutionVersion:
        warn(f'User configPPcustomsolution file is version {userConfigCustomSolutionVersion}, but the default file is ' +
            f'version {configCustomSolutionVersion}. Some settings may be missing; default values will be used. ' +
            f'To align the file version, delete configPPcustomsolution.py and run again, or execute reset.py ' +
            f'with python -m PlanetProfile.reset')
    if configGravityVersion != userConfigGravityVersion:
        warn(f'User configPPgravity file is version {userConfigGravityVersion}, but the default file is ' +
            f'version {configGravityVersion}. Some settings may be missing; default values will be used. ' +
            f'To align the file version, delete configPPgravity.py and run again, or execute reset.py ' +
            f'with python -m PlanetProfile.reset')
    if configMonteCarloVersion != userConfigMonteCarloVersion:
        warn(f'User configPPmontecarlo file is version {userConfigMonteCarloVersion}, but the default file is ' +
            f'version {configMonteCarloVersion}. Some settings may be missing; default values will be used. ' +
            f'To align the file version, delete configPPmontecarlo.py and run again, or execute reset.py ' +
            f'with python -m PlanetProfile.reset')
    if configInversionVersion != userConfigInversionVersion:
        warn(f'User configPPinversion file is version {userConfigInversionVersion}, but the default file is ' +
            f'version {configInversionVersion}. Some settings may be missing; default values will be used. ' +
            f'To align the file version, delete configPPinversion.py and run again, or execute reset.py ' +
            f'with python -m PlanetProfile.reset')
    userInductAssign = configInduct.inductAssign
    userTrajecAssign = configTrajec.trajecAssign
    userPlotAssign = configPlots.plotAssign
    userCustomSolutionAssign = configCustomSolution.customSolutionAssign
    userGravityAssign = configGravity.gravityAssign
    userMontecarloAssign = configMonteCarlo.montecarloAssign
    userInversionAssign = configInversion.inversionAssign

    userSigParams, userExcSpecParams, userInductParams, userTestBody = userInductAssign()
    userTrajecParams = userTrajecAssign()
    userColor, userStyle, userFigLbl, userFigSize, userFigMisc = userPlotAssign()
    userCustomSolutionParams = userCustomSolutionAssign()
    userGravityParams = userGravityAssign()
    userMonteCarloParams = userMontecarloAssign()
    userInversionParams = userInversionAssign()

    # Load user settings to allow for configuration
    for attr, value in userParams.__dict__.items():
        setattr(Params, attr, value)
    for attr, value in userExploreParams.__dict__.items():
        setattr(ExploreParams, attr, value)
    for attr, value in userSigParams.__dict__.items():
        setattr(SigParams, attr, value)
    for attr, value in userExcSpecParams.__dict__.items():
        setattr(ExcSpecParams, attr, value)
    for attr, value in userInductParams.__dict__.items():
        setattr(InductParams, attr, value)
    for attr, value in userGravityParams.__dict__.items():
        setattr(GravityParams, attr, value)
    for attr, value in userMonteCarloParams.__dict__.items():
        setattr(MonteCarloParams, attr, value)
    for attr, value in userInversionParams.__dict__.items():
        setattr(InversionParams, attr, value)

    for attr, value in userTrajecParams.__dict__.items():
        setattr(TrajecParams, attr, value)

    for attr, value in userColor.__dict__.items():
        setattr(Color, attr, value)
    for attr, value in userStyle.__dict__.items():
        setattr(Style, attr, value)
    for attr, value in userFigLbl.__dict__.items():
        setattr(FigLbl, attr, value)
    for attr, value in userFigSize.__dict__.items():
        setattr(FigSize, attr, value)
    for attr, value in userFigMisc.__dict__.items():
        setattr(FigMisc, attr, value)
    for attr, value in userCustomSolutionParams.__dict__.items():
        setattr(CustomSolutionParams, attr, value)

    # Execute necessary adjustments and add non-settings to Params objects
    Params.tStart_s = time.time()

    # Get RefProfile file names from the lists set in config file(s)
    Params.fNameRef = {comp:f'{comp}Ref.txt' for comp in Params.wRef_ppt.keys()}
    # Initialize array dicts for refprofiles
    Params.Pref_MPa = {}
    Params.rhoRef_kgm3 = {}
    Params.nRef = {}
    Params.nRefPts = {}

    # Parallel processing toggles
    if Params.DO_PARALLEL:
        Params.maxCores = mtp.cpu_count()
    else:
        Params.maxCores = 1
        log.info('DO_PARALLEL is False. Blocking parallel execution.')

    # Add Test body settings to InductParams
    inductOgramAttr = ['wMin', 'wMax', 'TbMin', 'TbMax', 'phiMin', 'phiMax', 'rhoMin',
                    'rhoMax', 'sigmaMin', 'sigmaMax', 'Dmin', 'Dmax', 'zbFixed_km']
    [getattr(InductParams, attr).update({'Test': getattr(InductParams, attr)[userTestBody]})
        for attr in inductOgramAttr]

    # Force calculations to be done for each oscillation to be plotted in inductograms
    for osc in InductParams.excSelectionPlot:
        if InductParams.excSelectionPlot[osc] and not InductParams.excSelectionCalc[osc]:
            InductParams.excSelectionCalc[osc] = True
        
    # Assign missing values with defaults in inductOgram settings
    Tnames = list(InductParams.cLevels['Test'].keys())
    zNames = ['Amp', 'Bx', 'By', 'Bz']
    for bodyname in _DefaultList:
        for attr in inductOgramAttr:
            if bodyname not in getattr(InductParams, attr).keys():
                getattr(InductParams, attr).update({'Test': getattr(InductParams, attr)[userTestBody]})
        
        if bodyname not in InductParams.cLevels.keys():
            InductParams.cLevels[bodyname] = {}
            for Tname in Tnames:
                InductParams.cLevels[bodyname][Tname] = {zName: None for zName in zNames}
        if bodyname not in InductParams.cfmt.keys():
            InductParams.cfmt[bodyname] = {}
            for Tname in Tnames:
                InductParams.cfmt[bodyname][Tname] = {zName: None for zName in zNames}

    # Assign extrapolation settings for Mixed Clathrates in dictionary
    mixedClathExtrap = ['MixedClathrateIh', 'MixedClathrateII', 'MixedClathrateIII', 'MixedClathrateV', 'MixedClathrateVI']
    icePhase = ['Ih', 'II', 'III', 'V', 'VI']
    for i, phase in enumerate(mixedClathExtrap):
        if phase not in Params.EXTRAP_ICE.keys():
            Params.EXTRAP_ICE[phase] = Params.EXTRAP_ICE['Clath'] and Params.EXTRAP_ICE[icePhase[i]]
    """ 
    Check CustomSolutionConfig inputs are valid, set file paths, and save global reference to object so Rkt file can use
    """
    # Ensure frezchem database is valid
    CustomSolutionParams.setPaths(_ROOT)
    for file_name in os.listdir(CustomSolutionParams.databasePath):
        if file_name == CustomSolutionParams.FREZCHEM_DATABASE:
            break
    else:
        log.warning(
            "Input frezchem database does not match any of the available saved files.\nCheck that the input is properly spelled and has .dat at end. Using default frezchem.dat file")
        CustomSolutionParams.FREZCHEM_DATABASE = "frezchem.dat"
    # Check the unit is 'g' or 'mol' (g - grams, mol - mols)
    if not CustomSolutionParams.SPECIES_CONCENTRATION_UNIT == "g" and not CustomSolutionParams.SPECIES_CONCENTRATION_UNIT == "mol":
        log.warning(
            "Input species concentration unit is not valid. Check that it is either g or mol. Using mol as default")
        CustomSolutionParams.SPECIES_CONCENTRATION_UNIT = "mol"
    # Load sub-settings into Params so they can be passed around together
    Params.Sig = SigParams
    Params.Induct = InductParams
    Params.MagSpectrum = ExcSpecParams
    Params.Explore = ExploreParams
    Params.Trajec = TrajecParams
    Params.CustomSolution = CustomSolutionParams
    Params.Gravity = GravityParams
    Params.MonteCarlo = MonteCarloParams
    Params.Inversion = InversionParams

loadUserSettings()