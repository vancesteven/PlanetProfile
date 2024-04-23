""" Load in default settings, then overwrite them with the user's settings in top-level dir. """
from warnings import warn
import os, sys
import time
import matplotlib
import spiceypy as spice
import logging
import multiprocessing as mtp
from functools import partial, partialmethod
from PlanetProfile import _DefaultList
import MoonMag.symmetry_funcs, MoonMag.asymmetry_funcs

# Fetch version numbers first to warn user about compatibility
from PlanetProfile.defaultConfig import configVersion
from configPP import configVersion as userConfigVersion

# Check config file versions and warn user if they differ
if configVersion != userConfigVersion:
    warn(f'User configPP file is version {userConfigVersion}, but the default file is ' +
                      f'version {configVersion}. Some settings may be missing; default values will be used. ' +
                      'To align the file version, delete configPP.py and run again, or execute reset.py ' +
                      'with python -m PlanetProfile.reset')

# Grab config settings first to load LSK for later adjustments with SPICE
from PlanetProfile.defaultConfig import configAssign
from configPP import configAssign as userConfigAssign
Params, ExploreParams = configAssign()
userParams, userExploreParams = userConfigAssign()

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

from PlanetProfile.MagneticInduction.defaultConfigInduct import configInductVersion
from configPPinduct import configInductVersion as userConfigInductVersion
from PlanetProfile.Plotting.defaultConfigPlots import configPlotsVersion
from configPPplots import configPlotsVersion as userConfigPlotsVersion

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

from PlanetProfile.MagneticInduction.defaultConfigInduct import inductAssign
from configPPinduct import inductAssign as userInductAssign
from PlanetProfile.TrajecAnalysis.defaultConfigTrajec import trajecAssign
from configPPtrajec import trajecAssign as userTrajecAssign
from PlanetProfile.Plotting.defaultConfigPlots import plotAssign
from configPPplots import plotAssign as userPlotAssign

SigParams, ExcSpecParams, InductParams, _ = inductAssign()
userSigParams, userExcSpecParams, userInductParams, userTestBody = userInductAssign()
TrajecParams = trajecAssign()
userTrajecParams = userTrajecAssign()
Color, Style, FigLbl, FigSize, FigMisc = plotAssign()
userColor, userStyle, userFigLbl, userFigSize, userFigMisc = userPlotAssign()

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

# Execute necessary adjustments and add non-settings to Params objects
Params.tStart_s = time.time()

# Get RefProfile file names from the lists set in config file(s)
Params.fNameRef = {comp:f'{comp}Ref.txt' for comp in Params.wRef_ppt.keys()}
# Initialize array dicts for refprofiles
Params.Pref_MPa = {}
Params.rhoRef_kgm3 = {}
Params.nRef = {}
Params.nRefPts = {}

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
    
# Set up message logging and apply verbosity level
log = logging.getLogger('PlanetProfile')
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

stream = logging.StreamHandler(sys.stdout)
stream.setFormatter(logging.Formatter(Params.printFmt))
log.setLevel(logLevel)
log.addHandler(stream)
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('PIL').setLevel(logging.WARNING)
logging.getLogger('MoonMag').setLevel(logLevelMoonMag)
logging.getLogger('lbftd').setLevel(logLevelLBF)
log.debug('Printing verbose runtime messages. Toggle with Params.VERBOSE in configPP.py.')

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

# Load sub-settings into Params so they can be passed around together
Params.Sig = SigParams
Params.Induct = InductParams
Params.MagSpectrum = ExcSpecParams
Params.Explore = ExploreParams
Params.Trajec = TrajecParams
