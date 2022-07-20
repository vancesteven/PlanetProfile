""" Load in default settings, then overwrite them with the user's settings in top-level dir. """
from warnings import warn
import os
import time
import matplotlib
import spiceypy as spice
import logging
import multiprocessing as mtp
from functools import partial, partialmethod
import MoonMag.symmetry_funcs, MoonMag.asymmetry_funcs

from PlanetProfile.defaultConfig import \
    Params, \
    ExploreParams, \
    configVersion
from configPP import \
    Params as userParams, \
    ExploreParams as userExploreParams, \
    configVersion as userConfigVersion
if hasattr(userParams,'spiceTLS') and hasattr(userParams,'spiceDir'):
    userLSK = os.path.join(userParams.spiceDir, userParams.spiceTLS)
    if not os.path.isfile(userLSK):
        raise FileNotFoundError(f'Leapseconds kernel was not found at {userLSK}. This likely means PlanetProfile ' +
                                f'has not been fully installed. Run the install script with the following command:\n' +
                                f'python -m PlanetProfile.install PPinstall')
    spice.furnsh(userLSK)
else:
    defLSK = os.path.join(Params.spiceDir, Params.spiceTLS)
    if not os.path.isfile(defLSK):
        raise FileNotFoundError(f'Leapseconds kernel was not found at {defLSK}. This likely means PlanetProfile ' +
                                f'has not been fully installed. Run the install script with the following command:\n' +
                                f'python -m PlanetProfile.install PPinstall')
    spice.furnsh(defLSK)
from PlanetProfile.MagneticInduction.defaultConfigInduct import \
    SigParams, \
    ExcSpecParams, \
    InductParams, \
    configInductVersion
from configPPinduct import \
    SigParams as userSigParams, \
    ExcSpecParams as userExcSpecParams, \
    InductParams as userInductParams, \
    testBody as userTestBody, \
    configInductVersion as userConfigInductVersion
from PlanetProfile.Plotting.defaultConfigPlots import \
    Color, \
    Style, \
    FigLbl, \
    FigSize, \
    FigMisc, \
    configPlotsVersion
from configPPplots import \
    Color as userColor, \
    Style as userStyle, \
    FigLbl as userFigLbl, \
    FigSize as userFigSize, \
    FigMisc as userFigMisc, \
    configPlotsVersion as userConfigPlotsVersion

# Check config file versions and warn user if they differ
if configVersion != userConfigVersion:
    warn(f'User configPP file is version {userConfigVersion}, but the default file is ' +
                      f'version {configVersion}. Some settings may be missing; default values will be used. ' +
                      'To align the file version, delete configPP.py and run again, or execute reset.py.')
if configInductVersion != userConfigInductVersion:
    warn(f'User configPPinduct file is version {userConfigInductVersion}, but the default file is ' +
                      f'version {configInductVersion}. Some settings may be missing; default values will be used. ' +
                      'To align the file version, delete configPPinduct.py and run again, or execute reset.py.')
if configPlotsVersion != userConfigPlotsVersion:
    warn(f'User configPPplots file is version {userConfigPlotsVersion}, but the default file is ' +
                      f'version {configPlotsVersion}. Some settings may be missing; default values will be used. ' +
                      'To align the file version, delete configPPplots.py and run again, or execute reset.py.')

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

stream = logging.StreamHandler()
stream.setFormatter(logging.Formatter(Params.printFmt))
log.setLevel(logLevel)
log.addHandler(stream)
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('PIL').setLevel(logging.WARNING)
logging.getLogger('MoonMag').setLevel(logLevelMoonMag)
log.debug('Printing verbose runtime messages. Toggle with Params.VERBOSE in configPP.py.')

# Parallel processing toggles
if Params.DO_PARALLEL:
    Params.maxCores = mtp.cpu_count()
else:
    Params.maxCores = 1
    log.info('DO_PARALLEL is False. Blocking parallel execution.')

# Add Test body settings to InductParams
[getattr(InductParams, attr).update({'Test': getattr(InductParams, attr)[userTestBody]})
    for attr in ['wMin', 'wMax', 'TbMin', 'TbMax', 'phiMin', 'phiMax', 'rhoMin',
                 'rhoMax', 'sigmaMin', 'sigmaMax', 'Dmin', 'Dmax', 'zbFixed_km']]

# Force calculations to be done for each oscillation to be plotted in inductograms
for osc in InductParams.excSelectionPlot:
    if InductParams.excSelectionPlot[osc] and not InductParams.excSelectionCalc[osc]:
        InductParams.excSelectionCalc[osc] = True

# Load induction and excitation spectrum settings into Params so they
# can be passed around together
Params.Sig = SigParams
Params.Induct = InductParams
Params.MagSpectrum = ExcSpecParams
Params.Explore = ExploreParams
