""" Load in default settings, then overwrite them with the user's settings in top-level dir. """
import logging as log
import multiprocessing as mtp
from functools import partial, partialmethod
import matplotlib
import time

from PlanetProfile.defaultConfig import \
    Params, \
    configVersion
from configPP import \
    Params as userParams, \
    configVersion as userConfigVersion
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
    raise UserWarning(f'User configPP file is version {userConfigVersion}, but the default file is ' +
                      f'version {configVersion}. Some settings may be missing; default values will be used. ' +
                      'To align the file version, delete configPP.py and run again, or execute reset.py.')
if configInductVersion != userConfigInductVersion:
    raise UserWarning(f'User configPPinduct file is version {userConfigInductVersion}, but the default file is ' +
                      f'version {configInductVersion}. Some settings may be missing; default values will be used. ' +
                      'To align the file version, delete configPPinduct.py and run again, or execute reset.py.')
if configPlotsVersion != userConfigPlotsVersion:
    raise UserWarning(f'User configPPplots file is version {userConfigPlotsVersion}, but the default file is ' +
                      f'version {configPlotsVersion}. Some settings may be missing; default values will be used. ' +
                      'To align the file version, delete configPPplots.py and run again, or execute reset.py.')

# Load user settings to allow for configuration
for attr, value in userParams.__dict__.items():
    setattr(Params, attr, value)
    
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

# Parallel processing toggles
if Params.DO_PARALLEL:
    Params.maxCores = mtp.cpu_count()
else:
    Params.maxCores = 1
    log.info('DO_PARALLEL is False. Blocking parallel execution.')
    
# Create parallel printout log level
log.PROFILE = log.WARN + 5
Params.logParallel = log.PROFILE + 0
log.addLevelName(log.PROFILE, 'PROFILE')
log.Logger.profile = partialmethod(log.Logger.log, log.PROFILE)
log.profile = partial(log.log, log.PROFILE)
if Params.VERBOSE:
    # Allow debug messages to be printed if VERBOSE is selected
    Params.logParallel -= 30
elif Params.QUIET:
    # Allow progress printout to be silenced if QUIET is selected
    Params.logParallel += 10
    
# Set up message logging and apply verbosity level
if Params.VERBOSE:
    logLevel = log.DEBUG
elif Params.QUIET:
    logLevel = log.WARN
else:
    logLevel = log.INFO
root = log.getLogger()
if root.handlers:
    for handler in root.handlers:
        root.removeHandler(handler)
log.basicConfig(level=logLevel, format=Params.printFmt)
log.getLogger('matplotlib').setLevel(log.WARNING)
log.debug('Printing verbose runtime messages. Toggle with Params.VERBOSE in configPP.py.')

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
