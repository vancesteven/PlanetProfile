"""
General runtime configuration parameters.
Overridden by any settings contained within PPBody.py files.
"""
import logging as log
from functools import partial, partialmethod
import os, time, numpy as np
import multiprocessing as mtp
from Utilities.defineStructs import ParamsStruct, Constants
from MagneticInduction.configInduct import InductParams, SigParams, ExcSpecParams

Params = ParamsStruct()
Params.configVersion = 1  # Integer number for config file version
Params.tStart_s = time.time()
Params.VERBOSE = True  # Provides extra runtime messages. Overrides QUIET below
Params.QUIET = False  # Hides all log messages except warnings and errors
Params.printFmt = '[%(levelname)s] %(message)s'  # Format for printing log messages
Params.DEBUG = False  # Special use
# The below flags allow or prevents extrapolation of EOS functions beyond the definition grid.
Params.EXTRAP_ICE = {'Ih':False, 'II':False, 'III':False, 'V':False, 'VI':False, 'Clath':False}
Params.EXTRAP_OCEAN = False
Params.EXTRAP_REF = True  # Allow refprofile extrapolation separate from normal ocean
Params.EXTRAP_SIL = False
Params.EXTRAP_Fe = False
Params.lookupInterpMethod = 'nearest'  # Interpolation method to use for EOS lookup tables. Options are 'nearest', 'linear', 'cubic'.

Params.CALC_NEW =         True  # Recalculate profiles? If not, read data from disk and re-plot.
Params.CALC_NEW_REF =     True  # Recalculate reference melting curve densities?
Params.CALC_NEW_INDUCT =  True  # Recalculate magnetic induction responses?
Params.CALC_SEISMIC =     True  # Calculate sound speeds and elastic moduli?
Params.CALC_CONDUCT =     True  # Calculate electrical conductivity?
Params.RUN_ALL_PROFILES = False  # Whether to run all PPBody.py files for the named body and plot together
Params.SPEC_FILE =        False  # Whether we are running a specific file or files
Params.COMPARE =          False  # Whether to plot each new run against other runs from the same body
Params.DO_PARALLEL =      True  # Whether to use multiprocessing module for parallel computation where applicable
Params.FORCE_EOS_RECALC = False  # Whether to reuse previously loaded EOS functions for multi-profile runs
Params.SKIP_INNER =       False  # Whether to skip past everything but ocean calculations after MoI matching (for large induction studies)
Params.NO_SAVEFILE =      False  # Whether to prevent printing run outputs to disk. Saves time and disk space for large induction studies.
Params.DISP_LAYERS =      True  # Whether to display layer depths and heat fluxes for user
Params.DISP_TABLE =       True  # Whether to print latex-formatted table
Params.DEPRECATED =       False  # Whether to allow deprecated code to run. Will often cause errors.

# Plot Settings
Params.SKIP_PLOTS = False  # Whether to skip creation of all plots
Params.PLOT_GRAVITY = False  # Whether to plot Gravity and Pressure
Params.PLOT_HYDROSPHERE = True  # Whether to plot Conductivity with Interior Properties (Hydrosphere)
Params.PLOT_REF = True  # Whether to plot reference melting curve densities on hydrosphere plot
Params.PLOT_SIGS = False  # Make a plot of conductivity as a function of radius to include on hydrosphere plot
Params.PLOT_TRADEOFF = False  # Whether to plot mantle properties tradeoff
Params.PLOT_WEDGE = True  # Whether to plot interior wedge diagram
Params.LEGEND = True  # Whether to include legends

# Magnetic induction plot settings
Params.DO_INDUCTOGRAM = False  # Whether to plot an inductogram for the body in question
Params.PLOT_FFT = True  # Whether to show plots of fourier space
Params.INDUCTOGRAM_IN_PROGRESS = False  # Whether we are currently working on constructing an inductogram
Params.COMBINE_BCOMPS = False  # Whether to plot Bx, By, Bz with phase all in one plot, or separate for each comp
Params.Sig = SigParams  # Load general induction settings
Params.Induct = InductParams  # Load inductogram settings
Params.MagSpectrum = ExcSpecParams  # Load excitation spectrum settings

# Reference profile settings
# Salinities of reference melting curves in ppt
Params.wRef_ppt = {'none':[0], 'PureH2O':[0],
                   'Seawater':[0, 0.5*Constants.stdSeawater_ppt, Constants.stdSeawater_ppt, 1.5*Constants.stdSeawater_ppt],
                   'MgSO4':[0, 33.3, 66.7, 100],
                   'NH3':[0, 10, 20],
                   'NaCl':[0, 17.5, 35]}
Params.fNameRef = {comp:os.path.join('Thermodynamics', 'RefProfiles', f'{comp}Ref.txt') for comp in Params.wRef_ppt.keys()}
# Initialize array dicts for refprofiles
Params.Pref_MPa = {}
Params.rhoRef_kgm3 = {}
Params.nRef = {}
Params.nRefPts = {}

# SPICE kernels
Params.spiceTLS = 'naif0012.tls'  # Leap-seconds kernel
Params.spicePCK = 'pck00010.tpc'  # Planetary Constants Kernel from SPICE in order to get body radii
Params.spiceJupiter = 'jup365.bsp'  # Generic kernel for Jupiter + Galilean moons
Params.spiceSaturn = 'sat427.bsp'  # Generic kernel for Saturn + large moons
Params.spiceUranus = 'ura111.bsp'  # Generic kernel for Uranus + large moons
Params.spiceNeptune = 'nep095.bsp'  # Generic kernel for Neptune + large moons
    
# Parallel processing
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


# The following step must always be done last in this file, to allow the user to override the settings above.

