"""
General runtime configuration parameters.
Overridden by any settings contained within PPBody.py files.
"""
import logging as log
from functools import partial, partialmethod
import os, shutil, time, numpy as np
from Utilities.defineStructs import ParamsStruct, Constants
import matplotlib as mpl
import matplotlib.pyplot as plt
import multiprocessing as mtp

Params = ParamsStruct()
Params.tStart_s = time.time()
Params.VERBOSE = False  # Provides extra runtime messages. Overrides QUIET below
Params.QUIET = False  # Hides all log messages except warnings and errors
Params.printFmt = '[%(levelname)s] %(message)s'  # Format for printing log messages
Params.DEBUG = False  # Special use
# The below flags allow or prevents extrapolation of EOS functions beyond the definition grid.
Params.EXTRAP_ICE = {'Ih':False, 'II':False, 'III':False, 'V':False, 'VI':False, 'Clath':False}
Params.EXTRAP_OCEAN = False
Params.EXTRAP_REF = True  # Allow refprofile extrapolation separate from normal ocean
Params.EXTRAP_SIL = False
Params.EXTRAP_Fe = False

Params.CALC_NEW =         True  # Recalculate profiles? If not, read data from disk and re-plot.
Params.CALC_NEW_REF =     True  # Recalculate reference melting curve densities?
Params.CALC_NEW_INDUCT =  True  # Calculate magnetic induction responses?
Params.CALC_SEISMIC =     True  # Calculate sound speeds and elastic moduli?
Params.CALC_CONDUCT =     True  # Calculate electrical conductivity?
Params.RUN_ALL_PROFILES = True  # Whether to run all PPBody.py files for the named body and plot together
Params.COMPARE =          True  # Whether to plot each new run against other runs from the same body
Params.DO_PARALLEL =      True  # Use multiprocessing module for parallel computation where applicable
Params.SKIP_INNER =       False  # Whether to skip past everything but ocean calculations after MoI matching (for large induction studies)
Params.NO_SAVEFILE =      False  # Whether to prevent printing run outputs to disk. Saves time and disk space for large induction studies.
Params.REDUCED_INDUCT =   True  # Whether to limit number of ocean layers for faster computation of layered induction
Params.INCLUDE_ASYM =     False  # Whether to include asymmetry in the induction conductivity profile based on J2 and C22 values
Params.DISP_LAYERS =      True  # Whether to display layer depths and heat fluxes for user
Params.DISP_TABLE =       False  # Whether to print latex-formatted table
Params.ALWAYS_SHOW_HP =   True  # Whether to force HP ices and clathrates to be shown in DISP_* outputs, even when none are present.
Params.ALWAYS_SHOW_PHI =  True  # Whether to force porosity printout in DISP_* outputs
Params.DEPRECATED =       False  # Whether to allow deprecated code to run. Will often cause errors.

# Plot Settings
Params.SKIP_PLOTS = False  # Whether to skip creation of all plots
Params.PLOT_GRAVITY = False  # Whether to plot Gravity and Pressure
Params.PLOT_HYDROSPHERE = True  # Whether to plot Conductivity with Interior Properties (Hydrosphere)
Params.PLOT_REF = True  # Whether to plot reference melting curve densities on hydrosphere plot
Params.PLOT_SIGS = False  # Make a plot of conductivity as a function of radius to include on hydrosphere plot
Params.PLOT_TRADEOFF = False  # Whether to plot mantle properties tradeoff
Params.PLOT_WEDGE = False  # Whether to plot interior wedge diagram
Params.LEGEND = True  # Whether to include legends

# General figure options
Params.figFormat = 'pdf'
Params.xtn = '.' + Params.figFormat  # Figure file extension. Good options are .eps, .pdf, and .png
Params.defaultFontName = 'STIXGeneral'  # Default font variables--STIX is what is used in Icarus journal submissions
Params.defaultFontCode = 'stix'  # Code name for default font needed in some function calls
Params.backupFont = 'Times New Roman'  # Backup font that looks similar to STIX that most users are likely to have
Params.LegendPosition = 'right'  # Where to place legends when forced
Params.refsInLegend = True  # Whether to include reference profiles in legend
plt.rcParams['font.family'] = 'serif'  # Choose serif font for figures to best match math mode variables in body text
plt.rcParams['font.serif'] = Params.defaultFontName  # Set plots to use the default font
Params.LineStyle = '-'  # Default line style to use on plots
Params.yLim = None  # y axis limits of hydrosphere density in "Conductivity with interior properties" plot

# Reference profiles
# Salinities of reference melting curves in ppt
Params.wRef_ppt = {'none':[0], 'PureH2O':[0],
                   'Seawater':[0, 0.5*Constants.stdSeawater_ppt, Constants.stdSeawater_ppt, 1.5*Constants.stdSeawater_ppt],
                   'MgSO4':[0, 33.3, 66.7, 100],
                   'NH3':[0, 10, 20],
                   'NaCl':[0, 17.5, 35]}
Params.fNameRef = {comp:os.path.join('Thermodynamics', 'RefProfiles', f'{comp}Ref.txt') for comp in Params.wRef_ppt.keys()}
# Style of lines showing reference melting curves of hydrosphere density plot
Params.refLS = {'none':None, 'PureH2O':'-', 'Seawater':':', 'MgSO4':'--', 'NH3':'--', 'NaCl':'--'}
Params.refLW = 0.75
Params.refColor = 'gray'
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

# Default figure sizes
Params.FigSize.vsP = (3,3)
Params.FigSize.vsR = (3,3)
Params.FigSize.vperm = (3,3)
Params.FigSize.vgsks = (3,3)
Params.FigSize.vseis = (3,3)
Params.FigSize.vhydro = (8,5)
Params.FigSize.vgrav = (6,5)
Params.FigSize.vmant = (6,6)
Params.FigSize.vcore = (6,6)
Params.FigSize.vpvt4 = (3,3)
Params.FigSize.vpvt6 = (3,3)
Params.FigSize.vwedg = (3,3)
Params.FigSize.induct = (6,6)

# Color selection
Params.cmap = 'inferno'
# Params.col_contSyn = cfg.cmap(floor(100*cc),:)
# Params.col_contOrb = cfg.cmap(floor( 10*cc),:)
# Params.col_contHrm = cfg.cmap(floor(200*cc),:)

# Params.col_Sw = summer(200)
Params.col_coldestSw = 'c'
# Params.col_midColdSw = cfg.col_Sw(25,:)
# Params.col_middestSw = cfg.col_Sw(50,:)
# Params.col_midWarmSw = cfg.col_Sw(75,:)
Params.col_warmestSw = '#b000ff'
Params.Sw_alt = (0, 175/255, 238/255)

#Params.col_MgSO4 = cool(133)
Params.col_coldestMgSO4 = 'b'
# Params.col_midColdMgSO4 = cfg.col_MgSO4(25,:)
# Params.col_middestMgSO4 = cfg.col_MgSO4(50,:)
# Params.col_midWarmMgSO4 = cfg.col_MgSO4(75,:)
Params.col_warmestMgSO4 = 'xkcd:purpley blue'

# Linestyle and linewidth options
Params.LS_Sw  = '-'  # linestyle for Seawater
Params.LS_Mg = '--'  # linestyle for MgSO4
Params.LS_sp =  ':'  # linestyle for special consideration models
Params.LW_sal = 3  # linewidth for higher salinity
Params.LW_dil = 1  # linewidth for dilute salinity
Params.LW_std = 2  # linewidth for standard salinity
Params.LW_sound = 1.5  # LineWidth for sound speed plots
Params.LW_seism = 1  # LineWidth for seismic plots (Attenuation)

# Wedge color options
Params.Colors.IonosphereTop = [1, 0, 1]
Params.Colors.Ionosphere = [1, 0, 1]
Params.Colors.IonosphereBot = [1, 0, 1]
Params.Colors.IceI = '#d3eefb'
Params.Colors.IceII = '#76b6ff'
Params.Colors.IceIII = '#a8deef'
Params.Colors.IceV = '#83d4f6'
Params.Colors.IceVI = '#cee5ea'
Params.Colors.Clath = '#86bcb8'
Params.Colors.OceanTop = [134/255, 149/255, 201/255] #'#4babdf'
Params.Colors.OceanBot = [45/255, 55/255, 100/255]
Params.Colors.Rock = [101/255, 46/255, 11/255]
Params.Colors.Core = [141/255, 122/255, 121/255]

# Magnetic induction calculation settings
Params.DO_INDUCTOGRAM = True  # Whether to plot an inductogram for the body in question
Params.DO_PER = True  # Convert frequency axes to periods
Params.DO_LEGEND = True  # Whether to force legends
Params.PLOT_FFT = True  # Whether to show plots of fourier space
Params.PLOT_CONTOURS = True  # Contours or surfaces
Params.PLOT_V2021 = True  # Mark the selected ocean/conductivity combos used in Vance et al. 2021
Params.wLims = None  # Minimum and maximum to use for frequency spectrum plots (magnetic induction)

Params.inductOtype = 'rho'  # Type of InductOgram plot to make. Options are "Tb", "phi", "rho", "sigma", where the first 3 are vs. salinity, and sigma is vs. thickness. Sigma/D plot is not self-consistent.
Params.excSelectionCalc = {'synodic': True, 'orbital': True, 'true anomaly': True,  'synodic harmonic': True}  # Which magnetic excitations to include in calculations
Params.excSelectionPlot = {'synodic': True, 'orbital': True, 'true anomaly': False, 'synodic harmonic': True}  # Which magnetic excitations to include in plotting
# Force calculations to be done for each oscillation to be plotted
for osc in Params.excSelectionPlot:
    if Params.excSelectionPlot[osc] and not Params.excSelectionCalc[osc]:
        Params.excSelectionCalc[osc] = True
Params.nwPts = 5  # Resolution for salinity values in ocean salinity vs. other plots
Params.wMin = {'Europa': np.log10(0.05 * Constants.stdSeawater_ppt)}
Params.wMax = {'Europa': np.log10(Constants.stdSeawater_ppt)}
Params.nTbPts = 60  # Resolution for Tb values in ocean salinity/Tb plots
Params.nphiPts = 60  # Resolution for phiRockMax values in ocean salinity/phiMax plots
Params.phiMin = {'Europa': np.log10(0.01)}
Params.phiMax = {'Europa': np.log10(0.75)}
Params.nrhoPts = 6  # Resolution for silicate density values in ocean salinity/rho plots
Params.rhoMin = {'Europa': 3500}
Params.rhoMax = {'Europa': 3700}
Params.nSigmaPts = 50  # Resolution for conductivity values in ocean conductivity/thickness plots
Params.sigmaMin = {'Europa': np.log10(1e-1)}
Params.sigmaMax = {'Europa': np.log10(1e1)}
Params.nDpts = 60  # Resolution for ocean thickness as for conductivity
Params.Dmin = {'Europa': np.log10(1e1)}
Params.Dmax = {'Europa': np.log10(200e1)}
Params.zbFixed_km = {'Europa': 20}
Params.EckhardtSolveMethod = 'RK45'  # Numerical solution method for scipy.integrate.solve_ivp. See https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
Params.rMinODE = 1e3  # Minimum radius to use for numerical solution. Cannot be zero because of singularity at the origin.
Params.oceanInterpMethod = 'linear'  # Interpolation method for determining ocean conductivities when REDUCED_INDUCT is True.
Params.nIntL = 5  # Number of ocean layers to use when REDUCED_INDUCT = 1
Params.nOmegaPts = 100  # Resolution in log frequency space for magnetic excitation spectra
Params.nOmegaFine = 1000  # Fine-spacing resolution for log frequency spectrum
#To be implemented- eventually need some ODE numerical solution parameters
#Params.opts_odeParams = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep', 2e3,'InitialStep',1e-2)
#Params.opts_odeLayers = odeset('RelTol',1e-8, 'AbsTol',1e-10,'MaxStep',10e3,'InitialStep',1e-2)

# Magnetic induction linestyle and linewidth settings
Params.LS_syn = '-'  # linestyle for the synodic period
Params.LS_orb = ':'  # linestyle for the orbital period
Params.LS_hrm = '-.'  # linestyle for the 2nd harmonic of the synodic period
Params.LW_syn = 2  # linewidth for the synodic period
Params.LW_orb = 2  # linewidth for the orbital period
Params.LW_hrm = 2  # linewidth for the 2nd harmonic of the synodic period

# Check if Latex executable is on the path so we can use backup options if Latex is not installed
if shutil.which('latex'):
    plt.rcParams['text.usetex'] = True  # Use Latex interpreter to render text on plots
    # Load in font package in Latex
    plt.rcParams['text.latex.preamble'] = f'\\usepackage{{{Params.defaultFontCode}}}' + \
        r'\usepackage[version=4]{mhchem}' + r'\usepackage{siunitx}' + r'\usepackage{upgreek}'
    Params.TEX_INSTALLED = True
else:
    print('A LaTeX installation was not found. Some plots may have fallback options in labels.')
    plt.rcParams['font.serif'] += ', ' + Params.backupFont  # Set plots to use the default font if installed, or a backup if not
    plt.rcParams['mathtext.fontset'] = Params.defaultFontCode
    Params.TEX_INSTALLED = False
    
# Parallel processing
if Params.DO_PARALLEL:
    Params.maxCores = mtp.cpu_count()
else:
    Params.maxCores = 1
    log.info('DO_PARALLEL is False. Blocking parallel execution.')
# Create parallel printout log level
Params.logParallel = log.WARN + 5
if Params.VERBOSE:
    # Allow warnings to be printed if VERBOSE is selected
    Params.logParallel -= 10
elif Params.QUIET:
    # Allow progress printout to be silenced if QUIET is selected
    Params.logParallel += 10
log.PROFILE = Params.logParallel
log.addLevelName(log.PROFILE, 'PROFILE')
log.Logger.profile = partialmethod(log.Logger.log, log.PROFILE)
log.profile = partial(log.log, log.PROFILE)


class InductionResults:
    def __init__(self):
        self.bodyname = None  # Name of body modeled.
        self.Amp = None  # Amplitude of dipole response (modulus of complex dipole response).
        self.phase = None  # (Positive) phase delay in degrees.
        self.Bix_nT = None  # Induced Bx dipole moments relative to body surface in nT for each excitation.
        self.Biy_nT = None  # Induced By dipole moments relative to body surface in nT for each excitation.
        self.Biz_nT = None  # Induced Bz dipole moments relative to body surface in nT for each excitation.
        self.Texc_hr = None  # Dict of excitation periods modeled.
        self.yName = None  # Name of variable along y axis. Options are "Tb", "phi", "rho", "sigma", where the first 3 are vs. salinity, and sigma is vs. thickness.
        self.x = None  # Values of salinity or conductivity used.
        self.y = None  # Values of yName used.
        self.sigmaMean_Sm = None  # Mean ocean conductivity. Used to map plots vs. salinity onto D/σ plots.
        self.sigmaTop_Sm = None  # Ocean top conductivity. Used to map plots vs. salinity onto D/σ plots.
        self.D_km = None  # Ocean layer thickness in km. Used to map plots vs. salinity onto D/σ plots.
        self.zb_km = None  # Upper ice shell thickness in km.


class ExcitationsList:
    def __init__(self):
        Texc_hr = {}
        # Approximate (.2f) periods to select from excitation spectrum for each body in hr
        Texc_hr['Europa'] = {'synodic':11.23, 'orbital':85.15, 'true anomaly':84.63, 'synodic harmonic':5.62}
        self.Texc_hr = Texc_hr
        self.Induction = InductionResults()

    def __call__(self, bodyname):
        return self.Texc_hr[bodyname]

Excitations = ExcitationsList()