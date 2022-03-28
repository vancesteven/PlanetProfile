"""
General runtime configuration parameters.
Overridden by any settings contained within PPBody.py files.
"""

import shutil
from Utilities.defineStructs import ParamsStruct
import matplotlib as mpl
import matplotlib.pyplot as plt


Params = ParamsStruct()
Params.VERBOSE = True  # Provides extra runtime messages. Overrides QUIET below
Params.QUIET = False  # Hides all log messages except warnings and errors
Params.printFmt = '[%(levelname)s] %(message)s'  # Format for printing log messages
Params.DEBUG = False  # Special use

Params.DO_PARALLEL = True  # Use multiprocessing module for parallel computation where applicable
Params.CALC_NEW =       True  # Recalculate profiles? If not, read data from disk and re-plot.
Params.CALC_NEW_REF =   True  # Recalculate reference phase curves?
Params.CALC_NEW_INDUC = True  # Calculate magnetic induction responses?
Params.CALC_SEISMIC =   True  # Calculate sound speeds?
Params.CALC_CONDUCT =   True  # Calculate electrical conductivity?
Params.SKIP_INNER =     False  # Whether to skip past everything but ocean calculations after MoI matching (for large induction studies)
Params.HOLD =           True  # Whether to overlay runs when possible
Params.REDUCED =        True  # Whether to limit number of ocean layers for faster computation of layered induction
Params.DISP_LAYERS =    False  # Whether to display layer depths and heat fluxes for user
Params.DISP_TABLES =    False  # Whether to print latex-formatted tables to Matlab command line
Params.DEPRECATED =     False  # Whether to allow deprecated code to run. Will often cause errors.

# Plot Settings
Params.SKIP_PLOTS = False  # Whether to skip creation of all plots
Params.PLOT_GRAVITY = True  # Whether to plot Gravity and Pressure
Params.PLOT_HYDROSPHERE = True  # Whether to plot Conductivity with Interior Properties (Hydrosphere)
Params.PLOT_TRADEOFF = True  # Whether to plot mantle properties tradeoff
Params.PLOT_WEDGE = True  # Whether to plot interior wedge diagram

# Magnetic induction calculation settings
Params.DO_EUR = True  # Whether to calculate induction responses for Europa
Params.DO_GAN = True  # Whether to calculate induction responses for Ganymede
Params.DO_CAL = True  # Whether to calculate induction responses for Callisto
Params.DO_ENC = True  # Whether to calculate induction responses for Enceladus
Params.DO_MIR = True  # Whether to calculate induction responses for Miranda
Params.DO_ARI = True  # Whether to calculate induction responses for Ariel
Params.DO_PER = True  # Convert frequency axes to periods
Params.DO_LEGEND = True  # Whether to force legends
Params.PLOT_FFT = True  # Whether to show plots of fourier space
Params.PLOT_CONTOURS = True  # Contours or surfaces
Params.PLOT_V2021S = True  # Mark the selected ocean/conductivity combos used in Vance et al. 2021

Params.interpMethod = 'nearest'  # Interpolation method. Options are nearest, linear, and cubic. Notably, the 'linear' method causes wiggles in magnetic contour plots in the Matlab version.
Params.npts_k = 50  # Resolution for conductivity values in ocean conductivity/thickness plots
Params.npts_D = 60  # Resolution for ocean thickness as for conductivity
Params.np_intp = 200  # Number of interpolation points to use for Eckhardt method induction calculations
Params.npts_w = 100  # Resolution in log frequency space for magnetic excitation spectra
Params.np_wfine = 1000  # Fine-spacing resolution for log frequency spectrum
Params.nIntL = 3  # Number of ocean layers to use when REDUCED = 1
#To be implemented- eventually need some ODE numerical solution parameters
#Params.opts_odeParams = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep', 2e3,'InitialStep',1e-2)
#Params.opts_odeLayers = odeset('RelTol',1e-8, 'AbsTol',1e-10,'MaxStep',10e3,'InitialStep',1e-2)

# General figure options
Params.defaultFontName = 'STIXGeneral'  # Default font variables--STIX is what is used in Icarus journal submissions
Params.defaultFontCode = 'stix'  # Code name for default font needed in some function calls
Params.backupFont = 'Times New Roman'  # Backup font that looks similar to STIX that most users are likely to have
Params.figFormat = 'pdf'
Params.xtn = '.' + Params.figFormat  # Figure file extension. Good options are .eps, .pdf, and .png
plt.rcParams['font.family'] = 'serif'  # Choose serif font for figures to best match math mode variables in body text
plt.rcParams['font.serif'] = Params.defaultFontName  # Set plots to use the default font
Params.PLOT_SIGS = False  # Make a plot of conductivity as a function of radius
Params.wlims = None  # Minimum and maximum to use for frequency spectrum plots (magnetic induction)
Params.LEGEND = False  # Whether to force legends to appear
Params.LegendPosition = None  # Where to place legends when forced
Params.yLim = None  # y axis limits of hydrosphere density in "Conductivity with interior properties" plot
Params.LineStyle = None  # Default line style to use on plots
Params.wRefLine_temporary = None  # Style of lines showing reference melting curves of hydrosphere density plot-should be done in config.py instead, delete this once implemented there
Params.wRef = None  # Salinities in ppt of reference melting curves
# Check if Latex executable is on the path so we can use backup options if Latex is not installed
if shutil.which('latex'):
    plt.rcParams['text.usetex'] = True  # Use Latex interpreter to render text on plots
    plt.rcParams['text.latex.preamble'] = f'\\usepackage{{{Params.defaultFontCode}}}'  # Load in font package in Latex
    Params.TEX_INSTALLED = True
else:
    print('A LaTeX installation was not found. Some plots may have fallback options in labels.')
    plt.rcParams['font.serif'] += ', ' + Params.backupFont  # Set plots to use the default font if installed, or a backup if not
    plt.rcParams['mathtext.fontset'] = Params.defaultFontCode
    Params.TEX_INSTALLED = False



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

# Linestyle options
#LS = linestyle
#LW = linewidth
# syn = synodic period (the time between, e.g. Jupiter showing the same longitude to Europa)
# orb = orbital period
# hrm = 2nd harmonic of the synodic period - the first excited harmonic of the synodic period magnetic oscillation
Params.LS_syn = '-'  # linestyle for the synodic period
Params.LS_orb = ':'  # linestyle for the orbital period
Params.LS_hrm = '-.'  # linestyle for the 2nd harmonic of the synodic period
Params.LS_Sw  = '-'  # linestyle for Saltwater
Params.LS_Mg = '--'  # linestyle for magnesium sulfate MgSO4
Params.LS_sp =  ':'
Params.LW_syn = 2  # linewidth for the synodic period
Params.LW_orb = 2  # linewidth for the orbital period
Params.LW_hrm = 2  # linewidth for the 2nd harmonic of the synodic period
Params.LW_sal = 3  # linewidth for higher salinity
Params.LW_dil = 1  # linewidth for dilute salinity
Params.LW_std = 2  # linewidth for standard salinity
Params.LW_sound = 1.5  # LineWidth for sound speed plots
Params.LW_seism = 1  # LineWidth for seismic plots (Attenuation)

# SPICE kernels
Params.spiceTLS = 'naif0012.tls'  # Leap-seconds kernel
Params.spicePCK = 'pck00010.tpc'  # Planetary Constants Kernel from SPICE in order to get body radii
Params.spiceJupiter = 'jup365.bsp'  # Generic kernel for Jupiter + Galilean moons
Params.spiceSaturn = 'sat427.bsp'  # Generic kernel for Saturn + large moons
Params.spiceUranus = 'ura111.bsp'  # Generic kernel for Uranus + large moons
Params.spiceNeptune = 'nep095.bsp'  # Generic kernel for Neptune + large moons

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


Params.FigSize.vsP = (3,3)
Params.FigSize.vsR = (3,3)
Params.FigSize.vperm = (3,3)
Params.FigSize.vgsks = (3,3)
Params.FigSize.vseis = (3,3)
Params.FigSize.vhydro = (6,5)
Params.FigSize.vgrav = (6,5)
Params.FigSize.vmant = (6,6)
Params.FigSize.vcore = (6,6)
Params.FigSize.vpvt4 = (3,3)
Params.FigSize.vpvt6 = (3,3)
Params.FigSize.vwedg = (3,3)
