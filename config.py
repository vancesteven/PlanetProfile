"""
General runtime configuration parameters.
Overridden by any settings contained within PPBody.py files.
"""

from Utilities.dataStructs import ParamsStruct
import matplotlib as mpl
import matplotlib.pyplot as plt


Params = ParamsStruct()
Params.VERBOSE = True #sets amount of printed runtime messages
if Params.VERBOSE: print('Printing verbose runtime messages. Toggle in config.py.')

Params.DO_PARALLEL = True # Use multiprocessing module for parallel computation where applicable
Params.CALC_NEW =       True  # Recalculate profiles? If not, read data from disk and re-plot.
Params.CALC_NEW_REF =   True  # Recalculate reference phase curves?
Params.CALC_NEW_SOUND = True  # Recalculate sound speeds?
Params.CALC_NEW_INDUC = True  # Recalculate magnetic induction responses?
Params.SKIP_PROFILES =  False  # Whether to skip past all PlanetProfile.m plotting and calculations
Params.NO_PLOTS =       False  # Suppress plot creation? If yes, just save profiles to disk.
Params.HOLD =           True  # Whether to overlay runs when possible
Params.CONDUCT =        True  # Calculate electrical conductivity
Params.REDUCED =        True  # Whether to limit number of ocean layers for faster computation of layered induction
Params.DISP_LAYERS =    False  # Whether to display layer depths and heat fluxes for user
Params.DISP_TABLES =    False  # Whether to print latex-formatted tables to Matlab command line
Params.DEPRECATED =     False  # Whether to allow deprecated code to run. Will often cause errors.

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

Params.intMethod = 'makima'  # Interpolation method. Certain ones can cause wiggles, notably 'linear'.
Params.npts_k = 50  # Resolution for conductivity values in ocean conductivity/thickness plots
Params.npts_D = 60  # Resolution for ocean thickness as for conductivity
Params.np_intp = 200  # Number of interpolation points to use for Eckhardt method induction calculations
Params.npts_w = 100  # Resolution in log frequency space for magnetic excitation spectra
Params.np_wfine = 1000  # Fine-spacing resolution for log frequency spectrum
Params.nIntL = 3  # Number of ocean layers to use when REDUCED = 1
#Params.opts_odeParams = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep', 2e3,'InitialStep',1e-2)
#Params.opts_odeLayers = odeset('RelTol',1e-8, 'AbsTol',1e-10,'MaxStep',10e3,'InitialStep',1e-2)

# General figure options
Params.dft_font = 'stix'  # default font variables- STIX is what is used in Icarus journal submissions
Params.dft_math = 'stix'  # default math font variables- STIX is what is used in Icarus journal submissions
Params.xtn = '.eps'  # figure file extension. Good options are .eps, .pdf, and .png


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
Params.spicePCK = 'pck00010.tpc'  # Planetary Constants Kernel from SPICE in order to get body radii
