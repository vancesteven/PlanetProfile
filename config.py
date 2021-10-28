"""
General runtime configuration parameters.
Overridden by any settings contained within PPBody.py files.
"""

from Utilities.dataStructs import ParamsStruct

Params = ParamsStruct()
Params.VERBOSE = True
if Params.VERBOSE: print('Printing verbose runtime messages. Toggle in config.py.')

Params.DO_PARALLEL = True
Params.CALC_NEW =       1  # Recalculate profiles? If not, read data from disk and re-plot.
Params.CALC_NEW_REF =   1  # Recalculate reference phase curves?
Params.CALC_NEW_SOUND = 1  # Recalculate sound speeds?
Params.CALC_NEW_INDUC = 1  # Recalculate magnetic induction responses?
Params.SKIP_PROFILES =  0  # Whether to skip past all PlanetProfile.m plotting and calculations
Params.NO_PLOTS =       0  # Suppress plot creation? If yes, just save profiles to disk.
Params.HOLD =           1  # Whether to overlay runs when possible
Params.CONDUCT =        1  # Calculate electrical conductivity
Params.REDUCED =        1  # Whether to limit number of ocean layers for faster computation of layered induction
Params.DISP_LAYERS =    0  # Whether to display layer depths and heat fluxes for user
Params.DISP_TABLES =    0  # Whether to print latex-formatted tables to Matlab command line
Params.DEPRECATED =     0  # Whether to allow deprecated code to run. Will often cause errors.

# Magnetic induction calculation settings
Params.DO_EUR = 1  # Whether to calculate induction responses for Europa
Params.DO_GAN = 1  # Ganymede
Params.DO_CAL = 1  # And so on...
Params.DO_ENC = 1
Params.DO_MIR = 1
Params.DO_ARI = 1
Params.DO_PER = 1  # Convert frequency axes to periods
Params.DO_LEGEND = 1  # Whether to force legends
Params.PLOT_FFT = 1  # Whether to show plots of fourier space
Params.PLOT_CONTOURS = 1  # Contours or surfaces
Params.PLOT_V2020S = 1  # Mark the selected ocean/conductivity combos used in Vance et al. 2020
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
Params.dft_font = 'stix'
Params.dft_math = 'stix'
#Params.interpreter = 'tex'
#Params.fig_fmt = '-depsc'
Params.xtn = '.eps'

# Color selection
Params.cmap = 'viridis'
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
Params.LS_syn = '-'
Params.LS_orb = ':'
Params.LS_hrm = '-.'
Params.LS_Sw  = '-'
Params.LS_Mg = '--'
Params.LS_sp =  ':'
Params.LW_syn = 2
Params.LW_orb = 2
Params.LW_hrm = 2
Params.LW_sal = 3
Params.LW_dil = 1
Params.LW_std = 2
Params.LW_sound = 1.5
Params.LW_seism = 1

# SPICE kernels
Params.spicePCK = 'pck00010.tpc'