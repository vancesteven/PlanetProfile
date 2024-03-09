"""
General runtime configuration parameters.
Overridden by any settings contained within PPBody.py files.
"""
import os
from PlanetProfile.Utilities.defineStructs import ParamsStruct, ExploreParamsStruct, Constants

configVersion = 15  # Integer number for config file version. Increment when new settings are added to the default config file.

def configAssign():
    Params = ParamsStruct()
    ExploreParams = ExploreParamsStruct()

    Params.VERBOSE =       False  # Provides extra runtime messages. Overrides QUIET below
    Params.QUIET =         False  # Hides all log messages except warnings and errors
    Params.QUIET_MOONMAG = True  # If True, sets MoonMag logging level to WARNING, otherwise uses the same as PlanetProfile.
    Params.QUIET_LBF = True  # If True, sets lbftd and mlbspline logging levels to ERROR, otherwise uses the same as PlanetProfile.
    Params.printFmt = '[%(levelname)s] %(message)s'  # Format for printing log messages
    # The below flags allow or prevents extrapolation of EOS functions beyond the definition grid.
    Params.EXTRAP_ICE = {'Ih':False, 'II':False, 'III':False, 'V':False, 'VI':False, 'Clath':False}
    Params.EXTRAP_OCEAN = False
    Params.EXTRAP_REF =   False  # Allow refprofile extrapolation separate from normal ocean
    Params.EXTRAP_SIL =   False
    Params.EXTRAP_Fe =    False
    Params.lookupInterpMethod = 'nearest'  # Interpolation method to use for EOS lookup tables. Options are 'nearest', 'linear', 'cubic'.
    Params.minPres_MPa = None  # Only applies to ice EOS! Applies a lower bound to how small the pressure step can be in loading the EOS. Avoids major slowdowns when chaining models for small and large bodies.
    Params.minTres_K = None  # Same as above. Set to None to allow default behavior for ice EOS resolution.

    Params.CALC_NEW =         True  # Recalculate profiles? If not, read data from disk and re-plot.
    Params.CALC_NEW_REF =     True  # Recalculate reference melting curve densities?
    Params.CALC_NEW_INDUCT =  True  # Recalculate magnetic induction responses?
    Params.CALC_NEW_ASYM =    False  # Recalculate asymmetric boundary plot(s)?
    Params.CALC_SEISMIC =     True  # Calculate sound speeds and elastic moduli?
    Params.CALC_CONDUCT =     True  # Calculate electrical conductivity?
    Params.CALC_VISCOSITY =   True  # Calculate viscosity for all layers as a post-processing step?
    Params.CALC_ASYM =        True  # Calculate induction with asymmetric shape?
    Params.RUN_ALL_PROFILES = False  # Whether to run all PPBody.py files for the named body and plot together
    Params.SPEC_FILE =        False  # Whether we are running a specific file or files
    Params.COMPARE =          False  # Whether to plot each new run against other runs from the same body
    Params.DO_PARALLEL =      True  # Whether to use multiprocessing module for parallel computation where applicable
    Params.threadLimit =      1000  # Upper limit to number of processors/threads for parallel computation
    Params.FORCE_EOS_RECALC = False  # Whether to reuse previously loaded EOS functions for multi-profile runs
    Params.SKIP_INNER =       False  # Whether to skip past everything but ocean calculations after MoI matching (for large induction studies)
    Params.NO_SAVEFILE =      False  # Whether to prevent printing run outputs to disk. Saves time and disk space for large induction studies.
    Params.DISP_LAYERS =      True  # Whether to display layer depths and heat fluxes for user
    Params.DISP_TABLE =       True  # Whether to print latex-formatted table
    Params.ALLOW_BROKEN_MODELS = False  # Whether to continue running models that don't match physical constraints (i.e. MoI), with many values set to nan. Currently only implemented for CONSTANT_INNER_DENSITY = True and only allows broken MoI matching. Broken Tb_K matching is also intended.
    Params.DEPRECATED =       False  # Whether to allow deprecated code to run. Will often cause errors.

    # Plot Settings
    Params.SKIP_PLOTS =       False  # Whether to skip creation of all plots
    Params.PLOT_GRAVITY =     True  # Whether to plot Gravity and Pressure
    Params.PLOT_HYDROSPHERE = True  # Whether to plot Conductivity with Interior Properties (Hydrosphere)
    Params.PLOT_REF =         True  # Whether to plot reference melting curve densities on hydrosphere plot
    Params.PLOT_SIGS =        True  # Whether to plot conductivities as a function of radius on hydrosphere plot if they have been calculated
    Params.PLOT_SOUNDS =      True  # Whether to plot sound speeds as a function of radius on hydrosphere plot if they have been calculated
    Params.PLOT_TRADEOFF =    True  # Whether to plot mantle properties tradeoff
    Params.PLOT_POROSITY =    True  # Whether to plot porosities in rock and/or ice for bodies that have it modeled
    Params.PLOT_SEISMIC =     True  # Whether to plot seismic quantities if they have been calculated
    Params.PLOT_VISCOSITY =   True  # Whether to plot viscosities
    Params.PLOT_WEDGE =       True  # Whether to plot interior wedge diagram
    Params.PLOT_HYDRO_PHASE = False  # Whether to plot phase diagram
    Params.PLOT_PVT_HYDRO =   False  # Whether to plot hydrosphere PT property plots
    Params.PLOT_PVT_INNER =   False  # Whether to plot silicate/core PT property plots
    Params.PLOT_BDIP =        False  # Whether to plot induced dipole surface strength in complex plane
    Params.PLOT_BSURF =       True  # Whether to plot induced field surface map
    Params.PLOT_ASYM =        False  # Whether to plot asymmetric boundary shape(s) when induced fields are calculated from them
    Params.PLOT_TRAJECS =     False  # Whether to plot spacecraft flyby trajectories for those involved in inversion
    Params.PLOT_BINVERSION =  False  # Whether to plot magnetic field vector components for flyby inversion
    Params.LEGEND =           True  # Whether to plot legends
    Params.TITLES =           True  # Whether to include a (sup)title on plots

    # Magnetic induction plot settings
    Params.DO_INDUCTOGRAM =          False  # Whether to evaluate and/or plot an inductogram for the body in question
    Params.INDUCTOGRAM_IN_PROGRESS = False  # Whether we are currently working on constructing an inductogram
    Params.COMBINE_BCOMPS =          False  # Whether to plot Bx, By, Bz with phase all in one plot, or separate for each comp -- same for Bdip components
    Params.PLOT_MAG_SPECTRUM_COMBO = False  # Whether to show plots of fourier space for magnetic induction, including induced spectrum
    Params.PLOT_MAG_SPECTRUM =       False  # Whether to show plot of fourier space of excitation moments for magnetic induction
    Params.tRangeCA_s =              120  # Range in seconds relative to named closest approach UTC datetime to search for the actual CA as identified by querying SPICE kernels

    # Parameter exploration plot settings
    Params.DO_EXPLOREOGRAM = False  # Whether to evaluate and/or plot an exploreogram for the body in question
    Params.SKIP_INDUCTION = False  # Whether to skip past induction calculations. Primarily intended to avoid duplicate calculations in exploreOgrams
    # Options for x/y variables: "xFeS", "rhoSilInput_kgm3", "wOcean_ppt", "Tb_K", "ionosTop_km", "sigmaIonos_Sm",
    # "silPhi_frac", "silPclosure_MPa", "icePhi_frac", "icePclosure_MPa", "Htidal_Wm3", "Qrad_Wkg", "qSurf_Wm2" (Do.NO_H2O only)
    ExploreParams.xName = 'wOcean_ppt'  # x variable over which to iterate for exploreograms. Options are as above.
    ExploreParams.yName = 'Tb_K'  # y variable over which to iterate for exploreograms. Options are as above.
    # Options for z variables: "CMR2mean", "D_km", "dzIceI_km", "dzIceI_km", "dzClath_km", "dzIceIII_km", "dzIceIIIund_km",
    # "dzIceV_km", "dzIceVund_km", "dzIceVI_km", "dzWetHPs_km", "eLid_km", "phiSeafloor_frac", "Rcore_km", "rhoSilMean_kgm3",
    # "sigmaMean_Sm", "silPhiCalc_frac", "zb_km", "zSeafloor_km", "qSurf_Wm2" (only if Do.NO_H2O is False).
    # New options must be added to ExplorationStruct attributes in Main (assign+save+reload) and in defineStructs, and
    # FigLbls.exploreDescrip, .<var>Label, and .axisLabels in defineStructs.
    ExploreParams.zName = ['CMR2calc', 'silPhiCalc_frac', 'phiSeafloor_frac', 'D_km', 'zb_km', 'dzWetHPs_km', 'rhoSilMean_kgm3', 'Pseafloor_MPa', 'zSeafloor_km', 'sigmaMean_Sm']  # heatmap/colorbar/z variable to plot for exploreograms. Options are as above; accepts a list.
    ExploreParams.xRange = [10.0, 100.0]  # [min, max] values for the x variable above
    ExploreParams.yRange = [249.0, 272.5]  # Same as above for y variable
    ExploreParams.nx = 30  # Number of points to use in linspace with above x range
    ExploreParams.ny = 24  # Same as above for y

    # Reference profile settings
    # Salinities of reference melting curves in ppt
    Params.wRef_ppt = {'none':[0], 'PureH2O':[0],
                       'Seawater':[0, 0.5*Constants.stdSeawater_ppt, Constants.stdSeawater_ppt, 1.5*Constants.stdSeawater_ppt],
                       'MgSO4':[0, 33.3, 66.7, 100],
                       'NH3':[0, 5, 10, 15],
                       'NaCl':[0, 17.5, 35]}
    Params.nRefRho = 50  # Number of values for plotting reference density curves (sets resolution)
    Params.PrefOverride_MPa = None  # Pressure setting to force refprofile recalc to go to a specific value instead of automatically using the first hydrosphere max

    # SPICE kernels to use
    Params.spiceDir = 'SPICE'
    Params.spiceTLS = 'naif0012.tls'  # Leap-seconds kernel
    Params.spicePCK = 'pck00010.tpc'  # Planetary Constants Kernel from SPICE in order to get body radii
    Params.spiceFK  = [  # Frames kernels used for converting between coordinate systems
        'custom_frames_v01.tf',  # Defines MOON_PHI_O (for Galilean moons), PSO, PSM, and PDSZ frames
    ]
    Params.spiceBSP = {
        'Earth': 'de430.bsp',  # Generic kernel for solar system planets
        'Jupiter': 'jup365.bsp',  # Generic kernel for Jupiter + Galilean moons
        'Saturn': 'sat441.bsp',  # Generic kernel for Saturn + large moons
        'Uranus': 'ura111.bsp',  # Generic kernel for Uranus + large moons
        'Neptune': 'nep097.bsp'  # Generic kernel for Neptune + Triton
    }

    return Params, ExploreParams
