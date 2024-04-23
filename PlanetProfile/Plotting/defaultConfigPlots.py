""" Default figure settings """
import numpy as np
import spiceypy as spice
from PlanetProfile.Utilities.defineStructs import ColorStruct, StyleStruct, \
    FigLblStruct, FigSizeStruct, FigMiscStruct

configPlotsVersion = 15  # Integer number for config file version. Increment when new settings are added to the default config file.

def plotAssign():
    Color = ColorStruct()
    Style = StyleStruct()
    FigLbl = FigLblStruct()
    FigSize = FigSizeStruct()
    FigMisc = FigMiscStruct()

    """ Figure color options """
    Color.cycler = None  # Color cycler to use for multi-line plots when colors are not important to specify individually. None uses a custom cycler. 'default' uses the default, which is not well-adapted for colorblind viewers.
    Color.Induction = {'synodic': 'blue', 'orbital': 'purple', 'true anomaly': 'green', 'synodic 2nd': 'goldenrod'}  # Colors for inductOgram plots
    Color.ref = 'gray'
    Color.geothermHydro = 'yellow'
    Color.geothermInner = 'white'
    Color.BdipInset = 'black'  # Color for inset box of surface induced dipole strength plots

    # Wedge diagram color options
    Color.none = '#FFFFFF00'
    Color.wedgeBd = 'black'
    Color.wedgeMarkRadii = 'black'
    Color.ionoCmapName = 'RdPu'
    Color.ionoTop = 0.0
    Color.ionoBot = 0.2
    Color.ionoN = 10
    Color.iceIcond = 'xkcd:ice blue'
    Color.iceIconv = 'xkcd:robin\'s egg blue'
    Color.iceII = 'xkcd:pale sky blue'
    Color.iceIII = 'xkcd:carolina blue'
    Color.iceV = 'xkcd:light periwinkle'
    Color.iceVI = '#91d1d4'
    Color.clathCond = 'xkcd:pastel blue'
    Color.clathConv = 'xkcd:dusty blue'
    Color.oceanCmapName = 'ocean_r'
    Color.oceanTop = 0.2  # Fraction of ocean colormap to start at
    Color.oceanBot = 0.45  # Fraction of ocean colormap to end at
    Color.oceanN = 20
    Color.silPorousCmapName = 'BrBG'
    Color.silPorousTop = 0.0
    Color.silPorousBot = 0.05
    Color.silPorousN = 5
    Color.silCondCmapName = 'BrBG'
    Color.silCondTop = 0.05 - 0.001
    Color.silCondBot = 0.05
    Color.silCondN = 1
    Color.silConvCmapName = 'gist_heat'
    Color.silConvTop = 0.5
    Color.silConvBot = 0.9
    Color.silConvN = 10
    Color.FeS = 'xkcd:puke'
    Color.Fe = '#2d3639'
    # Alternative color options for silicates
    Color.PALE_SILICATES = False  # Whether to use a lighter color scheme for silicate layers, or a more "orangey" saturated one
    Color.paleSilPorousCmapName = 'terrain'
    Color.paleSilPorousTop = 0.75
    Color.paleSilPorousBot = 0.78
    Color.paleSilCondCmapName = 'terrain'
    Color.paleSilCondTop = 0.78 - 0.001
    Color.paleSilCondBot = 0.78

    # PvThydro properties colormaps
    Color.PvThydroCmapName = 'Blues'
    Color.PvThydroHi = 1.0
    Color.PvThydroLo = 0.2
    Color.negPvThydroCmapName = 'Reds'
    Color.negPvThydroHi = 0.5
    Color.negPvThydroLo = 0.0

    # PvTPerpleX properties colormaps
    Color.PvTsilCmapName = 'copper'
    Color.PvTsilHi = 1.0
    Color.PvTsilLo = 0.0
    Color.PvTcoreCmapName = 'inferno'
    Color.PvTcoreHi = 1.0
    Color.PvTcoreLo = 0.0

    # Colormaps for inductogram phase space plots, hydrosphere plots, etc
    Color.cmapName = {
        'none': 'copper',
        'PureH2O': 'cividis',
        'Seawater': 'cool',
        'MgSO4': 'winter',
        'NH3': 'spring',
        'NaCl': 'summer',
        'Ice': 'coolwarm_r',
        'BmapPos': 'afmhot',
        'BmapNeg': 'afmhot_r',
        'BmapDiv': 'seismic',
        'asymDev': 'PuBu_r',
        'default': 'plasma'
    }
    # Select only a subset of the available colormap, if we choose to
    Color.cmapBounds = {
        'none': [0.0, 1.0],
        'PureH2O': [0.0, 0.8],
        'Seawater': [0.0, 1.0],
        'MgSO4': [0.0, 1.0],
        'NH3': [0.0, 1.0],
        'NaCl': [0.0, 1.0],
        'Ice': [0.2, 0.8],
        'BmapPos': [0.2, 1.0],
        'BmapNeg': [0.2, 1.0],
        'BmapDiv': [0.0, 1.0],
        'asymDev': [0.0, 1.0],
        'default': [0.0, 1.0]
    }
    # Set temperature bounds to use for colormap normalization
    Color.Tbounds_K = [255.0, 273.0]

    # Set upper bounds for max concentrations to use for darkening
    Color.saturation = {
        'none': 1.0,
        'PureH2O': 1.0,
        'Seawater': 304.0,
        'MgSO4': 282.0,
        'NH3': 100.0,
        'NaCl': 304.0
    }
    # Saturation & color brightness ("value" in HSV) values for salinity/conductivity axis bounds
    Color.fresh = [0.5, 1.0]
    Color.salty = [1.0, 0.5]

    # Assign colormaps to use the settings above
    Color.SetCmaps()

    # Color options for Fourier spectrum plots
    Color.BeiFT = {
        'x': 'blue',
        'y': 'black',
        'z': 'green'
    }
    Color.Ae1FT = 'purple'
    Color.TexcFT = 'red'

    # Color options for trajectory and CA plots
    Color.bodySurface = 'xkcd:grey'
    Color.CAdot = 'black'
    Color.thresh = 'blue'
    Color.CAline = 'black'
    Color.MAGdata = 'xkcd:pastel red'  # MAG data in trajectory plots
    Color.BcompsModelNet = 'xkcd:true green'  # Net magnetic field from models in trajectory plots
    Color.BcompsModelExc = 'xkcd:mustard yellow'  # Excitation field
    Color.BcompsModelInd = 'xkcd:royal blue'  # Induced field
    Color.BcompsModelPls = 'xkcd:dark magenta'  # Fields from plasma contributions


    """ Figure style options """
    Style.GRIDS = False  # Whether to plot grids
    Style.LS = {'none': None, 'PureH2O': ':', 'Seawater': '-', 'MgSO4': '--', 'NH3': '-.', 'NaCl': '-'}  # LineStyle options for hydrosphere plots
    Style.LWlims = [0.5, 2]  # Bounds of linewidths to use for salinity mapping
    Style.MW_hydro = 3  # Marker size for hydrosphere plot endpoint, in the factor by which to multiply the linewidth.
    Style.MS_hydro = 'o'  # Marker style for hydrosphere plot endpoint
    Style.LW_std = 1.5  # Standard linewidth to use when not mapping as above
    Style.LW_sound = 1  # Linewidth for sound speeds on hydrosphere plot
    Style.LW_geotherm = 1  # Linewidth for geotherm on PT plots
    Style.LS_geotherm = '-'  # Linestyle for geotherm on PT plots
    Style.LW_seis = 1  # Linewidth for seismic plots
    Style.LS_seis = {'KS': '-', 'GS': '--', 'VP': '-', 'VS': '--', 'QS': '-', 'P': '-', 'T': '--', 'rho': '-.'}
    Style.LS_ref = {'none': None, 'PureH2O': ':', 'Seawater': '-', 'MgSO4': '--', 'NH3': '-.', 'NaCl': '-'}  # Style for reference profiles
    Style.LW_ref = 0.75  # Linewidth for reference profiles
    Style.LS_Induction = {'synodic': '-', 'orbital': ':', 'true anomaly': ':', 'synodic 2nd': '--'}  # Style for inductOgram plots
    Style.LW_Induction = {'synodic': 1.5, 'orbital': 1.5, 'true anomaly': 1.5, 'synodic 2nd': 1.5}  # Widths for inductOgram plots
    Style.MW_Induction = 2  # Marker size to use for induction scatter plots
    Style.MS_Induction = 'o'  # Marker style for induction scatter plots

    # Wedge diagrams
    Style.wedgeAngle_deg = 25  # Angular size of wedge diagrams in degrees
    Style.LW_wedge = 0.125  # Linewidth in pt for minor boundaries in wedge diagrams
    Style.LW_wedgeMajor = 0.375  # Linewidth in pt for major layer boundaries in wedge diagrams
    Style.TS_ticks = 12  # Text size in pt for tick marks on radius scale
    Style.TS_desc = 14  # Text size in pt for model description and label
    Style.TS_super = 16  # Text size in pt for overall ("suptitle") label with multiple wedges
    Style.LS_markRadii = '--'  # Linestyle for radii mark line when toggled on
    Style.LW_markRadii = 0.375  # Linewidth for radii mark line when toggled on

    # Complex dipole plots
    Style.MW_dip = {'synodic': 5, 'orbital': 6, 'true anomaly': 6, 'synodic 2nd': 4.5}  # Marker size for each period in complex dipole plots
    Style.MS_dip = {'synodic': '*', 'orbital': 'o', 'true anomaly': 'P', 'synodic 2nd': 'h'}  # Marker style for each period in complex dipole plots
    Style.MAlims = [0, 1]  # Alpha channel (opacity) limits for markers
    Style.LS_BdipInset = '-'  # Linestyle for inset box
    Style.LW_BdipInset = 0.5  # Linewidth for inset box

    # Fourier spectrum plots
    Style.LS_FT = '-'  # Linestyle of Fourier spectrum plots
    Style.LW_FT = 0.75  # Linewidth for Ae1, Bx, By, Bz in Fourier spectrum plots
    Style.LS_TexcFT = '-'  # Linestyle for optional lines marking dominant excitations in Ae1 plot
    Style.LW_TexcFT = 0.5  # Linewidth for above

    # Trajectory and CA plots
    Style.LS_thresh = '-'  # Linestyle of MAG precision floor line
    Style.LW_thresh = 0.75  # Linewidth of MAG precision floor line
    Style.MS_CA = 'o'  # Marker style for closest approach points
    Style.MW_CA = 2  # Marker size for closest approach dots
    Style.LS_CA = '-'  # Linestyle for closest approach line
    Style.LW_CA = 0.5  # Linewidth for closest approach line
    Style.LS_MAGdata = '-'  # Linestyle of MAG data in trajectory plots
    Style.LW_MAGdata = 1  # Linewidth of MAG data in trajectory plots
    Style.LS_modelNet = '-'  # Linestyle of net model field in trajectory plots
    Style.LW_modelNet = 1  # Linewidth of net model field in trajectory plots
    Style.LS_modelExc = '--'  # Linestyle of excitation field
    Style.LW_modelExc = 1  # Linewidth of excitation field
    Style.LS_modelInd = '--'  # Linestyle of induced field
    Style.LW_modelInd = 1  # Linewidth of induced field
    Style.LS_modelPls = ':'  # Linestyle of fields from plasma contributions
    Style.LW_modelPls = 1  # Linewidth of fields from plasma contributions
    Style.LS_SCtrajec = {  # Linestyle of spacecraft trajectories
        'Cassini': '-',
        'Clipper': '-',
        'Galileo': '-',
        'Juno': '-',
        'JUICE': '-'
    }
    Style.LW_SCtrajec = 1  # Linewidth of spacecraft trajectories
    Style.MS_exit = 'o'  # Marker style for trajectory exit points
    Style.MW_exit = 10  # Marker size for trajectory exit dots


    """ Figure labels """
    # Note: Specific labels are set in PlanetProfile.Utilities.defineStructs.
    FigLbl.NEGATIVE_UNIT_POWERS = True  # Whether to use negative powers for units in latex tables, or instead a backslash.
    FigLbl.BODYNAME_IN_LABEL = True  # Whether to include the bodyname in legend labels for comparison plots including multiple bodies
    FigLbl.NAN_FOR_EMPTY = False  # Whether to use nan (or -) for empty layer parameters that were not calculated or not present.
    FigLbl.PFULL_IN_GPa = True  # Whether to plot P in GPa (or MPa) for full-body plots
    FigLbl.PHYDRO_IN_bar = False  # Whether to print P in bar (or MPa) for hydrosphere plots
    FigLbl.w_IN_WTPCT = False  # Whether to print salinities in wt% (or g/kg) in tables
    FigLbl.T_IN_C = False  # Whether to print T in deg C (or K) in plots
    FigLbl.x_IN_MOLPCT = True  # Whether to print silicate/core mass fractions in mol% (or fractional) in tables and plots
    FigLbl.qSURF_IN_mW = True  # Whether to print qSurf in mW/m^2 (or W/m^2)
    FigLbl.phi_IN_VOLPCT = False  # Whether to print porosity (phi) in vol% (or unitless volume fraction)
    FigLbl.PVT_CBAR_LABELS = False  # Whether to add short labels identifying silicate/core colorbars in PvT properties plots
    FigLbl.tCA_RELATIVE = True  # Whether to display trajectory x axes in time relative to closest approach or absolute times
    FigLbl.tCArelUnits = 'min'  # Units to use for times relative to closest approach in trajectory plots. Options are 'h', 'min', 's'.
    FigLbl.CAoffset = [0, 0.05]  # Offset for text of CA label from top-middle of marker lines. x units are axis units, y units are fractional of the line full height.
    FigLbl.AXES_INFO = True  # Whether to add explanatory info for IAU axis directions
    FigLbl.sciLimits = (-2, 4)  # Powers of 10 to use as limits on axis labels, e.g. [-2, 4] means anything < 0.01 or >= 10000 will use scientific notation.
    FigLbl.SetUnits()  # Make use of above toggles and assign labels


    """ Figure sizes """
    FigSize.vpore = (6, 6)
    FigSize.vperm = (6, 6)
    FigSize.vseis = (6, 6)
    FigSize.vhydro = (9, 5)
    FigSize.vgrav = (6, 5)
    FigSize.vmant = (6, 6)
    FigSize.vcore = (6, 6)
    FigSize.vvisc = (6, 6)
    FigSize.vpvt = (12, 6)
    FigSize.vwedg = (4.5, 4.5)
    FigSize.vphase = (5, 6)
    FigSize.explore = (6, 4)
    FigSize.phaseSpaceSolo = (6, 4)
    FigSize.phaseSpaceCombo = (9, 4)
    FigSize.induct = (8, 4)
    FigSize.inductCombo = (8, 8)
    FigSize.Bdip = (5, 3)
    FigSize.BdipCombo = (6, 9)
    FigSize.BdipSolo = (2.5, 3)
    FigSize.BdipSoloCombo = (3, 9)
    FigSize.MagFT = (6, 10)
    FigSize.MagFTexc = (5, 5)
    FigSize.MagSurf = (8, 5)
    FigSize.MagSurfCombo = (16, 5)
    FigSize.MagCA = (4, 4)
    FigSize.BtrajecCombo = (6, 9)
    FigSize.SCtrajecCombo = (6, 4)
    FigSize.SCtrajec3D = (6, 6)
    FigSize.asym = (8, 5)
    FigSize.apsidal = (6, 6)


    """ Miscellaneous figure options """
    # General figure options
    FigMisc.figFormat = 'pdf'
    FigMisc.dpi = 300  # Resolution in dots per inch for raster images (.png). Ignored for vector images (.pdf, .eps)
    FigMisc.xtn = '.' + FigMisc.figFormat  # Figure file extension. Good options are .eps, .pdf, and .png
    FigMisc.defaultFontName = 'STIXGeneral'  # Default font variables--STIX is what is used in Icarus journal submissions
    FigMisc.defaultFontCode = 'stix'  # Code name for default font needed in some function calls
    FigMisc.backupFont = 'Times New Roman'  # Backup font that looks similar to STIX that most users are likely to have
    FigMisc.FORCE_0_EDGES = True  # Sets the edge of plots with 0 radius, depth, pressure, etc. to be the edge of the axes, instead of including white space which is the default.

    # Hydrosphere plots
    FigMisc.LOG_SIG = False  # Whether to print conductivity plot on a log scale
    FigMisc.COMMON_ZMAX_SIG = False  # Whether to force conductivity plot to have the same maximum depth as other hydrosphere plots, or to let the bottom axis set automatically to zoom in on the ocean. Only has an effect for undersea HP ices.
    FigMisc.SHOW_ICE_CONDUCT = False  # Whether to force conductivity plot to include (usually arbitrarily small) conductivities in ice phases.
    FigMisc.SCALE_HYDRO_LW = True  # Whether to adjust thickness of lines on hydrosphere plot according to relative salinity
    FigMisc.MANUAL_HYDRO_COLORS = True  # Whether to set color of lines in hydrosphere according to melting temperature
    FigMisc.RELATIVE_Tb_K = True  # Whether to set colormap of lines based on relative comparison (or fixed settings in ColorStruct)
    FigMisc.lowSigCutoff_Sm = 1e-3  # Cutoff conductivity below which profiles will be excluded. Setting to None includes all profiles
    FigMisc.PminHydro_MPa = None  # Minimum pressure to use for hydrosphere and phase diagram PT plots in MPa. Set to None to use min of geotherm.
    FigMisc.TminHydro = 250  # Minimum temperature to display on hydrosphere plots
    FigMisc.PmaxHydro_MPa = None  # When set, maximum pressure to use for hydrosphere and phase diagram PT plots in MPa. Set to None to use max of geotherm.
    FigMisc.TmaxHydro_K = None  # When set, maximum temperature to use for hydrosphere and phase diagram PT plots in K. Set to None to use max of geotherm.
    FigMisc.PHASE_LABELS = True  # Whether to print phase labels on density plots
    FigLbl.TS_hydroLabels = 18  # Font size for hydrosphere phase labels in pt

    # Wedge diagrams
    FigMisc.IONOSPHERE_IN_WEDGE = False  # Whether to include specified ionosphere in wedge diagram
    FigMisc.WEDGE_ICE_TICKS = False  # Whether to print ticks for ice shell, which usually overlap with the body outer radius
    FigMisc.DRAW_IONOS_BOUND = False  # Whether to draw a boundary line around the ionosphere
    FigMisc.DRAW_CONVECTION_BOUND = False  # Whether to draw a boundary line between convecting and conducting regions
    FigMisc.DRAW_POROUS_BOUND = False  # Whether to draw a boundary line between porous and non-porous materials
    FigMisc.DRAW_FeS_BOUND = True  # Whether to draw a boundary line between Fe and FeS in the core
    FigMisc.minzbRratio_frac = 0.05  # Fraction of total body radius for ice shell and ocean thickness, above which ice shell ticks will automatically switch on (overrides WEDGE_ICE_TICKS)
    FigMisc.MARK_RADII = True  # Whether to add a marker line from radius labels to layer arc
    FigMisc.LABEL_RADII = False  # Whether to add a label to radius km numbers

    # Hydrosphere PT diagrams
    FigMisc.PT_RASTER = True  # Whether to rasterize gridded information in PT plots and phase diagrams. Dramatically speeds up figure creation time and reduces file size, but renders gridded data grainy upon zoom-in.
    FigMisc.nTphase = 240  # Number of temperature points to evaluate/plot for hydrosphere phase diagram. If None, use geotherm (Planet.P_MPa) for the phase diagram.
    FigMisc.nPphase = 360  # Number of pressure points to evaluate/plot for hydrosphere phase diagram. If None, use geotherm (Planet.T_K) for the phase diagram.
    FigMisc.nThydro = 160  # Number of temperature points to evaluate/plot for PT property plots
    FigMisc.nPhydro = 200  # Number of pressure points to evaluate/plot for PT property plots
    FigMisc.PminHydro_MPa = 0.1  # Minimum pressure to use for hydrosphere and phase diagram PT plots in MPa. Set to None to use min of geotherm.
    FigMisc.TminHydro_K = 220  # Minimum temperature to use for hydrosphere and phase diagram PT plots in K. Set to None to use min of geotherm.
    FigLbl.hydroPhaseSize = 14  # Font size of label for phase in phase diagram

    # Silicate/core PT diagrams
    FigMisc.nTgeo = 80  # Number of temperature points to evaluate/plot for PT property plots
    FigMisc.nPgeo = 100  # Number of pressure points to evaluate/plot for PT property plots
    FigMisc.nPgeoCore = 40  # Subset of nPgeo to use for core, if present
    FigMisc.PVT_INCLUDE_CORE = True  # Whether to include core as well as silicates in PT properties diagrams

    # Induced field surface strength plots
    FigMisc.FIXED_COLORBAR = True  # Whether to maintain the same colorbar scale for each evaluation time in B surface plots.
    FigMisc.FIXED_ALL_COMPS = False  # Whether to apply the above across x, y, and z or just within each component
    FigMisc.MAG_CBAR_SEPARATE = True  # Whether to use an independent colormap/norm from the above settings for the magnitude
    FigMisc.BdipZoomMult = 1.05  # Extra space to include around zoomed-in part, in fraction of largest value.
    FigMisc.SHOW_INSET = True  # Whether to show the inset box for the zoom-in plot, when applicable
    FigMisc.rMagEval_Rp = 1.0  # Fraction of body radius to use for surface over which PlotMagSurface is evaluated
    FigMisc.tMagLbl = None  # List of strings to use to describe the times listed in tMagEval_s
    FigMisc.tMagEval_s = spice.str2et('2000-01-01T11:58:55.816')  # Time in seconds past J2000 at which to evaluate magnetic field surface maps. spiceypy.str2et will convert a datetime string to the correct format. Accepts an array.
    FigMisc.LARGE_ADJUST = True  # Whether to make certain labels better for cramped spaces, including removing colorbars (True is more pared down, and overrides nLon/LatTicks below.)
    FigMisc.BASYM_WITH_SYM = True  # Whether to plot Basym plot and Bsym plot on the same figure
    FigMisc.vCompMagSurf = 'mag'  # Component to use for induced field surface strength plots. Options are ['x', 'y', 'z', 'mag', 'all'].
    FigMisc.mapTitleSize = None  # Override for title font size. Defaults for normal and LARGE_ADJUST are in Utilities.defineStructs.
    FigMisc.nPPGCmapRes = 180  # Number of points per great circle to use for map angular resolution
    FigMisc.DO_360 = False  # Whether to range longitudes from 0 to 360 or from -180 to +180
    FigMisc.nLatTicks = 7  # Number of ticks to mark on latitude axis
    FigMisc.nLonTicks = 9  # Number of ticks to mark on longitude axis
    FigMisc.nMagContours = 9  # Number of contour intervals to mark on magnetic plots
    FigMisc.nAsymContours = 7  # Number of contour intervals to mark on asymmetry maps
    FigMisc.latlonSize = 14  # Font size for lat/lon labels
    FigMisc.cLabelSize = 12  # Font size for contour labels
    FigMisc.cLabelPad = 1  # Padding in pt to use for contour labels
    FigMisc.vminMagSurf_nT = None  # Minimum value for colormap in magnetic field surface plots (None uses data min/max). Overwritten by FIXED_COLORBAR.
    FigMisc.vmaxMagSurf_nT = None  # Minimum value for colormap in magnetic field surface plots (None uses data min/max). Overwritten by FIXED_COLORBAR.
    FigMisc.vminMagSurfDiff_nT = None  # Minimum value for difference colormap in magnetic field surface plots (None uses data min/max). Overwritten by FIXED_COLORBAR.
    FigMisc.vmaxMagSurfDiff_nT = None  # Minimum value for difference colormap in magnetic field surface plots (None uses data min/max). Overwritten by FIXED_COLORBAR.
    FigMisc.vminMagSurfComp_nT = None  # Minimum value for model comparison colormap in magnetic field surface plots (None uses data min/max). Overwritten by FIXED_COLORBAR.
    FigMisc.vmaxMagSurfComp_nT = None  # Minimum value for model comparison colormap in magnetic field surface plots (None uses data min/max). Overwritten by FIXED_COLORBAR.

    # Magnetic field trajectory and CA plots
    FigMisc.trajLims = None  # Distance in body radii at which to cut off trajectory plots
    FigMisc.MARK_CA_B = True  # Whether to mark closest approach on B trajectory plots with a line and text
    FigMisc.MARK_CA_POS = True  # Whether to mark closest approach on trajectory plots with a line and text
    FigMisc.EXIT_ARROWS = True  # Whether to use arrows or other marker types to indicate exit points (selected in Style.MS_exit)
    FigMisc.CAlblSize = 12  # Size of text labels on CA points
    FigMisc.SHOW_MAG_THRESH = True  # Whether to show a line indicating the precision floor of a magnetometer
    FigMisc.thresh_nT = 0.0488  # Precision floor in nT for magnetometer to plot
    FigMisc.threshCenter = 100  # x coordinate to place the MAG floor label
    FigMisc.hCAmax_km = 500  # Maximum altitude to show on CA plot
    FigMisc.SHOW_EXCITATION = False  # Whether to show the background field in trajectory plots, separately from the net model field
    FigMisc.SHOW_INDUCED = False  # Whether to show induced field in trajectory plots, separately from the net model field
    FigMisc.SHOW_PLASMA = False  # Whether to show plasma contributions in trajectory plots, separately from the net model field

    # Inductogram phase space plots
    FigMisc.DARKEN_SALINITIES = False  # Whether to match hues to the colorbar, but darken points based on salinity, or to just use the colorbar colors.
    FigMisc.NORMALIZED_SALINITIES = False  # Whether to normalize salinities to absolute concentrations relative to the saturation limit for each salt
    FigMisc.NORMALIZED_TEMPERATURES = False  # Whether to normalize ocean mean temperatures to specified maxima and minima for the colormap
    # Inductograms
    FigMisc.MARK_INDUCT_BOUNDS = True  # Whether to draw a border around the models on sigma/D plot when combined
    FigMisc.PLOT_V2021 = False  # Whether to mark the selected ocean/conductivity combos used in Vance et al. 2021
    # Excitation spectra
    FigMisc.MAG_SPECTRA_PERIODS = True  # Plot against periods for magnetic spectra plots (or frequencies)
    FigMisc.MARK_TEXC = True  # Add lines marking the main excitation periods/frequencies on Ae1 plot
    FigMisc.MARK_BEXC_MAX = True  # Whether to annotate excitation spectrum plots with label for highest peak
    FigLbl.peakLblSize = 14  # Font size in pt for highest-peak annotation
    FigMisc.Tmin_hr = None  # Cutoff period to limit range of Fourier space plots

    # Legends
    FigMisc.REFS_IN_LEGEND = True  # Hydrosphere plot: Whether to include reference profiles in legend
    FigMisc.wedgeLegendPos = 'center right'  # Wedge diagram: Where in axes added at right to place legend
    FigMisc.legendFontSize = 'x-small'  # Font size to use in legends, set by rcParams.

    # Latex settings
    FigMisc.fontFamily = 'serif'
    FigMisc.latexPackages = [
        r'\usepackage[version=4]{mhchem}',
        r'\usepackage{siunitx}',
        r'\usepackage{upgreek}'
    ]
    FigMisc.SetLatex()
    FigMisc.SetLatLon()
    if not FigMisc.TEX_INSTALLED:
        FigLbl.StripLatex()

    # Table printout settings
    FigMisc.PRINT_BULK = True  # Whether to print bulk body properties, like mass and MoI
    FigMisc.ALWAYS_SHOW_HP = True  # Whether to force HP ices and clathrates to be shown in DISP_* outputs to the terminal, even when none are present.
    FigMisc.ALWAYS_SHOW_PHI = True  # Whether to force porosity printout in DISP_* outputs
    FigMisc.LATEX_VLINES = False  # Whether to include vertical lines at table edges and between entries. Some journals do not allow them.
    FigMisc.LATEX_HLINES = False  # Whether to print horizontal lines between table entries.
    FigMisc.HF_HLINES = True  # Whether to print horizontal lines at head and foot of latex tables
    FigMisc.COMP_ROW = True  # Whether to force composition into a row instead of printing a separate summary table for each ocean comp
    FigMisc.BODY_NAME_ROW = True  # Whether to print a row with body name in bold in summary table

    # Contour labels
    FigMisc.cLabelSize = 10  # Font size in pt for contour labels
    FigMisc.cLabelPad = 5  # Padding in pt to set beside contour labels
    FigMisc.cLegendOpacity = 1.0  # Opacity of legend backgrounds in contour plots.

    # Colorbar settings
    FigMisc.cbarTitleSize = 'small'  # Font size specifier for colorbar titles
    FigMisc.cbarFmt = '%.1f'  # Format string to use for colorbar units
    FigMisc.cbarSize = '5%'  # Description of the size of colorbar to use with make_axes_locatable (secondary colorbars in phase space plots)
    FigMisc.cbarSpace = 0.5  # Amount of whitespace in inches to use for colorbars
    FigMisc.extraPad = FigMisc.cbarSpace * 0.8  # Amount of extra padding to apply to secondary colorbars
    FigMisc.nCbarPts = 80  # Number of points to use for drawing colorbar gradient
    FigMisc.SetFontSizes()  # Assign sizes for fonts

    # Assign metadata for correct output format combination
    FigLbl.SetMeta(FigMisc.figFormat)

    return Color, Style, FigLbl, FigSize, FigMisc
