""" Default figure settings """
import numpy as np
import spiceypy as spice
from PlanetProfile.Utilities.defineStructs import ColorStruct, StyleStruct, \
    FigLblStruct, FigSizeStruct, FigMiscStruct

configPlotsVersion = 9  # Integer number for config file version. Increment when new settings are added to the default config file.
Color = ColorStruct()
Style = StyleStruct()
FigLbl = FigLblStruct()
FigSize = FigSizeStruct()
FigMisc = FigMiscStruct()

""" Figure color options """
Color.Induction = {'synodic': 'blue', 'orbital': 'purple', 'true anomaly': 'green', 'synodic harmonic': 'goldenrod'}  # Colors for inductOgram plots
Color.ref = 'gray'
Color.geotherm = 'white'
Color.BdipInset = 'black'  # Color for inset box of surface induced dipole strength plots

# Wedge diagram color options
Color.none = '#FFFFFF00'
Color.wedgeBd = 'black'
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

# PvT properties colormaps
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
    'PureH2O': [0.0, 1.0],
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
Color.CAdot = 'black'
Color.thresh = 'blue'


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
Style.LS_ref = {'none': None, 'PureH2O': '-', 'Seawater': ':', 'MgSO4': '--', 'NH3': '--', 'NaCl': '--'}  # Style for reference profiles
Style.LW_ref = 0.75  # Linewidth for reference profiles
Style.LS_Induction = {'synodic': '-', 'orbital': ':', 'true anomaly': ':', 'synodic harmonic': '--'}  # Style for inductOgram plots
Style.LW_Induction = {'synodic': 1.5, 'orbital': 1.5, 'true anomaly': 1.5, 'synodic harmonic': 1.5}  # Widths for inductOgram plots
Style.MW_Induction = 2  # Marker size to use for induction scatter plots
Style.MS_Induction = 'o'  # Marker style for induction scatter plots

# Wedge diagrams
Style.wedgeAngle_deg = 25  # Angular size of wedge diagrams in degrees
Style.LW_wedge = 0.125  # Linewidth in pt for minor boundaries in wedge diagrams
Style.LW_wedgeMajor = 0.375  # Linewidth in pt for major layer boundaries in wedge diagrams

# Complex dipole plots
Style.MW_dip = {'synodic': 5, 'orbital': 6, 'true anomaly': 6, 'synodic harmonic': 4.5}  # Marker size for each period in complex dipole plots
Style.MS_dip = {'synodic': '*', 'orbital': 'o', 'true anomaly': 'P', 'synodic harmonic': 'h'}  # Marker style for each period in complex dipole plots
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
Style.MW_CA = 5  # Marker size for closest approach dots


""" Figure labels """
# Note: Specific labels are set in PlanetProfile.Utilities.defineStructs.
FigLbl.NEGATIVE_UNIT_POWERS = True  # Whether to use negative powers for units in latex tables, or instead a backslash.
FigLbl.BODYNAME_IN_LABEL = True  # Whether to include the bodyname in legend labels for comparison plots including multiple bodies
FigLbl.NAN_FOR_EMPTY = False  # Whether to use nan (or -) for empty layer parameters that were not calculated or not present.
FigLbl.PFULL_IN_GPa = True  # Whether to plot P in GPa (or MPa) for full-body plots
FigLbl.PHYDRO_IN_bar = False  # Whether to print P in bar (or MPa) for hydrosphere plots
FigLbl.w_IN_WTPCT = False  # Whether to print salinities in wt% (or g/kg) in tables
FigLbl.T_IN_C = False  # Whether to print T in °C (or K) in plots
FigLbl.x_IN_MOLPCT = True  # Whether to print silicate/core mass fractions in mol% (or fractional) in tables and plots
FigLbl.qSURF_IN_mW = True  # Whether to print qSurf in mW/m^2 (or W/m^2)
FigLbl.phi_IN_VOLPCT = False  # Whether to print porosity (phi) in vol% (or unitless volume fraction)
FigLbl.PVT_CBAR_LABELS = False  # Whether to add short labels identifying silicate/core colorbars in PvT properties plots
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
FigSize.vpvt = (12, 6)
FigSize.vwedg = (4.5, 4.5)
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
FigSize.MagSurf = (8, 5)
FigSize.MagSurfCombo = (16, 5)
FigSize.MagCA = (4, 4)
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

# Hydrosphere plots
FigMisc.SCALE_HYDRO_LW = True  # Whether to adjust thickness of lines on hydrosphere plot according to relative salinity
FigMisc.MANUAL_HYDRO_COLORS = True  # Whether to set color of lines in hydrosphere according to melting temperature
FigMisc.RELATIVE_Tb_K = True  # Whether to set colormap of lines based on relative comparison (or fixed settings in ColorStruct)
FigMisc.TminHydro = 200  # Minimum temperature to display on hydrosphere plots

# Wedge diagrams
FigMisc.IONOSPHERE_IN_WEDGE = False  # Whether to include specified ionosphere in wedge diagram
FigMisc.WEDGE_ICE_TICKS = False  # Whether to print ticks for ice shell, which usually overlap with the body outer radius
FigMisc.DRAW_IONOS_BOUND = False  # Whether to draw a boundary line around the ionosphere
FigMisc.DRAW_CONVECTION_BOUND = False  # Whether to draw a boundary line between convecting and conducting regions
FigMisc.DRAW_POROUS_BOUND = False  # Whether to draw a boundary line between porous and non-porous materials
FigMisc.DRAW_FeS_BOUND = True  # Whether to draw a boundary line between Fe and FeS in the core
FigMisc.minzbRratio_frac = 0.05  # Fraction of total body radius for ice shell and ocean thickness, above which ice shell ticks will automatically switch on (overrides WEDGE_ICE_TICKS)

# Silicate/core PT diagrams
FigMisc.nTgeo = 80  # Number of temperature points to evaluate/plot for PT property plots
FigMisc.nPgeo = 100  # Number of pressure points to evaluate/plot for PT property plots
FigMisc.nPgeoCore = 40  # Subset of nPgeo to use for core, if present
FigMisc.PVT_INCLUDE_CORE = True  # Whether to include core as well as silicates in PT properties diagrams
        
# Induced field surface strength plots
FigMisc.BdipZoomMult = 1.05  # Extra space to include around zoomed-in part, in fraction of largest value.
FigMisc.SHOW_INSET = True  # Whether to show the inset box for the zoom-in plot, when applicable
FigMisc.rMagEval_Rp = 1.0  # Fraction of body radius to use for surface over which PlotMagSurface is evaluated
FigMisc.tMagLbl = None  # List of strings to use to describe the times listed in tMagEval_s
FigMisc.tMagEval_s = spice.str2et('2000-01-01T11:58:55.816')  # Time in seconds past J2000 at which to evaluate magnetic field surface maps. spiceypy.str2et will convert a datetime string to the correct format. Accepts an array.
FigMisc.LARGE_ADJUST = True  # Whether to make certain labels better for cramped spaces, including removing colorbars (True is more pared down, and overrides nLon/LatTicks below.)
FigMisc.BASYM_WITH_SYM = True  # Whether to plot Basym plot and Bsym plot on the same figure
FigMisc.vCompMagSurf = 'mag'  # Component to use for induced field surface strength plots. Options are ['x', 'y', 'z', 'mag'].
FigMisc.nPPGCmapRes = 180  # Number of points per great circle to use for map angular resolution
FigMisc.DO_360 = False  # Whether to range longitudes from 0 to 360 or from -180 to +180
FigMisc.nLatTicks = 7  # Number of ticks to mark on latitude axis
FigMisc.nLonTicks = 9  # Number of ticks to mark on longitude axis
FigMisc.nMagContours = 9  # Number of contour intervals to mark on magnetic plots
FigMisc.nAsymContours = 7  # Number of contour intervals to mark on asymmetry maps
FigMisc.latlonSize = 14  # Font size for lat/lon labels
FigMisc.cLabelSize = 12  # Font size for contour labels
FigMisc.cLabelPad = 1  # Padding in pt to use for contour labels
FigMisc.vminMagSurf_nT = None  # Minimum value for colormap in magnetic field surface plots (None uses data min/max)
FigMisc.vmaxMagSurf_nT = None  # Minimum value for colormap in magnetic field surface plots (None uses data min/max)
FigMisc.vminMagSurfDiff_nT = None  # Minimum value for difference colormap in magnetic field surface plots (None uses data min/max)
FigMisc.vmaxMagSurfDiff_nT = None  # Minimum value for difference colormap in magnetic field surface plots (None uses data min/max)

# Magnetic field trajectory and CA plots
FigMisc.CAlblSize = 12  # Size of text labels on CA points
FigMisc.SHOW_MAG_THRESH = True  # Whether to show a line indicating the precision floor of a magnetometer
FigMisc.thresh_nT = 0.0488  # Precision floor in nT for magnetometer to plot
FigMisc.threshCenter = 100  # x coordinate to place the MAG floor label
FigMisc.hCAmax_km = 500  # Maximum altitude to show on CA plot

# Inductogram phase space plots
FigMisc.DARKEN_SALINITIES = False  # Whether to match hues to the colorbar, but darken points based on salinity, or to just use the colorbar colors.
FigMisc.NORMALIZED_SALINITIES = False  # Whether to normalize salinities to absolute concentrations relative to the saturation limit for each salt
FigMisc.NORMALIZED_TEMPERATURES = False  # Whether to normalize ocean mean temperatures to specified maxima and minima for the colormap
# Inductograms
FigMisc.PLOT_V2021 = False  # Mark the selected ocean/conductivity combos used in Vance et al. 2021
# Excitation spectra
FigMisc.MAG_SPECTRA_PERIODS = True  # Plot against periods for magnetic spectra plots (or frequencies)
FigMisc.MARK_TEXC = True  # Add lines marking the main excitation periods/frequencies on Ae1 plot

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

# Contour labels
FigMisc.cLabelSize = 10  # Font size in pt for contour labels
FigMisc.cLabelPad = 5  # Padding in pt to set beside contour labels
FigMisc.cLegendOpacity = 1.0  # Opacity of legend backgrounds in contour plots.

# Colorbar settings
FigMisc.cbarTitleSize = 'small'  # Font size specifier for colorbar titles
FigMisc.cbarFmt = '%.1f'  # Format string to use for colorbar units
FigMisc.cbarSpace = 0.5  # Amount of whitespace in inches to use for colorbars
FigMisc.extraPad = FigMisc.cbarSpace * 0.8  # Amount of extra padding to apply to secondary colorbars
FigMisc.nCbarPts = 80  # Number of points to use for drawing colorbar gradient
FigMisc.SetFontSizes()  # Assign sizes for fonts
