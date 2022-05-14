""" Default figure settings """
import numpy as np
from PlanetProfile.Utilities.defineStructs import ColorStruct, StyleStruct, \
    FigLblStruct, FigSizeStruct, FigMiscStruct

configPlotsVersion = 1  # Integer number for config file version. Increment when new settings are added to the default config file.
Color = ColorStruct()
Style = StyleStruct()
FigLbl = FigLblStruct()
FigSize = FigSizeStruct()
FigMisc = FigMiscStruct()

""" Figure color options """
Color.Induction = {'synodic': 'blue', 'orbital': 'purple', 'true anomaly': 'green', 'synodic harmonic': 'goldenrod'}  # Colors for inductOgram plots
Color.ref = 'gray'

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
Color.innerCmapName = 'inferno'  # For plotting temperature profiles atop Perple_X data
# Alternative color options for silicates
Color.PALE_SILICATES = False  # Whether to use a lighter color scheme for silicate layers, or a more "orangey" saturated one
Color.paleSilPorousCmapName = 'terrain'
Color.paleSilPorousTop = 0.75
Color.paleSilPorousBot = 0.78
Color.paleSilCondCmapName = 'terrain'
Color.paleSilCondTop = 0.78 - 0.001
Color.paleSilCondBot = 0.78

# Colormaps for inductogram phase space plots, hydrosphere plots, etc
Color.cmapName = {
    'none': 'copper',
    'PureH2O': 'cividis',
    'Seawater': 'cool_r',
    'MgSO4': 'winter_r',
    'NH3': 'spring',
    'NaCl': 'summer',
    'Ice': 'coolwarm_r'
}
# Select only a subset of the available colormap, if we choose to
Color.cmapBounds = {
    'none': [0.0, 1.0],
    'PureH2O': [0.0, 1.0],
    'Seawater': [0.0, 1.0],
    'MgSO4': [0.0, 1.0],
    'NH3': [0.0, 1.0],
    'NaCl': [0.0, 1.0],
    'Ice': [0.2, 0.8]
}
# Set temperature bounds to use for colormap normalization
Color.Tbounds_K = [245.0, 300.0]

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


""" Figure style options """
Style.LS_dft = '-'  # Default line style to use on plots
Style.LS_Sw = '-'  # Linestyle for Seawater
Style.LS_Mg = '--'  # Linestyle for MgSO4
Style.LS_sp = ':'  # Linestyle for special consideration models
Style.LW_sal = 3  # Linewidth for higher salinity
Style.LW_dil = 1  # Linewidth for dilute salinity
Style.LW_std = 2  # Linewidth for standard salinity
Style.LW_sound = 1.5  # LineWidth for sound speed plots
Style.LW_seism = 1  # LineWidth for seismic plots (Attenuation)
Style.LS_ref = {'none': None, 'PureH2O': '-', 'Seawater': ':', 'MgSO4': '--', 'NH3': '--', 'NaCl': '--'}  # Style for reference profiles
Style.LW_ref = 0.75  # Linewidth for reference profiles
Style.LS_Induction = {'synodic': '-', 'orbital': ':', 'true anomaly': ':', 'synodic harmonic': '--'}  # Style for inductOgram plots
Style.LW_Induction = {'synodic': 1.5, 'orbital': 1.5, 'true anomaly': 1.5, 'synodic harmonic': 1.5}  # Widths for inductOgram plots
Style.MW_Induction = 2  # Marker size to use for induction scatter plots
Style.MS_Induction = 'o'  # Marker style for induction scatter plots

Style.wedgeAngle_deg = 25  # Angular size of wedge diagrams in degrees
Style.LW_wedge = 0.125  # Linewidth in pt for minor boundaries in wedge diagrams
Style.LW_wedgeMajor = 0.375  # Linewidth in pt for major layer boundaries in wedge diagrams


""" Figure labels """
# Note: Specific labels are set in PlanetProfile.Utilities.defineStructs.
FigLbl.NEGATIVE_UNIT_POWERS = True  # Whether to use negative powers for units in latex tables, or instead a backslash.
FigLbl.NAN_FOR_EMPTY = False  # Whether to use nan (or -) for empty layer parameters that were not calculated or not present.
FigLbl.w_IN_WTPCT = False  # Whether to print salinities in wt% (or g/kg) in tables
FigLbl.x_IN_WTPCT = True  # Whether to print silicate/core mass fractions in wt% (or g/kg) in tables
FigLbl.qSURF_IN_mW = True  # Whether to print qSurf in mW/m^2 (or W/m^2)
FigLbl.phi_IN_VOLPCT = False  # Whether to print porosity (phi) in vol% (or unitless volume fraction)
FigLbl.SetUnits()  # Make use of above toggles and assign labels


""" Figure sizes """
FigSize.vsP = (3, 3)
FigSize.vsR = (3, 3)
FigSize.vperm = (3, 3)
FigSize.vgsks = (3, 3)
FigSize.vseis = (3, 3)
FigSize.vhydro = (8, 5)
FigSize.vgrav = (6, 5)
FigSize.vmant = (6, 6)
FigSize.vcore = (6, 6)
FigSize.vpvt4 = (3, 3)
FigSize.vpvt6 = (3, 3)
FigSize.vwedg = (4.5, 4.5)
FigSize.phaseSpaceSolo = (6, 4)
FigSize.phaseSpaceCombo = (9, 4)
FigSize.induct = (8, 4)
FigSize.inductCombo = (8, 8)


""" Miscellaneous figure options """
# General figure options
FigMisc.figFormat = 'pdf'
FigMisc.dpi = 300  # Resolution in dots per inch for raster images (.png). Ignored for vector images (.pdf, .eps)
FigMisc.xtn = '.' + FigMisc.figFormat  # Figure file extension. Good options are .eps, .pdf, and .png
FigMisc.defaultFontName = 'STIXGeneral'  # Default font variables--STIX is what is used in Icarus journal submissions
FigMisc.defaultFontCode = 'stix'  # Code name for default font needed in some function calls
FigMisc.backupFont = 'Times New Roman'  # Backup font that looks similar to STIX that most users are likely to have

# Wedge diagrams
FigMisc.IONOSPHERE_IN_WEDGE = False  # Whether to include specified ionosphere in wedge diagram
FigMisc.WEDGE_ICE_TICKS = False  # Whether to print ticks for ice shell, which usually overlap with the body outer radius
FigMisc.DRAW_IONOS_BOUND = False  # Whether to draw a boundary line around the ionosphere
FigMisc.DRAW_CONVECTION_BOUND = False  # Whether to draw a boundary line between convecting and conducting regions
FigMisc.DRAW_POROUS_BOUND = False  # Whether to draw a boundary line between porous and non-porous materials
FigMisc.DRAW_FeS_BOUND = True  # Whether to draw a boundary line between Fe and FeS in the core

# Inductogram phase space plots
FigMisc.DARKEN_SALINITIES = False  # Whether to match hues to the colorbar, but darken points based on salinity, or to just use the colorbar colors.
FigMisc.NORMALIZED_SALINITIES = False  # Whether to normalize salinities to absolute concentrations relative to the saturation limit for each salt
FigMisc.NORMALIZED_TEMPERATURES = False  # Whether to normalize ocean mean temperatures to specified maxima and minima for the colormap
# Inductograms
FigMisc.PLOT_CONTOURS = True  # Contours or surfaces
FigMisc.PLOT_V2021 = True  # Mark the selected ocean/conductivity combos used in Vance et al. 2021
# Excitation spectra
FigMisc.DO_PER = True  # Convert frequency axes to periods for FFT plots

# Legends
FigMisc.LEGEND = True  # Whether to plot legends
FigMisc.REFS_IN_LEGEND = True  # Hydrosphere plot: Whether to include reference profiles in legend
FigMisc.hydroLegendBox = (0.15, 0.1, 0.33, 0.5)  # Hydrosphere plot: Bounding box for where to place legends. Values are x, y, dx, dy in fractions of the figure size, where x and y are for the bottom-left corner of the box.
FigMisc.hydroLegendPos = 'center left'  # Hydrosphere plot: Where to place legends within bounding box
FigMisc.wedgeLegendPos = 'center right'  # Wedge diagram: Where in axes added at right to place legend

# Latex settings
FigMisc.fontFamily = 'serif'
FigMisc.latexPackages = [
    r'\usepackage[version=4]{mhchem}',
    r'\usepackage{siunitx}',
    r'\usepackage{upgreek}'
]
FigMisc.SetLatex()

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
FigMisc.cbarSpace = 0.5  # Amount of whitespace in inches to use for colorbars
FigMisc.cbarSize = '5%'  # Description of the size of colorbar to use with make_axes_locatable
FigMisc.cbarHeight = 0.6  # Fraction of total figure height to use for colorbar size
FigMisc.cbarPad = 0.25  # Padding in pt to use for colorbars
FigMisc.extraPad = FigMisc.cbarSpace * 0.8  # Amount of extra padding to apply to secondary colorbars
FigMisc.cbarFmt = '%.1f'  # Format string to use for colorbar units
FigMisc.nCbarPts = 80  # Number of points to use for drawing colorbar gradient
FigMisc.cbarBottom = (1 - FigMisc.cbarHeight - FigMisc.cbarPad*2/72)/2  # Fraction of total figure height to use for bottom edge of colorbar
