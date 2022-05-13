""" Default figure settings """
import shutil
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import get_cmap
from matplotlib.colors import rgb_to_hsv, hsv_to_rgb
from scipy.interpolate import interp1d
import cmasher

""" Figure sizes """
class FigSizeStruct:
    def __init__(self):
        self.vsP = (3, 3)
        self.vsR = (3, 3)
        self.vperm = (3, 3)
        self.vgsks = (3, 3)
        self.vseis = (3, 3)
        self.vhydro = (8, 5)
        self.vgrav = (6, 5)
        self.vmant = (6, 6)
        self.vcore = (6, 6)
        self.vpvt4 = (3, 3)
        self.vpvt6 = (3, 3)
        self.vwedg = (4.5, 4.5)
        self.phaseSpaceSolo = (6, 4)
        self.phaseSpaceCombo = (9, 4)
        self.induct = (8, 4)
        self.inductCombo = (8, 8)


""" Figure color options """
class ColorsStruct:
    def __init__(self):
        self.Induction = {'synodic': 'blue', 'orbital': 'purple', 'true anomaly': 'green', 'synodic harmonic': 'goldenrod'}  # Colors for inductOgram plots
        self.ref = 'gray'
        self.innerCmap = get_cmap('inferno')

        # Wedge diagram color options
        self.none = '#FFFFFF00'
        self.wedgeBd = 'black'
        self.ionoCmapName = 'RdPu'
        self.ionoTop = 0.0
        self.ionoBot = 0.2
        self.ionoCmap = cmasher.get_sub_cmap(self.ionoCmapName, self.ionoTop, self.ionoBot)
        self.ionoN = 10
        self.iceIcond = 'xkcd:ice blue'
        self.iceIconv = 'xkcd:robin\'s egg blue'
        self.iceII = 'xkcd:pale sky blue'
        self.iceIII = 'xkcd:carolina blue'
        self.iceV = 'xkcd:light periwinkle'
        self.iceVI = '#91d1d4'
        self.clathCond = 'xkcd:pastel blue'
        self.clathConv = 'xkcd:dusty blue'
        self.oceanCmapName = 'ocean_r'
        self.oceanTop = 0.2  # Fraction of ocean colormap to start at
        self.oceanBot = 0.45  # Fraction of ocean colormap to end at
        self.oceanCmap = cmasher.get_sub_cmap(self.oceanCmapName, self.oceanTop, self.oceanBot)
        self.oceanN = 20
        # self.silPorousCmapName = 'terrain'
        # self.silPorousTop = 0.75
        # self.silPorousBot = 0.78
        self.silPorousCmapName = 'BrBG'
        self.silPorousTop = 0.0
        self.silPorousBot = 0.05
        self.silPorousCmap = cmasher.get_sub_cmap(self.silPorousCmapName, self.silPorousTop, self.silPorousBot)
        self.silPorousN = 5
        # self.silCondCmapName = 'terrain'
        # self.silCondTop = 0.78 - 0.001
        # self.silCondBot = 0.78
        self.silCondCmapName = 'BrBG'
        self.silCondTop = 0.05 - 0.001
        self.silCondBot = 0.05
        self.silCondCmap = cmasher.get_sub_cmap(self.silCondCmapName, self.silCondTop, self.silCondBot)
        self.silCondN = 1
        self.silConvCmapName = 'gist_heat'
        self.silConvTop = 0.5
        self.silConvBot = 0.9
        self.silConvCmap = cmasher.get_sub_cmap(self.silConvCmapName, self.silConvTop, self.silConvBot)
        self.silConvN = 10
        self.FeS = 'xkcd:puke'
        self.Fe = '#2d3639'

        self.cmapName = {
            'none': 'copper',
            'PureH2O': 'cividis',
            'Seawater': 'cool_r',
            'MgSO4': 'winter_r',
            'NH3': 'spring',
            'NaCl': 'summer',
            'Ice': 'coolwarm_r'
        }
        # Select only a subset of the available colormap, if we choose to
        self.cmapBounds = {
            'none': [0.0, 1.0],
            'PureH2O': [0.0, 1.0],
            'Seawater': [0.0, 1.0],
            'MgSO4': [0.0, 1.0],
            'NH3': [0.0, 1.0],
            'NaCl': [0.0, 1.0],
            'Ice': [0.2, 0.8]
        }
        # Use cmasher to return colormap objects that do the down-select for us
        self.cmap = {comp: cmasher.get_sub_cmap(cmap, self.cmapBounds[comp][0], self.cmapBounds[comp][1])
                     for comp, cmap in self.cmapName.items()}
        # Set temperature bounds to use for colormap normalization
        self.Tbounds_K = [245.0, 300.0]
        self.GetNormT = interp1d([self.Tbounds_K[0], self.Tbounds_K[1]], [0.0, 1.0],
                                 bounds_error=False, fill_value='extrapolate')
        # Set upper bounds for max concentrations
        self.saturation = {
            'none': 1.0,
            'PureH2O': 1.0,
            'Seawater': 304.0,
            'MgSO4': 282.0,
            'NH3': 100.0,
            'NaCl': 304.0
        }
        # Saturation & color brightness ("value" in HSV) values for salinity/conductivity axis bounds
        self.fresh = [0.5, 1.0]
        self.salty = [1.0, 0.5]
        self.getSat = interp1d([0.0, 1.0], [self.fresh[0], self.salty[0]], bounds_error=False, fill_value=self.salty[0])
        self.getVal = interp1d([0.0, 1.0], [self.fresh[1], self.salty[1]], bounds_error=False, fill_value=self.salty[1])

    def OceanCmap(self, comps, w_normFrac, Tmean_normFrac, DARKEN_SALINITIES=True):
        """ Get colormap RGBA vectors for each salinity/T combination.

            Args:
                comps (string, shape N): Ocean composition string for each point.
                w_normFrac (float, shape N): Normalized ocean salinities as a fraction
                    of the saturation/maximum concentration used for the colormap.
                Tmean_normFrac (float, shape N): Normalized mean ocean temperatures
                    as a fraction of the range of values to use for the colormap, where
                    0 is at the bottom of the range and 1 is at the top.

            Returns:
                cList (float, shape N x 3): RGB vectors as rows for each combination.
        """

        if DARKEN_SALINITIES:
            # Get the hue for each point by getting the colormap entry for the normalized
            # temperature, then stripping off the alpha channel, then converting to HSV,
            # then taking only the first value in the [H,S,V] output.
            hueMap = np.array([rgb_to_hsv(self.cmap[comp](T)[:3])[0] for comp, T in zip(comps, Tmean_normFrac)])
            # Get the saturation and value lists from the min/max bounds for the
            # normalized salinity
            satMap = self.getSat(w_normFrac)
            valMap = self.getVal(w_normFrac)

            cList = hsv_to_rgb(np.column_stack((hueMap, satMap, valMap)))
        else:
            # Just use standard evaluation of the colorbar
            cList = np.row_stack([self.cmap[comp](T) for comp, T in zip(comps, Tmean_normFrac)])

        return cList


""" Figure style options """
class StylesStruct:
    def __init__(self):
        self.LS_dft = '-'  # Default line style to use on plots
        self.LS_Sw = '-'  # Linestyle for Seawater
        self.LS_Mg = '--'  # Linestyle for MgSO4
        self.LS_sp = ':'  # Linestyle for special consideration models
        self.LW_sal = 3  # Linewidth for higher salinity
        self.LW_dil = 1  # Linewidth for dilute salinity
        self.LW_std = 2  # Linewidth for standard salinity
        self.LW_sound = 1.5  # LineWidth for sound speed plots
        self.LW_seism = 1  # LineWidth for seismic plots (Attenuation)
        self.LS_ref = {'none': None, 'PureH2O': '-', 'Seawater': ':', 'MgSO4': '--', 'NH3': '--', 'NaCl': '--'}  # Style for reference profiles
        self.LW_ref = 0.75  # Width for reference profiles
        self.LS_Induction = {'synodic': '-', 'orbital': ':', 'true anomaly': ':', 'synodic harmonic': '--'}  # Style for inductOgram plots
        self.LW_Induction = {'synodic': 1.5, 'orbital': 1.5, 'true anomaly': 1.5, 'synodic harmonic': 1.5}  # Widths for inductOgram plots
        self.MW_Induction = 2  # Marker size to use for induction scatter plots
        self.MS_Induction = 'o'  # Marker style for induction scatter plots

        self.wedgeAngle_deg = 25  # Angular size of wedge diagrams in degrees
        self.LW_wedge = 0.125  # Linewidth in pt for minor boundaries in wedge diagrams
        self.LW_wedgeMajor = 0.375  # Linewidth in pt for major layer boundaries in wedge diagrams


""" Miscellaneous figure options """
class MiscStruct:
    def __init__(self):
        # General figure options
        self.figFormat = 'pdf'
        self.dpi = 300  # Resolution in dots per inch for raster images (.png). Ignored for vector images (.pdf, .eps)
        self.xtn = '.' + self.figFormat  # Figure file extension. Good options are .eps, .pdf, and .png
        self.defaultFontName = 'STIXGeneral'  # Default font variables--STIX is what is used in Icarus journal submissions
        self.defaultFontCode = 'stix'  # Code name for default font needed in some function calls
        self.backupFont = 'Times New Roman'  # Backup font that looks similar to STIX that most users are likely to have
        self.IONOSPHERE_IN_WEDGE = False  # Whether to include specified ionosphere in wedge diagram
        self.DRAW_IONOS_BOUND = False  # Whether to draw a boundary line around the ionosphere
        self.DRAW_CONVECTION_BOUND = False  # Whether to draw a boundary line between convecting and conducting regions
        self.DRAW_POROUS_BOUND = False  # Whether to draw a boundary line between porous and non-porous materials
        self.DRAW_FeS_BOUND = True  # Whether to draw a boundary line between Fe and FeS in the core
        self.WEDGE_ICE_TICKS = False  # Whether to print ticks for ice shell, which usually overlap with the body outer radius

        self.DARKEN_SALINITIES = False  # For inductogram phase space plots, whether to match hues to the colorbar, but darken points based on salinity, or to just use the colorbar colors.
        self.NORMALIZED_SALINITIES = False  # For inductogram phase space plots, whether to normalize salinities to absolute concentrations relative to the saturation limit for each salt
        self.NORMALIZED_TEMPERATURES = False  # For inductogram phase space plots, whether to normalize ocean mean temperatures to specified maxima and minima for the colormap

        self.LEGEND = True  # Whether to plot legends
        self.hydroLegendBox = (0.15, 0.1, 0.33, 0.5)  # Hydrosphere plot: Bounding box for where to place legends. Values are x, y, dx, dy in fractions of the figure size, where x and y are for the bottom-left corner of the box.
        self.hydroLegendPos = 'center left'  # Hydrosphere plot: Where to place legends within bounding box
        self.refsInLegend = True  # Hydrosphere plot: Whether to include reference profiles in legend
        self.wedgeLegendPos = 'center right'  # Wedge diagram: Where in axes added at right to place legend

        plt.rcParams['font.family'] = 'serif'  # Choose serif font for figures to best match math mode variables in body text
        plt.rcParams['font.serif'] = self.defaultFontName  # Set plots to use the default font
        # Check if Latex executable is on the path so we can use backup options if Latex is not installed
        latexPackages = [
            r'\usepackage[version=4]{mhchem}',
            r'\usepackage{siunitx}',
            r'\usepackage{upgreek}'
        ]
        if shutil.which('latex'):
            plt.rcParams['text.usetex'] = True  # Use Latex interpreter to render text on plots
            # Load in font and option packages in Latex
            plt.rcParams['text.latex.preamble'] = f'\\usepackage{{{self.defaultFontCode}}}' \
                                                  + ''.join(latexPackages)
            self.TEX_INSTALLED = True
        else:
            print('A LaTeX installation was not found. Some plots may have fallback options in labels.')
            plt.rcParams['font.serif'] += ', ' + self.backupFont  # Set plots to use the default font if installed, or a backup if not
            plt.rcParams['mathtext.fontset'] = self.defaultFontCode
            self.TEX_INSTALLED = False

        packageCmds = r'\n        '.join(latexPackages)
        self.latexPreamble = f"""For table setup, add this to your document preamble:
        {packageCmds}
        \sisetup{{group-separator={{\,}}, group-minimum-digits={{5}}, group-digits={{integer}}}}
        """

        # Table printout settings
        self.PRINT_BULK = True  # Whether to print bulk body properties, like mass and MoI
        self.ALWAYS_SHOW_HP = True  # Whether to force HP ices and clathrates to be shown in DISP_* outputs to the terminal, even when none are present.
        self.ALWAYS_SHOW_PHI = True  # Whether to force porosity printout in DISP_* outputs
        self.NEGATIVE_UNIT_POWERS = True  # Whether to use negative powers for units in latex tables, or instead a backslash.
        self.NAN_FOR_EMPTY = False  # Whether to use nan (or -) for empty layer parameters that were not calculated or not present.
        self.w_IN_WTPCT = False  # Whether to print salinities in wt% (or g/kg) in tables
        self.qSURF_IN_mW = True  # Whether to print qSurf in mW/m^2 (or W/m^2)
        self.phi_IN_VOLPCT = False  # Whether to print porosity (phi) in vol% (or unitless volume fraction)
        self.LATEX_VLINES = False  # Whether to include vertical lines at table edges and between entries. Some journals do not allow them.
        self.LATEX_HLINES = False  # Whether to print horizontal lines between table entries.
        self.HF_HLINES = True  # Whether to print horizontal lines at head and foot of latex tables

        self.cLabelSize = 10  # Font size in pt for contour labels
        self.cLabelPad = 5  # Padding in pt to set beside contour labels
        self.cLegendOpacity = 1.0  # Opacity of legend backgrounds in contour plots.
        
        self.cbarSpace = 0.5  # Amount of whitespace in inches to use for colorbars
        self.cbarSize = '5%'  # Description of the size of colorbar to use with make_axes_locatable
        self.cbarHeight = 0.6  # Fraction of total figure height to use for colorbar size
        self.cbarPad = 0.25  # Padding in pt to use for colorbars
        self.extraPad = self.cbarSpace * 0.8  # Amount of extra padding to apply to secondary colorbars
        self.cbarFmt = '%.1f'  # Format string to use for colorbar units
        self.nCbarPts = 80  # Number of points to use for drawing colorbar gradient
        self.cbarBottom = (1 - self.cbarHeight - self.cbarPad*2/72)/2  # Fraction of total figure height to use for bottom edge of colorbar


class FigLabelStruct:
    def __init__(self):
        # Wedge diagram labels
        self.wedgeTitle = 'interior structure diagram'
        self.wedgeRadius = r'Radius ($\mathrm{km}$)'

        # Inductogram labels
        self.plotTitles = ['Amplitude $A$', '$B_x$ component', '$B_y$ component', '$B_z$ component']
        self.fLabels = ['Amp', 'Bx', 'By', 'Bz']
        self.compEnd = ''
        self.phaseTitle = r'Phase delay $\upphi$ ($^\circ$)'
        self.sigLabel = r'Mean conductivity $\overline{\sigma}$ ($\si{S/m}$)'
        self.Dlabel = r'Ocean thickness $D$ ($\si{km}$)'
        self.wLabel = r'Salinity $w$ ($\si{g/kg}$)'
        self.TbLabel = r'Ice bottom temp $T_b$ ($\si{K}$)'
        self.rhoLabel = r'Silicate density $\rho_\mathrm{sil}$ ($\si{kg/m^3}$)'
        self.phiLabel = r'Seafloor porosity $\phi_\mathrm{sil}$ (void frac)'
        self.iceThickLbl = r'Ice shell thickness ($\si{km}$)'
        self.oceanTempLbl = r'Mean ocean temp ($\si{K}$)'
        self.wScale = 'log'
        self.sigScale = 'log'
        self.Dscale = 'log'

        # Induction parameter-dependent settings
        self.phaseSpaceTitle = None
        self.inductionTitle = None
        self.inductCompareTitle = None
        self.sigLims = None
        self.Dlims = None
        self.legendTexc = None
        self.yLabelInduct = None
        self.yScaleInduct = None

    def singleComp(self, comp):
        # Set a tag to append to titles in the event all of what we're plotting
        # has a single composition, for additional clarity.
        self.compEnd = f', \ce{{{comp}}} ocean'

    def setInduction(self, bodyname, IndParams, Texc_h):
        # Set titles, labels, and axis settings pertaining to inductogram plots
        self.phaseSpaceTitle = f'\\textbf{{{bodyname} interior phase space}}'
        self.inductionTitle = f'\\textbf{{{bodyname} induction response{self.compEnd}}}'
        self.inductCompareTitle = f'\\textbf{{{bodyname} induction response on different axes{self.compEnd}}}'

        self.sigLims = [10**IndParams.sigmaMin[bodyname], 10**IndParams.sigmaMax[bodyname]]
        self.Dlims = [10**IndParams.Dmin[bodyname], 10**IndParams.Dmax[bodyname]]
        self.legendTexc = np.array([f'{T_h:.2f} h' for T_h in Texc_h])

        if IndParams.inductOtype != 'sigma':
            if IndParams.inductOtype == 'Tb':
                self.yLabelInduct = self.TbLabel
                self.yScaleInduct = 'linear'
            elif IndParams.inductOtype == 'rho':
                self.yLabelInduct = FigLbl.rhoLabel
                self.yScaleInduct = 'linear'
            elif IndParams.inductOtype == 'phi':
                self.yLabelInduct = FigLbl.phiLabel
                self.yScaleInduct = 'log'



FigSize = FigSizeStruct()
Color = ColorsStruct()
Style = StylesStruct()
FigMisc = MiscStruct()
FigLbl = FigLabelStruct()
