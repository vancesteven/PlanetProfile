""" Default figure settings """
import shutil
import matplotlib.pyplot as plt

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
        self.vwedg = (3, 3)
        self.induct = (8, 4)
        self.inductCombo = (8, 8)


""" Figure color options """
class ColorsStruct:
    def __init__(self):
        self.Induction = {'synodic': 'blue', 'orbital': 'purple', 'true anomaly': 'green', 'synodic harmonic': 'goldenrod'}  # Colors for inductOgram plots
        self.ref = 'gray'
        self.cmap = 'inferno'
        self.coldestSw = 'c'
        self.warmestSw = '#b000ff'
        self.Sw_alt = (0, 175/255, 238/255)
        self.coldestMgSO4 = 'b'
        self.warmestMgSO4 = 'xkcd:purpley blue'
        # Comments were copied from Matlab color options and need updating before use.
        # self.col_contSyn = cfg.cmap(floor(100*cc),:)
        # self.col_contOrb = cfg.cmap(floor( 10*cc),:)
        # self.col_contHrm = cfg.cmap(floor(200*cc),:)
        # self.col_Sw = summer(200)
        # self.col_midColdSw = cfg.col_Sw(25,:)
        # self.col_middestSw = cfg.col_Sw(50,:)
        # self.col_midWarmSw = cfg.col_Sw(75,:)
        # self.col_MgSO4 = cool(133)
        # self.col_midColdMgSO4 = cfg.col_MgSO4(25,:)
        # self.col_middestMgSO4 = cfg.col_MgSO4(50,:)
        # self.col_midWarmMgSO4 = cfg.col_MgSO4(75,:)

        # Wedge diagram color options
        self.IonosphereTop = [1, 0, 1]
        self.Ionosphere = [1, 0, 1]
        self.IonosphereBot = [1, 0, 1]
        self.IceI = '#d3eefb'
        self.IceII = '#76b6ff'
        self.IceIII = '#a8deef'
        self.IceV = '#83d4f6'
        self.IceVI = '#cee5ea'
        self.Clath = '#86bcb8'
        self.OceanTop = [134/255, 149/255, 201/255] #'#4babdf'
        self.OceanBot = [45/255, 55/255, 100/255]
        self.Rock = [101/255, 46/255, 11/255]
        self.Core = [141/255, 122/255, 121/255]


""" Figure style options """
class StylesStruct:
    def __init__(self):
        self.LS_dft = '-'  # Default line style to use on plots
        self.LS_Sw = '-'  # linestyle for Seawater
        self.LS_Mg = '--'  # linestyle for MgSO4
        self.LS_sp = ':'  # linestyle for special consideration models
        self.LW_sal = 3  # linewidth for higher salinity
        self.LW_dil = 1  # linewidth for dilute salinity
        self.LW_std = 2  # linewidth for standard salinity
        self.LW_sound = 1.5  # LineWidth for sound speed plots
        self.LW_seism = 1  # LineWidth for seismic plots (Attenuation)
        self.LS_ref = {'none': None, 'PureH2O': '-', 'Seawater': ':', 'MgSO4': '--', 'NH3': '--', 'NaCl': '--'}  # Style for reference profiles
        self.LW_ref = 0.75  # Width for reference profiles
        self.LS_Induction = {'synodic': '-', 'orbital': ':', 'true anomaly': ':', 'synodic harmonic': '--'}  # Style for inductOgram plots
        self.LW_Induction = {'synodic': 1.5, 'orbital': 1.5, 'true anomaly': 1.5, 'synodic harmonic': 1.5}  # Widths for inductOgram plots


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
        self.LEGEND = True  # Whether to plot legends
        self.LegendPosition = 'right'  # Where to place legends when forced
        self.refsInLegend = True  # Whether to include reference profiles in legend
        plt.rcParams['font.family'] = 'serif'  # Choose serif font for figures to best match math mode variables in body text
        plt.rcParams['font.serif'] = self.defaultFontName  # Set plots to use the default font

        # Check if Latex executable is on the path so we can use backup options if Latex is not installed
        if shutil.which('latex'):
            plt.rcParams['text.usetex'] = True  # Use Latex interpreter to render text on plots
            # Load in font package in Latex
            plt.rcParams['text.latex.preamble'] = f'\\usepackage{{{self.defaultFontCode}}}' + \
                r'\usepackage[version=4]{mhchem}' + r'\usepackage{siunitx}' + r'\usepackage{upgreek}'
            self.TEX_INSTALLED = True
        else:
            print('A LaTeX installation was not found. Some plots may have fallback options in labels.')
            plt.rcParams['font.serif'] += ', ' + self.backupFont  # Set plots to use the default font if installed, or a backup if not
            plt.rcParams['mathtext.fontset'] = self.defaultFontCode
            self.TEX_INSTALLED = False

        self.cLabelSize = 10  # Font size in pt for contour labels
        self.cLabelPad = 5  # Padding in pt to set beside contour labels
        self.cLegendOpacity = 1.0  # Opacity of legend backgrounds in contour plots.


FigSize = FigSizeStruct()
Color = ColorsStruct()
Style = StylesStruct()
FigMisc = MiscStruct()
