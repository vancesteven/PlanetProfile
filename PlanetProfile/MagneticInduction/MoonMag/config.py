""" Contains runtime options for the program eval_induced_field.py
    and dependent functions.
    Developed in Python 3.8 for "A perturbation method for evaluating the
    magnetic field induced from an arbitrary, asymmetric ocean world 
    analytically" by Styczinski et al.
    DOI: 10.1016/j.icarus.2021.114840
Author: M. J. Styczinski, mjstyczi@uw.edu """   

import logging, shutil

# @@@@@@@@@@@@@@@@@
#   Main settings
# @@@@@@@@@@@@@@@@@

ppgc = 540  # Points per great circle (sets longitudinal resolution for plots)

verbose = True
use_latex = True  # Whether to use Latex in rendering plot text
relative = False  # Whether to use relative (epsilon/R*chi_pq) or absolute (direct spherical harmonic coefficients) formulation for reading Ypq. Affects which files are needed/read to interpret the asymmetric interior structure.
synodic_period_only = False  # Whether to consider only the synodic period for induction
orbital_time_series = False  # Whether to plot a time series that considers only the orbital period
sub_planet_vert = False  # Whether to plot a vertical cut at the sub-planetary point (or nearby) for a snapshot in time
plot_diffs = True  # Whether to plot magnetic fields as differences arising due to asymmetry
gif_diff = True  # Whether to plot differences in animation frames or the absolute component

nprm_max_main = 1  # Highest degree in excitation field to use
eval_radius = 2  # Distance (in units of body radii) from body center to use for evaluating B. Overridden for Europa and Enceladus.

debug = False  # Special use debug flag
convert_depth_to_chipq = True  # Prints a file named interior/chi_pq_bodyname.txt for copying over relative harmonic coefficients to degree<p>_shapes_bodyname.txt files. Only used if relative = False.
output_Schmidt = False  # Record induced moment coefficients in Schmidt semi-normalized form in addition to fully normalized moments. Also reads in and plots fields using this normalization.

# @@@@@@@@@@@@@@@@@@@@@@@@
#   Calculation settings
# @@@@@@@@@@@@@@@@@@@@@@@@

do_parallel = True  # Whether to run certain calculations in parallel. Setting to False makes debugging much more straightforward.
digits_precision = 750  # Number of digits of precision for use in calculations involving Bessel functions. Calculations involve small differences between large numbers when conductivities are very large or very small.
                        # If you encounter divide-by-zero errors, try increasing this number.

default_eval_datetime = "2000-01-01T11:58:55.816"  # UTC datetime string for the time to evaluate the induced fields (and start animations from) in yyyy-mm-ddThh:mm:mm:ss.sss format

# The following are only used if synodic_period_only is True:
loclat = 0  # degrees latitude (IAU) for plotting a time series
loclon = 0  # degrees longitude (IAU) for plotting a time series
localt = 0  # km altitude for plotting a time series

# Model option strings
prevEuropa = 'prev'
TobieHigh = 'Tobie_high'
TobieLow = 'Tobie_low'
CochraneC3C9 = 'Cochrane_C3C9only'
CochraneAll = 'Cochrane_all'
DO_LARGE = 'large'
DO_GIF = 'do_gif'

# @@@@@@@@@@@@@@@@@@@
#   Plotting params
# @@@@@@@@@@@@@@@@@@@

t_pts = 200  # Number of points in 1-period time series (or spatial points in vertical cut)
n_frames = 100  # Number of animation frames to use, spaced evenly throughout one period of oscillation
vert_start = 1.0  # Units of body radii to start vertical cut
vert_stop = 1.0 + 2000/1561  # And where to stop it
vert_cut_hr = 0.7  # Time in hrs past J2000 to perform vertical cut
vert_cut_lat = 30  # Lat/lon to project upward from the surface
vert_cut_lon = 0

# Color and style options
c = [ "green", "black", "blue" ]
vc = "brown"
style1 = "solid"
style2 = "dashed"
cMark = "xkcd:grass green"
cMarkE = "xkcd:forest green"
szMark = 20
stMark = 'o'

# Figure formatting defaults
deft_figsize = (8,4)
lat_min = -90
lat_max = 90
do_360 = False  # Whether to plot from 0 to 360 (True) or from -180 to 180 (False).
no_title_text = False
save_vector = False  # Toggle for saving additional vector graphics (pdf)
pub_override = True  # Use fixed colorbar scale for both Europa models
clabel_pad = 5  # Whitespace to add adjacent to contour labels
fig_dpi = 200  # Default dpi to use for raster images when not passed as argument
fmt = 'pdf'  # Default figure format to use (file extension)
animFmt = 'png'  # Figure format to force for animations (PDFs can't be easily made into animation frames)

# Colorbar and tick formatting
cbar_pos = [0.90, 0.18, 0.02, 0.6]
cbar_adj = 5
deft_tsize = 16
deft_ticksize = 14

# PLOT FONTS
# Use "sans-serif" to match captions in Icarus preprints,
# and use "serif" to match captions in final print versions.
# serif selects stix fonts and sans-serif selects computer modern sans-serif (cmss).
font_choice = "serif"

# Set output message level
log = logging.getLogger('MoonMag')
if verbose:
    logLevel = logging.DEBUG
else:
    logLevel = logging.INFO
printFmt = '[%(levelname)s] %(message)s'
stream = logging.StreamHandler()
stream.setFormatter(logging.Formatter(printFmt))
log.setLevel(logLevel)
log.addHandler(stream)
logging.getLogger('matplotlib').setLevel(logging.WARNING)
logging.getLogger('PIL').setLevel(logging.WARNING)

if use_latex and not shutil.which('latex'):
    log.warning('A LaTeX installation was not found. LaTeX will not be used for plot labels.')
    use_latex = False
