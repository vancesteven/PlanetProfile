""" Contains functions for plotting topography and magnetic fields
    from near-spherical conductors.
    Developed in Python 3.8 for "A perturbation method for evaluating the
    magnetic field induced from an arbitrary, asymmetric ocean world 
    analytically" by Styczinski et al.
    DOI: 10.1016/j.icarus.2021.114840
Author: M. J. Styczinski, mjstyczi@uw.edu """

import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdt
import matplotlib.ticker as tick
from matplotlib.gridspec import GridSpec

from PlanetProfile.MagneticInduction.MoonMag import _interior
from PlanetProfile.MagneticInduction.MoonMag.config import *
import PlanetProfile.MagneticInduction.MoonMag.asymmetry_funcs as asym
import PlanetProfile.MagneticInduction.MoonMag.field_xyz as field

J2000 = np.datetime64("2000-01-01T11:58:55.816")

# Set up log messages
printFmt = '[%(levelname)s] %(message)s'
stream = logging.StreamHandler()
stream.setFormatter(logging.Formatter(printFmt))
log = logging.getLogger('MoonMag')
log.setLevel(logging.DEBUG)

mpl.rcParams.update({
    "text.usetex": use_latex,
    "font.family": font_choice,
    "text.latex.preamble": "\\usepackage{stix}"
})

def get_latlon(do_large):
    lenx = ppgc + 1
    leny = int(ppgc/2) + 1
    lon = np.linspace(-180, 180, lenx)
    lat = np.linspace(-90, 90, leny)
    phi = np.radians(lon)
    phi[:int(lenx/2)] = phi[:int(lenx/2)] + 2.0*np.pi
    tht = np.pi/2.0 - np.radians(lat)
    elon = np.linspace(0, 360, lenx)
    phi_std = np.radians(elon)

    if do_large:
        n_latticks = 5
        n_lonticks = 5
        lg_end = "_small"
    else:
        n_latticks = 7
        n_lonticks = 9
        lg_end = ""
    latticks = np.linspace(lat_min, lat_max, n_latticks, endpoint=True, dtype=np.int64)

    # Generate lat/lon formatters
    if do_360:
        phi = phi_std
        lon = elon

        lon_formatter = tick.FuncFormatter(lambda v, pos: east_formatted(v, EAST=True))
        lon_min = 0
        lon_max = 360
        lonticks = np.linspace(lon_min, lon_max, n_lonticks, dtype=np.int64)
    else:
        lon_formatter = tick.FuncFormatter(lambda v, pos: east_formatted(v, EAST=False))
        lon_min = -180
        lon_max = 180
        lonticks = np.linspace(lon_min, lon_max, n_lonticks, dtype=np.int64)

    return lon, lat, lon_min, lon_max, tht, phi, lenx, leny, lonticks, latticks, n_lonticks, n_latticks, lon_formatter, lg_end

lat_formatter = tick.FuncFormatter(lambda v, pos: lat_formatted(v))
con_formatter = tick.FuncFormatter(lambda v, pos:cformat(v))

"""
plotAsym()
    Calculates layer thicknesses gridded by lat,lon.
    Usage: plotAsym(`recalc`, `do_large`, `index=-2`, `cmp_index=-1`, `r_bds=None`, `asym_shape=None`,
                    `pvals=None`, `qvals=None`, `fpath=None`, `bodyname=None`, `descrip="Ice--ocean"`, `no_title=False`)
    Returns: None
    Parameters:
        recalc: boolean. Whether to recalculate gridded values, or attempt to read them in. They are printed to disk if recalc'd.
            All other inputs are optional, because if recalc=True most will not be needed.
        do_large: boolean. Whether to print larger versions of labels and strip off the colorbar, to better display
            in presentation slides or side-by-side figures.
        index: integer (-2). The boundary index to reference in asym_shape, corresponding to the perturbed boundary.
        cmp_index: integer (-1). The boundary index against which to compare the primary boundary index above.
        r_bds: float, shape(N) (None). Mean radius for each boundary in the conductivity model, in m.
        asym_shape: float, shape(n_bds,2,p_max+1,p_max+1) (None). A list of boundary shape parameters chi_pq for each boundary.
        pvals: integer, shape(Npq) (None). A linear list of p values for constructing (p,q) pairs in parallelized loops.
        qvals: integer, shape(Npq) (None). A linear list of q values for constructing (p,q) pairs in parallelized loops.
        fpath: string (None). Optional location to print image file to. Defaults to "figures/".
        bodyname: string (None). Optional body name to include in plot titles and file names.
        append: string(""). Optional string to append to filenames.
        descrip: string ("Ice--ocean"). Description of the asymmetric layer being plotted, for plot titles.
        no_title: boolean (False). If true, print figures without title text.
    """
def plotAsym(recalc, do_large, index=-2, cmp_index=-1, r_bds=None, asym_shape=None, pvals=None, qvals=None, inpath=None,
             fpath=None, outFname=None, bodyname=None, append="", descrip="Ice--ocean", no_title=False, lfig_dpi=None):
    lon, lat, lon_min, lon_max, tht, phi, lenx, leny, lonticks, latticks, n_lonticks, n_latticks, lon_formatter, lg_end = get_latlon(do_large)
    do_cbar = not do_large

    if fpath is None:
        fpath = "figures"
    if lfig_dpi is None:
        lfig_dpi = fig_dpi

    if bodyname is None:
        bfname = ""
        bstr = ""
    else:
        bfname = f"_{bodyname}"
        bstr = f"{bodyname}_"

    if inpath is None:
        inpath = _interior
    if outFname is None:
        outFname = f'asym_devs{bfname}{append}'
    datpath = os.path.join(inpath, f"{outFname}.dat")

    if recalc:
        # Convert from m to km
        r_mean = r_bds[index] / 1e3
        R_cmp = r_bds[cmp_index] / 1e3
        asym_shape_km = asym_shape / 1e3
        mean_thick = abs(R_cmp - r_mean)

        log.debug("Calculating topographic data, this may take some time.")
        # Evaluate the deviations to get r(theta, phi)
        surf = asym.get_rsurf(pvals, qvals, asym_shape_km[index, :, :, :], r_mean, tht, phi)

        thicks = np.abs(R_cmp - surf)

        # Save data to disk for reading in with recalc=False
        fout = open(datpath, "w")
        headinfo = "Each row runs from -90 deg to +90 deg latitude, i.e. lat = np.linspace(-90, 90, leny).\n"
        fout.write(headinfo)
        headlbls = "Layer descrip., leny (#lats), compared radius (km), mean layer thickness (km):\n"
        fout.write(headlbls)
        headnums = f"{descrip}\n{leny}\n{R_cmp}\n{abs(mean_thick)}\n"
        fout.write(headnums)
        colheader = "{:<10}, {:<24}\n".format("Lon (deg)", "Layer thicknesses (km)")
        fout.write(colheader)
        lon_fmt = "{:>10}"
        data_fmt = ", {:>24}"
        for i_lon in range(lenx):
            fout.write( lon_fmt.format(lon[i_lon]) )
            for i_lat in range(leny):
                fout.write(data_fmt.format(thicks[i_lat,i_lon]))
            fout.write( "\n" )
        fout.close()
        log.info(f"Data for asymmetric layer thicknesses written to file: {datpath}")

    else:
        fasym = open(datpath)
        _ = [fasym.readline() for _ in range(2)]
        descrip_old = fasym.readline()
        leny_str = fasym.readline()
        R_cmp_str = fasym.readline()
        mean_thick_str = fasym.readline()
        fasym.close()

        # Allow overwriting in the case of blown up figures
        if not do_large:
            descrip = descrip_old[:-1] # Strip trailing newline
        leny = int(leny_str)
        R_cmp = float(R_cmp_str)
        mean_thick = float(mean_thick_str)
        lat = np.linspace(-90, 90, leny)

        asym_contents = np.loadtxt(datpath, skiprows=7, unpack=False, delimiter=',')
        lon = asym_contents[:, 0]
        thicks = np.transpose(asym_contents[:, 1:])

    log.debug(f"Maximum surface deviation: {np.max(np.abs(thicks-mean_thick)):.3f} km")

    # Generate and format figures
    lon_adj = lon + 0
    if do_360:
        lon_formatter = tick.FuncFormatter(lambda v, pos: east_formatted(v, EAST=True))
        lon_adj[lon_adj==-180] = 180.1
        lon_adj[lon_adj<0] = lon_adj[lon_adj<0] + 360
        iSort = np.argsort(lon_adj)
    else:
        lon_formatter = tick.FuncFormatter(lambda v, pos: east_formatted(v, EAST=False))
        lon_adj[lon_adj==360] = -0.1
        lon_adj[lon_adj>=180] = lon_adj[lon_adj>=180] - 360
        iSort = np.argsort(lon_adj)
        
    lon_adj = lon_adj[iSort]
    thicks_adj = thicks[:,iSort]

    # Stitch endpoints together if we adjusted the plotting coords
    if not (np.max(lon_adj) == 180 or np.max(lon_adj) == 360):
        if do_360:
            lon_adj = np.append(lon_adj, 360)
        else:
            lon_adj = np.append(lon_adj, 180)
        thicks_adj = np.append(thicks_adj, thicks_adj[:,0:1], axis=1)

    if bodyname == "Miranda":
        levels = round(mean_thick) + np.array([-15, -11, -7, -3, 1, 5, 9])*49.810/20.852 # For matching plot from Hemingway and Mittal (2019)
    elif bodyname == "Enceladus":
        levels = round(mean_thick) + np.array([-15, -11, -7, -3, 1, 5, 9])
    else:
        levels = None

    fig, axes = plt.subplots(1, 1, figsize=deft_figsize)
    axes.set_xticks([])
    axes.set_yticks([])
    cbar_title = "Layer thickness $(\mathrm{km})$"
    cbar_top = ""
    if do_cbar:
        cbar_ax = fig.add_axes(cbar_pos)
    else:
        lonticks = np.linspace(lon_min, lon_max, 5, dtype=np.int64)
        latticks = np.linspace(lat_min, lat_max, 5, dtype=np.int64)

    themap = plt.axes()
    fig.subplots_adjust(left=0.07, right=0.88, wspace=0.15, hspace=0.05, top=0.90, bottom=0.07)

    # Apply plot formatting
    themap.set_xticks(lonticks)
    themap.set_yticks(latticks)
    themap.xaxis.set_major_formatter(lon_formatter)
    themap.yaxis.set_major_formatter(lat_formatter)
    if no_title:
        nt_adjust = 0.08
    else:
        nt_adjust = 0.0
    if do_large:
        fig.subplots_adjust(left=0.09, right=0.965, wspace=0.15, hspace=0.05, top=0.875+nt_adjust, bottom=0.1)
        tick_size = deft_ticksize * 10 / 7
        clabel_size = deft_ticksize * 8 / 7
    else:
        fig.subplots_adjust(left=0.07, right=0.89, wspace=0.15, hspace=0.05, top=0.90+nt_adjust, bottom=0.07)
        tick_size = deft_ticksize
        clabel_size = deft_ticksize * 6 / 7
    if do_cbar:
        ptitle = f"{bodyname}: {descrip} boundary topography, $\overline{{D}}={abs(mean_thick):.0f}\,\mathrm{{km}}$"
        tsize = deft_tsize
    else:
        ptitle = f"{descrip} thickness $(\mathrm{{km}})$, $\overline{{D}}={abs(mean_thick):.0f}\,\mathrm{{km}}$"
        tsize = deft_tsize * 1.5
    themap.tick_params(axis='both', which='major', labelsize=tick_size)
    if not no_title:
        fig.suptitle(ptitle, size=tsize)

    # Plot the data
    mesh = plt.pcolormesh(lon_adj, lat, thicks_adj, shading="auto", cmap="PuBu_r", rasterized=True)
    cont = plt.contour(lon_adj, lat, thicks_adj, levels=levels, colors="black")
    lbls = plt.clabel(cont, fmt="%1.0f", fontsize=clabel_size, inline_spacing=clabel_pad)

    # Finish formatting (as long as we're not making big-label plots)
    if do_cbar:
        cbar = fig.colorbar(mesh, cax=cbar_ax, label=cbar_title, ticks=levels, format="%1.0f")
        cbar.ax.set_title(cbar_top, pad=cbar_adj, size=12)

    # Save the figure
    topofig = f"{bstr}asym_contour"
    print_fname = os.path.join(fpath, f"{topofig}{append}{lg_end}")
    fig.savefig(f"{print_fname}.{fmt}", format=fmt, dpi=lfig_dpi,)
    log.info(f"Contour plot for asym bdy saved to: {print_fname}.{fmt}")

    return

#############################################

# For configuring longitudes from -180 to 180 or 0 to 360.
def east_formatted(longitude, EAST=True):
    fmt_string = u'{longitude:{num_format}}{degree}{hemisphere}'
    return fmt_string.format(longitude=abs(longitude), num_format='g',
                            hemisphere=lon_hemisphere(longitude, EAST=EAST),
                            degree=u'\u00B0')
def lon_hemisphere(longitude, EAST=True):
    if EAST:
        longitude = fix_lons(longitude)
    if longitude > 0:
        hemisphere = 'E'
    elif longitude < 0:
        hemisphere = 'W'
    elif longitude == 0:
        hemisphere = ''
    else:
        hemisphere = 'E'
    return hemisphere
def fix_lons(lons):
    fixed_lons = lons[lons!=360] % 360
    return fixed_lons
def cformat(field):
    fmt_string = u'{sign}{field:{num_format}}'
    return fmt_string.format(field=abs(field), sign=get_sign(field), num_format='g')
def get_sign(val):
    if val < 0:
        sign = u'\u2013'
    else:
        sign = ''
    return sign
def lat_formatted(latitude):
    fmt_string = u'{latitude:{num_format}}{degree}{hemisphere}'
    return fmt_string.format(latitude=abs(latitude), num_format='g',
                            hemisphere=lat_hemisphere(latitude),
                            degree=u'\u00B0')
def lat_hemisphere(latitude):
    if latitude == 0:
        hemisphere = ''
    elif latitude > 0:
        hemisphere = 'N'
    else:
        hemisphere = 'S'
    return hemisphere

#############################################

"""
plotMagSurf()
    Plots an equidistant cylindrical projection of magnetic field values over a boundary surface or sphere.
    Set asym_frac to None (or omit the argument) to plot the field on a sphere of radius r_surf_mean * R_P,
    where R_P is the planetary radius. If difference=True, spherically symmetric moments Binm_sph must be passed; then
    differences between the two are plotted over the desired surface.
    Usage: plotMagSurf(`n_peaks`, `Binm`, `nvals`, `mvals`, `do_large`, `Schmidt=False`, `r_surf_mean=1.0`, `asym_frac=None`,
                `pvals=None`, `qvals=None`, `difference=False`, `Binm_sph=None`, `nprmvals=None`, `mprmvals=None`,
                `fpath=None`, `bodyname=None`, `append=""`, `fend=""`, `tstr=""`, `component=None`, `absolute=False`,
                `no_title=False`)
    Returns: None
    Parameters:
        n_peaks: integer. Number of peaks in the frequency spectrum to iterate over.
        Binm: complex, shape(n_peaks,2,n_max+1,n_max+1) OR shape(n_peaks, (n_max+1)**2-1); if Schmidt=True, tuple of (gnm,hnm).
            Induced magnetic moments in nT.
        nvals: integer, shape(Nnm). A linear list of n values for constructing (n,m) pairs in parallelized loops.
        mvals: integer, shape(Nnm). A linear list of m values for constructing (n,m) pairs in parallelized loops.
        do_large: boolean. Whether to print larger versions of labels and strip off the colorbar, to better display
            in presentation slides or side-by-side figures.
        Schmidt: boolean (False). Whether input magnetic moments are in Schmidt semi-normalized form without Condon-Shortley
            phase. If False, moments must be in fully normalized form with the Condon-Shortley phase.
        r_surf_mean: float (1.0). Fraction of body radius for evaluation surface. For example, evaluating
            at the top of a 200 km ionosphere atop a 1000 km radius body has r_surf_mean of 1.2.
        asym_frac: float, shape(2,p_max+1,p_max+1) (None). asym_shape normalized to the desired boundary. The ratio of
            orthonormal harmonic coefficients in km to the mean radius in km that corresponds to r_surf_mean.
        pvals: integer, shape(Npq) (None). A linear list of p values for constructing (p,q) pairs in parallelized loops.
        qvals: integer, shape(Npq) (None). A linear list of q values for constructing (p,q) pairs in parallelized loops.
        difference: boolean (False). Whether to plot differences against the spherically symmetric case. Requires Binm_sph.
        Binm_sph: complex, shape(n_peaks,2,nprm_max+1,nprm_max+1) OR shape(n_peaks, (nprm_max+1)**2-1) (None).
            Optional magnetic moments for spherically symmetric case to compare against. Required if difference = True.
        nprmvals: integer, shape(Nnmprm) (None). A linear list of n' values for constructing (n',m') pairs in parallelized loops. Required if difference = True.
        mprmvals: integer, shape(Nnmprm) (None). A linear list of m' values for constructing (n',m') pairs in parallelized loops. Required if difference = True.
        fpath: string (None). Optional location to print image file to. Defaults to "figures/".
        bodyname: string (None). Optional body name to include in plot titles and file names.
        append: string (""). Optional string to append to filenames.
        fend: string ("").  Frame filename ending to append to keep animation frames in ascending order. Set to
            empty string to use standard labelling; otherwise images will get set up as animation frames.
        tstr: string (""). String to use to describe the frame's timestamp in the image title.
        component: string (None). Optional component to plot instead of magnitude. Accepts 'x', 'y', and 'z'.
        absolute: boolean (False). Optional flag to plot the absolute induced field in addition to differences.
        no_title: boolean (False). If true, print figures without title text.
        outFname: string (None). Optional string with full file path to override other infrastructure for determining output path.
        gnm, hnm: complex, shape(n_max+1,n_max+1). Schmidt semi-normalized magnetic moments. Passed as a tuple in Binm.
        marks: tuple of 2 arrays, each shape(j). Lat/lon pairs (array 1 is lat, array 2 is lon) in degrees to mark on the plot.
    """
def plotMagSurf(n_peaks, Binm, nvals, mvals, do_large, Schmidt=False, r_surf_mean=1.0, asym_frac=None, pvals=None, qvals=None,
                difference=False, Binm_sph=None, nprmvals=None, mprmvals=None, fpath=None, bodyname=None, append="", fend="",
                tstr="", component=None, absolute=False, no_title=False, outFname=None, title=None, sym_title=None,
                lfig_dpi=None, marks=None):
    lon, lat, lon_min, lon_max, tht, phi, lenx, leny, lonticks, latticks, n_lonticks, n_latticks, lon_formatter, lg_end = get_latlon(do_large)
    do_cbar = not do_large

    if Schmidt:
        gnm, hnm = Binm
        if difference:
            gnm_sph, hnm_sph = Binm_sph

    if outFname is None:
        if fpath is None:
            fpath = "figures"
    if lfig_dpi is None:
        lfig_dpi = fig_dpi
    if component is None:
        compstr = ""
        comptitlestr = " magnitude"
        substr = "{\mathrm{mag}"
        field_cmap = "afmhot"
    else:
        compstr = component
        comptitlestr = f" ${component}$ component"
        substr = f"{{{component}"
        field_cmap = "seismic"

    Re_Bx, Re_By, Re_Bz = (np.zeros((leny,lenx)) for _ in range(3))

    if difference:
        diffstr = " difference"
        if do_cbar:
            titlestr = " vs. symmetric"
        else:
            titlestr = diffstr
        fname_diff = f"_{compstr}diff"
        Re_Bx_sym, Re_By_sym, Re_Bz_sym = (np.zeros((leny,lenx)) for _ in range(3))
    else:
        diffstr = ""
        titlestr = ""
        fname_diff = compstr

    if bodyname is None:
        bodyname = "P"
    title_dist = f" at $r={r_surf_mean:.2f}R_{bodyname[0]}$"  # at J2000 epoch"
    titlestr = titlestr + title_dist

    if fend == "":
        log.debug("Creating magnetic plot, this may take some time.")
    fig, axes = plt.subplots(1, 1, figsize=deft_figsize)
    plt.clf()
    cbar_title = f"Magnetic field{diffstr} $(\mathrm{{nT}})$"
    cbar_top = ""
    if do_cbar:
        cbar_ax = fig.add_axes(cbar_pos)
    else:
        if do_360:
            lonticks = [0, 90, 180, 270, 360]
        else:
            lonticks = [-180, -90, 0, 90, 180]
        latticks = [-90, -45, 0, 45, 90]

    themap = plt.axes()
    themap.set_xticks([])
    themap.set_yticks([])

    if asym_frac is None:
        r_th_ph = r_surf_mean
    else:
        log.debug("Getting asymmetric surface shape...")
        r_th_ph = asym.get_rsurf(pvals,qvals,asym_frac, r_surf_mean,tht,phi)

    # Combine the moments from all periods of oscillation to find the instantaneous net induced moments
    if Schmidt:
        gnm_sum = np.sum(gnm, axis=0)
        hnm_sum = np.sum(hnm, axis=0)
        Binm_sum = (gnm_sum, hnm_sum)
        Bx, By, Bz = asym.getMagSurf(nvals,mvals,Binm_sum, r_th_ph,tht,phi, Schmidt=Schmidt)
    else:
        Binm_sum = np.sum(Binm, axis=0)
        Bx, By, Bz = asym.getMagSurf(nvals,mvals,Binm_sum, r_th_ph,tht,phi)
    Re_Bx = Re_Bx + np.real(Bx)
    Re_By = Re_By + np.real(By)
    Re_Bz = Re_Bz + np.real(Bz)

    if difference:
        if Schmidt:
            gnm_sph_sum = np.sum(gnm_sph, axis=0)
            hnm_sph_sum = np.sum(hnm_sph, axis=0)
            Binm_sph_sum = (gnm_sph_sum, hnm_sph_sum)
            Bx_sym, By_sym, Bz_sym = asym.getMagSurf(nprmvals,mprmvals,Binm_sph_sum, r_th_ph,tht,phi, Schmidt=Schmidt)
        else:
            Binm_sph_sum = np.sum(Binm_sph, axis=0)
            Bx_sym, By_sym, Bz_sym = asym.getMagSurf(nprmvals,mprmvals,Binm_sph_sum, r_th_ph,tht,phi)
        Re_Bx_sym = Re_Bx_sym + np.real(Bx_sym)
        Re_By_sym = Re_By_sym + np.real(By_sym)
        Re_Bz_sym = Re_Bz_sym + np.real(Bz_sym)

    Bmag = np.sqrt(Re_Bx**2 + Re_By**2 + Re_Bz**2)

    if difference:
        Bmag_sym = np.sqrt(Re_Bx_sym**2 + Re_By_sym**2 + Re_Bz_sym**2)

        Bx_diff = Re_Bx - Re_Bx_sym
        By_diff = Re_By - Re_By_sym
        Bz_diff = Re_Bz - Re_Bz_sym
        Bmag_diff = Bmag - Bmag_sym

        if component == "x":
            B_plot = Bx_diff
            sym_plot = Re_Bx_sym
            asym_plot = Re_Bx
        elif component == "y":
            B_plot = By_diff
            sym_plot = Re_By_sym
            asym_plot = Re_By
        elif component == "z":
            B_plot = Bz_diff
            sym_plot = Re_Bz_sym
            asym_plot = Re_Bz
        else:
            B_plot = Bmag_diff
            sym_plot = Bmag_sym
            asym_plot = Bmag
    else:
        if component == "x":
            B_plot = Re_Bx
        elif component == "y":
            B_plot = Re_By
        elif component == "z":
            B_plot = Re_Bz
        else:
            B_plot = Bmag

    # Set plot formatting
    themap.set_xticks(lonticks)
    themap.set_yticks(latticks)
    themap.xaxis.set_major_formatter(lon_formatter)
    themap.yaxis.set_major_formatter(lat_formatter)
    if no_title:
        nt_adjust = 0.08
    else:
        nt_adjust = 0.0
    if do_large:
        fig.subplots_adjust(left=0.09, right=0.959, wspace=0.15, hspace=0.05, top=0.875+nt_adjust, bottom=0.1)
        tick_size = deft_ticksize*10/7
        clabel_size = deft_ticksize*8/7
    else:
        fig.subplots_adjust(left=0.07, right=0.89, wspace=0.15, hspace=0.05, top=0.90+nt_adjust, bottom=0.07)
        tick_size = deft_ticksize
        clabel_size = deft_ticksize*6/7

    if do_cbar:
        ptitle = f"{bodyname} induced field{comptitlestr}{titlestr}"
        sym_ptitle = f"{bodyname} symmetric induced magnetic field magnitude"
        asym_ptitle = f"{bodyname} asymmetric induced field{comptitlestr}{title_dist}"
        abs_cbar_title = "Magnetic field ($\mathrm{nT}$)"
        tsize = deft_tsize
    else:
        ptitle = f"Induced field $B_{substr}}}${titlestr} $(\mathrm{{nT}})$"
        sym_ptitle = f"Induced field $B_\mathrm{{sym}}${title_dist} $(\mathrm{{nT}})$"
        asym_ptitle = f"Induced field $B_{substr},\mathrm{{asym}}}}${title_dist} $(\mathrm{{nT}})$"
        tsize = deft_tsize * 1.5

    if fend != "":
        if bodyname == "Miranda" and component is None:
            diff_cmap = "plasma"
        else:
            diff_cmap = "seismic"

        ptitle = f"{ptitle} at t={tstr} h"
        comp_adj = 0
        if bodyname == "Europa":
            if not difference:
                maxval = None
                minval = None
            elif append == "_prev":
                maxval = 2.4
                minval = -maxval
            else:
                maxval = 2.4
                minval = -maxval
        elif bodyname == "Miranda":
            maxval = 0.32
            if not difference:
                maxval = None
                minval = None
            elif component is None:
                if r_surf_mean < 1.2: maxval = 3.2
                minval = 0
            else:
                minval = -maxval
        elif bodyname == "Callisto":
            if not difference:
                maxval = None
                minval = None
            else:
                maxval = 0.08
                minval = -maxval
        elif bodyname == "Triton":
            if not difference:
                maxval = None
                minval = None
            else:
                maxval = 0.48
                minval = -maxval
        elif bodyname == "Enceladus":
            if component is None:
                maxval = 0.45
                minval = 0.05
            else:
                maxval = None
                minval = None
        else:
            maxval = None
            minval = None
        if maxval is not None:
            clevels = np.linspace(minval, maxval, 9+comp_adj, endpoint=True)
        else:
            clevels = None
    else:
        clevels = None
        # Ensure contrast in the event values all have the same sign
        if (B_plot > 0).all():
            diff_cmap = "plasma"
            minval = None
            maxval = None
        elif (B_plot < 0).all():
            diff_cmap = "cividis"
            minval = None
            maxval = None
        else:
            diff_cmap = "seismic"
            if (pub_override and (bodyname == "Europa" and append[:6] == "_Tobie")) and difference:
                maxval = 2.4
                minval = -maxval
                clevels = np.linspace(minval, maxval, 9, endpoint=True)
            else:
                maxval = np.max(np.abs(B_plot))
                minval = -maxval

    themap.tick_params(axis='both', which='major', labelsize=tick_size)
    if not no_title:
        if title is None:
            title = ptitle
        fig.suptitle(title, size=tsize)

    # Generate the plot

    if difference:
        plot_cmap = diff_cmap
    else:
        plot_cmap = field_cmap

    mesh = plt.pcolormesh(lon, lat, B_plot, shading="auto", cmap=plot_cmap, vmin=minval, vmax=maxval, rasterized=True)
    cont = plt.contour(lon, lat, B_plot, levels=clevels, colors="black")
    lbls = plt.clabel(cont, fmt=con_formatter, fontsize=clabel_size, inline_spacing=clabel_pad)
    # Mark points if passed
    if marks is not None:
        marked = plt.scatter(marks[1], marks[0], color=cMark, edgecolors=cMarkE, marker=stMark, s=szMark)

    # Show colorbar if there's room
    if do_cbar:
        cbar = fig.colorbar(mesh, cax=cbar_ax, label=cbar_title)
        cbar.ax.set_title(cbar_top, pad=cbar_adj, size=12)

    # Save the figure
    if bodyname is None or bodyname == "":
        topofig = f"field_asym{fname_diff}{lg_end}{fend}"
        sym_topofig = "field_sym"
        asym_topofig = "field_asym"
    else:
        topofig = f"{bodyname}_field_asym{append}{fname_diff}{lg_end}{fend}"
        sym_topofig = f"{bodyname}_field_sym"
        asym_topofig = f"{bodyname}_field_asym"

    lfmt = fmt + ''
    if outFname is None:
        if fend != "":
            fpath = os.path.join(fpath, "anim_frames")
            lfmt = animFmt

        outFname = os.path.join(fpath, f"{topofig}.{lfmt}")
    fig.savefig(outFname, format=lfmt, dpi=lfig_dpi)
    log.info(f"Contour plot for asym field saved to: {outFname}")

    # Plot the absolute induced field in addition to a difference
    if (absolute and difference) and fend=="":
        for old_cont in cont.collections:
            old_cont.remove()
        for old_lbls in lbls:
            old_lbls.remove()
        del cont
        del lbls

        if (asym_plot > 0).all() or (asym_plot > 0).all():
            minval = np.min(asym_plot)
            maxval = np.max(asym_plot)
        else:
            maxval = np.max(np.abs(asym_plot))
            minval = -maxval

        if component is None:
            abs_cmap = "afmhot"
        else:
            abs_cmap = "seismic"

        mesh.set_array(asym_plot)
        mesh.set_cmap(abs_cmap)
        mesh.set_clim(minval, maxval)
        cont = plt.contour(lon, lat, asym_plot, colors="black")
        lbls = plt.clabel(cont, fmt=con_formatter, fontsize=clabel_size, inline_spacing=clabel_pad)
        if marks is not None:
            marked = plt.scatter(marks[1], marks[0], color=cMark, edgecolors=cMarkE, marker=stMark, s=szMark)

        if do_cbar:
            cbar.set_label(abs_cbar_title)

        if not no_title:
            fig.suptitle(asym_ptitle, size=tsize)
        print_abs_asym_fname = os.path.join(fpath, f"{asym_topofig}{append}{lg_end}")
        fig.savefig(f"{print_abs_asym_fname}.{fmt}", format=fmt, dpi=lfig_dpi)
        log.info(f"Contour plot for absolute asym field saved to: {print_abs_asym_fname}.{fmt}")

        # Plot analogous plot for field from symmetric shape for comparison
        for old_cont in cont.collections:
            old_cont.remove()
        for old_lbls in lbls:
            old_lbls.remove()
        del cont
        del lbls

        if (sym_plot>0).all() or (sym_plot>0).all():
            minval=np.min(sym_plot)
            maxval=np.max(sym_plot)
        else:
            maxval = np.max(np.abs(sym_plot))
            minval = -maxval
        sym_clevels = np.linspace(minval, maxval, 10, endpoint=True)

        mesh.set_array(sym_plot)
        mesh.set_cmap(abs_cmap)
        mesh.set_clim(minval,maxval)
        cont = plt.contour(lon, lat, sym_plot, colors="black")
        lbls = plt.clabel(cont, fmt=con_formatter, fontsize=clabel_size, inline_spacing=clabel_pad)
        if marks is not None:
            marked = plt.scatter(marks[1], marks[0], color=cMark, edgecolors=cMarkE, marker=stMark, s=szMark)

        if do_cbar:
            cbar.set_label(abs_cbar_title)

        if not no_title:
            if fend != "":
                title_tstring = f" at t={tstr} h"
            else:
                title_tstring = ""
            if sym_title is None:
                sym_title = f"{sym_ptitle}{title_tstring}"
            fig.suptitle(sym_title, size=tsize)
        print_abs_sym_fname = os.path.join(fpath, f"{sym_topofig}{append}{lg_end}")
        fig.savefig(f"{print_abs_sym_fname}.{fmt}", format=fmt, dpi=lfig_dpi)
        log.info(f"Contour plot for absolute sym field saved to: {print_abs_sym_fname}.{fmt}")

    plt.close()
    return

#############################################


"""
plotTimeSeries()
    Usage: plotTimeSeries(`loc`, `Binm`, `Benm`, `t_start`, `T_hrs`, `nprm_max`, `n_max`, `nvals`, `mvals`, `n_pts=200`,
                          `component=None`, `Binm_sph=None`, `bodyname=None`, `append=""`, `fpath=None`)
    Returns: None
    Parameters:
        loc: float, shape(4). Location for time series in body-centered, Cartesian coordinates, with units of body radius.
        Binm: complex, shape(2, n_max+1, n_max+1) OR shape(Nnm). Induced moments to plot in the time series.
        Benm: complex, shape(2, nprm_max+1, nprm_max+1). Excitation moments to combine with the induced field to calculate net field.
        t_start: float. Time in seconds past J2000 at which to start the time series.
        T_hrs: float. Period in hours of magnetic oscillations.
        nprm_max: integer. Maximum degree n' in excitation moments.
        n_max: integer. Maximum degree n in induced moments.
        nvals: integer, shape(Nnm). A linear list of n values for constructing (n,m) pairs to evaluate moments.
        mvals: integer, shape(Nnm). A linear list of m values for constructing (n,m) pairs to evaluate moments.
        n_pts: integer (200). Number of points in the time series to plot, i.e. the resolution.
        component: string (None). Component of magnetic field to plot. Options are None (for magnitude), "x", "y", or "z".
        Binm_sph: complex, shape(2, nprm_max+1, nprm_max+1) OR shape(Nnmprm) (None). Optional induced moments for a spherically symmetric body,
            to plot for comparison purposes.
        bodyname: string (None). Name of the body being modeled for titles and filenames.
        append: string (""). Optional string to append to filenames.
        fpath: string (None). Optional path to figure save location. Passing None defaults to "figures/".  
    """
def plotTimeSeries(loc, Binm, Benm, t_start, T_hrs, nprm_max, n_max, nvals, mvals, n_pts=200, component=None, Binm_sph=None, bodyname=None, append="", fpath=None,
                   lfig_dpi=None):
    if fpath is None:
        fpath = "figures"
    if lfig_dpi is None:
        lfig_dpi = fig_dpi

    T_secs = T_hrs * 3600
    x = loc[0]
    y = loc[1]
    z = loc[2]
    r = loc[3]
    t_s = np.arange(t_start, t_start+T_secs, T_secs/n_pts)
    t_h = np.arange(0, T_hrs, T_hrs/n_pts)
    omega = 2*np.pi/T_secs

    if n_max > 10:
        n_max = 10
        log.warning(f'Evaluation of magnetic fields is supported only up to n={n_max}. n_max has been set to {n_max}.')
    Nnm = (n_max + 1) ** 2 - 1
    Nnmprm = (nprm_max + 1) ** 2 - 1

    # Linearize Binm values
    if np.size(np.shape(Binm)) > 2:
        lin_Binm = np.array([ Binm[int(mvals[iN]<0),nvals[iN],abs(mvals[iN])] for iN in range(Nnm) ])
    else:
        lin_Binm = Binm

    Bnet_x, Bnet_y, Bnet_z = (np.zeros(n_pts, dtype=np.complex128) for _ in range(3))
    # Linearize Binm_sph values
    if Binm_sph is not None:
        if np.size(np.shape(Binm_sph)) > 2:
            lin_Binm_sph = np.array([Binm_sph[int(mvals[iN] < 0), nvals[iN], abs(mvals[iN])] for iN in range(Nnmprm)])
        else:
            lin_Binm_sph = Binm_sph
        Bnet_x_sph, Bnet_y_sph, Bnet_z_sph = (np.zeros(n_pts, dtype=np.complex128) for _ in range(3))

    for iN in range(Nnm):
        n = nvals[iN]
        m = mvals[iN]
        if iN < Nnmprm:
            Be_x, Be_y, Be_z = field.eval_Be(n, m, Benm[int(m<0),n,abs(m)], x, y, z, r, omega=omega, t=t_s)
            Bnet_x += Be_x
            Bnet_y += Be_y
            Bnet_z += Be_z

            if Binm_sph is not None:
                Bnet_x_sph += Be_x
                Bnet_y_sph += Be_y
                Bnet_z_sph += Be_z

                Bi_x_sph, Bi_y_sph, Bi_z_sph = field.eval_Bi(n, m, lin_Binm_sph[iN], x, y, z, r, omega=omega, t=t_s)
                Bnet_x_sph += Bi_x_sph
                Bnet_y_sph += Bi_y_sph
                Bnet_z_sph += Bi_z_sph

        Bi_x, Bi_y, Bi_z = field.eval_Bi(n, m, lin_Binm[iN], x, y, z, r, omega=omega, t=t_s)
        Bnet_x += Bi_x
        Bnet_y += Bi_y
        Bnet_z += Bi_z

    fig, axes = plt.subplots(1, 1)

    # Legend labels
    asym_label = "asymmetric"
    sym_label = "symmetric"
    # Set plot labels
    if component is not None:
        substr = component
        compstr = f"${substr}$ component "
        coordstr = ", IAU coordinates"
        if component == "x":
            Bplot = np.real(Bnet_x)
            if Binm_sph is not None:
                Bplot_sph = np.real(Bnet_x_sph)
        elif component == "y":
            Bplot = np.real(Bnet_y)
            if Binm_sph is not None:
                Bplot_sph = np.real(Bnet_y_sph)
        elif component == "z":
            Bplot = np.real(Bnet_z)
            if Binm_sph is not None:
                Bplot_sph = np.real(Bnet_z_sph)
        else:
            raise ValueError(f"ERROR: Selected component is not supported: '{component}'. Please use None (for magnitude) or 'x', 'y', 'z'.")
    else:
        substr = "\mathrm{mag}"
        compstr = ""
        coordstr = ""
        Bplot = np.sqrt(np.real(Bnet_x)**2 + np.real(Bnet_y)**2 + np.real(Bnet_z)**2)
        if Binm_sph is not None:
            Bplot_sph = np.sqrt(np.real(Bnet_x_sph)**2 + np.real(Bnet_y_sph)**2 + np.real(Bnet_z_sph)**2)

    if t_start == 0:
        epstr = " J2000"
    else:
        epstr = ""

    if bodyname == "Europa" or bodyname == "Callisto":
        planetname = "jovian"
    elif bodyname == "Miranda":
        planetname = "uranian"
    elif bodyname == "Triton":
        planetname = "neptunian"
    elif bodyname == "Enceladus":
        planetname = "saturnian"
    else:
        planetname = "planetary"
    axes.set_title(f"{bodyname} net time-varying field {compstr}at sub-{planetname} point")
    axes.set_xlabel(f"Time after{epstr} epoch ($\mathrm{{hr}}$)")
    axes.set_ylabel(f"Net $B_{substr}$ ($\mathrm{{nT}}$){coordstr}")
    axes.grid()

    axes.plot(t_h, Bplot, color=c[0], label=asym_label)
    if Binm_sph is not None:
        axes.plot(t_h, Bplot_sph, color=c[1], label=sym_label)
        plt.legend(loc="best")

    #	Save and close
    fig_fname = os.path.join(fpath, f"{bodyname}_tSeries{append}")
    fig.savefig(f"{fig_fname}.{fmt}", format=fmt, dpi=lfig_dpi)
    plt.close()
    log.info(f"Time series plot saved to file: {fig_fname}.{fmt}")
    return

#############################################


"""
plotTrajec()
    Usage: plotTrajec(`x`,`y`,`z`,`r`,`t`, `Binm`, `Benm`, `peak_omegas`, `nprm_max`, `n_max`, `nvals`, `mvals`, `R_body=None`, 
                      `component=None`, `Binm_sph=None`, `net_field=False`, `bodyname=None`, `append=""`, `fpath=None`)
    Returns: None
    Parameters:
        t: float, shape(n_pts). Linear arrays of corresponding x,y,z, and r values in units of body radius and t values in seconds past J2000.
        Binm: complex, shape(n_peaks, 2, n_max+1, n_max+1) OR shape(Nnm). Induced moments at J2000.
        Benm: complex, shape(n_peaks, 2, nprm_max+1, nprm_max+1). Excitation moments at J2000.
        peak_omegas: float, shape(n_peaks). List of angular frequencies in rad/s for excitation moments.
        nprm_max: integer. Maximum degree n' in excitation moments.
        n_max: integer. Maximum degree n in induced moments.
        nvals: integer, shape(Nnm). A linear list of n values for constructing (n,m) pairs to evaluate moments.
        mvals: integer, shape(Nnm). A linear list of m values for constructing (n,m) pairs to evaluate moments.
        R_body: float (None). Body radius in km. Required if t is an array of a single value (if we're plotting a snapshot cut).
        difference: boolean (False). Whether to plot the difference between the asymmetric and symmetric fields. Requires Binm_sph.
        component: string (None). Component of magnetic field to plot. Options are None (for magnitude), "x", "y", or "z".
        Binm_sph: complex, shape(n_peaks, 2, nprm_max+1, nprm_max+1) OR shape(Nnmprm) (None). Optional induced moments for a spherically symmetric body,
            to plot for comparison purposes. Required if difference == True.
        bodyname: string (None). Name of the body being modeled for titles and filenames.
        net_field: boolean (False). Whether to plot the net magnetic field or just the induced field.
        append: string (""). Optional string to append to filenames.
        fpath: string (None). Optional path to figure save location. Passing None defaults to "figures/".  
    """
def plotTrajec(t, Bx, By, Bz, Bdat=None, bodyname=None, t_CA=None, append="", fpath=None, lfig_dpi=None):
    if fpath is None:
        fpath = "figures"
    if lfig_dpi is None:
        lfig_dpi = fig_dpi

    # Set plot labels
    fig, axes = plt.subplots(3, 1)

    fig.suptitle(bodyname + " net magnetic field, IAU coordinates")
    axes[-1].set_xlabel("Measurement time")
    axes[0].set_ylabel("$B_x (\mathrm{nT})$")
    axes[1].set_ylabel("$B_y (\mathrm{nT})$")
    axes[2].set_ylabel("$B_z (\mathrm{nT})$")
    axes[0].grid()
    axes[1].grid()
    axes[2].grid()
    datefmt = mdt.ConciseDateFormatter(mdt.AutoDateLocator())
    axes[0].xaxis.set_major_formatter(datefmt)
    axes[1].xaxis.set_major_formatter(datefmt)
    axes[2].xaxis.set_major_formatter(datefmt)

    Bxplot = np.real(Bx)
    Byplot = np.real(By)
    Bzplot = np.real(Bz)
    axes[0].plot(t, Bxplot, color=c[0])
    axes[1].plot(t, Byplot, color=c[0])
    axes[2].plot(t, Bzplot, color=c[0])
    # Legend labels
    if Bdat is not None:
        (BxDat, ByDat, BzDat) = Bdat
        axes[0].plot(t, BxDat, color=c[1])
        axes[1].plot(t, ByDat, color=c[1])
        axes[2].plot(t, BzDat, color=c[1])

    if t_CA is not None:
        axes[0].axvline(x=t_CA, color=c[2])
        axes[1].axvline(x=t_CA, color=c[2])
        axes[2].axvline(x=t_CA, color=c[2])
        topYmin, topYmax = axes[0].get_ylim()
        axes[0].text(t_CA, topYmax + (topYmax-topYmin)/20, "CA", ha="center")

    # Save and close
    fig_fname = os.path.join(fpath, f"{bodyname}-{append}")
    fig.savefig(f"{fig_fname}.{fmt}", format=fmt, dpi=lfig_dpi)
    plt.close()
    log.info(f"Trajectory plot saved to file: {fig_fname}.{fmt}")

    return

#############################################


"""
calcAndPlotTrajec()
    Usage: plotTrajec(`x`,`y`,`z`,`r`,`t`, `Binm`, `Benm`, `peak_omegas`, `nprm_max`, `n_max`, `nvals`, `mvals`, `R_body=None`, 
                      `component=None`, `Binm_sph=None`, `net_field=False`, `bodyname=None`, `append=""`, `fpath=None`)
    Returns: None
    Parameters:
        x,y,z,r,t: float, shape(n_pts). Linear arrays of corresponding x,y,z, and r values in units of body radius and t values in seconds past J2000.
        Binm: complex, shape(n_peaks, 2, n_max+1, n_max+1) OR shape(Nnm). Induced moments at J2000.
        Benm: complex, shape(n_peaks, 2, nprm_max+1, nprm_max+1). Excitation moments at J2000.
        peak_omegas: float, shape(n_peaks). List of angular frequencies in rad/s for excitation moments.
        nprm_max: integer. Maximum degree n' in excitation moments.
        n_max: integer. Maximum degree n in induced moments.
        nvals: integer, shape(Nnm). A linear list of n values for constructing (n,m) pairs to evaluate moments.
        mvals: integer, shape(Nnm). A linear list of m values for constructing (n,m) pairs to evaluate moments.
        R_body: float (None). Body radius in km. Required if t is an array of a single value (if we're plotting a snapshot cut).
        difference: boolean (False). Whether to plot the difference between the asymmetric and symmetric fields. Requires Binm_sph.
        component: string (None). Component of magnetic field to plot. Options are None (for magnitude), "x", "y", or "z".
        Binm_sph: complex, shape(n_peaks, 2, nprm_max+1, nprm_max+1) OR shape(Nnmprm) (None). Optional induced moments for a spherically symmetric body,
            to plot for comparison purposes. Required if difference == True.
        bodyname: string (None). Name of the body being modeled for titles and filenames.
        net_field: boolean (False). Whether to plot the net magnetic field or just the induced field.
        append: string (""). Optional string to append to filenames.
        fpath: string (None). Optional path to figure save location. Passing None defaults to "figures/".  
    """
def calcAndPlotTrajec(x,y,z,r,t, Binm, Benm, peak_omegas, nprm_max, n_max, nvals, mvals, R_body=None, difference=False,
               component=None, Binm_sph=None, net_field=False, bodyname=None, append="", fpath=None, lfig_dpi=None):
    if fpath is None:
        fpath = "figures"
    if lfig_dpi is None:
        lfig_dpi = fig_dpi

    if n_max > 10:
        n_max = 10
        log.warning(f'Evaluation of magnetic fields is supported only up to n={n_max}. n_max has been set to {n_max}.')

    Nnm = (n_max + 1) ** 2 - 1
    Nnmprm = (nprm_max + 1) ** 2 - 1
    n_peaks = np.size(peak_omegas)
    n_pts = np.maximum(np.size(t), np.size(r))

    # Linearize Binm values
    if np.size(np.shape(Binm)) > 2:
        lin_Binm = np.zeros((n_peaks,Nnm), dtype=np.complex128)
        for i_om in range(n_peaks):
            lin_Binm[i_om,:] = np.array([ Binm[i_om,int(mvals[iN]<0),nvals[iN],abs(mvals[iN])] for iN in range(Nnm) ])
    else:
        lin_Binm = Binm

    Bnet_x, Bnet_y, Bnet_z = (np.zeros(n_pts, dtype=np.complex128) for _ in range(3))
    if Binm_sph is not None:
        if np.size(np.shape(Binm_sph)) > 2:
            lin_Binm_sph = np.zeros((n_peaks,Nnmprm), dtype=np.complex128)
            for i_om in range(n_peaks):
                lin_Binm_sph[i_om,:] = np.array([ Binm_sph[i_om,int(mvals[iN]<0),nvals[iN],abs(mvals[iN])] for iN in range(Nnmprm) ])
        else:
            lin_Binm_sph = Binm_sph
        Bnet_x_sph, Bnet_y_sph, Bnet_z_sph = (np.zeros(n_pts, dtype=np.complex128) for _ in range(3))

    for i_om in range(n_peaks):
        for iN in range(Nnm):
            n = nvals[iN]
            m = mvals[iN]
            if iN < Nnmprm:
                if net_field and not difference:
                    Be_x, Be_y, Be_z = field.eval_Be(n, m, Benm[i_om,int(m<0),n,abs(m)], x, y, z, r, omega=peak_omegas[i_om], t=t)
                else:
                    Be_x, Be_y, Be_z = ( np.zeros(n_pts, dtype=np.complex128) for _ in range(3) )
                Bnet_x += Be_x
                Bnet_y += Be_y
                Bnet_z += Be_z

                if Binm_sph is not None:
                    Bnet_x_sph += Be_x
                    Bnet_y_sph += Be_y
                    Bnet_z_sph += Be_z

                    Bi_x_sph, Bi_y_sph, Bi_z_sph = field.eval_Bi(n, m, lin_Binm_sph[i_om,iN], x, y, z, r, omega=peak_omegas[i_om], t=t)
                    Bnet_x_sph += Bi_x_sph
                    Bnet_y_sph += Bi_y_sph
                    Bnet_z_sph += Bi_z_sph

            Bi_x, Bi_y, Bi_z = field.eval_Bi(n, m, lin_Binm[i_om,iN], x, y, z, r, omega=peak_omegas[i_om], t=t)
            Bnet_x += Bi_x
            Bnet_y += Bi_y
            Bnet_z += Bi_z

    # Legend labels
    asym_label = "asymmetric"
    sym_label = "symmetric"
    # Set plot labels
    if component is not None:
        if component == "all":
            fig, axes = plt.subplots(3, 1)
            Bxplot = np.real(Bnet_x)
            Byplot = np.real(Bnet_y)
            Bzplot = np.real(Bnet_z)
        else:
            fig, axes = plt.subplots(1, 1)
            substr = component
            compstr = "$" + substr + "$ component "
            coordstr = ", IAU coordinates"
            if component == "x":
                Bplot = np.real(Bnet_x)
                if Binm_sph is not None:
                    Bplot_sph = np.real(Bnet_x_sph)
            elif component == "y":
                Bplot = np.real(Bnet_y)
                if Binm_sph is not None:
                    Bplot_sph = np.real(Bnet_y_sph)
            elif component == "z":
                Bplot = np.real(Bnet_z)
                if Binm_sph is not None:
                    Bplot_sph = np.real(Bnet_z_sph)
            else:
                raise ValueError(f"ERROR: Selected component is not supported: '{component}'. Please use None (for magnitude) or 'x', 'y', 'z'.")
    else:
        component = "mag"
        fig, axes = plt.subplots(1, 1)
        substr = "\mathrm{mag}"
        compstr = ""
        coordstr = ""
        Bplot = np.sqrt(np.real(Bnet_x)**2 + np.real(Bnet_y)**2 + np.real(Bnet_z)**2)
        if Binm_sph is not None:
            Bplot_sph = np.sqrt(np.real(Bnet_x_sph)**2 + np.real(Bnet_y_sph)**2 + np.real(Bnet_z_sph)**2)

    if t[0] == 0:
        epstr = " J2000"
    else:
        epstr = ""

    if t[0] == t[-1]:
        plotx = (r - 1) * R_body
        xtitle = "Distance from surface ($\mathrm{km}$)"
        if vert_cut_lat >= 0:
            lat_ltr = "N"
        else:
            lat_ltr = "S"
        if vert_cut_lon >= 0:
            lon_ltr = "E"
        else:
            lon_ltr = "W"
        trajstr = f" above ${abs(vert_cut_lat)}^\circ\mathrm{{{lat_ltr}}}$, ${abs(vert_cut_lon)}^\circ\mathrm{{{lon_ltr}}}$"
        endstr = f", {vert_cut_hr}$\,\mathrm{{hr}}$ past J2000"
    else:
        plotx = t/3600
        xtitle = "Hours past J2000"
        trajstr = ""
        endstr = ""


    if net_field and not difference:
        netstr = " net time-varying field "
        ytitle = f"Net $B_{substr}$ ($\mathrm{{nT}}$){coordstr}"
    elif difference:
        netstr = " induced field difference "
        ytitle = f"$B_{substr}$ difference ($\mathrm{{nT}}$){coordstr}"
        Bplot -= Bplot_sph
    elif component != "all":
        netstr = " induced field "
        ytitle = f"Induced $B_{substr}$ ($\mathrm{{nT}}$){coordstr}"

    if component == "all":
        fig.suptitle(f"{bodyname} net time-varying field, ${bodyname[0]}\phi\Omega$ coordinates")
        axes[-1].set_xlabel(xtitle)
        axes[0].set_ylabel("$B_x (\mathrm{nT})$")
        axes[1].set_ylabel("$B_y (\mathrm{nT})$")
        axes[2].set_ylabel("$B_z (\mathrm{nT})$")
        axes[0].grid()
        axes[1].grid()
        axes[2].grid()

        axes[0].plot(plotx, Bxplot, color=c[0])
        axes[1].plot(plotx, Byplot, color=c[0])
        axes[2].plot(plotx, Bzplot, color=c[0])
    else:
        axes.set_title(f"{bodyname}{netstr}{trajstr}{endstr}")
        axes.set_xlabel(xtitle)
        axes.set_ylabel(ytitle)
        axes.grid()

        axes.plot(plotx, Bplot, color=c[0], label=asym_label)
        if Binm_sph is not None and not difference:
            axes.plot(plotx, Bplot_sph, color=c[1], label=sym_label)
            plt.legend(loc="best")

    # Save and close
    fig_fname = os.path.join(fpath, f"{bodyname}{append}")
    fig.savefig(f"{fig_fname}.{fmt}", format=fmt, dpi=lfig_dpi)
    plt.close()
    log.info(f"Trajectory plot saved to file: {fig_fname}.{fmt}")

    return

#############################################


"""
plotAfunctions()
    Plots kr spectra of complex response amplitudes Ae and At for a specified n, and Ad for n=1.
    This function exists primarily for debug purposes, to verify asymptotic dependence 
    of these important quantities in the case of a 1-layer body.
    Usage: plotAfunctions(`kr`, `n`, `Ae_mag`, `Ae_arg`, `At_mag`, `At_arg`, `Ad_mag`, `Ad_arg`, `AtAd_mag`, `AtAd_arg`)
    Returns: None
    Parameters:
        kr: mpc, shape(n_omegas). List of kr values against which to plot.
        n: integer. n value for which A functions have been calculated.
        Aes: mpc, shape(n_omegas). "Excitation" amplitude A_n^e.
        Ats: mpc, shape(n_omegas). "Tangential" amplitude A_n^t.
        Ads: mpc, shape(n_omegas). "Mixing" amplitude A_n'^\star.
        AtAds: mpc, shape(n_omegas). Ats and Ads multiplied together, which should be asymptotic to (1+0i).
    """
def plotAfunctions(kr, n, Ae_mag, Ae_arg, At_mag, At_arg, Ad_mag, Ad_arg, AtAd_mag, AtAd_arg, lfig_dpi=None):
    if lfig_dpi is None:
        lfig_dpi = fig_dpi

    c = ["black", "blue", "brown", "green"]
    style1 = "solid"
    style2 = "dashed"

    #	(This is Figure 1 in the paper)
    fig, axes = plt.subplots(1, 1, figsize=(5.5, 4))

    #	Set plot labels
    axes.set_title("1-layer response amplitudes $\mathcal{A}_n$")
    axes.set_xlabel("Ice--ocean boundary $|kr|$")
    axes.set_ylabel("Complex magnitude")

    phax = axes.secondary_yaxis('right', functions=(getphase, getmag))
    phax.set_ylabel("Phase delay (degrees)")

    nstr = str(n)
    AeM_label = "$|\mathcal{A}_" + nstr + "^e|$"
    AeA_label = "$-\mathrm{arg}(\mathcal{A}_" + nstr + "^e)$"
    AtM_label = "$|\mathcal{A}_" + nstr + "^t|$"
    AtA_label = "$-\mathrm{arg}(\mathcal{A}_" + nstr + "^t)$"
    AdM_label = "$|\mathcal{A}_1^\star|$"
    AdA_label = "$-\mathrm{arg}(\mathcal{A}_1^\star)$"
    AtAdM_label = "$|\mathcal{A}_1^t\mathcal{A}_1^\star|$"
    AtAdA_label = "$-\mathrm{arg}(\mathcal{A}_1^t\mathcal{A}_1^\star)$"

    piov2 = np.pi / 2.0
    Ae_arg = [ val / piov2 for val in Ae_arg ]
    At_arg = [ val / piov2 for val in At_arg ]
    Ad_arg = [ val / piov2 for val in Ad_arg ]
    AtAd_arg = [ val / piov2 for val in AtAd_arg ]

    axes.plot(kr, Ae_mag, color=c[0], linestyle=style1, label=AeM_label)
    axes.plot(kr, Ae_arg, color=c[0], linestyle=style2, label=AeA_label)
    axes.plot(kr, At_mag, color=c[1], linestyle=style1, label=AtM_label)
    axes.plot(kr, At_arg, color=c[1], linestyle=style2, label=AtA_label)
    axes.plot(kr, Ad_mag, color=c[2], linestyle=style1, label=AdM_label)
    axes.plot(kr, Ad_arg, color=c[2], linestyle=style2, label=AdA_label)
    axes.plot(kr, AtAd_mag, color=c[3], linestyle=style1, label=AtAdM_label)
    axes.plot(kr, AtAd_arg, color=c[3], linestyle=style2, label=AtAdA_label)

    #	Set bounds on axes:
    ylim = 1.02

    axes.set_xlim([np.min(kr), np.max(kr)])
    axes.set_ylim([0.0, ylim])

    #	Adjust decoration
    axes.grid()
    axes.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    #	Save and close
    plt.legend(loc="best")
    thefig = os.path.join("figures", f"A_functions.{fmt}")
    fig.savefig(thefig, format=fmt, dpi=lfig_dpi)
    plt.close()
    log.info(f"A functions plot printed to: {thefig}")

    return


def plotAes(T_h, Ae_mag, Ae_arg, Tinterest_h, BeRatio, intModels, n=1, fEnd="", lfig_dpi=None):
    if lfig_dpi is None:
        lfig_dpi = fig_dpi

    c = ["blue", "green", "brown", "black", "red"]
    style1 = "solid"
    style2 = "dashed"

    fig, axes = plt.subplots(1, 1, figsize=(5.5, 4))

    # Set plot labels
    axes.set_title(f"Complex response amplitudes $\mathcal{{A}}^e_{n}$")
    axes.set_xlabel("Period $(\mathrm{hr})$")
    axes.set_ylabel("Response amplitude $A$")

    plt.rcParams["text.latex.preamble"] = r"\usepackage{stix}\usepackage{upgreek}"
    AeM_label = f"$A$"
    AeA_label = "$\\upphi$"
    Tin_label = "$T_\mathrm{exc}$"

    phax = axes.secondary_yaxis('right', functions=(getphase, getmag))
    phax.set_ylabel("Phase delay $\\upphi$ (degrees)")

    lThick = 1.0
    for i,modelID in enumerate(intModels):
        axes.plot(T_h,  Ae_mag[modelID],         color=c[i%4], linestyle=style1, label=f"{AeM_label} {modelID}", linewidth=lThick)
        axes.plot(T_h, -Ae_arg[modelID]*2/np.pi, color=c[i%4], linestyle=style2, label=f"{AeA_label} {modelID}", linewidth=lThick)

    vThick = 0.5
    axes.axvline(Tinterest_h[0], color=c[4], alpha=BeRatio[0], zorder=-1, label=Tin_label, linewidth=vThick)
    for i in range(1, np.size(Tinterest_h)):
        axes.axvline(Tinterest_h[i], color=c[4], alpha=BeRatio[i], zorder=-1, linewidth=vThick)

    # Set bounds on axes:
    ylim = 1.02

    axes.set_xscale('log')
    axes.set_xlim([np.min(T_h), np.max(T_h)])
    axes.set_ylim([0.0, ylim])

    # Adjust decoration
    axes.grid()
    axes.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    # Save and close
    plt.legend(loc="upper right")
    thefig = os.path.join("figures", f"A1e{fEnd}.{fmt}")
    fig.savefig(thefig, format=fmt, dpi=lfig_dpi)
    plt.close()
    log.info(f"Ae plot printed to: {thefig}")

    return

#	For secondary axis labels
def getphase(A):
    return A*90.0
def getmag(phi):
    return phi/90.0

def SetMap(ax):
    from PlanetProfile.GetConfig import FigMisc

    ax.set_xticks(FigMisc.lonMapTicks_deg)
    ax.set_yticks(FigMisc.latMapTicks_deg)
    ax.tick_params(axis='both', which='major', labelsize=FigMisc.latlonSize)
    ax.xaxis.set_major_formatter(tick.FuncFormatter(FigMisc.LonMapFormatter))
    ax.yaxis.set_major_formatter(tick.FuncFormatter(FigMisc.LatMapFormatter))

    return

def healpixMap(nside, map, title, fName, levels=None, cmap=None, FIELD=True):
    """ This function requires additional imports and is currently intended only for testing
        and validation purposes. """
    from healpy.projector import CartesianProj
    from healpy.pixelfunc import vec2pix
    from PlanetProfile.GetConfig import FigLbl, FigSize, FigMisc

    lonMap_deg, latMap_deg, lon_min, lon_max, _, _, nLonMap, nLatMap, _, _, _, _, _, _ = get_latlon(False)
    cp = CartesianProj(flipconv='geo')
    cp.set_proj_plane_info(xsize=nLonMap, ysize=nLatMap,
        lonra=np.array([lon_min, lon_max]), latra=None)
    vecFunc = lambda vec: vec2pix(nside, vec[0], vec[1], vec[2])
    pix = vecFunc(cp.xy2vec(cp.ij2xy()))
    plotData = map[pix]

    fig = plt.figure(figsize=FigSize.asym)
    grid = GridSpec(1, 1)
    ax = fig.add_subplot(grid[0, 0])
    SetMap(ax)
    asymMap = ax.pcolormesh(lonMap_deg, latMap_deg, plotData,
                            shading='auto', cmap=cmap, rasterized=FigMisc.PT_RASTER)
    asymContours = ax.contour(lonMap_deg, latMap_deg, plotData,
                              levels=levels, colors='black')
    ax.clabel(asymContours, fmt=tick.FuncFormatter(FigMisc.Cformat),
              fontsize=FigMisc.cLabelSize, inline_spacing=FigMisc.cLabelPad)
    ax.set_title(title, size=FigMisc.mapTitleSize)
    ax.set_aspect(1)
    cax = ax.inset_axes([1 + 0.02, 0, 0.03, 1])
    cbarLbl = FigLbl.MagSurfCbarDiffLabel if FIELD else FigLbl.asymCbarLabel
    fig.colorbar(asymMap, ax=ax, cax=cax, label=cbarLbl)
    FigMisc.SetLatex()  # Shouldn't be necessary but rcParams seem to not be getting set
    plt.tight_layout()
    fig.savefig(fName, bbox_inches='tight', format='pdf', dpi=FigMisc.dpi, metadata=FigLbl.meta)
    log.debug(f'Map figure saved to file: {fName}.pdf')
    plt.close()
