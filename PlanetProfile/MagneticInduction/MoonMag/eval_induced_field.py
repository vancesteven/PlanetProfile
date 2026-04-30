""" Calculations for induced magnetic fields
    from near-spherical conductors and plots the results.
    Developed in Python 3.8 for "A perturbation method for evaluating the
    magnetic field induced from an arbitrary, asymmetric ocean world 
    analytically" by Styczinski et al.
    DOI: 10.1016/j.icarus.2021.114840
Author: M. J. Styczinski, mjstyczi@uw.edu """

import os
import numpy as np
from collections.abc import Iterable

import spiceypy as spice

from MoonMag import _interior, _induced, _excitation
from PlanetProfile.MagneticInduction.MoonMag.config import *
import PlanetProfile.MagneticInduction.MoonMag.symmetry_funcs as sym
import PlanetProfile.MagneticInduction.MoonMag.asymmetry_funcs as asym
import PlanetProfile.MagneticInduction.MoonMag.plotting_funcs as plots

def run_calcs(bname, comp, recalc, plot_field, plot_asym, synodic_only=False,
              inp_model_path=None, inp_Bi_path=None, inp_Be_path=None, eval_r=None,
              timePair=None, CAmarks=None, modelOpts=None):

    if eval_r is None:
        eval_r = eval_radius

    if timePair is not None:
        eval_datetime = timePair[1]
    else:
        eval_datetime = default_eval_datetime

    if modelOpts is None:
        modelOpts = []

    flyby_opt = ""  # Label for flyby-specific models

    absolute = (comp==None)   # Whether to plot the symmetric and asymmetric absolute magnitudes in addition to the desired component
    if comp is None:
        compstr = "_mag"
    else:
        compstr = f"_{comp}"

    if bname is None:
        bfname = ""
    else:
        bfname = f"_{bname}"
        bin_name = f"{bname}_"

    # IO file paths
    if inp_model_path is None:
        inp_model_path = _interior
    if inp_Bi_path is None:
        inp_Bi_path = _induced

    if 'BeDir' in modelOpts:
        inp_Be_dir = inp_model_path
    else:
        inp_Be_dir = None

    # p_max is highest degree in boundary shapes to use
    # bname_opt is a string to append to filenames to identify special cases
    bname_opt = ""
    p_max_main = 2
    Benm_model = None
    if bname == "Enceladus":
        p_max_main = 8
        eval_r = 1.1  # About 25 km altitude
    elif bname == "Miranda":
        p_max_main = 8
    elif bname == "Europa":
        if prevEuropa in modelOpts:
            p_max_main = 2
            bname_opt = '_prev'
            eval_r = 1.0
            synodic_only = True
        elif TobieHigh in modelOpts or TobieLow in modelOpts:
            p_max_main = 4
            if TobieHigh in modelOpts:
                bname_opt = '_Tobie_high'
            else:
                bname_opt = '_Tobie_low'
            eval_r = 1.016  # 25 km altitude, the planned CA for some Clipper flybys
    elif bname == "Callisto":
        p_max_main = 2
        if CochraneC3C9 in modelOpts or CochraneAll in modelOpts:
            Benm_model = 'Cochrane2024'
            if timePair is not None:
                flyby_opt = timePair[0]
            if CochraneC3C9 in modelOpts:
                bname_opt = f'_{CochraneC3C9}'
            else:
                bname_opt = f'_{CochraneAll}'

    if 'surface' in modelOpts:
        eval_r = 1.0

    # Highest degree of induced moments
    n_max_main = nprm_max_main + p_max_main

    # Set time to evaluate magnetic fields
    spiceTLS = 'naif0012.tls'
    spice.furnsh(os.path.join(_excitation, spiceTLS))

    #tsec = spice.str2et(dtime.now().isoformat())  # Evaluate as the bodies are right NOW
    #tsec = 0  # Evaluate AT J2000
    tsec = spice.str2et(eval_datetime)

    # Linear arrays of values to loop over for nprm,mprm,p,q,n,m
    nprmvals = [ nprm for nprm in range(1,nprm_max_main+1) for _    in range(-nprm,nprm+1) ]
    mprmvals = [ mprm for nprm in range(1,nprm_max_main+1) for mprm in range(-nprm,nprm+1) ]
    pvals = [ p for p in range(1,p_max_main+1) for _ in range(-p,p+1) ]
    qvals = [ q for p in range(1,p_max_main+1) for q in range(-p,p+1) ]
    nvals = [ n for n in range(1,n_max_main+1) for _ in range(-n,n+1) ]
    mvals = [ m for n in range(1,n_max_main+1) for m in range(-n,n+1) ]
    Nnmprm = len(nprmvals)
    Npq = len(pvals)
    Nnm = len(nvals)

    # Define scaling for specifying distances in terms of body radii
    # Note: r_io MUST be negative
    if bname is None:
        rscale_moments = 1.0
        rscale_asym = 1.0
        r_io = -3
    elif bname == "Callisto":
        rscale_moments = 1/2410.3/1e3
        rscale_asym = 1/2510300
        r_io = -1
    elif bname == "Ganymede":
        rscale_moments = 1/2634.1/1e3
        rscale_asym = 1/2484100
        r_io = -1
    elif bname == "Enceladus":
        rscale_moments = 1/252.1/1e3
        rscale_asym = 1/230035.1
        r_io = -2
    elif bname == "Europa":
        r_io = -2
        if prevEuropa in modelOpts:
            rscale_moments = 1/1560.0/1e3
            rscale_asym = 1/1537500
        else:
            rscale_moments = 1/1561.0/1e3
            if TobieHigh in modelOpts or TobieLow in modelOpts:
                if TobieHigh in modelOpts:
                    rscale_asym = 1/1538503.78421254
                else:
                    rscale_asym = 1/1538489.72081887
            else:
                rscale_asym = 1/(1561 - 30)/1e3
    elif bname == "Miranda":
        rscale_moments = 1/235.8/1e3
        rscale_asym = 1/185989.764241229
        if "surface" in modelOpts:
            r_io = -2
        else:
            r_io = -3
    elif bname == "Triton":
        rscale_moments = 1/1353.4/1e3
        rscale_asym = 1/1803400
        r_io = -1

    R = 1 / rscale_moments / 1e3

    if debug:
        log.debug("Getting Xi/Xid values")
        Xid = asym.get_all_Xid(nprm_max_main, p_max_main, n_max_main, nvals, mvals)
        log.debug("#######   Xid values   #######")
        asym.print_Xid_table(Xid, nprm_max_main, p_max_main, n_max_main)
        log.debug(f"\n#######   Xi values   #######")
        asym.print_Xi_table(nprm_max_main, p_max_main, n_max_main)

    if recalc:
        int_model = os.path.join(inp_model_path, f"interior_model_asym{bfname}{bname_opt}.txt")
        log.debug(f"Using interior model: {int_model}")

        r_bds, sigmas, bcdev = np.loadtxt(int_model, skiprows=1, unpack=True, delimiter=',')
        n_bds = np.size(r_bds)

        # Read in asymmetric shape information
        eps_scaled = r_bds * bcdev
        if relative:
            single_asym = None
        else:
            single_asym = r_io
        asym_shape, grav_shape = asym.read_shape(n_bds, p_max_main, rscale_asym, bodyname=bname, relative=relative, single_asym=single_asym, fpath=inp_model_path,
                                     eps_scaled=eps_scaled, r_bds=r_bds, r_io=r_io, append=bname_opt+flyby_opt, convert_depth_to_chipq=convert_depth_to_chipq)
        r_bds, sigmas, asym_shape = asym.validate(r_bds, sigmas, bcdev, asym_shape, p_max_main)

        # Read in Benm info
        if plot_field:
            peak_periods, Benm, B0, exc_names, _ = asym.read_Benm(nprm_max_main, p_max_main, bodyname=bname, synodic=synodic_only,
                                                               orbital=False, model=Benm_model, fName=inp_Be_path, fpath=inp_Be_dir)
            peak_omegas = 2*np.pi/(peak_periods*3600)
            if not isinstance(peak_omegas, Iterable):
                peak_periods = [peak_periods]
                peak_omegas = [peak_omegas]
            n_peaks = len(peak_omegas)

    else:
        r_bds = None
        asym_shape = None
        grav_shape = None

    if debug:
        fpath = os.path.join("interior", "Avals.dat")
        if recalc:
            Xid = None
            rscaling = None
            omegas = np.linspace(1e-3, 100, 250)  # These will actually be used as kr values
            Binm, Aes, Ats, Ads, krvals = asym.BiList(r_bds, sigmas, omegas, asym_shape, grav_shape, Benm, rscale_moments, nvals, mvals, p_max_main, writeout=False, nprm_max=nprm_max_main, bodyname=bname, debug=True)

            fout = open(fpath, "w")
            header = "{:<24}, {:<24}, {:<24}, {:<24}, {:<24}, {:<24}, {:<24}, {:<24}, {:<24}\n".format("|kr|,", "|A_1^e|,", "-arg(A_1^e),", "|A_1^t|,", "-arg(A_1^t),", "|A_1^\star|", "-arg(A_1^\star)", "|A_1^tA_1^\star|", "-arg(A_1^tA_1^\star)")
            fout.write(header)

            AtAd = Ats[:, 0]*Ads[:, 0]
            kr = [np.abs(krval) for krval in krvals]
            Ae_mag = [np.abs(val) for val in Aes[:, 0]]
            Ae_arg = [-np.angle(val) for val in Aes[:, 0]]
            At_mag = [np.abs(val) for val in Ats[:, 0]]
            At_arg = [-np.angle(val) for val in Ats[:, 0]]
            Ad_mag = [np.abs(val) for val in Ads[:, 0]]
            Ad_arg = [-np.angle(val) for val in Ads[:, 0]]
            AtAd_mag = [np.abs(val) for val in AtAd]
            AtAd_arg = [-np.angle(val) for val in AtAd]

            for i in range(len(omegas)):
                this_line = f"{kr[i]}, {Ae_mag[i]}, {Ae_arg[i]}, {At_mag[i]}, {At_arg[i]}, {Ad_mag[i]}, {Ad_arg[i]}, {AtAd_mag[i]}, {AtAd_arg[i]}\n"
                fout.write(this_line)
            fout.close()
            log.info(f"Data for A functions written to file: {fpath}")
        else:
            kr, Ae_mag, Ae_arg, At_mag, At_arg, Ad_mag, Ad_arg, AtAd_mag, AtAd_arg = np.loadtxt(fpath, skiprows=1, unpack=True, delimiter=',')

        plots.plotAfunctions(kr, 1, Ae_mag, Ae_arg, At_mag, At_arg, Ad_mag, Ad_arg, AtAd_mag, AtAd_arg)
        log.debug("Quitting.")
        quit()

    if plot_asym:
        if bname == "Triton" or bname == "Callisto":
            R_surface = -2
            r_io = -1
            descrip = "Ionosphere"
        else:
            R_surface = r_io + 1
            if DO_LARGE in modelOpts:
                descrip = "Ice shell"
            else:
                descrip = "Ice--ocean"
        plots.plotAsym(recalc, DO_LARGE in modelOpts, index=r_io, cmp_index=R_surface, r_bds=r_bds,
                       asym_shape=asym_shape, pvals=pvals, qvals=qvals, bodyname=bname,
                       append=bname_opt+flyby_opt, descrip=descrip, no_title=no_title_text)

    if plot_field:
        if recalc:
            # Calculate and print to data files
            Binm_sph = sym.BiList(r_bds, sigmas, peak_omegas, Benm, nprmvals, mprmvals, rscale_moments, n_max=nprm_max_main, bodyname=bname, append=bname_opt+flyby_opt, Schmidt=output_Schmidt)
            Binm = asym.BiList(r_bds, sigmas, peak_omegas, asym_shape, grav_shape, Benm, rscale_moments, nvals, mvals, p_max_main, nprm_max=nprm_max_main, bodyname=bname, append=bname_opt+flyby_opt, Schmidt=output_Schmidt)
            n_peaks = len(peak_omegas)
            if output_Schmidt:
                nprmvals = [nprm for nprm in range(1, nprm_max_main + 1) for _ in range(0, nprm + 1)]
                mprmvals = [mprm for nprm in range(1, nprm_max_main + 1) for mprm in range(0, nprm + 1)]
                nvals = [n for n in range(1, n_max_main + 1) for _ in range(0, n + 1)]
                mvals = [m for n in range(1, n_max_main + 1) for m in range(0, n + 1)]
                gnm, hnm = ( np.zeros((n_peaks,n_max_main+1,n_max_main+1), dtype=np.complex128) for _ in range(2) )
                gnm_sph, hnm_sph = ( np.zeros((n_peaks,nprm_max_main+1,nprm_max_main+1), dtype=np.complex128) for _ in range(2) )
                for i in range(n_peaks):
                    gnm[i,...], hnm[i,...] = asym.get_gh_from_Binm(n_max_main,Binm[i,...])
                    gnm_sph[i,...], hnm_sph[i,...] = asym.get_gh_from_Binm(nprm_max_main,Binm_sph[i,...])
        else:
            if output_Schmidt:
                nprmvals = [nprm for nprm in range(1, nprm_max_main + 1) for _ in range(0, nprm + 1)]
                mprmvals = [mprm for nprm in range(1, nprm_max_main + 1) for mprm in range(0, nprm + 1)]
                nvals = [n for n in range(1, n_max_main + 1) for _ in range(0, n + 1)]
                mvals = [m for n in range(1, n_max_main + 1) for m in range(0, n + 1)]
                T_hrs, n_asy, m_asy, lin_gnm_Re, lin_gnm_Im, lin_hnm_Re, lin_hnm_Im = np.loadtxt(
                    os.path.join(inp_Bi_path, f"{bin_name}ghnm_asym{bname_opt}{flyby_opt}.dat"),
                    skiprows=1, unpack=True, delimiter=',')
                T_hrs_sym, n_sph, m_sph, lin_gnm_sph_Re, lin_gnm_sph_Im, lin_hnm_sph_Re, lin_hnm_sph_Im = np.loadtxt(
                    os.path.join(inp_Bi_path, f"{bin_name}ghnm_sym{bname_opt}{flyby_opt}.dat"),
                    skiprows=1, unpack=True, delimiter=',')
            else:
                T_hrs, n_asy, m_asy, lin_Binm_Re, lin_Binm_Im = np.loadtxt(
                    os.path.join(inp_Bi_path, f"{bin_name}Binm_asym{bname_opt}{flyby_opt}.dat"),
                    skiprows=1, unpack=True, delimiter=',')
                T_hrs_sym, n_sph, m_sph, lin_Binm_sph_Re, lin_Binm_sph_Im = np.loadtxt(
                    os.path.join(inp_Bi_path, f"{bin_name}Binm_sym{bname_opt}{flyby_opt}.dat"),
                    skiprows=1, unpack=True, delimiter=',')

            peak_periods = np.unique(T_hrs)
            peak_omegas = 2*np.pi/(peak_periods*3600)
            n_peaks = len(peak_omegas)
            Nnm_asy = int(len(n_asy)/n_peaks)
            peak_periods_sym = np.unique(T_hrs_sym)

            try:
                same_spectrum = (np.equal(peak_periods, peak_periods_sym)).all()
            except:
                raise ValueError("Periods of excitation for Binm_sph and Binm did not match. Set recalc=True and run again to fix the problem.")
            if not same_spectrum:
                raise ValueError("Periods of excitation for Binm_sph and Binm did not match. Set recalc=True and run again to fix the problem.")
            Nnm_sph = int(len(n_sph)/n_peaks)

            if output_Schmidt:
                lin_gnm = lin_gnm_Re + 1j*lin_gnm_Im
                gnm = np.reshape(lin_gnm, (n_peaks, Nnm_asy))
                lin_hnm = lin_hnm_Re + 1j*lin_hnm_Im
                hnm = np.reshape(lin_hnm, (n_peaks, Nnm_asy))
                lin_gnm_sph = lin_gnm_sph_Re + 1j*lin_gnm_sph_Im
                gnm_sph = np.reshape(lin_gnm_sph, (n_peaks,Nnm_sph))
                lin_hnm_sph = lin_hnm_sph_Re + 1j*lin_hnm_sph_Im
                hnm_sph = np.reshape(lin_hnm_sph, (n_peaks,Nnm_sph))
            else:
                lin_Binm = lin_Binm_Re + 1j*lin_Binm_Im
                Binm = np.reshape(lin_Binm, (n_peaks,Nnm_asy))
                lin_Binm_sph = lin_Binm_sph_Re + 1j*lin_Binm_sph_Im
                Binm_sph = np.reshape(lin_Binm_sph, (n_peaks,Nnm_sph))

        plot_on_sphere = True
        if plot_on_sphere:
            asym_frac = None
        else:
            asym_frac = asym_shape[r_io,...] / r_bds[r_io]

        if output_Schmidt:
            gnm_sph_rot = gnm_sph * 1.0
            gnm_rot = gnm * 1.0
            hnm_sph_rot = hnm_sph * 1.0
            hnm_rot = hnm * 1.0
        else:
            Binm_sph_rot = Binm_sph * 1.0
            Binm_rot = Binm * 1.0

        if DO_GIF in modelOpts:
            log.debug(f"Making {n_frames} animation frames.")
            for iT in range(n_frames):
                iT_str = f"{iT:04}"
                t_hr = peak_periods[-1] * iT / n_frames
                tstr = f"{round(t_hr, 1):03}"
                tframe = tsec + t_hr*3600
                for i_om in range(n_peaks):
                    if output_Schmidt:
                        gnm_sph_rot[i_om,...] = gnm_sph[i_om,...] * np.exp(-1j * peak_omegas[i_om] * tframe)
                        gnm_rot[i_om,...]     = gnm[i_om,...]     * np.exp(-1j * peak_omegas[i_om] * tframe)
                        hnm_sph_rot[i_om,...] = hnm_sph[i_om,...] * np.exp(-1j * peak_omegas[i_om] * tframe)
                        hnm_rot[i_om,...]     = hnm[i_om,...]     * np.exp(-1j * peak_omegas[i_om] * tframe)
                    else:
                        Binm_sph_rot[i_om,...] = Binm_sph[i_om,...] * np.exp(-1j * peak_omegas[i_om] * tframe)
                        Binm_rot[i_om,...]     = Binm[i_om,...]     * np.exp(-1j * peak_omegas[i_om] * tframe)

                if output_Schmidt:
                    Binm_rot = (gnm_rot, hnm_rot)
                    Binm_sph_rot = (gnm_sph_rot, hnm_sph_rot)
                plots.plotMagSurf(n_peaks, Binm_rot, nvals, mvals, DO_LARGE in modelOpts, Schmidt=output_Schmidt, r_surf_mean=eval_r, asym_frac=asym_frac,
                                  pvals=pvals, qvals=qvals, difference=gif_diff, Binm_sph=Binm_sph_rot, nprmvals=nprmvals, mprmvals=mprmvals,
                                  bodyname=bname, append=bname_opt+flyby_opt, fend=iT_str, tstr=tstr, component=comp,
                                  absolute=True, no_title=False, marks=CAmarks)

            log.info("Animation frames printed to figures/anim_frames/ folder.\n" +
                     "Stack them into a gif with, e.g.:\n" +
                     "convert -delay 15 figures/anim_frames/Miranda_field_asym0*.png -loop 15 figures/anim_Miranda_asym.gif")
        else:
            for i_om in range(n_peaks):
                if output_Schmidt:
                    gnm_sph_rot[i_om,...] = gnm_sph[i_om,...] * np.exp(-1j * peak_omegas[i_om] * tsec)
                    gnm_rot[i_om,...]     = gnm[i_om,...]     * np.exp(-1j * peak_omegas[i_om] * tsec)
                    hnm_sph_rot[i_om,...] = hnm_sph[i_om,...] * np.exp(-1j * peak_omegas[i_om] * tsec)
                    hnm_rot[i_om,...]     = hnm[i_om,...]     * np.exp(-1j * peak_omegas[i_om] * tsec)
                else:
                    Binm_sph_rot[i_om,...] = Binm_sph[i_om,...] * np.exp(-1j * peak_omegas[i_om] * tsec)
                    Binm_rot[i_om,...]     = Binm[i_om,...]     * np.exp(-1j * peak_omegas[i_om] * tsec)

            # Plot symmetric moments/difference
            if output_Schmidt:
                Binm_rot = (gnm_rot, hnm_rot)
                Binm_sph_rot = (gnm_sph_rot, hnm_sph_rot)
            plots.plotMagSurf(n_peaks, Binm_rot, nvals, mvals, DO_LARGE in modelOpts, Schmidt=output_Schmidt, r_surf_mean=eval_r, asym_frac=asym_frac,
                            pvals=pvals, qvals=qvals, difference=plot_diffs, Binm_sph=Binm_sph_rot, nprmvals=nprmvals, mprmvals=mprmvals,
                            bodyname=bname, append=bname_opt+flyby_opt, component=comp, absolute=absolute, no_title=no_title_text,
                            marks=CAmarks)

        # Restrict plotting of certain diagnostic plots so we don't get spammed when we only want to do this for special conditions
        actually_plot_traces = bname == "Europa" and comp == "x" \
                               and (prevEuropa in modelOpts or TobieLow in modelOpts or TobieHigh in modelOpts)
        if (sub_planet_vert and actually_plot_traces) and not output_Schmidt:
            if not recalc:
                peak_periods, Benm, B0, exc_names, _ = asym.read_Benm(nprm_max_main, p_max_main, bodyname=bname, synodic=synodic_only, model=Benm_model, fname=inp_Be_path)
                int_model = os.path.join(inp_model_path, f"interior_model_asym{bfname}{bname_opt}{flyby_opt}.txt")
                r_bds, sigmas, bcdev = np.loadtxt(int_model, skiprows=1, unpack=True, delimiter=',')

            t_cut = vert_cut_hr * 3600
            vert_begin = np.maximum(vert_start, r_bds[-1]*rscale_moments)
            r = np.linspace(vert_begin, vert_stop, t_pts)
            x = r * np.cos(np.radians(vert_cut_lat)) * np.cos(np.radians(vert_cut_lon))
            y = r * np.cos(np.radians(vert_cut_lat)) * np.sin(np.radians(vert_cut_lon))
            z = r * np.sin(np.radians(vert_cut_lat))
            t = np.zeros(t_pts)
            t += t_cut
            plots.calcAndPlotTrajec(x,y,z,r,t, Binm, Benm, peak_omegas, nprm_max_main, n_max_main, nvals, mvals, 
                                    R_body=R, component=comp, difference=True, Binm_sph=Binm_sph, bodyname=bname, 
                                    append=bname_opt+flyby_opt+compstr)

        # Plot a time series at the sub-parent-planet point--(0 deg, 0 deg) in IAU coordinates.
        # Only tested for Europa, with no ionosphere.
        if (synodic_only and actually_plot_traces) and not output_Schmidt:
            if not recalc:
                synodic_period, Benm, B0, exc_names, _ = asym.read_Benm(nprm_max_main, p_max_main, bodyname=bname, synodic=True)
                peak_periods = [synodic_period]
            if not sub_planet_vert:
                int_model = os.path.join(inp_model_path, f"interior_model_asym{bfname}{bname_opt}{flyby_opt}.txt")
                r_bds, sigmas, bcdev = np.loadtxt(int_model, skiprows=1, unpack=True, delimiter=',')

            r = (localt + R) / R

            # Only valid if the evaluation point is outside the conducting region, so that B is subject to the Laplace equation.
            if r*R*1e3/r_bds[-1] >= 0.999:
                locx = r * np.cos(np.radians(loclat)) * np.cos(np.radians(loclon))
                locy = r * np.cos(np.radians(loclat)) * np.sin(np.radians(loclon))
                locz = r * np.sin(np.radians(loclat))
                loc = [ locx, locy, locz, r ] # Now IAU planetocentric in terms of body radii
                plots.plotTimeSeries(loc, Binm[0,...], Benm[0,...], tsec, peak_periods[0], nprm_max_main, n_max_main, nvals, mvals,
                                     n_pts=t_pts, component=comp, Binm_sph=Binm_sph[0,...], bodyname=bname, append=bname_opt+flyby_opt+compstr)

                if orbital_time_series:
                    if not recalc:
                        n_bds = np.size(r_bds)
                        eps_scaled = r_bds * bcdev
                        if relative:
                            single_asym = None
                        else:
                            single_asym = r_io
                        asym_shape, grav_shape = asym.read_shape(n_bds, p_max_main, rscale_asym, bodyname=bname, relative=relative, single_asym=single_asym,
                                                     eps_scaled=eps_scaled, r_bds=r_bds, r_io=r_io, append=bname_opt+flyby_opt, convert_depth_to_chipq=convert_depth_to_chipq)
                        r_bds, sigmas, asym_shape = asym.validate(r_bds, sigmas, bcdev, asym_shape, p_max_main)

                    orbital_period, Benm_orbital, B0, exc_names, _ = asym.read_Benm(nprm_max_main, p_max_main, bodyname=bname, orbital=True)
                    orbital_omega = 2*np.pi/(orbital_period*3600)
                    Binm_sph_orbital = sym.BiList(r_bds, sigmas, [orbital_omega], Benm_orbital, nprmvals, mprmvals, rscale_moments, n_max=nprm_max_main, bodyname=bname, append=bname_opt+flyby_opt)
                    Binm_orbital = asym.BiList(r_bds, sigmas, [orbital_omega], asym_shape, grav_shape, Benm_orbital, rscale_moments, nvals, mvals, p_max_main, nprm_max=nprm_max_main, bodyname=bname, append=bname_opt+flyby_opt)
                    plots.plotTimeSeries(loc, Binm_orbital[0, ...], Benm_orbital[0, ...], tsec, orbital_period, nprm_max_main, n_max_main, nvals, mvals, n_pts=t_pts, component=comp, Binm_sph=Binm_sph_orbital[0, ...], bodyname=bname, append=bname_opt+flyby_opt+compstr+"_orbital")
