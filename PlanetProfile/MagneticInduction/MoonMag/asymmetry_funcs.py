""" Contains functions for calculating induced magnetic fields
    from near-spherical conductors.
    Developed in Python 3.8 for "A perturbation method for evaluating the
    magnetic field induced from an arbitrary, asymmetric ocean world
    analytically" by Styczinski et al.
    DOI: 10.1016/j.icarus.2021.114840
Author: M. J. Styczinski, mjstyczi@uw.edu """

import os
import csv
import numpy as np
from scipy.io import savemat, loadmat
from collections.abc import Iterable
from numpy import sqrt
from math import floor
from scipy.special import factorial as ft, lpmv as legenp, sph_harm

import mpmath as mp
# mpmath is needed for enhanced precision to avoid
# divide-by-zero errors induced by underflow.

from PlanetProfile.MagneticInduction.MoonMag import _excitation, _interior, _induced
from PlanetProfile.MagneticInduction.MoonMag.config import *
from PlanetProfile.MagneticInduction.MoonMag.field_xyz import eval_Bi, eval_Bi_Schmidt

# Parallelization is through multiprocessing module
import multiprocessing as mtp
import platform
plat = platform.system()
if plat == 'Windows':
    mtpType = 'spawn'
else:
    mtpType = 'fork'
mtpContext = mtp.get_context(mtpType)
num_cores = mtp.cpu_count()

# Global variables and settings
# Set maximum precision for mpmath quantities
mp.mp.dps = digits_precision
# Numerical constants in high-precision mpmath format
zero = mp.mpf(0)
one = mp.mpf(1)
two = mp.mpf(2)
j = mp.mpc(1j)
mu_o = mp.mpf("4.0e-7")*mp.pi
sqrt2 = sqrt(2)
sqrt4pi = np.sqrt(4*np.pi)

# Get type objects for initializing numpy arrays of mpf/mpc
mpf_global = type(zero)
mpc_global = type(j)

#############################################

"""
read_shape()
    Get boundary shape parameters from disk.
    Usage: `asym_shape` = read_shape(`n_bds`, `p_max`, `rscale`, `bodyname=None`, `relative=False`, `eps_scaled=None`,
        `single_asym=None`, `concentric=True`, `fpath=None`, `r_bds=None`, `r_io=-2`, `append=""`, `convert_depth_to_chipq=False`)
    Returns:
        asym_shape: complex, shape(n_bds,2,p_max+1,p_max+1). A list of absolute boundary shape parameters chi_pq
            for each boundary in m, converted from values in km read from .txt files.
        grav_shape: complex, shape(n_bds,2,p_max+1,p_max+1). Analogous to asym_shape, but due to tides deforming the
            body surface. Kept separate from asym_shape to preserve ice shell thickness contours for comparison to prior work.
    Parameters:
        n_bds: integer. Number of conducting boundaries in the model (this is N).
        p_max: integer. Largest degree p in boundary shape.
        rscale: float. 1/R_body; used to scale absolute asymmetry values in making asym_shape arrays.
        bodyname: string (None). Body name as it appears in file names. If None, a generic file name will be searched for.
        relative: boolean (True). Optional argument to toggle whether chi_pq values are normalized
            to bcdev (True) or give absolute coefficients for the deviations (False).
        eps_scaled: float, shape(n_bds) (None). Maximum boundary deviation in km, for converting relative asym_shape values to absolute.
            Required if relative = True.
        single_asym: integer (None). If not None, apply asymmetry only to the boundary with this index.  
        concentric: boolean (True). Optional argument to toggle whether to use the same spherical harmonic
            coefficients for every boundary (concentric asymmetry). In this case, only a single file is read in.
            If concentric = True, r_bds is required.
        fpath: string (None). Optional location to search for .txt files. Defaults to "interior/".
        r_bds: float, shape(n_bds) (None). Boundary radii in m. Required if concentric = True.
        r_io: integer (-2). The index of r_bds corresponding to the ice-ocean interface, where asymmetry is expected to
            be its most pronounced. Defaults to -2, i.e. no ionosphere. Set to -3 for the ice-ocean boundary when a 
            symmetric ionosphere (one layer) is included.
        append: string (""). Optional string appended to file names.
        convert_depth_to_chipq: boolean (False). Optional flag to print relative chi_pq values to disk.
    """
def read_shape(n_bds, p_max, rscale, bodyname=None, relative=False, eps_scaled=None, single_asym=None, concentric=True,
               fpath=None, r_bds=None, r_io=-2, append="", convert_depth_to_chipq=False):
    if fpath is None:
        fpath = _interior
    if bodyname is None:
        bfname = ""
    else:
        bfname = f"_{bodyname}"

    # Initialize asymmetric shape array. The latter 3 indices are
    # for q positive (0) or negative (1), p, and |q|.
    asym_shape, grav_shape = ( np.zeros((n_bds,2,p_max+1,p_max+1),dtype=np.complex128) for _ in range(2) )

    if single_asym is not None:
        asym_model = os.path.join(fpath, f"depth_chi_pq_shape{bfname}{append}.txt")
        log.debug(f"Using asymmetry model: {asym_model} for layer index {single_asym}")
        shape_n = np.loadtxt(asym_model, skiprows=1, unpack=False, delimiter=',')
        scaled_rad = r_bds * rscale
        for p in range(1, p_max + 1):
            this_min = int(p * (p + 1) / 2)
            this_max = this_min + p + 1
            # Negate to convert from deviations of depth to radius
            Cpq = -shape_n[this_min:this_max, 2]
            Spq = -shape_n[this_min:this_max, 3]

            chi_pqs = get_chipq_from_CSpq_single(p, Cpq, Spq)
            asym_shape[single_asym,:,p,:p+1] = chi_pqs * 1e3 * scaled_rad[single_asym]

        if bodyname == "Miranda":
            # Scale ice shell thickness variations for Enceladus to Miranda
            asym_shape = asym_shape * 49.810/20.852
    elif relative:
        # Put scalar value into a list in case there's only a single boundary
        if not isinstance(eps_scaled, Iterable): eps_scaled = [ eps_scaled ]

        asym_model1 = os.path.join(fpath, "degree")

        for p in range(1,p_max+1):
            try:
                asym_model = f"{asym_model1}{p}_shapes{bfname}{append}.txt"
                shape_p = np.loadtxt(asym_model, skiprows=1, unpack=False, delimiter=',', dtype=np.complex128)
                log.debug(f"Using asymmetry model: {asym_model}")
            except:
                asym_model = f"{asym_model1}{p}_shapes.txt"
                shape_p = np.loadtxt(asym_model, skiprows=1, unpack=False, delimiter=',', dtype=np.complex128)
                log.debug(f"Using asymmetry model: {asym_model}")

            qcount = 0
            for q in range(-p, p+1):
                qsign = int(q<0)
                qabs = abs(q)
                for i in range(n_bds):
                    asym_shape[i,qsign,p,qabs] = shape_p[i,qcount]*eps_scaled[i]
                qcount += 1
    elif concentric:
        asym_model = os.path.join(fpath, f"depth_chi_pq_shape{bfname}{append}.txt")
        log.debug(f"Using asymmetry model: {asym_model} concentrically, scaled to a radius of {1/rscale:.2f}")
        shape_n = np.loadtxt(asym_model, skiprows=1, unpack=False, delimiter=',')
        scaled_rad = r_bds * rscale
        for p in range(1,p_max+1):
            this_min = int(p*(p+1)/2)
            this_max = this_min + p+1
            # Negate to convert from deviations of depth to radius
            Cpq = -shape_n[this_min:this_max,2]
            Spq = -shape_n[this_min:this_max,3]

            chi_pqs = get_chipq_from_CSpq_single(p,Cpq,Spq)
            for i_layer in range(n_bds):
                asym_shape[i_layer,:,p,:p+1] = chi_pqs*1e3 * scaled_rad[i_layer]

    else:
        # If we got here, none of the options were selected.
        # Apply asymmetry independently to each layer, and read in one file per layer.
        raise RuntimeError("Independently asymmetric layers are not implemented.")
        for n_layer in range(1, n_bds+1):
            shape_n = np.loadtxt(os.path.join(fpath, f"depth_chi_pq_shape{n_layer}_{bfname}.txt"), skiprows=1, unpack=False, delimiter=',')
            for p in range(1, p_max+1):
                """ Not yet implemented, but this is a starting point.
                this_min = int(p*(p+1)/2)
                this_max = this_min + p+1
                Cpq = -shape_n[this_min:this_max,2]
                Spq = -shape_n[this_min:this_max,3]

                chi_min = p**2 - 1
                chi_max = (p+1)**2 - 1
                asym_shape[n_layer,chi_min:chi_max] = get_chipq_from_CSpq_single(p,Cpq,Spq)"""

    if convert_depth_to_chipq and not relative:
        if bodyname == "Europa":
            if append == "_Tobie":
                eps_max = 2500
            elif append == "_prev":
                eps_max = 2500
            else:
                eps_max = 30000
        elif bodyname == "Callisto":
            eps_max = 100000
        elif bodyname == "Triton":
            eps_max = 100000
        elif bodyname == "Miranda":
            eps_max = 49810/20852 * 16377 # Scaling up from maximum deviation in Hemingway and Mittal (2019) Enceladus model to thicker Miranda ice shell
        else:
            eps_max = 16377
        # Whether to print lines suitable for copying to boundary shape files,
        # or to print for readability instead
        print_for_copy = True
        if print_for_copy:
            chi_head = f" p, chipq (-q, ..., q); bcdev={eps_max} m\n"
        else:
            chi_head = f" p, q, chipq_re, chipq_im, bcdev={eps_max} m\n"
        eps = asym_shape[r_io, ...] / eps_max
        asym_out = os.path.join(fpath, f"chi_pq_{bodyname}.txt")
        with open(asym_out, "w") as f_chi:
            f_chi.write(chi_head)
            for p in range(1, p_max+1):
                if print_for_copy:
                    this_line = " " + str(p)
                for q in range(-p, p+1):
                    qsign = int(q<0)
                    qabs = abs(q)
                    if print_for_copy:
                        if np.imag(eps[qsign, p, qabs]) >= 0:
                            btw_char = "+"
                        else:
                            btw_char = "-"
                        this_line = f"{this_line}, {np.real(eps[qsign, p, qabs])}{btw_char}{np.abs(np.imag(eps[qsign, p, qabs]))}j"
                    else:
                        this_line = f" {p}, {str(q).rjust(2)}, {np.real(eps[qsign, p, qabs])}{btw_char}{np.abs(np.imag(eps[qsign, p, qabs]))}\n"
                        f_chi.write(this_line)
                if print_for_copy:
                    f_chi.write(this_line+"\n")
        log.debug(f"Printed chi_pq values to {asym_out}")

    # If this file exists, model tidal perturbations. If not, return None.
    grav_model = os.path.join(fpath, f"gravity{bfname}{append}.txt")
    try:
        g_shape_n = np.loadtxt(grav_model, skiprows=1, unpack=False, delimiter=',')
    except:
        grav_shape = None
    else:
        log.debug(f"Using surface gravity shape: {grav_model}")
        for p in range(1, p_max+1):
            this_min = int(p * (p+1) / 2)
            this_max = this_min + p+1
            # Unlike above, we do not need to negate because deviations are already in radii
            gCpq = g_shape_n[this_min:this_max, 2]
            gSpq = g_shape_n[this_min:this_max, 3]

            g_chi_pqs = get_chipq_from_CSpq_single(p, gCpq, gSpq)
            for ii in range(n_bds + r_io + 1):
                # If we are modeling asymmetry in any surface, it is adjacent to a conducting layer
                if (asym_shape[ii, ...] != 0).any():
                    # Apply surface tidal perturbations to all asymmetric layers equally
                    # (to preserve thicknesses determined by thermodynamics)
                    grav_shape[ii, :, p, :p+1] = g_chi_pqs * 1e3

            # Always add perturbation to body surface
            grav_shape[n_bds+r_io+1, :, p, :p+1] = g_chi_pqs * 1e3
    return asym_shape, grav_shape

#############################################

"""
get_chipq_from_CSpq_single()
    Convert from real, 4pi-normalized harmonic coefficients with no Condon-Shortley phase
    (the common normalization in the geodesy community) to orthonormal harmonic coefficients having
    the C-S phase. Handles all values for a given p at once.
    Usage: `chipq` = get_chipq_from_CSpq_single(`p`,`Cpq`,`Spq`)
    Returns:
        chipq: complex, shape(2,p+1). chi_pq values for all q = [-p,p], organized such that chipq[int(q<0),abs(q)]
            returns the result for a particular q value. Orthonormal, with Condon-Shortley phase.
            Orthonormal here means the integral of |Ynm|^2 * dOmega over a unit sphere is 1 for all n and m.
    Parameters:
        p: integer. Degree of boundary shapes; results for all q values are returned for this p value.
        Cpq: float. Spherical harmonic coefficient that multiplies cos(m*phi) for p,q boundary. 4pi-normalized with no Condon-Shortley phase.
        Spq: float. Spherical harmonic coefficient that multiplies sin(m*phi) for p,q boundary. 4pi-normalized with no Condon-Shortley phase.
        CSchange: bool (True). If False, do NOT correct for Condon-Shortley phase, i.e. return gnm, hnm
            that contain the CS phase.
    """
def get_chipq_from_CSpq_single(p, Cpq, Spq, CSchange=True):
    chipq = np.zeros((2,p+1),dtype=np.complex128)

    norm = sqrt4pi / sqrt2

    for q in range(1,p+1):
        # Negative q (first index 1):
        chipq[1,q] = (Cpq[q] + 1j*Spq[q]) * norm
        # Positive q (first index 0):
        chipq[0,q] = (-1)**(q*CSchange) * (Cpq[q] - 1j*Spq[q]) * norm

    chipq[0,0] = Cpq[0] * sqrt4pi

    return chipq

#############################################

"""
get_chipq_from_CSpq()
    Convert from real, 4pi-normalized harmonic coefficients with no Condon-Shortley phase
    (the common normalization in the geodesy community) to orthonormal harmonic coefficients having
    the C-S phase. Handles all values at once.
    Usage: `chipq` = get_chipq_from_CSpq(`pmax`,`Cpq`,`Spq`)
    Returns:
        chipq: complex, shape(2,pmax+1,pmax+1). chi_pq values for all q = [-p,p] for p up to pmax,
            organized such that chipq[int(q<0),p,abs(q)] returns the result for a particular q value.
            Orthonormal, with Condon-Shortley phase. Orthonormal here means the integral of
            |Ynm|^2 * dOmega over a unit sphere is 1 for all n and m.
    Parameters:
        pmax: integer. Maximum degree of boundary shapes; results for all p and q values up to this pmax are returned.
        Cpq: float, shape(pmax+1,pmax+1). Spherical harmonic coefficient that multiplies cos(m*phi) for p,q boundary. 4pi-normalized with no Condon-Shortley phase.
        Spq: float, shape(pmax+1,pmax+1). Spherical harmonic coefficient that multiplies sin(m*phi) for p,q boundary. 4pi-normalized with no Condon-Shortley phase.
        CSchange: bool (True). If False, do NOT correct for Condon-Shortley phase, i.e. return gnm, hnm
            that contain the CS phase.
    """
def get_chipq_from_CSpq(pmax, Cpq, Spq, CSchange=True):
    chipq = np.zeros((2,pmax+1,pmax+1),dtype=np.complex128)

    norm = sqrt4pi / sqrt2

    for p in range(pmax+1):
        for q in range(1, p+1):
            # Negative q (first index 1):
            chipq[1,p,q] = (Cpq[p,q] + 1j*Spq[p,q]) * norm
            # Positive q (first index 0):
            chipq[0,p,q] = (-1)**(q*CSchange) * (Cpq[p,q] - 1j*Spq[p,q]) * norm

        chipq[0,p,0] = Cpq[p,0] * sqrt4pi

    return chipq

#############################################

"""
get_CSpq_from_chipq()
    Convert from orthonormal harmonic coefficients having the C-S phase to real, 4pi-normalized
    harmonic coefficients with no Condon-Shortley phase. Handles all values at once.
    Usage: `Cpq`, `Spq` = get_chipq_from_CSpq(`pmax`,`chipq`)
    Returns:
        Cpq: float, shape(pmax+1,pmax+1). Spherical harmonic coefficient that multiplies cos(m*phi) for p,q boundary. 4pi-normalized with no Condon-Shortley phase.
        Spq: float, shape(pmax+1,pmax+1). Spherical harmonic coefficient that multiplies sin(m*phi) for p,q boundary. 4pi-normalized with no Condon-Shortley phase.
    Parameters:
        pmax: integer. Maximum degree of boundary shapes; results for all p and q values up to this pmax are returned.
        chipq: complex, shape(2,pmax+1,pmax+1). chi_pq values for all q = [-p,p] for p up to pmax,
            organized such that chipq[int(q<0),p,abs(q)] returns the result for a particular q value.
            Orthonormal, with Condon-Shortley phase. Orthonormal here means the integral of
            |Ynm|^2 * dOmega over a unit sphere is 1 for all n and m.
        CSchange: bool (True). If False, do NOT correct for Condon-Shortley phase, i.e. return gnm, hnm
            that contain the CS phase.
        noRenorm: bool (False). If True, return Cpq and Spq as fully normalized, real harmonic coefficients.
    """
def get_CSpq_from_chipq(pmax, chipq, CSchange=True, noRenorm=False):
    Cpq, Spq = ( np.zeros((pmax+1,pmax+1)) for _ in range(2) )

    if noRenorm:
        norm = 1
        norm0 = 1
    else:
        norm = 1 / sqrt4pi / sqrt2
        norm0 = norm * sqrt2

    for q in range(1, pmax+1):
        CS = (-1)**(q*(not CSchange))
        # cos(m*phi) part:
        Cpq[:,q] = np.real(((-1)**q * chipq[0,:,q] + chipq[1,:,q]) * norm) * CS
        # sin(m*phi) part:
        Spq[:,q] = np.real(((-1)**q * chipq[0,:,q] - chipq[1,:,q]) * 1j * norm) * CS

    Cpq[:,0] = np.real(chipq[0,:,0] * norm0)

    return Cpq, Spq

#############################################

"""
get_apq_from_chipq()
    Convert complex orthonormal harmonic coefficients with the Condon-Shortley phase to 
    the format expected by healpy routines, i.e. a linear list with the m<0 values
    excluded. Handles all values at once.
    Usage: `apq` = get_apq_from_chipq(`pvals`, `qvals`, `chipq`)
    Returns:
        apq: complex, shape(pvals). Orthonormal spherical harmonic coefficients with the C-S phase for m>0 only. 
            Coefficients for m<0 are simply discarded; the values are inferred from what they must be to match the
            m>0 values in order to return a real-valued map from healpy.alm2map.
    Parameters:
        pvals, qvals: integer, shape(Npq). List of paired p and q indices for apq, as returned from healpy.sphtfunc.Alm.getlm(pmax).
        chipq: complex, shape(2,pmax+1,pmax+1). chi_pq values for all q = [-p,p] for p up to pmax,
            organized such that chipq[int(q<0),p,abs(q)] returns the result for a particular p,q pair.
            Orthonormal, with Condon-Shortley phase. Orthonormal here means the integral of
            |Ynm|^2 * dOmega over a unit sphere is 1 for all n and m.
    """
def get_apq_from_chipq(pvals, qvals, chipq):
    # Arrange into a linear list compatible with healpy evaluation routines
    apq = np.array([chipq[0,p,q] for p,q in zip(pvals, qvals)])

    return apq

#############################################

"""
get_chipq_from_apq()
    Convert complex orthonormal harmonic coefficients with the Condon-Shortley phase from 
    the format expected by healpy routines, i.e. a linear list with the m<0 values
    excluded, to the full set for all values of -p <= m <= p. Handles all values at once.
    Usage: `chipq` = get_chipq_from_apq(`pvals`, `qvals`, `apq`)
    Returns:
        chipq: complex, shape(2,pmax+1,pmax+1). chi_pq values for all q = [-p,p] for p up to pmax,
            organized such that chipq[int(q<0),p,abs(q)] returns the result for a particular p,q pair.
            Orthonormal, with Condon-Shortley phase. Orthonormal here means the integral of
            |Ynm|^2 * dOmega over a unit sphere is 1 for all n and m.
    Parameters:
        pvals, qvals: integer, shape(Npq). List of paired p and q indices for apq, as returned from healpy.sphtfunc.Alm.getlm(pmax).
        apq: complex, shape(pvals). Orthonormal spherical harmonic coefficients with the C-S phase for m>0 only. 
            Coefficients for m<0 are simply discarded; the values are inferred from what they must be to match the
            m>0 values in order to return a real-valued map from healpy.alm2map.
    """
def get_chipq_from_apq(pvals, qvals, apq):
    # Reconstruct all complex coefficients as needed to result in real-valued outputs
    pmax = np.max(pvals)
    chipq = np.zeros((2,pmax+1,pmax+1), dtype=np.complex128)
    for a,p,q in zip(apq, pvals, qvals):
        chipq[0,p,q] = a
        chipq[1,p,q] = np.conj(a) * (-1)**q

    return chipq

#############################################

"""
get_gh_from_Binm()
    Convert from orthonormal harmonic coefficients with the Condon-Shortley phase (common in physics)
    to Schmidt semi-normalized form without the C-S phase (common in geophysics).
    Handles all values for a given n at once.
    Usage: `gnm`, `hnm` = get_gh_from_Binm(`n`, `n_max`, `Binm`)
    Returns:
        gnm, hnm: complex, shape(n_max+1,n_max+1). g_nm and h_nm values for all m = [0,n].

            Schmidt normalization here means the integral of |Ynm|^2 * dOmega over a unit sphere is
            4pi/(2n+1) for all n and m. No Condon-Shortley phase.
    Parameters:
        n_max: integer. Maximum degree n of induced moments.
        Binm: complex, shape(2,n_max+1,n_max+1). Complex induced magnetic moments calculated using fully normalized spherical harmonic coefficients.
    """
def get_gh_from_Binm(n_max, Binm):
    gnm, hnm = ( np.zeros((n_max+1,n_max+1), dtype=np.complex128) for _ in range(2) )

    for n in range(1,n_max+1):
        norm = np.sqrt(2*n+1) / sqrt2 / sqrt4pi

        for m in range(1, n+1):
            # g terms:
            gnm[n,m] = ((-1)**m * Binm[0,n,m] + Binm[1,n,m]) * norm
            # h terms:
            hnm[n,m] = ((-1)**m * Binm[0,n,m] - Binm[1,n,m]) * 1j * norm

        gnm[n,0] = Binm[0,n,0] * norm * sqrt2

    return gnm, hnm

#############################################

def get_Binm_from_gh(n_max, gnm, hnm):
    """
    Convert from Schmidt semi-normalized form without the C-S phase (common in geophysics)
    to orthonormal harmonic coefficients with the Condon-Shortley phase (common in physics).
    Handles all values for a given n at once.
    Args:
        n_max: int. Maximum degree n of induced moments.
        gnm, hnm: complex, shape(n_max+1,n_max+1). g_nm and h_nm values for all m = [0,n].
            Schmidt normalization here means the integral of |Ynm|^2 * dOmega over a unit sphere is
            4pi/(2n+1) for all n and m. No Condon-Shortley phase.

    Returns:
        Binm: complex, shape(2,n_max+1,n_max+1). Complex induced magnetic moments for fully normalized spherical harmonics.
    """
    Binm = np.zeros((2,n_max+1,n_max+1), dtype=np.complex128)

    for n in range(1, n_max+1):
        norm = sqrt4pi / np.sqrt(2*n+1) / sqrt2

        for m in range(1, n+1):
            # m>0 terms:
            Binm[0,n,m] = (-1)**m * norm * (gnm[n,m] - 1j * hnm[n,m])
            # m<0 terms:
            Binm[1,n,m] =           norm * (gnm[n,m] + 1j * hnm[n,m])

        Binm[0,n,0] = norm * sqrt2 * gnm[n,0]

    return Binm

#############################################

"""
validate()
    Check inputs to be sure everything will be interpreted correctly.
    Usage: `r_bds`, `sigmas`, `omegas`, `asym_shape` = validate(`r_bds`, `sigmas`, `omegas`, `bcdev`, `asym_shape`, `p_max`)
    Returns:
        r_bds: float, shape(n_bds).
        sigmas: float, shape(n_bds).
        omegas: float, shape(n_peaks).
        asym_shape: complex, shape(n_bds,2,p_max+1,p_max+1).
    Parameters:
        r_bds: float, shape(n_bds).
        sigmas: float, shape(n_bds).
        omegas: float, shape(n_peaks).
        bcdev: float, shape(n_bds).
        asym_shape: complex, shape(n_bds,2,p_max+1,p_max+1).
        p_max: integer.
    """
def validate(r_bds, sigmas, bcdev, asym_shape, p_max):
    # Check lengths of model lists
    if np.shape(r_bds) != np.shape(sigmas):
        log.error(f"boundaries shape: {np.shape(r_bds)}")
        log.error(f"sigmas shape: {np.shape(sigmas)}")
        raise ValueError("The number of boundaries is not equal to the number of conductivities.")
    if np.shape(r_bds) != np.shape(bcdev):
        log.error(f"boundaries shape: {np.shape(r_bds)}")
        log.error(f"deviations shape: {np.shape(bcdev)}")
        raise ValueError("The number of boundaries is not equal to the number of deviations.")

    # Make sure interior model is iterable (it's not if there is only 1 boundary)
    if not isinstance(r_bds, Iterable):
        r_bds = [r_bds]
        sigmas = [sigmas]
        np.reshape(asym_shape, (1,2,p_max+1,p_max+1))

    # Double check array lengths all match up
    if np.shape(asym_shape) != (np.size(r_bds),2,p_max+1,p_max+1):
        if np.shape(r_bds) != np.shape(asym_shape)[0]:
            log.error(f"boundaries length: {np.shape(r_bds)}")
            log.error(f"deviations length: {np.shape(asym_shape)[0]}")
            raise ValueError("The number of boundaries is not equal to the number of deviation shapes.")
        else:
            raise ValueError("The number of deviation shapes is not equal to (p_max + 1)^2 - 1.")

    return r_bds, sigmas, asym_shape

#############################################

"""
jnx(), ynx(), jdx(), ydx()
    Spherical Bessel and Neumann functions for degree n and argument x.
    jdx and ydx are j_n^\star and y_n^\star in the above publication, and are
    equal to d(r*jnx(kr))/dr and d(r*ynx(kr))/dr, respectively.
    Usage: `eval` = jnx(`n`, `x`)
    Returns:
        eval: mpc, shape(len(x)). List of results of Bessel functions in mpmath complex (mpc) format. Always returns an array, possibly of length 1.
    Parameters:
        n: mpf. Degree of spherical Bessel function. Should be an integer, despite being of mpf type (must be passed as half-integer order to besselj).
        x: mpc, shape(n_bds). Complex argument of Bessel function (kr). MUST be a list or this will break. Make single values into lists with x = [x].
    """
def jnx(n,x):
    two = mp.mpf(2)
    jnx = np.array([ mp.besselj( n+mp.mpf("0.5"),xi ) * mp.sqrt(mp.pi/two/xi) for xi in x ])
    return jnx
def ynx(n,x):
    two = mp.mpf(2)
    ynx = np.array([ mp.bessely( n+mp.mpf("0.5"),xi ) * mp.sqrt(mp.pi/two/xi) for xi in x ])
    return ynx
def jdx(n,x):
    one = mp.mpf(1)
    jn = jnx(n,x)
    jP = jnx(n+one,x)
    jdx = np.array([ (n+one)*jn[i] - x[i]*jP[i] for i in range(np.size(x)) ])
    return jdx
def ydx(n,x):
    one = mp.mpf(1)
    yn = ynx(n,x)
    yP = ynx(n+one,x)
    ydx = np.array([ (n+one)*yn[i] - x[i]*yP[i] for i in range(np.size(x)) ])
    return ydx

#############################################

"""
read_Benm()
    Read in the complex frequency spectrum of excitation field oscillations for the strongest excitations.
    Usage: `peak_periods`, `Benm`, `B0`, `exc_names`, 'Bexyz' = read_Benm(`nprm_max`, `p_max`, `bodyname=None`, `fpath=None`, `synodic=False`, `orbital=False`)
    Returns:
        peak_periods: float, shape(n_peaks). Periods in hr of peak oscillations. Values read from files are assumed to be in hr.
        Benm: complex, shape(n_peaks,2,nprm_max+p_max+1,nprm_max+p_max+1). Excitation moments for each period in nT.
        B0: float, shape(3). Static background field. Used to reconstruct the net magnetic field for measurement comparisons.
        exc_names: string, shape(n_peaks). Descriptive names for peak oscillations.
    Parameters:
        nprm_max: integer. Maximum degree of excitation harmonics to input. n' values appear in file names, so appropriate files must
            all be present from n'=1 to n'=nprm_max.
        p_max: integer. Maximum degree of boundary shapes. Required to size Benm array correctly for fast computation of boundary conditions.
        bodyname: string (None). Body name as it appears in file names. If None, a generic file name will be searched for.
        fpath: string (None). Optional location to search for excitation moment files. Defaults to "excitation/".
        synodic: boolean (False). Option to consider only the strongest excitation period, for simplicity.
        orbital: boolean (False). Option to consider only the orbital period, analogous to the synodic period above.
        limit_osc: integer (None). Number of oscillation frequencies to limit calculations to, where 1 is just the strongest oscillation.
        model: string (None). Optional model code to append to file name.
        fName: string (None). Optional override of naming conventions to specify filename.  
    """
def read_Benm(nprm_max, p_max, bodyname=None, fpath=None, synodic=False, orbital=False, limit_osc=None,
              model=None, fName=None):
    if fpath is None:
        fpath = _excitation
    if model is None:
        modelStr = ""
    else:
        modelStr = f"_{model}"
    if bodyname is None:
        bfname = f"{modelStr}"
    else:
        bfname = f"_{bodyname}{modelStr}"

    if synodic and orbital:
        log.warning("Both 'synodic' and 'orbital' options passed to asymmetry_funcs.read_Benm. Only synodic will be used.")
        orbital = False

    if synodic:
        log.warning("Considering synodic period only for excitation.")
        if fName is None:
            Benm_moments = os.path.join(fpath, f"synodic_Be1xyz{bfname}.txt")
        else:
            limit_osc = 1
            Benm_moments = os.path.join(fpath, f"{fName}.txt")
    elif orbital:
        log.warning("Considering orbital period only for excitation.")
        if bodyname == "Europa":
            log.warning("Extra warning: This excitation is APPROXIMATE, in that two closely-spaced periods have been summed together and set to T = 85.2 h.")
        Benm_moments = os.path.join(fpath, f"orbital_Be1xyz{bfname}.txt")
    else:
        if fName is None:
            fName = f'Be1xyz{bfname}'
        Benm_moments = os.path.join(fpath, f"{fName}.txt")
    log.debug(f"Using excitation moments: {Benm_moments}")
    _, peak_per1, B0x, B0y, B0z, Bex_Re, Bex_Im, Bey_Re, Bey_Im, Bez_Re, Bez_Im \
        = np.genfromtxt(Benm_moments, skip_header=True, unpack=True, delimiter=',')
    n_peaks1_prelim = peak_per1.size

    with open(Benm_moments) as f_Benm:
        data = csv.reader(f_Benm)
        exc1_names = np.array([row[0] for row in data])

    Bex = Bex_Re + 1j * Bex_Im
    Bey = Bey_Re + 1j * Bey_Im
    Bez = Bez_Re + 1j * Bez_Im

    if limit_osc is None:
        limit_osc = 0
    elif limit_osc != 0:
        if n_peaks1_prelim <= limit_osc:
            log.debug(f"limit_osc is set to {limit_osc}, but {Benm_moments} contains only {n_peaks1_prelim} oscillation periods. Including all oscillations in calculations.")
            limit_osc = 0
        else:
            if limit_osc < 0: limit_osc = abs(limit_osc)
            log.debug(f"Limiting Benm to strongest {limit_osc} oscillations.")
    absBenm = np.sqrt(Bex_Re**2 + Bex_Im**2 + Bey_Re**2 + Bey_Im**2 + Bez_Re**2 + Bez_Im**2)
    if np.size(absBenm) > 1:
        iMax = np.argpartition(absBenm, -limit_osc)[-limit_osc:]
        peak_per1_lim = np.sort(peak_per1[iMax])
        iMaxSort = np.array([np.argwhere(peak_per1 == peak_per1_lim[i]) for i in range(np.size(peak_per1_lim))]).flatten()

        peak_per1 = peak_per1[iMaxSort]
        Bex = Bex[iMaxSort]
        Bey = Bey[iMaxSort]
        Bez = Bez[iMaxSort]
        exc1_names = exc1_names[iMaxSort]
        
    n_peaks1 = np.size(peak_per1)

    Benm1 = np.zeros((n_peaks1,2,nprm_max+p_max+1,nprm_max+p_max+1), dtype=np.complex128)

    A1 = sqrt(2*np.pi/3)

    Benm1[:,1,1,1] = -A1 * (Bex + 1j*Bey)
    Benm1[:,0,1,0] = -A1*sqrt(2) * Bez
    Benm1[:,0,1,1] =  A1 * (Bex - 1j*Bey)

    # Get static background field
    B0 = np.array([np.mean(B0x), np.mean(B0y), np.mean(B0z)])

    # Collect Bexyz to return for convenience
    Bexyz = (Bex, Bey, Bez)
    
    if nprm_max > 1:
        log.warning("n'=2 excitation moments are being considered, but the amplitudes are just a placeholder! Remove this warning when correct Be2xyz have been installed.")
        for nn in range(1,nprm_max+1):
            if nn == 2:
                Benm2_moments = f"Be2xyz{bfname}.txt"
                _, peak_per2, xBex_Re, xBex_Im, xBey_Re, xBey_Im, xBez_Re, xBez_Im, yBez_Re, yBez_Im, zBez_Re, zBez_Im    \
                    = np.genfromtxt(os.path.join(fpath, Benm2_moments), skip_header=True, unpack=True, delimiter=',')
                n_peaks2 = np.size(peak_per2)

                with open(Benm2_moments) as f_Benm2:
                    data = csv.reader(f_Benm2)
                    exc2_names = np.array([f'Be2{row[0]}' for row in data])

                Benm2 = np.zeros((n_peaks2,2,nprm_max+p_max+1,nprm_max+p_max+1), dtype=np.complex128)

                xBex = xBex_Re + 1j*xBex_Im
                xBey = xBey_Re + 1j*xBey_Im
                xBez = xBez_Re + 1j*xBez_Im
                yBez = yBez_Re + 1j*yBez_Im
                zBez = zBez_Re + 1j*zBez_Im

                ort8 = 1/sqrt(8)
                A2 = -sqrt(6*np.pi/5)

                Benm2[:,1,2,2] =  A2 * (xBex + 1j*xBey + ort8*xBez)
                Benm2[:,1,2,1] =  A2 * (xBez - 1j*yBez)
                Benm2[:,0,2,0] =  A2 * sqrt(3/2) * zBez
                Benm2[:,0,2,1] = -A2 * (xBez + 1j*yBez)
                Benm2[:,0,2,2] =  A2 * (xBex - 1j*xBey + ort8*xBez)

                Benm = np.zeros((n_peaks1+n_peaks2,2,nprm_max+p_max+1,nprm_max+p_max+1),dtype=np.complex128)

                Benm[:n_peaks1,:,:,:] = Benm1
                Benm[n_peaks1:,:,:,:] = Benm2
                peak_periods = np.append(peak_per1, peak_per2)
                exc_names = np.append(exc1_names, exc2_names)

            elif nn > 2:
                raise ValueError("n_max >= 3 is not implemented.")
    else:
        Benm = Benm1
        peak_periods = peak_per1
        exc_names = exc1_names

    if n_peaks1 == 1:
        peak_periods = np.array([peak_periods])

    return peak_periods, Benm, B0, exc_names, Bexyz

#############################################

"""
BiList()
    Evaluate the induced magnetic moments Binm for the given excitation moments Benm,
    for the given interior structure described by r_bds, sigmas, and asym_shape and for a list of omega values.
    Usage: `Binms` = BiList(`r_bds`, `sigmas`, `peak_omegas`, `asym_shape`, `Benm`, `nprmvals`, `mprmvals`, `rscale_moments`,
                         `pvals`, `qvals`, `nvals`, `mvals`, `nprm_max`, `p_max`, `writeout=True`, `path=None`,
                         `append=""`, `debug=False`)
    Returns:
        Binms: complex, shape(n_peaks,2,n_max+1,n_max+1). List of complex induced magnetic moments for the given parameters.
    Parameters:
        r_bds: float, shape(n_bds). List of mean radii for each boundary surface in m.
        sigmas: float, shape(n_bds). List of conductivity values for each layer underneath the corresponding boundary in S/m.
        peak_omegas: float, shape(n_peaks). List of angular frequencies in rad/s to evaluate. Typically only the peaks of
            the Fourier spectra.
        asym_shape_layers: complex, shape(n_bds,2,p_max+1,p_max+1). Absolute deviations for each boundary in m, ignoring
            gravitational perturbations.
        grav_shape: complex, shape(n_bds,2,p_max+1,p_max+1). Absolute deviations for each (conducting) boundary in m,
            due to gravitational perturbations.
        Benm: complex, shape(2,n_max+1,n_max+1). Excitation moments for each degree and order n',m'.
        rscale_moments: float. Equal to 1/R_body in units of 1/m. For proper scaling when conducting boundaries
            may or may not coincide with the body surface.
        nvals, mvals: integer, shape(Nnm). Linear arrays of paired n,m values for parallel computation.
        p_max: integer. Maximum degree p of boundary shapes.
        nprm_max: integer (1). Maximum degree n' of Benm to evaluate. n'=1 is uniform field vector.
        writeout: boolean (True). Whether to save computed values to disk for rapid replotting.
        path: string (None). Path relative to run directory to print output data to. Defaults to "<install_loc>/MoonMag/induced/".
        bodyname: string (None). Body name to include in writeout filename.
        verbose: boolean (True). Whether to print progress updates to the terminal.
        append: string (""). Optional string appended to default file names.
        debug: boolean (False). Special use flag.
        do_parallel: boolean (True). Toggle for running certain calculations in parallel.
        outFname: string (None). Output filename to use when writeout = True.
        outFnameS: string (None). As above, for output Gauss coefficients in the Schmidt normalization.
        Xid: complex (None, shape ...). Option to pass in Xid to avoid needing to reload from disk or recalculate.
    """
def BiList(r_bds, sigmas, peak_omegas, asym_shape_layers, grav_shape, Benm, rscale_moments, nvals, mvals, p_max, nprm_max=1, writeout=True, path=None, bodyname=None,
           verbose=True, append="", debug=False, do_parallel=True, Schmidt=False, outFname=None, outFnameS=None, Xid=None):

    # Clean inputs and initialize
    if not isinstance(peak_omegas, Iterable):
        peak_omegas = [peak_omegas]
    n_peaks = np.size(peak_omegas)
    Nnm = np.size(nvals)
    n_max = nprm_max + p_max
    Binms = np.zeros((n_peaks, 2, n_max+1, n_max+1), dtype=np.complex128)
    log.debug(f"Calculating asymmetric B_inm for {np.size(peak_omegas)} periods.")
    if writeout:
        if bodyname is None:
            bfname = ""
        else:
            bfname = f"{bodyname}_"
    if rscale_moments == 1.0:
        rscaling = rscale_moments
    else:
        rscaling  = rscale_moments * r_bds[-1]
    if grav_shape is None:
        asym_shape = asym_shape_layers
    else:
        asym_shape = asym_shape_layers + grav_shape

    # Get mixing coefficients
    if Xid is None:
        Xid = get_all_Xid(nprm_max, p_max, nprm_max+p_max, nvals, mvals, reload=True, do_parallel=do_parallel, fpath=path)

    if do_parallel and not debug:
        par_kw = {'nprm_max':nprm_max, 'verbose':verbose}
        # For each omega, evaluate Bi:
        pool = mtpContext.Pool(np.minimum(num_cores,n_peaks))
        par_result = [pool.apply_async( BinmResponse, (r_bds,sigmas,peak_omegas[i_om],asym_shape,Benm[i_om,...],Xid,p_max,rscaling), par_kw ) for i_om in range(n_peaks)]
        pool.close()
        pool.join()

        for i_om in range(n_peaks):
            Binms[i_om, ...] = par_result[i_om].get()
    else:
        if debug:
            Aes, Ats, Ads = [ np.zeros((n_peaks, nprm_max+p_max), dtype=np.complex128) for _ in range(3) ]
            krvals = np.zeros(n_peaks, dtype=np.complex128)
            for i_om in range(n_peaks):
                Binms[0, ...], Aes[i_om, :], Ats[i_om, :], Ads[i_om, :], krvals[i_om] = BinmResponse(r_bds, sigmas, peak_omegas[i_om], asym_shape, Benm[0, ...], Xid, p_max, rscaling, nprm_max=nprm_max, verbose=verbose, debug=debug)
                if (i_om+1) % 10 == 0: log.debug(f"{i_om + 1} of {n_peaks} complete.")
        else:
            for i_om in range(n_peaks):
                Binms[i_om, ...] = BinmResponse(r_bds, sigmas, peak_omegas[i_om], asym_shape, Benm[i_om, ...], Xid, p_max, rscaling, nprm_max=nprm_max, verbose=verbose)
                log.debug(f"{i_om + 1} of {n_peaks} complete.")

    if writeout:
        if path is None:
            path = _induced
        if outFname is None:
            outFname = f'{bfname}Binm_asym{append}'
        fpath = os.path.join(path, f"{outFname}.dat")
        fout = open(fpath, "w")
        header = "{:<13}, {:<4}, {:<4}, {:<24}, {:<24}\n".format("Period (hr) ", "n ", "m ", "Binm_Re (nT)", "Binm_Im (nT)")
        fout.write(header)
        for i in range(np.size(peak_omegas)):
            T_hrs = 2*np.pi/peak_omegas[i]/3600
            for i_nm in range(Nnm):
                sign = int(mvals[i_nm]<0)
                this_Binm = Binms[i,sign,nvals[i_nm],abs(mvals[i_nm])]
                fout.write( "{:<13}, {:<4}, {:<4}, {:<24}, {:<24}\n".format(round(T_hrs,5), nvals[i_nm], mvals[i_nm], np.real(this_Binm), np.imag(this_Binm)) )
        fout.close()
        log.info(f"Data for asymmetric Binm written to file: {fpath}")

        if Schmidt:
            if outFnameS is None:
                outFnameS = f'{bfname}ghnm_asym{append}'
            fpath = os.path.join(path, f"{outFnameS}.dat")
            fout = open(fpath, "w")
            header = "{:<13}, {:<4}, {:<4}, {:<24}, {:<24}, {:<24}, {:<24}\n".format("Period (hr) ", "n ", "m ", "g_nm_Re (nT)", "g_nm_Im (nT)", "h_nm_Re (nT)", "h_nm_Im (nT)")
            fout.write(header)
            for i in range(np.size(peak_omegas)):
                T_hrs = 2*np.pi/peak_omegas[i]/3600
                this_gnm, this_hnm = get_gh_from_Binm(n_max,Binms[i,...])
                for n in range(1,n_max+1):
                    for m in range(n+1):
                        fout.write( "{:<13}, {:<4}, {:<4}, {:<24}, {:<24}, {:<24}, {:<24}\n".format(round(T_hrs,5), n, m, np.real(this_gnm[n,m]), np.imag(this_gnm[n,m]), np.real(this_hnm[n,m]), np.imag(this_hnm[n,m])) )
            fout.close()
            log.info(f"Data for asymmetric, Schmidt semi-normalized g_nm and h_nm written to file: {fpath}")

    if debug:
        return Binms, Aes, Ats, Ads, krvals
    else:
        return Binms

#############################################

"""
BinmResponse()
    The induced magnetic moments Binm for all input degree n and order m,
    with a specified asymmetric boundary shape asym_shape and for the input
    Benm at a single frequency omega.
    Usage: `Binm` = BinmResponse(`r_bds`, `sigmas`, `omega`, `asym_shape`, `Benm`, `Xid`, `p_max`, `rscaling`, `nprm_max=1`, `verbose=True`, `debug=False`)
    Returns:
        Binm: complex, shape(2,n_max+1,n_max+1). The complex induced moments for the given Benm and asymmetric interior structure,
            evaluated for a single angular frequency of oscillation omega.
    Parameters:
        r_bds: float, shape(n_bds). Radii for each boundary surface in m.
        sigmas: float, shape(n_bds). Conductivity for each layer under the corresponding boundary in S/m.
        omega: float. Angular frequency of magnetic oscillations in rad/s.
        asym_shape: complex, shape(n_bds,2,p_max+1,p_max+1). A list of boundary shape parameters chi_pq for each boundary.
        Benm: complex, shape(2,n_max+1,n_max+1). Excitation moments for each n' and m', padded with zeros for n > n'.
        Xid: mpf, shape(2,nprm_max+1,nprm_max+1, 2,p_max+1,p_max+1, 2,n_max+1,n_max+1). Mixing coefficients
            resulting from multiplication of spherical harmonics. \Xi^{\star nm}_{n'm'pq} from the paper.
        p_max: integer. Largest degree p in interior model boundary shapes.
        rscaling: float. Ratio of outermost conducting boundary to body radius.
        nprm_max: integer (1). Maximum degree n' of Benm to evaluate. n'=1 is uniform field vector.
        verbose: boolean (True). Whether to print progress updates to the terminal.
        debug: boolean (False). Special use flag.
    """
def BinmResponse(r_bds, sigmas, omega, asym_shape, Benm, Xid, p_max, rscaling, nprm_max=1, verbose=True, debug=False):
    n_max = nprm_max + p_max
    Binm = np.zeros((2, n_max+1, n_max+1), dtype=np.complex128)
    n_bds = np.size(r_bds)
    n_inner = n_bds - 1

    ni_bds = n_bds - 1
    ni_inner = n_inner - 1
    lone = one
    czero = mp.mpc(0)

    # Get type object for initializing numpy arrays of mpc
    mpc = mpc_global

    # Convert input data for k values into mpmath floats
    # When we won't need to explicitly loop over indices, use numpy
    # arrays to enable vectorized computation.
    n = np.array([ mp.mpf(ni) for ni in range(1,n_max+1) ])
    if not isinstance(r_bds, Iterable):
        r_bds = [mp.mpf(r_bds)]
        sigmas = [mp.mpf(sigmas)]
    else:
        r_bds = [mp.mpf(r_i) for r_i in r_bds]
        sigmas = [mp.mpf(sig_i) for sig_i in sigmas]

    if not debug:
        omega = mp.mpf(omega)

        # Make lists of k values
        kr_ll = np.array([ mp.sqrt(j*omega*mu_o*sigmas[i_bdy])*r_bds[i_bdy] for i_bdy in range(n_inner) ])
        kr_ul = np.array([ mp.sqrt(j*omega*mu_o*sigmas[i_bdy+1])*r_bds[i_bdy] for i_bdy in range(n_inner) ])
        # ka value at outer boundary
        kr_ex = mp.sqrt(j*omega*mu_o*sigmas[-1])*r_bds[-1]
    else:
        n_bds = 1
        # For debug purposes, omega is passed in as a kr value
        kr_ex = mp.sqrt(j) * mp.mpf(omega)

    # Insert kr into Bessel functions for outer bdy only first
    jn_ex = np.zeros(n_max, dtype=mpc)
    jd_ex = np.zeros(n_max, dtype=mpc)
    yn_ex = np.zeros(n_max, dtype=mpc)
    yd_ex = np.zeros(n_max, dtype=mpc)
    for ni in range(n_max):
        jn_ex[ni] = jnx(n[ni],[kr_ex])[0]
        jd_ex[ni] = jdx(n[ni],[kr_ex])[0]
        yn_ex[ni] = ynx(n[ni],[kr_ex])[0]
        yd_ex[ni] = ydx(n[ni],[kr_ex])[0]

    # Calculate outer boundary transfer quantities.
    # All such quantities (here and more below) have Q appended---the "Q response" is analogous to \mathcal{A}.
    alphQ = 1/kr_ex
    betaQ = jd_ex - (n+one)*jn_ex
    gammQ = yd_ex - (n+one)*yn_ex
    deltQ = n*jn_ex + jd_ex
    epsiQ = n*yn_ex + yd_ex
    xi__Q = -kr_ex**2 * jn_ex[:nprm_max]
    rho_Q = -kr_ex**2 * yn_ex[:nprm_max]
    zetaQ = n*jn_ex
    eta_Q = n*yn_ex

    log.debug("    Initialized transfer calculations...")

    # If there's only one boundary, kr_ll and kr_ul are empty and we can immediately calculate A function values.
    if n_bds > 1:
        # Initialize Bessel function results
        jn_ll, jd_ll, yn_ll, yd_ll = ( np.zeros((n_max,n_inner), dtype=mpc) for _ in range(4) )
        jn_ul, jd_ul, yn_ul, yd_ul = ( np.zeros((n_max,n_inner), dtype=mpc) for _ in range(4) )

        # Continue inserting kr into Bessel functions
        for ni in range(n_max):
            jn_ll[ni,:] = jnx(n[ni],kr_ll)
            jd_ll[ni,:] = jdx(n[ni],kr_ll)
            yn_ll[ni,:] = ynx(n[ni],kr_ll)
            yd_ll[ni,:] = ydx(n[ni],kr_ll)

            jn_ul[ni,:] = jnx(n[ni],kr_ul)
            jd_ul[ni,:] = jdx(n[ni],kr_ul)
            yn_ul[ni,:] = ynx(n[ni],kr_ul)
            yd_ul[ni,:] = ydx(n[ni],kr_ul)

        # Calculate layer transfer quantities
        alph = np.zeros((n_max,n_inner), dtype=mpc)
        alph[:,:-1] = jn_ll[:,1:]*yd_ll[:,1:] - jd_ll[:,1:]*yn_ll[:,1:]
        alph[:,-1] = jn_ex*yd_ex - jd_ex*yn_ex
        beta = jn_ll*yd_ul - yn_ul*jd_ll
        gamm = yn_ll*yd_ul - yn_ul*yd_ll
        delt = jn_ul*jd_ll - jn_ll*jd_ul
        epsi = jn_ul*yd_ll - yn_ll*jd_ul

        # Make the spherically symmetric recursion:
        Lambd = eval_inner_recur_sym(n_max, n_bds, beta, gamm, delt, epsi)

        # Evaluate terms for the layer series:
        TransferQts1 = (jn_ll, yn_ll, jn_ul, yn_ul, jd_ll, yd_ll, kr_ex, jn_ex, yn_ex, alph, beta, gamm, delt, epsi, deltQ, epsiQ, Lambd)
        At, Kn, aBar = get_series_terms(n_max, nprm_max, n, n_bds, TransferQts1, r_bds[-1], Benm)

        TransferQts2 = (kr_ll, kr_ul, jn_ll, yn_ll, Lambd)

    else:
        Lambd = np.zeros((n_max,1), dtype=mpc)
        At = np.zeros((n_max,1), dtype=mpc)
        Kn = np.ones((n_max,1), dtype=mpc)
        aBar = np.zeros((2,nprm_max+1,nprm_max+1,1), dtype=mpc)  # This will not be used
        At[:,0] = jn_ex/deltQ
        TransferQts2 = (kr_ex, kr_ex, jn_ex, yn_ex, Lambd)  # These will not be used

    Ad = (xi__Q + Lambd[:nprm_max,ni_bds]*rho_Q) / (deltQ[:nprm_max] + Lambd[:nprm_max,ni_bds]*epsiQ[:nprm_max])

    log.debug("    Calculating series results...")

    # Finally, evaluate the asymmetry series terms:
    Delta = get_Deltanmi(n_max, p_max, nprm_max, n, n_bds, r_bds, Benm, asym_shape, Xid, aBar, Ad, TransferQts2)

    # And the series itself:
    AKDseries = np.zeros((2,n_max+1,n_max+1), dtype=mpc)
    for i in range(n_bds):
        for nn in range(1,n_max+1):
            ni = nn - 1
            AKDseries[:,nn,:] += At[ni,i] * Kn[ni,i] * Delta[:,nn,:,i]

    if verbose:
        T_hrs = round(2*np.pi/omega/3600, 2)
        log.debug(f"    Evaluating final products for T = {T_hrs}...")

    # Finally, zip it all together:
    Ae = cpx_div( (betaQ + Lambd[:,-1] * gammQ), (deltQ + Lambd[:,-1] * epsiQ) )
    for nn in range(1,n_max+1):
        ni = nn - 1
        nv = n[ni]
        Binm[:,nn,:] = (nv/(nv+lone) * Ae[ni]*Benm[:,nn,:] + nv*AKDseries[:,nn,:]) * rscaling**(nn+2)

    if debug:
        return Binm, Ae, At, Ad, kr_ex
    else:
        return Binm

#############################################

"""
get_series_terms()
    Calculate quantities that appear in sums over layers.
    Usage: `At`, `Kn`, `aBar` = get_series_terms(`n_max`, `nprm_max`, `n`, `n_bds`, `TransferQts`, `R`, `Benm`)
    Returns:
        At: mpc, shape(n_max,n_bds). Scaling amplitude determining the strength of induced field from each layer.
            This is \mathcal{A}_n^{t,i} in the paper.
        Kn: mpc, shape(n_max,n_bds). Transfer product propagating induced moments resulting
            from asymmetry from the asymmetric layer to the outermost boundary. This is K_n^i in the paper.
        aBar: mpc, shape(2,nprm_max+1,nprm_max+1,n_bds). Solutions for the Bessel function coefficients
            a^i_{n'm'} for the case of spherical symmetry. This is \overline{a}^\mathrm{i}_{n'm'} in the paper.
    Parameters:
        n_max: integer. Maximum degree n of induced moments.
        nprm_max: integer. Maximum degree n' of excitation moments.
        n: mpf, shape(n_max). Integer values of n in mpmath float format.
        n_bds: integer. Number of boundaries present in the interior model.
        TransferQts: Tuple of Bessel functions and various quantities using them, all of type mpc. See BinmResponse.
        R: float. Radius of the outermost conducting boundary in m.
        Benm: mpc, shape(2,nprm_max+1,nprm_max+1). Complex amplitudes of excitation moments for this period.
    """
def get_series_terms(n_max, nprm_max, n, n_bds, TransferQts, R, Benm):
    mpc = mpc_global
    lone = one
    ltwo = two
    n_inner = n_bds - 1
    jn_ll, yn_ll, jn_ul, yn_ul, jd_ll, yd_ll, kr_ex, jn_ex, yn_ex, alph, beta, gamm, delt, epsi, deltQ, epsiQ, Lambd = TransferQts

    At, Kn = ( np.ones((n_max,n_bds), dtype=mpc) for _ in range(2) )
    aBar = np.zeros((2,nprm_max+1,nprm_max+1,n_bds), dtype=mpc)

    denom = deltQ + Lambd[:,-1] * epsiQ

    aBarBase = -R*Benm[:,:nprm_max+1,:nprm_max+1]
    for nnp in range(1, nprm_max+1):
        nv = n[nnp-1]
        for mmp in range(-nnp,nnp+1):
            mpsign = int(mmp<0)
            mpabs = abs(mmp)
            aBar[mpsign,nnp,mpabs,-1] = cpx_div_val( aBarBase[mpsign,nnp,mpabs] * (ltwo*nv + lone)/(nv + lone), denom[nnp-1] )

    At[:,-1] = n * (jn_ex + Lambd[:,-1]*yn_ex) / denom
    for i in range(n_inner-1, -1, -1):
        At[:,i] = cpx_div( (jn_ll[:,i] + Lambd[:,i]*yn_ll[:,i]), denom)

        Kn[:,i] = cpx_div( (Kn[:,i+1] * alph[:,i]), (beta[:,i] + Lambd[:,i]*gamm[:,i]) )

        for nnp in range(1,nprm_max+1):
            npi = nnp -1
            for mmp in range(-nnp,nnp+1):
                mpsign = int(mmp < 0)
                mpabs = abs(mmp)
                aBar[mpsign,nnp,mpabs,i] = cpx_div_val( (aBar[mpsign,nnp,mpabs,i+1] * (jn_ul[npi,i] + Lambd[npi,i+1]*yn_ul[npi,i])), (jn_ll[npi,i] + Lambd[npi,i]*yn_ll[npi,i]) )

    return At, Kn, aBar

#############################################

"""
get_Deltanmi()
    Calculate Delta, which represents the amount of "mixing" from excitation moments into
        induced moments of other n,m due to asymmetry.
    Usage: `Delta` = get_Deltanmi(`n_max`, `p_max`, `nprm_max`, `n`, `n_bds`, `r_bds`, `Benm`, `asym_shape`, `Xid`, `aBar`, `Ad`, `TransferQts`)
    Returns:
        Delta: mpc, shape(2,n_max+1,n_max+1,n_bds). This is \Delta^\mathrm{i}_{nm} in the paper.
    Parameters:
        n_max: integer. Maximum degree n of induced moments.
        p_max: integer. Largest degree p in interior model boundary shapes.
        nprm_max: integer. Maximum degree n' of excitation moments.
        n: mpf, shape(n_max). Integer values of n in mpmath float format.
        n_bds: integer. Number of boundaries present in the interior model.
        r_bds: float, shape(n_bds). Radii for each boundary surface in m.
        Benm: complex, shape(2,n_max+1,n_max+1). Complex excitation moments for this period
            of oscillation. This is B^{e}_{nm} in the paper.
        asym_shape: complex, shape(n_bds,2,p_max+1,p_max+1). A list of boundary shape parameters chi_pq for each boundary.
        Xid: mpf, shape(2,nprm_max+1,nprm_max+1, 2,p_max+1,p_max+1, 2,n_max+1,n_max+1). Mixing coefficients
            resulting from multiplication of spherical harmonics. \Xi^{\star nm}_{n'm'pq} from the paper.
        aBar: mpc, shape(n_bds,2,n_max+1,n_max+1). Complex coefficients for the internal magnetic field from the
            spherically symmetric solutions. This is \overline{a}^\mathrm{i}_{n'm'} from the paper.
        Ad: mpc, shape(nprm_max). Complex amplitude mixing from the excitation harmonics into the 
            induced harmonics. This is \mathcal{A}^{\star}_{n'} from the paper.
    """
def get_Deltanmi(n_max, p_max, nprm_max, n, n_bds, r_bds, Benm, asym_shape, Xid, aBar, Ad, TransferQts):
    mpc = mpc_global
    lone = one
    ltwo = two
    Delta = np.zeros((2,n_max+1,n_max+1,n_bds), dtype=mpc)
    kr_ll, kr_ul, jn_ll, yn_ll, Lambd = TransferQts
    mixsum_init = np.zeros((2,nprm_max+1,nprm_max+1, 2,n_max+1,n_max+1, n_bds), dtype=mpc)

    for nn in range(1,n_max+1):
        ni = nn - 1
        nv = n[ni]
        for mm in range(-nn,nn+1):
            msign = int(mm<0)
            mabs = abs(mm)

            for nnp in range(1,nprm_max+1):
                npi = nnp - 1
                for mmp in range(-nnp,nnp+1):
                    mpsign = int(mmp<0)
                    mpabs = abs(mmp)

                    mixsum = mixsum_init + 0

                    for pp in range(1,p_max+1):
                        ppi = pp - 1
                        for qq in range(-pp,pp+1):
                            qsign = int(qq<0)
                            qabs = abs(qq)

                            mixsum[mpsign,nnp,mpabs, msign,nn,mabs, :] += asym_shape[:,qsign,pp,qabs]/r_bds[:] * \
                                                                          Xid[mpsign,nnp,mpabs, qsign,pp,qabs, msign,nn,mabs]

                    for i in range(n_bds-1):
                        Delta[msign,nn,mabs, i] += mixsum[mpsign,nnp,mpabs, msign,nn,mabs, i]/r_bds[-1] * \
                                                   aBar[mpsign,nnp,mpabs, i] * (jn_ll[npi,i] + Lambd[npi,i]*yn_ll[npi,i]) * (kr_ll[i]**2 - kr_ul[i]**2)

                    Delta[msign,nn,mabs,-1] += mixsum[mpsign,nnp,mpabs, msign,nn,mabs, -1] * \
                                               (ltwo*nv + lone)/(nv + lone) * Benm[mpsign,nnp,mpabs] * Ad[npi]

    return Delta

#############################################

"""
eval_inner_recur_sym()
    Evaluate the spherically symmetric inner recursion relations to obtain LambdBar, which is
    exactly the spherically symmetric solution to the boundary conditions. This quantity appears
    in the asymmetric recursions, and we need to evaluate it for all n (not just the excitations
    that result in the case of spherical symmetry).
    Usage: `LambdBar` = eval_inner_recur_sym(`n_max`, `n_bds`, `beta`, `gamm`, `delt`, `epsi`)
    Returns:
        LambdBar: mpc, shape(n_max,n_bds). This is \overline{\Lambda}^{l}_{n} in the paper.
    Parameters:
        n_max: integer. Maximum degree n of induced moments.
        n_bds: integer. Number of boundaries present in the interior model.
        beta, gamm, delt, epsi: mpc, shape(n_max,n_bds). Recursion quantities derived from
            Bessel functions. These are \beta, \gamma, \delta, and \epsilon from the paper.
            See BinmResponse for more details.
    """
def eval_inner_recur_sym(n_max, n_bds, beta, gamm, delt, epsi):
    Lambd = np.zeros((n_max,n_bds), dtype=mpc_global)
    for i in range(n_bds-1):
        Lambd[:,i+1] = cpx_div((delt[:,i] + Lambd[:,i]*epsi[:,i]) ,
                               (beta[:,i] + Lambd[:,i]*gamm[:,i]))

    return Lambd

# Avoids divide-by-zero errors when Lambd is very close to i (high conductivities)
def cpx_div(a, b):
    nvals = np.size(a)
    amag = [ mp.fabs(ai) for ai in a ]
    aarg = [ mp.atan2(mp.im(ai), mp.re(ai)) for ai in a ]
    bmag = [ mp.fabs(bi) for bi in b ]
    barg = [ mp.atan2(mp.im(bi), mp.re(bi)) for bi in b ]

    div_mag = [ amag[i] / bmag[i] for i in range(nvals) ]
    div_arg = [ aarg[i] - barg[i] for i in range(nvals) ]
    quotient = np.array([ div_mag[i] * mp.exp(j * div_arg[i]) for i in range(nvals) ])
    return quotient
def cpx_div_val(a, b):
    amag = mp.fabs(a)
    aarg = mp.atan2(mp.im(a), mp.re(a))
    bmag = mp.fabs(b)
    barg = mp.atan2(mp.im(b), mp.re(b))

    div_mag = amag / bmag
    div_arg = aarg - barg
    quotient = div_mag * mp.exp(j * div_arg)
    return quotient

#############################################
# NOTE THE NOTATION FOR n' AND n IS REVERSED IN THIS CODE BLOCK COMPARED TO MOST OF THE OTHER CODE
# This is for consistency with the paper and for shorter lines of code.
# In this block, n,m is the excitation harmonic and n',m' is the induced harmonic.
#############################################

"""
get_all_Xid()
    Returns mixing coefficients that result from multiplication of spherical harmonics.
    Usage: `Xid` = get_all_Xid(`n_max`, `p_max`, `nprm_max`, `nprmvals`, `mprmvals`)
    Returns:
        Xid: mpf, shape(2,n_max+1,n_max+1, 2,p_max+1,p_max+1, 2,nprm_max+1,nprm_max+1). Amount of Ynm that results
            from multiplying Ynm by Ypq. Values for all combinations of n,m,p,q,n',m' are returned.
    Parameters:
        n_max: integer. Maximum degree n of excitation moments.
        p_max: integer. Maximum degree p of boundary shapes.
        nprm_max: integer. Maximum degree nprm of induced moments. Typically nmax + pmax.
        nprmvals: integer, shape((nprm_max+1)**2-1). Flattened array of nprm values.
        mprmvals: integer, shape((nprm_max+1)**2-1). Flattened array of mprm values corresponding to each nrpm above.
    """
def get_all_Xid(n_max, p_max, nprm_max, nprmvals, mprmvals, do_parallel=True,
                writeout=True, reload=False, fpath=None, fname=None):
    if writeout or reload:
        if fpath is None:
            fpath = _induced
        if fname is None:
            fname = f"Xid_values_n{n_max}_p{p_max}_np{nprm_max}"

        fullFile = os.path.join(fpath, fname+".mat")
        fileExists = os.path.isfile(fullFile)
        if reload and not fileExists:
            log.warning(f'Xid file {fullFile} does not exist, but reload is True. Calculating instead.')
            reload = False
            writeout = True

    if reload:
        Xid = loadmat(fullFile)['Xid']
    else:
        Xid = np.zeros((2,n_max+1,n_max+1, 2,p_max+1,p_max+1, 2,nprm_max+1,nprm_max+1), dtype=mpf_global)
        Nnmprm = np.size(nprmvals)

        for n in range(1,n_max+1):
            for m in range(-n,n+1):
                msign = int(m<0)
                mabs = abs(m)

                for p in range(1,p_max+1):
                    for q in range(-p,p+1):
                        qsign = int(q<0)
                        qabs = abs(q)

                        if do_parallel:
                            pool = mtpContext.Pool(num_cores)
                            par_result = [pool.apply_async( calc_Xid, args=(n,m,p,q,nprmvals[iN],mprmvals[iN],nprm_max) ) for iN in range(Nnmprm)]
                            pool.close()
                            pool.join()

                        for iN in range(Nnmprm):
                            nprm = nprmvals[iN]
                            mpsign = int(mprmvals[iN]<0)
                            mpabs = abs(mprmvals[iN])
                            if do_parallel:
                                Xid[msign,n,mabs, qsign,p,qabs, mpsign,nprm,mpabs] = par_result[iN].get()
                            else:
                                Xid[msign, n, mabs, qsign, p, qabs, mpsign, nprm, mpabs] = calc_Xid(n,m,p,q,nprmvals[iN],mprmvals[iN],nprm_max)
        if writeout:
            savemat(fullFile, {'Xid': Xid.astype(np.float64)})
            log.debug(f'Saved Xid values to file: {fullFile}')

    return Xid

"""
print_Xid_table()
    Prints a table to the terminal of values for all Xid for debug purposes.
    Usage: print_Xid_table(`Xid`, `n_max`, `p_max`, `nprm_max`)
    Returns:
        None.
    Parameters:
        Xid: mpf, shape(2,n_max+1,n_max+1, 2,p_max+1,p_max+1, 2,nprm_max+1,nprm_max+1). Xid values to print.
        n_max: integer. Maximum degree n of excitation moments.
        p_max: integer. Maximum degree p of boundary shapes.
        nprm_max: integer. Maximum degree nprm of induced moments. Typically nmax + pmax.
    """
def print_Xid_table(Xid, n_max, p_max, nprm_max):
    # Print table header
    head_row = "    "
    for n in range(1, n_max+1):
        for m in range(-n, n+1):
            head_row = head_row + f"|    Y{n}{m}    ".ljust(25)
    head_row = f"{head_row}|"
    print(head_row)
    # Print shapes as rows
    for p in range(1, p_max+1):
        for q in range(-p, p+1):
            qsign = int(q<0)
            qabs = abs(q)
            row_str = f"S{p}{q}".ljust(4)

            # Print excitation moments as columns
            for n in range(1, n_max+1):
                for m in range(-n, n+1):
                    msign = int(m<0)
                    mabs = abs(m)
                    cell_str = "| "

                    # Print non-zero terms for each cell
                    for nprm in range(nprm_max, 0, -1):
                        mprm = m + q
                        mpsign = int(mprm<0)
                        mpabs = abs(mprm)

                        this_Xid = Xid[msign,n,mabs, qsign,p,qabs, mpsign,nprm,mpabs]
                        abs_Xi = abs(this_Xid)
                        round_Xi = round(abs(this_Xid),3)
                        if(abs_Xi > 1e-18):
                            if(cell_str != "| "):
                                if(int(this_Xid>0) or round_Xi==0.0):
                                    cell_str = f"{cell_str}+ "
                                else:
                                    cell_str = f"{cell_str}- "
                            else:
                                if(int(this_Xid>0) or round_Xi==0.0):
                                    cell_str = f"{cell_str} "
                                else:
                                    cell_str = f"{cell_str}-"
                            cell_str = f"{cell_str}{round_Xi}Y{nprm}{mprm} ".ljust(13)

                    row_str = f"{row_str}{cell_str.ljust(25)}"

            row_str = f"{row_str}|"
            print(row_str)

    return

"""
print_Xi_table()
    Prints a table to the terminal of values for all Xi for debug purposes.
    Usage: print_Xi_table(`n_max`, `p_max`, `nprm_max`)
    Returns:
        None.
    Parameters:
        n_max: integer. Maximum degree n of excitation moments.
        p_max: integer. Maximum degree p of boundary shapes.
        nprm_max: integer. Maximum degree nprm of induced moments. Typically nmax + pmax.
    """
def print_Xi_table(n_max, p_max, nprm_max):
    # Print table header
    head_row = "    "
    for n in range(1, n_max+1):
        for m in range(-n, n+1):
            head_row = f"{head_row}|    Y{n}{m}    ".ljust(25)
    head_row = f"{head_row}|"
    print(head_row)
    # Print shapes as rows
    for p in range(1, p_max+1):
        for q in range(-p, p+1):
            row_str = f"S{p}{q}".ljust(4)

            # Print excitation moments as columns
            for n in range(1, n_max+1):
                for m in range(-n, n+1):
                    cell_str = "| "

                    # Print non-zero terms for each cell
                    for nprm in range(nprm_max, 0, -1):
                        mprm = m + q

                        this_Xi = calc_Xi(n,m,p,q,nprm,mprm)
                        if abs(this_Xi) > 0.00001:
                            if cell_str != "| ":
                                if(int(this_Xi<0)):
                                    cell_str = f"{cell_str}- "
                                else:
                                    cell_str = f"{cell_str}+ "
                            else:
                                if(int(this_Xi<0)):
                                    cell_str = f"{cell_str}-"
                                else:
                                    cell_str = f"{cell_str} "
                            cell_str = f"{cell_str}{abs(this_Xi):.3f}Y{nprm}{mprm} ".ljust(13)

                    row_str = f"{row_str}{cell_str}".ljust(25)

            row_str = f"{row_str}|"
            print(row_str)

    return

"""
The remaining functions in this block are all steps and substeps in calculating Xid
    as in equations S91-S103 from the paper.
    """
def calc_Xid(n,m,p,q,nprm,mprm, nprm_max):
    this_Xid_num =        calc_Xiw(n, m, p, q, nprm, mprm) +  \
                    XidSeriesLower(n, m, p, q, nprm, mprm) +  \
                    XidSeriesUpper(n, m, p, q, nprm, mprm, nprm_max)
    if this_Xid_num != 0:
        Xid = this_Xid_num / XidDenom(nprm,mprm,nprm_max)
    else:
        Xid = zero
    return Xid

def calc_Xi(n, m, p, q, nprm, mprm):
    triang = (abs(n - p) <= nprm) and (nprm <= n + p)
    even = (n + p + nprm) % 2 == 0
    m_match = (m + q == mprm)

    if triang and even and m_match:
        nu = int((n + p - nprm)/2)
        cs_sign = (-1)**nu

        kap_min = max(0, 2*nu-(n+m), 2*nu-(p-q))
        kap_max = min(2*nu, n-m, p+q)

        norm = np.sqrt((2*n + 1) * (2*p + 1) * (2*nprm + 1) / 4 / np.pi)
        ft_sqrt = np.sqrt( ft(n+m)*ft(n-m)*ft(p+q)*ft(p-q)*ft(nprm+mprm)*ft(nprm-mprm) )
        ft_frac = ft(2*nu)/ft(nu) * ft(2*n-2*nu)/ft(n-nu) * ft(2*p-2*nu)/ft(p-nu) * ft(nprm+nu)/ft(2*nprm+1+2*nu)

        kap_sum = np.sum([ (-1)**kap / ( ft(kap)*ft(2*nu-kap) * ft(n-m-kap)*ft(n+m-(2*nu-kap)) *
            ft(p+q-kap)*ft(p-q-(2*nu-kap)) ) for kap in range(kap_min,kap_max+1) ])

        Xi = cs_sign * norm * ft_frac * ft_sqrt * kap_sum
    else:
        Xi = zero

    return Xi

def w_M(n,m):
    wm = mp.re( (n+1)*mp.sqrt(n**2 - m**2) / mp.sqrt((2*n-1) * (2*n+1)) )
    return wm

def w_P(n,m):
    wp = n*mp.re( mp.sqrt((n+1)**2 - m**2) / mp.sqrt((2*n+1)*(2*n+3)) )
    return wp

def calc_Xiw(n, m, p, q, nprm, mprm):
    mabs = abs(m)
    mpabs = abs(mprm)

    wMn = w_M(n,mabs)
    wPn = w_P(n,mabs)
    wMnp = w_M(nprm,mpabs)
    wPnp = w_P(nprm,mpabs)

    Xiw =  wMn * wMnp * calc_Xi(n-1, m, p, q, nprm-1, mprm)    \
         + wPn * wPnp * calc_Xi(n+1, m, p, q, nprm+1, mprm)    \
         - wMn * wPnp * calc_Xi(n-1, m, p, q, nprm+1, mprm)    \
         - wPn * wMnp * calc_Xi(n+1, m, p, q, nprm-1, mprm)

    return Xiw

def F_M(n, m, twokap):
    if n < 3 or twokap > n-abs(m)-2:
        return zero
    else:
        n_this = n-twokap
        n_down = n-(twokap+2)
        n_down2 = n-(twokap+4)
        if w_M(n_this,m)==0 or w_P(n_down,m)==0:
            return zero
        else:
            Fm = w_M(n_this,m)*w_P(n_down,m) /    \
            ( w_M(n_down,m)**2 + w_P(n_down,m)**2 -
                w_M(n_down,m)*w_P(n_down2,m)*F_M(n,m,twokap+2) )

    return Fm

def F_P(n, m, twokap, nprm_max):
    if twokap > nprm_max-n:
        return zero
    else:
        n_this = n+twokap
        n_up = n+(twokap+2)
        n_up2 = n+(twokap+4)
        Fp = w_P(n_this,m)*w_M(n_up,m) /    \
        ( w_M(n_up,m)**2 + w_P(n_up,m)**2 -
            w_P(n_up,m)*w_M(n_up2,m)*F_P(n,m,twokap+2,nprm_max) )

    return Fp

# Note that although this term is called SeriesLower, it contains mixing terms from higher degree n'' into the lower n'
def XidSeriesLower(n, m, p, q, nprm, mprm):
    if nprm >= 3:
        gmax_l = floor((nprm-abs(mprm))/2)
    else:
        gmax_l = 0
    thisLower = 0
    for g in range(1,gmax_l+1):
        thisProd = 1
        for kap in range(g):
            thisProd = thisProd * F_M(nprm,mprm,2*kap)

        thisLower = thisLower + calc_Xiw(n, m, p, q, nprm-2*g, mprm)*thisProd

    seriesLower = thisLower
    return seriesLower

# Note that although this term is called SeriesUpper, it contains mixing terms from lower degree n'' into the higher n'
def XidSeriesUpper(n, m, p, q, nprm, mprm, nprm_max):
    gmax_u = floor((nprm_max-nprm)/2)
    thisUpper = 0
    for g in range(1,gmax_u+1):
        thisProd = 1
        for kap in range(g):
            thisProd = thisProd * F_P(nprm,mprm,2*kap,nprm_max)

        thisUpper = thisUpper + calc_Xiw(n, m, p, q, nprm+2*g, mprm)*thisProd

    seriesUpper = thisUpper
    return seriesUpper

def XidDenom(n, m, nprm_max):
    XidDenom = w_M(n,m)**2 + w_P(n,m)**2 - w_M(n,m)*w_P(n-2,m)*F_M(n,m,0) - w_P(n,m)*w_M(n+2,m)*F_P(n,m,0,nprm_max)
    return XidDenom

#############################################
# -- END REVERSED n, n' NOTATION --
#############################################

"""
eval_dev()
    Evaluates the deviation from spherical symmetry as a function of lat/lon for a given degree p and order q of boundary shape.
    Values for chi_pq are coefficients that multiply the associated orthonormal spherical harmonic and then
    are added to the nominal radius. Returns same units as input chi_pq.
    Usage: `this_devs` = eval_dev(`p`, `q`, `chi_pq`, `ltht`, `lphi`, `lleny`, `llenx`)
    Returns:
        this_devs: float, shape(lleny,llenx). Boundary shape deviations due to this particular p,q combination.
    Parameters:
        p: integer. Degree of shape harmonic to be evaluated.
        q: integer. Order of shape harmonic to be evaluated.
        chi_pq: complex. Coefficient for spherical harmonics of degree and order p,q,
            such that r(theta, phi) = R + Sum[chi_pq * Ypq(theta, phi)].
        ltht: float, shape(lleny) or shape(l1D). Local copy of theta grid values (faster execution than referencing a global variable).
        lphi: float, shape(llenx) or shape(l1D). Local copy of phi grid values.
    """
def eval_dev(p, q, chi_pq, ltht, lphi, outShape):
    if chi_pq == 0:
        this_devs = np.zeros(outShape, dtype=np.float64)
    else:
        this_devs = np.real(chi_pq * sph_harm(q, p, lphi, ltht))
    log.debug(f"p,q = {p}{q} completed")
    return this_devs

#############################################

"""
get_rsurf()
    Calculates r(theta, phi) from r_mean and chi_pq, where r(theta, phi) = r_mean + Sum[chi_pq*Ypq(theta, phi)].
    Usage: get_rsurf(`Binm`, `n_max`, `fpath=None`, `difference=True`, `Binm_sph=None`)
    Returns:
        surf: float, shape(lleny,llenx). r(theta, phi) values corresponding to the input grid of
            ltht and lphi values. Same units as asym_shape and r_mean.
    Parameters:
        pvals, qvals: integer, shape(Npq). Linear arrays of paired p,q values for parallel computation.
        asym_shape: complex, shape(2,p_max+1,p_max+1). Absolute boundary deviations for all p,q values for the
            considered surface. Units must match r_mean.
        r_mean: float. Mean value of radius for the considered surface. Units must match asym_shape.
        ltht: float, shape(lleny). Local copy of theta grid values (faster execution than referencing a global variable).
        lphi: float, shape(lleny). Local copy of phi grid values.
    """
def get_rsurf(pvals,qvals,asym_shape, r_mean,ltht,lphi, do_parallel=True):
    Npq = np.size(pvals)
    lleny = np.size(ltht)
    llenx = np.size(lphi)
    if lleny == llenx:
        # Assume 1D arrays if both ltht and lphi are same length
        outShape = (lleny)
    else:
        outShape = (lleny,llenx)
        ltht, lphi = np.meshgrid(ltht, lphi, indexing='ij')
    devs = np.zeros(outShape)

    lin_bd_shape = np.array([ asym_shape[int(qvals[iN]<0),pvals[iN],abs(qvals[iN])] for iN in range(Npq) ])

    if do_parallel:
        pool = mtpContext.Pool(num_cores)
        par_result = [pool.apply_async( eval_dev, args=(pvals[iN],qvals[iN],lin_bd_shape[iN],ltht,lphi,outShape) ) for iN in range(Npq)]
        pool.close()
        pool.join()

        # Unpack the parallel processing results and sum them together
        for res in par_result:
            devs = devs + res.get()
    else:
        for iN in range(Npq):
            devs = devs + eval_dev(pvals[iN],qvals[iN],lin_bd_shape[iN],ltht,lphi,outShape)

    if r_mean is None:
        r_mean = np.abs(asym_shape[0,0,0]) / np.sqrt(4*np.pi)
    surf = devs + r_mean
    return surf

#############################################

"""
getMagSurf()
    Evaluates the induced magnetic field at the surface for all magnetic moments Binm.
    Usage: `Bx`, `By`, `Bz` = getMagSurf(`nvals`, `mvals`, `Binm`, `r_th_ph`, `ltht`, `lphi`, `nmax_plot=10`, `Schmidt=False`)
    Returns:
        Bx,By,Bz (each): complex, shape(lleny,llenx). A (lleny,llenx) array of field values due to the particular Binm values passed.
            Field values can be obtained for any future time by multiplying each Binm[i,:,:,:] by the corresponding e^-iωt factor.
            t=0 is defined to be the J2000 epoch.
    Parameters:
        nvals: integer, shape(Nnm). Linear list of magnetic moment degrees to be evaluated.
        mvals: integer, shape(Nnm). Linear list of magnetic moment orders to be evaluated, corresponding to nvals of the same index.
        Binm: complex, shape(2,n_max+1,n_max+1) OR shape(Nnm); if Schmidt=True, tuple of (gnm,hnm). Magnetic moment of degree and
            order n,m that can be indexed by the matched entries in nvals and mvals. Units match the output field.
        r_th_ph: float, shape(lleny,llenx). Meshgrid array of r(theta, phi) values corresponding to the following tht and phi values.
            r_th_ph has units of R_body, i.e. the physical surface is 1.0.
        ltht: float, shape(lleny). Array of theta values over which to evaluate the field components.
        lphi: float, shape(llenx). Array of phi values over which to evaluate the field components.
        nmax_plot: integer (4). Maximum value of n for evaluating magnetic fields. eval_Bi must have each n explicitly hard-coded,
            which has only been done up to n=4 because larger-degree moments are expected to have small contributions at altitude.
        Schmidt: boolean (False). Whether input magnetic moments are in Schmidt semi-normalized form without Condon-Shortley
            phase. If False, moments must be in fully normalized form with the Condon-Shortley phase.
        gnm, hnm: complex, shape(n_max+1,n_max+1). Schmidt semi-normalized magnetic moments. Passed as a tuple in Binm.
    """
def getMagSurf(nvals,mvals,Binm, r_th_ph,ltht,lphi, nmax_plot=10, Schmidt=False, do_parallel=True):
    if nmax_plot > 10:
        nmax_plot = 10
        log.warning(f"Evaluation of magnetic fields is supported only up to n={nmax_plot}. nmax_plot has been set to {nmax_plot}.")

    if Schmidt:
        Nnm = min( int((nmax_plot+1)*(nmax_plot+2)/2) - 1, np.size(nvals) )
        gnm, hnm = Binm
        if np.size(np.shape(gnm)) > 1:
            lin_gnm = np.array([gnm[nvals[iN], mvals[iN]] for iN in range(Nnm)])
            lin_hnm = np.array([hnm[nvals[iN], mvals[iN]] for iN in range(Nnm)])
        else:
            lin_gnm = gnm
            lin_hnm = hnm
    else:
        Nnm = min((nmax_plot + 1) ** 2 - 1, np.size(nvals))
        if np.size(np.shape(Binm)) > 2:
            lin_Binm = np.array([ Binm[int(mvals[iN]<0),nvals[iN],abs(mvals[iN])] for iN in range(Nnm) ])
        else:
            lin_Binm = Binm

    lleny = np.size(ltht)
    llenx = np.size(lphi)
    Bx, By, Bz = ( np.zeros((1,lleny*llenx), dtype=np.complex128) for _ in range(3) )

    # If we pass a single value for r(theta, phi) = R, evaluate over a sphere with that radius.
    # Otherwise, evaluate *at* the asymmetric 3D surface.
    if not isinstance(r_th_ph, Iterable):
        x = np.array([ r_th_ph*np.sin(thti)*np.cos(phii) for thti in ltht for phii in lphi ])
        y = np.array([ r_th_ph*np.sin(thti)*np.sin(phii) for thti in ltht for phii in lphi ])
        z = np.array([ r_th_ph*np.cos(thti) for thti in ltht for _ in lphi ])
        r = np.zeros((1,lleny*llenx)) + r_th_ph
    else:
        x = np.array([ r_th_ph[i_th,i_ph]*np.sin(ltht[i_th])*np.cos(lphi[i_ph]) for i_th in range(lleny) for i_ph in range(llenx) ])
        y = np.array([ r_th_ph[i_th,i_ph]*np.sin(ltht[i_th])*np.sin(lphi[i_ph]) for i_th in range(lleny) for i_ph in range(llenx) ])
        z = np.array([ r_th_ph[i_th,i_ph]*np.cos(ltht[i_th]) for i_th in range(lleny) for i_ph in range(llenx) ])
        r = np.reshape(r_th_ph, (1,lleny*llenx))

    if do_parallel:
        if Schmidt:
            pool = mtpContext.Pool(num_cores)
            par_result = [pool.apply_async( eval_Bi_Schmidt, args=(nvals[iN],mvals[iN],lin_gnm[iN],lin_hnm[iN], x,y,z,r) ) for iN in range(Nnm)]
            pool.close()
            pool.join()
        else:
            pool = mtpContext.Pool(num_cores)
            par_result = [pool.apply_async( eval_Bi, args=(nvals[iN],mvals[iN],lin_Binm[iN], x,y,z,r) ) for iN in range(Nnm)]
            pool.close()
            pool.join()
        # Unpack results from parallel processing and sum them
        for res in par_result:
            this_Bx, this_By, this_Bz = res.get()
            Bx = Bx + this_Bx
            By = By + this_By
            Bz = Bz + this_Bz
    else:
        for iN in range(Nnm):
            if Schmidt:
                this_Bx, this_By, this_Bz = eval_Bi_Schmidt(nvals[iN], mvals[iN], lin_gnm[iN], lin_hnm[iN], x, y, z, r)
            else:
                this_Bx, this_By, this_Bz = eval_Bi(nvals[iN],mvals[iN],lin_Binm[iN], x,y,z,r)
            Bx = Bx + this_Bx
            By = By + this_By
            Bz = Bz + this_Bz

    Bx = np.reshape(Bx, (lleny,llenx))
    By = np.reshape(By, (lleny,llenx))
    Bz = np.reshape(Bz, (lleny,llenx))
    return Bx, By, Bz


def norm4pi(n, m):
    """
    Calculate normalization factor for 4pi-normalized spherical harmonics.
    Args:
        n: int. Degree of spherical harmonic.
        m: int. Order of spherical harmonic.

    Returns:
        norm: float. Normalization factor, not including any possible Condon-Shortley phase.
    """
    norm = np.sqrt(2*n + 1)
    if m != 0:
        norm = norm * np.sqrt(2 * ft(n - abs(m)) / ft(n + abs(m)))
    return norm


def spherharm4pi(n, m, ltht, lphi):
    """
    Determines the 4pi-normalized surface spherical harmonic without the Condon-Shortley phase.
    Args:
        n: int. Degree of spherical harmonic to evaluate.
        m: int. Order of spherical harmonic to evaluate.
        ltht: float, shape(lleny). 1D list along theta at which to evaluate.
        lphi: float, shape(llenx). 1D list along phi at which to evaluate.

    Returns:
        cSnm, sSnm: float, shape(lleny, llenx). Real surface spherical harmonics in 4pi normalization
            for input n, m. cSnm is the cos of tesseral part and sSnm is the sin part.
    """
    norm = norm4pi(n, m) * (-1)**m  # Scipy lpmv function (legenp) includes Condon-Shortley phase; this removes it.
    Pnm = norm * legenp(m, n, np.cos(ltht))  # Scipy lpmv also has n and m reversed--m first.
    cmphi = np.cos(m*lphi)
    smphi = np.sin(m*lphi)
    cSnm = np.squeeze(np.array([cmphi * Pnmi for Pnmi in Pnm]))
    sSnm = np.squeeze(np.array([smphi * Pnmi for Pnmi in Pnm]))

    return cSnm + sSnm
