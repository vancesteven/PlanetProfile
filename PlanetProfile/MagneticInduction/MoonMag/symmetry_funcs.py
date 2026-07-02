""" Contains functions for calculating induced magnetic fields
    from spherical conductors. Outputs typically appear as mpc, which is
    the mpmath complex type.
    Developed in Python 3.8 for "A perturbation method for evaluating the
    magnetic field induced from an arbitrary, asymmetric ocean world 
    analytically" by Styczinski et al.
    DOI: 10.1016/j.icarus.2021.114840
Author: M. J. Styczinski, mjstyczi@uw.edu """

import os
import numpy as np
from collections.abc import Iterable
import mpmath as mp
# mpmath is needed for enhanced precision to avoid
# divide-by-zero errors induced by underflow.
import multiprocessing as mtp
import platform
plat = platform.system()
if plat == 'Windows':
    mtpType = 'spawn'
else:
    mtpType = 'fork'
mtpContext = mtp.get_context(mtpType)
num_cores = mtp.cpu_count()

from MoonMag import _induced
from MoonMag.config import *

# Global variables and settings
# Set maximum precision for mpmath quantities
mp.mp.dps = digits_precision
# Numerical constants in high-precision mpmath format
zero = mp.mpf("0")
one = mp.mpf("1")
j = mp.mpc(1j)
mu_o = mp.mpf("4.0e-7")*mp.pi
sqrt2 = np.sqrt(2)
sqrt4pi = np.sqrt(4*np.pi)

#############################################

"""
validate()
    Check inputs to be sure everything will be interpreted correctly.
    Usage: `r_bds`, `sigmas`, `omegas` = validate(`r_bds`, `sigmas`, `omegas`)
    Returns:
        r_bds: float, shape(N).
        sigmas: float, shape(N).
        omegas: float, shape(P).
    Parameters:
        r_bds: float, shape(N).
        sigmas: float, shape(N).
        omegas: float, shape(P).
    """
def validate(r_bds, sigmas, omegas):
    #    Check length of boundary radius and conductivity lists
    if np.shape(r_bds) != np.shape(sigmas):
        log.debug(f"boundaries shape: {np.shape(r_bds)}")
        log.debug(f"sigmas shape: {np.shape(sigmas)}")
        raise ValueError("The number of boundaries is not equal to the number of conductivities.")

    if not isinstance(r_bds, list):
        r_bds = [r_bds]
        sigmas = [sigmas]
    if not isinstance(omegas, list):
        omegas = [omegas]

    return r_bds, sigmas, omegas

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
        x: mpc, shape(N). Complex argument of Bessel function (kr). MUST be a list or this will break. Make single values into lists with x = [x].
    """
def jnx(n,x):
    two = mp.mpf("2")
    jnx = np.array([ mp.besselj( n+mp.mpf("0.5"),xi ) * mp.sqrt(mp.pi/two/xi) for xi in x ])
    return jnx
def ynx(n,x):
    two = mp.mpf("2")
    ynx = np.array([ mp.bessely( n+mp.mpf("0.5"),xi ) * mp.sqrt(mp.pi/two/xi) for xi in x ])
    return ynx
def jdx(n,x):
    one = mp.mpf("1")
    jn = jnx(n,x)
    jP = jnx(n+one,x)
    jdx = np.array([ (n+one)*jn[i] - x[i]*jP[i] for i in range(np.size(x)) ])
    return jdx
def ydx(n,x):
    one = mp.mpf("1")
    yn = ynx(n,x)
    yP = ynx(n+one,x)
    ydx = np.array([ (n+one)*yn[i] - x[i]*yP[i] for i in range(np.size(x)) ])
    return ydx

#############################################

"""
AeResponse()
    The complex response amplitude A^e for degree n excitation field,
    with spherical symmetry and for one value of omega.
    Usage: `Ae` = AeResponse(`r_bds`, `sigmas`, `omega`, `nn=1`)
    Returns:
        Ae: mpc. The complex amplitude response for the given n and interior structure.
    Parameters:
        r_bds: float, shape(N). Radii for each boundary surface in m.
        sigmas: float, shape(N). Conductivity for each layer under the corresponding boundary in S/m.
        omega: float. Angular frequency of magnetic oscillations in rad/s.
        rscaling: float. Ratio of outermost conducting boundary to body radius.
        nn: integer (1). Degree n of A_n^e to evaluate.
    """
def AeResponse(r_bds, sigmas, omega, rscaling, nn=1):
    n_bds = np.size(r_bds)

    # Convert input data into mpmath floats
    n = mp.mpf(nn)
    if not isinstance(r_bds, Iterable):
        r_bds = [ mp.mpf(r_bds) ]
        sigmas = [mp.mpf(sigmas)]
    else:
        r_bds = [mp.mpf(r_i) for r_i in r_bds]
        sigmas = [mp.mpf(sig_i) for sig_i in sigmas]

    omega = mp.mpf(omega)

    # Make lists of k values
    kr_ll = []
    kr_ul = []
    for i_bdy in range(n_bds-1):
        k_lower = mp.sqrt(j*omega*mu_o*sigmas[i_bdy])
        k_upper = mp.sqrt(j*omega*mu_o*sigmas[i_bdy+1])
        # kr values for lower/upper layers at each boundary
        kr_ll.append(k_lower*r_bds[i_bdy])
        kr_ul.append(k_upper*r_bds[i_bdy])
    # ka value at outer boundary
    kr_ex = [mp.sqrt(j*omega*mu_o*sigmas[-1])*r_bds[-1]]

    # Insert kr into Bessel functions for outer bdy only first
    jn_ex = jnx(nn,kr_ex)
    jd_ex = jdx(nn,kr_ex)
    yn_ex = ynx(nn,kr_ex)
    yd_ex = ydx(nn,kr_ex)

    #    Calculate outer boundary transfer quantities
    betaQ =  jd_ex - (n+one)*jn_ex
    gammaQ = yd_ex - (n+one)*yn_ex
    deltaQ = n*jn_ex + jd_ex
    epsQ =   n*yn_ex + yd_ex

    # If there's only one boundary, kr_ll and kr_ul are empty and we can immediately calculate Ae.
    if n_bds == 1:
        Ae = cpx_div(betaQ,deltaQ) * rscaling**(nn+2)
        return Ae

    # Continue inserting kr into Bessel functions
    jn_ll = jnx(nn,kr_ll)
    jd_ll = jdx(nn,kr_ll)
    yn_ll = ynx(nn,kr_ll[1:])
    yd_ll = ydx(nn,kr_ll[1:])
    yn_ll = np.insert(yn_ll, 0, zero)
    yd_ll = np.insert(yd_ll, 0, zero)

    jn_ul = jnx(nn,kr_ul)
    jd_ul = jdx(nn,kr_ul)
    yn_ul = ynx(nn,kr_ul)
    yd_ul = ydx(nn,kr_ul)

    # Initialize quantities for recursion relations
    Lambd = []
    alph = []
    beta = []
    gamm = []
    delt = []
    epsi = []
    C_rel = []
    D_rel = []

    # Calculate layer transfer quantities
    beta.append(jn_ll[0]*yd_ul[0] - yn_ul[0]*jd_ll[0])
    gamm.append(yn_ll[0]*yd_ul[0] - yn_ul[0]*yd_ll[0])
    delt.append(jn_ul[0]*jd_ll[0] - jn_ll[0]*jd_ul[0])
    epsi.append(jn_ul[0]*yd_ll[0] - yn_ll[0]*jd_ul[0])
    Lambd.append(cpx_div_val(delt[0],beta[0]))

    for i in range(1,n_bds-1):
        beta.append(jn_ll[i]*yd_ul[i] - yn_ul[i]*jd_ll[i])
        gamm.append(yn_ll[i]*yd_ul[i] - yn_ul[i]*yd_ll[i])
        delt.append(jn_ul[i]*jd_ll[i] - jn_ll[i]*jd_ul[i])
        epsi.append(jn_ul[i]*yd_ll[i] - yn_ll[i]*jd_ul[i])
        this_Lambd = cpx_div_val( delt[-1] + Lambd[-1]*epsi[-1] , beta[-1] + Lambd[-1]*gamm[-1] )
        Lambd.append(this_Lambd)

    # Finally, find Ae from outer BC, using special beta-eps we eval'd first.
    Ae = cpx_div((betaQ + Lambd[-1]*gammaQ), (deltaQ + Lambd[-1]*epsQ)) * rscaling**(nn+2)
    return Ae

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

"""
InducedAeList()
    The complex response amplitude A^e for degree n excitation field,
    with spherical symmetry and for a list of omega values.
    Usage: `Aes`, `AeM`, `AeA` = InducedAeList(`r_bds`, `sigmas`, `omegas`, `nn=1`, `writeout=False`, `path=None`)
    Returns:
        Aes: mpc, shape(P). List of complex amplitude responses for the given n and interior structure.
        AeM: mpc, shape(P). List of response magnitudes for the given n and interior structure.
        AeA: mpc, shape(P). List of phase arguments for the given n and interior structure.
    Parameters:
        r_bds: float, shape(N). Radii for each boundary surface in m.
        sigmas: float, shape(N). Conductivity for each layer under the corresponding boundary in S/m.
        omegas: float, shape(P). Angular frequencies of magnetic oscillations in rad/s.
        rscale_moments: float. Equal to 1/R_body in units of 1/m. For proper scaling when conducting boundaries
            may or may not coincide with the body surface.
        nn: integer (1). Degree n of A_n^e to evaluate.
        writeout: boolean (False). Whether to save a .txt file of values calculated.
        path: string (None). Path relative to run directory to print output data to. Defaults to './'.
        outFname: strgin (None). Output filename to use when writeout = False.
    """
def InducedAeList(r_bds, sigmas, omegas, rscale_moments, nn=1, writeout=False, path=None, append="", do_parallel=True,
                  outFname=None):
    if writeout:
        log.debug(f"Calculating A_e for {np.size(omegas)} omega values.")

    n_omegas = np.size(omegas)
    Aes = np.zeros(n_omegas, dtype=np.complex128)

    if rscale_moments == 1.0:
        rscaling = rscale_moments
    else:
        rscaling = rscale_moments * r_bds[-1]

    # For each omega, evaluate Ae:
    if do_parallel:
        pool = mtpContext.Pool(num_cores)
        par_result = [pool.apply_async(AeResponse, args=(r_bds, sigmas, omegas[i_om], rscaling), kwds={'nn':nn}) for i_om in range(n_omegas)]
        pool.close()
        pool.join()
        Aes = [par_result[i_om].get()[0] for i_om in range(n_omegas)]
    else:
        for i_om in range(n_omegas):
            Aes[i_om] = AeResponse(r_bds,sigmas,omegas[i_om],rscaling,nn=nn)[0]
            if (i_om*4) % n_omegas < 4: log.debug(f"{i_om + 1} of {n_omegas} complete.")

    Aes = np.array([ complex(val) for val in Aes ])
    AeM = np.abs(Aes)
    AeA = np.angle(Aes)

    if writeout:
        T_hrs = [ 2*np.pi/omega_i/3600 for omega_i in omegas ]

        if path is None:
            path = os.getcwd()+'/'
        if outFname is None:
            outFname = f'complexAes{append}'
        fpath = os.path.join(path, f"{outFname}.dat")
        fout = open(fpath, "w")
        header = "{:<13}, {:<24}, {:<24}\n".format("Period (hr)", "Ae.mag", "Ae.arg (rad)")
        fout.write(header)
        [ fout.write( f"{T_hrs[i]:13.5f}, {AeM[i]:24.12e}, {AeA[i]:24.12e}\n" ) for i in range(np.size(omegas)) ]
        fout.close()
        log.info(f"Data for Aes written to file: {fpath}")

    return Aes, AeM, AeA

#############################################

"""
BiList()
    The complex induced magnetic moments for a given excitation field,
    with spherical symmetry and for a list of omega values.
    Usage: `Binm` = BiList(`r_bds`, `sigmas`, `peak_omegas`, `Benm`, `nprmvals`, `mprmvals`, `rscale_moments`, `n_max=1`, `writeout=True`,
                            `path=None`, `bodyname=None`, `append=""`)
    Returns:
        Binm: complex, shape(n_peaks,2,n_max+1,n_max+1). List of complex amplitude responses for the given Benm, frequencies, and interior structure.
    Parameters:
        r_bds: float, shape(N). Radii for each boundary surface in m.
        sigmas: float, shape(N). Conductivity for each layer under the corresponding boundary in S/m.
        peak_omegas: float, shape(n_omegas). Angular frequencies of peak oscillations in rad/s.
        Benm: complex, shape(n_omegas,2,n_max+1,n_max+1). Complex excitation amplitudes for each peak oscillation.
        nprmvals: integer, shape(Nnmprm). A linear list of n' values for constructing (n',m') pairs in parallelized loops.
        mprmvals: integer, shape(Nnmprm). A linear list of m' values for constructing (n',m') pairs in parallelized loops.
        rscale_moments: float. Equal to 1/R_body in units of 1/m. For proper scaling when conducting boundaries
            may or may not coincide with the body surface.
        n_max: integer (1). Largest n' represented in the excitation moments.
        writeout: boolean (True). Whether to save computed values to disk for rapid replotting.
        path: string (None). Path relative to run directory to print output data to. Defaults to '<install_loc>/MoonMag/induced/'.
        bodyname: string (None). Body name to include in writeout filename.
        append: string (""). Optional string appended to default file names.
        outFname: string (None). Output filename to use when writeout = True.
        outFnameS: string (None). As above, for output Gauss coefficients in the Schmidt normalization.
    """
def BiList(r_bds, sigmas, peak_omegas, Benm, nprmvals, mprmvals, rscale_moments, n_max=1, writeout=True, path=None, 
           bodyname=None, append="", Schmidt=False, outFname=None, outFnameS=None):
    if writeout:
        log.debug(f"Calculating symmetric B_inm for {np.size(peak_omegas)} periods.")

        if bodyname is None:
            bfname = ""
        else:
            bfname = bodyname+"_"

    n_omegas = np.size(peak_omegas)
    Nnm = np.size(nprmvals)
    Binms = np.zeros((n_omegas,2,n_max+1,n_max+1), dtype=np.complex128)

    # For each degree n, evaluate all Ae:
    for ni in range(1,n_max+1):
        Aes, _, _ = InducedAeList(r_bds, sigmas, peak_omegas, rscale_moments, nn=ni)
        n = mp.mpf(ni)
        for i_om in range(n_omegas):
            Binms[i_om,:,ni,:] = n/(n+one)*Benm[i_om,:n_max+1,ni,:n_max+1]*Aes[i_om]

    if writeout:
        if path is None:
            path = _induced
        if outFname is None:
            outFname = f'{bfname}Binm_sym{append}'
        fpath = os.path.join(path, f"{outFname}.dat")
        fout = open(fpath, "w")
        header = "{:<13}, {:<4}, {:<4}, {:<24}, {:<24}\n".format("Period (hr) ", "n ", " m ", "Binm_Re (nT)", "Binm_Im (nT)")
        fout.write(header)
        for i in range(np.size(peak_omegas)):
            T_hrs = 2*np.pi/peak_omegas[i]/3600
            for i_nm in range(Nnm):
                sign = int(mprmvals[i_nm]<0)
                this_Binm = Binms[i,sign,nprmvals[i_nm],abs(mprmvals[i_nm])]
                fout.write( "{:<13}, {:<4}, {:<4}, {:<24}, {:<24}\n".format(round(T_hrs,5), nprmvals[i_nm], mprmvals[i_nm], np.real(this_Binm), np.imag(this_Binm)) )
        fout.close()
        log.info(f"Data for symmetric Binm written to file: {fpath}")

        if Schmidt:
            if outFnameS is None:
                outFnameS = f'{bfname}ghnm_sym{append}'
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
            log.info(f"Data for symmetric, Schmidt semi-normalized g_nm and h_nm written to file: {fpath}")

    return Binms

#############################################

"""
get_gh_from_Binm()
    Convert from orthonormal harmonic coefficients with the Condon-Shortley phase (common in physics)
    to Schmidt semi-normalized form without the C-S phase (common in geophysics).
    Handles all values for a given n at once.
    Usage: `gnm, hnm` = get_gh_from_Binm(`n`, `n_max`, `Binm`)
    Returns:
        gnm, hnm: complex, shape(n_max+1,n_max+1). g_nm and h_nm values for all m = [0,n].
            Schmidt normalization here means the integral of |Ynm|^2 * dOmega over a unit sphere is
            4pi/(2n+1) for all n and m. No Condon-Shortley phase.
    Parameters:
        n_max: integer. Maximum degree n of induced moments.
        Binm: complex, shape(2,n_max+1,n_max+1). Complex induced magnetic moments calculated using fully normalized spherical harmonic coefficients.
    """
def get_gh_from_Binm(n_max, Binm):
    gnm, hnm = ( np.zeros((n_max+1,n_max+1),dtype=np.complex128) for _ in range(2) )

    for n in range(1,n_max+1):
        norm = np.sqrt(2*n+1) / sqrt2 / sqrt4pi

        for m in range(1,n+1):
            # g terms:
            gnm[n,m] = ((-1)**m * Binm[0,n,m] + Binm[1,n,m]) * norm
            # h terms:
            hnm[n,m] = ((-1)**m * Binm[0,n,m] - Binm[1,n,m]) * 1j * norm

        gnm[n,0] = Binm[0,n,0] * norm * sqrt2

    return gnm, hnm

#############################################
