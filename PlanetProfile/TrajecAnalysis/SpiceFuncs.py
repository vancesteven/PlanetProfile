""" Functions for interfacing with SPICE kernels.
"""

import logging
import os.path
import numpy as np
import spiceypy as spice
from glob import glob as GetFilesFromPattern
from MoonMag.field_xyz import eval_Bi
from PlanetProfile.Utilities.defineStructs import Constants

# Parallel processing
import multiprocessing as mtp
mtpFork = mtp.get_context('fork')
# Assign logger
log = logging.getLogger('PlanetProfile')

def LoadKernels(Params, parent, scName):
    """ Load all SPICE kernels relevant to the task we intend.
        We load a TLS kernel in GetConfig so that we can use spiceypy.str2et
        in config settings, so we skip that one here.
    """

    if parent is not None and parent != 'None':
        log.debug(f'Loading all SPICE kernels for {scName}:')
        kernelList = [
            os.path.join(Params.spiceDir, Params.spicePCK),
            os.path.join(Params.spiceDir, Params.spiceBSP[parent])
        ] + [os.path.join(Params.spiceDir, FK) for FK in Params.spiceFK] \
          + GetFilesFromPattern(os.path.join(Params.spiceSC[scName], '*.bsp'))
        log.debug(', '.join(kernelList))
        for kernel in kernelList:
            if not os.path.isfile(kernel):
                log.error(f'Kernel file not found: {kernel}')
        spice.furnsh(kernelList)

    return


def BodyDist_km(spiceSCname, bodyname, ets, coord=None):
    """ Return distance from spacecraft to target body in km for each ephemeris time
        in ets.
    """

    spiceBody = bodyname.upper()
    if coord is None:
        coord = f'IAU_{spiceBody}'

    pos, _ = spice.spkpos(spiceSCname, ets, coord, 'NONE', spiceBody)
    x_km = pos[:, 0]
    y_km = pos[:, 1]
    z_km = pos[:, 2]
    r_km = np.sqrt(x_km ** 2 + y_km ** 2 + z_km ** 2)

    return x_km, y_km, z_km, r_km


def BiTrajec(Planet, Params, spiceSCname, ets):
    nExc = np.size(Planet.Magnetic.omegaExc_radps)
    nPts = np.size(ets)
    Bix_nT, Biy_nT, Biz_nT = (np.zeros((nExc, nPts), dtype=np.complex_) for _ in range(3))
    if Planet.Magnetic.BinmLin_nT is not None:
        Nnm = np.size(Planet.Magnetic.nLin[Planet.Magnetic.nLin <= 4])
        x_Rp, y_Rp, z_Rp, r_Rp = (xyz_km * 1e3 / Planet.Bulk.R_m for xyz_km in BodyDist_km(spiceSCname, Planet.bodyname, ets))
        nCores = np.min([Params.maxCores, Nnm, Params.threadLimit])

        for iExc, omega_radps in enumerate(Planet.Magnetic.omegaExc_radps):
            if Params.DO_PARALLEL:
                pool = mtpFork.Pool(nCores)
                par_result = [pool.apply_async(eval_Bi, args=(Planet.Magnetic.nLin[iN], Planet.Magnetic.mLin[iN],
                                                              Planet.Magnetic.BinmLin_nT[iExc,iN], x_Rp, y_Rp, z_Rp, r_Rp),
                                                        kwds={'omega': omega_radps, 't': ets})
                              for iN in range(Nnm)]
                pool.close()
                pool.join()
                # Unpack results from parallel processing and sum them
                for res in par_result:
                    this_Bx, this_By, this_Bz = res.get()
                    Bix_nT[iExc,:] = Bix_nT[iExc,:] + this_Bx
                    Biy_nT[iExc,:] = Biy_nT[iExc,:] + this_By
                    Biz_nT[iExc,:] = Biz_nT[iExc,:] + this_Bz
            else:
                for iN in range(Nnm):
                    this_Bx, this_By, this_Bz = eval_Bi(Planet.Magnetic.nLin[iN], Planet.Magnetic.mLin[iN],
                            Planet.Magnetic.BinmLin_nT[iExc,iN], x_Rp, y_Rp, z_Rp, r_Rp,
                            omega=omega_radps, t=ets)

                    Bix_nT[iExc,:] = Bix_nT[iExc,:] + this_Bx
                    Biy_nT[iExc,:] = Biy_nT[iExc,:] + this_By
                    Biz_nT[iExc,:] = Biz_nT[iExc,:] + this_Bz

    else:
        log.warning('Magnetic.BinmLin_nT has not been assigned. BiTrajec cannot be evaluated.')

    return np.real(Bix_nT), np.real(Biy_nT), np.real(Biz_nT)


def spiceCode(name):
    if name == 'Galileo':
        code, parent = (-77, 599)
    elif name == 'Cassini':
        code, parent = (-82, 699)
    elif name == 'Juno':
        code, parent = (-61, 599)

    elif name == 'Jupiter':
        code, parent = (599, 0)
    elif name == 'Saturn':
        code, parent = (699, 0)
    elif name == 'Uranus':
        code, parent = (799, 0)
    elif name == 'Neptune':
        code, parent = (899, 0)

    elif name == 'Io':
        code, parent = (501, 599)
    elif name == 'Europa':
        code, parent = (502, 599)
    elif name == 'Ganymede':
        code, parent = (503, 599)
    elif name == 'Callisto':
        code, parent = (504, 599)

    elif name == 'Mimas':
        code, parent = (601, 699)
    elif name == 'Enceladus':
        code, parent = (602, 699)
    elif name == 'Tethys':
        code, parent = (603, 699)
    elif name == 'Dione':
        code, parent = (604, 699)
    elif name == 'Rhea':
        code, parent = (605, 699)
    elif name == 'Titan':
        code, parent = (606, 699)

    elif name == 'Miranda':
        code, parent = (705, 799)  # Note the integer code out of radial sequence.
    elif name == 'Ariel':
        code, parent = (701, 799)
    elif name == 'Umbriel':
        code, parent = (702, 799)
    elif name == 'Titania':
        code, parent = (703, 799)
    elif name == 'Oberon':
        code, parent = (704, 799)

    elif name == 'Triton':
        code, parent = (801, 899)

    else:
        log.warning(f'Body name {name} did not match a defined spacrcraft, planet, or moon.')
        code, parent = (None, None)

    if parent == 0:
        parentName = 'Sun'
    elif parent == 599:
        parentName = 'Jupiter'
    elif parent == 699:
        parentName = 'Saturn'
    elif parent == 799:
        parentName = 'Uranus'
    elif parent == 899:
        parentName = 'Neptune'
    else:
        parentName = None

    return code, parent, parentName


parentGM = {  # Masses retrieved from https://ssd.jpl.nasa.gov/planets/phys_par.html . Units are in km^3/s^2, as required by SPICE routines.
    'Jupiter': 1898.1250e24 * Constants.G * 1e-9,
    'Saturn':   568.3170e24 * Constants.G * 1e-9,
    'Uranus':    86.8099e24 * Constants.G * 1e-9,
    'Neptune':  102.4092e24 * Constants.G * 1e-9
}


spiceSCname = {
    'Galileo': 'GALILEO ORBITER',
    'Cassini': 'CASSINI',
    'Juno': 'JUNO'
}
