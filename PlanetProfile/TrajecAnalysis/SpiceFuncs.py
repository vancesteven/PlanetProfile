""" Functions for interfacing with SPICE kernels.
"""

import logging
import os.path
import numpy as np
import spiceypy as spice
from collections.abc import Iterable
from PlanetProfile.Utilities.defineStructs import Constants, EOSlist
from PlanetProfile import _SPICE

# Parallel processing
import multiprocessing as mtp
import platform
plat = platform.system()
if plat == 'Windows':
    mtpType = 'spawn'
else:
    mtpType = 'fork'
mtpContext = mtp.get_context(mtpType)
# Assign logger
log = logging.getLogger('PlanetProfile')

def LoadKernels(Params, parent, scName):
    """ Load all SPICE kernels relevant to the task we intend.
        We load a TLS kernel in GetConfig so that we can use spiceypy.str2et
        in config settings, so we skip that one here.
    """

    # Only take action if we haven't yet loaded kernels
    spiceKey = f'SPICE{parent}{scName}'
    ALREADY_LOADED = spiceKey in EOSlist.loaded.keys()
    if not ALREADY_LOADED:
        # Load kernels
        if parent is not None and parent != 'None':
            log.debug(f'Loading all SPICE kernels for {scName}:')
            if np.size(parent) == 1:
                parentStr = parent
                extraParents = [None]
            else:
                parentStr = parent[0]
                extraParents = parent[1:]
            kernelList = [
                os.path.join(Params.spiceDir, Params.spicePCK),
                os.path.join(Params.spiceDir, Params.spiceBSP[parentStr])
            ] + [os.path.join(Params.spiceDir, FK) for FK in Params.spiceFK] \
              + [os.path.join(Params.spiceDir, scName, spkFile)
                 for spkFile in Params.Trajec.spiceSPK[scName]]
            kernelList = kernelList + [os.path.join(Params.spiceDir, Params.spiceBSP[xparent])
                                       for xparent in extraParents if xparent is not None]
            log.debug(', '.join(kernelList))
            ERRORS = False
            for kernel in kernelList:
                if not os.path.isfile(kernel):
                    log.error(f'Kernel file not found: {kernel}')
                    ERRORS = True

            if ERRORS:
                log.error(f'At least one kernel file was not found. Review the SPICE README file, copied below, for ' +
                          f'instructions on placing kernels where spiceypy can find them.')
                spiceREADME = open(os.path.join(_SPICE, 'README.md'), 'r').read()
                log.error(spiceREADME)
                raise FileNotFoundError('See above for missing SPICE kernel files.')

            spice.furnsh(kernelList)
            EOSlist.loaded[spiceKey] = kernelList

    return


def BodyDist_km(spiceSC, bodyname, ets, coord=None):
    """ Return distance from spacecraft to target body in km for each ephemeris time
        in ets.
    """

    spiceBody = bodyname.upper()
    if coord is None:
        coord = f'IAU_{spiceBody}'

    if not isinstance(ets, Iterable):
        etin = np.array([ets])
    else:
        etin = ets
    pos, _ = spice.spkpos(spiceSC, etin, coord, 'NONE', spiceBody)
    x_km = pos[:, 0]
    y_km = pos[:, 1]
    z_km = pos[:, 2]
    r_km = np.sqrt(x_km**2 + y_km**2 + z_km**2)

    return x_km, y_km, z_km, r_km


def BodyVel_kms(spiceSC, bodyname, ets, coord=None):
    """ Return relative velocity between spacecraft and target body in km/s for each ephemeris time
        in ets.
    """

    spiceBody = bodyname.upper()
    if coord is None:
        coord = f'IAU_{spiceBody}'

    if not isinstance(ets, Iterable):
        etin = np.array([ets])
    else:
        etin = ets

    state, _ = spice.spkezr(spiceSC, etin, coord, 'NONE', spiceBody)
    # The SPICE spkezr function returns a list of arrays instead of a 2D array. Convert to array
    state = np.array(state)
    vx_kms = state[:, 3]
    vy_kms = state[:, 4]
    vz_kms = state[:, 5]
    v_kms = np.sqrt(vx_kms**2 + vy_kms**2 + vz_kms**2)

    return vx_kms, vy_kms, vz_kms, v_kms


def RotateFrame(vec, ets, fromCoord, toCoord):
    """
    Rotate a vector from one frame to another for a list of ephemeris times.

    Parameters
    ----------
    vec : float, array_like, shape Nx3
        Vectors aligned to fromCoord bases corresponding to ephemeris times ets.
    ets : float, array_like, shape N or 1
        Ephemeris times at which to evaluate frame transformations for each vec. If only one value
        is passed, each vec will be transformed at this same ephemeris time.
    fromCoord : str
        Name of a frame implemented in SPICE via loaded kernels to which vec are aligned.
    toCoord : str
        Name of a frame available in loaded kernels to which to rotate vec.

    Returns
    -------
    outVec: float, array_like, shape Nx3
        Transformed vectors corresponding to vec in the toCoord frame.
    """

    if np.size(ets) == 1:
        rotMat = spice.pxform(fromCoord, toCoord, ets)
        outVec = np.array([spice.mxv(rotMat, veci) for veci in vec])

    else:
        outVec = np.array([
            spice.mxv(spice.pxform(fromCoord, toCoord, et), vec[i_et,:])
            for i_et, et in enumerate(ets)
        ])

    return outVec


def spiceCode(name):
    if name == 'Pioneer 11':
        code, parent = (-24, 0)
    elif name == 'Voyager 1':
        code, parent = (-31, 0)
    elif name == 'Voyager 2':
        code, parent = (-32, 0)
    elif name == 'Galileo':
        code, parent = (-77, 599)
    elif name == 'Cassini':
        code, parent = (-82, 699)
    elif name == 'Juno':
        code, parent = (-61, 599)
    elif name == 'JUICE':
        code, parent = (-28, 599)
    elif name == 'Clipper':
        code, parent = (-159, 599)

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
        log.warning(f'Body name {name} did not match a defined spacecraft, planet, or moon.')
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

# Masses below retrieved from https://ssd.jpl.nasa.gov/planets/phys_par.html.
# Units are in km^3/s^2, as required by SPICE routines.
parentGM = {
    'Jupiter': 1898.1250e24 * Constants.G * 1e-9,
    'Saturn':   568.3170e24 * Constants.G * 1e-9,
    'Uranus':    86.8099e24 * Constants.G * 1e-9,
    'Neptune':  102.4092e24 * Constants.G * 1e-9
}

# Code names recognized by SPICE for relevant spacecraft
spiceSCname = {
    'Cassini': 'CASSINI',
    'Clipper': 'EUROPA CLIPPER',
    'Galileo': 'GALILEO ORBITER',
    'Juno': 'JUNO',
    'JUICE': 'JUICE',
    'Voyager 1': 'VOYAGER 1',
    'Voyager 2': 'VOYAGER 2'
}

# System III coordinates in which magnetic field data are stored
spiceS3coords = {
    'Jupiter': 'IAU_JUPITER',
    'Saturn': 'IAU_SATURN',
    'Uranus': 'ULS',
    'Neptune': 'NLS',
}
