"""
Loads a selected spacecraft dataset from ASCII or CDF files and saves as a standardized HDF5 .mat
file format.
"""

import logging
import os
import numpy as np
import spiceypy as spice
from hdf5storage import savemat, loadmat
from PlanetProfile.Utilities.defineStructs import MAGdataStruct, Constants, ParentName
from PlanetProfile.TrajecAnalysis import _MAGdir, _scList, _MAGdataList
from PlanetProfile.TrajecAnalysis.SpiceFuncs import spiceCode, spiceSCname, spiceS3coords, \
    BodyDist_km, RotateFrame, LoadKernels
from PlanetProfile.TrajecAnalysis.FlybyEvents import GetFlybyCA
from PlanetProfile.TrajecAnalysis.MagneticFields import Bsph2Bxyz

# Assign logger
log = logging.getLogger('PlanetProfile')


def RefileName(targetBody, scName, MAGdir=None):
    """
    Returns a standardized file name for reformatted MAG data for the given target body and
    spacecraft.

    Parameters
    ----------
    targetBody : str
        Target body for which to return the file path.
    scName : str
        Spacecraft name for which to return file path.
    MAGdir : str, default='SpacecraftMAGdata'
        Directory where spacecraft data are stored.

    Returns
    -------
    fPath : str
        Full file path for the named spacecraft's reformatted MAG data.
    """

    if MAGdir is None:
        MAGdir = _MAGdir

    fName = f'{scName}MAG{targetBody}.mat'
    fPath = os.path.join(MAGdir, scName, fName)

    return fPath


def MAGtoHDF5(Params, scName, MAGdir=None):
    """
    Converts ASCII text files of spacecraft magnetometer measurements for the given target body to
    the IAU frame and saves as HDF5.

    Parameters
    ----------
    Params : ParamsStruct
        A ParamsStruct class object containing trajectory analysis settings and parameters in
        the TrajecParamsStruct subclass object Params.Trajec.
    targetBody : str
        Target body for which to convert measurements.
    scName : str
        Spacecraft name for which to convert measurements.
    MAGdir : str, default='SpacecraftMAGdata'
        Directory in which to search for spacecraft data.
    """

    if MAGdir is None:
        MAGdir = _MAGdir

    log.debug(f'Reformatting data for {scName} encounters of {Params.Trajec.targetBody}.')
    FlybyCA = GetFlybyCA()
    if Params.Trajec.targetBody in _MAGdataList[scName].keys():
        # We get here if we have a definition for file name format for relevant flyby B data
        pdsFiles = _MAGdataList[scName][Params.Trajec.targetBody]

        t_UTC, ets, BxS3_nT, ByS3_nT, BzS3_nT, BrS3_nT, BthS3_nT, BphiS3_nT = ({} for _ in range(8))
        if scName == 'Galileo':
            if Params.Trajec.EXPANDED_RANGE or Params.Trajec.targetBody == 'Jupiter':
                if Params.Trajec.EXPANDED_RANGE:
                    # Replace files to load with those for full orbit
                    pdsFiles = {fbID: _MAGdataList[scName]['Jupiter'][fbID]
                                for fbID in pdsFiles.keys()}
                for fbID, file in pdsFiles.items():
                    fpath = os.path.join(MAGdir, scName, 'Jupiter', file)
                    t_UTC[fbID], _, BrS3_nT[fbID], BthS3_nT[fbID], BphiS3_nT[fbID], _, _, _, _, _ \
                        = np.loadtxt(fpath, unpack=True, dtype='U23,U23,f,f,f,f,f,f,f,f')
                    log.debug(f'Loaded file {file}.')
            else:
                for fbID, file in pdsFiles.items():
                    fpath = os.path.join(MAGdir, scName, Params.Trajec.targetBody, file)
                    t_UTC[fbID], BrS3_nT[fbID], BthS3_nT[fbID], BphiS3_nT[fbID], _, _, _, _, _ \
                        = np.loadtxt(fpath, unpack=True, dtype='U23,f,f,f,f,f,f,f,f')
                    log.debug(f'Loaded file {file}.')

        elif scName == 'Cassini':
            for revID, file in pdsFiles.items():
                fpath = os.path.join(MAGdir, scName, Params.Trajec.targetBody, file)
                t_UTC[revID], BrS3_nT[revID], BthS3_nT[revID], BphiS3_nT[revID], _, _ \
                    = np.loadtxt(fpath, unpack=True, dtype='U21,f,f,f,f,d')
                log.debug(f'Loaded file {file}.')

        elif scName == 'Juno':
            for pjID, fileDict in pdsFiles.items():
                t_UTC[pjID] = []
                for file in fileDict.values():
                    fpath = os.path.join(MAGdir, scName, Params.Trajec.targetBody, file)
                    nHeadLines = 0
                    with open(fpath) as f:
                        headLine = f.readline()
                        while f'{headLine[2:6]}{headLine[7:10]}' != file[11:18]:
                            headLine = f.readline()
                            nHeadLines += 1

                    yyyy, doy, h, m, s, ms, _, BxS3_nT[pjID], ByS3_nT[pjID], BzS3_nT[pjID], _, _, \
                        _, _ = np.loadtxt(fpath, unpack=True, skiprows=nHeadLines,
                                     dtype='d,d,d,d,d,d,f,f,f,f,f,f,f,f')
                    t_UTC[pjID] += [f'{yyyy1:.0f}-{doy1:03.0f}//{h1:02.0f}:{m1:02.0f}:{s1:02.0f}.' +
                                    f'{ms1:03.0f}'
                                    for yyyy1, doy1, h1, m1, s1, ms1 in zip(yyyy, doy, h, m, s, ms)]
                    log.debug(f'Loaded file {file}.')

                t_UTC[pjID] = np.array(t_UTC[pjID])

        elif scName in _scList:
            raise ValueError(f'PDS data read-in has not yet been implemented for "{scName}".')
        else:
            raise ValueError(f'No data for "{scName}" found in {MAGdir}.')

        ets = {fbID: spice.str2et(t) for fbID, t in t_UTC.items()}

    else:
        # Generate trajectory data with zero B values for testing infrastructure for upcoming flybys
        pdsFiles = {fbID: 'N/A' for fbID in FlybyCA[scName].etCA.keys()}
        rangeHalf = Params.Trajec.etPredRange_s/2
        range_s = np.arange(-rangeHalf, rangeHalf, Params.Trajec.etStep_s)
        Bzero_nT = np.zeros_like(range_s)
        ets = {fbID: etCA + range_s
               for fbID, etCA in FlybyCA[scName].etCA[Params.Trajec.targetBody].items()}
        t_UTC = {fbID: spice.et2utc(fbets, 'ISOC', 3) for fbID, fbets in ets.items()}
        BxS3_nT, ByS3_nT, BzS3_nT = ({fbID: Bzero_nT for fbID in ets.keys()} for _ in range(3))
        BrS3_nT, BthS3_nT, BphiS3_nT = (BxS3_nT, ByS3_nT, BzS3_nT)

    if Params.Trajec.PLANETMAG_MODEL:
        ambModel = Params.Trajec.BextDefault[Params.Trajec.targetBody]
        BxS3amb_nT, ByS3amb_nT, BzS3amb_nT, BxIAUamb_nT, ByIAUamb_nT, BzIAUamb_nT \
            = ({} for _ in range(6))

        for fbID, pdsFile in pdsFiles.items():
            fdir = os.path.join(MAGdir, scName)
            if not Params.Trajec.EXPANDED_RANGE:
                fdir = os.path.join(fdir, Params.Trajec.targetBody)
            elif scName != 'Cassini':
                fdir = os.path.join(fdir, ParentName(Params.Trajec.targetBody))
            if isinstance(pdsFile, dict):
                # Juno pdsFiles are listed by the DOY in which they occur within each PJ
                pdsFileHere = list(pdsFile.values())[0]
            else:
                pdsFileHere = pdsFile + ''
            ambModelData = loadmat(os.path.join(fdir, f'{pdsFileHere[:-4]}{ambModel}.mat'))
            BxS3amb_nT[fbID] = np.squeeze(ambModelData['BxS3_nT'])
            ByS3amb_nT[fbID] = np.squeeze(ambModelData['ByS3_nT'])
            BzS3amb_nT[fbID] = np.squeeze(ambModelData['BzS3_nT'])

    else:
        ambModel, BxIAUamb_nT, ByIAUamb_nT, BzIAUamb_nT = (Constants.NA for _ in range(4))
        BxS3amb_nT, ByS3amb_nT, BzS3amb_nT = (None for _ in range(3))

    if Params.Trajec.EXPANDED_RANGE:
        span_s = Params.Trajec.etExpandRange_s / 2
        keep = {fbID: np.abs(fbets - FlybyCA[scName].etCA[Params.Trajec.targetBody][fbID]) <= span_s
                for fbID, fbets in ets.items()}

        t_UTC = {fbID: t_UTC[fbID][fbkeep] for fbID, fbkeep in keep.items()}
        ets =   {fbID: ets[fbID][fbkeep]   for fbID, fbkeep in keep.items()}
        if scName == 'Juno':
            BxS3_nT = {pjID: BxS3_nT[pjID][pjkeep] for pjID, pjkeep in keep.items()}
            ByS3_nT = {pjID: ByS3_nT[pjID][pjkeep] for pjID, pjkeep in keep.items()}
            BzS3_nT = {pjID: BzS3_nT[pjID][pjkeep] for pjID, pjkeep in keep.items()}
        else:
            BrS3_nT =   {fbID: BrS3_nT[fbID][fbkeep]   for fbID, fbkeep in keep.items()}
            BthS3_nT =  {fbID: BthS3_nT[fbID][fbkeep]  for fbID, fbkeep in keep.items()}
            BphiS3_nT = {fbID: BphiS3_nT[fbID][fbkeep] for fbID, fbkeep in keep.items()}

        if Params.Trajec.PLANETMAG_MODEL:
            BxS3amb_nT = {fbID: BxS3amb_nT[fbID][fbkeep] for fbID, fbkeep in keep.items()}
            ByS3amb_nT = {fbID: ByS3amb_nT[fbID][fbkeep] for fbID, fbkeep in keep.items()}
            BzS3amb_nT = {fbID: BzS3amb_nT[fbID][fbkeep] for fbID, fbkeep in keep.items()}

    # Adjust data to common format
    scCode, pCode, parent = spiceCode(scName)
    LoadKernels(Params, parent, scName)
    if scName in ['Galileo', 'Cassini']:
        # Convert from spherical coordinates to Cartesian for frame transformation
        for fbID, fbets in ets.items():
            x_km, y_km, z_km, r_km = BodyDist_km(spiceSCname[scName], parent,
                                                 fbets, coord=spiceS3coords[parent])
            thS3_rad = np.arccos(z_km/r_km)
            phiS3_rad = np.arctan2(y_km, x_km)
            BxS3_nT[fbID], ByS3_nT[fbID], BzS3_nT[fbID] \
                = Bsph2Bxyz(BrS3_nT[fbID], BthS3_nT[fbID], BphiS3_nT[fbID], thS3_rad, phiS3_rad)

    elif scName == 'Juno' and Params.Trajec.targetBody != 'Jupiter' \
         and not Params.Trajec.EXPANDED_RANGE:
        # Pare down data to flyby encounter itself
        Rp_km = spice.bodvcd(pCode, 'RADII', 3)[1][0]
        for pjID, pjets in ets.items():
            _, _, _, r_km = BodyDist_km(spiceSCname[scName], Params.Trajec.targetBody, pjets)
            keep = r_km <= (Rp_km * Params.Trajec.fbRange_Rp)
            t_UTC[pjID] = t_UTC[pjID][keep]
            ets[pjID] = ets[pjID][keep]
            BxS3_nT[pjID] = BxS3_nT[pjID][keep]
            ByS3_nT[pjID] = ByS3_nT[pjID][keep]
            BzS3_nT[pjID] = BzS3_nT[pjID][keep]

    # Convert measurements to moon body-fixed IAU frame
    BxIAU_nT, ByIAU_nT, BzIAU_nT = ({} for _ in range(3))
    IAUtarget = f'IAU_{Params.Trajec.targetBody.upper()}'
    for fbID, fbets in ets.items():
        BS3_nT = np.ascontiguousarray(np.vstack((BxS3_nT[fbID], ByS3_nT[fbID], BzS3_nT[fbID])).T)
        BIAU_nT = RotateFrame(BS3_nT, fbets, spiceS3coords[parent], IAUtarget)
        BxIAU_nT[fbID] = BIAU_nT[:, 0]
        ByIAU_nT[fbID] = BIAU_nT[:, 1]
        BzIAU_nT[fbID] = BIAU_nT[:, 2]

        if Params.Trajec.PLANETMAG_MODEL:
            BS3amb_nT = np.ascontiguousarray(np.vstack((BxS3amb_nT[fbID], ByS3amb_nT[fbID],
                                                        BzS3amb_nT[fbID])).T)
            BIAUamb_nT = RotateFrame(BS3amb_nT, fbets, spiceS3coords[parent], IAUtarget)
            BxIAUamb_nT[fbID] = BIAUamb_nT[:, 0]
            ByIAUamb_nT[fbID] = BIAUamb_nT[:, 1]
            BzIAUamb_nT[fbID] = BIAUamb_nT[:, 2]

    data = {
        'fbInclude': {scName: Params.Trajec.fbInclude[scName]},
        'pdsFiles': pdsFiles,
        't_UTC': t_UTC,
        'ets': ets,
        'BxS3_nT': BxS3_nT,
        'ByS3_nT': ByS3_nT,
        'BzS3_nT': BzS3_nT,
        'BxIAU_nT': BxIAU_nT,
        'ByIAU_nT': ByIAU_nT,
        'BzIAU_nT': BzIAU_nT,
        'ambModel': ambModel,
        'BxIAUamb_nT': BxIAUamb_nT,
        'ByIAUamb_nT': ByIAUamb_nT,
        'BzIAUamb_nT': BzIAUamb_nT
    }

    magFname = RefileName(Params.Trajec.targetBody, scName)
    if os.path.exists(magFname):
        os.remove(magFname)
    savemat(magFname, data)
    log.debug(f'Reformatted MAG data saved to file: {magFname}')

    return


def LoadMAG(Params, magFname, scName):
    """
    Load reformatted MAG data from disk into a MAGdata class object.

    Parameters
    ----------
    Params : ParamsStruct
        A ParamsStruct class object containing trajectory analysis settings and parameters in
        the TrajecParamsStruct subclass object Params.Trajec.
    magFname : str
        File name for HDF5-formatted data to load into MAGdata class object.

    Returns
    -------
    magData : MAGdata
        A MAGdata class object containing timestamps, locations, and magnetic field data in the
        parent planet's System III frame.
    """

    loadDict = loadmat(magFname)
    magData = MAGdataStruct(Params, magFname, loadDict)
    _, _, parent = spiceCode(Params.Trajec.targetBody)
    LoadKernels(Params, parent, scName)

    # Add spacecraft positions from SPICE kernels
    for fbID, ets in magData.ets.items():
        log.debug(f'Evaluating {np.size(ets)} {magData.scName} positions relative to ' +
                  f'{Params.Trajec.targetBody} for flyby {fbID}.')
        magData.x_km[fbID], magData.y_km[fbID], magData.z_km[fbID], magData.r_km[fbID] \
            = BodyDist_km(spiceSCname[scName], Params.Trajec.targetBody, ets)

    return magData
