import logging
import numpy as np
from PlanetProfile.TrajecAnalysis.SpiceFuncs import BodyDist_km
from PlanetProfile.Utilities.defineStructs import ModelDataStruct
from PlanetProfile.MagneticInduction.MagneticInduction import GetBexc
from MoonMag.field_xyz import eval_Bi as EvalBi

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


def CalcModel(Planet, Params, magData, modelData):

    # Calculate stuff
    modelData.fitProfileFname = Planet.saveFile
    modelData = CalcModelAmbient(Planet, Params, magData, modelData)
    modelData = CalcModelInduced(Planet, Params, modelData)
    modelData = CalcModelPlasma(Planet, Params, modelData)

    # Sum net fields
    modelData.BxAll_nT, modelData.ByAll_nT, modelData.BzAll_nT = (np.empty(0) for _ in range(3))
    for scName, ets in modelData.ets.items():
        for fbID in ets.keys():
            modelData.BxIAU_nT[scName][fbID] = modelData.BxIAUexc_nT[scName][fbID] + \
                                               modelData.BxIAUind_nT[scName][fbID] + \
                                               modelData.BxIAUpls_nT[scName][fbID]
            modelData.ByIAU_nT[scName][fbID] = modelData.ByIAUexc_nT[scName][fbID] + \
                                               modelData.ByIAUind_nT[scName][fbID] + \
                                               modelData.ByIAUpls_nT[scName][fbID]
            modelData.BzIAU_nT[scName][fbID] = modelData.BzIAUexc_nT[scName][fbID] + \
                                               modelData.BzIAUind_nT[scName][fbID] + \
                                               modelData.BzIAUpls_nT[scName][fbID]
            modelData.BxAll_nT = np.concatenate((modelData.BxAll_nT, modelData.BxIAU_nT[scName][fbID]))
            modelData.ByAll_nT = np.concatenate((modelData.ByAll_nT, modelData.ByIAU_nT[scName][fbID]))
            modelData.BzAll_nT = np.concatenate((modelData.BzAll_nT, modelData.BzIAU_nT[scName][fbID]))

    return modelData


def CalcModelAmbient(Planet, Params, magData, modelData):

    if Params.Trajec.PLANETMAG_MODEL:
        for scName, etsList in modelData.ets.items():
            for fbID, ets in etsList.items():
                modelData.BxIAUexc_nT[scName][fbID] = magData[scName].BxIAUamb_nT[fbID]
                modelData.ByIAUexc_nT[scName][fbID] = magData[scName].ByIAUamb_nT[fbID]
                modelData.BzIAUexc_nT[scName][fbID] = magData[scName].BzIAUamb_nT[fbID]
    else:
        modelData = CalcAmbientFromExcitation(Planet.Magnetic, modelData)

    return modelData


def CalcModelInduced(Planet, Params, modelData):

    Nnm = np.size(Planet.Magnetic.nLin)

    if isinstance(Planet.Magnetic.Benm_nT, dict):
        for scName, Binm_nT in Planet.Magnetic.Binm_nT.items():
            for fbID, ets in modelData.ets[scName].items():
                Bix_nT, Biy_nT, Biz_nT = (np.zeros(np.size(ets), dtype=np.complex_)
                                          for _ in range(3))

                # Sum over all excitations and moments
                for iExc in range(Planet.Magnetic.nExc[scName]):
                    for iN in range(Nnm):
                        this_Bx, this_By, this_Bz = EvalBi(Planet.Magnetic.nLin[iN],
                            Planet.Magnetic.mLin[iN], Planet.Magnetic.BinmLin_nT[scName][iExc,iN],
                            modelData.x_Rp[scName][fbID], modelData.y_Rp[scName][fbID],
                            modelData.z_Rp[scName][fbID], modelData.r_Rp[scName][fbID],
                            omega=Planet.Magnetic.omegaExc_radps[scName][iExc], t=ets)

                        Bix_nT = Bix_nT + this_Bx
                        Biy_nT = Biy_nT + this_By
                        Biz_nT = Biz_nT + this_Bz

                modelData.BxIAUind_nT[scName][fbID] = np.real(Bix_nT)
                modelData.ByIAUind_nT[scName][fbID] = np.real(Biy_nT)
                modelData.BzIAUind_nT[scName][fbID] = np.real(Biz_nT)

    else:
        if Params.Trajec.SCera is None:
            # We only get here if there is only one scName in scSelect.
            scName = Params.Trajec.scSelect[0]
        else:
            scName = Params.Trajec.SCera
        for fbID, ets in modelData.ets[scName].items():
            Bix_nT, Biy_nT, Biz_nT = (np.zeros(np.size(ets), dtype=np.complex_)
                                      for _ in range(3))

            # Sum over all excitations and moments
            for iExc in range(Planet.Magnetic.nExc):
                for iN in range(Nnm):
                    this_Bx, this_By, this_Bz = EvalBi(Planet.Magnetic.nLin[iN],
                        Planet.Magnetic.mLin[iN], Planet.Magnetic.BinmLin_nT[iExc,iN],
                        modelData.x_Rp[scName][fbID], modelData.y_Rp[scName][fbID],
                        modelData.z_Rp[scName][fbID], modelData.r_Rp[scName][fbID],
                        omega=Planet.Magnetic.omegaExc_radps[iExc], t=ets)

                    Bix_nT = Bix_nT + this_Bx
                    Biy_nT = Biy_nT + this_By
                    Biz_nT = Biz_nT + this_Bz

            modelData.BxIAUind_nT[scName][fbID] = np.real(Bix_nT)
            modelData.ByIAUind_nT[scName][fbID] = np.real(Biy_nT)
            modelData.BzIAUind_nT[scName][fbID] = np.real(Biz_nT)

    return modelData


def CalcModelPlasma(Planet, Params, modelData):

    if Params.Trajec.plasmaType == 'Alfven':
        raise ValueError(f'Params.Trajec.plasmaType "Alfven" not implemented yet.')

    else:
        if Params.Trajec.plasmaType != 'none':
            log.warning(f'Params.Trajec.plasmaType "{Params.Trajec.plasmaType}" not recognized. ' +
                        f'Defaulting to "none".')

        for scName in modelData.ets.keys():
            for fbID in modelData.ets[scName].keys():
                modelData.BxIAUpls_nT[scName][fbID] = np.zeros_like(modelData.BxIAUexc_nT[scName][fbID])
                modelData.ByIAUpls_nT[scName][fbID] = np.zeros_like(modelData.ByIAUexc_nT[scName][fbID])
                modelData.BzIAUpls_nT[scName][fbID] = np.zeros_like(modelData.BzIAUexc_nT[scName][fbID])

    return modelData


def InitModelData(Planet, Magnetic, Params, magData):

    R_km = Planet.Bulk.R_m / 1e3

    # Populate fields for calculations
    modelData = ModelDataStruct()
    for scName, data in magData.items():
        modelData.allFlybys = {scName: {fbID: f'{Params.Trajec.fbDescrip[scName]}{fbID}'
                                        for fbID in data.ets.keys()}
                               for scName, data in magData.items()}
        modelData.fbInclude[scName] = Params.Trajec.fbInclude[scName]
        modelData.t_UTC[scName] = data.t_UTC
        modelData.ets[scName] = data.ets

        modelData.BxIAUexc_nT[scName], modelData.ByIAUexc_nT[scName], \
            modelData.BzIAUexc_nT[scName], modelData.BxIAUind_nT[scName], \
            modelData.ByIAUind_nT[scName], modelData.BzIAUind_nT[scName], \
            modelData.BxIAUpls_nT[scName], modelData.ByIAUpls_nT[scName], \
            modelData.BzIAUpls_nT[scName], modelData.BxIAU_nT[scName], \
            modelData.ByIAU_nT[scName], modelData.BzIAU_nT[scName] = ({} for _ in range(12))

        # Get trajectory in terms of planetary radii
        modelData.x_Rp[scName] = {fbID: x_km/R_km for fbID, x_km in data.x_km.items()}
        modelData.y_Rp[scName] = {fbID: y_km/R_km for fbID, y_km in data.y_km.items()}
        modelData.z_Rp[scName] = {fbID: z_km/R_km for fbID, z_km in data.z_km.items()}
        modelData.r_Rp[scName] = {fbID: r_km/R_km for fbID, r_km in data.r_km.items()}

    modelData.etsAll = np.concatenate(([data.etsAll for data in magData.values()]))

    # Carry over MagneticSubStruct settings from body input file
    Magnetic.ionosBounds_m = Planet.Magnetic.ionosBounds_m
    Magnetic.sigmaIonosPedersen_Sm = Planet.Magnetic.sigmaIonosPedersen_Sm

    return modelData


def CalcAmbientFromExcitation(Magnetic, modelData):

    if isinstance(Magnetic.Benm_nT, dict):
        for scName, etsList in modelData.ets.items():
            Bex_nT, Bey_nT, Bez_nT = Magnetic.Bexyz_nT[scName]
            # Add extra axis for multiplying into phases
            Bex_nT, Bey_nT, Bez_nT = (Bex_nT[:,None], Bey_nT[:,None], Bez_nT[:,None])
            for fbID, ets in etsList.items():
                phases = np.exp(-1j * np.outer(Magnetic.omegaExc_radps[scName], ets))
                modelData.BxIAUexc_nT[scName][fbID] = np.real(np.sum(Bex_nT * phases, axis=0)) \
                                                      + Magnetic.B0_nT[scName][0]
                modelData.ByIAUexc_nT[scName][fbID] = np.real(np.sum(Bey_nT * phases, axis=0)) \
                                                      + Magnetic.B0_nT[scName][1]
                modelData.BzIAUexc_nT[scName][fbID] = np.real(np.sum(Bez_nT * phases, axis=0)) \
                                                      + Magnetic.B0_nT[scName][2]

    else:
        for scName, etsList in modelData.ets.items():
            Bex_nT, Bey_nT, Bez_nT = Magnetic.Bexyz_nT
            Bex_nT, Bey_nT, Bez_nT = (Bex_nT[:,None], Bey_nT[:,None], Bez_nT[:,None])
            for fbID, ets in etsList.items():
                phases = np.exp(-1j * np.outer(Magnetic.omegaExc_radps, ets))
                modelData.BxIAUexc_nT[scName][fbID] = np.real(np.sum(Bex_nT * phases, axis=0)) \
                                                      + Magnetic.B0_nT[0]
                modelData.ByIAUexc_nT[scName][fbID] = np.real(np.sum(Bey_nT * phases, axis=0)) \
                                                      + Magnetic.B0_nT[1]
                modelData.BzIAUexc_nT[scName][fbID] = np.real(np.sum(Bez_nT * phases, axis=0)) \
                                                      + Magnetic.B0_nT[2]

    return modelData


def SetupMagnetic(Planet, Params):

    Magnetic = Planet.Magnetic
    if Params.CALC_ASYM:
        # If we are doing asymmetry, limit to low-degree moments to compute in a reasonable time
        Magnetic.pMax = 2
    else:
        Magnetic.pMax = 0

    if Params.Trajec.SCera is None:
        SCeraIn = list(Params.Trajec.scSelect)
    else:
        SCeraIn = Params.Trajec.SCera
    if Params.Trajec.BextModel is None:
        BextModel = Params.Trajec.BextDefault[Params.Trajec.targetBody]
    else:
        BextModel = Params.Trajec.BextModel

    if np.size(SCeraIn) > 1:
        Texc_hr, omegaExc_radps, Benm_nT, B0_nT, Bexyz_nT, nExc = ({} for _ in range(6))
        for SCera in SCeraIn:
            Texc_hr[SCera], omegaExc_radps[SCera], Benm_nT[SCera], B0_nT[SCera], Bexyz_nT[SCera] \
                = GetBexc(Params.Trajec.targetBody, SCera, BextModel,
                          Params.Induct.excSelectionCalc, nprmMax=Magnetic.nprmMax,
                          pMax=Magnetic.pMax)
            nExc[SCera] = np.size(Texc_hr[SCera])

    else:
        if isinstance(SCeraIn, list):
            SCera = SCeraIn[0]
        else:
            SCera = SCeraIn
        Texc_hr, omegaExc_radps, Benm_nT, B0_nT, Bexyz_nT = GetBexc(Params.Trajec.targetBody, SCera,
            BextModel, Params.Induct.excSelectionCalc, nprmMax=Magnetic.nprmMax, pMax=Magnetic.pMax)
        nExc = np.size(Texc_hr)

    Magnetic.Texc_hr = Texc_hr
    Magnetic.omegaExc_radps = omegaExc_radps
    Magnetic.Benm_nT = Benm_nT
    Magnetic.B0_nT = B0_nT
    Magnetic.Bexyz_nT = Bexyz_nT
    Magnetic.nExc = nExc

    return Magnetic


def BiTrajecSingle(Planet, Params, spiceSC, ets):
    if isinstance(Planet.Magnetic.omegaExc_radps, dict):
        raise RuntimeError('BiTrajecSingle is implemented only for single-era excitation ' +
                           'evaluation.')
    nExc = np.size(Planet.Magnetic.omegaExc_radps)
    nPts = np.size(ets)
    Bix_nT, Biy_nT, Biz_nT = (np.zeros((nExc, nPts), dtype=np.complex_) for _ in range(3))
    if Planet.Magnetic.BinmLin_nT is not None:
        Nnm = np.size(Planet.Magnetic.nLin)
        x_Rp, y_Rp, z_Rp, r_Rp = (xyz_km * 1e3 / Planet.Bulk.R_m
                                  for xyz_km in BodyDist_km(spiceSC, Planet.bodyname, ets))
        nCores = np.min([Params.maxCores, Nnm, Params.threadLimit])

        for iExc, omega_radps in enumerate(Planet.Magnetic.omegaExc_radps):
            if Params.DO_PARALLEL:
                pool = mtpContext.Pool(nCores)
                par_result = [pool.apply_async(
                    EvalBi, args=(Planet.Magnetic.nLin[iN], Planet.Magnetic.mLin[iN],
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
                    this_Bx, this_By, this_Bz = EvalBi(Planet.Magnetic.nLin[iN],
                        Planet.Magnetic.mLin[iN], Planet.Magnetic.BinmLin_nT[iExc,iN],
                        x_Rp, y_Rp, z_Rp, r_Rp, omega=omega_radps, t=ets)

                    Bix_nT[iExc,:] = Bix_nT[iExc,:] + this_Bx
                    Biy_nT[iExc,:] = Biy_nT[iExc,:] + this_By
                    Biz_nT[iExc,:] = Biz_nT[iExc,:] + this_Bz

    else:
        log.warning('Magnetic.BinmLin_nT has not been assigned. BiTrajecSingle cannot be evaluated.')

    return np.real(Bix_nT), np.real(Biy_nT), np.real(Biz_nT)


def Bxyz2Bsph(Bx, By, Bz, theta, phi):
    """
    Convert vector components aligned to cartesian axes into vector components aligned to spherical
    coordinates.
    Source: Arfken, Weber, Harris, Mathematical Methods for Physicists, 7th ed, pg. 199 for the unit
    vectors.

    Parameters
    ----------
    Bx, By, Bz : float, array_like
    theta : float, array_like
    phi : float, array_like

    Returns
    -------
    Br, Bth, Bphi : float, array_like
    """

    Br   =  np.sin(theta) * np.cos(phi) * Bx + np.sin(theta) * np.sin(phi) * By + np.cos(theta) * Bz
    Bth  =  np.cos(theta) * np.cos(phi) * Bx + np.cos(theta) * np.sin(phi) * By - np.sin(theta) * Bz
    Bphi =                 -np.sin(phi) * Bx +                 np.cos(phi) * By

    return Br, Bth, Bphi


def Bsph2Bxyz(Br, Bth, Bphi, theta, phi):
    """
    Convert vector components aligned to spherical coordinates into vector components aligned to
    Cartesian axes.
    Source: Arfken, Weber, Harris, Mathematical Methods for Physicists, 7th ed, pg. 199 for the unit
    vectors.

    Parameters
    ----------
    Br, Bth, Bphi : float, array_like
    theta : float, array_like
    phi : float, array_like

    Returns
    -------
    Bx, By, Bz : float, array_like
    """
    Bx =  np.sin(theta) * np.cos(phi) * Br + np.cos(theta) * np.cos(phi) * Bth - np.sin(phi) * Bphi
    By =  np.sin(theta) * np.sin(phi) * Br + np.cos(theta) * np.sin(phi) * Bth + np.cos(phi) * Bphi
    Bz =                np.cos(theta) * Br               - np.sin(theta) * Bth

    return Bx, By, Bz
