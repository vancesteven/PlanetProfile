import logging
import numpy as np
from PlanetProfile.TrajecAnalysis.SpiceFuncs import RotateFrame, RotateFrameManual, spiceS3coords
from PlanetProfile.TrajecAnalysis.FlybyEvents import GetFlybyCA
from PlanetProfile.Utilities.defineStructs import Constants
from PlanetProfile.GetConfig import FigMisc

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

def AlfvenWingField(Planet, Params, magData, modelData):
    # Calculate shape, current distribution, and magnetic field from Alfven wings
    # after Khurana et al. (1997): https://doi.org/10.1029/97GL02507
    RionosTop_km = (Planet.Bulk.R_m + Planet.Magnetic.ionosBounds_m[
        -1]) / 1e3  # Radius of the top of the ionosphere, where Alfven wings intersect it
    rhoRings_km = RionosTop_km * np.cos(np.radians(Planet.Magnetic.latAlfvenIntersect_deg))
    zRings_km = RionosTop_km * np.sin(np.radians(Planet.Magnetic.latAlfvenIntersect_deg))
    phiStarts, DeltaPhi = np.linspace(0, 2 * np.pi, Params.Trajec.nWiresAlfven, endpoint=False,
                                      retstep=True)
    phiEnds = phiStarts + DeltaPhi
    phiMids = phiStarts + DeltaPhi/2

    S3frame = spiceS3coords[Planet.parent]
    IAUframe = f'IAU_{Planet.bodyname.upper()}'
    phiOframe = f'{Planet.bodyname.upper()}_PHI_OMEGA'

    for scName in modelData.ets.keys():
        modelData.wireAlfven_RP[scName], modelData.BxIAUpls_nT[scName], modelData.ByIAUpls_nT[scName], modelData.BzIAUpls_nT[scName], \
            modelData.IwireRings_A[scName], modelData.IwireWings_A[scName] = ({} for _ in range(6))

        for fbID in modelData.ets[scName].keys():
            log.debug(f'Calculating Alfven wing field for {scName} {Planet.bodyname} encounter {modelData.allFlybys[scName][fbID]}.')

            modelData.wireAlfven_RP[scName][fbID] = {}
            if fbID not in Planet.Magnetic.IoverImax[scName].keys():
                Planet.Magnetic.IoverImax[scName][fbID] = Planet.Magnetic.IoverImaxDefault
            if fbID not in Planet.Magnetic.nPlasmaAmbient_pcc[scName].keys():
                Planet.Magnetic.nPlasmaAmbient_pcc[scName][fbID] = Planet.Magnetic.nPlasmaAmbientDefault_pcc

            npts = np.size(modelData.ets[scName][fbID])
            # Get index of et nearest CA
            FlybyCA = GetFlybyCA()
            iCA = np.argmin(modelData.ets[scName][fbID] - FlybyCA[scName].etCA[Planet.bodyname][fbID])

            if Params.Trajec.FIXED_ALFVEN:
                # Get background field needed to find the direction of Alfven wing (often called "Alfven characteristic") at CA
                BxyzMoonPhiO_nT = np.squeeze(RotateFrame(np.vstack([
                    magData[scName].BxS3moon_nT[fbID][iCA],
                    magData[scName].ByS3moon_nT[fbID][iCA],
                    magData[scName].BzS3moon_nT[fbID][iCA]
                ]).T, modelData.ets[scName][fbID][iCA], S3frame, phiOframe))
                BmoonBGmag_nT = np.sqrt(BxyzMoonPhiO_nT[0]**2 + BxyzMoonPhiO_nT[1]**2 + BxyzMoonPhiO_nT[2]**2)
                uBxMoonPhiO = BxyzMoonPhiO_nT[0] / BmoonBGmag_nT
                uByMoonPhiO = BxyzMoonPhiO_nT[1] / BmoonBGmag_nT
                uBzMoonPhiO = BxyzMoonPhiO_nT[2] / BmoonBGmag_nT
                wireStartUpperPhiO_km, wireStartLowerPhiO_km, ringMidsUpperPhiO_km, ringMidsLowerPhiO_km, \
                    ringMidsEqPhiO_km \
                    = GetWingIntersectsFixed(Params.Trajec.nWiresAlfven, RionosTop_km, rhoRings_km,
                    zRings_km, phiStarts, phiMids, uBzMoonPhiO)


            else:
                # Get background field needed to find the direction of Alfven wing (often called "Alfven characteristic")
                BxyzMoonPhiO_nT = RotateFrame(np.ascontiguousarray(np.vstack((
                    magData[scName].BxS3moon_nT[fbID],
                    magData[scName].ByS3moon_nT[fbID],
                    magData[scName].BzS3moon_nT[fbID]
                )).T), modelData.ets[scName][fbID], S3frame, phiOframe)
                BmoonBGmag_nT = np.sqrt(BxyzMoonPhiO_nT[:,0]**2 + BxyzMoonPhiO_nT[:,1]**2 + BxyzMoonPhiO_nT[:,2]**2)
                uBxMoonPhiO = BxyzMoonPhiO_nT[:, 0] / BmoonBGmag_nT
                uByMoonPhiO = BxyzMoonPhiO_nT[:, 1] / BmoonBGmag_nT
                uBzMoonPhiO = BxyzMoonPhiO_nT[:, 2] / BmoonBGmag_nT
                wireStartUpperPhiO_km, wireStartLowerPhiO_km, ringMidsUpperPhiO_km, ringMidsLowerPhiO_km, \
                    ringMidsEqPhiO_km \
                    = GetWingIntersectsVariable(Params.Trajec.nWiresAlfven, RionosTop_km, rhoRings_km,
                    zRings_km, phiStarts, phiMids, uBzMoonPhiO[iCA])

            # Transform SC trajectory to PhiO frame
            xyzPhiO_km = RotateFrame(np.ascontiguousarray(np.vstack((
                magData[scName].x_km[fbID],
                magData[scName].y_km[fbID],
                magData[scName].z_km[fbID]
            )).T), modelData.ets[scName][fbID], IAUframe, phiOframe)
            rhoPhiO_km = np.sqrt(xyzPhiO_km[:, 0]**2 + xyzPhiO_km[:, 1]**2)

            # Get direction of each characteristic in PhiO frame
            rhoPlasmaAmbient_kgm3 = Planet.Magnetic.nPlasmaAmbient_pcc[scName][fbID] * 1e3 \
                                    * Planet.Magnetic.mPlasmaAmbientAvg_gmol / Constants.NAvo
            vAlfven_kms = BmoonBGmag_nT * 1e-12 / np.sqrt(Constants.mu0 * rhoPlasmaAmbient_kgm3)
            dirAlfvenxUpper = Planet.Magnetic.vPlasmaAmbient_kms - uBxMoonPhiO * vAlfven_kms
            dirAlfvenyUpper = -uByMoonPhiO * vAlfven_kms
            dirAlfvenzUpper = -uBzMoonPhiO * vAlfven_kms
            AlfvenUpperNorm = np.sqrt(dirAlfvenxUpper**2 + dirAlfvenyUpper**2 + dirAlfvenzUpper**2)
            dirAlfvenxLower = Planet.Magnetic.vPlasmaAmbient_kms + uBxMoonPhiO * vAlfven_kms
            dirAlfvenyLower = uByMoonPhiO * vAlfven_kms
            dirAlfvenzLower = uBzMoonPhiO * vAlfven_kms
            AlfvenLowerNorm = np.sqrt(dirAlfvenxLower**2 + dirAlfvenyLower**2 + dirAlfvenzLower**2)
            uAlfvenUpperPhiO = np.ascontiguousarray(np.vstack((
                dirAlfvenxUpper / AlfvenUpperNorm,
                dirAlfvenyUpper / AlfvenUpperNorm,
                dirAlfvenzUpper / AlfvenUpperNorm
            )).T)
            uAlfvenLowerPhiO = np.ascontiguousarray(np.vstack((
                dirAlfvenxLower / AlfvenLowerNorm,
                dirAlfvenyLower / AlfvenLowerNorm,
                dirAlfvenzLower / AlfvenLowerNorm
            )).T)

            # Transform coordinates to wire frames
            angUpper_rad = np.arccos(uAlfvenUpperPhiO[:, 2])
            angLower_rad = np.arccos(uAlfvenLowerPhiO[:, 2])
            rotAxisUpper = np.squeeze(np.cross([0, 0, 1], uAlfvenUpperPhiO))
            rotAxisLower = np.squeeze(np.cross([0, 0, 1], uAlfvenLowerPhiO))
            xyzWireUpper_km = np.array(
                [RotateFrameManual(xyzPhiO_km, rotAxisUpper, angUpper_rad, O=wireStart) for
                 wireStart in wireStartUpperPhiO_km])
            xyzWireLower_km = np.array(
                [RotateFrameManual(xyzPhiO_km, rotAxisLower, angLower_rad, O=wireStart) for
                 wireStart in wireStartLowerPhiO_km])

            # Calculate current in each wire
            Imax_A = 4 * Planet.Magnetic.vPlasmaAmbient_kms * RionosTop_km * 1e6 * np.sqrt(
                rhoPlasmaAmbient_kgm3 / Constants.mu0)
            Itotal_A = Imax_A * Planet.Magnetic.IoverImax[scName][fbID]
            Iwire_A = 1/2 * Itotal_A * np.sin(phiStarts)
            IwireRing_A = 1/2 * Itotal_A/3 * np.cos(phiStarts)
            modelData.IwireRings_A[scName][fbID] = IwireRing_A
            modelData.IwireWings_A[scName][fbID] = Iwire_A

            if Params.Trajec.FIXED_ALFVEN:
                # Calculate wire positions in PhiO frame for plotting purposes, but only when there's only
                # one Alfven wing to calculate per encounter
                ringStartEqPhiO_km = np.ascontiguousarray(np.vstack((
                    RionosTop_km * np.cos(phiStarts),
                    RionosTop_km * np.sin(phiStarts),
                    np.zeros(Params.Trajec.nWiresAlfven)
                )).T)
                wireEndUpperPhiO_km = np.squeeze(np.array(
                    [wireStart + uAlfvenUpperPhiO * Planet.Bulk.R_m / 1e3 * FigMisc.LAlfven_RP for
                     wireStart in wireStartUpperPhiO_km]))
                wireEndLowerPhiO_km = np.squeeze(np.array(
                    [wireStart + uAlfvenLowerPhiO * Planet.Bulk.R_m / 1e3 * FigMisc.LAlfven_RP for
                     wireStart in wireStartLowerPhiO_km]))
                R_km = Planet.Bulk.R_m / 1e3
                # Transform wire locations to IAU coordinates and store in modelData
                modelData.wireAlfven_RP[scName][fbID]['upperRing'] = np.squeeze(np.array(
                    [RotateFrame(wireStart[np.newaxis, :], modelData.ets[scName][fbID][iCA], phiOframe, IAUframe) for
                     wireStart in wireStartUpperPhiO_km])) / R_km
                modelData.wireAlfven_RP[scName][fbID]['lowerRing'] = np.squeeze(np.array(
                    [RotateFrame(wireStart[np.newaxis, :], modelData.ets[scName][fbID][iCA], phiOframe, IAUframe) for
                     wireStart in wireStartLowerPhiO_km])) / R_km
                modelData.wireAlfven_RP[scName][fbID]['eqRing'] = np.squeeze(np.array(
                    [RotateFrame(wireStart[np.newaxis, :], modelData.ets[scName][fbID][iCA], phiOframe, IAUframe) for
                     wireStart in ringStartEqPhiO_km])) / R_km
                modelData.wireAlfven_RP[scName][fbID]['upperEnd'] = np.squeeze(np.array(
                    [RotateFrame(wireStart[np.newaxis, :], modelData.ets[scName][fbID][iCA], phiOframe, IAUframe) for
                     wireStart in wireEndUpperPhiO_km])) / R_km
                modelData.wireAlfven_RP[scName][fbID]['lowerEnd'] = np.squeeze(np.array(
                    [RotateFrame(wireStart[np.newaxis, :], modelData.ets[scName][fbID][iCA], phiOframe, IAUframe) for
                     wireStart in wireEndLowerPhiO_km])) / R_km

            # Use Biot-Savart to determine Bxy at SC position in wing frames from wing currents
            # Start with catch for if we are very close to any wire
            rhoWireUpper_km = np.sqrt(xyzWireUpper_km[..., 0]**2 + xyzWireUpper_km[..., 1]**2)
            thetaWireUpper_rad = np.arctan(rhoWireUpper_km / xyzWireUpper_km[..., 2])
            rhoWireLower_km = np.sqrt(xyzWireLower_km[..., 0]**2 + xyzWireLower_km[..., 1]**2)
            thetaWireLower_rad = np.arctan2(rhoWireLower_km, xyzWireLower_km[..., 2])
            phiWireUpper_rad = np.arctan2(xyzWireUpper_km[..., 1], xyzWireUpper_km[..., 0])
            phiWireLower_rad = np.arctan2(xyzWireLower_km[..., 1], xyzWireLower_km[..., 0])

            # Get prefactors common to all wires
            # 1e6 factor in the below is for 1e-3 from m -> km in denominator and 1e-9 for T -> nT
            BwirePre_nTkm = Constants.mu0 * Iwire_A / 4/np.pi * 1e6
            # Initialize Bphi
            BphiWingUpper_nT = np.zeros_like(rhoWireUpper_km)
            BphiWingLower_nT = np.zeros_like(rhoWireLower_km)
            for i, Bpre in enumerate(BwirePre_nTkm):
                # Evaluate near-wire points
                jCloseToWireUpper = rhoWireUpper_km[i, :] < Planet.Magnetic.dWireAlfven_RP / 2
                jCloseToWireLower = rhoWireLower_km[i, :] < Planet.Magnetic.dWireAlfven_RP / 2
                BphiWingUpper_nT[i, jCloseToWireUpper] = Bpre \
                    * (1 + np.cos(thetaWireUpper_rad[i, jCloseToWireUpper])) * 8 \
                    * rhoWireUpper_km[i, jCloseToWireUpper] / Planet.Magnetic.dWireAlfven_RP ** 2
                BphiWingLower_nT[i, jCloseToWireLower] = Bpre \
                    * (1 + np.cos(thetaWireLower_rad[i, jCloseToWireLower])) * 8 \
                    * rhoWireLower_km[i, jCloseToWireLower] / Planet.Magnetic.dWireAlfven_RP ** 2
                # Evaluate remaining points
                jNotCloseUpper = np.logical_not(jCloseToWireUpper)
                jNotCloseLower = np.logical_not(jCloseToWireLower)
                BphiWingUpper_nT[i, jNotCloseUpper] = Bpre \
                    * (1 + np.cos(thetaWireUpper_rad[i, jNotCloseUpper])) / rhoWireUpper_km[i, jNotCloseUpper]
                BphiWingLower_nT[i, jNotCloseLower] = Bpre \
                    * (1 + np.cos(thetaWireLower_rad[i, jNotCloseLower])) / rhoWireLower_km[i, jNotCloseLower]

            # Get Bx, By in PhiO frame from wire frames
            BxWingUpper_nT = np.sum(BphiWingUpper_nT * np.cos(phiWireUpper_rad), axis=0)
            ByWingUpper_nT = np.sum(BphiWingUpper_nT * np.sin(phiWireUpper_rad), axis=0)
            BxWingLower_nT = np.sum(BphiWingLower_nT * np.cos(phiWireLower_rad), axis=0)
            ByWingLower_nT = np.sum(BphiWingLower_nT * np.sin(phiWireLower_rad), axis=0)

            # Use rotation-only condition to find Bz component in wing frames
            BzWingUpper_nT = -np.sqrt(BmoonBGmag_nT**2 - BxWingUpper_nT**2 - ByWingUpper_nT**2)
            BzWingLower_nT = -np.sqrt(BmoonBGmag_nT**2 - BxWingLower_nT**2 - ByWingLower_nT**2)

            # Transform back to PhiO frame
            BxyzWingUpper_nT = np.ascontiguousarray(np.vstack((
                BxWingUpper_nT,
                ByWingUpper_nT,
                BzWingUpper_nT
            )).T)
            BxyzWingLower_nT = np.ascontiguousarray(np.vstack((
                BxWingLower_nT,
                ByWingLower_nT,
                BzWingLower_nT
            )).T)
            BwingUpperPhiO_nT = RotateFrameManual(BxyzWingUpper_nT, rotAxisUpper, -angUpper_rad)
            BwingLowerPhiO_nT = RotateFrameManual(BxyzWingLower_nT, rotAxisLower, -angLower_rad)

            # Add contribution from ring currents using Biot-Savart in PhiO frame
            rEvalUpper_km = np.array([xyzPhiO_km - ringMidsUpper for ringMidsUpper in ringMidsUpperPhiO_km])
            rEvalLower_km = np.array([xyzPhiO_km - ringMidsLower for ringMidsLower in ringMidsLowerPhiO_km])
            rEvalMid_km   = np.array([xyzPhiO_km - ringMidsEq    for ringMidsEq    in ringMidsEqPhiO_km])
            rFactorUpper_km2 = np.sqrt(rEvalUpper_km[..., 0]**2 + rEvalUpper_km[..., 1]**2 \
                                       + rEvalUpper_km[..., 2]**2)**3 / rhoRings_km
            rFactorLower_km2 = np.sqrt(rEvalLower_km[..., 0]**2 + rEvalLower_km[..., 1]**2 \
                                       + rEvalLower_km[..., 2]**2)**3 / rhoRings_km
            rFactorMid_km2 = np.sqrt(rEvalMid_km[..., 0]**2 + rEvalMid_km[..., 1]**2 \
                                     + rEvalMid_km[..., 2]**2)**3 / RionosTop_km
            # 1e6 factor in the below is for 1e-3 from m -> km in denominator and 1e-9 for T -> nT
            B0rings_nTkm = Constants.mu0 / 4/np.pi * 1e6 * IwireRing_A
            xPhis = (DeltaPhi/2 + (np.sin(2*phiEnds) - np.sin(2*phiStarts)) / 4)
            yPhis = (np.cos(2*phiStarts) - np.cos(2*phiEnds)) / 4
            BxRings_nT = np.sum(np.array([
                B0rings_nTkm * (rEvalUpper_km[:, i, 2] - zRings_km) * xPhis / rFactorUpper_km2[:, i]
              + B0rings_nTkm * rEvalMid_km[:, i, 2] * xPhis / rFactorMid_km2[:, i]
              + B0rings_nTkm * (rEvalLower_km[:, i, 2] + zRings_km) * xPhis / rFactorLower_km2[:, i]
                for i in range(npts)]), axis=1)
            ByRings_nT = np.sum(B0rings_nTkm * np.array([
                (rEvalUpper_km[:, i, 2] - zRings_km) * yPhis / rFactorUpper_km2[:, i]
              + rEvalMid_km[:, i, 2] * yPhis / rFactorMid_km2[:, i]
              + (rEvalLower_km[:, i, 2] + zRings_km) * yPhis / rFactorLower_km2[:, i]
                for i in range(npts)]), axis=1)
            BzRings_nT = 2*np.pi * np.sum(B0rings_nTkm * np.array([
                (rhoRings_km - rho)    / rFactorUpper_km2[:, i]
              + (RionosTop_km - rho) / rFactorMid_km2[:, i]
              + (rhoRings_km - rho)    / rFactorLower_km2[:, i]
                for i, rho in enumerate(rhoPhiO_km)]), axis=1)

            # Sum wire current contributions and ring current contributions
            BxPhiO_nT = BwingUpperPhiO_nT[:, 0] + BxRings_nT + BwingLowerPhiO_nT[:, 0]
            ByPhiO_nT = BwingUpperPhiO_nT[:, 1] + ByRings_nT + BwingLowerPhiO_nT[:, 1]
            BzPhiO_nT = BwingUpperPhiO_nT[:, 2] + BzRings_nT + BwingLowerPhiO_nT[:, 2]

            # Transform back to IAU frame
            BxyzIAU_nT = RotateFrame(np.ascontiguousarray(np.vstack((
                BxPhiO_nT,
                ByPhiO_nT,
                BzPhiO_nT
            )).T), modelData.ets[scName][fbID], phiOframe, IAUframe)

            # Record in modelData
            modelData.BxIAUpls_nT[scName][fbID], modelData.ByIAUpls_nT[scName][fbID], modelData.BzIAUpls_nT[scName][fbID] = (
                BxyzIAU_nT[:, 0], BxyzIAU_nT[:, 1], BxyzIAU_nT[:, 2])

    return modelData


def GetWingIntersectsFixed(nWires, RionosTop_km, rhoRings_km, zRings_km, phiStarts, phiMids, uBzMoonPhiO):
    
    onePhis = np.ones(nWires)
    # Get the points of ionosphere--current wire intersect points in PhiO coordinates
    wireStartUpperPhiO_km = np.ascontiguousarray(np.vstack((
        rhoRings_km * np.cos(phiStarts),
        rhoRings_km * np.sin(phiStarts),
        zRings_km * -np.sign(uBzMoonPhiO) * onePhis
    )).T)
    wireStartLowerPhiO_km = wireStartUpperPhiO_km + np.zeros_like(wireStartUpperPhiO_km)
    wireStartLowerPhiO_km[:, 2] *= -1
    # Get eval locations for points in ring currents that close Alfven wing currents
    ringMidsUpperPhiO_km = np.ascontiguousarray(np.vstack((
        rhoRings_km * np.cos(phiMids),
        rhoRings_km * np.sin(phiMids),
        zRings_km * -np.sign(uBzMoonPhiO) * onePhis
    )).T)
    ringMidsLowerPhiO_km = ringMidsUpperPhiO_km + np.zeros_like(ringMidsUpperPhiO_km)
    ringMidsLowerPhiO_km[:, 2] *= -1
    ringMidsEqPhiO_km = np.ascontiguousarray(np.vstack((
        RionosTop_km * np.cos(phiMids),
        RionosTop_km * np.sin(phiMids),
        onePhis * 0
    )).T)

    return wireStartUpperPhiO_km, wireStartLowerPhiO_km, ringMidsUpperPhiO_km, ringMidsLowerPhiO_km, \
        ringMidsEqPhiO_km


def GetWingIntersectsVariable(nWires, RionosTop_km, rhoRings_km, zRings_km, phiStarts, phiMids, uBzMoonPhiO):
    
    onePhis = np.ones(nWires)
    # Get the points of ionosphere--current wire intersect points in PhiO coordinates
    wireStartUpperPhiO_km = np.ascontiguousarray(np.vstack((
        rhoRings_km * np.cos(phiStarts),
        rhoRings_km * np.sin(phiStarts),
        zRings_km * -np.sign(uBzMoonPhiO) * onePhis
    )).T)
    wireStartLowerPhiO_km = wireStartUpperPhiO_km + np.zeros_like(wireStartUpperPhiO_km)
    wireStartLowerPhiO_km[:, 2] *= -1
    # Get eval locations for points in ring currents that close Alfven wing currents
    ringMidsUpperPhiO_km = np.ascontiguousarray(np.vstack((
        rhoRings_km * np.cos(phiMids),
        rhoRings_km * np.sin(phiMids),
        zRings_km * -np.sign(uBzMoonPhiO) * onePhis
    )).T)
    ringMidsLowerPhiO_km = ringMidsUpperPhiO_km + np.zeros_like(ringMidsUpperPhiO_km)
    ringMidsLowerPhiO_km[:, 2] *= -1
    ringMidsEqPhiO_km = np.ascontiguousarray(np.vstack((
        RionosTop_km * np.cos(phiMids),
        RionosTop_km * np.sin(phiMids),
        onePhis * 0
    )).T)
    
    return wireStartUpperPhiO_km, wireStartLowerPhiO_km, ringMidsUpperPhiO_km, ringMidsLowerPhiO_km, \
        ringMidsEqPhiO_km
