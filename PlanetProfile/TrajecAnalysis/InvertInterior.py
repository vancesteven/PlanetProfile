"""
Perform a statistical inversion of interior structure parameters for a given set of assumptions from
magnetometer flyby data for selected spacecraft.
"""

import os
import logging
import importlib
import numpy as np
from copy import deepcopy
from hdf5storage import savemat, loadmat
from PlanetProfile.Main import WriteProfile, ReloadProfile, ExploreOgram, WriteExploreOgram, ReloadExploreOgram
from PlanetProfile.GetConfig import Params
from PlanetProfile.Utilities.SetupInit import SetupInversion
from PlanetProfile.Utilities.defineStructs import FitData, ModelDataStruct
from PlanetProfile.TrajecAnalysis.MagneticFields import InitModelData, SetupMagnetic, \
    CalcModelAndSumAll, CalcModelAmbient, CalcModelInduced, CalcModelPlasma
from PlanetProfile.MagneticInduction.MagneticInduction import GetBexc
from PlanetProfile.Utilities.SummaryTables import PrintTrajecFit, PrintTrajecTableLatex
from PlanetProfile.Plotting.TrajecPlots import PlotFlybys

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


def Invert(Params):
    """
    Invert for best-fit interior structure/constraints based on a set of input parameters.

    Parameters
    ----------
    Params : ParamsStruct
        A ParamsStruct class object containing trajectory analysis settings and parameters in
        the TrajecParamsStruct subclass object Params.Trajec.

    Returns
    -------
    Planet : PlanetStruct
        A PlanetStruct class object containing best-fit information.
    """

    Params.INVERSION_IN_PROGRESS = True

    # Load in MAG data from disk and prepare necessary Params
    Params, magData = SetupInversion(Params)

    if Params.CALC_NEW_INDUCT:
        # Fetch standard body file to get general settings, including radius and asymmetry
        InitPlanet = importlib.import_module(
            f'{Params.Trajec.targetBody}.PP{Params.Trajec.targetBody}').Planet
        InitPlanet.Magnetic = SetupMagnetic(InitPlanet, Params)

        # Initialize modelData object to be copied to each calc
        modelData = InitModelData(InitPlanet, Params, magData)

        # Run grid of models and compare against data
        PlanetGrid, InversionExploration = InvertOgram(InitPlanet, Params, magData, modelData)

        # Find best fit from among inversion exploration
        iFitPlanet = np.unravel_index(np.nanargmin(InversionExploration.RMSe), np.shape(InversionExploration.RMSe))
        FitPlanet = PlanetGrid[iFitPlanet]

        # Recalculate model evaluation for final selected model
        modelData = CalcModelAndSumAll(FitPlanet, Params, magData, modelData)
        modelData.fitProfileFname = FitPlanet.saveFile

        if not Params.NO_SAVEFILE:
            WriteProfile(FitPlanet, Params)
            WriteInversion(Params, modelData, InversionExploration)

    else:
        FitPlanet, Params, modelData, InversionExploration = ReloadInversion(Params)

    # Finished with calculations or reloading, turn flag off
    Params.INVERSION_IN_PROGRESS = False

    # Calculate and print details of selected fit model
    FitOutputs = FitCheck(Params, magData, modelData)

    # Plot results
    if not Params.SKIP_PLOTS:
        PlotFlybys(Params, magData, modelData)

    return FitPlanet, Params, FitOutputs, InversionExploration


def FitCheck(Params, magData, modelData):

    FitOutputs = FitData(Params, magData, modelData)
    PrintTrajecTableLatex(FitOutputs, Params)
    PrintTrajecFit(FitOutputs, Params)

    return FitOutputs


def InvertOgram(InitPlanet, Params, magData, modelData):

    # Run ExploreOgram over geophysical model ranges, returning raw grid
    PlanetGrid, InversionExploration, Params = ExploreOgram(InitPlanet.bodyname, Params,
                                                    RETURN_GRID=True, Magnetic=InitPlanet.Magnetic)

    # Save induction response info for each model
    invalidA = {scName: np.empty(InitPlanet.Magnetic.nExc[scName]) * np.nan for scName in InitPlanet.Magnetic.nExc.keys()}
    InversionExploration.Amp = np.array([[Planeti.Magnetic.Amp if Planeti.Magnetic.Amp is not None else invalidA for Planeti in line] for line in PlanetGrid])
    InversionExploration.phase = np.array([[Planeti.Magnetic.phase if Planeti.Magnetic.phase is not None else invalidA for Planeti in line] for line in PlanetGrid])

    # Get model data to copy over for when values are the same for each model
    modelData = CalcModelAmbient(PlanetGrid[0,0], Params, magData, modelData)
    modelData = CalcModelPlasma(PlanetGrid[0,0], Params, magData, modelData)

    # Evaluate fit for each model in grid
    InversionExploration = GridFitCalcs(PlanetGrid, InversionExploration, Params, magData, modelData)

    return PlanetGrid, InversionExploration


def GridFitCalcs(PlanetGrid, InversionExploration, Params, magData, modelDataIn):

    PlanetList1D = np.reshape(PlanetGrid, -1)
    FitFill = FitData(Params, magData, modelDataIn, VALID=False)
    FitDataList = np.array([deepcopy(FitFill) for _ in PlanetList1D], dtype=object)

    if Params.Trajec.plasmaType == 'none':
        FitFunc = FitModelInduced
    else:
        FitFunc = FitModelInducedAndPlasma

    if Params.DO_PARALLEL:
        # Prevent slowdowns from competing process spawning when #cores > #jobs
        nCores = np.min([Params.maxCores, np.prod(np.shape(PlanetList1D)), Params.threadLimit])
        pool = mtpContext.Pool(nCores)
        parResult = [pool.apply_async(FitFunc, (deepcopy(Planet), deepcopy(Params), deepcopy(magData), deepcopy(modelDataIn))) for Planet in
                     PlanetList1D]
        pool.close()
        pool.join()

        for i, result in enumerate(parResult):
            FitDataList[i] = result.get()

    else:
        log.profile('Running grid without parallel processing. This may take some time.')
        FitDataList = np.array(
            [FitFunc(deepcopy(Planet), deepcopy(Params), deepcopy(magData), deepcopy(modelDataIn))
             for Planet in PlanetList1D])

    FitDataGrid = np.reshape(FitDataList, np.shape(PlanetGrid))
    InversionExploration.RMSe = np.array([[FitDatai.RMSe['total'] for FitDatai in line] for line in FitDataGrid])
    InversionExploration.chiSquared = np.array([[FitDatai.chiSquared['total'] for FitDatai in line] for line in FitDataGrid])
    InversionExploration.stdDev = np.array([[FitDatai.stdDev['total'] for FitDatai in line] for line in FitDataGrid])
    InversionExploration.Rsquared = np.array([[FitDatai.Rsquared['total'] for FitDatai in line] for line in FitDataGrid])

    return InversionExploration


def FitModelInduced(Planet, Params, magData, modelData):
    # Ambient field and plasma field must be set in modelData or errors will result

    if Planet.Do.VALID:
        # Calculate stuff
        modelData.fitProfileFname = Planet.saveFile
        modelData = CalcModelInduced(Planet, Params, magData, modelData)

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

        fit = FitData(Params, magData, modelData)
    else:
        fit = FitData(Params, magData, modelData, VALID=False)

    return fit


def FitModelInducedAndPlasma(Planet, Params, magData, modelData):
    # Ambient field must be set in modelData or errors will result

    if Planet.Do.VALID:
        # Calculate stuff
        modelData.fitProfileFname = Planet.saveFile
        modelData = CalcModelInduced(Planet, Params, magData, modelData)
        modelData = CalcModelPlasma(Planet, Params, magData, modelData)

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

        fit = FitData(Params, magData, modelData)
    else:
        fit = FitData(Params, magData, modelData, VALID=False)

    return fit


def WriteInversion(Params, modelData, InversionExploration):

    outData = {
        'fitProfileFname': modelData.fitProfileFname,
        'allFlybys': modelData.allFlybys,
        'fbInclude': modelData.fbInclude,
        't_UTC': modelData.t_UTC,
        'ets': modelData.ets,
        'BxIAU_nT': modelData.BxIAU_nT,
        'ByIAU_nT': modelData.ByIAU_nT,
        'BzIAU_nT': modelData.BzIAU_nT,
        'BxIAUexc_nT': modelData.BxIAUexc_nT,
        'ByIAUexc_nT': modelData.ByIAUexc_nT,
        'BzIAUexc_nT': modelData.BzIAUexc_nT,
        'BxIAUind_nT': modelData.BxIAUind_nT,
        'ByIAUind_nT': modelData.ByIAUind_nT,
        'BzIAUind_nT': modelData.BzIAUind_nT,
        'BxIAUpls_nT': modelData.BxIAUpls_nT,
        'ByIAUpls_nT': modelData.ByIAUpls_nT,
        'BzIAUpls_nT': modelData.BzIAUpls_nT,
        'x_Rp': modelData.x_Rp,
        'y_Rp': modelData.y_Rp,
        'z_Rp': modelData.z_Rp,
        'r_Rp': modelData.r_Rp,
        'wireAlfven_RP': modelData.wireAlfven_RP
    }

    if os.path.exists(Params.Trajec.DataFiles.Btrajec):
        os.remove(Params.Trajec.DataFiles.Btrajec)
    savemat(Params.Trajec.DataFiles.Btrajec, outData)

    WriteExploreOgram(InversionExploration, Params, INVERSION=True)

    return


def ReloadInversion(Params):

    data = loadmat(Params.Trajec.DataFiles.Btrajec)
    modelData = ModelDataStruct(loadDict=data)

    # Concatenate net field data
    modelData.BxAll_nT, modelData.ByAll_nT, modelData.BzAll_nT = (np.empty(0) for _ in range(3))
    for scName, ets in modelData.ets.items():
        for fbID in ets.keys():
            modelData.BxAll_nT = np.concatenate((modelData.BxAll_nT, modelData.BxIAU_nT[scName][fbID]))
            modelData.ByAll_nT = np.concatenate((modelData.ByAll_nT, modelData.ByIAU_nT[scName][fbID]))
            modelData.BzAll_nT = np.concatenate((modelData.BzAll_nT, modelData.BzIAU_nT[scName][fbID]))

    Planet, Params = ReloadProfile(None, Params, fnameOverride=modelData.fitProfileFname)

    InversionExploration = ReloadExploreOgram(Planet.bodyname, Params, INVERSION=True)

    return Planet, Params, modelData, InversionExploration


if __name__ == '__main__':
    Planet, Params, FitOutputs, InversionExploration = Invert(Params)
    a=0
