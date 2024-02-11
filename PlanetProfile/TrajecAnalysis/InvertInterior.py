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
from PlanetProfile.Main import PlanetProfile, WriteProfile, ReloadProfile
from PlanetProfile.GetConfig import Params
from PlanetProfile.Utilities.SetupInit import SetupInversion
from PlanetProfile.Utilities.defineStructs import FitData, ModelDataStruct
from PlanetProfile.TrajecAnalysis.MagneticFields import CalcModel, InitModelData, SetupMagnetic
from PlanetProfile.Utilities.SummaryTables import PrintTrajecFit, PrintTrajecTableLatex
from PlanetProfile.Plotting.TrajecPlots import PlotFlybys

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

    # Load in MAG data from disk and prepare necessary Params
    Params, magData = SetupInversion(Params)

    if Params.CALC_NEW_INDUCT:
        # Fetch standard body file to get general settings, including radius and asymmetry
        InitPlanet = importlib.import_module(
            f'{Params.Trajec.targetBody}.PP{Params.Trajec.targetBody}').Planet
        Magnetic = SetupMagnetic(InitPlanet, Params)
        modelData = InitModelData(InitPlanet, Magnetic, Params, magData)

        # Initialize modelData object to be copied to each calc

        # Placeholder before we implement inversion
        FitPlanet = importlib.import_module(
            f'{Params.Trajec.targetBody}.PP{Params.Trajec.targetBody}').Planet
        FitPlanet.Magnetic = deepcopy(Magnetic)

        # Recalculate model evaluation for final selected model
        FitPlanet, Params = PlanetProfile(FitPlanet, Params)
        modelData = CalcModel(FitPlanet, Params, magData, modelData)

        if not Params.NO_SAVEFILE:
            WriteProfile(FitPlanet, Params)
            WriteInversion(Params, modelData)

    else:
        FitPlanet, Params, modelData = ReloadInversion(Params)

    # Finished with calculations or reloading, turn flag off
    Params.INVERSION_IN_PROGRESS = False

    # Calculate and print details of selected fit model
    FitOutputs = FitCheck(Params, magData, modelData)

    # Plot results
    if not Params.SKIP_PLOTS:
        PlotFlybys(Params, magData, modelData)

    return FitPlanet, Params, FitOutputs


def FitCheck(Params, magData, modelData):

    FitOutputs = FitData(Params, magData, modelData)
    PrintTrajecTableLatex(FitOutputs, Params)
    PrintTrajecFit(FitOutputs, Params)

    return FitOutputs


def WriteInversion(Params, modelData):

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
        'r_Rp': modelData.r_Rp
    }

    if os.path.exists(Params.Trajec.DataFiles.Btrajec):
        os.remove(Params.Trajec.DataFiles.Btrajec)
    savemat(Params.Trajec.DataFiles.Btrajec, outData)

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

    Planet = ReloadProfile(None, Params, fnameOverride=modelData.fitProfileFname)

    return Planet, Params, modelData


if __name__ == '__main__':
    Planet, Params, FitOutputs = Invert(Params)
