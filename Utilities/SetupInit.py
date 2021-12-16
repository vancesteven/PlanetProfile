""" Docstring explaining what we're doing here """

import os
import numpy as np
import Utilities.PPversion as PPver
from Utilities.dataStructs import DataFilesSubstruct, FigureFilesSubstruct

def SetupInit(Planet, Params):

    # Print version number
    verNum = PPver.verNum
    print('-- PlanetProfile v' + verNum + ' --')
    if verNum[-3:] == 'dev': print('This version is in development.')

    # Check dependency compatibility
    PPver.CheckSeaFreeze(PPver.seaFreezeCompatVer)  # SeaFreeze
    PPver.CheckGSW(PPver.gswCompatVer)  # Gibbs Seawater
    PPver.CheckTauP(PPver.taupCompatVer)  # TauP

    # Get filenames for saving/loading
    Params.DataFiles, Params.FigureFiles = SetupFilenames(Planet, Params)
    # Warn user if filename will round Tb_K value
    if Params.VERBOSE and round(Planet.Bulk.Tb_K, 3) != Planet.Bulk.Tb_K:
        print('WARNING: Planet.Tb_K has been rounded to generate saveFile name.')

    # Set number of steps for unused options to zero
    if not Planet.Do.Fe_CORE:
        Planet.Steps.nCore = 0
    if not Planet.Do.CLATHRATE:
        Planet.Steps.nClath = 0
        Planet.zClath_m = 0
    if not Planet.Do.BOTTOM_ICEIII and not Planet.Do.BOTTOM_ICEV:
        Planet.Steps.nIceIIILitho = 0
        Planet.Steps.nIceVLitho = 0
    elif not Planet.Do.BOTTOM_ICEV:
        Planet.Steps.nIceVLitho = 0
    if not Planet.Do.POROUS_ROCK:
        Planet.Sil.phiRockMax_frac = 0

    # Calculate bulk density from total mass and radius, and warn user if they specified density
    if Planet.Bulk.M_kg is None:
        Planet.Bulk.M_kg = Planet.Bulk.rho_kgm3 * (4/3*np.pi * Planet.Bulk.R_m**3)
    else:
        if Planet.Bulk.rho_kgm3 is not None and Params.VERBOSE:
            print('Both bulk mass and density were specified. Only one is required--density will be recalculated from bulk mass for consistency.')
        Planet.Bulk.rho_kgm3 = Planet.Bulk.M_kg / (4/3*np.pi * Planet.Bulk.R_m**3)

    # Preallocate layer physical quantity arrays
    Planet = SetupLayers(Planet)

    return Planet, Params


def SetupFilenames(Planet, Params):
    """
    """
    datPath = Planet.name
    figPath = os.path.join(Planet.name, 'figures')

    saveBase = Planet.name + 'Profile_'
    if Planet.Do.CLATHRATE: saveBase += 'Clathrates_'
    saveBase += Planet.Ocean.comp + '_' + str(round(Planet.Ocean.wOcean_ppt)) + 'WtPpt' \
        + '_Tb{:.3f}K'.format(Planet.Bulk.Tb_K)
    if Planet.Sil.mantleEOSName is not None: fName += Planet.Sil.mantleEOSname
    if Planet.Do.POROUS_ICE: fName += '_PorousIce'

    datBase = os.path.join(datPath, saveBase)
    DataFiles = DataFilesSubstruct(datBase)
    figBase = os.path.join(figPath, saveBase)
    FigureFiles = FigureFilesSubstruct(figBase, Params.xtn)

    #vcondFig.savefig(Params.figureFields.vcond, format = "png", dpi = 200)
    return DataFiles, FigureFiles


def SetupLayers(Planet):

    nOceanMax = int(Planet.Ocean.PHydroMax_MPa / Planet.Ocean.deltaP)
    Planet.Steps.nHydroMax = Planet.Steps.nClath + Planet.Steps.nIceI + Planet.Steps.nIceIIILitho + Planet.Steps.nIceVLitho + nOceanMax

    Planet.phase = np.zeros(Planet.Steps.nHydroMax, dtype=np.int_)
    Planet.P_MPa, Planet.T_K, Planet.r_m, Planet.rho_kgm3, \
        Planet.Cp_JkgK, Planet.alpha_pK, Planet.g_ms2, Planet.phi_frac, \
        Planet.sigma_Sm, Planet.z_m, Planet.MLayer_kg = \
        (np.zeros(Planet.Steps.nHydroMax) for _ in range(11))

    return Planet
