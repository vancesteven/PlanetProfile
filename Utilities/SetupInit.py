""" Docstring explaining what we're doing here """

import os
import numpy as np
import Utilities.PPversion as PPver
from Utilities.SwEOSChooser import SetupEOS

def SetupInit(Planet, Params):

    # Print version number
    verNum = PPver.verNum
    print('-- PlanetProfile v' + verNum + ' --')
    if verNum[-3:] == 'dev': print('This version is in development.')

    # Check dependency compatibility
    PPver.CheckSeaFreeze(PPver.seaFreezeCompatVer)  # SeaFreeze
    PPver.CheckGSW(PPver.gswCompatVer)  # Gibbs Seawater
    PPver.CheckTauP(PPver.taupCompatVer)  # TauP

    # Get ocean equation of state
    Planet.oceanEOS = SetupEOS(Planet.Ocean.comp)

    # Get filenames for saving/loading
    Params.dataFiles, Params.figureFiles = SetupFilenames(Planet, Params)
    # Warn user if filename will round Tb_K value
    if Params.VERBOSE and round(Planet.Bulk.Tb_K, 3) != Planet.Bulk.Tb_K:
        print('WARNING: Planet.Tb_K has been rounded to generate saveFile name.')

    # Set number of steps for unused options to zero
    if not Planet.Do.CLATHRATE:
        Planet.Steps.nClath = 0
        Planet.zClath_m = 0
    if not Planet.Do.BOTTOM_ICEIII and not Planet.Do.BOTTOM_ICEV:
        Planet.Steps.nIceIIILitho = 0
        Planet.Steps.nIceVLitho = 0
    elif not Planet.Do.BOTTOM_ICEV:
        Planet.Steps.nIceVLitho = 0
    if not Planet.Do.POROUS_ROCK:
        Planet.Sil.phiRockMax = 0

    # Initialize calculated quantities
    if not Planet.Do.EQUIL_Q:
        # Assign heat input into ocean from mantle to be ~radiogenic
        print('WARNING: QfromMantle_Wm2 is set to a value consistent only with Europa radiogenic heating.')
        Planet.Ocean.QfromMantle_Wm2 = 2.2e11/ 4/np.pi / Planet.Bulk.R_m**2

    # Preallocate layer physical quantity arrays
    Planet = SetupLayers(Planet)

    # Delete the next 6 lines after their calculation functions have been implemented!
    Planet.Steps.nHydro = Planet.Steps.nHydroMax + 0
    Planet.C2mean = 0.0
    Planet.Sil.RsilMean_m = 0.0
    Planet.Sil.RsilRange_m = 0.0
    Planet.Core.RFeMean_m = 0.0
    Planet.Core.RFeRange_m = 0.0

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
    dataFiles = dataFilesStruct(datBase)

    # Figure filename strings
    vsP = 'Porosity_vs_P'
    vsR = 'Porosity_vs_R'
    vperm = 'Permeability'
    vgsks = 'Gs_Ks'
    vseis = 'Seismic'
    vcond = 'Conductivity'
    vgrav = 'Gravity'
    vmant = 'MantleDens'
    vcore = 'CoreMantTrade'
    vpvt4 = 'PTx4'
    vpvt6 = 'PTx6'
    vwedg = 'Wedge'

    figBase = os.path.join(figPath, saveBase)
    figureFiles = figureFilesStruct(figBase, Params.xtn)

    #vcondFig.savefig(Params.figureFields.vcond, format = "png", dpi = 200)
    return dataFiles, figureFiles


def SetupLayers(Planet):

    nOceanMax = int(Planet.Bulk.PHydroMax_MPa / Planet.Bulk.deltaP)
    Planet.Steps.nHydroMax = Planet.Steps.nClath + Planet.Steps.nIceI + Planet.Steps.nIceIIILitho + Planet.Steps.nIceVLitho + nOceanMax

    Planet.phase = np.zeros(Planet.Steps.nHydroMax, dtype=np.int_)
    Planet.z_m, Planet.r_m, Planet.P_MPa, Planet.T_K, Planet.rho_kgm3, \
    Planet.Cp_JkgK, Planet.sigma_Sm, Planet.g_ms2, Planet.vFluid_kms, \
    Planet.alpha_pK = \
        (np.zeros(Planet.Steps.nHydroMax) for _ in range(10))

    return Planet

# Construct filenames for data, saving/reloading
class dataFilesStruct:
    def __init__(self, fName):
        self.saveFile = fName + '.txt'
        self.mantCoreFile = fName + '_mantleCore.txt'
        self.permFile = fName + '_mantlePerm.txt'

# Construct filenames for figures etc.
class figureFilesStruct:
    def __init__(self, fName, xtn):
        self.dummyFig = fName + xtn
