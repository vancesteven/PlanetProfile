""" Docstring explaining what we're doing here """

import os
import numpy as np
import Utilities.PPversion as PPver
from Utilities.SwEOSChooser import SetupEOS
from Utilities.dataStructs import Constants

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
    Planet.Steps.nHydroMax = Planet.Steps.nIceI + nOceanMax

    if Planet.Do.CLATHRATE:
        Planet.Steps.nHydroMax += Planet.Steps.nClath
    else:
        Planet.Steps.nClath = 0

    Planet.phase = np.zeros(Planet.Steps.nHydroMax, dtype=np.int_)
    Planet.z_m, Planet.r_m, Planet.P_MPa, Planet.T_K, Planet.rho_kgm3, \
    Planet.Cp_JkgK, Planet.sigma_Sm, Planet.g_ms2, Planet.vFluid_kms, \
    Planet.MAbove_kg, Planet.MBelow_kg = \
        (np.zeros(Planet.Steps.nHydroMax) for _ in range(11))

    Planet.phase = np.concatenate((np.zeros(Planet.Steps.nClath, dtype=np.int_) + 30, Planet.phase))  # Prepend clathrate phases
    Planet.phase[Planet.Steps.nClath:Planet.Steps.nClath + Planet.Steps.nIceI] = 1  # Set ice Ih phase
    Planet.r_m[0] = Planet.Bulk.R_m  # Set first layer to planetary surface radius
    Planet.g_ms2[0] = Constants.G * Planet.Bulk.M_kg / Planet.Bulk.R_m**2  # Set first layer gravity at surface
    Planet.T_K[0] = Planet.Bulk.Tsurf_K  # Set first layer surface temp
    Planet.P_MPa[0] = Planet.Bulk.Psurf_MPa  # Set first layer to surface pressure

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
