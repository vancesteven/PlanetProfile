""" Docstring explaining what we're doing here """

import numpy as np
import Utilities.PPversion as PPver
from Utilities.SwEOSChooser import SetupEOS

def SetupInit(Planet, Params, Constants):

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
    Params.dataFiles, Params.figureFiles = SetupFilenames(Planet)
    # Warn user if filename will round Tb_K value
    if Params.VERBOSE and round(Planet.Tb_K, 3) != Planet.Tb_K:
        print('WARNING: Planet.Tb_K has been rounded to generate saveFile name.')

    # Preallocate layer physical quantity arrays
    Planet.nStepsHydro, Layers = SetupLayers(Planet, Constants)

    return Planet, Params, Layers


def SetupFilenames(Planet):
    datPath = Planet.name + '/'
    figPath = Planet.name + '/figures/'

    saveBase = Planet.name + 'Profile_'
    if Planet.CLATHRATE: saveBase += 'Clathrates_'

    # Construct filenames for data, saving/reloading
    class dataFilesStruct:
        fName = datPath + saveBase + Planet.Ocean.comp + '_' + str(round(Planet.Ocean.wtOcean_ppt)) + 'WtPpt' \
            + '_Tb{:.3f}K'.format(Planet.Tb_K)
        if Planet.Silicate.mantleEOSName is not None: fName += Planet.Silicate.mantleEOSname
        if Planet.POROUS_ICE: fName += '_PorousIce'
        saveFile = fName + '.txt'
        saveFileMantle = fName + '_mantle.txt'

    dataFiles = dataFilesStruct()

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

    # Construct filenames for figures etc.
    class figureFilesStruct:
        pass
    figureFiles = figureFilesStruct()

    return dataFiles, figureFiles


class LayersStruct:
    def __init__(self, nTotal):
        self.nTotal = nTotal
        self.phase = np.zeros(nTotal, dtype=np.int_)
        self.T_K, self.P_MPa, self.rho_kgm3, self.z_m, \
        self.g_ms2, self.MAbove_kg, self.MBelow_kg, self.r_m \
            = (np.zeros(nTotal) for _ in range(8))


def SetupLayers(Planet, Constants):
    nStepsHydro = Planet.nStepsIceI + Planet.nStepsOcean

    Layers = LayersStruct(nStepsHydro)

    Layers.phase[:Planet.nStepsIceI] = 1 # Set ice Ih phase
    if Planet.nStepsClath is not None:
        nStepsHydro += Planet.nStepsClath
        Layers.phase = np.concatenate((np.zeros(Planet.nStepsClath) + 30, Layers.phase))  # Prepend clathrate phases

    Layers.r_m[0] = Planet.R_m # Set first layer to planetary surface
    Layers.g_ms2[0] = Constants.G * Planet.M_kg / Planet.R_m**2 # Set first layer to surface gravity
    Layers.T_K[0] = Planet.Tsurf_K # Set first layer to surface temp
    Layers.P_MPa[0] = Planet.Psurf_MPa # Set first layer to surface pressure

    return nStepsHydro, Layers
