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
        saveFile = datPath + saveBase + Planet.Ocean.comp + '_' + str(round(Planet.Ocean.wtOcean_ppt)) + 'WtPpt'
        if Planet.Silicate.mantleEOSName is not None: saveFile += Planet.Silicate.mantleEOSname
        if Planet.POROUS_ICE: saveFile += '_PorousIce'

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


def SetupLayers(Planet, Constants):
    nStepsHydro = Planet.nStepsIceI + Planet.nStepsOcean

    class LayersStruct:
        phase, T_K, P_MPa, rho_kgm3, z_m, g_ms2, MAbove_kg, MBelow_kg, r_m = (np.zeros(nStepsHydro) for _ in range(9))

    Layers = LayersStruct()

    Layers.phase[:Planet.nStepsIceI] = 1 # Set ice Ih phase
    if Planet.nStepsClath is not None:
        nStepsHydro += Planet.nStepsClath
        Layers.phase = np.concatenate((np.zeros(Planet.nStepsClath) + 30, Layers.phase))  # Prepend clathrate phases

    Layers.r_m[0] = Planet.R_m # Set first layer to planetary surface
    Layers.g_ms2[0] = Constants.G * Planet.M_kg / Planet.R_m**2 # Set first layer to surface gravity
    Layers.T_K[0] = Planet.Tsurf_K # Set first layer to surface temp
    Layers.P_MPa[0] = Planet.Psurf_MPa # Set first layer to surface pressure

    return nStepsHydro, Layers
