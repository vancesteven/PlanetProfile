import numpy as np
import config as cfg
#from Thermodynamics.FromLiterature.conductiveMantleTemperature import conductiveMantleTemperature
#from Thermodynamics.FromLiterature.ConvectionDeschampsSotin2001 import ConvectionDeschampsSotin2001
from MantleSizePlot import MantleSizePlot
from CoreSizePlot import CoreSizePlot
import MatToPy

def writeProfile(path,saveStr,header,data):
    with open(path+saveStr+".txt","w") as f:
        f.write(header+"\n")
        for line in data:
            f.write("\t".join( [ "{:3.5e}".format(val) for val in line] )+"\n")

def PlanetProfile(Planet, Seismic, Params):
    nTbs = len(Planet)
    savebase = Planet[0]['name'] + '/' + Planet[0]['name'] + 'Profile_'
    figbase = Planet[0]['name'] + '/figures/' + Planet[0]['name']

    mantleSizeR, mantleSizeRho = ( np.zeros((nTbs,Params['nsteps_mantle'])) for _ in range(2) )
    for iT in range(nTbs):
        if 'clathrate' in Planet[iT]: savebase += 'Clathrates_'

        thisMantleSizePath = savebase + cfg.vmant + str(Planet[iT]['Tb_K']) + '.csv'
        thisMantleSizeR, thisMantleSizeRho = np.loadtxt(thisMantleSizePath, skiprows=1, unpack=True, delimiter=",")
        mantleSizeR[iT,:len(thisMantleSizeR)] = thisMantleSizeR
        mantleSizeRho[iT,:len(thisMantleSizeRho)] = thisMantleSizeRho

    MantleSizePlot(mantleSizeRho, mantleSizeR, Planet, nTbs, figbase+cfg.vmant, show=False)
    outPlanet = Planet # Placeholder
    return outPlanet