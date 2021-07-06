import sys
import numpy as np
import importlib

import config as cfg
#from Thermodynamics.FromLiterature.conductiveMantleTemperature import conductiveMantleTemperature
#from Thermodynamics.FromLiterature.ConvectionDeschampsSotin2001 import ConvectionDeschampsSotin2001
from MantleSizePlot import MantleSizePlot
from CoreSizePlot import CoreSizePlot
import MatToPy
from Utilities.PPversion import verNum

def main():
    # Intro
    print("-- PlanetProfile v" + verNum + " --")
    if verNum[-3:] == "dev": print("This version is in development.")

    # Command line args
    if len(sys.argv) > 1:
        # Body name was passed as command line argument
        bodyname = sys.argv[1]
    else:
        # No command line argument, ask user which body to run
        bodyname = input("Please input body name: ")

    bodyname = bodyname.capitalize()
    body = importlib.import_module(bodyname+".PP"+bodyname)

    outPlanet = PlanetProfile(body.Planet, body.Seismic, body.Params)

def writeProfile(path,saveStr,header,data):
    with open(path+saveStr+".txt","w") as f:
        f.write(header+"\n")
        for line in data:
            f.write("\t".join( [ "{:3.5e}".format(val) for val in line] )+"\n")

def PlanetProfile(Planet, Seismic, Params):
    nTbs = len(Planet)
    nMantleInds = np.zeros(nTbs, dtype=np.int_)

    savebase = Planet[0]['name'] + '/' + Planet[0]['name'] + 'Profile_'
    figbase = Planet[0]['name'] + '/figures/' + Planet[0]['name']

    mantleSizeR, mantleSizeRho = ( np.zeros((nTbs,Params['nsteps_mantle'])) for _ in range(2) )
    for iT in range(nTbs):
        if 'clathrate' in Planet[iT]: savebase += 'Clathrates_'

        thisMantleSizePath = savebase + cfg.vmant + str(Planet[iT]['Tb_K']) + '.csv'
        thisMantleSizeR, thisMantleSizeRho = np.loadtxt(thisMantleSizePath, skiprows=1, unpack=True, delimiter=",")
        nMantleInds[iT] = len(thisMantleSizeR)
        mantleSizeR[iT,:nMantleInds[iT]] = thisMantleSizeR
        mantleSizeRho[iT,:nMantleInds[iT]] = thisMantleSizeRho

    MantleSizePlot(mantleSizeRho, mantleSizeR, Planet, nTbs, nMantleInds, figbase+cfg.vmant, show=False)
    outPlanet = Planet # Placeholder
    return outPlanet


if __name__ == '__main__': main()