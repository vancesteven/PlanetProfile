import sys
import numpy as np
import importlib

from Utilities.dataStructs import *
import config as cfg
#from Thermodynamics.FromLiterature.conductiveMantleTemperature import conductiveMantleTemperature
#from Thermodynamics.FromLiterature.ConvectionDeschampsSotin2001 import ConvectionDeschampsSotin2001
#from MantleSizePlot import MantleSizePlot
#from CoreSizePlot import CoreSizePlot
#import MatToPy
import Utilities.PPversion as PPver

""" MAIN RUN BLOCK """
def main():
    # Intro
    verNum = PPver.verNum
    print("-- PlanetProfile v" + verNum + " --")
    if verNum[-3:] == "dev": print("This version is in development.")

    # Command line args
    if len(sys.argv) > 1:
        # Body name was passed as command line argument
        bodyname = sys.argv[1]
    else:
        # No command line argument, ask user which body to run
        bodyname = input("Please input body name: ")
        if bodyname == "":
            print("No body name entered. Defaulting to Europa.")
            bodyname = "Europa"

    bodyname = bodyname.capitalize()
    body = importlib.import_module(bodyname+".PP"+bodyname)
    Planet = body.Planet
    Params = body.Params

    outPlanet = PlanetProfile(Planet, Params)

    return
    
""" END MAIN RUN BLOCK """

def PlanetProfile(Planet, Params):

    # File name bases, to which we add specifics to create full file names
    savebase = Planet.name + '/' + Planet.name + 'Profile_'
    figbase = Planet.name + '/figures/' + Planet.name
    # Attach extra identifiers in special cases, to distinguish from standard cases
    if hasattr(Planet,'clathrate'): savebase += 'Clathrates_'

    # Read in data from Matlab dump for testing purposes
    mantleSizePath = savebase + Params.lbls.vmant + str(Planet.Tb_K) + '.csv'
    #mantleSizeR, mantleSizeRho = np.loadtxt(mantleSizePath, skiprows=1, unpack=True, delimiter=",")
    #nMantleInds = len(MantleSizeR)
    #mantleSizePlot(Planet, Params, mantleSizeRho[:nMantleInds], mantleSizeR[:nMantleInds], figbase)

    outPlanet = Planet
    return outPlanet

def writeProfile(path,saveStr,header,data):
    with open(path+saveStr+".txt","w") as f:
        f.write(header+"\n")
        for line in data:
            f.write("\t".join( [ "{:3.5e}".format(val) for val in line] )+"\n")

    return

if __name__ == '__main__': main()