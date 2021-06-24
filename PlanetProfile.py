import numpy as np
#from Thermodynamics.FromLiterature.conductiveMantleTemperature import conductiveMantleTemperature
#from Thermodynamics.FromLiterature.ConvectionDeschampsSotin2001 import ConvectionDeschampsSotin2001
from MantlePlot import MantlePlot

def writeProfile(path,saveStr,header,data):
    with open(path+saveStr+".txt","w") as f:
        f.write(header+"\n")
        for line in data:
            f.write("\t".join( [ "{:3.5e}".format(val) for val in line] )+"\n")

def celltolist(inputCellPy):
    """
        Helper function for converting MATLAB nested cell to python 2d list.
    """

    out = []
    for row in inputCellPy:
        out.append([])
        for val in row:
            val = int(val)
            out[-1].append(val)

    return out