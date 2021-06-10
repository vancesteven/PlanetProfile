import numpy as np
from Thermodynamics.FromLiterature.conductiveMantleTemperature import conductiveMantleTemperature
from Thermodynamics.FromLiterature.ConvectionDeschampsSotin2001 import ConvectionDeschampsSotin2001

def writeProfile(path,saveStr,header,data):
    with open(path+saveStr+".txt","w") as f:
        f.write(header+"\n")
        for line in data:
            f.write("\t".join( [ "{:3.5e}".format(val) for val in line] )+"\n")