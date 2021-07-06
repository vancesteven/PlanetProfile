import scipy.interpolate as interpolate
import numpy as np
import os.path
import sys

def fnTfreeze_K():
    with open(os.path.join(sys.path[0], "L_Ice_MgSO4.txt"), "r") as f:
        f.readline() # Skip header info line
        PPg_line = f.readline()
        wwg_line = f.readline()
        Pmin, Pmax, Psteps = map(int, PPg_line.split(",")[1:])
        wtPctMin, wtPctMax, wtPctsteps = map(int, wwg_line.split(",")[1:])
    PPg = np.linspace(Pmin,Pmax,Psteps)
    wwg = np.linspace(wtPctMin,wtPctMax,wtPctsteps)
    TT = np.loadtxt( os.path.join(sys.path[0], "L_Ice_MgSO4.txt") , skiprows=3 , delimiter=',')
    output = interpolate.RegularGridInterpolator( (PPg,wwg) , TT )
    return output