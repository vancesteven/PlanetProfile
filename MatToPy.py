import numpy as np
from scipy import interpolate

def cellto2dlist(inputCell,datatype):
    """
        Helper function for converting MATLAB nested cell to python 2d list.

        Parameters
        ----------
        inputCell : a 1D cell array in MATLAB containing other 1D cell arrays
        datatype : string
            either "int" or "float" depending on application

        Returns
        -------
        out : a nested numpy array
    """

    nRows = len(inputCell)

    out = np.empty(shape=(nRows,),dtype=object)

    for i in range(nRows):
        rowLen = len(inputCell[i])
        if datatype=="int":
            out[i] = np.empty(shape=(rowLen,),dtype=int)
        if datatype=="float":
            out[i] = np.empty(shape=(rowLen,),dtype=float)
        for j in range(rowLen):
            out[i][j] = inputCell[i][j]

    return out

"""
    work in progress conversion things
def MATtoPyGriddedInterp(gridVectors, values):
    
        Helper function for converting MATLAB griddedInterpolant to python scipy.interpolate.RegularGridInterpolator

        Parameters
        ----------
        gridVectors : 1x2 cell in MATLAB, each cell containing a 1D list of floats
        values : double array in MATLAB

        Returns
        -------
        out : scipy.interpolate.RegularGridInterpolator
    

    out = interpolate.RegularGridInterpolator( (gridVectors[0],gridVectors[1]), values )

    return out

def PlanetStruct(Planet_MAT):
    print(Planet_MAT["Tb_K"])
"""