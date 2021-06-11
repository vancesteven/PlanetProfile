import warnings

def getK_Andersson2005(P,T,varstr,typestr):
    """
        Uses pressure or temperature to retrieve thermal conductivity of ice phase from table.

        Parameters
        ----------
        P : float
            pressure [MPa]
        T : float
            temperature [K]
        varstr : string
            ice phase label (e.g. Ih,VII,LDA)
        typestr : string
            "T" or "P" to choose which table to use

        Returns
        -------
        K : float
            thermal conductivity [W/m/K]
    """

    if typestr == "T":
        K = kT(T,varstr)
    elif typestr == "P":
        K = kP(P,varstr)
    
    return K

def kT(T,varstr):
    """
        Helper function for getK_Andersson2005 in case 'T'

        Parameters
        ----------
        T : float
            temperature [K]

        Returns
        -------
        K : float
            thermal conductivity [W/m/K]
    """

    valueDict = {}

    # see https://en.wikipedia.org/wiki/Ice#Phases for some more details

    valueDict['Ih'] = [630, 0.995, 0.1, 40, 180, 5]
    valueDict['Ic'] = [192, 0.813, 0.1, 150, 200, 6]
    valueDict['II'] = [695, 1.097, 240, 120, 240, 4]
    valueDict['III'] = [93.2, 0.822, 240, 180, 250, 4]
    # IV metastable, so not useful for relevant calculations
    valueDict['V'] = [38.0, 0.612, 530, 240, 270, 4]
    valueDict['VI'] = [50.9, 0.612, 1000, 135, 250, 4]
    valueDict['VII'] = [320, 0.821, 2400, 275, 300, 4]
    valueDict['VIII'] = [15, 700, 1.417, 2400, 240, 270, 4]
    valueDict['IX'] = [164, 0.879, 240, 120, 160, 4]
    # X forms at about 70GPa, so not useful for relevant calculations
    valueDict['XI'] = [994, 1.041, 0.1, 50, 72, 6]
    valueDict['XII'] = [80.1, 0.77, 100, 115, 140, 4] # also metastable

    # amorphous phases of ice
    valueDict['LDA']= [17.5, 0.55, 0.1, 80, 130, 4] # low density amorphous
    valueDict['HDA']= [0.528, 0.026, 0.1, 75, 120, 4] # high density amorphous
    valueDict['VHDA']= [0.257, 0.20, 1000, 145, 155, 6] # very high density amorphous

    if varstr in valueDict:
        coeffs = valueDict[varstr]
    else:
        warnings.warn( "Ice "+varstr+" does not have relevant data" )
        coeffs = [0, 0]

    Dwm = coeffs[0]
    xT = coeffs[1]
    K = Dwm * T**(-xT)

    return K

def kP(P,varstr):

    """
        Helper function for getK_Andersson2005 in case 'P'

        Parameters
        ----------
        P : float
            pressure [MPa]

        Returns
        -------
        K : float
            thermal conductivity [W/m/K]
    """

    valueDict = {}

    # see https://en.wikipedia.org/wiki/Ice#Phases for some more details

    valueDict['Ih'] = [1.60, -0.44, 0, 0.50, 130]
    valueDict['Ic'] = [1.28, -0.34, 0, 0.50, 130]
    valueDict['II'] = [1.25, 0.2, 0, 0.24, 120]
    valueDict['III'] = [-0.02, 0.2, 0.2, 0.35, 240]
    # IV metastable, so not useful for relevant calculations
    valueDict['V'] = [0.16, 0.2, 0.35, 0.60, 246]
    valueDict['VI'] = [0.37, 0.16, 0.7, 2.0, 246]
    valueDict['VII'] = [0.65, 0.2, 2.0, 2.4, 286]
    valueDict['VIII'] = [1.38, 0.2, 2.0, 2.4, 246]
    # IX ? ? ? ?
    # X forms at about 70GPa, so not useful for relevant calculations
    valueDict['XI'] = [2.67, 0.24, 0.08, 0.16, 58]
    valueDict['XII'] = [0.73, 0.09, 0.1, 0.8, 115]

    # amorphous phases of ice
    valueDict['LDA']= [0.17, -0.63, 0, 0.35, 130]
    valueDict['HDA130K18']= [-0.52, 0.19, 0, 0.50, 130]
    # VHDA ? ? ? ?

    if varstr in valueDict: # if desired ice is one we have values for
        coeffs = valueDict[varstr] # use those values
    else:
        warnings.warn( "Ice "+varstr+" does not have relevant data" )
        coeffs = [0, 0] # otherwise, return with no a

    Ek = coeffs[0]
    F = coeffs[1] # [GPam1]

    K = Ek + F*P*1e-3

    return K