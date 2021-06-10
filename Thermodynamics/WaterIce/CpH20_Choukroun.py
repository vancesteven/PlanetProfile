def CpH2O_Choukroun(P,T,ind):
    """
        Parameters:
        -----------
        P : float
            pressure [MPa] not currently used
        T : float
            temperature [K]
        ind : int
            index of phase

        Returns:
        --------
        Cp : float
            heat capacity [J/kh/K]
    """
    Cp = getCp2010(T,ind)

    return Cp

def getCp2010(T,ind):
    c0 = [4190, 74.11, 2200, 820, 700, 940, 2150]
    c1 = [9, 7.56, 0, 7, 7.56, 5.5, 3.19]

    return 