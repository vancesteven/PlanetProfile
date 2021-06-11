def Tm_p_Hirschmann2000(P):
    """
        Calculates temperature of peridotite solidus given pressure

        Parameters:
        -----------
        P : float
            pressure [GPa]

        Returns:
        --------
        T : float
            temperature [K]

        Notes:
        ------
        Uses recommended fit from Hirschmann 2000.
    """
    T = -5.140*P**2 + 132.899*P + 1120.661
    return T