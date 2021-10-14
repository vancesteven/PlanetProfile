
# Work in progress python file from June 2021 - Artyom
def Helgerud_sI(P,T):
    """
        Determines Vp,Vs,shear,K and rho

        Parameters
        ----------
        P : float
            pressure [MPa]
        T : float
            temperature [K]
        
        Returns
        -------
        out : dictionary with keys
            "Vp" : float
            "Vs" : float
            "shear" : float
            "K" : float
            "rho" : float
            "alpha" : float

        Notes
        -----
        From Helgerud 2009
        A correction was issued that changed formulas from a(T+b)(P+c) to aT+bP+c

        It is unclear where the formula for 'rho' was acquired
    """

    out = {}

    out["alpha"] = ( 2*3.5038*10**(-4)*T+0.2517 ) / (3.5038*10**(-4)*(T**2)+0.2517*T+1610.9373)
    
    Tc = T - 273.15 # temperature in celcius
    out["Vp"] = -1.84*Tc + 0.31*P + 3766
    out["Vs"] = -0.892*Tc - 0.1*P + 1957
    out["shear"] = -4.2e-3*Tc + 9e-5*P + 3.541
    out["K"] = -1.9e-2*Tc + 3.8e-3*P + 8.39
    out["rho"] = (-2.3815e-4*Tc + 1.1843e-4*P + 0.91801)*1e3

    return out

