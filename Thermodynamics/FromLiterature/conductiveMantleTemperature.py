import numpy as np

def conductiveMantleTemperature(r,rb,a,kr,rho,T0,Q,qb):
    """
        Finds temperatures of the mantle of planet at provided radii.

        Parameters
        ----------
        r : numpy array with shape (n,)
            List of radii [m]
        rb : float
            Radius at bottom of conducting mantle [m]
        a : float
            Radius at top of conducting mantle [m]
        kr : float
            Thermal conductivity of mantle [W/m]
        rho : float
            Density of mantle [kg/m^3]
        T0 : float
            Temperature at top of mantle [K]
        Q : float
            Internal tidal heat [J ??]
        qb : float
            Heat flow supplied to base [J ??]

        Returns
        -------
        T : float array-like with shape (n,)
            List of temperatures at given radii

        Notes
        -----
        As described in Cammarano et al. 2006 (using a result from Turcotte,Schubert 1982).
        In the paper, the formula is used for conduction of heat through an ice layer.
    """
    m = rho*4/3*np.pi*(a**3-rb**3) # mass of the mantle [kg]
    H = Q/m # internal tidal heating [J/kg]
    term1 = (rho * H) / (6 * kr) * ( a**2 - np.power(r,2) )
    term2 = ( (rho * H * rb**3) / (3*kr) - (qb * rb**2) / kr ) * (np.power(r,-1) - 1/a)
    # In Cammarano et al., term2 (erroneously) uses (1/a - 1/r) instead.
    # Switching the sign (as done here) makes the term positive, as desired.
    T = T0 + term1 + term2
    return T