import numpy as np
from Utilities.dataStructs import Constants

def ConductionClathLid():
    """ Thermodynamics calculations for a (thermally) conductive clathrate lid

        Args:
            ???
        Returns:
            ???
    """

    return


def ConvectionDeschampsSotin2001():
    """ Thermodynamics calculations for convection in the near-surface ice Ih layer

        Args:
            ???
        Returns:
            ???
    """

    return


def ConductiveTemperature(Ttop_K, rTop_m, rBot_m, kTherm_WmK, rho_kgm3, Qrad_Wkg, Htidal_Wm3, qb_Wm2=0):
    """ Thermal profile for purely thermally conductive layers, based on Turcotte and Schubert (1982),
        as described in Cammarano et al. (2006).

        Args:
            rTop_m, rBot_m (float, shape N): Radius at top and bottom of layer in m, respectively.
            kTherm_WmK (float): Thermal conductivity of layer in W/(mK).
            rho_kgm3 (float, shape N): Mass density of layer in kg/m^3.
            Ttop_K (float, shape N): Temperature at the top of the layer in K.
            Qrad_Wkg (float): Average radiogenic heating rate in W/kg.
            Htidal_Wm3 (float): Average tidal heating rate in W/m^3.
            qb_Wm2=0 (float): Heating rate from below in W/m^2.
        Returns:
            Tbot_K (float): Temperature at the bottom of the layer in K.
    """
    Qtot_Wkg = Qrad_Wkg + Htidal_Wm3 / rho_kgm3
    Tbot_K = Ttop_K + rho_kgm3*Qtot_Wkg / 6/kTherm_WmK * (rTop_m**2 - rBot_m**2) + \
             (rho_kgm3*Qtot_Wkg*rBot_m**3 / 3/kTherm_WmK - qb_Wm2*rBot_m**2 / kTherm_WmK) * (1/rTop_m - 1/rBot_m)
    return Tbot_K


def TsolidusHirschmann2000(P_MPa):
    """ Silicate melting temperature parameterization based on
        Hirschmann (2000): https://doi.org/10.1029/2000GC000070 .
    """
    P_GPa = P_MPa * 1e-3
    a = -5.104
    b = 132.899
    c = 1120.661
    T_K = a*P_GPa**2 + b*P_GPa + c + Constants.T0
    return T_K