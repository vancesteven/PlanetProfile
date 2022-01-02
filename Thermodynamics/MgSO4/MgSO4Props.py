import numpy as np

def MgSO4Props(P_MPa, T_K, wOcean_ppt):
    """ Determine density rho, heat capacity Cp, thermal expansivity alpha,
        and thermal conductivity kTherm as functions of pressure P and
        temperature T for dissolved MgSO4 at a concentration wOcean in
        ppt by mass.

        Args:
            P_MPa (float, shape N): Pressures in MPa
            T_K (float, shape M): Temperature in K
            wOcean_ppt (float): (Absolute) salinity of Seawater in ppt by mass (g/kg)
        Returns:
            rho_kgm3 (float, shape NxM): Mass density of liquid in kg/m^3
            Cp_JkgK (float, shape NxM): Isobaric heat capacity of liquid in J/(kg K)
            alpha_pK (float, shape NxM): Thermal expansivity of liquid in 1/K
            kTherm_WmK (float, shape NxM): Thermal conductivity of liquid in W/(m K)
    """
    nPs = np.size(P_MPa)
    nTs = np.size(T_K)

    return rho_kgm3, Cp_JkgK, alpha_pK, kTherm_WmK


def GetPhaseMgSO4(P_MPa, T_K, wOcean_ppt):
    """ Calculate phase of liquid/ice within the hydrosphere for an ocean with
        dissolved MgSO4, given a span of P_MPa, T_K, and w_ppt, based on models
        from Vance et al. 2014: https://doi.org/10.1016/j.pss.2014.03.011
    """
    fn_phase = 0

    return fn_phase


def Molal2ppt(b_molkg, m_gmol):
    """ Convert dissolved salt concentration from molality to ppt

        Args:
            b_molkg (float, shape N): Concentration(s) in mol/kg (molal)
            m_gmol (float, shape 1 or shape N): Molecular weight in g/mol
        Returns:
            w_ppt (float, shape N): Corresponding concentration(s) in ppt by mass
    """
    m_kgmol = m_gmol / 1e3
    w_ppt = b_molkg*m_kgmol / (1 + b_molkg*m_kgmol)

    return w_ppt
