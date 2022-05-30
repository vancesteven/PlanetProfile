import numpy as np
import logging
from gsw.freezing import t_freezing as gswTfreeze
from gsw.conversions import CT_from_t, C_from_SP as gswConduct_mScm
from gsw.density import rho, alpha as alphaCT, sound_speed as gswVP_ms
from gsw.energy import enthalpy
from gsw._wrapped_ufuncs import alpha_wrt_t_exact as alpha
from PlanetProfile.Thermodynamics.InnerEOS import ResetNearestExtrap
from PlanetProfile.Utilities.defineStructs import Constants

# Assign logger
log = logging.getLogger('PlanetProfile')

def SwProps(P_MPa, T_K, wOcean_ppt):
    """ Determine density rho, heat capacity Cp, thermal expansivity alpha,
        and thermal conductivity kTherm as functions of pressure P and
        temperature T for a Seawater composition at a concentration wOcean in
        ppt by mass. Uses the Gibbs Seawater (GSW) implementation of the
        Thermodynamic Equation of Seawater (TEOS-10).

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
    T_C = T_K - Constants.T0
    SP_dbar = MPa2seaPressure(P_MPa)
    CT_C = gswT2conservT(wOcean_ppt, T_C, SP_dbar)
    rho_kgm3 = gswDensity_kgm3(wOcean_ppt, CT_C, SP_dbar)

    dCTdT, Cp_JkgK = GetdCTdTanddHdT(wOcean_ppt, CT_C, SP_dbar, T_C)
    GSW_ALPHA_WRT_CT = False
    if GSW_ALPHA_WRT_CT:
        # Use only functions appearing in the main listing here
        log.warning('The Python GSW package does not yet contain a calculation for alpha with respect to ' +
                    'in-situ temperature in standard libraries, which we use as T_K. alpha_pK will be scaled ' +
                    'approximately by evaluating d(CT_C)/d(T_K) numerically and multiplying this by the result from GSW.')
        alpha_pK = gswExpansivityCT_pK(wOcean_ppt, CT_C, SP_dbar) * dCTdT
    else:
        alpha_pK = gswExpansivity_pK(wOcean_ppt, T_C, SP_dbar)
    kTherm_WmK = np.zeros_like(alpha_pK) + Constants.kThermWater_WmK  # Placeholder until we implement a self-consistent calculation

    return rho_kgm3, Cp_JkgK, alpha_pK, kTherm_WmK


def gswT2conservT(wOcean_ppt, T_C, SP_dbar, DO_1D=False):
    """ Wrapper for GSW function CT_from_t that can simultaneously handle
        P and T arrays.

        Optional args:
            DO_1D (bool): If True, treat T_C and SP_dbar as same-length arrays that
                represent the values of a depth profile and return a 1D array of CT_C
                values.
    """
    if DO_1D:
        CT_C = np.array([CT_from_t(wOcean_ppt, T_C[i], SP_dbar[i]) for i in range(np.size(SP_dbar))])
    else:
        CT_C = np.array([CT_from_t(wOcean_ppt, T_C, SPi_dbar) for SPi_dbar in SP_dbar])
    return CT_C


def gswEnthalpy_Jkg(wOcean_ppt, CT_C, SP_dbar):
    """ Wrapper for GSW function enthalpy that can simultaneously handle
        P and T arrays.
    """
    H_Jkg = np.array([enthalpy(wOcean_ppt, CT_C[i,:], SP_dbar[i]) for i in range(np.size(SP_dbar))])
    return H_Jkg


def gswDensity_kgm3(wOcean_ppt, CT_C, SP_dbar, DO_1D=False):
    """ Wrapper for GSW function rho that can simultaneously handle
        P and T arrays.

        Optional args:
            DO_1D (bool): If True, treat CT_C and SP_dbar as same-length arrays that
                represent the values of a depth profile and return a 1D array of rho
                values.
    """
    if DO_1D:
        rho_kgm3 = np.array([rho(wOcean_ppt, CT_C[i], SP_dbar[i]) for i in range(np.size(SP_dbar))])
    else:
        rho_kgm3 = np.array([rho(wOcean_ppt, CT_C[i,:], SP_dbar[i]) for i in range(np.size(SP_dbar))])
    return rho_kgm3


def gswExpansivityCT_pK(wOcean_ppt, CT_C, SP_dbar):
    """ Wrapper for GSW function alpha that can simultaneously handle
        P and T arrays.
    """
    alphawrtCT_pK = np.array([alpha(wOcean_ppt, CT_C[i,:], SP_dbar[i]) for i in range(np.size(SP_dbar))])
    return alphawrtCT_pK


def gswExpansivity_pK(wOcean_ppt, T_C, SP_dbar):
    """ Wrapper for GSW function alpha_wrt_t_exact that can simultaneously handle
        P and T arrays.
    """
    alphawrt_pK = np.array([alpha(wOcean_ppt, T_C, SPi_dbar) for SPi_dbar in SP_dbar])
    return alphawrt_pK


def GetdCTdTanddHdT(wOcean_ppt, CT_C, SP_dbar, T_C, dT=0.01):
    """ Evaluate the first partial derivative of conservative temperature CT in celsius
        with respect to in-situ temperature T in C.

        Args:
            wOcean_ppt (float): (Absolute) salinity of Seawater in ppt by mass (g/kg)
            CT_C (float, shape NxM): Conservative temperature values corresponding
                to input (P,T) arrays in celsius
            SP_dbar (float, shape N): Sea pressures in dbar
            T_C (float, shape M): In-situ temperatures in C
            dT (optional, float): Half the small change in temperature to use to get dCT/dT
                (as we are not dividing this quantity by 2 in operations, it is dT/2)
        Returns:
            dCTdT (float, shape NxM): First partial derivative of CT with respect to T
            dHdT (float, shape NxM): First partial derivative of specific enthalpy H with
                respect to T (equal to Cp_JkgK)
    """
    # Get CT a bit above and a bit below each of the evaluated points
    CTplus_C = gswT2conservT(wOcean_ppt, T_C+dT, SP_dbar)
    CTless_C = gswT2conservT(wOcean_ppt, T_C-dT, SP_dbar)
    # Get dCT, the difference in CT with respect to this small change in T dT
    dCTplus = CTplus_C - CT_C
    dCTless = CT_C - CTless_C
    # Take the average value between the numerically evaluated derivatives above and below
    dCTdT = np.mean([dCTplus/dT, dCTless/dT], axis=0)

    # Use the 3 sets of CT values to also evaluate H(T+/-dT)
    Hmid_Jkg = gswEnthalpy_Jkg(wOcean_ppt, CT_C, SP_dbar)
    Hplus_Jkg = gswEnthalpy_Jkg(wOcean_ppt, CTplus_C, SP_dbar)
    Hless_Jkg = gswEnthalpy_Jkg(wOcean_ppt, CTless_C, SP_dbar)
    # Get differences in H
    dHplus = Hplus_Jkg - Hmid_Jkg
    dHless = Hmid_Jkg - Hless_Jkg
    # Take the average value between the numerically evaluated derivatives above and below
    dHdT = np.mean([dHplus/dT, dHless/dT], axis=0)

    return dCTdT, dHdT


class SwPhase:
    def __init__(self, wOcean_ppt):
        self.w_ppt = wOcean_ppt

    def __call__(self, P_MPa, T_K):
        P_MPa = np.array(P_MPa)
        T_K = np.array(T_K)
        if(np.size(P_MPa)==0 or np.size(T_K)==0):
            # If input is empty, return empty array
            return np.array([])
        elif((np.size(P_MPa) != np.size(T_K)) and not (np.size(P_MPa)==1 or np.size(T_K)==1)):
            # If arrays are different lengths, they are probably meant to get a 2D output
            P_MPa, T_K = np.meshgrid(P_MPa, T_K, indexing='ij')

        # 1. Convert to "sea pressure" and T in celsius as needed for GSW input
        # 2. Subtract the freezing temperature from the input temperature
        # 3. Compare to zero -- if we are below the freezing temp, it's ice I, above, liquid
        # 4. Cast the above comparison (True if less than Tfreeze, False if greater) to int,
        #       so that we get 1 if we are below the freezing temp and 0 if above.
        return (((T_K - Constants.T0) - gswTfreeze(self.w_ppt, MPa2seaPressure(P_MPa), 0)) < 0).astype(np.int_)


class SwSeismic:
    def __init__(self, wOcean_ppt, EXTRAP):
        self.w_ppt = wOcean_ppt
        self.EXTRAP = EXTRAP

    def __call__(self, P_MPa, T_K):
        if not self.EXTRAP:
            # Set extrapolation boundaries to limits inferred from GSW
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, 0, 200, 250, 365)
        T_C = T_K - Constants.T0
        SP_dbar = MPa2seaPressure(P_MPa)
        CT_C = gswT2conservT(self.w_ppt, T_C, SP_dbar, DO_1D=True)
        VP_kms = gswVP_ms(self.w_ppt, CT_C, SP_dbar) * 1e-3  # 1e-3 to convert from m/s to km/s
        KS_GPa = gswDensity_kgm3(self.w_ppt, CT_C, SP_dbar, DO_1D=True) * VP_kms**2 * 1e-3  # 1e-3 because (km/s)^2 * (kg/m^3) gives units of MPa, so 1e-3 to convert to GPa
        return VP_kms, KS_GPa


class SwConduct:
    def __init__(self, wOcean_ppt):
        self.w_ppt = wOcean_ppt
        # For a Seawater composition, practical salinity SP on the PSS-78 scale
        # is directly proportional to absolute salinity in g/kg--see Eq. 1 of
        # https://www.teos-10.org/pubs/gsw/pdf/SAAR.pdf
        self.PracSalin = self.w_ppt * (35/Constants.stdSeawater_ppt)

    def __call__(self, P_MPa, T_K, grid=False):
        T_C = T_K - Constants.T0
        SP_dbar = MPa2seaPressure(P_MPa)
        # GSW C_from_SP function (gswConduct) returns mS/cm, so we
        # multiply by 0.1 to get S/m
        return gswConduct_mScm(self.PracSalin, T_C, SP_dbar) * 0.1


def MPa2seaPressure(P_MPa):
    """ Calculates "sea pressure" as needed for inputs to GSW functions.
        Sea pressure is defined as the pressure relative to the top of
        the ocean, equal to the absolute pressure less 10.1325 dbar.

        Args:
            P_MPa (float, shape N or NxM): Pressures in MPa
        Returns:
            SP_dbar (float, shape N or NxM): Sea pressures in dbar (1 dbar = 0.1 bar = 0.1 * bar2MPa MPa)
    """
    # Subtract off Earth atmospheric pressure
    SP_MPa = P_MPa - Constants.bar2MPa
    # Convert to dbar
    SP_dbar = SP_MPa * 10 / Constants.bar2MPa

    return SP_dbar
