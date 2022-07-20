import numpy as np
from scipy.optimize import root_scalar as GetZero
from collections.abc import Iterable
from PlanetProfile.Utilities.defineStructs import Constants

def ClathProps(Plin_MPa, Tlin_K):
    """ Evaluate methane clathrate physical properties using Helgerud et al. (2009): https://doi.org/10.1029/2009JB006451
        for density, Ning et al. (2015): https://doi.org/10.1039/C4CP04212C for thermal expansivity and heat capacity,
        and Waite et al. (2005): https://www.researchgate.net/profile/W-Waite/publication/252708287_Thermal_Property_Measurements_in_Tetrahydrofuran_THF_Hydrate_Between_-25_and_4deg_C_and_Their_Application_to_Methane_Hydrate/links/57b900ae08aedfe0ec94abd7/Thermal-Property-Measurements-in-Tetrahydrofuran-THF-Hydrate-Between-25-and-4deg-C-and-Their-Application-to-Methane-Hydrate.pdf
        for thermal conductivity.
        Range of validity:
            rho_kgm3: P from 30.5 to 97.7 MPa, T from -20 to 15 C
            Cp_JkgK: P at 20 MPa, T from 5 to 292 K
            alpha_pK: P at 0.1 MPa, T from 5 to 268 K (note that Ning et al. also give a parameterization for
                alpha_pK at 20 MPa that differs slightly. Since clathrates only appear at the surface of icy moons,
                we use just the 1 bar value for simplicity.
            kTherm_WmK: P from 13.8 to 24.8 MPa, T from -30 to 20 C


        Args:
            Plin_MPa (float, shape M): Pressures to evaluate in MPa
            Tlin_K (float, shape N): Temperatures to evaluate in K
        Returns:
            rho_kgm3 (float, shape MxN): Mass density in kg/m^3 for each P and T. rho_gcm3 = aT_C + bP_MPa + c
            Cp_JkgK (float, shape MxN): Heat capacity in J/(kg K) for each P and T. Cp_JkgK = aT_K + b
            alpha_pK (float, shape MxN): Thermal expansivity in 1/K for each P and T. alpha_pK = (2aT_K + b)/(aT_K^2 + bT_K + c)
            kTherm_WmK (float, shape MxN): Thermal conductivity in W/(m K) for each P and T. kTherm_WmK = c, a constant, over the specified range.
    """
    P_MPa, T_K = np.meshgrid(Plin_MPa, Tlin_K, indexing='ij')

    T_C = T_K - Constants.T0

    rho_kgm3 = (-2.3815e-4*T_C + 1.1843e-4*P_MPa + 0.92435) * 1e3
    Cp_JkgK = 3.19*T_K + 2150
    alpha_pK = (3.5697e-4*T_K + 0.2558)/(3.5697e-4*T_K**2 + 0.2558*T_K + 1612.8597)
    kTherm_WmK = np.zeros_like(P_MPa) + 0.5

    return rho_kgm3, Cp_JkgK, alpha_pK, kTherm_WmK


""" Dissociation temperatures for clathrates based on a degree-2 fit to
    the dissociation curve at low pressure and a logarithmic fit at higher pressures,
    for curves from Sloan (1998) as reported in Choukron et al. (2010):
    https://doi.org/10.1016/j.icarus.2009.08.011
"""
class ClathDissoc:
    def __init__(self, Tb_K, NAGASHIMA_CLATH_DISSOC=True, ALLOW_BROKEN_MODELS=False, DO_EXPLOREOGRAM=False):
        self.Tb_K = Tb_K
        self.NAGASHIMA = NAGASHIMA_CLATH_DISSOC
        self.ALLOW_BROKEN_MODELS = ALLOW_BROKEN_MODELS
        self.DO_EXPLOREOGRAM = DO_EXPLOREOGRAM
        if self.NAGASHIMA:
            self.LOWER = None
            self.Pends_MPa = None
        else:
            if self.Tb_K < 273:
                self.LOWER = True
                self.Pends_MPa = [0.0, 2.567]
            else:
                self.LOWER = False
                self.Pends_MPa = [2.567, Constants.PmaxLiquid_MPa]

    def Tdissoc_K(self, P_MPa):
        if isinstance(P_MPa, Iterable):
            P_MPa[P_MPa == 0] = Constants.Pmin_MPa
        else:
            if P_MPa == 0:
                P_MPa = Constants.Pmin_MPa

        if self.NAGASHIMA:
            Td_K = TclathDissocNagashima_K(P_MPa)
        else:
            if self.LOWER:
                # For use when P_MPa < 2.567 and T < 273
                Td_K = TclathDissocLower_K(P_MPa)
            else:
                # For use when P_MPa >= 2.567 and T >= 273
                Td_K = TclathDissocUpper_K(P_MPa)

        return Td_K

    def TbZero_K(self, P_MPa):
        return self.Tb_K - self.Tdissoc_K(P_MPa)

    def PbClath_MPa(self):
        if self.NAGASHIMA:
            PbClath_MPa = PclathDissocNagashima_MPa(self.Tb_K)
        else:
            try:
                PbClath_MPa = GetZero(self.TbZero_K, bracket=self.Pends_MPa).root
            except ValueError:
                msg = f'No Pb was found for clathrates with Tb_K = {self.Tb_K}.'
                if self.ALLOW_BROKEN_MODELS:
                    if self.DO_EXPLOREOGRAM:
                        log.info(msg)
                    else:
                        log.error(
                            msg + 'ALLOW_BROKEN_MODELS is True, so calculations will proceed, with many values set to nan.')
                    PbClath_MPa = np.nan
                else:
                    raise ValueError(msg)

        return PbClath_MPa

def TclathDissocLower_K(P_MPa):
    if isinstance(P_MPa, Iterable):
        P_MPa[P_MPa == 0] = Constants.Pmin_MPa
    else:
        if P_MPa == 0:
            P_MPa = Constants.Pmin_MPa
    return 212.33820985 + 43.37319252 * P_MPa - 7.83348412 * P_MPa ** 2
def TclathDissocUpper_K(P_MPa):
    if isinstance(P_MPa, Iterable):
        P_MPa[P_MPa == 0] = Constants.Pmin_MPa
    else:
        if P_MPa == 0:
            P_MPa = Constants.Pmin_MPa
    return -20.3058036 + 8.09637199 * np.log(P_MPa / 4.56717945e-16)
def TclathDissocNagashima_K(P_MPa):
    if isinstance(P_MPa, Iterable):
        P_MPa[P_MPa == 0] = Constants.Pmin_MPa
    else:
        if P_MPa == 0:
            P_MPa = Constants.Pmin_MPa
    return 2214.1 / (15.959 - np.log(1e3*P_MPa))
def PclathDissocNagashima_MPa(T_K):
    return 1e-3 * np.exp(15.959 - 2214.1 / T_K)

def ClathStableSloan1998(P_MPa, T_K):
    """ Returns a grid of boolean values to say whether fully occupied methane
        clathrates are stable at the given P,T conditions, based on the dissocation
        curves above.

        Args:
            P_MPa (float, shape N): Pressures in MPa
            T_K (float, shape M): Temperatures in K
        Returns:
            stable (int, shape NxM): Set to Constants.phaseClath (phase ID) if clathrates are stable
                at this P,T and 0 if they are not (if not stable, Ocean.EOS.fn_phase should be
                queried to determine the phase)
    """

    # Avoid log of 0 for surface pressure
    P_MPa[P_MPa==0] = Constants.Pmin_MPa
    # Get (P,T) pairs relevant for each portion of the dissociation curve
    Plow_MPa, Tlow_K = np.meshgrid(P_MPa[P_MPa < 2.567], T_K, indexing='ij')
    Pupp_MPa, Tupp_K = np.meshgrid(P_MPa[P_MPa >= 2.567], T_K, indexing='ij')
    # Get evaluation points for lower and upper dissociation curves
    TdissocLower_K = TclathDissocLower_K(Plow_MPa)
    TdissocUpper_K = TclathDissocUpper_K(Pupp_MPa)
    # Assign clathrate phase ID to (P,T) points below the dissociation curves and 0 otherwise
    stableLow = np.zeros_like(Tlow_K).astype(np.int_)
    stableUpp = np.zeros_like(Tupp_K).astype(np.int_)
    stableLow[Tlow_K < TdissocLower_K] = Constants.phaseClath
    stableUpp[Tupp_K < TdissocUpper_K] = Constants.phaseClath
    stable = np.concatenate((stableLow, stableUpp), axis=0)

    return stable


def ClathStableNagashima2017(P_MPa, T_K):
    """ Same as above but for extrapolation from Nagashima (2017) dissertation
        according to S. Nozaki (private communication with M. Styczinski)
    """
    P2D_MPa, T2D_K = np.meshgrid(P_MPa, T_K, indexing='ij')
    TclathDissoc_K = TclathDissocNagashima_K(P2D_MPa)
    stable = np.zeros_like(P2D_MPa).astype(np.int_)
    stable[T2D_K > TclathDissoc_K] = Constants.phaseClath

    return stable


class ClathSeismic:
    """ Calculate seismic velocities in clathrates based on Helgerud et al. (2009): https://doi.org/10.1029/2009JB006451
        Note that the original article had a correction for all of the tables' equations--
        the correction is linked in the above DOI.

        Args:
            P_MPa (float, shape N): Pressure values to evaluate in MPa.
            T_K (float, shape N): Temperature values to evaluate in K.
        Returns:
            VP_kms (float, shape N): P-wave seismic velocity in km/s.
            VS_kms (float, shape N): S-wave seismic velocity in km/s.
            KS_GPa (float, shape N): Bulk modulus in GPa.
            GS_GPa (float, shape N): Shear modulus in GPa.
    """
    def __call__(self, P_MPa, T_K):
        T_C = T_K - Constants.T0
        VP_kms = (-1.84*T_C + 0.31*P_MPa + 3766) * 1e-3
        VS_kms = (-0.892*T_C - 0.1*P_MPa + 1957) * 1e-3
        KS_GPa = -1.09e-2*T_C + 3.8e-3*P_MPa + 8.39
        GS_GPa = -4.2e-3*T_C + 9e-5*P_MPa + 3.541

        return VP_kms, VS_kms, KS_GPa, GS_GPa
