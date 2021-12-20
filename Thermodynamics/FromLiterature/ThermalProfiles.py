import numpy as np
from Utilities.dataStructs import Constants
from seafreeze import seafreeze as SeaFreeze
from Thermodynamics.HydroEOS import PhaseConv, kThermMelinder2007, kThermHobbs1974, \
    kThermIsobaricAnderssonIbari2005, kThermIsothermalAnderssonIbari2005

def ConductionClathLid():
    """ Thermodynamics calculations for a (thermally) conductive clathrate lid

        Args:
            ???
        Returns:
            ???
    """

    return


def ConvectionDeschampsSotin2001(Ttop_K, rTop_m, kTop_WmK, Tb_K, zb_m, gtop_ms2, Pmid_MPa, phase, EQUIL_Q):
    """ Thermodynamics calculations for convection in an ice layer
        based on Deschamps and Sotin (2001): https://doi.org/10.1029/2000JE001253
        Note that these authors solved for the scaling laws we apply in Cartesian
        geometry and for the case of zero internal heating within the ice shell.
        These are poor assumptions for thick ice shells and where significant
        tidal heating occurs in the ice shell. Strictly speaking, we use values
        for parameters such as the "mid" pressure that are not consistent with
        the equations of Deschamps and Sotin (2001), and several of their equations
        are poorly described, contain errors (e.g. unit problem in Eq. 2), or use
        physically inconsistent values, as in defining the "core" Rayleigh number Ra.
        A more robust model for convection in the ice shell should be found and
        used to replace this.

        Args:
            Ttop_K (float): Temperature at top of whole layer in K
            rTop_m (float): Radius of top of ice layer in m
            kTop_WmK (float): Thermal conductivity at top of ice layer in W/(mK)
            Tb_K (float): Assumed bottom/melting temperature in K
            zb_m (float): Thickness of the ice layer in m
            gtop_ms2 (float): Gravitational acceleration at layer top
            Pmid_MPa (float): Pressure at the "middle" of the convective region in MPa
            phase (int): Ice phase index
            EQUIL_Q (bool): Whether to set heat flux from interior to be consistent with that released
                through the convective profile according to Deschamps and Sotin (2001) (if True) or to
                set it to a value consistent with Ojakangas and Stevenson (1989) for the upper conductive
                lid portion (if False)
        Returns:
            Tconv_K (float): Temperature of "well-mixed", convective region in K
            etaConv_Pas (float): Viscosity of "well-mixed", convective region in Pa*s
            eLid_m (float): Thickness of the stagnant lid conductive layer in m
            deltaTBL_m (float): Thickness of the thermal boundary layer (TBL) between the ocean
                and convective region in m
            qbot_Wm2 (float): Heat flux at the bottom of the ice in W/m^2
    """
    # Numerical constants derived in Cartesian geometry from Deschamps and Sotin (2000) and used in
    # Deschamps and Sotin (2001) parameterization
    c1 = 1.43
    c2 = -0.03
    # Numerical constants appearing in equations
    A = Constants.Eact_kJmol[phase] * 1e3 / Constants.R / Tb_K
    B = Constants.Eact_kJmol[phase] * 1e3 / 2 / Constants.R / c1
    C = c2 * (Tb_K - Ttop_K)
    # Temperature and viscosity of the "well-mixed" convective region
    Tconv_K = B * (np.sqrt(1 + 2/B*(Tb_K - C)) - 1)
    etaConv_Pas = Constants.etaMelt_Pas[phase] * np.exp(A * (Tb_K/Tconv_K - 1))
    # Get physical properties of ice at the "middle" of the convective region
    seaOut = SeaFreeze(np.array([(Pmid_MPa, Tconv_K)], dtype='f,f').astype(object), PhaseConv(phase))
    alphaMid_pK = seaOut.alpha[0]
    rhoMid_kgm3 = seaOut.rho[0]
    CpMid_JkgK = seaOut.Cp[0]
    kMid_WmK = kThermIsobaricAnderssonIbari2005(Tconv_K, phase)
    # Rayleigh number of whole ice layer, derived using viscosity of convective region
    Ra = alphaMid_pK * CpMid_JkgK * rhoMid_kgm3**2 * gtop_ms2 * (Tb_K - Ttop_K) * zb_m**3 / etaConv_Pas / kMid_WmK
    # Rayleigh number of lower thermal boundary layer, from parameterization results of Deschamps and Sotin (2000)
    Radelta = 0.28 * Ra**0.21
    # Thickness of lower thermal boundary layer
    deltaTBL_m = (etaConv_Pas * kMid_WmK * Radelta /
                  alphaMid_pK / CpMid_JkgK / rhoMid_kgm3**2 / gtop_ms2 / (Tb_K - Tconv_K))**(1/3)
    # Heat flux entering the bottom of the ice layer
    qbot_Wm2 = kMid_WmK * (Tb_K - Tconv_K) / deltaTBL_m
    # Heat flux leaving the top of the ice layer (adjusted for spherical geometry compared to Deschamps and Sotin, 2001)
    qtop_Wm2 = (rTop_m - zb_m)**2 / rTop_m**2 * qbot_Wm2
    # Thickness of conductive stagnant lid
    # This matches the Matlab, but based on Deschamps and Sotin (2001), who assume a fixed thermal conductivity
    # throughout the ice shell, the thickness of the conductive lid should probably use kTop_WmK, because that
    # will be what determines the conductive thermal profile based on the heat flux through the lid.
    eLid_m = kMid_WmK * (Tconv_K - Ttop_K) / qtop_Wm2
    #eLid_m = kTop_WmK * (Tconv_K - Ttop_K) / qtop_Wm2

    # If the Rayleigh number is less than some critical value, convection does not occur.
    if(Ra < Constants.RaCrit):
        print('Rayleigh number of ' + str(round(Ra)) + ' in the ice shell is less than the critical value ' +
              'of ' + str(round(Constants.RaCrit)) + '. Only conduction will be modeled in the ice shell.')
        # Set conductive layer thicknesses to whole shell thickness to force a whole-layer conductive profile
        eLid_m = zb_m
        deltaTBL_m = zb_m
        Tconv_K = Tb_K

    if (Ra < Constants.RaCrit or not EQUIL_Q):
        # Set heat flux to be equal to that escaping the conductive lid
        # according to Ojakangas and Stevenson (1989): https://doi.org/10.1016/0019-1035(89)90052-3
        Qbot_W = kMid_WmK * Ttop_K / eLid_m * np.log(Tconv_K/Ttop_K) * 4*np.pi * (rTop_m - eLid_m)**2

        # The below matches the Matlab, but this is not what Deschamps and Sotin (2001) do
        # and it is not consistent with the results of Andersson and Inaba (2005).
        # Dcond = np.array([np.nan, 632, 418, 242, np.nan, 328, 183])
        # qbot_Wm2 = Dcond[phase] * np.log(Tb_K/Ttop_K) / zb_m
    else:
        # Convection is occurring
        # Total heat flow entering the bottom of the ice shell
        Qbot_W = qbot_Wm2 * 4*np.pi * (rTop_m - zb_m)**2

    return Tconv_K, etaConv_Pas, eLid_m, deltaTBL_m, Qbot_W, Ra


def ConductiveTemperature(Ttop_K, rTop_m, rBot_m, kTherm_WmK, rho_kgm3, Qrad_Wkg, Htidal_Wm3, qTop_Wm2):
    """ Thermal profile for purely thermally conductive layers, based on Turcotte and Schubert (1982),
        equation 4.40: T = -rho*H/6/k * r^2 + c1/r + c2, where c1 and c2 are integration constants
        found through boundary conditions, rho is mass density of the conductive layer in kg/m^3,
        H is internal heating in W/kg, and k is thermal conductivity in W/m/K.
        The main equations we use here are developed similar to equation 2 of
        Cammarano et al. (2006): https://doi.org/10.1029/2006JE002710, but we parameterize in terms of
        the heat flux leaving the top of the layer instead of entering the bottom.

        Args:
            Ttop_K (float, shape N): Temperature at the top of the layer in K.
            rTop_m, rBot_m (float, shape N): Radius at top and bottom of layer in m, respectively.
            kTherm_WmK (float, shape N): Thermal conductivity of layer in W/(mK).
            rho_kgm3 (float, shape N): Mass density of layer in kg/m^3.
            Qrad_Wkg (float): Average radiogenic heating rate in W/kg.
            Htidal_Wm3 (float): Average tidal heating rate in W/m^3.
            qTop_Wm2 (float): Heat flux leaving the top of the layer in W/m^2.
        Returns:
            Tbot_K (float): Temperature at the bottom of the layer in K.
            qBot_Wm2 (float): Heat flux entering the bottom of the layer in W/m^2.
    """
    # Calculate needed values from inputs
    Qtot_Wkg = Qrad_Wkg + Htidal_Wm3 / rho_kgm3
    c1 = qTop_Wm2 * rTop_m**2 / 2/kTherm_WmK - rho_kgm3 * Qtot_Wkg / 6/kTherm_WmK * rTop_m**3
    # Find the temperature at the bottom of the layer
    Tbot_K = Ttop_K + rho_kgm3*Qtot_Wkg / 6/kTherm_WmK * (rTop_m**2 - rBot_m**2) + c1 * (1/rBot_m - 1/rTop_m)
    # Find the heat flux into the bottom of the layer
    qBot_Wm2 = rho_kgm3 * Qtot_Wkg / 3 * rBot_m + 2*kTherm_WmK / rBot_m**2 * c1

    return Tbot_K, qBot_Wm2

