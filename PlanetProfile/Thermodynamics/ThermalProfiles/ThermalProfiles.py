import numpy as np
import logging
from PlanetProfile.Utilities.defineStructs import Constants
from seafreeze import seafreeze as SeaFreeze
from PlanetProfile.Thermodynamics.HydroEOS import PhaseConv, GetTfreeze, GetPfreeze, kThermMelinder2007, \
    kThermHobbs1974, kThermIsobaricAnderssonIbari2005, kThermIsothermalAnderssonIbari2005

# Assign logger
log = logging.getLogger('PlanetProfile')

def ConvectionDeschampsSotin2001(Ttop_K, rTop_m, kTop_WmK, Tb_K, zb_m, gtop_ms2, Pmid_MPa,
                                 oceanEOS, iceEOS, phaseBot, EQUIL_Q):
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
            kTop_WmK (float): Thermal conductivity at top of ice layer in W/(m K)
            Tb_K (float): Assumed bottom temperature in K
            zb_m (float): Thickness of the ice layer in m
            gtop_ms2 (float): Gravitational acceleration at layer top
            Pmid_MPa (float): Pressure at the "middle" of the convective region in MPa
            oceanEOS (OceanEOSStruct): Interpolator functions for evaluating the ocean EOS
            iceEOS (IceEOSStruct): Interpolator functions for evaluating the ice EOS
            phaseBot (int): Ice phase index at the bottom of the layer
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
    # Get phase of convecting region from passed iceEOS
    phaseMid = iceEOS.phaseID
    # Numerical constants derived in Cartesian geometry from Deschamps and Sotin (2000) and used in
    # Deschamps and Sotin (2001) parameterization
    c1 = 1.43
    c2 = -0.03
    # Numerical constants appearing in equations
    A = Constants.Eact_kJmol[phaseMid] * 1e3 / Constants.R / Tb_K
    B = Constants.Eact_kJmol[phaseMid] * 1e3 / 2 / Constants.R / c1
    C = c2 * (Tb_K - Ttop_K)
    # Temperature and viscosity of the "well-mixed" convective region
    Tconv_K = B * (np.sqrt(1 + 2/B*(Tb_K - C)) - 1)
    if(Tconv_K < Ttop_K):
        Tconv_K = Ttop_K
        log.debug(f'Convecting temperature for ice {PhaseConv(phaseBot)} is less than the ' +
                   'temperature at the top of the layer. Tconv has been set equal to Ttop and ' +
                   'no conductive lid will be modeled.')
    if phaseMid != oceanEOS.fn_phase(Pmid_MPa, Tconv_K):
        if abs(phaseMid) != Constants.phaseClath:
            iceType = f'ice {PhaseConv(phaseMid)}'
            suggestion = 'Try adjusting Tb_K values to achieve a possible configuration.'
        else:
            iceType = f'clathrate'
            suggestion = 'This likely means whole-shell convection of clathrates is inconsistent ' + \
                         'with the low pressure found to be consistent with Tb_K. Try *increasing* Tb_K. ' + \
                         'Higher pressures are needed to keep warmer clathrates stable, so a higher Tb_K ' + \
                         'leads to a thicker ice shell prediction, which thickens the upper thermal ' + \
                         'boundary layer as needed for the temperature to be below the instability threshold ' + \
                         'throughout the ice shell.'
        log.warning(f'Convecting temperature of {iceType} exceeds a phase transition. ' + suggestion)
        oldPmid_MPa = Pmid_MPa + 0.0
        Pmid_MPa = GetPfreeze(oceanEOS, phaseMid, Tconv_K, UNDERPLATE=False)
        log.warning(f'Pmid_MPa has been adjusted upward from {oldPmid_MPa} to {Pmid_MPa} to compensate.')

    # Get melting temperature for calculating viscosity relative to this temp
    if phaseMid == 1:
        Tmelt_K = GetTfreeze(oceanEOS, Pmid_MPa, Tconv_K, TfreezeRange_K=275-Tconv_K)
    else:
        Tmelt_K = GetTfreeze(oceanEOS, Pmid_MPa, Tconv_K)
    etaConv_Pas = Constants.etaMelt_Pas[phaseMid] * np.exp(A * (Tmelt_K/Tconv_K - 1))
    # Get physical properties of ice at the "middle" of the convective region
    rhoMid_kgm3 = iceEOS.fn_rho_kgm3(Pmid_MPa, Tconv_K)
    CpMid_JkgK = iceEOS.fn_Cp_JkgK(Pmid_MPa, Tconv_K)
    alphaMid_pK = iceEOS.fn_alpha_pK(Pmid_MPa, Tconv_K)
    kMid_WmK = iceEOS.fn_kTherm_WmK(Pmid_MPa, Tconv_K)
    if iceEOS.POROUS:
        log.warning('Porosity corrections are not applied in calculating the Rayleigh number for convection models.')
    # Rayleigh number of whole ice layer, derived using viscosity of convective region
    Ra = alphaMid_pK * CpMid_JkgK * rhoMid_kgm3**2 * gtop_ms2 * (Tb_K - Ttop_K) * zb_m**3 / etaConv_Pas / kMid_WmK
    # Rayleigh number of lower thermal boundary layer, from parameterization results of Deschamps and Sotin (2000)
    Radelta = 0.28 * Ra**0.21
    # Thickness of lower thermal boundary layer
    deltaTBL_m = (etaConv_Pas * kMid_WmK * Radelta /
                  alphaMid_pK / CpMid_JkgK / rhoMid_kgm3**2 / gtop_ms2 / (Tb_K - Tconv_K))**(1/3)
    # Heat flux entering the bottom of the ice layer
    qBot_Wm2 = kMid_WmK * (Tb_K - Tconv_K) / deltaTBL_m
    # Heat flux leaving the top of the ice layer (adjusted for spherical geometry compared to Deschamps and Sotin, 2001)
    qTop_Wm2 = (rTop_m - zb_m)**2 / rTop_m**2 * qBot_Wm2
    # Thickness of conductive stagnant lid
    # This matches the Matlab, but based on Deschamps and Sotin (2001), who assume a fixed thermal conductivity
    # throughout the ice shell, the thickness of the conductive lid should probably use kTop_WmK, because that
    # will be what determines the conductive thermal profile based on the heat flux through the lid.
    #eLid_m = kMid_WmK * (Tconv_K - Ttop_K) / qTop_Wm2
    eLid_m = kTop_WmK * (Tconv_K - Ttop_K) / qTop_Wm2

    # If the Rayleigh number is less than some critical value, convection does not occur.
    RaCrit = GetRaCrit(Constants.Eact_kJmol[phaseBot], Tb_K, Ttop_K, Tconv_K)
    if(Ra < RaCrit):
        log.debug(f'Rayleigh number of {Ra:.3e} in the surface ice {PhaseConv(phaseBot)} ' +
                  f'layer is less than the critical value of {RaCrit:.3e}. ' +
                   'Only conduction will be modeled in this layer.')
        # Set conductive layer thicknesses to whole shell thickness to force a whole-layer conductive profile
        eLid_m = zb_m
        deltaTBL_m = 0.0
        Tconv_K = Ttop_K

    if not EQUIL_Q:
        # Set heat flux to be equal to that passing the conductive lid
        # according to Ojakangas and Stevenson (1989): https://doi.org/10.1016/0019-1035(89)90052-3
        #Qbot_W = kMid_WmK * Ttop_K / eLid_m * np.log(Tconv_K/Ttop_K) * 4*np.pi * (rTop_m - eLid_m)**2
        # Doing lower TBL instead because eLid_m gets set to 0 if Tconv > Ttop
        qBot_Wm2 = kMid_WmK * Tconv_K / deltaTBL_m * np.log(Tb_K/Tconv_K)

        # The commented out lines match the Matlab, but this is not what Deschamps
        # and Sotin (2001) do, it is not consistent with the results of Andersson
        # and Inaba (2005), and seems to be an incorrect evaluation of Eq. 2.7
        # from Ojakangas and Stevenson.
        #Dcond = np.array([np.nan, 632, 418, 242, np.nan, 328, 183])
        #qbot_Wm2 = Dcond[phase] * np.log(Tb_K/Ttop_K) / zb_m

    Dconv_m = zb_m - eLid_m - deltaTBL_m
    Qbot_W = qBot_Wm2 * 4*np.pi * (rTop_m - zb_m)**2

    return Tconv_K, etaConv_Pas, eLid_m, Dconv_m, deltaTBL_m, Qbot_W, Ra, RaCrit


def ConductiveTemperature(Ttop_K, rTop_m, rBot_m, kTherm_WmK, rhoRad_kgm3, Qrad_Wkg, Htidal_Wm3, qTop_Wm2):
    """ Thermal profile for purely thermally conductive layers, based on Turcotte and Schubert (1982),
        equation 4.40: T = -rho*H/6/k * r^2 + c1/r + c2, where c1 and c2 are integration constants
        found through boundary conditions, rho is mass density of the conductive layer in kg/m^3,
        H is internal heating in W/kg, and k is thermal conductivity in W/m/K.
        The main equations we use here are developed similar to equation 2 of
        Cammarano et al. (2006): https://doi.org/10.1029/2006JE002710, but we parameterize in terms of
        the heat flux leaving the top of the layer instead of entering the bottom. We also configure
        the internal heating so as to account for porosity by passing rho as the mass density of just
        the material contributing to radiogenic heat, i.e. silicates.

        Args:
            Ttop_K (float, shape N): Temperature at the top of the layer in K.
            rTop_m, rBot_m (float, shape N): Radius at top and bottom of layer in m, respectively.
            kTherm_WmK (float, shape N): Overall thermal conductivity of layer in W/(m K).
            rhoRad_kgm3 (float, shape N): Mass density of radiogenic material in layer in kg/m^3.
            Qrad_Wkg (float): Average radiogenic heating rate in W/kg.
            Htidal_Wm3 (float, shape N): Average tidal heating rate of the layer in W/m^3.
            qTop_Wm2 (float): Heat flux leaving the top of the layer in W/m^2.
        Returns:
            Tbot_K (float): Temperature at the bottom of the layer in K.
            qBot_Wm2 (float): Heat flux entering the bottom of the layer in W/m^2.
    """
    # Calculate needed values from inputs
    Htot_Wm3 = Qrad_Wkg * rhoRad_kgm3 + Htidal_Wm3
    c1 = qTop_Wm2 * rTop_m**2 / 2/kTherm_WmK - Htot_Wm3 / 6/kTherm_WmK * rTop_m**3
    # Find the temperature at the bottom of the layer
    Tbot_K = Ttop_K + Htot_Wm3 / 6/kTherm_WmK * (rTop_m**2 - rBot_m**2) + c1 * (1/rBot_m - 1/rTop_m)
    # The below calc is suspect. It seems to report too-high values for qBot.
    # Find the heat flux into the bottom of the layer
    #qBot_Wm2 = Htot_Wm3 / 3 * rBot_m + 2*kTherm_WmK / rBot_m**2 * c1
    # Find the approximate heat flux into the bottom of the layer
    qBot_Wm2 = kTherm_WmK * (Tbot_K - Ttop_K) / (rTop_m - rBot_m)

    return Tbot_K, qBot_Wm2


def ConductiveTemperatureActual(Ttop_K, rTop_m, rBot_m, kTherm_WmK, rhoRad_kgm3, Qrad_Wkg, Htidal_Wm3, qTop_Wm2):
    """ Thermal profile for purely thermally conductive layers, based on Turcotte and Schubert (1982),
        equation 4.40: T = -rho*H/6/k * r^2 + c1/r + c2, where c1 and c2 are integration constants
        found through boundary conditions, rho is mass density of the conductive layer in kg/m^3,
        H is internal heating in W/kg, and k is thermal conductivity in W/m/K.
        The main equations we use here are developed similar to equation 2 of
        Cammarano et al. (2006): https://doi.org/10.1029/2006JE002710, but we parameterize in terms of
        the heat flux leaving the top of the layer instead of entering the bottom. We also configure
        the internal heating so as to account for porosity by passing rho as the mass density of just
        the material contributing to radiogenic heat, i.e. silicates.

        Args:
            Ttop_K (float, shape N): Temperature at the top of the layer in K.
            rTop_m, rBot_m (float, shape N): Radius at top and bottom of layer in m, respectively.
            kTherm_WmK (float, shape N): Overall thermal conductivity of layer in W/(m K).
            rhoRad_kgm3 (float, shape N): Mass density of radiogenic material in layer in kg/m^3.
            Qrad_Wkg (float): Average radiogenic heating rate in W/kg.
            Htidal_Wm3 (float, shape N): Average tidal heating rate of the layer in W/m^3.
            qTop_Wm2 (float): Heat flux leaving the top of the layer in W/m^2.
        Returns:
            Tbot_K (float): Temperature at the bottom of the layer in K.
            qBot_Wm2 (float): Heat flux entering the bottom of the layer in W/m^2.
    """
    # Calculate needed values from inputs
    Htot_Wm3 = Qrad_Wkg * rhoRad_kgm3 + Htidal_Wm3
    c1 = qTop_Wm2 * rTop_m**2 / 2/kTherm_WmK - Htot_Wm3 / 6/kTherm_WmK * rTop_m**3
    # Find the temperature at the bottom of the layer
    Tbot_K = Ttop_K + Htot_Wm3 / 6/kTherm_WmK * (rTop_m**2 - rBot_m**2) + c1 * (1/rBot_m - 1/rTop_m)
    # The below calc is suspect. It seems to report too-high values for qBot.
    # Find the heat flux into the bottom of the layer
    qBot_Wm2 = Htot_Wm3 / 3 * rBot_m + 2*kTherm_WmK / rBot_m**2 * c1
    # Find the approximate heat flux into the bottom of the layer
    #qBot_Wm2 = kTherm_WmK * (Tbot_K - Ttop_K) / (rTop_m - rBot_m)

    return Tbot_K, qBot_Wm2


def GetRaCrit(Eact_kJmol, Tb_K, Ttop_K, Tconv_K):
    """ Calculates the critical Rayleigh number, above which convection is permitted
        based on Solomatov (1995) and reported as Eq. 3 of Hammond et al. (2016):
        https://doi.org/10.1002/2016GL069220

        Args:
            Eact_kJmol (float): Activation energy for diffusion in kJ/mol
            Tb_K (float): Temperature at bottom of ice layer in K
            Ttop_K (float): Temperature at top of ice layer in K
            Tconv_K (float): Average temperature of convecting layer in K
        Returns:
            RaCrit (float): Critical Rayleigh number, above which convection is permitted (dimensionless)
    """
    RaCrit = 20.9 * (Eact_kJmol*1e3 * (Tb_K - Ttop_K) / Constants.R / Tconv_K**2)**4

    return RaCrit


def GetPbConduct(Ttop_K, Tb_K, rTop_m, Ptop_MPa, gTop_ms2, qTop_Wm2, EOS, rRes_m=1e2, Qrad_Wkg=0, Htidal_Wm3=0):
    """ Find the pressure associated with the bottom temperature of a conductive layer,
        given the top temperature, outer radius, heat flux, and EOS.

        Args:
            Ttop_K (float): Temperature at the top of the conductive layer in K
            Tb_K (float): Temperature at the bottom of the conductive layer in K
            rTop_m (float): Radius of the top of the layer in m
            Ptop_MPa (float): Pressure at the top of the layer in MPa
            gTop_ms2 (float): Gravity at the top of the layer in m/s^2
            qTop_Wm2 (float): Heat flux leaving the top of the layer in W/m^2
            EOS (EOSStruct): Ice, Ocean, or Perple_X EOS struct that can be queried for
                density (fn_rho_kgm3) and thermal conductivity (fn_kTherm_WmK) as
                functions of P_MPa and T_K
            rRes_m = 1e2 (float): Radial resolution in m for calculating depth profile
            Htidal_Wm3 = 0 (float): Volumetric (tidal) heating rate in W/m^3
            Qrad_Wkg = 0 (float): Radiogenic heating rate in W/kg
    """
    Tbot_K = Ttop_K
    Pb_MPa = Ptop_MPa
    thisqTop_Wm2 = qTop_Wm2 + 0
    i = 0
    while Tbot_K < Tb_K:
        thisrTop_m = rTop_m - i*rRes_m
        rBot_m = thisrTop_m - rRes_m
        rho_kgm3 = EOS.fn_rho_kgm3(Pb_MPa, Tbot_K)
        kTherm_WmK = EOS.fn_kTherm_WmK(Pb_MPa, Tbot_K)
        Tbot_K, thisqTop_Wm2 = ConductiveTemperatureActual(Tbot_K, thisrTop_m, rBot_m,
                            kTherm_WmK, rho_kgm3, Qrad_Wkg, Htidal_Wm3, thisqTop_Wm2)
        MLayer_kg = 4/3 * np.pi * (thisrTop_m**3 - rBot_m**3) * rho_kgm3
        gTop_ms2 = (gTop_ms2 * thisrTop_m**2 - Constants.G * MLayer_kg) / rBot_m**2
        Pb_MPa += MLayer_kg * gTop_ms2 / (4*np.pi*rBot_m**2) / 1e6
        i += 1

    return Pb_MPa

