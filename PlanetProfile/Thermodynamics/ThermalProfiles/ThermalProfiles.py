import numpy as np
import logging
from PlanetProfile.Utilities.defineStructs import Constants
from PlanetProfile.Thermodynamics.HydroEOS import GetTfreeze, GetPfreeze, kThermMelinder2007, \
    kThermHobbs1974, GetOceanEOS
from PlanetProfile.Utilities.Indexing import PhaseConv, MixedPhaseSeparator, PhaseInv

# Assign logger
log = logging.getLogger('PlanetProfile')

def ConvectionPetricca2024(Tmelt_K, Ttop_K, zb_m, Eact_kJmol, etaMelt_Pas, rhoMid_kgm3, alphaMid_pK, KthermMid_WmK, CpMid_JkgK, gtop_ms2, phaseBot, rTop_m, Htidal_Wm3=0):
    """ Thermodynamics calculations for convection in an ice layer
        based on Petricca et al. (2024): https://doi.org/10.1016/j.icarus.2024.116120
        Petricca presents a model for convection in an ice layer based on the approximation
        of convective temperature in the ice shell presented by Howell et al. (2021), approximation of
        viscosity in the convective layer given by the Arrhenius equation, and
        the assumption that the ice behavior is described by diffusion creep.
    """
    # Calculate approximate convective temperature using Howell et al. (2021)
    Tconv_K = (np.sqrt((4 * Tmelt_K * Constants.R / (Eact_kJmol * 1e3)) + 1) - 1) / (2 * Constants.R / (Eact_kJmol * 1e3))
    
    
    # Calculate viscosity by Arrhenius equation
    etaConv_Pas = etaMelt_Pas * np.exp(Eact_kJmol * 1e3 * (Tmelt_K / Tconv_K - 1) / Constants.R / Tmelt_K)
    
    # Calculate thermal diffusivity
    alphaMid_m2s = KthermMid_WmK / (rhoMid_kgm3 * CpMid_JkgK)
    
    # Calculate Rayleigh number
    Ra = alphaMid_pK * rhoMid_kgm3 * gtop_ms2 * (Tmelt_K - Ttop_K) * zb_m**3 / etaConv_Pas / alphaMid_m2s
    
    # Calcualte critical Rayleigh number
    RaCrit = GetRaCrit(Eact_kJmol, Tmelt_K, Ttop_K, Tconv_K)
    
    # Calculate curvature f, which is the ratio of the ice shell thickness to the radius of the planet
    f = 1 - zb_m / rTop_m
    # Calculate gamma, which describes the temperature dependence of viscosity
    gamma = Eact_kJmol * 1e3 * (Tmelt_K - Ttop_K) / (Constants.R * Tconv_K**2)
    # Calculate Nusselt number
    Nu = (1.46 * Ra**(0.27)) / ((gamma ** 1.21)*(f ** 0.78))
    # Calculate deltaTv
    deltaTv = (Tmelt_K - Ttop_K) / gamma
    # Calculate heat flux throughout entire ice shell
    qSurface_Wm2 = KthermMid_WmK * (Tmelt_K - Ttop_K) / (zb_m * f)
    # Calculate heat flux through the convective region
    qBot_Wm2 = 1.46 * Ra ** 0.27 / (f**1.78) * (deltaTv / (Tmelt_K - Ttop_K)) ** 1.21 * qSurface_Wm2
    # Add tidal heating contribution to heat flux (H_tidal * layer_thickness)
    if Htidal_Wm3 > 0:
        qBot_Wm2 = qBot_Wm2 + Htidal_Wm3 * zb_m

    # Calculate thickness of conductive portion
    eLid_m = KthermMid_WmK * (Tconv_K - Ttop_K) / qBot_Wm2
    
    if eLid_m > zb_m:
        log.warning('Conductive portion of ice shell is thicker than the ice shell. Only conduction will be modeled in this layer.')
        eLid_m = zb_m
        Tconv_K = Ttop_K
        Dconv_m = 0
    else:
        Dconv_m = zb_m - eLid_m

    # This method does not calculate the thickness of the lower thermalboundary layer, so we set it to 0
    deltaTBL_m = 0.0
    
    Qbot_W = qBot_Wm2 * 4*np.pi * (rTop_m - zb_m)**2
    return Tconv_K, etaConv_Pas, eLid_m, Dconv_m, deltaTBL_m, Qbot_W, Ra, RaCrit
        
    

def ConvectionKalousova2018(Ttop_K, rTop_m, kTop_WmK, Tb_K, zb_m, gtop_ms2, Pmid_MPa,
                            oceanEOS, iceEOS, phaseBot, EQUIL_Q, Eact_kJmol, qBot_Wm2=None,
                            Htidal_Wm3=0):
    """ Thermodynamics calculations for convection in HP ice layers (III, V, VI)
        based on Kalousova & Sotin (2018): https://doi.org/10.1029/2018GL078889

        Implements two-phase convection model for HP ice layers with potential partial melting.
        Key features:
        1. Temperature at ice-ocean interface follows melting curve (computed before calling this)
        2. Checks if modified Rayleigh number exceeds critical value for temperate layer formation
        3. If Ra* > Ra*c, a temperate (partially molten) layer exists at silicate/ice interface
        4. Returns convection parameters consistent with Deschamps & Sotin interface

        Args:
            Ttop_K (float): Temperature at top of layer (from ocean, follows melting curve) in K
            rTop_m (float): Radius of top of ice layer in m
            kTop_WmK (float): Thermal conductivity at top of ice layer in W/(m K)
            Tb_K (float): Temperature at bottom of layer (at silicate interface) in K
            zb_m (float): Thickness of the ice layer in m
            gtop_ms2 (float): Gravitational acceleration at layer top
            Pmid_MPa (float): Pressure at the "middle" of the convective region in MPa
            oceanEOS (OceanEOSStruct): Interpolator functions for evaluating the ocean EOS
            iceEOS (IceEOSStruct): Interpolator functions for evaluating the ice EOS
            phaseBot (int): Ice phase index at the bottom of the layer
            EQUIL_Q (bool): Whether to set heat flux consistent with convective profile
            Eact_kJmol (dict): Activation energies per phase in kJ/mol
            qBot_Wm2 (float, optional): Heat flux from silicate layer in W/m². If None, will be computed.

        Returns:
            Tconv_K (float): Interior temperature (follows melting curve in convecting region) in K
            etaConv_Pas (float): Reference viscosity at melting temperature in Pa*s
            eLid_m (float): Thickness of the top temperate layer (if present) in m
            Dconv_m (float): Thickness of convecting interior in m
            deltaTBL_m (float): Thickness of lower thermal boundary layer in m
            Qbot_W (float): Total heat flux at bottom of layer in W
            Ra (float): Modified Rayleigh number Ra*
            RaCrit (float): Critical Rayleigh number for temperate layer formation

        Note: Melt fraction is stored separately in Planet.meltFractionIII/V/VI
    """
    # Get phase info
    phaseMid = iceEOS.phaseID
    phaseMidString = PhaseConv(phaseMid)

    # Get reference viscosity at melting temperature
    # Use etaMelt_Pas from Constants for this phase
    if phaseMid > Constants.phaseClath and phaseMid < Constants.phaseClath + 10:
        # Mixed clathrate phase - use averaging
        stringClath, stringIce = MixedPhaseSeparator(PhaseConv(phaseMid))
        phaseClath, phaseIce = PhaseInv(stringClath), PhaseInv(stringIce)
        etaMelt_Pas = iceEOS.fn_averageValuesAccordingtoRule(Constants.etaMelt_Pas[phaseClath],
                                                               Constants.etaMelt_Pas[phaseIce],
                                                               'Carahan2004Averaging')
    else:
        etaMelt_Pas = Constants.etaMelt_Pas[phaseMid]

    # Get physical properties at mid-layer pressure and melting temperature
    # Interior temperature T^i follows the melting curve
    Tconv_K = Ttop_K  # Interior is at melting temperature (vigorous convection)
    rhoMid_kgm3 = iceEOS.fn_rho_kgm3(Pmid_MPa, Tconv_K)
    CpMid_JkgK = iceEOS.fn_Cp_JkgK(Pmid_MPa, Tconv_K)
    alphaMid_pK = iceEOS.fn_alpha_pK(Pmid_MPa, Tconv_K)
    kMid_WmK = iceEOS.fn_kTherm_WmK(Pmid_MPa, Tconv_K)

    # Temperature contrast across the convective layer
    # ΔTc = T_melt(silicate interface) - T_interior
    # T_interior follows melting curve, so ΔTc = Tb_K - Tconv_K
    DeltaTc_K = Tb_K - Tconv_K

    if DeltaTc_K <= 0:
        # No thermal contrast - no convection
        log.debug(f'No thermal contrast in HP ice {phaseMidString} layer (ΔT = {DeltaTc_K:.2f} K). '
                  'Convection absent.')
        eLid_m = 0.0
        Dconv_m = zb_m
        deltaTBL_m = 0.0
        Qbot_W = kMid_WmK * abs(DeltaTc_K) / zb_m * 4*np.pi * (rTop_m - zb_m)**2
        Ra = 0.0
        RaCrit = np.inf
        return Tconv_K, etaMelt_Pas, eLid_m, Dconv_m, deltaTBL_m, Qbot_W, Ra, RaCrit

    # Calculate far-field Rayleigh number (Kalousova & Sotin 2018, Eq. 2;
    # Sotin & Labrosse 1999): Ra* = α ρ g ΔTc Hc³ / (κ μ₀)
    # Equivalently in terms of conductivity: Ra* = α ρ² g Cp ΔTc Hc³ / (k μ₀)
    Hc_m = zb_m  # Layer thickness
    Ra_star = alphaMid_pK * CpMid_JkgK * rhoMid_kgm3**2 * gtop_ms2 * DeltaTc_K * Hc_m**3 / \
              (etaMelt_Pas * kMid_WmK)

    # Calculate heat flux if not provided
    if qBot_Wm2 is None:
        # Use conductive estimate
        qBot_Wm2 = kMid_WmK * DeltaTc_K / Hc_m

    # Add tidal heating contribution to heat flux (H_tidal * layer_thickness)
    if Htidal_Wm3 > 0:
        qBot_Wm2 = qBot_Wm2 + Htidal_Wm3 * Hc_m

    # Convert to mW/m²
    qs_mWm2 = qBot_Wm2 * 1e3

    # Calculate critical Rayleigh number for temperate layer formation (Eq. 7)
    # Ra*c = 19.965 × 10³ × (qs [mW/m²])^3.690
    RaCrit = 19.965e3 * (qs_mWm2**3.690)

    # Check if temperate layer forms
    if Ra_star > RaCrit:
        # Supercritical: temperate layer present - calculate thickness (Eq. 9)
        # Ht[km] = (0.145 × 10⁻³ × qs[mW/m²] + 0.015) × μ₀[Pa·s]^0.21
        Ht_km = (0.145e-3 * qs_mWm2 + 0.015) * (etaMelt_Pas**0.21)
        eLid_m = Ht_km * 1e3
        log.debug(f'HP ice {phaseMidString} temperate layer present: '
                  f'Ra* = {Ra_star:.3e} > Ra*c = {RaCrit:.3e}, '
                  f'Ht = {Ht_km:.3f} km')

        # Calculate TBL thickness using scaling (Eq. 4)
        # delta_i = 2.746 (Ra*)^(-0.271)
        # delta = delta_i * Hc
        delta_i = 2.746 * (Ra_star**(-0.271))
        deltaTBL_m = delta_i * Hc_m

        # Convecting layer thickness
        Dconv_m = zb_m - eLid_m - deltaTBL_m

        if Dconv_m < 0:
            # Boundary layers thicker than total layer - no convection possible
            log.warning(f'HP ice {phaseMidString} boundary layers ({(eLid_m + deltaTBL_m)/1e3:.1f} km) '
                        f'exceed layer thickness ({zb_m/1e3:.1f} km). Setting to conductive profile.')
            eLid_m = 0.0
            Dconv_m = 0.0
            deltaTBL_m = 0.0
    else:
        # Subcritical: entire layer is conductive, no convection
        eLid_m = zb_m
        Dconv_m = 0.0
        deltaTBL_m = 0.0
        log.debug(f'HP ice {phaseMidString} subcritical: '
                  f'Ra* = {Ra_star:.3e} < Ra*c = {RaCrit:.3e}. '
                  f'Entire {zb_m/1e3:.1f} km layer is conductive.')

    # Calculate heat flux
    Qbot_W = qBot_Wm2 * 4*np.pi * (rTop_m - zb_m)**2

    # Return convection parameters
    # etaConv is the reference viscosity at melting temperature
    return Tconv_K, etaMelt_Pas, eLid_m, Dconv_m, deltaTBL_m, Qbot_W, Ra_star, RaCrit


def ConvectionDeschampsSotin2001(Ttop_K, rTop_m, kTop_WmK, Tb_K, zb_m, gtop_ms2, Pmid_MPa,
                                 oceanEOS, iceEOS, phaseBot, EQUIL_Q, Eact_kJmol, Htidal_Wm3=0):
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
    if Tb_K < Ttop_K:
        raise ValueError('Tb_K is less than Ttop_K, which will result in a negative Rayleigh' +
                         'number. Try adjusting Tb_K, TbIII_K, TbV_K, etc.')

    # Get phase of convecting region from passed iceEOS
    phaseMid = iceEOS.phaseID
    phaseMidString = PhaseConv(phaseMid)
    # Numerical constants derived in Cartesian geometry from Deschamps and Sotin (2000) and used in
    # Deschamps and Sotin (2001) parameterization
    c1 = 1.43
    c2 = -0.03
    if not np.isnan(Eact_kJmol[phaseMidString]):
        # If we specify Eact_kJmol in Planet, then we should use those values, otherwise use constants
        A = Eact_kJmol[phaseMidString] * 1e3 / Constants.R / Tb_K
        B = Eact_kJmol[phaseMidString] * 1e3 / 2 / Constants.R / c1
    else:
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
    if phaseMid != oceanEOS.fn_phase(Pmid_MPa, Tconv_K) and Tconv_K > oceanEOS.Tmin:
        if not (abs(phaseMid) >= Constants.phaseClath and abs(phaseMid) < Constants.phaseClath + 10):
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
        Pmid_MPa = GetPfreeze(oceanEOS, phaseMid, Tconv_K, UNDERPLATE=False, HPNOOCEAN=False)
        log.warning(f'Pmid_MPa has been adjusted upward from {oldPmid_MPa} to {Pmid_MPa} to compensate.')
        if not np.isfinite(Pmid_MPa):
            log.warning(f'Pmid_MPa is not finite after adjustment. Returning conductive profile.')
            eLid_m = zb_m
            Dconv_m = 0.0
            deltaTBL_m = 0.0
            Qbot_W = kTop_WmK * abs(Tb_K - Ttop_K) / zb_m * 4*np.pi * (rTop_m - zb_m)**2
            return Tconv_K, Constants.etaMelt_Pas[phaseMid], eLid_m, Dconv_m, deltaTBL_m, Qbot_W, 0.0, np.inf

    # Get melting temperature for calculating viscosity relative to this temp
    Pmelt_MPa = np.arange(Pmid_MPa - 0.05*3, Pmid_MPa+0.05*3, 0.05)
    if phaseMid >= Constants.phaseClath and phaseMid < Constants.phaseClath + 10:
        meltEOS = oceanEOS
        Tupper_K = Tb_K
    else:
        Tupper_K = 274
        meltEOS = GetOceanEOS('PureH2O', 0.0, Pmelt_MPa,
                               np.arange(Tconv_K, 274.0, 0.05), None,
                               phaseType='calc', MELT=True)
    Tmelt_K = GetTfreeze(meltEOS, Pmid_MPa, Tconv_K, TfreezeRange_K=Tupper_K-Tconv_K)
    if phaseMid > Constants.phaseClath and phaseMid < Constants.phaseClath + 10:
        stringClath, stringIce = MixedPhaseSeparator(PhaseConv(phaseMid))
        phaseClath, phaseIce = PhaseInv(stringClath), PhaseInv(stringIce)
        etaMelt_Pas = iceEOS.fn_averageValuesAccordingtoRule(Constants.etaMelt_Pas[phaseClath], Constants.etaMelt_Pas[phaseIce], 'Carahan2004Averaging')
    else:
        etaMelt_Pas = Constants.etaMelt_Pas[phaseMid]
    etaConv_Pas = etaMelt_Pas * np.exp(A * (Tmelt_K/Tconv_K - 1))
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
    # Add tidal heating contribution to heat flux (H_tidal * layer_thickness)
    if Htidal_Wm3 > 0:
        qBot_Wm2 = qBot_Wm2 + Htidal_Wm3 * zb_m
    # Heat flux leaving the top of the ice layer (adjusted for spherical geometry compared to Deschamps and Sotin, 2001)
    qTop_Wm2 = (rTop_m - zb_m)**2 / rTop_m**2 * qBot_Wm2
    # Thickness of conductive stagnant lid
    # This matches the Matlab, but based on Deschamps and Sotin (2001), who assume a fixed thermal conductivity
    # throughout the ice shell, the thickness of the conductive lid should probably use kTop_WmK, because that
    # will be what determines the conductive thermal profile based on the heat flux through the lid.
    #eLid_m = kMid_WmK * (Tconv_K - Ttop_K) / qTop_Wm2
    eLid_m = kTop_WmK * (Tconv_K - Ttop_K) / qTop_Wm2
    phaseBotString = PhaseConv(phaseBot)
    if not np.isnan(Eact_kJmol[phaseBotString]):
        # Again, if we specify Eact_kJmol[phaseBot] in Planet, then we should use those values, otherwise use constants
        RaCrit = GetRaCrit(Eact_kJmol[phaseBotString], Tb_K, Ttop_K, Tconv_K)
    else:
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
        qBot_Wm2 = kMid_WmK * Tconv_K / eLid_m * np.log(Tb_K/Tconv_K)

        # The commented out lines match the Matlab, but this is not what Deschamps
        # and Sotin (2001) do, it is not consistent with the results of Andersson
        # and Inaba (2005), and seems to be an incorrect evaluation of Eq. 2.7
        # from Ojakangas and Stevenson.
        #Dcond = np.array([np.nan, 632, 418, 242, np.nan, 328, 183])
        #qbot_Wm2 = Dcond[phase] * np.log(Tb_K/Ttop_K) / zb_m

    Dconv_m = zb_m - eLid_m - deltaTBL_m
    Qbot_W = qBot_Wm2 * 4*np.pi * (rTop_m - zb_m)**2

    return Tconv_K, etaConv_Pas, eLid_m, Dconv_m, deltaTBL_m, Qbot_W, Ra, RaCrit


def ConvectionYao2014(Ttop_K, rTop_m, kTop_WmK, Tb_K, zb_m, gtop_ms2, Pmid_MPa,
                      oceanEOS, iceEOS, phaseBot, EQUIL_Q, Eact_kJmol, Htidal_Wm3=0):
    """ Stagnant-lid convection in spherical ice shells based on Yao et al. (2014):
        https://doi.org/10.1002/2014JE004653
        Uses 3D spherical geometry scaling laws with curvature parameter f.
        Inherently uses Arrhenius viscosity (Frank-Kamenetskii approximation).

        Args:
            Ttop_K (float): Temperature at top of ice layer in K
            rTop_m (float): Radius of top of ice layer in m
            kTop_WmK (float): Thermal conductivity at top of ice layer in W/(m K)
            Tb_K (float): Bottom temperature (melting point) in K
            zb_m (float): Thickness of the ice layer in m
            gtop_ms2 (float): Gravitational acceleration at layer top in m/s^2
            Pmid_MPa (float): Pressure at the middle of the convective region in MPa
            oceanEOS: Ocean EOS interpolator
            iceEOS: Ice EOS interpolator
            phaseBot (int): Ice phase index at the bottom of the layer
            EQUIL_Q (bool): Whether to use equilibrium heat flux formulation
            Eact_kJmol (dict): Activation energies per phase in kJ/mol
            Htidal_Wm3 (float): Volumetric tidal heating rate in W/m^3
        Returns:
            Tconv_K (float): Convective interior temperature in K
            etaConv_Pas (float): Viscosity of convective region in Pa*s
            eLid_m (float): Stagnant lid thickness in m
            Dconv_m (float): Convecting layer thickness in m
            deltaTBL_m (float): Bottom thermal boundary layer thickness in m
            Qbot_W (float): Total heat flux at bottom of ice in W
            Ra (float): Rayleigh number (viscous-temperature based)
            RaCrit (float): Critical Rayleigh number
    """
    if Tb_K < Ttop_K:
        raise ValueError('Tb_K is less than Ttop_K in ConvectionYao2014.')

    phaseMid = iceEOS.phaseID
    phaseMidString = PhaseConv(phaseMid)

    # Curvature parameter: f = R_base / R_top (Eq. 4 of Yao et al. 2014)
    rBot_m = rTop_m - zb_m
    f = rBot_m / rTop_m

    # Activation energy for the relevant phase
    if not np.isnan(Eact_kJmol[phaseMidString]):
        Eact_Jmol = Eact_kJmol[phaseMidString] * 1e3
    else:
        Eact_Jmol = Constants.Eact_kJmol[phaseMid] * 1e3

    # Frank-Kamenetskii viscosity contrast parameter (Eq. 6)
    gamma = Eact_Jmol / (Constants.R * Tb_K)

    # Solve for convective interior temperature iteratively (Eqs. 22, 31, 32)
    # Non-dimensional: theta_m = 1 - alpha_T / (gamma * f^beta)
    # with alpha_T = 1.23, beta = 1.5
    alpha_T = 1.23
    beta = 1.5
    DeltaT = Tb_K - Ttop_K
    theta_m = 1.0 - alpha_T / (gamma * f**beta)
    theta_m = np.clip(theta_m, 0.01, 0.99)
    Tconv_K = Ttop_K + theta_m * DeltaT

    if Tconv_K < Ttop_K:
        Tconv_K = Ttop_K
    if Tconv_K > Tb_K:
        Tconv_K = Tb_K - 0.1

    # Get melting temperature for Arrhenius viscosity reference
    Pmelt_MPa = np.arange(Pmid_MPa - 0.05*3, Pmid_MPa + 0.05*3, 0.05)
    Tupper_K = 274.0
    meltEOS = GetOceanEOS('PureH2O', 0.0, Pmelt_MPa,
                          np.arange(Tconv_K, 274.0, 0.05), None,
                          phaseType='calc', MELT=True)
    Tmelt_K = GetTfreeze(meltEOS, Pmid_MPa, Tconv_K, TfreezeRange_K=Tupper_K - Tconv_K)

    # Arrhenius viscosity at convective temperature (Eq. 30)
    etaMelt_Pas = Constants.etaMelt_Pas[phaseMid]
    etaConv_Pas = etaMelt_Pas * np.exp(Eact_Jmol / Constants.R * (1.0/Tconv_K - 1.0/Tmelt_K))

    # Material properties at mid-layer conditions
    rhoMid_kgm3 = iceEOS.fn_rho_kgm3(Pmid_MPa, Tconv_K)
    CpMid_JkgK = iceEOS.fn_Cp_JkgK(Pmid_MPa, Tconv_K)
    alphaMid_pK = iceEOS.fn_alpha_pK(Pmid_MPa, Tconv_K)
    kMid_WmK = iceEOS.fn_kTherm_WmK(Pmid_MPa, Tconv_K)
    kappaMid_m2s = kMid_WmK / (rhoMid_kgm3 * CpMid_JkgK)

    # Viscous temperature scale (Eq. 32)
    DeltaT_v = Constants.R * Tconv_K**2 / Eact_Jmol

    # Full-ΔT Rayleigh number for convection onset check (same definition as DS2001)
    Ra = (alphaMid_pK * CpMid_JkgK * rhoMid_kgm3**2 * gtop_ms2 * DeltaT * zb_m**3
          / (etaConv_Pas * kMid_WmK))

    # Viscous-temperature Rayleigh number for heat flux scaling (Eq. 34)
    Ra_m = (alphaMid_pK * rhoMid_kgm3 * gtop_ms2 * DeltaT_v * zb_m**3
            / (etaConv_Pas * kappaMid_m2s))

    # Bottom heat flux from scaling law (Eq. 35)
    # Phi_bot = a_F * Ra_m^b / f^d * (DeltaT_v/DeltaT)^c * Phi_c
    a_F = 1.46
    b = 0.27
    c = 1.21
    d = 1.78
    Phi_c = kMid_WmK * DeltaT / zb_m
    qBot_Wm2 = a_F * Ra_m**b / f**d * (DeltaT_v / DeltaT)**c * Phi_c

    # Add tidal heating contribution
    if Htidal_Wm3 > 0:
        qBot_Wm2 = qBot_Wm2 + Htidal_Wm3 * zb_m

    # Nusselt number (Eq. 27)
    Nu = qBot_Wm2 * zb_m / (kMid_WmK * DeltaT)

    # Stagnant lid thickness
    if Nu > 1.0:
        eLid_m = zb_m / Nu
    else:
        eLid_m = zb_m

    # Bottom thermal boundary layer thickness
    if qBot_Wm2 > 0:
        deltaTBL_m = kMid_WmK * (Tb_K - Tconv_K) / qBot_Wm2
    else:
        deltaTBL_m = 0.0

    # Critical Rayleigh number check
    phaseBotString = PhaseConv(phaseBot)
    if not np.isnan(Eact_kJmol[phaseBotString]):
        RaCrit = GetRaCrit(Eact_kJmol[phaseBotString], Tb_K, Ttop_K, Tconv_K)
    else:
        RaCrit = GetRaCrit(Constants.Eact_kJmol[phaseBot], Tb_K, Ttop_K, Tconv_K)

    if Ra < RaCrit:
        log.debug(f'Rayleigh number {Ra:.3e} < critical {RaCrit:.3e} in Yao2014 model. '
                  'Only conduction will be modeled.')
        eLid_m = zb_m
        deltaTBL_m = 0.0
        Tconv_K = Ttop_K

    Dconv_m = zb_m - eLid_m - deltaTBL_m
    if Dconv_m < 0:
        Dconv_m = 0.0
        eLid_m = zb_m
        deltaTBL_m = 0.0

    # Total heat flux through bottom spherical surface
    Qbot_W = qBot_Wm2 * 4 * np.pi * rBot_m**2

    return Tconv_K, etaConv_Pas, eLid_m, Dconv_m, deltaTBL_m, Qbot_W, Ra, RaCrit


def ConductiveTemperature(Ttop_K, rTop_m, rBot_m, kTherm_WmK, rhoRad_kgm3, Qrad_Wkg, Htidal_Wm3, qTop_Wm2):
    """ Thermal profile for purely thermally conductive layers, based on Turcotte and
        Schubert (2002) equation 4.40.  The main equations here were developed similar
        to equation 2 of Cammarano et al. (2006): https://doi.org/10.1029/2006JE002710,
        parameterized in terms of the heat flux leaving the top of the layer rather
        than entering the bottom.

        *** KNOWN c1 FACTOR-OF-2 DISCREPANCY — see `ConductiveTemperatureCorrect` ***

        The c1 formula below has `/2` and `/6` in the denominators, which do not
        match a direct derivation from T&S 4.40 (which would give `/1` and `/3`).
        The thin-layer planar limit of this function gives Delta_T = qTop * Delta_R / (2k),
        half of Fourier's law qTop * Delta_R / k.  `ConductiveTemperatureActual` below
        carries a compensating *2 prefactor in its spherical qBot term that cancels
        the /2 in c1 — so its qBot return is correct, even though Tbot is not.

        This function is now used by `SilRecursionSolid` / `SilRecursionPorous` only
        for **shell** silicate bodies (Fe_CORE=True or CONSTANT_INNER_DENSITY=True).
        For those cases the c1/r term is well-defined (rBot > 0) and the legacy
        halved-c1 form is preserved to avoid disturbing existing test-suite behavior.
        Solid-sphere silicate bodies (no Fe core, no CONSTANT_INNER_DENSITY) now use
        `ConductiveTemperatureCorrect` with an overridden qTop = Htot*rTop/3 (which
        forces c1 = 0 and yields the closed-form profile finite at r=0); see
        `Geophysical.py` for that path.

        For the narrow use-case where we need the true T&S 4.40 behavior (namely,
        `GetPbConduct` which integrates downward from the top of the clathrate layer
        until reaching a target Tbot), use `ConductiveTemperatureCorrect` below.

        Args:
            Ttop_K (float, shape N): Temperature at the top of the layer in K.
            rTop_m, rBot_m (float, shape N): Radius at top and bottom of layer in m.
            kTherm_WmK (float, shape N): Thermal conductivity of layer in W/(m K).
            rhoRad_kgm3 (float, shape N): Mass density of the radiogenic material in kg/m^3.
            Qrad_Wkg (float): Average radiogenic heating rate in W/kg.
            Htidal_Wm3 (float, shape N): Tidal heating rate of the layer in W/m^3.
            qTop_Wm2 (float): Heat flux leaving the top of the layer in W/m^2.
        Returns:
            Tbot_K (float): Temperature at the bottom of the layer in K.
            qBot_Wm2 (float): Heat flux entering the bottom of the layer in W/m^2.
    """
    # Calculate needed values from inputs
    Htot_Wm3 = Qrad_Wkg * rhoRad_kgm3 + Htidal_Wm3
    # NOTE: /2 and /6 are known to differ from a strict T&S 4.40 derivation; see
    # module docstring above and `ConductiveTemperatureCorrect`.  Do not change
    # without addressing the silicate BC regression.
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
    """ Spherical-form counterpart of `ConductiveTemperature` — same c1 as above, but
        returns qBot from the exact spherical formula with a compensating *2 prefactor
        that cancels the /2 in c1.  That makes qBot correct (matches Fourier's law in
        the thin-layer limit) even though the Tbot formula is not.

        *** SEE `ConductiveTemperature` DOCSTRING FOR THE c1 DISCREPANCY ***

        Args:
            (same as ConductiveTemperature)
        Returns:
            (same as ConductiveTemperature)
    """
    # Calculate needed values from inputs
    Htot_Wm3 = Qrad_Wkg * rhoRad_kgm3 + Htidal_Wm3
    c1 = qTop_Wm2 * rTop_m**2 / 2/kTherm_WmK - Htot_Wm3 / 6/kTherm_WmK * rTop_m**3
    # Find the temperature at the bottom of the layer
    Tbot_K = Ttop_K + Htot_Wm3 / 6/kTherm_WmK * (rTop_m**2 - rBot_m**2) + c1 * (1/rBot_m - 1/rTop_m)
    # Spherical qBot: the *2 prefactor here cancels the /2 in c1, yielding the correct
    # Fourier flux.  Do not remove the *2 without also removing the /2 from c1.
    qBot_Wm2 = Htot_Wm3 / 3 * rBot_m + 2*kTherm_WmK / rBot_m**2 * c1

    return Tbot_K, qBot_Wm2


def ConductiveTemperatureCorrect(Ttop_K, rTop_m, rBot_m, kTherm_WmK, rhoRad_kgm3, Qrad_Wkg, Htidal_Wm3, qTop_Wm2):
    """ Strict T&S 4.40 spherical-conduction step — the mathematically correct version
        of `ConductiveTemperatureActual`.

        T(r) = -Htot/(6 k) * r^2 + c1/r + c2

        Derivation of c1 from the boundary condition q(rTop) = qTop, with
        q(r) = -k dT/dr = Htot * r / 3 + k * c1 / r^2:

            c1 = qTop * rTop^2 / k  -  Htot * rTop^3 / (3 k)

        Planar limit (rBot → rTop): Delta_T = qTop * Delta_R / k (Fourier's law).

        This function is used by `GetPbConduct` to integrate downward from the top of
        a clathrate layer until T reaches a target value.  Using this (correct) form
        makes the resulting clathrate layer thickness match `Bulk.clathMaxThick_m`.

        The older `ConductiveTemperature` / `ConductiveTemperatureActual` above use
        halved c1; they are retained with that behavior for **shell** silicate bodies
        (Fe core or CONSTANT_INNER_DENSITY) where the c1/r term is well-defined.
        Solid-sphere silicate bodies were reworked in 2026-05 to use this function
        with an overridden qTop = Htot*rTop/3 that forces c1 = 0 and yields the
        closed-form profile finite at r=0 (see `SilRecursionSolid` /
        `SilRecursionPorous` in Geophysical.py).

        Args:
            (same as ConductiveTemperature)
        Returns:
            (same as ConductiveTemperature — Tbot_K, qBot_Wm2 both from T&S 4.40)
    """
    Htot_Wm3 = Qrad_Wkg * rhoRad_kgm3 + Htidal_Wm3
    # Strict T&S 4.40 c1 from q(rTop) = qTop boundary condition.
    c1 = qTop_Wm2 * rTop_m**2 / kTherm_WmK - Htot_Wm3 * rTop_m**3 / (3 * kTherm_WmK)
    Tbot_K = Ttop_K + Htot_Wm3 / (6 * kTherm_WmK) * (rTop_m**2 - rBot_m**2) + c1 * (1/rBot_m - 1/rTop_m)
    # Exact spherical qBot from q(r) = Htot * r / 3 + k c1 / r^2 at rBot.
    qBot_Wm2 = Htot_Wm3 * rBot_m / 3 + kTherm_WmK * c1 / rBot_m**2

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
        # Use the T&S 4.40 correct form here (not the legacy halved-c1
        # `ConductiveTemperatureActual`).  The correct form makes the GetPbConduct
        # integration terminate at the physically expected depth, so the clathrate
        # layer thickness matches `Bulk.clathMaxThick_m` (was 2x too deep prior).
        Tbot_K, thisqTop_Wm2 = ConductiveTemperatureCorrect(Tbot_K, thisrTop_m, rBot_m,
                            kTherm_WmK, rho_kgm3, Qrad_Wkg, Htidal_Wm3, thisqTop_Wm2)
        MLayer_kg = 4/3 * np.pi * (thisrTop_m**3 - rBot_m**3) * rho_kgm3
        gTop_ms2 = (gTop_ms2 * thisrTop_m**2 - Constants.G * MLayer_kg) / rBot_m**2
        Pb_MPa += MLayer_kg * gTop_ms2 / (4*np.pi*rBot_m**2) / 1e6
        i += 1

    return Pb_MPa
