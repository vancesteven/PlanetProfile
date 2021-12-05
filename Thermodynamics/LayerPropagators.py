import numpy as np
from collections.abc import Iterable
from Utilities.dataStructs import Constants
from Thermodynamics.Hydrosphere import GetIceThermo, GetPfreeze, GetTfreeze, GetPfreezeHP

def IceLayers(Planet):
    """ Layer propagation from surface downward through the ice using geophysics.
        Iteratively sets up the thermal profile (the density and temperature)
        of the layer with each pressure step for all ice layers including
        ice Ih, ice III, ice V by calling the necessary subfunctions

        Requires Planet attributes:
            Do.BOTTOM_ICEIII
            Do.BOTTOM_ICEV
            Do.CLATHRATE
            Steps.nClath
            more in functions
        Assigns Planet attributes:
            Steps.nSurfIce, phase, r_m, z_m, g_ms2, T_K, P_MPa, rho_kgm3, Cp_JkgK, alpha_pK, PbI_MPa, Pb_MPa
    """
    Planet.Steps.nIbottom = Planet.Steps.nClath + Planet.Steps.nIceI
    Planet.Steps.nIIIbottom = Planet.Steps.nIbottom + Planet.Steps.nIceIIILitho
    Planet.Steps.nSurfIce = Planet.Steps.nIIIbottom + Planet.Steps.nIceVLitho
    # Assign phase values for near-surface ice Ih
    Planet.phase[:Planet.Steps.nClath] = 30  # Clathrate layers
    Planet.phase[Planet.Steps.nClath:Planet.Steps.nIbottom] = 1  # Ice Ih layers
    Planet.phase[Planet.Steps.nIbottom:Planet.Steps.nIIIbottom] = 3  # Ice III layers
    Planet.phase[Planet.Steps.nIIIbottom:Planet.Steps.nSurfIce] = 5  # Ice V layers

    # Layer property initialization for surface
    Planet.z_m[0] = 0.0  # Set first layer depth to zero (layer properties correspond to outer radius)
    Planet.r_m[0] = Planet.Bulk.R_m  # Set first layer to planetary surface radius
    Planet.g_ms2[0] = Constants.G * Planet.Bulk.M_kg / Planet.Bulk.R_m**2  # Set first layer gravity at surface
    Planet.T_K[0] = Planet.Bulk.Tsurf_K  # Set first layer surface temp
    Planet.P_MPa[0] = Planet.Bulk.Psurf_MPa  # Set first layer to surface pressure

    # Surface ice layer propagators -- if present, clathrates are nearest the surface
    # so we do clathrate layers first.
    if Planet.Do.CLATHRATE:
        ClathrateLayers(Planet)
    else:
        Planet.PbClath_MPa = Planet.Bulk.Psurf_MPa

    # Get the pressure consistent with the bottom of the ice Ih layer that is
    # consistent with the choice of Tb_K we suppose for this model
    Planet.PbI_MPa = GetPfreeze(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, Planet.Bulk.Tb_K,
        PfreezeLower_MPa=Planet.PfreezeLower_MPa, PfreezeUpper_MPa=Planet.PfreezeUpper_MPa, PfreezeRes_MPa=Planet.PfreezeRes_MPa)
    # Set linear P and adiabatic T in ice I layers
    PIceI_MPa = np.linspace(Planet.PbClath_MPa, Planet.PbI_MPa, Planet.Steps.nIceI)
    TIceI_K = Planet.Bulk.Tb_K**(PIceI_MPa/Planet.PbI_MPa) * Planet.T_K[Planet.Steps.nClath]**(1 - PIceI_MPa/Planet.PbI_MPa)
    Planet.P_MPa[Planet.Steps.nClath:Planet.Steps.nIbottom] = PIceI_MPa
    Planet.T_K[Planet.Steps.nClath:Planet.Steps.nIbottom] = TIceI_K

    # Additional adiabats in ice III and/or V underplate layers -- for thick, cold ice shells
    if Planet.Do.BOTTOM_ICEV:
        Planet.PbIII_MPa = GetPfreezeHP(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, Planet.Bulk.TbIII_K, 3)
        Planet.Pb_MPa = GetPfreezeHP(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, Planet.Bulk.TbV_K, 5)
        # Set linear P and adiabatic T in ice III layers
        PIceIII_MPa = np.linspace(Planet.P_MPa[Planet.Steps.nIbottom-1], Planet.PbIII_MPa, Planet.Steps.nIceIIILitho)
        TIceIII_K = Planet.Bulk.TbIII_K**(PIceIII_MPa/Planet.PbIII_MPa) * \
                    Planet.T_K[Planet.Steps.nIbottom-1]**(1 - PIceIII_MPa/Planet.PbIII_MPa)
        Planet.P_MPa[Planet.Steps.nIbottom:Planet.Steps.nIIIbottom] = PIceIII_MPa
        Planet.T_K[Planet.Steps.nIbottom:Planet.Steps.nIIIbottom] = TIceIII_K
        # Now ice V layers too
        PIceV_MPa = np.linspace(Planet.P_MPa[Planet.Steps.nIIIbottom-1], Planet.Pb_MPa, Planet.Steps.nIceVLitho)
        TIceV_K = Planet.Bulk.TbV_K**(PIceV_MPa/Planet.Pb_MPa) * \
                    Planet.T_K[Planet.Steps.nIIIbottom-1]**(1 - PIceV_MPa/Planet.Pb_MPa)
        Planet.P_MPa[Planet.Steps.nIIIbottom:Planet.Steps.nSurfIce] = PIceV_MPa
        Planet.T_K[Planet.Steps.nIIIbottom:Planet.Steps.nSurfIce] = TIceV_K
    elif Planet.Do.BOTTOM_ICEIII:
        Planet.Pb_MPa = GetPfreezeHP(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, Planet.Bulk.TbIII_K, 3)
        # Set linear P and adiabatic T in ice III layers
        PIceIII_MPa = np.linspace(Planet.P_MPa[Planet.Steps.nIbottom-1], Planet.PbIII_MPa, Planet.Steps.nIceIIILitho)
        TIceIII_K = Planet.Bulk.TbIII_K**(PIceIII_MPa/Planet.PbIII_MPa) * \
                    Planet.T_K[Planet.Steps.nIbottom-1]**(1 - PIceIII_MPa/Planet.PbIII_MPa)
        Planet.P_MPa[Planet.Steps.nIbottom:Planet.Steps.nIIIbottom] = PIceIII_MPa
        Planet.T_K[Planet.Steps.nIbottom:Planet.Steps.nIIIbottom] = TIceIII_K
    else:
        Planet.Pb_MPa = Planet.PbI_MPa

    # Evaluate thermodynamic properties of uppermost ice with SeaFreeze
    rho_surfIce_kgm3, Cp_surfIce_JkgK, alpha_surfIce_pK = GetIceThermo(Planet.P_MPa[:Planet.Steps.nSurfIce],
        Planet.T_K[:Planet.Steps.nSurfIce], Planet.phase[:Planet.Steps.nSurfIce])
    Planet.rho_kgm3[:Planet.Steps.nSurfIce] = rho_surfIce_kgm3
    Planet.Cp_JkgK[:Planet.Steps.nSurfIce] = Cp_surfIce_JkgK
    Planet.alpha_pK[:Planet.Steps.nSurfIce] = alpha_surfIce_pK
    # Calculate remaining physical properties of uppermost ice
    M_above_kg = 0
    for i in range(1, Planet.Steps.nSurfIce):
        # Increment depth based on change in pressure, combined with gravity and density
        Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
        Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
        M_above_kg += 4/3*np.pi * rho_surfIce_kgm3[i] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
        M_below_kg = Planet.Bulk.M_kg - M_above_kg
        Planet.g_ms2[i] = Constants.G * M_below_kg / Planet.r_m[i]**2
    Planet.zb_km = Planet.z_m[Planet.Steps.nSurfIce-1] / 1e3

    return Planet


def ClathrateLayers(Planet):
    """ For ice shells insulated by a layer of clathrate at the surface
        Calculates state variables of the layer with each pressure step

        Requires Planet attributes:
            Planet.Steps.nClath
        Assigns Planet attributes:
            phase
    """

    return Planet


def OceanLayers(Planet):
    """ Geophysical and thermodynamic calculations for ocean layer
        Calculates state variables of the layer with each pressure step

        Requires Planet attributes:
            ???
        Assigns Planet attributes:
            ???
    """

    return Planet


def InnerLayers(Planet):
    """ Geophysical and thermodynamic calculations for silicate and core layers
        Calculates state variables of the layer with each pressure step

        Requires Planet attributes:
            ???
        Assigns Planet attributes:
            ???
    """

    Planet = SilicateLayers(Planet)

    if Planet.Do.Fe_CORE:
        Planet = IronCoreLayers(Planet)

    return Planet


def SilicateLayers(Planet):
    """ Geophysical and thermodynamic calculations for silicate layers
        Calculates state variables of the layer with each pressure step

        Requires Planet attributes:
            ???
        Assigns Planet attributes:
            ???
    """

    Planet.Steps.indSil = 0
    Planet.Steps.nTotal = Planet.Steps.indSil + Planet.Steps.nSil
    return Planet


def IronCoreLayers(Planet):
    """ Geophysical and thermodynamic calculations for core layers
        Calculates state variables of the layer with each pressure step

        Requires Planet attributes:
            ???
        Assigns Planet attributes:
            ???
    """

    Planet.Steps.nTotal += Planet.Steps.nCore
    return Planet