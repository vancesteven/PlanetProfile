import numpy as np
from collections.abc import Iterable
from Utilities.dataStructs import Constants
from Thermodynamics.HydroEOS import GetIceThermo, GetPfreeze, GetTfreeze, GetPfreezeHP, FluidEOS, GetPhase, OceanEOSStruct
from Thermodynamics.FromLiterature.ThermalProfiles import ConductionClathLid, ConvectionDeschampsSotin2001, \
    ConductiveTemperature, kThermIsobaricAnderssonIbari2005
from Thermodynamics.InnerEOS import PerplexEOSStruct

def IceLayers(Planet, Params):
    """ Layer propagation from surface downward through the ice using geophysics.
        Iteratively sets up the thermal profile (the density and temperature)
        of the layer with each pressure step for all ice layers including
        ice Ih, ice III, ice V by calling the necessary subfunctions

        Assigns Planet attributes:
            Steps.nSurfIce, phase, r_m, z_m, g_ms2, T_K, P_MPa, rho_kgm3, Cp_JkgK, alpha_pK, MLayer_kg, PbI_MPa, Pb_MPa
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

    # Get the pressure consistent with the bottom of the ice Ih layer that is
    # consistent with the choice of Tb_K we suppose for this model
    if Params.VERBOSE: print('Finding uppermost ice melting pressure...')
    Planet.PbI_MPa = GetPfreeze(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, Planet.Bulk.Tb_K,
        PfreezeLower_MPa=Planet.PfreezeLower_MPa, PfreezeUpper_MPa=Planet.PfreezeUpper_MPa, PfreezeRes_MPa=Planet.PfreezeRes_MPa,
        Pguess=None)

    # Now, we want to check for a convective profile. First, we need to get zb_km, so we need to suppose
    # a whole-layer conductive profile. The densities will change slightly, so we depart from self-consistency
    # here. Repeated applications of IceConvect will get closer to self-consistency.

    # Surface ice layer propagators -- if present, clathrates are nearest the surface
    # so we do clathrate layers first.
    if Planet.Do.CLATHRATE:
        if Params.VERBOSE: print('Evaluating clathrate layers.')
        ClathrateLayers(Planet, Params)
    else:
        Planet.PbClath_MPa = Planet.Bulk.Psurf_MPa
        Planet.T_K[Planet.Steps.nClath] = Planet.Bulk.Tsurf_K
        Planet.zClath_m = 0

    # Set linear P and adiabatic T in ice I layers
    PIceI_MPa = np.linspace(Planet.PbClath_MPa, Planet.PbI_MPa, Planet.Steps.nIceI+1)
    TIceI_K = Planet.Bulk.Tb_K**(PIceI_MPa[:-1]/Planet.PbI_MPa) * Planet.T_K[Planet.Steps.nClath]**(1 - PIceI_MPa[:-1]/Planet.PbI_MPa)
    Planet.P_MPa[Planet.Steps.nClath:Planet.Steps.nIbottom] = PIceI_MPa[:-1]
    Planet.T_K[Planet.Steps.nClath:Planet.Steps.nIbottom] = TIceI_K

    # Evaluate thermodynamic properties of uppermost ice with SeaFreeze
    rhoIceI_kgm3, CpIceI_JkgK, alphaIceI_pK = GetIceThermo(PIceI_MPa[:-1], TIceI_K,
                            Planet.phase[Planet.Steps.nClath:Planet.Steps.nIbottom])
    Planet.rho_kgm3[Planet.Steps.nClath:Planet.Steps.nIbottom] = rhoIceI_kgm3
    Planet.Cp_JkgK[Planet.Steps.nClath:Planet.Steps.nIbottom] = CpIceI_JkgK
    Planet.alpha_pK[Planet.Steps.nClath:Planet.Steps.nIbottom] = alphaIceI_pK

    # Calculate remaining physical properties of upper ice I
    thisMAbove_kg = 0
    for i in range(1, Planet.Steps.nIbottom):
        # Increment depth based on change in pressure, combined with gravity and density
        Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
        Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
        Planet.MLayer_kg[i-1] = 4/3*np.pi * rhoIceI_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
        thisMAbove_kg += Planet.MLayer_kg[i-1]
        thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
        Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2
        if Params.VERBOSE: print('il: ' + str(i) +
                                 '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                 '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                 '; phase: ' + str(Planet.phase[i]))

    if Params.VERBOSE: print('Upper ice initial conductive profile complete.')

    if not Planet.Do.NO_ICEI_CONVECTION:
        # Record zb_m to see if it gets adjusted significantly
        zbOld_m = Planet.z_m[Planet.Steps.nIbottom-1] + 0.0
        # Now check for convective region and get dimensions if present
        Planet = IceConvect(Planet, Params)
        # Run IceConvect a second time if zb_km changed by more than a set tolerance
        if(np.abs(Planet.z_m[Planet.Steps.nIbottom] - zbOld_m)/Planet.z_m[Planet.Steps.nIbottom] > Planet.Bulk.zbChangeTol_frac):
            if Params.VERBOSE: print('The bottom depth of surface ice I changed by ' +
                str(round((Planet.z_m[Planet.Steps.nIbottom] - zbOld_m)/1e3,2)) + ' km from IceConvect, which is greater than ' +
                str(Planet.Bulk.zbChangeTol_frac * 100) + '%. running IceConvect a second time...')
            Planet = IceConvect(Planet, Params)
    else:
        if Params.VERBOSE: print('NO_ICEI_CONVECTION is True -- skipping convection calculations.')
        # Find the surface heat flux from the conductive profile. This assumes there is no tidal heating!
        qSurf_Wm2 = (Planet.T_K[1] - Planet.Bulk.Tsurf_K) / (Planet.Bulk.R_m - Planet.r_m[1]) * Planet.kTherm_WmK[0]
        Planet.Ocean.QfromMantle_W = qSurf_Wm2 * 4*np.pi * Planet.Bulk.R_m**2

    # Additional adiabats in ice III and/or V underplate layers -- for thick, cold ice shells
    if Planet.Do.BOTTOM_ICEV:
        if Params.VERBOSE: print('...with ice III and V underplating...')
        Planet = IceVUnderplate(Planet, Params)
    elif Planet.Do.BOTTOM_ICEIII:
        if Params.VERBOSE: print('...with ice III underplating...')
        Planet = IceIIIUnderplate(Planet, Params)
    else:
        Planet.Pb_MPa = Planet.PbI_MPa

    if Params.VERBOSE: print('Upper ice melting pressure: '+str(Planet.Pb_MPa))

    Planet.zb_km = Planet.z_m[Planet.Steps.nSurfIce-1] / 1e3

    return Planet


def ClathrateLayers(Planet, Params):
    """ For ice shells insulated by a layer of clathrate at the surface
        Calculates state variables of the layer with each pressure step

        Requires Planet attributes:
            Planet.Steps.nClath
        Assigns Planet attributes:
            phase
    """

    if Params.VERBOSE: print('Applying clathrate lid conduction.')
    Planet = ConductionClathLid()

    return Planet


def IceConvect(Planet, Params):
    """ Apply convection models from literature to determine thermal profile and
        state variables for possibly convecting ice layers

        Assigns Planet attributes:
            Tconv_K, etaConv_Pas, eLid_m, deltaTBL_m, QfromMantle_W, all physical layer arrays
    """

    if Params.VERBOSE: print('Applying solid-state convection to surface ice I based on Deschamps and Sotin (2001).')
    zbI_m = Planet.z_m[Planet.Steps.nIbottom-1]
    # Get "middle" pressure
    Pmid_MPa = (Planet.PbI_MPa - Planet.PbClath_MPa) / 2
    # Get thermal conductivity at surface
    if Planet.Do.CLATHRATE:
        Planet.kTherm_WmK[0] = Constants.kClath_WmK
    else:
        Planet.kTherm_WmK[0] = kThermIsobaricAnderssonIbari2005(Planet.Bulk.Tsurf_K, 1)

    # Run calculations to get convection layer parameters
    Planet.Tconv_K, Planet.etaConv_Pas, Planet.eLid_m, Planet.deltaTBL_m, qbot_Wm2, Planet.RaConvect = \
        ConvectionDeschampsSotin2001(Planet.Bulk.Tsurf_K, Planet.Bulk.R_m, Planet.kTherm_WmK[0],
                                     Planet.Bulk.Tb_K, zbI_m, Planet.g_ms2[0], Pmid_MPa, 1)

    # Convert bottom heat flux into total heat flow rate so it's normalized to the area it passes through
    Planet.Ocean.QfromMantle_W = 4*np.pi * (Planet.Bulk.R_m - zbI_m)**2 * qbot_Wm2

    if(Planet.zClath_m > Planet.eLid_m):
        Planet.Bulk.clathMaxDepth_m = Planet.eLid_m
        print('WARNING: Clathrate lid thickness was greater than the conductive lid thickness.' +
              'Planet.Bulk.clathMaxDepth_m has been reduced to be equal to the conductive lid thickness.')

    # Check for whole-lid conduction
    if(zbI_m <= Planet.eLid_m + Planet.deltaTBL_m):
        if Params.VERBOSE: print('Ice shell thickness is less than that of the thermal boundary layers, convection is absent.\n' +
                                 'Applying whole-shell conductive profile.')

        # Recalculate heat flux, as it will be too high for conduction-only:
        qSurf_Wm2 = (Planet.T_K[1] - Planet.Bulk.Tsurf_K) / (Planet.Bulk.R_m - Planet.r_m[1]) * Planet.kTherm_WmK[0]
        Planet.Ocean.QfromMantle_W = qSurf_Wm2 * 4*np.pi * Planet.Bulk.R_m**2

        # We again need to do clathrate layers first if there are any
        if Planet.Do.CLATHRATE:
            if Params.VERBOSE: print('Evaluating clathrate layers.')
            ClathrateLayers(Planet, Params)
        else:
            thisMAbove_kg = 0

        Planet.z_m[Planet.Steps.nClath] = Planet.zClath_m + 0.0
        Planet.r_m[Planet.Steps.nClath] = Planet.Bulk.R_m - Planet.z_m[Planet.Steps.nClath]
        thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
        Planet.g_ms2[Planet.Steps.nClath] = Constants.G * thisMBelow_kg / Planet.r_m[Planet.Steps.nClath]**2
        qTop_Wm2 = Planet.Ocean.QfromMantle_W / 4/np.pi/Planet.r_m[Planet.Steps.nClath]**2

        # Evaluate thermodynamic properties at top of ice I with SeaFreeze
        Planet.rho_kgm3[Planet.Steps.nClath], Planet.Cp_JkgK[Planet.Steps.nClath], Planet.alpha_pK[Planet.Steps.nClath] \
            = GetIceThermo([Planet.P_MPa[Planet.Steps.nClath]], [Planet.T_K[Planet.Steps.nClath]], Planet.phase[Planet.Steps.nClath])

        for i in range(Planet.Steps.nClath+1, Planet.Steps.nIbottom):
            Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
            Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
            Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
            thisMAbove_kg += Planet.MLayer_kg[i-1]
            thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
            Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2
            Planet.kTherm_WmK[i] = Planet.kTherm_WmK[i-1]  # Placeholder until a self-consistent determination is implemented

            Planet.T_K[i], qTop_Wm2 = ConductiveTemperature(Planet.T_K[i-1], Planet.r_m[i-1], Planet.r_m[i],
                                                            Planet.kTherm_WmK[i], Planet.rho_kgm3[i], 0, 0, qTop_Wm2)
            Planet.rho_kgm3[i], Planet.Cp_JkgK[i], Planet.alpha_pK[i] = GetIceThermo([Planet.P_MPa[i]], [Planet.T_K[i]],
                                                                                      Planet.phase[i])
            if Params.VERBOSE: print('il: ' + str(i) +
                                     '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                     '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                     '; phase: ' + str(Planet.phase[i]))

        # Force final temperature value to be the assigned value
        Planet.T_K[Planet.Steps.nIbottom] = Planet.Bulk.Tb_K

    else:
        # Model conduction in the stagnant lid
        # We now know the temperature difference and thickness of the upper thermal
        # boundary layer. We first need to reconfigure the profile and set a linear
        # depth dependence, so we can use ConductiveTemperature to propagate things
        nConduct = int(Planet.Steps.nIbottom * Planet.Steps.eTBL_frac)
        if(Planet.Steps.nClath > nConduct):
            raise ValueError('Planet.Steps.nClath is greater than the fraction of surface ice layers used for the ' +
                             'upper conductive portion of the shell. Either increase Planet.Steps.eTBL_frac or ' +
                             'decrease Planet.Steps.nClath and re-run to solve this problem.')
        Planet.z_m[:nConduct] = np.linspace(0, Planet.eLid_m, nConduct)
        Planet.r_m[:nConduct] = Planet.Bulk.R_m - Planet.z_m[:nConduct]

        if Planet.Do.CLATHRATE:
            if Params.VERBOSE: print('Evaluating clathrate layers in stagnant lid.')
            ClathrateLayers(Planet, Params)
        else:
            thisMAbove_kg = 0

        thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
        Planet.g_ms2[Planet.Steps.nClath] = Constants.G * thisMBelow_kg / Planet.r_m[Planet.Steps.nClath]**2
        qTop_Wm2 = Planet.Ocean.QfromMantle_W / 4/np.pi/Planet.r_m[Planet.Steps.nClath]**2

        if Params.VERBOSE: print('Modeling conduction in stagnant lid...')
        for i in range(Planet.Steps.nClath+1, nConduct):
            Planet.rho_kgm3[i-1], Planet.Cp_JkgK[i-1], Planet.alpha_pK[i-1] \
                = GetIceThermo([Planet.P_MPa[i-1]], [Planet.T_K[i-1]], Planet.phase[i-1])
            Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
            thisMAbove_kg += Planet.MLayer_kg[i-1]
            thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
            Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2
            Planet.P_MPa[i] = Planet.P_MPa[i-1] + Planet.MLayer_kg[i-1] * Planet.g_ms2[i-1] / (4*np.pi*Planet.r_m[i-1]**2) / 1e6
            Planet.kTherm_WmK[i] = Planet.kTherm_WmK[i-1]  # Placeholder until a self-consistent determination is implemented

            Planet.T_K[i], qTop_Wm2 = ConductiveTemperature(Planet.T_K[i-1], Planet.r_m[i-1], Planet.r_m[i],
                                                            Planet.kTherm_WmK[i], Planet.rho_kgm3[i], 0, 0, qTop_Wm2)
            if Params.VERBOSE: print('il: ' + str(i) +
                                     '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                     '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                     '; phase: ' + str(Planet.phase[i]))

        # Force final temperature value to be the determined value
        Planet.T_K[nConduct] = Planet.Tconv_K
        Planet.rho_kgm3[nConduct], Planet.Cp_JkgK[nConduct], Planet.alpha_pK[nConduct] \
            = GetIceThermo([Planet.P_MPa[nConduct]], [Planet.T_K[nConduct]], Planet.phase[nConduct])

        # Next, we need the pressure difference across the lower thermal boundary layer.
        # We use the "old" radius and gravity here, so we depart from self-consistency.
        # Repeated iterations will get closer to self-consistent.
        nLowerTBL = int(Planet.Steps.nIbottom * Planet.Steps.deltaTBL_frac)
        rhoTBL_kgm3, _, _ = GetIceThermo([Planet.PbI_MPa], [Planet.Bulk.Tb_K], 1)
        gTBL_ms2 = Planet.g_ms2[Planet.Steps.nIbottom]
        rTBLbot_m = Planet.Bulk.R_m - zbI_m
        rTBLtop_m = rTBLbot_m + Planet.deltaTBL_m
        DeltaPlowerTBL_MPa = rhoTBL_kgm3[0] * gTBL_ms2 / 3 * (rTBLtop_m**3 - rTBLbot_m**3) / rTBLbot_m**2 / 1e6
        PtopLowerTBL_MPa = Planet.PbI_MPa + DeltaPlowerTBL_MPa

        # Now, we set the convecting layer to be isothermal and assign a linear
        # increase in pressure across the layers
        nConvect = Planet.Steps.nIbottom - nConduct - nLowerTBL
        if(nConvect <= 0): raise ValueError('Planet.Steps.eTBL_frac + Planet.Steps.deltaTBL_frac exceeds 1.')
        Planet.T_K[nConduct:nConduct+nConvect] = Planet.Tconv_K
        Planet.P_MPa[nConduct:nConduct+nConvect+1] = np.linspace(Planet.P_MPa[nConduct], PtopLowerTBL_MPa, nConvect+1)
        # Evaluate thermodynamic properties with SeaFreeze
        rhoConvect_kgm3, CpConvect_JkgK, alphaConvect_pK \
            = GetIceThermo(Planet.P_MPa[nConduct:nConduct+nConvect], Planet.T_K[nConduct:nConduct+nConvect],
                           Planet.phase[nConduct:nConduct+nConvect])
        Planet.rho_kgm3[nConduct:nConduct+nConvect] = rhoConvect_kgm3
        Planet.Cp_JkgK[nConduct:nConduct+nConvect] = CpConvect_JkgK
        Planet.alpha_pK[nConduct:nConduct+nConvect] = alphaConvect_pK

        # Calculate remaining physical properties of convecting layer and find depths
        if Params.VERBOSE: print('Modeling isothermal convecting layer...')
        for i in range(nConduct, nConduct + nConvect):
            Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
            Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
            Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
            thisMAbove_kg += Planet.MLayer_kg[i-1]
            thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
            Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2
            Planet.kTherm_WmK[i] = Planet.kTherm_WmK[i-1]  # Placeholder until a self-consistent determination is implemented
            if Params.VERBOSE: print('il: ' + str(i) +
                                     '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                     '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                     '; phase: ' + str(Planet.phase[i]))

        # Finally, we get properties for the lower TBL layers and add them to the profile
        PlowerTBL_MPa = np.linspace(Planet.P_MPa[nConduct+nConvect], Planet.PbI_MPa, nLowerTBL+1)
        Planet.P_MPa[nConduct+nConvect:Planet.Steps.nIbottom] = PlowerTBL_MPa[:-1]
        TlowerTBL_K = Planet.Bulk.Tb_K**(PlowerTBL_MPa[:-1]/Planet.PbI_MPa) \
                      * Planet.T_K[nConduct+nConvect]**(1 - PlowerTBL_MPa[:-1]/Planet.PbI_MPa)
        # Evaluate thermodynamic properties of remaining ice with SeaFreeze
        rhoLowerTBL_kgm3, CpLowerTBL_JkgK, alphaLowerTBL_pK = GetIceThermo(PlowerTBL_MPa[:-1], TlowerTBL_K,
                   Planet.phase[nConduct+nConvect:Planet.Steps.nIbottom])
        Planet.rho_kgm3[nConduct+nConvect:Planet.Steps.nIbottom] = rhoLowerTBL_kgm3
        Planet.Cp_JkgK[nConduct+nConvect:Planet.Steps.nIbottom] = CpLowerTBL_JkgK
        Planet.alpha_pK[nConduct+nConvect:Planet.Steps.nIbottom] = alphaLowerTBL_pK
        # Get remaining physical properties of final layer

        if Params.VERBOSE: print('Modeling conduction in lower thermal boundary layer...')
        for i in range(nConduct+nConvect+1, Planet.Steps.nIbottom):
            # Increment depth based on change in pressure, combined with gravity and density
            Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
            Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
            Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
            thisMAbove_kg += Planet.MLayer_kg[i-1]
            thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
            Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2
            if Params.VERBOSE: print('il: ' + str(i) +
                                     '; P_MPa: ' + str(round(Planet.P_MPa[i],3)) +
                                     '; T_K: ' + str(round(Planet.T_K[i],3)) +
                                     '; phase: ' + str(Planet.phase[i]))

    if Params.VERBOSE: print('Convection calculations complete.')

    return Planet


def IceVUnderplate(Planet, Params):
    """ Conductive profile calculations for ice III and ice V layers between
        the ocean and surface ice I layer.

        Assigns Planet attributes:
            PbIII_MPa, PbV_MPa, Pb_MPa, P_MPa, T_K
    """
    Planet.PbIII_MPa = GetPfreezeHP(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, Planet.Bulk.TbIII_K, 3)
    Planet.Pb_MPa = GetPfreezeHP(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, Planet.Bulk.TbV_K, 5)
    # Set linear P and adiabatic T in ice III layers
    if(Planet.PbIII_MPa < Planet.P_MPa[Planet.Steps.nIbottom-1]):
        raise ValueError('Ice III bottom melting pressure is less than the pressure at the bottom of the ' +
                         'ice I shell. This means TbIII_K has been set too high to be consistent with Tb_K.')
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

    return Planet


def IceIIIUnderplate(Planet, Params):
    """ Conductive profile calculations for ice III layers between
        the ocean and surface ice I layer.

        Assigns Planet attributes:
            PbIII_MPa, Pb_MPa, P_MPa, T_K
    """

    Planet.Pb_MPa = GetPfreezeHP(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, Planet.Bulk.TbIII_K, 3)
    # Set linear P and adiabatic T in ice III layers
    PIceIII_MPa = np.linspace(Planet.P_MPa[Planet.Steps.nIbottom-1], Planet.PbIII_MPa, Planet.Steps.nIceIIILitho)
    TIceIII_K = Planet.Bulk.TbIII_K**(PIceIII_MPa/Planet.PbIII_MPa) * \
                Planet.T_K[Planet.Steps.nIbottom-1]**(1 - PIceIII_MPa/Planet.PbIII_MPa)
    Planet.P_MPa[Planet.Steps.nIbottom:Planet.Steps.nIIIbottom] = PIceIII_MPa
    Planet.T_K[Planet.Steps.nIbottom:Planet.Steps.nIIIbottom] = TIceIII_K
    return Planet


def OceanLayers(Planet, Params):
    """ Geophysical and thermodynamic calculations for ocean layer
        Calculates state variables of the layer with each pressure step

        Assigns Planet attributes:
            phase, r_m, z_m, g_ms2, T_K, P_MPa, rho_kgm3, Cp_JkgK, alpha_pK, MLayer_kg
    """
    if Params.VERBOSE: print('Evaluating ocean layers.')

    # Confirm that we haven't made mistakes in phase assignment in IceLayers()
    if Planet.phase[Planet.Steps.nSurfIce] != 0:
        raise ValueError('Phase of first "ocean" layer is not zero.')

    # Start ocean pressure a tiny bit above the melting pressure so that GetPhase doesn't get confused
    POcean_MPa = np.arange(Planet.Pb_MPa + 0.1*Planet.Ocean.deltaP, Planet.Ocean.PHydroMax_MPa, Planet.Ocean.deltaP)
    Planet.Steps.nOceanMax = np.size(POcean_MPa)

    # Initialize remaining local arrays
    TOcean_K, rhoOcean_kgm3, CpOcean_JkgK, alphaOcean_pK = (np.zeros(Planet.Steps.nOceanMax) for _ in range(4))
    TOcean_K = np.insert(TOcean_K, 0, Planet.T_K[Planet.Steps.nSurfIce-1])

    # Get ocean EOS functions
    TOceanLin_K = np.arange(Planet.Bulk.Tb_K, Planet.Ocean.THydroMax_K, Planet.Ocean.deltaT)
    Planet.Ocean.EOS = OceanEOSStruct(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, POcean_MPa, TOceanLin_K)

    for i in range(Planet.Steps.nOceanMax):
        Planet.phase[Planet.Steps.nSurfIce+i] = \
            GetPhase(Planet.Ocean.EOS, Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, POcean_MPa[i], TOcean_K[i])
        if Params.VERBOSE: print('il: '+str(Planet.Steps.nSurfIce+i) +
                               '; P_MPa: '+str(round(POcean_MPa[i],3)) +
                               '; T_K: '+str(round(TOcean_K[i],3)) +
                               '; phase: '+str(Planet.phase[Planet.Steps.nSurfIce+i]))
        if Planet.phase[Planet.Steps.nSurfIce+i] == 0:
            # Liquid water layers -- get fluid properties for the present layer but with the
            # overlaying layer's temperature
            rhoOcean_kgm3[i] = Planet.Ocean.EOS.fn_rho_kgm3(POcean_MPa[i], TOcean_K[i])
            CpOcean_JkgK[i] = Planet.Ocean.EOS.fn_Cp_JkgK(POcean_MPa[i], TOcean_K[i])
            alphaOcean_pK[i] = Planet.Ocean.EOS.fn_alpha_pK(POcean_MPa[i], TOcean_K[i])
            # Now use the present layer's properties to calculate an adiabatic thermal profile for layers below
            TOcean_K[i+1] = TOcean_K[i] + alphaOcean_pK[i] * TOcean_K[i] / \
                            CpOcean_JkgK[i] / rhoOcean_kgm3[i] * Planet.Ocean.deltaP*1e6
        else:
            # Undersea high-pressure ices -- we use GetTfreeze here to propagate the layer temperatures.
            # This is based on an assumption that the undersea HP ices are vigorously mixed by
            # two-phase convection, such that each layer is in local equilibrium with the liquid,
            # meaning each layer's temperature is equal to the melting temperature.
            rhoOcean_kgm3[i], CpOcean_JkgK[i], alphaOcean_pK[i] = \
                GetIceThermo([POcean_MPa[i]], [TOcean_K[i]], [Planet.phase[Planet.Steps.nSurfIce+i]])
            TOcean_K[i+1] = GetTfreeze(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, POcean_MPa[i], TOcean_K[i])

    TOcean_K = np.delete(TOcean_K, 0)

    Planet.P_MPa[Planet.Steps.nSurfIce:Planet.Steps.nSurfIce + Planet.Steps.nOceanMax] = POcean_MPa
    Planet.T_K[Planet.Steps.nSurfIce:Planet.Steps.nSurfIce + Planet.Steps.nOceanMax] = TOcean_K
    Planet.rho_kgm3[Planet.Steps.nSurfIce:Planet.Steps.nSurfIce + Planet.Steps.nOceanMax] = rhoOcean_kgm3
    Planet.Cp_JkgK[Planet.Steps.nSurfIce:Planet.Steps.nSurfIce + Planet.Steps.nOceanMax] = CpOcean_JkgK
    Planet.alpha_pK[Planet.Steps.nSurfIce:Planet.Steps.nSurfIce + Planet.Steps.nOceanMax] = alphaOcean_pK

    MAbove_kg = np.sum(Planet.MLayer_kg[:Planet.Steps.nSurfIce])
    for i in range(Planet.Steps.nSurfIce, Planet.Steps.nSurfIce + Planet.Steps.nOceanMax):
        Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1]) * 1e6 / Planet.g_ms2[i-1] / \
                        Planet.rho_kgm3[i-1]
        Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
        Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
        MAbove_kg += Planet.MLayer_kg[i-1]
        MBelow_kg = Planet.Bulk.M_kg - MAbove_kg
        Planet.g_ms2[i] = Constants.G * MBelow_kg / Planet.r_m[i]**2

    if Params.VERBOSE: print('Ocean layers complete.')

    return Planet


def InnerLayers(Planet, Params):
    """ Geophysical and thermodynamic calculations for silicate and core layers
        Calculates state variables of the layer with each pressure step

        Assigns Planet attributes:
            Steps.nTotal, all layer arrays
    """
    if Planet.Do.CONSTANT_INNER_DENSITY or Params.SKIP_INNER:
        # It will be better to calculate the layers always using the EOS, but no-core models
        # require more inputs to be physically reasonable.
        if not Planet.Do.CONSTANT_INNER_DENSITY:
            # Force this flag on in case Params.SKIP_INNER is True and this flag is not,
            # to avoid problems in SilicateLayers:
            Planet.Do.CONSTANT_INNER_DENSITY = True
            print('WARNING: Do.CONSTANT_INNER_DENSITY forced on based on implementation in SilicateLayers.')
        Planet, mantleProps, coreProps = CalcMoIConstantRho(Planet, Params)
    else:
        Planet, mantleProps, coreProps = CalcMoIWithEOS(Planet, Params)

    if Planet.Steps.nHydro <= Planet.Steps.nSurfIce: print('WARNING: For these run settings, the hydrosphere is entirely frozen.')
    Planet.Steps.nTotal = Planet.Steps.nHydro + Planet.Steps.nSil + Planet.Steps.nCore

    if Params.VERBOSE: print('Evaluating remaining quantities for layer arrays...')
    # Extend Planet layer arrays to make space for silicate and possible core layers
    extend = np.zeros(Planet.Steps.nSil + Planet.Steps.nCore)
    Planet.P_MPa = np.concatenate((Planet.P_MPa[:Planet.Steps.nHydro], extend))
    Planet.T_K = np.concatenate((Planet.T_K[:Planet.Steps.nHydro], extend))
    Planet.r_m = np.concatenate((Planet.r_m[:Planet.Steps.nHydro], extend))
    Planet.phase = np.concatenate((Planet.phase[:Planet.Steps.nHydro], extend.astype(np.int_)))
    Planet.rho_kgm3 = np.concatenate((Planet.rho_kgm3[:Planet.Steps.nHydro], extend))
    Planet.Cp_JkgK = np.concatenate((Planet.Cp_JkgK[:Planet.Steps.nHydro], extend))
    Planet.alpha_pK = np.concatenate((Planet.alpha_pK[:Planet.Steps.nHydro], extend))
    Planet.g_ms2 = np.concatenate((Planet.g_ms2[:Planet.Steps.nHydro], extend))
    Planet.phi_frac = np.concatenate((Planet.phi_frac[:Planet.Steps.nHydro], extend))
    Planet.sigma_Sm = np.concatenate((Planet.sigma_Sm[:Planet.Steps.nHydro], extend))
    Planet.z_m = np.concatenate((Planet.z_m[:Planet.Steps.nHydro], extend))

    # Assign phase values for silicates and core
    Planet.phase[Planet.Steps.nHydro:Planet.Steps.nHydro + Planet.Steps.nSil] = 50
    Planet.phase[Planet.Steps.nHydro + Planet.Steps.nSil:Planet.Steps.nTotal] = 100

    # Unpack results from MoI calculations
    iOS = Planet.Steps.nHydro
    iSC = Planet.Steps.nHydro + Planet.Steps.nSil
    Planet.P_MPa[iOS:iSC], Planet.T_K[iOS:iSC], Planet.r_m[iOS:iSC], Planet.rho_kgm3[iOS:iSC], \
    Planet.g_ms2[iOS:iSC], Planet.phi_frac[iOS:iSC] = mantleProps

    iCC = Planet.Steps.nTotal
    if Planet.Do.Fe_CORE:
        # Unpack results from MoI calculations
        Planet.P_MPa[iSC:iCC], Planet.T_K[iSC:iCC], Planet.r_m[iSC:iCC], Planet.rho_kgm3[iSC:iCC], \
        Planet.g_ms2[iSC:iCC], Planet.Cp_JkgK[iSC:iCC], Planet.alpha_pK[iSC:iCC] = coreProps

    Planet.z_m[iOS:iCC] = Planet.Bulk.R_m - Planet.r_m[iOS:iCC]

    # Check for any negative temperature gradient (indicates non-equilibrium conditions)
    gradTneg = np.where(np.diff(Planet.T_K) < 0)
    if np.any(gradTneg): print('WARNING: Negative temperature gradient starting at index ' + str(gradTneg[0][0])
                               + '. This indicates that internal heating parameters Qrad_Wkg and/or '
                               + 'Htidal_Wm3 are set too high to be consistent with the heat flux '
                               + 'into the ocean. The current configuration represents a '
                               + 'non-equilibrium state.')

    return Planet


def CalcMoIConstantRho(Planet, Params):
    """ Find the relative sizes of silicate, core, and hydrosphere layers that are
        consistent with the measured moment of inertia, based on calculated hydrosphere
        properties and assumptions about the silicate and possible core layers.

        Assigns Planet attributes:
            CMR2mean, Sil.RsilMean_m, Sil.RsilRange_m, Core.RFeMean_m, Core.RFeRange_m, Steps.nHydro
    """
    if Params.VERBOSE: print('Finding MoI consistent with measured value for constant-density inner layers...')
    # Get MR^2 -- we will need to divide each C by this later.
    MR2 = Planet.Bulk.M_kg * Planet.Bulk.R_m**2

    # Get final number of layers modeled in "overshoot" hydrosphere
    nHydroActual = Planet.Steps.nSurfIce + Planet.Steps.nOceanMax
    # Find contribution to axial moment of inertia C from each ocean layer
    dC_H2O = 8*np.pi/15 * Planet.rho_kgm3[:-1] * (Planet.r_m[:-1]**5 - Planet.r_m[1:]**5)
    # Find total mass contained above each hydrosphere layer
    MAbove_kg = np.array([np.sum(Planet.MLayer_kg[:i]) for i in range(nHydroActual)])
    # Find volume of a full sphere of silicate corresponding to each valid layer
    VsilSphere_m3 = 4/3*np.pi * Planet.r_m[Planet.Steps.iSilStart:]**3

    if Planet.Do.Fe_CORE:
        # Find core bulk density based on assumed sulfide content (Eq 10 of Vance et al., 2014)
        rhoCore_kgm3 = Planet.Core.rhoFeS_kgm3 * Planet.Core.rhoFe_kgm3 / \
            (Planet.Core.xFeS * (Planet.Core.rhoFe_kgm3 - Planet.Core.rhoFeS_kgm3) + Planet.Core.rhoFeS_kgm3)
        # Calculate core volume for a silicate layer with outer radius equal to bottom of each hydrosphere layer
        # and inner radius equal to the core radius
        VCore_m3 = np.array([(Planet.Bulk.M_kg - MAbove_kg[i] - VsilSphere_m3[i-Planet.Steps.iSilStart]*
                             Planet.Sil.rhoSilWithCore_kgm3) / (rhoCore_kgm3 - Planet.Sil.rhoSilWithCore_kgm3)
                             for i in range(Planet.Steps.iSilStart, nHydroActual-1)])
        # Find values for which the silicate radius is too large
        nTooBig = next((i[0] for i, val in np.ndenumerate(VCore_m3) if val>0))
        # Calculate corresponding core radii based on above density
        rCore_m = (VCore_m3[nTooBig:]*3/4/np.pi)**(1/3)
        # Assign fixed density to an array for dual-use code looking for compatible C/MR^2
        rhoSil_kgm3 = np.ones_like(rCore_m) * Planet.Sil.rhoSilWithCore_kgm3
    else:
        # Find silicate density consistent with observed bulk mass for each radius
        rhoSil_kgm3 = np.array([(Planet.Bulk.M_kg - MAbove_kg[i]) / VsilSphere_m3[i-Planet.Steps.iSilStart]
                              for i in range(Planet.Steps.iSilStart, nHydroActual-1)])
        # Density of silicates is scaled to fit the total mass, so there is no nTooBig in this case.
        nTooBig = 0
        # Set core radius and density to zero so calculations can proceed
        rCore_m = np.zeros(nHydroActual-1 - Planet.Steps.iSilStart)
        rhoCore_kgm3 = 0

    # Calculate C for a mantle extending up to each hydrosphere layer in turn
    C = np.zeros(nHydroActual - 1)
    C[Planet.Steps.iSilStart + nTooBig:] = [np.sum(dC_H2O[:i + Planet.Steps.iSilStart + nTooBig + 1]) +
            8*np.pi/15 * rhoSil_kgm3[i] * (Planet.r_m[i + Planet.Steps.iSilStart + nTooBig]**5 - rCore_m[i]**5) +
            8*np.pi/15 * rhoCore_kgm3 * rCore_m[i]**5
            for i in range(nHydroActual - Planet.Steps.iSilStart - nTooBig - 1)]
    CMR2 = C / MR2

    CMR2inds = [i[0] for i, valCMR2 in np.ndenumerate(CMR2)
                 if valCMR2 > Planet.Bulk.Cmeasured - Planet.Bulk.Cuncertainty
                and valCMR2 < Planet.Bulk.Cmeasured + Planet.Bulk.Cuncertainty]

    if len(CMR2inds) == 0:
        raise ValueError('No MoI found matching ' +
                        ' C/MR^2 = ' + str(round(Planet.Bulk.Cmeasured, 3)) + '±' + str(round(Planet.Bulk.Cuncertainty,3)) + '.\n ' +
                         'Min: ' + str(round(min(CMR2[CMR2>0]),3)) + ', Max: ' + str(round(max(CMR2),3)) + '.\n ' +
                         'Try increasing PHydroMax_MPa or adjusting properties of silicates and core.')

    # Find the C/MR^2 value most closely matching the measured value
    CMR2diff = np.abs(CMR2[CMR2inds] - Planet.Bulk.Cmeasured)
    # Get index of closest match in CMR2inds
    iCMR2inds = np.argmin(CMR2diff)
    # Find Planet array index corresponding to closest matching value
    iCMR2 = CMR2inds[iCMR2inds]
    iCMR2inner = iCMR2 - Planet.Steps.iSilStart - nTooBig
    CMR2indsInner = [ind - Planet.Steps.iSilStart - nTooBig for ind in CMR2inds]
    # Record the best-match C/MR^2 value
    Planet.CMR2mean = CMR2[iCMR2]
    # Record interior sizes
    Planet.Sil.rhoTrade_kgm3 = rhoSil_kgm3[CMR2indsInner]
    Planet.Sil.Rmean_m = Planet.r_m[iCMR2]
    Planet.Sil.Rtrade_m = Planet.r_m[CMR2inds]
    Planet.Sil.Rrange_m = Planet.Sil.Rtrade_m[0] - Planet.Sil.Rtrade_m[-1]
    Planet.Core.Rmean_m = rCore_m[iCMR2inner]
    Planet.Core.Rtrade_m = rCore_m[CMR2indsInner]
    Planet.Core.Rrange_m = Planet.Core.Rtrade_m[-1] - Planet.Core.Rtrade_m[0]
    # Now we finally know how many layers there are in the hydrosphere
    Planet.Steps.nHydro = iCMR2 - 1
    # Number of steps in the silicate layer is fixed for the constant-density approach
    Planet.Steps.nSil = Planet.Steps.nSilMax

    if not Params.SKIP_INNER:
        # Load in Perple_X table for silicate properties
        if Params.VERBOSE: print('Loading silicate Perple_X table...')
        Planet.Sil.EOS = PerplexEOSStruct(Planet.Sil.mantleEOS, EOSinterpMethod=Params.interpMethod)
        # Propagate the silicate EOS from each hydrosphere layer to the center of the body
        nSilTooBig, nProfiles, Psil_MPa, Tsil_K, rSil_m, rhoSilEOS_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, phiSil_frac \
            = SilicateLayers(Planet, Params)

        # Fill core/mantle trade arrays and set mean values consistent with MoI
        MtotSil_kg = np.sum(MLayerSil_kg)
        Planet.Sil.rhoMean_kgm3 = MtotSil_kg / (4/3*np.pi * (rSil_m[0,0]**3 - rSil_m[0,-1]**3))

        if Planet.Do.Fe_CORE:
            # Load in Perple_X table for core properties
            if Params.VERBOSE: print('Loading core Perple_X table...')
            Planet.Core.EOS = PerplexEOSStruct(Planet.Core.coreEOS, EOSinterpMethod=Params.interpMethod)
            # Propagate the silicate EOS from each hydrosphere layer to the center of the body
            _, Pcore_MPa, Tcore_K, rCoreEOS_m, rhoCore_kgm3, MLayerCore_kg, gCore_ms2, CpCore_JkgK, alphaCore_pK = \
                IronCoreLayers(Planet, Params,
                               nSilTooBig, nProfiles, Psil_MPa, Tsil_K, rSil_m, MAboveSil_kg, gSil_ms2)

            MtotCore_kg = np.sum(MLayerCore_kg)
            Planet.Core.rhoMean_kgm3 = MtotCore_kg / VCore_m3[iCMR2inner]

            coreProps = (Pcore_MPa, Tcore_K, rCoreEOS_m[0,:-1], rhoCore_kgm3, gCore_ms2, CpCore_JkgK, alphaCore_pK)
        else:
            MtotCore_kg = 0
            Planet.Core.rhoMean_kgm3 = 0
            Planet.Core.Rtrade_m = np.zeros_like(Planet.Sil.Rtrade_m)
            Planet.Core.Rrange_m = 0
            coreProps = ()

        if Params.VERBOSE:
            Mtot_kg = np.sum(Planet.MLayer_kg[:iCMR2]) + MtotSil_kg + MtotCore_kg
            print('Found matching MoI of ' + str(round(Planet.CMR2mean,3)) +
                ' (C/MR^2 = ' + str(round(Planet.Bulk.Cmeasured,3)) + '±' + str(round(Planet.Bulk.Cuncertainty,3)) + ') for ' +
                  'rho_sil = ' + str(round(Planet.Sil.rhoMean_kgm3)) + ' kg/m^3, ' +
                  'R_sil = ' + str(round(Planet.Sil.Rmean_m / Planet.Bulk.R_m,2)) + ' R, ' +
                  'R_core = ' + str(round(Planet.Core.Rmean_m / Planet.Bulk.R_m,2)) + ' R, ' +
                  'M_tot = ' + str(round(Mtot_kg/Planet.Bulk.M_kg,4)) + ' M_' + Planet.name[0] + '.')
            print('WARNING: Because silicate and core properties were determined from the EOS after finding their ' +
                  'sizes by assuming constant densities, the body mass may not match the measured value.')

    else:
        Psil_MPa, Tsil_K, rhoSilEOS_kgm3, gSil_ms2, phiSil_frac = (np.zeros(Planet.Steps.nSil) for _ in range(5))
        rSil_m = np.zeros((1, Planet.Steps.nSil+1))
        coreProps = (np.zeros(Planet.Steps.nCore) for _ in range(7))
        Planet.Sil.rhoMean_kgm3 = 0
        Planet.Core.rhoMean_kgm3 = 0

    mantleProps = (Psil_MPa, Tsil_K, rSil_m[0,:-1], rhoSilEOS_kgm3, gSil_ms2, phiSil_frac)

    return Planet, mantleProps, coreProps


def CalcMoIWithEOS(Planet, Params):
    """ Find the relative sizes of silicate, core, and hydrosphere layers that are
        consistent with the measured moment of inertia, based on calculated hydrosphere
        properties and EOS data for assumed mantle and core compositions output by Perple_X.

        Assigns Planet attributes:
            CMR2mean, Sil.RsilMean_m, Sil.RsilRange_m, Core.RFeMean_m, Core.RFeRange_m, Steps.nHydro, Steps.nSil,
            all layer arrays
    """
    if Params.VERBOSE: print('Finding MoI consistent with measured value...')
    # Get MR^2 -- we will need to divide each C by this later.
    MR2 = Planet.Bulk.M_kg * Planet.Bulk.R_m**2

    # Load in Perple_X table for silicate properties
    if Params.VERBOSE: print('Loading silicate Perple_X table...')
    Planet.Sil.EOS = PerplexEOSStruct(Planet.Sil.mantleEOS, EOSinterpMethod=Params.interpMethod)
    # Propagate the silicate EOS from each hydrosphere layer to the center of the body
    nSilTooBig, nProfiles, Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, phiSil_frac \
        = SilicateLayers(Planet, Params)

    # Find contribution to axial moment of inertia C from each ocean layer
    dCfromH2O = 8*np.pi/15 * Planet.rho_kgm3[:-1] * (Planet.r_m[:-1]**5 - Planet.r_m[1:]**5)
    # Same for silicate layers
    dCfromSil = 8*np.pi/15 * rhoSil_kgm3 * (rSil_m[:,:-1]**5 - rSil_m[:,1:]**5)

    if Planet.Do.Fe_CORE:
        # Load in Perple_X table for core properties
        if Params.VERBOSE: print('Loading core Perple_X table...')
        Planet.Core.EOS = PerplexEOSStruct(Planet.Core.coreEOS, EOSinterpMethod=Params.interpMethod)
        # Propagate the silicate EOS from each hydrosphere layer to the center of the body
        nSilFinal, Pcore_MPa, Tcore_K, rCore_m, rhoCore_kgm3, MLayerCore_kg, gCore_ms2, CpCore_JkgK, alphaCore_pK = \
            IronCoreLayers(Planet, Params,
                           nSilTooBig, nProfiles, Psil_MPa, Tsil_K, rSil_m, MAboveSil_kg, gSil_ms2)

        dCfromCore = 8*np.pi/15 * rhoCore_kgm3 * (rCore_m[:,:-1]**5 - rCore_m[:,1:]**5)
    else:
        # To implement silicates-only here, we need some other variable besides silicate size to tweak
        # that can permit us to match both the body mass and MoI. Composition or porosity would each work,
        # but there's no example to draw from for tweaking these to match the MoI yet.
        coreProps = None
        nSilFinal = Planet.Steps.nSilMax
        raise ValueError('Handling for MoI matching with silicate EOS and no core is not implemented yet. Set Planet.Do.CONSTANT_INNER_DENSITY to True and run again.')

    # Calculate C for a mantle extending up to each hydrosphere layer in turn
    C = np.zeros(nProfiles + Planet.Steps.iSilStart)
    C[Planet.Steps.iSilStart + nSilTooBig:] = [np.sum(dCfromH2O[:i + Planet.Steps.iSilStart + 1]) + \
            np.sum(dCfromSil[i,:nSilFinal[i]]) + \
            np.sum(dCfromCore[i - nSilTooBig,:]) \
            for i in range(nSilTooBig, nProfiles)]
    CMR2 = C / MR2

    CMR2inds = [i[0] for i, valCMR2 in np.ndenumerate(CMR2)
                 if valCMR2 > Planet.Bulk.Cmeasured - Planet.Bulk.Cuncertainty
                and valCMR2 < Planet.Bulk.Cmeasured + Planet.Bulk.Cuncertainty]

    if len(CMR2inds) == 0:
        raise ValueError('No MoI found matching ' +
                        ' C/MR^2 = ' + str(round(Planet.Bulk.Cmeasured, 3)) + '±' + str(round(Planet.Bulk.Cuncertainty,3)) + '.\n ' +
                         'Min: ' + str(round(min(CMR2[CMR2>0]),3)) + ', Max: ' + str(round(max(CMR2),3)) + '.\n ' +
                         'Try increasing PHydroMax_MPa or adjusting properties of silicates and core.')

    # Find the C/MR^2 value most closely matching the measured value
    CMR2diff = np.abs(CMR2[CMR2inds] - Planet.Bulk.Cmeasured)
    # Get index of closest match in CMR2inds
    iCMR2inds = np.argmin(CMR2diff)
    # Find Planet array index corresponding to closest matching value
    iCMR2 = CMR2inds[iCMR2inds]
    # Get indices for inner layer arrays
    iCMR2sil = iCMR2 - Planet.Steps.iSilStart
    iCMR2core = iCMR2sil - nSilTooBig
    CMR2indsSil = [ind - Planet.Steps.iSilStart for ind in CMR2inds]
    CMR2indsCore = [ind - nSilTooBig for ind in CMR2indsSil]
    # Record the best-match C/MR^2 value
    Planet.CMR2mean = CMR2[iCMR2]
    # Now we finally know how many layers there are in the hydrosphere and silicates
    Planet.Steps.nHydro = iCMR2 - 1
    Planet.Steps.nSil = nSilFinal[iCMR2sil]

    # Fill core/mantle trade arrays and set mean values consistent with MoI
    MtotSil_kg = np.sum(MLayerSil_kg[iCMR2sil,:nSilFinal[iCMR2sil]])
    Planet.Sil.rhoMean_kgm3 = MtotSil_kg / (4/3*np.pi * (rSil_m[iCMR2sil,0]**3 - rSil_m[iCMR2sil,nSilFinal[iCMR2sil]-1]**3))
    Planet.Sil.rhoTrade_kgm3 = np.array([np.sum(MLayerSil_kg[i,:nSilFinal[i]]) / (4/3*np.pi * (rSil_m[i,0]**3 - rSil_m[i,nSilFinal[i]-1]**3)) for i in CMR2indsSil])
    Planet.Sil.Rmean_m = Planet.r_m[iCMR2]
    Planet.Sil.Rtrade_m = Planet.r_m[CMR2inds]
    Planet.Sil.Rrange_m = Planet.Sil.Rtrade_m[0] - Planet.Sil.Rtrade_m[-1]
    if Planet.Do.Fe_CORE:
        MtotCore_kg = np.sum(MLayerCore_kg[iCMR2core,:])
        Planet.Core.rhoMean_kgm3 = MtotCore_kg / (4/3*np.pi * rCore_m[iCMR2core,0]**3)
        Planet.Core.Rmean_m = rCore_m[iCMR2core,0]
        Planet.Core.Rtrade_m = rCore_m[CMR2indsCore,0]
        Planet.Core.Rrange_m = Planet.Core.Rtrade_m[-1] - Planet.Core.Rtrade_m[0]

        # Package up core properties for returning
        coreProps = (Pcore_MPa[iCMR2core,:], Tcore_K[iCMR2core,:], rCore_m[iCMR2core,:-1],
                     rhoCore_kgm3[iCMR2core,:], gCore_ms2[iCMR2core,:], CpCore_JkgK[iCMR2core,:],
                     alphaCore_pK[iCMR2core,:])
    else:
        MtotCore_kg = 0
        Planet.Core.rhoMean_kgm3 = 0
        Planet.Core.Rmean_m = 0
        Planet.Core.Rtrade_m = 0
        Planet.Core.Rrange_m = 0
        coreProps = ()

    if Params.VERBOSE:
        Mtot_kg = np.sum(Planet.MLayer_kg[:iCMR2]) + MtotSil_kg + MtotCore_kg
        print('Found matching MoI of ' + str(round(Planet.CMR2mean,3)) +
            ' (C/MR^2 = ' + str(round(Planet.Bulk.Cmeasured,3)) + '±' + str(round(Planet.Bulk.Cuncertainty,3)) + ') for ' +
              'rho_sil = ' + str(round(Planet.Sil.rhoMean_kgm3)) + ' kg/m^3, ' +
              'R_sil = ' + str(round(Planet.Sil.Rmean_m / Planet.Bulk.R_m,2)) + ' R, ' +
              'R_core = ' + str(round(Planet.Core.Rmean_m / Planet.Bulk.R_m,2)) + ' R, ' +
              'M_tot = ' + str(round(Mtot_kg/Planet.Bulk.M_kg,4)) + ' M_' + Planet.name[0] + '.')

    mantleProps = (Psil_MPa[iCMR2sil,:nSilFinal[iCMR2sil]], Tsil_K[iCMR2sil,:nSilFinal[iCMR2sil]],
                   rSil_m[iCMR2sil,:nSilFinal[iCMR2sil]], rhoSil_kgm3[iCMR2sil,:nSilFinal[iCMR2sil]],
                   gSil_ms2[iCMR2sil,:nSilFinal[iCMR2sil]], phiSil_frac[iCMR2sil,:nSilFinal[iCMR2sil]])

    return Planet, mantleProps, coreProps


def SilicateLayers(Planet, Params):
    """ Determines properties of silicate layers based on input Perple_X table
        and seafloor properties.

        Returns:
            nSilTooBig (int): Number of silicate profiles that have a mass that exceeds the body mass
            nProfiles (int): Number of silicate profiles considered
            Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, phiSil_frac
                (float, shape nHydroMax-2): State variables needed to determine underlying core properties
                 to proceed with MoI calculations.
    """
    if Planet.Do.CONSTANT_INNER_DENSITY:
        nProfiles = 1
        profRange = [Planet.Steps.nHydro - Planet.Steps.iSilStart + 1]
    else:
        nProfiles = Planet.Steps.nSurfIce + Planet.Steps.nOceanMax - Planet.Steps.iSilStart - 1
        profRange = range(nProfiles)
    # Initialize output arrays and working arrays
    Psil_MPa, Tsil_K, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg, gSil_ms2, phiSil_frac, \
        = (np.zeros((nProfiles, Planet.Steps.nSilMax)) for _ in range(8))
    # Check if we set the core radius to 0 or a found C/MR^2 value (for constant-density approach)
    if Planet.Core.Rmean_m is not None:
        rSilEnd_m = Planet.Core.Rmean_m
    else:
        rSilEnd_m = 0

    rSil_m = np.array([np.linspace(Planet.r_m[i+Planet.Steps.iSilStart], rSilEnd_m, Planet.Steps.nSilMax+1) for i in profRange])
    Psil_MPa[:,0] = [Planet.P_MPa[i+Planet.Steps.iSilStart] for i in profRange]
    Tsil_K[:,0] = [Planet.T_K[i+Planet.Steps.iSilStart] for i in profRange]
    rhoSil_kgm3[:,0] = [Planet.Sil.EOS.fn_rho_kgm3(Psil_MPa[i,0], Tsil_K[i,0]) for i in range(nProfiles)]
    kThermSil_WmK[:,0] = [Planet.Sil.kTherm_WmK for i in range(nProfiles)]
    MLayerSil_kg[:,0] = [rhoSil_kgm3[i,0] * 4/3*np.pi*(rSil_m[i,0]**3 - rSil_m[i,1]**3) for i in range(nProfiles)]
    gSil_ms2[:,0] = [Planet.g_ms2[i+Planet.Steps.iSilStart] for i in profRange]
    if Planet.Do.POROUS_ROCK:
        print('POROUS_ROCK not implemented yet. Only the top silicate layer will have porosity set.')
        phiSil_frac[:,0] = [Planet.Sil.phiRockMax_frac for _ in profRange]

    MHydro_kg = np.array([np.sum(Planet.MLayer_kg[:i]) for i in range(Planet.Steps.iSilStart, Planet.Steps.iSilStart + nProfiles)])
    # Initialize MAbove_kg to 0th silicate layer, so that the hydrosphere mass is equal to the mass above the silicates.
    MAboveSil_kg[:,0] = MHydro_kg + 0.0
    # Initialize Qtop_WmK, the heat flux leaving the top of each layer.
    # For the top layer, this is equal to the heat flux warming the hydrosphere.
    qTop_WmK = Planet.Ocean.QfromMantle_W / 4/np.pi/rSil_m[:,0]**2

    if Params.VERBOSE: print('Propagating silicate EOS for each possible mantle size...')
    for j in range(1, Planet.Steps.nSilMax):
        MAboveSil_kg[:,j] = MAboveSil_kg[:,j-1] + MLayerSil_kg[:,j-1]
        Psil_MPa[:,j] = Psil_MPa[:,j-1] + 1e-6 * MLayerSil_kg[:,j-1] * gSil_ms2[:,j-1] / rSil_m[:,j]**2
        Tsil_K[:,j], qTop_WmK = ConductiveTemperature(Tsil_K[:,j-1], rSil_m[:,j-1], rSil_m[:,j],
                    kThermSil_WmK[:,j-1], rhoSil_kgm3[:,j-1], Planet.Sil.Qrad_Wkg, Planet.Sil.Htidal_Wm3,
                    qTop_WmK)
        rhoSil_kgm3[:,j] = [Planet.Sil.EOS.fn_rho_kgm3(Psil_MPa[i,j], Tsil_K[i,j]) for i in range(nProfiles)]
        kThermSil_WmK[:,j] = kThermSil_WmK[:,j-1]  # Placeholder until we implement variable thermal conductivity
        MLayerSil_kg[:,j] = rhoSil_kgm3[:,j] * 4/3*np.pi*(rSil_m[:,j]**3 - rSil_m[:,j+1]**3)
        # Calculate gravity using absolute values, as we will use MAboveSil to check for exceeding body mass later.
        gSil_ms2[:,j] = Constants.G * np.abs(Planet.Bulk.M_kg - MAboveSil_kg[:,j]) / rSil_m[:,j]**2

    if Planet.Do.CONSTANT_INNER_DENSITY:
        # Set to zero for later calculations if we already found the desired C/MR^2 match
        nSilTooBig = 0
    else:
        # Get total mass for each possible silicate layer size
        Mtot_kg = MLayerSil_kg[:,-1] + MAboveSil_kg[:,-1]
        # Find silicate radii for which the total mass is too high so we can exclude them
        try:
            nSilTooBig = next(i[0] for i, val in np.ndenumerate(Mtot_kg) if val <= Planet.Bulk.M_kg)
        except:
            raise RuntimeError('No core/mantle combination was found matching total mass and moment of inertia. ' +
                               'Minimum mass: ' + str(round(np.min(Mtot_kg/Planet.Bulk.M_kg),3)) + ' M_' + Planet.name[0] +
                               '. Try adjusting run settings that affect mantle density, like silicate composition ' +
                               'and heat flux settings.')

    return nSilTooBig, nProfiles, Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, \
           MLayerSil_kg, MAboveSil_kg, gSil_ms2, phiSil_frac


def IronCoreLayers(Planet, Params,
                   nSilTooBig, nProfiles, Psil_MPa, Tsil_K, rSil_m, MAboveSil_kg, gSil_ms2):
    """ Determines properties of core layers based on input Perple_X table
        and seafloor properties.

        Args:
            nSilTooBig (int): Number of silicate profiles to skip past due to masses exceeding body mass.
            Psil_MPa, Tsil_K, rSil_m, MLayerSil_kg, MAboveSil_kg, gSil_ms2 (float, shape NxM): Outputs from
                SilicateLayers with layer properties for each silicate region size possibility.
        Returns:
            nSilFinal (int): Index in silicate profiles of core with a total mass just under body mass.
            Pcore_MPa, Tcore_K, rCore_m, rhoCore_kgm3, MLayerCore_kg, gCore_ms2 (float, shape Planet.Steps.nCore):
                Core properties needed to determine MoI.
    """
    # Initialize output arrays and working arrays
    Pcore_MPa, Tcore_K, rhoCore_kgm3, MLayerCore_kg, gCore_ms2, CpCore_JkgK, alphaCore_pK = \
        (np.zeros((nProfiles-nSilTooBig, Planet.Steps.nCore)) for _ in range(7))
    rCore_m = np.zeros((nProfiles-nSilTooBig, Planet.Steps.nCore+1))
    # Initialize matching indices as -1 as a flag for unfilled values
    iCoreMatch, nSilFinal = (-1 * np.ones(nProfiles).astype(np.int_) for _ in range(2))

    # Calculate maximum core size based on minimum plausible density setting
    MCore_kg = Planet.Bulk.M_kg - MAboveSil_kg
    MCore_kg[:nSilTooBig,:] = 0

    if Planet.Do.CONSTANT_INNER_DENSITY:
        iCoreStart = [-1]
        silEnd = 0
    else:
        rCoreMax_m = (MCore_kg/Planet.Core.rhoMin_kgm3 * 3/4/np.pi)**(1/3)
        # Find first silicate layer smaller than the max core radius
        iCoreStart = [next(j[0] for j,val in np.ndenumerate(rSil_m[i,:]) if val < rCoreMax_m[i,j])
                      for i in range(nSilTooBig, nProfiles)]
        silEnd = Planet.Steps.nSilMax

    if Params.VERBOSE: print('Evaluating core EOS for possible configurations...')
    for iProf in range(nSilTooBig, nProfiles):
        # Get index for profile starting at 0 for first valid one (past nSilTooBig)
        iValid = iProf - nSilTooBig
        # Get index for which silicate layer gets replaced by core layers there and below
        thisCoreStart = iCoreStart[iValid]
        # Get number of remaining silicate layers to iterate over for possible core configs
        nSilRemain = silEnd - thisCoreStart
        #(Re-)initialize placeholder arrays for each core possibility (they change length for each)
        thisPcore_MPa, thisTcore_K, thisrhoCore_kgm3, thisMLayerCore_kg, thisgCore_ms2, thisCpCore_JkgK, \
        thisalphaCore_pK = (np.zeros((nSilRemain, Planet.Steps.nCore)) for _ in range(7))

        # Set starting core values for all possibilities to be equal to silicates at this transition radius
        thisrCore_m = np.array([np.linspace(rSil_m[iProf,thisCoreStart+j], 0, Planet.Steps.nCore+1) for j in range(nSilRemain)])
        thisPcore_MPa[:,0] = [Psil_MPa[iProf,thisCoreStart+j] for j in range(nSilRemain)]
        thisTcore_K[:,0] = [Tsil_K[iProf,thisCoreStart+j] for j in range(nSilRemain)]
        thisrhoCore_kgm3[:,0] = [Planet.Core.EOS.fn_rho_kgm3(thisPcore_MPa[j,0], thisTcore_K[j,0]) for j in range(nSilRemain)]
        thisCpCore_JkgK[:,0] = [Planet.Core.EOS.fn_Cp_JkgK(thisPcore_MPa[j,0], thisTcore_K[j,0]) for j in range(nSilRemain)]
        thisalphaCore_pK[:,0] = [Planet.Core.EOS.fn_alpha_pK(thisPcore_MPa[j,0], thisTcore_K[j,0]) for j in range(nSilRemain)]
        thisMLayerCore_kg[:,0] = [thisrhoCore_kgm3[j,0] * 4/3*np.pi*(thisrCore_m[j,0]**3 - thisrCore_m[j,1]**3) for j in range(nSilRemain)]
        thisgCore_ms2[:,0] = [gSil_ms2[iProf,thisCoreStart+j] for j in range(nSilRemain)]
        MAbove_kg = np.array([MAboveSil_kg[iProf,thisCoreStart+j] for j in range(nSilRemain)])

        for k in range(1, Planet.Steps.nCore):
            MAbove_kg += thisMLayerCore_kg[:,k-1]
            thisDeltaP = 1e-6 * thisMLayerCore_kg[:,k-1] * thisgCore_ms2[:,k-1] / thisrCore_m[:,k]**2
            thisPcore_MPa[:,k] = thisPcore_MPa[:,k-1] + thisDeltaP
            thisTcore_K[:,k] = thisTcore_K[:,k-1] + thisalphaCore_pK[:,k-1]*thisTcore_K[:,k] / \
                           thisCpCore_JkgK[:,k-1] / thisrhoCore_kgm3[:,k-1] * thisDeltaP*1e6
            thisrhoCore_kgm3[:,k] = [Planet.Core.EOS.fn_rho_kgm3(thisPcore_MPa[i,k], thisTcore_K[i,k]) for i in range(nSilRemain)]
            thisCpCore_JkgK[:,k] = [Planet.Core.EOS.fn_Cp_JkgK(thisPcore_MPa[i,k], thisTcore_K[i,k]) for i in range(nSilRemain)]
            thisalphaCore_pK[:,k] = [Planet.Core.EOS.fn_alpha_pK(thisPcore_MPa[i,k], thisTcore_K[i,k]) for i in range(nSilRemain)]
            thisMLayerCore_kg[:,k] = thisrhoCore_kgm3[:,k] * 4/3*np.pi*(thisrCore_m[:,k]**3 - thisrCore_m[:,k+1]**3)
            # Approximate gravity as linear to avoid blowing up for total mass less than body mass (accurate for constant density only)
            thisgCore_ms2[:,k] = thisgCore_ms2[:,0] * thisrCore_m[:,k] / thisrCore_m[:,0]

        if not Planet.Do.CONSTANT_INNER_DENSITY:
            # Find the first core profile that has a mass just below the body mass
            Mtot_kg = MAbove_kg + thisMLayerCore_kg[:,-1]
            iCoreMatch[iProf] = next(ii[0] for ii,val in np.ndenumerate(Mtot_kg) if val < Planet.Bulk.M_kg)
            nSilFinal[iProf] = iCoreStart[iValid] + iCoreMatch[iProf]
            if Params.VERBOSE: print('Core match for iProf = ' + str(iProf) +
                                     ' with Steps.nSil = ' + str(nSilFinal[iProf]) +
                                     ' and M = ' + str(round(Mtot_kg[iCoreMatch[iProf]]/Planet.Bulk.M_kg,4)) + ' M_P.')

        # Assign the values for the core profile with matching total mass to output arrays
        Pcore_MPa[iValid,:] = thisPcore_MPa[iCoreMatch[iProf],:]
        Tcore_K[iValid,:] = thisTcore_K[iCoreMatch[iProf],:]
        rCore_m[iValid,:] = thisrCore_m[iCoreMatch[iProf],:]
        rhoCore_kgm3[iValid,:] = thisrhoCore_kgm3[iCoreMatch[iProf],:]
        MLayerCore_kg[iValid,:] = thisMLayerCore_kg[iCoreMatch[iProf],:]
        gCore_ms2[iValid,:] = thisgCore_ms2[iCoreMatch[iProf],:]
        CpCore_JkgK[iValid,:] = thisCpCore_JkgK[iCoreMatch[iProf],:]
        alphaCore_pK[iValid,:] = thisalphaCore_pK[iCoreMatch[iProf],:]

    return nSilFinal, Pcore_MPa, Tcore_K, rCore_m, rhoCore_kgm3, MLayerCore_kg, gCore_ms2, CpCore_JkgK, alphaCore_pK
