import numpy as np
from collections.abc import Iterable
from Utilities.dataStructs import Constants
from Thermodynamics.HydroEOS import GetIceThermo, GetPfreeze, GetTfreeze, GetPfreezeHP, FluidEOS, GetPhase
from Thermodynamics.FromLiterature.IceShell import ConductionClathLid, ConvectionDeschampsSotin2001
from Thermodynamics.InnerEOS import PerplexEOSStruct, MantleEOS

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

    # Surface ice layer propagators -- if present, clathrates are nearest the surface
    # so we do clathrate layers first.
    if Planet.Do.CLATHRATE:
        if Params.VERBOSE: print('Evaluating clathrate layers.')
        ClathrateLayers(Planet)
    else:
        Planet.PbClath_MPa = Planet.Bulk.Psurf_MPa

    # Get the pressure consistent with the bottom of the ice Ih layer that is
    # consistent with the choice of Tb_K we suppose for this model
    if Params.VERBOSE: print('Finding uppermost ice melting pressure...')
    Planet.PbI_MPa = GetPfreeze(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, Planet.Bulk.Tb_K,
        PfreezeLower_MPa=Planet.PfreezeLower_MPa, PfreezeUpper_MPa=Planet.PfreezeUpper_MPa, PfreezeRes_MPa=Planet.PfreezeRes_MPa,
        Pguess=None)
    # Set linear P and adiabatic T in ice I layers
    PIceI_MPa = np.linspace(Planet.PbClath_MPa, Planet.PbI_MPa, Planet.Steps.nIceI)
    TIceI_K = Planet.Bulk.Tb_K**(PIceI_MPa/Planet.PbI_MPa) * Planet.T_K[Planet.Steps.nClath]**(1 - PIceI_MPa/Planet.PbI_MPa)
    Planet.P_MPa[Planet.Steps.nClath:Planet.Steps.nIbottom] = PIceI_MPa
    Planet.T_K[Planet.Steps.nClath:Planet.Steps.nIbottom] = TIceI_K

    # Additional adiabats in ice III and/or V underplate layers -- for thick, cold ice shells
    if Planet.Do.BOTTOM_ICEV:
        if Params.VERBOSE: print('...with ice III and V underplating...')
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
        if Params.VERBOSE: print('...with ice III underplating...')
        Planet.Pb_MPa = GetPfreezeHP(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, Planet.Bulk.TbIII_K, 3)
        # Set linear P and adiabatic T in ice III layers
        PIceIII_MPa = np.linspace(Planet.P_MPa[Planet.Steps.nIbottom-1], Planet.PbIII_MPa, Planet.Steps.nIceIIILitho)
        TIceIII_K = Planet.Bulk.TbIII_K**(PIceIII_MPa/Planet.PbIII_MPa) * \
                    Planet.T_K[Planet.Steps.nIbottom-1]**(1 - PIceIII_MPa/Planet.PbIII_MPa)
        Planet.P_MPa[Planet.Steps.nIbottom:Planet.Steps.nIIIbottom] = PIceIII_MPa
        Planet.T_K[Planet.Steps.nIbottom:Planet.Steps.nIIIbottom] = TIceIII_K
    else:
        Planet.Pb_MPa = Planet.PbI_MPa

    if Params.VERBOSE: print('Upper ice melting pressure: '+str(Planet.Pb_MPa))

    # Evaluate thermodynamic properties of uppermost ice with SeaFreeze
    rhoSurfIce_kgm3, CpSurfIce_JkgK, alphaSurfIce_pK = GetIceThermo(Planet.P_MPa[:Planet.Steps.nSurfIce],
        Planet.T_K[:Planet.Steps.nSurfIce], Planet.phase[:Planet.Steps.nSurfIce])
    Planet.rho_kgm3[:Planet.Steps.nSurfIce] = rhoSurfIce_kgm3
    Planet.Cp_JkgK[:Planet.Steps.nSurfIce] = CpSurfIce_JkgK
    Planet.alpha_pK[:Planet.Steps.nSurfIce] = alphaSurfIce_pK

    # Calculate remaining physical properties of uppermost ice
    thisMAbove_kg = 0
    for i in range(1, Planet.Steps.nSurfIce):
        # Increment depth based on change in pressure, combined with gravity and density
        Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
        Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
        Planet.MLayer_kg[i-1] = 4/3*np.pi * rhoSurfIce_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
        thisMAbove_kg += Planet.MLayer_kg[i-1]
        thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
        Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2
        if Params.VERBOSE: print('il: ' + str(i) +
                                 '; P_MPa: ' + str(Planet.P_MPa[i]) +
                                 '; T_K: ' + str(Planet.T_K[i]) +
                                 '; phase: ' + str(Planet.phase[i]))
    Planet.zb_km = Planet.z_m[Planet.Steps.nSurfIce-1] / 1e3

    if Params.VERBOSE: print('Upper ice complete.')

    return Planet


def ClathrateLayers(Planet, Params):
    """ For ice shells insulated by a layer of clathrate at the surface
        Calculates state variables of the layer with each pressure step

        Requires Planet attributes:
            Planet.Steps.nClath
        Assigns Planet attributes:
            phase
    """

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

    POcean_MPa = np.arange(Planet.Pb_MPa+Planet.Ocean.deltaP, Planet.Ocean.PHydroMax_MPa, Planet.Ocean.deltaP)
    Planet.Steps.nOceanMax = np.size(POcean_MPa)

    # Initialize remaining local arrays
    TOcean_K, rhoOcean_kgm3, CpOcean_JkgK, alphaOcean_pK = (np.zeros(Planet.Steps.nOceanMax) for _ in range(4))
    TOcean_K = np.insert(TOcean_K, 0, Planet.T_K[Planet.Steps.nSurfIce-1])

    for i in range(Planet.Steps.nOceanMax):
        Planet.phase[Planet.Steps.nSurfIce+i] = \
            GetPhase(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, POcean_MPa[i], TOcean_K[i])
        if Params.VERBOSE: print('il: '+str(Planet.Steps.nSurfIce+i) +
                               '; P_MPa: '+str(POcean_MPa[i]) +
                               '; T_K: '+str(TOcean_K[i]) +
                               '; phase: '+str(Planet.phase[Planet.Steps.nSurfIce+i]))
        if Planet.phase[Planet.Steps.nSurfIce+i] == 0:
            # Liquid water layers -- get fluid properties for the present layer but with the
            # overlaying layer's temperature
            rhoOcean_kgm3[i], CpOcean_JkgK[i], alphaOcean_pK[i] = \
                FluidEOS(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, [POcean_MPa[i]], [TOcean_K[i]])
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


def HydroConvect(Planet, Params):
    """ Apply convection models from literature to adjust thermal profile and
        state variables for hydrosphere layers

        Assigns Planet attributes:
            ???
    """
    if Planet.Do.CLATHRATE:
        if Params.VERBOSE: print('Applying clathrate lid conduction.')
        Planet = ConductionClathLid()

    if Params.VERBOSE: print('Applying solid state convection to the ice shell based on Deschamps and Sotin (2001).')
    ConvectionDeschampsSotin2001()

    if not Planet.Do.EQUIL_Q:
        # Assign heat input into ocean from mantle to be ~radiogenic
        print('WARNING: QfromMantle_Wm2 is set to a value consistent only with Europa radiogenic heating.')
        Planet.Ocean.QfromMantle_Wm2 = 2.2e11 / 4/np.pi / Planet.Bulk.R_m**2

    return Planet


def InnerLayers(Planet, Params):
    """ Geophysical and thermodynamic calculations for silicate and core layers
        Calculates state variables of the layer with each pressure step

        Assigns Planet attributes:
            Steps.nTotal, phase
    """
    Planet = FindSeafloorMoI(Planet, Params)

    if Planet.Steps.nHydro <= Planet.Steps.nSurfIce: print('WARNING: For these run settings, the hydrosphere is entirely frozen.')
    Planet.Steps.nTotal = Planet.Steps.nHydro + Planet.Steps.nSil + Planet.Steps.nCore

    # Assign phase values for silicates and core
    Planet.phase[Planet.Steps.nHydro:Planet.Steps.nHydro + Planet.Steps.nSil] = 50
    Planet.phase[Planet.Steps.nHydro + Planet.Steps.nSil:Planet.Steps.nTotal] = 100

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

    Planet = SilicateLayers(Planet, Params)

    if Planet.Do.Fe_CORE:
        Planet = IronCoreLayers(Planet, Params)

    return Planet


def FindSeafloorMoI(Planet, Params):
    """ Find the relative sizes of silicate, core, and hydrosphere layers that are
        consistent with the measured moment of inertia, based on calculated hydrosphere
        properties and assumptions about the silicate and possible core layers.

        Assigns Planet attributes:
            CMR2mean, Sil.RsilMean_m, Sil.RsilRange_m, Core.RFeMean_m, Core.RFeRange_m, Steps.nHydro
    """
    if Params.VERBOSE: print('Finding MoI consistent with measured value...')
    MR2 = Planet.Bulk.M_kg * Planet.Bulk.R_m**2
    nHydroActual = Planet.Steps.nSurfIce + Planet.Steps.nOceanMax
    # Find contribution to axial moment of inertia C from each ocean layer
    dC_H2O = 8*np.pi/15 * Planet.rho_kgm3[:-1] * (Planet.r_m[:-1]**5 - Planet.r_m[1:]**5)
    # Find total mass contained above each hydrosphere layer
    MAbove_kg = np.array([np.sum(Planet.MLayer_kg[:i]) for i in range(nHydroActual)])
    # Initialize arrays for axial moment of inertia C and other necessary quantities
    C, rCore_m = (np.zeros(nHydroActual-1) for _ in range(2))
    # Find volume of a full sphere of silicate corresponding to each valid layer
    VsilSphere_m3 = 4/3*np.pi * Planet.r_m[1:]**3

    if Planet.Do.Fe_CORE:
        # If there is a core, we assume fixed densities in the silicates and core to
        # make the problem of finding a matching MoI tractable.

        # Find core bulk density based on assumed sulfide content (Eq 10 of Vance et al., 2014)
        rhoCore_kgm3 = Planet.Core.rhoFeS_kgm3 * Planet.Core.rhoFe_kgm3 / \
            (Planet.Core.xFeS * (Planet.Core.rhoFe_kgm3 - Planet.Core.rhoFeS_kgm3) + Planet.Core.rhoFeS_kgm3)
        # Calculate core volume for a silicate layer with outer radius equal to bottom of each hydrosphere layer
        # and inner radius equal to the core radius
        VCore_m3 = np.array([(Planet.Bulk.M_kg - MAbove_kg[i] - VsilSphere_m3[i]*Planet.Sil.rhoSilWithCore_kgm3) /
                             (rhoCore_kgm3 - Planet.Sil.rhoSilWithCore_kgm3) for i in range(nHydroActual-1)])
        # Find values for which the silicate radius is too large
        nTooBig = next((i[0] for i, val in np.ndenumerate(VCore_m3) if val>0))
        # Calculate corresponding core radii based on above density
        rCore_m[nTooBig:] = (VCore_m3[nTooBig:]*3/4/np.pi)**(1/3)
        # Assign fixed density to an array for dual-use code looking for compatible C/MR^2
        rhoSil_kg = np.ones_like(rCore_m) * Planet.Sil.rhoSilWithCore_kgm3
    else:
        # When there is no core, we treat the silicate density as a free parameter to
        # find a matching solution for the MoI.

        # Find silicate density consistent with observed bulk mass for each radius
        rhoSil_kg = np.array([(Planet.Bulk.M_kg - MAbove_kg[i]) / VsilSphere_m3[i] for i in range(nHydroActual-1)])
        # Density of silicates is scaled to fit the total mass, so there is no nTooBig in this case.
        nTooBig = 0
        # Set core density to zero so calculations can proceed
        rhoCore_kgm3 = 0

    silEOS = PerplexEOSStruct(Planet.Sil.mantleEOS, EOSinterpMethod=Params.interpMethod)
    coreEOS = PerplexEOSStruct(Planet.Core.coreEOS)
    silStuff = MantleEOS(silEOS, Planet.P_MPa[Planet.Steps.nHydro], Planet.T_K[Planet.Steps.nHydro],
                         Planet.r_m[Planet.Steps.nHydro], rCore_m)

    # Calculate C for a mantle extending up to each hydrosphere layer in turn
    C[nTooBig:] = [np.sum(dC_H2O[:i+1]) + \
            8*np.pi/15 * rhoSil_kg[i] * (Planet.r_m[i]**5 - rCore_m[i]**5) + \
            8*np.pi/15 * rhoCore_kgm3 * rCore_m[i]**5 \
            for i in range(nTooBig, nHydroActual-1)]
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
    # Record the best-match C/MR^2 value
    Planet.CMR2mean = CMR2[iCMR2]
    # Now we finally know how many layers there are in the hydrosphere
    Planet.Steps.nHydro = iCMR2 - 1

    # Fill core/mantle trade arrays and set mean values consistent with MoI
    Planet.Sil.rhoMean_kgm3 = rhoSil_kg[iCMR2]  # PLACEHOLDER until mantle EOS searching is implemented for C/MR^2 calculations
    Planet.Sil.rhoTrade_kgm3 = rhoSil_kg[CMR2inds]
    Planet.Sil.Rmean_m = Planet.r_m[iCMR2]
    Planet.Sil.Rtrade_m = Planet.r_m[CMR2inds]
    Planet.Sil.Rrange_m = Planet.Sil.Rtrade_m[0] - Planet.Sil.Rtrade_m[-1]
    Planet.Core.rhoMean_kgm3 = rhoCore_kgm3  # PLACEHOLDER until core EOS searching is implemented for C/MR^2 calculations
    Planet.Core.Rmean_m = rCore_m[iCMR2]
    Planet.Core.Rtrade_m = rCore_m[CMR2inds]
    Planet.Core.Rrange_m = Planet.Core.Rtrade_m[0] - Planet.Core.Rtrade_m[-1]

    if Params.VERBOSE: print('Found matching MoI of ' + str(round(Planet.CMR2mean,3)) +
                           ' (C/MR^2 = ' + str(round(Planet.Bulk.Cmeasured,3)) + '±' + str(round(Planet.Bulk.Cuncertainty,3)) + ') for ' +
                             'rho_sil = ' + str(round(Planet.Sil.rhoMean_kgm3)) + ' kg/m^3, ' +
                             'R_sil = ' + str(round(Planet.Sil.Rmean_m / Planet.Bulk.R_m,2)) + 'R, ' +
                             'R_core = ' + str(round(Planet.Core.Rmean_m / Planet.Bulk.R_m,2)) + 'R.')

    return Planet


def SilicateLayers(Planet, Params):
    """ Geophysical and thermodynamic calculations for silicate layers
        Calculates state variables of the layer with each pressure step

        Assigns Planet attributes:
            phase
    """

    return Planet


def IronCoreLayers(Planet, Params):
    """ Geophysical and thermodynamic calculations for core layers
        Calculates state variables of the layer with each pressure step

        Assigns Planet attributes:
            phase, g_ms2
    """

    nCoreStart = Planet.Steps.nHydro + Planet.Steps.nSil
    nCoreEnd = Planet.Steps.nTotal
    # Find gravity at outer core radius
    Planet.g_ms2[nCoreStart] = Constants.G * (Planet.Bulk.M_kg - np.sum(Planet.MLayer_kg[:nCoreStart])) \
                               / Planet.r_m[nCoreStart]**2
    # Set gravity to simply be linear in the core (exact for a uniform-density core, approximate if
    # density variations are only slight in the core. Doing so avoids negative gravity problems.
    Planet.g_ms2[nCoreStart:nCoreEnd] = Planet.g_ms2[nCoreStart] / Planet.r_m[nCoreStart] * Planet.r_m[nCoreStart:nCoreEnd]

    return Planet
