import numpy as np
import logging
from PlanetProfile.Thermodynamics.HydroEOS import PhaseConv
from PlanetProfile.Thermodynamics.ThermalProfiles.ThermalProfiles import ConductiveTemperature
from PlanetProfile.Utilities.defineStructs import Constants

# Assign logger
log = logging.getLogger('PlanetProfile')

def EvalLayerProperties(Planet, Params, iStart, iEnd, EOS, P_MPa, T_K):
    """ Evaluate EOS functions for a given pressure and temperature range.

        Args:
            iStart, iEnd (int): Layer array indices to start and end the calculations.
                P_MPa[0] should correspond to Planet.P_MPa[iStart], etc.
            EOS (EOSStruct): Query-able ice, ocean, sil, or core EOS containing functions
                for the listed layer properties to be calculated.
            P_MPa, T_K (float, shape (iEnd-iStart)): List of (P,T) pairs to evaluate.
        Assigns Planet attributes:
            rhoMatrix_kgm3, Cp_JkgK, alpha_pK, kTherm_WmK
    """

    Planet.rhoMatrix_kgm3[iStart:iEnd] = EOS.fn_rho_kgm3(  P_MPa, T_K)
    Planet.Cp_JkgK[iStart:iEnd] =        EOS.fn_Cp_JkgK(   P_MPa, T_K)
    Planet.alpha_pK[iStart:iEnd] =       EOS.fn_alpha_pK(  P_MPa, T_K)
    Planet.kTherm_WmK[iStart:iEnd] =     EOS.fn_kTherm_WmK(P_MPa, T_K)

    return Planet


def PorosityCorrectionVacIce(Planet, Params, iStart, iEnd, EOS, P_MPa, T_K):
    """ Correct layer properties retrieved with EvalLayerProperties according
        to the porosity of the material. This function assumes pores are empty.

        Args:
            iStart, iEnd (int): Layer array indices to start and end the calculations.
                P_MPa[0] should correspond to Planet.P_MPa[iStart], etc.
            EOS (EOSStruct): Query-able ice EOS containing functions for the listed
                layer properties to be calculated, including fn_phi_frac.
            P_MPa, T_K (float, shape (iEnd-iStart)): List of (P,T) pairs to evaluate.
        Assigns Planet attributes:
            phi_frac, rho_kgm3, Cp_JkgK, alpha_pK, kTherm_WmK
    """

    Planet.phi_frac[iStart:iEnd] = EOS.fn_phi_frac(P_MPa, T_K)
    # Only evaluate at places where porosity is non-negligible, so that we
    # can limit calling EOS functions with invalid parameters, which can give NaNs
    iEval = np.where(Planet.phi_frac[iStart:iEnd] >= Planet.Ocean.phiMin_frac)[0]
    iSolid = np.where(Planet.phi_frac[iStart:iEnd] < Planet.Ocean.phiMin_frac)[0]

    # Assign matrix densities as total densities for those layers with porosities
    # low enough to be below our modeling threshold
    Planet.rho_kgm3[iStart:iEnd][iSolid] = Planet.rhoMatrix_kgm3[iStart:iEnd][iSolid]
    # Correct each quantity according to J rules (see Thermodynamics.InnerEOS)
    Planet.rho_kgm3[iStart:iEnd][iEval] = EOS.fn_porosCorrect(
        Planet.rhoMatrix_kgm3[iStart:iEnd][iEval], 0,
        Planet.phi_frac[iStart:iEnd][iEval], Planet.Ocean.Jrho)
    Planet.Cp_JkgK[iStart:iEnd][iEval] = EOS.fn_porosCorrect(
        Planet.Cp_JkgK[iStart:iEnd][iEval] * Planet.rhoMatrix_kgm3[iStart:iEnd][iEval], 0,
        Planet.phi_frac[iStart:iEnd][iEval], Planet.Ocean.JCp) / Planet.rho_kgm3[iStart:iEnd][iEval]
    Planet.alpha_pK[iStart:iEnd][iEval] = EOS.fn_porosCorrect(
        Planet.alpha_pK[iStart:iEnd][iEval], 0,
        Planet.phi_frac[iStart:iEnd][iEval], Planet.Ocean.Jalpha)
    Planet.kTherm_WmK[iStart:iEnd][iEval] = EOS.fn_porosCorrect(
        Planet.kTherm_WmK[iStart:iEnd][iEval], 0,
        Planet.phi_frac[iStart:iEnd][iEval], Planet.Ocean.JkTherm)

    return Planet


def PorosityCorrectionFilledIce(Planet, Params, iStart, iEnd, EOS, EOSpore):
    """ Correct ice layer properties retrieved with EvalLayerProperties according
        to the porosity of the material. This function assumes pores contain
        ocean fluids.

        Args:
            iStart, iEnd (int): Layer array indices to start and end the calculations.
                P_MPa[0] should correspond to Planet.P_MPa[iStart], etc.
            EOS (EOSStruct): Query-able ice EOS for the matrix material
                containing functions for the listed layer properties to be calculated,
                including fn_phi_frac.
            EOSpore (EOSStruct): Ocean EOS to evaluate for pore material.
        Assigns Planet attributes:
            Ppore_MPa, phi_frac, rho_kgm3, Cp_JkgK, alpha_pK, kTherm_WmK
        Returns:
            phasePore
    """

    alphaPeff = Planet.Ocean.alphaPeff[EOS.phaseStr]
    DeltaPpore_MPa = 0.0
    if iStart == 0:
        # Record last entry in Ppore, and set it as needed to ensure
        # Ppore[i-1] reference will work correctly in loop
        PporeSave = Planet.Ppore_MPa[-1] + 0.0
        Planet.Ppore_MPa[-1] = Planet.P_MPa[0] + 0.0
    elif Planet.Ppore_MPa[iStart-1] == 0:
        # Make sure Ppore[i-1] has been assigned, and set it so Ppore[iStart] = P_MPa[iStart],
        # because we are at the top of a porous layer, above which is liquid or vacuum
        Planet.Ppore_MPa[iStart-1] = Planet.P_MPa[iStart] + 0.0

    for i in range(iStart, iEnd):
        Planet.Ppore_MPa[i] = Planet.Ppore_MPa[i-1] + DeltaPpore_MPa
        Peff_MPa = Planet.P_MPa[i] - alphaPeff * Planet.Ppore_MPa[i]
        Planet.phi_frac[i] = EOS.fn_phi_frac(Peff_MPa, Planet.T_K[i])
        if Planet.phi_frac[i] >= Planet.Ocean.phiMin_frac:
            # In ices, assume pores are filled only with liquid.
            rhoPore_kgm3 =   EOSpore.fn_rho_kgm3(  Planet.Ppore_MPa[i], Planet.T_K[i])
            CpPore_JkgK =    EOSpore.fn_Cp_JkgK(   Planet.Ppore_MPa[i], Planet.T_K[i])
            alphaPore_pK =   EOSpore.fn_alpha_pK(  Planet.Ppore_MPa[i], Planet.T_K[i])
            kThermPore_WmK = EOSpore.fn_kTherm_WmK(Planet.Ppore_MPa[i], Planet.T_K[i])

            # Correct each quantity according to J rules (see Thermodynamics.InnerEOS)
            Planet.rho_kgm3[i] = EOS.fn_porosCorrect(Planet.rhoMatrix_kgm3[i], rhoPore_kgm3,
                                                       Planet.phi_frac[i], Planet.Ocean.Jrho)
            Planet.Cp_JkgK[i] = EOS.fn_porosCorrect(Planet.Cp_JkgK[i] * Planet.rhoMatrix_kgm3[i], CpPore_JkgK * rhoPore_kgm3,
                                                       Planet.phi_frac[i], Planet.Ocean.JCp) / Planet.rho_kgm3[i]
            Planet.alpha_pK[i] = EOS.fn_porosCorrect(Planet.alpha_pK[i], alphaPore_pK,
                                                       Planet.phi_frac[i], Planet.Ocean.Jalpha)
            Planet.kTherm_WmK[i] = EOS.fn_porosCorrect(Planet.kTherm_WmK[i], kThermPore_WmK,
                                                       Planet.phi_frac[i], Planet.Ocean.JkTherm)
            # Get approximate pore pressure change, using the previous depth change. This value
            # should be close to the correct value -- we do it this way to avoid array indexing
            # issues at the start/end points and the need to calculate the depth change twice
            DeltaPpore_MPa = 1e-6 * rhoPore_kgm3 * Planet.g_ms2[i] * (Planet.z_m[i] - Planet.z_m[i-1])
        else:
            Planet.phi_frac[i] = 0.0
            Planet.rho_kgm3[i] = Planet.rhoMatrix_kgm3[i] + 0.0
            Planet.Ppore_MPa[i] = Planet.P_MPa[i] + 0.0
            DeltaPpore_MPa = Planet.P_MPa[i] - Planet.P_MPa[i-1]

    # Restore final Ppore entry if it was beyond those affected by the loop,
    # because we used it as the -1 index to start the loop
    if iStart == 0 and iEnd < np.size(Planet.Ppore_MPa):
        Planet.Ppore_MPa[-1] = PporeSave

    return Planet, phasePore


def PropagateConduction(Planet, Params, iStart, iEnd):
    """ Use P, T, and layer properties as determined from conductive layer profile to evaluate
        an EOS to get the layer thicknesses, gravity, and masses.

        Args:
            iStart, iEnd (int): Layer array indices corresponding to the last evaluated
                layer and the end of the conductive profile (e.g. material transition),
                respectively.
        Assigns Planet attributes:
            z_m, r_m, MLayer_kg, g_ms2
    """

    # Add a catch in case we call with invalid indices, which is convenient
    # after convection calculations when no convection is happening
    if iStart < iEnd:
        thisMAbove_kg = np.sum(Planet.MLayer_kg[:iStart])
        for i in range(iStart+1, iEnd+1):
            # Increment depth based on change in pressure, combined with gravity and density
            Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1]) * 1e6 / Planet.g_ms2[i-1] / \
                            Planet.rho_kgm3[i-1]
            Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
            Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1] ** 3 - Planet.r_m[i] ** 3)
            thisMAbove_kg += Planet.MLayer_kg[i-1]
            thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
            Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i] ** 2
            log.debug(f'il: {i:d}; P_MPa: {Planet.P_MPa[i]:.3f}; T_K: {Planet.T_K[i]:.3f}; phase: {Planet.phase[i]:d}')

    return Planet


def PropagateAdiabaticSolid(Planet, Params, iStart, iEnd, EOS):
    """ Use layer-top values and assumption of an adiabatic thermal profile
        to evaluate the conditions at the bottom of each layer in the specified
        zone (iStart to iEnd). This function assumes no porosity.

        Args:
            iStart, iEnd (int): Layer array indices corresponding to the last evaluated
                layer and the end of the conductive profile (e.g. material transition),
                respectively.
            EOS (EOSStruct): Ice, ocean, sil, or core EOS to query for layer properties.
        Assigns Planet attributes:
            z_m, r_m, MLayer_kg, g_ms2, rhoMatrix_kgm3, Cp_JkgK, alpha_pK, kTherm_WmK,
            rho_kgm3
    """

    thisMAbove_kg = np.sum(Planet.MLayer_kg[:iStart-1])
    for i in range(iStart, iEnd):
        # Increment depth based on change in pressure, combined with gravity and density
        Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1]) * 1e6 / Planet.g_ms2[i-1] / \
                        Planet.rhoMatrix_kgm3[i-1]
        # Convert depth to radius
        Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
        # Calculate layer mass
        Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rhoMatrix_kgm3[i-1] * (Planet.r_m[i-1] ** 3 - Planet.r_m[i] ** 3)
        thisMAbove_kg += Planet.MLayer_kg[i-1]
        thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
        # Use remaining mass below in Gauss's law for gravity to get g at the top of this layer
        Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i] ** 2

        # Propagate adiabatic thermal profile
        Planet.T_K[i] = Planet.T_K[i-1] + Planet.T_K[i-1] * Planet.alpha_pK[i-1] / \
                        Planet.Cp_JkgK[i-1] / Planet.rhoMatrix_kgm3[i-1] * (Planet.P_MPa[i] - Planet.P_MPa[i-1]) * 1e6
        # Now use P and T for this layer to get physical properties
        Planet.rhoMatrix_kgm3[i] = EOS.fn_rho_kgm3(  Planet.P_MPa[i], Planet.T_K[i])
        Planet.Cp_JkgK[i] =        EOS.fn_Cp_JkgK(   Planet.P_MPa[i], Planet.T_K[i])
        Planet.alpha_pK[i] =       EOS.fn_alpha_pK(  Planet.P_MPa[i], Planet.T_K[i])
        Planet.kTherm_WmK[i] =     EOS.fn_kTherm_WmK(Planet.P_MPa[i], Planet.T_K[i])

        log.debug(f'il: {i:d}; P_MPa: {Planet.P_MPa[i]:.3f}; T_K: {Planet.T_K[i]:.3f}; phase: {Planet.phase[i]:d}')

    Planet.rho_kgm3[iStart:iEnd] = Planet.rhoMatrix_kgm3[iStart:iEnd] + 0.0

    return Planet


def PropagateAdiabaticPorousVacIce(Planet, Params, iStart, iEnd, EOS):
    """ Use layer-top values and assumption of an adiabatic thermal profile in ices
        to evaluate the conditions at the bottom of each layer in the specified
        zone (iStart to iEnd). This function assumes pores are empty.

        Args:
            iStart, iEnd (int): Layer array indices corresponding to the last evaluated
                layer and the end of the conductive profile (e.g. material transition),
                respectively.
            EOS (EOSStruct): Ice EOS to query for layer properties.
        Assigns Planet attributes:
            z_m, r_m, MLayer_kg, g_ms2, rhoMatrix_kgm3, Cp_JkgK, alpha_pK, kTherm_WmK,
            phi_frac, rho_kgm3
    """

    if iStart == 0:
        raise RuntimeError('Adiabatic calculations rely on overlying layer properties to begin recursion.')
    thisMAbove_kg = np.sum(Planet.MLayer_kg[:iStart-1])

    for i in range(iStart, iEnd):
        # Increment depth based on change in pressure, combined with gravity and density
        Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
        # Convert depth to radius
        Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
        # Calculate layer mass
        Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
        thisMAbove_kg += Planet.MLayer_kg[i-1]
        thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
        # Use remaining mass below in Gauss's law for gravity to get g at the top of this layer
        Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2

        # Propagate adiabatic thermal profile
        Planet.T_K[i] = Planet.T_K[i-1] + Planet.T_K[i-1] * Planet.alpha_pK[i-1] / \
                        Planet.Cp_JkgK[i-1] / Planet.rho_kgm3[i-1] * (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6
        # Now use P and T for this layer to get physical properties
        Planet.rhoMatrix_kgm3[i] = EOS.fn_rho_kgm3(  Planet.P_MPa[i], Planet.T_K[i])
        Planet.Cp_JkgK[i] =        EOS.fn_Cp_JkgK(   Planet.P_MPa[i], Planet.T_K[i])
        Planet.alpha_pK[i] =       EOS.fn_alpha_pK(  Planet.P_MPa[i], Planet.T_K[i])
        Planet.kTherm_WmK[i] =     EOS.fn_kTherm_WmK(Planet.P_MPa[i], Planet.T_K[i])

        # Correct for porosity in ice I layers, assuming pores are empty.
        Planet.phi_frac[i] = EOS.fn_phi_frac(Planet.P_MPa[i], Planet.T_K[i])
        # Only model porosity if above some threshold amount, to save on computation when
        # the porosity is negligible.
        if Planet.phi_frac[i] >= Planet.Ocean.phiMin_frac:
            Planet.rho_kgm3[i] = EOS.fn_porosCorrect(Planet.rhoMatrix_kgm3[i], 0,
                                                       Planet.phi_frac[i], Planet.Ocean.Jrho)
            Planet.Cp_JkgK[i] = EOS.fn_porosCorrect(Planet.Cp_JkgK[i] * Planet.rhoMatrix_kgm3[i], 0,
                                                       Planet.phi_frac[i], Planet.Ocean.JCp) / Planet.rho_kgm3[i]
            Planet.alpha_pK[i] = EOS.fn_porosCorrect(Planet.alpha_pK[i], 0,
                                                       Planet.phi_frac[i], Planet.Ocean.Jalpha)
            Planet.kTherm_WmK[i] = EOS.fn_porosCorrect(Planet.kTherm_WmK[i], 0,
                                                       Planet.phi_frac[i], Planet.Ocean.JkTherm)
        else:
            Planet.phi_frac[i] = 0.0
            Planet.rho_kgm3[i] = Planet.rhoMatrix_kgm3[i] + 0.0

        log.debug(f'il: {i:d}; P_MPa: {Planet.P_MPa[i]:.3f}; T_K: {Planet.T_K[i]:.3f}; phase: {Planet.phase[i]:d}')

    return Planet


def PropagateAdiabaticPorousFilledIce(Planet, Params, iStart, iEnd, EOS, EOSpore):
    """ Use layer-top values and assumption of an adiabatic thermal profile in ices
        to evaluate the conditions at the bottom of each layer in the specified
        zone (iStart to iEnd). This function assumes pores are filled with ocean
        fluids.

        Args:
            iStart, iEnd (int): Layer array indices corresponding to the last evaluated
                layer and the end of the conductive profile (e.g. material transition),
                respectively.
            EOS (EOSStruct): Ice EOS to query for layer properties.
            EOSpore (EOSStruct): Ocean EOS to query for pore material properties.
        Assigns Planet attributes:
            Ppore_MPa, z_m, r_m, MLayer_kg, g_ms2, rhoMatrix_kgm3, Cp_JkgK, alpha_pK, kTherm_WmK,
            phi_frac, rho_kgm3
        Returns:
            phasePore
    """

    # Get the coupling between pore pressure and effective pressure (how effectively
    # the pressure of pore materials resists their closure) for this material, so
    # we don't have to keep evaluating it as we do here.
    alphaPeff = Planet.Ocean.alphaPeff[EOS.phaseStr]

    # Assign phasePore all zeros. Since we assume pores contain liquid, we won't
    # change any of these values. This is for forward compatibility, in case the
    # assumption of pores containing ocean fluid is relaxed.
    phasePore = np.zeros(iEnd - iStart, dtype=np.int_)
    # Initialize DeltaPpore at 0 so we can use it to set the top pore pressure
    # equal to the pressure of the overlying material. In this, we're assuming
    # iStart corresponds to a layer for which there is communication between
    # the overlying material and pore space, i.e. that overlying ocean fluids
    # contribute the layer-top pressure for both matrix and pore materials.
    DeltaPpore_MPa = 0.0

    if iStart == 0:
        raise RuntimeError('Adiabatic calculations rely on overlying layer properties to begin recursion.')
    # Initialize overlying mass
    thisMAbove_kg = np.sum(Planet.MLayer_kg[:iStart-1])
    # Ensure the top-layer pore pressure will be incremented up to the matching overburden pressure
    Planet.Ppore_MPa[iStart-1] = Planet.P_MPa[iStart] + 0.0

    # Begin adiabatic profile propagation
    for i in range(iStart, iEnd):
        # Increment depth based on change in pressure, combined with gravity and density
        Planet.z_m[i] = Planet.z_m[i-1] + (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6 / Planet.g_ms2[i-1] / Planet.rho_kgm3[i-1]
        # Convert depth to radius
        Planet.r_m[i] = Planet.Bulk.R_m - Planet.z_m[i]
        # Calculate layer mass
        Planet.MLayer_kg[i-1] = 4/3*np.pi * Planet.rho_kgm3[i-1] * (Planet.r_m[i-1]**3 - Planet.r_m[i]**3)
        thisMAbove_kg += Planet.MLayer_kg[i-1]
        thisMBelow_kg = Planet.Bulk.M_kg - thisMAbove_kg
        # Use remaining mass below in Gauss's law for gravity to get g at the top of this layer
        Planet.g_ms2[i] = Constants.G * thisMBelow_kg / Planet.r_m[i]**2

        # Propagate adiabatic thermal profile
        Planet.T_K[i] = Planet.T_K[i-1] + Planet.T_K[i-1] * Planet.alpha_pK[i-1] / \
                        Planet.Cp_JkgK[i-1] / Planet.rho_kgm3[i-1] * (Planet.P_MPa[i] - Planet.P_MPa[i-1])*1e6
        # Now use P and T for this layer to get physical properties
        Planet.rhoMatrix_kgm3[i] = EOS.fn_rho_kgm3(  Planet.P_MPa[i], Planet.T_K[i])
        Planet.Cp_JkgK[i] =        EOS.fn_Cp_JkgK(   Planet.P_MPa[i], Planet.T_K[i])
        Planet.alpha_pK[i] =       EOS.fn_alpha_pK(  Planet.P_MPa[i], Planet.T_K[i])
        Planet.kTherm_WmK[i] =     EOS.fn_kTherm_WmK(Planet.P_MPa[i], Planet.T_K[i])

        # Correct for porosity in ice I layers
        # First, increment pore pressure based on hydrostatic pressure from
        # overlying pore material, assuming pore connectivity.
        Planet.Ppore_MPa[i] = Planet.Ppore_MPa[i-1] + DeltaPpore_MPa
        # Next, use pore pressure coupling constant (alpha) to find effective
        # pore closure pressure
        Peff_MPa = Planet.P_MPa[i] - alphaPeff * Planet.Ppore_MPa[i]
        # Use Peff to find porosity (alpha = 0 for Do.P_EFFECTIVE = False)
        Planet.phi_frac[i] = EOS.fn_phi_frac(Peff_MPa, Planet.T_K[i])

        # Only model porosity if above some threshold amount, to save on computation when
        # the porosity is negligible.
        if Planet.phi_frac[i] >= Planet.Ocean.phiMin_frac:
            # In ices, assume pores are filled only with liquid. Find liquid physical properties
            rhoPore_kgm3 =   EOSpore.fn_rho_kgm3(  Planet.Ppore_MPa[i], Planet.T_K[i])
            CpPore_JkgK =    EOSpore.fn_Cp_JkgK(   Planet.Ppore_MPa[i], Planet.T_K[i])
            alphaPore_pK =   EOSpore.fn_alpha_pK(  Planet.Ppore_MPa[i], Planet.T_K[i])
            kThermPore_WmK = EOSpore.fn_kTherm_WmK(Planet.Ppore_MPa[i], Planet.T_K[i])

            # Combine them with matrix physical properties
            Planet.rho_kgm3[i] = EOS.fn_porosCorrect(Planet.rhoMatrix_kgm3[i], rhoPore_kgm3,
                                                       Planet.phi_frac[i], Planet.Ocean.Jrho)
            Planet.Cp_JkgK[i] = EOS.fn_porosCorrect(Planet.Cp_JkgK[i] * Planet.rhoMatrix_kgm3[i], CpPore_JkgK * rhoPore_kgm3,
                                                       Planet.phi_frac[i], Planet.Ocean.JCp) / Planet.rho_kgm3[i]
            Planet.alpha_pK[i] = EOS.fn_porosCorrect(Planet.alpha_pK[i], alphaPore_pK,
                                                       Planet.phi_frac[i], Planet.Ocean.Jalpha)
            Planet.kTherm_WmK[i] = EOS.fn_porosCorrect(Planet.kTherm_WmK[i], kThermPore_WmK,
                                                       Planet.phi_frac[i], Planet.Ocean.JkTherm)
            # Get approximate pore pressure change, using the previous depth change. This value
            # should be close to the correct value -- we do it this way to avoid array indexing
            # issues at the start/end points and the need to calculate the depth change twice
            DeltaPpore_MPa = 1e-6 * rhoPore_kgm3 * Planet.g_ms2[i] * (Planet.z_m[i] - Planet.z_m[i-1])
        else:
            # Zero out porosity and ignore it if it's below threshold (this also takes care of
            # negative porosity values that can happen from our RectBivariateSpline implementation).
            Planet.phi_frac[i] = 0.0
            Planet.rho_kgm3[i] = Planet.rhoMatrix_kgm3[i] + 0.0
            Planet.Ppore_MPa[i] = Planet.P_MPa[i] + 0.0
            DeltaPpore_MPa = Planet.P_MPa[i] - Planet.P_MPa[i-1]

        log.debug(f'il: {i:d}; P_MPa: {Planet.P_MPa[i]:.3f}; T_K: {Planet.T_K[i]:.3f}; phase: {Planet.phase[i]:d}')

    return Planet, phasePore


def PropagateConductionProfilesSolid(Planet, Params, nProfiles, profRange, rSilEnd_m):
    """ Same as PropagateConductionSolid, but for silicate layers, for which we calculate
        a number of conductive profiles simultaneously as part of MoI matching.

        Args:
            nProfiles (int): Number of parallel layer profiles to calculate for silicates.
                Determined by the number of hydrosphere layers greater than Steps.iSilStart,
                because we use each hydrosphere layer as a starting radius option for the
                silicates.
            profRange (int range, Iterable): Usually range(nProfiles), but included for
                compatibility with waterless bodies using the same calculations.
            rSilEnd_m (float): Inner radius of silicates, if known. Usually 0, but nonzero
                if we already found the core radius, i.e. by assuming constant silicate and
                core densities using Do.CONSTANT_INNER_DENSITY = True.
        Returns:
            Psil_MPa (float, shape (nProfiles, Planet.Steps.nSilMax)): Pressures in silicate layers in MPa.
                Test-case 2D array to evaluate -- we pick the closest mass match along one dimension
                and the pressures correspond to the silicate layers for that match along the other
                dimension.
            rSil_m (float, shape (nProfiles, Planet.Steps.nSilMax+1)): Layer outer radii in m. Note
                that this array contains one extra value along the profile axis, i.e. an extra layer
                radius for each profile. This is for convenience in taking differences between values,
                and because r = 0 is not at the top of any layer. The last value does not get assigned
                to Planet.r_m.
            Tsil_K, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg, gSil_kg, phiSil_frac,
                HtidalSil_Wm3: Same as Psil_MPa but for other layer physical properties.
            MHydro_kg (float, shape nProfiles): Total mass of hydrosphere above each silicate profile in kg.
            KSsil_GPa, GSsil_GPa (float, shape (nProfiles, Planet.Steps.nSilMax)): Bulk and shear moduli,
                respectively, of combined porous layer properties in GPa.
            Ppore_kg, rhoPore_kgm3 (float, shape (nProfiles, Planet.Steps.nSilMax)): Pore properties.
            phasePore (int, (nProfiles, Planet.Steps.nSilMax)): Liquid/ice phase indices of pore materials.
            qTop_Wm2 (float, shape nProfiles): Heat flux leaving the top of the mantle for each profile in W/m^2.
    """

    # Initialize output arrays and working arrays
    Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg, \
    MHydro_kg, gSil_ms2, phiSil_frac, HtidalSil_Wm3, KSsil_GPa, GSsil_GPa, \
    Ppore_MPa, rhoPore_kgm3, phasePore, qTop_Wm2, fn_g_ms2 \
        = InitSil(Planet, Params, nProfiles, profRange, rSilEnd_m)

    # Calculate initial values based on matrix properties since we're not modeling porosity
    MLayerSil_kg[:,0] = rhoSil_kgm3[:,0] * 4/3*np.pi*(rSil_m[:,0]**3 - rSil_m[:,1]**3)
    HtidalSil_Wm3[:,0] = Planet.Sil.fn_Htidal_Wm3(rhoSil_kgm3[:,0], gSil_ms2[:,0], KSsil_GPa, GSsil_GPa)

    # Layer recursion across profiles
    Psil_MPa, Tsil_K, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, \
    HtidalSil_Wm3, kThermSil_WmK = SilRecursionSolid(Planet, Params,
        Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg,
        gSil_ms2, HtidalSil_Wm3, KSsil_GPa, GSsil_GPa, qTop_Wm2, fn_g_ms2)

    return Planet, Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, \
           gSil_ms2, phiSil_frac, HtidalSil_Wm3, kThermSil_WmK, Ppore_MPa, rhoSil_kgm3, \
           rhoPore_kgm3, phasePore


def PropagateConductionProfilesPorous(Planet, Params, nProfiles, profRange, rSilEnd_m):
    """ See PropagateConductionProfilesSolid for variable descriptions.
        Generally, Sil corrsponds to the rock matrix, Pore corresponds to the pore
        material, and Tot is the combined physical properties of both.
    """

    # Initialize output arrays and working arrays in common with non-porous modeling
    # with top-layer values
    Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg, \
    MHydro_kg, gSil_ms2, phiSil_frac, HtidalSil_Wm3, KSsil_GPa, GSsil_GPa, \
    Ppore_MPa, rhoPore_kgm3, phasePore, qTop_Wm2, fn_g_ms2 \
        = InitSil(Planet, Params, nProfiles, profRange, rSilEnd_m)

    # Initialize porosity-specific arrays
    kThermPore_WmK, KSpore_GPa, GSpore_GPa, DeltaPpore_MPa, \
    kThermTot_WmK, KStot_GPa, GStot_GPa, rhoTot_kgm3, \
    MLayerSil_kg[:,0], HtidalSil_Wm3[:,0], rhoPore_kgm3[:,0], phasePore[:,0], \
    Ppore_MPa[:,0], phiSil_frac[:,0] = InitPorous(Planet, Params, nProfiles,
        rSil_m[:,0], rSil_m[:,1], Psil_MPa[:,0], Tsil_K[:,0], gSil_ms2[:,0],
        kThermSil_WmK[:,0], rhoSil_kgm3[:,0], KSsil_GPa[:,0], GSsil_GPa[:,0])

    # Layer recursion across profiles
    Psil_MPa, Tsil_K, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, \
    HtidalSil_Wm3, kThermSil_WmK, rhoTot_kgm3, phiSil_frac, kThermTot_WmK, \
    Ppore_MPa, rhoPore_kgm3, phasePore = SilRecursionPorous(Planet, Params,
        Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg,
        gSil_ms2, HtidalSil_Wm3, qTop_Wm2, fn_g_ms2,
        KSsil_GPa, GSsil_GPa, phiSil_frac, Ppore_MPa, phasePore,
        rhoPore_kgm3, kThermPore_WmK, KSpore_GPa, GSpore_GPa, DeltaPpore_MPa,
        rhoTot_kgm3, kThermTot_WmK, KStot_GPa, GStot_GPa)

    return Planet, Psil_MPa, Tsil_K, rSil_m, rhoTot_kgm3, MLayerSil_kg, MAboveSil_kg, \
           gSil_ms2, phiSil_frac, HtidalSil_Wm3, kThermTot_WmK, Ppore_MPa, rhoSil_kgm3, \
           rhoPore_kgm3, phasePore


def InitSil(Planet, Params, nProfiles, profRange, rSilEnd_m):
    """ See PropagateConductionProfilesSolid for variable definitions.
        Note that we will later truncate phasePore to be 1D along the
        MoI- and mass-matching profile. The main purpose of this function
        is to assign the top-layer values for the silicate profiles to
        start the layer propagation.

        Returns:
            fn_g_ms2 (func, args: MAboveLayer_kg, rLayerBot_m): Gravity function to use for silicates. Linear if
                there is no core, and same as for the hydrosphere (Gauss's law for gravity) if not. The linear
                approximation is valid if the interior has constant density, which is clearly not the case,
                but it avoids some dramatic blow-up problems in intermediate steps and isn't super far off reality.
    """

    Psil_MPa, Tsil_K, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg, gSil_ms2, \
    phiSil_frac, HtidalSil_Wm3, Ppore_MPa, rhoPore_kgm3, KSsil_GPa, GSsil_GPa \
        = (np.zeros((nProfiles, Planet.Steps.nSilMax)) for _ in range(13))
    phasePore = np.zeros((nProfiles, Planet.Steps.nSilMax), dtype=np.int_)

    # Assign ocean-bottom options as silicate-top values
    rSil_m = np.array([np.linspace(Planet.r_m[i+Planet.Steps.iSilStart], rSilEnd_m, Planet.Steps.nSilMax+1) for i in profRange])
    Psil_MPa[:,0] = [Planet.P_MPa[i+Planet.Steps.iSilStart] for i in profRange]
    Tsil_K[:,0] = [Planet.T_K[i+Planet.Steps.iSilStart] for i in profRange]
    rhoSil_kgm3[:,0] = Planet.Sil.EOS.fn_rho_kgm3(Psil_MPa[:,0], Tsil_K[:,0])
    kThermSil_WmK[:,0] = Planet.Sil.EOS.fn_kTherm_WmK(Psil_MPa[:,0], Tsil_K[:,0])
    KSsil_GPa[:,0] = Planet.Sil.EOS.fn_KS_GPa(Psil_MPa[:,0], Tsil_K[:,0])
    GSsil_GPa[:,0] = Planet.Sil.EOS.fn_GS_GPa(Psil_MPa[:,0], Tsil_K[:,0])
    gSil_ms2[:,0] = [Planet.g_ms2[i+Planet.Steps.iSilStart] for i in profRange]

    MHydro_kg = np.array([np.sum(Planet.MLayer_kg[:i]) for i in range(Planet.Steps.iSilStart, Planet.Steps.iSilStart + nProfiles)])
    # Initialize MAbove_kg to 0th silicate layer, so that the hydrosphere mass is equal to the mass above the silicates.
    MAboveSil_kg[:,0] = MHydro_kg + 0.0
    # Initialize qTop_WmK, the heat flux leaving the top of each layer.
    # For the top layer, this is equal to the heat flux warming the hydrosphere.
    qTop_Wm2 = Planet.Ocean.QfromMantle_W / 4/np.pi/rSil_m[:,0]**2

    # Adjust gravity model behavior in case these layers extend to the origin.
    # This implementation may need to be changed to account for varying silicate
    # physical properties with depth, esp. porosity, that will affect the validity
    # of the constant-density approximation giving linear gravity.
    if rSilEnd_m < 1 and not Planet.Do.Fe_CORE:
        # Constant density approximation yields usable gravity values that
        # introduce less inaccuracy than gravity values exploding because of
        # bulk mass misfit from MoI searching
        fn_g_ms2 = lambda MAboveLayer_kg, rLayerBot_m: \
            gSil_ms2[:,0] * rLayerBot_m / rSil_m[:,0]
    else:
        fn_g_ms2 = lambda MAboveLayer_kg, rLayerBot_m: \
            Constants.G * np.abs(Planet.Bulk.M_kg - MAboveLayer_kg) / rLayerBot_m**2

    return Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg, \
           MHydro_kg, gSil_ms2, phiSil_frac, HtidalSil_Wm3, KSsil_GPa, GSsil_GPa, \
           Ppore_MPa, rhoPore_kgm3, phasePore, qTop_Wm2, fn_g_ms2


def InitPorous(Planet, Params, nProfiles, rSil_m0, rSil_m1, Psil_MPa0, Tsil_K0, gSil_ms20, kThermSil_WmK0,
    rhoSil_kgm30, KSsil_GPa0, GSsil_GPa0):
    """ Initialize the layer propagation variables needed for porosity, which are mostly
        not generated or used when porosity is not modeled.

        Args:
            nProfiles (int): Number of hydrosphere layers we use as outer silicate radii.
            rSil_m0, Psil_MPa0, Tsil_K0, etc. (float, shape nProfiles): 0-index (top) layer
                properties of the rock matrix for each profile.
            rSil_m1 (float, shape nProfiles): 1-index (next one under top) layer radius.
        Returns same as InitSil, plus:
            kThermPore_WmK, KSpore_GPa, GSpore_GPa, DeltaPpore_MPa, kThermTot_WmK, KStot_GPa,
                GStot_GPa, rhoTot_kgm3 (float, shape (nProfiles, Planet.Steps.nSilMax)):
                Pore and combined porous material property arrays.
    """

    kThermPore_WmK, KSpore_GPa, GSpore_GPa, DeltaPpore_MPa, \
        kThermTot_WmK, KStot_GPa, GStot_GPa, rhoTot_kgm3 \
        = (np.zeros((nProfiles, Planet.Steps.nSilMax)) for _ in range(8))
    # Set uppermost pore pressure equal to the uppermost matrix material
    Ppore_MPa0 = Psil_MPa0 + 0.0
    # Get effective pore closure pressure for top layer
    Peff_MPa = Psil_MPa0 - Planet.Sil.alphaPeff * Ppore_MPa0
    # Get porosity of top layer
    phiSil_frac0 = Planet.Sil.fn_phi_frac(Peff_MPa, Tsil_K0)
    # Get pore fluid properties
    # First check for HP ice phases
    phasePore0 = Planet.Sil.poreEOS.fn_phase(Ppore_MPa0, Tsil_K0)
    liqP = phasePore0 == 0
    # If all the pores are filled with liquid, do everything together
    if np.all(liqP):
        rhoPore_kgm30 = Planet.Sil.poreEOS.fn_rho_kgm3(Ppore_MPa0, Tsil_K0)
        kThermPore_WmK[:,0] = Planet.Sil.poreEOS.fn_kTherm_WmK(Ppore_MPa0, Tsil_K0)
        _, KSpore_GPa[:,0] = Planet.Sil.poreEOS.fn_Seismic(Ppore_MPa0, Tsil_K0)
        GSpore_GPa[:,0] = np.zeros_like(KSpore_GPa[:,0])
        DeltaPpore_MPa[:,0] = 1e-6 * rhoPore_kgm30 * gSil_ms20 * (rSil_m0 - rSil_m1)
    else:
        rhoPore_kgm30 = np.zeros_like(Psil_MPa0)
        # For the pores filled with liquid, use ocean liquid EOS
        rhoPore_kgm30[liqP] = Planet.Sil.poreEOS.fn_rho_kgm3(Ppore_MPa0[liqP], Tsil_K0[liqP])
        kThermPore_WmK[liqP,0] = Planet.Sil.poreEOS.fn_kTherm_WmK(Ppore_MPa0[liqP], Tsil_K0[liqP])
        _, KSpore_GPa[liqP,0] = Planet.Sil.poreEOS.fn_Seismic(Ppore_MPa0[liqP], Tsil_K0[liqP])
        GSpore_GPa[liqP,0] = np.zeros_like(KSpore_GPa[liqP,0])
        DeltaPpore_MPa[liqP,0] = 1e-6 * rhoPore_kgm30[liqP] * gSil_ms20[liqP] * (rSil_m0[liqP] - rSil_m1[liqP])
        phases = np.unique(phasePore0)
        # For each HP ice phase represented, perform the same calculations
        for phase in phases[phases != 0]:
            if phase == 1:
                # We only get here if the MoI is consistent with a fully frozen ocean
                thisIceEOS = Planet.Ocean.surfIceEOS['Ih']
            else:
                thisIceEOS = Planet.Ocean.iceEOS[PhaseConv(phase)]
            iP = np.where(phasePore0 == phase)[0]
            rhoPore_kgm30[iP] = thisIceEOS.fn_rho_kgm3(Ppore_MPa0[iP], Tsil_K0[iP])
            kThermPore_WmK[iP,0] = thisIceEOS.fn_kTherm_WmK(Ppore_MPa0[iP], Tsil_K0[iP])
            GSpore_GPa[iP,0], KSpore_GPa[iP,0], _, _ = thisIceEOS.fn_Seismic(Ppore_MPa0[iP], Tsil_K0[iP])
            # Keep DeltaPpore calculation separate for liquid and ice in case of future updates to adjust matrix coupling/pressure scaling
            DeltaPpore_MPa[iP,0] = 1e-6 * rhoPore_kgm30[iP] * gSil_ms20[iP] * (rSil_m0[iP] - rSil_m1[iP])
    # Combine properties using rules defined in EOS functions for porosity
    rhoTot_kgm3[:,0] = Planet.Sil.EOS.fn_porosCorrect(rhoSil_kgm30, rhoPore_kgm30, phiSil_frac0, Planet.Sil.Jrho)
    kThermTot_WmK[:,0] = Planet.Sil.EOS.fn_porosCorrect(kThermSil_WmK0, kThermPore_WmK[:,0], phiSil_frac0, Planet.Sil.JkTherm)
    KStot_GPa[:,0] = Planet.Sil.EOS.fn_porosCorrect(KSsil_GPa0, KSpore_GPa[:,0], phiSil_frac0, Planet.Sil.JKS)
    GStot_GPa[:,0] = Planet.Sil.EOS.fn_porosCorrect(GSsil_GPa0, GSpore_GPa[:,0], phiSil_frac0, Planet.Sil.JGS)

    MLayerSil_kg0 = rhoTot_kgm3[:,0] * 4/3*np.pi*(rSil_m0**3 - rSil_m1**3)
    HtidalSil_Wm30 = Planet.Sil.fn_Htidal_Wm3(rhoTot_kgm3[:,0], gSil_ms20, KStot_GPa[:,0], GStot_GPa[:,0])

    return kThermPore_WmK, KSpore_GPa, GSpore_GPa, DeltaPpore_MPa, \
           kThermTot_WmK, KStot_GPa, GStot_GPa, rhoTot_kgm3, \
           MLayerSil_kg0, HtidalSil_Wm30, rhoPore_kgm30, phasePore0, \
           Ppore_MPa0, phiSil_frac0


def SilRecursionSolid(Planet, Params,
    Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg,
    gSil_ms2, HtidalSil_Wm3, KSsil_GPa, GSsil_GPa, qTop_Wm2, fn_g_ms2):
    """ Propagate a conductive profile down through all silicate profiles.
        In this case, we have different "knowns" than in the ice shell. There,
        we assumed the bottom melting temperature, which gave us the bottom pressure
        from the melting curve. Here, we know the radii from propagating the ice shell
        and ocean down from the surface, so we use a linear radial profile now instead
        of a linear pressure profile. This function assumes no porosity in the rocks.

        Args: See PropagateConductionProfilesSolid
        Returns:
            Psil_MPa, Tsil_K, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2,
                HtidalSil_Wm3, kThermSil_WmK (float, shape (nProfiles, Planet.Steps.nSilMax)):
                Layer properties for silicates.
    """

    # Start propagation at index 1, because we already calculated the 0-index values to
    # get us started here.
    for j in range(1, Planet.Steps.nSilMax):
        # Increment overlying mass using initialization calculation
        MAboveSil_kg[:,j] = MAboveSil_kg[:,j-1] + MLayerSil_kg[:,j-1]
        # Step pressure according to the local layer gravity and overlying mass increase
        Psil_MPa[:,j] = Psil_MPa[:,j-1] + 1e-6 * MLayerSil_kg[:,j-1] * gSil_ms2[:,j-1] / (4*np.pi*rSil_m[:,j]**2)
        # Use Fourier's law and the heat flux consistent with that through the ice shell and ocean,
        # corrected for the difference in radius and adding radiogenic + volumetric heating
        # to determine the temperature change across the layer
        Tsil_K[:,j], qTop_Wm2 = ConductiveTemperature(Tsil_K[:,j-1], rSil_m[:,j-1], rSil_m[:,j],
                    kThermSil_WmK[:,j-1], rhoSil_kgm3[:,j-1], Planet.Sil.Qrad_Wkg, HtidalSil_Wm3[:,j-1],
                    qTop_Wm2)
        rhoSil_kgm3[:,j] = Planet.Sil.EOS.fn_rho_kgm3(Psil_MPa[:,j], Tsil_K[:,j])
        kThermSil_WmK[:,j] = Planet.Sil.EOS.fn_kTherm_WmK(Psil_MPa[:,j], Tsil_K[:,j])
        # Get KS and GS now as they are needed for Htidal calculation;
        # we will calculate them again later along with other seismic calcs
        KSsil_GPa[:,j] = Planet.Sil.EOS.fn_KS_GPa(Psil_MPa[:,j], Tsil_K[:,j])
        GSsil_GPa[:,j] = Planet.Sil.EOS.fn_GS_GPa(Psil_MPa[:,j], Tsil_K[:,j])
        # Calculate gravity using absolute values, as we will use MAboveSil to check for exceeding body mass later.
        gSil_ms2[:,j] = fn_g_ms2(MAboveSil_kg[:,j], rSil_m[:,j])

        MLayerSil_kg[:,j] = rhoSil_kgm3[:,j] * 4/3*np.pi*(rSil_m[:,j]**3 - rSil_m[:,j+1]**3)
        HtidalSil_Wm3[:,j] = Planet.Sil.fn_Htidal_Wm3(rhoSil_kgm3[:,j], gSil_ms2[:,j], KSsil_GPa[:,j], GSsil_GPa[:,j])

    return Psil_MPa, Tsil_K, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, \
           HtidalSil_Wm3, kThermSil_WmK


def SilRecursionPorous(Planet, Params,
    Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, kThermSil_WmK, MLayerSil_kg, MAboveSil_kg,
    gSil_ms2, HtidalSil_Wm3, qTop_Wm2, fn_g_ms2,
    KSsil_GPa, GSsil_GPa, phiSil_frac, Ppore_MPa, phasePore,
    rhoPore_kgm3, kThermPore_WmK, KSpore_GPa, GSpore_GPa, DeltaPpore_MPa,
    rhoTot_kgm3, kThermTot_WmK, KStot_GPa, GStot_GPa):
    """ This is a far-more-complicated version of SilRecursionSolid that accounts
        for porosity in the rock, assuming pores are entirely filled with ocean
        fluid. The ocean EOS is queried using the pore pressure and layer temp to
        determine the phase of the pore-filling material before the combined
        properties are calculated.
    """

    # Start at index 1 because we already did index 0 to get started here.
    for j in range(1, Planet.Steps.nSilMax):
        # Increment overlying mass
        MAboveSil_kg[:,j] = MAboveSil_kg[:,j-1] + MLayerSil_kg[:,j-1]
        # Increment pressure based on layer properties of next layer up and size difference
        Psil_MPa[:,j] = Psil_MPa[:,j-1] + 1e-6 * MLayerSil_kg[:,j-1] * gSil_ms2[:,j-1] / (4*np.pi*rSil_m[:,j]**2)
        # Apply Fourier's law using heat flux out the top, radiogenic, and tidal heating to find
        # temperature change and heat flux into the bottom of the layer
        Tsil_K[:,j], qTop_Wm2 = ConductiveTemperature(Tsil_K[:,j-1], rSil_m[:,j-1], rSil_m[:,j],
                    kThermSil_WmK[:,j-1], rhoSil_kgm3[:,j-1], Planet.Sil.Qrad_Wkg, HtidalSil_Wm3[:,j-1],
                    qTop_Wm2)
        # Get matrix material physical properties
        rhoSil_kgm3[:,j] = Planet.Sil.EOS.fn_rho_kgm3(Psil_MPa[:,j], Tsil_K[:,j])
        kThermSil_WmK[:,j] = Planet.Sil.EOS.fn_kTherm_WmK(Psil_MPa[:,j], Tsil_K[:,j])
        # Get KS and GS now as they are needed for Htidal calculation;
        # we will calculate them again later along with other seismic calcs
        KSsil_GPa[:,j] = Planet.Sil.EOS.fn_KS_GPa(Psil_MPa[:,j], Tsil_K[:,j])
        GSsil_GPa[:,j] = Planet.Sil.EOS.fn_GS_GPa(Psil_MPa[:,j], Tsil_K[:,j])
        # Calculate gravity using absolute values, as we will use MAboveSil to check for exceeding body mass later.
        gSil_ms2[:,j] = fn_g_ms2(MAboveSil_kg[:,j], rSil_m[:,j])

        # Adjust for properties of pore material
        # First, increment pore pressure, assuming hydrostatic pressure in pores, i.e.
        # pores are fully connected and rock matrix supports them so that only overlying
        # ocean pressure increases the pore pressure
        Ppore_MPa[:,j] = Ppore_MPa[:,j-1] + DeltaPpore_MPa[:,j-1]
        # Get effective pore closure pressure using coupling constant (alpha)
        Peff_MPa = Psil_MPa[:,j] - Planet.Sil.alphaPeff * Ppore_MPa[:,j]
        # Get porosity of rock matrix based on these considerations
        phiSil_frac[:,j] = Planet.Sil.fn_phi_frac(Peff_MPa, Tsil_K[:,j])
        # Get pore fluid properties
        # First check for HP ice phases
        phasePore[:,j] = Planet.Sil.poreEOS.fn_phase(Ppore_MPa[:,j], Tsil_K[:,j])
        liqP = phasePore[:,j] == 0
        # If all the pores are filled with liquid, do everything together
        if np.all(liqP):
            rhoPore_kgm3[:,j] = Planet.Sil.poreEOS.fn_rho_kgm3(Ppore_MPa[:,j], Tsil_K[:,j])
            kThermPore_WmK[:,j] = Planet.Sil.poreEOS.fn_kTherm_WmK(Ppore_MPa[:,j], Tsil_K[:,j])
            _, KSpore_GPa[:,j] = Planet.Sil.poreEOS.fn_Seismic(Ppore_MPa[:,j], Tsil_K[:,j])
            GSpore_GPa[:,j] = np.zeros_like(KSpore_GPa[:,j])
            DeltaPpore_MPa[:,j] = 1e-6 * rhoPore_kgm3[:,j] * gSil_ms2[:,j] * (rSil_m[:,j] - rSil_m[:,j+1])
        else:
            # For the pores filled with liquid, use ocean liquid EOS
            rhoPore_kgm3[liqP,j] = Planet.Sil.poreEOS.fn_rho_kgm3(Ppore_MPa[liqP,j], Tsil_K[liqP,j])
            kThermPore_WmK[liqP,j] = Planet.Sil.poreEOS.fn_kTherm_WmK(Ppore_MPa[liqP,j], Tsil_K[liqP,j])
            _, KSpore_GPa[liqP,j] = Planet.Sil.poreEOS.fn_Seismic(Ppore_MPa[liqP,j], Tsil_K[liqP,j])
            GSpore_GPa[liqP,j] = np.zeros_like(KSpore_GPa[liqP,j])
            DeltaPpore_MPa[liqP,j] = 1e-6 * rhoPore_kgm3[liqP,j] * gSil_ms2[liqP,j] * (rSil_m[liqP,j] - rSil_m[liqP,j+1])
            phases = np.unique(phasePore[:,j])
            # For each HP ice phase represented, perform the same calculations
            for phase in phases[phases != 0]:
                # We never need to load iceEOS for ice Ih, so we need special handling to get
                # pore-filling ice Ih properties
                if phase == 1:
                    # We only get here if the MoI is consistent with a fully frozen ocean
                    thisIceEOS = Planet.Ocean.surfIceEOS['Ih']
                else:
                    thisIceEOS = Planet.Ocean.iceEOS[PhaseConv(phase)]
                # Get indices where this ice phase is present
                iP = np.where(phasePore[:,j] == phase)[0]
                # Evaluate pore material properties for this phase
                rhoPore_kgm3[iP,j] = thisIceEOS.fn_rho_kgm3(Ppore_MPa[iP,j], Tsil_K[iP,j])
                kThermPore_WmK[iP,j] = thisIceEOS.fn_kTherm_WmK(Ppore_MPa[iP,j], Tsil_K[iP,j])
                GSpore_GPa[iP,j], KSpore_GPa[iP,j], _, _ = thisIceEOS.fn_Seismic(Ppore_MPa[iP,j], Tsil_K[iP,j])
                # Keep DeltaPpore calculation separate for liquid and ice in case of future updates to adjust matrix coupling/pressure scaling
                DeltaPpore_MPa[iP,j] = 1e-6 * rhoPore_kgm3[iP,j] * gSil_ms2[iP,j] * (rSil_m[iP,j] - rSil_m[iP,j+1])
        # Combine properties using rule defined in EOS function for porosity as in
        # Yu et al. (2016): http://dx.doi.org/10.1016/j.jrmge.2015.07.004
        rhoTot_kgm3[:,j] = Planet.Sil.EOS.fn_porosCorrect(rhoSil_kgm3[:,j], rhoPore_kgm3[:,j], phiSil_frac[:,j], Planet.Sil.Jrho)
        kThermTot_WmK[:,j] = Planet.Sil.EOS.fn_porosCorrect(kThermSil_WmK[:,j], kThermPore_WmK[:,j], phiSil_frac[:,j], Planet.Sil.JkTherm)
        KStot_GPa[:,j] = Planet.Sil.EOS.fn_porosCorrect(KSsil_GPa[:,j], KSpore_GPa[:,j], phiSil_frac[:,j], Planet.Sil.JKS)
        GStot_GPa[:,j] = Planet.Sil.EOS.fn_porosCorrect(GSsil_GPa[:,j], GSpore_GPa[:,j], phiSil_frac[:,j], Planet.Sil.JGS)

        # Finally, get layer total mass and tidal heating rate
        MLayerSil_kg[:,j] = rhoTot_kgm3[:,j] * 4/3*np.pi*(rSil_m[:,j]**3 - rSil_m[:,j+1]**3)
        HtidalSil_Wm3[:,j] = Planet.Sil.fn_Htidal_Wm3(rhoTot_kgm3[:,j], gSil_ms2[:,j], KStot_GPa[:,j], GStot_GPa[:,j])

    return Psil_MPa, Tsil_K, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, \
           HtidalSil_Wm3, kThermSil_WmK, rhoTot_kgm3, phiSil_frac, kThermTot_WmK, \
           Ppore_MPa, rhoPore_kgm3, phasePore


