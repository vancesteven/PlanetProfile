import numpy as np
import logging

# Assign logger
log = logging.getLogger('PlanetProfile')

def IronCoreLayers(Planet, Params,
                   indsSilValid, nSilTooBig, nProfiles, Psil_MPa, Tsil_K, rSil_m, MAboveSil_kg, gSil_ms2):
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
    Pcore_MPa, Tcore_K, rhoCore_kgm3, MLayerCore_kg, gCore_ms2, CpCore_JkgK, alphaCore_pK, kThermCore_WmK = \
        (np.zeros((nProfiles-nSilTooBig, Planet.Steps.nCore)) for _ in range(8))
    rCore_m = np.zeros((nProfiles-nSilTooBig, Planet.Steps.nCore+1))
    # Initialize matching indices as -1 as a flag for unfilled values
    iCoreMatch, nSilFinal = (-1 * np.ones(nProfiles).astype(np.int_) for _ in range(2))

    if Planet.Do.CONSTANT_INNER_DENSITY:
        iCoreStart = [-1]
        silEnd = 0
        indsSilValid = [0]
    else:
        # Calculate maximum core size based on minimum plausible density setting
        MCore_kg = Planet.Bulk.M_kg - MAboveSil_kg[indsSilValid,:]
        rCoreMax_m = (MCore_kg/Planet.Core.rhoMin_kgm3 * 3/4/np.pi)**(1/3)
        # Find first silicate layer smaller than the max core radius
        iCoreStart = [next(j[0] for j,val in np.ndenumerate(rSil_m[i,:-1]) if val < rCoreMax_m[i-nSilTooBig,j])
                      for i in indsSilValid]
        silEnd = Planet.Steps.nSilMax

    log.debug('Evaluating core EOS for possible configurations...')
    for iValid in range(nProfiles - nSilTooBig):
        # Get profile index among silicate layers from number of profiles
        iProf = indsSilValid[iValid]
        # Get index for which silicate layer gets replaced by core layers there and below
        thisCoreStart = iCoreStart[iValid]
        # Get number of remaining silicate layers to iterate over for possible core configs
        nSilRemain = silEnd - thisCoreStart
        # (Re-)initialize placeholder arrays for each core possibility (they change length for each)
        thisPcore_MPa, thisTcore_K, thisrhoCore_kgm3, thisMLayerCore_kg, thisgCore_ms2, thisCpCore_JkgK, \
        thisalphaCore_pK = (np.zeros((nSilRemain, Planet.Steps.nCore)) for _ in range(7))

        # Set starting core values for all possibilities to be equal to silicates at this transition radius
        thisrCore_m = np.array([np.linspace(rSil_m[iProf,thisCoreStart+j], 0, Planet.Steps.nCore+1) for j in range(nSilRemain)])
        thisPcore_MPa[:,0] = [Psil_MPa[iProf,thisCoreStart+j] for j in range(nSilRemain)]
        thisTcore_K[:,0] = [Tsil_K[iProf,thisCoreStart+j] for j in range(nSilRemain)]
        thisrhoCore_kgm3[:,0] = Planet.Core.EOS.fn_rho_kgm3(thisPcore_MPa[:nSilRemain,0], thisTcore_K[:nSilRemain,0])
        thisCpCore_JkgK[:,0] = Planet.Core.EOS.fn_Cp_JkgK(thisPcore_MPa[:nSilRemain,0], thisTcore_K[:nSilRemain,0])
        thisalphaCore_pK[:,0] = Planet.Core.EOS.fn_alpha_pK(thisPcore_MPa[:nSilRemain,0], thisTcore_K[:nSilRemain,0])
        thisMLayerCore_kg[:,0] = [thisrhoCore_kgm3[j,0] * 4/3*np.pi*(thisrCore_m[j,0]**3 - thisrCore_m[j,1]**3) for j in range(nSilRemain)]
        thisgCore_ms2[:,0] = [gSil_ms2[iProf,thisCoreStart+j] for j in range(nSilRemain)]
        MAbove_kg = np.array([MAboveSil_kg[iProf,thisCoreStart+j] for j in range(nSilRemain)])

        for k in range(1, Planet.Steps.nCore):
            MAbove_kg += thisMLayerCore_kg[:,k-1]
            thisDeltaP = 1e-6 * thisMLayerCore_kg[:,k-1] * thisgCore_ms2[:,k-1] / (4*np.pi*thisrCore_m[:,k]**2)
            thisPcore_MPa[:,k] = thisPcore_MPa[:,k-1] + thisDeltaP
            thisTcore_K[:,k] = thisTcore_K[:,k-1] + thisalphaCore_pK[:,k-1]*thisTcore_K[:,k] / \
                           thisCpCore_JkgK[:,k-1] / thisrhoCore_kgm3[:,k-1] * thisDeltaP*1e6
            thisrhoCore_kgm3[:,k] = Planet.Core.EOS.fn_rho_kgm3(thisPcore_MPa[:nSilRemain,k], thisTcore_K[:nSilRemain,k])
            thisCpCore_JkgK[:,k] = Planet.Core.EOS.fn_Cp_JkgK(thisPcore_MPa[:nSilRemain,k], thisTcore_K[:nSilRemain,k])
            thisalphaCore_pK[:,k] = Planet.Core.EOS.fn_alpha_pK(thisPcore_MPa[:nSilRemain,k], thisTcore_K[:nSilRemain,k])
            thisMLayerCore_kg[:,k] = thisrhoCore_kgm3[:,k] * 4/3*np.pi*(thisrCore_m[:,k]**3 - thisrCore_m[:,k+1]**3)
            # Approximate gravity as linear to avoid blowing up for total mass less than body mass (accurate for constant density only)
            thisgCore_ms2[:,k] = thisgCore_ms2[:,0] * thisrCore_m[:,k] / thisrCore_m[:,0]

        if not Planet.Do.CONSTANT_INNER_DENSITY:
            # Find the first core profile that has a mass just below the body mass
            Mtot_kg = MAbove_kg + thisMLayerCore_kg[:,-1]
            iCoreMatch[iProf] = next(ii[0] for ii,val in np.ndenumerate(Mtot_kg) if val < Planet.Bulk.M_kg)
            nSilFinal[iProf] = iCoreStart[iValid] + iCoreMatch[iProf]
            log.debug(f'Core match for iProf = {iProf:d} with Steps.nSil = {nSilFinal[iProf]:d} ' +
                      f'and M = {Mtot_kg[iCoreMatch[iProf]]/Planet.Bulk.M_kg:.4f} M_{Planet.name[0]}.')

        # Assign the values for the core profile with matching total mass to output arrays
        Pcore_MPa[iValid,:] = thisPcore_MPa[iCoreMatch[iProf],:]
        Tcore_K[iValid,:] = thisTcore_K[iCoreMatch[iProf],:]
        rCore_m[iValid,:] = thisrCore_m[iCoreMatch[iProf],:]
        rhoCore_kgm3[iValid,:] = thisrhoCore_kgm3[iCoreMatch[iProf],:]
        MLayerCore_kg[iValid,:] = thisMLayerCore_kg[iCoreMatch[iProf],:]
        gCore_ms2[iValid,:] = thisgCore_ms2[iCoreMatch[iProf],:]
        CpCore_JkgK[iValid,:] = thisCpCore_JkgK[iCoreMatch[iProf],:]
        alphaCore_pK[iValid,:] = thisalphaCore_pK[iCoreMatch[iProf],:]
        kThermCore_WmK[iValid,:] = Planet.Core.EOS.fn_kTherm_WmK(Pcore_MPa[iValid,:], Tcore_K[iValid,:])

    return nSilFinal, Pcore_MPa, Tcore_K, rCore_m, rhoCore_kgm3, MLayerCore_kg, gCore_ms2, CpCore_JkgK, alphaCore_pK, \
        kThermCore_WmK
