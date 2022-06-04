import numpy as np
import logging
from PlanetProfile.Utilities.defineStructs import Constants
from PlanetProfile.Thermodynamics.Geophysical import PropagateConductionProfilesSolid, \
    PropagateConductionProfilesPorous

# Assign logger
log = logging.getLogger('PlanetProfile')

def SilicateLayers(Planet, Params):
    """ Determines properties of silicate layers based on input Perple_X table
        and seafloor properties, for only non-porous silicates.

        Returns:
            nSilTooBig (int): Number of silicate profiles that have a mass that exceeds the body mass
            nProfiles (int): Number of silicate profiles considered
            Psil_MPa, Tsil_K, rSil_m, rhoSil_kgm3, MLayerSil_kg, MAboveSil_kg, gSil_ms2, phiSil_frac,
                HtidalSil_Wm3, kThermTot_WmK, Ppore_MPa, rhoSil_kgm3, rhoPore_kgm3
                (float, shape nHydroMax-2): State variables corresponding to the silicate layers for
                each physically possible configuration based on the inputs.
            phasePore (int, shape (nProfiles, Planet.Steps.nSil)): Liquid/ice phase of pore material.
                 Later, we truncate this array to shape Planet.Steps.nSil with the MoI- and mass-
                 matching profile.
    """
    if Planet.Do.CONSTANT_INNER_DENSITY or Planet.Do.NO_H2O:
        # If CONSTANT_INNER_DENSITY is True, we have already done the C/MR^2 calculations
        # and now we are just evaluating the EOS for the winning silicate layer set
        # Similarly, if Do.NO_H2O is True, we have 0 ice layers so there is just 1 profile to run
        nProfiles = 1
        if Planet.Do.NO_H2O:
            Planet.Steps.iSilStart = 0
            Planet.Steps.nHydro = 0
        profRange = [Planet.Steps.nHydro - Planet.Steps.iSilStart]
    else:
        nProfiles = Planet.Steps.nSurfIce - Planet.Steps.iSilStart + Planet.Steps.nOceanMax - 1
        profRange = range(nProfiles)

    # Check if we set the core radius to 0, or a found C/MR^2 value (for constant-density approach)
    if Planet.Core.Rset_m is not None:
        rSilEnd_m = Planet.Core.Rset_m
    else:
        rSilEnd_m = 0

    # Finally, we're ready to perform the layer propagation within silicates
    if Planet.Do.POROUS_ROCK:
        Planet, Psil_MPa, Tsil_K, rSil_m, rhoTot_kgm3, MLayerSil_kg, MAboveSil_kg, \
        gSil_ms2, phiSil_frac, HtidalSil_Wm3, kThermTot_WmK, Ppore_MPa, rhoSil_kgm3, \
        rhoPore_kgm3, phasePore \
            = PropagateConductionProfilesPorous(Planet, Params, nProfiles, profRange, rSilEnd_m)
    else:
        Planet, Psil_MPa, Tsil_K, rSil_m, rhoTot_kgm3, MLayerSil_kg, MAboveSil_kg, \
        gSil_ms2, phiSil_frac, HtidalSil_Wm3, kThermTot_WmK, Ppore_MPa, rhoSil_kgm3, \
        rhoPore_kgm3, phasePore \
            = PropagateConductionProfilesSolid(Planet, Params, nProfiles, profRange, rSilEnd_m)

    # Perform validity checks on outputs and package for return
    if Planet.Do.CONSTANT_INNER_DENSITY:
        # Include all indices for later calculations if we already found the desired C/MR^2 match
        indsSilValid = profRange
    else:
        # Get total mass for each possible silicate layer size
        Mtot_kg = MLayerSil_kg[:,-1] + MAboveSil_kg[:,-1]
        # Find silicate radii for which the total mass is too high so we can exclude them
        indsSilValid = np.where(Mtot_kg <= Planet.Bulk.M_kg)[0]
        if Planet.Do.Fe_CORE:
            if(np.size(indsSilValid) == 0):
                msg = 'No silicate mantle size had less than the total body mass.\n' + \
                     f'Min mass: {np.min(Mtot_kg/Planet.Bulk.M_kg):.3f} M_{Planet.name[0]}, ' + \
                     f'max mass: {np.max(Mtot_kg/Planet.Bulk.M_kg):.3f} M_{Planet.name[0]}. '
                suggestion = 'Try adjusting run settings that affect mantle density, like silicate composition ' + \
                             'and heat flux settings.'
                if Params.ALLOW_BROKEN_MODELS:
                    if Params.DO_EXPLOREOGRAM:
                        log.info(msg)
                    else:
                        log.error(msg + suggestion + ' Params.ALLOW_BROKEN_MODELS is True, so calculations ' +
                                  'will proceed with many values set to nan.')
                else:
                    raise RuntimeError(msg + suggestion)
                indsSilValid = range(0)
        elif np.all(Mtot_kg < Planet.Bulk.M_kg):
            Mdiff_frac = 1 - Mtot_kg / Planet.Bulk.M_kg
            MdiffThresh = 0.05
            if np.min(Mdiff_frac) > MdiffThresh:
                raise RuntimeError(f'All masses for some SilicateLayers solutions are more than {100*MdiffThresh:d}% ' +
                                   'less than the total body mass. This likely means the ice shell is too thick to be ' +
                                   'consistent with the value of Steps.iSilStart -- try to increase Bulk.Tb_K or ' +
                                   'decrease Steps.iSilStart.')
            else:
                log.warning(f'All masses for some SilicateLayers solutions are more than {100*MdiffThresh:d}% less ' +
                            f'than the total body mass, though the closest is only {100*np.min(Mdiff_frac):d}% less.')
        else:
            # Record the first entry in the list of models with a total mass lower than the bulk mass --
            # this is *the* match for this Htidal coupling
            if np.size(indsSilValid) != 0:
                indsSilValid = indsSilValid[0]
            # Mark this model as invalid if it has negative temps
            if Tsil_K[indsSilValid, -1] < 0:
                indsSilValid = range(0)

        if np.any(Tsil_K[indsSilValid,:] < 0):
            raise RuntimeError('Negative temperatures encountered in silicates. This likely indicates Qrad_Wkg + Htidal_Wm3 ' +
                           'is too high to be consistent with the heat flow through the ice shell.')

    return indsSilValid, nProfiles, Psil_MPa, Tsil_K, rSil_m, rhoTot_kgm3, \
           MLayerSil_kg, MAboveSil_kg, gSil_ms2, phiSil_frac, HtidalSil_Wm3, kThermTot_WmK, \
           Ppore_MPa, rhoSil_kgm3, rhoPore_kgm3, phasePore
