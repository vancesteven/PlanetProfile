import numpy as np
import logging
import scipy.interpolate as spi
from PlanetProfile.Thermodynamics.HydroEOS import GetOceanEOS
from PlanetProfile.Utilities.Indexing import GetPhaseIndices
from PlanetProfile.Utilities.defineStructs import Constants, EOSlist, Timing
import time
# Assign logger
log = logging.getLogger('PlanetProfile')
def LiquidOceanPropsCalcs(Planet, Params):
    """ Calculate/assign aqueous-specific properties for each layer

        Assigns Planet attributes:
            Ocean.Bulk_pHs
            Ocean.aqueousSpeciesAmount_mol
            Ocean.aqueousSpecies
            Ocean.
    """
    Timing.setFunctionTime(time.time())
    # Only perform calculations if this is a valid profile
    setNaN = False
    if (Planet.Do.VALID or (Params.ALLOW_BROKEN_MODELS and Planet.Do.STILL_CALCULATE_BROKEN_PROPERTIES)) and not Planet.Do.NON_SELF_CONSISTENT:
        if Params.CALC_OCEAN_PROPS:
        
            # Identify indices of liquid phases
            indsLiq, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund, indsV, indsVund, indsVI, indsVIund, \
                indsClath, indsClathWet, indsMixedClathrateIh, indsMixedClathrateII, indsMixedClathrateIII, indsMixedClathrateV, indsMixedClathrateVI, \
                indsMixedClathrateIhwet, indsMixedClathrateIIund, indsMixedClathrateIIIund, indsMixedClathrateVund, indsMixedClathrateVIund, \
                indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, indsSilV, indsSilVI, \
                indsFe = GetPhaseIndices(Planet.phase)

            if not Planet.Do.NO_OCEAN and Planet.Ocean.EOS.key not in EOSlist.loaded.keys():
                POcean_MPa = np.arange(Planet.PfreezeLower_MPa, Planet.Ocean.PHydroMax_MPa, Planet.Ocean.deltaP)
                TOcean_K = np.arange(Planet.Bulk.Tb_K, Planet.Ocean.THydroMax_K, Planet.Ocean.deltaT)
                Planet.Ocean.EOS = GetOceanEOS(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, POcean_MPa, TOcean_K,
                                Planet.Ocean.MgSO4elecType, rhoType=Planet.Ocean.MgSO4rhoType,
                                scalingType=Planet.Ocean.MgSO4scalingType, FORCE_NEW=Params.FORCE_EOS_RECALC,
                                phaseType=Planet.Ocean.phaseType, EXTRAP=Params.EXTRAP_OCEAN,
                                sigmaFixed_Sm=Planet.Ocean.sigmaFixed_Sm, propsStepReductionFactor=Planet.Ocean.propsStepReductionFactor)
            # Check if we have liquid phases
            if np.size(indsLiq) != 0:
                if 'CustomSolution' in Planet.Ocean.comp:
                    # If so, then get pH and speciation of ocean
                    Planet.Ocean.Reaction = setupReactionSubstruct(Planet.Ocean.Reaction)
                    Planet.Ocean.Bulk_pHs, Planet.Ocean.aqueousSpeciesAmount_mol, Planet.Ocean.aqueousSpecies, Planet.Ocean.affinity_kJ = (
                        Planet.Ocean.EOS.fn_species(Planet.P_MPa[indsLiq], Planet.T_K[indsLiq], reactionSubstruct = Planet.Ocean.Reaction))
                else:
                    Planet.Ocean.Bulk_pHs, Planet.Ocean.aqueousSpeciesAmount_mol, Planet.Ocean.aqueousSpecies = (
                        Planet.Ocean.EOS.fn_species(Planet.P_MPa[indsLiq], Planet.T_K[indsLiq]))
                    Planet.Ocean.affinity_kJ = np.repeat(np.nan, len(indsLiq))
                Planet.Ocean.Mean_pH = np.mean(Planet.Ocean.Bulk_pHs)
                Planet.Ocean.pHSeafloor = Planet.Ocean.Bulk_pHs[-1]
                Planet.Ocean.pHTop = Planet.Ocean.Bulk_pHs[0]
                Planet.Ocean.affinitySeafloor_kJ = Planet.Ocean.affinity_kJ[-1]
                Planet.Ocean.affinityTop_kJ = Planet.Ocean.affinity_kJ[0]
                Planet.Ocean.affinityMean_kJ = np.mean(Planet.Ocean.affinity_kJ)
            else:
                setNaN = True
        else:
            setNaN = True
    else:
        setNaN = True
    if setNaN:
        Planet.Ocean.Bulk_pHs, Planet.Ocean.Mean_pH, Planet.Ocean.pHSeafloor, Planet.Ocean.pHTop, \
            Planet.Ocean.aqueousSpeciesAmount_mol, Planet.Ocean.aqueousSpecies, Planet.Ocean.affinity_kJ, \
                Planet.Ocean.affinitySeafloor_kJ, Planet.Ocean.affinityTop_kJ, Planet.Ocean.affinityMean_kJ = np.repeat(np.nan, 10)
    Timing.printFunctionTimeDifference('LiquidOceanPropsCalcs()', time.time())
    return Planet
def setupReactionSubstruct(reactionSubstruct):
    if reactionSubstruct.reaction is None:
        reactionSubstruct.reaction = 'NaN'
    if reactionSubstruct.reaction != 'NaN':
        reactionSubstruct.parsed_reaction = reaction_parser(reactionSubstruct.reaction)
        for species in reactionSubstruct.parsed_reaction["allSpecies"]:
            if reactionSubstruct.useReferenceSpecies:
                reactionSubstruct.disequilibriumConcentrations[species] = None
                referenceSpecies = reactionSubstruct.referenceSpecies
                if reactionSubstruct.useH2ORatio:
                    referenceRatioToH2O = reactionSubstruct.mixingRatioToH2O[referenceSpecies]
                    if species in reactionSubstruct.mixingRatioToH2O.keys():
                        reactionSubstruct.disequilibriumConcentrations[species] = reactionSubstruct.mixingRatioToH2O[species] / referenceRatioToH2O
                else:
                    if species in reactionSubstruct.relativeRatioToReferenceSpecies.keys():
                        reactionSubstruct.disequilibriumConcentrations[species] = reactionSubstruct.disequilibriumConcentrations[referenceSpecies]
    else:
        reactionSubstruct.disequilibriumConcentrations = {}
        reactionSubstruct.useReferenceSpecies = False
        reactionSubstruct.useH2ORatio = False
        reactionSubstruct.referenceSpecies = 'NaN'
    return reactionSubstruct

def reaction_parser(reaction):
    """
        Parse a chemical reaction string into reactants, products, and optional disequilibrium species.

        Parameters:
        reaction_str (str): The chemical reaction string (e.g., "CO2 + 4 H2(aq) -> CH4(aq) + 2 H2O(aq)").

        Returns:
        dict: Parsed reaction with reactants, products, and optional disequilibrium species.
        """
    reaction_parts = reaction.split("->")
    reactants_str, products_str = reaction_parts[0], reaction_parts[1]

    def parse_side(side_str):
        species_dict = {}
        components = side_str.split("+")
        for component in components:
            component = component.strip()
            if " " in component:
                coeff, species = component.split(" ", 1)
                species_dict[species.strip()] = float(coeff)
            else:
                species_dict[component.strip()] = 1.0
        return species_dict

    reactants = parse_side(reactants_str)
    products = parse_side(products_str)


    return {"reactants": reactants, "products": products, "allSpecies": reactants.keys() | products.keys()}
    
def WriteLiquidOceanProps(Planet, Params):
    """ Write out liquid ocean property calculations to disk """
    if not Planet.Ocean.Reaction.useReferenceSpecies:
        reactionSpeciesDescription = 'Reaction Species Relative Ratio at Disequilibrium'
    else:
        reactionSpeciesDescription = f'Reaction Species Relative Ratio to {Planet.Ocean.Reaction.referenceSpecies} at Disequilibrium'
        
    headerLines = [
        f'Significant Species in Ocean = ' + f'{"; ".join(Planet.Ocean.aqueousSpecies)}',
        f'Reaction Considered in Ocean = {Planet.Ocean.Reaction.reaction}',
        f'Use Reference Species = {Planet.Ocean.Reaction.useReferenceSpecies}',
        f'Reference Species = {Planet.Ocean.Reaction.referenceSpecies}',
        f'{reactionSpeciesDescription} = {Planet.Ocean.Reaction.disequilibriumConcentrations}',
        ]

    colHeaders = ([' P (MPa)'.ljust(24),
                  'T (K)'.ljust(24),
                   'Bulk pH'.ljust(24), 'Affinity (kJ)'.ljust(24)]
                  + [f'{SpeciesCol.ljust(23)}' for SpeciesCol in Planet.Ocean.aqueousSpecies])
    indsLiq, indsIceI, indsIceIwet, indsIceII, indsIceIIund, indsIceIII, indsIceIIIund, indsIceV, indsIceVund, \
               indsIceVI, indsIceVIund, indsClath, indsClathWet, indsMixedClathrateIh, indsMixedClathrateII, indsMixedClathrateIII, indsMixedClathrateV, indsMixedClathrateVI, \
               indsMixedClathrateIhwet, indsMixedClathrateIIund, indsMixedClathrateIIIund, indsMixedClathrateVund, indsMixedClathrateVIund, \
               indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, \
               indsSilV, indsSilVI, indsFe = GetPhaseIndices(
        Planet.phase)
    Params.nHeadLines = np.size(headerLines) + 2
    headerLines = np.insert(headerLines, 0, f'nHeadLines = {Params.nHeadLines:d}')
    # Write out data from core/mantle trade
    with open(Params.DataFiles.oceanPropsFile, 'w') as f:
        f.write('\n  '.join(headerLines) + '\n')
        f.write(' '.join(colHeaders) + '\n')
        for i, idx in enumerate(indsLiq):
            line = ' '.join([
                f'{Planet.P_MPa[idx]:24.17e}',
                f'{Planet.T_K[idx]:24.17e}',
                f'{Planet.Ocean.Bulk_pHs[i]:24.17e}',
                f'{Planet.Ocean.affinity_kJ[i]:24.17e}'])
            for j in range(len(Planet.Ocean.aqueousSpecies)):
                line = line + f'{Planet.Ocean.aqueousSpeciesAmount_mol[j][i]:24.17e} '
            line = line + '\n'
            f.write(line)
    log.info(f'Ocean specific properties saved to file: {Params.DataFiles.oceanPropsFile}')
    return
