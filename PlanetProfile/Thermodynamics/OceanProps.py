import numpy as np
import logging
import scipy.interpolate as spi
from PlanetProfile.Thermodynamics.HydroEOS import GetOceanEOS
from PlanetProfile.Utilities.Indexing import GetPhaseIndices
from PlanetProfile.Utilities.defineStructs import Constants, EOSlist
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
    # Only perform calculations if this is a valid profile
    if Planet.Do.VALID:
        # Identify indices of liquid phases
        indsLiq, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund, indsV, indsVund, indsVI, indsVIund, \
            indsClath, indsClathWet, indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, indsSilV, indsSilVI, \
            indsFe = GetPhaseIndices(Planet.phase)

        if not Planet.Do.NO_OCEAN and Planet.Ocean.EOS.key not in EOSlist.loaded.keys():
            POcean_MPa = np.arange(Planet.PfreezeLower_MPa, Planet.Ocean.PHydroMax_MPa, Planet.Ocean.deltaP)
            TOcean_K = np.arange(Planet.Bulk.Tb_K, Planet.Ocean.THydroMax_K, Planet.Ocean.deltaT)
            Planet.Ocean.EOS = GetOceanEOS(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, POcean_MPa, TOcean_K,
                               Planet.Ocean.MgSO4elecType, rhoType=Planet.Ocean.MgSO4rhoType,
                               scalingType=Planet.Ocean.MgSO4scalingType, FORCE_NEW=Params.FORCE_EOS_RECALC,
                               phaseType=Planet.Ocean.phaseType, EXTRAP=Params.EXTRAP_OCEAN,
                               sigmaFixed_Sm=Planet.Ocean.sigmaFixed_Sm)
        # Check if we have liquid phases
        if np.size(indsLiq) != 0:
            # If so, then get pH and speciation of ocean
            Planet.Ocean.Bulk_pHs, Planet.Ocean.aqueousSpeciesAmount_mol, Planet.Ocean.aqueousSpecies = (
                Planet.Ocean.EOS.fn_species(Planet.P_MPa[indsLiq], Planet.T_K[indsLiq]))
            Planet.Ocean.Mean_pH = np.mean(Planet.Ocean.Bulk_pHs)
            if "CustomSolution" in Planet.Ocean.comp and Planet.Ocean.reaction is not None:
                Planet.Ocean.affinity_kJ = Planet.Ocean.EOS.fn_rxn_affinity(Planet.P_MPa[indsLiq], Planet.T_K[indsLiq], Planet.Ocean.reaction, Planet.Ocean.reactionDisequilibriumConcentrations)
                Planet.Ocean.affinityMean_kJ = np.mean(Planet.Ocean.affinity_kJ)
                Planet.Ocean.affinitySeafloor_kJ = Planet.Ocean.affinity_kJ[-1]
            else:
                Planet.Ocean.affinity_kJ = (np.zeros(np.size(indsLiq))) * np.nan
                Planet.Ocean.affinitySeafloor_kJ = np.nan
                Planet.Ocean.affinityMean_kJ = np.nan
                Planet.Ocean.reaction = 'NaN'
                Planet.Ocean.reactionDisequilibriumConcentrations = 'NaN'
        else:
            Planet.Ocean.Bulk_pHs, Planet.Ocean.Mean_pH, Planet.Ocean.aqueousSpeciesAmount_mol, Planet.Ocean.aqueousSpecies, Planet.Ocean.affinity_kJ, Planet.Ocean.affinitySeafloor_kJ, Planet.Ocean.affinityMean_kJ = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
            Planet.Ocean.reaction = 'NaN'
            Planet.Ocean.reactionDisequilibriumConcentrations = 'NaN'
    else:
        Planet.Ocean.Bulk_pHs, Planet.Ocean.Mean_pH, Planet.Ocean.aqueousSpeciesAmount_mol, Planet.Ocean.aqueousSpecies, Planet.Ocean.affinity_kJ, Planet.Ocean.affinitySeafloor_kJ, Planet.Ocean.affinityMean_kJ = np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
        Planet.Ocean.reaction = 'NaN'
        Planet.Ocean.reactionDisequilibriumConcentrations = 'NaN'
    return Planet

def WriteLiquidOceanProps(Planet, Params):
    """ Write out liquid ocean property calculations to disk """
    headerLines = [
        f'Significant Species in Ocean = ' + f'{", ".join(Planet.Ocean.aqueousSpecies)}',
        f'Reaction Considered in Ocean = {Planet.Ocean.reaction}',
        f'Concentration of Reaction Species at Disequilibrium = {Planet.Ocean.reactionDisequilibriumConcentrations}'
        ]

    colHeaders = ([' P (MPa)'.ljust(24),
                  'T (K)'.ljust(24),
                   'Bulk pH'.ljust(24), 'Affinity (kJ)'.ljust(24)]
                  + [f'{SpeciesCol.ljust(23)}' for SpeciesCol in Planet.Ocean.aqueousSpecies])
    indsLiq, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund, indsV, indsVund, indsVI, indsVIund, indsClath, indsClathWet, indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, indsSilV, indsSilVI, indsFe = GetPhaseIndices(
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
                line = line + f'{Planet.Ocean.aqueousSpeciesAmount_mol[j][i]:24.17e}'
            line = line + '\n'
            f.write(line)
    log.info(f'Ocean specific properties saved to file: {Params.DataFiles.oceanPropsFile}')
    return
