import numpy as np
import logging
import scipy.interpolate as spi
from PlanetProfile.Thermodynamics.HydroEOS import GetPlanetOceanEOS
from PlanetProfile.Utilities.Indexing import GetPhaseIndices
from PlanetProfile.Utilities.defineStructs import Constants, EOSlist, Timing
from PlanetProfile.Utilities.ResultsIO import ensure_parent_dir
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
    if (Planet.Do.VALID or (Params.ALLOW_BROKEN_MODELS and Planet.Do.STILL_CALCULATE_BROKEN_PROPERTIES)):
        if Params.CALC_OCEAN_PROPS:
            if Planet.Ocean.reactionEquation is None:
                Planet.Ocean.reactionEquation = 'none'
        
            # Identify indices of liquid phases
            indsLiq, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund, indsV, indsVund, indsVI, indsVIund, \
                indsClath, indsClathWet, indsMixedClathrateIh, indsMixedClathrateII, indsMixedClathrateIII, indsMixedClathrateV, indsMixedClathrateVI, \
                indsMixedClathrateIhwet, indsMixedClathrateIIund, indsMixedClathrateIIIund, indsMixedClathrateVund, indsMixedClathrateVIund, \
                indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, indsSilV, indsSilVI, \
                indsFe = GetPhaseIndices(Planet.phase)

            if not Planet.Do.NO_OCEAN and Planet.Ocean.EOS.key not in EOSlist.loaded.keys():
                POcean_MPa = np.arange(Planet.PfreezeLower_MPa, Planet.Ocean.PHydroMax_MPa, Planet.Ocean.deltaP)
                TOcean_K = np.arange(Planet.Bulk.Tb_K, Planet.Ocean.THydroMax_K, Planet.Ocean.deltaT)
                Planet.Ocean.EOS = GetPlanetOceanEOS(Planet, Params, POcean_MPa, TOcean_K)
            # Check if we have liquid phases
            if np.size(indsLiq) != 0:
                if 'CustomSolution' in Planet.Ocean.comp:
                    Planet.Ocean.Bulk_pHs, Planet.Ocean.aqueousSpeciesAmount_mol, Planet.Ocean.aqueousSpecies, Planet.Ocean.equilibriumReactionConstant = (
                        Planet.Ocean.EOS.fn_species(Planet.P_MPa[indsLiq], Planet.T_K[indsLiq], reactionEquation = Planet.Ocean.reactionEquation))
                else:
                    Planet.Ocean.Bulk_pHs, Planet.Ocean.aqueousSpeciesAmount_mol, Planet.Ocean.aqueousSpecies = (
                        Planet.Ocean.EOS.fn_species(Planet.P_MPa[indsLiq], Planet.T_K[indsLiq]))
                    Planet.Ocean.equilibriumReactionConstant = np.repeat(np.nan, len(indsLiq))
                Planet.Ocean.Mean_pH = np.mean(Planet.Ocean.Bulk_pHs)
                Planet.Ocean.pHSeafloor = Planet.Ocean.Bulk_pHs[-1]
                Planet.Ocean.pHTop = Planet.Ocean.Bulk_pHs[0]
            else:
                setNaN = True
        else:
            setNaN = True
    else:
        setNaN = True
    if setNaN:
        Planet.Ocean.reactionEquation = 'none'
        Planet.Ocean.Bulk_pHs, Planet.Ocean.Mean_pH, Planet.Ocean.pHSeafloor, Planet.Ocean.pHTop, \
            Planet.Ocean.aqueousSpeciesAmount_mol, Planet.Ocean.aqueousSpecies, Planet.Ocean.equilibriumReactionConstant = np.repeat(np.nan, 7)
    Timing.printFunctionTimeDifference('LiquidOceanPropsCalcs()', time.time())
    return Planet
    
def WriteLiquidOceanProps(Planet, Params):
    """ Write out liquid ocean property calculations to disk """
    headerLines = [
        f'Significant Species in Ocean = ' + f'{"; ".join(Planet.Ocean.aqueousSpecies)}',
        f'Reaction Considered in Ocean = {Planet.Ocean.reactionEquation}',
        ]

    colHeaders = ([' P (MPa)'.ljust(24),
                  'T (K)'.ljust(24),
                   'Bulk pH'.ljust(24), 'Equilibrium Reaction Constant'.ljust(24)]
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
    ensure_parent_dir(Params.DataFiles.oceanPropsFile)
    with open(Params.DataFiles.oceanPropsFile, 'w') as f:
        f.write('\n  '.join(headerLines) + '\n')
        f.write(' '.join(colHeaders) + '\n')
        for i, idx in enumerate(indsLiq):
            line = ' '.join([
                f'{Planet.P_MPa[idx]:24.17e}',
                f'{Planet.T_K[idx]:24.17e}',
                f'{Planet.Ocean.Bulk_pHs[i]:24.17e}',
                f'{Planet.Ocean.equilibriumReactionConstant[i]:24.17e}'])
            for j in range(len(Planet.Ocean.aqueousSpecies)):
                line = line + f'{Planet.Ocean.aqueousSpeciesAmount_mol[j][i]:24.17e} '
            line = line + '\n'
            f.write(line)
    log.info(f'Ocean specific properties saved to file: {Params.DataFiles.oceanPropsFile}')
    return
