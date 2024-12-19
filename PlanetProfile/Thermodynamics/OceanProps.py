import numpy as np
import logging
import scipy.interpolate as spi
from PlanetProfile.Thermodynamics.HydroEOS import GetOceanEOS
from PlanetProfile.Utilities.Indexing import GetPhaseIndices
from PlanetProfile.Utilities.defineStructs import Constants, EOSlist
# Assign logger
log = logging.getLogger('PlanetProfile')
def LiquidOceanPropsCalc(Planet, Params):
    """ Calculate/assign aqueous-specific properties for each layer

        Assigns Planet attributes:
            pHs
            aqueousSpeciesAmount_mol
            aqueousSpecies
    """
    # Only perform calculations if this is a valid profile
    if Planet.Do.VALID:
        # Identify which indices correspond to which phases
        indsLiq, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund, indsV, indsVund, indsVI, indsVIund, \
            indsClath, indsClathWet, indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, indsSilV, indsSilVI, \
            indsFe = GetPhaseIndices(Planet.phase)
        if not Planet.Do.NO_H2O and Planet.Ocean.EOS.key not in EOSlist.loaded.keys():
            POcean_MPa = np.arange(Planet.PfreezeLower_MPa, Planet.Ocean.PHydroMax_MPa, Planet.Ocean.deltaP)
            TOcean_K = np.arange(Planet.Bulk.Tb_K, Planet.Ocean.THydroMax_K, Planet.Ocean.deltaT)
            Planet.Ocean.EOS = GetOceanEOS(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, POcean_MPa, TOcean_K,
                               Planet.Ocean.MgSO4elecType, rhoType=Planet.Ocean.MgSO4rhoType,
                               scalingType=Planet.Ocean.MgSO4scalingType, FORCE_NEW=Params.FORCE_EOS_RECALC,
                               phaseType=Planet.Ocean.phaseType, EXTRAP=Params.EXTRAP_OCEAN,
                               sigmaFixed_Sm=Planet.Ocean.sigmaFixed_Sm)
        if np.size(indsLiq) != 0:
            Planet.Ocean.pHs, Planet.Ocean.aqueousSpeciesAmount_mol, Planet.Ocean.aqueousSpecies = (
                Planet.Ocean.EOS.fn_species(Planet.P_MPa[indsLiq], Planet.T_K[indsLiq]))
        else:
            Planet.Ocean.pHs, Planet.Ocean.aqueousSpeciesAmount_mol, Planet.Ocean.aqueousSpecies = np.nan, np.nan, np.nan
    else:
        Planet.Ocean.pHs, Planet.Ocean.aqueousSpeciesAmount_mol, Planet.Ocean.aqueousSpecies = np.nan, np.nan, np.nan
    return Planet

def WriteLiquidOceanProps(Planet, Params):
    """ Write out liquid ocean property calculations to disk """
    colHeaders = [' P (MPa)'.ljust(24),
                  'T (K)'.ljust(24),
                   'pH'.ljust(24)] + [f'{SpeciesCol.ljust(23)}' for SpeciesCol in Planet.Ocean.aqueousSpecies]
    indsLiq, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund, indsV, indsVund, indsVI, indsVIund, indsClath, indsClathWet, indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, indsSilV, indsSilVI, indsFe = GetPhaseIndices(
        Planet.phase)
    # Write out data from core/mantle trade
    with open(Params.DataFiles.oceanPropsFile, 'w') as f:
        f.write(f'Significant Species in Ocean = ' + f'{", ".join(Planet.Ocean.aqueousSpecies)}' + '\n')
        f.write(' '.join(colHeaders) + '\n')
        for i, idx in enumerate(indsLiq):
            line = ' '.join([
                f'{Planet.P_MPa[idx]:24.17e}',
                f'{Planet.T_K[idx]:24.17e}',
                f'{Planet.Ocean.pHs[i]:24.17e}'])
            for j in range(len(Planet.Ocean.aqueousSpecies)):
                line = line + f'{Planet.Ocean.aqueousSpeciesAmount_mol[i][j]:24.17e}'
            line = line + '\n'
            f.write(line)
    log.info(f'Ocean specific properties saved to file: {Params.DataFiles.oceanPropsFile}')
    return