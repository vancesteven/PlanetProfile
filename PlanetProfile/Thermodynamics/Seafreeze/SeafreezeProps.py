import numpy as np
import seafreeze.seafreeze as sfz
# Use direct seafreeze import for efficiency
from seafreeze.seafreeze import defpath, _get_tdvs, _is_scatter, phases as seafreeze_phases
from mlbspline import load
from sqlalchemy.sql import false
from PlanetProfile.Utilities.defineStructs import Constants, EOSlist
from itertools import repeat
def IceSeaFreezeProps(PTgrid, phaseName):
    # Get boundarys of P_MPa and T_K grid to generate a unique tag
    sfz_Pmin = np.min(PTgrid[0])
    sfz_Pmax = np.max(PTgrid[0])
    sfz_Tmin = np.min(PTgrid[1])
    sfz_Tmax = np.max(PTgrid[1])
    # Ensure space is linear so that we can use the already loaded EOS if it exists
    PTgrid[0] = np.linspace(sfz_Pmin, sfz_Pmax, PTgrid[0].size)
    PTgrid[1] = np.linspace(sfz_Tmin, sfz_Tmax, PTgrid[1].size)
    if PTgrid[0].size == 1:
        sfz_deltaP = 0
    else:
        sfz_deltaP = np.round(np.mean(np.diff(PTgrid[0])), 3)
    sfz_deltaT = np.round(np.mean(np.diff(PTgrid[1])), 3)
    # Generate unique tags for the seafreeze lookup tables so we only have to generate them once
    seafreezeRange = (sfz_Pmin, sfz_Pmax, sfz_Tmin, sfz_Tmax, sfz_deltaP, sfz_deltaT)
    alreadyLoaded = False
    if phaseName in EOSlist.loaded.keys():
        if EOSlist.ranges[phaseName] == seafreezeRange:
            alreadyLoaded = True
            sfzOut = EOSlist.loaded[phaseName]
        else:
            alreadyLoadedPmin = EOSlist.ranges[phaseName][0]
            alreadyLoadedPmax = EOSlist.ranges[phaseName][1]
            alreadyLoadedTmin = EOSlist.ranges[phaseName][2]
            alreadyLoadedTmax = EOSlist.ranges[phaseName][3]
            alreadyLoadedDeltaP = EOSlist.ranges[phaseName][4]
            alreadyLoadedDeltaT = EOSlist.ranges[phaseName][5]
            # Check if we have already loaded a higher fidelity grid
            noP = sfz_Pmin < alreadyLoadedPmin or sfz_Pmax > alreadyLoadedPmax
            noT = sfz_Tmin < alreadyLoadedTmin or sfz_Tmax > alreadyLoadedTmax
            noDeltaP = sfz_deltaP < alreadyLoadedDeltaP
            noDeltaT = sfz_deltaT < alreadyLoadedDeltaT
            if not noP and not noT and not noDeltaP and not noDeltaT:
                alreadyLoaded = True
                sfzOut = EOSlist.loaded[phaseName]
                PTgrid[0] = np.linspace(alreadyLoadedPmin, alreadyLoadedPmax, sfzOut.G.shape[0])
                PTgrid[1] = np.linspace(alreadyLoadedTmin, alreadyLoadedTmax, sfzOut.G.shape[1])
    
    if not alreadyLoaded:
        sfzOut = sfz.seafreeze(PTgrid, phaseName)
        EOSlist.loaded[phaseName] = sfzOut
        EOSlist.ranges[phaseName] = seafreezeRange
    return sfzOut, PTgrid[0], PTgrid[1]
    
def GenerateSeafreezeChemicalPotentials(P_MPa, T_K, doPureWater = False):
    """
    This function generates the minimum ice chemical potential and the most stable ice phase for a given PT grid.
    It also generates the pure water chemical potential if requested. Functiosn are saved to EOSlist and the tags to access them are returned.
    This is purposely done to avoid recopying large numpy grids, which can be expensive and take up large amounts of temporary memory space, especially for high fidelity grids.
    """
    # Put P_MPa and T_K grid points into format compatible with seafreeze
    evalPts_sfz = np.array([P_MPa, T_K], dtype=object)
    # Get boundarys of P_MPa and T_K grid to generate a unique tag
    sfz_Pmin = np.min(P_MPa)
    sfz_Pmax = np.max(P_MPa)
    sfz_Tmin = np.min(T_K)
    sfz_Tmax = np.max(T_K)
    if P_MPa.size == 1:
        sfz_deltaP = 0
    else:
        sfz_deltaP = np.round(np.mean(np.diff(P_MPa)), 3)
    sfz_deltaT = np.round(np.mean(np.diff(T_K)), 3)
    # Generate unique tags for the seafreeze lookup tables so we only have to generate them once
    seafreezeRange = f'Pmin_{sfz_Pmin}_Pmax_{sfz_Pmax}_Tmin_{sfz_Tmin}_Tmax_{sfz_Tmax}_deltaP_{sfz_deltaP}_deltaT_{sfz_deltaT}'
    seafreezeMuTag = f'mu_J_mol_{seafreezeRange}'
    seafreezePureWaterMuTag = f'mu_J_mol_pure_water_{seafreezeRange}'
    seafreezeIcePhaseTag = f'icePhase_{seafreezeRange}'
    # Check if we have already loaded the minimum ice chemical potential grid and its associated most stable phase grid
    if seafreezeMuTag in EOSlist.loaded.keys() and seafreezeIcePhaseTag in EOSlist.loaded.keys():
        pass
    else:
        # Check if we need to calculate the pure water chemical potential
        calcPureWater = True if doPureWater and seafreezePureWaterMuTag not in EOSlist.loaded.keys() else False
        # Check if we need to calculate the high-pressure ice chemical potentials
        doHPIces = sfz_Pmax > Constants.PminHPices_MPa 
        # Determine number of ice phases to consider
        if doHPIces:
            num_ice_phases = len(Constants.seafreeze_ice_phases)
        else:
            num_ice_phases = 1
        # Create array to store chemical potentials of pure water and for all ice phases
        muWater = np.full((P_MPa.size, T_K.size), np.nan)
        muIces = np.full((P_MPa.size, T_K.size, num_ice_phases), np.nan)
        # Calculate chemical potential for each ice phase
        for phase, name in Constants.seafreeze_ice_phases.items():
            # Skip high-pressure ices if not needed or skip pure water if not needed
            if (phase > 1 and not doHPIces) or (phase == 0 and not calcPureWater):
                continue
            # Calculate chemical potential using seafreeze
            if P_MPa.size == 1 or T_K.size == 1:
                # Handle single value arrays
                sfz_PT = np.array([P_MPa, T_K], dtype=object)
                sfzMu = sfz.getProp(sfz_PT, name).G * Constants.m_gmol['H2O'] / 1000
            else:
                phasedesc = seafreeze_phases[name]
                sp = load.loadSpline(defpath, phasedesc.sp_name)
                isscatter = _is_scatter(evalPts_sfz)
                sfzMu = _get_tdvs(sp, evalPts_sfz, isscatter, 'G').G * Constants.m_gmol['H2O'] / 1000
            # Save pure water chemical potential separately
            if phase == 0:
                muWater = sfzMu
            else:
                # Save chemical potential for other ice phases to 3d array to compare
                muIces[:, :, phase-1] = sfzMu
        
        # Find most stable ice phase and minimum chemical potential
        minIce_mu_Jmol_all = np.full((P_MPa.size, T_K.size), np.nan)
        stableIcePhase = np.zeros((P_MPa.size, T_K.size), dtype=np.uint8)
        if doHPIces:
            # Find minimum chemical potential across all ice phases
            minIce_mu_Jmol_all = np.nanmin(muIces, axis=-1)
            # Only process valid (non-NaN) points
            valid_mask = ~np.isnan(minIce_mu_Jmol_all)
            # Store minimum chemical potential and corresponding phase
            stableIcePhase[valid_mask] = np.nanargmin(muIces[valid_mask], axis=-1) + 1
        else:
            # Only ice Ih is stable at low pressures
            minIce_mu_Jmol_all = muIces[:, :, 0]
            stableIcePhase = np.ones((P_MPa.size, T_K.size), dtype=np.uint8)
        # Save to EOSlist so we do not have to recalculate them next time (good for generating many EOS objects with same PT grid)
        EOSlist.loaded[seafreezeMuTag] = minIce_mu_Jmol_all
        EOSlist.loaded[seafreezeIcePhaseTag] = stableIcePhase
        EOSlist.loaded[seafreezePureWaterMuTag] = muWater
    return seafreezeMuTag, seafreezeIcePhaseTag, seafreezePureWaterMuTag