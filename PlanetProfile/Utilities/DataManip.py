"""
DataManip: Miscellaneous functions for adjusting data that don't belong in
a more specific location.
"""

import numpy as np
from PlanetProfile.Utilities.defineStructs import EOSlist
from PlanetProfile.Utilities.defineStructs import Constants
import logging

# Assign logger
log = logging.getLogger('PlanetProfile')

def ResetNearestExtrap(var1, var2, min1, max1, min2, max2):
    """ Choose the nearest value when EOS function inputs are outside
        the generation domain.

        var1 (float, shape N): First variable, typically P_MPa.
        var2 (float, shape M): First variable, typically T_K.
        min1, max1, min2, max2 (float): Domain boundaries.
    """
    outVar1 = var1 + 0.0
    outVar2 = var2 + 0.0
    if np.size(var1) == 1:
        outVar1 = np.array(var1)
    if np.size(var2) == 1:
        outVar2 = np.array(var2)

    outVar1[var1 < min1] = min1
    outVar1[var1 > max1] = max1
    outVar2[var2 < min2] = min2
    outVar2[var2 > max2] = max2

    return outVar1, outVar2


class ReturnZeros:
    """ Returns an array or tuple of arrays of zeros, for functions of properties
        not modeled that still work with querying routines. We have to run things
        this way and not with a lambda because anonymous functions can't be used
        in parallel processing.
    """
    def __init__(self, nVar):
        self.nVar = nVar

    def __call__(self, P, T, grid=False):
        if grid:
            out = np.zeros((np.size(P), np.size(T)))
        else:
            nPs = np.size(P)
            nTs = np.size(T)
            if self.nVar > 1:
                out = (np.zeros(np.maximum(nPs, nTs)) for _ in range(self.nVar))
            else:
                out = np.zeros(np.maximum(nPs, nTs))
        return out


class ReturnConstantPTw:
    """ Returns an array or tuple of arrays of a constant value, for functions of
        properties not modeled that still work with querying routines. We have to
        run things this way and not with a lambda because anonymous functions
        can't be used in parallel processing.
    """
    def __init__(self, const=None, nVar=1):
        if const is None:
            self.const = 0
        else:
            self.const = const
        self.nVar = nVar

    def __call__(self, P, T, w, grid=False, KEEP_wDIM=False):
        nPs = np.size(P)
        nTs = np.size(T)
        nws = np.size(w)
        if grid or nPs != nTs:
            if KEEP_wDIM:
                out = np.zeros((nPs, nTs, nws)) + self.const
            else:
                out = np.squeeze(np.zeros((nPs, nTs, nws)) + self.const)
        else:
            if self.nVar > 1:
                out = (np.zeros(np.max([nPs, nTs, nws])) for _ in range(self.nVar))
            else:
                out = np.zeros(np.max([nPs, nTs, nws]))
        return out

class ReturnConstantSpecies:
    """ Returns 3 arrays of constant values, for functions of
        properties not modeled that still work with querying routines. We have to
        run things this way and not with a lambda because anonymous functions
        can't be used in parallel processing.
        reference_ppt: Only applicable to Seawater. Gives the reference ppt of the speciation in the associated
        Standard Seawater dictionary and its species
    """
    def __init__(self, ppt, reference_ppt, pH, species):
        if reference_ppt is None:
            speciation_factor = ppt
        else:
            speciation_factor = ppt/reference_ppt # Relative ratio between ppt of HydroEOS and reference ppt in
            # Constant dictionary
        if pH is None:
            self.pH = np.nan
        else:
            self.pH = pH
        if species is None:
            self.species_names = np.array([])
            self.speciation = np.array([])
        else:
            self.species_names = np.array(list(species.keys()))
            self.speciation = np.array(list(species.values())) * speciation_factor
        # Add 55.5mol/kg of water to speciation and species_name
        self.species_names = np.append(self.species_names,'H2O(aq)')
        self.speciation = np.append(self.speciation, 1/Constants.m_gmol['H2O']*1000)

    def __call__(self, P, T, grid = False):
        nPs = np.size(P)
        nTs = np.size(T)
        pH = (np.zeros(nPs)) + self.pH

        speciation = np.zeros((len(self.species_names), nPs))

        # Populate the speciation array
        for row_index, value in enumerate(self.speciation):
            speciation[row_index] = value  # Fill the entire row with the corresponding speciation value
        species_names = np.array(self.species_names)
        return pH, speciation, species_names




class EOSwrapper:
    """ Lightweight wrapper for accessing EOS functions stored in the EOSlist dict. """

    def __init__(self, key):
        self.key = key

        # Assign only those attributes we reference in functions
        if EOSlist.loaded[self.key].EOStype == 'ice':
            self.phaseID = EOSlist.loaded[self.key].phaseID
            self.POROUS = EOSlist.loaded[self.key].POROUS
            self.Tmin = EOSlist.loaded[self.key].Tmin
        elif EOSlist.loaded[self.key].EOStype == 'ocean':
            self.Pmin = EOSlist.loaded[self.key].Pmin
            self.Pmax = EOSlist.loaded[self.key].Pmax
            self.propsPmax = EOSlist.loaded[self.key].propsPmax
            self.Tmin = EOSlist.loaded[self.key].Tmin
            self.Tmax = EOSlist.loaded[self.key].Tmax
            self.deltaP = EOSlist.loaded[self.key].deltaP
            self.deltaT = EOSlist.loaded[self.key].deltaT
            self.EOSdeltaP = EOSlist.loaded[self.key].EOSdeltaP
            self.EOSdeltaT = EOSlist.loaded[self.key].EOSdeltaT
            self.comp = EOSlist.loaded[self.key].comp
            self.w_ppt = EOSlist.loaded[self.key].w_ppt
        elif EOSlist.loaded[self.key].EOStype == 'inner':
            self.comp = EOSlist.loaded[self.key].comp

    def fn_phase(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_phase(P_MPa, T_K, grid=grid)
    def fn_rho_kgm3(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_rho_kgm3(P_MPa, T_K, grid=grid)
    def fn_Cp_JkgK(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_Cp_JkgK(P_MPa, T_K, grid=grid)
    def fn_alpha_pK(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_alpha_pK(P_MPa, T_K, grid=grid)
    def fn_kTherm_WmK(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_kTherm_WmK(P_MPa, T_K, grid=grid)
    def fn_VP_kms(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_VP_kms(P_MPa, T_K, grid=grid)
    def fn_VS_kms(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_VS_kms(P_MPa, T_K, grid=grid)
    def fn_KS_GPa(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_KS_GPa(P_MPa, T_K, grid=grid)
    def fn_GS_GPa(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_GS_GPa(P_MPa, T_K, grid=grid)
    def fn_phi_frac(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_phi_frac(P_MPa, T_K, grid=grid)
    def fn_porosCorrect(self, propBulk, propPore, phi, J):
        return EOSlist.loaded[self.key].fn_porosCorrect(propBulk, propPore, phi, J)
    def fn_sigma_Sm(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_sigma_Sm(P_MPa, T_K, grid=grid)
    def fn_eta_Pas(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_eta_Pas(P_MPa, T_K, grid=grid)
    def fn_Seismic(self, P_MPa, T_K, grid=False):
        return EOSlist.loaded[self.key].fn_Seismic(P_MPa, T_K, grid=grid)
    def fn_species(self, P_MPa, T_K, grid = False):
        return EOSlist.loaded[self.key].fn_species(P_MPa, T_K, grid=grid)
    def fn_rxn_affinity(self, P_MPa, T_K, reaction, concentrations, grid = False):
        return EOSlist.loaded[self.key].fn_rxn_affinity(P_MPa, T_K, reaction, concentrations, grid=grid)
    def fn_averageValuesAccordingtoRule(self, prop1, prop2, rule):
        return EOSlist.loaded[self.key].fn_averageValuesAccordingtoRule(prop1, prop2, rule)
