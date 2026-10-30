import os
import numpy as np
import logging
from scipy.io import loadmat
from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator, interp1d as Interp1D, griddata as GridData
from PlanetProfile import _ROOT
from PlanetProfile.Utilities.defineStructs import Constants, EOSlist, ConstantPropsStruct
from PlanetProfile.Utilities.DataManip import ResetNearestExtrap, ReturnZeros, EOSwrapper, ReturnConstantSpecies
from PlanetProfile.Utilities.Indexing import PhaseInv
# Assign logger
log = logging.getLogger('PlanetProfile')

def GetConstantEOS(constantProperties, EOStype, EOScomp = None):
    constantEOS = ConstantEOSStruct(constantProperties, EOStype, EOScomp)
    if constantEOS.ALREADY_LOADED:
        log.debug(f'{EOStype}{EOScomp} constant EOS already loaded. Reusing existing EOS.')
        constantEOS = EOSlist.loaded[constantEOS.EOSlabel]

    # Ensure each EOSlabel is included in EOSlist, in case we have reused EOSs with
    # e.g. a smaller range that can reuse the larger-range already-loaded EOS.
    if constantEOS.EOSlabel not in EOSlist.loaded.keys():
        EOSlist.loaded[constantEOS.EOSlabel] = constantEOS
        EOSlist.ranges[constantEOS.EOSlabel] = constantEOS.rangeLabel

    constantEOSwrapper = EOSwrapper(constantEOS.EOSlabel)
    
    return constantEOSwrapper

class ConstantEOSStruct:
    def __init__(self, constantProperties, EOStype, EOScomp = None):
        self.EOSlabel = getConstantEOSLabel(constantProperties, EOStype, EOScomp)

        if self.EOSlabel in EOSlist.loaded.keys():
            self.ALREADY_LOADED = True
        else:
            self.ALREADY_LOADED = False
            self.EOStype = EOStype
            self.comp = 'constant'
            self.EXTRAP = True
            
            if EOStype == 'ice':
                self.phaseStr = EOScomp
                self.phaseID = PhaseInv(EOScomp)
                self.POROUS = False
            elif EOStype == 'ocean':
                self.w_ppt = 0
            elif EOStype == 'inner':
                self.phaseStr = EOScomp
                self.phaseID = PhaseInv(EOScomp)
            else:
                raise ValueError(f'Invalid constant EOStype: {EOStype}')

            # Set valid ranges to large values since we are not using them
            self.Pmin = 0
            self.Pmax = 1e8
            self.Tmin = 0
            self.Tmax = 1e8
            self.deltaT = 1e5
            self.deltaP = 1e5
            self.EOSTmin = self.Tmin
            self.EOSTmax = self.Tmax
            self.propsPmax = self.Pmax
            self.EOSdeltaP = self.deltaP
            self.EOSdeltaT = self.deltaT
            self.rangeLabel = f'{self.Pmin:.2f},{self.Pmax:.2f},{self.deltaP:.2e},' + \
                f'{self.Tmin},{self.Tmax},{self.deltaT}'

            # Set up functions for the constant properties
            self.ufn_rho_kgm3  = returnVal(constantProperties.rho_kgm3)
            self.ufn_Cp_JkgK   = returnVal(constantProperties.Cp_JkgK)
            self.ufn_alpha_pK  = returnVal(constantProperties.alpha_pK)
            self.ufn_kTherm_WmK = returnVal(constantProperties.kTherm_WmK)
            self.ufn_VP_kms    = returnVal(constantProperties.VP_kms)
            self.ufn_VS_kms    = returnVal(constantProperties.VS_kms)
            self.ufn_KS_GPa    = returnVal(constantProperties.KS_GPa)
            self.ufn_GS_GPa    = returnVal(constantProperties.GS_GPa)
            self.ufn_phi_frac  = ReturnZeros(1)
            self.ufn_sigma_Sm  = returnVal(constantProperties.sigma_Sm)
            self.ufn_species   = ReturnConstantSpecies(0, None, None, None)
            # Set up function for the viscosity which may consider a transition temperature or not
            if not isinstance(constantProperties.eta_Pas, list):
                self.ufn_eta_Pas = returnVal(constantProperties.eta_Pas)
            else:
                self.ufn_eta_Pas = returnValWithThreshold(
                    constantProperties.eta_Pas[0],
                    constantProperties.eta_Pas[1],
                    constantProperties.TviscTrans_K)
    
            self.ufn_Seismic = ReturnSeismic(self.ufn_VP_kms, self.ufn_VS_kms, self.ufn_KS_GPa, self.ufn_GS_GPa)

            # Set up phase function for ocean and ice
            if self.EOStype == 'ocean':
                self.ufn_phase = returnStepPhase2D(constantProperties._phase_Pb_MPa, constantProperties._phase_Tb_K)
            else:
                self.ufn_phase = returnVal(self.phaseID)


    def fn_porosCorrect(self, propBulk, propPore, phi, J):
        # Combine pore fluid properties with matrix properties in accordance with
        # Yu et al. (2016): http://dx.doi.org/10.1016/j.jrmge.2015.07.004
        return (propBulk**J * (1 - phi) + propPore**J * phi) ** (1/J)

    def fn_phase(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_phase(P_MPa, T_K, grid=grid)

    def fn_rho_kgm3(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_rho_kgm3(P_MPa, T_K, grid=grid)

    def fn_Cp_JkgK(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_Cp_JkgK(P_MPa, T_K, grid=grid)

    def fn_alpha_pK(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_alpha_pK(P_MPa, T_K, grid=grid)

    def fn_kTherm_WmK(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_kTherm_WmK(P_MPa, T_K, grid=grid)
    def fn_Seismic(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_Seismic(P_MPa, T_K, grid=grid)

    def fn_phi_frac(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_phi_frac(P_MPa, T_K, grid=grid)

    def fn_eta_Pas(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_eta_Pas(P_MPa, T_K, grid=grid)
    
    def fn_sigma_Sm(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_sigma_Sm(P_MPa, T_K, grid=grid)
    
    def fn_species(self, P_MPa, T_K, grid=False, reactionEquation=None):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_species(P_MPa, T_K, grid=grid, reactionEquation=reactionEquation)

    def updateConvectionViscosity(self, etaConv_Pas, Tconv_K):
        if isinstance(self.ufn_eta_Pas, returnValWithThreshold):
            self.ufn_eta_Pas.updateThreshold(etaConv_Pas, Tconv_K)
        else:
            pass

class returnVal:
    def __init__(self, val):
        self.val = val

    def __call__(self, P, T, grid=False):
        if grid:
            P, _ = np.meshgrid(P, T, indexing='ij')
        return (np.ones_like(P) * self.val)


class returnValWithThreshold:
    def __init__(self, val_below, val_above, threshold_K):
        self.val_below   = val_below
        self.val_above   = val_above
        self.threshold_K = threshold_K

    def __call__(self, P, T, grid=False):
        if grid:
            P, T_grid = np.meshgrid(P, T, indexing='ij')
            result = np.where(T_grid < self.threshold_K, self.val_below, self.val_above)
        else:
            result = np.where(T < self.threshold_K, self.val_below, self.val_above)
        return result
    
    def updateThreshold(self, valAbove, threshold_K):
        self.val_above = valAbove
        self.threshold_K = threshold_K


class ReturnMultipleVal:
    def __init__(self, vals):
        self.vals = vals

    def __call__(self, P, T, grid=False):
        if grid:
            P, _ = np.meshgrid(P, T, indexing='ij')
        return tuple(np.ones_like(P) * val for val in self.vals)


class ReturnMultipleValWithThreshold:
    def __init__(self, vals_below, vals_above, threshold_K):
        self.vals_below  = vals_below
        self.vals_above  = vals_above
        self.threshold_K = threshold_K

    def __call__(self, P, T, grid=False):
        if grid:
            P, T_grid = np.meshgrid(P, T, indexing='ij')
            result = np.where(T_grid < self.threshold_K, self.vals_below, self.vals_above)
        else:
            result = np.where(T < self.threshold_K, self.vals_below, self.vals_above)
        return result
    

class ReturnSeismic:
    def __init__(self, VP_kms, VS_kms, KS_GPa, GS_GPa):
        self.VP_kms = VP_kms
        self.VS_kms = VS_kms
        self.KS_GPa = KS_GPa
        self.GS_GPa = GS_GPa

    def __call__(self, P, T, grid=False):
        return (self.VP_kms(P, T, grid=grid), self.VS_kms(P, T, grid=grid),
                self.KS_GPa(P, T, grid=grid), self.GS_GPa(P, T, grid=grid))


class returnStepPhase2D:
    """2D step phase for the constant-property melt EOS.

    Returns Ice Ih (phase=1) at P <= Pb_MPa AND T <= Tb_K, liquid (0) otherwise.
    """

    def __init__(self, Pb_MPa, Tb_K):
        self.Pb_MPa = Pb_MPa
        self.Tb_K   = Tb_K
        
    def __call__(self, P, T, grid=False):
        if grid:
            P, T = np.meshgrid(P, T, indexing='ij')
        return np.where((P <= self.Pb_MPa) & (T <= self.Tb_K), 1, 0).astype(int)

def getConstantEOSLabel(constantProperties, EOStype, EOScomp):
    key = f'constant{EOStype}{EOScomp}_rho{constantProperties.rho_kgm3:.2f}_Cp{constantProperties.Cp_JkgK:.2f}_alpha{constantProperties.alpha_pK:.2f}' + \
            f'_kTherm{constantProperties.kTherm_WmK:.2f}_VP{constantProperties.VP_kms:.2f}_VS{constantProperties.VS_kms:.2f}_KS{constantProperties.KS_GPa:.2f}_GS{constantProperties.GS_GPa:.2f}' + \
            f'_eta{constantProperties.eta_Pas}_TviscTrans{constantProperties.TviscTrans_K}'
    if EOStype == 'ocean':
        key += f'_Pb{constantProperties._phase_Pb_MPa:.2f}_Tb{constantProperties._phase_Tb_K:.2f}'
    return key