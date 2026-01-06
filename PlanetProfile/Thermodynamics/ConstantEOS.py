import os
import numpy as np
import logging
from scipy.io import loadmat
from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator, interp1d as Interp1D, griddata as GridData
from PlanetProfile import _ROOT
from PlanetProfile.Utilities.defineStructs import Constants, EOSlist
from PlanetProfile.Utilities.DataManip import ResetNearestExtrap, ReturnZeros, EOSwrapper


class  ConstantEOSStruct:
    def __init__(self, constantProperties, TviscTrans_K=None, EOStype = None, innerComp = None):
        if EOStype == 'inner':
            self.EOStype = 'inner'
            self.comp = innerComp
            self.EOSlabel = f'constant{EOStype}comp{innerComp}constantProps{constantProperties}'
        elif EOStype == 'ocean':
            self.EOStype = 'ocean'
            self.EOSlabel = f'constant{EOStype}constantProps{constantProperties}'
        else:
            raise ValueError(f'Invalid EOStype: {EOStype}')
        self.EOSlabel = f'constant{EOStype}constantProps{constantProperties}'
        if self.EOSlabel in EOSlist.loaded.keys():
            self.ALREADY_LOADED = True
        else:
            self.ALREADY_LOADED = False

        if not self.ALREADY_LOADED:
            self.TviscTrans_K = TviscTrans_K
            self.EXTRAP = True
            self.Pmin, self.EOSPmin = 0, 0
            self.Pmax, self.EOSPmax = 10000, 10000
            self.Tmin, self.EOSTmin = 0, 0
            self.Tmax, self.EOSTmax = 10000, 10000
            
            if EOStype == 'inner':
                if constantProperties['kTherm_WmK'] is None:
                    if self.comp == 'core':
                        constantProperties['kTherm_WmK'] = Constants.kThermFe_WmK
                    elif self.comp == 'sil':
                        constantProperties['kTherm_WmK'] = Constants.kThermSil_WmK
                if constantProperties['eta_Pas'] is None:
                    if self.comp == 'core':
                        constantProperties['eta_Pas'] = Constants.etaFeSolid_Pas
                    elif self.comp == 'sil':
                        constantProperties['eta_Pas'] = Constants.etaRock_Pas
                if TviscTrans_K is None:
                    if self.comp == 'core':
                        self.TviscTrans_K = Constants.TviscFe_K
                    elif self.comp == 'sil':
                        self.TviscTrans_K = Constants.TviscRock_K
                if constantProperties['GS_GPa'] is None:
                    if self.comp == 'core':
                        constantProperties['GS_GPa'] = Constants.GS_GPa[Constants.phaseFe]
                    elif self.comp == 'sil':
                        constantProperties['GS_GPa'] = Constants.GS_GPa[Constants.phaseSil]
                if constantProperties['VP_kms'] is None:
                    if self.comp == 'core':
                        constantProperties['VP_kms'] = np.nan
                    elif self.comp == 'sil':
                        constantProperties['VP_kms'] = np.nan
                if constantProperties['VS_kms'] is None:
                    if self.comp == 'core':
                        constantProperties['VS_kms'] = np.nan
                    elif self.comp == 'sil':
                        constantProperties['VS_kms'] = np.nan
                
                if constantProperties['KS_GPa'] is None:
                    if self.comp == 'core':
                        constantProperties['KS_GPa'] = np.nan
                    elif self.comp == 'sil':
                        constantProperties['KS_GPa'] = np.nan
                if constantProperties['sigma_Sm'] is None:
                    if self.comp == 'core':
                        constantProperties['sigma_Sm'] = np.nan
                    elif self.comp == 'sil':
                        constantProperties['sigma_Sm'] = np.nan
            
            self.ufn_rho_kgm3 = returnVal(constantProperties['rho_kgm3'])
            self.ufn_Cp_JkgK = returnVal(constantProperties['Cp_JkgK'])
            self.ufn_alpha_pK = returnVal(constantProperties['alpha_pK'])
            self.ufn_kTherm_WmK = returnVal(constantProperties['kTherm_WmK'])
            self.ufn_VP_kms = returnVal(constantProperties['VP_kms'])
            self.ufn_VS_kms = returnVal(constantProperties['VS_kms'])
            self.ufn_KS_GPa = returnVal(constantProperties['KS_GPa'])
            self.ufn_GS_GPa = returnVal(constantProperties['GS_GPa'])
            self.ufn_phi_frac = ReturnZeros(1)
            self.ufn_sigma_Sm = returnVal(constantProperties['sigma_Sm'])
            if not isinstance(constantProperties['eta_Pas'], list):
                self.ufn_eta_Pas = returnVal(constantProperties['eta_Pas'])
            else:
                self.ufn_eta_Pas = returnValWithThreshold(constantProperties['eta_Pas'][0], constantProperties['eta_Pas'][1], self.TviscTrans_K)
            self.EOSdeltaP = None
            self.EOSdeltaT = None
            self.propsPmax = 0
            

    
    def fn_porosCorrect(self, propBulk, propPore, phi, J):
        # Combine pore fluid properties with matrix properties in accordance with
        # Yu et al. (2016): http://dx.doi.org/10.1016/j.jrmge.2015.07.004
        return (propBulk**J * (1 - phi) + propPore**J * phi) ** (1/J)
    def fn_phase(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_phase(P_MPa, T_K, grid=grid)
    def fn_rho_kgm3(self, P_MPa, T_K, grid=False):
        # Limit extrapolation to use nearest value from evaluated fit if desired
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
    def fn_VP_kms(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_VP_kms(P_MPa, T_K, grid=grid)
    def fn_VS_kms(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_VS_kms(P_MPa, T_K, grid=grid)
    def fn_KS_GPa(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_KS_GPa(P_MPa, T_K, grid=grid)
    def fn_GS_GPa(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_GS_GPa(P_MPa, T_K, grid=grid)
    def fn_phi_frac(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_phi_frac(P_MPa, T_K, grid=grid)
    def fn_eta_Pas(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            P_MPa, T_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
        return self.ufn_eta_Pas(P_MPa, T_K, grid=grid)

    
    

class returnVal:
    def __init__(self, val):
        self.val = val
    def __call__(self, P, T, grid=False):
        if grid:
            P, _ = np.meshgrid(P, T, indexing='ij')
        return (np.ones_like(P) * self.val)

class returnValWithThreshold:
    def __init__(self, val_below, val_above, threshold_K):
        self.val_below = val_below
        self.val_above = val_above
        self.threshold_K = threshold_K
    
    def __call__(self, P, T, grid=False):
        if grid:
            P, T_grid = np.meshgrid(P, T, indexing='ij')
            result = np.where(T_grid < self.threshold_K, self.val_below, self.val_above)
        else:
            result = np.where(T < self.threshold_K, self.val_below, self.val_above)
        return result

class ReturnMultipleVal:
    def __init__(self, vals):
        self.vals = vals
    def __call__(self, P, T, grid=False):
        if grid:
            P, _ = np.meshgrid(P, T, indexing='ij')
        # Return a tuple of arrays, each filled with the corresponding value from vals
        return tuple(np.ones_like(P) * val for val in self.vals)
class ReturnMultipleValWithThreshold:
    def __init__(self, vals_below, vals_above, threshold_K):
        self.vals_below = vals_below
        self.vals_above = vals_above
        self.threshold_K = threshold_K
    def __call__(self, P, T, grid=False):
        if grid:
            P, T_grid = np.meshgrid(P, T, indexing='ij')
            result = np.where(T_grid < self.threshold_K, self.vals_below, self.vals_above)
        else:
            result = np.where(T < self.threshold_K, self.vals_below, self.vals_above)
        return result