import os
import numpy as np
import logging
from scipy.io import loadmat
from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator, interp1d as Interp1D, griddata as GridData
from PlanetProfile import _ROOT
from PlanetProfile.Utilities.defineStructs import Constants, EOSlist

# Assign logger
log = logging.getLogger('PlanetProfile')

def GetInnerEOS(EOSfname, EOSinterpMethod='nearest', nHeaders=13, Fe_EOS=False, kThermConst_WmK=None,
                HtidalConst_Wm3=0, porosType=None, phiTop_frac=0, Pclosure_MPa=350, phiMin_frac=None,
                EXTRAP=False, wFeCore_ppt=None, wScore_ppt=None):
    innerEOS = PerplexEOSStruct(EOSfname, EOSinterpMethod=EOSinterpMethod, nHeaders=nHeaders,
                                Fe_EOS=Fe_EOS, kThermConst_WmK=kThermConst_WmK,
                                HtidalConst_Wm3=HtidalConst_Wm3, porosType=porosType,
                                phiTop_frac=phiTop_frac, Pclosure_MPa=Pclosure_MPa,
                                phiMin_frac=phiMin_frac, EXTRAP=EXTRAP,
                                wFeCore_ppt=wFeCore_ppt, wScore_ppt=wScore_ppt)
    if innerEOS.ALREADY_LOADED:
        log.debug(f'{innerEOS.comp} EOS already loaded. Reusing existing EOS.')
        innerEOS = EOSlist.loaded[innerEOS.EOSlabel]

    # Ensure each EOSlabel is included in EOSlist, in case we have reused EOSs with
    # e.g. a smaller range that can reuse the larger-range already-loaded EOS.
    if innerEOS.EOSlabel not in EOSlist.loaded.keys():
        EOSlist.loaded[innerEOS.EOSlabel] = innerEOS

    innerEOSwrapper = EOSwrapper(innerEOS.EOSlabel)

    return innerEOSwrapper

class PerplexEOSStruct:
    """ Loads in Perple_X table and creates interpolators that can be called to
        obtain silicate/core properties as functions of P and T.
    """
    def __init__(self, EOSfname, EOSinterpMethod='nearest', nHeaders=13, Fe_EOS=False, kThermConst_WmK=None,
                 HtidalConst_Wm3=0, porosType=None, phiTop_frac=0, Pclosure_MPa=350, phiMin_frac=None,
                 EXTRAP=False, wFeCore_ppt=None, wScore_ppt=None):
        self.comp = EOSfname[:-4]
        self.EOSlabel = f'{self.comp}{EOSinterpMethod}{kThermConst_WmK}{HtidalConst_Wm3}{porosType}' + \
                        f'{phiTop_frac}{Pclosure_MPa}{phiMin_frac}{EXTRAP}{wFeCore_ppt}{wScore_ppt}'
        if self.EOSlabel in EOSlist.loaded.keys():
            self.ALREADY_LOADED = True
        else:
            self.ALREADY_LOADED = False

        if not self.ALREADY_LOADED:
            if Fe_EOS:
                descrip = ' core composition.'
            else:
                descrip = f' silicate composition with Htidal_Wm3 = {HtidalConst_Wm3}, porosType ' + \
                          f'{porosType}, phiTop_frac = {phiTop_frac}, Pclosure_MPa = {Pclosure_MPa}, ' + \
                          f'and EXTRAP = {EXTRAP}.'
            log.debug(f'Loading Perplex EOS for {self.comp}{descrip}')
            self.dir = os.path.join(_ROOT, 'Thermodynamics', 'EOStables', 'Perple_X')
            self.fpath = os.path.join(self.dir, EOSfname)
            self.Fe_EOS = Fe_EOS
            self.EXTRAP = EXTRAP
            self.EOStype = 'inner'
            self.porosType = porosType
            self.phiTop_frac = phiTop_frac
            self.Pclosure_MPa = Pclosure_MPa
            self.phiMin_frac = phiMin_frac
            self.wFeCore_ppt = wFeCore_ppt
            self.wScore_ppt = wScore_ppt

            if '3D_EOS' in self.fpath:
                if wFeCore_ppt is None and wScore_ppt is None:
                    log.warning(f'3D EOS {self.fpath} is being loaded, but wFeCore_ppt and wScore_ppt ' +
                                f'are both None -- at least one must be specified to select the table. ' +
                                f'Default value of {Constants.wFeDef_ppt} will be used.')
                    wFeCore_ppt = Constants.wFeDef_ppt
                elif wFeCore_ppt is not None and wScore_ppt is not None:
                    log.warning('Both wFeCore_ppt and wScore_ppt were specified. wScore_ppt will be ignored.')
                if wFeCore_ppt is None:
                    wFeCore_ppt = 1e3 - wScore_ppt
                tableKey = f'{self.fpath}_wFe{wFeCore_ppt}'
            else:
                tableKey = self.fpath

            if tableKey in EOSlist.loaded.keys():
                log.debug('Reloading previously imported Perplex table.')
                P_MPa, T_K, self.ufn_rho_kgm3, self.ufn_VP_kms, self.ufn_VS_kms, \
                self.ufn_KS_GPa, self.ufn_GS_GPa, self.ufn_Cp_JkgK, self.ufn_alpha_pK \
                    = EOSlist.loaded[tableKey]
                self.Pmin, self.Pmax, self.Tmin, self.Tmax, self.deltaP, self.deltaT \
                    = EOSlist.ranges[tableKey]
            else:
                if '3D_EOS' in self.fpath:
                    EOS3D = loadmat(self.fpath)
                    wFe_ppt = EOS3D['wFe_ppt'][0]
                    P1D_MPa = EOS3D['P_MPa'][0]
                    T1D_K = EOS3D['T_K'][0]
                    intPts = (wFe_ppt, P1D_MPa, T1D_K)
                    Peval_MPa, Teval_K = np.meshgrid(P1D_MPa, T1D_K, indexing='ij')
                    Peval_MPa = Peval_MPa.flatten()
                    Teval_K = Teval_K.flatten()
                    evalPts = np.array([wFeCore_ppt * np.ones_like(Peval_MPa), Peval_MPa, Teval_K]).T
                    rho_kgm3 = RegularGridInterpolator(intPts, EOS3D['rho_kgm3'], method='linear')(evalPts)
                    VP_kms = RegularGridInterpolator(intPts, EOS3D['VP_kms'], method='linear')(evalPts)
                    VS_kms = RegularGridInterpolator(intPts, EOS3D['VS_kms'], method='linear')(evalPts)
                    Cp_JkgK = RegularGridInterpolator(intPts, EOS3D['Cp_JkgK'], method='linear')(evalPts)
                    alpha_pK = RegularGridInterpolator(intPts, EOS3D['alpha_pK'], method='linear')(evalPts)
                    KS_GPa = RegularGridInterpolator(intPts, EOS3D['KS_GPa'], method='linear')(evalPts)
                    GS_GPa = RegularGridInterpolator(intPts, EOS3D['GS_GPa'], method='linear')(evalPts)

                    tableShape = (np.size(P1D_MPa), -1)
                    rho_kgm3 = np.reshape(rho_kgm3, tableShape)
                    VP_kms = np.reshape(VP_kms, tableShape)
                    VS_kms = np.reshape(VS_kms, tableShape)
                    Cp_JkgK = np.reshape(Cp_JkgK, tableShape)
                    alpha_pK = np.reshape(alpha_pK, tableShape)
                    KS_GPa = np.reshape(KS_GPa, tableShape)
                    GS_GPa = np.reshape(GS_GPa, tableShape)

                    self.Pmin = np.min(P1D_MPa)
                    self.Pmax = np.max(P1D_MPa)
                    self.Tmin = np.min(T1D_K)
                    self.Tmax = np.max(T1D_K)
                    self.deltaP = P1D_MPa[1] - P1D_MPa[0]
                    self.deltaT = T1D_K[1] - T1D_K[0]
                else:
                    # Load in Perple_X data. Note that all P, KS, and GS are stored as bar
                    firstPT, secondPT, rho_kgm3, VP_kms, VS_kms, Cp_Jm3K, alpha_pK, KS_bar, GS_bar \
                        = np.loadtxt(self.fpath, skiprows=nHeaders, unpack=True)
                    # We don't know yet whether P or T is in the first column, which increments the fastest.
                    # The second column increments once for each time the first column runs through the whole range.
                    dim2 = next(i[0] for i, val in np.ndenumerate(secondPT) if val != secondPT[0])
                    # Check the column header to see if P or T is printed first
                    with open(self.fpath) as f:
                        [f.readline() for _ in range(nHeaders-1)]
                        colHeaderLine = f.readline().strip()
                    # Get necessary values to generate P, T arrays
                    if colHeaderLine[0] == 'P':
                        P_FIRST = True
                        Plin_MPa = firstPT * Constants.bar2MPa  # Pressures saved as bar
                        Tlin_K = secondPT
                        lenP = dim2
                    elif colHeaderLine[0] == 'T':
                        P_FIRST = False
                        Plin_MPa = secondPT * Constants.bar2MPa
                        Tlin_K = firstPT
                        lenP = int(len(Plin_MPa)/dim2)
                    else:
                        raise ValueError(f'Perple_X table {EOSfname} does not have T or P in the first column.')

                    # Make 1D arrays of P and T that just span the axes (unlike the 2D meshes of the dependent variables)
                    self.Pmin = Plin_MPa[0]
                    self.Pmax = Plin_MPa[-1]
                    self.Tmin = Tlin_K[0]
                    self.Tmax = Tlin_K[-1]
                    lenT = int(len(Tlin_K) / lenP)
                    P1D_MPa, self.deltaP = np.linspace(self.Pmin, self.Pmax, lenP, retstep=True)
                    T1D_K, self.deltaT = np.linspace(self.Tmin, self.Tmax, lenT, retstep=True)

                    # Set unphysical values to NaN so they will be caught by the next step (for all but alpha, which can cross zero)
                    rho_kgm3[rho_kgm3 <= 0] = np.nan
                    VP_kms[VP_kms <= 0] = np.nan
                    VS_kms[VS_kms <= 0] = np.nan
                    Cp_Jm3K[Cp_Jm3K <= 0] = np.nan
                    KS_bar[KS_bar <= 0] = np.nan
                    GS_bar[GS_bar <= 0] = np.nan

                    # Interpolate dependent variables where the values are NaN
                    PTpts = (Plin_MPa, Tlin_K)
                    thisVarValid = np.isfinite(rho_kgm3)
                    if not np.all(thisVarValid):
                        thisValidPs = Plin_MPa[np.where(thisVarValid)]
                        thisValidTs = Tlin_K[np.where(thisVarValid)]
                        rho_kgm3 = GridData((thisValidPs, thisValidTs), rho_kgm3[thisVarValid], PTpts, method=EOSinterpMethod)
                    thisVarValid = np.isfinite(VP_kms)
                    if not np.all(thisVarValid):
                        thisValidPs = Plin_MPa[np.where(thisVarValid)]
                        thisValidTs = Tlin_K[np.where(thisVarValid)]
                        VP_kms = GridData((thisValidPs, thisValidTs), VP_kms[thisVarValid], PTpts, method=EOSinterpMethod)
                    thisVarValid = np.isfinite(VS_kms)
                    if not np.all(thisVarValid):
                        thisValidPs = Plin_MPa[np.where(thisVarValid)]
                        thisValidTs = Tlin_K[np.where(thisVarValid)]
                        VS_kms = GridData((thisValidPs, thisValidTs), VS_kms[thisVarValid], PTpts, method=EOSinterpMethod)
                    thisVarValid = np.isfinite(Cp_Jm3K)
                    if not np.all(thisVarValid):
                        thisValidPs = Plin_MPa[np.where(thisVarValid)]
                        thisValidTs = Tlin_K[np.where(thisVarValid)]
                        # Note that Perple_X tables store heat capacities as J/m^3/K, so we convert to J/kg/K by dividing by rho in a moment
                        Cp_Jm3K = GridData((thisValidPs, thisValidTs), Cp_Jm3K[thisVarValid], PTpts, method=EOSinterpMethod)
                    thisVarValid = np.isfinite(alpha_pK)
                    if not np.all(thisVarValid):
                        thisValidPs = Plin_MPa[np.where(thisVarValid)]
                        thisValidTs = Tlin_K[np.where(thisVarValid)]
                        alpha_pK = GridData((thisValidPs, thisValidTs), alpha_pK[thisVarValid], PTpts, method=EOSinterpMethod)
                    thisVarValid = np.isfinite(KS_bar)
                    if not np.all(thisVarValid):
                        thisValidPs = Plin_MPa[np.where(thisVarValid)]
                        thisValidTs = Tlin_K[np.where(thisVarValid)]
                        KS_bar = GridData((thisValidPs, thisValidTs), KS_bar[thisVarValid], PTpts, method=EOSinterpMethod)
                    thisVarValid = np.isfinite(GS_bar)
                    if not np.all(thisVarValid):
                        thisValidPs = Plin_MPa[np.where(thisVarValid)]
                        thisValidTs = Tlin_K[np.where(thisVarValid)]
                        GS_bar = GridData((thisValidPs, thisValidTs), GS_bar[thisVarValid], PTpts, method=EOSinterpMethod)

                    # Check that NaN removal worked correctly
                    errNaNstart = 'Failed to interpolate over NaNs in PerplexEOS '
                    errNaNend = ' values. The NaN gap may be too large to use this Perple_X output.'
                    if np.any(np.isnan(rho_kgm3)): raise RuntimeError(errNaNstart + 'rho' + errNaNend)
                    if np.any(np.isnan(VP_kms)): raise RuntimeError(errNaNstart + 'VP' + errNaNend)
                    if np.any(np.isnan(VS_kms)): raise RuntimeError(errNaNstart + 'VS' + errNaNend)
                    if np.any(np.isnan(Cp_Jm3K)): raise RuntimeError(errNaNstart + 'Cp' + errNaNend)
                    if np.any(np.isnan(alpha_pK)): raise RuntimeError(errNaNstart + 'alpha' + errNaNend)
                    if np.any(np.isnan(KS_bar)): raise RuntimeError(errNaNstart + 'KS' + errNaNend)
                    if np.any(np.isnan(GS_bar)): raise RuntimeError(errNaNstart + 'GS' + errNaNend)

                    # Now make 2D grids of values.
                    rho_kgm3 = np.reshape(rho_kgm3, (-1,dim2))
                    VP_kms = np.reshape(VP_kms, (-1,dim2))
                    VS_kms = np.reshape(VS_kms, (-1,dim2))
                    Cp_JkgK = np.reshape(Cp_Jm3K, (-1,dim2)) / rho_kgm3
                    alpha_pK = np.reshape(alpha_pK, (-1,dim2))
                    KS_GPa = np.reshape(KS_bar, (-1,dim2)) * Constants.bar2GPa
                    GS_GPa = np.reshape(GS_bar, (-1,dim2)) * Constants.bar2GPa

                    if P_FIRST:
                        # Transpose 2D meshes if P is the first column.
                        rho_kgm3 = rho_kgm3.T
                        VP_kms = VP_kms.T
                        VS_kms = VS_kms.T
                        Cp_JkgK = Cp_JkgK.T
                        alpha_pK = alpha_pK.T
                        KS_GPa = KS_GPa.T
                        GS_GPa = GS_GPa.T

                P_MPa = P1D_MPa
                T_K = T1D_K

                # Assign temporary functions we will wrap with porosity if modeled
                self.ufn_rho_kgm3 = RectBivariateSpline(P1D_MPa, T1D_K, rho_kgm3)
                self.ufn_VP_kms = RectBivariateSpline(P1D_MPa, T1D_K, VP_kms)
                self.ufn_VS_kms = RectBivariateSpline(P1D_MPa, T1D_K, VS_kms)
                self.ufn_KS_GPa = RectBivariateSpline(P1D_MPa, T1D_K, KS_GPa)
                self.ufn_GS_GPa = RectBivariateSpline(P1D_MPa, T1D_K, GS_GPa)
                self.ufn_Cp_JkgK = RectBivariateSpline(P1D_MPa, T1D_K, Cp_JkgK)
                self.ufn_alpha_pK = RectBivariateSpline(P1D_MPa, T1D_K, alpha_pK)

                EOSlist.loaded[tableKey] = (P_MPa, T_K, self.ufn_rho_kgm3, self.ufn_VP_kms, self.ufn_VS_kms,
                                              self.ufn_KS_GPa, self.ufn_GS_GPa, self.ufn_Cp_JkgK, self.ufn_alpha_pK)
                EOSlist.ranges[tableKey] = (self.Pmin, self.Pmax, self.Tmin, self.Tmax, self.deltaP, self.deltaT)

            self.rangeLabel = f'{self.Pmin},{self.Pmax},{self.deltaP},' + \
                              f'{self.Tmin},{self.Tmax},{self.deltaT}'

            # Assign thermal conductivity
            # (currently a placeholder)
            if kThermConst_WmK is None:
                if Fe_EOS:
                    kThermConst_WmK = Constants.kThermFe_WmK
                else:
                    kThermConst_WmK = Constants.kThermSil_WmK
            kTherm_WmK = np.zeros((np.size(P_MPa), np.size(T_K))) + kThermConst_WmK  # Placeholder until a self-consistent determination is implemented
            self.ufn_kTherm_WmK = RectBivariateSpline(P_MPa, T_K, kTherm_WmK)

            # Assign tidal heating function
            # (currently a placeholder)
            Htidal_Wm3 = np.zeros_like(kTherm_WmK) + HtidalConst_Wm3  # Placeholder until a self-consistent determination is implemented
            self.fn_Htidal_Wm3 = GetHtidalFunc(HtidalConst_Wm3)

            # Assign porosity model function, if applicable
            if Fe_EOS:
                # No porosity modeled
                self.ufn_phi_frac = ReturnZeros(1)

            elif porosType is None or porosType == 'none':
                # No porosity modeled, but need a dummy field for cross-compatibility
                self.ufn_phi_frac = ReturnZeros(1)

            else:
                if porosType == 'Vitovtova2014' or porosType == 'Chen2020':
                    # Pressure-depth lookup using the Preliminary Earth Reference Model,
                    # Dziewonski and Anderson (1981): https://doi.org/10.1016/0031-9201(81)90046-7
                    PREMzPfile = os.path.join(_ROOT, 'Thermodynamics', 'EOSdata', 'PREMtable.txt')
                    zPREM_km, PPREM_kbar = np.loadtxt(PREMzPfile, skiprows=2, unpack=True, delimiter=',')
                    PPREM_MPa = PPREM_kbar * 1e3 * Constants.bar2MPa
                    self.PREMlookup = Interp1D(PPREM_MPa, zPREM_km)
                else:
                    self.PREMlookup = None
                self.ufn_phi_frac = GetphiFunc(self.porosType, self.phiTop_frac, self.Pclosure_MPa,
                                              self.PREMlookup, P_MPa, T_K)

            # Store complete EOSStruct in global list of loaded EOSs
            EOSlist.loaded[self.EOSlabel] = self
            EOSlist.ranges[self.EOSlabel] = self.rangeLabel

    def fn_porosCorrect(self, propBulk, propPore, phi, J):
        # Combine pore fluid properties with matrix properties in accordance with
        # Yu et al. (2016): http://dx.doi.org/10.1016/j.jrmge.2015.07.004
        return (propBulk**J * (1 - phi) + propPore**J * phi) ** (1/J)

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


class EOSwrapper:
    """ Lightweight wrapper for accessing EOS functions stored in the EOSlist dict. """

    def __init__(self, key):
        self.key = key
        # Assign only those attributes we reference in functions
        if EOSlist.loaded[self.key].EOStype == 'ice':
            self.phaseID = EOSlist.loaded[self.key].phaseID
            self.POROUS = EOSlist.loaded[self.key].POROUS
        elif EOSlist.loaded[self.key].EOStype == 'ocean':
            self.Pmax = EOSlist.loaded[self.key].Pmax
        elif EOSlist.loaded[self.key].EOStype == 'inner':
            self.comp = EOSlist.loaded[self.key].comp

    def fn_phase(self, P_MPa, T_K):
        return EOSlist.loaded[self.key].fn_phase(P_MPa, T_K)
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
    def fn_sigma_Sm(self, P_MPa, T_K):
        return EOSlist.loaded[self.key].fn_sigma_Sm(P_MPa, T_K)
    def fn_Seismic(self, P_MPa, T_K):
        return EOSlist.loaded[self.key].fn_Seismic(P_MPa, T_K)


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
        nPs = np.size(P)
        nTs = np.size(T)
        if self.nVar > 1:
            out = (np.zeros(np.maximum(nPs, nTs)) for _ in range(self.nVar))
        else:
            out = np.zeros(np.maximum(nPs, nTs))
        return out


def TsolidusHirschmann2000(P_MPa):
    """ Silicate melting temperature parameterization based on
        Hirschmann (2000): https://doi.org/10.1029/2000GC000070 .

        Args:
            P_MPa (float, shape N): Pressure values in MPa to evaluate
        Returns:
            Tsolidus_K (float, shape N): Solidus temperature for each P value
    """
    P_GPa = P_MPa * 1e-3
    a = -5.104
    b = 132.899
    c = 1120.661
    Tsolidus_K = a*P_GPa**2 + b*P_GPa + c + Constants.T0
    return Tsolidus_K


class GetHtidalFunc:
    def __init__(self, HtidalConst_Wm3):
        # Tidal heating as a function of density, gravity, and bulk and shear moduli
        # Currently a placeholder returning a constant value, awaiting a self-
        # consistent calculation.
        self.HtidalConst_Wm3 = HtidalConst_Wm3

    def __call__(self, rho_kgm3, g_ms2, KS_GPa, GS_GPa):
        return self.HtidalConst_Wm3


def GetphiFunc(porosType, phiTop_frac, Pclosure_MPa, PREMlookup, P1D_MPa, T1D_K):
    # Porosity models as functions of P and T, set by variable vacuum maximum
    # porosity, pore closure pressure, or model coefficients.

    # Get 2D grids of P, T values for constructing T-independent functions
    # of porosity
    P2D_MPa, T2D_MPa = np.meshgrid(P1D_MPa, T1D_K, indexing='ij')

    if porosType == 'Han2014':
        # Han et al. (2014) porosity model: https://doi.org/10.1002/2014GL059378
        # This is an exponential model based on a surface value and a pore closure
        # pressure. A more robust model might have a closure pressure that is
        # self-consistently determined as a function of temperature, instead of
        # using a value prescribed by the user.
        if phiTop_frac == 0 or Pclosure_MPa == 0:
            phi_frac = np.zeros_like(P2D_MPa)
        else:
            c = 6.15
            phi_frac = phiTop_frac * np.exp(-c * P2D_MPa / Pclosure_MPa)
    elif porosType == 'Vitovtova2014' or porosType == 'Chen2020':
        # Earth rock porosity parameterizations referenced to ocean worlds
        # via PREM lookup
        # Create an array of repeated (T-independent) arrays of PREM-equivalent depths
        zPREMequiv_km = np.tile(PREMlookup(P1D_MPa), (np.size(T1D_K), 1)).T
        if porosType == 'Vitovtova2014':
            # Vitovtova et al. (2014) porosity model: https://doi.org/10.1134/S1069351314040181
            # This is "generalized depth trend" model from Vitovtova et al., 2014
            # for porosity of Earth crust. Note that the equation reported in the text
            # does not match that in the figure, and the equation from the figure
            # matches the plot, so we assume the equation from the figure is
            # the correct one.
            phi_frac = 10 ** (-0.65 - 0.16 * zPREMequiv_km + 0.0019 * zPREMequiv_km ** 2)
            # Because the second-order term has a positive coefficient, porosity
            # will increase above about ~42 km PREM depth, but the porosity should
            # continue to decrease. For physical consistency, we set porosity to
            # zero beyond the minimum point.
            zPhiMin_km = 0.16 / (2 * 0.0019)
            phi_frac[zPREMequiv_km > zPhiMin_km] = 0.0
        else:
            # Chen et al. (2020) porosity model: https://doi.org/10.1007/s10040-020-02214-x
            # 1/(1 + z)^n model based on a surface value and an empirical parameterization
            # relevant to Earth. Here we use the "oceanic crust" fitting parameters.
            m = 0.008
            n = 89.53
            phi0 = 0.678
            phi_frac = phi0 / (1 + m * zPREMequiv_km) ** n
    else:
        # Invalid porosity type
        raise ValueError(f'Porosity type "{porosType}" is not supported.')

    # Create unchanging function for the expected porosity
    fn_phi_frac = RectBivariateSpline(P1D_MPa, T1D_K, phi_frac)
    return fn_phi_frac


class GetphiCalc:
    def __init__(self, phiMax_frac, fn_phiEOS_frac, phiMin_frac):
        self.phiMin_frac = phiMin_frac
        self.phiMax_frac = phiMax_frac
        self.multFactor = 1.0
        if type(fn_phiEOS_frac(0,0)) == ReturnZeros:
            # For some reason, binding a ReturnZeros func as a method requires an extra dummy call
            self.fn_phiEOS_frac = fn_phiEOS_frac(0, 0)
        else:
            self.fn_phiEOS_frac = fn_phiEOS_frac

    def update(self, newPhiMax_frac):
        self.multFactor = newPhiMax_frac / self.phiMax_frac

    def __call__(self, P_MPa, T_K, grid=False):
        if type(self.fn_phiEOS_frac) == ReturnZeros:
            phi_frac = self.fn_phiEOS_frac(P_MPa, T_K, grid=grid)
        else:
            phi_frac = self.multFactor * self.fn_phiEOS_frac(P_MPa, T_K, grid=grid)
            if np.size(P_MPa) == 1:
                if phi_frac < self.phiMin_frac:
                    phi_frac = 0
            else:
                phi_frac[phi_frac < self.phiMin_frac] = 0

        return phi_frac
