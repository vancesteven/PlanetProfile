import os
import numpy as np
from scipy.interpolate import RectBivariateSpline, interp1d as Interp1D, griddata as GridData
from Utilities.defineStructs import Constants

class PerplexEOSStruct:
    """ Loads in Perple_X table and creates interpolators that can be called to
        obtain silicate/core properties as functions of P and T.
    """
    def __init__(self, EOSfname, EOSinterpMethod='nearest', nHeaders=13, Fe_EOS=False, kThermConst_WmK=None,
                 HtidalConst_Wm3=0, porosType=None, phiTop_frac=0, Pclosure_MPa=350):
        self.comp = EOSfname[:-4]
        self.dir = os.path.join('Thermodynamics', 'Perple_X', 'output_data')
        self.fpath = os.path.join(self.dir, EOSfname)

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
        self.Pmin_MPa = Plin_MPa[0]
        self.Pmax_MPa = Plin_MPa[-1]
        self.Tmin_K = Tlin_K[0]
        self.Tmax_K = Tlin_K[-1]
        lenT = int(len(Tlin_K) / lenP)
        P1D_MPa = np.linspace(self.Pmin_MPa, self.Pmax_MPa, lenP)
        T1D_K = np.linspace(self.Tmin_K, self.Tmax_K, lenT)

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
        self.rho_kgm3 = np.reshape(rho_kgm3, (-1,dim2))
        self.VP_kms = np.reshape(VP_kms, (-1,dim2))
        self.VS_kms = np.reshape(VS_kms, (-1,dim2))
        self.Cp_JkgK = np.reshape(Cp_Jm3K, (-1,dim2)) / self.rho_kgm3
        self.alpha_pK = np.reshape(alpha_pK, (-1,dim2))
        self.KS_GPa = np.reshape(KS_bar, (-1,dim2)) * Constants.bar2GPa
        self.GS_GPa = np.reshape(GS_bar, (-1,dim2)) * Constants.bar2GPa

        if P_FIRST:
            # Transpose 2D meshes if P is the first column.
            self.rho_kgm3 = self.rho_kgm3.T
            self.VP_kms = self.VP_kms.T
            self.VS_kms = self.VS_kms.T
            self.Cp_JkgK = self.Cp_JkgK.T
            self.alpha_pK = self.alpha_pK.T
            self.KS_GPa = self.KS_GPa.T
            self.GS_GPa = self.GS_GPa.T

        if kThermConst_WmK is None:
            if Fe_EOS:
                kThermConst_WmK = Constants.kThermFe_WmK
            else:
                kThermConst_WmK = Constants.kThermSil_WmK
        self.kTherm_WmK = np.zeros_like(self.alpha_pK) + kThermConst_WmK  # Placeholder until a self-consistent determination is implemented
        self.Htidal_Wm3 = np.zeros_like(self.alpha_pK) + HtidalConst_Wm3  # Placeholder until a self-consistent determination is implemented

        self.P_MPa = P1D_MPa
        self.T_K = T1D_K
        # Assign temporary functions we will wrap with porosity if modeled
        self.fn_rho_kgm3 = RectBivariateSpline(P1D_MPa, T1D_K, self.rho_kgm3)
        self.fn_VP_kms = RectBivariateSpline(P1D_MPa, T1D_K, self.VP_kms)
        self.fn_VS_kms = RectBivariateSpline(P1D_MPa, T1D_K, self.VS_kms)
        self.fn_KS_GPa = RectBivariateSpline(P1D_MPa, T1D_K, self.KS_GPa)
        self.fn_GS_GPa = RectBivariateSpline(P1D_MPa, T1D_K, self.GS_GPa)
        self.fn_kTherm_WmK = RectBivariateSpline(P1D_MPa, T1D_K, self.kTherm_WmK)
        self.fn_Cp_JkgK = RectBivariateSpline(P1D_MPa, T1D_K, self.Cp_JkgK)
        self.fn_alpha_pK = RectBivariateSpline(P1D_MPa, T1D_K, self.alpha_pK)

        # Assign tidal heating function
        # (currently a placeholder)
        self.fn_Htidal_Wm3 = GetHtidalFunc(HtidalConst_Wm3)

        # Assign porosity model function, if applicable
        if Fe_EOS:
            # No porosity modeled, and no need for dummy field
            pass

        elif porosType is None or porosType == 'none':
            # No porosity modeled, but need a dummy field for cross-compatibility
            self.fn_phi_frac = lambda P, T: np.zeros_like(P)

        else:
            if porosType == 'Vitovtova2014' or porosType == 'Chen2020':
                # Pressure-depth lookup using the Preliminary Earth Reference Model,
                # Dziewonski and Anderson (1981): https://doi.org/10.1016/0031-9201(81)90046-7
                PREMzPfile = os.path.join('Thermodynamics', 'FromLiterature', 'PREMtable.txt')
                zPREM_km, PPREM_kbar = np.loadtxt(PREMzPfile, skiprows=2, unpack=True, delimiter=',')
                PPREM_MPa = PPREM_kbar * 1e3 * Constants.bar2MPa
                self.PREMlookup = Interp1D(PPREM_MPa, zPREM_km)
            else:
                self.PREMlookup = None
            self.fn_phi_frac = GetphiFunc(porosType, phiTop_frac, Pclosure_MPa, self.PREMlookup, P1D_MPa, T1D_K)

            # Combine pore fluid properties with matrix properties in accordance with
            # Yu et al. (2016): http://dx.doi.org/10.1016/j.jrmge.2015.07.004
            self.fn_porosCorrect = lambda propBulk, propPore, phi, J: (propBulk**J * (1 - phi)
                                                 + propPore**J * phi)**(1/J)


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


def GetHtidalFunc(HtidalConst_Wm3):
    # Tidal heating as a function of density, gravity, and bulk and shear moduli
    # Currently a placeholder returning a constant value, awaiting a self-
    # consistent calculation.
    fn_Htidal_Wm3 = lambda rho, g, KS, GS: HtidalConst_Wm3
    return fn_Htidal_Wm3


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
            # Set porosity to be zero when it is vanishingly small, i.e.
            # at pressures greater than double the pore closure pressure.
            phi_frac[P2D_MPa > 2 * Pclosure_MPa] = 0.0
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
