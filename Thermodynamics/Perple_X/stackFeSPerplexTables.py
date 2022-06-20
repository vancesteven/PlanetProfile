import os
import numpy as np
from glob import glob
from scipy.io import savemat
from scipy.interpolate import RectBivariateSpline, interp1d, griddata as GridData

bar2MPa = 1.01325e-1
bar2GPa = 1.01325e-4

class PerplexEOS:
    def __init__(self, fpath, nHeaders, EOSinterpMethod='nearest', nP=None, nT=None):
        # Get Fe wt% and S wt% from fpath
        fname = os.path.basename(fpath)
        # fname format: 'Fe###_S@@@_1.tab', where ### is iron wt% and @@@ is sulfur wt%
        self.wFe_ppt = float(fpath[2:5]) * 10
        self.wS_ppt = float(fpath[7:10]) * 10

        # Load in Perple_X data. Note that all P, KS, and GS are stored as bar
        firstPT, secondPT, rho_kgm3, VP_kms, VS_kms, Cp_Jm3K, alpha_pK, KS_bar, GS_bar \
            = np.loadtxt(fpath, skiprows=nHeaders, unpack=True)
        # We don't know yet whether P or T is in the first column, which increments the fastest.
        # The second column increments once for each time the first column runs through the whole range.
        dim2 = next(i[0] for i, val in np.ndenumerate(secondPT) if val != secondPT[0])
        # Check the column header to see if P or T is printed first
        with open(fpath) as f:
            [f.readline() for _ in range(nHeaders - 1)]
            colHeaderLine = f.readline().strip()
        # Get necessary values to generate P, T arrays
        if colHeaderLine[0] == 'P':
            P_FIRST = True
            Plin_MPa = firstPT * bar2MPa  # Pressures saved as bar
            Tlin_K = secondPT
            lenP = dim2
        elif colHeaderLine[0] == 'T':
            P_FIRST = False
            Plin_MPa = secondPT * bar2MPa
            Tlin_K = firstPT
            lenP = int(len(Plin_MPa) / dim2)
        else:
            raise ValueError(f'Perple_X table {fpath} does not have T or P in the first column.')

        # Make 1D arrays of P and T that just span the axes (unlike the 2D meshes of the dependent variables)
        self.Pmin = Plin_MPa[0]
        self.Pmax = Plin_MPa[-1]
        self.Tmin = Tlin_K[0]
        self.Tmax = Tlin_K[-1]
        lenT = int(len(Tlin_K) / lenP)
        P1D_MPa, self.deltaP = np.linspace(self.Pmin, self.Pmax, lenP, retstep=True)
        T1D_K, self.deltaT = np.linspace(self.Tmin, self.Tmax, lenT, retstep=True)

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
        self.rho_kgm3 = np.reshape(rho_kgm3, (-1, dim2))
        self.VP_kms = np.reshape(VP_kms, (-1, dim2))
        self.VS_kms = np.reshape(VS_kms, (-1, dim2))
        self.Cp_JkgK = np.reshape(Cp_Jm3K, (-1, dim2)) / self.rho_kgm3
        self.alpha_pK = np.reshape(alpha_pK, (-1, dim2))
        self.KS_GPa = np.reshape(KS_bar, (-1, dim2)) * bar2GPa
        self.GS_GPa = np.reshape(GS_bar, (-1, dim2)) * bar2GPa

        if P_FIRST:
            # Transpose 2D meshes if P is the first column.
            self.rho_kgm3 = self.rho_kgm3.T
            self.VP_kms = self.VP_kms.T
            self.VS_kms = self.VS_kms.T
            self.Cp_JkgK = self.Cp_JkgK.T
            self.alpha_pK = self.alpha_pK.T
            self.KS_GPa = self.KS_GPa.T
            self.GS_GPa = self.GS_GPa.T

        self.P_MPa = P1D_MPa
        self.T_K = T1D_K

        # Shrink output tables to desired amount if number of P or T values is specified
        if nP is not None or nT is not None:
            fn_rho_kgm3 = RectBivariateSpline(P1D_MPa, T1D_K, self.rho_kgm3)
            fn_VP_kms = RectBivariateSpline(P1D_MPa, T1D_K, self.VP_kms)
            fn_VS_kms = RectBivariateSpline(P1D_MPa, T1D_K, self.VS_kms)
            fn_KS_GPa = RectBivariateSpline(P1D_MPa, T1D_K, self.KS_GPa)
            fn_GS_GPa = RectBivariateSpline(P1D_MPa, T1D_K, self.GS_GPa)
            fn_Cp_JkgK = RectBivariateSpline(P1D_MPa, T1D_K, self.Cp_JkgK)
            fn_alpha_pK = RectBivariateSpline(P1D_MPa, T1D_K, self.alpha_pK)

            if nP is not None:
                self.P_MPa, self.deltaP = np.linspace(self.Pmin, self.Pmax, nP, retstep=True)
            if nT is not None:
                self.T_K, self.deltaT = np.linspace(self.Tmin, self.Tmax, nT, retstep=True)

            self.rho_kgm3 = fn_rho_kgm3(self.P_MPa, self.T_K, grid=True)
            self.VP_kms = fn_VP_kms(self.P_MPa, self.T_K, grid=True)
            self.VS_kms = fn_VS_kms(self.P_MPa, self.T_K, grid=True)
            self.KS_GPa = fn_KS_GPa(self.P_MPa, self.T_K, grid=True)
            self.GS_GPa = fn_GS_GPa(self.P_MPa, self.T_K, grid=True)
            self.Cp_JkgK = fn_Cp_JkgK(self.P_MPa, self.T_K, grid=True)
            self.alpha_pK = fn_alpha_pK(self.P_MPa, self.T_K, grid=True)


# Assume we're in the same directory as the Perplex tables and get file list
tableFiles = np.sort(glob('Fe*.tab'))
nTables = np.size(tableFiles)
tables = np.empty(nTables, dtype=object)
for i, tableFile in enumerate(tableFiles):
    print(f'Loading Perplex table: {tableFile}, {i+1} of {nTables}.')
    tables[i] = PerplexEOS(tableFile, 13, nP=325, nT=325)

# Save .mat file with 3D grids for each variable
outFname = 'Fe-S_3D_EOS.mat'
saveDict = {
    'Perplex version': '6.6.6',
    'wFe_ppt': np.array([table.wFe_ppt for table in tables]),
    'wS_ppt': np.array([table.wS_ppt for table in tables]),
    'P_MPa': tables[0].P_MPa,
    'T_K': tables[0].T_K,
    'rho_kgm3': np.array([table.rho_kgm3 for table in tables]),
    'VP_kms': np.array([table.VP_kms for table in tables]),
    'VS_kms': np.array([table.VS_kms for table in tables]),
    'Cp_JkgK': np.array([table.Cp_JkgK for table in tables]),
    'alpha_pK': np.array([table.alpha_pK for table in tables]),
    'KS_GPa': np.array([table.KS_GPa for table in tables]),
    'GS_GPa': np.array([table.GS_GPa for table in tables])
}

savemat(outFname, saveDict)
