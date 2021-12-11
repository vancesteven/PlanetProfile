import numpy as np
import scipy.interpolate as spi
from Utilities.dataStructs import Constants

class PerplexEOSStruct:
    def __init__(self, EOSfname, nHeaders=13):
        self.comp = EOSfname[:-4]
        self.dir = 'Thermodynamics/Perple_X/output_data/'
        self.fpath = self.dir + EOSfname

        # Load in Perple_X data. Note that all P, KS, and GS are stored as bar
        firstPT, secondPT, rho_kgm3, VP_kms, VS_kms, Cp_JkgK, alpha_pK, KS_bar, GS_bar \
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
            raise ValueError('Perple_X table ' + EOSfname + ' does not have T or P in the first column.')

        # Make 1D arrays of P and T that just span the axes (unlike the 2D meshes of the dependent variables)
        P1D_MPa = np.linspace(Plin_MPa[0], Plin_MPa[-1], lenP)
        T1D_K = np.linspace(Tlin_K[0], Tlin_K[-1], int(len(Tlin_K)/lenP))

        # Interpolate dependent variables where the values are NaN
        PTpts = (Plin_MPa, Tlin_K)
        thisNaNyn = np.isnan(rho_kgm3)
        if np.any(thisNaNyn):
            thisNaNPs = Plin_MPa[np.where(thisNaNyn)]
            thisNaNTs = Tlin_K[np.where(thisNaNyn)]
            rho_kgm3[np.where(thisNaNyn)] = spi.griddata(PTpts, rho_kgm3, (thisNaNPs, thisNaNTs), method='cubic')
        thisNaNyn = np.isnan(VP_kms)
        if np.any(thisNaNyn):
            thisNaNPs = Plin_MPa[np.where(thisNaNyn)]
            thisNaNTs = Tlin_K[np.where(thisNaNyn)]
            VP_kms[np.where(thisNaNyn)] = spi.griddata(PTpts, VP_kms, (thisNaNPs, thisNaNTs), method='cubic')
        thisNaNyn = np.isnan(VS_kms)
        if np.any(thisNaNyn):
            thisNaNPs = Plin_MPa[np.where(thisNaNyn)]
            thisNaNTs = Tlin_K[np.where(thisNaNyn)]
            VS_kms[np.where(thisNaNyn)] = spi.griddata(PTpts, VS_kms, (thisNaNPs, thisNaNTs), method='cubic')
        thisNaNyn = np.isnan(Cp_JkgK)
        if np.any(thisNaNyn):
            thisNaNPs = Plin_MPa[np.where(thisNaNyn)]
            thisNaNTs = Tlin_K[np.where(thisNaNyn)]
            Cp_JkgK[np.where(thisNaNyn)] = spi.griddata(PTpts, Cp_JkgK, (thisNaNPs, thisNaNTs), method='cubic')
        thisNaNyn = np.isnan(alpha_pK)
        if np.any(thisNaNyn):
            thisNaNPs = Plin_MPa[np.where(thisNaNyn)]
            thisNaNTs = Tlin_K[np.where(thisNaNyn)]
            alpha_pK[np.where(thisNaNyn)] = spi.griddata(PTpts, alpha_pK, (thisNaNPs, thisNaNTs), method='cubic')
        thisNaNyn = np.isnan(KS_bar)
        if np.any(thisNaNyn):
            thisNaNPs = Plin_MPa[np.where(thisNaNyn)]
            thisNaNTs = Tlin_K[np.where(thisNaNyn)]
            KS_bar[np.where(thisNaNyn)] = spi.griddata(PTpts, KS_bar, (thisNaNPs, thisNaNTs), method='cubic')
        thisNaNyn = np.isnan(GS_bar)
        if np.any(thisNaNyn):
            thisNaNPs = Plin_MPa[np.where(thisNaNyn)]
            thisNaNTs = Tlin_K[np.where(thisNaNyn)]
            GS_bar[np.where(thisNaNyn)] = spi.griddata(PTpts, GS_bar, (thisNaNPs, thisNaNTs), method='cubic')

        # Check that NaN removal worked correctly
        errNaNstart = 'Failed to interpolate over NaNs in PerplexEOS '
        errNaNend = ' values. The NaN gap may be too large to use this Perple_X output.'
        if np.any(np.isnan(rho_kgm3)): raise RuntimeError(errNaNstart + 'rho' +errNaNend)
        if np.any(np.isnan(VP_kms)): raise RuntimeError(errNaNstart + 'VP' +errNaNend)
        if np.any(np.isnan(VS_kms)): raise RuntimeError(errNaNstart + 'VS' +errNaNend)
        if np.any(np.isnan(Cp_JkgK)): raise RuntimeError(errNaNstart + 'Cp' +errNaNend)
        if np.any(np.isnan(alpha_pK)): raise RuntimeError(errNaNstart + 'alpha' +errNaNend)
        if np.any(np.isnan(KS_bar)): raise RuntimeError(errNaNstart + 'KS' +errNaNend)
        if np.any(np.isnan(GS_bar)): raise RuntimeError(errNaNstart + 'GS' +errNaNend)

        # Now make 2D grids of values.
        rho_kgm3 = np.reshape(rho_kgm3, (-1,dim2))
        VP_kms = np.reshape(VP_kms, (-1,dim2))
        VS_kms = np.reshape(VS_kms, (-1,dim2))
        Cp_JkgK = np.reshape(Cp_JkgK, (-1,dim2))
        alpha_pK = np.reshape(alpha_pK, (-1,dim2))
        KS_GPa = np.reshape(KS_bar, (-1,dim2)) * Constants.bar2GPa
        GS_GPa = np.reshape(GS_bar, (-1,dim2)) * Constants.bar2GPa

        if P_FIRST:
            # Transpose 2D meshes if P is the first column
            rho_kgm3 = rho_kgm3.T
            VP_kms = VP_kms.T
            VS_kms = VS_kms.T
            Cp_JkgK = Cp_JkgK.T
            alpha_pK = alpha_pK.T
            KS_GPa = KS_GPa.T
            GS_GPa = GS_GPa.T

        self.fn_rho_kgm3 = spi.RectBivariateSpline(P1D_MPa, T1D_K, rho_kgm3)
        self.fn_VP_kms = spi.RectBivariateSpline(P1D_MPa, T1D_K, VP_kms)
        self.fn_VS_kms = spi.RectBivariateSpline(P1D_MPa, T1D_K, VS_kms)
        self.fn_Cp_JkgK = spi.RectBivariateSpline(P1D_MPa, T1D_K, Cp_JkgK)
        self.fn_alpha_pK = spi.RectBivariateSpline(P1D_MPa, T1D_K, alpha_pK)
        self.fn_KS_GPa = spi.RectBivariateSpline(P1D_MPa, T1D_K, KS_GPa)
        self.fn_GS_GPa = spi.RectBivariateSpline(P1D_MPa, T1D_K, GS_GPa)

def MantleEOS(PerplexEOS, Pstart, Tstart, rStart, rEnd):
    fn_rhoSil_kgm3 = PerplexEOS.fn_rho_kgm3
    fn_VPSil_kms = PerplexEOS.fn_VP_kms

    return