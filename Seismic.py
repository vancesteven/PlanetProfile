import numpy as np
from Thermodynamics.FromLiterature.ThermalProfiles import TsolidusHirschmann2000

def SeismicCalcs(Planet, Params):
    """ Calculation of seismic properties, including wave speeds

        Assigns Planet attributes:
            Seismic.VP_kms, Seismic.VS_kms, Seismic.QS, Seismic.KS_GPa, Seismic.GS_GPa
    """
    # Initialize arrays
    Planet.Seismic.VP_kms, Planet.Seismic.VS_kms, Planet.Seismic.QS, Planet.Seismic.KS_GPa, \
        Planet.Seismic.GS_GPa = (np.zeros(Planet.Steps.nTotal) for _ in range(5))

    # Evaluate silicate EOS for seismic properties
    iOS = Planet.Steps.nHydro
    iSC = Planet.Steps.nHydro + Planet.Steps.nSil
    Planet.Cp_JkgK[iOS:iSC] = [Planet.Sil.EOS.fn_Cp_JkgK(Planet.P_MPa[i], Planet.T_K[i]) for i in range(iOS, iSC)]
    Planet.alpha_pK[iOS:iSC] = [Planet.Sil.EOS.fn_alpha_pK(Planet.P_MPa[i], Planet.T_K[i]) for i in range(iOS, iSC)]
    Planet.Seismic.VP_kms[iOS:iSC] = [Planet.Sil.EOS.fn_VP_kms(Planet.P_MPa[i], Planet.T_K[i]) for i in range(iOS, iSC)]
    Planet.Seismic.VS_kms[iOS:iSC] = [Planet.Sil.EOS.fn_VS_kms(Planet.P_MPa[i], Planet.T_K[i]) for i in range(iOS, iSC)]
    Planet.Seismic.KS_GPa[iOS:iSC] = [Planet.Sil.EOS.fn_KS_GPa(Planet.P_MPa[i], Planet.T_K[i]) for i in range(iOS, iSC)]
    Planet.Seismic.GS_GPa[iOS:iSC] = [Planet.Sil.EOS.fn_GS_GPa(Planet.P_MPa[i], Planet.T_K[i]) for i in range(iOS, iSC)]
    TsilMelt = TsolidusHirschmann2000(Planet.P_MPa[iOS:iSC])
    Planet.Seismic.QS[iOS:iSC] = Planet.Seismic.BSil * np.exp(
        Planet.Seismic.gammaSil * Planet.Seismic.gSil * TsilMelt / Planet.T_K[iOS:iSC])

    # Evaluate core EOS for seismic properties
    iCC = Planet.Steps.nTotal
    Planet.Seismic.VP_kms[iSC:iCC] = [Planet.Core.EOS.fn_VP_kms(Planet.P_MPa[i], Planet.T_K[i]) for i in range(iSC, iCC)]
    Planet.Seismic.VS_kms[iSC:iCC] = [Planet.Core.EOS.fn_VS_kms(Planet.P_MPa[i], Planet.T_K[i]) for i in range(iSC, iCC)]
    Planet.Seismic.KS_GPa[iSC:iCC] = [Planet.Core.EOS.fn_KS_GPa(Planet.P_MPa[i], Planet.T_K[i]) for i in range(iSC, iCC)]
    Planet.Seismic.GS_GPa[iSC:iCC] = [Planet.Core.EOS.fn_GS_GPa(Planet.P_MPa[i], Planet.T_K[i]) for i in range(iSC, iCC)]
    if Planet.Seismic.QScore is not None: Planet.Seismic.QS[iSC:iCC] = Planet.Seismic.QScore

    return Planet