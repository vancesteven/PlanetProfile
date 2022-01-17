import numpy as np
from Thermodynamics.FromLiterature.InnerEOS import TsolidusHirschmann2000
from Thermodynamics.FromLiterature.HydroEOS import GetPhaseIndices
from Utilities.dataStructs import Constants

def SeismicCalcs(Planet, Params):
    """ Calculation of seismic properties, including wave speeds

        Assigns Planet attributes:
            Seismic.VP_kms, Seismic.VS_kms, Seismic.QS, Seismic.KS_GPa, Seismic.GS_GPa
    """
    # Initialize arrays
    Planet.Seismic.VP_kms, Planet.Seismic.VS_kms, Planet.Seismic.QS, Planet.Seismic.KS_GPa, \
        Planet.Seismic.GS_GPa = (np.zeros(Planet.Steps.nTotal) for _ in range(5))

    if Params.CALC_SEISMIC:

        indsLiquid, indsIceI, indsIceII, indsIceIII, indsIceV, indsIceVI, indsClath, indsSil, indsFe = GetPhaseIndices(Planet.phase)

        if Planet.Do.CLATHRATE:
            # Get seismic properties of surface clathrate layer
            Planet.Seismic.VP_kms[indsClath], Planet.Seismic.VS_kms[indsClath], \
            Planet.Seismic.KS_GPa[indsClath], Planet.Seismic.GS_GPa[indsClath] \
                = Planet.Ocean.surfIceEOS['Clath'].fn_Seismic(Planet.P_MPa[indsClath], Planet.T_K[indsClath])

            Hclath = Planet.Seismic.gClath * np.max(Planet.T_K[indsClath])
            Planet.Seismic.QS[indsClath] = Planet.Seismic.BClath * np.exp(
                Planet.Seismic.gammaClath * Hclath / Planet.T_K[indsClath])

        # Get seismic properties of all ice layers
        if np.size(indsIceI) != 0:
            Planet.Seismic.VP_kms[indsIceI], Planet.Seismic.VS_kms[indsIceI], Planet.Seismic.KS_GPa[indsIceI], \
            Planet.Seismic.GS_GPa[indsIceI] = Planet.Ocean.surfIceEOS['Ih'].fn_Seismic(Planet.P_MPa[indsIceI], Planet.T_K[indsIceI])
            HiceI = Planet.Seismic.gIceI * Planet.Bulk.Tb_K
            Planet.Seismic.QS[indsIceI] = Planet.Seismic.BIceI * np.exp(
                Planet.Seismic.gammaIceI * HiceI / Planet.T_K[indsIceI])
        if np.size(indsIceII) != 0:
            surfIce = indsIceII < Planet.Steps.nSurfIce
            if np.any(surfIce):
                Planet.Seismic.VP_kms[indsIceII[surfIce]], Planet.Seismic.VS_kms[indsIceII[surfIce]], \
                Planet.Seismic.KS_GPa[indsIceII[surfIce]], Planet.Seismic.GS_GPa[indsIceII[surfIce]] \
                    = Planet.Ocean.surfIceEOS['II'].fn_Seismic(
                            Planet.P_MPa[indsIceII[surfIce]], Planet.T_K[indsIceII[surfIce]])
            if not np.all(surfIce):
                lowerIce = np.logical_not(surfIce)
                Planet.Seismic.VP_kms[indsIceII[lowerIce]], Planet.Seismic.VS_kms[indsIceII[lowerIce]], \
                Planet.Seismic.KS_GPa[indsIceII[lowerIce]], Planet.Seismic.GS_GPa[indsIceII[lowerIce]] \
                    = Planet.Ocean.iceEOS['II'].fn_Seismic(
                    Planet.P_MPa[indsIceII[lowerIce]], Planet.T_K[indsIceII[lowerIce]])
            HiceII = Planet.Seismic.gIceII * np.max(Planet.T_K[indsIceII])
            Planet.Seismic.QS[indsIceII] = Planet.Seismic.BIceII * np.exp(
                Planet.Seismic.gammaIceII * HiceII / Planet.T_K[indsIceII])
        if np.size(indsIceIII) != 0:
            surfIce = indsIceIII < Planet.Steps.nSurfIce
            if np.any(surfIce):
                Planet.Seismic.VP_kms[indsIceIII[surfIce]], Planet.Seismic.VS_kms[indsIceIII[surfIce]], \
                Planet.Seismic.KS_GPa[indsIceIII[surfIce]], Planet.Seismic.GS_GPa[indsIceIII[surfIce]] \
                    = Planet.Ocean.surfIceEOS['III'].fn_Seismic(
                            Planet.P_MPa[indsIceIII[surfIce]], Planet.T_K[indsIceIII[surfIce]])
            if not np.all(surfIce):
                lowerIce = np.logical_not(surfIce)
                Planet.Seismic.VP_kms[indsIceIII[lowerIce]], Planet.Seismic.VS_kms[indsIceIII[lowerIce]], \
                Planet.Seismic.KS_GPa[indsIceIII[lowerIce]], Planet.Seismic.GS_GPa[indsIceIII[lowerIce]] \
                    = Planet.Ocean.iceEOS['III'].fn_Seismic(
                    Planet.P_MPa[indsIceIII[lowerIce]], Planet.T_K[indsIceIII[lowerIce]])
            HiceIII = Planet.Seismic.gIceIII * np.max(Planet.T_K[indsIceIII])
            Planet.Seismic.QS[indsIceIII] = Planet.Seismic.BIceIII * np.exp(
                Planet.Seismic.gammaIceIII * HiceIII / Planet.T_K[indsIceIII])
        if np.size(indsIceV) != 0:
            surfIce = indsIceV < Planet.Steps.nSurfIce
            if np.any(surfIce):
                Planet.Seismic.VP_kms[indsIceV[surfIce]], Planet.Seismic.VS_kms[indsIceV[surfIce]], \
                Planet.Seismic.KS_GPa[indsIceV[surfIce]], Planet.Seismic.GS_GPa[indsIceV[surfIce]] \
                    = Planet.Ocean.surfIceEOS['V'].fn_Seismic(
                            Planet.P_MPa[indsIceV[surfIce]], Planet.T_K[indsIceV[surfIce]])
            if not np.all(surfIce):
                lowerIce = np.logical_not(surfIce)
                Planet.Seismic.VP_kms[indsIceV[lowerIce]], Planet.Seismic.VS_kms[indsIceV[lowerIce]], \
                Planet.Seismic.KS_GPa[indsIceV[lowerIce]], Planet.Seismic.GS_GPa[indsIceV[lowerIce]] \
                    = Planet.Ocean.iceEOS['V'].fn_Seismic(
                    Planet.P_MPa[indsIceV[lowerIce]], Planet.T_K[indsIceV[lowerIce]])
            HiceV = Planet.Seismic.gIceV * np.max(Planet.T_K[indsIceV])
            Planet.Seismic.QS[indsIceV] = Planet.Seismic.BIceV * np.exp(
                Planet.Seismic.gammaIceV * HiceV / Planet.T_K[indsIceV])
        if np.size(indsIceVI) != 0:
            surfIce = indsIceVI < Planet.Steps.nSurfIce
            if np.any(surfIce):
                Planet.Seismic.VP_kms[indsIceVI[surfIce]], Planet.Seismic.VS_kms[indsIceVI[surfIce]], \
                Planet.Seismic.KS_GPa[indsIceVI[surfIce]], Planet.Seismic.GS_GPa[indsIceVI[surfIce]] \
                    = Planet.Ocean.surfIceEOS['VI'].fn_Seismic(
                            Planet.P_MPa[indsIceVI[surfIce]], Planet.T_K[indsIceVI[surfIce]])
            if not np.all(surfIce):
                lowerIce = np.logical_not(surfIce)
                Planet.Seismic.VP_kms[indsIceVI[lowerIce]], Planet.Seismic.VS_kms[indsIceVI[lowerIce]], \
                Planet.Seismic.KS_GPa[indsIceVI[lowerIce]], Planet.Seismic.GS_GPa[indsIceVI[lowerIce]] \
                    = Planet.Ocean.iceEOS['VI'].fn_Seismic(
                    Planet.P_MPa[indsIceVI[lowerIce]], Planet.T_K[indsIceVI[lowerIce]])
            HiceVI = Planet.Seismic.gIceVI * np.max(Planet.T_K[indsIceVI])
            Planet.Seismic.QS[indsIceVI] = Planet.Seismic.BIceVI * np.exp(
                Planet.Seismic.gammaIceVI * HiceVI / Planet.T_K[indsIceVI])
            
        # Get seismic properties of ocean layers
        if np.size(indsLiquid) != 0:
            Planet.Seismic.VP_kms[indsLiquid], Planet.Seismic.KS_GPa[indsLiquid] \
                = Planet.Ocean.EOS.fn_Seismic(Planet.P_MPa[indsLiquid], Planet.T_K[indsLiquid])
            # VS_kms, QS, and GS_GPa were all initialized as zero and are all zero for these
            # layers, so we just leave them as-is.

        if not Params.SKIP_INNER:
            # Evaluate silicate EOS for seismic properties
            Planet.Cp_JkgK[indsSil] = Planet.Sil.EOS.fn_Cp_JkgK(Planet.P_MPa[indsSil], Planet.T_K[indsSil], grid=False)
            Planet.alpha_pK[indsSil] = Planet.Sil.EOS.fn_alpha_pK(Planet.P_MPa[indsSil], Planet.T_K[indsSil], grid=False)
            Planet.Seismic.VP_kms[indsSil] = Planet.Sil.EOS.fn_VP_kms(Planet.P_MPa[indsSil], Planet.T_K[indsSil], grid=False)
            Planet.Seismic.VS_kms[indsSil] = Planet.Sil.EOS.fn_VS_kms(Planet.P_MPa[indsSil], Planet.T_K[indsSil], grid=False)
            Planet.Seismic.KS_GPa[indsSil] = Planet.Sil.EOS.fn_KS_GPa(Planet.P_MPa[indsSil], Planet.T_K[indsSil], grid=False)
            Planet.Seismic.GS_GPa[indsSil] = Planet.Sil.EOS.fn_GS_GPa(Planet.P_MPa[indsSil], Planet.T_K[indsSil], grid=False)
            Hsil = Planet.Seismic.gSil * TsolidusHirschmann2000(Planet.P_MPa[indsSil])
            Planet.Seismic.QS[indsSil] = Planet.Seismic.BSil * np.exp(
                Planet.Seismic.gammaSil * Hsil / Planet.T_K[indsSil])

            if Planet.Do.Fe_CORE:
                # Evaluate core EOS for seismic properties
                Planet.kTherm_WmK[indsFe] = Planet.Core.EOS.fn_kTherm_WmK(Planet.P_MPa[indsFe], Planet.T_K[indsFe], grid=False)
                Planet.Seismic.VP_kms[indsFe] = Planet.Core.EOS.fn_VP_kms(Planet.P_MPa[indsFe], Planet.T_K[indsFe], grid=False)
                Planet.Seismic.VS_kms[indsFe] = Planet.Core.EOS.fn_VS_kms(Planet.P_MPa[indsFe], Planet.T_K[indsFe], grid=False)
                Planet.Seismic.KS_GPa[indsFe] = Planet.Core.EOS.fn_KS_GPa(Planet.P_MPa[indsFe], Planet.T_K[indsFe], grid=False)
                Planet.Seismic.GS_GPa[indsFe] = Planet.Core.EOS.fn_GS_GPa(Planet.P_MPa[indsFe], Planet.T_K[indsFe], grid=False)
                if Planet.Seismic.QScore is not None:
                    Planet.Seismic.QS[indsFe] = Planet.Seismic.QScore
                else:
                    Planet.Seismic.QS[indsFe] = Constants.QScore

    return Planet
