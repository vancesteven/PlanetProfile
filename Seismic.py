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

            # Correct for porous clathrates, assuming pore space may be treated as vacuum
            if Planet.Do.POROUS_ICE and Planet.Ocean.phiMax_frac['Clath'] is not None:
                Planet.Seismic.VP_kms[indsClath] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(
                    Planet.Seismic.VP_kms[indsClath], 0, Planet.phi_frac[indsClath], Planet.Ocean.JVP)
                Planet.Seismic.VS_kms[indsClath] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(
                    Planet.Seismic.VS_kms[indsClath], 0, Planet.phi_frac[indsClath], Planet.Ocean.JVS)
                Planet.Seismic.KS_GPa[indsClath] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(
                    Planet.Seismic.KS_GPa[indsClath], 0, Planet.phi_frac[indsClath], Planet.Ocean.JKS)
                Planet.Seismic.GS_GPa[indsClath] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(
                    Planet.Seismic.GS_GPa[indsClath], 0, Planet.phi_frac[indsClath], Planet.Ocean.JGS)

            Hclath = Planet.Seismic.gClath * np.max(Planet.T_K[indsClath])
            Planet.Seismic.QS[indsClath] = Planet.Seismic.BClath * np.exp(
                Planet.Seismic.gammaClath * Hclath / Planet.T_K[indsClath])

        # Get seismic properties of all ice layers
        if np.size(indsIceI) != 0:
            Planet.Seismic.VP_kms[indsIceI], Planet.Seismic.VS_kms[indsIceI], Planet.Seismic.KS_GPa[indsIceI], \
            Planet.Seismic.GS_GPa[indsIceI] = Planet.Ocean.surfIceEOS['Ih'].fn_Seismic(Planet.P_MPa[indsIceI], Planet.T_K[indsIceI])

            # Correct for porous ice, assuming pore space may be treated as vacuum
            if Planet.Do.POROUS_ICE:
                Planet.Seismic.VP_kms[indsIceI] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
                    Planet.Seismic.VP_kms[indsIceI], 0, Planet.phi_frac[indsIceI], Planet.Ocean.JVP)
                Planet.Seismic.VS_kms[indsIceI] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
                    Planet.Seismic.VS_kms[indsIceI], 0, Planet.phi_frac[indsIceI], Planet.Ocean.JVS)
                Planet.Seismic.KS_GPa[indsIceI] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
                    Planet.Seismic.KS_GPa[indsIceI], 0, Planet.phi_frac[indsIceI], Planet.Ocean.JKS)
                Planet.Seismic.GS_GPa[indsIceI] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
                    Planet.Seismic.GS_GPa[indsIceI], 0, Planet.phi_frac[indsIceI], Planet.Ocean.JGS)

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
            # Get Cp and alpha here, because we didn't calculate them earlier since we didn't need them
            # in calculating a conductive profile in the silicates and it would contribute extra,
            # unnecessary computational overhead there.
            Planet.Cp_JkgK[indsSil] = Planet.Sil.EOS.fn_Cp_JkgK(Planet.P_MPa[indsSil], Planet.T_K[indsSil], grid=False)
            Planet.alpha_pK[indsSil] = Planet.Sil.EOS.fn_alpha_pK(Planet.P_MPa[indsSil], Planet.T_K[indsSil], grid=False)
            # Evaluate silicate EOS for seismic properties
            Planet.Seismic.VP_kms[indsSil] = Planet.Sil.EOS.fn_VP_kms(Planet.P_MPa[indsSil], Planet.T_K[indsSil], grid=False)
            Planet.Seismic.VS_kms[indsSil] = Planet.Sil.EOS.fn_VS_kms(Planet.P_MPa[indsSil], Planet.T_K[indsSil], grid=False)
            Planet.Seismic.KS_GPa[indsSil] = Planet.Sil.EOS.fn_KS_GPa(Planet.P_MPa[indsSil], Planet.T_K[indsSil], grid=False)
            Planet.Seismic.GS_GPa[indsSil] = Planet.Sil.EOS.fn_GS_GPa(Planet.P_MPa[indsSil], Planet.T_K[indsSil], grid=False)
            if Planet.Do.POROUS_ROCK:
                VPpore_kms, KSpore_GPa = Planet.Ocean.EOS.fn_Seismic(Planet.Ppore_MPa[indsSil], Planet.T_K[indsSil])
                VSpore_kms, GSpore_GPa = (np.zeros_like(VPpore_kms) for _ in range(2))
                CpPore_JkgK = Planet.Ocean.EOS.fn_Cp_JkgK(Planet.Ppore_MPa[indsSil], Planet.T_K[indsSil], grid=False)
                alphaPore_pK = Planet.Ocean.EOS.fn_alpha_pK(Planet.Ppore_MPa[indsSil], Planet.T_K[indsSil], grid=False)
                Planet.Seismic.VP_kms[indsSil] = Planet.Sil.EOS.fn_porosCorrect(Planet.Seismic.VP_kms[indsSil], VPpore_kms, Planet.phi_frac[indsSil], Planet.Sil.JVP)
                Planet.Seismic.VS_kms[indsSil] = Planet.Sil.EOS.fn_porosCorrect(Planet.Seismic.VS_kms[indsSil], VSpore_kms, Planet.phi_frac[indsSil], Planet.Sil.JVS)
                Planet.Seismic.KS_GPa[indsSil] = Planet.Sil.EOS.fn_porosCorrect(Planet.Seismic.KS_GPa[indsSil], KSpore_GPa, Planet.phi_frac[indsSil], Planet.Sil.JKS)
                Planet.Seismic.GS_GPa[indsSil] = Planet.Sil.EOS.fn_porosCorrect(Planet.Seismic.GS_GPa[indsSil], GSpore_GPa, Planet.phi_frac[indsSil], Planet.Sil.JGS)
                Planet.Cp_JkgK[indsSil] = Planet.Sil.EOS.fn_porosCorrect(Planet.Cp_JkgK[indsSil]*Planet.rhoMatrix_kgm3[indsSil],
                        CpPore_JkgK*Planet.rhoPore_kgm3[indsSil], Planet.phi_frac[indsSil], Planet.Sil.JCp) / Planet.rho_kgm3[indsSil]
                Planet.alpha_pK[indsSil] = Planet.Sil.EOS.fn_porosCorrect(Planet.alpha_pK[indsSil], alphaPore_pK, Planet.phi_frac[indsSil], Planet.Sil.Jalpha)

            Hsil = Planet.Seismic.gSil * TsolidusHirschmann2000(Planet.P_MPa[indsSil])
            Planet.Seismic.QS[indsSil] = Planet.Seismic.BSil * np.exp(
                Planet.Seismic.gammaSil * Hsil / Planet.T_K[indsSil])

            if Planet.Do.Fe_CORE:
                # Evaluate core EOS for seismic properties
                Planet.Seismic.VP_kms[indsFe] = Planet.Core.EOS.fn_VP_kms(Planet.P_MPa[indsFe], Planet.T_K[indsFe], grid=False)
                Planet.Seismic.VS_kms[indsFe] = Planet.Core.EOS.fn_VS_kms(Planet.P_MPa[indsFe], Planet.T_K[indsFe], grid=False)
                Planet.Seismic.KS_GPa[indsFe] = Planet.Core.EOS.fn_KS_GPa(Planet.P_MPa[indsFe], Planet.T_K[indsFe], grid=False)
                Planet.Seismic.GS_GPa[indsFe] = Planet.Core.EOS.fn_GS_GPa(Planet.P_MPa[indsFe], Planet.T_K[indsFe], grid=False)
                if Planet.Seismic.QScore is not None:
                    Planet.Seismic.QS[indsFe] = Planet.Seismic.QScore
                else:
                    Planet.Seismic.QS[indsFe] = Constants.QScore

    return Planet
