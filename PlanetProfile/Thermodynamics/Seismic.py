import numpy as np
from PlanetProfile.Thermodynamics.HydroEOS import GetPhaseIndices
from PlanetProfile.Thermodynamics.InnerEOS import TsolidusHirschmann2000
from PlanetProfile.Utilities.defineStructs import Constants

def SeismicCalcs(Planet, Params):
    """ Calculation of seismic properties, including wave speeds

        Assigns Planet attributes:
            Seismic.VP_kms, Seismic.VS_kms, Seismic.QS, Seismic.KS_GPa, Seismic.GS_GPa
    """
    # Initialize arrays
    Planet.Seismic.VP_kms, Planet.Seismic.VS_kms, Planet.Seismic.QS, Planet.Seismic.KS_GPa, \
        Planet.Seismic.GS_GPa = (np.zeros(Planet.Steps.nTotal) for _ in range(5))

    if Params.CALC_SEISMIC and Planet.Do.VALID:

        indsLiq, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund, indsV, indsVund, indsVI, indsVIund, \
            indsClath, indsClathWet, indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, indsSilV, indsSilVI, \
            indsFe = GetPhaseIndices(Planet.phase)


        # Get seismic properties of all ice layers, starting with dry phases
        indsAllI = np.concatenate((indsI, indsIwet))
        if np.size(indsAllI) != 0:
            Planet.Seismic.VP_kms[indsAllI], Planet.Seismic.VS_kms[indsAllI], Planet.Seismic.KS_GPa[indsAllI], \
            Planet.Seismic.GS_GPa[indsAllI] = Planet.Ocean.surfIceEOS['Ih'].fn_Seismic(Planet.P_MPa[indsAllI], Planet.T_K[indsAllI])
            HiceI = Planet.Seismic.gIceI * Planet.Bulk.Tb_K
            Planet.Seismic.QS[indsAllI] = Planet.Seismic.BIceI * np.exp(
                Planet.Seismic.gammaIceI * HiceI / Planet.T_K[indsAllI])

        indsAllClath = np.concatenate((indsClath, indsClathWet))
        if np.size(indsAllClath) != 0:
            Planet.Seismic.VP_kms[indsAllClath], Planet.Seismic.VS_kms[indsAllClath], \
            Planet.Seismic.KS_GPa[indsAllClath], Planet.Seismic.GS_GPa[indsAllClath] \
                = Planet.Ocean.surfIceEOS['Clath'].fn_Seismic(Planet.P_MPa[indsAllClath], Planet.T_K[indsAllClath])
            Hclath = Planet.Seismic.gClath * np.max(Planet.T_K[indsAllClath])
            Planet.Seismic.QS[indsAllClath] = Planet.Seismic.BClath * np.exp(
                Planet.Seismic.gammaClath * Hclath / Planet.T_K[indsAllClath])

        indsAllII = np.concatenate((indsIIund, indsII))
        if np.size(indsAllII) != 0:
            if np.size(indsIIund) != 0:
                Planet.Seismic.VP_kms[indsIIund], Planet.Seismic.VS_kms[indsIIund], \
                Planet.Seismic.KS_GPa[indsIIund], Planet.Seismic.GS_GPa[indsIIund] \
                    = Planet.Ocean.surfIceEOS['II'].fn_Seismic(Planet.P_MPa[indsIIund], Planet.T_K[indsIIund])
            if np.size(indsII) != 0:
                Planet.Seismic.VP_kms[indsII], Planet.Seismic.VS_kms[indsII], \
                Planet.Seismic.KS_GPa[indsII], Planet.Seismic.GS_GPa[indsII] \
                    = Planet.Ocean.iceEOS['II'].fn_Seismic(Planet.P_MPa[indsII], Planet.T_K[indsII])
            HiceII = Planet.Seismic.gIceII * np.max(Planet.T_K[indsAllII])
            Planet.Seismic.QS[indsAllII] = Planet.Seismic.BIceII * np.exp(
                Planet.Seismic.gammaIceII * HiceII / Planet.T_K[indsAllII])

        indsAllIII = np.concatenate((indsIIIund, indsIII))
        if np.size(indsAllIII) != 0:
            if np.size(indsIIIund) != 0:
                Planet.Seismic.VP_kms[indsIIIund], Planet.Seismic.VS_kms[indsIIIund], \
                Planet.Seismic.KS_GPa[indsIIIund], Planet.Seismic.GS_GPa[indsIIIund] \
                    = Planet.Ocean.surfIceEOS['III'].fn_Seismic(Planet.P_MPa[indsIIIund], Planet.T_K[indsIIIund])
            if np.size(indsIII) != 0:
                Planet.Seismic.VP_kms[indsIII], Planet.Seismic.VS_kms[indsIII], \
                Planet.Seismic.KS_GPa[indsIII], Planet.Seismic.GS_GPa[indsIII] \
                    = Planet.Ocean.iceEOS['III'].fn_Seismic(Planet.P_MPa[indsIII], Planet.T_K[indsIII])
            HiceIII = Planet.Seismic.gIceIII * np.max(Planet.T_K[indsAllIII])
            Planet.Seismic.QS[indsAllIII] = Planet.Seismic.BIceIII * np.exp(
                Planet.Seismic.gammaIceIII * HiceIII / Planet.T_K[indsAllIII])

        indsAllV = np.concatenate((indsVund, indsV))
        if np.size(indsAllV) != 0:
            if np.size(indsVund) != 0:
                Planet.Seismic.VP_kms[indsVund], Planet.Seismic.VS_kms[indsVund], \
                Planet.Seismic.KS_GPa[indsVund], Planet.Seismic.GS_GPa[indsVund] \
                    = Planet.Ocean.surfIceEOS['V'].fn_Seismic(Planet.P_MPa[indsVund], Planet.T_K[indsVund])
            if np.size(indsV) != 0:
                Planet.Seismic.VP_kms[indsV], Planet.Seismic.VS_kms[indsV], \
                Planet.Seismic.KS_GPa[indsV], Planet.Seismic.GS_GPa[indsV] \
                    = Planet.Ocean.iceEOS['V'].fn_Seismic(Planet.P_MPa[indsV], Planet.T_K[indsV])
            HiceV = Planet.Seismic.gIceV * np.max(Planet.T_K[indsAllV])
            Planet.Seismic.QS[indsAllV] = Planet.Seismic.BIceV * np.exp(
                Planet.Seismic.gammaIceV * HiceV / Planet.T_K[indsAllV])

        indsAllVI = np.concatenate((indsVIund, indsVI))
        if np.size(indsAllVI) != 0:
            if np.size(indsVIund) != 0:
                Planet.Seismic.VP_kms[indsVIund], Planet.Seismic.VS_kms[indsVIund], \
                Planet.Seismic.KS_GPa[indsVIund], Planet.Seismic.GS_GPa[indsVIund] \
                    = Planet.Ocean.surfIceEOS['VI'].fn_Seismic(Planet.P_MPa[indsVIund], Planet.T_K[indsVIund])
            if np.size(indsVI) != 0:
                Planet.Seismic.VP_kms[indsVI], Planet.Seismic.VS_kms[indsVI], \
                Planet.Seismic.KS_GPa[indsVI], Planet.Seismic.GS_GPa[indsVI] \
                    = Planet.Ocean.iceEOS['VI'].fn_Seismic(Planet.P_MPa[indsVI], Planet.T_K[indsVI])
            HiceVI = Planet.Seismic.gIceVI * np.max(Planet.T_K[indsAllVI])
            Planet.Seismic.QS[indsAllVI] = Planet.Seismic.BIceVI * np.exp(
                Planet.Seismic.gammaIceVI * HiceVI / Planet.T_K[indsAllVI])

        if Planet.Do.POROUS_ICE:
            Planet = CalcSeisPorIce(Planet, Params, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund, indsV,
                                    indsVund, indsVI, indsVIund, indsClath, indsClathWet)
            
        # Get seismic properties of ocean layers
        if np.size(indsLiq) != 0:
            Planet.Seismic.VP_kms[indsLiq], Planet.Seismic.KS_GPa[indsLiq] \
                = Planet.Ocean.EOS.fn_Seismic(Planet.P_MPa[indsLiq], Planet.T_K[indsLiq])
            # VS_kms, QS, and GS_GPa were all initialized as zero and are all zero for these
            # layers, so we just leave them as-is.

        if not Params.SKIP_INNER:
            # Get Cp and alpha here, because we didn't calculate them earlier since we didn't need them
            # in calculating a conductive profile in the silicates and it would contribute extra,
            # unnecessary computational overhead there.
            Planet.Cp_JkgK[indsSil] = Planet.Sil.EOS.fn_Cp_JkgK(Planet.P_MPa[indsSil], Planet.T_K[indsSil])
            Planet.alpha_pK[indsSil] = Planet.Sil.EOS.fn_alpha_pK(Planet.P_MPa[indsSil], Planet.T_K[indsSil])
            # Evaluate silicate EOS for seismic properties
            Planet.Seismic.VP_kms[indsSil] = Planet.Sil.EOS.fn_VP_kms(Planet.P_MPa[indsSil], Planet.T_K[indsSil])
            Planet.Seismic.VS_kms[indsSil] = Planet.Sil.EOS.fn_VS_kms(Planet.P_MPa[indsSil], Planet.T_K[indsSil])
            Planet.Seismic.KS_GPa[indsSil] = Planet.Sil.EOS.fn_KS_GPa(Planet.P_MPa[indsSil], Planet.T_K[indsSil])
            Planet.Seismic.GS_GPa[indsSil] = Planet.Sil.EOS.fn_GS_GPa(Planet.P_MPa[indsSil], Planet.T_K[indsSil])
            Hsil = Planet.Seismic.gSil * TsolidusHirschmann2000(Planet.P_MPa[indsSil])
            Planet.Seismic.QS[indsSil] = Planet.Seismic.BSil * np.exp(
                Planet.Seismic.gammaSil * Hsil / Planet.T_K[indsSil])

            if Planet.Do.POROUS_ROCK:
                Planet = CalcSeisPorRock(Planet, Params, indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, indsSilV, indsSilVI)

            if Planet.Do.Fe_CORE:
                # Evaluate core EOS for seismic properties
                Planet.Seismic.VP_kms[indsFe] = Planet.Core.EOS.fn_VP_kms(Planet.P_MPa[indsFe], Planet.T_K[indsFe])
                Planet.Seismic.VS_kms[indsFe] = Planet.Core.EOS.fn_VS_kms(Planet.P_MPa[indsFe], Planet.T_K[indsFe])
                Planet.Seismic.KS_GPa[indsFe] = Planet.Core.EOS.fn_KS_GPa(Planet.P_MPa[indsFe], Planet.T_K[indsFe])
                Planet.Seismic.GS_GPa[indsFe] = Planet.Core.EOS.fn_GS_GPa(Planet.P_MPa[indsFe], Planet.T_K[indsFe])
                if Planet.Seismic.QScore is not None:
                    Planet.Seismic.QS[indsFe] = Planet.Seismic.QScore
                else:
                    Planet.Seismic.QS[indsFe] = Constants.QScore

    return Planet


def CalcSeisPorRock(Planet, Params, indsSil, indsSilLiq, indsSilI, indsSilII, indsSilIII, indsSilV, indsSilVI):

    VPpore_kms, VSpore_kms, KSpore_GPa, GSpore_GPa, CpPore_JkgK, alphaPore_pK \
        = (np.zeros_like(Planet.Ppore_MPa) for _ in range(6))

    if np.size(indsSilLiq) != 0:
        VPpore_kms[indsSilLiq], KSpore_GPa[indsSilLiq] \
            = Planet.Ocean.EOS.fn_Seismic(Planet.Ppore_MPa[indsSilLiq], Planet.T_K[indsSilLiq])
        # VS and GS are already initialized to 0, so we can just leave them be.
        CpPore_JkgK[indsSilLiq] = Planet.Ocean.EOS.fn_Cp_JkgK(Planet.Ppore_MPa[indsSilLiq], Planet.T_K[indsSilLiq])
        alphaPore_pK[indsSilLiq] = Planet.Ocean.EOS.fn_alpha_pK(Planet.Ppore_MPa[indsSilLiq], Planet.T_K[indsSilLiq])

    if np.size(indsSilI) != 0:
        VPpore_kms[indsSilI], VSpore_kms[indsSilI], KSpore_GPa[indsSilI], GSpore_GPa[indsSilI] \
            = Planet.Ocean.surfIceEOS['Ih'].fn_Seismic(Planet.Ppore_MPa[indsSilI], Planet.T_K[indsSilI])
        CpPore_JkgK[indsSilI] = Planet.Ocean.surfIceEOS['Ih'].fn_Cp_JkgK(Planet.Ppore_MPa[indsSilI], Planet.T_K[indsSilI])
        alphaPore_pK[indsSilI] = Planet.Ocean.surfIceEOS['Ih'].fn_alpha_pK(Planet.Ppore_MPa[indsSilI], Planet.T_K[indsSilI])

    if np.size(indsSilII) != 0:
        VPpore_kms[indsSilII], VSpore_kms[indsSilII], KSpore_GPa[indsSilII], GSpore_GPa[indsSilII] \
            = Planet.Ocean.iceEOS['II'].fn_Seismic(Planet.Ppore_MPa[indsSilII], Planet.T_K[indsSilII])
        CpPore_JkgK[indsSilII] = Planet.Ocean.iceEOS['II'].fn_Cp_JkgK(Planet.Ppore_MPa[indsSilII], Planet.T_K[indsSilII])
        alphaPore_pK[indsSilII] = Planet.Ocean.iceEOS['II'].fn_alpha_pK(Planet.Ppore_MPa[indsSilII], Planet.T_K[indsSilII])

    if np.size(indsSilIII) != 0:
        VPpore_kms[indsSilIII], VSpore_kms[indsSilIII], KSpore_GPa[indsSilIII], GSpore_GPa[indsSilIII] \
            = Planet.Ocean.iceEOS['III'].fn_Seismic(Planet.Ppore_MPa[indsSilIII], Planet.T_K[indsSilIII])
        CpPore_JkgK[indsSilIII] = Planet.Ocean.iceEOS['III'].fn_Cp_JkgK(Planet.Ppore_MPa[indsSilIII], Planet.T_K[indsSilIII])
        alphaPore_pK[indsSilIII] = Planet.Ocean.iceEOS['III'].fn_alpha_pK(Planet.Ppore_MPa[indsSilIII], Planet.T_K[indsSilIII])

    if np.size(indsSilV) != 0:
        VPpore_kms[indsSilV], VSpore_kms[indsSilV], KSpore_GPa[indsSilV], GSpore_GPa[indsSilV] \
            = Planet.Ocean.iceEOS['V'].fn_Seismic(Planet.Ppore_MPa[indsSilV], Planet.T_K[indsSilV])
        CpPore_JkgK[indsSilV] = Planet.Ocean.iceEOS['V'].fn_Cp_JkgK(Planet.Ppore_MPa[indsSilV], Planet.T_K[indsSilV])
        alphaPore_pK[indsSilV] = Planet.Ocean.iceEOS['V'].fn_alpha_pK(Planet.Ppore_MPa[indsSilV], Planet.T_K[indsSilV])

    if np.size(indsSilVI) != 0:
        VPpore_kms[indsSilVI], VSpore_kms[indsSilVI], KSpore_GPa[indsSilVI], GSpore_GPa[indsSilVI] \
            = Planet.Ocean.iceEOS['VI'].fn_Seismic(Planet.Ppore_MPa[indsSilVI], Planet.T_K[indsSilVI])
        CpPore_JkgK[indsSilVI] = Planet.Ocean.iceEOS['VI'].fn_Cp_JkgK(Planet.Ppore_MPa[indsSilVI], Planet.T_K[indsSilVI])
        alphaPore_pK[indsSilVI] = Planet.Ocean.iceEOS['VI'].fn_alpha_pK(Planet.Ppore_MPa[indsSilVI], Planet.T_K[indsSilVI])

    # Finally, combine the pore properties with the matrix properties across all phases
    Planet.Seismic.VP_kms[indsSil] = Planet.Sil.EOS.fn_porosCorrect(
        Planet.Seismic.VP_kms[indsSil], VPpore_kms[indsSil], Planet.phi_frac[indsSil], Planet.Sil.JVP)
    Planet.Seismic.VS_kms[indsSil] = Planet.Sil.EOS.fn_porosCorrect(
        Planet.Seismic.VS_kms[indsSil], VSpore_kms[indsSil], Planet.phi_frac[indsSil], Planet.Sil.JVS)
    Planet.Seismic.KS_GPa[indsSil] = Planet.Sil.EOS.fn_porosCorrect(
        Planet.Seismic.KS_GPa[indsSil], KSpore_GPa[indsSil], Planet.phi_frac[indsSil], Planet.Sil.JKS)
    Planet.Seismic.GS_GPa[indsSil] = Planet.Sil.EOS.fn_porosCorrect(
        Planet.Seismic.GS_GPa[indsSil], GSpore_GPa[indsSil], Planet.phi_frac[indsSil], Planet.Sil.JGS)
    Planet.Cp_JkgK[indsSil] = Planet.Sil.EOS.fn_porosCorrect(Planet.Cp_JkgK[indsSil] * Planet.rhoMatrix_kgm3[indsSil],
                                                             CpPore_JkgK[indsSil] * Planet.rhoPore_kgm3[indsSil],
                                                             Planet.phi_frac[indsSil], Planet.Sil.JCp) / \
                                                             Planet.rho_kgm3[indsSil]
    Planet.alpha_pK[indsSil] = Planet.Sil.EOS.fn_porosCorrect(Planet.alpha_pK[indsSil], alphaPore_pK[indsSil],
                                                              Planet.phi_frac[indsSil], Planet.Sil.Jalpha)

    return Planet


def CalcSeisPorIce(Planet, Params, indsI, indsIwet, indsII, indsIIund, indsIII, indsIIIund, indsV,
                                   indsVund, indsVI, indsVIund, indsClath, indsClathWet):

    # Make corrections for dry phases first
    if np.size(indsI) != 0:
        Planet.Seismic.VP_kms[indsI] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.Seismic.VP_kms[indsI], 0, Planet.phi_frac[indsI], Planet.Ocean.JVP)
        Planet.Seismic.VS_kms[indsI] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.Seismic.VS_kms[indsI], 0, Planet.phi_frac[indsI], Planet.Ocean.JVS)
        Planet.Seismic.KS_GPa[indsI] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.Seismic.KS_GPa[indsI], 0, Planet.phi_frac[indsI], Planet.Ocean.JKS)
        Planet.Seismic.GS_GPa[indsI] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.Seismic.GS_GPa[indsI], 0, Planet.phi_frac[indsI], Planet.Ocean.JGS)

    if np.size(indsClath) != 0:
        Planet.Seismic.VP_kms[indsClath] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(
            Planet.Seismic.VP_kms[indsClath], 0, Planet.phi_frac[indsClath], Planet.Ocean.JVP)
        Planet.Seismic.VS_kms[indsClath] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(
            Planet.Seismic.VS_kms[indsClath], 0, Planet.phi_frac[indsClath], Planet.Ocean.JVS)
        Planet.Seismic.KS_GPa[indsClath] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(
            Planet.Seismic.KS_GPa[indsClath], 0, Planet.phi_frac[indsClath], Planet.Ocean.JKS)
        Planet.Seismic.GS_GPa[indsClath] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(
            Planet.Seismic.GS_GPa[indsClath], 0, Planet.phi_frac[indsClath], Planet.Ocean.JGS)

    if np.size(indsIIund) != 0:
        Planet.Seismic.VP_kms[indsIIund] = Planet.Ocean.surfIceEOS['II'].fn_porosCorrect(
            Planet.Seismic.VP_kms[indsIIund], 0, Planet.phi_frac[indsIIund], Planet.Ocean.JVP)
        Planet.Seismic.VS_kms[indsIIund] = Planet.Ocean.surfIceEOS['II'].fn_porosCorrect(
            Planet.Seismic.VS_kms[indsIIund], 0, Planet.phi_frac[indsIIund], Planet.Ocean.JVS)
        Planet.Seismic.KS_GPa[indsIIund] = Planet.Ocean.surfIceEOS['II'].fn_porosCorrect(
            Planet.Seismic.KS_GPa[indsIIund], 0, Planet.phi_frac[indsIIund], Planet.Ocean.JKS)
        Planet.Seismic.GS_GPa[indsIIund] = Planet.Ocean.surfIceEOS['II'].fn_porosCorrect(
            Planet.Seismic.GS_GPa[indsIIund], 0, Planet.phi_frac[indsIIund], Planet.Ocean.JGS)

    if np.size(indsIIIund) != 0:
        Planet.Seismic.VP_kms[indsIIIund] = Planet.Ocean.surfIceEOS['III'].fn_porosCorrect(
            Planet.Seismic.VP_kms[indsIIIund], 0, Planet.phi_frac[indsIIIund], Planet.Ocean.JVP)
        Planet.Seismic.VS_kms[indsIIIund] = Planet.Ocean.surfIceEOS['III'].fn_porosCorrect(
            Planet.Seismic.VS_kms[indsIIIund], 0, Planet.phi_frac[indsIIIund], Planet.Ocean.JVS)
        Planet.Seismic.KS_GPa[indsIIIund] = Planet.Ocean.surfIceEOS['III'].fn_porosCorrect(
            Planet.Seismic.KS_GPa[indsIIIund], 0, Planet.phi_frac[indsIIIund], Planet.Ocean.JKS)
        Planet.Seismic.GS_GPa[indsIIIund] = Planet.Ocean.surfIceEOS['III'].fn_porosCorrect(
            Planet.Seismic.GS_GPa[indsIIIund], 0, Planet.phi_frac[indsIIIund], Planet.Ocean.JGS)

    if np.size(indsVund) != 0:
        Planet.Seismic.VP_kms[indsVund] = Planet.Ocean.surfIceEOS['V'].fn_porosCorrect(
            Planet.Seismic.VP_kms[indsVund], 0, Planet.phi_frac[indsVund], Planet.Ocean.JVP)
        Planet.Seismic.VS_kms[indsVund] = Planet.Ocean.surfIceEOS['V'].fn_porosCorrect(
            Planet.Seismic.VS_kms[indsVund], 0, Planet.phi_frac[indsVund], Planet.Ocean.JVS)
        Planet.Seismic.KS_GPa[indsVund] = Planet.Ocean.surfIceEOS['V'].fn_porosCorrect(
            Planet.Seismic.KS_GPa[indsVund], 0, Planet.phi_frac[indsVund], Planet.Ocean.JKS)
        Planet.Seismic.GS_GPa[indsVund] = Planet.Ocean.surfIceEOS['V'].fn_porosCorrect(
            Planet.Seismic.GS_GPa[indsVund], 0, Planet.phi_frac[indsVund], Planet.Ocean.JGS)

    if np.size(indsVIund) != 0:
        Planet.Seismic.VP_kms[indsVIund] = Planet.Ocean.surfIceEOS['VI'].fn_porosCorrect(
            Planet.Seismic.VP_kms[indsVIund], 0, Planet.phi_frac[indsVIund], Planet.Ocean.JVP)
        Planet.Seismic.VS_kms[indsVIund] = Planet.Ocean.surfIceEOS['VI'].fn_porosCorrect(
            Planet.Seismic.VS_kms[indsVIund], 0, Planet.phi_frac[indsVIund], Planet.Ocean.JVS)
        Planet.Seismic.KS_GPa[indsVIund] = Planet.Ocean.surfIceEOS['VI'].fn_porosCorrect(
            Planet.Seismic.KS_GPa[indsVIund], 0, Planet.phi_frac[indsVIund], Planet.Ocean.JKS)
        Planet.Seismic.GS_GPa[indsVIund] = Planet.Ocean.surfIceEOS['VI'].fn_porosCorrect(
            Planet.Seismic.GS_GPa[indsVIund], 0, Planet.phi_frac[indsVIund], Planet.Ocean.JGS)

    # Now wet phases. First, fluid properties need to be calculated for each, then these
    # are used to correct for the phase mixture.
    if np.size(indsIwet) != 0:
        VPpore_kms, KSpore_GPa = Planet.Ocean.EOS.fn_Seismic(Planet.Ppore_MPa[indsIwet], Planet.T_K[indsIwet])
        Planet.Seismic.VP_kms[indsIwet] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.Seismic.VP_kms[indsIwet], VPpore_kms, Planet.phi_frac[indsIwet], Planet.Ocean.JVP)
        Planet.Seismic.VS_kms[indsIwet] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.Seismic.VS_kms[indsIwet], 0, Planet.phi_frac[indsIwet], Planet.Ocean.JVS)
        Planet.Seismic.KS_GPa[indsIwet] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.Seismic.KS_GPa[indsIwet], KSpore_GPa, Planet.phi_frac[indsIwet], Planet.Ocean.JKS)
        Planet.Seismic.GS_GPa[indsIwet] = Planet.Ocean.surfIceEOS['Ih'].fn_porosCorrect(
            Planet.Seismic.GS_GPa[indsIwet], 0, Planet.phi_frac[indsIwet], Planet.Ocean.JGS)

    if np.size(indsClathWet) != 0:
        VPpore_kms, KSpore_GPa = Planet.Ocean.EOS.fn_Seismic(Planet.Ppore_MPa[indsClathWet], Planet.T_K[indsClathWet])
        Planet.Seismic.VP_kms[indsClathWet] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(
            Planet.Seismic.VP_kms[indsClathWet], VPpore_kms, Planet.phi_frac[indsClathWet], Planet.Ocean.JVP)
        Planet.Seismic.VS_kms[indsClathWet] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(
            Planet.Seismic.VS_kms[indsClathWet], 0, Planet.phi_frac[indsClathWet], Planet.Ocean.JVS)
        Planet.Seismic.KS_GPa[indsClathWet] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(
            Planet.Seismic.KS_GPa[indsClathWet], KSpore_GPa, Planet.phi_frac[indsClathWet], Planet.Ocean.JKS)
        Planet.Seismic.GS_GPa[indsClathWet] = Planet.Ocean.surfIceEOS['Clath'].fn_porosCorrect(
            Planet.Seismic.GS_GPa[indsClathWet], 0, Planet.phi_frac[indsClathWet], Planet.Ocean.JGS)

    if np.size(indsII) != 0:
        VPpore_kms, KSpore_GPa = Planet.Ocean.EOS.fn_Seismic(Planet.Ppore_MPa[indsII], Planet.T_K[indsII])
        Planet.Seismic.VP_kms[indsII] = Planet.Ocean.iceEOS['II'].fn_porosCorrect(
            Planet.Seismic.VP_kms[indsII], VPpore_kms, Planet.phi_frac[indsII], Planet.Ocean.JVP)
        Planet.Seismic.VS_kms[indsII] = Planet.Ocean.iceEOS['II'].fn_porosCorrect(
            Planet.Seismic.VS_kms[indsII], 0, Planet.phi_frac[indsII], Planet.Ocean.JVS)
        Planet.Seismic.KS_GPa[indsII] = Planet.Ocean.iceEOS['II'].fn_porosCorrect(
            Planet.Seismic.KS_GPa[indsII], KSpore_GPa, Planet.phi_frac[indsII], Planet.Ocean.JKS)
        Planet.Seismic.GS_GPa[indsII] = Planet.Ocean.iceEOS['II'].fn_porosCorrect(
            Planet.Seismic.GS_GPa[indsII], 0, Planet.phi_frac[indsII], Planet.Ocean.JGS)

    if np.size(indsIII) != 0:
        VPpore_kms, KSpore_GPa = Planet.Ocean.EOS.fn_Seismic(Planet.Ppore_MPa[indsIII], Planet.T_K[indsIII])
        Planet.Seismic.VP_kms[indsIII] = Planet.Ocean.iceEOS['III'].fn_porosCorrect(
            Planet.Seismic.VP_kms[indsIII], VPpore_kms, Planet.phi_frac[indsIII], Planet.Ocean.JVP)
        Planet.Seismic.VS_kms[indsIII] = Planet.Ocean.iceEOS['III'].fn_porosCorrect(
            Planet.Seismic.VS_kms[indsIII], 0, Planet.phi_frac[indsIII], Planet.Ocean.JVS)
        Planet.Seismic.KS_GPa[indsIII] = Planet.Ocean.iceEOS['III'].fn_porosCorrect(
            Planet.Seismic.KS_GPa[indsIII], KSpore_GPa, Planet.phi_frac[indsIII], Planet.Ocean.JKS)
        Planet.Seismic.GS_GPa[indsIII] = Planet.Ocean.iceEOS['III'].fn_porosCorrect(
            Planet.Seismic.GS_GPa[indsIII], 0, Planet.phi_frac[indsIII], Planet.Ocean.JGS)

    if np.size(indsV) != 0:
        VPpore_kms, KSpore_GPa = Planet.Ocean.EOS.fn_Seismic(Planet.Ppore_MPa[indsV], Planet.T_K[indsV])
        Planet.Seismic.VP_kms[indsV] = Planet.Ocean.iceEOS['V'].fn_porosCorrect(
            Planet.Seismic.VP_kms[indsV], VPpore_kms, Planet.phi_frac[indsV], Planet.Ocean.JVP)
        Planet.Seismic.VS_kms[indsV] = Planet.Ocean.iceEOS['V'].fn_porosCorrect(
            Planet.Seismic.VS_kms[indsV], 0, Planet.phi_frac[indsV], Planet.Ocean.JVS)
        Planet.Seismic.KS_GPa[indsV] = Planet.Ocean.iceEOS['V'].fn_porosCorrect(
            Planet.Seismic.KS_GPa[indsV], KSpore_GPa, Planet.phi_frac[indsV], Planet.Ocean.JKS)
        Planet.Seismic.GS_GPa[indsV] = Planet.Ocean.iceEOS['V'].fn_porosCorrect(
            Planet.Seismic.GS_GPa[indsV], 0, Planet.phi_frac[indsV], Planet.Ocean.JGS)

    if np.size(indsVI) != 0:
        VPpore_kms, KSpore_GPa = Planet.Ocean.EOS.fn_Seismic(Planet.Ppore_MPa[indsVI], Planet.T_K[indsVI])
        Planet.Seismic.VP_kms[indsVI] = Planet.Ocean.iceEOS['VI'].fn_porosCorrect(
            Planet.Seismic.VP_kms[indsVI], VPpore_kms, Planet.phi_frac[indsVI], Planet.Ocean.JVP)
        Planet.Seismic.VS_kms[indsVI] = Planet.Ocean.iceEOS['VI'].fn_porosCorrect(
            Planet.Seismic.VS_kms[indsVI], 0, Planet.phi_frac[indsVI], Planet.Ocean.JVS)
        Planet.Seismic.KS_GPa[indsVI] = Planet.Ocean.iceEOS['VI'].fn_porosCorrect(
            Planet.Seismic.KS_GPa[indsVI], KSpore_GPa, Planet.phi_frac[indsVI], Planet.Ocean.JKS)
        Planet.Seismic.GS_GPa[indsVI] = Planet.Ocean.iceEOS['VI'].fn_porosCorrect(
            Planet.Seismic.GS_GPa[indsVI], 0, Planet.phi_frac[indsVI], Planet.Ocean.JGS)


    return Planet
