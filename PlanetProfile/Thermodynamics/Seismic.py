import numpy as np
from scipy.interpolate import interp1d
from PlanetProfile.Thermodynamics.HydroEOS import GetIceEOS
from PlanetProfile.Utilities.Indexing import GetPhaseIndices
from PlanetProfile.Thermodynamics.InnerEOS import TsolidusHirschmann2000
from PlanetProfile.Utilities.defineStructs import Constants, EOSlist
from PlanetProfile.Utilities.PPversion import ppVerNum
import logging

# Assign logger
log = logging.getLogger('PlanetProfile')

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
            # Get ice EOS if not currently loaded
            icePhase = 'Ih'
            if Planet.Ocean.surfIceEOS[icePhase].key not in EOSlist.loaded.keys():
                PIce_MPa = np.linspace(Planet.P_MPa[indsAllI][0], Planet.P_MPa[indsAllI][-1] + Planet.Ocean.deltaP * 3, np.maximum(np.size(indsAllI), 4))
                TIce_K = np.linspace(Planet.T_K[indsAllI][0], Planet.T_K[indsAllI][-1] + Planet.Ocean.deltaT * 3, np.maximum(np.size(indsAllI), 4))
                Planet.Ocean.surfIceEOS[icePhase] = GetIceEOS(PIce_MPa, TIce_K, icePhase,
                                                              porosType=Planet.Ocean.porosType[icePhase],
                                                              phiTop_frac=Planet.Ocean.phiMax_frac[icePhase],
                                                              Pclosure_MPa=Planet.Ocean.Pclosure_MPa[icePhase],
                                                              phiMin_frac=Planet.Ocean.phiMin_frac,
                                                              EXTRAP=Params.EXTRAP_ICE[icePhase],
                                                              ICEIh_DIFFERENT=Planet.Do.ICEIh_DIFFERENT)

            Planet.Seismic.VP_kms[indsAllI], Planet.Seismic.VS_kms[indsAllI], Planet.Seismic.KS_GPa[indsAllI], \
            Planet.Seismic.GS_GPa[indsAllI] = Planet.Ocean.surfIceEOS['Ih'].fn_Seismic(Planet.P_MPa[indsAllI], Planet.T_K[indsAllI])
            HiceI = Planet.Seismic.gIceI * Planet.Bulk.Tb_K
            Planet.Seismic.QS[indsAllI] = Planet.Seismic.BIceI * np.exp(
                Planet.Seismic.gammaIceI * HiceI / Planet.T_K[indsAllI])

        indsAllClath = np.concatenate((indsClath, indsClathWet))
        if np.size(indsAllClath) != 0:
            # Get ice EOS if not currently loaded
            icePhase = 'Clath'
            if Planet.Ocean.surfIceEOS[icePhase].key not in EOSlist.loaded.keys():
                PIce_MPa = np.linspace(Planet.P_MPa[indsAllClath][0], Planet.P_MPa[indsAllClath][-1] + Planet.Ocean.deltaP * 3, np.maximum(np.size(indsAllClath), 4))
                TIce_K = np.linspace(Planet.T_K[indsAllClath][0], Planet.T_K[indsAllClath][-1] + Planet.Ocean.deltaT * 3, np.maximum(np.size(indsAllClath), 4))
                Planet.Ocean.surfIceEOS[icePhase] = GetIceEOS(PIce_MPa, TIce_K, icePhase,
                                                              porosType=Planet.Ocean.porosType[icePhase],
                                                              phiTop_frac=Planet.Ocean.phiMax_frac[icePhase],
                                                              Pclosure_MPa=Planet.Ocean.Pclosure_MPa[icePhase],
                                                              phiMin_frac=Planet.Ocean.phiMin_frac,
                                                              EXTRAP=Params.EXTRAP_ICE[icePhase],
                                                              ClathDissoc=Planet.Ocean.ClathDissoc)

            Planet.Seismic.VP_kms[indsAllClath], Planet.Seismic.VS_kms[indsAllClath], \
            Planet.Seismic.KS_GPa[indsAllClath], Planet.Seismic.GS_GPa[indsAllClath] \
                = Planet.Ocean.surfIceEOS['Clath'].fn_Seismic(Planet.P_MPa[indsAllClath], Planet.T_K[indsAllClath])
            Hclath = Planet.Seismic.gClath * np.max(Planet.T_K[indsAllClath])
            Planet.Seismic.QS[indsAllClath] = Planet.Seismic.BClath * np.exp(
                Planet.Seismic.gammaClath * Hclath / Planet.T_K[indsAllClath])

        indsAllII = np.concatenate((indsIIund, indsII))
        if np.size(indsAllII) != 0:
            icePhase = 'II'
            if np.size(indsIIund) != 0:
                # Get ice EOS if not currently loaded
                if Planet.Ocean.surfIceEOS[icePhase].key not in EOSlist.loaded.keys():
                    PIce_MPa = np.linspace(Planet.P_MPa[indsIIund][0], Planet.P_MPa[indsIIund][-1] + Planet.Ocean.deltaP * 3, np.maximum(np.size(indsIIund), 4))
                    TIce_K = np.linspace(Planet.T_K[indsIIund][0], Planet.T_K[indsIIund][-1] + Planet.Ocean.deltaT * 3, np.maximum(np.size(indsIIund), 4))
                    Planet.Ocean.surfIceEOS[icePhase] = GetIceEOS(PIce_MPa, TIce_K, icePhase,
                                                                  porosType=Planet.Ocean.porosType[icePhase],
                                                                  phiTop_frac=Planet.Ocean.phiMax_frac[icePhase],
                                                                  Pclosure_MPa=Planet.Ocean.Pclosure_MPa[icePhase],
                                                                  phiMin_frac=Planet.Ocean.phiMin_frac,
                                                                  EXTRAP=Params.EXTRAP_ICE[icePhase])

                Planet.Seismic.VP_kms[indsIIund], Planet.Seismic.VS_kms[indsIIund], \
                Planet.Seismic.KS_GPa[indsIIund], Planet.Seismic.GS_GPa[indsIIund] \
                    = Planet.Ocean.surfIceEOS['II'].fn_Seismic(Planet.P_MPa[indsIIund], Planet.T_K[indsIIund])
            if np.size(indsII) != 0:
                # Get ice EOS if not currently loaded
                if Planet.Ocean.iceEOS[icePhase].key not in EOSlist.loaded.keys():
                    PIce_MPa = np.linspace(Planet.P_MPa[indsII][0], Planet.P_MPa[indsII][-1] + Planet.Ocean.deltaP * 3, np.maximum(np.size(indsII), 4))
                    TIce_K = np.linspace(Planet.T_K[indsII][0], Planet.T_K[indsII][-1] + Planet.Ocean.deltaT * 3, np.maximum(np.size(indsII), 4))
                    Planet.Ocean.iceEOS[icePhase] = GetIceEOS(PIce_MPa, TIce_K, icePhase,
                                                              porosType=Planet.Ocean.porosType[icePhase],
                                                              phiTop_frac=Planet.Ocean.phiMax_frac[icePhase],
                                                              Pclosure_MPa=Planet.Ocean.Pclosure_MPa[icePhase],
                                                              phiMin_frac=Planet.Ocean.phiMin_frac,
                                                              EXTRAP=Params.EXTRAP_ICE[icePhase])

                Planet.Seismic.VP_kms[indsII], Planet.Seismic.VS_kms[indsII], \
                Planet.Seismic.KS_GPa[indsII], Planet.Seismic.GS_GPa[indsII] \
                    = Planet.Ocean.iceEOS['II'].fn_Seismic(Planet.P_MPa[indsII], Planet.T_K[indsII])
            HiceII = Planet.Seismic.gIceII * np.max(Planet.T_K[indsAllII])
            Planet.Seismic.QS[indsAllII] = Planet.Seismic.BIceII * np.exp(
                Planet.Seismic.gammaIceII * HiceII / Planet.T_K[indsAllII])

        indsAllIII = np.concatenate((indsIIIund, indsIII))
        if np.size(indsAllIII) != 0:
            icePhase = 'III'
            if np.size(indsIIIund) != 0:
                # Get ice EOS if not currently loaded
                if Planet.Ocean.surfIceEOS[icePhase].key not in EOSlist.loaded.keys():
                    PIce_MPa = np.linspace(Planet.P_MPa[indsIIIund][0], Planet.P_MPa[indsIIIund][-1] + Planet.Ocean.deltaP * 3, np.maximum(np.size(indsIIIund), 4))
                    TIce_K = np.linspace(Planet.T_K[indsIIIund][0], Planet.T_K[indsIIIund][-1] + Planet.Ocean.deltaT * 3, np.maximum(np.size(indsIIIund), 4))
                    Planet.Ocean.surfIceEOS[icePhase] = GetIceEOS(PIce_MPa, TIce_K, icePhase,
                                                                  porosType=Planet.Ocean.porosType[icePhase],
                                                                  phiTop_frac=Planet.Ocean.phiMax_frac[icePhase],
                                                                  Pclosure_MPa=Planet.Ocean.Pclosure_MPa[icePhase],
                                                                  phiMin_frac=Planet.Ocean.phiMin_frac,
                                                                  EXTRAP=Params.EXTRAP_ICE[icePhase])

                Planet.Seismic.VP_kms[indsIIIund], Planet.Seismic.VS_kms[indsIIIund], \
                Planet.Seismic.KS_GPa[indsIIIund], Planet.Seismic.GS_GPa[indsIIIund] \
                    = Planet.Ocean.surfIceEOS['III'].fn_Seismic(Planet.P_MPa[indsIIIund], Planet.T_K[indsIIIund])
            if np.size(indsIII) != 0:
                # Get ice EOS if not currently loaded
                if Planet.Ocean.iceEOS[icePhase].key not in EOSlist.loaded.keys():
                    PIce_MPa = np.linspace(Planet.P_MPa[indsIII][0], Planet.P_MPa[indsIII][-1] + Planet.Ocean.deltaP * 3, np.maximum(np.size(indsIII), 4))
                    TIce_K = np.linspace(Planet.T_K[indsIII][0], Planet.T_K[indsIII][-1] + Planet.Ocean.deltaT * 3, np.maximum(np.size(indsIII), 4))
                    Planet.Ocean.iceEOS[icePhase] = GetIceEOS(PIce_MPa, TIce_K, icePhase,
                                                              porosType=Planet.Ocean.porosType[icePhase],
                                                              phiTop_frac=Planet.Ocean.phiMax_frac[icePhase],
                                                              Pclosure_MPa=Planet.Ocean.Pclosure_MPa[icePhase],
                                                              phiMin_frac=Planet.Ocean.phiMin_frac,
                                                              EXTRAP=Params.EXTRAP_ICE[icePhase])

                Planet.Seismic.VP_kms[indsIII], Planet.Seismic.VS_kms[indsIII], \
                Planet.Seismic.KS_GPa[indsIII], Planet.Seismic.GS_GPa[indsIII] \
                    = Planet.Ocean.iceEOS['III'].fn_Seismic(Planet.P_MPa[indsIII], Planet.T_K[indsIII])
            HiceIII = Planet.Seismic.gIceIII * np.max(Planet.T_K[indsAllIII])
            Planet.Seismic.QS[indsAllIII] = Planet.Seismic.BIceIII * np.exp(
                Planet.Seismic.gammaIceIII * HiceIII / Planet.T_K[indsAllIII])

        indsAllV = np.concatenate((indsVund, indsV))
        if np.size(indsAllV) != 0:
            icePhase ='V'
            if np.size(indsVund) != 0:
                # Get ice EOS if not currently loaded
                if Planet.Ocean.surfIceEOS[icePhase].key not in EOSlist.loaded.keys():
                    PIce_MPa = np.linspace(Planet.P_MPa[indsVund][0], Planet.P_MPa[indsVund][-1] + Planet.Ocean.deltaP * 3, np.maximum(np.size(indsVund), 4))
                    TIce_K = np.linspace(Planet.T_K[indsVund][0], Planet.T_K[indsVund][-1] + Planet.Ocean.deltaT * 3, np.maximum(np.size(indsVund), 4))
                    Planet.Ocean.surfIceEOS[icePhase] = GetIceEOS(PIce_MPa, TIce_K, icePhase,
                                                                  porosType=Planet.Ocean.porosType[icePhase],
                                                                  phiTop_frac=Planet.Ocean.phiMax_frac[icePhase],
                                                                  Pclosure_MPa=Planet.Ocean.Pclosure_MPa[icePhase],
                                                                  phiMin_frac=Planet.Ocean.phiMin_frac,
                                                                  EXTRAP=Params.EXTRAP_ICE[icePhase])

                Planet.Seismic.VP_kms[indsVund], Planet.Seismic.VS_kms[indsVund], \
                Planet.Seismic.KS_GPa[indsVund], Planet.Seismic.GS_GPa[indsVund] \
                    = Planet.Ocean.surfIceEOS['V'].fn_Seismic(Planet.P_MPa[indsVund], Planet.T_K[indsVund])
            if np.size(indsV) != 0:
                # Get ice EOS if not currently loaded
                if Planet.Ocean.iceEOS[icePhase].key not in EOSlist.loaded.keys():
                    PIce_MPa = np.linspace(Planet.P_MPa[indsV][0], Planet.P_MPa[indsV][-1] + Planet.Ocean.deltaP * 3, np.maximum(np.size(indsV), 4))
                    TIce_K = np.linspace(Planet.T_K[indsV][0], Planet.T_K[indsV][-1] + Planet.Ocean.deltaT * 3, np.maximum(np.size(indsV), 4))
                    Planet.Ocean.iceEOS[icePhase] = GetIceEOS(PIce_MPa, TIce_K, icePhase,
                                                              porosType=Planet.Ocean.porosType[icePhase],
                                                              phiTop_frac=Planet.Ocean.phiMax_frac[icePhase],
                                                              Pclosure_MPa=Planet.Ocean.Pclosure_MPa[icePhase],
                                                              phiMin_frac=Planet.Ocean.phiMin_frac,
                                                              EXTRAP=Params.EXTRAP_ICE[icePhase])

                Planet.Seismic.VP_kms[indsV], Planet.Seismic.VS_kms[indsV], \
                Planet.Seismic.KS_GPa[indsV], Planet.Seismic.GS_GPa[indsV] \
                    = Planet.Ocean.iceEOS['V'].fn_Seismic(Planet.P_MPa[indsV], Planet.T_K[indsV])
            HiceV = Planet.Seismic.gIceV * np.max(Planet.T_K[indsAllV])
            Planet.Seismic.QS[indsAllV] = Planet.Seismic.BIceV * np.exp(
                Planet.Seismic.gammaIceV * HiceV / Planet.T_K[indsAllV])

        indsAllVI = np.concatenate((indsVIund, indsVI))
        if np.size(indsAllVI) != 0:
            icePhase = 'VI'
            if np.size(indsVIund) != 0:
                # Get ice EOS if not currently loaded
                if Planet.Ocean.surfIceEOS[icePhase].key not in EOSlist.loaded.keys():
                    PIce_MPa = np.linspace(Planet.P_MPa[indsVIund][0], Planet.P_MPa[indsVIund][-1] + Planet.Ocean.deltaP * 3, np.maximum(np.size(indsVIund), 4))
                    TIce_K = np.linspace(Planet.T_K[indsVIund][0], Planet.T_K[indsVIund][-1] + Planet.Ocean.deltaT * 3, np.maximum(np.size(indsVIund), 4))
                    Planet.Ocean.surfIceEOS[icePhase] = GetIceEOS(PIce_MPa, TIce_K, icePhase,
                                                                  porosType=Planet.Ocean.porosType[icePhase],
                                                                  phiTop_frac=Planet.Ocean.phiMax_frac[icePhase],
                                                                  Pclosure_MPa=Planet.Ocean.Pclosure_MPa[icePhase],
                                                                  phiMin_frac=Planet.Ocean.phiMin_frac,
                                                                  EXTRAP=Params.EXTRAP_ICE[icePhase])

                Planet.Seismic.VP_kms[indsVIund], Planet.Seismic.VS_kms[indsVIund], \
                Planet.Seismic.KS_GPa[indsVIund], Planet.Seismic.GS_GPa[indsVIund] \
                    = Planet.Ocean.surfIceEOS['VI'].fn_Seismic(Planet.P_MPa[indsVIund], Planet.T_K[indsVIund])
            if np.size(indsVI) != 0:
                # Get ice EOS if not currently loaded
                if Planet.Ocean.iceEOS[icePhase].key not in EOSlist.loaded.keys():
                    PIce_MPa = np.linspace(Planet.P_MPa[indsVI][0], Planet.P_MPa[indsVI][-1] + Planet.Ocean.deltaP * 3, np.maximum(np.size(indsVI), 4))
                    TIce_K = np.linspace(Planet.T_K[indsVI][0], Planet.T_K[indsVI][-1] + Planet.Ocean.deltaT * 3, np.maximum(np.size(indsVI), 4))
                    Planet.Ocean.iceEOS[icePhase] = GetIceEOS(PIce_MPa, TIce_K, icePhase,
                                                              porosType=Planet.Ocean.porosType[icePhase],
                                                              phiTop_frac=Planet.Ocean.phiMax_frac[icePhase],
                                                              Pclosure_MPa=Planet.Ocean.Pclosure_MPa[icePhase],
                                                              phiMin_frac=Planet.Ocean.phiMin_frac,
                                                              EXTRAP=Params.EXTRAP_ICE[icePhase])

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

    else:
        Planet.Ocean.GSmeanIIIwet_GPa = np.nan
        Planet.Ocean.GSmeanVwet_GPa = np.nan
        Planet.Ocean.GSmeanVI_GPa = np.nan
        Planet.Sil.GSmean_GPa = np.nan
        Planet.Core.GSmeanFe_GPa = np.nan
        Planet.Core.GSmeanFeS_GPa = np.nan

    if np.any(Planet.Seismic.QS > Planet.Seismic.QSmax):
        log.debug(f'Resetting unnecessarily high QS values to max value: {Planet.Seismic.QSmax}')
        Planet.Seismic.QS[Planet.Seismic.QS > Planet.Seismic.QSmax] = Planet.Seismic.QSmax

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


def WriteSeismic(Planet, Params):
    # Print outputs usable by seismic wave post-processing software

    """ @@@@@@
        AxiSEM
        @@@@@@
    """
    # Extend final layer to center
    rAxi_m = Planet.r_m + 0
    rhoAxi_kgm3 = np.append(Planet.rho_kgm3, Planet.rho_kgm3[-1])
    VPAxi_ms = np.append(Planet.Seismic.VP_kms*1e3, Planet.Seismic.VP_kms[-1]*1e3)
    VSAxi_ms = np.append(Planet.Seismic.VS_kms*1e3, Planet.Seismic.VS_kms[-1]*1e3)
    QSAxi = np.append(Planet.Seismic.QS, Planet.Seismic.QS[-1])
    Qkappa = Planet.Seismic.Qkappa * np.ones_like(rAxi_m)  # Qkappa is currently simply set to a constant.
    # Add entries for discontinuities -- note that atmosphere layers must be added manually.
    nDisc = 0
    for i, phase in enumerate(Planet.phase[1:]):
        if phase != Planet.phase[i]:
            # Check that we aren't just changing from one pore fluid phase to another
            if not ((phase >= Constants.phaseSil and phase < Constants.phaseSil + 10) and
                (Planet.phase[i] >= Constants.phaseSil and Planet.phase[i] < Constants.phaseSil + 10)):
                # Insert a duplicate entry at the bottom of each bulk layer
                rAxi_m = np.insert(rAxi_m, i+nDisc+1, rAxi_m[i+nDisc+1])
                rhoAxi_kgm3 = np.insert(rhoAxi_kgm3, i+nDisc+1, rhoAxi_kgm3[i+nDisc])
                VPAxi_ms = np.insert(VPAxi_ms, i+nDisc+1, VPAxi_ms[i+nDisc])
                VSAxi_ms = np.insert(VSAxi_ms, i+nDisc+1, VSAxi_ms[i+nDisc])
                QSAxi = np.insert(QSAxi, i+nDisc+1, QSAxi[i+nDisc])
                Qkappa = np.insert(Qkappa, i+nDisc+1, Qkappa[i+nDisc])
                nDisc += 1

    # Construct header
    leadWSAxi = 2
    headerLines = [
        f'# Model for {Planet.name}',
        f'# Created with PlanetProfile v{ppVerNum}:',
        f'# {Planet.label}',
        f'NAME         {Planet.name.lower()}',
        f'ANELASTIC    T',
        f'ANISOTROPIC  F',
        f'UNITS        m',
        f'COLUMNS     ' + ' '.join([
            'radius',
            'rho',
            'vpv',
            'vsv',
            'qka',
            'qmu'
        ])
    ]
    with open(Params.DataFiles.AxiSEMfile,'w') as f:
        f.write('\n'.join(headerLines) + '\n')
        for i in range(np.size(rAxi_m)):
            f.write(' '*leadWSAxi + ' '.join([
                f'{rAxi_m[i]:.1f}',
                f'{rhoAxi_kgm3[i]:.2f}',
                f'{VPAxi_ms[i]:.2f}',
                f'{VSAxi_ms[i]:.2f}',
                f'{Qkappa[i]:.1f}',
                f'{QSAxi[i]:.1f}'
            ]) + '\n')


    """ @@@@@@@@@@@@@@@
        minEOS velmodel
        @@@@@@@@@@@@@@@
    """
    # Reconfigure data into minEOS-appropriate format
    rminEOSpre_m = np.arange(0, Planet.r_m[0] + Planet.Seismic.minEOS_rRes_m, Planet.Seismic.minEOS_rRes_m)
    zminEOSpre_m = Planet.Bulk.R_m - rminEOSpre_m
    # Save copies we can then edit
    z_m, rho_kgm3, VP_ms, VS_ms, KS_GPa, QS = (Planet.z_m, Planet.rho_kgm3, Planet.Seismic.VP_kms*1e3,
                                               Planet.Seismic.VS_kms*1e3, Planet.Seismic.KS_GPa, Planet.Seismic.QS)
    # Set interpolation kwargs to extrapolate to finer resolution
    intArgs = {'kind': 'cubic', 'bounds_error': False, 'fill_value': 'extrapolate'}

    # Adjust layer boundaries at top and bottom of ocean to avoid "mushy" layers
    # that combine properties in the interpolation.
    if Planet.Do.NO_H2O or not np.any(Planet.phase == 0):
        # Simple case -- no ocean, just interpolate.
        rminEOS_m = rminEOSpre_m
        zminEOS_m = zminEOSpre_m
        nminEOS = np.size(rminEOS_m)
        iOceanTop, iOceanBot = (nminEOS, nminEOS)
        # Interpolate at desired output radii/depths
        rhominEOS_kgm3 = interp1d(z_m[:-1], rho_kgm3, **intArgs)(zminEOS_m)
        VPminEOS_ms =    interp1d(z_m[:-1], VP_ms,    **intArgs)(zminEOS_m)
        VSminEOS_ms =    interp1d(z_m[:-1], VS_ms,    **intArgs)(zminEOS_m)
        QkappaminEOS =   interp1d(z_m[:-1], KS_GPa,   **intArgs)(zminEOS_m)
        QmuminEOS =      interp1d(z_m[:-1], QS,       **intArgs)(zminEOS_m)
    else:
        # This only identifies the top and bottom of the contiguous ocean layer.
        # If the ocean is separated into further layers between deep-sea HP ices,
        # additional detail will be required in creating minEOS input files to
        # handle the multiple liquid layers.
        iObotPP = next(i for i,phase in enumerate(Planet.phase[:Planet.Steps.nHydro]) if phase == 0 and phase != Planet.phase[i+1])
        iInner = np.where(rminEOSpre_m <= Planet.r_m[iObotPP+1])[0]
        iShell = np.where(rminEOSpre_m > Planet.r_m[Planet.Steps.nSurfIce])[0]
        iOcean = np.hstack((iInner[-1], np.where(np.logical_and(rminEOSpre_m > Planet.r_m[iObotPP+1],
                                                         rminEOSpre_m <= Planet.r_m[Planet.Steps.nSurfIce]))[0], iShell[0]))
        iOceanBot = iOcean[0]
        iOceanTop = iOcean[-1]

        zInner_m, zOcean_m, zShell_m = (zminEOSpre_m[iInner], zminEOSpre_m[iOcean], zminEOSpre_m[iShell])
        rminEOS_m = np.hstack((rminEOSpre_m[iInner], rminEOSpre_m[iOcean], rminEOSpre_m[iShell]))
        nminEOS = np.size(rminEOS_m)

        # Interpolate layers below ocean
        rhoInner_kgm3 = interp1d(z_m[iObotPP+1:-1], rho_kgm3[iObotPP+1:], **intArgs)(zInner_m)
        VPinner_ms = interp1d(z_m[iObotPP+1:-1], VP_ms[iObotPP+1:],       **intArgs)(zInner_m)
        VSinner_ms = interp1d(z_m[iObotPP+1:-1], VS_ms[iObotPP+1:],       **intArgs)(zInner_m)
        QkappaInner = interp1d(z_m[iObotPP+1:-1], KS_GPa[iObotPP+1:],     **intArgs)(zInner_m)
        QmuInner = interp1d(z_m[iObotPP+1:-1], QS[iObotPP+1:],            **intArgs)(zInner_m)

        # Interpolate layers within ocean
        rhoOcean_kgm3 = interp1d(z_m[Planet.Steps.nSurfIce:iObotPP+1], rho_kgm3[Planet.Steps.nSurfIce:iObotPP+1], **intArgs)(zOcean_m)
        VPocean_ms = interp1d(z_m[Planet.Steps.nSurfIce:iObotPP+1], VP_ms[Planet.Steps.nSurfIce:iObotPP+1],       **intArgs)(zOcean_m)
        VSocean_ms = interp1d(z_m[Planet.Steps.nSurfIce:iObotPP+1], VS_ms[Planet.Steps.nSurfIce:iObotPP+1],       **intArgs)(zOcean_m)
        QkappaOcean = interp1d(z_m[Planet.Steps.nSurfIce:iObotPP+1], KS_GPa[Planet.Steps.nSurfIce:iObotPP+1],     **intArgs)(zOcean_m)
        QmuOcean = interp1d(z_m[Planet.Steps.nSurfIce:iObotPP+1], QS[Planet.Steps.nSurfIce:iObotPP+1],            **intArgs)(zOcean_m)

        # Interpolate layers above ocean
        rhoShell_kgm3 = interp1d(z_m[:Planet.Steps.nSurfIce], rho_kgm3[:Planet.Steps.nSurfIce], **intArgs)(zShell_m)
        VPshell_ms = interp1d(z_m[:Planet.Steps.nSurfIce], VP_ms[:Planet.Steps.nSurfIce],       **intArgs)(zShell_m)
        VSshell_ms = interp1d(z_m[:Planet.Steps.nSurfIce], VS_ms[:Planet.Steps.nSurfIce],       **intArgs)(zShell_m)
        QkappaShell = interp1d(z_m[:Planet.Steps.nSurfIce], KS_GPa[:Planet.Steps.nSurfIce],     **intArgs)(zShell_m)
        QmuShell = interp1d(z_m[:Planet.Steps.nSurfIce], QS[:Planet.Steps.nSurfIce],            **intArgs)(zShell_m)

        # Finally, stitch together into output arrays
        rhominEOS_kgm3 = np.hstack((rhoInner_kgm3, rhoOcean_kgm3, rhoShell_kgm3))
        VPminEOS_ms = np.hstack((VPinner_ms, VPocean_ms, VPshell_ms))
        VSminEOS_ms = np.hstack((VSinner_ms, VSocean_ms, VSshell_ms))
        QkappaminEOS = np.hstack((QkappaInner, QkappaOcean, QkappaShell))
        QmuminEOS = np.hstack((QmuInner, QmuOcean, QmuShell))

        # Adjust ocean bottom and top indices since we added extra layers at boundaries
        iOceanBot += 1
        iOceanTop += 2

    # PlanetProfile calculates bulk modulus, not Qkappa; set values to 99999.0 to flag them
    QkappaminEOS[QkappaminEOS > -1] = 99999.0
    # Set max Qmu (QS) value to 99999.0
    QmuminEOS[np.logical_or(QmuminEOS > 99999, QmuminEOS < 0)] = 99999.0

    leadWSminEOS = 1
    with open(Params.DataFiles.minEOSvelFile,'w') as f:
        # Header
        f.write(f'velmodel\n' +
                f'   1    1.00000  1\n' +
                f'{nminEOS} {iOceanBot} {iOceanTop} 0\n')

        # Write data
        for i in range(nminEOS):
            f.write(' ' * leadWSminEOS + ' '.join([
                f'{rminEOS_m[i]:7.0f}',
                f'{rhominEOS_kgm3[i]:8.2f}',
                f'{VPminEOS_ms[i]:8.2f}',
                f'{VSminEOS_ms[i]:8.2f}',
                f'{QkappaminEOS[i]:8.1f}',
                f'{QmuminEOS[i]:8.1f}',
                f'{VPminEOS_ms[i]:8.2f}',
                f'{VSminEOS_ms[i]:8.2f}',
                f' 1.0000'
            ]) + '\n')

    """ @@@@@@@@@@@@@@@@@@
        minEOS Yannos file
        @@@@@@@@@@@@@@@@@@
    """
    # Write Yannos file
    yannosCfg = '# Input body model:\n' + \
                'velmodel\n' + \
                '# Output information file:\n' + \
                'per_model\n' + \
                '# Prefix for output eigenfunction files\n' + \
                'fct_model\n' + \
                '# Type code (3=spheroidal, 2=toroidal, 0=both)\n' + \
                f'{Planet.Seismic.minEOS_mode}\n' + \
                '# Precision (2) and switch off gravity perturbation frequency\n' + \
                f'{Planet.Seismic.minEOSyanPrec} {Planet.Seismic.minEOSyanPrec} {Planet.Seismic.fGravCutoff_mHz}\n' + \
                '# lmin, lmax, fmin, fmax, nmax:\n' + \
                f'{Planet.Seismic.minEOS_lmin} {Planet.Seismic.minEOS_lmax} {Planet.Seismic.minEOS_fmin_mHz} ' + \
                f'{Planet.Seismic.minEOS_fmax_mHz} {Planet.Seismic.minEOS_nmax}\n' + \
                '#### If you need to output only some radius layers of eigenfunction\n' + \
                '# number of layer for output: \n' + \
                '1\n' + \
                '# radius layers\n' + \
                f'00. {Planet.Bulk.R_m:.0f}.\n' + \
                '#### Computation flags\n' + \
                '# force fmin (usually F):\n' + \
                'F\n' + \
                '# Cancel Gravity (usually F)\n' + \
                'F\n' + \
                '# Never use start level (usually F):\n' + \
                'F\n' + \
                '# Use T ref (usually T):\n' + \
                'T\n' + \
                '# Check modes (usually F):\n' + \
                'F\n' + \
                '# use_remedy awkward modes (usually T):\n' + \
                'T\n' + \
                '# rescue (try to do something if missing a mode, usually T):\n' + \
                'T\n' + \
                '# restart (to restart in case of crash during computation)\n' + \
                'F\n' + \
                '# force_systematic_search\n' + \
                'F\n' + \
                '# keep_bad_modes\n' + \
                'F\n' + \
                '# modout_format (ipg, ucb, olm; ipg is standard)\n' + \
                'ipg\n' + \
                '# seuil_ray\n' + \
                '.0001\n' + \
                '# l_startlevel\n' + \
                '000'

    with open(Params.DataFiles.minEOSyanFile,'w') as f:
        f.write(yannosCfg)

    return
