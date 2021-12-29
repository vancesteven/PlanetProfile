import numpy as np
from Thermodynamics.FromLiterature.InnerEOS import TsolidusHirschmann2000
from Thermodynamics.FromLiterature.HydroEOS import GetPhaseIndices
from seafreeze import seafreeze as SeaFreeze
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
                = ClathrateSeismicHelgerud2009(Planet.P_MPa[indsClath], Planet.T_K[indsClath])

            Hclath = Planet.Seismic.gClath * np.max(Planet.T_K[indsClath])
            Planet.Seismic.QS[indsClath] = Planet.Seismic.BClath * np.exp(
                Planet.Seismic.gammaClath * Hclath / Planet.T_K[indsClath])

        # Get seismic properties of all ice layers
        if np.size(indsIceI) != 0:
            Planet.Seismic.VP_kms[indsIceI], Planet.Seismic.VS_kms[indsIceI], Planet.Seismic.KS_GPa[indsIceI], \
            Planet.Seismic.GS_GPa[indsIceI] = SeaFreezeSeismicIce(Planet.P_MPa[indsIceI], Planet.T_K[indsIceI], 1)
            HiceI = Planet.Seismic.gIceI * Planet.Bulk.Tb_K
            Planet.Seismic.QS[indsIceI] = Planet.Seismic.BIceI * np.exp(
                Planet.Seismic.gammaIceI * HiceI / Planet.T_K[indsIceI])
        if np.size(indsIceII) != 0:
            Planet.Seismic.VP_kms[indsIceII], Planet.Seismic.VS_kms[indsIceII], Planet.Seismic.KS_GPa[indsIceII], \
            Planet.Seismic.GS_GPa[indsIceII] = SeaFreezeSeismicIce(Planet.P_MPa[indsIceII], Planet.T_K[indsIceII], 2)
            HiceII = Planet.Seismic.gIceII * np.max(Planet.T_K[indsIceII])
            Planet.Seismic.QS[indsIceII] = Planet.Seismic.BIceII * np.exp(
                Planet.Seismic.gammaIceII * HiceII / Planet.T_K[indsIceII])
        if np.size(indsIceIII) != 0:
            Planet.Seismic.VP_kms[indsIceIII], Planet.Seismic.VS_kms[indsIceIII], Planet.Seismic.KS_GPa[indsIceIII], \
            Planet.Seismic.GS_GPa[indsIceIII] = SeaFreezeSeismicIce(Planet.P_MPa[indsIceIII], Planet.T_K[indsIceIII], 3)
            HiceIII = Planet.Seismic.gIceIII * np.max(Planet.T_K[indsIceIII])
            Planet.Seismic.QS[indsIceIII] = Planet.Seismic.BIceIII * np.exp(
                Planet.Seismic.gammaIceIII * HiceIII / Planet.T_K[indsIceIII])
        if np.size(indsIceV) != 0:
            Planet.Seismic.VP_kms[indsIceV], Planet.Seismic.VS_kms[indsIceV], Planet.Seismic.KS_GPa[indsIceV], \
            Planet.Seismic.GS_GPa[indsIceV] = SeaFreezeSeismicIce(Planet.P_MPa[indsIceV], Planet.T_K[indsIceV], 5)
            HiceV = Planet.Seismic.gIceV * np.max(Planet.T_K[indsIceV])
            Planet.Seismic.QS[indsIceV] = Planet.Seismic.BIceV * np.exp(
                Planet.Seismic.gammaIceV * HiceV / Planet.T_K[indsIceV])
        if np.size(indsIceVI) != 0:
            Planet.Seismic.VP_kms[indsIceVI], Planet.Seismic.VS_kms[indsIceVI], Planet.Seismic.KS_GPa[indsIceVI], \
            Planet.Seismic.GS_GPa[indsIceVI] = SeaFreezeSeismicIce(Planet.P_MPa[indsIceVI], Planet.T_K[indsIceVI], 6)
            HiceVI = Planet.Seismic.gIceVI * np.max(Planet.T_K[indsIceVI])
            Planet.Seismic.QS[indsIceVI] = Planet.Seismic.BIceVI * np.exp(
                Planet.Seismic.gammaIceVI * HiceVI / Planet.T_K[indsIceVI])
            
        # Get seismic properties of ocean layers
        if np.size(indsLiquid) != 0:
            Planet.Seismic.VP_kms[indsLiquid], Planet.Seismic.VS_kms[indsLiquid], Planet.Seismic.KS_GPa[indsLiquid], \
            Planet.Seismic.GS_GPa[indsLiquid] = SeaFreezeSeismicFluid(Planet.P_MPa[indsLiquid], Planet.T_K[indsLiquid],
                                                                    Planet.Ocean.comp, Planet.Ocean.wOcean_ppt)
            Planet.Seismic.QS[indsLiquid] = np.zeros(np.size(indsLiquid))

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


def ClathrateSeismicHelgerud2009(P_MPa, T_K):
    """ Calculate seismic velocities in clathrates based on Helgerud et al. (2009): https://doi.org/10.1029/2009JB006451
        Note that the original article had a correction for all of the tables' equations--
        the correction is linked in the above DOI.

        Args:
            P_MPa (float, shape N): Pressure values to evaluate in MPa.
            T_K (float, shape N): Temperature values to evaluate in K.
        Returns:
            VP_kms (float, shape N): P-wave seismic velocity in km/s.
            VS_kms (float, shape N): S-wave seismic velocity in km/s.
            KS_GPa (float, shape N): Bulk modulus in GPa.
            GS_GPa (float, shape N): Shear modulus in GPa.
    """

    T_C = T_K - Constants.T0
    VP_kms = (-1.84*T_C + 0.31*P_MPa + 3766) * 1e-3
    VS_kms = (-0.892*T_C - 0.1*P_MPa + 1957) * 1e-3
    KS_GPa = -1.09e-2*T_C + 3.8e-3*P_MPa + 8.39
    GS_GPa = -4.2e-3*T_C + 9e-5*P_MPa + 3.541

    return VP_kms, VS_kms, KS_GPa, GS_GPa


def SeaFreezeSeismicIce(P_MPa, T_K, phase):
    """ Calculate seismic properties in non-clathrate ices using SeaFreeze.
        Only calculates one phase at a time, because no advantage is offered
        in terms of the number of calls to SeaFreeze for packaging the data
        into larger chunks for this function.
        Note that all values returned from SeaFreeze for velocities are m/s and all moduli are
        in MPa, so everything must be scaled by a factor 1e-3.

        Args:
            P_MPa (float, shape N): Pressure values to evaluate in MPa.
            T_K (float, shape N): Temperature values to evaluate in K.
            phase (int): Phase of ice for all PT pairs.
        Returns:
            VP_kms (float, shape N): P-wave seismic velocity in km/s.
            VS_kms (float, shape N): S-wave seismic velocity in km/s.
            KS_GPa (float, shape N): Bulk modulus in GPa.
            GS_GPa (float, shape N): Shear modulus in GPa.
    """
    # Arrange input data into (P,T) value pair tuples compatible with SeaFreeze
    PTpairs = np.array([(P_MPa[i], T_K[i]) for i in range(np.size(P_MPa))], dtype='f,f').astype(object)

    if phase == 1:
        seaOut = SeaFreeze(PTpairs, 'Ih')
        VP_kms = seaOut.Vp * 1e-3
        VS_kms = seaOut.Vs * 1e-3
        KS_GPa = seaOut.Ks * 1e-3
        GS_GPa = seaOut.shear * 1e-3
    elif phase == 2:
        seaOut = SeaFreeze(PTpairs, 'II')
        VP_kms = seaOut.Vp * 1e-3
        VS_kms = seaOut.Vs * 1e-3
        KS_GPa = seaOut.Ks * 1e-3
        GS_GPa = seaOut.shear * 1e-3
    elif phase == 3:
        seaOut = SeaFreeze(PTpairs, 'III')
        VP_kms = seaOut.Vp * 1e-3
        VS_kms = seaOut.Vs * 1e-3
        KS_GPa = seaOut.Ks * 1e-3
        GS_GPa = seaOut.shear * 1e-3
    elif phase == 5:
        seaOut = SeaFreeze(PTpairs, 'V')
        VP_kms = seaOut.Vp * 1e-3
        VS_kms = seaOut.Vs * 1e-3
        KS_GPa = seaOut.Ks * 1e-3
        GS_GPa = seaOut.shear * 1e-3
    elif phase == 6:
        seaOut = SeaFreeze(PTpairs, 'VI')
        VP_kms = seaOut.Vp * 1e-3
        VS_kms = seaOut.Vs * 1e-3
        KS_GPa = seaOut.Ks * 1e-3
        GS_GPa = seaOut.shear * 1e-3
    else:
        raise ValueError('SeaFreezeSeismicIce was passed a phase that does not correspond to ice.')

    return VP_kms, VS_kms, KS_GPa, GS_GPa


def SeaFreezeSeismicFluid(P_MPa, T_K, compstr, w_ppt):
    """ Calculate seismic properties of ocean layers using SeaFreeze.
        Note that all values returned from SeaFreeze for velocities are m/s and all moduli are
        in MPa, so everything must be scaled by a factor 1e-3.

        Args:
            P_MPa (float, shape N): Pressure of the fluid in MPa
            T_K (float, shape N): Temperature of the fluid in K
            compstr (string): Composition of dissolved salt
            w_ppt (float): Salinity of fluid in ppt
        Returns:
            VP_kms (float, shape N): P-wave seismic velocity in km/s.
            VS_kms (float, shape N): S-wave seismic velocity in km/s.
            KS_GPa (float, shape N): Bulk modulus in GPa.
            GS_GPa (float, shape N): Shear modulus in GPa.
    """
    # Arrange input data into (P,T) value pair tuples compatible with SeaFreeze
    PTpairs = np.array([(P_MPa[i], T_K[i]) for i in range(np.size(P_MPa))], dtype='f,f').astype(object)

    if w_ppt == 0:
        seaOut = SeaFreeze(PTpairs, 'water1')
        VP_kms = seaOut.vel * 1e-3
        VS_kms = np.zeros_like(VP_kms)
        KS_GPa = seaOut.Ks * 1e-3
        GS_GPa = np.zeros_like(VP_kms)
    elif compstr == 'Seawater':
        raise ValueError('Unable to get seismic properties of ocean. Seawater is not implemented yet.')
    elif compstr == 'NH3':
        raise ValueError('Unable to get seismic properties of ocean. NH3 is not implemented yet.')
    elif compstr == 'MgSO4':
        raise ValueError('Unable to get seismic properties of ocean. MgSO4 is not implemented yet.')
    elif compstr == 'NaCl':
        raise ValueError('Unable to get seismic properties of ocean. NaCl is not implemented yet.')
    else:
        raise ValueError('Unable to get seismic properties of ocean. compstr="'+compstr+'" but options are Seawater, NH3, MgSO4, and NaCl.')

    return VP_kms, VS_kms, KS_GPa, GS_GPa