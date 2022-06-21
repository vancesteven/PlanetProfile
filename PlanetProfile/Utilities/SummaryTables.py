import numpy as np
import logging
from PlanetProfile.GetConfig import FigMisc, FigLbl
from PlanetProfile.Thermodynamics.HydroEOS import PhaseConv
from PlanetProfile.Utilities.defineStructs import Constants

# Assign logger
log = logging.getLogger('PlanetProfile')

def GetLayerMeans(PlanetList, Params):
    """ For calculating layer means we didn't need at any other point in our analysis,
        but that we might want to consider in comparing/analyzing profiles.
    """
    # Record if all models pertain to the same body
    Params.ALL_ONE_BODY = np.all([Planet.bodyname == PlanetList[0].bodyname for Planet in PlanetList])
    # Similar for waterless bodies
    Params.ALL_NO_H2O = np.all([Planet.Do.NO_H2O for Planet in PlanetList])

    # Get flags to pass on regarding types of layers we have
    Params.yesPorousRock = np.any([Planet.Do.POROUS_ROCK for Planet in PlanetList])
    Params.yesPorousIce = np.any([Planet.Do.POROUS_ICE for Planet in PlanetList])
    Params.yesClath = np.any([np.any(abs(Planet.phase) == Constants.phaseClath) for Planet in PlanetList])
    Params.yesIceIII = np.any([np.any(abs(Planet.phase) == 3) for Planet in PlanetList])
    Params.yesIceVund = np.any([np.any(Planet.phase == -5) for Planet in PlanetList])
    Params.yesIceV = np.any([np.any(Planet.phase == 5) for Planet in PlanetList])
    Params.yesIceVI = np.any([np.any(abs(Planet.phase) == 6) for Planet in PlanetList])
    Params.yesIceVandVI = np.any([np.any(Planet.phase == 5) or np.any(Planet.phase == 6) for Planet in PlanetList])

    for Planet in PlanetList:
        # Mean values are set to nan by default. Set relevant values here.

        if Planet.Do.VALID:
            # Inner layer means
            iSil = np.logical_and(Planet.phase >= Constants.phaseSil, Planet.phase < Constants.phaseSil + 10)
            if Planet.Do.POROUS_ROCK:
                Planet.dzSilPorous_km = (Planet.Sil.Rmean_m -
                    np.max(Planet.r_m[:-1][np.logical_and(iSil, Planet.phi_frac < Planet.Sil.phiMin_frac)],
                           initial=Planet.Core.Rmean_m)) / 1e3
            else:
                Planet.dzSilPorous_km = 0.0
            if Params.CALC_CONDUCT:
                Planet.Sil.sigmaMean_Sm = np.mean(Planet.sigma_Sm[iSil])
            # Get mean shear modulus in silicates
            if Params.CALC_SEISMIC:
                Planet.Sil.GSmean_GPa = np.mean(Planet.Seismic.GS_GPa[iSil])

            iFe = Planet.phase >= Constants.phaseFe
            if np.sum(iFe) > 0:
                iPureFe = np.logical_or(Planet.phase >= Constants.phaseFe, Planet.phase < Constants.phaseFeS)
                iFeS = np.logical_or(Planet.phase >= Constants.phaseFeS, Planet.phase < Constants.phaseFeS + 10)
                if np.sum(iPureFe) > 0:
                    Planet.Core.rhoMeanFe_kgm3 = np.mean(Planet.rho_kgm3[iPureFe])
                if np.sum(iFeS) > 0:
                    Planet.Core.rhoMeanFeS_kgm3 = np.mean(Planet.rho_kgm3[iFeS])
                    Planet.dzFeS_km = (np.max(Planet.r_m[:-1][iFeS]) - np.max(Planet.r_m[:-1][iPureFe], initial=0)) / 1e3
                else:
                    # We need this to be zero, not nan, for some wedge calculations
                    Planet.dzFeS_km = 0.0
                if Params.CALC_CONDUCT:
                    Planet.Core.sigmaMean_Sm = np.mean(Planet.sigma_Sm[iFe])
                # Get mean shear modulus in core
                if Params.CALC_SEISMIC:
                    Planet.Core.GSmean_GPa = np.mean(Planet.Seismic.GS_GPa[iFe])
                if np.sum(iPureFe) > 0:
                    Planet.Core.GSmeanFe_GPa = np.mean(Planet.Seismic.GS_GPa[iPureFe])
                if np.sum(iFeS) > 0:
                    Planet.Core.GSmeanFeS_GPa = np.mean(Planet.Seismic.GS_GPa[iFeS])

            # Hydrosphere layer means
            if not Planet.Do.NO_H2O:
                iCond = Planet.z_m[:-1] < Planet.eLid_m
                iConv = np.logical_and(Planet.z_m[:-1] >= Planet.eLid_m, Planet.z_m[:-1] < Planet.zb_km*1e3)
                iCondI = abs(Planet.phase[iCond]) == 1
                iCondClath = abs(Planet.phase[iCond]) == Constants.phaseClath
                iConvI = abs(Planet.phase[iConv]) == 1
                iConvClath = abs(Planet.phase[iConv]) == Constants.phaseClath
                if np.any(iCondI):
                    Planet.Ocean.rhoCondMean_kgm3['Ih'] = np.sum(Planet.MLayer_kg[iCond][iCondI]) / np.sum(Planet.VLayer_m3[iCond][iCondI])
                    # Get mean conductivity, ignoring spherical effects
                    if Params.CALC_CONDUCT:
                        Planet.Ocean.sigmaCondMean_Sm['Ih'] = np.mean(Planet.sigma_Sm[iCond][iCondI])
                    # Get mean shear modulus, ignoring spherical effects
                    if Params.CALC_SEISMIC:
                        Planet.Ocean.GScondMean_GPa['Ih'] = np.mean(Planet.Seismic.GS_GPa[iCond][iCondI])
                if np.any(iConvI):
                    Planet.Ocean.rhoConvMean_kgm3['Ih'] = np.sum(Planet.MLayer_kg[iConv][iConvI]) / np.sum(Planet.VLayer_m3[iConv][iConvI])
                    if Params.CALC_CONDUCT:
                        Planet.Ocean.sigmaConvMean_Sm['Ih'] = np.mean(Planet.sigma_Sm[iConv][iConvI])
                    if Params.CALC_SEISMIC:
                        Planet.Ocean.GSconvMean_GPa['Ih'] = np.mean(Planet.Seismic.GS_GPa[iConv][iConvI])
                if np.any(iCondClath):
                    Planet.Ocean.rhoCondMean_kgm3['Clath'] = np.sum(Planet.MLayer_kg[iCond][iCondClath]) / np.sum(Planet.VLayer_m3[iCond][iCondClath])
                    if Params.CALC_CONDUCT:
                        Planet.Ocean.sigmaCondMean_Sm['Clath'] = np.mean(Planet.sigma_Sm[iCond][iCondClath])
                    if Params.CALC_SEISMIC:
                        Planet.Ocean.GScondMean_GPa['Clath'] = np.mean(Planet.Seismic.GS_GPa[iCond][iCondClath])
                if np.any(iConvClath):
                    Planet.Ocean.rhoConvMean_kgm3['Clath'] = np.sum(Planet.MLayer_kg[iConv][iConvClath]) / np.sum(Planet.VLayer_m3[iConv][iConvClath])
                    if Params.CALC_CONDUCT:
                        Planet.Ocean.sigmaConvMean_Sm['Clath'] = np.mean(Planet.sigma_Sm[iConv][iConvClath])
                    if Params.CALC_SEISMIC:
                        Planet.Ocean.GSconvMean_GPa['Clath'] = np.mean(Planet.Seismic.GS_GPa[iConv][iConvClath])

                if Planet.Do.BOTTOM_ICEIII or Planet.Do.BOTTOM_ICEV:
                    iCondIII = np.logical_and(Planet.z_m[:-1] >= Planet.z_m[Planet.Steps.nIbottom],
                                                           Planet.z_m[:-1] < Planet.z_m[Planet.Steps.nIbottom] + Planet.eLidIII_m)
                    iConvIII = np.logical_and(Planet.z_m[:-1] >= Planet.z_m[Planet.Steps.nIbottom] + Planet.eLidIII_m,
                                                           Planet.z_m[:-1] < Planet.z_m[Planet.Steps.nIIIbottom])
                    Planet.Ocean.rhoCondMean_kgm3['III'] = np.sum(Planet.MLayer_kg[iCondIII]) / np.sum(Planet.VLayer_m3[iCondIII])
                    Planet.Ocean.rhoConvMean_kgm3['III'] = np.sum(Planet.MLayer_kg[iConvIII]) / np.sum(Planet.VLayer_m3[iConvIII])
                    if Params.CALC_CONDUCT:
                        Planet.Ocean.sigmaCondMean_Sm['III'] = np.mean(Planet.sigma_Sm[iCondIII])
                        if np.sum(iConvIII) > 0:
                            Planet.Ocean.sigmaConvMean_Sm['III'] = np.mean(Planet.sigma_Sm[iConvIII])
                    # Get mean shear moduli
                    if Params.CALC_SEISMIC:
                        Planet.Ocean.GScondMean_GPa['III'] = np.mean(Planet.Seismic.GS_GPa[iCondIII])
                        if np.sum(iConvIII) > 0:
                            Planet.Ocean.GSconvMean_GPa['III'] = np.mean(Planet.Seismic.GS_GPa[iConvIII])
                if Planet.Do.BOTTOM_ICEV:
                    iCondV = np.logical_and(Planet.z_m[:-1] >= Planet.z_m[Planet.Steps.nIIIbottom],
                                                         Planet.z_m[:-1] < Planet.z_m[Planet.Steps.nIIIbottom] + Planet.eLidV_m)
                    iConvV = np.logical_and(Planet.z_m[:-1] >= Planet.z_m[Planet.Steps.nIIIbottom] + Planet.eLidV_m,
                                                         Planet.z_m[:-1] < Planet.z_m[Planet.Steps.nSurfIce])
                    Planet.Ocean.rhoCondMean_kgm3['V'] = np.sum(Planet.MLayer_kg[iCondV]) / np.sum(Planet.VLayer_m3[iCondV])
                    Planet.Ocean.rhoConvMean_kgm3['V'] = np.sum(Planet.MLayer_kg[iConvV]) / np.sum(Planet.VLayer_m3[iConvV])
                    if Params.CALC_CONDUCT:
                        Planet.Ocean.sigmaCondMean_Sm['V'] = np.mean(Planet.sigma_Sm[iCondV])
                        if np.sum(iConvV) > 0:
                            Planet.Ocean.sigmaConvMean_Sm['V'] = np.mean(Planet.sigma_Sm[iConvV])
                    # Get mean shear moduli
                    if Params.CALC_SEISMIC:
                        Planet.Ocean.GScondMean_GPa['V'] = np.mean(Planet.Seismic.GS_GPa[iCondV])
                        if np.sum(iConvV) > 0:
                            Planet.Ocean.GSconvMean_GPa['V'] = np.mean(Planet.Seismic.GS_GPa[iConvV])

                # Non-underplate ice layer sizes
                if np.any(abs(Planet.phase) == 1):
                    Planet.zIceI_m = np.min(Planet.z_m[:-1][abs(Planet.phase) == 1])
                    Planet.dzIceI_km = (next(z_m for i, z_m in enumerate(Planet.z_m[:-1])
                                            if i > np.where(abs(Planet.phase) == 1)[0][0]
                                               and not abs(Planet.phase[i]) == 1) - Planet.zIceI_m) / 1e3
                else:
                    Planet.zIceI_m = np.nan
                    Planet.dzIceI_km = np.nan
                if np.any(abs(Planet.phase) == Constants.phaseClath):
                    # Note that this differs from Planet.zClath_m, which is used to set the thickness/depth of the BOTTOM
                    # of the clathrate lid in the "top" clathrate model.
                    Planet.zClath_km = np.min(Planet.z_m[:-1][abs(Planet.phase) == Constants.phaseClath])/1e3
                    Planet.dzClath_km = np.max(Planet.z_m[:-1][abs(Planet.phase) == Constants.phaseClath])/1e3 \
                                       - Planet.zClath_km
                else:
                    Planet.zClath_km = np.nan
                    Planet.dzClath_km = np.nan
                if np.any(abs(Planet.phase) == 3):
                    Planet.zIceIII_m = np.min(Planet.z_m[:-1][abs(Planet.phase) == 3])
                    Planet.dzIceIII_km = (Planet.z_m[Planet.Steps.nIIIbottom] - Planet.z_m[Planet.Steps.nIbottom])/1e3
                else:
                    Planet.zIceIII_m = np.nan
                    Planet.dzIceIII_km = np.nan
                if np.any(Planet.phase == -5):
                    Planet.zIceVund_m = np.min(Planet.z_m[:-1][Planet.phase == -5])
                    Planet.dzIceVund_km = Planet.zb_km - Planet.z_m[Planet.Steps.nIIIbottom]/1e3
                else:
                    Planet.zIceVund_m = np.nan
                    Planet.dzIceVund_km = np.nan
                if np.any(Planet.phase == 5):
                    Planet.zIceV_m = np.min(Planet.z_m[:-1][Planet.phase == 5])
                    Planet.dzIceV_km = (next(z_m for i, z_m in enumerate(Planet.z_m[:-1])
                                            if i > np.where(Planet.phase == 5)[0][0]
                                               and not Planet.phase[i] in [0, 5]) - Planet.zIceV_m) / 1e3
                    Planet.Ocean.rhoMeanVwet_kgm3 = np.mean(Planet.rho_kgm3[Planet.phase == 5])
                    if Params.CALC_CONDUCT:
                        Planet.Ocean.sigmaMeanVwet_Sm = np.mean(Planet.sigma_Sm[Planet.phase == 5])
                    if Params.CALC_SEISMIC:
                        Planet.Ocean.GSmeanVwet_GPa = np.mean(Planet.Seismic.GS_GPa[Planet.phase == 5])
                else:
                    Planet.zIceV_m = np.nan
                    Planet.dzIceV_km = np.nan
                if np.any(abs(Planet.phase) == 6):
                    Planet.zIceVI_m = np.min(Planet.z_m[:-1][abs(Planet.phase) == 6])
                    Planet.dzIceVI_km = (next(z_m for i, z_m in enumerate(Planet.z_m[:-1])
                                            if i > np.where(abs(Planet.phase) == 6)[0][0]
                                               and not Planet.phase[i] in [-6, 0, 6]) - Planet.zIceVI_m) / 1e3
                    Planet.Ocean.rhoMeanVI_kgm3 = np.mean(Planet.rho_kgm3[Planet.phase == 6])
                    if Params.CALC_CONDUCT:
                        Planet.Ocean.sigmaMeanVI_Sm = np.mean(Planet.sigma_Sm[Planet.phase == 6])
                    if Params.CALC_SEISMIC:
                        Planet.Ocean.GSmeanVI_GPa = np.mean(Planet.Seismic.GS_GPa[abs(Planet.phase) == 6])
                else:
                    Planet.zIceVI_m = np.nan
                    Planet.dzIceVI_km = np.nan
                if np.any(Planet.phase == 5) or np.any(Planet.phase == 6):
                    Planet.dzIceVandVI_km = (next(z_m for i, z_m in enumerate(Planet.z_m[:-1])
                                            if i > np.where(np.logical_or(Planet.phase == 5, Planet.phase == 6))[0][0]
                                               and not Planet.phase[i] in [0, 5, 6]) - Planet.zIceVI_m) / 1e3
                else:
                    Planet.dzIceVandVI_km = np.nan

    return PlanetList, Params


def PrintGeneralSummary(PlanetList, Params):
    """ Print out all of the bulk calculation outputs the user
        is likely to want for understanding the model and/or
        comparisons.
    """

    # Construct strings for printing #
    # Newline+indent constant for optional lines
    endl = '\n    '
    models = 'Models printed:\n' + '\n'.join([Planet.saveLabel for Planet in PlanetList])
    bodyMass = f'Body mass (kg): {PlanetList[0].Bulk.M_kg:.5e}'
    computedMass = 'Computed mass (kg): ' + ', '.join([f'{Planet.Mtot_kg:.5e}' for Planet in PlanetList])
    inputCMR2 = f'Input C/MR^2: {PlanetList[0].Bulk.Cmeasured} ± {PlanetList[0].Bulk.Cuncertainty}'
    if (np.any([Planet.Bulk.Cmeasured != PlanetList[0].Bulk.Cmeasured for Planet in PlanetList])
        or np.any([Planet.Bulk.Cuncertainty != PlanetList[0].Bulk.Cuncertainty for Planet in PlanetList])):
        log.warning('One or more moment of inertia parameters do not match the first profile. Only ' +
                    'the first profile\'s parameters will be printed.')
    belowCMR2 =    '            (-)  ' + ', '.join([f'{Planet.CMR2less:.5f}' for Planet in PlanetList])
    computedCMR2 = 'Computed C/MR^2: ' + ', '.join([f'{Planet.CMR2mean:.5f}' for Planet in PlanetList])
    aboveCMR2 =    '            (+)  ' + ', '.join([f'{Planet.CMR2more:.5f}' for Planet in PlanetList])
    zUpper = 'z_b (km): ' + ', '.join([f'{Planet.zb_km:.1f}' for Planet in PlanetList])
    # Note the below-surface-ice "wet hydrosphere" thickness separately from the liquid ocean thickness
    wetHydro = 'Total wet hydrosphere (km): ' + ', '.join([f'{(Planet.z_m[Planet.Steps.nHydro]/1e3 - Planet.zb_km):.1f}' for Planet in PlanetList])
    oceanThick = 'Ocean thickness D (km): ' + ', '.join([f'{Planet.D_km:.1f}' for Planet in PlanetList])
    oceanSigma = 'Mean ocean conductivity σ (S/m): ' + ', '.join([f'{Planet.Ocean.sigmaMean_Sm:.2f}' for Planet in PlanetList])
    oceanDensity = f'Mean ocean density (kg/m^3): ' + ', '.join([f'{Planet.Ocean.rhoMean_kgm3:.1f}' for Planet in PlanetList])
    zIceI = 'z(km) ice I: ' + ', '.join([f'{Planet.zIceI_m/1e3:.1f}' for Planet in PlanetList])

    # Optional output strings
    poreSigma, phiRockMax, phiIceMax, zClath, zIceIII, zIceVund, zIceV, zIceVI, dzClath, dzIceIII, dzIceVund, dzIceV, \
    dzIceVI, dzIceVandVI, RaIII, RaCritIII, eLidIII, DconvIII, deltaTBLIII, RaV, RaCritV, eLidV, DconvV, deltaTBLV \
        = ('' for _ in range(24))

    # Porosity
    if FigMisc.ALWAYS_SHOW_PHI or Params.yesPorousRock:
        poreSigma = f'{endl}Mean pore σ (S/m): ' + ', '.join([f'{Planet.Sil.sigmaPoreMean_Sm:.2f}' for Planet in PlanetList])
        phiRockMax = f'{endl}{endl}Max ϕsil (%): ' + ', '.join([f'{Planet.phi_frac[Planet.Steps.nHydro]*100:.1f}' for Planet in PlanetList])
    if FigMisc.ALWAYS_SHOW_PHI or Params.yesPorousIce:
        phiIceMax = f'{endl}Max ϕice (%): ' + ', '.join([f'{Planet.phi_frac[0]*100:.1f}' for Planet in PlanetList])

    # Temperatures and heat flux
    TbK = 'T_b (K): ' + ', '.join([f'{Planet.Bulk.Tb_K}' for Planet in PlanetList])
    TsilK = 'T_silTop (K): ' + ', '.join([f'{Planet.T_K[Planet.Steps.nHydro]:.1f}' for Planet in PlanetList])
    TCMBK = 'T_CMB (K): ' + ', '.join([f'{Planet.T_K[Planet.Steps.nHydro + Planet.Steps.nSil - 1]:.1f}' for Planet in PlanetList])
    TcenterK = 'T_center (K): ' + ', '.join([f'{Planet.T_K[Planet.Steps.nTotal - 1]:.1f}' for Planet in PlanetList])
    qOcean = 'qOcBot (mW/m^2): ' + ', '.join([f'{1e3*Planet.Ocean.QfromMantle_W/(4*np.pi*Planet.Sil.Rmean_m**2):.1f}' for Planet in PlanetList])
    qSurf = 'qSurf (mW/m^2): ' + ', '.join([f'{1e3*Planet.Ocean.QfromMantle_W/(4*np.pi*Planet.Bulk.R_m**2):.1f}' for Planet in PlanetList])

    # Only print HP ice values if they are present or forced on
    if FigMisc.ALWAYS_SHOW_HP or Params.yesClath:
        zClath = f'{endl}z(km) clath: ' + ', '.join([f'{Planet.zClath_km:.1f}' for Planet in PlanetList])
        dzClath = f'{endl}dz(km) clath: ' + ', '.join([f'{Planet.dzClath_km:.1f}' for Planet in PlanetList])
    if FigMisc.ALWAYS_SHOW_HP or Params.yesIceIII:
        zIceIII = f'{endl}z(km) ice III: ' + ', '.join([f'{Planet.zIceIII_m/1e3:.1f}' for Planet in PlanetList])
        dzIceIII = f'{endl}dz(km) ice III: ' + ', '.join([f'{Planet.dzIceIII_km:.1f}' for Planet in PlanetList])
    if FigMisc.ALWAYS_SHOW_HP or Params.yesIceVund:
        zIceVund = f'{endl}z(km) ice V (und): ' + ', '.join([f'{Planet.zIceVund_m/1e3:.1f}' for Planet in PlanetList])
        dzIceVund = f'{endl}dz(km) ice V (und): ' + ', '.join([f'{Planet.dzIceVund_km:.1f}' for Planet in PlanetList])
    if FigMisc.ALWAYS_SHOW_HP or Params.yesIceV:
        zIceV = f'{endl}z(km) ice V (wet): ' + ', '.join([f'{Planet.zIceV_m/1e3:.1f}' for Planet in PlanetList])
        dzIceV = f'{endl}dz(km) ice V (wet): ' + ', '.join([f'{Planet.dzIceV_km:.1f}' for Planet in PlanetList])
    if FigMisc.ALWAYS_SHOW_HP or Params.yesIceVI:
        zIceVI = f'{endl}z(km) ice VI: ' + ', '.join([f'{Planet.zIceVI_m/1e3:.1f}' for Planet in PlanetList])
        dzIceVI = f'{endl}dz(km) ice VI: ' + ', '.join([f'{Planet.dzIceVI_km:.1f}' for Planet in PlanetList])
    if FigMisc.ALWAYS_SHOW_HP or Params.yesIceVandVI:
        dzIceVandVI = f'{endl}dz(km) ice V (wet) + ice VI: ' + ', '.join([f'{Planet.dzIceVandVI_km:.1f}' for Planet in PlanetList])

    # Convection parameters for ice I/clathrate shell
    RaI = f'Ice shell Rayleigh number Ra: ' + ', '.join([f'{Planet.RaConvect:.2e}' for Planet in PlanetList])
    RaCritI = f'Critical Rayleigh number Ra_crit: ' + ', '.join([f'{Planet.RaCrit:.2e}' for Planet in PlanetList])
    eLidI = f'Conductive lid thickness e (km): ' + ', '.join([f'{Planet.eLid_m/1e3:.2f}' for Planet in PlanetList])
    DconvI = f'Convecting layer thickness D_conv (km): ' + ', '.join([f'{Planet.Dconv_m/1e3:.2f}' for Planet in PlanetList])
    deltaTBLI = f'Lower TBL thickness δ (km): ' + ', '.join([f'{Planet.deltaTBL_m/1e3:.2f}' for Planet in PlanetList])

    # Body mass fractions
    MfracH2O = f'Mass fraction H2O (%): ' + ', '.join([f'{100*Planet.MH2O_kg/Planet.Mtot_kg:.1f}' for Planet in PlanetList])
    MfracSalt = f'Mass fraction solutes (%): ' + ', '.join([f'{100*Planet.Msalt_kg/Planet.Mtot_kg:.1f}' for Planet in PlanetList])
    MfracSil = f'Mass fraction silicates (%): ' + ', '.join([f'{100*Planet.Mrock_kg/Planet.Mtot_kg:.1f}' for Planet in PlanetList])
    MfracFe = f'Mass fraction iron core (%): ' + ', '.join([f'{100*Planet.Mcore_kg/Planet.Mtot_kg:.1f}' for Planet in PlanetList])

    # Range of interior sizes/densities within ±1σ of the MoI match
    minRFe =  '               (-)  ' + ', '.join([f'{np.min(Planet.Core.Rtrade_m)/1e3:.1f}' for Planet in PlanetList])
    bestRFe = 'Best fit R_Fe (km): ' + ', '.join([f'{Planet.Core.Rmean_m/1e3:.1f}' for Planet in PlanetList])
    maxRFe =  '               (+)  ' + ', '.join([f'{np.max(Planet.Core.Rtrade_m)/1e3:.1f}' for Planet in PlanetList])
    minRsil =  '                (-)  ' + ', '.join([f'{np.min(Planet.Sil.Rtrade_m)/1e3:.1f}' for Planet in PlanetList])
    bestRsil = 'Best fit R_sil (km): ' + ', '.join([f'{Planet.Sil.Rmean_m/1e3:.1f}' for Planet in PlanetList])
    maxRsil =  '                (+)  ' + ', '.join([f'{np.max(Planet.Sil.Rtrade_m)/1e3:.1f}' for Planet in PlanetList])
    minrhoSil =  '                    (-)  ' + ', '.join([f'{np.min(Planet.Sil.rhoTrade_kgm3):.1f}' for Planet in PlanetList])
    bestrhoSil = 'Best fit ρ_sil (kg/m^3): ' + ', '.join([f'{Planet.Sil.rhoMean_kgm3:.1f}' for Planet in PlanetList])
    maxrhoSil =  '                    (+)  ' + ', '.join([f'{np.max(Planet.Sil.rhoTrade_kgm3):.1f}' for Planet in PlanetList])

    log.info(f"""
    {models}
    {bodyMass}
    {computedMass}
    {inputCMR2}
    {belowCMR2}
    {computedCMR2}
    {aboveCMR2}
    {zUpper}
    {oceanThick}
    {oceanDensity}
    {TbK}
    {TsilK}
    {TCMBK}
    {TcenterK}
    {qOcean}
    {qSurf}
    {oceanSigma}{poreSigma}{phiRockMax}{phiIceMax}
    
    {zIceI}{zClath}{zIceIII}{zIceVund}{zIceV}{zIceVI}{dzClath}{dzIceIII}{dzIceVund}{dzIceV}{dzIceVI}{dzIceVandVI}
    {wetHydro}
    {RaI}
    {RaCritI}
    {eLidI}
    {DconvI}
    {deltaTBLI}
    
    {MfracH2O}
    {MfracSalt}
    {MfracSil}
    {MfracFe}
    
    {minRFe}
    {bestRFe}
    {maxRFe}
    {minRsil}
    {bestRsil}
    {maxRsil}
    {minrhoSil}
    {bestrhoSil}
    {maxrhoSil}
    """)
    return


def PrintLayerSummaryLatex(PlanetList, Params):
    """ Construct multiple latex tables for summarizing layer structures
        for each model.
    """

    # Table horizontal division constant
    tab = ' & '
    hline = r'\hline'
    # Table vertical division constant
    endl = r' \\ '
    if FigMisc.LATEX_HLINES:
        endl = endl + hline
    newline = '\n            '
    zero = r'\num{0.0}'
    # Vertical lines for table, if present
    if FigMisc.LATEX_VLINES:
        v = ' | '
    else:
        v = ' '
    # Table end
    tClose = r'\end{tabular}'
    if FigMisc.HF_HLINES and not FigMisc.LATEX_HLINES:
        tClose = hline + tClose

    # Layer table labels
    columns = [r'\textbf{Layer}', r'\textbf{Radius ($\si{km}$)}', r'\textbf{Density ($\si{' + FigLbl.rhoUnits + '}$)}',
               r'\textbf{Thickness ($\si{km}$)}', r'\textbf{Shear modulus ($\si{GPa}$)}',
                r'\textbf{Conductivity ($\si{' + FigLbl.sigUnits + '}$)}']
    header = '\n' + tab.join(columns) + endl
    if FigMisc.HF_HLINES:
        header = hline + header
    tOpen = r'\begin{tabular}{' + v + v.join(['l' for _ in columns]) + v + '}'
    condIceLbl = 'Conductive ice Ih'
    convIceLbl = 'Convective ice Ih'
    condIceIIIlbl = 'Conductive ice III'
    convIceIIIlbl = 'Convective ice III'
    condIceVlbl = 'Conductive ice V'
    convIceVlbl = 'Convective ice V'
    condClathLbl = 'Conductive \ce{CH4} clathrates'
    convClathLbl = 'Convective \ce{CH4} clathrates'
    wetIceVlbl = 'Ice V'
    iceVIlbl = 'Ice VI'
    oceanLbl = 'ocean'
    silLbl = 'Mantle'
    coreLbl = 'Core'
    notPresent = [FigLbl.NA, FigLbl.NA, zero, FigLbl.NA, FigLbl.NA]
    emptyCondIce = newline + tab.join(np.append(condIceLbl, notPresent)) + endl
    emptyConvIce = newline + tab.join(np.append(convIceLbl, notPresent)) + endl
    emptyConvClath = newline + tab.join(np.append(convClathLbl, notPresent)) + endl

    log.info(FigMisc.latexPreamble)

    log.info('Layer tables:')
    for Planet in PlanetList:
        title = f'\section*{{{Planet.name}}}'
        if Planet.Do.NO_H2O:
            surfIceLayers = ''
            HPiceLayers = ''
        else:
            if Planet.Do.CLATHRATE:
                if Planet.Bulk.clathType == 'top' or Planet.Bulk.clathType == 'whole':
                    condClathLayers = newline + tab.join([condClathLbl,
                                                          f'\\num{{{Planet.Bulk.R_m/1e3:.1f}}}',
                                                          f'\\num{{{Planet.Ocean.rhoCondMean_kgm3["Clath"]:.0f}}}',
                                                          f'\\num{{{np.minimum(Planet.dzClath_km, Planet.eLid_m/1e3):.1f}}}',
                                                          f'\\num{{{Planet.Ocean.GScondMean_GPa["Clath"]:.1f}}}',
                                                          f'\\num{{{Planet.Ocean.sigmaCondMean_Sm["Clath"]:.1e}}}']) + endl
                    # For "top" clathrate model, clathrates are limited to the conductive lid
                    # of the ice shell, so if there are any convecting layers they will not be
                    # clathrates.
                    if Planet.Bulk.clathType == 'top':
                        convClathLayers = emptyConvClath
                        if Planet.zIceI_m < Planet.eLid_m:
                            condIceLayers = newline + tab.join([condIceLbl,
                                                                f'\\num{{{(Planet.Bulk.R_m - Planet.zIceI_m)/1e3:.1f}}}',
                                                                f'\\num{{{Planet.Ocean.rhoCondMean_kgm3["Ih"]:.0f}}}',
                                                                f'\\num{{{Planet.dzIceI_km - (Planet.Dconv_m + Planet.deltaTBL_m)/1e3}}}',
                                                                f'\\num{{{Planet.Ocean.GScondMean_GPa["Ih"]:.1f}}}',
                                                                f'\\num{{{Planet.Ocean.sigmaCondMean_Sm["Ih"]:.1e}}}']) + endl
                        else:
                            condIceLayers = emptyCondIce
                        if Planet.Dconv_m > 0:
                            convIceLayers = newline + tab.join([convIceLbl,
                                                                f'\\num{{{(Planet.Bulk.R_m - Planet.eLid_m)/1e3:.1f}}}',
                                                                f'\\num{{{Planet.Ocean.rhoConvMean_kgm3["Ih"]:.0f}}}',
                                                                f'\\num{{{(Planet.Dconv_m + Planet.deltaTBL_m)/1e3:.1f}}}',
                                                                f'\\num{{{Planet.Ocean.GSconvMean_GPa["Ih"]:.1f}}}',
                                                                f'\\num{{{Planet.Ocean.sigmaConvMean_Sm["Ih"]:.1e}}}']) + endl
                        else:
                            convIceLayers = emptyConvIce
                        surfIceLayers = condClathLayers + convClathLayers + condIceLayers + convIceLayers
                    else:
                        convClathLayers = newline + tab.join([convClathLbl,
                                                              f'\\num{{{(Planet.Bulk.R_m - Planet.eLid_m)/1e3:.1f}}}',
                                                              f'\\num{{{Planet.Ocean.rhoConvMean_kgm3["Clath"]:.0f}}}',
                                                              f'\\num{{{(Planet.Dconv_m + Planet.deltaTBL_m)/1e3:.1f}}}',
                                                              f'\\num{{{Planet.Ocean.GSconvMean_GPa["Clath"]:.1f}}}',
                                                              f'\\num{{{Planet.Ocean.sigmaCondMean_Sm["Clath"]:.1e}}}']) + endl
                        surfIceLayers = condClathLayers + convClathLayers
    
                else:
                    condIceLayers = newline + tab.join([condIceLbl,
                                                        f'\\num{{{Planet.Bulk.R_m/1e3:.1f}}}',
                                                        f'\\num{{{Planet.Ocean.rhoCondMean_kgm3["Ih"]:.0f}}}',
                                                        f'\\num{{{Planet.dzIceI_km:.1f}}}',
                                                        f'\\num{{{Planet.Ocean.GScondMean_GPa["Ih"]:.1f}}}',
                                                        f'\\num{{{Planet.Ocean.sigmaCondMean_Sm["Ih"]:.1e}}}']) + endl
                    # For "bottom" clathrate model, whole shell is assumed conductive because of the expected
                    # shallow temperature gradient due to the very strong insulating effect of the clathrate layer.
                    # This results in an assumption of no convection in the ice shell because of the small Rayleigh
                    # number.
                    convIceLayers = emptyConvIce
                    condClathLayers = newline + tab.join([condClathLbl,
                                                          f'\\num{{{Planet.Bulk.R_m - Planet.zClath_km:.1f}}}',
                                                          f'\\num{{{Planet.Ocean.rhoCondMean_kgm3["Clath"]:.0f}}}',
                                                          f'\\num{{{Planet.dzClath_km:.1f}}}',
                                                          f'\\num{{{Planet.Ocean.GScondMean_GPa["Clath"]:.1f}}}',
                                                          f'\\num{{{Planet.Ocean.sigmaCondMean_Sm["Clath"]:.1e}}}']) + endl
                    convClathLayers = emptyConvClath
                    surfIceLayers = condIceLayers + convIceLayers + condClathLayers + convClathLayers
            else:
                # No clathrates in upper ice
                if Planet.eLid_m > 0:
                    condIceLayers = newline + tab.join([condIceLbl,
                                                        f'\\num{{{Planet.Bulk.R_m/1e3:.1f}}}',
                                                        f'\\num{{{Planet.Ocean.rhoCondMean_kgm3["Ih"]:.0f}}}',
                                                        f'\\num{{{Planet.eLid_m/1e3:.1f}}}',
                                                        f'\\num{{{Planet.Ocean.GScondMean_GPa["Ih"]:.1f}}}',
                                                        f'\\num{{{Planet.Ocean.sigmaCondMean_Sm["Ih"]:.1e}}}']) + endl
                else:
                    condIceLayers = emptyCondIce
                if Planet.Dconv_m > 0:
                    convIceLayers = newline + tab.join([convIceLbl,
                                                        f'\\num{{{(Planet.Bulk.R_m - Planet.eLid_m)/1e3:.1f}}}',
                                                        f'\\num{{{Planet.Ocean.rhoConvMean_kgm3["Ih"]:.0f}}}',
                                                        f'\\num{{{(Planet.Dconv_m + Planet.deltaTBL_m)/1e3:.1f}}}',
                                                        f'\\num{{{Planet.Ocean.GSconvMean_GPa["Ih"]:.1f}}}',
                                                        f'\\num{{{Planet.Ocean.sigmaConvMean_Sm["Ih"]:.1e}}}']) + endl
                else:
                    convIceLayers = emptyConvIce
                surfIceLayers = condIceLayers + convIceLayers
                    
            # Underplating HP ices
            if Planet.Do.BOTTOM_ICEIII or Planet.Do.BOTTOM_ICEV:
                if Planet.eLidIII_m > 0:
                    condIceIIIlayers = newline + tab.join([condIceIIIlbl,
                                                           f'\\num{{{(Planet.Bulk.R_m - Planet.zIceIII_m)/1e3:.1f}}}',
                                                           f'\\num{{{Planet.Ocean.rhoCondMean_kgm3["III"]:.0f}}}',
                                                           f'\\num{{{Planet.eLidIII_m/1e3:.1f}}}',
                                                           f'\\num{{{Planet.Ocean.GScondMean_GPa["III"]:.1f}}}',
                                                           f'\\num{{{Planet.Ocean.sigmaCondMean_Sm["III"]:.1e}}}']) + endl
                else:
                    condIceIIIlayers = newline + tab.join(np.append(condIceIIIlbl, notPresent)) + endl
                if Planet.DconvIII_m > 0:
                    convIceIIIlayers = newline + tab.join([convIceIIIlbl, 
                                                           f'\\num{{{(Planet.Bulk.R_m - Planet.zIceIII_m - Planet.eLidIII_m)/1e3:.1f}}}',
                                                           f'\\num{{{Planet.Ocean.rhoConvMean_kgm3["III"]:.0f}}}',
                                                           f'\\num{{{(Planet.DconvIII_m + Planet.deltaTBLIII_m)/1e3:.1f}}}',
                                                           f'\\num{{{Planet.Ocean.GSconvMean_GPa["III"]:.1f}}}',
                                                           f'\\num{{{Planet.Ocean.sigmaConvMean_Sm["III"]:.1e}}}']) + endl
                else:
                    convIceIIIlayers = newline + tab.join(np.append(convIceIIIlbl, notPresent)) + endl
                surfIceLayers = surfIceLayers + condIceIIIlayers + convIceIIIlayers
            if Planet.Do.BOTTOM_ICEV:                
                if Planet.eLidV_m > 0:
                    condIceVlayers = newline + tab.join([condIceVlbl,
                                                         f'\\num{{{(Planet.Bulk.R_m - Planet.zIceVund_m)/1e3:.1f}}}',
                                                         f'\\num{{{Planet.Ocean.rhoCondMean_kgm3["V"]:.0f}}}',
                                                         f'\\num{{{Planet.eLidV_m/1e3:.1f}}}',
                                                         f'\\num{{{Planet.Ocean.GScondMean_GPa["V"]:.1f}}}',
                                                         f'\\num{{{Planet.Ocean.sigmaCondMean_Sm["V"]:.1e}}}']) + endl
                else:
                    condIceVlayers = newline + tab.join(np.append(condIceVlbl, notPresent)) + endl
                if Planet.DconvV_m > 0:
                    convIceVlayers = newline + tab.join([convIceVlbl, 
                                                         f'\\num{{{(Planet.Bulk.R_m - Planet.zIceVund_m - Planet.eLidV_m)/1e3:.1f}}}',
                                                         f'\\num{{{Planet.Ocean.rhoConvMean_kgm3["V"]:.0f}}}',
                                                         f'\\num{{{(Planet.DconvV_m + Planet.deltaTBLV_m)/1e3:.1f}}}',
                                                         f'\\num{{{Planet.Ocean.GSconvMean_GPa["V"]:.1f}}}',
                                                         f'\\num{{{Planet.Ocean.sigmaConvMean_Sm["V"]:.1e}}}']) + endl
                else:
                    convIceVlayers = newline + tab.join(np.append(convIceVlbl, notPresent)) + endl
                surfIceLayers = surfIceLayers + condIceVlayers + convIceVlayers

            # In-ocean HP ices
            HPiceLayers = ''
            if np.any(Planet.phase == 5):
                wetIceVlayers = newline + tab.join([wetIceVlbl,
                                                    f'\\num{{{(Planet.Bulk.R_m - Planet.zIceV_m)/1e3:.1f}}}',
                                                    f'\\num{{{Planet.Ocean.rhoMeanVwet_kgm3:.0f}}}',
                                                    f'\\num{{{Planet.dzIceV_km:.1f}}}',
                                                    f'\\num{{{Planet.Ocean.GSmeanVwet_GPa:.1f}}}',
                                                    f'\\num{{{Planet.Ocean.sigmaMeanVwet_Sm:.1e}}}']) + endl
                HPiceLayers = HPiceLayers + wetIceVlayers
            if np.any(abs(Planet.phase) == 6):
                iceVIlayers = newline + tab.join([iceVIlbl,
                                                  f'\\num{{{(Planet.Bulk.R_m - Planet.zIceVI_m)/1e3:.1f}}}',
                                                  f'\\num{{{Planet.Ocean.rhoMeanVI_kgm3:.0f}}}',
                                                  f'\\num{{{Planet.dzIceVI_km:.1f}}}',
                                                  f'\\num{{{Planet.Ocean.GSmeanVI_GPa:.1f}}}',
                                                  f'\\num{{{Planet.Ocean.sigmaMeanVI_Sm:.1e}}}']) + endl
                HPiceLayers = HPiceLayers + iceVIlayers
                
        # Ocean layers
        if Planet.Do.NO_H2O:
            salt = 'No'
            rhoOcean = FigLbl.NA
            GSocean = FigLbl.NA
            sigmaOcean = FigLbl.NA
        else:
            rhoOcean = f'\\num{{{Planet.Ocean.rhoMean_kgm3:.0f}}}'
            GSocean = zero
            if Planet.Ocean.sigmaMean_Sm > 0.01:
                sigmaOcean = f'\\num{{{Planet.Ocean.sigmaMean_Sm:.2f}}}'
            else:
                sigmaOcean = f'\\num{{{Planet.Ocean.sigmaMean_Sm:.2e}}}'
            if Planet.Ocean.wOcean_ppt == 0:
                salt = r'Pure~\ce{H2O}'
            else:
                salt = f'$\SI{{{Planet.Ocean.wOcean_ppt*FigLbl.wMult:.1f}}}{{{FigLbl.wUnits}\,\ce{{{Planet.Ocean.comp}}}}}$'

        oceanLayers = tab.join([f'{salt} {oceanLbl}',
                                f'\\num{{{Planet.Bulk.R_m/1e3 - Planet.zb_km:.1f}}}',
                                rhoOcean,
                                f'\\num{{{Planet.D_km:.1f}}}',
                                GSocean,
                                sigmaOcean]) + endl
        silLayers = tab.join([silLbl,
                              f'\\num{{{Planet.Sil.Rmean_m/1e3:.1f}}}',
                              f'\\num{{{Planet.Sil.rhoMean_kgm3:.0f}}}',
                              f'\\num{{{(Planet.Sil.Rmean_m - Planet.Core.Rmean_m)/1e3:.1f}}}',
                              f'\\num{{{Planet.Sil.GSmean_GPa:.1f}}}',
                              f'\\num{{{Planet.Sil.sigmaMean_Sm:.1e}}}']) + endl
        if Planet.Do.Fe_CORE:
            coreLayers = newline + tab.join([coreLbl,
                                   f'\\num{{{Planet.Core.Rmean_m/1e3:.1f}}}',
                                   f'\\num{{{Planet.Core.rhoMean_kgm3:.0f}}}',
                                   f'\\num{{{Planet.Core.Rmean_m/1e3:.1f}}}',
                                   f'\\num{{{Planet.Core.GSmean_GPa:.1f}}}',
                                   f'\\num{{{Planet.Core.sigmaMean_Sm:.1e}}}']) + endl
        else:
            coreLayers = ''
        log.info(f"""{Planet.saveLabel}
        {title}
        {tOpen}
        {header}{surfIceLayers}
            {oceanLayers}{HPiceLayers}
            {silLayers}{coreLayers}
        {tClose}
        """)

    return


def PrintLayerTableLatex(PlanetList, Params):
    """ Construct a Latex table for collating comparisons between models.
    """

    # Table horizontal division constant
    tab = ' & '
    # Table vertical division constant
    if FigMisc.LATEX_HLINES:
        endl = r' \\ \hline'
    else:
        endl = r' \\'
    newline = '\n            '
    # Vertical lines for table, if present
    if FigMisc.LATEX_VLINES:
        v = ' | '
    else:
        v = ' '
    
    # Subscript strings
    ice = r'\mathrm{ice}'
    Ih = r'\mathrm{Ih}'
    III = r'\mathrm{III}'
    Vund = r'\mathrm{V,under}'
    Vwet = r'\mathrm{V}'
    VI = r'\mathrm{VI}'
    clath = r'\mathrm{clath}'
    mixed = r'\mathrm{Ih+clath}'
    ocean = r'\mathrm{ocean}'
    rock = r'\mathrm{rock}'
    rockMean = r'\mathrm{rock,mean}'
    model = r'\mathrm{model}'
    surf = r'\mathrm{surf}'
    con = r'\mathrm{con}'
    core = r'\mathrm{core}'

    # Header line
    if FigMisc.HF_HLINES:
        header = '\hline'
    else:
        header = r'\ '
    # Table end
    tClose = header + r'\end{tabular}'

    nModels = np.size(PlanetList)
    boolDclath_km, boolDIII_km, boolDVund_km, boolDVwet_km, boolDVI_km, boolRcore_km, \
    boolTop, boolBottom, boolWhole, boolphiIceMax_frac, boolphiRockMax_frac \
        = (np.zeros(nModels, dtype=np.bool_) for _ in range(11))

    strMmeas_kg, strMcalc_kg, strCmeas, strCcalc, strRsurf_km, strrhoRock_kgm3, strTb_K, \
    strqSurf_Wm2, strqCon_Wm2, stretaI_Pas, strDIh_km, strDclath_km, strDIII_km, \
    strDVund_km, strD_km, strDVwet_km, strDVI_km, strsigOcean_Sm, strRrock_km, strRcore_km,\
    strphiIceMax_frac, strphiRockMax_frac \
        = (np.full(nModels, FigLbl.NA, dtype='<U100') for _ in range(22))

    # Organize values into strings for tabulating
    for i, Planet in enumerate(PlanetList):
        Cp = Planet.CMR2more - Planet.CMR2mean
        Cm = Planet.CMR2mean - Planet.CMR2less
        if not np.isnan(Cp) or not np.isnan(Cm):
            CMR2pm = f'\substack{{+{Cp:.5f} \\\\ -{Cm:.5f}}}'
        else:
            CMR2pm = ''

        strMmeas_kg[i] = f'$\\num{{{Planet.Bulk.M_kg:.4e}}}$'
        strMcalc_kg[i] = f'$\\num{{{Planet.Mtot_kg:.4e}}}$'
        strCmeas[i] = f'${Planet.Bulk.Cmeasured:.5f}\pm{Planet.Bulk.Cuncertainty:.5f}$'
        strCcalc[i] = f'${Planet.CMR2mean:.5f}{CMR2pm}$'
        strRsurf_km[i] = f'$\\num{{{Planet.Bulk.R_m/1e3:.1f}}}$'
        strrhoRock_kgm3[i] = f'$\\num{{{Planet.Sil.rhoMean_kgm3:.0f}}}$'
        strTb_K[i] = f'$\\num{{{Planet.Bulk.Tb_K}}}$'
        strqSurf_Wm2[i] = f'$\\num{{{Planet.qSurf_Wm2*FigLbl.qMult:.1f}}}$'
        if not Planet.Do.NO_H2O:
            strqCon_Wm2[i] = f'$\\num{{{Planet.qCon_Wm2*FigLbl.qMult:.1f}}}$'
            stretaI_Pas[i] = f'$\\num{{{Planet.etaConv_Pas:.2e}}}$'
            strDIh_km[i] = f'$\\num{{{Planet.dzIceI_km:.1f}}}$'
            strD_km[i] = f'$\\num{{{Planet.D_km:.1f}}}$'
            strsigOcean_Sm[i] = f'$\\num{{{Planet.Ocean.sigmaMean_Sm:.1f}}}$'
            if Planet.Do.CLATHRATE:
                strDclath_km[i] = f'$\\num{{{Planet.dzClath_km:.1f}}}$'
                if Planet.Bulk.clathType == 'top':
                    boolTop[i] = True
                elif Planet.Bulk.clathType == 'bottom':
                    boolBottom[i] = True
                elif Planet.Bulk.clathType == 'whole':
                    boolWhole[i] = True
                else:
                    raise ValueError(f'Bulk.clathType "{Planet.Bulk.clathType}" not recognized.')

            # HP ices
            if not np.isnan(Planet.dzIceIII_km) and not round(Planet.dzIceIII_km) == 0:
                strDIII_km[i] = f'$\\num{{{Planet.dzIceIII_km:.1f}}}$'
                boolDIII_km[i] = True
            if not np.isnan(Planet.dzIceVund_km) and not round(Planet.dzIceVund_km) == 0:
                strDVund_km[i] = f'$\\num{{{Planet.dzIceVund_km:.1f}}}$'
                boolDVund_km[i] = True
            if not np.isnan(Planet.dzIceV_km) and not round(Planet.dzIceV_km) == 0:
                strDVwet_km[i] = f'$\\num{{{Planet.dzIceV_km:.1f}}}$'
                boolDVwet_km[i] = True
            if not np.isnan(Planet.dzIceVI_km) and not round(Planet.dzIceVI_km) == 0:
                strDVI_km[i] = f'$\\num{{{Planet.dzIceVI_km:.1f}}}$'
                boolDVI_km[i] = True

        strRrock_km[i] = f'$\\num{{{Planet.Sil.Rmean_m/1e3:.1f}}}$'
        if Planet.Do.Fe_CORE:
            strRcore_km[i] = f'$\\num{{{Planet.Core.Rmean_m/1e3:.1f}}}$'
            boolRcore_km[i] = True

        if Planet.Do.POROUS_ICE:
            strphiIceMax_frac[i] = f'$\\num{{{Planet.Ocean.phiMax_frac["Ih"]*FigLbl.phiMult:.2f}}}$'
            boolphiIceMax_frac[i] = True
        if Planet.Do.POROUS_ROCK:
            strphiRockMax_frac[i] = f'$\\num{{{Planet.Sil.phiRockMax_frac*FigLbl.phiMult:.2f}}}$'
            boolphiRockMax_frac[i] = True

    # Begin constructing and printing tables
    log.info(FigMisc.latexPreamble)
    log.info('Comparison tables:')

    for thisComp in np.unique([Planet.Ocean.comp for Planet in PlanetList]):
        for thisw in np.unique([Planet.Ocean.wOcean_ppt for Planet in PlanetList if Planet.Ocean.comp == thisComp]):
            thisSubset = [Planet.Ocean.comp == thisComp and Planet.Ocean.wOcean_ppt == thisw
                          for Planet in PlanetList]
            # Table begin
            tOpen = r'\begin{tabular}{' + v + v.join(['c', 'c'] + ['c' for _ in PlanetList[thisSubset]]) + v + '}'
            if thisComp == 'PureH2O':
                compStr = r'Pure~\ce{H2O}'
                wStr = ''
                Tsub = 'b'
            elif thisComp == 'none':
                compStr = r'No~\ce{H2O}'
                wStr = ''
                Tsub = surf
            else:
                compStr = f'\ce{{{thisComp}}}'
                wStr = f'$\SI{{{thisw*FigLbl.wMult:.1f}}}{{{FigLbl.wUnits}}}$'
                Tsub = 'b'

            Tb_K = f'{tab}$T_{Tsub}~(\si{{K}})${tab}' + tab.join(strTb_K[thisSubset]) + endl
            rhoRock = f'{tab}$\\rho_{rockMean}~(\si{{{FigLbl.rhoUnits}}})${tab}' + tab.join(strrhoRock_kgm3[thisSubset]) + endl
            if FigMisc.PRINT_BULK:
                Mmeas = newline + f'{tab}$M~(\si{{kg}})${tab}' + tab.join(strMmeas_kg[thisSubset]) + endl
                Mcalc = newline + f'{wStr}{tab}$M_{model}~(\si{{kg}})${tab}' + tab.join(strMcalc_kg[thisSubset]) + endl
                Cmeas = newline + f'{tab}$C/MR^2${tab}' + tab.join(strCmeas[thisSubset]) + endl
                Ccalc = newline + f'{tab}$C_{model}/MR^2${tab}' + tab.join(strCcalc[thisSubset]) + endl
                rhoRock = newline + rhoRock
                Rsurf = newline + f'{tab}$R_{surf}~(\si{{km}})${tab}' + tab.join(strRsurf_km[thisSubset]) + endl
            else:
                Mmeas, Mcalc, Cmeas, Ccalc, Rsurf = ('' for _ in range(5))
                Tb_K = wStr + Tb_K

            qSurf = f'{tab}$q_{surf}~(\si{{{FigLbl.fluxUnits}}})${tab}' + tab.join(strqSurf_Wm2[thisSubset]) + endl
            if thisComp != 'none':
                qCon = newline + f'{tab}$q_{con}~(\si{{{FigLbl.fluxUnits}}})${tab}' + tab.join(strqCon_Wm2[thisSubset]) + endl
                etaI = newline + f'{tab}$\eta_{con}~(\si{{Pa\,s}})${tab}' + tab.join(stretaI_Pas[thisSubset]) + endl
                Docean = newline + f'{tab}$D_{ocean}~(\si{{km}})${tab}' + tab.join(strD_km[thisSubset]) + endl
                sigOcean = newline + f'{tab}$\overline{{\sigma}}_{ocean}~(\si{{{FigLbl.sigUnits}}})${tab}' + tab.join(strsigOcean_Sm[thisSubset]) + endl

                # Surface ices
                if np.any(np.logical_or(boolTop[thisSubset], boolWhole[thisSubset], boolWhole[thisSubset])):
                    if np.all(boolTop[thisSubset]):
                        DI = newline + f'{tab}$D_{clath}~(\si{{km}})${tab}' + tab.join(strDclath_km[thisSubset]) + endl \
                            + newline + f'{tab}$D_{Ih}~(\si{{km}})${tab}' + tab.join(strDIh_km[thisSubset]) + endl
                    elif np.all(boolBottom[thisSubset]):
                        DI = newline + f'{tab}$D_{Ih}~(\si{{km}})${tab}' + tab.join(strDIh_km[thisSubset]) + endl \
                            + newline + f'{tab}$D_{clath}~(\si{{km}})${tab}' + tab.join(strDclath_km[thisSubset]) + endl
                    elif np.all(boolWhole[thisSubset]):
                        DI = newline + f'{tab}$D_{clath}~(\si{{km}})${tab}' + tab.join(strDclath_km[thisSubset]) + endl
                    else:
                        DI = newline + f'{tab}$D_{mixed}~(\si{{km}})${tab}' + tab.join([
                            f'$\\num{{{Planet.dzIceI_km + Planet.dzClath_km:.1f}}}$'
                            for Planet in PlanetList[thisSubset]]) + endl
                else:
                    DI = newline + f'{tab}$D_{Ih}~(\si{{km}})${tab}' + tab.join(strDIh_km[thisSubset]) + endl

                # HP ices
                if np.any(boolDIII_km[thisSubset]):
                    DIII = newline + f'{tab}$D_{III}~(\si{{km}})${tab}' + tab.join(strDIII_km[thisSubset]) + endl
                else:
                    DIII = ''
                if np.any(boolDVund_km[thisSubset]):
                    DVund = newline + f'{tab}$D_{Vund}~(\si{{km}})${tab}' + tab.join(strDVund_km[thisSubset]) + endl
                else:
                    DVund = ''
                if np.any(boolDVwet_km[thisSubset]):
                    DVwet = newline + f'{tab}$D_{Vwet}~(\si{{km}})${tab}' + tab.join(strDVwet_km[thisSubset]) + endl
                else:
                    DVwet = ''
                if np.any(boolDVI_km[thisSubset]):
                    DVI = newline + f'{tab}$D_{VI}~(\si{{km}})${tab}' + tab.join(strDVI_km[thisSubset]) + endl
                else:
                    DVI = ''
                if np.any(boolphiIceMax_frac[thisSubset]):
                    phiIce = newline + f'{tab}$\phi_{ice}${FigLbl.phiUnitsParen}{tab}' + tab.join(strphiIceMax_frac[thisSubset]) + endl
                else:
                    phiIce = ''
            else:
                qCon, etaI, DI, DIII, DVund, Docean, DVwet, DVI, sigOcean, phiIce = ('' for _ in range(10))

            Rrock = f'{tab}$R_{rock}~(\si{{km}})${tab}' + tab.join(strRrock_km[thisSubset]) + endl
            if np.any(boolRcore_km[thisSubset]):
                Rcore = f'{tab}$R_{core}~(\si{{km}})${tab}' + tab.join(strRcore_km[thisSubset]) + endl
            else:
                Rcore = ''
            if np.any(boolphiRockMax_frac[thisSubset]):
                phiRock = newline + f'{tab}$\phi_{rock}${FigLbl.phiUnitsParen}{tab}' + tab.join(strphiRockMax_frac[thisSubset]) + endl
            else:
                phiRock = ''

            log.info(f"""{tOpen}
                {header}
                    {compStr}{Mmeas}{Mcalc}{Cmeas}{Ccalc}{rhoRock}
                    {Tb_K}
                    {qSurf}{qCon}{etaI}{DI}{DIII}{DVund}{Docean}{DVwet}{DVI}{sigOcean}{Rsurf}
                    {Rrock}{Rcore}{phiIce}{phiRock}
                {tClose}
                """)

    return
