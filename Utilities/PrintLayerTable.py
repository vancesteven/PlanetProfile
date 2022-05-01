import numpy as np
import logging as log
from Utilities.defineStructs import Constants

def GetLayerMeans(PlanetList, Params):
    # For calculating layer means we didn't need at any other point in our analysis,
    # but that we might want to consider in comparing/analyzing profiles.

    for Planet in PlanetList:

        if Planet.Do.NO_H2O:
            Planet.Ocean.rhoCondImean_kgm3 = np.nan
            Planet.Ocean.rhoConvImean_kgm3 = np.nan
        else:
            iCond = Planet.z_m[:-1] < Planet.eLid_m
            iConv = np.logical_and(Planet.z_m[:-1] >= Planet.eLid_m, Planet.z_m[:-1] < Planet.zb_km*1e3)
            Planet.Ocean.rhoCondImean_kgm3 = np.sum(Planet.MLayer_kg[iCond]) / np.sum(Planet.VLayer_m3[iCond])
            Planet.Ocean.rhoConvImean_kgm3 = np.sum(Planet.MLayer_kg[iConv]) / np.sum(Planet.VLayer_m3[iConv])
            # Get mean ice layer conductivities, neglecting spherical effects
            Planet.Ocean.sigmaCondImean_Sm = np.mean(Planet.sigma_Sm[iCond])
            if np.sum(iConv) > 0:
                Planet.Ocean.sigmaConvImean_Sm = np.mean(Planet.sigma_Sm[iConv])
            else:
                Planet.Ocean.sigmaConvImean_Sm = np.nan

            if Planet.Do.BOTTOM_ICEIII or Planet.Do.BOTTOM_ICEV:
                iCondIII = np.logical_and(Planet.z_m[:-1] >= Planet.z_m[Planet.Steps.nIbottom],
                                                       Planet.z_m[:-1] < Planet.z_m[Planet.Steps.nIbottom] + Planet.eLidIII_m)
                iConvIII = np.logical_and(Planet.z_m[:-1] >= Planet.z_m[Planet.Steps.nIbottom] + Planet.eLidIII_m,
                                                       Planet.z_m[:-1] < Planet.z_m[Planet.Steps.nIIIbottom])
                Planet.Ocean.rhoCondIIImean_kgm3 = np.sum(Planet.MLayer_kg[iCondIII]) / np.sum(VLayer_m3[iCondIII])
                Planet.Ocean.rhoConvIIImean_kgm3 = np.sum(Planet.MLayer_kg[iConvIII]) / np.sum(VLayer_m3[iConvIII])
                Planet.Ocean.sigmaCondIIImean_Sm = np.mean(Planet.sigma_Sm[iCondIII])
                if np.sum(iConvIII) > 0:
                    Planet.Ocean.sigmaConvIIImean_Sm = np.mean(Planet.sigma_Sm[iConvIII])
                else:
                    Planet.Ocean.sigmaConvIIImean_Sm = np.nan
            if Planet.Do.BOTTOM_ICEV:
                iCondV = np.logical_and(Planet.z_m[:-1] >= Planet.z_m[Planet.Steps.nIIIbottom],
                                                     Planet.z_m[:-1] < Planet.z_m[Planet.Steps.nIIIbottom] + Planet.eLidV_m)
                iConvV = np.logical_and(Planet.z_m[:-1] >= Planet.z_m[Planet.Steps.nIIIbottom] + Planet.eLidV_m,
                                                     Planet.z_m[:-1] < Planet.z_m[Planet.Steps.nSurfIce])
                Planet.Ocean.rhoCondVmean_kgm3 = np.sum(Planet.MLayer_kg[iCondV]) / np.sum(VLayer_m3[iCondV])
                Planet.Ocean.rhoConvVmean_kgm3 = np.sum(Planet.MLayer_kg[iConvV]) / np.sum(VLayer_m3[iConvV])
                if np.sum(iConvV) > 0:
                    Planet.Ocean.sigmaConvVmean_Sm = np.mean(Planet.sigma_Sm[iConvV])
                else:
                    Planet.Ocean.sigmaConvVmean_Sm = np.nan

            Planet.Sil.sigmaMean_Sm = np.mean(Planet.sigma_Sm[np.logical_and(Planet.phase >= Constants.phaseSil, Planet.phase < Constants.phaseSil + 10)])
            if np.sum(Planet.phase >= Constants.phaseFe) > 0:
                Planet.Core.sigmaMean_Sm = np.mean(Planet.sigma_Sm[Planet.phase >= Constants.phaseFe])
            else:
                Planet.Core.sigmaMean_Sm = np.nan

    return PlanetList

def PrintLayerTable(PlanetList, Params):

    # Construct strings for printing
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
    zIceI = 'z(km) ice I: ' + ', '.join([f'{Planet.z_m[Planet.Steps.nIbottom]/1e3:.1f}' for Planet in PlanetList])

    # Optional output strings
    poreSigma, phiRockMax, phiIceMax, zClath, zIceIII, zIceV, zIceVI, dzClath, dzIceIII, dzIceV, dzIceVI, \
    dzIceVandVI, RaIII, RaCritIII, eLidIII, DconvIII, deltaTBLIII, RaV, RaCritV, eLidV, DconvV, deltaTBLV \
        = ('' for _ in range(22))

    # Porosity
    yesPorousRock = np.any([Planet.Do.POROUS_ROCK for Planet in PlanetList])
    yesPorousIce = np.any([Planet.Do.POROUS_ICE for Planet in PlanetList])
    if Params.ALWAYS_SHOW_PHI or yesPorousRock:
        poreSigma = f'{endl}Mean pore σ (S/m): ' + ', '.join([f'{Planet.Sil.sigmaPoreMean_Sm:.2f}' for Planet in PlanetList])
        phiRockMax = f'{endl}{endl}Max ϕsil (%): ' + ', '.join([f'{Planet.phi_frac[Planet.Steps.nHydro]*100:.1f}' for Planet in PlanetList])
    if Params.ALWAYS_SHOW_PHI or yesPorousIce:
        phiIceMax = f'{endl}Max ϕice (%): ' + ', '.join([f'{Planet.phi_frac[0]*100:.1f}' for Planet in PlanetList])

    # Temperatures and heat flux
    TbK = 'T_b (K): ' + ', '.join([f'{Planet.Bulk.Tb_K}' for Planet in PlanetList])
    TsilK = 'T_sil (K): ' + ', '.join([f'{Planet.T_K[Planet.Steps.nHydro]:.1f}' for Planet in PlanetList])
    TCMBK = 'T_CMB (K): ' + ', '.join([f'{Planet.T_K[Planet.Steps.nHydro + Planet.Steps.nSil]:.1f}' for Planet in PlanetList])
    TcenterK = 'T_center (K): ' + ', '.join([f'{Planet.T_K[Planet.Steps.nTotal-1]:.1f}' for Planet in PlanetList])
    qOcean = 'qOcBot (mW/m^2): ' + ', '.join([f'{1e3*Planet.Ocean.QfromMantle_W/(4*np.pi*Planet.Sil.Rmean_m**2):.1f}' for Planet in PlanetList])
    qSurf = 'qSurf (mW/m^2): ' + ', '.join([f'{1e3*Planet.Ocean.QfromMantle_W/(4*np.pi*Planet.Bulk.R_m**2):.1f}' for Planet in PlanetList])

    # Only print HP ice values if they are present or forced on
    yesClath = np.any([np.any(abs(Planet.phase) == Constants.phaseClath) for Planet in PlanetList])
    yesIceIII = np.any([np.any(abs(Planet.phase) == 3) for Planet in PlanetList])
    yesIceVund = np.any([np.any(Planet.phase == -5) for Planet in PlanetList])
    yesIceVI = np.any([np.any(abs(Planet.phase) == 6) for Planet in PlanetList])
    yesIceVandVI = np.any([np.any(Planet.phase == 5) or np.any(Planet.phase == 6) for Planet in PlanetList])
    if Params.ALWAYS_SHOW_HP or yesClath:
        zClath = f'{endl}z(km) clath: ' + ', '.join([f'{Planet.zClath_m/1e3:.1f}' for Planet in PlanetList])
        dzClath = f'{endl}dz(km) clath: ' + ', '.join([f'{(np.max(Planet.z_m[:-1][abs(Planet.phase)==Constants.phaseClath], initial=0) - np.min(Planet.z_m[:-1][abs(Planet.phase)==Constants.phaseClath], initial=0))/1e3:.1f}' for Planet in PlanetList])
    if Params.ALWAYS_SHOW_HP or yesIceIII:
        zIceIII = f'{endl}z(km) ice III: ' + ', '.join([f'{np.max(Planet.z_m[:-1][abs(Planet.phase) == 3], initial=0)/1e3:.1f}' for Planet in PlanetList])
        dzIceIII = f'{endl}dz(km) ice III: ' + ', '.join([
            f'{(Planet.z_m[Planet.Steps.nIIIbottom] - Planet.z_m[Planet.Steps.nIbottom + Planet.Steps.nClath])/1e3:.1f}'
            for Planet in PlanetList])
    if Params.ALWAYS_SHOW_HP or yesIceVund:
        zIceV = f'{endl}z(km) ice V (und): ' + ', '.join([f'{np.max(Planet.z_m[:-1][Planet.phase == -5], initial=0)/1e3:.1f}' for Planet in PlanetList])
        dzIceV = f'{endl}dz(km) ice V (und): ' + ', '.join([
            f'{Planet.zb_km - Planet.z_m[Planet.Steps.nIIIbottom]/1e3:.1f}' for Planet in PlanetList])
    if Params.ALWAYS_SHOW_HP or yesIceVI:
        zIceVI = f'{endl}z(km) ice VI: ' + ', '.join([f'{np.max(Planet.z_m[:-1][abs(Planet.phase) == 6], initial=0)/1e3:.1f}' for Planet in PlanetList])
        dzIceVI = f'{endl}dz(km) ice VI: ' + ', '.join([
            f'{(np.max(Planet.z_m[:-1][abs(Planet.phase)==6], initial=0) - np.min(Planet.z_m[:-1][abs(Planet.phase)==6], initial=0))/1e3:.1f}'
            for Planet in PlanetList])
    if Params.ALWAYS_SHOW_HP or yesIceVandVI:
        dzIceVandVI = f'{endl}dz(km) ice V (wet) + ice VI: ' + ', '.join([
            f'{Planet.z_m[Planet.Steps.nHydro]/1e3 - Planet.D_km - Planet.zb_km:.1f}'
            for Planet in PlanetList])
    RaI = f'Ice shell Rayleigh number Ra: ' + ', '.join([f'{Planet.RaConvect:.2e}' for Planet in PlanetList])
    RaCritI = f'Critical Rayleigh number Ra_crit: ' + ', '.join([f'{Planet.RaCrit:.2e}' for Planet in PlanetList])
    eLidI = f'Conductive lid thickness e (km): ' + ', '.join([f'{Planet.eLid_m/1e3:.2f}' for Planet in PlanetList])
    DconvI = f'Convecting layer thickness D_conv (km): ' + ', '.join([f'{Planet.Dconv_m/1e3:.2f}' for Planet in PlanetList])
    deltaTBLI = f'Lower TBL thickness δ (km): ' + ', '.join([f'{Planet.deltaTBL_m/1e3:.2f}' for Planet in PlanetList])

    MfracH2O = f'Mass fraction H2O (%): ' + ', '.join([f'{100*Planet.MH2O_kg/Planet.Mtot_kg:.1f}' for Planet in PlanetList])
    MfracSalt = f'Mass fraction solutes (%): ' + ', '.join([f'{100*Planet.Msalt_kg/Planet.Mtot_kg:.1f}' for Planet in PlanetList])
    MfracSil = f'Mass fraction silicates (%): ' + ', '.join([f'{100*Planet.Mrock_kg/Planet.Mtot_kg:.1f}' for Planet in PlanetList])
    MfracFe = f'Mass fraction iron core (%): ' + ', '.join([f'{100*Planet.Mcore_kg/Planet.Mtot_kg:.1f}' for Planet in PlanetList])

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
    
    {zIceI}{zClath}{zIceIII}{zIceV}{zIceVI}{dzClath}{dzIceIII}{dzIceV}{dzIceVI}{dzIceVandVI}
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


def PrintLayerTableLatex(PlanetList, Params):

    # Construct a Latex table for collating comparisons
    # Table horizontal division constant
    tab = ' & '
    # Table vertical division constant
    endl = r' \\ \hline'
    # Vertical lines for table, if present
    if Params.LATEX_VLINES:
        v = ' | '
    else:
        v = ' '
    # Table begin
    tOpen = r'\begin{tabular}{' + v + v.join(['l'] + ['c' for _ in PlanetList]) + v + '}'
    # Table end
    tClose = r'\end{tabular}'
    # Header line
    header = ''

    # Layer table labels
    columns = [r'\textbf{Layer}', r'\textbf{Radius ($\si{km}$)}', r'\textbf{Density ($\si{kg/m^3}$)}',
                              r'\textbf{Thickness ($\si{km}$)}', r'\textbf{Conductivity ($\si{S/m}$)}']
    headerLayers = '\hline\n' + tab.join(columns) + endl
    tOpenLayers = r'\begin{tabular}{' + v + v.join(['l' for _ in columns]) + v + '}'
    condIceLbl = 'Conductive ice'
    convIceLbl = 'Convective ice'
    oceanLbl = 'ocean'
    silLbl = 'Mantle'
    coreLbl = 'Core'

    log.info('Layer tables:')
    for Planet in PlanetList:
        condIceLayers = tab.join([condIceLbl, f'{Planet.Bulk.R_m/1e3:.1f}', f'{Planet.Ocean.rhoCondImean_kgm3:.0f}',
                            f'{Planet.eLid_m/1e3:.1f}', f'{Planet.Ocean.sigmaCondImean_Sm:.1e}']) + endl
        convIceLayers = tab.join([convIceLbl, f'{(Planet.Bulk.R_m - Planet.eLid_m)/1e3:.1f}', f'{Planet.Ocean.rhoConvImean_kgm3:.0f}',
                            f'{Planet.Dconv_m/1e3:.1f}', f'{Planet.Ocean.sigmaConvImean_Sm:.1e}']) + endl
        if Planet.Do.NO_H2O:
            salt = 'No'
        elif Planet.Ocean.wOcean_ppt == 0:
            salt = r'Pure~\ce{H2O}'
        else:
            salt = f'$\SI{{{Planet.Ocean.wOcean_ppt:.1f}}}{{g/kg\,\ce{{{Planet.Ocean.comp}}}}}$'
        oceanLayers = tab.join([f'{salt} {oceanLbl}', f'{Planet.Bulk.R_m/1e3 - Planet.zb_km:.1f}', f'{Planet.Ocean.rhoMean_kgm3:.0f}',
                            f'{Planet.D_km:.1f}', f'{Planet.Ocean.sigmaMean_Sm:.2f}']) + endl
        silLayers = tab.join([silLbl, f'{Planet.Sil.Rmean_m/1e3:.1f}', f'{Planet.Sil.rhoMean_kgm3:.0f}',
                            f'{(Planet.Sil.Rmean_m - Planet.Core.Rmean_m)/1e3:.1f}', f'{Planet.Sil.sigmaMean_Sm:.1e}']) + endl
        coreLayers = tab.join([coreLbl, f'{Planet.Core.Rmean_m/1e3:.1f}', f'{Planet.Core.rhoMean_kgm3:.0f}',
                            f'{Planet.Core.Rmean_m/1e3:.1f}', f'{Planet.Core.sigmaMean_Sm:.1e}']) + endl
        log.info(f"""{Planet.saveLabel}
        {tOpenLayers}
        {headerLayers}
            {condIceLayers}
            {convIceLayers}
            {oceanLayers}
            {silLayers}
            {coreLayers}
        {tClose}
        """)

    log.info(f"""Comparison table:
    {tOpen}
    {header}
    {tClose}
    """)
    return
