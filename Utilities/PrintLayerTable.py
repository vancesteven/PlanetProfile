import numpy as np
import logging as log
from Utilities.defineStructs import Constants

def PrintLayerTable(PlanetList, Params):

    # Construct strings for printing
    # Create newline+indent constant for optional lines
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
    zIceI = 'z(km) ice I: ' + ', '.join([f'{Planet.z_m[Planet.Steps.nIbottom]/1e3:.1f}' for Planet in PlanetList])

    # Optional output strings
    poreSigma, phiRockMax, phiIceMax, zClath, zIceIII, zIceV, zIceVI, dzClath, dzIceIII, dzIceV, dzIceVI, \
    dzIceVandVI \
        = ('' for _ in range(12))

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
    {TbK}
    {TsilK}
    {TCMBK}
    {TcenterK}
    {qOcean}
    {qSurf}
    {oceanSigma}{poreSigma}{phiRockMax}{phiIceMax}
    
    {zIceI}{zClath}{zIceIII}{zIceV}{zIceVI}{dzClath}{dzIceIII}{dzIceV}{dzIceVI}{dzIceVandVI}
    {wetHydro}
    
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
