import numpy as np
import logging as log

def PrintLayerTable(PlanetList, Params):

    # Construct strings for printing
    models = 'Models printed:\n' + '\n'.join([Planet.saveLabel for Planet in PlanetList])
    bodyMass = f'Body mass (kg): {PlanetList[0].Bulk.M_kg:.5e}'
    computedMass = 'Computed mass (kg): ' + ', '.join([f'{Planet.Mtot_kg:.5e}' for Planet in PlanetList])
    inputCMR2 = f'Input C/MR^2: {PlanetList[0].Bulk.Cmeasured} ± {PlanetList[0].Bulk.Cuncertainty}'
    belowCMR2 =    '            (-)  ' + ', '.join([f'{Planet.CMR2less:.5f}' for Planet in PlanetList])
    computedCMR2 = 'Computed C/MR^2: ' + ', '.join([f'{Planet.CMR2mean:.5f}' for Planet in PlanetList])
    aboveCMR2 =    '            (+)  ' + ', '.join([f'{Planet.CMR2more:.5f}' for Planet in PlanetList])
    TbK = 'Tb (K): ' + ', '.join([f'{Planet.Bulk.Tb_K}' for Planet in PlanetList])
    zUpper = 'zb (km): ' + ', '.join([f'{Planet.zb_km}' for Planet in PlanetList])
    oceanThk = 'Ocean thickness D (km): ' + ', '.join([f'{Planet.D_km:.1f}' for Planet in PlanetList])
    zIceI = 'z(km) ice I: ' + ', '.join([f'{Planet.z_m[Planet.Steps.nIbottom]/1e3:.1f}' for Planet in PlanetList])
    zClath = 'z(km) clath: ' + ', '.join([f'{Planet.zClath_m/1e3:.1f}' for Planet in PlanetList])
    zIceIII = 'z(km) ice III: ' + ', '.join([f'{np.max(Planet.z_m[np.abs(Planet.phase) == 3], initial=0)/1e3:.1f}' for Planet in PlanetList])
    zIceV = 'z(km) ice V: ' + ', '.join([f'{np.max(Planet.z_m[np.abs(Planet.phase) == 5], initial=0)/1e3:.1f}' for Planet in PlanetList])
    zIceVI = 'z(km) ice VI: ' + ', '.join([f'{np.max(Planet.z_m[np.abs(Planet.phase) == 6], initial=0)/1e3:.1f}' for Planet in PlanetList])
    wetHydro = 'Total wet hydrosphere: ' + ', '.join([f'{(Planet.z_m[Planet.Steps.nHydro]/1e3 - Planet.zb_km):.1f}' for Planet in PlanetList])
    dzIceIII = 'dz(km) ice III: ' + ', '.join([
        f'{(Planet.z_m[Planet.Steps.nIIIbottom] - Planet.z_m[Planet.Steps.nIbottom + Planet.Steps.nClath])/1e3:.1f}'
        for Planet in PlanetList])
    dzIceV = 'dz(km) ice V: ' + ', '.join([
        f'{(Planet.z_m[Planet.Steps.nSurfIce] - Planet.z_m[Planet.Steps.nIIIbottom])/1e3:.1f}' for Planet in PlanetList])
    dzIceVI = 'dz(km) ice VI: ' + ', '.join([
        f'{Planet.z_m[Planet.Steps.nHydro]/1e3 - Planet.D_km - Planet.z_m[Planet.Steps.nSurfIce]/1e3:.1f}'
        for Planet in PlanetList])
    dzIceVandVI = 'dz(km) ice V + ice VI: ' + ', '.join([
        f'{Planet.z_m[Planet.Steps.nHydro]/1e3 - Planet.D_km - Planet.z_m[Planet.Steps.nIIIbottom]/1e3:.1f}'
        for Planet in PlanetList])

    log.info(f"""
        {models}
        {bodyMass}
        {computedMass}
        {inputCMR2}
        {belowCMR2}
        {computedCMR2}
        {aboveCMR2}
        {TbK}
        {zUpper}
        {oceanThk}
        {zIceI}
        {zClath}
        {zIceIII}
        {zIceV}
        {zIceVI}
        {wetHydro}
        {dzIceIII}
        {dzIceV}
        {dzIceVI}
        {dzIceVandVI}
    """)
    return
