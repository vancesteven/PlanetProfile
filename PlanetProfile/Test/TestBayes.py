from PlanetProfile.Main import InitBayes, PlanetProfile, UpdateRun
import logging

# Assign logger
log = logging.getLogger('PlanetProfile')

def TestBayes(bodyname):

    fEnd = 'Bayes'
    # Load in parameters set in base run config for this body
    Planet, Params = InitBayes(bodyname, fEnd)

    # Optional run of the contents of Europa/PPEuropaBayes.py
    log.profile(f'Running initial config for {bodyname} UpdateRun.')
    Planet, Params = PlanetProfile(Planet, Params)

    # Starting outputs
    Texc_hr = Planet.Magnetic.Texc_hr
    omegaExc_radps = Planet.Magnetic.Texc_hr
    nLin = Planet.Magnetic.nLin
    mLin = Planet.Magnetic.mLin
    Binm_nT = Planet.Magnetic.BinmLin_nT
    log.profile(f'Starting w_ppt: {Planet.Ocean.wOcean_ppt}. Binm_nT: {Binm_nT}')

    # Update to several new settings for ocean salinity and run each
    for w_ppt in [3, 10, 20]:
        updateTo = {'wOcean_ppt': w_ppt}
        Planet, Params = UpdateRun(Planet, Params, changes=updateTo)
        # Updated induced moments
        Binm_nT = Planet.Magnetic.BinmLin_nT
        log.profile(f'New w_ppt: {w_ppt}. Binm_nT dipoles: {Binm_nT[:,:3]}')

    # Test all the options written in
    testDict = {
        'wOcean_ppt': 30,
        'Tb_K': 268.5,
        'rhoSilInput_kgm3': 3450,
        'ionosTop_km': 50,
        'sigmaIonos_Sm': 1e-4,
        'silPhi_frac': 0.3,
        'silPclosure_MPa': 340,
        'icePhi_frac': 0.7,
        'icePclosure_MPa': 21,
        'Htidal_Wm3': 2e-10,
        'Qrad_Wkg': 1.4e-14,
        'qSurf_Wm2': 135e-3,
        'compOcean': 'MgSO4',
        'compSil': 'CM_hydrous_differentiated_Ganymede_Core100Fe000S_excluding_fluid_properties.tab',
        'wFeCore_ppt': 1000
    }

    # Loop over each option setting, one at a time. Multiple can be set at once but
    # this setup is just to check that they all work correctly.
    for key, value in testDict.items():
        updateTo = {key: value}
        Planet, Params = UpdateRun(Planet, Params, changes=updateTo)
        # Updated induced moments
        Binm_nT = Planet.Magnetic.BinmLin_nT
        log.profile(f'New {key}: {value}. Binm_nT dipoles: {Binm_nT[:,:3]}')

    return Planet, Params

if __name__ == '__main__':
    _, _ = TestBayes('Europa')
