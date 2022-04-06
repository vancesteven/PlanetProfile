import numpy as np
from Thermodynamics.FromLiterature.HydroEOS import OceanEOSStruct, GetTfreeze

def CalcRefProfiles(Planet, Params):

    # Fetch the values we need and initialize
    wList = Params.wRef_ppt[Planet.Ocean.comp]
    Params.nRef = np.size(wList)
    Params.nRefPts = Planet.Steps.nRefRho + 0
    Params.rhoRef_kgm3 = np.zeros((Params.nRef, Params.nRefPts))
    Params.Pref_MPa = np.linspace(0, Planet.Ocean.PHydroMax_MPa, Params.nRefPts)
    Tref_K = np.arange(220, 320, 1)
    for i,w_ppt in enumerate(wList):
        EOSref = OceanEOSStruct(Planet.Ocean.comp, w_ppt, Params.Pref_MPa, Tref_K, None)
        Tfreeze_K = np.array([GetTfreeze(EOSref, P_MPa, Tref_K[0], TfreezeRange_K=100) for P_MPa in Params.Pref_MPa])
        Params.rhoRef_kgm3[i,:] = EOSref.fn_rho_kgm3(Params.Pref_MPa, Tfreeze_K, grid=False)

    # Save to disk for quick reloading
    with open(Params.fnameRef, 'w') as f:
        f.write(f'This file contains melting curve densities for one or more "{Planet.Ocean.comp}" salinity values.\n')
        wListStr = ''
        colHeader = f'P (MPa)'.ljust(24)
        for w_ppt in wList:
            wListStr = wListStr + f' {w_ppt:.3f},'
            colHeader = ' '.join([colHeader, f'rho_{w_ppt:.3f} (kg/m3)'.ljust(24)])

        f.write(f'  w_ppt = {wListStr[1:-1]}\n')
        f.write(colHeader + '\n')

        for i in range(Planet.Steps.nRefRho):
            line = f'{Params.Pref_MPa[i]:24.17e}'
            for j in range(Params.nRef):
                line = ' '.join([line, f'{Params.rhoRef_kgm3[j,i]:24.17e}'])
            f.write(line + '\n')

    return Params


def ReloadRefProfiles(Planet, Params):

    with open(Params.fnameRef) as f:
        _ = f.readline()
        Params.wRef_ppt[Planet.Ocean.comp] = np.array(f.readline().split('=')[-1].split(',')).astype(np.float_)
    PrhoRef = np.loadtxt(Params.fnameRef, skiprows=3, unpack=False)
    PrhoRef = PrhoRef.T
    Params.nRef = np.size(Params.wRef_ppt[Planet.Ocean.comp])
    Params.nRefPts = np.shape(PrhoRef)[1]

    Params.Pref_MPa = PrhoRef[0,:]
    Params.rhoRef_kgm3 = PrhoRef[1:,:]

    return Params
