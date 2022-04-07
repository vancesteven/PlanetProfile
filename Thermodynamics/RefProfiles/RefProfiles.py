import numpy as np
import logging as log
from Thermodynamics.FromLiterature.HydroEOS import OceanEOSStruct, GetTfreeze

def CalcRefProfiles(PlanetList, Params):

    comps = np.unique([Planet.Ocean.comp for Planet in PlanetList])
    newRef = {comp: True for comp in comps}

    for Planet in PlanetList:
        if newRef[Planet.Ocean.comp]:
            log.info(f'Calculating reference profiles for {Planet.Ocean.comp}.')

            # Fetch the values we need and initialize
            wList = Params.wRef_ppt[Planet.Ocean.comp]
            Params.nRef[Planet.Ocean.comp] = np.size(wList)
            Params.nRefPts[Planet.Ocean.comp] = Planet.Steps.nRefRho + 0
            Params.rhoRef_kgm3[Planet.Ocean.comp] = np.zeros((Params.nRef[Planet.Ocean.comp], Params.nRefPts[Planet.Ocean.comp]))
            Params.Pref_MPa[Planet.Ocean.comp] = np.linspace(0, Planet.Ocean.PHydroMax_MPa, Params.nRefPts[Planet.Ocean.comp])
            Tref_K = np.arange(220, 450, 1)
            for i,w_ppt in enumerate(wList):
                EOSref = OceanEOSStruct(Planet.Ocean.comp, w_ppt, Params.Pref_MPa[Planet.Ocean.comp], Tref_K, Planet.Ocean.MgSO4elecType)
                Tfreeze_K = np.array([GetTfreeze(EOSref, P_MPa, Tref_K[0], TfreezeRange_K=230) for P_MPa in Params.Pref_MPa[Planet.Ocean.comp]])
                Params.rhoRef_kgm3[Planet.Ocean.comp][i,:] = EOSref.fn_rho_kgm3(Params.Pref_MPa[Planet.Ocean.comp], Tfreeze_K, grid=False)

            # Save to disk for quick reloading
            with open(Params.fNameRef[Planet.Ocean.comp], 'w') as f:
                f.write(f'This file contains melting curve densities for one or more "{Planet.Ocean.comp}" salinity values.\n')
                wListStr = ''
                colHeader = f'P (MPa)'.ljust(24)
                for w_ppt in wList:
                    wListStr = wListStr + f' {w_ppt:.3f},'
                    colHeader = ' '.join([colHeader, f'rho_{w_ppt:.3f} (kg/m3)'.ljust(24)])

                f.write(f'  w_ppt = {wListStr[1:-1]}\n')
                f.write(colHeader + '\n')

                for i in range(Planet.Steps.nRefRho):
                    line = f'{Params.Pref_MPa[Planet.Ocean.comp][i]:24.17e}'
                    for j in range(Params.nRef[Planet.Ocean.comp]):
                        line = ' '.join([line, f'{Params.rhoRef_kgm3[Planet.Ocean.comp][j,i]:24.17e}'])
                    f.write(line + '\n')

            newRef[Planet.Ocean.comp] = False

    return Params


def ReloadRefProfiles(PlanetList, Params):

    comps = np.unique([Planet.Ocean.comp for Planet in PlanetList])
    newRef = {comp: True for comp in comps}

    for Planet in PlanetList:
        if newRef[Planet.Ocean.comp]:

            with open(Params.fNameRef[Planet.Ocean.comp]) as f:
                _ = f.readline()
                Params.wRef_ppt[Planet.Ocean.comp] = np.array(f.readline().split('=')[-1].split(',')).astype(np.float_)
            try:
                PrhoRef = np.loadtxt(Params.fNameRef[Planet.Ocean.comp], skiprows=3, unpack=False)
            except:
                raise ValueError(f'Reference melting curves for {Planet.Ocean.comp} have not been generated. '
                                  'Re-run with CALC_NEW_REF = True in config.py')
            PrhoRef = PrhoRef.T
            Params.nRef[Planet.Ocean.comp] = np.size(Params.wRef_ppt[Planet.Ocean.comp])
            Params.nRefPts[Planet.Ocean.comp] = np.shape(PrhoRef)[1]

            Params.Pref_MPa[Planet.Ocean.comp] = PrhoRef[0,:]
            Params.rhoRef_kgm3[Planet.Ocean.comp] = PrhoRef[1:,:]
            newRef[Planet.Ocean.comp] = False

    return Params
