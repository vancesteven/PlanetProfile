import os
import numpy as np
import logging
from PlanetProfile import _ROOT
from PlanetProfile.Thermodynamics.HydroEOS import GetOceanEOS, GetTfreeze
from PlanetProfile.Utilities.defineStructs import EOSlist

# Assign logger
log = logging.getLogger('PlanetProfile')

def CalcRefProfiles(PlanetList, Params):

    comps = np.unique([Planet.Ocean.comp for Planet in PlanetList])
    newRef = {comp: True for comp in comps}
    maxPmax = np.max([Planet.P_MPa[Planet.Steps.nHydro-1] for Planet in PlanetList])

    for Planet in PlanetList:
        if newRef[Planet.Ocean.comp] and Planet.Ocean.comp != 'none' and not Planet.Do.VARIABLE_COMP_OCEAN:
            wList = Params.wRef_ppt[Planet.Ocean.comp]

            thisRefLabel = f'{Planet.Ocean.comp}' + ','.join([f'{w_ppt}' for w_ppt in wList])
            thisRefRange = maxPmax
            if thisRefLabel in EOSlist.loaded.keys() and thisRefRange <= EOSlist.ranges[thisRefLabel]:
                log.debug(f'Reference profiles for {Planet.Ocean.comp} already loaded. Reusing existing.')
                Params.Pref_MPa[Planet.Ocean.comp], Params.rhoRef_kgm3[Planet.Ocean.comp] = EOSlist.loaded[thisRefLabel]
                newRef[Planet.Ocean.comp] = False
            else:
                log.info(f'Calculating reference profiles for {Planet.Ocean.comp} at {{' + ','.join([f'{w_ppt}' for w_ppt in wList]) + '} ppt.')

                # Fetch the values we need and initialize
                Params.nRef[Planet.Ocean.comp] = np.size(wList)
                Params.nRefPts[Planet.Ocean.comp] = Params.nRefRho + 0
                Params.rhoRef_kgm3[Planet.Ocean.comp] = np.zeros((Params.nRef[Planet.Ocean.comp], Params.nRefPts[Planet.Ocean.comp]))
                if Params.PrefOverride_MPa is None:
                    if Planet.Ocean.EOS is None:
                        PmaxEOS = maxPmax
                    else:
                        PmaxEOS = Planet.Ocean.EOS.Pmax
                    Pmax = np.minimum(maxPmax, PmaxEOS)
                else:
                    Pmax = Params.PrefOverride_MPa
                Params.Pref_MPa[Planet.Ocean.comp] = np.linspace(0.1, Pmax, Params.nRefPts[Planet.Ocean.comp])
                Tref_K = np.arange(220, 450, 0.05)
                for i, w_ppt in enumerate(wList):
                    EOSref = GetOceanEOS(Planet.Ocean.comp, w_ppt, Params.Pref_MPa[Planet.Ocean.comp], Tref_K, Planet.Ocean.MgSO4elecType,
                            rhoType=Planet.Ocean.MgSO4rhoType, scalingType=Planet.Ocean.MgSO4scalingType, phaseType='lookup',
                            EXTRAP=Params.EXTRAP_REF, FORCE_NEW=Params.FORCE_EOS_RECALC, MELT=True,
                            VARIABLE_COMP=Planet.Do.VARIABLE_COMP_OCEAN, Pstratified_MPa=Planet.Ocean.Pstratified_MPa,
                            wStratified_ppt=Planet.Ocean.wStratified_ppt,compStratified=Planet.Ocean.compStratified,
                            CONTINUOUS_SALINITY=Planet.Do.CONTINUOUS_SALINITY)
                    if EOSref.propsPmax < Pmax or EOSref.Pmax < Pmax:
                        Params.Pref_MPa[Planet.Ocean.comp] = np.linspace(Params.Pref_MPa[Planet.Ocean.comp][0], np.minimum(EOSref.propsPmax, EOSref.Pmax),
                                                                         Params.nRefPts[Planet.Ocean.comp])
                    try:
                        Tfreeze_K = np.array([GetTfreeze(EOSref, P_MPa, Tref_K[0], TfreezeRange_K=230) for P_MPa in Params.Pref_MPa[Planet.Ocean.comp]])
                    except:
                        raise RuntimeError(f'Unable to calculate reference melting curve for {Planet.Ocean.comp} with ' +
                                           f'maximum Pref_MPa = {Params.Pref_MPa[Planet.Ocean.comp][-1]}. Try to recalculate ' +
                                           'with new models by setting BOTH Params.CALC_NEW and Params.CALC_NEW_REF to True ' +
                                           'in configPP.py.')
                    Params.rhoRef_kgm3[Planet.Ocean.comp][i,:] = EOSref.fn_rho_kgm3(Params.Pref_MPa[Planet.Ocean.comp], Tfreeze_K)

                # Save to disk for quick reloading
                with open(os.path.join(_ROOT, 'Thermodynamics', 'RefProfiles', Params.fNameRef[Planet.Ocean.comp]), 'w') as f:
                    f.write(f'This file contains melting curve densities for one or more "{Planet.Ocean.comp}" salinity values.\n')
                    wListStr = ''
                    colHeader = f'P (MPa)'.ljust(24)
                    for w_ppt in wList:
                        wListStr = wListStr + f' {w_ppt:.3f},'
                        colHeader = ' '.join([colHeader, f'rho_{w_ppt:.3f} (kg/m3)'.ljust(24)])

                    f.write(f'  w_ppt = {wListStr[1:-1]}\n')
                    f.write(colHeader + '\n')

                    for i in range(Params.nRefRho):
                        line = f'{Params.Pref_MPa[Planet.Ocean.comp][i]:24.17e}'
                        for j in range(Params.nRef[Planet.Ocean.comp]):
                            line = ' '.join([line, f'{Params.rhoRef_kgm3[Planet.Ocean.comp][j,i]:24.17e}'])
                        f.write(line + '\n')

                EOSlist.loaded[thisRefLabel] = Params.Pref_MPa[Planet.Ocean.comp], Params.rhoRef_kgm3[Planet.Ocean.comp]
                EOSlist.ranges[thisRefLabel] = Pmax
                newRef[Planet.Ocean.comp] = False

        elif Planet.Do.VARIABLE_COMP_OCEAN:
            log.error('VARIABLE_COMP_OCEAN not yet implemented for plotting refProfiles.')
            Params.nRef[Planet.Ocean.comp] = 0
            Params.Pref_MPa[Planet.Ocean.comp] = np.nan * np.empty(Params.nRefRho)
            Params.rhoRef_kgm3[Planet.Ocean.comp] = np.nan * np.empty(Params.nRefRho)

    return Params


def ReloadRefProfiles(PlanetList, Params):

    comps = np.unique([Planet.Ocean.comp for Planet in PlanetList])
    newRef = {comp: True for comp in comps}

    for Planet in PlanetList:
        if newRef[Planet.Ocean.comp] and Planet.Ocean.comp != 'none':

            fNameRefReload = os.path.join(_ROOT, 'Thermodynamics', 'RefProfiles', Params.fNameRef[Planet.Ocean.comp])
            if not os.path.isfile(fNameRefReload):
                raise RuntimeError(f'CALC_NEW_REF is set to False, but a reference profile for {Planet.Ocean.comp} ' +
                                   'was not found. Try running again with CALC_NEW_REF set to True in configPP.py.')
            with open(fNameRefReload) as f:
                _ = f.readline()
                Params.wRef_ppt[Planet.Ocean.comp] = np.array(f.readline().split('=')[-1].split(',')).astype(np.float_)
            PrhoRef = np.loadtxt(fNameRefReload, skiprows=3, unpack=False)
            PrhoRef = PrhoRef.T
            Params.nRef[Planet.Ocean.comp] = np.size(Params.wRef_ppt[Planet.Ocean.comp])
            Params.nRefPts[Planet.Ocean.comp] = np.shape(PrhoRef)[1]

            Params.Pref_MPa[Planet.Ocean.comp] = PrhoRef[0,:]
            Params.rhoRef_kgm3[Planet.Ocean.comp] = PrhoRef[1:,:]
            newRef[Planet.Ocean.comp] = False

    return Params
