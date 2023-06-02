import os
import numpy as np
import logging
from PlanetProfile import _ROOT
from PlanetProfile.Thermodynamics.HydroEOS import GetOceanEOS, GetTfreeze
from PlanetProfile.Utilities.defineStructs import EOSlist, Constants

# Assign logger
log = logging.getLogger('PlanetProfile')

def CalcRefProfiles(PlanetList, Params):

    if np.any([Planet.Ocean.comp == Constants.varCompStr for Planet in PlanetList]):
        comps = [Planet.Ocean.comp if Planet.Ocean.comp != Constants.varCompStr else np.unique(Planet.Ocean.compStratified) for Planet in PlanetList]
        # Flatten list in for stratified oceans
        comps = np.unique([comp for items in comps for comp in items])
    else:
        comps = np.unique([Planet.Ocean.comp for Planet in PlanetList])
    newRef = {comp: True for comp in comps}
    maxPmax = {comp: np.max([Planet.P_MPa[Planet.Steps.nHydro-1] for Planet in PlanetList if Planet.Ocean.comp in [comp, Constants.varCompStr]]) for comp in comps}

    for comp in comps:
        if newRef[comp] and comp != 'none':
            wList = Params.wRef_ppt[comp]

            thisRefLabel = f'{comp}' + ','.join([f'{w_ppt}' for w_ppt in wList])
            thisRefRange = maxPmax[comp]
            if thisRefLabel in EOSlist.loaded.keys() and thisRefRange <= EOSlist.ranges[thisRefLabel]:
                log.debug(f'Reference profiles for {comp} already loaded. Reusing existing.')
                Params.Pref_MPa[comp], Params.rhoRef_kgm3[comp] = EOSlist.loaded[thisRefLabel]
                newRef[comp] = False
            else:
                log.info(f'Calculating reference profiles for {comp} at {{' + ','.join([f'{w_ppt}' for w_ppt in wList]) + '} ppt.')

                # Fetch the values we need and initialize
                Params.nRef[comp] = np.size(wList)
                Params.nRefPts[comp] = Params.nRefRho + 0
                Params.rhoRef_kgm3[comp] = np.zeros((Params.nRef[comp], Params.nRefPts[comp]))
                Pmax = maxPmax[comp] if Params.PrefOverride_MPa is None else Params.PrefOverride_MPa
                Params.Pref_MPa[comp] = np.linspace(0.1, Pmax, Params.nRefPts[comp])
                Tref_K = np.arange(220, 450, 0.05)
                for i, w_ppt in enumerate(wList):
                    EOSref = GetOceanEOS(comp, w_ppt, Params.Pref_MPa[comp], Tref_K, PlanetList[0].Ocean.MgSO4elecType,
                            rhoType=PlanetList[0].Ocean.MgSO4rhoType, scalingType=PlanetList[0].Ocean.MgSO4scalingType, phaseType='lookup',
                            EXTRAP=Params.EXTRAP_REF, FORCE_NEW=Params.FORCE_EOS_RECALC, MELT=True, VARIABLE_COMP=False)
                    if EOSref.propsPmax < Pmax or EOSref.Pmax < Pmax:
                        Params.Pref_MPa[comp] = np.linspace(Params.Pref_MPa[comp][0], np.minimum(EOSref.propsPmax, EOSref.Pmax),
                                                                         Params.nRefPts[comp])
                    try:
                        Tfreeze_K = np.array([GetTfreeze(EOSref, P_MPa, Tref_K[0], TfreezeRange_K=230) for P_MPa in Params.Pref_MPa[comp]])
                    except:
                        raise RuntimeError(f'Unable to calculate reference melting curve for {comp} with ' +
                                           f'maximum Pref_MPa = {Params.Pref_MPa[comp][-1]}. Try to recalculate ' +
                                           'with new models by setting BOTH Params.CALC_NEW and Params.CALC_NEW_REF to True ' +
                                           'in configPP.py.')
                    Params.rhoRef_kgm3[comp][i,:] = EOSref.fn_rho_kgm3(Params.Pref_MPa[comp], Tfreeze_K)

                # Save to disk for quick reloading
                with open(os.path.join(_ROOT, 'Thermodynamics', 'RefProfiles', Params.fNameRef[comp]), 'w') as f:
                    f.write(f'This file contains melting curve densities for one or more "{comp}" salinity values.\n')
                    wListStr = ''
                    colHeader = f'P (MPa)'.ljust(24)
                    for w_ppt in wList:
                        wListStr = wListStr + f' {w_ppt:.3f},'
                        colHeader = ' '.join([colHeader, f'rho_{w_ppt:.3f} (kg/m3)'.ljust(24)])

                    f.write(f'  w_ppt = {wListStr[1:-1]}\n')
                    f.write(colHeader + '\n')

                    for i in range(Params.nRefRho):
                        line = f'{Params.Pref_MPa[comp][i]:24.17e}'
                        for j in range(Params.nRef[comp]):
                            line = ' '.join([line, f'{Params.rhoRef_kgm3[comp][j,i]:24.17e}'])
                        f.write(line + '\n')

                EOSlist.loaded[thisRefLabel] = Params.Pref_MPa[comp], Params.rhoRef_kgm3[comp]
                EOSlist.ranges[thisRefLabel] = Pmax
                newRef[comp] = False

    return Params


def ReloadRefProfiles(PlanetList, Params):

    if np.any([Planet.Ocean.comp == Constants.varCompStr for Planet in PlanetList]):
        comps = [Planet.Ocean.comp if Planet.Ocean.comp != Constants.varCompStr else np.unique(Planet.Ocean.compStratified) for Planet in PlanetList]
        # Flatten list in for stratified oceans
        comps = np.unique([comp for items in comps for comp in items])
    else:
        comps = np.unique([Planet.Ocean.comp for Planet in PlanetList])
    newRef = {comp: True for comp in comps}

    for comp in comps:
        if newRef[comp] and comp != 'none':

            fNameRefReload = os.path.join(_ROOT, 'Thermodynamics', 'RefProfiles', Params.fNameRef[comp])
            if not os.path.isfile(fNameRefReload):
                raise RuntimeError(f'CALC_NEW_REF is set to False, but a reference profile for {comp} ' +
                                   'was not found. Try running again with CALC_NEW_REF set to True in configPP.py.')
            with open(fNameRefReload) as f:
                _ = f.readline()
                Params.wRef_ppt[comp] = np.array(f.readline().split('=')[-1].split(',')).astype(np.float_)
            PrhoRef = np.loadtxt(fNameRefReload, skiprows=3, unpack=False)
            PrhoRef = PrhoRef.T
            Params.nRef[comp] = np.size(Params.wRef_ppt[comp])
            Params.nRefPts[comp] = np.shape(PrhoRef)[1]

            Params.Pref_MPa[comp] = PrhoRef[0,:]
            Params.rhoRef_kgm3[comp] = PrhoRef[1:,:]
            newRef[comp] = False

    return Params
