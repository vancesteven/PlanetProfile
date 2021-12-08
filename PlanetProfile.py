# Import necessary Python modules
import sys
import numpy as np
import importlib
from distutils.util import strtobool
from os.path import isfile

# Import all function definitions for this file
from Utilities.SetupInit import SetupInit, SetupFilenames
from Thermodynamics.LayerPropagators import IceLayers, OceanLayers, HydroConvect, InnerLayers
from Thermodynamics.Electrical import ElecConduct
from Seismic import SeismicCalcs

#from Thermodynamics.FromLiterature.conductiveMantleTemperature import conductiveMantleTemperature


""" MAIN RUN BLOCK """
def main():

    # Command line args
    if len(sys.argv) > 1:
        # Body name was passed as command line argument
        bodyname = sys.argv[1]
    else:
        # No command line argument, ask user which body to run
        bodyname = input('Please input body name: ')
        if bodyname == '':
            print('No body name entered. Defaulting to Europa.')
            bodyname = 'Europa'

    bodyname = bodyname.capitalize()
    body = importlib.import_module(bodyname+'.PP'+bodyname)
    Planet = body.Planet
    Params = body.Params
    if Params.VERBOSE: print('Body name: ' + Planet.name)

    if Params.DEBUG:
        # Compare to test Matlab output
        Planet, Params = LoadBody(bodyname)
        CompareProfile(Planet, Params, 'Europa/EuropaProfile_Seawater_0WtPpt_Tb269.800K_mat.txt')
    else:
        Planet = PlanetProfile(Planet, Params)

    return

""" END MAIN RUN BLOCK """


def PlanetProfile(Planet, Params):

    if Params.CALC_NEW:
        # Initialize
        Planet, Params = SetupInit(Planet, Params)
        Planet = IceLayers(Planet, Params)
        Planet = OceanLayers(Planet, Params)
        Planet = HydroConvect(Planet, Params)
        Planet = InnerLayers(Planet, Params)
        Planet = ElecConduct(Planet, Params)
        Planet = SeismicCalcs(Planet, Params)

        # Save data after modeling
        WriteProfile(Planet, Params)
    else:
        # Reload previous run
        Planet, Params = ReloadProfile(Planet, Params)

    if not Params.SKIP_PROFILES:
        # Plotting functions
        pass

    if Params.VERBOSE: print('Run complete!')
    return Planet


def WriteProfile(Planet, Params):
    """ Write out all profile calculations to disk """
    Params.nHeadLines = 26  # Increment as new header lines are added
    with open(Params.DataFiles.saveFile,'w') as f:
        # Print number of header lines first so we can skip the rest on read-in if we want to
        f.write('  nHeadLines = ' + str(Params.nHeadLines) + '\n')
        f.write('  Ocean salt = ' + str(Planet.Ocean.comp) + '\n')
        f.write('  Iron core = ' + str(Planet.Do.Fe_CORE) + '\n')
        f.write('  Salinity(ppt) = ' + str(Planet.Ocean.wOcean_ppt) + '\n')
        f.write('  Tb_K = ' + str(Planet.Bulk.Tb_K) + '\n')
        f.write('  zb_km = ' + str(Planet.zb_km) + '\n')
        f.write('  zClath_m = ' + str(Planet.zClath_m) + '\n')
        f.write('  Pb_MPa = ' + str(Planet.Pb_MPa) + '\n')
        f.write('  PbI_MPa = ' + str(Planet.PbI_MPa) + '\n')
        f.write('  deltaP = ' + str(Planet.Ocean.deltaP) + '\n')
        f.write('  CMR2mean = ' + str(Planet.CMR2mean) + '\n')
        f.write('  QfromMantle_Wm2 = ' + str(Planet.Ocean.QfromMantle_Wm2) + '\n')
        f.write('  phiRockMax = ' + str(Planet.Sil.phiRockMax) + '\n')
        f.write('  RsilMean_m = ' + str(Planet.Sil.Rmean_m) + '\n')
        f.write('  RFeMean_m = ' + str(Planet.Core.Rmean_m) + '\n')
        f.write('  RsilRange_m = ' + str(Planet.Sil.Rrange_m) + '\n')
        f.write('  RFeRange_m = ' + str(Planet.Core.Rrange_m) + '\n')
        f.write('  rhoFe_kgm3 = ' + str(Planet.Core.rhoFe_kgm3) + '\n')
        f.write('  Steps.nClath = ' + str(Planet.Steps.nClath) + '\n')
        f.write('  Steps.nIceI = ' + str(Planet.Steps.nIceI) + '\n')
        f.write('  Steps.nIceIIILitho = ' + str(Planet.Steps.nIceIIILitho) + '\n')
        f.write('  Steps.nIceVLitho = ' + str(Planet.Steps.nIceVLitho) + '\n')
        f.write('  Steps.nHydro = ' + str(Planet.Steps.nHydro) + '\n')
        f.write('  Steps.nSil = ' + str(Planet.Steps.nSil) + '\n')
        f.write('  Steps.nCore = ' + str(Planet.Steps.nCore) + '\n')
        f.write(' '.join(['P (MPa)'.ljust(24), \
                          'T (K)'.ljust(24), \
                          'r (m)'.ljust(24), \
                          'phase ID'.ljust(8), \
                          'rho (kg/m3)'.ljust(24), \
                          'Cp (J/kg/K)'.ljust(24), \
                          'alpha (1/K)'.ljust(24), \
                          'g (m/s2)'.ljust(24), \
                          'phi (void/solid frac)'.ljust(24), \
                          'sigma (S/m)'.ljust(24), \
                          'VP (km/s)'.ljust(24), \
                          'VS (km/s)'.ljust(24), \
                          'QS'.ljust(24), \
                          'KS (GPa)'.ljust(24), \
                          'GS (GPa)']) + '\n')
        # Now print the columnar data
        for i in range(Planet.Steps.nTotal):
            line = \
                '{:24.17e} '.format(Planet.P_MPa[i]) + \
                '{:24.17e} '.format(Planet.T_K[i]) + \
                '{:24.17e} '.format(Planet.r_m[i]) + \
                '{:8d} '.format(Planet.phase[i]) + \
                '{:24.17e} '.format(Planet.rho_kgm3[i]) + \
                '{:24.17e} '.format(Planet.Cp_JkgK[i]) + \
                '{:24.17e} '.format(Planet.alpha_pK[i]) + \
                '{:24.17e} '.format(Planet.g_ms2[i]) + \
                '{:24.17e} '.format(Planet.phi_frac[i]) + \
                '{:24.17e} '.format(Planet.sigma_Sm[i]) + \
                '{:24.17e} '.format(Planet.Seismic.VP_kms[i]) + \
                '{:24.17e} '.format(Planet.Seismic.VS_kms[i]) + \
                '{:24.17e} '.format(Planet.Seismic.QS[i]) + \
                '{:24.17e} '.format(Planet.Seismic.KS_GPa[i]) + \
                '{:24.17e}\n '.format(Planet.Seismic.GS_GPa[i])
            f.write(line)

    print('Profile saved to file: ' + Params.DataFiles.saveFile)
    return


def ReloadProfile(Planet, Params, fnameOverride=None):
    """ Reload previously saved PlanetProfile run from disk """

    if fnameOverride is not None:
        Params.DataFiles.saveFile = fnameOverride
    else:
        Params.DataFiles, Params.FigureFiles = SetupFilenames(Planet, Params)
    print('Reloading previously saved run from file: ' + Params.DataFiles.saveFile)
    if Params.VERBOSE: print('WARNING: Steps.n settings from PP' + Planet.name + '.py will be ignored.')
    if not isfile(Params.DataFiles.saveFile):
        raise ValueError('Params.CALC_NEW is set to False in config.py but the reload file was not found.\n' \
                        +'Re-run with CALC_NEW set to True to generate the profile.')

    with open(Params.DataFiles.saveFile) as f:
        # Get number of header lines to read in from (and skip for columnar data)
        Params.nHeadLines = int(f.readline().split('=')[-1])
        # Get dissolved salt supposed for ocean (present in filename, but this is intended for future-proofing when we move to a database lookup)
        Planet.Ocean.comp = f.readline().split('=')[-1].strip()
        # Get whether iron core is modeled
        Planet.Do.Fe_CORE = bool(strtobool(f.readline().split('=')[-1].strip()))
        # Get float values from header
        Planet.Ocean.wOcean_ppt, Planet.Bulk.Tb_K, Planet.zb_km, Planet.zClath_m, \
        Planet.Pb_MPa, Planet.PbI_MPa, Planet.Ocean.deltaP, Planet.CMR2mean, Planet.Ocean.QfromMantle_Wm2, \
        Planet.Sil.phiRockMax, Planet.Sil.Rmean_m, Planet.Core.Rmean_m, Planet.Sil.Rrange_m, \
        Planet.Core.Rrange_m, Planet.Core.rhoFe_kgm3 \
            = (float(f.readline().split('=')[-1]) for _ in range(15))
        # Get integer values from header (nSteps values)
        Planet.Steps.nClath, Planet.Steps.nIceI, \
        Planet.Steps.nIceIIILitho, Planet.Steps.nIceVLitho, \
        Planet.Steps.nHydro, Planet.Steps.nSil, Planet.Steps.nCore \
            = (int(f.readline().split('=')[-1]) for _ in range(7))

    Planet.Steps.nTotal = Planet.Steps.nHydro + Planet.Steps.nSil + Planet.Steps.nCore
    # Read in columnar data that follows header lines -- full-body
    Planet.P_MPa, Planet.T_K, Planet.r_m, Planet.phase, Planet.rho_kgm3, Planet.Cp_JkgK, Planet.alpha_pK, \
    Planet.g_ms2, Planet.phi_frac, Planet.sigma_Sm, Planet.Seismic.VP_kms, Planet.Seismic.VS_kms, Planet.Seismic.QS, \
    Planet.Seismic.KS_GPa, Planet.Seismic.GS_GPa \
        = np.loadtxt(Params.DataFiles.saveFile, skiprows=Params.nHeadLines, unpack=True)
    Planet.z_m = Planet.Bulk.R_m - Planet.r_m
    Planet.phase = Planet.phase.astype(np.int_)

    # Read in data for core/mantle trade
    Planet.Sil.Rtrade_m, Planet.Core.Rtrade_m, Planet.Sil.rhoTrade_kgm3, \
        = np.loadtxt(Params.DataFiles.mantCoreFile, skiprows=1, unpack=True)

    # Read in data for mantle permeability properties
    Planet.Sil.perm1, Planet.Sil.perm2, Planet.Sil.perm3, Planet.Sil.perm4, \
    Planet.Sil.perm5 \
        = np.loadtxt(Params.DataFiles.permFile, skiprows=1, unpack=True)

    return Planet, Params


def LoadBody(bodyname, fnameOverride=None):
    """ Loads the settings in PPBody.py to reload a previous run. """
    # NOTE: Seems to overwrite all existing Planet instances with the newly read in values. Do not use to load in a second profile!
    # Use copy.deepcopy and ReloadProfile instead.
    bodyname = bodyname.capitalize()
    body = importlib.import_module(bodyname + '.PP' + bodyname)
    Planet = body.Planet
    Params = body.Params
    Planet, Params = ReloadProfile(Planet, Params, fnameOverride=fnameOverride)

    return Planet, Params


def CompareProfile(Planet, Params, fname2, tol=0.01, tiny=1e-6):
    """ Checks saved data in a named file against the current Planet run
        to say whether the contents of the files differ.

        Args:
             Planet (PlanetStruct): The current run's completed data
             Params (ParamsStruct): Params for the current run
             fname2 (string): The filename for save data to compare against
             tol (float): Tolerance in % difference to consider values the same
                in comparing layer arrays
        Returns:
            None
    """
    import copy

    print('Comparing current run with ' + fname2 + '...')
    Planet2 = copy.deepcopy(Planet)
    Params2 = copy.deepcopy(Params)
    Planet2, Params2 = ReloadProfile(Planet2, Params2, fnameOverride=fname2)

    # Avoid divide-by-zero errors
    if Planet.Ocean.wOcean_ppt == 0: Planet.Ocean.wOcean_ppt = tiny
    if Planet.Bulk.Tb_K == 0: Planet.Bulk.Tb_K = tiny
    if Planet.zb_km == 0: Planet.zb_km = tiny
    if Planet.zClath_m == 0: Planet.zClath_m = tiny
    if Planet.Pb_MPa == 0: Planet.Pb_MPa = tiny
    if Planet.PbI_MPa == 0: Planet.PbI_MPa = tiny
    if Planet.Ocean.deltaP == 0: Planet.Ocean.deltaP = tiny
    if Planet.CMR2mean == 0: Planet.CMR2mean = tiny
    if Planet.Ocean.QfromMantle_Wm2 == 0: Planet.Ocean.QfromMantle_Wm2 = tiny
    if Planet.Sil.phiRockMax == 0: Planet.Sil.phiRockMax = tiny
    if Planet.Sil.Rmean_m == 0: Planet.Sil.Rmean_m = tiny
    if Planet.Core.Rmean_m == 0: Planet.Core.Rmean_m = tiny
    if Planet.Sil.Rrange_m == 0: Planet.Sil.Rrange_m = tiny
    if Planet.Core.Rrange_m == 0: Planet.Core.Rrange_m = tiny
    if Planet.Core.rhoFe_kgm3 == 0: Planet.Core.rhoFe_kgm3 = tiny

    # Avoid false positives for zero values
    if Planet2.Ocean.wOcean_ppt == 0: Planet2.Ocean.wOcean_ppt = tiny
    if Planet2.Bulk.Tb_K == 0: Planet2.Bulk.Tb_K = tiny
    if Planet2.zb_km == 0: Planet2.zb_km = tiny
    if Planet2.zClath_m == 0: Planet2.zClath_m = tiny
    if Planet2.Pb_MPa == 0: Planet2.Pb_MPa = tiny
    if Planet2.PbI_MPa == 0: Planet2.PbI_MPa = tiny
    if Planet2.Ocean.deltaP == 0: Planet2.Ocean.deltaP = tiny
    if Planet2.CMR2mean == 0: Planet2.CMR2mean = tiny
    if Planet2.Ocean.QfromMantle_Wm2 == 0: Planet2.Ocean.QfromMantle_Wm2 = tiny
    if Planet2.Sil.phiRockMax == 0: Planet2.Sil.phiRockMax = tiny
    if Planet2.Sil.Rmean_m == 0: Planet2.Sil.Rmean_m = tiny
    if Planet2.Core.Rmean_m == 0: Planet2.Core.Rmean_m = tiny
    if Planet2.Sil.Rrange_m == 0: Planet2.Sil.Rrange_m = tiny
    if Planet2.Core.Rrange_m == 0: Planet2.Core.Rrange_m = tiny
    if Planet2.Core.rhoFe_kgm3 == 0: Planet2.Core.rhoFe_kgm3 = tiny

    # Compare header info
    same_nHeadLines = Params.nHeadLines == Params2.nHeadLines
    same_comp = Planet.Ocean.comp == Planet2.Ocean.comp
    same_Fe_CORE = Planet.Do.Fe_CORE and Planet2.Do.Fe_CORE
    same_wOcean_ppt = (Planet.Ocean.wOcean_ppt - Planet2.Ocean.wOcean_ppt) / Planet.Ocean.wOcean_ppt < tol
    same_Tb_K = (Planet.Bulk.Tb_K - Planet2.Bulk.Tb_K) / Planet.Bulk.Tb_K < tol
    same_zb_km = (Planet.zb_km > Planet2.zb_km) / Planet.zb_km < tol
    same_zClath_m = (Planet.zClath_m - Planet2.zClath_m) / Planet.zClath_m < tol
    same_Pb_MPa = (Planet.Pb_MPa - Planet2.Pb_MPa) / Planet.Pb_MPa < tol
    same_PbI_MPa = (Planet.PbI_MPa - Planet2.PbI_MPa) / Planet.PbI_MPa < tol
    same_deltaP = (Planet.Ocean.deltaP - Planet2.Ocean.deltaP) / Planet.Ocean.deltaP < tol
    same_CMR2mean = (Planet.CMR2mean - Planet2.CMR2mean) / Planet.CMR2mean < tol
    same_QfromMantle_Wm2 = (Planet.Ocean.QfromMantle_Wm2 - Planet2.Ocean.QfromMantle_Wm2) / Planet.Ocean.QfromMantle_Wm2 < tol
    same_phiRockMax = (Planet.Sil.phiRockMax - Planet2.Sil.phiRockMax) / Planet.Sil.phiRockMax < tol
    same_silRmean_m = (Planet.Sil.Rmean_m - Planet2.Sil.Rmean_m) / Planet.Sil.Rmean_m < tol
    same_coreRmean_m = (Planet.Core.Rmean_m - Planet2.Core.Rmean_m) / Planet.Core.Rmean_m < tol
    same_silRrange_m = (Planet.Sil.Rrange_m - Planet2.Sil.Rrange_m) / Planet.Sil.Rrange_m < tol
    same_coreRrange_m = (Planet.Core.Rrange_m - Planet2.Core.Rrange_m) / Planet.Core.Rrange_m < tol
    same_rhoFe_kgm3 = (Planet.Core.rhoFe_kgm3 - Planet2.Core.rhoFe_kgm3) / Planet.Core.rhoFe_kgm3 < tol
    same_nClath = Planet.Steps.nClath == Planet2.Steps.nClath
    same_nIceI = Planet.Steps.nIceI == Planet2.Steps.nIceI
    same_nIceIIILitho = Planet.Steps.nIceIIILitho == Planet2.Steps.nIceIIILitho
    same_nIceVLitho = Planet.Steps.nIceVLitho == Planet2.Steps.nIceVLitho
    same_nHydro = Planet.Steps.nHydro == Planet2.Steps.nHydro
    same_nSil = Planet.Steps.nSil == Planet2.Steps.nSil
    same_nCore = Planet.Steps.nCore == Planet2.Steps.nCore

    same_steps = same_nClath and same_nIceI and same_nIceIIILitho and same_nIceVLitho and same_nHydro and same_nSil \
        and same_nCore
    headers_match = same_nHeadLines and same_comp and same_Fe_CORE and same_wOcean_ppt and same_Tb_K and same_zb_km \
        and same_zClath_m and same_Pb_MPa and same_PbI_MPa and same_deltaP and same_CMR2mean and same_QfromMantle_Wm2 \
        and same_phiRockMax and same_silRmean_m and same_coreRmean_m and same_silRrange_m and same_coreRrange_m \
        and same_rhoFe_kgm3 and same_steps

    if not headers_match:
        if not same_nHeadLines: print('nHeadLines differs. ' + str(Params.nHeadLines) + ' | ' + str(Params2.nHeadLines))
        if not same_comp: print('Ocean.comp differs. ' + str(Planet.Ocean.comp) + ' | ' + str(Planet2.Ocean.comp))
        if not same_Fe_CORE: print('Do.Fe_CORE differs. ' + str(Planet.Do.Fe_CORE) + ' | ' + str(Planet2.Do.Fe_CORE))
        if not same_wOcean_ppt: print('Ocean.wOcean_ppt differs. ' + str(Planet.Ocean.wOcean_ppt) + ' | ' + str(Planet2.Ocean.wOcean_ppt))
        if not same_Tb_K: print('Bulk.Tb_K differs. ' + str(Planet.Bulk.Tb_K) + ' | ' + str(Planet2.Bulk.Tb_K))
        if not same_zb_km: print('zb_km differs. ' + str(Planet.zb_km) + ' | ' + str(Planet2.zb_km))
        if not same_zClath_m: print('zClath_m differs. ' + str(Planet.zClath_m) + ' | ' + str(Planet2.zClath_m))
        if not same_Pb_MPa: print('Pb_MPa differs. ' + str(Planet.Pb_MPa) + ' | ' + str(Planet2.Pb_MPa))
        if not same_PbI_MPa: print('PbI_MPa differs. ' + str(Planet.PbI_MPa) + ' | ' + str(Planet2.PbI_MPa))
        if not same_deltaP: print('Ocean.deltaP differs. ' + str(Planet.Ocean.deltaP) + ' | ' + str(Planet2.Ocean.deltaP))
        if not same_CMR2mean: print('CMR2mean differs. ' + str(Planet.CMR2mean) + ' | ' + str(Planet2.CMR2mean))
        if not same_QfromMantle_Wm2: print('Ocean.QfromMantle_Wm2 differs. ' + str(Planet.Ocean.QfromMantle_Wm2) + ' | ' + str(Planet2.Ocean.QfromMantle_Wm2))
        if not same_phiRockMax: print('Sil.phiRockMax differs. ' + str(Planet.Sil.phiRockMax) + ' | ' + str(Planet2.Sil.phiRockMax))
        if not same_silRmean_m: print('Sil.Rmean_m differs. ' + str(Planet.Sil.Rmean_m) + ' | ' + str(Planet2.Sil.Rmean_m))
        if not same_coreRmean_m: print('Core.Rmean_m differs. ' + str(Planet.Core.Rmean_m) + ' | ' + str(Planet2.Core.Rmean_m))
        if not same_silRrange_m: print('Sil.Rrange_m differs. ' + str(Planet.Sil.Rrange_m) + ' | ' + str(Planet2.Sil.Rrange_m))
        if not same_coreRrange_m: print('Core.Rrange_m differs. ' + str(Planet.Core.Rrange_m) + ' | ' + str(Planet2.Core.Rrange_m))
        if not same_rhoFe_kgm3: print('Core.rhoFe_kgm3 differs. ' + str(Planet.Core.rhoFe_kgm3) + ' | ' + str(Planet2.Core.rhoFe_kgm3))
        if not same_nClath: print('Steps.nClath differs. ' + str(Planet.Steps.nClath) + ' | ' + str(Planet2.Steps.nClath))
        if not same_nIceI: print('Steps.nIceI differs. ' + str(Planet.Steps.nIceI) + ' | ' + str(Planet2.Steps.nIceI))
        if not same_nIceIIILitho: print('Steps.nIceIIILitho differs. ' + str(Planet.Steps.nIceIIILitho) + ' | ' + str(Planet2.Steps.nIceIIILitho))
        if not same_nIceVLitho: print('Steps.nIceVLitho differs. ' + str(Planet.Steps.nIceVLitho) + ' | ' + str(Planet2.Steps.nIceVLitho))
        if not same_nHydro: print('Steps.nHydro differs. ' + str(Planet.Steps.nHydro) + ' | ' + str(Planet2.Steps.nHydro))
        if not same_nSil: print('Steps.nSil differs. ' + str(Planet.Steps.nSil) + ' | ' + str(Planet2.Steps.nSil))
        if not same_nCore: print('Steps.nCore differs. ' + str(Planet.Steps.nCore) + ' | ' + str(Planet2.Steps.nCore))
    else:
        print('All header values match!')

    if same_steps:
        print('Some arrays differ. The first index more than ' + str(round(tol*100,2)) + '% different will be printed.')
        iIO = Planet.Steps.nSurfIce
        iOS = Planet.Steps.nHydro
        iSC = Planet.Steps.nHydro+Planet.Steps.nSil
        iCC = Planet.Steps.nTotal

        # Avoid divide-by-zero errors
        Planet.P_MPa[Planet.P_MPa==0] = tiny
        Planet.T_K[Planet.T_K==0] = tiny
        Planet.r_m[Planet.r_m==0] = tiny
        Planet.rho_kgm3[Planet.rho_kgm3==0] = tiny
        Planet.Cp_JkgK[Planet.Cp_JkgK==0] = tiny
        Planet.alpha_pK[Planet.alpha_pK==0] = tiny
        Planet.g_ms2[Planet.g_ms2==0] = tiny
        Planet.phi_frac[Planet.phi_frac==0] = tiny
        Planet.sigma_Sm[Planet.sigma_Sm==0] = tiny
        Planet.Seismic.VP_kms[Planet.Seismic.VP_kms==0] = tiny
        Planet.Seismic.VS_kms[Planet.Seismic.VS_kms==0] = tiny
        Planet.Seismic.QS[Planet.Seismic.QS==0] = tiny
        Planet.Seismic.KS_GPa[Planet.Seismic.KS_GPa==0] = tiny
        Planet.Seismic.GS_GPa[Planet.Seismic.GS_GPa==0] = tiny

        # Avoid false positives from zero values
        Planet2.P_MPa[Planet2.P_MPa==0] = tiny
        Planet2.T_K[Planet2.T_K==0] = tiny
        Planet2.r_m[Planet2.r_m==0] = tiny
        Planet2.rho_kgm3[Planet2.rho_kgm3==0] = tiny
        Planet2.Cp_JkgK[Planet2.Cp_JkgK==0] = tiny
        Planet2.alpha_pK[Planet2.alpha_pK==0] = tiny
        Planet2.g_ms2[Planet2.g_ms2==0] = tiny
        Planet2.phi_frac[Planet2.phi_frac==0] = tiny
        Planet2.sigma_Sm[Planet2.sigma_Sm==0] = tiny
        Planet2.Seismic.VP_kms[Planet2.Seismic.VP_kms==0] = tiny
        Planet2.Seismic.VS_kms[Planet2.Seismic.VS_kms==0] = tiny
        Planet2.Seismic.QS[Planet2.Seismic.QS==0] = tiny
        Planet2.Seismic.KS_GPa[Planet2.Seismic.KS_GPa==0] = tiny
        Planet2.Seismic.GS_GPa[Planet2.Seismic.GS_GPa==0] = tiny

        # Compare layer calculations
        diff_P_MPa = [i[0] for i, val in np.ndenumerate(Planet.P_MPa) if abs(val-Planet2.P_MPa[i[0]])/val>tol]
        diff_T_K = [i[0] for i, val in np.ndenumerate(Planet.T_K) if abs(val-Planet2.T_K[i[0]])/val>tol]
        diff_r_m = [i[0] for i, val in np.ndenumerate(Planet.r_m) if abs(val-Planet2.r_m[i[0]])/val>tol]
        diff_phase = [i[0] for i, val in np.ndenumerate(Planet.phase) if Planet2.phase[i[0]]!=val]
        diff_rho_kgm3 = [i[0] for i, val in np.ndenumerate(Planet.rho_kgm3) if abs(val-Planet2.rho_kgm3[i[0]])/val>tol]
        diff_Cp_JkgK = [i[0] for i, val in np.ndenumerate(Planet.Cp_JkgK) if abs(val-Planet2.Cp_JkgK[i[0]])/val>tol]
        diff_alpha_pK = [i[0] for i, val in np.ndenumerate(Planet.alpha_pK) if abs(val-Planet2.alpha_pK[i[0]])/val>tol]
        diff_g_ms2 = [i[0] for i, val in np.ndenumerate(Planet.g_ms2) if abs(val-Planet2.g_ms2[i[0]])/val>tol]
        diff_phi_frac = [i[0] for i, val in np.ndenumerate(Planet.phi_frac) if abs(val-Planet2.phi_frac[i[0]])/val>tol]
        diff_sigma_Sm = [i[0] for i, val in np.ndenumerate(Planet.sigma_Sm) if abs(val-Planet2.sigma_Sm[i[0]])/val>tol]
        diff_VP_kms = [i[0] for i, val in np.ndenumerate(Planet.Seismic.VP_kms) if abs(val-Planet2.Seismic.VP_kms[i[0]])/val>tol]
        diff_VS_kms = [i[0] for i, val in np.ndenumerate(Planet.Seismic.VS_kms) if abs(val-Planet2.Seismic.VS_kms[i[0]])/val>tol]
        diff_QS = [i[0] for i, val in np.ndenumerate(Planet.Seismic.QS) if abs(val-Planet2.Seismic.QS[i[0]])/val>tol]
        diff_KS_GPa = [i[0] for i, val in np.ndenumerate(Planet.Seismic.KS_GPa) if abs(val-Planet2.Seismic.KS_GPa[i[0]])/val>tol]
        diff_GS_GPa = [i[0] for i, val in np.ndenumerate(Planet.Seismic.GS_GPa) if abs(val-Planet2.Seismic.GS_GPa[i[0]])/val>tol]

        same_P_MPa = len(diff_P_MPa) == 0
        same_T_K = len(diff_T_K) == 0
        same_r_m = len(diff_r_m) == 0
        same_phase = len(diff_phase) == 0
        same_rho_kgm3 = len(diff_rho_kgm3) == 0
        same_Cp_JkgK = len(diff_Cp_JkgK) == 0
        same_alpha_pK = len(diff_alpha_pK) == 0
        same_g_ms2 = len(diff_g_ms2) == 0
        same_phi_frac = len(diff_phi_frac) == 0
        same_sigma_Sm = len(diff_sigma_Sm) == 0
        same_VP_kms = len(diff_VP_kms) == 0
        same_VS_kms = len(diff_VS_kms) == 0
        same_QS = len(diff_QS) == 0
        same_KS_GPa = len(diff_KS_GPa) == 0
        same_GS_GPa = len(diff_GS_GPa) == 0

        layers_match = same_P_MPa and same_T_K and same_r_m and same_phase and same_rho_kgm3 and same_Cp_JkgK \
            and same_alpha_pK and same_g_ms2 and same_phi_frac and same_sigma_Sm and same_VP_kms and same_VS_kms and \
            same_QS and same_KS_GPa and same_GS_GPa
        all_match = headers_match and layers_match

        if not layers_match:
            if not same_P_MPa: print('P_MPa differs in position ' + str(diff_P_MPa[0]) + ': ' + str(Planet.P_MPa[diff_P_MPa[0]]) + ' | ' + str(Planet2.P_MPa[diff_P_MPa[0]]))
            if not same_T_K: print('T_K differs in position ' + str(diff_T_K[0]) + ': ' + str(Planet.T_K[diff_T_K[0]]) + ' | ' + str(Planet2.T_K[diff_T_K[0]]))
            if not same_r_m: print('r_m differs in position ' + str(diff_r_m[0]) + ': ' + str(Planet.r_m[diff_r_m[0]]) + ' | ' + str(Planet2.r_m[diff_r_m[0]]))
            if not same_phase: print('phase differs in position ' + str(diff_phase[0]) + ': ' + str(Planet.phase[diff_phase[0]]) + ' | ' + str(Planet2.phase[diff_phase[0]]))
            if not same_rho_kgm3: print('rho_kgm3 differs in position ' + str(diff_rho_kgm3[0]) + ': ' + str(Planet.rho_kgm3[diff_rho_kgm3[0]]) + ' | ' + str(Planet2.rho_kgm3[diff_rho_kgm3[0]]))
            if not same_Cp_JkgK: print('Cp_JkgK differs in position ' + str(diff_Cp_JkgK[0]) + ': ' + str(Planet.Cp_JkgK[diff_Cp_JkgK[0]]) + ' | ' + str(Planet2.Cp_JkgK[diff_Cp_JkgK[0]]))
            if not same_alpha_pK: print('alpha_pK differs in position ' + str(diff_alpha_pK[0]) + ': ' + str(Planet.alpha_pK[diff_alpha_pK[0]]) + ' | ' + str(Planet2.alpha_pK[diff_alpha_pK[0]]))
            if not same_g_ms2: print('g_ms2 differs in position ' + str(diff_g_ms2[0]) + ': ' + str(Planet.g_ms2[diff_g_ms2[0]]) + ' | ' + str(Planet2.g_ms2[diff_g_ms2[0]]))
            if not same_phi_frac: print('phi_frac differs in position ' + str(diff_phi_frac[0]) + ': ' + str(Planet.phi_frac[diff_phi_frac[0]]) + ' | ' + str(Planet2.phi_frac[diff_phi_frac[0]]))
            if not same_sigma_Sm: print('sigma_Sm differs in position ' + str(diff_sigma_Sm[0]) + ': ' + str(Planet.sigma_Sm[diff_sigma_Sm[0]]) + ' | ' + str(Planet2.sigma_Sm[diff_sigma_Sm[0]]))
            if not same_VP_kms: print('VP_kms differs in position ' + str(diff_VP_kms[0]) + ': ' + str(Planet.Seismic.VP_kms[diff_VP_kms[0]]) + ' | ' + str(Planet2.Seismic.VP_kms[diff_VP_kms[0]]))
            if not same_VS_kms: print('VS_kms differs in position ' + str(diff_VS_kms[0]) + ': ' + str(Planet.Seismic.VS_kms[diff_VS_kms[0]]) + ' | ' + str(Planet2.Seismic.VS_kms[diff_VS_kms[0]]))
            if not same_QS: print('QS differs in position ' + str(diff_QS[0]) + ': ' + str(Planet.Seismic.QS[diff_QS[0]]) + ' | ' + str(Planet2.Seismic.QS[diff_QS[0]]))
            if not same_KS_GPa: print('KS_GPa differs in position ' + str(diff_KS_GPa[0]) + ': ' + str(Planet.Seismic.KS_GPa[diff_KS_GPa[0]]) + ' | ' + str(Planet2.Seismic.KS_GPa[diff_KS_GPa[0]]))
            if not same_GS_GPa: print('GS_GPa differs in position ' + str(diff_GS_GPa[0]) + ': ' + str(Planet.Seismic.GS_GPa[diff_GS_GPa[0]]) + ' | ' + str(Planet2.Seismic.GS_GPa[diff_GS_GPa[0]]))
        elif not all_match:
            print('All layer calculations match!')
        else:
            print('All values match! Both profiles are identical.')
    else:
        print('Unable to compare layer calculations, as Steps values do not match.')

    return


if __name__ == '__main__': main()
