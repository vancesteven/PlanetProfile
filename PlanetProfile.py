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
    with open(Params.dataFiles.saveFile,'w') as f:
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
        f.write('  RsilMean_m = ' + str(Planet.Sil.RsilMean_m) + '\n')
        f.write('  RFeMean_m = ' + str(Planet.Core.RFeMean_m) + '\n')
        f.write('  RsilRange_m = ' + str(Planet.Sil.RsilRange_m) + '\n')
        f.write('  RFeRange_m = ' + str(Planet.Core.RFeRange_m) + '\n')
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
                          'VP (km/s)'.ljust(24), \
                          'VS (km/s)'.ljust(24), \
                          'QS'.ljust(24), \
                          'KS (GPa)'.ljust(24), \
                          'GS (GPa)'.ljust(24), \
                          'sigma (S/m)']) + '\n')
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
                '{:24.17e} '.format(Planet.Seismic.VP_kms[i]) + \
                '{:24.17e} '.format(Planet.Seismic.VS_kms[i]) + \
                '{:24.17e} '.format(Planet.Seismic.QS[i]) + \
                '{:24.17e} '.format(Planet.Seismic.KS_GPa[i]) + \
                '{:24.17e} '.format(Planet.Seismic.GS_GPa[i]) + \
                '{:24.17e}\n'.format(Planet.sigma_Sm[i])
            f.write(line)

    print('Profile saved to file: ' + Params.dataFiles.saveFile)
    return


def ReloadProfile(Planet, Params):
    """ Reload previously saved PlanetProfile run from disk """
    Params.dataFiles, Params.figureFiles = SetupFilenames(Planet, Params)
    print('Reloading previously saved run from file: ' + Params.dataFiles.saveFile)
    if Params.VERBOSE: print('WARNING: Steps.n settings from PP' + Planet.name + '.py will be ignored.')
    if not isfile(Params.dataFiles.saveFile):
        raise ValueError('Params.CALC_NEW is set to False in config.py but the reload file was not found.\n' \
                        +'Re-run with CALC_NEW set to True to generate the profile.')

    with open(Params.dataFiles.saveFile) as f:
        # Get number of header lines to read in from (and skip for columnar data)
        Params.nHeadLines = int(f.readline().split('=')[-1])
        # Get dissolved salt supposed for ocean (present in filename, but this is intended for future-proofing when we move to a database lookup)
        Planet.Ocean.comp = f.readline().split('=')[-1].strip()
        # Get whether iron core is modeled
        Planet.Do.Fe_CORE = bool(strtobool(f.readline().split('=')[-1].strip()))
        # Get float values from header
        Planet.Ocean.wOcean_ppt, Planet.Bulk.Tb_K, Planet.zb_km, Planet.zClath_m, \
        Planet.Pb_MPa, Planet.PbI_MPa, Planet.Ocean.deltaP, Planet.CMR2mean, Planet.Ocean.QfromMantle_Wm2, \
        Planet.Sil.phiRockMax, Planet.Sil.RsilMean_m, Planet.Core.RFeMean_m, Planet.Sil.RsilRange_m, \
        Planet.Core.RFeRange_m, Planet.Core.rhoFe_kgm3 \
            = (float(f.readline().split('=')[-1]) for _ in range(15))
        # Get integer values from header (nSteps values)
        Planet.Steps.nClath, Planet.Steps.nIceI, \
        Planet.Steps.nIceIIILitho, Planet.Steps.nIceVLitho, \
        Planet.Steps.nHydro, Planet.Steps.nSil, Planet.Steps.nCore \
            = (int(f.readline().split('=')[-1]) for _ in range(7))

    Planet.Steps.nTotal = Planet.Steps.nHydro + Planet.Steps.nSil + Planet.Steps.nCore
    # Read in columnar data that follows header lines -- full-body
    Planet.P_MPa, Planet.T_K, Planet.r_m, Planet.phase, Planet.rho_kgm3, Planet.Cp_JkgK, Planet.alpha_pK, \
    Planet.g_ms2, Planet.phi_frac, Planet.Seismic.VP_kms, Planet.Seismic.VS_kms, Planet.Seismic.QS, \
    Planet.Seismic.KS_GPa, Planet.Seismic.GS_GPa, Planet.sigma_Sm \
        = np.loadtxt(Params.dataFiles.saveFile, skiprows=Params.nHeadLines, unpack=True)
    Planet.z_m = Planet.Bulk.R_m - Planet.r_m
    Planet.phase = Planet.phase.astype(np.int_)

    # Read in data for core/mantle trade
    Planet.Sil.RsilTrade_m, Planet.Core.RFeTrade_m, Planet.Sil.rhoSilTrade_kgm3, \
        = np.loadtxt(Params.dataFiles.mantCoreFile, skiprows=1, unpack=True)

    # Read in data for mantle permeability properties
    Planet.Sil.perm1, Planet.Sil.perm2, Planet.Sil.perm3, Planet.Sil.perm4, \
    Planet.Sil.perm5 \
        = np.loadtxt(Params.dataFiles.permFile, skiprows=1, unpack=True)

    return Planet, Params


def CompareProfile(Planet, Params, fname2):
    """ Checks saved data in a named file against the current Planet run
        to say whether the contents of the files differ.

        Args:
             Planet (PlanetStruct): The current run's completed data
             Params (ParamsStruct): Params for the current run
             fname2 (string): The filename for save data to compare against
        Returns:
            None
    """
    import copy.deepcopy as copy

    Params2 = copy(Params)
    Params2.dataFiles.saveFile = fname2

    Planet2, _ = ReloadProfile(Planet, Params2)

    # Compare header info
    same_nHeadLines = Params.nHeadLines == Params2.nHeadLines
    same_comp = Planet.Ocean.comp == Planet2.Ocean.comp
    same_Fe_CORE = Planet.Do.Fe_CORE and Planet2.Do.Fe_CORE
    same_wOcean_ppt = Planet.Ocean.wOcean_ppt == Planet2.Ocean.wOcean_ppt
    same_Tb_K = Planet.Bulk.Tb_K == Planet2.Bulk.Tb_K
    same_zb_km = Planet.zb_km == Planet2.zb_km
    same_zClath_m = Planet.zClath_m == Planet2.zClath_m
    same_Pb_MPa = Planet.Pb_MPa == Planet2.Pb_MPa
    same_PbI_MPa = Planet.PbI_MPa == Planet2.PbI_MPa
    same_deltaP = Planet.Ocean.deltaP == Planet2.Ocean.deltaP
    same_CMR2mean = Planet.CMR2mean == Planet2.CMR2mean
    same_QfromMantle_Wm2 = Planet.Ocean.QfromMantle_Wm2 == Planet2.Ocean.QfromMantle_Wm2
    same_phiRockMax = Planet.Sil.phiRockMax == Planet2.Sil.phiRockMax
    same_RsilMean_m = Planet.Sil.RsilMean_m == Planet2.Sil.RsilMean_m
    same_RFeMean_m = Planet.Core.RFeMean_m == Planet2.Core.RFeMean_m
    same_RsilRange_m = Planet.Sil.RsilRange_m == Planet2.Sil.RsilRange_m
    same_RFeRange_m = Planet.Core.RFeRange_m == Planet2.Core.RFeRange_m
    same_rhoFe_kgm3 = Planet.Core.rhoFe_kgm3 == Planet2.Core.rhoFe_kgm3
    same_nClath = Planet.Steps.nClath == Planet2.Steps.nClath
    same_nIceI = Planet.Steps.nIceI == Planet2.Steps.nIceI
    same_nIceIIILitho = Planet.Steps.nIceIIILitho == Planet2.Steps.nIceIIILitho
    same_nIceVLitho = Planet.Steps.nIceVLitho == Planet2.Steps.nIceVLitho
    same_nHydro = Planet.Steps.nHydro == Planet2.Steps.nHydro
    same_nSil = Planet.Steps.nSil == Planet2.Steps.nSil
    same_nCore = Planet.Steps.nCore == Planet2.Steps.nCore

    if not same_nHeadLines:
        print('nHeadLines differs. ' + str(Params.nHeadLines) + ' | ' + str(Params2.nHeadLines))

    return


if __name__ == '__main__': main()
