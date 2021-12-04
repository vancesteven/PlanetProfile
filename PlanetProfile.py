# Import necessary Python modules
import sys
import numpy as np
import importlib
from distutils.util import strtobool
from os.path import isfile

# Import all function definitions for this file
from Utilities.SetupInit import SetupInit, SetupFilenames
from Thermodynamics.LayerPropagators import IceLayers, OceanLayers, InnerLayers, PlanetDepths

#from Thermodynamics.FromLiterature.conductiveMantleTemperature import conductiveMantleTemperature
#from Thermodynamics.FromLiterature.ConvectionDeschampsSotin2001 import ConvectionDeschampsSotin2001


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
        Planet = IceLayers(Planet)
        Planet = OceanLayers(Planet)
        Planet = PlanetDepths(Planet)
        Planet = InnerLayers(Planet)

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
        f.write('  deltaP = ' + str(Planet.Bulk.deltaP) + '\n')
        f.write('  C2mean = ' + str(Planet.C2mean) + '\n')
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
        f.write(' '.join(['z (km)'.ljust(24), \
                           'r (m)'.ljust(24), \
                           'P (MPa)'.ljust(24), \
                           'T (K)'.ljust(24), \
                           'phase ID'.ljust(8), \
                           'rho (kg/m3)'.ljust(24), \
                           'Cp (J/kg/K)'.ljust(24), \
                           'sigma (S/m)'.ljust(24), \
                           'g (m/s2)'.ljust(24), \
                           'vFluid (km/s)']) + '\n')
        print('WARNING: Reload data write-out only includes Steps up to nHydroMax but should be nTotal.' \
             +'This is a placeholder while LayerPropagators is being developed.')
        for i in range(Planet.Steps.nHydroMax):
            line = \
                '{:24.17e} '.format(Planet.z_m[i]/1e3) + \
                '{:24.17e} '.format(Planet.r_m[i]) + \
                '{:24.17e} '.format(Planet.P_MPa[i]) + \
                '{:24.17e} '.format(Planet.T_K[i]) + \
                '{:8d} '.format(Planet.phase[i]) + \
                '{:24.17e} '.format(Planet.rho_kgm3[i]) + \
                '{:24.17e} '.format(Planet.Cp_JkgK[i]) + \
                '{:24.17e} '.format(Planet.sigma_Sm[i]) + \
                '{:24.17e} '.format(Planet.g_ms2[i]) + \
                '{:24.17e}\n'.format(Planet.vFluid_kms[i])
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
        Planet.Pb_MPa, Planet.PbI_MPa, Planet.Bulk.deltaP, Planet.C2mean, Planet.Ocean.QfromMantle_Wm2, \
        Planet.Sil.phiRockMax, Planet.Sil.RsilMean_m, Planet.Core.RFeMean_m, Planet.Sil.RsilRange_m, \
        Planet.Core.RFeRange_m, Planet.Core.rhoFe_kgm3 \
            = (float(f.readline().split('=')[-1]) for _ in range(15))
        # Get integer values from header (nSteps values)
        Planet.Steps.nClath, Planet.Steps.nIceI, \
        Planet.Steps.nIceIIILitho, Planet.Steps.nIceVLitho, \
        Planet.Steps.nHydro, Planet.Steps.nSil, Planet.Steps.nCore \
            = (int(f.readline().split('=')[-1]) for _ in range(7))

    Planet.Steps.nTotal = Planet.Steps.nHydro + Planet.Steps.nSil + Planet.Steps.nCore

    # Read in columnar data that follows header lines -- ocean
    Planet.z_m, Planet.r_m, Planet.P_MPa, Planet.T_K, Planet.phase, Planet.rho_kgm3, Planet.Cp_JkgK, \
    Planet.sigma_Sm, Planet.g_ms2, Planet.vFluid_kms \
        = np.loadtxt(Params.dataFiles.saveFile, skiprows=Params.nHeadLines, unpack=True)
    Planet.z_m *= 1e3  # Stored to disk as km
    Planet.phase = Planet.phase.astype(np.int_)

    # Read in in data for core/mantle trade
    Planet.Sil.RsilTrade_m, Planet.Core.RFeTrade_m, Planet.Sil.rhoSilTrade_kgm3, \
        = np.loadtxt(Params.dataFiles.mantCoreFile, skiprows=1, unpack=True)

    # Read in data for mantle permeability properties
    Planet.Sil.perm1, Planet.Sil.perm2, Planet.Sil.perm3, Planet.Sil.perm4, \
    Planet.Sil.perm5 \
        = np.loadtxt(Params.dataFiles.permFile, skiprows=1, unpack=True)

    return Planet, Params


if __name__ == '__main__': main()
