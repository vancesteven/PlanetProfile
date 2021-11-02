# Import necessary Python modules
import sys
import numpy as np
import importlib
from distutils.util import strtobool

# Import all function definitions for this file
from Utilities.SetupInit import SetupInit, SetupFilenames
from Thermodynamics.LayerPropagators import IceLayers, OceanLayers, InnerLayers, PlanetDepths

#from Thermodynamics.FromLiterature.conductiveMantleTemperature import conductiveMantleTemperature
#from Thermodynamics.FromLiterature.ConvectionDeschampsSotin2001 import ConvectionDeschampsSotin2001


"""
Naming Conventions
CamelCase with Capital at start for FunctionNames
ALLCAPS for SWITCHES or USER_ORIENTED_README_FILES
camelCase with lowercase at start for variableNames or filesContainingOnlyVariables
snake_case for output_file_names
FilesContainingFunctions named the same as FunctionName
"""


"""
Step 1:
setup function that check version numbers and makes sure the version is compatible with
seafreeze version, leave configs as is for now, tauP leave for now as well

not indexing over iT in python
separate temps are separate PlanetProfile calls

porosity can be later

Step 2:
catch all function or file for default settings/setup
check if need porosity-> class fields turned on or off

Step 3
check fields function? or create class?

Class full of optional checks that check to see if fields are nonzero (i.e. want to do those calculations)
checks not of existence but if nonzero (ie if not None)
hardcode all fields that expect planets to have
set fields to None as default- PPBody file
if POROUS_ROCK is None- don't run the POROUS_ROCK code

Step 4 : Establish EOS
create plots funciton not needed in python- plots drawn to file instead
can skip plot making and implement printing only

EOS state switch- swEOS_chooser-> chooses between which dissolved salt you are looking at

establish EOS/ OceanEOS that gives property need-> SetupEOS
pass planet struct to EOS, have it return or grab needed EOS for core, water, ice, etc.
if no core, core EOS set to none
flag for if want iron core or not


Step 5: Begin creating outputs
those named with strings/ Str
output function
SetupFilenames


setup step numbers somewhere else, not necessarily in icelayer
Step 6: IceLayer
funtion that itereatively sets up the thermal profile
establishes surface temp and Pressure, then find the temp at pressure at each step
specific heat
clathrate function: ClathrateLayer- function within IceLayer, called as opetion when clathrate flag is set to yes
iceIlayer function setup
IceIIILayer


step7:
OceanLayer function
very important
inputs: number of steps, presure at bottom of icelayer
each layer; cal P plus delta p
check that p seafloor is big but not too big
phase vector: keeps track of what ice type ->outputs global pressure, temp. density,
returns index of which ice types we are
outputs are what are written to disk
else if related to stability of EOS model





PLANET PROFILE MAIN BODY
SetupClathrates  190-205
SetupLayers


"""


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
    header = 'Header lines'
    with open(Params.dataFiles.saveFile,'w') as f:
        f.write(header+'\n')
        for i in range(Planet.Steps.nTotal):
            #f.write('\t'.join( [ '{:3.5e}'.format(val) for val in line] )+'\n')
            f.write('Test\n')

    print('Profile saved to file: ' + Params.dataFiles.saveFile)
    return


def ReloadProfile(Planet, Params):
    """ Reload previously saved PlanetProfile run from disk """
    Params.dataFiles, Params.figureFiles = SetupFilenames(Planet, Params)
    print('Reloading previously saved run from file: ' + Params.dataFiles.saveFile)
    if Params.VERBOSE: print('WARNING: nSteps settings from PP' + Planet.name + '.py will be ignored.')

    with open(Params.dataFiles.saveFile) as f:
        # Get number of header lines to read in from (and skip for columnar data)
        Params.nHeadLines = int(f.readline().split('=')[-1])
        # Get dissolved salt supposed for ocean (present in filename, but this is intended for future-proofing when we move to a database lookup)
        Planet.Ocean.comp = f.readline().split('=')[-1].strip()
        # Get whether iron core is modeled
        Planet.Core.FeCore = bool(strtobool(f.readline().split('=')[-1].strip()))
        # Get float values from header
        Planet.Ocean.wtOcean_ppt, Planet.Tb_K, Planet.Zb_km, Planet.zClath_m, Planet.Pb_MPa, \
        Planet.PbI_MPa, Planet.deltaP, Planet.C2mean, Planet.Sil.Qmantle_Wm2, \
        Planet.Sil.Qb, Planet.Sil.Q_Wm2, Planet.Sil.R_sil_mean_m, \
        Planet.Core.R_Fe_mean_m, Planet.Sil.R_sil_range_m, Planet.Core.R_Fe_range_m, \
        Planet.rho_kgm3 \
            = (float(f.readline().split('=')[-1]) for _ in range(16))
        # Get integer values from header (nSteps values)
        Planet.nStepsIceI, Planet.nStepsOcean, Planet.nStepsHydro, Planet.nStepsTotal, \
        Planet.indSil, Planet.nStepsMantle, Planet.nStepsCore, \
        Planet.nIceIIILitho, Planet.nIceVLitho \
            = (int(f.readline().split('=')[-1]) for _ in range(9))
    Planet.nStepsClath = Planet.nStepsHydro - Planet.nStepsIceI - Planet.nStepsOcean - Planet.nIceIIILitho - Planet.nIceVLitho

    # Read in columnar data that follows header lines -- ocean
    Planet.z_m, Planet.r_m, Planet.P_MPa, Planet.T_K, Planet.phase, Planet.rho_kgm3, Planet.g_ms2, \
    Planet.Cp, Planet.vfluid_kms \
        = np.loadtxt(Params.dataFiles.saveFile, skiprows=Params.nHeadLines, unpack=True)
    Planet.z_m *= 1e3  # Stored as km
    Planet.phase = Planet.phase.astype(np.int_)

    # Read in in data for core/mantle trade
    Planet.Sil.R_sil_trade_m, Planet.Core.R_Fe_trade_m, Planet.Sil.rho_sil_trade_kgm3, \
        = np.loadtxt(Params.dataFiles.mantCoreFile, skiprows=1, unpack=True)

    # Read in data for mantle physical properties
    Planet.Sil.r_mantle_m, Planet.Sil.P_mantle_MPa, Planet.Sil.porosDens, \
    Planet.Sil.perm1, Planet.Sil.perm2, Planet.Sil.perm3, Planet.Sil.perm4, \
    Planet.Sil.perm5 \
        = np.loadtxt(Params.dataFiles.mantPropsFile, skiprows=1, unpack=True)

    # Read in full-body layer values
    # Currently in PlanetProfile.m these values are labeled differently from those above,
    # but I think these are the ones we will ultimately want to keep & use so I'm overwriting
    Planet.r_m, Planet.P_MPa, Planet.T_K, Planet.rho_kgm3, Planet.sig_Sm, Planet.phase, \
        Planet.g_ms2 \
        = np.loadtxt(Params.dataFiles.fullLayersFile, skiprows=1, unpack=True)
    Planet.phase = Planet.phase.astype(np.int_)

    return Planet, Params


if __name__ == '__main__': main()
