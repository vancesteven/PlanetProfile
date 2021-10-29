# Import necessary Python modules
import sys
import numpy as np
import importlib

# Import all function definitions for this file
from Utilities.SetupInit import SetupInit, SetupFilenames, LayersStruct
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
setup funciton that check version numbers and makes sure the version is compatible with
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
if PORUS_ROCK is None- don't run the POROUS_ROCK code

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
    Constants = body.Constants
    if Params.VERBOSE: print('Body name: ' + Planet.name)

    Planet = PlanetProfile(Planet, Params, Constants)

    return

""" END MAIN RUN BLOCK """

def PlanetProfile(Planet, Params, Constants):

    if Params.CALC_NEW:
        # Initialize
        Planet, Params, Layers = SetupInit(Planet, Params, Constants)
        Layers = IceLayers(Planet, Layers, Constants)
        Layers = OceanLayers(Planet, Layers, Constants)
        Layers = PlanetDepths(Planet, Layers, Constants)
        Planet, Layers = InnerLayers(Planet, Layers, Constants)
        WriteProfile(Planet, Params, Layers)
    else:
        # Reload previous run
        Planet, Params, Layers = ReloadProfile(Planet, Params)

    if not Params.SKIP_PROFILES:
        # Plotting functions
        pass

    if Params.VERBOSE: print('Run complete!')
    return Planet


def WriteProfile(Planet, Params, Layers):
    header = 'Header lines'
    with open(Params.dataFiles.saveFile,'w') as f:
        f.write(header+'\n')
        for i in range(Planet.nStepsTotal):
            #f.write('\t'.join( [ '{:3.5e}'.format(val) for val in line] )+'\n')
            f.write('Test\n')

    print('Profile saved to file: ' + Params.dataFiles.saveFile)
    return


def ReloadProfile(Planet, Params):
    """ Reload previously saved PlanetProfile run from disk """
    Params.dataFiles, Params.figureFiles = SetupFilenames(Planet)
    print('Reloading previously saved run from file: ' + Params.dataFiles.saveFile)
    if Params.VERBOSE: print('WARNING: nSteps settings from PP' + Planet.name + '.py will be ignored.')

    with open(Params.dataFiles.saveFile) as f:
        # Get number of header lines to read in from (and skip for columnar data)
        Params.nHeadLines = int(f.readline().split('=')[-1])
        # Get float values from header
        Planet.Tb_K, Planet.Zb_km, Planet.Pb_MPa, Planet.PbI_MPa, Planet.deltaP \
            = (float(f.readline().split('=')[-1]) for _ in range(5))
        # Get integer values from header (nSteps values)
        Planet.nStepsIceI, Planet.nStepsOcean, Planet.nStepsHydro, Planet.nIceIIILitho, Planet.nIceVLitho \
            = (int(f.readline().split('=')[-1]) for _ in range(5))
    Planet.nStepsClath = Planet.nStepsHydro - Planet.nStepsIceI - Planet.nStepsOcean - Planet.nIceIIILitho - Planet.nIceVLitho

    # Read in columnar data that follows header lines
    Layers = LayersStruct(Planet.nStepsHydro)
    Layers.z_m, Layers.P_MPa, Layers.T_K, Layers.phase, Layers.rho_kgm3, Layers.g_ms2, Layers.Cp \
        = np.loadtxt(Params.dataFiles.saveFile, skiprows=Params.nHeadLines, unpack=True)


    return Planet, Params, Layers


if __name__ == '__main__': main()
