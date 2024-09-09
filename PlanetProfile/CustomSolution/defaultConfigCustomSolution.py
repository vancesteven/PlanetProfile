""" Default custom ocean solution settings """
from PlanetProfile.Utilities.defineStructs import CustomSolutionParamsStruct

configCustomSolutionVersion = 2 # Integer number for config file version. Increment when new settings are added to the default config file.

def customSolutionAssign():
    CustomSolutionParams = CustomSolutionParamsStruct()

    # Frezchem database to use - frezchem.dat, frezchemNH3.dat (has ammonia species), frezchemSiCH4.dat (has CH4 and many minerals)
    CustomSolutionParams.FREZCHEM_DATABASE = "frezchemNH3.dat"

    # Unit of species input (can be "g" for grams, or "mol" for mols
    CustomSolutionParams.SPECIES_CONCENTRATION_UNIT = "mol"

    # Pressure and temperature increment to use when generating EOS for custom solution. Smaller increments increases the accuracy of grid but increases initial runtime when generating grid
    CustomSolutionParams.EOS_deltaP = 1.0
    CustomSolutionParams.EOS_deltaT = 0.5

    return CustomSolutionParams
