""" Default custom ocean solution settings """
from PlanetProfile.Utilities.defineStructs import CustomSolutionStruct

configCustomSolutionVersion = 1 # Integer number for config file version. Increment when new settings are added to the default config file.

def customSolutionAssign():
    CustomSolutionParams = CustomSolutionStruct()
    # Frezchem database to use - frezchem.dat, frezchemNH3.dat (has ammonia species), frezchemSiCH4.dat (has CH4 and many minerals)
    CustomSolutionParams.FREZCHEM_DATABASE = "frezchemNH3.dat"
    # Unit of species input (can be "g" for grams, or "mol" for mols
    CustomSolutionParams.SPECIES_CONCENTRATION_UNIT = "mol"

    return CustomSolutionParams
