""" Default custom ocean solution settings """
from PlanetProfile.Utilities.defineStructs import CustomSolutionParamsStruct

configCustomSolutionVersion = 4 # Integer number for config file version. Increment when new settings are added to the default config file.

def customSolutionAssign():
    CustomSolutionParams = CustomSolutionParamsStruct()

    # Frezchem database to use - frezchem.dat, frezchemNH3.dat (has ammonia and CO2 species), frezchemSiCH4.dat (has CH4 and many minerals)
    CustomSolutionParams.FREZCHEM_DATABASE = "frezchemNH3.dat"
    CustomSolutionParams.SUPCRT_DATABASE = "supcrt16"

    # Unit of species input (can be "g" for grams, or "mol" for mols
    CustomSolutionParams.SPECIES_CONCENTRATION_UNIT = "mol"

    # Pressure and temperature increment to use when generating EOS for custom solution. Smaller increments increases the accuracy of grid but increases initial runtime when generating grid
    CustomSolutionParams.EOS_deltaP = 2.0
    CustomSolutionParams.EOS_deltaT = 2.0

    # Consider solid phases when calculating thermodynamic properties
    CustomSolutionParams.SOLID_PHASES = True
    # Only valid if SOLID_PHASES is True. Specify the Solid Phases to consider - specifying only primary phases can speed up runtime. None considers all possible phases available in database
    CustomSolutionParams.SOLID_PHASES_TO_CONSIDER = ['Carbonates', 'Sulfates'] # Can specify minerals specifically, or add keywords that include 'Carbonates', 'Sulfates'

    return CustomSolutionParams
