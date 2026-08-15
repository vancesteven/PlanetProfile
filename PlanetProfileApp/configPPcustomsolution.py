""" Default custom ocean solution settings """
from PlanetProfile.Utilities.defineStructs import CustomSolutionParamsStruct

configCustomSolutionVersion = 6 # Integer number for config file version. Increment when new settings are added to the default config file.

def customSolutionAssign():
    CustomSolutionParams = CustomSolutionParamsStruct()

    # Frezchem database to use - frezchem.dat, frezchemNH3.dat (has ammonia and CO2 species), frezchemSiCH4.dat (has CH4 and many minerals)
    CustomSolutionParams.FREZCHEM_DATABASE = "frezchemNH3.dat"
    # Supcrt database to use - supcrt16, supcrt16-organics (has organic species) are primary options
    CustomSolutionParams.SUPCRT_DATABASE = "supcrt16"

    # Unit of species input (can be "g" for grams, or "mol" for mols
    CustomSolutionParams.SPECIES_CONCENTRATION_UNIT = "mol"

    # Pressure and temperature increment to use when generating EOS for custom solution. Smaller increments increases the accuracy of grid but increases initial runtime when generating grid
    CustomSolutionParams.EOS_deltaP = 2.0
    CustomSolutionParams.EOS_deltaT = 2.0

    # Consider solid phases when calculating thermodynamic properties
    CustomSolutionParams.SOLID_PHASES = True
    # Only valid if SOLID_PHASES is True. Specify the Solid Phases to consider - specifying only primary phases can speed up runtime. 'All' considers all possible phases available in database
    CustomSolutionParams.SOLID_PHASES_TO_CONSIDER = ['Carbonates', 'Sulfates'] # Can specify minerals specifically, or add keywords that include 'Carbonates', 'Sulfates' that are defined in Constants.py

    # Have PlanetProfile remove species that Frezchem does not have in database so we can consider more diverse speciated compositions in the ocean thermodynamcis
    # This removes self-consistency between the phase equilibria (up to 200MPa) and ocean thermodynamics since frezchem is calculating liquid-IceI equilibria of a less speciated chemistry
    # Example usage is when user wants to consider an ocean with Fe aqueous species, which Frezchem does not have. So we remove any Fe species from the Frezchem system, but they are still considered in Supcrt system.
    CustomSolutionParams.REMOVE_SPECIES_NA_IN_FREZCHEM = False

    # Maximum number of iterations to allow when calculating equilibrium before throwing convergence error
    # For complicated systems, increase this number by an order of magnitude (i.e. 2000) - dramatically incresaes runtime though
    CustomSolutionParams.maxIterations = 200

    return CustomSolutionParams
