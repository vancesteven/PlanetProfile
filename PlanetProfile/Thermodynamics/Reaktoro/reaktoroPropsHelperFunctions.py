import copy
import reaktoro as rkt
import numpy as np
import logging
import pandas as pd
import os
import pickle as pickle
from seafreeze import seafreeze as sfz
from scipy.interpolate import BSpline, splrep
from PlanetProfile.Utilities.defineStructs import Constants
from PlanetProfile.Utilities.DataManip import ResetNearestExtrap
from scipy import interpolate
log = logging.getLogger('PlanetProfile')


def SpeciesParser(species_string_with_ratios):
    '''
    Converts the provided String of species and their molar ratios into formats necessary for Reaktoro. Namely, creates
    a String of all the species in the list and a dictionary with 'active' species that are added to solution (the observer species are
    automatically generated to be 1e-16 moles in the solution by Reaktoro). It also returns the w_ppt of the solution. If any of the species do not exist in the database Reaktoro is implementing
    (namely frezchem), this method raises an error that the species does not exist.

     Parameters
     ----------
     species_string_with_ratios: String of all the species that should be considered in aqueous phase and their corresponding molar ratios.
        For example, "Cl-: 19.076, Na+: 5.002, Ca2+: 0.0"
     Returns
     -------
     aqueous_species_string: String that has all species names that should be considered in aqueous phase
     speciation_ratio_mol_kg: Dictionary of active species and the values of their molar ratio (mol/kg of water)

    '''
    # Initialize the Phreeqc database with frezchem
    db = rkt.PhreeqcDatabase("frezchem.dat")
    # Create a new string that will hold all the aqueous species in a format compatible with Reaktoro
    aqueous_species_string = ""
    # Create a new dictionary that will hold all the aqueous species with a specified amount to add in a format compatible with Reaktoro
    speciation_ratio_mol_kg = {}
    # Go through each species and corresponding ratio_mol_kg and add to corresponding lists
    for species_with_ratio in species_string_with_ratios.split(", "):
        species, ratio_mol_kg = species_with_ratio.split(": ")
        # Ensure that species is in frezchem database and if not then raise error
        try:
            db.species(species)
        except:
            raise ValueError(f'{species} does not exist in the Phreeqc database. Check that it is entered correctly in the Planet.ocean.species')
        # Add species to string
        aqueous_species_string = aqueous_species_string + species + " "
        # Check if the species is active (amount > 0 mol) and if so, add it to the dictionary
        if (float(ratio_mol_kg) > 0):
            speciation_ratio_mol_kg[species] = float(ratio_mol_kg)
    # Check if water is in the aqueous species string and dictionary and if not, add it, ensuring to update the weight to be 1kg
    if not "H2O" in aqueous_species_string:
        aqueous_species_string = aqueous_species_string + "H2O "
    # Ensure H2O amount is a mol equivalent of 1kg
    speciation_ratio_mol_kg.update({"H2O": float(1/rkt.waterMolarMass)})
    # Return the species string and dictionary (remove the trailing white space from the String as well with rstrip())
    return aqueous_species_string.rstrip(" "), speciation_ratio_mol_kg



def TemperatureCorrectionSplineGenerator():
    current_directory = os.path.dirname(__file__)
    folder = 'Reaktoro_Saved_Files'
    filename = 'rkt_temperature_correction.pkl'
    file_path = os.path.join(current_directory, folder, filename)
    if os.path.exists(file_path):
        with open(file_path, 'rb') as file:
            dictionary = pickle.load(file)
    else:
        dictionary = dictionary_pressure_correction_generator(file_path)
    t, c, k = splrep(dictionary["P_MPa"], dictionary["Difference_in_T_freezings_K"])
    spline = BSpline(t, c, k, extrapolate=True)
    return spline


def PhreeqcGenerator(aqueous_species_list, speciation_ratio_mol_kg, database_name):
    """ Create a Phreeqc Reaktoro System with the solid and liquid phase whose relevant species are determined by the provided aqueous_species_list.
        Works for both core10.dat and frezchem.dat.
    Args:
        aqueous_species_list: aqueous species in reaction. Should be formatted in one long string with a space in between each species
     speciation_ratio_mol_kg: the ratio of species in the aqueous solution in mol/kg of water. Should be a dictionary
     with the species as the key and its ratio as its value.
        database_name: frezchem or core10
    Returns:
        db, system, state, conditions, solver, props, ice_name: Relevant reaktoro objects
    """
    # Initialize the database
    db = rkt.PhreeqcDatabase(database_name)
    # Prescribe the solution
    solution = rkt.AqueousPhase(aqueous_species_list)
    solution.setActivityModel(rkt.chain(rkt.ActivityModelPitzer(), rkt.ActivityModelPhreeqcIonicStrengthPressureCorrection()))
    # Obtain all related solid phases
    solids = rkt.MineralPhases()
    # Initialize the system
    system = rkt.ChemicalSystem(db, solution, solids)
    # Create constraints on equilibrium - pressure and temperature
    specs = rkt.EquilibriumSpecs(system)
    specs.pressure()
    specs.temperature()
    # Create a solver object
    solver = rkt.EquilibriumSolver(specs)
    # Create a chemical state and its associated properties
    state = rkt.ChemicalState(system)
    props = rkt.ChemicalProps(state)
    # Populate the state with the prescribed species at the given ratios
    for ion, ratio in speciation_ratio_mol_kg.items():
        state.add(ion, ratio, "mol")
    # Create a conditions object
    conditions = rkt.EquilibriumConditions(specs)
    # Obtain ice name
    if database_name == "frezchem.dat":
        ice_name = "Ice(s)"
    else:
        ice_name = "Ice"
    # Return the Reaktoro objects that user will need to interact with
    return db, system, state, conditions, solver, props, ice_name, database_name


def PhreeqcGeneratorForChemicalConstraint(aqueous_species_list, speciation_ratio_mol_kg, database):
    """ Create a Phreeqc Reaktoro System with the solid and liquid phase whose relevant species are determined by the provided aqueous_species_list.
        Works for both core10.dat and frezchem.dat.
        THIS IS DIFFERENT IN THAT IT ASSUMES TEMPERATURE IS UNKNOWN AND SPECIFIES CHEMICAL CONSTRAINT AT EQUILIBIRUM.
    Args:
        aqueous_species_list: aqueous species in reaction. Should be formatted in one long string with a space in between each species
     speciation_ratio_mol_kg: the ratio of species in the aqueous solution in mol/kg of water. Should be a dictionary
     with the species as the key and its ratio as its value.
        database: frezchem or core10
    Returns:
        db, system, state, conditions, solver, props: Relevant reaktoro objects
    """
    # Initialize the database
    db = rkt.PhreeqcDatabase(database)
    # Prescribe the solution
    solution = rkt.AqueousPhase(aqueous_species_list)
    solution.setActivityModel(rkt.chain(rkt.ActivityModelPitzer(), rkt.ActivityModelPhreeqcIonicStrengthPressureCorrection()))
    # Obtain all related solid phases
    solids = rkt.MineralPhases()
    # Initialize the system
    system = rkt.ChemicalSystem(db, solution, solids)
    # Create constraints on equilibrium - pressure and temperature
    specs = rkt.EquilibriumSpecs(system)
    specs.pressure()
    specs.unknownTemperature()
    # Create equilibrium constraint on the phase amount of ices to significant threshold
    # This constraint will allow Reaktoro to query for the pressure at which ice begins to form at the prescribed pressure
    idx_ice_phase = specs.addInput("IP")
    ices_phase_constraint = rkt.EquationConstraint()
    ices_phase_constraint.id = "icePhaseAmountConstraint"
    ices_phase_constraint.fn = lambda props, w: ices_phases_amount_mol(props) - w[idx_ice_phase]
    specs.addConstraint(ices_phase_constraint)
    # Create a solver object
    solver = rkt.EquilibriumSolver(specs)
    # Create a chemical state and its associated properties
    state = rkt.ChemicalState(system)
    props = rkt.ChemicalProps(state)
    # Populate the state with the prescribed species at the given ratios
    for ion, ratio in speciation_ratio_mol_kg.items():
        state.add(ion, ratio, "mol")
    # Create a conditions object
    conditions = rkt.EquilibriumConditions(specs)
    # Return the Reaktoro objects that user will need to interact with
    return db, system, state, conditions, solver, props


def SupcrtGenerator(aqueous_species_list, speciation_ratio_mol_kg, database):
    """ Create a Supcrt Reaktoro System with the solid and liquid phase whose relevant species are determined by the provided aqueous_species_list.
    Args:
    aqueous_species_list: aqueous species in reaction. Should be formatted in one long string with a space in between each species
     speciation_ratio_mol_kg: the ratio of species in the aqueous solution in mol/kg of water. Should be a dictionary
     with the species as the key and its ratio as its value.
    database: Supcrt database to use
    Returns:
        db, system, state, conditions, solver, props, ice_name: Relevant reaktoro objects
    """
    # Initialize the database
    db = rkt.SupcrtDatabase(database)
    # Prescribe the solution
    solution = rkt.AqueousPhase(aqueous_species_list)
    # Obtain all related solid phases
    solids = rkt.MineralPhases()
    # Initialize the system
    system = rkt.ChemicalSystem(db, solution, solids)
    # Create constraints on equilibrium - pressure and temperature
    specs = rkt.EquilibriumSpecs(system)
    specs.pressure()
    specs.temperature()
    # Create a solver object
    solver = rkt.EquilibriumSolver(specs)
    # Create a chemical state and its associated properties
    state = rkt.ChemicalState(system)
    props = rkt.ChemicalProps(state)
    # Populate the state with the prescribed species at the given ratios
    for ion, ratio in speciation_ratio_mol_kg.items():
        state.add(ion, ratio, "mol")
    # Create a conditions object
    conditions = rkt.EquilibriumConditions(specs)
    # Return the Reaktoro objects that user will need to interact with
    return db, system, state, conditions, solver, props


def reset_state(system, speciation_ratio_mol_kg):
    """ Returns a new Reaktoro state for given system populated with provided species.

    Args:
        system: Reaktoro system
        speciation_ratio_mol_kg: the ratio of species in the aqueous solution in mol/kg of water. Should be a dictionary
    Returns:
        state: Reaktoro state populated with provided species
    """
    state = rkt.ChemicalState(system)
    for ion, ratio in speciation_ratio_mol_kg.items():
        state.add(ion, ratio, "mol")
    return state


def ices_phases_amount_mol(props: rkt.ChemicalProps):
    '''
    Helper equilibirum constraint function to constrain the total amount of moles of all ices in the current state using its associated properties. Function is used in
    the ice constraint for rkt_p_freeze().
    '''
    # Get name of H2O ice (either Ice(s) for frezchem or Ice for core10.dat)
    ice_name = ""
    try:
        # Check if ice is labeled as "Ice(s)" in database (core10.dat labeling)
        props.system().database().species("Ice(s)")
        ice_name = "Ice(s)"
    except Exception as e:
        # If exception is thrown, then we are using frezchem.dat database and ice_name should be "Ice"
        ice_name = "Ice"
    # Return the amount of solid H2O in the state, given by moles
    ice_chem_potential = props.speciesChemicalPotential(ice_name)
    water_chem_potential = props.speciesChemicalPotential("H2O")
    return ice_chem_potential - water_chem_potential


def species_convertor_compatible_with_supcrt(aqueous_species_string, speciation_ratio_mol_kg):
    """
    Converts aqueous species string and speciation ratio dictionary into formats compatible with supcrt. Namely, in phreeqc the liquid phase of H2O
    is labeled "H2O", whereas in supcrt it requires "H2O(aq)". Thus, converts "H2O" in the string and speciation ratio dictionary to "H2O(aq)".
    Importantly, since speciation_ratio_mol_kg is a dictionary, we must make a deep copy before editing so as to not disturb the original dictionary, which
    will still be used by the phreeqc database in the phase change function.
    Args:
        aqueous_species_string: String that has all species names that should be considered in aqueous phase
        speciation_ratio_mol_kg: Dictionary of active species and the values of their molar ratio (mol/kg of water)
    Returns:
        aqueous_species_string: Adapted string that has all species names in format compatible with supcrt
        supcrt_speciation_ratio_mol_kg: Deep copy of dictionary that has species names in format compatible with supcrt

    """
    # Since python passes dictionary by reference, need to make deep copy to preserve original dictionary
    # Check if "H2O" is in the string (and thus dictionary), indicating it is not compatible with supcrt
    deep_copy_ratio_mol_kg = copy.deepcopy(speciation_ratio_mol_kg)
    if "H2O" in aqueous_species_string:
        # Change "H2O" to "H2O(aq) in the string
        aqueous_species_string = aqueous_species_string.replace("H2O", "H2O(aq)")
        # Now change the "H2O" key in the dictionary to "H2O(aq)"
        deep_copy_ratio_mol_kg["H2O(aq)"] = deep_copy_ratio_mol_kg.pop("H2O")
        # Return the string and adapted dictionary
    if "CO2" in aqueous_species_string:
        # Change "H2O" to "H2O(aq) in the string
        aqueous_species_string = aqueous_species_string.replace("CO2", "CO2(aq)")
        # Now change the "H2O" key in the dictionary to "H2O(aq)"
        deep_copy_ratio_mol_kg["CO2(aq)"] = deep_copy_ratio_mol_kg.pop("CO2")
    return aqueous_species_string, deep_copy_ratio_mol_kg


def McClevskyIonParser(aqueous_species_list, speciation_ratio_mol_kg):
    ions = {}
    for species, mol_kg in speciation_ratio_mol_kg.items():
        # Check that the species is an ion
        if "+" in species or "-" in species:
            # Rewrite the species into a format that McCLevsky can handle
            # Namely, change the + to a _p or the - to a _m
            if "+" in species:
                # If there is no number after the +, then we must append _p1, not just _p
                if species.endswith("+"):
                    species = species + "1"
                # CHange the - to a _p
                species = species.replace('+', '_p')
            elif "-" in species:
                # If there is no number after the -, then we must append a 1 to signify it is a -1 charge
                if species.endswith("-"):
                    species = species + "1"
                # Change the - to a _m
                species = species.replace('-', '_m')
            # Format the mol amount into a dictionary
            mol_dictionary = {'mols': mol_kg}
            # Append the new species and mol dictionary to the ion dictionary
            ions[species] = mol_dictionary
        else:
            # If the species is not an ion, then it is not relevant to electrical conductivity so do not add to ion dictionary
            pass
    # Return ions list
    return ions



def database_to_use(T_K):
    """ Decides between frezchem.dat and core10.dat database to use according to the temperature array.
        As per Comparison of thermodynamic data files for PHREEQC (Lu et. al 2022):
            frezchem.dat: −73–25°C at 1 bar
            core10.dat: 0.01–100°C at 1bar, 100–300°C along P_Sat (water vapor saturation pressure)
        Currently, frezchem will be utilized for 200K<T<298.15T, and dat10 for T>=298.15.
    Args:
        T_K (shape N): The temperatures in increasing order to determine between frezchem.dat and core10.dat
    Returns:
        index: The starting index of the T_K at which core10.dat should be used. All indices before will be assumed to use frezchem.dat
            If the index is size N, this indicates that all temperatures are within frezchem limits.
    """
    # The temperature to switch between frezchem and core10
    switch_T_K = 298.15
    if T_K.ndim == 1:
        index = np.searchsorted(T_K, switch_T_K, side = 'left')
    elif T_K.ndim == 2:
        row_of_temperatures = T_K[0]
        index = np.searchsorted(row_of_temperatures, switch_T_K, side='left')
    return index


def pressure_constraint(P_MPa, aqueous_species_list, speciation_ratio_mol_kg, database, dP = -1):
    """ Find the pressure constraint at which Reaktoro can find equilibrium for the given speciation and database. Starts at P_MPa and checks if rkt can find equilibrium
        with a temperature of 273K. If it cannot, then adjusts P_MPa by dP and tries again, continuing this process until a compatible P_MPa is found.
        If we are looking for upper pressure cosntraint (dP is negative), then we stop also if it reaches <= 1.1 MPa, which Rkt should be compatible with.
    Args:
        P_MPa: Initial pressure constraint in MPa
        aqueous_species_list: String that has all species names that should be considered in aqueous phase
        speciation_ratio_mol_kg: Dictionary of active species and the values of their molar ratio (mol/kg of water)
        database: Database to find pressure constraint for
        dP: The amount to change P_MPa by if equilibrium is not achieved. Defaults to -1 (for upper constraint)

    Returns:
        P_MPa: New pressure constraint in MPa for given database.

    """
    # Disable chemical convergence warnings that Reaktoro raises. We handle these internally instead and throw more specific warnings when they appear.
    rkt.Warnings.disable(906)
    # Initilialize the database, either being supcrt of phreeqc
    if "supcrt" in database:
        # Since supcrt labels "H2O(aq)", we need to adjust our species list and dictionary accordingly
        aqueous_species_list, speciation_ratio_mol_kg = species_convertor_compatible_with_supcrt(
            aqueous_species_list, speciation_ratio_mol_kg)
        db, system, state, conditions, solver, props = SupcrtGenerator(aqueous_species_list, speciation_ratio_mol_kg, database)
    else:
        db, system, state, conditions, solver, props, ice_name, database_name = PhreeqcGenerator(aqueous_species_list, speciation_ratio_mol_kg, database)
    # Establish pressure constraint of 273 K
    conditions.temperature(273, "K")
    # Tracker variable for equilibrium being found
    success = False
    # While we have not found chemical equilibrium, keep iterating through
    while not success:
        # Establish pressure constraint at P_MPa
        conditions.pressure(P_MPa, "MPa")
        # Solve the equilibrium problem
        result = solver.solve(state, conditions)
        # Check if the equilibrium problem succeeded
        if result.succeeded():
            # If so, change success to true
            success = True
        # To prevent an infinite loop, added this statement which will stop when P_MPa is <= 1.1 and dP < 0, meanign we are looking for upper pressure constraint
        elif P_MPa <= 1.1 and dP < 0:
            log.debug("Infinite loop prevention stopped for pressure_constraint. Should likely not be occurring")
            success = True
        # Otherwise, equilibrium is not achieved so we should reset the state and adjust P_MPa by dP, and reattempt the equilibrium problem
        else:
            state = reset_state(system, speciation_ratio_mol_kg)
            P_MPa += dP
    # Return the adjusted P_MPa
    return P_MPa


def temperature_constraint(T_K, aqueous_species_list, speciation_ratio_mol_kg, database, dT = 1):
    """ Find the tempearture constraint at which Reaktoro can find equilibrium for the given speciation and database. Starts at T_K and checks if rkt can find equilibrium
        with a pressure of 0.1 MPa. If it cannot, then adjusts T_K by dT and tries again, continuing this process until a compatible T_K is found.
        For finding the lower temperature constraint if it reaches >= 273K we stop, which Rkt should be compatible with.
    Args:
        T_K: Initial temperature constraint in K
        aqueous_species_list: String that has all species names that should be considered in aqueous phase
        speciation_ratio_mol_kg: Dictionary of active species and the values of their molar ratio (mol/kg of water)
        database: Database to find temperature constraint for
        dT: The amount to change T_K by if equilibrium is not achieved. Defaults to 1 (for lower constraint)

    Returns:
        T_K: New temperature constraint in T_K for given database.

    """
    # Disable chemical convergence warnings that Reaktoro raises. We handle these internally instead and throw more specific warnings when they appear.
    rkt.Warnings.disable(906)
    # Initilialize the database, either being supcrt of phreeqc
    if "supcrt" in database:
        # Since supcrt labels "H2O(aq)", we need to adjust our species list and dictionary accordingly
        aqueous_species_list, speciation_ratio_mol_kg = species_convertor_compatible_with_supcrt(
            aqueous_species_list, speciation_ratio_mol_kg)
        db, system, state, conditions, solver, props = SupcrtGenerator(aqueous_species_list, speciation_ratio_mol_kg, database)
    else:
        db, system, state, conditions, solver, props, ice_name, database_name = PhreeqcGenerator(aqueous_species_list, speciation_ratio_mol_kg, database)
    # Establish pressure constraint of 1bar
    conditions.pressure(0.1, "MPa")
    # Tracker variable for equilibrium being found
    success = False
    # While we have not found chemical equilibrium, keep iterating through
    while not success:
        # Establish temperature constraint at T_K
        conditions.temperature(T_K, "K")
        # Solve the equilibrium problem
        result = solver.solve(state, conditions)
        # Check if the equilibrium problem succeeded
        if result.succeeded():
            # If so, change success to true
            success = True
        # If it failed but T_K >= 273 and dT > 0, which indicates we are finding lower constraint, then set success to true since RKt should be compatible
        elif T_K >= 273 and dT > 0:
            log.debug("Infinite loop prevention stopped for temperature_constraint. Should likely not be occurring")
            success = True
        # Otherwise, equilibrium is not achieved so we should reset the state and adjust T_K by dT, and reattempt the equilibrium problem
        else:
            T_K += dT
            state = reset_state(system, speciation_ratio_mol_kg)
    # Return the adjusted T_K
    return T_K

def interpolation_2d_seismic(zero_indices, sound_speeds, densities):
    initial_array_lists = [sound_speeds, densities]
    updated_array_lists = []
    for array in initial_array_lists:
        # Create mask of non-zero values
        nonzero_mask = (array != 0)
        # Interpolate using scipy interp2d function
        x, y = np.meshgrid(np.arange(array.shape[1]), np.arange(array.shape[0]))
        x_interp = x[nonzero_mask]
        y_interp = y[nonzero_mask]
        z_interp = array[nonzero_mask]
        updated_array_lists.append(interpolate.griddata((x_interp, y_interp), z_interp, (x, y)))
    sound_speeds = updated_array_lists[0]
    densities = updated_array_lists[1]
    log.warning("Performed 2d linear interpolation on missing thermodynamic properties")
    return sound_speeds, densities
def interpolation_2d(zero_indices, rho_kgm3, Cp_JKgK, alpha_pK):
    """ Utilized as a helper function for thermodynamic properties calculation. Performs a 2d interpolation on any values that are non-zero in the
        provided arrays, allowing the zero values to be interpolated.
    Args:
        zero_indices: Array of indices where there are zeros
        rho_kgm3 (float, shape NxM): Mass density of liquid in kg/m^3
        Cp_JkgK (float, shape NxM): Isobaric heat capacity of liquid in J/(kg K)
        alpha_pK (float, shape NxM): Thermal expansivity of liquid in 1/K
    Returns:
        rho_kgm3 (float, shape NxM): Mass density of liquid in kg/m^3 whose zero values have been replaced by interpolated values
        Cp_JkgK (float, shape NxM): Isobaric heat capacity of liquid in J/(kg K) whose zero values have been replaced by interpolated values
        alpha_pK (float, shape NxM): Thermal expansivity of liquid in 1/K whose zero values have been replaced by interpolated values

    """
    initial_array_lists = [rho_kgm3, Cp_JKgK, alpha_pK]
    updated_array_lists = []
    for array in initial_array_lists:
        # Create mask of non-zero values
        nonzero_mask = (array != 0)
        # Interpolate using scipy interp2d function
        x, y = np.meshgrid(np.arange(array.shape[1]), np.arange(array.shape[0]))
        x_interp = x[nonzero_mask]
        y_interp = y[nonzero_mask]
        z_interp = array[nonzero_mask]
        updated_array_lists.append(interpolate.griddata((x_interp, y_interp), z_interp, (x, y)))
    rho_kgm3 = updated_array_lists[0]
    Cp_JKgK = updated_array_lists[1]
    alpha_pK = updated_array_lists[2]
    log.warning("Performed 2d linear interpolation on missing thermodynamic properties")
    return rho_kgm3, Cp_JKgK, alpha_pK


def interpolation_1d(zero_indices, array1, array2):
    array_lists = [array1, array2]
    for array in array_lists:
        # Create mask of non-zero values
        nonzero_mask = (array != 0)
        # Interpolate using scipy interp2d function
        x_interp = np.arange(len(array))
        x_interp_nonzero = x_interp[nonzero_mask]
        y_interp_nonzero = array[nonzero_mask]
        f = interpolate.interp1d(x_interp_nonzero, y_interp_nonzero, kind = 'linear', fill_value = 'extrapolate')
        for idx in zero_indices:
            if array[idx] == 0:
                array[idx] = f(idx)
    log.warning("Performed 1d linear interpolation on missing seismic properties")
    return array1, array2


def panda_df_generator(system):
    columns = ["P (MPa)", "T (K)", "pH", "Amount", "Charge", "SolidAmount", "LiquidAmount"]
    species = []
    for speciesItem in system.species():
        species.append(speciesItem.name())
    # Dictionary that maps certain column names that require using species list to Reaktoro code
    translation_dictionary_for_species_list = {"Amount": ["amount" + name for name in species]}
    df_column_names = []
    for column in columns:
        if column in translation_dictionary_for_species_list:
            df_column_names += translation_dictionary_for_species_list[column]
        else:
            df_column_names += [column]
    df = pd.DataFrame(columns=df_column_names)
    return df


def dictionary_pressure_correction_generator(file_path_to_save):
    log.warning(f"Reaktoro's current implementation of Frezchem does not properly correct for the effect of pressure. In order to remedy this correction,\n" +
    f"we create an internal correction spline that is generated by comparing the difference in freezing temperatures between a pure solution of water of Reaktoro and Seafreeze over a range of pressures\n."
                + f"We then apply this spline to adjust the calculated equilibrium temperatures for Reaktoro.")

    eos_P_MPa = np.linspace(0.1, 100, 150)
    eos_T_K = np.linspace(260, 274, 300)
    eos_PT = np.array([eos_P_MPa, eos_T_K], dtype=object)
    sfz_phases = sfz.whichphase(eos_PT)
    diff = np.diff(sfz_phases)
    sfz_T_freezing = []
    for row in diff:
        i = np.where(row == -1.0)[0]
        sfz_T_freezing.append(((eos_T_K[i] + eos_T_K[i+1])/2)[0])


    pure_water_speciation = "H+: 1e-7, OH-: 1e-7"
    aqueous_species_list, speciation_ratio_mol_kg = SpeciesParser(pure_water_speciation)
    frezchem = PhreeqcGeneratorForChemicalConstraint(aqueous_species_list, speciation_ratio_mol_kg, "frezchem.dat")
    rkt_T_freezing = rkt_t_freeze_without_pressure_correction(aqueous_species_list, speciation_ratio_mol_kg, eos_P_MPa, frezchem, 100, 0, {})

    sfz_T_freezing = np.array(sfz_T_freezing)
    rkt_T_freezing = np.array(rkt_T_freezing, dtype = np.float64)

    difference_in_T_freezings =  sfz_T_freezing - rkt_T_freezing

    folder = 'Reaktoro_Saved_Files'

    # Ensure the directory exists
    current_directory = os.path.dirname(__file__)
    folder_path = os.path.join(current_directory, folder)
    os.makedirs(folder_path, exist_ok=True)

    dictionary_to_save = {"P_MPa": eos_P_MPa, "Difference_in_T_freezings_K": difference_in_T_freezings}
    with open(file_path_to_save, 'wb') as f:
        pickle.dump(dictionary_to_save, f)
    return dictionary_to_save

def rkt_t_freeze_without_pressure_correction(aqueous_species_list, speciation_ratio_mol_kg, P_MPa, frezchem, PMax_MPa, freezing_temperature_spline, calculated_freezing_temperatures, TMin_K = 250, TMax_K = 300, significant_threshold = 0.1):
    """
     Calculates the temperature at which the prescribed aqueous solution freezes. Utilizes the reaktoro framework to
     constrain the equilibrium position at the prescribed pressure, the lower and upper limits of temperature (in K),
     and the total amount of ice at a significant threshold of 1e-14, therefore calculating and returning the
     temperature (within the range) at which ice begins to form.

     Parameters
     ----------
     aqueous_species_list: aqueous species in reaction. Should be formatted in one long string with a space in between each species
     speciation_ratio_mol_kg: the ratio of species in the aqueous solution in mol/kg of water. Should be a dictionary
     with the species as the key and its ratio as its value.
     P_MPa: the desired equilibrium freezing pressure(s).
     TMin_K: the lower limit of temperature that Reaktoro should query over
     TMax_K: the upper limit of temperature that Reaktoro should query over
     database: the Phreeqc database to be using in this calculation (either frezchem.dat or core10.dat)
     significant_threshold: the amount of moles of ice present for H2O to be considered in solid phase. Default is 1e-14 moles.

     Returns
     -------
     t_freezing_K: the temperature at which the solution begins to freeze.
     """
    # Disable chemical convergence warnings that Reaktoro raises. We handle these internally instead and throw more specific warnings when they appear.
    # rkt.Warnings.disable(906)
    # Create list that holds boolean values of whether ice is present
    freezing_temperatures = []
    db, system, state, conditions, solver, props = frezchem
    state = rkt.ChemicalState(state)
    # Create an iterator to go through P_MPa
    it = np.nditer([P_MPa])
    conditions.set("IP", significant_threshold)
    conditions.setLowerBoundTemperature(TMin_K, "K")
    conditions.setUpperBoundTemperature(TMax_K, "K")
    for P in it:
        P = float(P)
        # Adjust Pressure
        if P in calculated_freezing_temperatures.keys():
            equilibrium_temperature = calculated_freezing_temperatures[P]
            freezing_temperatures.append(equilibrium_temperature)
            continue
        # Check that pressure is below PMax_MPa, the maximum pressure constraint calculated previously
        if P <= PMax_MPa:
            # Specify equilibrium constraints
            conditions.pressure(P, "MPa")
            # Solve the equilibrium problem
            result = solver.solve(state, conditions)
            # Update the properties
            props.update(state)
            # Obtain the equilibrium temperature
            equilibrium_temperature = props.temperature()
            # Check if the result succeeded
            if result.succeeded():
                state = reset_state(system, speciation_ratio_mol_kg)
            # If the result failed, we will use spline
            else:
                log.warning(f"While attempting to find bottom freezing temperature for pressure of {P_MPa} MPa, \n"
                            +f"Reaktoro was unable to find a temperature within range of {TMin_K} K and {TMax_K}.\n"
                            +f"Instead, we will use a spline of freezing temperatures that we generated for this EOS to find the associated freezing temperature.")
                equilibrium_temperature = freezing_temperature_spline(P)
                state = reset_state(system, speciation_ratio_mol_kg)
        # If Pressure is >= PMax, Reaktoro is at its limits of computation and will likely fail. Thus, we utilize our spline approximation of freezing temperatures
        # over a range of pressures that we are confident that Frezchem works at and find the associated freezing temperature at the given pressure.
        # If the temperature we are querying over is less than that freezing temperature, then we assume its ice, otherwise we assume its liquid
        else:
            equilibrium_temperature = freezing_temperature_spline(P)
        freezing_temperatures.append(equilibrium_temperature)
        calculated_freezing_temperatures[P] = equilibrium_temperature
    # Return the equilibrium temperature
    return np.array(freezing_temperatures)