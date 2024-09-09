import copy
import reaktoro as rkt
import numpy as np
import logging
import pandas as pd
import os
import pickle as pickle
from seafreeze import seafreeze as sfz
from scipy.interpolate import BSpline, splrep
from scipy import interpolate
log = logging.getLogger('PlanetProfile')


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


def PhreeqcGenerator(aqueous_species_list, speciation_ratio_per_kg, species_unit, database_file):
    """ Create a Phreeqc Reaktoro System with the solid and liquid phase whose relevant species are determined by the provided aqueous_species_list.
        Works for both core10.dat and frezchem.dat.
    Args:
        aqueous_species_list: aqueous species in reaction. Should be formatted in one long string with a space in between each species
        speciation_ratio_per_kg: the ratio of species in the aqueous solution in per kg of water. Should be a dictionary
            with the species as the key and its ratio as its value.
        species_unit: "mol" or "g" that species ratio is in
        database_file: file path of Phreeqc database
    Returns:
        db, system, state, conditions, solver, props, ice_name: Relevant reaktoro objects
    """
    # Initialize the database
    db = rkt.PhreeqcDatabase.fromFile(database_file)
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
    for ion, ratio in speciation_ratio_per_kg.items():
        state.add(ion, ratio, species_unit)
    # Create a conditions object
    conditions = rkt.EquilibriumConditions(specs)
    # Obtain ice name
    if "frezchem" in database_file:
        ice_name = "Ice(s)"
    else:
        ice_name = "Ice"
    # Return the Reaktoro objects that user will need to interact with
    return db, system, state, conditions, solver, props, ice_name, database_file


def PhreeqcGeneratorForChemicalConstraint(aqueous_species_list, speciation_ratio_per_kg, species_unit, database_file):
    """ Create a Phreeqc Reaktoro System with the solid and liquid phase whose relevant species are determined by the provided aqueous_species_list.
        Works for both core10.dat and frezchem.dat.
        THIS IS DIFFERENT IN THAT IT ASSUMES TEMPERATURE IS UNKNOWN AND SPECIFIES CHEMICAL CONSTRAINT AT EQUILIBIRUM.
    Args:
        aqueous_species_list: aqueous species in reaction. Should be formatted in one long string with a space in between each species
        speciation_ratio_per_kg: the ratio of species in the aqueous solution per kg of water. Should be a dictionary
            with the species as the key and its ratio as its value.
        species_unit: "mol" or "g" that species ratio is in
        database_file: file path of Phreeqc database
    Returns:
        db, system, state, conditions, solver, props: Relevant reaktoro objects
    """
    # Initialize the database
    db = rkt.PhreeqcDatabase.fromFile(database_file)
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
    for ion, ratio in speciation_ratio_per_kg.items():
        state.add(ion, ratio, species_unit)
    # Create a conditions object
    conditions = rkt.EquilibriumConditions(specs)
    # Return the Reaktoro objects that user will need to interact with
    return db, system, state, conditions, solver, props


def SupcrtGenerator(aqueous_species_list, speciation_ratio_per_kg, species_unit, database):
    """ Create a Supcrt Reaktoro System with the solid and liquid phase whose relevant species are determined by the provided aqueous_species_list.
    Args:
    aqueous_species_list: aqueous species in reaction. Should be formatted in one long string with a space in between each species
    speciation_ratio_per_kg: the ratio of species in the aqueous solution per kg of water. Should be a dictionary
        with the species as the key and its ratio as its value.
    species_unit: "mol" or "g" that species ratio is in
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
    for ion, ratio in speciation_ratio_per_kg.items():
        state.add(ion, ratio, species_unit)
    # Create a conditions object
    conditions = rkt.EquilibriumConditions(specs)
    # Return the Reaktoro objects that user will need to interact with
    return db, system, state, conditions, solver, props


def reset_state(system, speciation_ratio_per_kg, species_unit):
    """ Returns a new Reaktoro state for given system populated with provided species.

    Args:
        system: Reaktoro system
        speciation_ratio_per_kg: the ratio of species in the aqueous solution in mol/kg of water. Should be a dictionary
        species_unit: "mol" or "g" that species ratio is in
    Returns:
        state: Reaktoro state populated with provided species
    """
    state = rkt.ChemicalState(system)
    for ion, ratio in speciation_ratio_per_kg.items():
        state.add(ion, ratio, species_unit)
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


def species_convertor_compatible_with_supcrt(aqueous_species_string, speciation_ratio_per_kg, Phreeqc_to_Supcrt_names):
    """
    Converts aqueous species string and speciation ratio dictionary into formats compatible with supcrt. Namely, in phreeqc the liquid phase of H2O
    is labeled "H2O", whereas in supcrt it requires "H2O(aq)". Thus, converts "H2O" in the string and speciation ratio dictionary to "H2O(aq)".
    Importantly, since speciation_ratio_mol_kg is a dictionary, we must make a deep copy before editing so as to not disturb the original dictionary, which
    will still be used by the phreeqc database in the phase change function.
    Args:
        aqueous_species_string: String that has all species names that should be considered in aqueous phase
        speciation_ratio_per_kg: Dictionary of active species and the values of their molar ratio (mol/kg of water)
        Phreeqc_to_Supcrt_names: Dictionary of Phreeqc names that must be converted to Supcrt for compatibility
    Returns:
        aqueous_species_string: Adapted string that has all species names in format compatible with supcrt
        supcrt_speciation_ratio_per_kg: Deep copy of dictionary that has species names in format compatible with supcrt
    """
    # Since python passes dictionary by reference, need to make deep copy to preserve original dictionary
    supcrt_speciation_ratio_per_kg = copy.deepcopy(speciation_ratio_per_kg)
    # Dictionary of known values that are different between Phreeqc and supcrt
    for phreeqc_name, supcrt_name in Phreeqc_to_Supcrt_names.items():
        # Check if name is in the string (and thus dictionary), indicating it is not compatible with supcrt
        if phreeqc_name in aqueous_species_string:
            # Change label
            aqueous_species_string = aqueous_species_string.replace(phreeqc_name, supcrt_name)
            # Now change the phreeqc key in the dictionary to supcrt key
            supcrt_speciation_ratio_per_kg[supcrt_name] = supcrt_speciation_ratio_per_kg.pop(phreeqc_name)
    # Return the string and adapted dictionary
    return aqueous_species_string, supcrt_speciation_ratio_per_kg


def McClevskyIonParser(speciation_ratio_per_kg, species_unit):
    """
    Parse through provided species list and convert to format compatible with McClevsky
    Args:
        speciation_ratio_per_kg: the ratio of species in the aqueous solution in mol/kg of water. Should be a dictionary
        species_unit: "mol" or "g" that species ratio is in
    Returns:
        Ion list compatible with McClevsky
    """
    ions = {}
    for species, per_kg in speciation_ratio_per_kg.items():
        # Check that the species is an ion
        if "+" in species or "-" in species:
            # Rewrite the species into a format that McCLevsky can handle
            # First, convert speciation ratio into mols to be compatible with McClevskyIonParser
            mol_amount = per_kg
            if species_unit == "g":
                formula = rkt.ChemicalFormula(species)
                molar_mass_kg_mol = formula.molarMass()
                mol_amount = (1/molar_mass_kg_mol)/1000*per_kg
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
            mol_dictionary = {'mols': mol_amount}
            # Append the new species and mol dictionary to the ion dictionary
            ions[species] = mol_dictionary
        else:
            # If the species is not an ion, then it is not relevant to electrical conductivity so do not add to ion dictionary
            pass
    # Return ions list
    return ions

def interpolation_2d(P_MPa, arrays):
    """ Utilized as a helper function for thermodynamic properties calculation. Performs a 2d interpolation on any values that are NaN in the
        provided arrays, allowing the zero values to be interpolated."""
    interpolated_arrays = []
    P_MPa = P_MPa[:, 0]
    for array in arrays:
        interpolated_array = np.copy(array)
        for col in range(array.shape[1]):
            column_data = array[:, col]
            nan_mask = np.isnan(column_data)
            if np.any(nan_mask):
                x_known = P_MPa[~nan_mask]
                y_known = column_data[~nan_mask]
                spline = interpolate.make_interp_spline(x_known, y_known, k = 2)
                interpolated_values = spline(P_MPa)
                interpolated_array[:, col] = interpolated_values
        interpolated_arrays.append(interpolated_array)
    return tuple(interpolated_arrays)


def interpolation_1d(P_MPa, arrays):
    interpolated_arrays = []
    for array in arrays:
        # Create mask for known values (not NaN)
        nan_mask = np.isnan(array)
        # Extract known points and values
        x_known = P_MPa[~nan_mask]
        y_known = array[~nan_mask]
        spline = interpolate.make_interp_spline(x_known, y_known, k=2)
        # Perform the interpolation
        interpolated_results = spline(P_MPa)
        interpolated_arrays.append(interpolated_results)
    return tuple(interpolated_arrays)


def panda_df_generator(system):
    columns = ["P (MPa)", "T (K)", "pH", "Amount", "Charge", "SolidTotal", "LiquidTotal"]
    species = []
    for speciesItem in system.species():
        species.append(speciesItem.aggregateState().name + "Amount" + speciesItem.name())
    # Dictionary that maps certain column names that require using species list to Reaktoro code
    translation_dictionary_for_species_list = {"Amount": species}
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

    aqueous_species_list = 'H+ OH- H2O'
    speciation_ratio_mol_kg = {'H+': 1e-7, 'OH-': 1e-7, 'H2O': float(1/rkt.waterMolarMass)}
    frezchem = PhreeqcGeneratorForChemicalConstraint(aqueous_species_list, speciation_ratio_mol_kg, "mol", "frezchem.dat")
    rkt_T_freezing = rkt_t_freeze_without_pressure_correction(speciation_ratio_mol_kg, eos_P_MPa, frezchem, 100, 0, {})
    sfz_T_freezing = np.array(sfz_T_freezing)
    rkt_T_freezing = np.array(rkt_T_freezing, dtype = np.float64)
    difference_in_T_freezings =  sfz_T_freezing - rkt_T_freezing


    # Ensure the directory exists
    folder = 'Reaktoro_Saved_Files'
    current_directory = os.path.dirname(__file__)
    folder_path = os.path.join(current_directory, folder)
    os.makedirs(folder_path, exist_ok=True)

    dictionary_to_save = {"P_MPa": eos_P_MPa, "Difference_in_T_freezings_K": difference_in_T_freezings}
    with open(file_path_to_save, 'wb') as f:
        pickle.dump(dictionary_to_save, f)
    return dictionary_to_save


def rkt_t_freeze_without_pressure_correction(speciation_ratio_mol_kg, P_MPa, frezchem, PMax_MPa, freezing_temperature_spline, calculated_freezing_temperatures, TMin_K = 250, TMax_K = 300, significant_threshold = 0.1):
    """
     Calculates the temperature at which the prescribed aqueous solution freezes. Utilizes the reaktoro framework to
     constrain the equilibrium position at the prescribed pressure, the lower and upper limits of temperature (in K),
     and the total amount of ice at a significant threshold of 1e-14, therefore calculating and returning the
     temperature (within the range) at which ice begins to form.

     Parameters
     ----------
     speciation_ratio_mol_kg: the ratio of species in the aqueous solution in mol/kg of water. Should be a dictionary
     with the species as the key and its ratio as its value.
     P_MPa: the desired equilibrium freezing pressure(s).
     TMin_K: the lower limit of temperature that Reaktoro should query over
     TMax_K: the upper limit of temperature that Reaktoro should query over
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
                state = reset_state(system, speciation_ratio_mol_kg, "mol")
            # If the result failed, we will use spline
            else:
                log.warning(f"While attempting to find bottom freezing temperature for pressure of {P_MPa} MPa, \n"
                            +f"Reaktoro was unable to find a temperature within range of {TMin_K} K and {TMax_K}.\n"
                            +f"Instead, we will use a spline of freezing temperatures that we generated for this EOS to find the associated freezing temperature.")
                equilibrium_temperature = freezing_temperature_spline(P)
                state = reset_state(system, speciation_ratio_mol_kg, "mol")
        # If Pressure is >= PMax, Reaktoro is at its limits of computation and will likely fail. Thus, we utilize our spline approximation of freezing temperatures
        # over a range of pressures that we are confident that Frezchem works at and find the associated freezing temperature at the given pressure.
        # If the temperature we are querying over is less than that freezing temperature, then we assume its ice, otherwise we assume its liquid
        else:
            equilibrium_temperature = freezing_temperature_spline(P)
        freezing_temperatures.append(equilibrium_temperature)
        calculated_freezing_temperatures[P] = equilibrium_temperature
    # Return the equilibrium temperature
    return np.array(freezing_temperatures)


"""
NOT USED FUNCTIONS


def pressure_constraint(P_MPa, aqueous_species_list, speciation_ratio_mol_kg, database, dP = -1):
    Find the pressure constraint at which Reaktoro can find equilibrium for the given speciation and database. Starts at P_MPa and checks if rkt can find equilibrium
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
"""