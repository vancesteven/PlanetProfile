import copy
import reaktoro as rkt
import numpy as np
import logging
from PlanetProfile.Utilities.defineStructs import Constants
from PlanetProfile.Utilities.DataManip import ResetNearestExtrap
from scipy import interpolate
from PlanetProfile.Thermodynamics.Reaktoro.sigmaElectricMcCleskey2012 import elecCondMcCleskey2012
log = logging.getLogger('PlanetProfile')


def PhreeqcGenerator(aqueous_species_list, speciation_ratio_mol_kg, database):
    """ Create a Phreeqc Reaktoro System with the solid and liquid phase whose relevant species are determined by the provided aqueous_species_list.
        Works for both core10.dat and frezchem.dat.
    Args:
        aqueous_species_list: aqueous species in reaction. Should be formatted in one long string with a space in between each species
     speciation_ratio_mol_kg: the ratio of species in the aqueous solution in mol/kg of water. Should be a dictionary
     with the species as the key and its ratio as its value.
        database: frezchem or core10
    Returns:
        db, system, state, conditions, solver, props, ice_name: Relevant reaktoro objects
    """
    # Initialize the database
    db = rkt.PhreeqcDatabase(database)
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
    # Obtain ice name
    if database == "frezchem.dat":
        ice_name = "Ice(s)"
    else:
        ice_name = "Ice"
    # Return the Reaktoro objects that user will need to interact with
    return db, system, state, conditions, solver, props, ice_name


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
    return props.speciesAmount(ice_name)


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
    # Check if "H2O" is in the string (and thus dictionary), indicating it is not compatible with supcrt
    if "H2O" in aqueous_species_string:
        # Change "H2O" to "H2O(aq) in the string
        aqueous_species_list = aqueous_species_string.replace("H2O", "H2O(aq)")
        # Since python passes dictionary by reference, need to make deep copy to preserve original dictionary
        supcrt_speciation_ratio_mol_kg = copy.deepcopy(speciation_ratio_mol_kg)
        # Now change the "H2O" key in the dictionary to "H2O(aq)"
        supcrt_speciation_ratio_mol_kg["H2O(aq)"] = supcrt_speciation_ratio_mol_kg.pop("H2O")
        # Return the string and adapted dictionary
        return aqueous_species_list, supcrt_speciation_ratio_mol_kg
    else:
        # If "H2O" is not in the string (and thus dictionary), then just return them back
        return aqueous_species_string, speciation_ratio_mol_kg


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
    index = np.searchsorted(T_K, switch_T_K, side = 'left')
    frezchem_indices = range(0, index)
    core_indices = range(index, np.size(T_K))
    return frezchem_indices, core_indices


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
        db, system, state, conditions, solver, props, ice_name = PhreeqcGenerator(aqueous_species_list, speciation_ratio_mol_kg, database)
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
        db, system, state, conditions, solver, props, ice_name = PhreeqcGenerator(aqueous_species_list, speciation_ratio_mol_kg, database)
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
            state = reset_state(system, aqueous_species_list)
    # Return the adjusted T_K
    return T_K


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
    array_lists = [rho_kgm3, Cp_JKgK, alpha_pK]
    for array in array_lists:
        # Create mask of non-zero values
        nonzero_mask = (array != 0)
        # Interpolate using scipy interp2d function
        x, y = np.meshgrid(np.arange(array.shape[1]), np.arange(array.shape[0]))
        x_interp = x[nonzero_mask]
        y_interp = y[nonzero_mask]
        z_interp = array[nonzero_mask]
        f = interpolate.interp2d(x_interp, y_interp, z_interp, z_interp, kind='linear')
        for idx in zero_indices:
            if array[idx[0], idx[1]] == 0:
                array[idx[0], idx[1]] = f(idx[1], idx[0])
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
    log.warning("Performed 1d linear interpolation on missing seoismic properties")
    return array1, array2
