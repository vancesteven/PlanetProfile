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


def SupcrtGenerator(aqueous_species_list, speciation_ratio_per_kg, species_unit, database, consider_solid_phases, solid_phases_to_consider):
    """ Create a Supcrt Reaktoro System with the solid and liquid phase whose relevant species are determined by the provided aqueous_species_list.
    Args:
    aqueous_species_list: aqueous species in reaction. Should be formatted in one long string with a space in between each species
    speciation_ratio_per_kg: the ratio of species in the aqueous solution per kg of water. Should be a dictionary
        with the species as the key and its ratio as its value.
    species_unit: "mol" or "g" that species ratio is in
    database: Supcrt database to use
    consider_solid_phases: whether or not to consider solid phases in calculations
    solid_phases_to_consider: the solid phases to consider
    Returns:
        db, system, state, conditions, solver, props, ice_name: Relevant reaktoro objects
    """
    # Initialize the database
    db = rkt.SupcrtDatabase(database)
    # Prescribe the solution
    solution = rkt.AqueousPhase(aqueous_species_list)
    # If we are considering solid phases, create a solid phase
    if consider_solid_phases:
        if solid_phases_to_consider is None:
            solids = rkt.MineralPhases()
        else:
            # If we are specifying what phases to consider, note that we should only consider the phases that are relevant to the speices in the solution
            # I.e. don't consider all clathrates if they aren't possible to form (greatly decreases runtime if we do consider irrelevant phases)
            final_phases_to_consider = ''
            solids_tester = rkt.MineralPhases()
            system_tester = rkt.ChemicalSystem(db, solution, solids_tester)
            RelevantPhases = system_tester.species().withAggregateState(rkt.AggregateState.Solid)
            for solid in solid_phases_to_consider:
                if RelevantPhases.findWithName(solid) < RelevantPhases.size():
                    final_phases_to_consider = final_phases_to_consider + f' {solid}'
            solids = rkt.MineralPhases(final_phases_to_consider)
        system = rkt.ChemicalSystem(db, solution, solids)
    else:
        system = rkt.ChemicalSystem(db, solution)
    # Obtain all related solid phases
    # solids = rkt.MineralPhases()
    # Obtain all related gas phases
    # gases = rkt.GaseousPhase()
    # Initialize the system
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
    Importantly, since speciation_ratio_per_kg is a dictionary, we must make a deep copy before editing so as to not disturb the original dictionary, which
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


def freezing_temperature_correction_calculator():
    eos_P_MPa = np.linspace(0.1, 100, 150)
    eos_T_K = np.linspace(260, 280, 500)
    eos_PT = np.array([eos_P_MPa, eos_T_K], dtype=object)
    sfz_phases = sfz.whichphase(eos_PT)
    diff = np.diff(sfz_phases)
    sfz_T_freezing = []
    for row in diff:
        i = np.where(row == -1.0)[0]
        sfz_T_freezing.append(((eos_T_K[i] + eos_T_K[i+1])/2)[0])
    sfz_T_freezing = np.array(sfz_T_freezing)

    # Obtain frezchem freezing temperatures
    rkt_T_freezing = []
    aqueous_species_list = 'H+ OH- H2O'
    speciation_ratio_mol_kg = {'H2O': float(1/rkt.waterMolarMass)}
    frezchem_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Databases', 'frezchem.dat')
    frezchem = PhreeqcGeneratorForChemicalConstraint(aqueous_species_list, speciation_ratio_mol_kg, "mol", frezchem_file_path)
    db, system, state, conditions, solver, props = frezchem
    state = rkt.ChemicalState(state)
    # Create an iterator to go through P_MPa
    it = np.nditer([eos_P_MPa])
    conditions.set("IP", 0.1)
    conditions.setLowerBoundTemperature(240, "K")
    conditions.setUpperBoundTemperature(280, "K")
    for P in it:
        P = float(P)
        conditions.pressure(P, "MPa")
        # Solve the equilibrium problem
        result = solver.solve(state, conditions)
        # Update the properties
        props.update(state)
        # Obtain the equilibrium temperature
        rkt_T_freezing.append(float(props.temperature()))
        # Reset the state
        state = reset_state(system, speciation_ratio_mol_kg, "mol")
    rkt_T_freezing = np.array(rkt_T_freezing)

    # Find difference in freezing temperatures
    difference_in_T_freezings = sfz_T_freezing-rkt_T_freezing
    return eos_P_MPa, difference_in_T_freezings

def SupcrtH2OChemicalPotentialCorrectionSplineGenerator():
    """ Obtain the Supcrt16 aqueous H2O chemical potential differnce of pure water for Supcrt and Seafreeze across a 2-D grid of pressure and temperatures."""
    eos_P_MPa = np.linspace(0.1, 500, 500)
    eos_T_K = np.linspace(240, 400, 500)
    PT = np.array([eos_P_MPa, eos_T_K])
    out = sfz.getProp(PT, 'water1')
    sfz_chem_potential = out.G * rkt.waterMolarMass # Multiply Gibbs free energy by molar mass of H2O to get chemical potential of pure water

    rkt_chemical_potential = []
    aqueous_species_list = 'H+ OH- H2O(aq)'
    speciation_ratio_mol_kg = {'H2O(aq)': float(1/rkt.waterMolarMass)}
    supcrt = SupcrtGenerator(aqueous_species_list, speciation_ratio_mol_kg, "mol", "supcrt16")
    db, system, state, conditions, solver, props = supcrt
    P_MPa, T_K = np.meshgrid(eos_P_MPa, eos_T_K, indexing='ij')
    # Create a nditer iterator
    it = np.nditer([P_MPa, T_K], flags=['multi_index'])
    # Go through each P, T combination
    for P, T in it:
        P = float(P)
        T = float(T)
        conditions.temperature(T, "K")
        # Establish equilibrium pressure constraint value
        conditions.pressure(P, "MPa")
        # Solve the equilibrium problem
        result = solver.solve(state, conditions)
        # Update the properties
        props.update(state)
        # Check if the equilibrium problem succeeded
        if result.succeeded():
            rkt_chemical_potential.append(float(props.speciesChemicalPotential('H2O(aq)')))
        else:
            print("HELLO")
    rkt_chemical_potential = np.array(rkt_chemical_potential).reshape(P_MPa.shape)

    # Find difference in chemical potentials
    difference_in_chemical_potential = rkt_chemical_potential-sfz_chem_potential

    return eos_P_MPa, eos_T_K, difference_in_chemical_potential





"""
NOT USED FUNCTIONS


def pressure_constraint(P_MPa, aqueous_species_list, speciation_ratio_per_kg, database, dP = -1):
    Find the pressure constraint at which Reaktoro can find equilibrium for the given speciation and database. Starts at P_MPa and checks if rkt can find equilibrium
        with a temperature of 273K. If it cannot, then adjusts P_MPa by dP and tries again, continuing this process until a compatible P_MPa is found.
        If we are looking for upper pressure cosntraint (dP is negative), then we stop also if it reaches <= 1.1 MPa, which Rkt should be compatible with.
    Args:
        P_MPa: Initial pressure constraint in MPa
        aqueous_species_list: String that has all species names that should be considered in aqueous phase
        speciation_ratio_per_kg: Dictionary of active species and the values of their molar ratio (mol/kg of water)
        database: Database to find pressure constraint for
        dP: The amount to change P_MPa by if equilibrium is not achieved. Defaults to -1 (for upper constraint)

    Returns:
        P_MPa: New pressure constraint in MPa for given database.

    # Disable chemical convergence warnings that Reaktoro raises. We handle these internally instead and throw more specific warnings when they appear.
    rkt.Warnings.disable(906)
    # Initilialize the database, either being supcrt of phreeqc
    if "supcrt" in database:
        # Since supcrt labels "H2O(aq)", we need to adjust our species list and dictionary accordingly
        aqueous_species_list, speciation_ratio_per_kg = species_convertor_compatible_with_supcrt(
            aqueous_species_list, speciation_ratio_per_kg)
        db, system, state, conditions, solver, props = SupcrtGenerator(aqueous_species_list, speciation_ratio_per_kg, database)
    else:
        db, system, state, conditions, solver, props, ice_name, database_name = PhreeqcGenerator(aqueous_species_list, speciation_ratio_per_kg, database)
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
            state = reset_state(system, speciation_ratio_per_kg)
            P_MPa += dP
    # Return the adjusted P_MPa
    return P_MPa
"""