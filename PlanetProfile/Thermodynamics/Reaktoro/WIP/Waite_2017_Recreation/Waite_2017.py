from reaktoro import *
import matplotlib.pyplot as plt
import numpy as np
import csv

# Load the thermodynamic database
db = SupcrtDatabase("supcrt16-organics")

# Define constants of CH4 and H2 concentrations from CO2 ratios, as defined in supplementary material
ratio_CH4_CO2 = 0.4
ratio_H2_CO2 = 1.6

# Model setups
P_MPa = [0.1]  # Pressure range in MPa
T_K = [273]  # Temperature range in K

# Enable LaTeX rendering and add mhchem to the preamble for plotting
stix = 'stix'
mhchem = 'mhchem'
siunitx = 'siunitx'
plt.rcParams.update({
    'text.usetex': True,
    'text.latex.preamble': f'\\usepackage{{{stix}}}\\usepackage{{{mhchem}}}\\usepackage{{{siunitx}}}',
    'font.family': 'sans-serif',  # Changed from 'serif' to 'sans-serif'
    'font.sans-serif': 'Arial',    # Set to Arial
    'axes.labelsize': 16,      # Axis label font size
    'xtick.labelsize': 14,     # X tick label font size
    'ytick.labelsize': 14,     # Y tick label font size
    'axes.titlesize': 18,
        'font.weight': 'bold',           # Makes all text bold
    'axes.labelweight': 'bold',      # Makes axis labels bold
    'axes.titleweight': 'bold', 
})

def TableS11_recreation():
    """
    Function to recreate the supplementary table S11 from Waite 2017 paper.
    This function simulates equilibrium calculations at pH 9 and 11 based on table S11 modeling setup.
    Outputs are saved to a text file.
    """
    # Define the range of pH values to be considered
    pH = [9, 11]

    # Define the initial concentrations of aqueous species (in mol/kg H2O)
    aqueous_species = {"CO3-2": 0.015, "H+": 0, "OH-": 0, "Cl-": 0.1, "HCO3-": 0.015, "CO2(aq)": 0, "H2O(aq)": 55.51}

    # Generate the SUPCRT system for equilibrium calculations
    db, system, initial_state, conditions, solver, props = generate_supcrt(aqueous_species=aqueous_species,
        use_activity_model=True)

    # Clone the initial state to modify during calculations
    state = initial_state.clone()

    # Open a file to save results
    with open("TableS11_Output.txt", "w") as file:
        # Iterate through each pressure value
        for P in P_MPa:
            conditions.pressure(P, "MPa")  # Set pressure condition

            # Iterate through each temperature value
            for T in T_K:
                conditions.temperature(T, "K")  # Set temperature condition

                # Iterate through each pH value
                for ph in pH:
                    conditions.pH(ph)  # Set pH condition

                    # Update state with current conditions
                    state.setPressure(P, "MPa")
                    state.setTemperature(T, "K")

                    # Solve the equilibrium problem
                    result = solver.solve(state, conditions)

                    if result.succeeded():
                        # If solution is successful, update properties
                        props.update(state)
                        aprops = AqueousProps(props)

                        # Retrieve molality of CO2(aq)
                        CO2_molal = float(aprops.speciesMolality("CO2(aq)"))

                        # Compute molality of CH4(aq) and H2(aq) based on defined ratios
                        CH4_molal = CO2_molal * ratio_CH4_CO2
                        H2_molal = CO2_molal * ratio_H2_CO2

                        # Save the results to file
                        file.write(f"At P = {P} MPa, T = {T} K, pH = {ph}, [CO2(aq)] is {CO2_molal}, "
                                   f"[CH4(aq)] is {CH4_molal}, [H2(aq)] is {H2_molal}\n")
                    else:
                        # Save an error message to file if the solver fails
                        file.write(f"ERROR at P = {P} MPa and T = {T} K\n")

                    # Reset state to the initial conditions before next iteration
                    state = initial_state.clone()
                    props.update(state)

    print("Output saved to 'TableS11_Output.txt'")


def Fig4_recreation():
    """
    Function to recreate Figure 4 from Waite 2017 paper.
    This function calculates reaction quotients for methanogenesis at different pH and H2/H2O mixing ratios, based on Figure 4 modeling setup.
    """
    # Define the pH range for analysis
    pH_list = [8, 9, 10, 11, 12, 13]

    # Define H2/H2O ratios in the system (x-axis range)
    ratio_H2_H2O = [0.00001, 0.0001, 0.001, 0.005, 0.009, 0.014, 0.1]

    # Define the CO2/H2O mixing ratio from Table 1
    CO2_H2O_mixing_ratio = 0.0055

    # Convert H2/H2O ratios to H2/CO2 ratios by dividing by the CO2_H2O_mixing_ratio
    ratio_H2_CO2 = [ratio / CO2_H2O_mixing_ratio for ratio in ratio_H2_H2O]

    # Define the methanogenesis reaction: CO2(aq) + 4H2(aq) -> CH4(aq) + 2H2O(aq)
    reaction = {"products": {"Methane(aq)": 1, "H2O(aq)": 2}, "reactants": {"CO2(aq)": 1, "H2(aq)": 4}}

    # Define function to calculate CO2 molality based on equation S18
    CO2_mol_equation = lambda ph: -0.1213 * (ph ** 2) + 0.9832 * ph - 3.1741

    # Define initial concentrations of gaseous and aqueous species
    gaseous_species = {"H2(g)": 0, "CO2(g)": 0, "CH4(g)": 0}
    aqueous_species = {"CO3-2": 0.03, "Na+": 0, "H2(aq)": 0, "H+": 0, "OH-": 0, "Cl-": 0.1, "HCO3-": 0.0, "CO2(aq)": 0,
                       "H2O(aq)": 55.51}

    # Initialize a dictionary to store affinity data for each pH level
    pH_affinity_data = {ph: [] for ph in pH_list}

    # Iterate through pressure and temperature conditions
    for P in P_MPa:
        for T in T_K:
            for idx, ph in enumerate(pH_list):
                # Calculate the molality of CO2 based on pH
                CO2_mol = 10 ** CO2_mol_equation(ph)
                aqueous_species["CO2(aq)"] = CO2_mol

                # Calculate the molality of CH4 based on CO2 amount and reaction ratio
                CH4_mol = CO2_mol * ratio_CH4_CO2
                aqueous_species["Methane(aq)"] = CH4_mol

                # Iterate through each H2/CO2 ratio
                for index, ratio in enumerate(ratio_H2_CO2):
                    corresponding_ratio_H2_H2O = ratio_H2_H2O[index]  # Get associated H2/H2O ratio

                    # Calculate the molality of H2
                    H2_mol = CO2_mol * ratio
                    aqueous_species["H2(aq)"] = H2_mol

                    # Generate SUPCRT database for equilibrium calculations
                    db, system, state, conditions, solver, props = generate_supcrt(aqueous_species, gaseous_species=gaseous_species)

                    # Set experimental conditions
                    conditions.pressure(P, "MPa")
                    conditions.temperature(T, "K")
                    conditions.pH(ph)
                    state.setPressure(P, "MPa")
                    state.setTemperature(T, "K")
                    props.update(state)

                    # Calculate disequilibrium reaction quotient (Q) based on current properties
                    Q = calculate_reaction_quotient(props, reaction)

                    # Solve equilibrium problem
                    result = solver.solve(state, conditions)
                    if result.succeeded():
                        # Update properties post-solution
                        props.update(state)

                        # Calculate equilibrium constant (K) and affinity (A)
                        K = calculate_reaction_quotient(props, reaction)
                        R = 8.31446  # Gas constant in J/(mol·K)
                        A = 2.3026 * R * T * (np.log10(K) - np.log10(Q)) / 1000  # Affinity in kJ

                        # Store affinity data for plotting
                        pH_affinity_data[ph].append((corresponding_ratio_H2_H2O, A))
                    else:
                        print(f"ERROR at P = {P} MPa and T = {T} K")

    # Plot the calculated affinity values against H2/H2O ratios
    plt.figure(figsize=(8.48, 6.727))
    for ph in pH_list:
        ratios = [data[0] for data in pH_affinity_data[ph]]
        affinities = [data[1] for data in pH_affinity_data[ph]]
        plt.plot(ratios, affinities, label=f'pH = {ph}', color='blue')
        plt.text(ratios[-1], affinities[-1], f'pH {ph}', fontsize=10, color='blue', ha='left', va='center')

    # Configure plot labels and scales
    plt.ylabel(r'Apparent Affinity (A) $[kJ \ (mol \ CH_4)^{-1}]$')
    plt.xlabel(r'$H_2/H_2O$ Mixing Ratio in Plume Gas')
    plt.xscale('log')

    # Add observed region indicators
    plt.axvline(x=0.005, color='orange', linestyle='-', label='Observed region', zorder=5)
    plt.axvline(x=0.014, color='orange', linestyle='-', label='Observed region', zorder=5)
    plt.text(0.0085, 160, 'Observed', ha='center', fontsize=14, color='orange')
    plt.axhline(y=0, color='red', linestyle='--')
    plt.text(0.0001, 5, 'Equilibrium', color='red', fontsize=16, ha='left')

    # Highlight the observed region between pH 9 and pH 11
    ratios_pH9 = np.array([data[0] for data in pH_affinity_data[9]])
    affinity_pH9 = np.array([data[1] for data in pH_affinity_data[9]])
    affinity_pH11 = np.array([data[1] for data in pH_affinity_data[11]])
    mask = (ratios_pH9 >= 0.005) & (ratios_pH9 <= 0.014)
    plt.fill_between(ratios_pH9[mask], affinity_pH9[mask], affinity_pH11[mask], color='skyblue', alpha=1,
                     label='Observed region', zorder=3)

    # Adjust plot limits and ticks
    plt.ylim(-40, 170)
    plt.xlim(0.00001, 0.1)
    plt.yticks(np.arange(-40, 161, 40))

    # Remove unnecessary spines for cleaner visualization
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Save and display the plot
    plt.savefig("Fig4_Recreation.png", dpi=600, bbox_inches='tight')
    plt.show()


def calculate_reaction_quotient(prop, reaction):
    """
    Function to calculate the reaction quotient (Q) for a given chemical reaction.
    The reaction quotient is computed as the ratio of the activities of products to reactants,
    each raised to the power of their respective stoichiometric coefficients.

    :param prop: Object containing species activity data.
    :param reaction: Dictionary defining reactants and products with their stoichiometric coefficients.
    :return: Reaction quotient (Q).
    """
    # Initialize numerator and denominator for Q
    Q_numerator = 1.0
    Q_denominator = 1.0

    # Compute the numerator by multiplying the activities of the products
    # raised to their respective stoichiometric coefficients
    for species, coefficient in reaction["products"].items():
        speciesActivity = float(prop.speciesActivity(species))
        Q_numerator *= speciesActivity ** coefficient

    # Compute the denominator by multiplying the activities of the reactants
    # raised to their respective stoichiometric coefficients
    for species, coefficient in reaction["reactants"].items():
        speciesActivity = float(prop.speciesActivity(species))
        Q_denominator *= speciesActivity ** coefficient

    # Calculate the reaction quotient Q as the ratio of numerator to denominator
    Q = Q_numerator / Q_denominator
    return Q


def generate_supcrt(aqueous_species, gaseous_species=None, use_activity_model=False):
    """
    Function to generate a SUPCRT-based chemical system for equilibrium calculations.
    This function defines the aqueous and gaseous phases, sets up equilibrium conditions,
    and initializes the solver for the system.

    :param aqueous_species: Dictionary of aqueous species with their initial molalities.
    :param gaseous_species: (Optional) Dictionary of gaseous species with their initial molalities.
    :param use_activity_model: Boolean indicating whether to use an activity model for the solution.
    :return: Tuple containing the database, system, state, conditions, solver, and properties.
    """
    # Convert aqueous species keys to a space-separated string for defining the aqueous phase
    aqueous_species_str = " ".join(aqueous_species.keys())

    # Define the chemical system with an aqueous phase
    solution = AqueousPhase(aqueous_species_str)

    # Apply an activity model if specified
    if use_activity_model:
        solution.setActivityModel(chain(ActivityModelPitzer(), ActivityModelPhreeqcIonicStrengthPressureCorrection()))

    # Define gaseous phase if gaseous species are provided
    if gaseous_species is not None:
        gaseous_species_str = " ".join(gaseous_species.keys())
        gases = GaseousPhase(gaseous_species_str)
        system = ChemicalSystem(db, solution, gases)
    else:
        system = ChemicalSystem(db, solution)

    # Define the equilibrium state
    state = ChemicalState(system)

    # Add gaseous species to the state if provided
    if gaseous_species is not None:
        for species, amount in gaseous_species.items():
            state.add(species, amount, "mol")

    # Add aqueous species to the state
    for species, amount in aqueous_species.items():
        state.add(species, amount, "mol")

    # Set equilibrium specifications including temperature, pressure, and pH
    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.pH()
    conditions = EquilibriumConditions(specs)

    # Calculate the aqueous phase properties
    props = ChemicalProps(state)

    # Initialize solver for equilibrium calculations
    solver = EquilibriumSolver(specs)

    # Update properties for the initial state
    props.update(state)

    # Return the necessary components for further equilibrium calculations
    return db, system, state, conditions, solver, props

if __name__ == '__main__':
    TableS11_recreation()
    Fig4_recreation()
