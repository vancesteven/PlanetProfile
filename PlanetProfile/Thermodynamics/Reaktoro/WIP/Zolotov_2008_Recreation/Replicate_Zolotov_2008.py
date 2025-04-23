from reaktoro import *
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import LogFormatterExponent
from scipy.io import savemat
from adjustText import adjust_text

# Define ROOT (for loading files)
_ROOT = os.path.dirname(os.path.abspath(__file__))

# Load the Core11 database into Reaktoro
path_file = os.path.join(_ROOT, "core11_Diab_2023.dat")
db = PhreeqcDatabase.fromFile(path_file)

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

# Modeling setup
P_MPa = [137.5]  # Pressure in MPa
T_K = [273.2]  # Temperature in K
mass_g = 1000  # Mass of CI Chondrite to add in grams
# W/R ratio to use
water_rock_ratio = 1
# Elements to have in chemical system
elements = "O H C S Cl Na K Mg Fe Ca Si Al"

def Replicate_Zolotov_H2():
    """
    Function to run to replicate the Zolotov H2 fugacity model using the Core11 database and plot the results.
    """
    ci_chondrite_composition = {'aqueous species': {}, 'mineral species': {'CI-Chondrite-Zolotov-TenPercent': mass_g}} # Define the initial state's aqueous and mineral species

    # Log H2 fugacities to query over
    log_H2_fugacity = np.linspace(-12, -3, 100)

    # Call function to generate the Reaktoro chemical system and all associated objects we need
    system, initial_state, conditions, solver, restrictions = generate_chemical_system(db, ci_chondrite_composition, True, DO_speciate=True, elements=elements, DO_PH=False, DO_CO2 = False, DO_H2 = True)

    #Prepare lists for all the information we care about to know at equilibrium (pH, elemental concentration, solid species, etc.)
    pH_list = []
    element_names = np.array([element.symbol() for element in system.elements()])
    element_list = [[] for _ in range(element_names.size)]
    aqueous_species_names = np.array([species.name() for species in system.species() if species.aggregateState() == AggregateState.Aqueous])
    aqueous_species_list = [[] for _ in range(aqueous_species_names.size)]
    solid_species_names = np.array([species.name() for species in system.species() if species.aggregateState() == AggregateState.Solid])
    solid_species_mol_list = [[] for _ in range(solid_species_names.size)]
    solid_phases_names = np.array([phase.name() for phase in system.phases() if phase.aggregateState() == AggregateState.Solid])
    solid_phases_volume_list = [[] for _ in range(solid_phases_names.size)]

    # Add mass_g * water_rock_ratio (1kg) to the state
    initial_state.add("H2O", mass_g * water_rock_ratio, "g")

    # The Zolotov model separately adds Cl- as it assumes 'The majority of Cl was extracted from accreted rocks' rather than the CI chondrite, so we also add a determined Cl- amount
    initial_state.add("Cl-", 0.12, "mol")

    # More reaktoro wrangling that includes cloning our initial state and defining a properties object
    state = initial_state.clone()
    props = ChemicalProps(state)

    # Define the Pressure and Temperature conditions
    for P in P_MPa:
        conditions.pressure(P, "MPa")
        for T in T_K:
            conditions.temperature(T, "K")
            # Go through each H2 fugacity in the array
            for H2 in log_H2_fugacity:
                conditions.fugacity("H2(g)", 10**H2, "bar")
                # Set the initial state's equilibrium pressure and temperature
                state.setPressure(P, "MPa")
                state.setTemperature(T, "K")

                # Here we increase the number of iterations Reaktoro will do while looking for thermodynamic equilibrium
                options = EquilibriumOptions()
                options.optima.maxiters = 20000
                solver.setOptions(options)

                # Solve the equilibrium problem
                result = solver.solve(state, conditions)

                # If we did not succeed, then let's try again from the initial state
                if not result.succeeded():
                    state = initial_state.clone()
                    result = solver.solve(state, conditions)
                # If we did succeed, then save all the equilibrium data to the previously created lists
                if result.succeeded():
                    # Update the properties
                    props.update(state)
                    aprops = AqueousProps(props)
                    # Save all the equilibrium data to the previously created lists
                    pH = float(aprops.pH())
                    pH_list.append(float(aprops.pH()))
                    for k, species_molality in enumerate(aprops.elementMolalities()):
                        element_list[k].append(float(f'{float(species_molality):.3e}'))
                    for k, species_molality in enumerate(aprops.speciesMolalities()):
                        aqueous_species_list[k].append(float(f'{float(species_molality):.3e}'))
                    for k, species in enumerate(solid_species_names):
                        solid_species_mol_list[k].append(float(props.speciesAmount(species)))
                    for k, species in enumerate(solid_phases_names):
                        solid_phases_volume_list[k].append(float(props.phaseProps(species).volume())*100**3)
                # If we did not succeed again, then we should print error and fill our lists accordingly with empty data
                else:
                    aprops = AqueousProps(props)
                    pH_list.append(np.nan)
                    for k, species_molality in enumerate(aprops.elementMolalities()):
                        element_list[k].append(np.nan)
                    for k, species_molality in enumerate(aprops.speciesMolalities()):
                        aqueous_species_list[k].append(np.nan)
                    for k, species in enumerate(solid_species_names):
                        solid_species_mol_list[k].append(np.nan)
                    for k, species in enumerate(solid_phases_names):
                        solid_phases_volume_list[k].append(np.nan)
                    props.update(state)
                    # Print error
                    print(f"ERROR at P = {P} MPa and T = {T} K, H2 = {H2}")

    # Convert lists to arrays
    pH_array = np.array(pH_list)
    aqueous_species_array_molal = np.array(aqueous_species_list)
    element_species_array_molal = np.array(element_list)
    solid_species_array_mol = np.array(solid_species_mol_list)
    solid_phases_volume_array = np.array(solid_phases_volume_list)

    # Define the axis limits of the plots to generate
    species_lim = (1e-6, 1e0)
    pH_lim = (6, 12.5)

    # Save the data if desired
    # save_aqueous_species_mat(aqueous_species_names, aqueous_species_array_molal)

    # Generate plots
    figname = os.path.join(_ROOT, 'Zolotov_H2_Replication')
    generate_Zolotov_plots(data=(element_names, aqueous_species_names, solid_phases_names, pH_array, aqueous_species_array_molal, element_species_array_molal, solid_species_array_mol, solid_phases_volume_array), x=log_H2_fugacity, species_lim = species_lim, pH_lim = pH_lim, figname=figname)




def Replicate_Zolotov_Core11_Melwani_CO2():
    """
    Function to run to replicate the Zolotov CO2 fugacity model using the Core11 database and plot the results.
    """
    ci_chondrite_composition = {'aqueous species': {}, 'mineral species': {'CI-Chondrite-Zolotov': mass_g}} # Define the initial state's aqueous and mineral species
    # Define the log_H2_fugacity
    log_H2_fugacity = -10

    # Log CO2 fugacities to query over
    log_CO2_fugacity = np.linspace(-6, 2, 100)

    # Call function to generate the Reaktoro chemical system and all associated objects we need
    system, initial_state, conditions, solver, restrictions = generate_chemical_system(db, ci_chondrite_composition, True, DO_speciate=True, elements=elements, DO_PH=False)

    #Prepare lists for pH and species amounts
    pH_list = []
    element_names = np.array([element.symbol() for element in system.elements()])
    element_list = [[] for _ in range(element_names.size)]
    aqueous_species_names = np.array([species.name() for species in system.species() if species.aggregateState() == AggregateState.Aqueous])
    aqueous_species_list = [[] for _ in range(aqueous_species_names.size)]
    solid_species_names = np.array([species.name() for species in system.species() if species.aggregateState() == AggregateState.Solid])
    solid_species_mol_list = [[] for _ in range(solid_species_names.size)]
    solid_phases_names = np.array([phase.name() for phase in system.phases() if phase.aggregateState() == AggregateState.Solid])
    solid_phases_volume_list = [[] for _ in range(solid_phases_names.size)]
    initial_state.add("H2O", mass_g * water_rock_ratio, "g")
    initial_state.add("Cl-", 0.12, "mol")
    state = initial_state.clone()

    # Create props object
    props = ChemicalProps(state)
    # Go through each pressure and temperature combination
    for P in P_MPa:
        conditions.pressure(P, "MPa")
        for T in T_K:
            for CO2 in log_CO2_fugacity:
                # set the conditions
                conditions.temperature(T, "K")
                state.setPressure(P, "MPa")
                state.setTemperature(T, "K")
                conditions.fugacity("H2(g)", 10**log_H2_fugacity, "bar")
                conditions.fugacity("CO2(g)", 10**CO2, "bar")
                options = EquilibriumOptions()
                options.optima.maxiters = 20000
                solver.setOptions(options)

                # Solve the equilibrium problem
                result = solver.solve(state, conditions)
                if not result.succeeded():
                    state = initial_state.clone()
                    result = solver.solve(state, conditions)
                if result.succeeded():
                    # Update the properties
                    props.update(state)
                    aprops = AqueousProps(props)

                    pH_list.append(float(aprops.pH()))
                    for k, species_molality in enumerate(aprops.elementMolalities()):
                        element_list[k].append(float(species_molality))
                    for k, species_molality in enumerate(aprops.speciesMolalities()):
                        aqueous_species_list[k].append(float(species_molality))
                    for k, species in enumerate(solid_species_names):
                        solid_species_mol_list[k].append(float(props.speciesAmount(species)))
                    for k, species in enumerate(solid_phases_names):
                        solid_phases_volume_list[k].append(float(props.phaseProps(species).volume())*100**3)
                else:
                    aprops = AqueousProps(props)
                    pH_list.append(np.nan)
                    for k, species_molality in enumerate(aprops.elementMolalities()):
                        element_list[k].append(np.nan)
                    for k, species_molality in enumerate(aprops.speciesMolalities()):
                        aqueous_species_list[k].append(np.nan)
                    for k, species in enumerate(solid_species_names):
                        solid_species_mol_list[k].append(np.nan)
                    for k, species in enumerate(solid_phases_names):
                        solid_phases_volume_list[k].append(np.nan)
                    props.update(state)
                    # Print error
                    print(f"ERROR at P = {P} MPa and T = {T} K, CO2 = {CO2}")
    # Convert lists to arrays
    pH_array = np.array(pH_list)
    aqueous_species_array_molal = np.array(aqueous_species_list)
    element_species_array_molal = np.array(element_list)
    solid_species_array_mol = np.array(solid_species_mol_list)
    solid_phases_volume_array = np.array(solid_phases_volume_list)

    # Define the axis limits of the plots to generate
    species_lim = (1e-6, 1e1)
    pH_lim = (4, 8)

    # Save the data if desired
    # save_aqueous_species_mat(aqueous_species_names, aqueous_species_array_molal)

    # Generate plots
    figname = os.path.join(_ROOT, 'Zolotov_CO2_Replication')
    generate_Zolotov_plots(data=(element_names, aqueous_species_names, solid_phases_names, pH_array, aqueous_species_array_molal, element_species_array_molal, solid_species_array_mol, solid_phases_volume_array), x=log_CO2_fugacity, species_lim = species_lim, pH_lim = pH_lim, figname = figname, exclude_species={'HS-+H2S'})




def generate_Zolotov_plots(data, x, species_lim, pH_lim, figname, exclude_species=None):
    """
    Function to generate plots, modeled after Zolotov's
    :param data: data to plot
    :param x: x range to plot
    :param species_lim: species axis limits
    :param pH_lim: pH axis limits
    :param figname: output filename
    :param exclude_species: set of species names to exclude from plotting (optional)
    """
    if exclude_species is None:
        exclude_species = set()
    # Define the elements and aqeuous species we want to plot as lines
    lines = {
     ChemicalSpecies('tab:cyan', 'Fe', 'Fe'),
     ChemicalSpecies('tab:orange', 'Ca', 'Ca'),
    ChemicalSpecies('tab:green', 'K', 'K'),
     ChemicalSpecies('tab:red', 'Na', 'Na'),
     ChemicalSpecies('tab:red', 'Cl', 'Cl'),
    ChemicalSpecies('tab:purple', 'SO4^{2-}', 'SO4'),
        ChemicalSpecies('tab:olive', 'Mg', 'Mg'),
     ChemicalSpecies('tab:blue', 'Si', 'Si'),
    ChemicalSpecies('tab:olive', 'C', 'C'),
        ChemicalSpecies('tab:olive', 'HS^{-}+H_{2}S', 'HS-+H2S')
    }
    # Filter out excluded species
    desired_species_to_plot = [s for s in lines if s.species not in exclude_species]
    # Get data
    element_names, aqueous_species_names, solid_phases_names, pH_array, aqueous_species_array_molal, element_species_array_molal, solid_species_array_mol, solid_phases_volume_array = data

    # Set up the plot
    # Create subplots (one on top and one below)
    fig = plt.figure(figsize=(8,10))
    grid = GridSpec(7, 1)
    ax1 = fig.add_subplot(grid[0:3, 0])
    ax2 = fig.add_subplot(grid[3, 0])
    ax3 = fig.add_subplot(grid[4:7, 0])
    texts_ax1 = []
    texts_ax3 = []
    for species in desired_species_to_plot:
        if species.species == 'SO4':
            index_species = np.where(np.char.find(aqueous_species_names, species.species) != -1)
        elif species.species == 'HS-+H2S':
            index_species = np.where((np.char.find(aqueous_species_names, 'HS-') != -1) | (np.char.find(aqueous_species_names, 'H2S') != -1))

        else:
            index_species = np.where(aqueous_species_names == species.species)
        index_element = np.where(element_names == species.species)
        if index_species[0].size > 0:
            Plot = True
            species_molality = np.zeros(np.shape(aqueous_species_array_molal[index_species[0][0]]))
            for index in index_species[0]:
                species_molality += aqueous_species_array_molal[index]
        elif index_element[0].size > 0:
            Plot = True
            species_molality = element_species_array_molal[index_element[0][0]]
        else:
            Plot = False

        if Plot:
            # Find the index of the maximum y value
            max_index = np.argmax(species_molality)
            line,  = ax1.plot(x, species_molality, label=species.species, color = species.color)
            text_obj = ax1.text(x[max_index], species_molality[max_index], rf'$\ce{{{species.mhchem_name}}}$', ha='left',
                       va='top', color=line.get_color(), fontsize=14)
            texts_ax1.append(text_obj)
    # Add labels and title
    ax1.set_yscale('log')
    # Apply LogFormatterExponent to show only exponents
    formatter = LogFormatterExponent(base=10)
    ax1.yaxis.set_major_formatter(formatter)
    ax1.set_xticklabels([]) # hide labels but keep ticks
    y_limit_min, y_limit_max = species_lim  # Adjust this factor as needed; 10 is an example value for the log scale
    # Set y-limits for aqueous species plots
    ax1.set_ylim(bottom=y_limit_min, top = y_limit_max)
    ax1.set_ylabel(r'Log Mole (kg $\ce{H2O^{-1}}$)')
    xmax = np.max(x)
    xmin = np.min(x)
    ax2.set_xlim(xmin, xmax)
    ax1.set_xlim(xmin, xmax)
    ax2.set_ylim(6, 12.5)
    ax2.set_yticks([6, 8, 10, 12])  # Set desired tick locations
    ax2.set_yticklabels(['6', '8', '10', '12'])  # Align text to left
    ax2.set_xticklabels([]) # hide labels but keep ticks

    # Plot the lower plot: pH vs log CO2 fugacity
    ax2.plot(x, pH_array, color='black', linestyle='-', label='pH')
    ax2.set_ylabel('pH')
    ax3.set_xlabel(r'Log $\text{H}_2$ Fugacity')
    desired_solid_species_to_plot = solid_phases_names
    for i, species_name in enumerate(desired_solid_species_to_plot):
        species_molality = solid_phases_volume_array[i]
        Plot = False
        if np.any(species_molality > 1e-14):
            Plot = True
        # Plot the phase plot
        if Plot:
            # Find the index of the maximum y value
            max_index = np.argmax(species_molality)
            line, = ax3.plot(x, species_molality, label=species_name)
            text_obj = ax3.text(x[max_index], species_molality[max_index], species_name, ha='left', va='top',
                     color=line.get_color(), fontsize=10)
            texts_ax3.append(text_obj)
    y_limit_min, y_limit_max = pH_lim
    ax2.set_ylim(bottom=y_limit_min, top = y_limit_max)
    ax3.set_ylabel(r"Vol. of minerals, $\mathrm{cm^{3}}/\mathrm{kg}$ rock")
    # Set x-axis limits and ticks for ax3
    x_start = float(np.min(x))
    x_end = float(np.max(x))
    ax3.set_xlim(x_start, x_end)
    # Generate ticks every 2 units from min(x) up to max(x) (not exceeding max)
    tick_start = int(np.ceil(x_start))
    tick_end = int(np.floor(x_end))
    ticks = list(range(tick_start, tick_end + 1, 2))
    if ticks and ticks[-1] > x_end:
        ticks = ticks[:-1]
    ax3.set_xticks(ticks)
    # After plotting everything, synchronize x-ticks across all axes
    ax1.set_xticks(ticks)
    ax2.set_xticks(ticks)
    # Adjust text to prevent overlap and keep within bounds
    adjust_text(texts_ax1, ax=ax1, only_move={'points':'y', 'text':'xy'}, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
    adjust_text(texts_ax3, ax=ax3, only_move={'points':'y', 'text':'xy'}, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
    plt.tight_layout()
    plt.show()

    # Save figure
    fig.savefig(figname, dpi = 600)






def generate_chemical_system(db, mineral_species, use_activity_model=False, DO_speciate=False,
                             elements=None, DO_PH=True, DO_CO2 = True, DO_H2 = True):
    # Define the aqueous phase by speciating the elements
    solution = AqueousPhase(speciate(elements))


    # Define the solid solution phases in the system
    MgNaCaKSaponite = MineralPhase("Saponite-Mg-Mg Saponite-Mg-Na Saponite-Mg-K Saponite-Mg-Ca")
    MgNaCaKSaponite.setName("Mg-Na-Ca-K-saponite")
    Montmor = MineralPhase("Montmor-Ca Montmor-K Montmor-Mg Montmor-Na")
    Montmor.setName("Na-Mg-K-Ca-montmorillonite")

    # Define the mineral phases
    sepiolite = MineralPhase("Sepiolite")
    Chrysotile = MineralPhase("Chrysotile")
    goethite = MineralPhase("Goethite")
    amorphous_silica = MineralPhase("SiO2(am)")
    dolomite = MineralPhase("Dolomite")
    calcite = MineralPhase("Calcite")
    andradite = MineralPhase("Andradite")
    Sylvite = MineralPhase("Sylvite")
    gypsum = MineralPhase("Gypsum")
    magnetite = MineralPhase("Magnetite")
    pyrrhotite = MineralPhase("Pyrrhotite")
    magnesite = MineralPhase("Magnesite")
    pyrite = MineralPhase("Pyrite")
    daphnite = MineralPhase('Daphnite-14A')
    siderite = MineralPhase("Siderite")
    minerals = MineralPhases(" ".join(mineral_species['mineral species'].keys()))

    # Define the gas phases by speciating the elements
    gases = GaseousPhase(speciate(elements))

    # Use pitzer and peng robinson activity models
    if use_activity_model:
        solution.set(chain(ActivityModelPitzer(), ActivityModelPhreeqcIonicStrengthPressureCorrection()))
        gases.set(ActivityModelPengRobinson())

    # Create the system
    system = ChemicalSystem(db, solution, gases, minerals, MgNaCaKSaponite, goethite, amorphous_silica, gypsum, daphnite, magnesite, siderite, Montmor, pyrite, dolomite, calcite, Chrysotile, Sylvite, andradite, pyrrhotite, magnetite, sepiolite)

    # Define the state
    state = ChemicalState(system)

    # Add the mineral species to the state
    for species, amount in mineral_species['mineral species'].items():
        state.add(species, amount, "g")

    # Set equilibrium specifications needed for Reaktoro
    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    if DO_H2:
        specs.fugacity("H2(g)")
    if DO_CO2:
        specs.fugacity("CO2(g)")
    if DO_PH:
        specs.pH()
    conditions = EquilibriumConditions(specs)
    restrictions = EquilibriumRestrictions(system)

    # Initialize solver
    solver = EquilibriumSolver(specs)

    # Return the system and objects needed to interact with the system
    return system, state, conditions, solver, restrictions


class ChemicalSpecies:
    def __init__(self, color, mhchem_name, species):
        """
        Initialize the ChemicalSpecies object with color, mhchem_name, and species.

        Parameters:
        color (str): The color associated with the species.
        mhchem_name (str): The name of the species in mhchem format.
        species (str): The species name.
        """
        self.color = color
        self.mhchem_name = mhchem_name
        self.species = species

    def __repr__(self):
        return f"ChemicalSpecies(color='{self.color}', mhchem_name='{self.mhchem_name}', species='{self.species}')"


def save_aqueous_species_mat(aqueous_species_names, aqueous_species_array_molal):
    """
    Function to save the aqueous species to .mat file of strings formatted for compatibility with PlanetProfile
    :param aqueous_species_names: name of aqueous species
    :param aqueous_species_array_molal: data of molal (mol/kg) of each aqueous species
    :return:
    """
    # Save to list compatible with PlanetProfile
    m_strings = [f"CustomSolution{j}" for j in range(aqueous_species_array_molal.shape[1])]
    # Transpose to iterate over conditions (columns)
    m_strings = [f"{m_strings[j]} = " + ", ".join(
        f"{species}: {aqueous_species_array_molal[i, j]}".rstrip(', ') for i, species in
        enumerate(aqueous_species_names)) for j in range(aqueous_species_array_molal.shape[1])]
    # Save to mat file
    dict = {"Data": np.array(m_strings)}
    savemat('/Users/Changster101/PlanetProfile/Europa/xRangeData.mat', dict)

if __name__ == "__main__":
    Replicate_Zolotov_H2()
    Replicate_Zolotov_Core11_Melwani_CO2()
