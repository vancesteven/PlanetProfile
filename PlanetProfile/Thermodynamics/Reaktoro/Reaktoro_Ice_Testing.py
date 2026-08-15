import reaktoro as rkt
from reaktoroProps import SpeciesParser
import numpy as np
import pandas as pd
import reaktplot as plot
import math

def reaktoro_species(aqueous_species_list, speciation_ratio_mol_kg, P_MPa, T_K):
    # Initialize the database
    db = rkt.PhreeqcDatabase("frezchem.dat")
    # Prescribe the solution
    solution = rkt.AqueousPhase(aqueous_species_list)
    solution.setActivityModel(rkt.chain(rkt.ActivityModelPitzer(), rkt.ActivityModelPhreeqcIonicStrengthPressureCorrection()))
    # Obtain all related solid phases
    solids = rkt.MineralPhases()
    # Initialize the system
    system = rkt.ChemicalSystem(db, solution, solids)
    # Create dataframe and state
    df = panda_df_generator(system)
    state, conditions, specs = state_generator(system, aqueous_species_list, speciation_ratio_mol_kg)
    amounts_initial = []
    for speciesItem in system.species():
        amounts_initial.append(state.speciesAmount(speciesItem.name()))
    for reaction in system.reactions():
        print(reaction.name())
    for P in P_MPa:
        for T in T_K:
            if P > 92:
                print("Here")
            conditions.pressure(P, "MPa")
            conditions.temperature(T, "K")
            conditions.pH(7.0)
            solver = rkt.EquilibriumSolver(specs)
            restrictions = rkt.EquilibriumRestrictions(system)
            restrictions.cannotReact("H+")
            restrictions.cannotReact("OH-")
            result = solver.solve(state, conditions, restrictions)
            props = rkt.ChemicalProps(state)
            if not result.succeeded():
                print(f"{P} MPa and {T} K error")
                state = reset_state(system, speciation_ratio_mol_kg)
                continue
            amounts_final = []
            for speciesItem in system.species():
                amounts_final.append(state.speciesAmount(speciesItem.name()))
            chemical_potential = []
            for speciesItem in system.species():
                chemical_potential.append(props.speciesChemicalPotential(speciesItem.name()))
            amounts = (np.array(amounts_final) - np.array(amounts_initial))
            aprops = rkt.AqueousProps(state)
            pH = [np.array(float(aprops.pH()))]
            charge = [np.array(float(state.charge()))]
            chemPotentialDifference = [np.array(float(props.speciesChemicalPotential("H2O")-props.speciesChemicalPotential("Ice(s)")))]
            volume = [np.array(float(props.volume()))]
            internalEnergy = [np.array(float(props.internalEnergy()))]
            df.loc[len(df)] = np.concatenate(([P], [T], amounts, charge, chemical_potential, chemPotentialDifference, pH, volume, internalEnergy))
            state = reset_state(system, speciation_ratio_mol_kg)
    return df

def reset_state(system, speciation_ratio_mol_kg):
    state = rkt.ChemicalState(system)
    for ion, ratio in speciation_ratio_mol_kg.items():
        state.add(ion, ratio, "mol")
    return state





def panda_df_generator(system):
    species = []
    for speciesItem in system.species():
        species.append(speciesItem.name())
    columns = ["P"] + ["T"] + ["amount" + name for name in species] + ["Charge"] + ["ChemPotential" + name for name in species] + ["ChemPotentialDifference"] + ["pH"] + ["Volume"] + ["Internal Energy"]
    df = pd.DataFrame(columns=columns)
    return df

def state_generator(system, aqueous_species_list, speciation_ratio_mol_kg):
    # Create specs
    specs = rkt.EquilibriumSpecs(system)
    specs.pressure()
    specs.temperature()
    specs.pH()
    # Create conditions
    conditions = rkt.EquilibriumConditions(specs)
    # Create state
    state = rkt.ChemicalState(system)
    for ion, ratio in speciation_ratio_mol_kg.items():
        state.add(ion, ratio, "mol")
    return state, conditions, specs

def reaktoro_plot_species(df, query):
    if query == "P":
        pressures = df.P.unique()
        for pressure in pressures:
            pressure_df = df[df['P'] == pressure]
            fig = plot.Figure()
            fig.xaxisTitle('Temperature')
            fig.yaxisTitle(f'Amounts vs temperature at {pressure} MPa')
            column_list = list(pressure_df.columns)
            species_list = [species for species in column_list if 'amount' in species]

            for species in species_list:
                fig.drawLine(pressure_df["T"], pressure_df[species], name=species)
            fig.show()
    else:
        temperatures = df["T"].unique()
        for temperature in temperatures:
            pressure_df = df[df['T'] == temperature]
            fig = plot.Figure()
            fig.xaxisTitle('Pressure')
            fig.yaxisTitle(f'Amounts vs pressure at {temperature} K')
            column_list = list(pressure_df.columns)
            species_list = [species for species in column_list if 'amount' in species]

            for species in species_list:
                fig.drawLine(pressure_df["P"], pressure_df[species], name=species)
            fig.show()

def reaktoro_plot_chem_potential(df, query):
    if query == "P":
        pressures = df.P.unique()
        for pressure in pressures:
            pressure_df = df[df['P'] == pressure]
            fig = plot.Figure()
            fig.xaxisTitle('Temperature')
            fig.yaxisTitle(f'Chem potential vs temperature at {pressure} MPa')
            fig.drawLine(pressure_df["T"], pressure_df["ChemPotentialDifference"], name = "Delta in J/mol")
            fig.show()
    else:
        temperatures = df["T"].unique()
        for temperature in temperatures:
            pressure_df = df[df['T'] == temperature]
            fig = plot.Figure()
            fig.xaxisTitle('Pressure')
            fig.yaxisTitle(f'Chem potential vs pressure at {temperature} K')
            species_list = list(pressure_df.columns)
            fig.drawLine(pressure_df["P"], pressure_df["ChemPotentialDifference"], name = "Delta in J/mol")
            fig.show()
def reaktoro_plot_charge(df, query):
    if query == "P":
        pressures = df.P.unique()
        for pressure in pressures:
            pressure_df = df[df['P'] == pressure]
            fig = plot.Figure()
            fig.xaxisTitle('Temperature')
            fig.yaxisTitle(f'Charge vs temperature at {pressure} MPa')
            species_list = list(pressure_df.columns)
            fig.drawLine(pressure_df["T"], pressure_df['Charge'], name="Charge")
            fig.show()
    else:
        temperatures = df["T"].unique()
        for temperature in temperatures:
            pressure_df = df[df['T'] == temperature]
            fig = plot.Figure()
            fig.xaxisTitle('Temperature')
            fig.yaxisTitle(f'Charge vs pressure at {temperature} K')
            species_list = list(pressure_df.columns)
            fig.drawLine(pressure_df["P"], pressure_df['Charge'], name="Charge")
            fig.show()


if __name__ == "__main__":
    # speciation_ratio = "Cl-: 0.5657647, Na+: 0.4860597, Mg+2: 0.0547421, Ca+2: 0.0106568, K+: 0.0105797, SO4-2: 0.0292643"
    speciation_ratio = "H+: 1e-7, OH-: 1e-7"
    P_MPa = np.linspace(7, 200, 100)
    T_K = np.linspace(263, 263, 1)


    aqueous_species_string, speciation_ratio_mol_kg = SpeciesParser(speciation_ratio)
    df = reaktoro_species(aqueous_species_string, speciation_ratio_mol_kg, P_MPa, T_K)
    reaktoro_plot_species(df, "T")
    reaktoro_plot_chem_potential(df, "T")

