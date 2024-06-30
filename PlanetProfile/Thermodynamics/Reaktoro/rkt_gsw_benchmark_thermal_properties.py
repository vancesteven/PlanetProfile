import reaktoro as rkt
import numpy as np
import matplotlib.pyplot as plt
from PlanetProfile.Thermodynamics.HydroEOS import GetOceanEOS
def thermal_properties_supcrt(aqueous_species_list, speciation_ratio_mol_kg, database, T_K, P):
    """ Gets thermal properties for given aqueuous solution in provided temperature ranger
      """
    # Create lists of thermodynamic properties that will be of length NxM
    rho_kgm3 = []
    Cp_JKgK = []
    alpha_pK = []
    kTherm_WmK = []
    P_MPa = np.atleast_1d(P)
    # Go through each (P,T) combination, iterating through each P in P_MPa for a T in T_K
    # To increase efficiency, set up the chemical problem before iterating through each P
    for T in T_K:
        # Initilialize the database
        db = rkt.SupcrtDatabase(database)
        # Prescribe the solution and solids
        solution = rkt.AqueousPhase(aqueous_species_list)
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
        # Establish equilibrium constraint values
        conditions = rkt.EquilibriumConditions(specs)
        conditions.temperature(T, "K")
        # Now that we have created a Reaktoro instance for a given temperature, go through each pressure
        for P in P_MPa:
            # Establish equilibrium pressure constraint value
            conditions.pressure(P, "MPa")
            # Populate the state with the prescribed species at the given ratios
            for ion, ratio in speciation_ratio_mol_kg.items():
                state.add(ion, ratio, "mol")
            # Solve the equilibrium problem
            result = solver.solve(state, conditions)
            # Update the properties
            props.update(state)
            # Obtain the relevant aqueous only phase props
            aqueous_props = props.phaseProps("AqueousPhase")
            # Obtain the thermodynamic properties of the aqueous phase
            rho_kg_m3 = aqueous_props.density()
            Cp_J_kg_K = aqueous_props.specificHeatCapacityConstP()
            # To calculate thermal coefficient, we will multiply inverse of specific volume by its
            # partial derivative with respect to temperature
            specific_volume_m3_kg = aqueous_props.specificVolume()
            dSpecificVolumedT = aqueous_props.specificVolumeT()
            thermalExpansivity_1_K = 1 / float(specific_volume_m3_kg) * float(dSpecificVolumedT)


            # Append the values to the associated array
            rho_kgm3.append(float(rho_kg_m3))
            Cp_JKgK.append(float(Cp_J_kg_K))
            alpha_pK.append(thermalExpansivity_1_K)
            # Reset the state
            state = rkt.ChemicalState(system)
    # Turn the list into an array and transpose it into a format identical to what PP does, where pressure is in row and temperature is in columns
    rho_kgm3 = np.array(rho_kgm3).reshape(len(T_K), len(P_MPa)).T
    Cp_JKgK = np.array(Cp_JKgK).reshape(len(T_K), len(P_MPa)).T
    alpha_pK = np.array(alpha_pK).reshape(len(T_K), len(P_MPa)).T
    return rho_kgm3, Cp_JKgK, alpha_pK


def thermal_properties_phreeqc(aqueous_species_list, speciation_ratio_mol_kg, database, T_K):
    """ Gets thermal properties for given aqueuous solution in range from 0C to 25C
    """
    # Create lists of thermodynamic properties that will be of length NxM
    rho_kgm3 = []
    Cp_JKgK = []
    alpha_pK = []
    kTherm_WmK = []
    P_MPa = np.linspace(1, 1, 1)
    # Go through each (P,T) combination, iterating through each P in P_MPa for a T in T_K
    # To increase efficiency, set up the chemical problem before iterating through each P
    for T in T_K:
        # Initilialize the database
        db = rkt.PhreeqcDatabase(database)
        # Prescribe the solution
        solution = rkt.AqueousPhase(aqueous_species_list)
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
        # Establish equilibrium constraint values
        conditions = rkt.EquilibriumConditions(specs)
        conditions.temperature(T, "K")
        # Now that we have created a Reaktoro instance for a given temperature, go through each pressure
        for P in P_MPa:
            # Establish equilibrium pressure constraint value
            conditions.pressure(P, "MPa")
            # Populate the state with the prescribed species at the given ratios
            for ion, ratio in speciation_ratio_mol_kg.items():
                state.add(ion, ratio, "mol")
            # Solve the equilibrium problem
            result = solver.solve(state, conditions)
            # Update the properties
            props.update(state)
            # Obtain the relevant aqueous phase props
            aqueous_props = props.phaseProps("AqueousPhase")
            density_kg_m3 = aqueous_props.density()
            rho_J_kg_K = aqueous_props.specificHeatCapacityConstP()
            rho_kgm3.append(float(density_kg_m3))
            Cp_JKgK.append(float(rho_J_kg_K))
            H2O = db.species("H2O")
            # To obtain thermal expansivity, we will multiply recipricol of specific volume by temperature derivative of the specific volume
            specific_volume_m3_kg = float(aqueous_props.specificVolume())
            # To obtain the temperature derivative of specific volume, we will do a finite differiation calculation
            # Namely, get T a bit above and below the evaluated T
            dT = 0.01
            deltaTemperatures = [T - dT, T + dT]
            deltaSpecificVolumes = []
            # Reset the state and resolve the equilibrium at these delta temperatures
            for delta in deltaTemperatures:
                state = rkt.ChemicalState(system)
                for ion, ratio in speciation_ratio_mol_kg.items():
                    state.add(ion, ratio, "mol")
                conditions.temperature(delta)
                result = solver.solve(state, conditions)
                props.update(state)
                aqueous_props = props.phaseProps("AqueousPhase")
                deltaSpecificVolumes.append(float(aqueous_props.specificVolume()))
            # Now obtain the change in specific volume
            dSpecificVolumePlus = specific_volume_m3_kg - deltaSpecificVolumes[0]
            dSpecificVolumeLess = deltaSpecificVolumes[1] - specific_volume_m3_kg
            # Finally, find the average value between the numerically evaluated derivatives above and below
            dSpecificVolumedT = np.mean([dSpecificVolumePlus / dT, dSpecificVolumeLess / dT], axis=0)
            # Now calcuate the thermal expansivity using the finitely calculated dSpecificVolumeT
            thermalExpansivity_1_K = 1 / specific_volume_m3_kg * dSpecificVolumedT
            alpha_pK.append(thermalExpansivity_1_K)
            # Reset the state
            state = rkt.ChemicalState(system)
    # arrange the array
    rho_kgm3 = np.array(rho_kgm3).reshape(len(T_K), len(P_MPa)).T
    Cp_JKgK = np.array(Cp_JKgK).reshape(len(T_K), len(P_MPa)).T
    alpha_pK = np.array(alpha_pK).reshape(len(T_K), len(P_MPa)).T
    return rho_kgm3, Cp_JKgK, alpha_pK


def best_fit_line_plot(x1, y1, x2, y2, x3, y3, title1, title2, title3):
    '''
    Plot data in scatterplot with best fit line
    '''
    fig = plt.figure()
    # Plot #1
    slope, intercept = np.polyfit(x1, y1, 1)
    best_fit_line = slope * x1 + intercept
    ax1 = fig.add_subplot(131)
    ax1.scatter(x1, y1, color='blue', label='Data points')
    ax1.plot(x1, best_fit_line, color='red', label='Best fit line')
    ax1.axhline(0, color='green', linestyle='--')
    ax1.set_title(title1)
    # Plot #2
    slope, intercept = np.polyfit(x2, y2, 1)
    best_fit_line = slope * x2 + intercept
    ax2 = fig.add_subplot(132)
    ax2.scatter(x2, y2, color='blue', label='Data points')
    ax2.plot(x2, best_fit_line, color='red', label='Best fit line')
    ax2.axhline(0, color='green', linestyle='--')
    ax2.set_title(title2)
    #Plot #3
    slope, intercept = np.polyfit(x3, y3, 1)
    best_fit_line = slope * x3 + intercept
    ax3 = fig.add_subplot(133)
    ax3.scatter(x3, y3, color='blue', label='Data points')
    ax3.plot(x3, best_fit_line, color='red', label='Best fit line')
    ax3.axhline(0, color='green', linestyle='--')
    ax3.set_title(title3)



def core10_versus_frezchem():
    species_list = "H2O Cl- Na+ Mg+2 Ca+2 K+ SO4-2"
    speciation_ratio = {"H2O": 55.508, "Cl-": 0.5657647, "Na+": 0.4860597, "Mg+2": 0.0547421, "Ca+2": 0.0106568, "K+": 0.0105797, "SO4-2": 0.0292643}
    T_K = np.linspace(273, 300, 26)
    core10 = thermal_properties_phreeqc(species_list, speciation_ratio, "core10.dat", T_K)
    frezchem = thermal_properties_phreeqc(species_list, speciation_ratio, "frezchem.dat", T_K)
    rho_difference = (frezchem[0]-core10[0])[0]
    Cp_difference = (frezchem[1]-core10[1])[0]
    alpha_difference = (frezchem[2]-core10[2])[0]

    rho_weight_difference = rho_difference/(((frezchem[0]+core10[0])/2)[0])
    Cp_weight_difference = Cp_difference/(((frezchem[1]+core10[1])/2)[0])
    alpha_weight_difference = alpha_difference/(((frezchem[2]+core10[2])/2)[0])



    plt.show()

def supcrt_versus_core10():
    T_K = np.linspace(273, 300, 26)
    species_list = "H2O(aq) Cl- Na+ Mg+2 Ca+2 K+ SO4-2"
    speciation_ratio = {"H2O(aq)": 55.508, "Cl-": 0.5657647, "Na+": 0.4860597, "Mg+2": 0.0547421, "Ca+2": 0.0106568, "K+": 0.0105797, "SO4-2": 0.0292643}
    supcrt = thermal_properties_supcrt(species_list, speciation_ratio, "supcrt16", T_K)
    species_list = "H2O Cl- Na+ Mg+2 Ca+2 K+ SO4-2"
    speciation_ratio = {"H2O": 55.56, "Cl-": 0.5657647, "Na+": 0.4860597, "Mg+2": 0.0547421, "Ca+2": 0.0106568, "K+": 0.0105797, "SO4-2": 0.0292643}
    core10 = thermal_properties_phreeqc(species_list, speciation_ratio, "core10.dat", T_K)

    rho_difference = (supcrt[0]-core10[0])[0]
    Cp_difference = (supcrt[1]-core10[1])[0]
    alpha_difference = (supcrt[2]-core10[2])[0]

    rho_weight_difference = rho_difference/(((supcrt[0]+core10[0])/2)[0])
    Cp_weight_difference = Cp_difference/(((supcrt[1]+core10[1])/2)[0])
    alpha_weight_difference = alpha_difference/(((supcrt[2]+core10[2])/2)[0])
    title1 = "Rho weight difference % vs temperature"
    title2 = "Cp weight difference % vs temperature"
    title3 = "Alpha weight difference % vs temperature"
    best_fit_line_plot(T_K, rho_weight_difference*100, T_K, Cp_weight_difference*100, T_K, alpha_weight_difference*100, title1, title2, title3)
    plt.show()
def supcr_versus_gsw():
    # Temperature range to query over
    T_K = np.linspace(265, 268, 25)
    # Pressure value to compare (in MPa)
    P = 1
    # Species list (comparable to GSW)
    species_list = "H2O(aq) Cl- Na+ Mg+2 Ca+2 K+ SO4-2"
    speciation_ratio = {"H2O(aq)": 55.508, "Cl-": 0.5657647, "Na+": 0.4860597, "Mg+2": 0.0547421, "Ca+2": 0.0106568, "K+": 0.0105797, "SO4-2": 0.0292643}
    # Get the thermal properties from supcrt, where supcrt[0] has rho, supcrt[1] has Cp, and supcrt[2] has alpha
    supcrt = thermal_properties_supcrt(species_list, speciation_ratio, "supcrt16", T_K, P)

    # Set up the GSW implementation in PlanetProfile
    # Need to prescribe a Pressure range to obtain EOS - should be around the range of pressure that we want to query
    P_MPa = np.linspace(0.5, 1.4, 11)
    gsw_EOS = GetOceanEOS("Seawater", 35.155, P_MPa, T_K,  'Vance2018', speciation_ratio,'Millero', 'Vance2018', 'lookup')
    # Obtain the thermodyanamic properties at P for range of T_K
    rho = gsw_EOS.fn_rho_kgm3(P, T_K, grid = True)
    Cp = gsw_EOS.fn_Cp_JkgK(P, T_K, grid = True)
    alpha_1_K = gsw_EOS.fn_alpha_pK(P, T_K, grid = True)

    # Compare the thermodynamic properties betweeen Reaktoro and GSW and compute weight difference %
    # Have to use [0] at the end since it is in the form of a 2d array
    rho_difference = (supcrt[0] - rho)[0]
    rho_weight_difference = (rho_difference / (((supcrt[0] + rho) / 2)[0]))*100
    Cp_difference = (supcrt[1] - Cp)[0]
    Cp_weight_difference = (Cp_difference / (((supcrt[1] + Cp) / 2)[0])) * 100
    alpha_difference = (supcrt[2] - alpha_1_K)[0]
    alpha_weight_difference = (alpha_difference/(((supcrt[2]+alpha_1_K)/2)[0]))*100
    # Set up graph titles and plot in one figure
    title1 = "Rho weight difference % vs temperature"
    title2 = "Cp weight difference % vs temperature"
    title3 = "Alpha weight difference % vs temperature"
    best_fit_line_plot(T_K, rho_weight_difference, T_K, Cp_weight_difference, T_K,
                       alpha_weight_difference, title1, title2, title3)
    # Print the average difference between alpha of Rkt and GSW to show the difference is very small
    print(f"Average difference between alpha of RKt and GSW {np.mean(alpha_difference)}")

    rkt_coefficient = supcrt[2]/(supcrt[0]*supcrt[1])
    gsw_coefficient = alpha_1_K/(rho*Cp)
    weight_difference = (rkt_coefficient-gsw_coefficient)/(((rkt_coefficient+gsw_coefficient))/2)*100
    print(np.mean(weight_difference))

    plt.show()


if __name__ == "__main__":
    supcr_versus_gsw()

