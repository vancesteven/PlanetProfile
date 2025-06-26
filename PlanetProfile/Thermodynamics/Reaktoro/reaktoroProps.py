import time
from PlanetProfile.Thermodynamics.Reaktoro.sigmaElectricMcCleskey2012 import elecCondMcCleskey2012
from PlanetProfile.Thermodynamics.Reaktoro.reaktoroPropsHelperFunctions import *
from PlanetProfile import _ROOT
from PlanetProfile.Utilities.defineStructs import Constants, EOSlist
from PlanetProfile.GetConfig import CustomSolutionParams
from hdf5storage import loadmat, savemat
from PlanetProfile.Utilities.DataManip import ResetNearestExtrap, ReturnConstantPTw
from collections.abc import Iterable
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import root_scalar as GetZero
from itertools import repeat
import hashlib

# Assign logger
log = logging.getLogger('PlanetProfile')
# Get FileLock, which is important to prevent race conditions in parallel computing when saving EOS
from multiprocessing import Lock
FileLock = Lock()

""" 
Check CustomSolutionConfig inputs are valid, set file paths, and save global reference to object so Rkt file can use
"""
# Ensure frezchem database is valid
CustomSolutionParams.setPaths(_ROOT)
for file_name in os.listdir(CustomSolutionParams.databasePath):
    if file_name == CustomSolutionParams.FREZCHEM_DATABASE:
        break
else:
    log.warning(
        "Input frezchem database does not match any of the available saved files.\nCheck that the input is properly spelled and has .dat at end. Using default frezchem.dat file")
    CustomSolutionParams.FREZCHEM_DATABASE = "frezchem.dat"
# Check the unit is 'g' or 'mol' (g - grams, mol - mols)
if not CustomSolutionParams.SPECIES_CONCENTRATION_UNIT == "g" and not CustomSolutionParams.SPECIES_CONCENTRATION_UNIT == "mol":
    log.warning(
        "Input species concentration unit is not valid. Check that it is either g or mol. Using mol as default")
    CustomSolutionParams.SPECIES_CONCENTRATION_UNIT = "mol"
    

def MolalConverter(ocean_species_string_g_kg):
    """
    Convert the ocean species string from g/kg of H2O to mol/kg of H2O
    """
    if CustomSolutionParams.SPECIES_CONCENTRATION_UNIT == "g":
        # Get CustomSolutionLabel
        CustomSolutionLabel = ocean_species_string_g_kg.split('=')[0]
        # Get the species string
        species_string_with_ratios =  ocean_species_string_g_kg.split('=')[1].strip()
        # Convert the species string into a string of species and its corresponding dictionary of ratios
        aqueous_species, speciation_g_kg = SpeciesFormatter(species_string_with_ratios )
        # Construct a new string that is in mol/kg
        ocean_species_string_mol_kg = f'{CustomSolutionLabel} = '
        # Go through each species in dictionary and convert ratio from g/kg to mol/kg
        for species, ratio_g_kg in speciation_g_kg.items():
            species_molar_mass_g_mol = rkt.Species(species).molarMass() * 1000
            species_mol_kg = ratio_g_kg * (1 / species_molar_mass_g_mol)
            ocean_species_string_mol_kg = ocean_species_string_mol_kg + f'{species}: {species_mol_kg}, '
        # Return adjusted string
        return ocean_species_string_mol_kg.strip(', ')


def wpptCalculator(species_string_with_ratios_mol_kg):
    """
    Calculate the total solute mass in grams of the species string and return. Assumes that the species string is in mol/kg
     Parameters
     ----------
     species_string_with_ratios_mol_kg: String of all the species that should be considered in aqueous phase and their corresponding molal ratios.
        For example, "Cl-: 19.076, Na+: 5.002, Ca2+: 0.0"
     Returns
     -------
     total_solute_mass_g: total solute mass of solution in grams
    """
    # Convert the species string into a string of species and its corresponding dictionary of ratios
    aqueous_species_string, speciation_mol_kg = SpeciesFormatter(species_string_with_ratios_mol_kg)
    # Sum up total amount of mass in grams
    total_solute_mass_g = 0
    # Go through each species in dictionary and calculate total mass in grams
    for species, ratio_mol_kg in speciation_mol_kg.items():
        # Don't add H2O mass
        if species != "H2O":
            try:
                species_molar_mass_g_mol = rkt.Species(species).molarMass() * 1000
            except:
                try:
                    db = rkt.SupcrtDatabase(CustomSolutionParams.SUPCRT_DATABASE)
                    species_molar_mass_g_mol = rkt.Species(str(db.species(species).formula())).molarMass() * 1000
                except:
                    raise ValueError(f'Species {species} not found in Reaktoro database. Check that the species is spelled correctly and is in the database.')
            species_gram_kg = species_molar_mass_g_mol * ratio_mol_kg
            total_solute_mass_g += species_gram_kg
    return total_solute_mass_g


def SpeciesParser(species_string_with_ratios_mol_kg, w_ppt):
    '''
    Converts the provided String of species and their molal ratios into formats necessary for Reaktoro. Namely, creates
    a String of all the species in the list and a dictionary with 'active' species that are added to solution (the observer species are
    automatically generated to be 1e-16 moles in the solution by Reaktoro). It also adjusts speciation to match the w_ppt of the solution, if specified.

     Parameters
     ----------
     species_string_with_ratios_mol_kg: String of all the species that should be considered in aqueous phase and their corresponding molal ratios.
        For example, "Cl-: 19.076, Na+: 5.002, Ca2+: 0.0"
     Returns
     -------
     aqueous_species_string: String that has all species names that should be considered in aqueous phase
     speciation_ratio_mol_per_kg: Dictionary of active species and the values of their molal ratio (mol/kg of water) that matches input w_ppt
     EOS_lookup_label: EOS table lookup label
    '''
    # Get CustomSolutionLabel
    CustomSolutionLabel = species_string_with_ratios_mol_kg.split('=')[0].strip()
    # Get the species string
    species_string_with_ratios = species_string_with_ratios_mol_kg.split('=')[1].strip()
    # Convert the species string into a string of species and its corresponding dictionary of ratios
    aqueous_species_string, speciation_mol_kg = SpeciesFormatter(species_string_with_ratios)
    # Sum up total amount of mass in grams
    total_solute_mass_g = wpptCalculator(species_string_with_ratios)
    # If w_ppt is specified >= 0, then adjust ocean_species_string so it corresponds to input w_ppt
    # However, if w_ppt = total_solute_mass_g, this is the case where we want to use Planet's input speciation so we don't adjust speciation
    if w_ppt is not None and w_ppt > 0 and w_ppt != total_solute_mass_g:
        log.info(
            f'Changing ocean comp to match input w_ppt of {w_ppt}.')
        # Find the common ratio to get to w_ppt
        common_ratio = w_ppt / total_solute_mass_g
        # Multiply every value in dictionary by common_ratio
        final_speciation_mol_kg = {key: value * common_ratio for key, value in speciation_mol_kg.items()
                                             if key != "H2O"}
    # In this case, we have pure water so let's simplify strings and speciation to nothing, which will be filled with H+, OH-, H2O below
    elif w_ppt == 0:
        aqueous_species_string = ''
        final_speciation_mol_kg = {}
    # Otherwise, we assume the input species is the desired input and we return the calculated w_ppt
    else:
        final_speciation_mol_kg = speciation_mol_kg
        log.info(
            f'Calculated w_ppt is {total_solute_mass_g}.')
    # Now, let's ensure we have H+, OH- and H2O in speciation
    # Add H+ and OH- to species string so we can track pH
    if not "H+" in final_speciation_mol_kg:
        aqueous_species_string = aqueous_species_string + " H+"
        final_speciation_mol_kg['H+'] = 0
    if not "OH-" in final_speciation_mol_kg:
        aqueous_species_string = aqueous_species_string + " OH-"
        final_speciation_mol_kg['OH-'] = 0
    # Check if water is in the aqueous species string and dictionary and if not, add it, ensuring to update the weight to be 1kg
    if not "H2O" in final_speciation_mol_kg:
        aqueous_species_string = aqueous_species_string + " H2O"
        final_speciation_mol_kg["H2O"] =  float(1 / rkt.waterMolarMass)
    # Construct an EOS lookup label
    EOS_lookup_label = "_".join(f"{key}-{value}" for key, value in final_speciation_mol_kg.items())
    # If we are considering solid phases, then get all relevant solid phases
    if CustomSolutionParams.SOLID_PHASES:
        if CustomSolutionParams.SOLID_PHASES_TO_CONSIDER == 'All':
            solid_phases = 'All'
        else:
            solid_phases = []
            for solid_phase in CustomSolutionParams.SOLID_PHASES_TO_CONSIDER:
                if solid_phase in Constants.SolidPhases:
                    solid_phases = solid_phases + Constants.SolidPhases[solid_phase]
                else:
                    solid_phases.append(solid_phase)
        # Get only the solid phases that are relevant to the system - reduces runtime by not considering solids that would not appear in system
        db = rkt.SupcrtDatabase(CustomSolutionParams.SUPCRT_DATABASE)
        supcrt_aqueous_species_string, supcrt_speciation_ratio_mol_kg = species_convertor_compatible_with_supcrt(db,
            aqueous_species_string, final_speciation_mol_kg, Constants.PhreeqcToSupcrtNames)
        solid_phases_to_consider = RelevantSolidSpecies(db, supcrt_aqueous_species_string, solid_phases)
        # Append solids to EOS lookup label
        EOS_lookup_label = f'{EOS_lookup_label}_{"_".join(solid_phases_to_consider.split())}'
    else:
        solid_phases_to_consider = None
    # Return aqueous species, solid species, and EOS lookup label
    return aqueous_species_string, final_speciation_mol_kg, solid_phases_to_consider, EOS_lookup_label


def SpeciesFormatter(species_string_with_ratios):
    '''
    Converts provided String of species and ratios into formats necessary for Reaktoro. Namely, creates
    a String of all the species in the list and a dictionary with 'active' species that are added to solution (the observer species are
    automatically generated to be 1e-16 moles in the solution by Reaktoro).

     Parameters
     ----------
     species_string_with_ratios: String of all the species that should be considered in aqueous phase and their corresponding ratios.
        For example, "Cl-: 19.076, Na+: 5.002, Ca2+: 0.0"
     Returns
     -------
     aqueous_species_string: String that has all species names that should be considered in aqueous phase
     speciation_ratio_mol_per_kg: Dictionary of active species and the values of their molar ratio (mol/kg of water)
    '''
    # Dictionary of speciation
    speciation_ratio_mol_per_kg = {}
    # Go through each species and corresponding ratio_mol_per_kg and add to corresponding lists
    for species_with_ratio in species_string_with_ratios.split(", "):
        species, ratio_mol_per_kg = species_with_ratio.split(": ")
        speciation_ratio_mol_per_kg[species] = float(ratio_mol_per_kg)
    # String of speciation
    aqueous_species_str = " ".join(speciation_ratio_mol_per_kg.keys())
    # Return string and species
    return aqueous_species_str, speciation_ratio_mol_per_kg

def checkSpeciesCompatibleWithFrezchem(aqueous_species_string, speciation_ratio_mol_kg, frezchemPath):
    """
    Check that the input species are available in Frezchem. If not, then remove the species fromthe string and dictionary. Return final string compatible with Frezchem.
    Only called when CustomSolutionParams.REMOVE_SPECIES_NA_IN_FREZCHEM is enabled.
    """
    if CustomSolutionParams.REMOVE_SPECIES_NA_IN_FREZCHEM:
        frezchem_aqueous_species_string = ''
        frezchem_speciation_ratio_mol_per_kg = {}
        # Initialize the database
        db = rkt.PhreeqcDatabase.fromFile(frezchemPath)
        for species_name in speciation_ratio_mol_kg:
            try:
                db.species(species_name)
                # If error not thrown, then species exists in database so we add to final species string
                frezchem_aqueous_species_string = frezchem_aqueous_species_string + f'{species_name} '
                frezchem_speciation_ratio_mol_per_kg[species_name] = speciation_ratio_mol_kg[species_name]
            except:
                # If error thrown, then species does not exist in database so we don't add to final species string
                log.debug(f'Removing {species_name} from consideration in frezchem database')
                pass
        return frezchem_aqueous_species_string, frezchem_speciation_ratio_mol_per_kg
    else:
        return aqueous_species_string, speciation_ratio_mol_kg



class EOSLookupTableLoader():
    def __init__(self, aqueous_species_string, speciation_ratio_mol_per_kg, ocean_solid_species, EOS_lookup_label):
        self.name = EOS_lookup_label
        self.fHashName = hashlib.md5(EOS_lookup_label.encode()).hexdigest()
        self.EOSname = f'{EOS_lookup_label}' # Name used as unique identifier to save to EOS list - See SaveEOSToDisk for use of EOSname
        self.props_fLookup = os.path.join(CustomSolutionParams.rktPath, 'CustomSolutionLookupTables',
                                          f'props_{self.fHashName}.pkl.gz')
        self.phase_fLookup = os.path.join(CustomSolutionParams.rktPath, 'CustomSolutionLookupTables',
                                          f'phase_{self.fHashName}.pkl.gz')
        if os.path.exists(self.props_fLookup) and os.path.exists(self.phase_fLookup):
            log.debug('EOS table saved to disk for ocean composition with:\n'
                      f'The following aqueous species: {aqueous_species_string}.\n'
                      f'The following solid species: {ocean_solid_species}')
            self.alreadySavedToDisk = True
            self.props_EOS = load_dict_from_pkl(self.props_fLookup)
            self.phase_EOS = load_dict_from_pkl(self.phase_fLookup)
            TRkt_K = self.props_EOS['T_K']
            PRkt_MPa = self.props_EOS['P_MPa']
            self.Pmin = np.min(PRkt_MPa)
            self.Pmax = np.max(PRkt_MPa)
            self.EOSdeltaP = np.maximum(np.round(np.mean(np.diff(PRkt_MPa)), 2), 0.001)
            self.Tmin = np.min(TRkt_K)
            self.Tmax = np.max(TRkt_K)
            self.EOSdeltaT = np.maximum(np.round(np.mean(np.diff(TRkt_K)), 2), 0.001)
        else:
            log.warning(
                f'EOSlookup table with label {EOS_lookup_label} does not exist, so we will generate a property EOS.\n'
                f'Namely, we will query Reaktoro with input pressures to find associated properties and save to disk.\n'
                f'Considering the following aqueous species: {aqueous_species_string}.\n'
                f'Considering the following solid species: {ocean_solid_species}')
            self.alreadySavedToDisk = False
            # Obtain properties from Reaktoro
            start_time = time.time()
            # Define lower and upper pressure limits and create pressure range to query over
            Pinitialmin = Constants.RktPmin_MPa
            Pinitialmax = Constants.FrezchemPmax_MPa
            deltaP = CustomSolutionParams.EOS_deltaP
            nPs = round((Pinitialmax - Pinitialmin) / deltaP)
            PRkt_MPa = np.linspace(Pinitialmin, Pinitialmax, nPs)

            # Define Frezchem and Supcrt systems
            # When defining the Frezchem system, let's check if we want to remove self-consistency between aqueous species and those of frezchem system (see REMOVE_SPECIES_NA_IN_FREZCHEM for details)
            if CustomSolutionParams.REMOVE_SPECIES_NA_IN_FREZCHEM:
                log.warning(f'CustomSolutionParams.REMOVE_SPECIES_NA_IN_FREZCHEM is enabled, meaning we will check to see if any species in input \n'
                            f'ocean string are not compatible with frezchem and remove it from consideration in iceIh-liquid equilibria consideration.\n'
                            f'This will remove self-consistency between ocean thermodynamics and Ih phase equilibria')
                frezchem_aqueous_species_string, frezchem_speciation_ratio_mol_per_kg = checkSpeciesCompatibleWithFrezchem(aqueous_species_string, speciation_ratio_mol_per_kg, CustomSolutionParams.frezchemPath)
                Frezchem_System = PhreeqcGeneratorForChemicalConstraint(frezchem_aqueous_species_string, frezchem_speciation_ratio_mol_per_kg,
                                                                        "mol", CustomSolutionParams.frezchemPath, CustomSolutionParams.maxIterations)
            else:
                Frezchem_System = PhreeqcGeneratorForChemicalConstraint(aqueous_species_string,
                                                                        speciation_ratio_mol_per_kg, "mol",
                                                                        CustomSolutionParams.frezchemPath, CustomSolutionParams.maxIterations)

            Supcrt_System = SupcrtGenerator(aqueous_species_string, speciation_ratio_mol_per_kg, "mol",
                                            CustomSolutionParams.SUPCRT_DATABASE, ocean_solid_species,
                                            Constants.PhreeqcToSupcrtNames, CustomSolutionParams.maxIterations)
            # Get freezing temperatures for pressure range calculated by Frezchem
            TFreezing_K = self.RktFreezingTemperatureFinder(Frezchem_System, PRkt_MPa)
            # Get Seafreeze correction
            SeafreezePureWaterCorrectorFunctions = SeafreezePureWaterCorrector()
            # Correct data
            TFreezing_K = TFreezing_K - SeafreezePureWaterCorrectorFunctions.fn_freezing_temp_correction(PRkt_MPa)

            # Load properties into dictionary and save to a .mat file
            RktPhase_Dictionary = {'P_MPa': PRkt_MPa, 'TFreezing_K': TFreezing_K}

            self.Tmin, self.Tmax, self.Pmin, Pmax_supcrt = self.PropsConstraintFinder(Supcrt_System)
            self.EOSdeltaP = CustomSolutionParams.EOS_deltaP
            self.EOSdeltaT = CustomSolutionParams.EOS_deltaT
            self.Pmax = Constants.EOSPmax_MPa
            log.warning(f'Input EOS will span from [Pmin, Pmax] = {self.Pmin}, {Pmax_supcrt} MPa by a step of {self.EOSdeltaP} and [Tmin, Tmax] = {self.Tmin}, {self.Tmax} K by a step of {self.EOSdeltaT}.\n'
                        f'This data will be linearly extrapolated up to {self.Pmax} MPa and then have a pure water correction applied to its data')
            nTs = round((self.Tmax - self.Tmin) / self.EOSdeltaT)
            nPs = round((Pmax_supcrt - self.Pmin) / self.EOSdeltaP)
            nPs_extra = round(
                (self.Pmax - (Pmax_supcrt + Constants.EOSdeltaP_For_Extrapolation)) / Constants.EOSdeltaP_For_Extrapolation)
            TRkt_K = np.linspace(self.Tmin, self.Tmax, nTs)
            PRkt_MPa = np.linspace(self.Pmin, Pmax_supcrt, nPs)
            PExtrap_MPa = np.linspace(Pmax_supcrt + Constants.EOSdeltaP_For_Extrapolation, self.Pmax, nPs_extra)
            P_MPa = np.concatenate((PRkt_MPa, PExtrap_MPa))

            rho_kgm3, Cp_JKgK, alpha_pK, VP_kms, mu_J_mol = self.RktProps(Supcrt_System, PRkt_MPa, TRkt_K)
            # Calculate bulk modulus from sound speed
            KS_GPa = rho_kgm3 * VP_kms ** 2 * 1e-3

            # Load data into linear functions
            ufn_rho_kgm3 = RegularGridInterpolator((PRkt_MPa, TRkt_K), rho_kgm3, method='linear', bounds_error=False,
                                                   fill_value=None)
            ufn_Cp_JkgK = RegularGridInterpolator((PRkt_MPa, TRkt_K), Cp_JKgK, method='linear', bounds_error=False,
                                                  fill_value=None)
            ufn_alpha_pK = RegularGridInterpolator((PRkt_MPa, TRkt_K), alpha_pK, method='linear', bounds_error=False,
                                                   fill_value=None)
            ufn_VP_kms = RegularGridInterpolator((PRkt_MPa, TRkt_K), VP_kms, method='linear', bounds_error=False,
                                                 fill_value=None)
            ufn_KS_GPa = RegularGridInterpolator((PRkt_MPa, TRkt_K), KS_GPa, method='linear', bounds_error=False,
                                                 fill_value=None)
            ufn_mu_J_mol = RegularGridInterpolator((PRkt_MPa, TRkt_K), mu_J_mol, method='linear', bounds_error=False,
                                                   fill_value=None)

            # Extrapolate up to Constants.EOSPmax_MPa
            P_MPa_extrap_mesh, T_K_extrap_mesh = np.meshgrid(PExtrap_MPa, TRkt_K, indexing='ij')
            extrapEvalPts = np.column_stack((P_MPa_extrap_mesh.ravel(), T_K_extrap_mesh.ravel()))
            rho_kgm3 = np.concatenate((rho_kgm3, ufn_rho_kgm3(extrapEvalPts).reshape(P_MPa_extrap_mesh.shape)))
            Cp_JKgK = np.concatenate((Cp_JKgK, ufn_Cp_JkgK(extrapEvalPts).reshape(P_MPa_extrap_mesh.shape)))
            alpha_pK = np.concatenate((alpha_pK, ufn_alpha_pK(extrapEvalPts).reshape(P_MPa_extrap_mesh.shape)))
            VP_kms = np.concatenate((VP_kms, ufn_VP_kms(extrapEvalPts).reshape(P_MPa_extrap_mesh.shape)))
            KS_GPa = np.concatenate((KS_GPa, ufn_KS_GPa(extrapEvalPts).reshape(P_MPa_extrap_mesh.shape)))
            mu_J_mol = np.concatenate((mu_J_mol, ufn_mu_J_mol(extrapEvalPts).reshape(P_MPa_extrap_mesh.shape)))

            # Get Seafreeze correction
            P_MPa_mesh, T_K_mesh = np.meshgrid(P_MPa, TRkt_K, indexing='ij')
            evalPts = np.column_stack((P_MPa_mesh.ravel(), T_K_mesh.ravel()))
            SeafreezePureWaterCorrectorFunctions = SeafreezePureWaterCorrector()

            # Correct data
            rho_kgm3 = (rho_kgm3 - SeafreezePureWaterCorrectorFunctions.fn_rho_correction(evalPts).reshape(
                P_MPa_mesh.shape)).astype(np.float32)

            Cp_JKgK = (Cp_JKgK - SeafreezePureWaterCorrectorFunctions.fn_Cp_correction(evalPts).reshape(
                P_MPa_mesh.shape)).astype(np.float32)

            alpha_pK = (alpha_pK - SeafreezePureWaterCorrectorFunctions.fn_alpha_correction(evalPts).reshape(
                P_MPa_mesh.shape)).astype(np.float32)

            VP_kms = (VP_kms - SeafreezePureWaterCorrectorFunctions.fn_VP_correction(evalPts).reshape(
                P_MPa_mesh.shape)).astype(np.float32)

            KS_GPa = (KS_GPa - SeafreezePureWaterCorrectorFunctions.fn_KS_correction(evalPts).reshape(
                P_MPa_mesh.shape)).astype(np.float32)

            mu_J_mol = (mu_J_mol - SeafreezePureWaterCorrectorFunctions.fn_mu_correction(evalPts).reshape(
                P_MPa_mesh.shape)).astype(np.float32)

            reaktoro_calculation_times = time.time() - start_time
            log.warning(f'EOS grid generation took {reaktoro_calculation_times} seconds.')

            # Load properties into a dictionary and save to a .pkl file
            RktProps_Dictionary = {'T_K': TRkt_K, 'P_MPa': P_MPa, 'rho': rho_kgm3, 'Cp': Cp_JKgK, 'alpha': alpha_pK,
                                   'VP': VP_kms, 'KS': KS_GPa, 'mu': mu_J_mol}

            self.props_EOS = RktProps_Dictionary
            self.phase_EOS = RktPhase_Dictionary
            """
            We use the multiprocessing FileLock here to prevent race conditions where multiple processes try to write to save file at same time.
            This is mainly applicable to Parallel computing in Explorogram or Inductograms. Might not be the most efficient solution, so should be revisited.
            """
            # Acquire the lock before writing
            with FileLock:
                save_dict_to_pkl(self.props_EOS, self.props_fLookup)
                save_dict_to_pkl(self.phase_EOS, self.phase_fLookup)


    def RktFreezingTemperatureFinder(self, Frezchem_System, P_MPa, TMin_K=220, TMax_K=300, significant_threshold=0.1):
        """
         Calculates the temperature at which the prescribed aqueous solution freezes. Utilizes the reaktoro framework to
         constrain the equilibrium position at the prescribed pressure and the chemical potential difference between ice and liquid water at 0.1,
          therefore calculating and returning the temperature (within the range) at which ice begins to form.

         Parameters
         ----------
         speciation_ratio_mol_per_kg: the ratio of species in the aqueous solution in mol/kg of water. Should be a dictionary
         with the species as the key and its ratio as its value.
         P_MPa: the desired equilibrium freezing pressure(s).
         TMin_K: the lower limit of temperature that Reaktoro should query over
         TMax_K: the upper limit of temperature that Reaktoro should query over
         significant_threshold: the amount of moles of ice present for H2O to be considered in solid phase. Default is 1e-14 moles.

         Returns
         -------
         t_freezing_K: the temperature at which the solution begins to freeze.
         P_MPa_adjusted: adjusted pressure range that removes values that did not converge
         """
        # Disable chemical convergence warnings that Reaktoro raises. We handle these internally instead and throw more specific warnings when they appear.
        rkt.Warnings.disable(906)
        # Create freezing temperatures list and indices of pressures to remove, if necessary
        freezing_temperatures = []
        # Create frezchem system
        db, system, initial_state, conditions, solver, props = Frezchem_System
        # Set conditions
        conditions.set("IP", significant_threshold)
        conditions.setLowerBoundTemperature(TMin_K, "K")
        conditions.setUpperBoundTemperature(TMax_K, "K")
        state = initial_state.clone()
        for index, P in enumerate(P_MPa):
            P = float(P)
            conditions.pressure(P, "MPa")
            # Solve the equilibrium problem with warm start
            result = solver.solve(state, conditions)
            if not result.succeeded():
                state = initial_state.clone()
                result = solver.solve(state, conditions)
            # Update the properties
            props.update(state)
            # Obtain the equilibrium temperature
            equilibrium_temperature = props.temperature()
            # Check if the result succeeded
            if result.succeeded():
                freezing_temperatures.append(float(equilibrium_temperature))
            # If the result failed, then do not include this in list of freezing temperatures and remove from pressure list to avoid its use in creating spline
            else:
                log.debug(
                    f'Failed to find freezing temperature at pressure {P} MPa when generating phase data from Frezchem. Will extrapolate this value.')
                freezing_temperatures.append(np.nan)
                state = initial_state.clone()
        freezing_temperatures = np.array(freezing_temperatures)
        # Make interpolator for values that could converge
        P_MPa_calculated = P_MPa[~np.isnan(freezing_temperatures)]
        freezing_temperatures_calculated = freezing_temperatures[~np.isnan(freezing_temperatures)]
        fn_frezchem_phaseRGI = RegularGridInterpolator((P_MPa_calculated,), freezing_temperatures_calculated,
                                                       method='linear', bounds_error=False, fill_value=None)
        # Interpolate the missing values
        P_MPa_extrapolate = P_MPa[np.isnan(freezing_temperatures)]
        freezing_temperatures_extrapolate = fn_frezchem_phaseRGI(P_MPa_extrapolate)
        freezing_temperatures[np.isnan(freezing_temperatures)] = freezing_temperatures_extrapolate
        # Return freezing temperatures
        return freezing_temperatures

    def PropsConstraintFinder(self, Supcrt_System):
        """ Finds the constraints of pressure and temperature that are compatible with both phreeqc and supcrt databases for the given solution speciation.
            ~ The lower temperature limit will be dynamically set by supcrt as equilibrium for this database can only be found starting at temperature of ~240K.
            ~ The upper temperature is statically set to 400K
            ~ The lower pressure limit will by statistically set to 0.1 MPa
            ~ The upper pressure limit will be set to 200 MPa, the maximum pressure before high pressure ices

            Returns adjusted pressure and temperature ranges with the found constraints. This approach is similar to MgSO4 solution, however we take advantage
            of RKt's dynamic ability to find chemical equilibrium different depending on the solution, rather than hard coding in constraints.
        Args:
            Supcrt_System: Supcrt system tuple with following items: db, system, state, conditions, solver, props

        Returns:
            P_MPa (AT MOST shape N): Adjusted array that only has values from initial P_MPa that are compatible with RKt.
            P (AT MOST shape M): Adjusted array that only has values from initial T_K that are compatible with RKt.
        """
        # Dynamically find lower Tmin_K that supcrt can converge at
        temperatureChange = lambda T: 0.5 - temperature_constraint(T, Supcrt_System)
        if temperatureChange(Constants.SupcrtTmin_K) < 0:
            log.debug(f"Will set lower EOS temperature boundary to {Constants.SupcrtTmin_K} K.")
            Tmin_K = Constants.SupcrtTmin_K
        else:
            try:
                Tmin_K = GetZero(temperatureChange, bracket=[Constants.SupcrtTmin_K, Constants.T0]).root
            except:
                log.warning(
                    f"Supcrt could not converge at equilibrium for minimum EOS temperatures of {Constants.SupcrtTmin_K} to {273} K."
                    f"We will still attempt to generate EOS for aqueous composition, but note that extrapolation will be occurring."
                    f"Try to simplify aqueous composition so convergence is possible at these lower temperatures.")
                Tmin_K = Constants.T0

        Tmax_K = Constants.SupcrtTmax_K
        Pmin_MPa = Constants.RktPmin_MPa

        # Dynamically find upper Tmax_K that supcrt can converge at
        pressureChange = lambda P: 0.5 - pressure_constraint(P, Supcrt_System)
        if pressureChange(Constants.SupcrtPmax_MPa) < 0:
            log.debug(f"Will set upper EOS pressure boundary to {Constants.SupcrtPmax_MPa} MPa.")
            Pmax_MPa = Constants.SupcrtPmax_MPa
        else:
            try:
                Pmax_MPa = GetZero(pressureChange, bracket=[Constants.PminHPices_MPa, Constants.SupcrtPmax_MPa]).root
            except:
                log.warning(
                    f"Supcrt could not converge at equilibrium for maxiumum pressures of {Constants.PminHPices_MPa} to {Constants.SupcrtPmax_MPa} MPa."
                    f"We will still attempt to generate EOS for aqueous composition, but note that extrapolation will be occurring."
                    f"Try to simplify aqueous composition so convergence is possible at these higher pressures.")
                Pmax_MPa = Constants.PminHPices_MPa
        return Tmin_K, Tmax_K, Pmin_MPa, Pmax_MPa

    def RktProps(self, Supcrt_System, P_MPa, T_K):
        """ Determine density rho, heat capacity Cp, thermal expansivity alpha,
            and thermal conductivity kTherm as functions of pressure P and
            temperature T for the provided solution species list and corresponding molarity ratios.
            Implements Reaktoro Supcrt database to find these thermal properties at equilibrium for the prescribed pressure and temperature.
            Importantly, if any thermal properties cannot be found due to equilibrium divergence, we perform a linear interpolation on these values, however
            this outcome is HIGHLY UNLIKELY since we have found constraints that are compatible with Rkt, and user is warned if this does occur.
        Args:
            Supcrt_System: Supcrt system tuple with following items: db, system, state, conditions, solver, props
            P_MPa (float, shape N): Pressures in MPa
            T_K (float, shape M): Temperature in K
        Returns:
            rho_kgm3 (float, shape NxM): Mass density of liquid in kg/m^3
            Cp_JkgK (float, shape NxM): Isobaric heat capacity of liquid in J/(kg K)
            alpha_pK (float, shape NxM): Thermal expansivity of liquid in 1/K
            kTherm_WmK (float, shape NxM): Thermal conductivity of liquid in W/(m K)
        """
        # Create lists of thermodynamic properties that will be of length NxM
        rho_kgm3 = []
        Cp_JKgK = []
        alpha_pK = []
        Vp_kms = []
        mu_J_mol = []
        # Disable chemical convergence warnings that Reaktoro raises. We handle these internally instead and throw more specific warnings when they appear.
        rkt.Warnings.disable(906)
        # Create Reaktoro objects
        db, system, initial_state, conditions, solver, props = Supcrt_System
        state = initial_state.clone()
        # Go through each pressure, temperature combination
        for P in P_MPa:
            # Establish equilibrium pressure constraint value
            conditions.pressure(P, "MPa")
            for T in T_K:
                conditions.temperature(T, "K")
                # Solve the equilibrium problem
                result = solver.solve(state, conditions)
                # Check if the equilibrium problem failed from this warm start approach
                if not result.succeeded():
                    # Cold restart the state and try again
                    state = initial_state.clone()
                    result = solver.solve(state, conditions)
                if result.succeeded():
                    # Update the properties
                    props.update(state)
                    # If it did succeed, obtain the relevant aqueous only phase props
                    aqueous_props = props.phaseProps("AqueousPhase")
                    # Obtain the thermodynamic properties of the aqueous phase
                    rho_kg_m3 = aqueous_props.density()
                    Cp_J_kg_K = aqueous_props.specificHeatCapacityConstP()

                    # Obtain chemical potential of H2O(aq) and correct with correction spline
                    aqueous_mu_J_mol = props.speciesChemicalPotential("H2O(aq)")
                    aqueous_mu_J_mol = float(aqueous_mu_J_mol)

                    # To calculate thermal coefficient, we will multiply inverse of specific volume by its
                    # partial derivative with respect to temperature
                    specific_volume_m3_kg = aqueous_props.specificVolume()
                    dSpecificVolumedT = aqueous_props.specificVolumeT()
                    thermalExpansivity_1_K = 1 / float(specific_volume_m3_kg) * float(dSpecificVolumedT)

                    # To calculate sound speed, obtain specific volume and pressure derivative of specific volume
                    specific_volume_m3_kg = float(aqueous_props.specificVolume())
                    pressure_derivative_specific_volume = float(aqueous_props.specificVolumeP())
                    # Calculate commpressibility factor, multiplying negative inverse of specific volume by pressure derivative of specific volume
                    k_compressibility = -(1 / specific_volume_m3_kg) * pressure_derivative_specific_volume
                    # Calculate sound speed, multiplying k by density and taking the product to the -0.5
                    c_m_s = (rho_kg_m3 * k_compressibility) ** (-0.5)
                    # Convert to km/s
                    c_km_s = c_m_s / 1000

                    # Append the values to the associated array
                    rho_kgm3.append(float(rho_kg_m3))
                    Cp_JKgK.append(float(Cp_J_kg_K))
                    alpha_pK.append(float(thermalExpansivity_1_K))
                    Vp_kms.append(float(c_km_s))
                    mu_J_mol.append(float(aqueous_mu_J_mol))
                # Otherwise, the equilibrium problem failed so we need to handle it accordingly
                # THIS IS VERY UNLIKELY SINCE WE HAVE ESTABLISHED CONSTRAINTS OF TEMPERATURE AND PRESSURE THAT SHOULD WORK WITH SUPCRT
                else:
                    # Log to the user that the computation was unsuccessful and that missed value will be extrapolated
                    log.debug(
                        f"Unsuccessful computation at: {props.pressure() / 1e+6} MPa and {props.temperature()} K.\n"
                        f"The temperature and pressure may be out of bounds, thus we will linearly extrapolate this value from previous ones using a 1d interopolation")
                    # Reset the state
                    state = initial_state.clone()
                    # Append zeros, which will be modified later
                    rho_kgm3.append(np.nan)
                    Cp_JKgK.append(np.nan)
                    alpha_pK.append(np.nan)
                    Vp_kms.append(np.nan)
                    mu_J_mol.append(np.nan)
            # After going through each entire temperature range for given pressure, reset state
            state = initial_state.clone()
        # Turn the list into an array of shape NxM
        rho_kgm3 = np.array(rho_kgm3).reshape((P_MPa.size, -1))
        Cp_JKgK = np.array(Cp_JKgK).reshape((P_MPa.size, -1))
        alpha_pK = np.array(alpha_pK).reshape((P_MPa.size, -1))
        Vp_kms = np.array(Vp_kms).reshape((P_MPa.size, -1))
        mu_J_mol = np.array(mu_J_mol).reshape((P_MPa.size, -1))

        # In the case that any properties could not be calculated, we must linearly interpolate these
        # THIS IS HIGHLY UNLIKELY SINCE WE HAVE FOUND CONSTRAINTS COMPATIBLE WITH RKT, BUT JUST IN CASE THIS IS IMPLEMENTED (has not been rigorously tested)
        if np.sum(np.isnan(rho_kgm3) | np.isinf(rho_kgm3)) > 0:
            rho_kgm3, Cp_JKgK, alpha_pK, Vp_kms, mu_J_mol = interpolation_2d(P_MPa,
                                                                             [rho_kgm3, Cp_JKgK, alpha_pK, Vp_kms,
                                                                              mu_J_mol])
        # Return the thermodynamic properties
        return rho_kgm3, Cp_JKgK, alpha_pK, Vp_kms, mu_J_mol

def RktProps(EOSLookupTable, P_MPa, T_K, EXTRAP):
    fn_RktProps = RktPropsLookup(EOSLookupTable)
    EOS_deltaP = fn_RktProps.EOSdeltaP
    EOS_deltaT = fn_RktProps.EOSdeltaT
    if not EXTRAP:
        newP_MPa, newT_K = ResetNearestExtrap(P_MPa, T_K, fn_RktProps.Pmin, fn_RktProps.Pmax,
                                              fn_RktProps.Tmin, fn_RktProps.Tmax)
        if (not np.all(newP_MPa == P_MPa)) or (not np.all(newT_K == T_K)):
            log.warning('Extrapolation is disabled for ocean fluids, and input EOS P and/or T ' +
                        'extend beyond the EOS properties lookup table limits of [Pmin, Pmax] = ' +
                        f'{fn_RktProps.Pmin}, {fn_RktProps.Pmax} MPa and [Tmin, Tmax] = ' +
                        f'{fn_RktProps.Tmin}, {fn_RktProps.Tmax} K.')
            P_MPa = np.unique(newP_MPa)
            T_K = np.unique(newT_K)
        # Ensure that P_MPa and T_K have at least 5 values in order to make rectbivariatespline
        if np.size(T_K) < 5:
            T_K = np.linspace(T_K[0], T_K[-1] + EOS_deltaT, 5)
            EOS_deltaT = np.maximum(np.round(np.mean(np.diff(T_K)), 2), 0.001)
        if np.size(P_MPa) < 5:
            P_MPa = np.linspace(P_MPa[0], P_MPa[-1] + EOS_deltaP, 5)
            EOS_deltaP = np.maximum(np.round(np.mean(np.diff(P_MPa)), 2), 0.001)

    evalPts = fn_RktProps.fn_evalPts(P_MPa, T_K)
    nPs = np.size(P_MPa)
    # Interpolate the input data to get the values corresponding to the current ocean comp,
    # then get the property values for the input (P,T) pairs and reshape to how they need
    # to be formatted for use in the ocean EOS.
    rho_kgm3 = np.reshape(fn_RktProps.fn_rho_kgm3(evalPts), (nPs, -1))
    Cp_JkgK = np.reshape(fn_RktProps.fn_Cp_JkgK(evalPts), (nPs, -1))
    alpha_pK = np.reshape(fn_RktProps.fn_alpha_pK(evalPts), (nPs, -1))
    kTherm_WmK = fn_RktProps.fn_kTherm_WmK(P_MPa, T_K, 0, grid =True)

    return P_MPa, T_K, rho_kgm3, Cp_JkgK, alpha_pK, kTherm_WmK, EOS_deltaP, EOS_deltaT


class RktPropsLookup:
    def __init__(self, EOSLookupTable):
        self.fLookup = f'{EOSLookupTable.name}_Props'
        if self.fLookup in EOSlist.loaded.keys():
            log.debug(f'EOS properties lookup table with label {self.fLookup} already loaded.')
            self.fn_rho_kgm3, self.fn_Cp_JkgK, self.fn_alpha_pK, self.fn_kTherm_WmK, self.fn_VP_kms, self.fn_KS_GPa, self.fn_mu_J_mol, self.fn_evalPts = \
            EOSlist.loaded[self.fLookup]
            self.Pmin, self.Pmax, self.EOSdeltaP, self.Tmin, self.Tmax, self.EOSdeltaT = EOSlist.ranges[self.fLookup]
        else:
            fRktProps = EOSLookupTable.props_EOS
            TRkt_K = fRktProps['T_K']
            PRkt_MPa = fRktProps['P_MPa']
            self.fn_rho_kgm3 = RegularGridInterpolator((PRkt_MPa, TRkt_K), fRktProps['rho'], method='linear',
                                                       bounds_error=False, fill_value=None)
            self.fn_Cp_JkgK = RegularGridInterpolator((PRkt_MPa, TRkt_K), fRktProps['Cp'], method='linear',
                                                      bounds_error=False, fill_value=None)
            self.fn_alpha_pK = RegularGridInterpolator((PRkt_MPa, TRkt_K), fRktProps['alpha'], method='linear',
                                                       bounds_error=False, fill_value=None)
            self.fn_VP_kms = RegularGridInterpolator((PRkt_MPa, TRkt_K), fRktProps['VP'], method='linear',
                                                     bounds_error=False, fill_value=None)
            self.fn_KS_GPa = RegularGridInterpolator((PRkt_MPa, TRkt_K), fRktProps['KS'], method='linear',
                                                     bounds_error=False, fill_value=None)
            self.fn_mu_J_mol = RegularGridInterpolator((PRkt_MPa, TRkt_K), fRktProps['mu'], method='linear',
                                                       bounds_error=False, fill_value=None)
            self.fn_kTherm_WmK = ReturnConstantPTw(const=Constants.kThermWater_WmK)

            self.Pmin = EOSLookupTable.Pmin
            self.Pmax = EOSLookupTable.Pmax
            self.EOSdeltaP = EOSLookupTable.EOSdeltaP
            self.Tmin = EOSLookupTable.Tmin
            self.Tmax = EOSLookupTable.Tmax
            self.EOSdeltaT = EOSLookupTable.EOSdeltaT

            self.fn_kTherm_WmK = ReturnConstantPTw(const=Constants.kThermWater_WmK)
            # Save functions to EOSlist so they can be referenced in future
            EOSlist.loaded[self.fLookup] = (
                self.fn_rho_kgm3, self.fn_Cp_JkgK, self.fn_alpha_pK, self.fn_kTherm_WmK, self.fn_VP_kms, self.fn_KS_GPa,
                self.fn_mu_J_mol, self.fn_evalPts)
            EOSlist.ranges[self.fLookup] = (self.Pmin, self.Pmax, self.EOSdeltaP, self.Tmin, self.Tmax, self.EOSdeltaT)

    def fn_evalPts(self, Pin_MPa, Tin_K):
        P_MPa = ensureArray(Pin_MPa)
        T_K = ensureArray(Tin_K)
        P_Mesh, T_Mesh = np.meshgrid(P_MPa, T_K, indexing='ij')
        out = np.column_stack((P_Mesh.ravel(), T_Mesh.ravel()))
        return np.array(out)


class SeafreezePureWaterCorrector:
    def __init__(self):
        # Obtain internal aqueous chemical potential correction spline (generating one if necessary)
        spline_path = os.path.join(CustomSolutionParams.rktPath, 'Reaktoro_Saved_Files',
                                   'seafreeze_pure_water_correction.mat')
        if os.path.exists(spline_path):
            correction_data = loadmat(spline_path)
            supcrt_P_MPa = correction_data['Supcrt_P_MPa']
            supcrt_T_K = correction_data['Supcrt_T_K']
            delta_rho = correction_data['density_difference']
            delta_Cp = correction_data['isobaric_heat_capacity_difference']
            delta_alpha = correction_data['thermal_expansivity_difference']
            delta_VP = correction_data['sound_speed_difference']
            delta_KS = correction_data['bulk_modulus_difference']
            delta_mu = correction_data['chemical_potential_difference']

            frezchem_P_MPa = correction_data['frezchem_P_MPa']
            delta_freezing_temp = correction_data['freezing_temperature_difference']
        # This file will be uploaded to Github, so user should not need to generate this data
        else:
            log.warning('The Seafreeze pure water correction data is not available on disk. We will generate this data and save to disk so it can be called in the future.')
            frezchem_P_MPa, delta_freezing_temp = self.FrezchemCorrectionData()
            supcrt_P_MPa, supcrt_T_K, delta_rho, delta_Cp, delta_alpha, delta_VP, delta_KS, delta_mu = self.SupcrtSeafreezeCorrectionData()
            SupcrtSeafreezeCorrectionData = {'Supcrt_P_MPa': supcrt_P_MPa, 'Supcrt_T_K': supcrt_T_K,
                                                                  'density_difference': delta_rho,
                                             'isobaric_heat_capacity_difference': delta_Cp,
                                             'thermal_expansivity_difference': delta_alpha,
                                             'sound_speed_difference': delta_VP,
                                             'bulk_modulus_difference': delta_KS,
                                             'chemical_potential_difference': delta_mu,
                                             'frezchem_P_MPa': frezchem_P_MPa,
                                             'freezing_temperature_difference': delta_freezing_temp}
            savemat(spline_path, SupcrtSeafreezeCorrectionData)
        self.fn_rho_correction = RegularGridInterpolator((supcrt_P_MPa, supcrt_T_K), delta_rho, method = 'linear', bounds_error = False, fill_value = None)
        self.fn_Cp_correction = RegularGridInterpolator((supcrt_P_MPa, supcrt_T_K), delta_Cp, method='linear', bounds_error=False,
                                                        fill_value=None)
        self.fn_alpha_correction = RegularGridInterpolator((supcrt_P_MPa, supcrt_T_K), delta_alpha, method='linear',
                                                           bounds_error=False, fill_value=None)
        self.fn_VP_correction = RegularGridInterpolator((supcrt_P_MPa, supcrt_T_K), delta_VP, method='linear', bounds_error=False,
                                                        fill_value=None)
        self.fn_KS_correction = RegularGridInterpolator((supcrt_P_MPa, supcrt_T_K), delta_KS, method='linear', bounds_error=False,
                                                        fill_value=None)
        self.fn_mu_correction = RegularGridInterpolator((supcrt_P_MPa, supcrt_T_K), delta_mu, method='linear', bounds_error=False,
                                                        fill_value=None)
        self.fn_freezing_temp_correction = RegularGridInterpolator((frezchem_P_MPa,), delta_freezing_temp,
                                                           method='linear', bounds_error=False, fill_value=None)
    def FrezchemCorrectionData(self):
        def whichphaseChooser(P, T):
            PT = np.empty((1,), dtype='object')
            PT[0] = (P, T)
            return sfz.whichphase(PT)[0]

        eos_P_MPa = np.linspace(0.1, 200, 200)
        sfz_T_freezing = []
        for P in eos_P_MPa:
            phaseChange = lambda T: 0.5 - (1 - int(whichphaseChooser(P, T) > 0))
            Tfreeze_K = GetZero(phaseChange, bracket=[240, 280]).root
            sfz_T_freezing.append(Tfreeze_K)
        sfz_T_freezing = np.array(sfz_T_freezing)

        # Obtain frezchem freezing temperatures
        rkt_T_freezing = []
        aqueous_species_list = 'H+ OH- H2O'
        speciation_ratio_mol_kg = {'H2O': float(1 / rkt.waterMolarMass)}
        frezchem_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Databases', 'frezchem.dat')
        frezchem = PhreeqcGeneratorForChemicalConstraint(aqueous_species_list, speciation_ratio_mol_kg, "mol",
                                                         frezchem_file_path, CustomSolutionParams.maxIterations)
        db, system, initial_state, conditions, solver, props = frezchem
        state = initial_state.clone()
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
            if not result.succeeded():
                print("HELLO")
            # Update the properties
            props.update(state)
            # Obtain the equilibrium temperature
            rkt_T_freezing.append(float(props.temperature()))
            # Reset the state
            state = initial_state.clone()
        rkt_T_freezing = np.array(rkt_T_freezing)

        # Find difference in freezing temperatures
        difference_in_T_freezings = rkt_T_freezing - sfz_T_freezing
        return eos_P_MPa, difference_in_T_freezings
    def SupcrtSeafreezeCorrectionData(self):
        """ Obtain the aqueous H2O thermodynamic differnce of pure water for Supcrt and Seafreeze across a 2-D grid of pressure and temperatures."""
        P_MPa = np.linspace(Constants.RktPmin_MPa, Constants.EOSPmax_MPa, 500)
        T_K = np.linspace(Constants.SupcrtTmin_K, Constants.SupcrtTmax_K, 200)
        P_MPa_supcrt = np.linspace(Constants.RktPmin_MPa, Constants.SupcrtPmax_MPa, 500)

        PT = np.array([P_MPa, T_K], dtype = object)
        out = sfz.getProp(PT, 'water1')
        sfz_rho = out.rho
        sfz_Cp = out.Cp
        sfz_alpha = out.alpha
        sfz_VP = out.vel * 1e-3
        sfz_KS = out.Ks * 1e-3
        sfz_mu = out.G * rkt.waterMolarMass # Multiply Gibbs free energy by molar mass of H2O to get chemical potential of pure water

        rkt_rho = []
        rkt_Cp = []
        rkt_alpha = []
        rkt_VP = []
        rkt_KS = []
        rkt_mu = []
        aqueous_species_list = 'H+ OH- H2O(aq)'
        speciation_ratio_mol_kg = {'H2O(aq)': float(1 / rkt.waterMolarMass)}
        supcrt = SupcrtGenerator(aqueous_species_list, speciation_ratio_mol_kg, "mol", "supcrt16", None, Constants.PhreeqcToSupcrtNames, CustomSolutionParams.maxIterations)
        db, system, initial_state, conditions, solver, props = supcrt
        state = initial_state.clone()
        P_MPa_mesh, T_K_mesh = np.meshgrid(P_MPa_supcrt, T_K, indexing='ij')
        # Create a nditer iterator
        it = np.nditer([P_MPa_mesh, T_K_mesh], flags=['multi_index'])
        # Go through each P, T combination
        for P, T in it:
            P = float(P)
            T = float(T)
            conditions.temperature(T, "K")
            # Establish equilibrium pressure constraint value
            conditions.pressure(P, "MPa")
            # Solve the equilibrium problem
            result = solver.solve(state, conditions)
            if result.succeeded():
                # Update the properties
                props.update(state)
                # Get the aqueous phase props
                aqueousProps = props.phaseProps("AqueousPhase")
                # Check if the equilibrium problem succeeded
                rho_kgm3 = float(aqueousProps.density())
                Cp_JKgK = float(aqueousProps.specificHeatCapacityConstP())
                mu_J_mol = float(props.speciesChemicalPotential("H2O(aq)"))

                # To calculate thermal coefficient, we will multiply inverse of specific volume by its
                # partial derivative with respect to temperature
                specific_volume_m3_kg = aqueousProps.specificVolume()
                dSpecificVolumedT = aqueousProps.specificVolumeT()
                thermalExpansivity_1_K = 1 / float(specific_volume_m3_kg) * float(dSpecificVolumedT)

                # To calculate sound speed, obtain specific volume and pressure derivative of specific volume
                specific_volume_m3_kg = float(aqueousProps.specificVolume())
                pressure_derivative_specific_volume = float(aqueousProps.specificVolumeP())
                # Calculate commpressibility factor, multiplying negative inverse of specific volume by pressure derivative of specific volume
                k_compressibility = -(1 / specific_volume_m3_kg) * pressure_derivative_specific_volume
                # Calculate sound speed, multiplying k by density and taking the product to the -0.5
                c_m_s = (rho_kgm3 * k_compressibility) ** (-0.5)
                # Convert to km/s
                c_km_s = c_m_s / 1000
                # Calculate KS
                KS_GPa = rho_kgm3 * float(c_km_s) ** 2 * 1e-3
            else:
                rho_kgm3 = np.nan
                Cp_JKgK = np.nan
                mu_J_mol = np.nan
                thermalExpansivity_1_K = np.nan
                c_km_s = np.nan
                KS_GPa = np.nan


            rkt_rho.append(rho_kgm3)
            rkt_Cp.append(Cp_JKgK)
            rkt_mu.append(mu_J_mol)
            rkt_alpha.append(thermalExpansivity_1_K)
            rkt_VP.append(float(c_km_s))
            rkt_KS.append(float(KS_GPa))
            state = initial_state.clone()
        # Turn the list into a 2d array
        rkt_rho = np.array(rkt_rho).reshape((P_MPa_supcrt.size, -1))
        rkt_Cp = np.array(rkt_Cp).reshape((P_MPa_supcrt.size, -1))
        rkt_alpha = np.array(rkt_alpha).reshape((P_MPa_supcrt.size, -1))
        rkt_VP = np.array(rkt_VP).reshape((P_MPa_supcrt.size, -1))
        rkt_KS = np.array(rkt_KS).reshape((P_MPa_supcrt.size, -1))
        rkt_mu = np.array(rkt_mu).reshape((P_MPa_supcrt.size, -1))

        # Create RegularGridInterpolators and cubic extrapolate up to 2000MPa
        fn_rkt_rho = RegularGridInterpolator((P_MPa_supcrt, T_K), rkt_rho, method='linear', bounds_error=False,
                                         fill_value=None)
        fn_rkt_Cp = RegularGridInterpolator((P_MPa_supcrt, T_K), rkt_Cp, method='linear', bounds_error=False,
                                         fill_value=None)
        fn_rkt_alpha = RegularGridInterpolator((P_MPa_supcrt, T_K), rkt_alpha, method='linear', bounds_error=False,
                                            fill_value=None)
        fn_rkt_VP = RegularGridInterpolator((P_MPa_supcrt, T_K), rkt_VP, method='linear', bounds_error=False,
                                         fill_value=None)
        fn_rkt_KS = RegularGridInterpolator((P_MPa_supcrt, T_K), rkt_KS, method='linear', bounds_error=False,
                                         fill_value=None)
        fn_rkt_mu = RegularGridInterpolator((P_MPa_supcrt, T_K), rkt_mu, method='linear', bounds_error=False,
                                         fill_value=None)

        # Get data points up to 2000MPa
        P_MPa_mesh, T_K_mesh = np.meshgrid(P_MPa, T_K, indexing='ij')
        evalPts = np.column_stack((P_MPa_mesh.ravel(), T_K_mesh.ravel()))
        rkt_rho = fn_rkt_rho(evalPts).reshape(P_MPa_mesh.shape)
        rkt_Cp = fn_rkt_Cp(evalPts).reshape(P_MPa_mesh.shape)
        rkt_alpha = fn_rkt_alpha(evalPts).reshape(P_MPa_mesh.shape)
        rkt_VP = fn_rkt_VP(evalPts).reshape(P_MPa_mesh.shape)
        rkt_KS = fn_rkt_KS(evalPts).reshape(P_MPa_mesh.shape)
        rkt_mu = fn_rkt_mu(evalPts).reshape(P_MPa_mesh.shape)

        # Find difference between rkt and sfz
        delta_rho = rkt_rho - sfz_rho
        delta_Cp = rkt_Cp - sfz_Cp
        delta_alpha = rkt_alpha - sfz_alpha
        delta_VP = rkt_VP - sfz_VP
        delta_KS = rkt_KS - sfz_KS
        delta_mu = rkt_mu - sfz_mu

        # Return
        return P_MPa, T_K, delta_rho, delta_Cp, delta_alpha, delta_VP, delta_KS, delta_mu




class RktSeismic:
    def __init__(self, EOSLookupTable, EXTRAP):
        self.EXTRAP = EXTRAP
        self.WARNED = False  # Track whether user has been warned about extrapolation
        # Gets the RktPropsLookup again. This should be quick as we have already loaded it into EOSlist using RktProps called before
        fn_RktProps = RktPropsLookup(EOSLookupTable)
        # Get pressure and temperature limits
        self.Pmin = fn_RktProps.Pmin
        self.Pmax = fn_RktProps.Pmax
        self.Tmin = fn_RktProps.Tmin
        self.Tmax = fn_RktProps.Tmax
        # Reassign the functions so they can be referenced when object is called
        self.fn_VP_kms = fn_RktProps.fn_VP_kms
        self.fn_KS_GPa = fn_RktProps.fn_KS_GPa

    def __call__(self, P_MPa, T_K, grid=False):
        if not self.EXTRAP:
            newP_MPa, newT_K = ResetNearestExtrap(P_MPa, T_K, self.Pmin, self.Pmax, self.Tmin, self.Tmax)
            if (not np.all(newP_MPa == P_MPa)) or (not np.all(newT_K == T_K)):
                if not self.WARNED:
                    log.warning(
                        'Extrapolation is disabled for ocean fluids, and input EOS P and/or T for seismic calculations ' +
                        f'extend beyond the EOS properties lookup table limits of\n[Pmin, Pmax] = ' +
                        f'{self.Pmin}, {self.Pmax} MPa and [Tmin, Tmax] = ' +
                        f'{self.Tmin}, {self.Tmax} K. Will reset the inputs to stay within the ranges.')
                    self.WARNED = True
                P_MPa = newP_MPa
                T_K = newT_K
        if grid:
            evalPts = tuple(np.meshgrid(P_MPa, T_K, indexing='ij'))
        else:
            evalPts = np.column_stack((P_MPa, T_K))
        VP_kms = np.squeeze(self.fn_VP_kms(evalPts))
        KS_GPa = np.squeeze(self.fn_KS_GPa(evalPts))
        return VP_kms, KS_GPa

    def fn_evalPts(self, Pin_MPa, Tin_K):
        P_MPa = ensureArray(Pin_MPa)
        T_K = ensureArray(Tin_K)
        out = [[P, T] for P in P_MPa for T in T_K]
        return np.array(out)

class RktPhaseLookup:
    def __init__(self, EOSLookupTable, P_MPa, T_K, EOS_deltaP, EOS_deltaT):
        self.fLookup = f'{EOSLookupTable.name}_Phase'
        if self.fLookup in EOSlist.loaded.keys():
            log.debug(f'EOS phase lookup table with label {self.fLookup} already loaded.')
            self.fn_frezchem_phaseRGI, self.fn_mu_J_mol = EOSlist.loaded[self.fLookup]
            self.Pmin, self.Pmax, self.deltaP, self.deltaT = EOSlist.ranges[self.fLookup]
        else:
            fRktPhase = EOSLookupTable.phase_EOS
            PRkt_MPa = fRktPhase['P_MPa']
            self.Pmin = np.min(PRkt_MPa)
            self.deltaP = np.maximum(np.round(np.mean(np.diff(PRkt_MPa)), 2), 0.001)
            self.fn_frezchem_phaseRGI = RegularGridInterpolator((PRkt_MPa,), fRktPhase['TFreezing_K'],
                                                       method='linear', bounds_error=False, fill_value=None)

            # Gets the RktPropsLookup again. This should be quick as we have already loaded it into EOSlist using RktProps called before
            fn_RktProps = RktPropsLookup(EOSLookupTable)

            # Get the temperature limits
            self.Pmax = fn_RktProps.Pmax
            self.Tmin = fn_RktProps.Tmin
            self.Tmax = fn_RktProps.Tmax
            self.deltaT = fn_RktProps.EOSdeltaT
            # Reassign the functions so they can be referenced when object is called
            self.fn_mu_J_mol = fn_RktProps.fn_mu_J_mol

            EOSlist.loaded[self.fLookup] = self.fn_frezchem_phaseRGI, self.fn_mu_J_mol
            EOSlist.ranges[self.fLookup] = (self.Pmin, self.Pmax, self.deltaP, self.deltaT)
        # Generate pressure and temperature arrays that extend to limits of EOS pressure and temperature but
        # use the fidelity of the temperature and pressure steps of EOS (i.e. deltaT and deltaP)
        self.deltaP = np.min([self.deltaP, EOS_deltaP])
        self.deltaT = np.min([self.deltaT, EOS_deltaT])
        self.P_MPa_to_query = np.arange(P_MPa[0], P_MPa[-1], self.deltaP)
        self.T_K_to_query = np.arange(T_K[0], T_K[-1], self.deltaT)
        self.phase_lookup_grid = self.phase_lookup_grid_generator(self.P_MPa_to_query, self.T_K_to_query, self.fn_frezchem_phaseRGI, self.fn_mu_J_mol)
        self.fn_phase = RegularGridInterpolator((self.P_MPa_to_query, self.T_K_to_query), self.phase_lookup_grid, method='nearest', bounds_error=False, fill_value=None)

    def __call__(self, P_MPa, T_K,  grid=False):
        if grid:
            P_MPa, T_K = np.meshgrid(P_MPa, T_K, indexing='ij')
        return (self.fn_phase((P_MPa, T_K)))

    def phase_lookup_grid_generator(self, P_MPa, T_K, freezing_temperature_function_below_200_MPa, mu_function_above_200_MPa):
        P_MPa_below_200_MPa_index = np.searchsorted(P_MPa, 200, side = 'left')
        phases = []
        if P_MPa_below_200_MPa_index > 0:
            P_MPa_below_200_MPa = P_MPa[0:P_MPa_below_200_MPa_index]
            freezing_temperatures = freezing_temperature_function_below_200_MPa(P_MPa_below_200_MPa)
            freezing_temperatures, T_K_pts = np.meshgrid(freezing_temperatures, T_K, indexing='ij')
            phases = phases + (T_K_pts < freezing_temperatures).astype(np.int_).tolist()

        if P_MPa_below_200_MPa_index < P_MPa.size:
            P_MPa_above_200_MPa = P_MPa[P_MPa_below_200_MPa_index:]
            evalPts_RGI = tuple(np.meshgrid(P_MPa_above_200_MPa, T_K))
            # Make sure we're always passing a properly formatted array to seafreeze
            # Create a 2-element array with pressure and temperature arrays
            evalPts_sfz = np.array([P_MPa_above_200_MPa, T_K], dtype=object)
            ptsh = (P_MPa_above_200_MPa.size, T_K.size)
            max_phase_num = max([p for p in Constants.seafreeze_ice_phases.keys()])
            comp = np.full(ptsh + (max_phase_num + 1,), np.nan)
            
            for phase, name in Constants.seafreeze_ice_phases.items():
                if phase == 0:
                    mu_J_mol = mu_function_above_200_MPa(evalPts_RGI).T
                else:
                    # Ensure we handle single value arrays properly
                    if P_MPa_above_200_MPa.size == 1 or T_K.size == 1:
                        # Create a proper grid for seafreeze to work with
                        sfz_P, sfz_T = np.meshgrid(P_MPa_above_200_MPa, T_K, indexing='ij')
                        sfz_PT = np.array([sfz_P.flatten(), sfz_T.flatten()], dtype=object)
                        mu_J_mol = sfz.getProp(sfz_PT, name).G * Constants.m_gmol['H2O'] / 1000
                        mu_J_mol = mu_J_mol.reshape(sfz_P.shape)
                    else:
                        try:
                            mu_J_mol = sfz.getProp(evalPts_sfz, name).G * Constants.m_gmol['H2O'] / 1000
                        except:
                            from seafreeze.seafreeze import defpath, _get_tdvs, _is_scatter
                            from seafreeze.seafreeze import phases as seafreeze_phases
                            from mlbspline import load
                            phasedesc = seafreeze_phases[name]
                            sp = load.loadSpline(defpath, phasedesc.sp_name)
                            # Calc density and isentropic bulk modulus
                            isscatter = _is_scatter(evalPts_sfz)
                            tdvs = _get_tdvs(sp, evalPts_sfz, isscatter)
                            mu_J_mol = tdvs.G * Constants.m_gmol['H2O'] / 1000
                            #raise ValueError(f"Error in seafreeze calculation for phase {name}. Check the input values {mu_J_mol}.")
                sl = tuple(repeat(slice(None), 2)) + (phase,)
                comp[sl] = np.squeeze(mu_J_mol)
            all_nan_sl = np.all(np.isnan(comp), -1)  # Find slices where all values are nan along the innermost axis
            out_phase = np.full(ptsh, np.nan)
            out_phase[~all_nan_sl] = np.nanargmin(comp[~all_nan_sl], -1)
            phases = phases + out_phase.tolist()
        phases = np.array(phases).reshape(P_MPa.size, T_K.size)
        return phases



def ensureArray(var):
    if isinstance(var, Iterable):
        return var
    else:
        return np.array([var])


class RktConduct():
    def __init__(self, aqueous_species_list, speciation_ratio_mol_kg, ocean_solid_species, fn_species):
        """
        Initialize the RKtConduct() object and parse the aqueous species list into a format compatible with elecCondMcClevskey2012()
        """
        # Convert H2O label to H2O(aq) label for compatability with Supcrt database
        db = rkt.SupcrtDatabase(CustomSolutionParams.SUPCRT_DATABASE)
        self.aqueous_species_list, self.speciation_ratio_mol_per_kg = species_convertor_compatible_with_supcrt(db, aqueous_species_list, speciation_ratio_mol_kg, Constants.PhreeqcToSupcrtNames)
        self.ocean_solid_phases = ocean_solid_species
        self.fn_species = fn_species
        # Get reference to dictionary that holds already calculated speciations in form of key (P_MPa, T_K)
        self.calculated_speciations = fn_species.calculated_speciations

    def __call__(self, P_MPa, T_K, grid=False):
        """
        Finds electrical conductivity for given ions. Currently does not adjust for pressure.

        Args:
            P_MPa (array): Array of pressures
            T_K (float or array): Array of temperatures
            grid (bool): Whether or not to convert to coordinate grid (optional)

        Returns:
            VP_kms (float, Shape N): Corresponding sound speeds in km/s
            KS_GPa (float, Shape N): Corresponding bulk modulus in GPa
        """
        # Ensure P_MPa and T_K are numpy arrays
        P_MPa, T_K = np.array(P_MPa), np.array(T_K)

        # Store the original P_MPa and T_K for the cache key
        original_P_MPa, original_T_K = P_MPa, T_K

        # If grid is needed or if inputs are not 1D arrays, create a grid
        if grid or (P_MPa.ndim != 1 or T_K.ndim != 1):
            P_MPa, T_K = np.meshgrid(P_MPa, T_K, indexing='ij')
        elif P_MPa.size == 0 or T_K.size == 0:
            # Return empty array if input is empty
            return np.array([])

        # Check if arrays are mismatched but might be intended for grid
        if P_MPa.size != T_K.size and not (P_MPa.size == 1 or T_K.size == 1):
            P_MPa, T_K = np.meshgrid(P_MPa, T_K, indexing='ij')
            grid = True  # Indicate that grid should be used

        # Check if speciation data has already been calculated for this grid
        key = (tuple(P_MPa.ravel()), tuple(T_K.ravel()))  # Use original arrays for the key
        if key not in self.calculated_speciations:
            # Calculate speciation if not found
            self.calculated_speciations[key] = self.fn_species(original_P_MPa, original_T_K, grid=grid)
        pH, speciation, species_names = self.calculated_speciations[key]

        # Convert species names to compatible format
        McClevsky_speciation = self.McClevskyIonParser(speciation, species_names, self.speciation_ratio_mol_per_kg)
        # McCleskey function requires Celcius units
        T_C = T_K - Constants.T0
        return elecCondMcCleskey2012(P_MPa, T_C, McClevsky_speciation)

    def McClevskyIonParser(self, speciation_array, species_names_array, speciation_ratio_mol_kg):
        """
        Parse through provided species list and convert to format compatible with McClevsky
        Args:
            speciation_array: array of species names
            species_array: array of species names
        Returns:
            Ion list compatible with McClevsky
        """
        ions = {}
        for index, species in np.ndenumerate(species_names_array):
            # Check that the species is an ion
            if "+" in species or "-" in species:
                # Rewrite the species into a format that McCLevsky can handle
                # First, convert speciation ratio into mols to be compatible with McClevskyIonParser
                mol_amount = speciation_array[index]
                # If any values are NaN, then we replace with constant speciation
                nan_indices = np.isnan(mol_amount)
                # Interpolate NaN values using the other values
                if np.any(nan_indices):
                    mol_amount[nan_indices] = np.interp(np.flatnonzero(nan_indices), np.flatnonzero(~nan_indices), mol_amount[~nan_indices])
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


class RktHydroSpecies():
    def __init__(self, aqueous_species_list, speciation_ratio_mol_kg, ocean_solid_phases):
        self.aqueous_species_list = aqueous_species_list
        self.speciation_ratio_mol_kg = speciation_ratio_mol_kg
        self.ocean_solid_phases = ocean_solid_phases
        # Create a dictionary of calculated speciations that will hold calculated speciations for input P_MPa and T_K
        self.calculated_speciations = {}

    def __call__(self, P_MPa, T_K, grid=False):
        """ Calculates speciation of composition at provided pressure and temperature using Supcrt. Notably,
        we have to reset the pressure since we cannot calculate equilibrium above 500MPa using Supcrt."""
        # Reset P_MPa so that it does not extend above 500MPa, since supcrt cannot go above this pressure
        newP_MPa, newT_K = ResetNearestExtrap(P_MPa, T_K, P_MPa[0], Constants.SupcrtPmax_MPa, T_K[0], T_K[-1])
        if (not np.all(newP_MPa == P_MPa)) or (not np.all(newT_K == T_K)):
            log.warning(
                'Supcrt can only accurately calculate hydrosphere species up to 500MPa, so we will reset hydro species function call up to the ocean depth that correlates with 500MPa.' + f'{newP_MPa[0]}, {newP_MPa[-1]} MPa and [Tmin, Tmax] = ' + f'{newT_K[0]}, {newT_K[-1]} K. Will reset the inputs to stay within the ranges.')
            P_MPa = newP_MPa
            T_K = newT_K

        if grid:
            P_MPa, T_K = np.meshgrid(P_MPa, T_K, indexing='ij')
        else:
            P_MPa = np.array(P_MPa)
            T_K = np.array(T_K)
            if (np.size(P_MPa) == 0 or np.size(T_K) == 0):
                # If input is empty, return empty array
                return np.array([])
            elif ((np.size(P_MPa) != np.size(T_K)) and not (np.size(P_MPa) == 1 or np.size(T_K) == 1)):
                # If arrays are different lengths, they are probably meant to get a 2D output
                P_MPa, T_K = np.meshgrid(P_MPa, T_K, indexing='ij')
                grid = True

        # If we are doing grid, we must first unravel meshgrid and then reformat output
        if grid:
            P_MPa_flat = P_MPa.ravel()
            T_K_flat = T_K.ravel()
            pH, species, species_names = self.species_at_equilibrium(P_MPa_flat, T_K_flat)
            # Reshape species to (num_species, P_MPa.size, T_K.size)
            num_species = species_names.size
            species = species.reshape((num_species, P_MPa.shape[0], P_MPa.shape[1]))
            pH = pH.reshape(P_MPa.shape)
        else:
            pH, species, species_names = self.species_at_equilibrium(P_MPa, T_K)
        # Let's save the speciation in the dictionary (which we will reference in RktConduct to reduce runtime)
        self.calculated_speciations[(tuple(P_MPa.ravel()), tuple(T_K.ravel()))] = pH, species, species_names
        return pH, species, species_names

    def species_at_equilibrium(self, P_MPa, T_K):
        """
        Go through P_MPa and T_K  and calculate equilibrium speciation of aqueous and solid species, as well as pH.
        Return species above
        """
        # Keep track of time it takes to do calculation
        start_time = time.time()
        # Establish supcrt generator
        db, system, initial_state, conditions, solver, props = SupcrtGenerator(self.aqueous_species_list, self.speciation_ratio_mol_kg,
                                      "mol", CustomSolutionParams.SUPCRT_DATABASE, self.ocean_solid_phases, Constants.PhreeqcToSupcrtNames, CustomSolutionParams.maxIterations)
        state = initial_state.clone()
        # Prepare lists for pH and species amounts
        pH_list = []
        species_list = [[] for _ in range(len(system.species()))]
        species_names = np.array([species.name() for species in system.species()])  # Extract species names
        for P, T in zip(P_MPa, T_K):
            conditions.pressure(P, "MPa")
            conditions.temperature(T, "K")
            # Solve the equilibrium problem using the hot-start approach
            result = solver.solve(state, conditions)
            if not result.succeeded():
                # Attempt a cold start
                state = initial_state.clone()
                result = solver.solve(state, conditions)
            if result.succeeded():
                # Update props and extract data
                props.update(state)
                aprops = rkt.AqueousProps(props)
                pH_list.append(float(aprops.pH()))
                for k, species in enumerate(system.species()):
                    species_list[k].append(float(state.speciesAmount(species.name())))
            else:
                # If we fail to find equilibrium, let's just append the pH from last successful attempt
                log.warning(f"Failed to find equilibrium at {P} MPa and {T} K. Filling with NaN.")
                pH_list.append(np.nan)
                for k in range(len(system.species())):
                    species_list[k].append(np.nan)
                # Reset after each temperature
                state = initial_state.clone()
        # Convert lists to arrays
        pH_array = np.array(pH_list)
        species_array = np.array(species_list)

        # Log time it took to calculate speciation
        end_time = time.time()
        log.debug(f'{end_time-start_time} seconds to calculate hydrosphere species')

        # Return the filtered results
        return pH_array, species_array, species_names

class RktRxnAffinity():
    def __init__(self, aqueous_species_list, speciation_ratio_mol_per_kg, ocean_solid_species):
        # Convert H2O label to H2O(aq) label for compatability with Supcrt database
        db = rkt.SupcrtDatabase(CustomSolutionParams.SUPCRT_DATABASE)
        self. aqueous_species_list, self.speciation_ratio_mol_per_kg = species_convertor_compatible_with_supcrt(db, aqueous_species_list, speciation_ratio_mol_per_kg, Constants.PhreeqcToSupcrtNames)
        self.ocean_solid_species = ocean_solid_species

    def __call__(self, P_MPa, T_K, reaction, concentrations, grid=False):
        """ Calculates affinity of reaction, whose species are at prescribed concentrations at disequilibrium, at provided pressure and temperature
        using Supcrt. Notably,we have to reset the pressure since we cannot calculate equilibrium above 500MPa using Supcrt."""
        # Reset P_MPa so that it does not extend above 500MPa, since supcrt cannot go above this pressure
        newP_MPa, newT_K = ResetNearestExtrap(P_MPa, T_K, P_MPa[0], Constants.SupcrtPmax_MPa, T_K[0], T_K[-1])
        if (not np.all(newP_MPa == P_MPa)) or (not np.all(newT_K == T_K)):
            log.warning(
                'Supcrt can only accurate calculate hydrosphere species up to 500MPa, so we will reset reaction affinity '
                'function call up to the ocean depth that correlates with 500MPa.' +
                f'{newP_MPa[0]}, {newP_MPa[-1]} MPa and [Tmin, Tmax] = ' +
                f'{newT_K[0]}, {newT_K[-1]} K. Will reset the inputs to stay within the ranges.')
            P_MPa = newP_MPa
            T_K = newT_K
        if grid:
            P_MPa, T_K = np.meshgrid(P_MPa, T_K, indexing='ij')
        else:
            P_MPa = np.array(P_MPa)
            T_K = np.array(T_K)
            if (np.size(P_MPa) == 0 or np.size(T_K) == 0):
                # If input is empty, return empty array
                return np.array([])
            elif ((np.size(P_MPa) != np.size(T_K)) and not (np.size(P_MPa) == 1 or np.size(T_K) == 1)):
                # If arrays are different lengths, they are probably meant to get a 2D output
                P_MPa, T_K = np.meshgrid(P_MPa, T_K, indexing='ij')
        parsed_reaction = self.reaction_parser(reaction)



        return self.reaction_affinity(parsed_reaction, concentrations, P_MPa, T_K)

    def reaction_affinity(self, reaction, rxn_disequilibrium_concentrations, P_MPa, T_K):
        # First, parse through the rxn_disequilibrium_concentrations and convert to a dictionary of species and their concentrations
        for species, concentration in rxn_disequilibrium_concentrations.items():
            if concentration is None:
                concentration = self.speciation_ratio_mol_per_kg[species]
                rxn_disequilibrium_concentrations[species] = concentration
            if type(concentration) == float:
                rxn_disequilibrium_concentrations[species] = concentration
            elif type(concentration) == dict:
                species_reference = concentration['reference species']
                lambda_function = concentration['equation']
                concentration = lambda_function(self.speciation_ratio_mol_per_kg[species_reference])
                rxn_disequilibrium_concentrations[species] = float(concentration)
        # Update speciation_ratio_mol_per_kg dictionary with disequilibrium concentrations
        speciation_ratio_mol_per_kg = {**self.speciation_ratio_mol_per_kg, **rxn_disequilibrium_concentrations}

        aqueous_species_list = " ".join(speciation_ratio_mol_per_kg.keys())

        # Keep track of time it takes to do calculation
        start_time = time.time()
        # Establish supcrt generator
        db, system, initial_state, conditions, solver, props = SupcrtGenerator(aqueous_species_list,
                                                                       speciation_ratio_mol_per_kg,
                                                                       "mol",
                                                                       CustomSolutionParams.SUPCRT_DATABASE, self.ocean_solid_species, Constants.PhreeqcToSupcrtNames, CustomSolutionParams.maxIterations)
        affinity_kJ = []
        # Create a copy of the state
        state = initial_state.clone()
        # Go through each P and T
        for P, T in zip(P_MPa, T_K):
            # Set conditions
            conditions.temperature(T, "K")
            conditions.pressure(P, "MPa")
            state.setPressure(P, "MPa")
            state.setTemperature(T, "K")
            # Let's always update props with the initial state that has disequilibrium speciation
            props.update(initial_state)
            # Calculate Q disequilibrium constant
            Q = self.calculate_reaction_quotient(props, reaction)
            # Solve the equilibrium problem
            result = solver.solve(state, conditions)
            if not result.succeeded():
                # Attempt a cold start
                state = initial_state.clone()
                result = solver.solve(state, conditions)
            if result.succeeded():
                # Update the properties
                props.update(state)
                # Calculate K equilibrium constant
                K = self.calculate_reaction_quotient(props, reaction)

                # Calculate affinity
                R = 8.31446
                A = 2.3026 * R * T * (np.log10(K) - np.log10(Q)) / 1000  # Affinity in kJ
                # Store the affinity (A)
                affinity_kJ.append(A)
            else:
                log.warning(f"Failed to find equilibrium at {P} MPa and {T} K. Filling with NaN.")
                affinity_kJ.append(np.nan)
                state = initial_state.clone()
        # Convert lists to arrays
        affinity_kJ = np.array(affinity_kJ)

        # Log time it took to calculate speciation
        end_time = time.time()
        log.debug(f'{end_time - start_time} seconds to calculate affinity of reaction')

        # Return the filtered results
        return affinity_kJ

    def calculate_reaction_quotient(self, prop, reaction):

        # Initialize numerator and denominator for Q
        Q_numerator = 1.0
        Q_denominator = 1.0

        # Multiply activities raised to their stoichiometric coefficients for products
        for species, coefficient in reaction["products"].items():
            speciesActivity = float(prop.speciesActivity(species))
            Q_numerator *= speciesActivity ** coefficient
        # Multiply activities raised to their stoichiometric coefficients for reactants
        for species, coefficient in reaction["reactants"].items():
            speciesActivity = float(prop.speciesActivity(species))
            Q_denominator *= speciesActivity ** coefficient

        # Calculate the reaction quotient Q
        Q = Q_numerator / Q_denominator
        return Q

    def reaction_parser(self, reaction):
        """
           Parse a chemical reaction string into reactants, products, and optional disequilibrium species.

           Parameters:
           reaction_str (str): The chemical reaction string (e.g., "CO2 + 4 H2(aq) -> CH4(aq) + 2 H2O(aq)").

           Returns:
           dict: Parsed reaction with reactants, products, and optional disequilibrium species.
           """
        reaction_parts = reaction.split("->")
        reactants_str, products_str = reaction_parts[0], reaction_parts[1]

        def parse_side(side_str):
            species_dict = {}
            components = side_str.split("+")
            for component in components:
                component = component.strip()
                if " " in component:
                    coeff, species = component.split(" ", 1)
                    species_dict[species.strip()] = float(coeff)
                else:
                    species_dict[component.strip()] = 1.0
            return species_dict

        reactants = parse_side(reactants_str)
        products = parse_side(products_str)


        return {"reactants": reactants, "products": products}

def temperature_constraint(T_K, System):
    """ Find the pressure constraint at which Reaktoro can find equilibrium for the given speciation and database. Checks if rkt can find equilibrium
        with a pressure of 0.1 MPa at T_K temperature. If it cannot, then returns 0. If it can, then returns 1.
    Args:
        T_K: Initial temperature constraint in K
        System: System tuple with following items: db, system, state, conditions, solver, props

    Returns:
        1 if equilibrium found
        0 if equilibrium not found
    """
    # Initialize the database
    db, system, state, conditions, solver, props = System
    initial_state = state.clone()
    # Establish pressure constraint of 1bar
    conditions.pressure(1, "MPa")
    conditions.temperature(T_K, "K")
    # Solve the equilibrium problem
    result = solver.solve(initial_state, conditions)
    # Check if the equilibrium problem succeeded
    if result.succeeded():
        return 1
    else:
        return 0


def pressure_constraint(P_MPa, System):
    """ Find the pressure constraint at which Reaktoro can find equilibrium for the given speciation and database. Checks if rkt can find equilibrium at 273K at P_MPa pressure.
    If it cannot, then returns 0. If it can, then returns 1.
    Args:
        P_MPa: Pressure to find equilibrium: P_MPa
        System: System tuple with following items: db, system, state, conditions, solver, props

    Returns:
        1 if equilibrium found
        0 if equilibrium not found
    """
    # Disable chemical convergence warnings that Reaktoro raises. We handle these internally instead and throw more specific warnings when they appear.
    rkt.Warnings.disable(906)
    # Initilialize the database
    db, system, state, conditions, solver, props = System
    initial_state = state.clone()
    # Establish pressure constraint of 1bar
    conditions.pressure(P_MPa, "MPa")
    conditions.temperature(Constants.T0, "K")
    # Solve the equilibrium problem
    result = solver.solve(initial_state, conditions)
    # Check if the equilibrium problem succeeded
    if result.succeeded():
        return 1
    else:
        return 0




class RktSeismicOnDemand():
    def __init__(self, aqueous_species_list, speciation_ratio_mol_kg, TMin_K, TMax_K, PMin_MPa, PMax_MPa):
        """
        Initialize the RKtSeismic() object with the minimum and maximum T (in K) and P (in MPa) and speciation info.
        """
        self.aqueous_species_list = aqueous_species_list
        self.speciation_ratio_mol_kg = speciation_ratio_mol_kg
        self.TMin_K = TMin_K
        self.TMax_K = TMax_K
        self.PMin_MPa = PMin_MPa
        self.PMax_MPa = PMax_MPa

    def __call__(self, P_MPa, T_K, grid=False):
        """
        Finds the sound speed in km/s and bulk modulus in GPa for the input of P_MPa and T_K
        Args:
            P_MPa (float, Shape N): Array of pressures
            T_K (float, Shape N): Array of temperatures
            grid: Whether or not they need to be made into a coordinate grid
        Returns:
            VP_kms (float, Shape N): Corresponding sound speeds in km/s
            KS_GPa (float, Shape N): Corresponding bulk modulus in GPa
        """
        if grid:
            P_MPa, T_K = np.meshgrid(P_MPa, T_K, indexing='ij')
        else:
            P_MPa = np.array(P_MPa)
            T_K = np.array(T_K)
            if (np.size(P_MPa) == 0 or np.size(T_K) == 0):
                # If input is empty, return empty array
                return np.array([])
            elif ((np.size(P_MPa) != np.size(T_K)) and not (np.size(P_MPa) == 1 or np.size(T_K) == 1)):
                # If arrays are different lengths, they are probably meant to get a 2D output
                P_MPa, T_K = np.meshgrid(P_MPa, T_K, indexing='ij')
        # Finds the sound speed and associated densities for the given pressure and temperature inputs
        VP_kms, rho_kgm3 = self.seismic_calculations(self.aqueous_species_list, self.speciation_ratio_mol_kg, P_MPa,
                                                     T_K)
        # Calculates the bulk modulus from its relaiontship with sound speed and density
        KS_GPa = rho_kgm3 * VP_kms ** 2 * 1e-3  # 1e-3 because (km/s)^2 * (kg/m^3) gives units of MPa, so 1e-3 to convert to GPa
        # Return sound speed and bulk modulus
        return VP_kms, KS_GPa

    def seismic_calculations(self, aqueous_species_list, speciation_ratio_mol_kg, P_MPa, T_K):
        """
        Calculates the sound speed and densities of the aqueous phase of the solution for the input pressure and temperatures,
        utilizing the thermodynamic properties provided by the supcrt database at equilibrium.
        Args:
            aqueous_species_list: aqueous species in reaction. Should be formatted in one long string with a space in between each species
            speciation_ratio_mol_kg: the ratio of species in the aqueous solution in mol/kg of water. Should be a dictionary with the species as the key and its ratio as its value.
            P_MPa (float, shape N): the desired equilibrium freezing pressure(s) in an array of size N.
            T_K (float, shape N): the desired equilibrium freezing temperature(s) in an array of size N.
        Returns:
            sound_speeds (float, shape N): Associated aqueous sound speeds in km/s
            densities (float, shape N): Associated aqueous densities in km/s
        """

        # Disable chemical convergence warnings that Reaktoro raises. We handle these internally instead and throw more specific warnings when they appear.
        rkt.Warnings.disable(906)
        # Create list that will become array that will hold sound speeds
        sound_speeds = []
        # Create list that will become array taht holds corresponding aqueous densities
        densities = []
        # Create Reaktoro objects
        db, system, initial_state, conditions, solver, props = SupcrtGenerator(aqueous_species_list,
                                                                       speciation_ratio_mol_kg,
                                                                       "mol",
                                                                       CustomSolutionParams.SUPCRT_DATABASE, CustomSolutionParams.SOLID_PHASES_TO_CONSIDER, Constants.PhreeqcToSupcrtNames, CustomSolutionParams.maxIterations)
        state = initial_state.clone()
        # Create an iterator to go through P_MPa and T_K
        it = np.nditer([P_MPa, T_K], flags=['multi_index'])
        # Go through each (P,T) combination, where we assume that P_MPa and T_K are same size and correspond to one another
        for P, T in it:
            P = float(P)
            T = float(T)
            conditions.temperature(T, "K")
            conditions.pressure(P, "MPa")
            # Solve the equilibrium problem
            result = solver.solve(state, conditions)
            # Update the properties
            props.update(state)
            # Check if the equilibrium problem succeeded
            if result.succeeded():
                # If it did succeed, obtain the relevant aqueous only phase props
                aqueous_props = props.phaseProps("AqueousPhase")
                # Obtain density, specific volume, and pressure derivative of specific volume
                density_kg_m3 = float(aqueous_props.density())
                specific_volume_m3_kg = float(aqueous_props.specificVolume())
                pressure_derivative_specific_volume = float(aqueous_props.specificVolumeP())
                # Calculate commpressibility factor, multiplying negative inverse of specific volume by pressure derivative of specific volume
                k_compressibility = -(1 / specific_volume_m3_kg) * pressure_derivative_specific_volume
                # Calculate sound speed, multiplying k by density and taking the product to the -0.5
                c_m_s = (density_kg_m3 * k_compressibility) ** (-0.5)
                # Convert to km/s
                c_km_s = c_m_s / 1000
                # Append to lists
                sound_speeds.append(c_km_s)
                densities.append(density_kg_m3)
                # Reset the state
                state = initial_state.clone()
            # If equilibrium did not succeed, then we will linearly interpolate these values
            # THIS IS VERY UNLIKELY SINCE WE ARE CALCULATING SEISMIC VALUES FOR OCEAN LAYER WHICH HAS ALREADY BEEN CALCULATED BY RKT, SO VALUES SHOULD BE WITHIN RANGE OF RKT
            else:
                # For now, just append np.nan and we will handle later
                sound_speeds.append(np.nan)
                densities.append(np.nan)
                # Reset the state
                state = initial_state.clone()
        # Convert lists to arrays
        sound_speeds = np.array(sound_speeds).reshape(P_MPa.shape)
        densities = np.array(densities).reshape(P_MPa.shape)
        # Check if any values are np.nan, which indicates they were not found in chemical equilibrium
        # THIS IS VERY UNLIKELY SINCE RKT SHOULD BE ABLE TO HANDLE INPUT P AND T RANGE, HOWEVER FOR SELF-CONSISTENCY WE ADD THIS IN CASE
        if np.sum(np.isnan(sound_speeds)) > 0:
            sound_speeds, densities = interpolation_1d(P_MPa, [sound_speeds, densities])
            log.warning("Performed 1d linear interpolation on missing seismic properties")
        return sound_speeds, densities


class RktPhaseOnDemand:
    """
    Class that can find the phase of a given speciation over a range of temperatures and pressures
    """

    def __init__(self, aqueous_species_list, speciation_ratio_mol_kg):
        """
        Initialize the RKtPhaseOnDemand object that will find the phase (liquid or Ice 1) for a input pressure and
        temperature using the Frezchem database. Uses a temperature correction function previously calculated in comparsion with Seafreeze.
        """
        self.aqueous_species_list = aqueous_species_list
        self.speciation_ratio_mol_kg = speciation_ratio_mol_kg
        # Create both frezchem and core Reaktoro systems that can be utilized later on
        self.frezchem = PhreeqcGeneratorForChemicalConstraint(self.aqueous_species_list, self.speciation_ratio_mol_kg,
                                                              "mol",
                                                              CustomSolutionParams.frezchemPath, CustomSolutionParams.maxIterations)
        # Obtain internal temperature correction spline
        # self.temperature_correction_spline = FrezchemFreezingTemperatureCorrectionSplineGenerator()
        self.spline_for_pressures_above_100_MPa, self.PMax_MPa = self.Frezchem_Spline_Generator(
            self.aqueous_species_list, self.speciation_ratio_mol_kg, self.temperature_correction_spline)
        self.calculated_freezing_temperatures = {}

    def __call__(self, P_MPa, T_K, grid=False):
        """
        Call the ice_freezing function for the given input P_MPa and T_K coordinates.
        Importantly, ice_frrezing assumes P_MPa and T_K are the same size and correspond to one another (coordinate pairs),
        so we get them into that format if they are not already.
        Args:
            P_MPa: Array of pressures
            T_K: Array of temperatures
            grid: Whether or not they need to be made into a coordinate grid
        Returns:
            The associated phases of the solution for the input pressure and temperatures
        """
        if not grid:
            np.array(P_MPa)
            np.array(T_K)
        freezing_temperatures = self.rkt_t_freeze(self.aqueous_species_list, self.speciation_ratio_mol_kg, P_MPa,
                                                  self.frezchem, self.PMax_MPa, self.spline_for_pressures_above_100_MPa,
                                                  self.calculated_freezing_temperatures,
                                                  self.temperature_correction_spline)
        if grid:
            freezing_temperatures, T_K = np.meshgrid(freezing_temperatures, T_K, indexing='ij')
        return (T_K < freezing_temperatures).astype(np.int_)

    def Frezchem_Spline_Generator(self, aqueous_species_list, speciation_ratio_mol_kg, temperature_correction_spline,
                                  data_points=30,
                                  significant_threshold=0.1):
        P_MPa = np.linspace(0.1, 100, data_points)
        # Disable chemical convergence warnings that Reaktoro raises. We handle these internally instead and throw more specific warnings when they appear.
        rkt.Warnings.disable(906)
        TMin_K = 220
        TMax_K = 300
        # Create Reaktoro objects
        db, system, initial_state, conditions, solver, props = PhreeqcGeneratorForChemicalConstraint(aqueous_species_list,
                                                                                             speciation_ratio_mol_kg,
                                                                                             "mol",
                                                                                             CustomSolutionParams.frezchemPath, CustomSolutionParams.maxIterations)
        # Create freezing temperatures list and indices of pressures to remove, if necessary
        freezing_temperatures = []
        indices_to_remove = []
        for index, P in enumerate(P_MPa):
            state = initial_state.clone()
            # Specify equilibrium constraints
            conditions.pressure(P, "MPa")
            conditions.set("IP", significant_threshold)
            conditions.setLowerBoundTemperature(TMin_K, "K")
            conditions.setUpperBoundTemperature(TMax_K, "K")
            # Solve the equilibrium problem
            result = solver.solve(state, conditions)
            # Update the properties
            props.update(state)
            # Obtain the equilibrium temperature
            equilibrium_temperature = props.temperature()
            # Correct equilibrium temperature with correction spline
            corrected_temperature = float(equilibrium_temperature) + temperature_correction_spline(P)
            # Check if the result succeeded
            if result.succeeded():
                freezing_temperatures.append(corrected_temperature)
            # If the result failed, then do not include this in list of freezing temperatures and remove from pressure list to avoid its use in creating spline
            else:
                indices_to_remove.append(index)
        # Remove pressures that failed to find equilibrium
        if len(indices_to_remove) > 5:
            log.warning(
                "Reaktoro had a difficult time finding convergence of species composition to generate spline over range of [0.1, 100] MPa. Be warned that results are being extrapolated far beyond Frezchem's ability.")
        P_MPa = np.delete(P_MPa, indices_to_remove)
        # Convert freezing_temperatures to array
        freezing_temperatures = np.array(freezing_temperatures)
        # Create b-spline
        spline = interpolate.make_interp_spline(P_MPa, freezing_temperatures, k=2)
        # Find highest pressure at which Frezchem could find freezing temperature for given composition
        max_pressure = P_MPa[-1]
        # Return spline and max pressure
        return spline, max_pressure

    def rkt_t_freeze(self, aqueous_species_list, speciation_ratio_mol_kg, P_MPa, frezchem, PMax_MPa,
                     freezing_temperature_spline, calculated_freezing_temperatures, temperature_correction_spline,
                     TMin_K=220, TMax_K=300,
                     significant_threshold=0.1):
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
         significant_threshold: the amount of moles of ice present for H2O to be considered in solid phase. Default is 1e-14 moles.

         Returns
         -------
         t_freezing_K: the temperature at which the solution begins to freeze.
         """
        # Disable chemical convergence warnings that Reaktoro raises. We handle these internally instead and throw more specific warnings when they appear.
        # rkt.Warnings.disable(906)
        # Create list that holds boolean values of whether ice is present
        freezing_temperatures = []
        db, system, initial_state, conditions, solver, props = frezchem
        state = initial_state.clone()
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
                # Correct the equilibrium temperature using the temperature correction spline
                equilibrium_temperature = equilibrium_temperature + temperature_correction_spline(P)
                # Check if the result succeeded
                if result.succeeded():
                    state = initial_state.clone()
                # If the result failed, we will use spline
                else:
                    log.warning(f"While attempting to find bottom freezing temperature for pressure of {P} MPa, \n"
                                + f"Reaktoro was unable to find a temperature within range of {TMin_K} K and {TMax_K}.\n"
                                + f"Instead, we will use a spline of freezing temperatures that we generated for this EOS to find the associated freezing temperature.")
                    equilibrium_temperature = freezing_temperature_spline(P)
                    state = initial_state.clone()
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
def ice_freezing(aqueous_species_list, speciation_ratio_mol_per_kg, P_MPa, T_K, PMax_MPa, TMin_K, TMax_K, freezing_temperatures_spline, frezchem_database, core_database, significant_threshold = 0.1):

     Calculates whether ice is freezing at the given temperature and pressure. Utilizes the reaktoro framework to
     constrain the equilibrium position at the prescribed pressure and temperature with the given composition,
     and determines if ice has formed at a significant threshold of 1e-14, therefore calculating the phase of the solution
     and returning true if the phase is ice and false if the phase is liquid.

     Parameters
     ----------
     aqueous_species_list: aqueous species in reaction. Should be formatted in one long string with a space in between each species
     speciation_ratio_mol_per_kg: the ratio of species in the aqueous solution in mol/kg of water. Should be a dictionary
     with the species as the key and its ratio as its value.
     P_MPa (float, shape N): the desired equilibrium freezing pressure(s) in an array of size N.
     T_K (float, shape N): the desired equilibrium freezing temperature(s) in an array of size N.
     significant_threshold: the amount of moles of ice present for H2O to be considered in solid phase. Default is 1e-14 moles.

     Returns
     -------
     ice_present (boolean, shape N): an array of true and falses that indicate whether for prescribed P_MPa and T_K coordinates
        if there is ice present at significant threshold (default of 1e-14)

    # Disable chemical convergence warnings that Reaktoro raises. We handle these internally instead and throw more specific warnings when they appear.
    rkt.Warnings.disable(906)
    # Create list that holds boolean values of whether ice is present
    ice_present = []
    # Create an iterator to go through P_MPa and T_K
    it = np.nditer([P_MPa, T_K], flags =['multi_index'])
    for P, T in it:
        P = float(P)
        T = float(T)
        # Check that pressure is below PMax_MPa, the maximum pressure constraint calculated previously
        if P <= PMax_MPa:
            # Determine whether to use frezchem or core database based on temperature
            if T < 298.15:
                db, system, state, conditions, solver, props, ice_name, database_name = frezchem_database
            else:
                db, system, state, conditions, solver, props, ice_name, database_name = core_database
            # Specify equilibrium conditions
            conditions.pressure(P, "MPa")
            conditions.temperature(T, "K")
            # Solve the equilibrium problem
            result = solver.solve(state, conditions)
            # Update the properties
            props.update(state)
            # Check if the equilibrium problem succeeded
            if result.succeeded():
                # If so, obtain the total amount of solid H2O at equilibrium in moles
                ice_amount = props.speciesAmount(ice_name)
                # If the amount of solid H2O at equilibrium is > the significant threshold, then append true to the list, otherwise append false
                ice_present.append(ice_amount > significant_threshold)
                # Reset the state
                state = reset_state(system, speciation_ratio_mol_per_kg, "mol")
            # Otherwise, the equilibrium problem failed so we need to handle it accordingly
            else:
                # Warn the user about the failed equilibrium (good for debugging)
                log.debug(
                    f"Unsuccessful computation at: {props.pressure() / 1e+6} MPa and {props.temperature()} K.\n"
                    f"The temperature and pressure may be out of bounds, but sometimes Reaktoro fails to converge on specific values even within bounds.\n"
                    f"Will take an alternative approach of rerunning the equilibrium calculation by determining the temperature at which ice begins to forms at the given pressure and composition,\n"
                    f"and if this temperature is below the freezing temperature then we will assume the state is solid, and if the temperature is above the freezing temperature then we will assume the state is liquid.")
                # Instead, we will find the associated freezing temperature with the prescribed pressure and compare that freezing temperature to T
                freezing_temperature = rkt_t_freeze(aqueous_species_list, speciation_ratio_mol_per_kg, P, TMin_K,
                                                    TMax_K, database_name, significant_threshold)
                # If T is <= the freezing temperature for the given pressure, then we will assume ice is also present at this temperature and append True to list
                # Otherwise, if T > the freezing temperature, we will assume ice is not present and append False to list
                ice_present.append(T <= freezing_temperature)  # If T <= freezing_temperature, then we assume ice is present
                # Reset the state
                state = reset_state(system, speciation_ratio_mol_per_kg, "mol")
        # If Pressure is >= PMax, Reaktoro is at its limits of computation and will likely fail. Thus, we utilize our spline approximation of freezing temperatures
        # over a range of pressures that we are confident that Frezchem works at and find the associated freezing temperature at the given pressure.
        # If the temperature we are querying over is less than that freezing temperature, then we assume its ice, otherwise we assume its liquid
        else:
            TFreezing_At_P = freezing_temperatures_spline(P)
            ice_present.append(T < TFreezing_At_P)
    # Convert the ice_present list into an array (of shape N) and return
    ice_present = np.reshape(ice_present, P_MPa.shape)
    return ice_present
    
    

 """
