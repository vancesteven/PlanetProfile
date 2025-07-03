""" File for functions relevant to Custom Solution implementation into PlanetProfile"""
from PlanetProfile.GetConfig import Color, Style, FigMisc
from PlanetProfile.Thermodynamics.Reaktoro.reaktoroProps import MolalConverter, wpptCalculator, SpeciesParser, EOSLookupTableLoader
import numpy as np
import logging
import re
# Assign logger
log = logging.getLogger('PlanetProfile')

def SetupCustomSolution(Planet, Params):
    """
    Configure a Planet's ocean comp and wOcean_ppt based on settings specified in PPCustomSolution.py
    """
    if 'CustomSolution' in Planet.Ocean.comp:
        log.info('Setting up Planet for compatability with CustomSolution implementation. This includes converting string to mols if setting is in grams,\n'
                  'calculating wOcean_ppt if None is passed in, and setting up plot settings and generating EOS data.')
        # Ensure that ocean composition is in molal
        if Params.CustomSolution.SPECIES_CONCENTRATION_UNIT == 'g':
            Planet.Ocean.comp = MolalConverter(Planet.Ocean.comp)
        # Calculate w_ppt for Planet Ocean comp if not specified
        if Planet.Ocean.wOcean_ppt is None or Planet.Ocean.wOcean_ppt < 0:
            # Flag that we are not using wOcean_ppt as independent parameter - used in file name generation
            Planet.Do.USE_WOCEAN_PPT = False
            Planet.Ocean.wOcean_ppt = wpptCalculator(Planet.Ocean.comp.split('=')[1].strip())
        Params = SetupCustomSolutionPlotSettings(np.array(Planet.Ocean.comp), Params)
        SetupCustomSolutionEOS(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt)
    return Planet, Params


def SetupCustomSolutionEOS(CustomSolutionComp, wOcean_ppt):
    """
    Generate/Load a Planet's custom solution EOS. We generate and load here so when GetOceanEOS is called, the file is already loaded.
    """
    # Parse out the species list and ratio into a format compatible with Reaktoro and create a CustomSolution EOS label
    aqueous_species_string, speciation_ratio_mol_kg, ocean_solid_phases, EOS_lookup_label = SpeciesParser(
        CustomSolutionComp, wOcean_ppt)
    # Call function to generate EOS table
    EOSLookupTableLoader(aqueous_species_string, speciation_ratio_mol_kg, ocean_solid_phases, EOS_lookup_label)


def SetupCustomSolutionPlotSettings(PlanetOceanArray, Params):
    """ Setup Color and Linestyle Settings. Namely, We must set iterate through Ocean comp List and add each ocean comp to list of Color and Linestyles
        for CustomSolution
    """

    CustomSolutionOceanComps = [CustomOceanComp for CustomOceanComp in PlanetOceanArray.flatten() if 'CustomSolution' in CustomOceanComp]
    for CustomSolutionOceanComp in CustomSolutionOceanComps:
        # Force native string since numpy strings can be problematic (i.e. '<str_, len() = 344>')
        CustomSolutionOceanComp = str(CustomSolutionOceanComp)
        # Here we need to add the Planets CustomSolution composition to some parameter dictionaries for plotting purposes, which we must do dynamically since input can be anything
        # Add wRef_ppts - namely, we will add the Planet.Ocean.wOcean_ppt and any wRef_ppt in CustomSolution
        if CustomSolutionOceanComp not in Color.cmapName:
            if FigMisc.CustomSolutionSingleCmap:
                Color.cmapName[CustomSolutionOceanComp] = Color.CustomSolutionCmapNames[0]
                
                # Count existing CustomSolution comps in cmapName to determine bounds
                custom_comps = [comp for comp in Color.cmapName.keys() if 'CustomSolution' in comp]
                n_custom_comps = len(custom_comps)
                
                # Set cmapBounds based on number of CustomSolution comps
                for i, comp in enumerate(custom_comps):
                    lower_bound = i / n_custom_comps
                    upper_bound = (i + 1) / n_custom_comps
                    Color.cmapBounds[comp] = [lower_bound, upper_bound]
            else:
                Color.cmapName[CustomSolutionOceanComp] = Color.CustomSolutionCmapNames.pop(0)
                Color.CustomSolutionCmapNames.append(Color.cmapName[CustomSolutionOceanComp])
                Color.cmapBounds[CustomSolutionOceanComp] = Color.cmapBounds["CustomSolution"]
            Color.saturation[CustomSolutionOceanComp] = Color.saturation["CustomSolution"]
            Color.SetCmaps()
            Style.LS[CustomSolutionOceanComp] = Style.LS["CustomSolution"]
            Style.LS_ref[CustomSolutionOceanComp] = Style.LS_ref["CustomSolution"]
            Params.wRef_ppt[CustomSolutionOceanComp] = Params.wRef_ppt["CustomSolution"]
            Params.fNameRef[CustomSolutionOceanComp] = f'{CustomSolutionOceanComp}Ref.txt'
    return Params

def strip_latex_formatting_from_CustomSolutionLabel(s):
    # Remove $...$
    s = re.sub(r'\$+', '', s)
    # Remove \mathrm{...}
    s = re.sub(r'\\mathrm\{(.*?)\}', r'\1', s)
    # Remove other latex commands like \textbf{...}
    s = re.sub(r'\\[a-zA-Z]+\{(.*?)\}', r'\1', s)
    # Remove remaining braces
    s = s.replace('{', '').replace('}', '')
    return s