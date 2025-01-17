""" File for functions relevant to Custom Solution"""
from PlanetProfile.GetConfig import Color, Style, FigLbl, FigMisc
from PlanetProfile.Thermodynamics.Reaktoro.reaktoroProps import MolalConverter, wpptCalculator
import numpy as np
from PlanetProfile.Utilities.defineStructs import Constants, EOSlist
from PlanetProfile.Thermodynamics.Reaktoro.reaktoroProps import EOSLookupTableLoader

def SetupCustomSolution(Planet, Params):
    # Ensure that ocean composition is in molal
    if Params.CustomSolution.SPECIES_CONCENTRATION_UNIT == 'g':
        Planet.Ocean.comp = MolalConverter(Planet.Ocean.comp)
    # Calculate w_ppt for Planet Ocean comp if not specified
    if Planet.Ocean.wOcean_ppt is None or Planet.Ocean.wOcean_ppt < 0:
        # Flag that we are not using wOcean_ppt as independent parameter - used in file name generation
        Planet.Do.USE_WOCEAN_PPT = False
        Planet.Ocean.wOcean_ppt = wpptCalculator(Planet.Ocean.comp.split('=')[1].strip())
    SetupCustomSolutionPlotSettings(np.array(Planet.Ocean.comp), Params)
    return Planet, Params


def SaveEOSToDisk(EOSList):
    """
    Save EOS to disk, so we can load from disk rather than having to re-generate next time.
    This currently applies to CustomSolution, where we generate EOS during Profile run and save to EOSList
    """
    for EOS in EOSlist.loaded['CustomSolutionEOS']:
        # If we have an EOSLookupTableLoader, then we need to save to disk
        if isinstance(EOSlist.loaded['CustomSolutionEOS'][EOS], EOSLookupTableLoader):
            EOSlist.loaded['CustomSolutionEOS'][EOS].saveEOSToDisk()



def SetupCustomSolutionPlotSettings(PlanetOceanArray, Params):
    """ Setup Color and Linestyle Settings. Namely, We must set iterate through Ocean comp List and add each ocean comp to list of Color and Linestyles
        for CustomSolution
    """

    CustomSolutionOceanComps = [CustomOceanComp for CustomOceanComp in PlanetOceanArray.flatten() if 'CustomSolution' in CustomOceanComp]
    for CustomSolutionOceanComp in CustomSolutionOceanComps:
        # Here we need to add the Planets CustomSolution composition to some parameter dictionaries for plotting purposes, which we must do dynamically since input can be anything
        # Add wRef_ppts - namely, we will add the Planet.Ocean.wOcean_ppt and any wRef_ppt in CustomSolution
        if CustomSolutionOceanComp not in Color.cmapName:
            Color.cmapName[CustomSolutionOceanComp] = Color.CustomSolutionCmapNames.pop(0)
            Color.CustomSolutionCmapNames.append(Color.cmapName[CustomSolutionOceanComp])
            Color.cmapBounds[CustomSolutionOceanComp] = Color.cmapBounds["CustomSolution"]
            Color.saturation[CustomSolutionOceanComp] = Color.saturation["CustomSolution"]
            Color.SetCmaps()
            Style.LS[CustomSolutionOceanComp] = Style.LS["CustomSolution"]
            Style.LS_ref[CustomSolutionOceanComp] = Style.LS_ref["CustomSolution"]
            Params.wRef_ppt[CustomSolutionOceanComp] = Params.wRef_ppt["CustomSolution"]
            Params.fNameRef[CustomSolutionOceanComp] = f'{CustomSolutionOceanComp}Ref.txt'