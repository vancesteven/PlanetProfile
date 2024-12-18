from PlanetProfile.GetConfig import Color, Style, FigLbl, FigSize, FigMisc
from PlanetProfile.Thermodynamics.Reaktoro.reaktoroProps import FileSetupFromConfig


def ReaktoroConfigAdjustments(Planet, Params, ROOT):
    Params.wRef_ppt[Planet.Ocean.comp] = Params.wRef_ppt["CustomSolution"]
    Params.fNameRef[Planet.Ocean.comp] = f'{Planet.Ocean.comp}Ref.txt'

    Color.cmapName[Planet.Ocean.comp] = Color.cmapName["CustomSolution"]
    Color.cmapBounds[Planet.Ocean.comp] = Color.cmapBounds["CustomSolution"]
    Color.saturation[Planet.Ocean.comp] = Color.saturation["CustomSolution"]
    Color.SetCmaps()
    Style.LS[Planet.Ocean.comp] = Style.LS["CustomSolution"]
    Style.LS_ref[Planet.Ocean.comp] = Style.LS_ref["CustomSolution"]
    Params.CustomSolution.setPaths(ROOT)
    FileSetupFromConfig(Params.CustomSolution)

    return Params
