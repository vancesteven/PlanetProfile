from Utilities.dataStructs import Constants

def IceLayers(Planet, Layers):
    """ Layer propagation from surface downward through ice using geophysics """
    Layers = IceILayers(Planet, Layers)

    if Planet.BOTTOM_ICEV:
        Layers = IceVUnderplateLayers(Planet, Layers)
    elif Planet.BOTTOM_ICEIII:
        Layers = IceIIIUnderplateLayers(Planet, Layers)
        Planet.nIceVLitho = 0
    else:
        Planet.nIceIIILitho = 0
        Planet.nIceVLitho = 0

    return Layers


def IceILayers(Planet, Layers):
    """ Geophysical and thermodynamic calculations for outermost ice layer """

    return Layers


def IceIIIUnderplateLayers(Planet, Layers):
    """ For cold, thick ice shells, model ice III under ice I layer """

    return Layers


def IceVUnderplateLayers(Planet, Layers):
    """ For cold, thick ice shells, model ice V and ice III under ice I layer """

    return Layers


def OceanLayers(Planet, Layers):
    """ Geophysical and thermodynamic calculations for ocean layer """

    return Layers


def PlanetDepths(Planet, Layers):
    """ Convert from organization by radius into organization by depth """

    return Layers


def InnerLayers(Planet, Layers):
    """ Geophysical and thermodynamic calculations for silicate and core layers """

    nStepsSilicate, Layers = SilicateLayers(Planet, Layers)
    Planet.nStepsTotal = Planet.nStepsHydro + nStepsSilicate

    if Planet.Core.FeCORE:
        Layers = IronCoreLayers(Planet, Layers)
        Planet.nStepsTotal += Planet.nStepsCore

    return Planet, Layers


def SilicateLayers(Planet, Layers):
    """ Geophysical and thermodynamic calculations for silicate layers """
    nStepsSilicate = 0

    return nStepsSilicate, Layers


def IronCoreLayers(Planet, Layers):
    """ Geophysical and thermodynamic calculations for core layers """

    return Layers