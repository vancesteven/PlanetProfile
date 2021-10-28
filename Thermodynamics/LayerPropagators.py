def IceLayers(Planet, Layers, Constants):
    """ Layer propagation from surface downward through ice using geophysics """
    Layers = IceILayers(Planet, Layers, Constants)

    if Planet.BOTTOM_ICEIII:
        Layers = IceIIIUnderplateLayers(Planet, Layers, Constants)
    elif Planet.BOTTOM_ICEV:
        Layers = IceVUnderplateLayers(Planet, Layers, Constants)

    return Layers


def IceILayers(Planet, Layers, Constants):
    """ Geophysical and thermodynamic calculations for outermost ice layer """

    return Layers


def IceIIIUnderplateLayers(Planet, Layers, Constants):
    """ For cold, thick ice shells, model ice III under ice I layer """

    return Layers


def IceVUnderplateLayers(Planet, Layers, Constants):
    """ For cold, thick ice shells, model ice V and ice III under ice I layer """

    return Layers


def OceanLayers(Planet, Layers, Constants):
    """ Geophysical and thermodynamic calculations for ocean layer """

    return Layers


def PlanetDepths(Planet, Layers, Constants):
    """ Convert from organization by radius into organization by depth """

    return Layers


def InnerLayers(Planet, Layers, Constants):
    """ Geophysical and thermodynamic calculations for silicate and core layers """

    nStepsSilicate, Layers = SilicateLayers(Planet, Layers, Constants)
    Planet.nStepsTotal = Planet.nStepsHydro + nStepsSilicate

    if Planet.Core.FeCORE:
        Layers = IronCoreLayers(Planet, Layers, Constants)
        Planet.nStepsTotal += Planet.nStepsCore

    return Planet, Layers


def SilicateLayers(Planet, Layers, Constants):
    """ Geophysical and thermodynamic calculations for silicate layers """
    nStepsSilicate = 0

    return nStepsSilicate, Layers


def IronCoreLayers(Planet, Layers, Constants):
    """ Geophysical and thermodynamic calculations for core layers """

    return Layers