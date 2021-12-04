from Utilities.dataStructs import Constants
from Thermodynamics.FluidFuncs import GetPfreeze, GetTfreeze

def IceLayers(Planet):
    """ Layer propagation from surface downward through the ice using geophysics.
        Iteratively sets up the thermal profile (the density and temperature)
        of the layer with each pressure step for all ice layers including
        ice Ih, ice III, ice V by calling the necessary subfunctions

        Requires Planet attributes:
            ???
        Assigns Planet attributes:
            ???
    """

    Planet = IceILayers(Planet)

    if Planet.Do.BOTTOM_ICEV:
        Planet = IceVUnderplateLayers(Planet)
    elif Planet.Do.BOTTOM_ICEIII:
        Planet = IceIIIUnderplateLayers(Planet)
        Planet.Steps.nIceVLitho = 0
    else:
        Planet.Steps.nIceIIILitho = 0
        Planet.Steps.nIceVLitho = 0

    return Planet


def IceILayers(Planet):
    """ Geophysical and thermodynamic calculations for outermost ice layer
        Calculates the density and temperature of the layer with each pressure step

        Requires Planet attributes:
            Tb_K
            Ocean.comp
            Ocean.wOcean_ppt
        Assigns Planet attributes:
            ???
    """

    Planet.Pb_MPa = GetPfreeze(Planet.Ocean.comp, Planet.Ocean.wOcean_ppt, Planet.Bulk.Tb_K)

    if Planet.Steps.nClath is not None:
        ii = Planet.Steps.nClath + 1
    else:
        ii = 2


    #Planet.PbI_MPa = Planet.Pb_MPa
    return Planet


def IceIIIUnderplateLayers(Planet):
    """ For cold, thick ice shells, model ice III under ice I layer
        Calculates the density and temperature of the layer with each pressure step

        Requires Planet attributes:
            ???
        Assigns Planet attributes:
            ???
    """

    return Planet


def IceVUnderplateLayers(Planet):
    """ For cold, thick ice shells, model ice V and ice III under ice I layer
        Calculates the density and temperature of the layer with each pressure step

        Requires Planet attributes:
            ???
        Assigns Planet attributes:
            ???
    """

    return Planet


def OceanLayers(Planet):
    """ Geophysical and thermodynamic calculations for ocean layer
        Calculates the density and temperature of the layer with each pressure step

        Requires Planet attributes:
            ???
        Assigns Planet attributes:
            ???
    """

    return Planet


def PlanetDepths(Planet):
    """ Convert from organization by radius into organization by depth
        Calculates the density and temperature of the layer with each pressure step

        Requires Planet attributes:
            ???
        Assigns Planet attributes:
            ???
    """

    return Planet


def InnerLayers(Planet):
    """ Geophysical and thermodynamic calculations for silicate and core layers
        Calculates the density and temperature of the layer with each pressure step

        Requires Planet attributes:
            ???
        Assigns Planet attributes:
            ???
    """

    Planet = SilicateLayers(Planet)

    if Planet.Do.FeCORE:
        Planet = IronCoreLayers(Planet)

    return Planet


def SilicateLayers(Planet):
    """ Geophysical and thermodynamic calculations for silicate layers
        Calculates the density and temperature of the layer with each pressure step

        Requires Planet attributes:
            ???
        Assigns Planet attributes:
            ???
    """

    Planet.Steps.indSil = 0
    Planet.Steps.nTotal = Planet.Steps.indSil + Planet.Steps.nSil
    return Planet


def IronCoreLayers(Planet):
    """ Geophysical and thermodynamic calculations for core layers
        Calculates the density and temperature of the layer with each pressure step

        Requires Planet attributes:
            ???
        Assigns Planet attributes:
            ???
    """

    Planet.Steps.nTotal += Planet.Steps.nCore
    return Planet
