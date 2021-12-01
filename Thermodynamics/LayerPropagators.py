from Utilities.dataStructs import Constants
from Thermodynamics.SwEOSChooser import GetModgswPressureFreezingCT as p_freezing

def IceLayers(Planet):
    """ Layer propagation from surface downward through the ice using geophysics.
        Iteratively sets up the thermal profile (the density and temperature)
        of the layer with each pressure step for all ice layers including
        ice Ih, ice III, ice V by calling the necessary subfunctions
        Args:
            Planet(object): instance of the Planet class
            Layers (array) : empty array to be filled in with layer data
            Constants (struct) : structure of constants needed
        Returns:
            Layers : (array) density and temperature of the layer with each pressure step

        Examples:
    """

    Planet.IceILayers()

    Planet = IceILayers(Planet)
# I am a bit cnfused how we are trying ot implement this- it looks like we are overwriting the layers array
#with just the bottom layers
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

        Args:
            Planet(object): instance of the Planet class
            Layers (array) : empty array to be filled in with layer data
            Constants (struct) : structure of constants needed
        Returns:
            Layers : (array) density and temperature of the layer with each pressure step
    """

    Planet.Pb_MPa = getPfreeze(Planet.Tb_K, wo , Planet.Ocean.comp)
    deltaP = Planet.Pb_MPa/(Planet.Steps.nClath + Planet.Steps.n_IceI - 1)

    if Planet.Steps.nClath is not None:
        ii = Planet.Steps.nClath + 1
    else:
        ii = 2

    for il = [ii :(Planet.Steps.n_IceI + Planet.Steps.nClath)]:  # propagate P,T,rho to the bottom of the ice
        Planet.P_MPa(il) = Planet.P_MPa(il-1) + deltaP
        Planet.T_K(il) = (Planet.Bulk.Tb_K**(Planet.P_MPa(il) / Planet.Pb_MPa) * (Planet.Bulk.Tsurf_K**(1-Planet.P_MPa(il) / Planet.Pb_MPa))
        Planet.rho_kgm3(il) = getRhoIce(Planet.P_MPa(il), Planet.T_K(il), 1)
        [Planet.Cp(il), Planet.alpha_perK(il)] = getCpIce(Planet.P_MPa(il), Planet.T_K(il), Planet.phase(il))

    Planet.Steps.nIceIIILitho = 0
    Planet.Steps.nIceVLitho = 0
    #Planet.PbI_MPa = Planet.Pb_MPa
    return Planet


def IceIIIUnderplateLayers(Planet):
    """ For cold, thick ice shells, model ice III under ice I layer
        Calculates the density and temperature of the layer with each pressure step

        Args:
            Planet(object): instance of the Planet class
            Layers (array) : empty array to be filled in with layer data
            Constants (struct) : structure of constants needed
        Returns:
            Layers : (array) density and temperature of the layer with each pressure step

        Examples:
    """

    return Planet


def IceVUnderplateLayers(Planet):
    """ For cold, thick ice shells, model ice V and ice III under ice I layer
        Calculates the density and temperature of the layer with each pressure step

        Args:
            Planet(object): instance of the Planet class
            Layers (array) : empty array to be filled in with layer data
            Constants (struct) : structure of constants needed
        Returns:
            Layers : (array) density and temperature of the layer with each pressure step

        Examples:
    """

    return Planet


def OceanLayers(Planet):
    """ Geophysical and thermodynamic calculations for ocean layer
        Calculates the density and temperature of the layer with each pressure step

        Args:
            Planet(object): instance of the Planet class
            Layers (array) : empty array to be filled in with layer data
            Constants (struct) : structure of constants needed
        Returns:
            Layers : (array) density and temperature of the layer with each pressure step

        Examples:
    """

    return Planet


def PlanetDepths(Planet):
    """ Convert from organization by radius into organization by depth
        Calculates the density and temperature of the layer with each pressure step

        Args:
            Planet(object): instance of the Planet class
            Layers (array) : empty array to be filled in with layer data
            Constants (struct) : structure of constants needed
        Returns:
            Layers : (array) density and temperature of the layer with each pressure step

        Examples:
    """

    return Planet


def InnerLayers(Planet):
    """ Geophysical and thermodynamic calculations for silicate and core layers
        Calculates the density and temperature of the layer with each pressure step

        Args:
            Planet(object): instance of the Planet class
            Layers (array) : empty array to be filled in with layer data
            Constants (struct) : structure of constants needed
        Returns:
            Layers : (array) density and temperature of the layer with each pressure step

        Examples:
    """

    Planet = SilicateLayers(Planet)

    if Planet.Do.FeCORE:
        Planet = IronCoreLayers(Planet)

    return Planet


def SilicateLayers(Planet):
    """ Geophysical and thermodynamic calculations for silicate layers
        Calculates the density and temperature of the layer with each pressure step

        Args:
            Planet(object): instance of the Planet class
            Layers (array) : empty array to be filled in with layer data
            Constants (struct) : structure of constants needed
        Returns:
            Layers : (array) density and temperature of the layer with each pressure step

        Examples:
    """

    Planet.Steps.indSil = 0
    Planet.Steps.nTotal = Planet.Steps.indSil + Planet.Steps.nSil
    return Planet


def IronCoreLayers(Planet):
    """ Geophysical and thermodynamic calculations for core layers
        Calculates the density and temperature of the layer with each pressure step

        Args:
            Planet(object): instance of the Planet class
            Layers (array) : empty array to be filled in with layer data
            Constants (struct) : structure of constants needed
        Returns:
            Layers : (array) density and temperature of the layer with each pressure step

        Examples:
    """

    Planet.Steps.nTotal += Planet.Steps.nCore
    return Planet




def getPfreeze(Planet):
    """ Returns the pressure at which the ocean freezes based on temperature, salinity, and ocean composition

        Args:
            T_K (float) : temperature in Kelvin of the layer
            w_o (float) : salinity in weight percent
            Planet(object): instance of the Planet class
        Returns:
            Layers : (array) density and temperature of the layer with each pressure step

        Examples:
    """
    if Planet.Ocean.comp == "MgSO4":
        import sympy
        Pfreeze_MPa = scipy.optimize.root_scalar(L_IceMgSO4, args=(Planet.P, Planet.T_K, Planet.Ocean.wOcean_ppt, 1), bracket=[0, 250],  method='bisect')
        return Pfreeze_MPa


    elif Planet.Ocean.comp == "Seawater":
        Pfreeze_MPa = 0.1 * p_freezing(w_o, T_K)
        return Pfreeze_MPa

    elif Planet.Ocean.comp == "NaCl":
        pass
    elif Planet.Ocean.comp == "NH3":
        pass
#Note to self: Pfreeze is a case-by-case function, there are also PfreezeIII and PfreezeV in the MATLAB code currently


'''
def getIcePhaseMgSO4(P,T,w):

    # phase = getIcePhaseMgSO4(P,T,w)
    # P in MPa
    # T in K
    # w in Wt%
    xH2O = 1 - 1  / (1 + (1 - 0.01 * w)  / (0.01.*w) * 120.3686 / 18.0142)

dmu = deltaMuILiqMgSO4(P,T,xH2O);

phase = find(dmu == min(dmu));

if dmu(phase) > 0
    phase = 0;
end

phase(phase==5) = 6;
phase(phase==4) = 5;

def deltaMuILiqMgSO4(P,T,x):
R = 8.314;
dmu = MuIceLiqH2O(P,T) - R./0.018.*T.*(activity(P,T,x,-1.8e6,150,1.45e-4,-12,246)*x);% this line contains the coefficients from Vance et al. 2014 for MgSO4

function y = activity(P,T,x,w0,w1,w2,w3,To)
R = 8.314./0.018;
y = (((1-x).^2).*w0.*(1+w1.*tanh(w2.*P)).*(1+w3./(T-To)^2)./(R.*T));
'''
