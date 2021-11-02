from Utilities.dataStructs import Constants

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

        Examples:
    """


    ''' for il=ii:(n_iceI(iT)+n_clath(iT)) % propagate P,T,rho to the bottom of the ice % THIS CAN BE COMPUTED AS A VECTOR INSTEAD.
            P_MPa(iT,il) = P_MPa(iT,il-1) + deltaP;
            T_K(iT,il) = (Planet.Tb_K(iT).^(P_MPa(iT,il)./Pb_MPa(iT))).*(Planet.Tsurf_K.^(1-P_MPa(iT,il)./Pb_MPa(iT)));
            rho_kgm3(iT,il) = getRhoIce(P_MPa(iT,il),T_K(iT,il),1); % I THINK THIS CALL CAN ACCEPT A VECTOR
            [Cp(iT,il) alpha_K(iT,il)]= getCpIce(P_MPa(iT,il),T_K(iT,il),phase(iT,il)) ; % I THINK THIS CALL CAN ACCEPT A VECTOR

            %rho_kgm3(iT,il) = getRhoIce(P_MPa(iT,il),T_K(iT,il),1);
        end


        nIceIIILithosphere=0;
        nIceVLithosphere=0;
        PbI_MPa = Pb_MPa;
    catch
        disp('PlanetProfile failed to get Pb! Here''s why:')
        disp(lasterr)
        disp('Maybe it''s okay? If execution stopped, then probably not. Try looking where getPfreeze is called.')
    '''


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
