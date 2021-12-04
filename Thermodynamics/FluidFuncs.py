from seafreeze import seafreeze as SeaFreeze

def FluidEOS(compstr, w_ppt, P_MPa, T_K):
    """ Returns mass density, heat capacity, and thermal expansivity based on thermodynamics
     from SeaFreeze and input pressure, temperature, salinity, and composition

        Args:
            compstr (string): Composition of dissolved salt
            w_ppt (float): Salinity of fluid in ppt
            P_MPa (float, shape N): Pressure of the fluid in MPa
            T_K (float, shape N): Temperature of the fluid in K
        Returns:
            rho_kgm3 (float, shape N): Density of fluid in kg/m^3
            Cp_JkgK (float, shape N): Heat capacity at constant pressure of fluid in J/kg/K
            alpha_pK (float, shape N): Thermal expansivity of fluid in K^-1
            vFluid_kms (float, shape N): Sound speeds in fluid in km/s
    """
    if w_ppt != 0: raise ValueError('SeaFreeze only applies to pure water so far. Set Planet.Bulk.wOcean_ppt = 0.')

    # Arrange input data into (P,T) value pair tuples compatible with SeaFreeze
    nPTs = np.size(P_MPa)
    PTpairs = np.empty((nPTs,), dtype=object)
    for i in range(nPTs):
        PTpairs[i] = (P_MPa[i], T_K[i])

    if compstr == 'Seawater':
        raise ValueError('Unable to set FluidEOS. Only NH3 is implemented so far.')
    elif compstr == 'NH3':
        seaOutEOS = SeaFreeze(PTpairs, 'water1')
    elif compstr == 'MgSO4':
        raise ValueError('Unable to set FluidEOS. Only NH3 is implemented so far.')
    elif compstr == 'NaCl':
        raise ValueError('Unable to set FluidEOS. Only NH3 is implemented so far.')
    else:
        raise ValueError('Unable to set FluidEOS. compstr="'+compstr+'" but options are Seawater, NH3, MgSO4, and NaCl.')

    rho_kgm3 = seaOutEOS.rho
    Cp_JkgK = seaOutEOS.Cp
    alpha_pK = seaOutEOS.alpha
    vFluid_kms = seaOutEOS.vel/1e3
    return rho_kgm3, Cp_JkgK, alpha_pK, vFluid_kms


def GetPfreeze(compstr, w_ppt, T_K):
    """ Returns the pressure at which the fluid freezes based on temperature, salinity, and composition

        Args:
            compstr (string): Composition of dissolved salt
            w_ppt (float): Salinity of fluid in ppt
            T_K (float): Temperature of the fluid in K
        Returns:
            Layers : (array) density and temperature of the layer with each pressure step
    """
    if compstr == 'Seawater':
        pass#seaOut = SeaFreeze([P_MPa, T_K], 'Ih')
    elif compstr == 'NH3':
        raise ValueError('Unable to GetPfreeze. NH3 is not implemented yet.')
    elif compstr == 'MgSO4':
        raise ValueError('Unable to GetPfreeze. MgSO4 is not implemented yet.')
    elif compstr == 'NaCl':
        raise ValueError('Unable to GetPfreeze. NaCl is not implemented yet.')
    else:
        raise ValueError('Unable to GetPfreeze. compstr="'+compstr+'" but options are Seawater, NH3, MgSO4, and NaCl.')

    Pfreeze_MPa = 100
    return Pfreeze_MPa


def GetTfreeze(compstr, w_ppt, P_MPa, T_K):
    """ Returns the pressure at which the fluid freezes based on temperature, salinity, and composition

        Args:
            compstr (string): Composition of dissolved salt
            w_ppt (float): Salinity of fluid in ppt
            P_MPa (float): Pressure of the fluid in MPa
            T_K (float): Temperature of the fluid in K
        Returns:
            Layers : (array) density and temperature of the layer with each pressure step
    """
    if compstr == 'Seawater':
        Tfreeze_K = SeaFreeze([P_MPa, T_K], 'Ih')
    elif compstr == 'NH3':
        raise ValueError('Unable to GetPfreeze. NH3 is not implemented yet.')
    elif compstr == 'MgSO4':
        raise ValueError('Unable to GetPfreeze. MgSO4 is not implemented yet.')
    elif compstr == 'NaCl':
        raise ValueError('Unable to GetPfreeze. NaCl is not implemented yet.')
    else:
        raise ValueError('Unable to GetPfreeze. compstr="'+compstr+'" but options are Seawater, NH3, MgSO4, and NaCl.')

    return Tfreeze_K

