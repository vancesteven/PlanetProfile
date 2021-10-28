"""swEOS_chooser.py-> chooses between which dissolved salt you are looking at
Using Planet.Ocean.comp and parameters for each function, returns or grab needed EOS for core, water, ice, etc.
Establish EOS/ OceanEOS that gives property need-> SetupEOS
"""


import gsw

#if Planet.Ocean.comp == 'NH3':
#    raise ValueError(['NH3 is not currently implemented.'])

#if Planet.Ocean.comp == 'NaCl':
#    raise ValueError(['NaCl is not currently implemented.'])
#Note to self: define OceanComposition variable and replace Planet-profile specific defs

def GetModgswAdiabaticLapseRate(planetOceanComp, SA, P, T = None, CT = None):
    """Calculates the adiabatic lapse rate of sea water

        Args:
            planetOceanComp(string): seawater composition, pass in as Planet.Ocean.comp
            SA (MxN array-like, float): absolute salinity [g/kg]
            T (MxN array-like, float) : in situ temperature in degrees celsius [o^C]
            CT (MxN array-like, float): Conservative temperature [o^C (ITS-90)]
            P(1x1, Mx1, 1xN, MxN array-like, float) : pressure in absolute bars [dbar]
        Returns:
            adiabatic_lapse_rate : (array_like) [K/Pa]

        Examples:
    """

    if planetOceanComp == 'Seawater':
        if TA is not None:
            return gsw.adiabatic_lapse_rate_t_exact(SA,T-T0,10*(P-P0/1e5))*1e5
        if CT is not None:
            return gsw.adiabatic_lapse_rate_from_CT(SA, CT, P)
        raise ValueError('Adiabatic Lapse Rate needs either T or CT.')
    elif planetOceanComp == 'MgSO4':
        pass



def GetModgswAlpha(planetOceanComp, SA, P, T = None, CT = None, PT = None):
    """Calculates Alpha, the thermal expansion coefficient of seawater with
       respect to in situ temperature.

        Args:
            planetOceanComp(string): seawater composition, pass in as Planet.Ocean.comp
            SA (MxN array-like, float): absolute salinity [g/kg]
            T (MxNarray-like, float) : in situ temperature in degrees celsius [o^C]
            P(1x1, Mx1, 1xN, MxN array-like, float) : pressure in absolute bars [dbar]
            CT (MxN array-like, float): conservative temperature [o^C (ITS-90)]
            PT : potential temperature, with a reference pressure of zero [o^C (ITS-90)]
        Returns:
            alpha(array-like):  thermal expansion coefficient [1/K]

            Examples:
    """
    if planetOceanComp == 'Seawater':
        if T is not None:
            return gsw.alpha_wrt_t_exact(SA, T, P)
        if CT is not None:
            return gsw.alpha_wrt_CT_t_exact(SA, CT, P)
        if PT is not None:
            return gsw.alpha_wrt_pt_t_exact(SA, T, P)
    elif planetOceanComp == 'MgSO4':
        pass

def GetModgswBeta(planetOceanComp, SA, P, T = None, CT = None, PT = None):
    """Calculates the saline (i.e. haline) contraction coefficient beta of seawater
       at constant temperature.

       Args:
        planetOceanComp(string): seawater composition, pass in as Planet.Ocean.comp
        SA (MxN array-like, float): absolute salinity [g/kg]
        T (MxNarray-like, float) : in situ temperature in degrees celsius [o^C]
        P(1x1, Mx1, 1xN, MxN array-like, float) : pressure in absolute bars [dbar]
        CT (MxN array-like, float): conservative temperature [o^C (ITS-90)]
        PT : potential temperature, with a reference pressure of 0 dbar [o^C (ITS-90)]
        Returns:
        beta(array-like):  saline (i.e. haline) contraction coefficient [kg/g]

            Examples:
    """
    if planetOceanComp == 'Seawater':
        if T is not None:
            return gsw.beta_const_t_exact(SA, T, P)
        if CT is not None:
            return gsw.beta_const_CT_t_exact(SA, CT, P)
        if PT is not None:
            return gsw.beta_const_pt_t_exact(SA, T, P)
    elif planetOceanComp == 'MgSO4':
        pass



def GetModgswCFromSA(planetOceanComp, SA, T, P):
    """Calculates conductivity, C, from (SP, t, p) using PSS-78 in the range
    2 < SP < 42.

    Args:
    planetOceanComp(string): seawater composition, pass in as Planet.Ocean.comp
    SA (MxN array-like, float): absolute salinity [g/kg]
    T (MxNarray-like, float) : temperature in degrees celsius [o^C]
    P(1x1, Mx1, 1xN, MxN array-like, float) : pressure in absolute bars [dbar]
    Returns:
    C(array-like): conductivity [S/m]

    Examples:
    """
    if planetOceanComp == 'Seawater':
        #First, we check if we are in correct range of SP
        SP = SP_from_SA(SA, 10.1325, 0, 0) # set lat lon to 0
        if (SP<2) or (SP>42):
            raise ValueError('Practical Salinity SP must be between 2 and 42')
        else:
            return 0.1*gsw.C_from_SP(SP, T-gsw_T0, 10*(P-gsw_P0/1e5))


    #NOTE: C_from_SA not defined on own in gsw, but C_from_SP is
    elif planetOceanComp == 'MgSO4':
        pass

def GetModgswCPtExact(planetOceanComp, SA, T, P):
    """Calculates the isobaric heat capacity of seawater.

       Args:
        planetOceanComp(string): seawater composition, pass in as Planet.Ocean.comp
        SA (MxN array-like, float): absolute salinity [g/kg]
        T (MxNarray-like, float) : in situ temperature in degrees celsius [o^C]
        P(1x1, Mx1, 1xN, MxN array-like, float) : pressure in absolute bars [dbar]
        Returns:
        CP(array-like): isobaric heat capacity of seawater. [J/(kg*K)]

            Examples:
    """
    if planetOceanComp == 'Seawater':
        return gsw.cp_t_exact(SA, T, P)
    elif planetOceanComp == 'MgSO4':
        pass


def GetModgswRhotExact(planetOceanComp, SA, T, P):
    """Calculates in situ density of seawater from absolute salinity
        and in situ temperature.

       Args:
        planetOceanComp(string): seawater composition, pass in as Planet.Ocean.comp
        SA (MxN array-like, float): absolute salinity [g/kg]
        T (MxNarray-like, float) : in situ temperature in degrees celsius [o^C]
        P(1x1, Mx1, 1xN, MxN array-like, float) : pressure in absolute bars [dbar]
        Returns:
        rho_t_exact(array-like):  in situ density of seawater [kg/m^3]

            Examples:
    """
    if planetOceanComp == 'Seawater':
        return gsw.rho_t_exact(SA,T-gsw_T0,10*(P-gsw_P0/1e5))
    elif planetOceanComp == 'MgSO4':
        pass
def LatentHeatMelting(SA,P):
    """   NOT YET DEFINED IN PYTHON-GSW :(

    Args:
        SA (MxN array-like, float): absolute salinity [g/kg]
        T (MxNarray-like, float) : temperature in degrees celsius [o^C (ITS-90)]
        P(1x1, Mx1, 1xN, MxN array-like, float) : pressure in absolute bars [dbar]
        Returns:
        L_m(list?float?): ???

    Examples:
    """
    #L_m = need to find this in gsw
    #return  L_m
    pass

def GetModgswPotRhotExact(planetOceanComp, SA, T, P, P_ref = 0):
    """Calculates potential density of seawater.

       Args:
        planetOceanComp(string): seawater composition, pass in as Planet.Ocean.comp
        SA (MxN array-like, float): absolute salinity [g/kg]
        T (MxNarray-like, float) : in situ temperature in degrees celsius [o^C]
        P(1x1, Mx1, 1xN, MxN array-like, float) : pressure in absolute bars [dbar]
        P_ref(int, float, optional) :  reference pressure, default = 0 [dbar]
        Returns:
        pot_rho(array-like):  potential density of seawater [kg/m^3]

            Examples:
    """
    if planetOceanComp == "Seawater":
        return gsw.pot_rho_t_exact(SA,T-gsw_T0,10*(P-gsw_P0/1e5),10*(P_ref-gsw_P0/1e5));

    elif planetOceanComp == "MgSO4":
        pass


def GetModgswPTFromT(planetOceanComp, SA, T, P, P_ref = 0):
    """Calculates potential temperature with the general reference pressure,
       P_ref, from in situ temperature.

       Args:
        planetOceanComp(string): seawater composition, pass in as Planet.Ocean.comp
        SA (MxN array-like, float): absolute salinity [g/kg]
        T (MxNarray-like, float) : in situ temperature in degrees celsius [o^C]
        P(1x1, Mx1, 1xN, MxN array-like, float) : pressure in absolute bars [dbar]
        P_ref(int, float, optional) :  reference pressure, default = 0 [dbar]
        Returns:
        pt(array-like):  potential temperature [oC (ITS-90)]

            Examples:
    """
    if planetOceanComp == "Seawater":
        return T0 + gsw.pt_from_t(SA,T-T0,10*(P-P0/1e5),10*(P_ref-P0/1e5))

    elif planetOceanComp == "MgSO4":
        pass

def GetModgswTFreezing(planetOceanComp, SA, T, saturation_fraction = 1):
    """ Calculates the in-situ temperature at which seawater freezes.

       Args:
        planetOceanComp(string): seawater composition, pass in as Planet.Ocean.comp
        SA (MxN array-like, float): absolute salinity [g/kg]
        T (MxNarray-like, float) : in situ temperature in degrees celsius [o^C]
        saturation_fraction: fraction between 0 and 1.  The saturation fraction of
                          dissolved air in seawater.  Default is 0 or
                          completely saturated.
        Returns:
        t_freezing(array-like):  in-situ temperature at which seawater freezes in K

        Examples:
    """

    if planetOceanComp == "Seawater":
        return  T0 + gsw.t_freezing(SA,10*(P-P0/1e5))
    elif planetOceanComp == "MgSO4":
        pass

def GetModgswPressureFreezingCT(planetOceanComp, SA, T):
    """ NEED BETTER DESCRIPTION HERE

       Args:
        planetOceanComp(string): seawater composition, pass in as Planet.Ocean.comp
        SA (MxN array-like, float): absolute salinity [g/kg]
        P(1x1, Mx1, 1xN, MxN array-like, float) : sea pressure in absolute bars [dbar]

        Returns:
        P_f_CT(array-like):  Pressure ?????? [dbar]

        Examples:
    """

    if planetOceanComp == "Seawater":
        import sympy
        return 0.1*(sympy.solve(T - T0 - t_freezing(SA,P))) + P0/1e5;
    elif planetOceanComp == "MgSO4":
        pass


def GetModgswSoundSpeedTExact(planetOceanComp, SA, T, P):
    """ Calculates the speed of sound (c) in seawater.
    The speed of sound in seawater :math:`c` is given by:
    .. math::
        c(SA, t, p) = \sqrt{ \partial P  / \partial \rho |_{SA,\eta}} =
                      \sqrt{(\rho\kappa)^{-1}} =
                      g_P \sqrt{g_{TT}/(g^2_{TP} - g_{TT}g_{PP})}
    Note that in these expressions, since sound speed is in m s :sup`-1` and
    density has units of kg m :sup:`-3` it follows that the pressure of the
    partial derivatives must be in Pa and the isentropic compressibility
    :math:`kappa` must have units of Pa :sup:`-1`. The sound speed c produced
    by both the SIA and the GSW software libraries (appendices M and N) has
    units of m s :sup:`-1`.

       Args:
        planetOceanComp(string): seawater composition, pass in as Planet.Ocean.comp
        SA (MxN array-like, float): absolute salinity [g/kg]
        T (MxNarray-like, float) : in situ temperature in degrees celsius [o^C]
        P(1x1, Mx1, 1xN, MxN array-like, float) : sea pressure in absolute bars [dbar]

        Returns:
        c(array-like): speed of sound in seawater [m/s]

        Examples:
    """
    if planetOceanComp == "Seawater":
        return gsw.sound_speed_t_exact(SA, t, p)
    elif planetOceanComp == "MgSO4":
        pass

def GetModgswSoundSpeedIce(planetOceanComp, T, P):
    """ Calculates the speed of sound (c) in ice

       Args:
        planetOceanComp(string): seawater composition, pass in as Planet.Ocean.comp
        T (MxNarray-like, float) : in situ temperature in degrees celsius [o^C]
        P(1x1, Mx1, 1xN, MxN array-like, float) : sea pressure in absolute bars [dbar]

        Returns:
        c(array-like): speed of sound in ice [m/s]

        Examples:
    """
    '''if planetOceanComp == "Seawater":
        return T - T0+ 10*(P-gsw_P0/1e5)
        #NOTE: THIS SEEMS TO BE IMPLEMENTED INCORRECTLY IN MATLAB AND WE NEED TO DOUBLE CHECK THE EQUATION ABOVE
    elif planetOceanComp == "MgSO4":
        pass'''


def GetVelSIceIh(planetOceanComp, T, P):
    """ Calculates the speed of sound (c) in iceIh

       Args:
        planetOceanComp(string): seawater composition, pass in as Planet.Ocean.comp
        T (MxNarray-like, float) : in situ temperature in degrees celsius [o^C]
        P(1x1, Mx1, 1xN, MxN array-like, float) : sea pressure in absolute bars [dbar]

        Returns:
        c(array-like): speed of sound in iceIh [m/s]

        Examples:
    """
    if planetOceanComp == 'Seawater':
        pass
'''import math
        vP = sound_speed_ice(T-T0, 10*(P-P0/1e5)) #[m/s]
        gi_tt = gsw_gibbs_ice(2,0,T,P);
        gi_tp = gsw_gibbs_ice(1,1,T,P);

kappa_ice = (gi_tp.*gi_tp - gi_tt.*gsw_gibbs_ice(0,2,t,p))./ ...
                  (gsw_gibbs_ice(0,1,t,p).*gi_tt);

        KS = kappa_ice(T-T0,10*(P-P0/1e5)) #[1/Pa]
        V = specvol_ice(T-T0,10*(P-P0/1e5))# [m^3/kg]
        return math.sqrt(0.75*(vP.^2-V./KS)) # [m/s]"""
'''
    #if planetOceanComp == "MgSO4":
        #pass


def SetupEOS(oceanComp):
    print('SetupEOS not implemented yet')