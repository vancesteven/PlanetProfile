import numpy as np
from seafreeze import seafreeze
from getK_Andersson2005 import getK_Andersson2005
from Helgerud_sI import Helgerud_sI
import warnings

def ConvectionDeschampsSotin2001( Ttop , Tbot , Pmid , h , g , waterPhase ):
    """
        Determines solid state convection for ices.

        Parameters
        ----------
        Ttop : float
            surface temperature [K]
        Tbot : float
            bottom temperature [K]
        Pmid : float
            half the pressure at bottom of ice (presumably pressure at middle) [MPa]
        h : float
            thickness of ice layer [m]
        g : float
            gravity [m/s^2]
        waterPhase : int
            index to use
            1 -> Ocean
            2 -> Pure water and Ice I
            ... 6 -> Pure water and Ice VI
            16 -> Clathrates

        Returns
        -------
        Q :
            heat flux at bottom of ice layer [W/m^2]
        deltaTBL
            thickness of (lower) thermal bounder layer [m]
        eTBL
            thickness of the thermal lithosphere [m]
        Tc
            temperature of well mixed interior [K]
        rhoIce
            density of ice [kg/m^3]
        alphaIce
            coefficient of thermal expansion [1/K]
        CpIce
            specific heat [J/kg/K]
        kIce
            thermal conductivity [W/m/K]
        nu
            viscocity [Pa*s]
        CONVECTION_FLAG : boolean


        Notes
        -----
        Based on Deschamps and Sotin 2001,
        Thermal Convection in Icy Satellites, J. Geophys. Res. 106(E3):5107-5121 (DS2001)
    """

    if waterPhase==1:
        raise ValueError("Solid state convection isn't computed for liquid water (ind=1)")
    
    waterPhase -= 1 # python has 0-indexing, so index must be lowered to keep consistent

    Rg = 8.314 # ideal gas constant J/mol/K

    # defining lists with index respective to relevant value
    varstrs = ['water','Ih','II','III','V','VI'] # name
    # DS2001 value for Ih ; mean values for II and III ; high T value for VI
    E = [0, 60e3, np.mean([98, 55])*1e3, np.mean([103, 151])*1e3, 136e3, 110e3] # energy [Joules/mol]
    nu0 = [0, 1e14, 1e18, 5e12, 5e14, 5e14] # viscocity at melting point [Pa*s]
    Dcond = [0, 632, 418, 242, 328, 183] # ???
    
    # parameters relevant to Eq. 9 in DS2001
    # found by numerical experiments in DS2000
    c1 = 1.43
    c2 = -0.03
    
    DeltaT = Tbot - Ttop # difference in temperature across layer [K]
    #Tlith = Ttop + 0.3*DeltaT # lithosphere isotherm, not used in code

    if waterPhase < 29: # not clathrates ??? seems inconsistent with header definitions of indices
        B = E[waterPhase] / 2 / Rg / c1 # [K]
        C = c2*DeltaT # [K]
        Tc = B*( np.sqrt( 1 + 2/B*(Tbot-C) ) - 1 ) # [K] DS2001 Eq. 18
        A = E[waterPhase] / Rg / Tbot # [dimensionless]
        nu = nu0[waterPhase] * np.exp(A*(Tbot/Tc-1)) # DS2001 Eq. 11

        PT = np.array([[Pmid],[Tc]])
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore",category=np.VisibleDeprecationWarning) # seafreeze uses Numpy in a deprecated manner
            SF = seafreeze( PT , varstrs[waterPhase])

        rhoIce = SF.rho[0][0] # [kg/m^3] density of ice
        alphaIce = SF.alpha[0][0] # [1/K] coefficient of thermal expansion
        CpIce = SF.Cp[0][0] # [J/kg/K] specific heat

        kIce = getK_Andersson2005(Pmid,Tc,varstrs[waterPhase],'T')

        # how was this critical value selected
        Ra_crit = 1e5 # critical Rayleigh number (convection may happen if actual Rayleigh number is above this value)
    else:
        #E_clath = 90000
        #nu0_clath = (1e14)*20
        kIce = 0.5 # [W/m/K]
        #Vcm = 19

        B = E / 2 / Rg / c1 # [K]
        C = c2*DeltaT # [K]
        Tc = B*( np.sqrt( 1 + 2/B*(Tbot-C) ) - 1 ) # [K] DS2001 Eq. 18
        A = E / Rg / Tbot # [dimensionless]
        nu = nu0 * np.exp(A*(Tbot/Tc-1)) # DS2001 Eq. 11

        SF = Helgerud_sI(Pmid,Tc)

        rhoIce = SF["rho"] # [kg/m^3] density of ice
        alphaIce = SF["alpha"] # [1/K] coefficient of thermal expansion
        CpIce = 3.19*Tc + 2150 # [J/kg/K] specific heat (Eq.13 Ning et al 2015)
        #Q_crit = kIce*(Tc-Ttop)/h # critical heat ???
        #TBL_crit = kIce*(Tbot-Tc)/Q_crit # [m]
    
    Kappa = kIce / rhoIce / CpIce # [m^2 / s ???]

    # Rayleigh number
    Ra = alphaIce*rhoIce*g*DeltaT*(h**3)/Kappa/nu # DS2001 Eq.4 (Kalosuova 2017 uses nu0 instead of nu)

    if Ra > Ra_crit:
        CONVECTION_FLAG = True
        Ra_del = 0.28*Ra**0.21 # DS2001 Eq. 8

        deltaTBL = (nu*Kappa/alphaIce/rhoIce/g/(Tbot-Tc)*Ra_del)**(1/3) # DS2001 Eq. 19
        Q = kIce*(Tbot-Tc)/deltaTBL # DS2001 Eq. 20
        eTBL = kIce*(Tc-Ttop)/Q # DS2001 Eq. 21

        if eTBL > h:
            CONVECTION_FLAG = False
            Tc = DeltaT
            deltaTBL = 0
            eTBL = 0
            if waterPhase < 7:
                Q = Dcond[waterPhase]*np.log(Tbot/Ttop)/h
            else:
                Q = kIce*DeltaT/h
    else:
        CONVECTION_FLAG = False
        Tc = DeltaT
        deltaTBL = 0
        eTBL = 0
        if waterPhase < 7:
            Q = Dcond[waterPhase]*np.log(Tbot/Ttop)/h
        else:
            Q = kIce*DeltaT/h

    return Q,deltaTBL,eTBL,Tc,rhoIce,alphaIce,CpIce,kIce,nu,CONVECTION_FLAG