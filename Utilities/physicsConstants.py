#Physics Constants



c = 2.99792458e8 #[m/s], speed of light
G = 6.67408e-11 # [m^3/(kg*s^2)], Gravitational constant
g = 9.80665  # [m/s^2], gravitational acceleration (at Earth's surface)
e = 1.6021766208e-19 #[C], The charge of an electron
h = 6.626070040e-34 #[J*s] Planck's constant
hb = h/(2*pi) #[J*s], H-bar
uB = 5.7883818012e-11 #[MeV/T] Bohr Magneton
E0 = 8.854187817e-12 #[F/m] The permittivity of free space, E0
u0 = 12.566370614e-7 # [N/A^2] The permeability of free space, u0 exact
me = 9.10938356e-31 #[kg], mass of the electron
mp = 1.672621898e-27 #[kg], mass of proton
NA = 6.022140857e23 #[1/mol], Avogadro's number
kB = 1.38064852e-23 # [J/K], Boltzmann's constant
eV = 1.6021766208e-19 #[J], Electron-Volt Joule conversion
k = 1/(4*pi*E0) #[N*m^2/C^2], Coulomb electrostatic constant
a_0 = 4*pi*E0*hb^2/(me*e^2) #[m], Bohr Radius



#T0 and P0 as defined by paper below
"""
AUTHOR:
Trevor McDougall and Paul Barker                    [ help@teos-10.org ]

VERSION NUMBER: 3.05 (27th January 2015)

REFERENCES:
IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
seawater - 2010: Calculation and use of thermodynamic properties.
Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
See section 2.2, appendix A.2 and Table D.1 of this TEOS-10 Manual.

The software is available from http://www.TEOS-10.org
"""
T0 = 273.15 # gsw_T0     Celcius zero point
#The Celcius zero point; 273.15 K.  That is T = t + T0 where T is the
#Absolute Temperature (in degrees K) and t is temperature in degrees C.'''
P0 = 101325
#Returns the  Absolute Pressure of one standard atmosphere in Pa, 1atm = 101325 Pa
