r""" Functions for calculating magnetic fields
    in planetocentric Cartesian coordinates from spherical magnetic moments.
    Developed in Python 3.8 for "A perturbation method for evaluating the
    magnetic field induced from an arbitrary, asymmetric ocean world 
    analytically" by Styczinski et al.
    DOI: 10.1016/j.icarus.2021.114840
Author: M. J. Styczinski, mjstyczi@uw.edu
"""

import numpy as np
from math import *

"""
eval_Be()
    Evaluates the excitation field at the given coordinates due to a particular magnetic moment Benm.
    Usage: `Bx`, `By`, `Bz` = eval_Be(`n`, `m`, `Benm`, `x`, `y`, `z`, `r`, `omega=None`, `t=None`)
    Returns:
        Bx, By, Bz (each): complex, ndarray shape(Nvals). A linear array of field values due to these particular n,m values,
            at each of the specified points. Returns a time sequence if omega and t are passed.
    Parameters:
        n: integer. Degree of magnetic moment to be evaluated.
        m: integer. Order of magnetic moment to be evaluated.
        Benm: complex. Excitation moment of degree and order n,m. Units match the output field.
        x,y,z,r: float, shape(Nvals). Linear arrays of corresponding x,y, and z values. r is not needed
            but requiring it as an argument makes the function call identical to eval_Bi. If omega and t
            are passed, these quantities are the trajectory locations.
        omega: float (None). Optional oscillation frequency in rads/s for evaluating time series. Requires t to be passed as well.
        t: float, shape(Nvals) (None). Optional time values in TDB seconds since J2000 epoch. Required if omega is passed.
    """
def eval_Be(n,m,Benm, x,y,z,r, omega=None, t=None):

    if omega is None:
        timeRot = 1.0
    else:
        timeRot = np.exp(-1j * omega * t)

    if n == 1:
        # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   Uniform field components
        # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A1 = sqrt(3/8/np.pi)
        B_base = A1*Benm

        if m == -1:
            Bx = -B_base
            By = 1j*B_base
            Bz = 0
        elif m == 0:
            Bx = 0
            By = 0
            Bz = -sqrt(2)*B_base
        elif m == 1:
            Bx = B_base
            By = 1j*B_base
            Bz = 0
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Be, n={n} and m is not between -n and n.")

    elif n==2:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   Linear field components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A2 = sqrt(15/32/np.pi)
        B_base = 2*A2*Benm

        if m==-2:
            Bx = B_base * (-x + 1j*y)
            By = B_base * (y + 1j*x)
            Bz = 0
        elif m==-1:
            Bx = B_base * -z
            By = B_base * 1j*z
            Bz = B_base * (-x + 1j*y)
        elif m==0:
            rt2o3 = sqrt(2/3)
            Bx = B_base * rt2o3*x
            By = B_base * rt2o3*y
            Bz = B_base * -2*rt2o3*z
        elif m==1:
            Bx = B_base * z
            By = B_base * 1j*z
            Bz = B_base * (x + 1j*y)
        elif m==2:
            Bx = B_base * (-x - 1j*y)
            By = B_base * (y - 1j*x)
            Bz = 0
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Be, n={n} and m is not between -n and n.")

    else:
        print(" n = ", n)
        raise ValueError("In field_xyz.eval_Be, n>2 but only n=1 to n=2 are supported.")

    Bx = Bx * timeRot
    By = By * timeRot
    Bz = Bz * timeRot

    return Bx, By, Bz
    
#############################################

"""
eval_Be_Schmidt()
    Evaluates the excitation field at the given coordinates due to a particular Schmidt semi-normalized
    magnetic moment identified by Gauss coefficients Gnm, Hnm.
    Usage: `Bx`, `By`, `Bz` = eval_Be(`n`, `m`, `Gnm`, `Hnm`, `x`, `y`, `z`, `r`, `omega=None`, `t=None`)
    Returns:
        Bx, By, Bz (each): complex, ndarray shape(Nvals). A linear array of field values due to these particular n,m values,
            at each of the specified points. Returns a time sequence if omega and t are passed.
    Parameters:
        n: integer. Degree of magnetic moment to be evaluated.
        m: integer. Order of magnetic moment to be evaluated.
        Gnm, Hnm: complex. Gauss coefficient of degree and order n,m. Units match the output field.
        x,y,z,r: float, shape(Nvals). Linear arrays of corresponding x,y, and z values. r is not needed
            but requiring it as an argument makes the function call identical to eval_Bi_Schmidt. If omega and t
            are passed, these quantities are the trajectory locations.
        omega: float (None). Optional oscillation frequency in rads/s for evaluating time series. Requires t to be passed as well.
        t: float, shape(Nvals) (None). Optional time values in TDB seconds since J2000 epoch. Required if omega is passed.
    """
def eval_Be_Schmidt(n,m,Gnm,Hnm, x,y,z,r, omega=None, t=None):

    if omega is None:
        timeRot = 1.0
    else:
        timeRot = np.exp(-1j * omega * t)

    if n == 1:
        # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   Uniform field components
        # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        if m == 0:
            Bx = 0j
            By = 0j
            Bz = -Gnm
        elif m == 1:
            Bx = -Gnm
            By = -Hnm
            Bz = 0j
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Be_Schmidt, n={n} and m is not between 0 and n.")

    elif n==2:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   Linear field components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        rt3 = sqrt(3)
        A2 = rt3/2
        Be = -2*A2

        if m==0:
            Bx = Be * Gnm*-x/rt3
            By = Be * Gnm*-y/rt3
            Bz = Be * Gnm*2*z/rt3
        elif m==1:
            Bx = Be * Gnm*z
            By = Be * Hnm*z
            Bz = Be * (Gnm*x + Hnm*y)
        elif m==2:
            Bx = Be * (Gnm*x + Hnm*y)
            By = Be * (Gnm*-y + Hnm*x)
            Bz = 0j
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Be_Schmidt, n={n} and m is not between 0 and n.")

    else:
        print(" n = ", n)
        raise ValueError("In field_xyz.eval_Be_Schmidt, n>2 but only n=1 to n=2 are supported.")

    Bx = Bx * timeRot
    By = By * timeRot
    Bz = Bz * timeRot

    return Bx, By, Bz

#############################################

"""
eval_Bi()
    Evaluates the induced magnetic field at the given coordinates due to a particular magnetic moment Binm.
    Usage: `Bx`, `By`, `Bz` = eval_Bi(`n`, `m`, `Binm`, `x`, `y`, `z`, `r`, `omega=None`, `t=None`)
    Returns:
        Bx, By, Bz (each): complex, ndarray shape(Nvals). A linear array of field values due to these particular n,m values,
            at each of the specified points. Returns a time sequence if omega and t are passed.
    Parameters:
        n: integer. Degree of magnetic moment to be evaluated.
        m: integer. Order of magnetic moment to be evaluated.
        Binm: complex. Magnetic moment of degree and order n,m. Units match the output field.
        x,y,z,r: float, shape(Nvals). Linear arrays of corresponding x,y, and z values. r is redundant but
            saves on computation time to avoid recalculating on every call to this function. If omega and t
            are passed, these quantities are the trajectory locations.
        omega: float (None). Optional oscillation frequency in rads/s for evaluating time series. Requires t to be passed as well.
        t: float, shape(Nvals) (None). Optional time values in TDB seconds since J2000 epoch. Required if omega is passed.
    """
def eval_Bi(n,m,Binm, x,y,z,r, omega=None, t=None):

    if omega is None:
        timeRot = 1.0
    else:
        timeRot = np.exp(-1j*omega*t)

    if n==1:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   Dipole field components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A1 = sqrt(3/8/np.pi)
        Bx = A1*Binm/r**5
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==-1:
            Bx = Bx * ( (2*x**2 - y**2 - z**2) + 1j*(-3*x*y) )
            By = By * ( (3*x*y) + 1j*(x**2 - 2*y**2 + z**2) )
            Bz = Bz * ( (3*x*z) + 1j*(-3*y*z) )
        elif m==0:
            rt2 = sqrt(2)
            Bx = Bx * 3*rt2*(x*z)
            By = By * 3*rt2*(y*z)
            Bz = Bz * -rt2*(x**2 + y**2 - 2*z**2)
        elif m==1:
            Bx = Bx * ( (-2*x**2 + y**2 + z**2) + 1j*(-3*x*y) )
            By = By * ( (-3*x*y) + 1j*(x**2 - 2*y**2 + z**2) )
            Bz = Bz * ( (-3*x*z) + 1j*(-3*y*z) )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi, n={n} and m is not between -n and n.")

    elif n==2:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   Quadrupole moment components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A2 = sqrt(15/32/np.pi)
        Bx = A2*Binm/r**7
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==-2:
            Bx = Bx * ( x*(3*x**2 - 7*y**2 - 2*z**2) + 1j*2*y*(-4*x**2 + y**2 + z**2) )
            By = By * ( y*(7*x**2 - 3*y**2 + 2*z**2) + 1j*2*x*(x**2 - 4*y**2 + z**2) )
            Bz = Bz * ( 5*(x**2 - y**2)*z + 1j*(-10*x*y*z) )
        elif m==-1:
            Bx = Bx * ( (8*x**2*z - 2*z*(y**2 + z**2)) + 1j*(-10*x*y*z) )
            By = By * ( (10*x*y*z) + 1j*2*z*(x**2 - 4*y**2 + z**2) )
            Bz = Bz * ( -2*x*(x**2 + y**2 - 4*z**2) + 1j*2*y*(x**2 + y**2 - 4*z**2) )
        elif m==0:
            rt6 = sqrt(6)
            Bx = Bx * -rt6*x*(x**2 + y**2 - 4*z**2)
            By = By * -rt6*y*(x**2 + y**2 - 4*z**2)
            Bz = Bz * rt6*z*(-3*x**2 - 3*y**2 + 2*z**2)
        elif m==1:
            Bx = Bx * ( 2*z*(-4*x**2 + y**2 + z**2) + 1j*(-10*x*y*z) )
            By = By * ( (-10*x*y*z) + 1j*2*z*(x**2 - 4*y**2 + z**2) )
            Bz = Bz * ( 2*x*(x**2 + y**2 - 4*z**2) + 1j*2*y*(x**2 + y**2 - 4*z**2) )
        elif m==2:
            Bx = Bx * ( x*(3*x**2 - 7*y**2 - 2*z**2) + 1j*(8*x**2*y - 2*y*(y**2 + z**2)) )
            By = By * ( y*(7*x**2 - 3*y**2 + 2*z**2) + 1j*-2*x*(x**2 - 4*y**2 + z**2) )
            Bz = Bz * ( 5*(x**2 - y**2)*z + 1j*(10*x*y*z) )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi, n={n} and m is not between -n and n.")

    elif n==3:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   Octupole moment components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A3 = sqrt(7/64/np.pi)
        Bx = A3*Binm/r**9
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==-3:
            rt5 = sqrt(5)
            Bx = Bx * rt5*( (4*x**4 + 3*y**2*(y**2 + z**2) - 3*x**2*(7*y**2 + z**2)) + 1j*x*y*(-15*x**2 + 13*y**2 + 6*z**2) )
            By = By * rt5*( x*y*(13*x**2 - 15*y**2 + 6*z**2) + 1j*(3*x**4 + 4*y**4 - 3*y**2*z**2 + 3*x**2*(-7*y**2 + z**2)) )
            Bz = Bz * 7*rt5*( x*(x**2 - 3*y**2)*z + 1j*y*(-3*x**2 + y**2)*z )
        elif m==-2:
            rt30 = sqrt(30)
            Bx = Bx * rt30*( x*z*(5*x**2 - 9*y**2 - 2*z**2) + 1j*2*y*z*(-6*x**2 + y**2 + z**2) )
            By = By * rt30*( y*z*(9*x**2 - 5*y**2 + 2*z**2) + 1j*2*x*z*(x**2 - 6*y**2 + z**2) )
            Bz = Bz * rt30*( -(x**2 - y**2)*(x**2 + y**2 - 6*z**2) + 1j*2*x*y*(x**2 + y**2 - 6*z**2) )
        elif m==-1:
            rt3 = sqrt(3)
            Bx = Bx * rt3*( (-4*x**4 + y**4 - 3*y**2*z**2 - 4*z**4 - 3*x**2*(y**2 - 9*z**2)) + 1j*5*x*y*(x**2 + y**2 - 6*z**2) )
            By = By * rt3*( -5*x*y*(x**2 + y**2 - 6*z**2) + 1j*(-x**4 + 4*y**4 - 27*y**2*z**2 + 4*z**4 + 3*x**2*(y**2 + z**2)) )
            Bz = Bz * rt3*( -5*x*z*(3*x**2 + 3*y**2 - 4*z**2) + 1j*5*y*z*(3*x**2 + 3*y**2 - 4*z**2) )
        elif m==0:
            Bx = Bx * -10*x*z*(3*x**2 + 3*y**2 - 4*z**2)
            By = By * -10*y*z*(3*x**2 + 3*y**2 - 4*z**2)
            Bz = Bz * 2*(3*x**4 + 3*y**4 - 24*y**2*z**2 + 8*z**4 + 6*x**2*(y**2 - 4*z**2))
        elif m==1:
            rt3 = sqrt(3)
            Bx = Bx * rt3*( (4*x**4 - y**4 + 3*y**2*z**2 + 4*z**4 + 3*x**2*(y**2 - 9*z**2)) + 1j*5*x*y*(x**2 + y**2 - 6*z**2) )
            By = By * rt3*( 5*x*y*(x**2 + y**2 - 6*z**2) + 1j*(-x**4 + 4*y**4 - 27*y**2*z**2 + 4*z**4 + 3*x**2*(y**2 + z**2)) )
            Bz = Bz * 5*rt3*( x*z*(3*x**2 + 3*y**2 - 4*z**2) + 1j*y*z*(3*x**2 + 3*y**2 - 4*z**2) )
        elif m==2:
            rt30 = sqrt(30)
            Bx = Bx * rt30*( x*z*(5*x**2 - 9*y**2 - 2*z**2) + 1j*-2*y*z*(-6*x**2 + y**2 + z**2) )
            By = By * rt30*( y*z*(9*x**2 - 5*y**2 + 2*z**2) + 1j*-2*x*z*(x**2 - 6*y**2 + z**2) )
            Bz = Bz * rt30*( -(x**2 - y**2)*(x**2 + y**2 - 6*z**2) + 1j*-2*x*y*(x**2 + y**2 - 6*z**2) )
        elif m==3:
            rt5 = sqrt(5)
            Bx = Bx * rt5*( (-4*x**4 - 3*y**2*(y**2 + z**2) + 3*x**2*(7*y**2 + z**2)) + 1j*x*y*(-15*x**2 + 13*y**2 + 6*z**2) )
            By = By * rt5*( x*y*(-13*x**2 + 15*y**2 - 6*z**2) + 1j*(3*x**4 + 4*y**4 - 3*y**2*z**2 + 3*x**2*(-7*y**2 + z**2)) )
            Bz = Bz * 7*rt5*( -x*(x**2 - 3*y**2)*z + 1j*y*(-3*x**2 + y**2)*z )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi, n={n} and m is not between -n and n.")

    elif n==4:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   Hexadecapole moment components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A4 = sqrt(1/192/np.pi)
        Bx = A4*Binm/r**11
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==-4:
            rt2 = sqrt(2)
            rt105 = sqrt(105)
            Bx = Bx * 3*rt105*( 1/2/rt2*x*(5*x**4 - 2*x**2*(23*y**2 + 2*z**2) + 3*y**2*(7*y**2 + 4*z**2)) + 1j*-rt2*y*(6*x**4 + y**2*(y**2 + z**2) - x**2*(11*y**2 + 3*z**2)) )
            By = By * 3*rt105*( 1/2/rt2*y*(21*x**4 + 5*y**4 - 4*y**2*z**2 + x**2*(-46*y**2 + 12*z**2)) + 1j*rt2*x*(x**4 + 6*y**4 - 3*y**2*z**2 + x**2*(-11*y**2 + z**2)) )
            Bz = Bz * 27*rt105*( 1/2/rt2*(x**4 - 6*x**2*y**2 + y**4)*z + 1j*rt2*x*y*(-x**2 + y**2)*z )
        elif m==-3:
            rt105 = sqrt(105)
            Bx = Bx * 9*rt105*( z*(2*x**4 + y**2*(y**2 + z**2) - x**2*(9*y**2 + z**2)) + 1j*x*y*z*(-7*x**2 + 5*y**2 + 2*z**2) )
            By = By * 9*rt105*( x*y*z*(5*x**2 - 7*y**2 + 2*z**2) + 1j*z*(x**4 + 2*y**4 - y**2*z**2 + x**2*(-9*y**2 + z**2)) )
            Bz = Bz * -3*rt105*( x*(x**2 - 3*y**2)*(x**2 + y**2 - 8*z**2) + 1j*y*(-3*x**2 + y**2)*(x**2 + y**2 - 8*z**2) )
        elif m==-2:
            rt2 = sqrt(2)
            rt15 = sqrt(15)
            Bx = Bx * 3*rt15*( -1/rt2*x*(5*x**4 - 9*y**4 + 66*y**2*z**2 + 12*z**4 - 2*x**2*(2*y**2 + 23*z**2)) + 1j*rt2*y*(6*x**4 - y**4 + 5*y**2*z**2 + 6*z**4 + x**2*(5*y**2 - 51*z**2)) )
            By = By * 3*rt15*( 1/rt2*y*(-9*x**4 + 5*y**4 - 46*y**2*z**2 + 12*z**4 + x**2*(-4*y**2 + 66*z**2)) + 1j*-rt2*x*(x**4 - 6*y**4 + 51*y**2*z**2 - 6*z**4 - 5*x**2*(y**2 + z**2)) )
            Bz = Bz * 63*rt15*( -1/rt2*(x**2 - y**2)*z*(x**2 + y**2 - 2*z**2) + 1j*rt2*x*y*z*(x**2 + y**2 - 2*z**2) )
        elif m==-1:
            rt15 = sqrt(15)
            Bx = Bx * 3*rt15*( -z*(18*x**4 - 3*y**4 + y**2*z**2 + 4*z**4 + x**2*(15*y**2 - 41*z**2)) + 1j*21*x*y*z*(x**2 + y**2 - 2*z**2) )
            By = By * 3*rt15*( -21*x*y*z*(x**2 + y**2 - 2*z**2) + 1j*z*(-3*x**4 + 18*y**4 - 41*y**2*z**2 + 4*z**4 + x**2*(15*y**2 + z**2)) )
            Bz = Bz * 9*rt15*( x*(x**4 + y**4 - 12*y**2*z**2 + 8*z**4 + 2*x**2*(y**2 - 6*z**2)) + 1j*-y*(x**4 + y**4 - 12*y**2*z**2 + 8*z**4 + 2*x**2*(y**2 - 6*z**2)) )
        elif m==0:
            rt3 = sqrt(3)
            Bx = Bx * 45/2*rt3*x*(x**4 + y**4 - 12*y**2*z**2 + 8*z**4 + 2*x**2*(y**2 - 6*z**2))
            By = By * 45/2*rt3*y*(x**4 + y**4 - 12*y**2*z**2 + 8*z**4 + 2*x**2*(y**2 - 6*z**2))
            Bz = Bz * 15/2*rt3*z*(15*x**4 + 15*y**4 - 40*y**2*z**2 + 8*z**4 + 10*x**2*(3*y**2 - 4*z**2))
        elif m==1:
            rt15 = sqrt(15)
            Bx = Bx * 3*rt15*( z*(18*x**4 - 3*y**4 + y**2*z**2 + 4*z**4 + x**2*(15*y**2 - 41*z**2)) + 1j*21*x*y*z*(x**2 + y**2 - 2*z**2) )
            By = By * 3*rt15*( 21*x*y*z*(x**2 + y**2 - 2*z**2) + 1j*z*(-3*x**4 + 18*y**4 - 41*y**2*z**2 + 4*z**4 + x**2*(15*y**2 + z**2)) )
            Bz = Bz * -9*rt15*( x*(x**4 + y**4 - 12*y**2*z**2 + 8*z**4 + 2*x**2*(y**2 - 6*z**2)) + 1j*y*(x**4 + y**4 - 12*y**2*z**2 + 8*z**4 + 2*x**2*(y**2 - 6*z**2)) )
        elif m==2:
            rt2 = sqrt(2)
            rt15 = sqrt(15)
            Bx = Bx * 3*rt15*( -1/rt2*x*(5*x**4 - 9*y**4 + 66*y**2*z**2 + 12*z**4 - 2*x**2*(2*y**2 + 23*z**2)) + 1j*rt2*y*(-6*x**4 + y**4 - 5*y**2*z**2 - 6*z**4 + x**2*(-5*y**2 + 51*z**2)) )
            By = By * 3*rt15*( 1/rt2*y*(-9*x**4 + 5*y**4 - 46*y**2*z**2 + 12*z**4 + x**2*(-4*y**2 + 66*z**2)) + 1j*rt2*x*(x**4 - 6*y**4 + 51*y**2*z**2 - 6*z**4 - 5*x**2*(y**2 + z**2)) )
            Bz = Bz * -63*rt15*( 1/rt2*(x**2 - y**2)*z*(x**2 + y**2 - 2*z**2) + 1j*rt2*x*y*z*(x**2 + y**2 - 2*z**2) )
        elif m==3:
            rt105 = sqrt(105)
            Bx = Bx * 9*rt105*( -z*(2*x**4 + y**2*(y**2 + z**2) - x**2*(9*y**2 + z**2)) + 1j*x*y*z*(-7*x**2 + 5*y**2 + 2*z**2) )
            By = By * 9*rt105*( -x*y*z*(5*x**2 - 7*y**2 + 2*z**2) + 1j*z*(x**4 + 2*y**4 - y**2*z**2 + x**2*(-9*y**2 + z**2)) )
            Bz = Bz * 3*rt105*( x*(x**2 - 3*y**2)*(x**2 + y**2 - 8*z**2) + 1j*-y*(-3*x**2 + y**2)*(x**2 + y**2 - 8*z**2) )
        elif m==4:
            rt2 = sqrt(2)
            rt105 = sqrt(105)
            Bx = Bx * 3*rt105*( 1/2/rt2*x*(5*x**4 - 2*x**2*(23*y**2 + 2*z**2) + 3*y**2*(7*y**2 + 4*z**2)) + 1j*rt2*y*(6*x**4 + y**2*(y**2 + z**2) - x**2*(11*y**2 + 3*z**2)) )
            By = By * 3*rt105*( 1/2/rt2*y*(21*x**4 + 5*y**4 - 4*y**2*z**2 + x**2*(-46*y**2 + 12*z**2)) + 1j*-rt2*x*(x**4 + 6*y**4 - 3*y**2*z**2 + x**2*(-11*y**2 + z**2)) )
            Bz = Bz * 27*rt105*( 1/2/rt2*(x**4 - 6*x**2*y**2 + y**4)*z + 1j*rt2*x*y*(x**2 - y**2)*z )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi, n={n} and m is not between -n and n.")

    elif n==5:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   n=5 moment components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A5 = sqrt(11/16/np.pi)
        Bx = A5*Binm/r**13
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==-5:
            rt7 = sqrt(7)
            Bx = Bx * 3/8*rt7*( (6*x**6 - 5*y**4*(y**2 + z**2) - 5*x**4*(17*y**2 + z**2) + 10*x**2*(8*y**4 + 3*y**2*z**2)) + 1j*-x*y*(35*x**4 + 31*y**4 + 20*y**2*z**2 - 10*x**2*(11*y**2 + 2*z**2)) )
            By = By * 3/8*rt7*( x*y*(31*x**4 + 5*y**2*(7*y**2 - 4*z**2) + x**2*(-110*y**2 + 20*z**2)) + 1j*(5*x**6 - 6*y**6 + 5*y**4*z**2 + x**4*(-80*y**2 + 5*z**2) + 5*x**2*(17*y**4 - 6*y**2*z**2)) )
            Bz = Bz * 33/8*rt7*( x*(x**4 - 10*x**2*y**2 + 5*y**4)*z + 1j*-y*(5*x**4 - 10*x**2*y**2 + y**4)*z )
        elif m==-4:
            rt35o2 = sqrt(35/2)
            Bx = Bx * 3*rt35o2*( 1/4*x*z*(7*x**4 + 23*y**4 + 12*y**2*z**2 - 2*x**2*(29*y**2 + 2*z**2)) + 1j*-y*z*(8*x**4 + y**2*(y**2 + z**2) - x**2*(13*y**2 + 3*z**2)) )
            By = By * 3*rt35o2*( 1/4*y*z*(23*x**4 + 7*y**4 - 4*y**2*z**2 + x**2*(-58*y**2 + 12*z**2)) + 1j*x*z*(x**4 + 8*y**4 - 3*y**2*z**2 + x**2*(-13*y**2 + z**2)) )
            Bz = Bz * 3*rt35o2*( -1/4*(x**4 - 6*x**2*y**2 + y**4)*(x**2 + y**2 - 10*z**2) + 1j*x*y*(x**2 - y**2)*(x**2 + y**2 - 10*z**2) )
        elif m==-3:
            rt35 = sqrt(35)
            Bx = Bx * 3/8*rt35*( -(2*x**6 + y**6 - 7*y**4*z**2 - 8*y**2*z**4 - x**4*(7*y**2 + 23*z**2) + x**2*(-8*y**4 + 90*y**2*z**2 + 8*z**4)) + 1j*x*y*(7*x**4 - 5*y**4 + 44*y**2*z**2 + 16*z**4 + 2*x**2*(y**2 - 38*z**2)) )
            By = By * -3/8*rt35*( x*y*(5*x**4 - 7*y**4 + 76*y**2*z**2 - 16*z**4 - 2*x**2*(y**2 + 22*z**2)) + 1j*(x**6 + 2*y**6 - 23*y**4*z**2 + 8*y**2*z**4 - x**4*(8*y**2 + 7*z**2) + x**2*(-7*y**4 + 90*y**2*z**2 - 8*z**4)) )
            Bz = Bz * -9/8*rt35*( x*(x**2 - 3*y**2)*z*(3*x**2 + 3*y**2 - 8*z**2) + 1j*y*(-3*x**2 + y**2)*z*(3*x**2 + 3*y**2 - 8*z**2) )
        elif m==-2:
            rt105o2 = sqrt(105/2)
            Bx = Bx * rt105o2*( 1/2*x*z*(-7*x**4 + 11*y**4 - 26*y**2*z**2 - 4*z**4 + x**2*(4*y**2 + 22*z**2)) + 1j*y*z*(8*x**4 - y**4 + y**2*z**2 + 2*z**4 + x**2*(7*y**2 - 23*z**2)) )
            By = By * rt105o2*( 1/2*y*z*(-11*x**4 + 7*y**4 - 22*y**2*z**2 + 4*z**4 + x**2*(-4*y**2 + 26*z**2)) + 1j*x*z*(-x**4 + 8*y**4 - 23*y**2*z**2 + 2*z**4 + x**2*(7*y**2 + z**2)) )
            Bz = Bz * rt105o2*( 1/2*(x**2 - y**2)*(x**4 + y**4 - 16*y**2*z**2 + 16*z**4 + 2*x**2*(y**2 - 8*z**2)) + 1j*-x*y*(x**4 + y**4 - 16*y**2*z**2 + 16*z**4 + 2*x**2*(y**2 - 8*z**2)) )
        elif m==-1:
            rt15o2 = sqrt(15/2)
            Bx = Bx * 1/4*rt15o2*( (6*x**6 - y**6 + 11*y**4*z**2 + 4*y**2*z**4 - 8*z**6 + x**4*(11*y**2 - 101*z**2) + 2*x**2*(2*y**4 - 45*y**2*z**2 + 58*z**4)) + 1j*-7*x*y*(x**4 + y**4 - 16*y**2*z**2 + 16*z**4 + 2*x**2*(y**2 - 8*z**2)) )
            By = By * 1/4*rt15o2*( 7*x*y*(x**4 + y**4 - 16*y**2*z**2 + 16*z**4 + 2*x**2*(y**2 - 8*z**2)) + 1j*(x**6 - 6*y**6 + 101*y**4*z**2 - 116*y**2*z**4 + 8*z**6 - x**4*(4*y**2 + 11*z**2) + x**2*(-11*y**4 + 90*y**2*z**2 - 4*z**4)) )
            Bz = Bz * 7/4*rt15o2*( x*z*(5*x**4 + 5*y**4 - 20*y**2*z**2 + 8*z**4 + 10*x**2*(y**2 - 2*z**2)) + 1j*-y*z*(5*x**4 + 5*y**4 - 20*y**2*z**2 + 8*z**4 + 10*x**2*(y**2 - 2*z**2)) )
        elif m==0:
            Bx = Bx * 21/4*x*z*(5*x**4 + 5*y**4 - 20*y**2*z**2 + 8*z**4 + 10*x**2*(y**2 - 2*z**2))
            By = By * 21/4*y*z*(5*x**4 + 5*y**4 - 20*y**2*z**2 + 8*z**4 + 10*x**2*(y**2 - 2*z**2))
            Bz = Bz * -3/4*(5*x**6 + 5*y**6 - 90*y**4*z**2 + 120*y**2*z**4 - 16*z**6 + 15*x**4*(y**2 - 6*z**2) + 15*x**2*(y**4 - 12*y**2*z**2 + 8*z**4))
        elif m==1:
            rt15o2 = sqrt(15/2)
            Bx = Bx * 1/4*rt15o2*( (-6*x**6 + y**6 - 11*y**4*z**2 - 4*y**2*z**4 + 8*z**6 + x**4*(-11*y**2 + 101*z**2) - 2*x**2*(2*y**4 - 45*y**2*z**2 + 58*z**4)) + 1j*-7*x*y*(x**4 + y**4 - 16*y**2*z**2 + 16*z**4 + 2*x**2*(y**2 - 8*z**2)) )
            By = By * 1/4*rt15o2*( -7*x*y*(x**4 + y**4 - 16*y**2*z**2 + 16*z**4 + 2*x**2*(y**2 - 8*z**2)) + 1j*(x**6 - 6*y**6 + 101*y**4*z**2 - 116*y**2*z**4 + 8*z**6 - x**4*(4*y**2 + 11*z**2) + x**2*(-11*y**4 + 90*y**2*z**2 - 4*z**4)) )
            Bz = Bz * -7/4*rt15o2*( x*z*(5*x**4 + 5*y**4 - 20*y**2*z**2 + 8*z**4 + 10*x**2*(y**2 - 2*z**2)) + 1j*y*z*(5*x**4 + 5*y**4 - 20*y**2*z**2 + 8*z**4 + 10*x**2*(y**2 - 2*z**2)) )
        elif m==2:
            rt105o2 = sqrt(105/2)
            Bx = Bx * rt105o2*( 1/2*x*z*(-7*x**4 + 11*y**4 - 26*y**2*z**2 - 4*z**4 + x**2*(4*y**2 + 22*z**2)) + 1j*y*z*(-8*x**4 + y**4 - y**2*z**2 - 2*z**4 + x**2*(-7*y**2 + 23*z**2)) )
            By = By * rt105o2*( 1/2*y*z*(-11*x**4 + 7*y**4 - 22*y**2*z**2 + 4*z**4 + x**2*(-4*y**2 + 26*z**2)) + 1j*x*z*(x**4 - 8*y**4 + 23*y**2*z**2 - 2*z**4 - x**2*(7*y**2 + z**2)) )
            Bz = Bz * rt105o2*( 1/2*(x**2 - y**2)*(x**4 + y**4 - 16*y**2*z**2 + 16*z**4 + 2*x**2*(y**2 - 8*z**2)) + 1j*x*y*(x**4 + y**4 - 16*y**2*z**2 + 16*z**4 + 2*x**2*(y**2 - 8*z**2)) )
        elif m==3:
            rt35 = sqrt(35)
            Bx = Bx * 3/8*rt35*( (2*x**6 + y**6 - 7*y**4*z**2 - 8*y**2*z**4 - x**4*(7*y**2 + 23*z**2) + x**2*(-8*y**4 + 90*y**2*z**2 + 8*z**4)) + 1j*x*y*(7*x**4 - 5*y**4 + 44*y**2*z**2 + 16*z**4 + 2*x**2*(y**2 - 38*z**2)) )
            By = By * 3/8*rt35*( x*y*(5*x**4 - 7*y**4 + 76*y**2*z**2 - 16*z**4 - 2*x**2*(y**2 + 22*z**2)) + 1j*-(x**6 + 2*y**6 - 23*y**4*z**2 + 8*y**2*z**4 - x**4*(8*y**2 + 7*z**2) + x**2*(-7*y**4 + 90*y**2*z**2 - 8*z**4)) )
            Bz = Bz * 9/8*rt35*( x*(x**2 - 3*y**2)*z*(3*x**2 + 3*y**2 - 8*z**2) + 1j*-y*(-3*x**2 + y**2)*z*(3*x**2 + 3*y**2 - 8*z**2) )
        elif m==4:
            rt35o2 = sqrt(35/2)
            Bx = Bx * 3*rt35o2*( 1/4*x*z*(7*x**4 + 23*y**4 + 12*y**2*z**2 - 2*x**2*(29*y**2 + 2*z**2)) + 1j*y*z*(8*x**4 + y**2*(y**2 + z**2) - x**2*(13*y**2 + 3*z**2)) )
            By = By * 3*rt35o2*( 1/4*y*z*(23*x**4 + 7*y**4 - 4*y**2*z**2 + x**2*(-58*y**2 + 12*z**2)) + 1j*-x*z*(x**4 + 8*y**4 - 3*y**2*z**2 + x**2*(-13*y**2 + z**2)) )
            Bz = Bz * -3*rt35o2*( 1/4*(x**4 - 6*x**2*y**2 + y**4)*(x**2 + y**2 - 10*z**2) + 1j*x*y*(x**2 - y**2)*(x**2 + y**2 - 10*z**2) )
        elif m==5:
            rt7 = sqrt(7)
            Bx = Bx * -3/8*rt7*( (6*x**6 - 5*y**4*(y**2 + z**2) - 5*x**4*(17*y**2 + z**2) + 10*x**2*(8*y**4 + 3*y**2*z**2)) + 1j*x*y*(35*x**4 + 31*y**4 + 20*y**2*z**2 - 10*x**2*(11*y**2 + 2*z**2)) )
            By = By * 3/8*rt7*( -x*y*(31*x**4 + 5*y**2*(7*y**2 - 4*z**2) + x**2*(-110*y**2 + 20*z**2)) + 1j*(5*x**6 - 6*y**6 + 5*y**4*z**2 + x**4*(-80*y**2 + 5*z**2) + 5*x**2*(17*y**4 - 6*y**2*z**2)) )
            Bz = Bz * -33/8*rt7*( x*(x**4 - 10*x**2*y**2 + 5*y**4)*z + 1j*y*(5*x**4 - 10*x**2*y**2 + y**4)*z )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi, n={n} and m is not between -n and n.")

    elif n==6:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   n=6 moment components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A6 = sqrt(91/2048/np.pi)
        Bx = A6*Binm/r**15
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==-6:
            rt2 = sqrt(2)
            rt33 = sqrt(33)
            Bx = Bx * rt33*( 1/rt2*x*(7*x**6 - 43*y**6 - 30*y**4*z**2 - 3*x**4*(47*y**2 + 2*z**2) + 15*x**2*(15*y**4 + 4*y**2*z**2)) + 1j*rt2*y*(-24*x**6 + 3*y**4*(y**2 + z**2) + 5*x**4*(23*y**2 + 3*z**2) - 6*x**2*(11*y**4 + 5*y**2*z**2)) )
            By = By * rt33*( 1/rt2*y*(43*x**6 - 7*y**6 + 6*y**4*z**2 + x**4*(-225*y**2 + 30*z**2) + 3*x**2*(47*y**4 - 20*y**2*z**2)) + 1j*rt2*x*(3*x**6 - 24*y**6 + 15*y**4*z**2 + x**4*(-66*y**2 + 3*z**2) + 5*x**2*(23*y**4 - 6*y**2*z**2)) )
            Bz = Bz * 13*rt33*( 1/rt2*(x**6 - 15*x**4*y**2 + 15*x**2*y**4 - y**6)*z + 1j*-rt2*x*y*(3*x**4 - 10*x**2*y**2 + 3*y**4)*z )
        elif m==-5:
            rt22 = sqrt(22)
            Bx = Bx * 3*rt22*( z*(8*x**6 - 5*y**4*(y**2 + z**2) + 30*x**2*y**2*(3*y**2 + z**2) - 5*x**4*(21*y**2 + z**2)) + 1j*-x*y*z*(45*x**4 + 33*y**4 + 20*y**2*z**2 - 10*x**2*(13*y**2 + 2*z**2)) )
            By = By * 3*rt22*( x*y*z*(33*x**4 + 5*y**2*(9*y**2 - 4*z**2) + x**2*(-130*y**2 + 20*z**2)) + 1j*z*(5*x**6 - 8*y**6 + 5*y**4*z**2 + x**4*(-90*y**2 + 5*z**2) + 15*x**2*(7*y**4 - 2*y**2*z**2)) )
            Bz = Bz * 3*rt22*( -x*(x**4 - 10*x**2*y**2 + 5*y**4)*(x**2 + y**2 - 12*z**2) + 1j*y*(5*x**4 - 10*x**2*y**2 + y**4)*(x**2 + y**2 - 12*z**2) )
        elif m==-4:
            Bx = Bx * 3*( -x*(7*x**6 + 23*y**6 - 240*y**4*z**2 - 120*y**2*z**4 - 3*x**4*(17*y**2 + 32*z**2) + x**2*(-35*y**4 + 720*y**2*z**2 + 40*z**4)) + 1j*4*y*(8*x**6 + y**6 - 9*y**4*z**2 - 10*y**2*z**4 - 5*x**4*(y**2 + 21*z**2) - 6*x**2*(2*y**4 - 25*y**2*z**2 - 5*z**4)) )
            By = By * -3*( y*(23*x**6 + 7*y**6 - 96*y**4*z**2 + 40*y**2*z**4 - 5*x**4*(7*y**2 + 48*z**2) - 3*x**2*(17*y**4 - 240*y**2*z**2 + 40*z**4)) + 1j*4*x*(x**6 + 8*y**6 - 105*y**4*z**2 + 30*y**2*z**4 - 3*x**4*(4*y**2 + 3*z**2) - 5*x**2*(y**4 - 30*y**2*z**2 + 2*z**4)) )
            Bz = Bz * 33*( -(x**4 - 6*x**2*y**2 + y**4)*z*(3*x**2 + 3*y**2 - 10*z**2) + 1j*4*x*y*(x**2 - y**2)*z*(3*x**2 + 3*y**2 - 10*z**2) )
        elif m==-3:
            rt30 = sqrt(30)
            Bx = Bx * rt30*( z*(-24*x**6 - 9*y**6 + 15*y**4*z**2 + 24*y**2*z**4 + x**4*(75*y**2 + 95*z**2) + 6*x**2*(15*y**4 - 55*y**2*z**2 - 4*z**4)) + 1j*x*y*z*(81*x**4 - 51*y**4 + 140*y**2*z**2 + 48*z**4 + 30*x**2*(y**2 - 10*z**2)) )
            By = By * rt30*( x*y*z*(-51*x**4 + 81*y**4 - 300*y**2*z**2 + 48*z**4 + 10*x**2*(3*y**2 + 14*z**2)) + 1j*z*(-9*x**6 - 24*y**6 + 95*y**4*z**2 - 24*y**2*z**4 + 15*x**4*(6*y**2 + z**2) + 3*x**2*(25*y**4 - 110*y**2*z**2 + 8*z**4)) )
            Bz = Bz * rt30*( x*(x**2 - 3*y**2)*(3*x**4 + 3*y**4 - 60*y**2*z**2 + 80*z**4 + 6*x**2*(y**2 - 10*z**2)) + 1j*y*(-3*x**2 + y**2)*(3*x**4 + 3*y**4 - 60*y**2*z**2 + 80*z**4 + 6*x**2*(y**2 - 10*z**2)) )
        elif m==-2:
            rt2 = sqrt(2)
            rt15 = sqrt(15)
            Bx = Bx * rt15*( 1/rt2*x*(7*x**6 - 11*y**6 + 210*y**4*z**2 - 240*y**2*z**4 - 32*z**6 + 3*x**4*(y**2 - 50*z**2) - 15*x**2*(y**4 - 4*y**2*z**2 - 16*z**4)) + 1j*rt2*y*(-8*x**6 + y**6 - 15*y**4*z**2 + 16*z**6 - 15*x**4*(y**2 - 11*z**2) - 6*x**2*(y**4 - 25*y**2*z**2 + 40*z**4)) )
            By = By * rt15*( 1/rt2*y*(11*x**6 - 7*y**6 + 150*y**4*z**2 - 240*y**2*z**4 + 32*z**6 + 15*x**4*(y**2 - 14*z**2) - 3*x**2*(y**4 + 20*y**2*z**2 - 80*z**4)) + 1j*rt2*x*(x**6 - 8*y**6 + 165*y**4*z**2 - 240*y**2*z**4 + 16*z**6 - 3*x**4*(2*y**2 + 5*z**2) - 15*x**2*(y**4 - 10*y**2*z**2)) )
            Bz = Bz * 3*rt15*( 1/rt2*(x**2 - y**2)*z*(15*x**4 + 15*y**4 - 80*y**2*z**2 + 48*z**4 + 10*x**2*(3*y**2 - 8*z**2)) + 1j*-rt2*x*y*z*(15*x**4 + 15*y**4 - 80*y**2*z**2 + 48*z**4 + 10*x**2*(3*y**2 - 8*z**2)) )
        elif m==-1:
            rt3 = sqrt(3)
            Bx = Bx * 2*rt3*( z*(40*x**6 - 5*y**6 + 15*y**4*z**2 + 12*y**2*z**4 - 8*z**6 + 75*x**4*(y**2 - 3*z**2) + 6*x**2*(5*y**4 - 35*y**2*z**2 + 26*z**4)) + 1j*-3*x*y*z*(15*x**4 + 15*y**4 - 80*y**2*z**2 + 48*z**4 + 10*x**2*(3*y**2 - 8*z**2)) )
            By = By * 2*rt3*( 3*x*y*z*(15*x**4 + 15*y**4 - 80*y**2*z**2 + 48*z**4 + 10*x**2*(3*y**2 - 8*z**2)) + 1j*z*(5*x**6 - 40*y**6 + 225*y**4*z**2 - 156*y**2*z**4 + 8*z**6 - 15*x**4*(2*y**2 + z**2) - 3*x**2*(25*y**4 - 70*y**2*z**2 + 4*z**4)) )
            Bz = Bz * 2*rt3*( -x*(5*x**6 + 5*y**6 - 120*y**4*z**2 + 240*y**2*z**4 - 64*z**6 + 15*x**4*(y**2 - 8*z**2) + 15*x**2*(y**4 - 16*y**2*z**2 + 16*z**4)) + 1j*y*(5*x**6 + 5*y**6 - 120*y**4*z**2 + 240*y**2*z**4 - 64*z**6 + 15*x**4*(y**2 - 8*z**2) + 15*x**2*(y**4 - 16*y**2*z**2 + 16*z**4)) )
        elif m==0:
            rt14 = sqrt(14)
            Bx = Bx * -rt14*x*(5*x**6 + 5*y**6 - 120*y**4*z**2 + 240*y**2*z**4 - 64*z**6 + 15*x**4*(y**2 - 8*z**2) + 15*x**2*(y**4 - 16*y**2*z**2 + 16*z**4))
            By = By * -rt14*y*(5*x**6 + 5*y**6 - 120*y**4*z**2 + 240*y**2*z**4 - 64*z**6 + 15*x**4*(y**2 - 8*z**2) + 15*x**2*(y**4 - 16*y**2*z**2 + 16*z**4))
            Bz = Bz * rt14*z*(-35*x**6 - 35*y**6 + 210*y**4*z**2 - 168*y**2*z**4 + 16*z**6 - 105*x**4*(y**2 - 2*z**2) - 21*x**2*(5*y**4 - 20*y**2*z**2 + 8*z**4))
        elif m==1:
            rt3 = sqrt(3)
            Bx = Bx * -2*rt3*( z*(40*x**6 - 5*y**6 + 15*y**4*z**2 + 12*y**2*z**4 - 8*z**6 + 75*x**4*(y**2 - 3*z**2) + 6*x**2*(5*y**4 - 35*y**2*z**2 + 26*z**4)) + 1j*3*x*y*z*(15*x**4 + 15*y**4 - 80*y**2*z**2 + 48*z**4 + 10*x**2*(3*y**2 - 8*z**2)) )
            By = By * 2*rt3*( -3*x*y*z*(15*x**4 + 15*y**4 - 80*y**2*z**2 + 48*z**4 + 10*x**2*(3*y**2 - 8*z**2)) + 1j*z*(5*x**6 - 40*y**6 + 225*y**4*z**2 - 156*y**2*z**4 + 8*z**6 - 15*x**4*(2*y**2 + z**2) - 3*x**2*(25*y**4 - 70*y**2*z**2 + 4*z**4)) )
            Bz = Bz * 2*rt3*( x*(5*x**6 + 5*y**6 - 120*y**4*z**2 + 240*y**2*z**4 - 64*z**6 + 15*x**4*(y**2 - 8*z**2) + 15*x**2*(y**4 - 16*y**2*z**2 + 16*z**4)) + 1j*y*(5*x**6 + 5*y**6 - 120*y**4*z**2 + 240*y**2*z**4 - 64*z**6 + 15*x**4*(y**2 - 8*z**2) + 15*x**2*(y**4 - 16*y**2*z**2 + 16*z**4)) )
        elif m==2:
            rt2 = sqrt(2)
            rt15 = sqrt(15)
            Bx = Bx * rt15*( 1/rt2*x*(7*x**6 - 11*y**6 + 210*y**4*z**2 - 240*y**2*z**4 - 32*z**6 + 3*x**4*(y**2 - 50*z**2) - 15*x**2*(y**4 - 4*y**2*z**2 - 16*z**4)) + 1j*rt2*y*(8*x**6 - y**6 + 15*y**4*z**2 - 16*z**6 + 15*x**4*(y**2 - 11*z**2) + 6*x**2*(y**4 - 25*y**2*z**2 + 40*z**4)) )
            By = By * rt15*( 1/rt2*y*(11*x**6 - 7*y**6 + 150*y**4*z**2 - 240*y**2*z**4 + 32*z**6 + 15*x**4*(y**2 - 14*z**2) - 3*x**2*(y**4 + 20*y**2*z**2 - 80*z**4)) + 1j*rt2*x*(-x**6 + 8*y**6 - 165*y**4*z**2 + 240*y**2*z**4 - 16*z**6 + 3*x**4*(2*y**2 + 5*z**2) + 15*x**2*(y**4 - 10*y**2*z**2)) )
            Bz = Bz * 3*rt15*( 1/rt2*(x**2 - y**2)*z*(15*x**4 + 15*y**4 - 80*y**2*z**2 + 48*z**4 + 10*x**2*(3*y**2 - 8*z**2)) + 1j*rt2*x*y*z*(15*x**4 + 15*y**4 - 80*y**2*z**2 + 48*z**4 + 10*x**2*(3*y**2 - 8*z**2)) )
        elif m==3:
            rt30 = sqrt(30)
            Bx = Bx * rt30*( z*(24*x**6 - 5*x**4*(15*y**2 + 19*z**2) + 3*y**2*(3*y**4 - 5*y**2*z**2 - 8*z**4) + x**2*(-90*y**4 + 330*y**2*z**2 + 24*z**4)) + 1j*x*y*z*(81*x**4 - 51*y**4 + 140*y**2*z**2 + 48*z**4 + 30*x**2*(y**2 - 10*z**2)) )
            By = By * rt30*( x*y*z*(51*x**4 - 81*y**4 + 300*y**2*z**2 - 48*z**4 - 10*x**2*(3*y**2 + 14*z**2)) + 1j*z*(-9*x**6 - 24*y**6 + 95*y**4*z**2 - 24*y**2*z**4 + 15*x**4*(6*y**2 + z**2) + 3*x**2*(25*y**4 - 110*y**2*z**2 + 8*z**4)) )
            Bz = Bz * rt30*( -x*(x**2 - 3*y**2)*(3*x**4 + 3*y**4 - 60*y**2*z**2 + 80*z**4 + 6*x**2*(y**2 - 10*z**2)) + 1j*y*(-3*x**2 + y**2)*(3*x**4 + 3*y**4 - 60*y**2*z**2 + 80*z**4 + 6*x**2*(y**2 - 10*z**2)) )
        elif m==4:
            Bx = Bx * -3*( x*(7*x**6 + 23*y**6 - 240*y**4*z**2 - 120*y**2*z**4 - 3*x**4*(17*y**2 + 32*z**2) + x**2*(-35*y**4 + 720*y**2*z**2 + 40*z**4)) + 1j*4*y*(8*x**6 + y**6 - 9*y**4*z**2 - 10*y**2*z**4 - 5*x**4*(y**2 + 21*z**2) - 6*x**2*(2*y**4 - 25*y**2*z**2 - 5*z**4)) )
            By = By * 3*( -y*(23*x**6 + 7*y**6 - 96*y**4*z**2 + 40*y**2*z**4 - 5*x**4*(7*y**2 + 48*z**2) - 3*x**2*(17*y**4 - 240*y**2*z**2 + 40*z**4)) + 1j*4*x*(x**6 + 8*y**6 - 105*y**4*z**2 + 30*y**2*z**4 - 3*x**4*(4*y**2 + 3*z**2) - 5*x**2*(y**4 - 30*y**2*z**2 + 2*z**4)) )
            Bz = Bz * -33*( (x**4 - 6*x**2*y**2 + y**4)*z*(3*x**2 + 3*y**2 - 10*z**2) + 1j*4*x*y*(x**2 - y**2)*z*(3*x**2 + 3*y**2 - 10*z**2) )
        elif m==5:
            rt22 = sqrt(22)
            Bx = Bx * -3*rt22*( z*(8*x**6 - 5*y**4*(y**2 + z**2) + 30*x**2*y**2*(3*y**2 + z**2) - 5*x**4*(21*y**2 + z**2)) + 1j*x*y*z*(45*x**4 + 33*y**4 + 20*y**2*z**2 - 10*x**2*(13*y**2 + 2*z**2)) )
            By = By * 3*rt22*( -x*y*z*(33*x**4 + 5*y**2*(9*y**2 - 4*z**2) + x**2*(-130*y**2 + 20*z**2)) + 1j*z*(5*x**6 - 8*y**6 + 5*y**4*z**2 + x**4*(-90*y**2 + 5*z**2) + 15*x**2*(7*y**4 - 2*y**2*z**2)) )
            Bz = Bz * 3*rt22*( x*(x**4 - 10*x**2*y**2 + 5*y**4)*(x**2 + y**2 - 12*z**2) + 1j*y*(5*x**4 - 10*x**2*y**2 + y**4)*(x**2 + y**2 - 12*z**2) )
        elif m==6:
            rt2 = sqrt(2)
            rt33 = sqrt(33)
            Bx = Bx * rt33*( 1/rt2*x*(7*x**6 - 43*y**6 - 30*y**4*z**2 - 3*x**4*(47*y**2 + 2*z**2) + 15*x**2*(15*y**4 + 4*y**2*z**2)) + 1j*rt2*y*(24*x**6 - 3*y**4*(y**2 + z**2) - 5*x**4*(23*y**2 + 3*z**2) + 6*x**2*(11*y**4 + 5*y**2*z**2)) )
            By = By * rt33*( 1/rt2*y*(43*x**6 - 7*y**6 + 6*y**4*z**2 + x**4*(-225*y**2 + 30*z**2) + 3*x**2*(47*y**4 - 20*y**2*z**2)) + 1j*rt2*x*(-3*x**6 + 3*y**4*(8*y**2 - 5*z**2) + x**4*(66*y**2 - 3*z**2) - 5*x**2*(23*y**4 - 6*y**2*z**2)) )
            Bz = Bz * 13*rt33*( 1/rt2*(x**6 - 15*x**4*y**2 + 15*x**2*y**4 - y**6)*z + 1j*rt2*x*y*(3*x**4 - 10*x**2*y**2 + 3*y**4)*z )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi, n={n} and m is not between -n and n.")

    elif n==7:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   n=7 moment components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A7 = sqrt(15/16/np.pi)
        Bx = A7*Binm/r**17
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==-7:
            rt429o2 = sqrt(429/2)
            Bx = Bx * 1/16*rt429o2*( (8*x**8 + 7*y**6*(y**2 + z**2) + 105*x**4*y**2*(5*y**2 + z**2) - 7*x**6*(31*y**2 + z**2) - 7*x**2*(29*y**6 + 15*y**4*z**2)) + 1j*x*y*(-63*x**6 + 57*y**6 + 42*y**4*z**2 + 7*x**4*(61*y**2 + 6*z**2) - 7*x**2*(59*y**4 + 20*y**2*z**2)) )
            By = By * 1/16*rt429o2*( x*y*(57*x**6 - 63*y**6 + 42*y**4*z**2 - 7*x**4*(59*y**2 - 6*z**2) + 7*x**2*(61*y**4 - 20*y**2*z**2)) + 1j*(7*x**8 + 8*y**8 - 7*y**6*z**2 + 7*x**6*(-29*y**2 + z**2) + 105*x**4*(5*y**4 - y**2*z**2) - 7*x**2*(31*y**6 - 15*y**4*z**2)) )
            Bz = Bz * 15/16*rt429o2*( x*(x**6 - 21*x**4*y**2 + 35*x**2*y**4 - 7*y**6)*z + 1j*y*(-7*x**6 + 35*x**4*y**2 - 21*x**2*y**4 + y**6)*z )
        elif m==-6:
            rt3003 = sqrt(3003)
            Bx = Bx * 3/8*rt3003*( 1/2*x*z*(3*x**6 - x**4*(57*y**2 + 2*z**2) + 5*x**2*(17*y**4 + 4*y**2*z**2) - 5*(3*y**6 + 2*y**4*z**2)) + 1j*y*z*(-10*x**6 + y**4*(y**2 + z**2) + 5*x**4*(9*y**2 + z**2) - 2*x**2*(12*y**4 + 5*y**2*z**2)) )
            By = By * 3/8*rt3003*( 1/2*y*z*(15*x**6 - 3*y**6 + 2*y**4*z**2 + x**4*(-85*y**2 + 10*z**2) + x**2*(57*y**4 - 20*y**2*z**2)) + 1j*x*z*(x**6 + x**4*(-24*y**2 + z**2) + 5*y**4*(-2*y**2 + z**2) + 5*x**2*(9*y**4 - 2*y**2*z**2)) )
            Bz = Bz * 1/8*rt3003*( -1/2*(x**6 - 15*x**4*y**2 + 15*x**2*y**4 - y**6)*(x**2 + y**2 - 14*z**2) + 1j*x*y*(3*x**4 - 10*x**2*y**2 + 3*y**4)*(x**2 + y**2 - 14*z**2) )
        elif m==-5:
            rt231o2 = sqrt(231/2)
            Bx = Bx * 1/16*rt231o2*( (-8*x**8 + x**6*(97*y**2 + 127*z**2) + 5*y**4*(y**4 - 11*y**2*z**2 - 12*z**4) + 15*x**4*(y**4 - 103*y**2*z**2 - 4*z**4) + x**2*(-85*y**6 + 1185*y**4*z**2 + 360*y**2*z**4)) + 1j*x*y*(45*x**6 + 33*y**6 - 402*y**4*z**2 - 240*y**2*z**4 - 5*x**4*(17*y**2 + 138*z**2) + x**2*(-97*y**4 + 1820*y**2*z**2 + 240*z**4)) )
            By = By * 1/16*rt231o2*( x*y*(-33*x**6 + x**4*(97*y**2 + 402*z**2) - 15*y**2*(3*y**4 - 46*y**2*z**2 + 16*z**4) + 5*x**2*(17*y**4 - 364*y**2*z**2 + 48*z**4)) + 1j*(-5*x**8 + 8*y**8 - 127*y**6*z**2 + 60*y**4*z**4 + x**6*(85*y**2 + 55*z**2) - 15*x**4*(y**4 + 79*y**2*z**2 - 4*z**4) + x**2*(-97*y**6 + 1545*y**4*z**2 - 360*y**2*z**4)) )
            Bz = Bz * 39/16*rt231o2*( -x*(x**4 - 10*x**2*y**2 + 5*y**4)*z*(x**2 + y**2 - 4*z**2) + 1j*y*(5*x**4 - 10*x**2*y**2 + y**4)*z*(x**2 + y**2 - 4*z**2) )
        elif m==-4:
            rt231o2 = sqrt(231/2)
            Bx = Bx * 1/2*rt231o2*( 1/4*x*z*(-27*x**6 + x**4*(183*y**2 + 128*z**2) + 5*x**2*(27*y**4 - 176*y**2*z**2 - 8*z**4) + 15*y**2*(-5*y**4 + 16*y**2*z**2 + 8*z**4)) + 1j*y*z*(30*x**6 + 3*y**6 - 7*y**4*z**2 - 10*y**2*z**4 - 15*x**4*(y**2 + 9*z**2) + x**2*(-42*y**4 + 170*y**2*z**2 + 30*z**4)) )
            By = By * 1/2*rt231o2*( 1/4*y*z*(-75*x**6 - 27*y**6 + 128*y**4*z**2 - 40*y**2*z**4 + 15*x**4*(9*y**2 + 16*z**2) + x**2*(183*y**4 - 880*y**2*z**2 + 120*z**4)) + 1j*x*z*(-3*x**6 + 7*x**4*(6*y**2 + z**2) + 5*x**2*(3*y**4 - 34*y**2*z**2 + 2*z**4) - 15*(2*y**6 - 9*y**4*z**2 + 2*y**2*z**4)) )
            Bz = Bz * 3/2*rt231o2*( 1/4*(x**4 - 6*x**2*y**2 + y**4)*(x**4 + y**4 - 24*y**2*z**2 + 40*z**4 + 2*x**2*(y**2 - 12*z**2)) + 1j*-x*y*(x**2 - y**2)*(x**4 + y**4 - 24*y**2*z**2 + 40*z**4 + 2*x**2*(y**2 - 12*z**2)) )
        elif m==-3:
            rt21o2 = sqrt(21/2)
            Bx = Bx * 3/16*rt21o2*( (8*x**8 + 3*y**8 - 57*y**6*z**2 + 20*y**4*z**4 + 80*y**2*z**6 - x**6*(17*y**2 + 207*z**2) + x**4*(-55*y**4 + 585*y**2*z**2 + 420*z**4) - x**2*(27*y**6 - 735*y**4*z**2 + 1320*y**2*z**4 + 80*z**6)) + 1j*-x*y*(27*x**6 - 17*y**6 + 378*y**4*z**2 - 480*y**2*z**4 - 160*z**6 + x**4*(37*y**2 - 678*z**2) + x**2*(-7*y**4 - 300*y**2*z**2 + 1280*z**4)) )
            By = By * 3/16*rt21o2*( x*y*(17*x**6 - 27*y**6 + 678*y**4*z**2 - 1280*y**2*z**4 + 160*z**6 + 7*x**4*(y**2 - 54*z**2) + x**2*(-37*y**4 + 300*y**2*z**2 + 480*z**4)) + 1j*(3*x**8 + 8*y**8 - 207*y**6*z**2 + 420*y**4*z**4 - 80*y**2*z**6 - 3*x**6*(9*y**2 + 19*z**2) + x**4*(-55*y**4 + 735*y**2*z**2 + 20*z**4) + x**2*(-17*y**6 + 585*y**4*z**2 - 1320*y**2*z**4 + 80*z**6)) )
            Bz = Bz * 55/16*rt21o2*( x*(x**2 - 3*y**2)*z*(3*x**4 + 3*y**4 - 20*y**2*z**2 + 16*z**4 + x**2*(6*y**2 - 20*z**2)) + 1j*y*(-3*x**2 + y**2)*z*(3*x**4 + 3*y**4 - 20*y**2*z**2 + 16*z**4 + x**2*(6*y**2 - 20*z**2)) )
        elif m==-2:
            rt21 = sqrt(21)
            Bx = Bx * 1/8*rt21*( 1/2*x*z*(135*x**6 + x**4*(75*y**2 - 970*z**2) + x**2*(-255*y**4 + 260*y**2*z**2 + 944*z**4) - 3*(65*y**6 - 410*y**4*z**2 + 272*y**2*z**4 + 32*z**6)) + 1j*y*z*(-150*x**6 + 15*y**6 - 65*y**4*z**2 - 32*y**2*z**4 + 48*z**6 - 15*x**4*(19*y**2 - 69*z**2) - 2*x**2*(60*y**4 - 485*y**2*z**2 + 456*z**4)) )
            By = By * 1/8*rt21*( 1/2*y*z*(195*x**6 - 135*y**6 + 970*y**4*z**2 - 944*y**2*z**4 + 96*z**6 + 15*x**4*(17*y**2 - 82*z**2) + x**2*(-75*y**4 - 260*y**2*z**2 + 816*z**4)) + 1j*x*z*(15*x**6 - 150*y**6 + 1035*y**4*z**2 - 912*y**2*z**4 + 48*z**6 - 5*x**4*(24*y**2 + 13*z**2) + x**2*(-285*y**4 + 970*y**2*z**2 - 32*z**4)) )
            Bz = Bz * 15/8*rt21*( -1/2*(x**2 - y**2)*(x**6 + y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6 + 3*x**4*(y**2 - 10*z**2) + x**2*(3*y**4 - 60*y**2*z**2 + 80*z**4)) + 1j*x*y*(x**6 + y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6 + 3*x**4*(y**2 - 10*z**2) + x**2*(3*y**4 - 60*y**2*z**2 + 80*z**4)) )
        elif m==-1:
            rt7o2 = sqrt(7/2)
            Bx = Bx * 1/16*rt7o2*( -(40*x**8 - 5*y**8 + 115*y**6*z**2 - 120*y**4*z**4 - 176*y**2*z**6 + 64*z**8 + 5*x**6*(23*y**2 - 247*z**2) + 15*x**4*(7*y**4 - 157*y**2*z**2 + 232*z**4) + x**2*(25*y**6 - 1005*y**4*z**2 + 3360*y**2*z**4 - 1616*z**6)) + 1j*45*x*y*(x**6 + y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6 + 3*x**4*(y**2 - 10*z**2) + x**2*(3*y**4 - 60*y**2*z**2 + 80*z**4)) )
            By = By * 1/16*rt7o2*( -45*x*y*(x**6 + y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6 + 3*x**4*(y**2 - 10*z**2) + x**2*(3*y**4 - 60*y**2*z**2 + 80*z**4)) + 1j*(-5*x**8 + 40*y**8 - 1235*y**6*z**2 + 3480*y**4*z**4 - 1616*y**2*z**6 + 64*z**8 + 5*x**6*(5*y**2 + 23*z**2) + 15*x**4*(7*y**4 - 67*y**2*z**2 - 8*z**4) + x**2*(115*y**6 - 2355*y**4*z**2 + 3360*y**2*z**4 - 176*z**6)) )
            Bz = Bz * 9/16*rt7o2*( -x*z*(35*x**6 + 35*y**6 - 280*y**4*z**2 + 336*y**2*z**4 - 64*z**6 + 35*x**4*(3*y**2 - 8*z**2) + 7*x**2*(15*y**4 - 80*y**2*z**2 + 48*z**4)) + 1j*y*z*(35*x**6 + 35*y**6 - 280*y**4*z**2 + 336*y**2*z**4 - 64*z**6 + 35*x**4*(3*y**2 - 8*z**2) + 7*x**2*(15*y**4 - 80*y**2*z**2 + 48*z**4)) )
        elif m==0:
            Bx = Bx * -9/8*x*z*(35*x**6 + 35*y**6 - 280*y**4*z**2 + 336*y**2*z**4 - 64*z**6 + 35*x**4*(3*y**2 - 8*z**2) + 7*x**2*(15*y**4 - 80*y**2*z**2 + 48*z**4))
            By = By * -9/8*y*z*(35*x**6 + 35*y**6 - 280*y**4*z**2 + 336*y**2*z**4 - 64*z**6 + 35*x**4*(3*y**2 - 8*z**2) + 7*x**2*(15*y**4 - 80*y**2*z**2 + 48*z**4))
            Bz = Bz * (35/8*(x**8 + y**8) - 140*y**6*z**2 + 420*y**4*z**4 - 224*y**2*z**6 + 16*z**8 + 35/2*x**6*(y**2 - 8*z**2) + 105/4*x**4*(y**4 - 16*y**2*z**2 + 16*z**4) + 7/2*x**2*(5*y**6 - 120*y**4*z**2 + 240*y**2*z**4 - 64*z**6))
        elif m==1:
            rt7o2 = sqrt(7/2)
            Bx = Bx * 1/16*rt7o2*( (40*x**8 - 5*y**8 + 115*y**6*z**2 - 120*y**4*z**4 - 176*y**2*z**6 + 64*z**8 + 5*x**6*(23*y**2 - 247*z**2) + 15*x**4*(7*y**4 - 157*y**2*z**2 + 232*z**4) + x**2*(25*y**6 - 1005*y**4*z**2 + 3360*y**2*z**4 - 1616*z**6)) + 1j*45*x*y*(x**6 + y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6 + 3*x**4*(y**2 - 10*z**2) + x**2*(3*y**4 - 60*y**2*z**2 + 80*z**4)) )
            By = By * 1/16*rt7o2*( 45*x*y*(x**6 + y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6 + 3*x**4*(y**2 - 10*z**2) + x**2*(3*y**4 - 60*y**2*z**2 + 80*z**4)) + 1j*(-5*x**8 + 40*y**8 - 1235*y**6*z**2 + 3480*y**4*z**4 - 1616*y**2*z**6 + 64*z**8 + 5*x**6*(5*y**2 + 23*z**2) + 15*x**4*(7*y**4 - 67*y**2*z**2 - 8*z**4) + x**2*(115*y**6 - 2355*y**4*z**2 + 3360*y**2*z**4 - 176*z**6)) )
            Bz = Bz * 9/16*rt7o2*( x*z*(35*x**6 + 35*y**6 - 280*y**4*z**2 + 336*y**2*z**4 - 64*z**6 + 35*x**4*(3*y**2 - 8*z**2) + 7*x**2*(15*y**4 - 80*y**2*z**2 + 48*z**4)) + 1j*y*z*(35*x**6 + 35*y**6 - 280*y**4*z**2 + 336*y**2*z**4 - 64*z**6 + 35*x**4*(3*y**2 - 8*z**2) + 7*x**2*(15*y**4 - 80*y**2*z**2 + 48*z**4)) )
        elif m==2:
            rt21 = sqrt(21)
            Bx = Bx * 1/8*rt21*( 1/2*x*z*(135*x**6 + x**4*(75*y**2 - 970*z**2) + x**2*(-255*y**4 + 260*y**2*z**2 + 944*z**4) - 3*(65*y**6 - 410*y**4*z**2 + 272*y**2*z**4 + 32*z**6)) + 1j*y*z*(150*x**6 - 15*y**6 + 65*y**4*z**2 + 32*y**2*z**4 - 48*z**6 + 15*x**4*(19*y**2 - 69*z**2) + 2*x**2*(60*y**4 - 485*y**2*z**2 + 456*z**4)) )
            By = By * 1/8*rt21*( 1/2*y*z*(195*x**6 - 135*y**6 + 970*y**4*z**2 - 944*y**2*z**4 + 96*z**6 + 15*x**4*(17*y**2 - 82*z**2) + x**2*(-75*y**4 - 260*y**2*z**2 + 816*z**4)) + 1j*x*z*(-15*x**6 + 5*x**4*(24*y**2 + 13*z**2) + x**2*(285*y**4 - 970*y**2*z**2 + 32*z**4) + 3*(50*y**6 - 345*y**4*z**2 + 304*y**2*z**4 - 16*z**6)) )
            Bz = Bz * -15/8*rt21*( 1/2*(x**2 - y**2)*(x**6 + y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6 + 3*x**4*(y**2 - 10*z**2) + x**2*(3*y**4 - 60*y**2*z**2 + 80*z**4)) + 1j*x*y*(x**6 + y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6 + 3*x**4*(y**2 - 10*z**2) + x**2*(3*y**4 - 60*y**2*z**2 + 80*z**4)) )
        elif m==3:
            rt21o2 = sqrt(21/2)
            Bx = Bx * -3/16*rt21o2*( (8*x**8 + 3*y**8 - 57*y**6*z**2 + 20*y**4*z**4 + 80*y**2*z**6 - x**6*(17*y**2 + 207*z**2) + x**4*(-55*y**4 + 585*y**2*z**2 + 420*z**4) - x**2*(27*y**6 - 735*y**4*z**2 + 1320*y**2*z**4 + 80*z**6)) + 1j*x*y*(27*x**6 - 17*y**6 + 378*y**4*z**2 - 480*y**2*z**4 - 160*z**6 + x**4*(37*y**2 - 678*z**2) + x**2*(-7*y**4 - 300*y**2*z**2 + 1280*z**4)) )
            By = By * 3/16*rt21o2*( -x*y*(17*x**6 - 27*y**6 + 678*y**4*z**2 - 1280*y**2*z**4 + 160*z**6 + 7*x**4*(y**2 - 54*z**2) + x**2*(-37*y**4 + 300*y**2*z**2 + 480*z**4)) + 1j*(3*x**8 + 8*y**8 - 207*y**6*z**2 + 420*y**4*z**4 - 80*y**2*z**6 - 3*x**6*(9*y**2 + 19*z**2) + x**4*(-55*y**4 + 735*y**2*z**2 + 20*z**4) + x**2*(-17*y**6 + 585*y**4*z**2 - 1320*y**2*z**4 + 80*z**6)) )
            Bz = Bz * 55/16*rt21o2*( -x*(x**2 - 3*y**2)*z*(3*x**4 + 3*y**4 - 20*y**2*z**2 + 16*z**4 + x**2*(6*y**2 - 20*z**2)) + 1j*y*(-3*x**2 + y**2)*z*(3*x**4 + 3*y**4 - 20*y**2*z**2 + 16*z**4 + x**2*(6*y**2 - 20*z**2)) )
        elif m==4:
            rt231o2 = sqrt(231/2)
            Bx = Bx * 1/2*rt231o2*( 1/4*x*z*(-27*x**6 + x**4*(183*y**2 + 128*z**2) + 5*x**2*(27*y**4 - 176*y**2*z**2 - 8*z**4) + 15*y**2*(-5*y**4 + 16*y**2*z**2 + 8*z**4)) + 1j*y*z*(-30*x**6 - 3*y**6 + 7*y**4*z**2 + 10*y**2*z**4 + 15*x**4*(y**2 + 9*z**2) + 2*x**2*(21*y**4 - 85*y**2*z**2 - 15*z**4)) )
            By = By * 1/2*rt231o2*( 1/4*y*z*(-75*x**6 - 27*y**6 + 128*y**4*z**2 - 40*y**2*z**4 + 15*x**4*(9*y**2 + 16*z**2) + x**2*(183*y**4 - 880*y**2*z**2 + 120*z**4)) + 1j*x*z*(3*x**6 - 7*x**4*(6*y**2 + z**2) - 5*x**2*(3*y**4 - 34*y**2*z**2 + 2*z**4) + 15*(2*y**6 - 9*y**4*z**2 + 2*y**2*z**4)) )
            Bz = Bz * 3/2*rt231o2*( 1/4*(x**4 - 6*x**2*y**2 + y**4)*(x**4 + y**4 - 24*y**2*z**2 + 40*z**4 + 2*x**2*(y**2 - 12*z**2)) + 1j*x*y*(x**2 - y**2)*(x**4 + y**4 - 24*y**2*z**2 + 40*z**4 + 2*x**2*(y**2 - 12*z**2)) )
        elif m==5:
            rt231o2 = sqrt(231/2)
            Bx = Bx * 1/16*rt231o2*( (8*x**8 - 5*y**8 + 55*y**6*z**2 + 60*y**4*z**4 - x**6*(97*y**2 + 127*z**2) - 15*x**4*(y**4 - 103*y**2*z**2 - 4*z**4) + 5*x**2*(17*y**6 - 237*y**4*z**2 - 72*y**2*z**4)) + 1j*x*y*(45*x**6 + 33*y**6 - 402*y**4*z**2 - 240*y**2*z**4 - 5*x**4*(17*y**2 + 138*z**2) + x**2*(-97*y**4 + 1820*y**2*z**2 + 240*z**4)) )
            By = By * 1/16*rt231o2*( x*y*(33*x**6 - x**4*(97*y**2 + 402*z**2) + 15*y**2*(3*y**4 - 46*y**2*z**2 + 16*z**4) - 5*x**2*(17*y**4 - 364*y**2*z**2 + 48*z**4)) + 1j*(-5*x**8 + 8*y**8 - 127*y**6*z**2 + 60*y**4*z**4 + x**6*(85*y**2 + 55*z**2) - 15*x**4*(y**4 + 79*y**2*z**2 - 4*z**4) + x**2*(-97*y**6 + 1545*y**4*z**2 - 360*y**2*z**4)) )
            Bz = Bz * 39/16*rt231o2*( x*(x**4 - 10*x**2*y**2 + 5*y**4)*z*(x**2 + y**2 - 4*z**2) + 1j*y*(5*x**4 - 10*x**2*y**2 + y**4)*z*(x**2 + y**2 - 4*z**2) )
        elif m==6:
            rt3003 = sqrt(3003)
            Bx = Bx * 3/8*rt3003*( 1/2*x*z*(3*x**6 - x**4*(57*y**2 + 2*z**2) + 5*x**2*(17*y**4 + 4*y**2*z**2) - 5*(3*y**6 + 2*y**4*z**2)) + 1j*-y*z*(-10*x**6 + y**4*(y**2 + z**2) + 5*x**4*(9*y**2 + z**2) - 2*x**2*(12*y**4 + 5*y**2*z**2)) )
            By = By * 3/8*rt3003*( 1/2*y*z*(15*x**6 - 3*y**6 + 2*y**4*z**2 + x**4*(-85*y**2 + 10*z**2) + x**2*(57*y**4 - 20*y**2*z**2)) + 1j*-x*z*(x**6 + x**4*(-24*y**2 + z**2) + 5*y**4*(-2*y**2 + z**2) + 5*x**2*(9*y**4 - 2*y**2*z**2)) )
            Bz = Bz * -1/8*rt3003*( 1/2*(x**6 - 15*x**4*y**2 + 15*x**2*y**4 - y**6)*(x**2 + y**2 - 14*z**2) + 1j*x*y*(3*x**4 - 10*x**2*y**2 + 3*y**4)*(x**2 + y**2 - 14*z**2) )
        elif m==7:
            rt429o2 = sqrt(429/2)
            Bx = Bx * 1/16*rt429o2*( (-8*x**8 - 7*y**6*(y**2 + z**2) - 105*x**4*y**2*(5*y**2 + z**2) + 7*x**6*(31*y**2 + z**2) + 7*x**2*(29*y**6 + 15*y**4*z**2)) + 1j*x*y*(-63*x**6 + 57*y**6 + 42*y**4*z**2 + 7*x**4*(61*y**2 + 6*z**2) - 7*x**2*(59*y**4 + 20*y**2*z**2)) )
            By = By * 1/16*rt429o2*( x*y*(-57*x**6 + 63*y**6 - 42*y**4*z**2 + 7*x**4*(59*y**2 - 6*z**2) - 7*x**2*(61*y**4 - 20*y**2*z**2)) + 1j*(7*x**8 + 8*y**8 - 7*y**6*z**2 + 7*x**6*(-29*y**2 + z**2) + 105*x**4*(5*y**4 - y**2*z**2) - 7*x**2*(31*y**6 - 15*y**4*z**2)) )
            Bz = Bz * 15/16*rt429o2*( -x*(x**6 - 21*x**4*y**2 + 35*x**2*y**4 - 7*y**6)*z + 1j*y*(-7*x**6 + 35*x**4*y**2 - 21*x**2*y**4 + y**6)*z )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi, n={n} and m is not between -n and n.")

    elif n==8:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   n=8 moment components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A8 = 3/256*sqrt(17/np.pi)
        Bx = A8*Binm/r**19
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==-8:
            rt2 = sqrt(2)
            rt715 = sqrt(715)
            Bx = Bx * rt715*( 1/rt2*x*(9*x**8 + 73*y**8 + 56*y**6*z**2 - 4*x**6*(79*y**2 + 2*z**2) + 14*x**4*(77*y**4 + 12*y**2*z**2) - 140*x**2*(5*y**6 + 2*y**4*z**2)) + 1j*-4*rt2*y*(10*x**8 + y**6*(y**2 + z**2) - 7*x**6*(13*y**2 + z**2) + 7*x**4*(19*y**4 + 5*y**2*z**2) - x**2*(37*y**6 + 21*y**4*z**2)) )
            By = By * rt715*( 1/rt2*y*(73*x**8 + 9*y**8 - 8*y**6*z**2 + x**6*(-700*y**2 + 56*z**2) + 14*x**4*(77*y**4 - 20*y**2*z**2) - 4*x**2*(79*y**6 - 42*y**4*z**2)) + 1j*4*rt2*x*(x**8 + 10*y**8 - 7*y**6*z**2 + x**6*(-37*y**2 + z**2) + 7*x**4*(19*y**4 - 3*y**2*z**2) + x**2*(-91*y**6 + 35*y**4*z**2)) )
            Bz = Bz * 17*rt715*( 1/rt2*(x**8 - 28*x**6*y**2 + 70*x**4*y**4 - 28*x**2*y**6 + y**8)*z + 1j*-4*rt2*x*y*(x**6 - 7*x**4*y**2 + 7*x**2*y**4 - y**6)*z )
        elif m==-7:
            rt1430 = sqrt(1430)
            Bx = Bx * 2*rt1430*( z*(10*x**8 + 7*y**6*(y**2 + z**2) - 7*x**6*(37*y**2 + z**2) + 35*x**4*(17*y**4 + 3*y**2*z**2) - 7*x**2*(31*y**6 + 15*y**4*z**2)) + 1j*x*y*z*(-77*x**6 + 59*y**6 + 42*y**4*z**2 + 7*x**4*(71*y**2 + 6*z**2) - 35*x**2*(13*y**4 + 4*y**2*z**2)) )
            By = By * 2*rt1430*( x*y*z*(59*x**6 - 77*y**6 + 42*y**4*z**2 - 7*x**4*(65*y**2 - 6*z**2) + 7*x**2*(71*y**4 - 20*y**2*z**2)) + 1j*z*(7*x**8 + 10*y**8 - 7*y**6*z**2 + 7*x**6*(-31*y**2 + z**2) + 35*x**4*(17*y**4 - 3*y**2*z**2) - 7*x**2*(37*y**6 - 15*y**4*z**2)) )
            Bz = Bz * -2*rt1430*( x*(x**6 - 21*x**4*y**2 + 35*x**2*y**4 - 7*y**6)*(x**2 + y**2 - 16*z**2) + 1j*y*(-7*x**6 + 35*x**4*y**2 - 21*x**2*y**4 + y**6)*(x**2 + y**2 - 16*z**2) )
        elif m==-6:
            rt429 = sqrt(429)
            Bx = Bx * -2*rt429*( x*(3*x**8 - 54*x**6*(y**2 + z**2) + 14*x**4*(2*y**4 + 69*y**2*z**2 + 2*z**4) + 5*y**4*(-3*y**4 + 42*y**2*z**2 + 28*z**4) + 70*x**2*(y**6 - 19*y**4*z**2 - 4*y**2*z**4)) + 1j*2*y*(-10*x**8 + y**8 - 13*y**6*z**2 - 14*y**4*z**4 + 35*x**6*(y**2 + 5*z**2) + 7*x**4*(3*y**4 - 105*y**2*z**2 - 10*z**4) + x**2*(-23*y**6 + 357*y**4*z**2 + 140*y**2*z**4)) )
            By = By * 2*rt429*( y*(-15*x**8 + 3*y**8 - 54*y**6*z**2 + 28*y**4*z**4 + 70*x**6*(y**2 + 3*z**2) + 14*x**4*(2*y**4 - 95*y**2*z**2 + 10*z**4) + x**2*(-54*y**6 + 966*y**4*z**2 - 280*y**2*z**4)) + 1j*-2*x*(x**8 - x**6*(23*y**2 + 13*z**2) + 7*x**4*(3*y**4 + 51*y**2*z**2 - 2*z**4) - 5*y**4*(2*y**4 - 35*y**2*z**2 + 14*z**4) + 35*x**2*(y**6 - 21*y**4*z**2 + 4*y**2*z**4)) )
            Bz = Bz * 10*rt429*( -(x**6 - 15*x**4*y**2 + 15*x**2*y**4 - y**6)*z*(3*x**2 + 3*y**2 - 14*z**2) + 1j*2*x*y*(3*x**4 - 10*x**2*y**2 + 3*y**4)*z*(3*x**2 + 3*y**2 - 14*z**2) )
        elif m==-5:
            rt2002 = sqrt(2002)
            Bx = Bx * 10*rt2002*( -z*(2*x**8 - y**8 + 3*y**6*z**2 + 4*y**4*z**4 - x**6*(23*y**2 + 11*z**2) + x**4*(-5*y**4 + 125*y**2*z**2 + 4*z**4) + x**2*(19*y**6 - 85*y**4*z**2 - 24*y**2*z**4)) + 1j*x*y*z*(11*x**6 + 7*y**6 - 26*y**4*z**2 - 16*y**2*z**4 - x**4*(19*y**2 + 58*z**2) + x**2*(-23*y**4 + 140*y**2*z**2 + 16*z**4)) )
            By = By * -10*rt2002*( x*y*z*(7*x**6 + 11*y**6 - 58*y**4*z**2 + 16*y**2*z**4 - x**4*(23*y**2 + 26*z**2) + x**2*(-19*y**4 + 140*y**2*z**2 - 16*z**4)) + 1j*z*(x**8 - 2*y**8 + 11*y**6*z**2 - 4*y**4*z**4 - x**6*(19*y**2 + 3*z**2) + x**4*(5*y**4 + 85*y**2*z**2 - 4*z**4) + x**2*(23*y**6 - 125*y**4*z**2 + 24*y**2*z**4)) )
            Bz = Bz * 2*rt2002*( x*(x**4 - 10*x**2*y**2 + 5*y**4)*(x**4 + y**4 - 28*y**2*z**2 + 56*z**4 + 2*x**2*(y**2 - 14*z**2)) + 1j*-y*(5*x**4 - 10*x**2*y**2 + y**4)*(x**4 + y**4 - 28*y**2*z**2 + 56*z**4 + 2*x**2*(y**2 - 14*z**2)) )
        elif m==-4:
            rt154 = sqrt(154)
            Bx = Bx * rt154*( x*(9*x**8 - 4*x**6*(13*y**2 + 68*z**2) + x**4*(-106*y**4 + 1728*y**2*z**2 + 664*z**4) - 20*x**2*(y**6 - 68*y**4*z**2 + 212*y**2*z**4 + 8*z**6) + 5*y**2*(5*y**6 - 128*y**4*z**2 + 184*y**2*z**4 + 96*z**6)) + 1j*-4*y*(10*x**8 + y**8 - 23*y**6*z**2 + 16*y**4*z**4 + 40*y**2*z**6 + 5*x**6*(y**2 - 59*z**2) + x**4*(-19*y**4 + 115*y**2*z**2 + 680*z**4) - x**2*(13*y**6 - 387*y**4*z**2 + 760*y**2*z**4 + 120*z**6)) )
            By = By * rt154*( y*(25*x**8 + 9*y**8 - 272*y**6*z**2 + 664*y**4*z**4 - 160*y**2*z**6 - 20*x**6*(y**2 + 32*z**2) + x**4*(-106*y**4 + 1360*y**2*z**2 + 920*z**4) + x**2*(-52*y**6 + 1728*y**4*z**2 - 4240*y**2*z**4 + 480*z**6)) + 1j*4*x*(x**8 - x**6*(13*y**2 + 23*z**2) + x**4*(-19*y**4 + 387*y**2*z**2 + 16*z**4) + 5*x**2*(y**6 + 23*y**4*z**2 - 152*y**2*z**4 + 8*z**6) + 5*(2*y**8 - 59*y**6*z**2 + 136*y**4*z**4 - 24*y**2*z**6)) )
            Bz = Bz * 65*rt154*( (x**4 - 6*x**2*y**2 + y**4)*z*(x**4 + y**4 - 8*y**2*z**2 + 8*z**4 + 2*x**2*(y**2 - 4*z**2)) + 1j*-4*x*y*(x**2 - y**2)*z*(x**4 + y**4 - 8*y**2*z**2 + 8*z**4 + 2*x**2*(y**2 - 4*z**2)) )
        elif m==-3:
            rt2310 = sqrt(2310)
            Bx = Bx * 2*rt2310*( -z*(-10*x**8 - 3*y**8 + 17*y**6*z**2 + 4*y**4*z**4 - 16*y**2*z**6 + x**6*(19*y**2 + 87*z**2) + x**4*(65*y**4 - 225*y**2*z**2 - 108*z**4) + x**2*(33*y**6 - 295*y**4*z**2 + 312*y**2*z**4 + 16*z**6)) + 1j*x*y*z*(-33*x**6 + 19*y**6 - 138*y**4*z**2 + 96*y**2*z**4 + 32*z**6 + x**4*(-47*y**2 + 278*z**2) + 5*x**2*(y**4 + 28*y**2*z**2 - 64*z**4)) )
            By = By * 2*rt2310*( x*y*z*(19*x**6 - 33*y**6 + 278*y**4*z**2 - 320*y**2*z**4 + 32*z**6 + x**4*(5*y**2 - 138*z**2) + x**2*(-47*y**4 + 140*y**2*z**2 + 96*z**4)) + 1j*z*(3*x**8 + 10*y**8 - 87*y**6*z**2 + 108*y**4*z**4 - 16*y**2*z**6 - x**6*(33*y**2 + 17*z**2) + x**4*(-65*y**4 + 295*y**2*z**2 - 4*z**4) + x**2*(-19*y**6 + 225*y**4*z**2 - 312*y**2*z**4 + 16*z**6)) )
            Bz = Bz * 2*rt2310*( -x*(x**2 - 3*y**2)*(x**6 + y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6 + 3*x**4*(y**2 - 12*z**2) + 3*x**2*(y**4 - 24*y**2*z**2 + 40*z**4)) + 1j*y*(3*x**2 - y**2)*(x**6 + y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6 + 3*x**4*(y**2 - 12*z**2) + 3*x**2*(y**4 - 24*y**2*z**2 + 40*z**4)) )
        elif m==-2:
            rt35 = sqrt(35)
            Bx = Bx * 2*rt35*( -x*(9*x**8 - 13*y**8 + 454*y**6*z**2 - 1420*y**4*z**4 + 608*y**2*z**6 + 64*z**8 + 2*x**6*(7*y**2 - 169*z**2) - 2*x**4*(6*y**4 + 111*y**2*z**2 - 610*z**4) - 10*x**2*(3*y**6 - 57*y**4*z**2 + 20*y**2*z**4 + 80*z**6)) + 1j*2*y*(10*x**8 - y**8 + 29*y**6*z**2 - 50*y**4*z**4 - 48*y**2*z**6 + 32*z**8 + x**6*(29*y**2 - 367*z**2) + x**4*(27*y**4 - 705*y**2*z**2 + 1270*z**4) + x**2*(7*y**6 - 309*y**4*z**2 + 1220*y**2*z**4 - 752*z**6)) )
            By = By * -2*rt35*( y*(13*x**8 - 9*y**8 + 338*y**6*z**2 - 1220*y**4*z**4 + 800*y**2*z**6 - 64*z**8 + x**6*(30*y**2 - 454*z**2) + 2*x**4*(6*y**4 - 285*y**2*z**2 + 710*z**4) + x**2*(-14*y**6 + 222*y**4*z**2 + 200*y**2*z**4 - 608*z**6)) + 1j*2*x*(x**8 - 10*y**8 + 367*y**6*z**2 - 1270*y**4*z**4 + 752*y**2*z**6 - 32*z**8 - x**6*(7*y**2 + 29*z**2) + x**4*(-27*y**4 + 309*y**2*z**2 + 50*z**4) + x**2*(-29*y**6 + 705*y**4*z**2 - 1220*y**2*z**4 + 48*z**6)) )
            Bz = Bz * 22*rt35*( -(x**2 - y**2)*z*(7*x**6 + 7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6 + 7*x**4*(3*y**2 - 10*z**2) + 7*x**2*(3*y**4 - 20*y**2*z**2 + 16*z**4)) + 1j*2*x*y*z*(7*x**6 + 7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6 + 7*x**4*(3*y**2 - 10*z**2) + 7*x**2*(3*y**4 - 20*y**2*z**2 + 16*z**4)) )
        elif m==-1:
            rt2 = sqrt(2)
            Bx = Bx * 2*rt2*( -z*(350*x**8 - 35*y**8 + 245*y**6*z**2 - 56*y**4*z**4 - 272*y**2*z**6 + 64*z**8 + 35*x**6*(29*y**2 - 103*z**2) + 7*x**4*(135*y**4 - 995*y**2*z**2 + 872*z**4) + x**2*(245*y**6 - 3115*y**4*z**2 + 6048*y**2*z**4 - 2032*z**6)) + 1j*55*x*y*z*(7*x**6 + 7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6 + 7*x**4*(3*y**2 - 10*z**2) + 7*x**2*(3*y**4 - 20*y**2*z**2 + 16*z**4)) )
            By = By * -2*rt2*( 55*x*y*z*(7*x**6 + 7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6 + 7*x**4*(3*y**2 - 10*z**2) + 7*x**2*(3*y**4 - 20*y**2*z**2 + 16*z**4)) + 1j*z*(35*x**8 - 350*y**8 + 3605*y**6*z**2 - 6104*y**4*z**4 + 2032*y**2*z**6 - 64*z**8 - 245*x**6*(y**2 + z**2) + x**4*(-945*y**4 + 3115*y**2*z**2 + 56*z**4) + x**2*(-1015*y**6 + 6965*y**4*z**2 - 6048*y**2*z**4 + 272*z**6)) )
            Bz = Bz * 10*rt2*( x*(7*x**8 + 7*y**8 - 280*y**6*z**2 + 1120*y**4*z**4 - 896*y**2*z**6 + 128*z**8 + 28*x**6*(y**2 - 10*z**2) + 14*x**4*(3*y**4 - 60*y**2*z**2 + 80*z**4) + 28*x**2*(y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6)) + 1j*-y*(7*x**8 + 7*y**8 - 280*y**6*z**2 + 1120*y**4*z**4 - 896*y**2*z**6 + 128*z**8 + 28*x**6*(y**2 - 10*z**2) + 14*x**4*(3*y**4 - 60*y**2*z**2 + 80*z**4) + 28*x**2*(y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6)) )
        elif m==0:
            Bx = Bx * 15*x*(7*x**8 + 7*y**8 - 280*y**6*z**2 + 1120*y**4*z**4 - 896*y**2*z**6 + 128*z**8 + 28*x**6*(y**2 - 10*z**2) + 14*x**4*(3*y**4 - 60*y**2*z**2 + 80*z**4) + 28*x**2*(y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6))
            By = By * 15*y*(7*x**8 + 7*y**8 - 280*y**6*z**2 + 1120*y**4*z**4 - 896*y**2*z**6 + 128*z**8 + 28*x**6*(y**2 - 10*z**2) + 14*x**4*(3*y**4 - 60*y**2*z**2 + 80*z**4) + 28*x**2*(y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6))
            Bz = Bz * 3*z*(315*x**8 + 315*y**8 - 3360*y**6*z**2 + 6048*y**4*z**4 - 2304*y**2*z**6 + 128*z**8 + 420*x**6*(3*y**2 - 8*z**2) + 126*x**4*(15*y**4 - 80*y**2*z**2 + 48*z**4) + 36*x**2*(35*y**6 - 280*y**4*z**2 + 336*y**2*z**4 - 64*z**6))
        elif m==1:
            rt2 = sqrt(2)
            Bx = Bx * 2*rt2*( z*(350*x**8 - 35*y**8 + 245*y**6*z**2 - 56*y**4*z**4 - 272*y**2*z**6 + 64*z**8 + 35*x**6*(29*y**2 - 103*z**2) + 7*x**4*(135*y**4 - 995*y**2*z**2 + 872*z**4) + x**2*(245*y**6 - 3115*y**4*z**2 + 6048*y**2*z**4 - 2032*z**6)) + 1j*55*x*y*z*(7*x**6 + 7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6 + 7*x**4*(3*y**2 - 10*z**2) + 7*x**2*(3*y**4 - 20*y**2*z**2 + 16*z**4)) )
            By = By * 2*rt2*( 55*x*y*z*(7*x**6 + 7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6 + 7*x**4*(3*y**2 - 10*z**2) + 7*x**2*(3*y**4 - 20*y**2*z**2 + 16*z**4)) + 1j*-z*(35*x**8 - 350*y**8 + 3605*y**6*z**2 - 6104*y**4*z**4 + 2032*y**2*z**6 - 64*z**8 - 245*x**6*(y**2 + z**2) + x**4*(-945*y**4 + 3115*y**2*z**2 + 56*z**4) + x**2*(-1015*y**6 + 6965*y**4*z**2 - 6048*y**2*z**4 + 272*z**6)) )
            Bz = Bz * -10*rt2*( x*(7*x**8 + 7*y**8 - 280*y**6*z**2 + 1120*y**4*z**4 - 896*y**2*z**6 + 128*z**8 + 28*x**6*(y**2 - 10*z**2) + 14*x**4*(3*y**4 - 60*y**2*z**2 + 80*z**4) + 28*x**2*(y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6)) + 1j*y*(7*x**8 + 7*y**8 - 280*y**6*z**2 + 1120*y**4*z**4 - 896*y**2*z**6 + 128*z**8 + 28*x**6*(y**2 - 10*z**2) + 14*x**4*(3*y**4 - 60*y**2*z**2 + 80*z**4) + 28*x**2*(y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6)) )
        elif m==2:
            rt35 = sqrt(35)
            Bx = Bx * -2*rt35*( x*(9*x**8 - 13*y**8 + 454*y**6*z**2 - 1420*y**4*z**4 + 608*y**2*z**6 + 64*z**8 + 2*x**6*(7*y**2 - 169*z**2) - 2*x**4*(6*y**4 + 111*y**2*z**2 - 610*z**4) - 10*x**2*(3*y**6 - 57*y**4*z**2 + 20*y**2*z**4 + 80*z**6)) + 1j*2*y*(10*x**8 - y**8 + 29*y**6*z**2 - 50*y**4*z**4 - 48*y**2*z**6 + 32*z**8 + x**6*(29*y**2 - 367*z**2) + x**4*(27*y**4 - 705*y**2*z**2 + 1270*z**4) + x**2*(7*y**6 - 309*y**4*z**2 + 1220*y**2*z**4 - 752*z**6)) )
            By = By * 2*rt35*( -y*(13*x**8 - 9*y**8 + 338*y**6*z**2 - 1220*y**4*z**4 + 800*y**2*z**6 - 64*z**8 + x**6*(30*y**2 - 454*z**2) + 2*x**4*(6*y**4 - 285*y**2*z**2 + 710*z**4) + x**2*(-14*y**6 + 222*y**4*z**2 + 200*y**2*z**4 - 608*z**6)) + 1j*2*x*(x**8 - 10*y**8 + 367*y**6*z**2 - 1270*y**4*z**4 + 752*y**2*z**6 - 32*z**8 - x**6*(7*y**2 + 29*z**2) + x**4*(-27*y**4 + 309*y**2*z**2 + 50*z**4) + x**2*(-29*y**6 + 705*y**4*z**2 - 1220*y**2*z**4 + 48*z**6)) )
            Bz = Bz * -22*rt35*( (x**2 - y**2)*z*(7*x**6 + 7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6 + 7*x**4*(3*y**2 - 10*z**2) + 7*x**2*(3*y**4 - 20*y**2*z**2 + 16*z**4)) + 1j*2*x*y*z*(7*x**6 + 7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6 + 7*x**4*(3*y**2 - 10*z**2) + 7*x**2*(3*y**4 - 20*y**2*z**2 + 16*z**4)) )
        elif m==3:
            rt2310 = sqrt(2310)
            Bx = Bx * 2*rt2310*( z*(-10*x**8 - 3*y**8 + 17*y**6*z**2 + 4*y**4*z**4 - 16*y**2*z**6 + x**6*(19*y**2 + 87*z**2) + x**4*(65*y**4 - 225*y**2*z**2 - 108*z**4) + x**2*(33*y**6 - 295*y**4*z**2 + 312*y**2*z**4 + 16*z**6)) + 1j*x*y*z*(-33*x**6 + 19*y**6 - 138*y**4*z**2 + 96*y**2*z**4 + 32*z**6 + x**4*(-47*y**2 + 278*z**2) + 5*x**2*(y**4 + 28*y**2*z**2 - 64*z**4)) )
            By = By * 2*rt2310*( -x*y*z*(19*x**6 - 33*y**6 + 278*y**4*z**2 - 320*y**2*z**4 + 32*z**6 + x**4*(5*y**2 - 138*z**2) + x**2*(-47*y**4 + 140*y**2*z**2 + 96*z**4)) + 1j*z*(3*x**8 + 10*y**8 - 87*y**6*z**2 + 108*y**4*z**4 - 16*y**2*z**6 - x**6*(33*y**2 + 17*z**2) + x**4*(-65*y**4 + 295*y**2*z**2 - 4*z**4) + x**2*(-19*y**6 + 225*y**4*z**2 - 312*y**2*z**4 + 16*z**6)) )
            Bz = Bz * 2*rt2310*( x*(x**2 - 3*y**2)*(x**6 + y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6 + 3*x**4*(y**2 - 12*z**2) + 3*x**2*(y**4 - 24*y**2*z**2 + 40*z**4)) + 1j*y*(3*x**2 - y**2)*(x**6 + y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6 + 3*x**4*(y**2 - 12*z**2) + 3*x**2*(y**4 - 24*y**2*z**2 + 40*z**4)) )
        elif m==4:
            rt154 = sqrt(154)
            Bx = Bx * rt154*( x*(9*x**8 - 4*x**6*(13*y**2 + 68*z**2) + x**4*(-106*y**4 + 1728*y**2*z**2 + 664*z**4) - 20*x**2*(y**6 - 68*y**4*z**2 + 212*y**2*z**4 + 8*z**6) + 5*y**2*(5*y**6 - 128*y**4*z**2 + 184*y**2*z**4 + 96*z**6)) + 1j*4*y*(10*x**8 + y**8 - 23*y**6*z**2 + 16*y**4*z**4 + 40*y**2*z**6 + 5*x**6*(y**2 - 59*z**2) + x**4*(-19*y**4 + 115*y**2*z**2 + 680*z**4) - x**2*(13*y**6 - 387*y**4*z**2 + 760*y**2*z**4 + 120*z**6)) )
            By = By * rt154*( y*(25*x**8 + 9*y**8 - 272*y**6*z**2 + 664*y**4*z**4 - 160*y**2*z**6 - 20*x**6*(y**2 + 32*z**2) + x**4*(-106*y**4 + 1360*y**2*z**2 + 920*z**4) + x**2*(-52*y**6 + 1728*y**4*z**2 - 4240*y**2*z**4 + 480*z**6)) + 1j*-4*x*(x**8 - x**6*(13*y**2 + 23*z**2) + x**4*(-19*y**4 + 387*y**2*z**2 + 16*z**4) + 5*x**2*(y**6 + 23*y**4*z**2 - 152*y**2*z**4 + 8*z**6) + 5*(2*y**8 - 59*y**6*z**2 + 136*y**4*z**4 - 24*y**2*z**6)) )
            Bz = Bz * 65*rt154*( (x**4 - 6*x**2*y**2 + y**4)*z*(x**4 + y**4 - 8*y**2*z**2 + 8*z**4 + 2*x**2*(y**2 - 4*z**2)) + 1j*4*x*y*(x**2 - y**2)*z*(x**4 + y**4 - 8*y**2*z**2 + 8*z**4 + 2*x**2*(y**2 - 4*z**2)) )
        elif m==5:
            rt2002 = sqrt(2002)
            Bx = Bx * 10*rt2002*( z*(2*x**8 - y**8 + 3*y**6*z**2 + 4*y**4*z**4 - x**6*(23*y**2 + 11*z**2) + x**4*(-5*y**4 + 125*y**2*z**2 + 4*z**4) + x**2*(19*y**6 - 85*y**4*z**2 - 24*y**2*z**4)) + 1j*x*y*z*(11*x**6 + 7*y**6 - 26*y**4*z**2 - 16*y**2*z**4 - x**4*(19*y**2 + 58*z**2) + x**2*(-23*y**4 + 140*y**2*z**2 + 16*z**4)) )
            By = By * 10*rt2002*( x*y*z*(7*x**6 + 11*y**6 - 58*y**4*z**2 + 16*y**2*z**4 - x**4*(23*y**2 + 26*z**2) + x**2*(-19*y**4 + 140*y**2*z**2 - 16*z**4)) + 1j*-z*(x**8 - 2*y**8 + 11*y**6*z**2 - 4*y**4*z**4 - x**6*(19*y**2 + 3*z**2) + x**4*(5*y**4 + 85*y**2*z**2 - 4*z**4) + x**2*(23*y**6 - 125*y**4*z**2 + 24*y**2*z**4)) )
            Bz = Bz * -2*rt2002*( x*(x**4 - 10*x**2*y**2 + 5*y**4)*(x**4 + y**4 - 28*y**2*z**2 + 56*z**4 + 2*x**2*(y**2 - 14*z**2)) + 1j*y*(5*x**4 - 10*x**2*y**2 + y**4)*(x**4 + y**4 - 28*y**2*z**2 + 56*z**4 + 2*x**2*(y**2 - 14*z**2)) )
        elif m==6:
            rt429 = sqrt(429)
            Bx = Bx * 2*rt429*( -x*(3*x**8 - 54*x**6*(y**2 + z**2) + 14*x**4*(2*y**4 + 69*y**2*z**2 + 2*z**4) + 5*y**4*(-3*y**4 + 42*y**2*z**2 + 28*z**4) + 70*x**2*(y**6 - 19*y**4*z**2 - 4*y**2*z**4)) + 1j*2*y*(-10*x**8 + y**8 - 13*y**6*z**2 - 14*y**4*z**4 + 35*x**6*(y**2 + 5*z**2) + 7*x**4*(3*y**4 - 105*y**2*z**2 - 10*z**4) + x**2*(-23*y**6 + 357*y**4*z**2 + 140*y**2*z**4)) )
            By = By * 2*rt429*( y*(-15*x**8 + 3*y**8 - 54*y**6*z**2 + 28*y**4*z**4 + 70*x**6*(y**2 + 3*z**2) + 14*x**4*(2*y**4 - 95*y**2*z**2 + 10*z**4) + x**2*(-54*y**6 + 966*y**4*z**2 - 280*y**2*z**4)) + 1j*2*x*(x**8 - x**6*(23*y**2 + 13*z**2) + 7*x**4*(3*y**4 + 51*y**2*z**2 - 2*z**4) - 5*y**4*(2*y**4 - 35*y**2*z**2 + 14*z**4) + 35*x**2*(y**6 - 21*y**4*z**2 + 4*y**2*z**4)) )
            Bz = Bz * -10*rt429*( (x**6 - 15*x**4*y**2 + 15*x**2*y**4 - y**6)*z*(3*x**2 + 3*y**2 - 14*z**2) + 1j*2*x*y*(3*x**4 - 10*x**2*y**2 + 3*y**4)*z*(3*x**2 + 3*y**2 - 14*z**2) )
        elif m==7:
            rt1430 = sqrt(1430)
            Bx = Bx * 2*rt1430*( -z*(10*x**8 + 7*y**6*(y**2 + z**2) - 7*x**6*(37*y**2 + z**2) + 35*x**4*(17*y**4 + 3*y**2*z**2) - 7*x**2*(31*y**6 + 15*y**4*z**2)) + 1j*x*y*z*(-77*x**6 + 59*y**6 + 42*y**4*z**2 + 7*x**4*(71*y**2 + 6*z**2) - 35*x**2*(13*y**4 + 4*y**2*z**2)) )
            By = By * 2*rt1430*( -x*y*z*(59*x**6 - 77*y**6 + 42*y**4*z**2 - 7*x**4*(65*y**2 - 6*z**2) + 7*x**2*(71*y**4 - 20*y**2*z**2)) + 1j*z*(7*x**8 + 10*y**8 - 7*y**6*z**2 + 7*x**6*(-31*y**2 + z**2) + 35*x**4*(17*y**4 - 3*y**2*z**2) - 7*x**2*(37*y**6 - 15*y**4*z**2)) )
            Bz = Bz * 2*rt1430*( x*(x**6 - 21*x**4*y**2 + 35*x**2*y**4 - 7*y**6)*(x**2 + y**2 - 16*z**2) + 1j*-y*(-7*x**6 + 35*x**4*y**2 - 21*x**2*y**4 + y**6)*(x**2 + y**2 - 16*z**2) )
        elif m==8:
            rt2 = sqrt(2)
            rt715 = sqrt(715)
            Bx = Bx * rt715*( 1/rt2*x*(9*x**8 + 73*y**8 + 56*y**6*z**2 - 4*x**6*(79*y**2 + 2*z**2) + 14*x**4*(77*y**4 + 12*y**2*z**2) - 140*x**2*(5*y**6 + 2*y**4*z**2)) + 1j*4*rt2*y*(10*x**8 + y**6*(y**2 + z**2) - 7*x**6*(13*y**2 + z**2) + 7*x**4*(19*y**4 + 5*y**2*z**2) - x**2*(37*y**6 + 21*y**4*z**2)) )
            By = By * rt715*( 1/rt2*y*(73*x**8 + 9*y**8 - 8*y**6*z**2 + x**6*(-700*y**2 + 56*z**2) + 14*x**4*(77*y**4 - 20*y**2*z**2) - 4*x**2*(79*y**6 - 42*y**4*z**2)) + 1j*-4*rt2*x*(x**8 + 10*y**8 - 7*y**6*z**2 + x**6*(-37*y**2 + z**2) + 7*x**4*(19*y**4 - 3*y**2*z**2) + x**2*(-91*y**6 + 35*y**4*z**2)) )
            Bz = Bz * 17*rt715*( 1/rt2*(x**8 - 28*x**6*y**2 + 70*x**4*y**4 - 28*x**2*y**6 + y**8)*z + 1j*4*rt2*x*y*(x**6 - 7*x**4*y**2 + 7*x**2*y**4 - y**6)*z )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi, n={n} and m is not between -n and n.")

    elif n==9:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   n=9 moment components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A9 = sqrt(95/2/np.pi)/256
        Bx = A9*Binm/r**21
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==-9:
            rt2431o2 = sqrt(2431/2)
            Bx = Bx * rt2431o2*( (10*x**10 - 9*y**8*(y**2 + z**2) + 252*x**6*y**2*(8*y**2 + z**2) - 9*x**8*(49*y**2 + z**2) - 42*x**4*(47*y**6 + 15*y**4*z**2) + 18*x**2*(23*y**8 + 14*y**6*z**2)) + 1j*x*y*(-99*x**8 - 91*y**8 - 72*y**6*z**2 + 12*x**6*(97*y**2 + 6*z**2) - 126*x**4*(19*y**4 + 4*y**2*z**2) + 36*x**2*(31*y**6 + 14*y**4*z**2)) )
            By = By * rt2431o2*( x*y*(91*x**8 + 99*y**8 - 72*y**6*z**2 - 36*x**6*(31*y**2 - 2*z**2) + 126*x**4*(19*y**4 - 4*y**2*z**2) - 12*x**2*(97*y**6 - 42*y**4*z**2)) + 1j*(9*x**10 - 10*y**10 + 9*y**8*z**2 + 9*x**8*(-46*y**2 + z**2) + 42*x**6*(47*y**4 - 6*y**2*z**2) - 126*x**4*(16*y**6 - 5*y**4*z**2) + 63*x**2*(7*y**8 - 4*y**6*z**2)) )
            Bz = Bz * 19*rt2431o2*( x*(x**8 - 36*x**6*y**2 + 126*x**4*y**4 - 84*x**2*y**6 + 9*y**8)*z + 1j*-y*(9*x**8 - 84*x**6*y**2 + 126*x**4*y**4 - 36*x**2*y**6 + y**8)*z )
        elif m==-8:
            rt2431 = sqrt(2431)
            Bx = Bx * 3*rt2431*( x*z*(11*x**8 + 75*y**8 + 56*y**6*z**2 - 4*x**6*(93*y**2 + 2*z**2) + 42*x**4*(29*y**4 + 4*y**2*z**2) - 28*x**2*(27*y**6 + 10*y**4*z**2)) + 1j*-8*y*z*(12*x**8 + y**6*(y**2 + z**2) - 7*x**6*(15*y**2 + z**2) + 7*x**4*(21*y**4 + 5*y**2*z**2) - 3*x**2*(13*y**6 + 7*y**4*z**2)) )
            By = By * 3*rt2431*( y*z*(75*x**8 + 11*y**8 - 8*y**6*z**2 + x**6*(-756*y**2 + 56*z**2) + 14*x**4*(87*y**4 - 20*y**2*z**2) + x**2*(-372*y**6 + 168*y**4*z**2)) + 1j*8*x*z*(x**8 + 12*y**8 - 7*y**6*z**2 + x**6*(-39*y**2 + z**2) + 21*x**4*(7*y**4 - y**2*z**2) - 35*x**2*(3*y**6 - y**4*z**2)) )
            Bz = Bz * 3*rt2431*( -(x**8 - 28*x**6*y**2 + 70*x**4*y**4 - 28*x**2*y**6 + y**8)*(x**2 + y**2 - 18*z**2) + 1j*8*x*y*(x**6 - 7*x**4*y**2 + 7*x**2*y**4 - y**6)*(x**2 + y**2 - 18*z**2) )
        elif m==-7:
            rt143o2 = sqrt(143/2)
            Bx = Bx * 3*rt143o2*( -(10*x**10 - 3*x**8*(83*y**2 + 67*z**2) + 7*y**6*(y**4 - 15*y**2*z**2 - 16*z**4) + 28*x**6*(12*y**4 + 177*y**2*z**2 + 4*z**4) + 42*x**4*(9*y**6 - 255*y**4*z**2 - 40*y**2*z**4) - 42*x**2*(5*y**8 - 86*y**6*z**2 - 40*y**4*z**4)) + 1j*x*y*(77*x**8 - 59*y**8 + 936*y**6*z**2 + 672*y**4*z**4 - 84*x**6*(5*y**2 + 18*z**2) - 42*x**4*(y**4 - 220*y**2*z**2 - 16*z**4) + 4*x**2*(99*y**6 - 1974*y**4*z**2 - 560*y**2*z**4)) )
            By = By * -3*rt143o2*( x*y*(59*x**8 - 36*x**6*(11*y**2 + 26*z**2) + 42*x**4*(y**4 + 188*y**2*z**2 - 16*z**4) - 7*y**4*(11*y**4 - 216*y**2*z**2 + 96*z**4) + 140*x**2*(3*y**6 - 66*y**4*z**2 + 16*y**2*z**4)) + 1j*(7*x**10 - 105*x**8*(2*y**2 + z**2) + 14*x**6*(27*y**4 + 258*y**2*z**2 - 8*z**4) + y**6*(10*y**4 - 201*y**2*z**2 + 112*z**4) + 42*x**4*(8*y**6 - 255*y**4*z**2 + 40*y**2*z**4) - 3*x**2*(83*y**8 - 1652*y**6*z**2 + 560*y**4*z**4)) )
            Bz = Bz * -51*rt143o2*( x*(x**6 - 21*x**4*y**2 + 35*x**2*y**4 - 7*y**6)*z*(3*x**2 + 3*y**2 - 16*z**2) + 1j*y*(-7*x**6 + 35*x**4*y**2 - 21*x**2*y**4 + y**6)*z*(3*x**2 + 3*y**2 - 16*z**2) )
        elif m==-6:
            rt858 = sqrt(858)
            Bx = Bx * 2*rt858*( -x*z*(33*x**8 - 141*y**8 + 602*y**6*z**2 + 420*y**4*z**4 - 2*x**6*(285*y**2 + 103*z**2) + 42*x**4*(6*y**4 + 83*y**2*z**2 + 2*z**4) + 42*x**2*(17*y**6 - 105*y**4*z**2 - 20*y**2*z**4)) + 1j*2*y*z*(108*x**8 - 9*y**8 + 33*y**6*z**2 + 42*y**4*z**4 - 21*x**6*(17*y**2 + 31*z**2) - 7*x**4*(33*y**4 - 365*y**2*z**2 - 30*z**4) + 3*x**2*(75*y**6 - 371*y**4*z**2 - 140*y**2*z**4)) )
            By = By * 2*rt858*( y*z*(-141*x**8 + 33*y**8 - 206*y**6*z**2 + 84*y**4*z**4 + 14*x**6*(51*y**2 + 43*z**2) + 42*x**4*(6*y**4 - 105*y**2*z**2 + 10*z**4) - 6*x**2*(95*y**6 - 581*y**4*z**2 + 140*y**2*z**4)) + 1j*-2*x*z*(9*x**8 - 3*x**6*(75*y**2 + 11*z**2) + 21*x**4*(11*y**4 + 53*y**2*z**2 - 2*z**4) - 3*y**4*(36*y**4 - 217*y**2*z**2 + 70*z**4) + 7*x**2*(51*y**6 - 365*y**4*z**2 + 60*y**2*z**4)) )
            Bz = Bz * 2*rt858*( (x**6 - 15*x**4*y**2 + 15*x**2*y**4 - y**6)*(3*x**4 + 3*y**4 - 96*y**2*z**2 + 224*z**4 + 6*x**2*(y**2 - 16*z**2)) + 1j*-2*x*y*(3*x**4 - 10*x**2*y**2 + 3*y**4)*(3*x**4 + 3*y**4 - 96*y**2*z**2 + 224*z**4 + 6*x**2*(y**2 - 16*z**2)) )
        elif m==-5:
            rt1430 = sqrt(1430)
            Bx = Bx * 3*rt1430*( (2*x**10 - 3*x**8*(7*y**2 + 23*z**2) - 28*x**6*(y**4 - 27*y**2*z**2 - 7*z**4) + 14*x**4*(y**6 + 15*y**4*z**2 - 150*y**2*z**4 - 4*z**6) - y**4*(y**6 - 27*y**4*z**2 + 28*y**2*z**4 + 56*z**6) + 6*x**2*(3*y**8 - 98*y**6*z**2 + 210*y**4*z**4 + 56*y**2*z**6)) + 1j*-x*y*(11*x**8 + 7*y**8 - 204*y**6*z**2 + 336*y**4*z**4 + 224*y**2*z**6 - 4*x**6*(2*y**2 + 93*z**2) - 42*x**4*(y**4 - 14*y**2*z**2 - 24*z**4) - 4*x**2*(4*y**6 - 189*y**4*z**2 + 560*y**2*z**4 + 56*z**6)) )
            By = By * 3*rt1430*( x*y*(7*x**8 + 11*y**8 - 372*y**6*z**2 + 1008*y**4*z**4 - 224*y**2*z**6 - 4*x**6*(4*y**2 + 51*z**2) - 42*x**4*(y**4 - 18*y**2*z**2 - 8*z**4) - 4*x**2*(2*y**6 - 147*y**4*z**2 + 560*y**2*z**4 - 56*z**6)) + 1j*(x**10 - 2*y**10 + 69*y**8*z**2 - 196*y**6*z**4 + 56*y**4*z**6 - 9*x**8*(2*y**2 + 3*z**2) - 14*x**6*(y**4 - 42*y**2*z**2 - 2*z**4) + 14*x**4*(2*y**6 - 15*y**4*z**2 - 90*y**2*z**4 + 4*z**6) + 21*x**2*(y**8 - 36*y**6*z**2 + 100*y**4*z**4 - 16*y**2*z**6)) )
            Bz = Bz * 3*rt1430*( x*(x**4 - 10*x**2*y**2 + 5*y**4)*z*(15*x**4 + 15*y**4 - 140*y**2*z**2 + 168*z**4 + 10*x**2*(3*y**2 - 14*z**2)) + 1j*-y*(5*x**4 - 10*x**2*y**2 + y**4)*z*(15*x**4 + 15*y**4 - 140*y**2*z**2 + 168*z**4 + 10*x**2*(3*y**2 - 14*z**2)) )
        elif m==-4:
            rt1001 = sqrt(1001)
            Bx = Bx * 6*rt1001*( x*z*(11*x**8 + 27*y**8 - 224*y**6*z**2 + 168*y**4*z**4 + 96*y**2*z**6 - 4*x**6*(15*y**2 + 28*z**2) - 42*x**4*(3*y**4 - 16*y**2*z**2 - 4*z**4) - 4*x**2*(7*y**6 - 140*y**4*z**2 + 252*y**2*z**4 + 8*z**6)) + 1j*-4*y*z*(12*x**8 + y**8 - 7*y**6*z**2 + 8*y**2*z**6 + 7*x**6*(y**2 - 17*z**2) - 7*x**4*(3*y**4 - 5*y**2*z**2 - 24*z**4) - 3*x**2*(5*y**6 - 49*y**4*z**2 + 56*y**2*z**4 + 8*z**6)) )
            By = By * 6*rt1001*( y*z*(27*x**8 + 11*y**8 - 112*y**6*z**2 + 168*y**4*z**4 - 32*y**2*z**6 - 28*x**6*(y**2 + 8*z**2) - 14*x**4*(9*y**4 - 40*y**2*z**2 - 12*z**4) - 12*x**2*(5*y**6 - 56*y**4*z**2 + 84*y**2*z**4 - 8*z**6)) + 1j*4*x*z*(x**8 + 12*y**8 - 119*y**6*z**2 + 168*y**4*z**4 - 24*y**2*z**6 - x**6*(15*y**2 + 7*z**2) - 21*x**4*(y**4 - 7*y**2*z**2) + x**2*(7*y**6 + 35*y**4*z**2 - 168*y**2*z**4 + 8*z**6)) )
            Bz = Bz * 6*rt1001*( -(x**4 - 6*x**2*y**2 + y**4)*(x**6 + y**6 - 42*y**4*z**2 + 168*y**2*z**4 - 112*z**6 + 3*x**4*(y**2 - 14*z**2) + 3*x**2*(y**4 - 28*y**2*z**2 + 56*z**4)) + 1j*4*x*y*(x**2 - y**2)*(x**6 + y**6 - 42*y**4*z**2 + 168*y**2*z**4 - 112*z**6 + 3*x**4*(y**2 - 14*z**2) + 3*x**2*(y**4 - 28*y**2*z**2 + 56*z**4)) )
        elif m==-3:
            rt462 = sqrt(462)
            Bx = Bx * rt462*( (-10*x**10 + 9*x**8*(y**2 + 49*z**2) + 84*x**6*(y**4 - 9*y**2*z**2 - 23*z**4) + 14*x**4*(7*y**6 - 195*y**4*z**2 + 330*y**2*z**4 + 116*z**6) - 3*y**2*(y**8 - 35*y**6*z**2 + 84*y**4*z**4 + 56*y**2*z**6 - 64*z**8) + 6*x**2*(5*y**8 - 238*y**6*z**2 + 1050*y**4*z**4 - 728*y**2*z**6 - 32*z**8)) + 1j*x*y*(33*x**8 - 19*y**8 + 756*y**6*z**2 - 2688*y**4*z**4 + 1120*y**2*z**6 + 384*z**8 + 4*x**6*(20*y**2 - 357*z**2) + 42*x**4*(y**4 - 50*y**2*z**2 + 144*z**4) - 12*x**2*(2*y**6 - 7*y**4*z**2 - 280*y**2*z**4 + 392*z**6)) )
            By = By * rt462*( x*y*(-19*x**8 + 33*y**8 - 1428*y**6*z**2 + 6048*y**4*z**4 - 4704*y**2*z**6 + 384*z**8 + x**6*(-24*y**2 + 756*z**2) + 42*x**4*(y**4 + 2*y**2*z**2 - 64*z**4) + 20*x**2*(4*y**6 - 105*y**4*z**2 + 168*y**2*z**4 + 56*z**6)) + 1j*(-3*x**10 - 10*y**10 + 441*y**8*z**2 - 1932*y**6*z**4 + 1624*y**4*z**6 - 192*y**2*z**8 + 15*x**8*(2*y**2 + 7*z**2) + 14*x**6*(7*y**4 - 102*y**2*z**2 - 18*z**4) + 42*x**4*(2*y**6 - 65*y**4*z**2 + 150*y**2*z**4 - 4*z**6) + 3*x**2*(3*y**8 - 252*y**6*z**2 + 1540*y**4*z**4 - 1456*y**2*z**6 + 64*z**8)) )
            Bz = Bz * 13*rt462*( -x*(x**2 - 3*y**2)*z*(7*x**6 + 7*y**6 - 84*y**4*z**2 + 168*y**2*z**4 - 64*z**6 + 21*x**4*(y**2 - 4*z**2) + 21*x**2*(y**4 - 8*y**2*z**2 + 8*z**4)) + 1j*y*(3*x**2 - y**2)*z*(7*x**6 + 7*y**6 - 84*y**4*z**2 + 168*y**2*z**4 - 64*z**6 + 21*x**4*(y**2 - 4*z**2) + 21*x**2*(y**4 - 8*y**2*z**2 + 8*z**4)) )
        elif m==-2:
            rt22 = sqrt(22)
            Bx = Bx * 6*rt22*( -x*z*(77*x**8 - 105*y**8 + 1218*y**6*z**2 - 2268*y**4*z**4 + 672*y**2*z**6 + 64*z**8 + 42*x**6*(3*y**2 - 23*z**2) - 42*x**4*(2*y**4 + 17*y**2*z**2 - 50*z**4) - 2*x**2*(119*y**6 - 735*y**4*z**2 + 84*y**2*z**4 + 496*z**6)) + 1j*2*y*z*(84*x**8 - 7*y**8 + 63*y**6*z**2 - 42*y**4*z**4 - 80*y**2*z**6 + 32*z**8 + 49*x**6*(5*y**2 - 21*z**2) + 21*x**4*(11*y**4 - 95*y**2*z**2 + 102*z**4) + 3*x**2*(21*y**6 - 301*y**4*z**2 + 700*y**2*z**4 - 304*z**6)) )
            By = By * -6*rt22*( y*z*(105*x**8 - 77*y**8 + 966*y**6*z**2 - 2100*y**4*z**4 + 992*y**2*z**6 - 64*z**8 + 14*x**6*(17*y**2 - 87*z**2) + 42*x**4*(2*y**4 - 35*y**2*z**2 + 54*z**4) - 42*x**2*(3*y**6 - 17*y**4*z**2 - 4*y**2*z**4 + 16*z**6)) + 1j*2*x*z*(7*x**8 - 84*y**8 + 1029*y**6*z**2 - 2142*y**4*z**4 + 912*y**2*z**6 - 32*z**8 - 63*x**6*(y**2 + z**2) + x**4*(-231*y**4 + 903*y**2*z**2 + 42*z**4) - 5*x**2*(49*y**6 - 399*y**4*z**2 + 420*y**2*z**4 - 16*z**6)) )
            Bz = Bz * 6*rt22*( (x**2 - y**2)*(7*x**8 + 7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8 + 28*x**6*(y**2 - 12*z**2) + 42*x**4*(y**4 - 24*y**2*z**2 + 40*z**4) + 28*x**2*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6)) + 1j*-2*x*y*(7*x**8 + 7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8 + 28*x**6*(y**2 - 12*z**2) + 42*x**4*(y**4 - 24*y**2*z**2 + 40*z**4) + 28*x**2*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6)) )
        elif m==-1:
            Bx = Bx * 3*( (70*x**10 - 7*y**10 + 273*y**8*z**2 - 840*y**6*z**4 - 224*y**4*z**6 + 768*y**2*z**8 - 128*z**10 + 21*x**8*(13*y**2 - 163*z**2) + 196*x**6*(2*y**4 - 51*y**2*z**2 + 90*z**4) + 14*x**4*(17*y**6 - 675*y**4*z**2 + 2460*y**2*z**4 - 1424*z**6) + 6*x**2*(7*y**8 - 434*y**6*z**2 + 2660*y**4*z**4 - 3360*y**2*z**6 + 832*z**8)) + 1j*-11*x*y*(7*x**8 + 7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8 + 28*x**6*(y**2 - 12*z**2) + 42*x**4*(y**4 - 24*y**2*z**2 + 40*z**4) + 28*x**2*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6)) )
            By = By * 3*( 11*x*y*(7*x**8 + 7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8 + 28*x**6*(y**2 - 12*z**2) + 42*x**4*(y**4 - 24*y**2*z**2 + 40*z**4) + 28*x**2*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6)) + 1j*(7*x**10 - 70*y**10 + 3423*y**8*z**2 - 17640*y**6*z**4 + 19936*y**4*z**6 - 4992*y**2*z**8 + 128*z**10 - 21*x**8*(2*y**2 + 13*z**2) + x**6*(-238*y**4 + 2604*y**2*z**2 + 840*z**4) - 14*x**4*(28*y**6 - 675*y**4*z**2 + 1140*y**2*z**4 - 16*z**6) - 3*x**2*(91*y**8 - 3332*y**6*z**2 + 11480*y**4*z**4 - 6720*y**2*z**6 + 256*z**8)) )
            Bz = Bz * 33*( x*z*(63*x**8 + 63*y**8 - 840*y**6*z**2 + 2016*y**4*z**4 - 1152*y**2*z**6 + 128*z**8 + 84*x**6*(3*y**2 - 10*z**2) + 126*x**4*(3*y**4 - 20*y**2*z**2 + 16*z**4) + 36*x**2*(7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6)) + 1j*-y*z*(63*x**8 + 63*y**8 - 840*y**6*z**2 + 2016*y**4*z**4 - 1152*y**2*z**6 + 128*z**8 + 84*x**6*(3*y**2 - 10*z**2) + 126*x**4*(3*y**4 - 20*y**2*z**2 + 16*z**4) + 36*x**2*(7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6)) )
        elif m==0:
            rt10 = sqrt(10)
            Bx = Bx * 11*rt10*x*z*(63*x**8 + 63*y**8 - 840*y**6*z**2 + 2016*y**4*z**4 - 1152*y**2*z**6 + 128*z**8 + 84*x**6*(3*y**2 - 10*z**2) + 126*x**4*(3*y**4 - 20*y**2*z**2 + 16*z**4) + 36*x**2*(7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6))
            By = By * 11*rt10*y*z*(63*x**8 + 63*y**8 - 840*y**6*z**2 + 2016*y**4*z**4 - 1152*y**2*z**6 + 128*z**8 + 84*x**6*(3*y**2 - 10*z**2) + 126*x**4*(3*y**4 - 20*y**2*z**2 + 16*z**4) + 36*x**2*(7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6))
            Bz = Bz * -rt10*(63*x**10 + 63*y**10 - 3150*y**8*z**2 + 16800*y**6*z**4 - 20160*y**4*z**6 + 5760*y**2*z**8 - 256*z**10 + 315*x**8*(y**2 - 10*z**2) + 210*x**6*(3*y**4 - 60*y**2*z**2 + 80*z**4) + 630*x**4*(y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6) + 45*x**2*(7*y**8 - 280*y**6*z**2 + 1120*y**4*z**4 - 896*y**2*z**6 + 128*z**8))
        elif m==1:
            Bx = Bx * -3*( (70*x**10 - 7*y**10 + 273*y**8*z**2 - 840*y**6*z**4 - 224*y**4*z**6 + 768*y**2*z**8 - 128*z**10 + 21*x**8*(13*y**2 - 163*z**2) + 196*x**6*(2*y**4 - 51*y**2*z**2 + 90*z**4) + 14*x**4*(17*y**6 - 675*y**4*z**2 + 2460*y**2*z**4 - 1424*z**6) + 6*x**2*(7*y**8 - 434*y**6*z**2 + 2660*y**4*z**4 - 3360*y**2*z**6 + 832*z**8)) + 1j*11*x*y*(7*x**8 + 7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8 + 28*x**6*(y**2 - 12*z**2) + 42*x**4*(y**4 - 24*y**2*z**2 + 40*z**4) + 28*x**2*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6)) )
            By = By * 3*( -11*x*y*(7*x**8 + 7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8 + 28*x**6*(y**2 - 12*z**2) + 42*x**4*(y**4 - 24*y**2*z**2 + 40*z**4) + 28*x**2*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6)) + 1j*(7*x**10 - 70*y**10 + 3423*y**8*z**2 - 17640*y**6*z**4 + 19936*y**4*z**6 - 4992*y**2*z**8 + 128*z**10 - 21*x**8*(2*y**2 + 13*z**2) + x**6*(-238*y**4 + 2604*y**2*z**2 + 840*z**4) - 14*x**4*(28*y**6 - 675*y**4*z**2 + 1140*y**2*z**4 - 16*z**6) - 3*x**2*(91*y**8 - 3332*y**6*z**2 + 11480*y**4*z**4 - 6720*y**2*z**6 + 256*z**8)) )
            Bz = Bz * -33*( x*z*(63*x**8 + 63*y**8 - 840*y**6*z**2 + 2016*y**4*z**4 - 1152*y**2*z**6 + 128*z**8 + 84*x**6*(3*y**2 - 10*z**2) + 126*x**4*(3*y**4 - 20*y**2*z**2 + 16*z**4) + 36*x**2*(7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6)) + 1j*y*z*(63*x**8 + 63*y**8 - 840*y**6*z**2 + 2016*y**4*z**4 - 1152*y**2*z**6 + 128*z**8 + 84*x**6*(3*y**2 - 10*z**2) + 126*x**4*(3*y**4 - 20*y**2*z**2 + 16*z**4) + 36*x**2*(7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6)) )
        elif m==2:
            rt22 = sqrt(22)
            Bx = Bx * -6*rt22*( x*z*(77*x**8 - 105*y**8 + 1218*y**6*z**2 - 2268*y**4*z**4 + 672*y**2*z**6 + 64*z**8 + 42*x**6*(3*y**2 - 23*z**2) - 42*x**4*(2*y**4 + 17*y**2*z**2 - 50*z**4) - 2*x**2*(119*y**6 - 735*y**4*z**2 + 84*y**2*z**4 + 496*z**6)) + 1j*2*y*z*(84*x**8 - 7*y**8 + 63*y**6*z**2 - 42*y**4*z**4 - 80*y**2*z**6 + 32*z**8 + 49*x**6*(5*y**2 - 21*z**2) + 21*x**4*(11*y**4 - 95*y**2*z**2 + 102*z**4) + 3*x**2*(21*y**6 - 301*y**4*z**2 + 700*y**2*z**4 - 304*z**6)) )
            By = By * 6*rt22*( -y*z*(105*x**8 - 77*y**8 + 966*y**6*z**2 - 2100*y**4*z**4 + 992*y**2*z**6 - 64*z**8 + 14*x**6*(17*y**2 - 87*z**2) + 42*x**4*(2*y**4 - 35*y**2*z**2 + 54*z**4) - 42*x**2*(3*y**6 - 17*y**4*z**2 - 4*y**2*z**4 + 16*z**6)) + 1j*2*x*z*(7*x**8 - 84*y**8 + 1029*y**6*z**2 - 2142*y**4*z**4 + 912*y**2*z**6 - 32*z**8 - 63*x**6*(y**2 + z**2) + x**4*(-231*y**4 + 903*y**2*z**2 + 42*z**4) - 5*x**2*(49*y**6 - 399*y**4*z**2 + 420*y**2*z**4 - 16*z**6)) )
            Bz = Bz * 6*rt22*( (x**2 - y**2)*(7*x**8 + 7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8 + 28*x**6*(y**2 - 12*z**2) + 42*x**4*(y**4 - 24*y**2*z**2 + 40*z**4) + 28*x**2*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6)) + 1j*2*x*y*(7*x**8 + 7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8 + 28*x**6*(y**2 - 12*z**2) + 42*x**4*(y**4 - 24*y**2*z**2 + 40*z**4) + 28*x**2*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6)) )
        elif m==3:
            rt462 = sqrt(462)
            Bx = Bx * rt462*( (10*x**10 - 9*x**8*(y**2 + 49*z**2) - 84*x**6*(y**4 - 9*y**2*z**2 - 23*z**4) - 14*x**4*(7*y**6 - 195*y**4*z**2 + 330*y**2*z**4 + 116*z**6) + 3*y**2*(y**8 - 35*y**6*z**2 + 84*y**4*z**4 + 56*y**2*z**6 - 64*z**8) - 6*x**2*(5*y**8 - 238*y**6*z**2 + 1050*y**4*z**4 - 728*y**2*z**6 - 32*z**8)) + 1j*x*y*(33*x**8 - 19*y**8 + 756*y**6*z**2 - 2688*y**4*z**4 + 1120*y**2*z**6 + 384*z**8 + 4*x**6*(20*y**2 - 357*z**2) + 42*x**4*(y**4 - 50*y**2*z**2 + 144*z**4) - 12*x**2*(2*y**6 - 7*y**4*z**2 - 280*y**2*z**4 + 392*z**6)) )
            By = By * rt462*( x*y*(19*x**8 - 33*y**8 + 1428*y**6*z**2 - 6048*y**4*z**4 + 4704*y**2*z**6 - 384*z**8 + 12*x**6*(2*y**2 - 63*z**2) - 42*x**4*(y**4 + 2*y**2*z**2 - 64*z**4) - 20*x**2*(4*y**6 - 105*y**4*z**2 + 168*y**2*z**4 + 56*z**6)) + 1j*(-3*x**10 - 10*y**10 + 441*y**8*z**2 - 1932*y**6*z**4 + 1624*y**4*z**6 - 192*y**2*z**8 + 15*x**8*(2*y**2 + 7*z**2) + 14*x**6*(7*y**4 - 102*y**2*z**2 - 18*z**4) + 42*x**4*(2*y**6 - 65*y**4*z**2 + 150*y**2*z**4 - 4*z**6) + 3*x**2*(3*y**8 - 252*y**6*z**2 + 1540*y**4*z**4 - 1456*y**2*z**6 + 64*z**8)) )
            Bz = Bz * 13*rt462*( x*(x**2 - 3*y**2)*z*(7*x**6 + 7*y**6 - 84*y**4*z**2 + 168*y**2*z**4 - 64*z**6 + 21*x**4*(y**2 - 4*z**2) + 21*x**2*(y**4 - 8*y**2*z**2 + 8*z**4)) + 1j*y*(3*x**2 - y**2)*z*(7*x**6 + 7*y**6 - 84*y**4*z**2 + 168*y**2*z**4 - 64*z**6 + 21*x**4*(y**2 - 4*z**2) + 21*x**2*(y**4 - 8*y**2*z**2 + 8*z**4)) )
        elif m==4:
            rt1001 = sqrt(1001)
            Bx = Bx * 6*rt1001*( x*z*(11*x**8 + 27*y**8 - 224*y**6*z**2 + 168*y**4*z**4 + 96*y**2*z**6 - 4*x**6*(15*y**2 + 28*z**2) - 42*x**4*(3*y**4 - 16*y**2*z**2 - 4*z**4) - 4*x**2*(7*y**6 - 140*y**4*z**2 + 252*y**2*z**4 + 8*z**6)) + 1j*4*y*z*(12*x**8 + y**8 - 7*y**6*z**2 + 8*y**2*z**6 + 7*x**6*(y**2 - 17*z**2) - 7*x**4*(3*y**4 - 5*y**2*z**2 - 24*z**4) - 3*x**2*(5*y**6 - 49*y**4*z**2 + 56*y**2*z**4 + 8*z**6)) )
            By = By * 6*rt1001*( y*z*(27*x**8 + 11*y**8 - 112*y**6*z**2 + 168*y**4*z**4 - 32*y**2*z**6 - 28*x**6*(y**2 + 8*z**2) - 14*x**4*(9*y**4 - 40*y**2*z**2 - 12*z**4) - 12*x**2*(5*y**6 - 56*y**4*z**2 + 84*y**2*z**4 - 8*z**6)) + 1j*-4*x*z*(x**8 + 12*y**8 - 119*y**6*z**2 + 168*y**4*z**4 - 24*y**2*z**6 - x**6*(15*y**2 + 7*z**2) - 21*x**4*(y**4 - 7*y**2*z**2) + x**2*(7*y**6 + 35*y**4*z**2 - 168*y**2*z**4 + 8*z**6)) )
            Bz = Bz * -6*rt1001*( (x**4 - 6*x**2*y**2 + y**4)*(x**6 + y**6 - 42*y**4*z**2 + 168*y**2*z**4 - 112*z**6 + 3*x**4*(y**2 - 14*z**2) + 3*x**2*(y**4 - 28*y**2*z**2 + 56*z**4)) + 1j*4*x*y*(x**2 - y**2)*(x**6 + y**6 - 42*y**4*z**2 + 168*y**2*z**4 - 112*z**6 + 3*x**4*(y**2 - 14*z**2) + 3*x**2*(y**4 - 28*y**2*z**2 + 56*z**4)) )
        elif m==5:
            rt1430 = sqrt(1430)
            Bx = Bx * -3*rt1430*( (2*x**10 - 3*x**8*(7*y**2 + 23*z**2) - 28*x**6*(y**4 - 27*y**2*z**2 - 7*z**4) + 14*x**4*(y**6 + 15*y**4*z**2 - 150*y**2*z**4 - 4*z**6) - y**4*(y**6 - 27*y**4*z**2 + 28*y**2*z**4 + 56*z**6) + 6*x**2*(3*y**8 - 98*y**6*z**2 + 210*y**4*z**4 + 56*y**2*z**6)) + 1j*x*y*(11*x**8 + 7*y**8 - 204*y**6*z**2 + 336*y**4*z**4 + 224*y**2*z**6 - 4*x**6*(2*y**2 + 93*z**2) - 42*x**4*(y**4 - 14*y**2*z**2 - 24*z**4) - 4*x**2*(4*y**6 - 189*y**4*z**2 + 560*y**2*z**4 + 56*z**6)) )
            By = By * 3*rt1430*( -x*y*(7*x**8 + 11*y**8 - 372*y**6*z**2 + 1008*y**4*z**4 - 224*y**2*z**6 - 4*x**6*(4*y**2 + 51*z**2) - 42*x**4*(y**4 - 18*y**2*z**2 - 8*z**4) - 4*x**2*(2*y**6 - 147*y**4*z**2 + 560*y**2*z**4 - 56*z**6)) + 1j*(x**10 - 2*y**10 + 69*y**8*z**2 - 196*y**6*z**4 + 56*y**4*z**6 - 9*x**8*(2*y**2 + 3*z**2) - 14*x**6*(y**4 - 42*y**2*z**2 - 2*z**4) + 14*x**4*(2*y**6 - 15*y**4*z**2 - 90*y**2*z**4 + 4*z**6) + 21*x**2*(y**8 - 36*y**6*z**2 + 100*y**4*z**4 - 16*y**2*z**6)) )
            Bz = Bz * -3*rt1430*( x*(x**4 - 10*x**2*y**2 + 5*y**4)*z*(15*x**4 + 15*y**4 - 140*y**2*z**2 + 168*z**4 + 10*x**2*(3*y**2 - 14*z**2)) + 1j*y*(5*x**4 - 10*x**2*y**2 + y**4)*z*(15*x**4 + 15*y**4 - 140*y**2*z**2 + 168*z**4 + 10*x**2*(3*y**2 - 14*z**2)) )
        elif m==6:
            rt858 = sqrt(858)
            Bx = Bx * -2*rt858*( x*z*(33*x**8 - 141*y**8 + 602*y**6*z**2 + 420*y**4*z**4 - 2*x**6*(285*y**2 + 103*z**2) + 42*x**4*(6*y**4 + 83*y**2*z**2 + 2*z**4) + 42*x**2*(17*y**6 - 105*y**4*z**2 - 20*y**2*z**4)) + 1j*2*y*z*(108*x**8 - 9*y**8 + 33*y**6*z**2 + 42*y**4*z**4 - 21*x**6*(17*y**2 + 31*z**2) - 7*x**4*(33*y**4 - 365*y**2*z**2 - 30*z**4) + 3*x**2*(75*y**6 - 371*y**4*z**2 - 140*y**2*z**4)) )
            By = By * 2*rt858*( y*z*(-141*x**8 + 33*y**8 - 206*y**6*z**2 + 84*y**4*z**4 + 14*x**6*(51*y**2 + 43*z**2) + 42*x**4*(6*y**4 - 105*y**2*z**2 + 10*z**4) - 6*x**2*(95*y**6 - 581*y**4*z**2 + 140*y**2*z**4)) + 1j*2*x*z*(9*x**8 - 3*x**6*(75*y**2 + 11*z**2) + 21*x**4*(11*y**4 + 53*y**2*z**2 - 2*z**4) - 3*y**4*(36*y**4 - 217*y**2*z**2 + 70*z**4) + 7*x**2*(51*y**6 - 365*y**4*z**2 + 60*y**2*z**4)) )
            Bz = Bz * 2*rt858*( (x**6 - 15*x**4*y**2 + 15*x**2*y**4 - y**6)*(3*x**4 + 3*y**4 - 96*y**2*z**2 + 224*z**4 + 6*x**2*(y**2 - 16*z**2)) + 1j*2*x*y*(3*x**4 - 10*x**2*y**2 + 3*y**4)*(3*x**4 + 3*y**4 - 96*y**2*z**2 + 224*z**4 + 6*x**2*(y**2 - 16*z**2)) )
        elif m==7:
            rt143o2 = sqrt(143/2)
            Bx = Bx * 3*rt143o2*( (10*x**10 - 3*x**8*(83*y**2 + 67*z**2) + 7*y**6*(y**4 - 15*y**2*z**2 - 16*z**4) + 28*x**6*(12*y**4 + 177*y**2*z**2 + 4*z**4) + 42*x**4*(9*y**6 - 255*y**4*z**2 - 40*y**2*z**4) - 42*x**2*(5*y**8 - 86*y**6*z**2 - 40*y**4*z**4)) + 1j*x*y*(77*x**8 - 59*y**8 + 936*y**6*z**2 + 672*y**4*z**4 - 84*x**6*(5*y**2 + 18*z**2) - 42*x**4*(y**4 - 220*y**2*z**2 - 16*z**4) + 4*x**2*(99*y**6 - 1974*y**4*z**2 - 560*y**2*z**4)) )
            By = By * 3*rt143o2*( x*y*(59*x**8 - 36*x**6*(11*y**2 + 26*z**2) + 42*x**4*(y**4 + 188*y**2*z**2 - 16*z**4) - 7*y**4*(11*y**4 - 216*y**2*z**2 + 96*z**4) + 140*x**2*(3*y**6 - 66*y**4*z**2 + 16*y**2*z**4)) + 1j*-(7*x**10 - 105*x**8*(2*y**2 + z**2) + 14*x**6*(27*y**4 + 258*y**2*z**2 - 8*z**4) + y**6*(10*y**4 - 201*y**2*z**2 + 112*z**4) + 42*x**4*(8*y**6 - 255*y**4*z**2 + 40*y**2*z**4) - 3*x**2*(83*y**8 - 1652*y**6*z**2 + 560*y**4*z**4)) )
            Bz = Bz * 51*rt143o2*( x*(x**6 - 21*x**4*y**2 + 35*x**2*y**4 - 7*y**6)*z*(3*x**2 + 3*y**2 - 16*z**2) + 1j*-y*(-7*x**6 + 35*x**4*y**2 - 21*x**2*y**4 + y**6)*z*(3*x**2 + 3*y**2 - 16*z**2) )
        elif m==8:
            rt2431 = sqrt(2431)
            Bx = Bx * 3*rt2431*( x*z*(11*x**8 + 75*y**8 + 56*y**6*z**2 - 4*x**6*(93*y**2 + 2*z**2) + 42*x**4*(29*y**4 + 4*y**2*z**2) - 28*x**2*(27*y**6 + 10*y**4*z**2)) + 1j*8*y*z*(12*x**8 + y**6*(y**2 + z**2) - 7*x**6*(15*y**2 + z**2) + 7*x**4*(21*y**4 + 5*y**2*z**2) - 3*x**2*(13*y**6 + 7*y**4*z**2)) )
            By = By * 3*rt2431*( y*z*(75*x**8 + 11*y**8 - 8*y**6*z**2 + x**6*(-756*y**2 + 56*z**2) + 14*x**4*(87*y**4 - 20*y**2*z**2) + x**2*(-372*y**6 + 168*y**4*z**2)) + 1j*-8*x*z*(x**8 + 12*y**8 - 7*y**6*z**2 + x**6*(-39*y**2 + z**2) + 21*x**4*(7*y**4 - y**2*z**2) - 35*x**2*(3*y**6 - y**4*z**2)) )
            Bz = Bz * -3*rt2431*( (x**8 - 28*x**6*y**2 + 70*x**4*y**4 - 28*x**2*y**6 + y**8)*(x**2 + y**2 - 18*z**2) + 1j*8*x*y*(x**6 - 7*x**4*y**2 + 7*x**2*y**4 - y**6)*(x**2 + y**2 - 18*z**2) )
        elif m==9:
            rt2431o2 = sqrt(2431/2)
            Bx = Bx * rt2431o2*( (-10*x**10 + 9*y**8*(y**2 + z**2) - 252*x**6*y**2*(8*y**2 + z**2) + 9*x**8*(49*y**2 + z**2) + 42*x**4*(47*y**6 + 15*y**4*z**2) - 18*x**2*(23*y**8 + 14*y**6*z**2)) + 1j*x*y*(-99*x**8 - 91*y**8 - 72*y**6*z**2 + 12*x**6*(97*y**2 + 6*z**2) - 126*x**4*(19*y**4 + 4*y**2*z**2) + 36*x**2*(31*y**6 + 14*y**4*z**2)) )
            By = By * rt2431o2*( x*y*(-91*x**8 - 99*y**8 + 72*y**6*z**2 + 36*x**6*(31*y**2 - 2*z**2) - 126*x**4*(19*y**4 - 4*y**2*z**2) + 12*x**2*(97*y**6 - 42*y**4*z**2)) + 1j*(9*x**10 - 10*y**10 + 9*y**8*z**2 + 9*x**8*(-46*y**2 + z**2) + 42*x**6*(47*y**4 - 6*y**2*z**2) - 126*x**4*(16*y**6 - 5*y**4*z**2) + 63*x**2*(7*y**8 - 4*y**6*z**2)) )
            Bz = Bz * -19*rt2431o2*( x*(x**8 - 36*x**6*y**2 + 126*x**4*y**4 - 84*x**2*y**6 + 9*y**8)*z + 1j*y*(9*x**8 - 84*x**6*y**2 + 126*x**4*y**4 - 36*x**2*y**6 + y**8)*z )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi, n={n} and m is not between -n and n.")

    elif n==10:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   n=10 moment components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A10 = sqrt(1155/2/np.pi)/64
        Bx = A10*Binm/r**23
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==-10:
            rt4199o10 = sqrt(4199/10)
            Bx = Bx * 1/4*rt4199o10*( 1/2*x*(11*x**10 - 5*x**8*(119*y**2 + 2*z**2) - 3*y**8*(37*y**2 + 30*z**2) + 90*x**6*(39*y**4 + 4*y**2*z**2) - 210*x**4*(23*y**6 + 6*y**4*z**2) + 15*x**2*(113*y**8 + 56*y**6*z**2)) + 1j*y*(-60*x**10 + 5*y**8*(y**2 + z**2) + 15*x**8*(59*y**2 + 3*z**2) - 84*x**6*(29*y**4 + 5*y**2*z**2) + 90*x**4*(19*y**6 + 7*y**4*z**2) - 20*x**2*(14*y**8 + 9*y**6*z**2)) )
            By = By * 1/4*rt4199o10*( 1/2*y*(111*x**10 - 11*y**10 + 10*y**8*z**2 + x**8*(-1695*y**2 + 90*z**2) + 210*x**6*(23*y**4 - 4*y**2*z**2) - 90*x**4*(39*y**6 - 14*y**4*z**2) + 5*x**2*(119*y**8 - 72*y**6*z**2)) + 1j*x*(5*x**10 - 60*y**10 + 45*y**8*z**2 + 5*x**8*(-56*y**2 + z**2) + 90*x**6*(19*y**4 - 2*y**2*z**2) - 42*x**4*(58*y**6 - 15*y**4*z**2) + 15*x**2*(59*y**8 - 28*y**6*z**2)) )
            Bz = Bz * 21/4*rt4199o10*( 1/2*(x**10 - 45*x**8*y**2 + 210*x**6*y**4 - 210*x**4*y**6 + 45*x**2*y**8 - y**10)*z + 1j*-x*y*(5*x**8 - 60*x**6*y**2 + 126*x**4*y**4 - 60*x**2*y**6 + 5*y**8)*z )
        elif m==-9:
            rt4199o2 = sqrt(4199/2)
            Bx = Bx * 3/4*rt4199o2*( z*(4*x**10 - 3*y**8*(y**2 + z**2) + 84*x**6*y**2*(9*y**2 + z**2) - 3*x**8*(57*y**2 + z**2) - 42*x**4*(17*y**6 + 5*y**4*z**2) + 12*x**2*(12*y**8 + 7*y**6*z**2)) + 1j*-x*y*z*(39*x**8 + 31*y**8 + 24*y**6*z**2 - 12*x**6*(37*y**2 + 2*z**2) + 42*x**4*(21*y**4 + 4*y**2*z**2) - 12*x**2*(33*y**6 + 14*y**4*z**2)) )
            By = By * 3/4*rt4199o2*( x*y*z*(31*x**8 + 39*y**8 - 24*y**6*z**2 + x**6*(-396*y**2 + 24*z**2) + 42*x**4*(21*y**4 - 4*y**2*z**2) + x**2*(-444*y**6 + 168*y**4*z**2)) + 1j*z*(3*x**10 - 4*y**10 + 3*y**8*z**2 + 3*x**8*(-48*y**2 + z**2) + 42*x**6*(17*y**4 - 2*y**2*z**2) - 42*x**4*(18*y**6 - 5*y**4*z**2) + 3*x**2*(57*y**8 - 28*y**6*z**2)) )
            Bz = Bz * 1/4*rt4199o2*( -x*(x**8 - 36*x**6*y**2 + 126*x**4*y**4 - 84*x**2*y**6 + 9*y**8)*(x**2 + y**2 - 20*z**2) + 1j*y*(9*x**8 - 84*x**6*y**2 + 126*x**4*y**4 - 36*x**2*y**6 + y**8)*(x**2 + y**2 - 20*z**2) )
        elif m==-8:
            rt221 = sqrt(221)
            Bx = Bx * rt221*( 1/8*x*(-11*x**10 + x**8*(361*y**2 + 244*z**2) - 18*x**6*(47*y**4 + 440*y**2*z**2 + 8*z**4) + 3*y**6*(-25*y**4 + 444*y**2*z**2 + 336*z**4) - 42*x**4*(11*y**6 - 588*y**4*z**2 - 72*y**2*z**4) + 3*x**2*(227*y**8 - 4816*y**6*z**2 - 1680*y**4*z**4)) + 1j*y*(12*x**10 - 3*x**8*(31*y**2 + 87*z**2) + y**6*(y**4 - 17*y**2*z**2 - 18*z**4) + 42*x**6*(y**4 + 52*y**2*z**2 + 3*z**4) + 18*x**4*(6*y**6 - 161*y**4*z**2 - 35*y**2*z**4) + x**2*(-38*y**8 + 720*y**6*z**2 + 378*y**4*z**4)) )
            By = By * rt221*( 1/8*y*(-75*x**10 - 11*y**10 + 244*y**8*z**2 - 144*y**6*z**4 + 3*x**8*(227*y**2 + 444*z**2) - 42*x**6*(11*y**4 + 344*y**2*z**2 - 24*z**4) - 18*x**4*(47*y**6 - 1372*y**4*z**2 + 280*y**2*z**4) + x**2*(361*y**8 - 7920*y**6*z**2 + 3024*y**4*z**4)) + 1j*-x*(x**10 - x**8*(38*y**2 + 17*z**2) + 18*x**6*(6*y**4 + 40*y**2*z**2 - z**4) + 3*y**6*(4*y**4 - 87*y**2*z**2 + 42*z**4) + 42*x**4*(y**6 - 69*y**4*z**2 + 9*y**2*z**4) - 3*x**2*(31*y**8 - 728*y**6*z**2 + 210*y**4*z**4)) )
            Bz = Bz * 57*rt221*( -1/8*(x**8 - 28*x**6*y**2 + 70*x**4*y**4 - 28*x**2*y**6 + y**8)*z*(x**2 + y**2 - 6*z**2) + 1j*x*y*(x**6 - 7*x**4*y**2 + 7*x**2*y**4 - y**6)*z*(x**2 + y**2 - 6*z**2) )
        elif m==-7:
            rt663o2 = sqrt(663/2)
            Bx = Bx * 1/4*rt663o2*( z*(-36*x**10 + x**8*(867*y**2 + 251*z**2) - 28*x**6*(39*y**4 + 211*y**2*z**2 + 4*z**4) + 7*y**6*(-3*y**4 + 13*y**2*z**2 + 16*z**4) - 42*x**4*(31*y**6 - 285*y**4*z**2 - 40*y**2*z**4) + 28*x**2*(24*y**8 - 131*y**6*z**2 - 60*y**4*z**4)) + 1j*x*y*z*(273*x**8 - 183*y**8 + 888*y**6*z**2 + 672*y**4*z**4 - 84*x**6*(17*y**2 + 22*z**2) - 14*x**4*(15*y**4 - 764*y**2*z**2 - 48*z**4) + 4*x**2*(327*y**6 - 2114*y**4*z**2 - 560*y**2*z**4)) )
            By = By * 1/4*rt663o2*( -x*y*z*(183*x**8 - 12*x**6*(109*y**2 + 74*z**2) + 14*x**4*(15*y**4 + 604*y**2*z**2 - 48*z**4) + 28*x**2*(51*y**6 - 382*y**4*z**2 + 80*y**2*z**4) - 21*(13*y**8 - 88*y**6*z**2 + 32*y**4*z**4)) + 1j*z*(-21*x**10 - 36*y**10 + 251*y**8*z**2 - 112*y**6*z**4 + 7*x**8*(96*y**2 + 13*z**2) - 14*x**6*(93*y**4 + 262*y**2*z**2 - 8*z**4) - 42*x**4*(26*y**6 - 285*y**4*z**2 + 40*y**2*z**4) + x**2*(867*y**8 - 5908*y**6*z**2 + 1680*y**4*z**4)) )
            Bz = Bz * 3/4*rt663o2*( x*(x**6 - 21*x**4*y**2 + 35*x**2*y**4 - 7*y**6)*(x**4 + y**4 - 36*y**2*z**2 + 96*z**4 + 2*x**2*(y**2 - 18*z**2)) + 1j*y*(-7*x**6 + 35*x**4*y**2 - 21*x**2*y**4 + y**6)*(x**4 + y**4 - 36*y**2*z**2 + 96*z**4 + 2*x**2*(y**2 - 18*z**2)) )
        elif m==-6:
            rt39o2 = sqrt(39/2)
            Bx = Bx * 3/4*rt39o2*( 1/2*x*(11*x**10 - x**8*(179*y**2 + 426*z**2) + x**6*(-106*y**4 + 7080*y**2*z**2 + 1376*z**4) + 14*x**4*(23*y**6 - 186*y**4*z**2 - 1584*y**2*z**4 - 32*z**6) - y**4*(47*y**6 - 1542*y**4*z**2 + 2912*y**2*z**4 + 2240*z**6) + x**2*(191*y**8 - 8568*y**6*z**2 + 25760*y**4*z**4 + 4480*y**2*z**6)) + 1j*y*(-36*x**10 + x**8*(83*y**2 + 1371*z**2) + 28*x**6*(7*y**4 - 153*y**2*z**2 - 152*z**4) + y**4*(3*y**6 - 93*y**4*z**2 + 128*y**2*z**4 + 224*z**6) + 2*x**4*(y**6 - 1491*y**4*z**2 + 7840*y**2*z**4 + 560*z**6) - 4*x**2*(18*y**8 - 645*y**6*z**2 + 1512*y**4*z**4 + 560*y**2*z**6)) )
            By = By * 3/4*rt39o2*( 1/2*y*(47*x**10 - 11*y**10 + 426*y**8*z**2 - 1376*y**6*z**4 + 448*y**4*z**6 - x**8*(191*y**2 + 1542*z**2) + x**6*(-322*y**4 + 8568*y**2*z**2 + 2912*z**4) + 2*x**4*(53*y**6 + 1302*y**4*z**2 - 12880*y**2*z**4 + 1120*z**6) + x**2*(179*y**8 - 7080*y**6*z**2 + 22176*y**4*z**4 - 4480*y**2*z**6)) + 1j*x*(3*x**10 - 36*y**10 + 1371*y**8*z**2 - 4256*y**6*z**4 + 1120*y**4*z**6 - 3*x**8*(24*y**2 + 31*z**2) + 2*x**6*(y**4 + 1290*y**2*z**2 + 64*z**4) + 14*x**4*(14*y**6 - 213*y**4*z**2 - 432*y**2*z**4 + 16*z**6) + x**2*(83*y**8 - 4284*y**6*z**2 + 15680*y**4*z**4 - 2240*y**2*z**6)) )
            Bz = Bz * 17/4*rt39o2*( 1/2*(x**6 - 15*x**4*y**2 + 15*x**2*y**4 - y**6)*z*(15*x**4 + 15*y**4 - 160*y**2*z**2 + 224*z**4 + 10*x**2*(3*y**2 - 16*z**2)) + 1j*-x*y*(3*x**4 - 10*x**2*y**2 + 3*y**4)*z*(15*x**4 + 15*y**4 - 160*y**2*z**2 + 224*z**4 + 10*x**2*(3*y**2 - 16*z**2)) )
        elif m==-5:
            rt39o10 = sqrt(39/10)
            Bx = Bx * -1/2*rt39o10*( z*(-180*x**10 + 5*x**8*(363*y**2 + 419*z**2) + 28*x**6*(90*y**4 - 785*y**2*z**2 - 131*z**4) - 210*x**4*(5*y**6 + 35*y**4*z**2 - 178*y**2*z**4 - 4*z**6) + 5*y**4*(15*y**6 - 125*y**4*z**2 + 28*y**2*z**4 + 168*z**6) - 20*x**2*(75*y**8 - 805*y**6*z**2 + 987*y**4*z**4 + 252*y**2*z**6)) + 1j*x*y*z*(975*x**8 + 555*y**8 - 5220*y**6*z**2 + 4368*y**4*z**4 + 3360*y**2*z**6 - 300*x**6*(2*y**2 + 37*z**2) - 70*x**4*(51*y**4 - 230*y**2*z**2 - 264*z**4) - 20*x**2*(72*y**6 - 1099*y**4*z**2 + 1904*y**2*z**4 + 168*z**6)) )
            By = By * 1/2*rt39o10*( x*y*z*(555*x**8 - 180*x**6*(8*y**2 + 29*z**2) + x**4*(-3570*y**4 + 21980*y**2*z**2 + 4368*z**4) - 20*x**2*(30*y**6 - 805*y**4*z**2 + 1904*y**2*z**4 - 168*z**6) + 15*(65*y**8 - 740*y**6*z**2 + 1232*y**4*z**4 - 224*y**2*z**6)) + 1j*z*(75*x**10 - 180*y**10 + 2095*y**8*z**2 - 3668*y**6*z**4 + 840*y**4*z**6 - 125*x**8*(12*y**2 + 5*z**2) - 70*x**6*(15*y**4 - 230*y**2*z**2 - 2*z**4) + 210*x**4*(12*y**6 - 35*y**4*z**2 - 94*y**2*z**4 + 4*z**6) + 5*x**2*(363*y**8 - 4396*y**6*z**2 + 7476*y**4*z**4 - 1008*y**2*z**6)) )
            Bz = Bz * 3/2*rt39o10*( -x*(x**4 - 10*x**2*y**2 + 5*y**4)*(5*x**6 + 5*y**6 - 240*y**4*z**2 + 1120*y**2*z**4 - 896*z**6 + 15*x**4*(y**2 - 16*z**2) + 5*x**2*(3*y**4 - 96*y**2*z**2 + 224*z**4)) + 1j*y*(5*x**4 - 10*x**2*y**2 + y**4)*(5*x**6 + 5*y**6 - 240*y**4*z**2 + 1120*y**2*z**4 - 896*z**6 + 15*x**4*(y**2 - 16*z**2) + 5*x**2*(3*y**4 - 96*y**2*z**2 + 224*z**4)) )
        elif m==-4:
            rt39 = sqrt(39)
            Bx = Bx * rt39*( -1/4*x*(11*x**10 - x**8*(49*y**2 + 556*z**2) - 6*x**6*(31*y**4 - 480*y**2*z**2 - 476*z**4) - 14*x**4*(11*y**6 - 444*y**4*z**2 + 1164*y**2*z**4 + 208*z**6) + 3*y**2*(9*y**8 - 404*y**6*z**2 + 1624*y**4*z**4 - 672*y**2*z**6 - 448*z**8) + x**2*(-y**8 + 1568*y**6*z**2 - 14280*y**4*z**4 + 16576*y**2*z**6 + 448*z**8)) + 1j*y*(12*x**10 + x**8*(19*y**2 - 597*z**2) - 14*x**6*(y**4 + 28*y**2*z**2 - 213*z**4) - 6*x**4*(6*y**6 - 161*y**4*z**2 + 105*y**2*z**4 + 476*z**6) + y**2*(y**8 - 41*y**6*z**2 + 126*y**4*z**4 + 56*y**2*z**6 - 112*z**8) + x**2*(-14*y**8 + 720*y**6*z**2 - 3486*y**4*z**4 + 2576*y**2*z**6 + 336*z**8)) )
            By = By * -rt39*( 1/4*y*(27*x**10 - x**8*(y**2 + 1212*z**2) + x**6*(-154*y**4 + 1568*y**2*z**2 + 4872*z**4) - 6*x**4*(31*y**6 - 1036*y**4*z**2 + 2380*y**2*z**4 + 336*z**6) + x**2*(-49*y**8 + 2880*y**6*z**2 - 16296*y**4*z**4 + 16576*y**2*z**6 - 1344*z**8) + y**2*(11*y**8 - 556*y**6*z**2 + 2856*y**4*z**4 - 2912*y**2*z**6 + 448*z**8)) + 1j*x*(x**10 - x**8*(14*y**2 + 41*z**2) - 18*x**6*(2*y**4 - 40*y**2*z**2 - 7*z**4) - 14*x**4*(y**6 - 69*y**4*z**2 + 249*y**2*z**4 - 4*z**6) + x**2*(19*y**8 - 392*y**6*z**2 - 630*y**4*z**4 + 2576*y**2*z**6 - 112*z**8) + 3*y**2*(4*y**8 - 199*y**6*z**2 + 994*y**4*z**4 - 952*y**2*z**6 + 112*z**8)) )
            Bz = Bz * 21*rt39*( -1/4*(x**4 - 6*x**2*y**2 + y**4)*z*(5*x**6 + 5*y**6 - 70*y**4*z**2 + 168*y**2*z**4 - 80*z**6 + 5*x**4*(3*y**2 - 14*z**2) + x**2*(15*y**4 - 140*y**2*z**2 + 168*z**4)) + 1j*x*y*(x**2 - y**2)*z*(5*x**6 + 5*y**6 - 70*y**4*z**2 + 168*y**2*z**4 - 80*z**6 + 5*x**4*(3*y**2 - 14*z**2) + x**2*(15*y**4 - 140*y**2*z**2 + 168*z**4)) )
        elif m==-3:
            rt39o2 = sqrt(39/2)
            Bx = Bx * 3/2*rt39o2*( -z*(28*x**10 - 7*x**8*(3*y**2 + 59*z**2) - 28*x**6*(8*y**4 - 23*y**2*z**2 - 39*z**4) - 2*x**4*(133*y**6 - 1225*y**4*z**2 + 1218*y**2*z**4 + 332*z**6) + y**2*(7*y**8 - 77*y**6*z**2 + 84*y**4*z**4 + 104*y**2*z**6 - 64*z**8) - 4*x**2*(21*y**8 - 329*y**6*z**2 + 861*y**4*z**4 - 420*y**2*z**6 - 16*z**8)) + 1j*x*y*z*(91*x**8 - 49*y**8 + 644*y**6*z**2 - 1344*y**4*z**4 + 352*y**2*z**6 + 128*z**8 + 28*x**6*(8*y**2 - 47*z**2) + 14*x**4*(9*y**4 - 142*y**2*z**2 + 240*z**4) - 4*x**2*(14*y**6 + 7*y**4*z**2 - 504*y**2*z**4 + 472*z**6)) )
            By = By * 3/2*rt39o2*( -x*y*z*(49*x**8 - 91*y**8 + 1316*y**6*z**2 - 3360*y**4*z**4 + 1888*y**2*z**6 - 128*z**8 + 28*x**6*(2*y**2 - 23*z**2) - 14*x**4*(9*y**4 - 2*y**2*z**2 - 96*z**4) - 4*x**2*(56*y**6 - 497*y**4*z**2 + 504*y**2*z**4 + 88*z**6)) + 1j*z*(-7*x**10 - 28*y**10 + 413*y**8*z**2 - 1092*y**6*z**4 + 664*y**4*z**6 - 64*y**2*z**8 + x**8*(84*y**2 + 77*z**2) + 14*x**6*(19*y**4 - 94*y**2*z**2 - 6*z**4) + 2*x**4*(112*y**6 - 1225*y**4*z**2 + 1722*y**2*z**4 - 52*z**6) + x**2*(21*y**8 - 644*y**6*z**2 + 2436*y**4*z**4 - 1680*y**2*z**6 + 64*z**8)) )
            Bz = Bz * 7/2*rt39o2*( x*(x**2 - 3*y**2)*(x**8 + y**8 - 56*y**6*z**2 + 336*y**4*z**4 - 448*y**2*z**6 + 128*z**8 + 4*x**6*(y**2 - 14*z**2) + 6*x**4*(y**4 - 28*y**2*z**2 + 56*z**4) + 4*x**2*(y**6 - 42*y**4*z**2 + 168*y**2*z**4 - 112*z**6)) + 1j*y*(-3*x**2 + y**2)*(x**8 + y**8 - 56*y**6*z**2 + 336*y**4*z**4 - 448*y**2*z**6 + 128*z**8 + 4*x**6*(y**2 - 14*z**2) + 6*x**4*(y**4 - 28*y**2*z**2 + 56*z**4) + 4*x**2*(y**6 - 42*y**4*z**2 + 168*y**2*z**4 - 112*z**6)) )
        elif m==-2:
            rt3 = sqrt(3)
            Bx = Bx * 1/4*rt3*( 1/2*x*(77*x**10 + 7*x**8*(29*y**2 - 634*z**2) + 42*x**6*(y**4 - 180*y**2*z**2 + 664*z**4) - 14*x**4*(23*y**6 - 282*y**4*z**2 - 1608*y**2*z**4 + 2896*z**6) + x**2*(-343*y**8 + 12824*y**6*z**2 - 38640*y**4*z**4 + 448*y**2*z**6 + 14464*z**8) - 3*(35*y**10 - 1918*y**8*z**2 + 11088*y**6*z**4 - 13664*y**4*z**6 + 2944*y**2*z**8 + 256*z**10)) + 1j*y*(-84*x**10 + 7*y**10 - 329*y**8*z**2 + 1344*y**6*z**4 - 112*y**4*z**6 - 1408*y**2*z**8 + 384*z**10 - 7*x**8*(47*y**2 - 681*z**2) - 28*x**6*(17*y**4 - 499*y**2*z**2 + 1044*z**4) - 42*x**4*(7*y**6 - 317*y**4*z**2 + 1360*y**2*z**4 - 968*z**6) - 4*x**2*(14*y**8 - 945*y**6*z**2 + 6636*y**4*z**4 - 10136*y**2*z**6 + 3264*z**8)) )
            By = By * 1/4*rt3*( 1/2*y*(105*x**10 - 77*y**10 + 4438*y**8*z**2 - 27888*y**6*z**4 + 40544*y**4*z**6 - 14464*y**2*z**8 + 768*z**10 + 7*x**8*(49*y**2 - 822*z**2) + 14*x**6*(23*y**4 - 916*y**2*z**2 + 2376*z**4) - 42*x**4*(y**6 + 94*y**4*z**2 - 920*y**2*z**4 + 976*z**6) + x**2*(-203*y**8 + 7560*y**6*z**2 - 22512*y**4*z**4 - 448*y**2*z**6 + 8832*z**8)) + 1j*x*(7*x**10 - 84*y**10 + 4767*y**8*z**2 - 29232*y**6*z**4 + 40656*y**4*z**6 - 13056*y**2*z**8 + 384*z**10 - 7*x**8*(8*y**2 + 47*z**2) - 42*x**6*(7*y**4 - 90*y**2*z**2 - 32*z**4) - 14*x**4*(34*y**6 - 951*y**4*z**2 + 1896*y**2*z**4 + 8*z**6) + x**2*(-329*y**8 + 13972*y**6*z**2 - 57120*y**4*z**4 + 40544*y**2*z**6 - 1408*z**8)) )
            Bz = Bz * 39/4*rt3*( 1/2*(x**2 - y**2)*z*(21*x**8 + 21*y**8 - 336*y**6*z**2 + 1008*y**4*z**4 - 768*y**2*z**6 + 128*z**8 + 84*x**6*(y**2 - 4*z**2) + 126*x**4*(y**4 - 8*y**2*z**2 + 8*z**4) + 12*x**2*(7*y**6 - 84*y**4*z**2 + 168*y**2*z**4 - 64*z**6)) + 1j*-x*y*z*(21*x**8 + 21*y**8 - 336*y**6*z**2 + 1008*y**4*z**4 - 768*y**2*z**6 + 128*z**8 + 84*x**6*(y**2 - 4*z**2) + 126*x**4*(y**4 - 8*y**2*z**2 + 8*z**4) + 12*x**2*(7*y**6 - 84*y**4*z**2 + 168*y**2*z**4 - 64*z**6)) )
        elif m==-1:
            Bx = Bx * 1/4*( z*(756*x**10 - 63*y**10 + 777*y**8*z**2 - 1176*y**6*z**4 - 864*y**4*z**6 + 1024*y**2*z**8 - 128*z**10 + 21*x**8*(141*y**2 - 587*z**2) + 84*x**6*(51*y**4 - 431*y**2*z**2 + 454*z**4) + 18*x**4*(147*y**6 - 1925*y**4*z**2 + 4172*y**2*z**4 - 1712*z**6) + 4*x**2*(126*y**8 - 2499*y**6*z**2 + 8946*y**4*z**4 - 7920*y**2*z**6 + 1504*z**8)) + 1j*-39*x*y*z*(21*x**8 + 21*y**8 - 336*y**6*z**2 + 1008*y**4*z**4 - 768*y**2*z**6 + 128*z**8 + 84*x**6*(y**2 - 4*z**2) + 126*x**4*(y**4 - 8*y**2*z**2 + 8*z**4) + 12*x**2*(7*y**6 - 84*y**4*z**2 + 168*y**2*z**4 - 64*z**6)) )
            By = By * 1/4*( 39*x*y*z*(21*x**8 + 21*y**8 - 336*y**6*z**2 + 1008*y**4*z**4 - 768*y**2*z**6 + 128*z**8 + 84*x**6*(y**2 - 4*z**2) + 126*x**4*(y**4 - 8*y**2*z**2 + 8*z**4) + 12*x**2*(7*y**6 - 84*y**4*z**2 + 168*y**2*z**4 - 64*z**6)) + 1j*z*(63*x**10 - 756*y**10 + 12327*y**8*z**2 - 38136*y**6*z**4 + 30816*y**4*z**6 - 6016*y**2*z**8 + 128*z**10 - 21*x**8*(24*y**2 + 37*z**2) - 294*x**6*(9*y**4 - 34*y**2*z**2 - 4*z**4) - 18*x**4*(238*y**6 - 1925*y**4*z**2 + 1988*y**2*z**4 - 48*z**6) + x**2*(-2961*y**8 + 36204*y**6*z**2 - 75096*y**4*z**4 + 31680*y**2*z**6 - 1024*z**8)) )
            Bz = Bz * 3/4*( -x*(21*x**10 + 21*y**10 - 1260*y**8*z**2 + 8400*y**6*z**4 - 13440*y**4*z**6 + 5760*y**2*z**8 - 512*z**10 + 105*x**8*(y**2 - 12*z**2) + 210*x**6*(y**4 - 24*y**2*z**2 + 40*z**4) + 210*x**4*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6) + 15*x**2*(7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8)) + 1j*y*(21*x**10 + 21*y**10 - 1260*y**8*z**2 + 8400*y**6*z**4 - 13440*y**4*z**6 + 5760*y**2*z**8 - 512*z**10 + 105*x**8*(y**2 - 12*z**2) + 210*x**6*(y**4 - 24*y**2*z**2 + 40*z**4) + 210*x**4*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6) + 15*x**2*(7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8)) )
        elif m==0:
            rt11o10 = sqrt(11/10)
            Bx = Bx * -3/4*rt11o10*x*(21*x**10 + 21*y**10 - 1260*y**8*z**2 + 8400*y**6*z**4 - 13440*y**4*z**6 + 5760*y**2*z**8 - 512*z**10 + 105*x**8*(y**2 - 12*z**2) + 210*x**6*(y**4 - 24*y**2*z**2 + 40*z**4) + 210*x**4*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6) + 15*x**2*(7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8))
            By = By * -3/4*rt11o10*y*(21*x**10 + 21*y**10 - 1260*y**8*z**2 + 8400*y**6*z**4 - 13440*y**4*z**6 + 5760*y**2*z**8 - 512*z**10 + 105*x**8*(y**2 - 12*z**2) + 210*x**6*(y**4 - 24*y**2*z**2 + 40*z**4) + 210*x**4*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6) + 15*x**2*(7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8))
            Bz = Bz * 1/4*rt11o10*z*(-693*x**10 - 693*y**10 + 11550*y**8*z**2 - 36960*y**6*z**4 + 31680*y**4*z**6 - 7040*y**2*z**8 + 256*z**10 - 1155*x**8*(3*y**2 - 10*z**2) - 2310*x**6*(3*y**4 - 20*y**2*z**2 + 16*z**4) - 990*x**4*(7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6) - 55*x**2*(63*y**8 - 840*y**6*z**2 + 2016*y**4*z**4 - 1152*y**2*z**6 + 128*z**8))
        elif m==1:
            Bx = Bx * 1/4*( z*(-756*x**10 + 63*y**10 - 777*y**8*z**2 + 1176*y**6*z**4 + 864*y**4*z**6 - 1024*y**2*z**8 + 128*z**10 - 21*x**8*(141*y**2 - 587*z**2) - 84*x**6*(51*y**4 - 431*y**2*z**2 + 454*z**4) - 18*x**4*(147*y**6 - 1925*y**4*z**2 + 4172*y**2*z**4 - 1712*z**6) - 4*x**2*(126*y**8 - 2499*y**6*z**2 + 8946*y**4*z**4 - 7920*y**2*z**6 + 1504*z**8)) + 1j*-39*x*y*z*(21*x**8 + 21*y**8 - 336*y**6*z**2 + 1008*y**4*z**4 - 768*y**2*z**6 + 128*z**8 + 84*x**6*(y**2 - 4*z**2) + 126*x**4*(y**4 - 8*y**2*z**2 + 8*z**4) + 12*x**2*(7*y**6 - 84*y**4*z**2 + 168*y**2*z**4 - 64*z**6)) )
            By = By * 1/4*( -39*x*y*z*(21*x**8 + 21*y**8 - 336*y**6*z**2 + 1008*y**4*z**4 - 768*y**2*z**6 + 128*z**8 + 84*x**6*(y**2 - 4*z**2) + 126*x**4*(y**4 - 8*y**2*z**2 + 8*z**4) + 12*x**2*(7*y**6 - 84*y**4*z**2 + 168*y**2*z**4 - 64*z**6)) + 1j*z*(63*x**10 - 756*y**10 + 12327*y**8*z**2 - 38136*y**6*z**4 + 30816*y**4*z**6 - 6016*y**2*z**8 + 128*z**10 - 21*x**8*(24*y**2 + 37*z**2) - 294*x**6*(9*y**4 - 34*y**2*z**2 - 4*z**4) - 18*x**4*(238*y**6 - 1925*y**4*z**2 + 1988*y**2*z**4 - 48*z**6) + x**2*(-2961*y**8 + 36204*y**6*z**2 - 75096*y**4*z**4 + 31680*y**2*z**6 - 1024*z**8)) )
            Bz = Bz * 3/4*( x*(21*x**10 + 21*y**10 - 1260*y**8*z**2 + 8400*y**6*z**4 - 13440*y**4*z**6 + 5760*y**2*z**8 - 512*z**10 + 105*x**8*(y**2 - 12*z**2) + 210*x**6*(y**4 - 24*y**2*z**2 + 40*z**4) + 210*x**4*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6) + 15*x**2*(7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8)) + 1j*y*(21*x**10 + 21*y**10 - 1260*y**8*z**2 + 8400*y**6*z**4 - 13440*y**4*z**6 + 5760*y**2*z**8 - 512*z**10 + 105*x**8*(y**2 - 12*z**2) + 210*x**6*(y**4 - 24*y**2*z**2 + 40*z**4) + 210*x**4*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6) + 15*x**2*(7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8)) )
        elif m==2:
            rt3 = sqrt(3)
            Bx = Bx * 1/4*rt3*( 1/2*x*(77*x**10 + 7*x**8*(29*y**2 - 634*z**2) + 42*x**6*(y**4 - 180*y**2*z**2 + 664*z**4) - 14*x**4*(23*y**6 - 282*y**4*z**2 - 1608*y**2*z**4 + 2896*z**6) + x**2*(-343*y**8 + 12824*y**6*z**2 - 38640*y**4*z**4 + 448*y**2*z**6 + 14464*z**8) - 3*(35*y**10 - 1918*y**8*z**2 + 11088*y**6*z**4 - 13664*y**4*z**6 + 2944*y**2*z**8 + 256*z**10)) + 1j*-y*(-84*x**10 + 7*y**10 - 329*y**8*z**2 + 1344*y**6*z**4 - 112*y**4*z**6 - 1408*y**2*z**8 + 384*z**10 - 7*x**8*(47*y**2 - 681*z**2) - 28*x**6*(17*y**4 - 499*y**2*z**2 + 1044*z**4) - 42*x**4*(7*y**6 - 317*y**4*z**2 + 1360*y**2*z**4 - 968*z**6) - 4*x**2*(14*y**8 - 945*y**6*z**2 + 6636*y**4*z**4 - 10136*y**2*z**6 + 3264*z**8)) )
            By = By * 1/4*rt3*( 1/2*y*(105*x**10 - 77*y**10 + 4438*y**8*z**2 - 27888*y**6*z**4 + 40544*y**4*z**6 - 14464*y**2*z**8 + 768*z**10 + 7*x**8*(49*y**2 - 822*z**2) + 14*x**6*(23*y**4 - 916*y**2*z**2 + 2376*z**4) - 42*x**4*(y**6 + 94*y**4*z**2 - 920*y**2*z**4 + 976*z**6) + x**2*(-203*y**8 + 7560*y**6*z**2 - 22512*y**4*z**4 - 448*y**2*z**6 + 8832*z**8)) + 1j*-x*(7*x**10 - 84*y**10 + 4767*y**8*z**2 - 29232*y**6*z**4 + 40656*y**4*z**6 - 13056*y**2*z**8 + 384*z**10 - 7*x**8*(8*y**2 + 47*z**2) - 42*x**6*(7*y**4 - 90*y**2*z**2 - 32*z**4) - 14*x**4*(34*y**6 - 951*y**4*z**2 + 1896*y**2*z**4 + 8*z**6) + x**2*(-329*y**8 + 13972*y**6*z**2 - 57120*y**4*z**4 + 40544*y**2*z**6 - 1408*z**8)) )
            Bz = Bz * 39/4*rt3*( 1/2*(x**2 - y**2)*z*(21*x**8 + 21*y**8 - 336*y**6*z**2 + 1008*y**4*z**4 - 768*y**2*z**6 + 128*z**8 + 84*x**6*(y**2 - 4*z**2) + 126*x**4*(y**4 - 8*y**2*z**2 + 8*z**4) + 12*x**2*(7*y**6 - 84*y**4*z**2 + 168*y**2*z**4 - 64*z**6)) + 1j*x*y*z*(21*x**8 + 21*y**8 - 336*y**6*z**2 + 1008*y**4*z**4 - 768*y**2*z**6 + 128*z**8 + 84*x**6*(y**2 - 4*z**2) + 126*x**4*(y**4 - 8*y**2*z**2 + 8*z**4) + 12*x**2*(7*y**6 - 84*y**4*z**2 + 168*y**2*z**4 - 64*z**6)) )
        elif m==3:
            rt39o2 = sqrt(39/2)
            Bx = Bx * 3/2*rt39o2*( z*(28*x**10 - 7*x**8*(3*y**2 + 59*z**2) - 28*x**6*(8*y**4 - 23*y**2*z**2 - 39*z**4) - 2*x**4*(133*y**6 - 1225*y**4*z**2 + 1218*y**2*z**4 + 332*z**6) + y**2*(7*y**8 - 77*y**6*z**2 + 84*y**4*z**4 + 104*y**2*z**6 - 64*z**8) - 4*x**2*(21*y**8 - 329*y**6*z**2 + 861*y**4*z**4 - 420*y**2*z**6 - 16*z**8)) + 1j*x*y*z*(91*x**8 - 49*y**8 + 644*y**6*z**2 - 1344*y**4*z**4 + 352*y**2*z**6 + 128*z**8 + 28*x**6*(8*y**2 - 47*z**2) + 14*x**4*(9*y**4 - 142*y**2*z**2 + 240*z**4) - 4*x**2*(14*y**6 + 7*y**4*z**2 - 504*y**2*z**4 + 472*z**6)) )
            By = By * 3/2*rt39o2*( x*y*z*(49*x**8 - 91*y**8 + 1316*y**6*z**2 - 3360*y**4*z**4 + 1888*y**2*z**6 - 128*z**8 + 28*x**6*(2*y**2 - 23*z**2) - 14*x**4*(9*y**4 - 2*y**2*z**2 - 96*z**4) - 4*x**2*(56*y**6 - 497*y**4*z**2 + 504*y**2*z**4 + 88*z**6)) + 1j*z*(-7*x**10 - 28*y**10 + 413*y**8*z**2 - 1092*y**6*z**4 + 664*y**4*z**6 - 64*y**2*z**8 + x**8*(84*y**2 + 77*z**2) + 14*x**6*(19*y**4 - 94*y**2*z**2 - 6*z**4) + 2*x**4*(112*y**6 - 1225*y**4*z**2 + 1722*y**2*z**4 - 52*z**6) + x**2*(21*y**8 - 644*y**6*z**2 + 2436*y**4*z**4 - 1680*y**2*z**6 + 64*z**8)) )
            Bz = Bz * 7/2*rt39o2*( -x*(x**2 - 3*y**2)*(x**8 + y**8 - 56*y**6*z**2 + 336*y**4*z**4 - 448*y**2*z**6 + 128*z**8 + 4*x**6*(y**2 - 14*z**2) + 6*x**4*(y**4 - 28*y**2*z**2 + 56*z**4) + 4*x**2*(y**6 - 42*y**4*z**2 + 168*y**2*z**4 - 112*z**6)) + 1j*y*(-3*x**2 + y**2)*(x**8 + y**8 - 56*y**6*z**2 + 336*y**4*z**4 - 448*y**2*z**6 + 128*z**8 + 4*x**6*(y**2 - 14*z**2) + 6*x**4*(y**4 - 28*y**2*z**2 + 56*z**4) + 4*x**2*(y**6 - 42*y**4*z**2 + 168*y**2*z**4 - 112*z**6)) )
        elif m==4:
            rt39 = sqrt(39)
            Bx = Bx * -rt39*( 1/4*x*(11*x**10 - x**8*(49*y**2 + 556*z**2) - 6*x**6*(31*y**4 - 480*y**2*z**2 - 476*z**4) - 14*x**4*(11*y**6 - 444*y**4*z**2 + 1164*y**2*z**4 + 208*z**6) + 3*y**2*(9*y**8 - 404*y**6*z**2 + 1624*y**4*z**4 - 672*y**2*z**6 - 448*z**8) + x**2*(-y**8 + 1568*y**6*z**2 - 14280*y**4*z**4 + 16576*y**2*z**6 + 448*z**8)) + 1j*y*(12*x**10 + x**8*(19*y**2 - 597*z**2) - 14*x**6*(y**4 + 28*y**2*z**2 - 213*z**4) - 6*x**4*(6*y**6 - 161*y**4*z**2 + 105*y**2*z**4 + 476*z**6) + y**2*(y**8 - 41*y**6*z**2 + 126*y**4*z**4 + 56*y**2*z**6 - 112*z**8) + x**2*(-14*y**8 + 720*y**6*z**2 - 3486*y**4*z**4 + 2576*y**2*z**6 + 336*z**8)) )
            By = By * rt39*( -1/4*y*(27*x**10 - x**8*(y**2 + 1212*z**2) + x**6*(-154*y**4 + 1568*y**2*z**2 + 4872*z**4) - 6*x**4*(31*y**6 - 1036*y**4*z**2 + 2380*y**2*z**4 + 336*z**6) + x**2*(-49*y**8 + 2880*y**6*z**2 - 16296*y**4*z**4 + 16576*y**2*z**6 - 1344*z**8) + y**2*(11*y**8 - 556*y**6*z**2 + 2856*y**4*z**4 - 2912*y**2*z**6 + 448*z**8)) + 1j*x*(x**10 - x**8*(14*y**2 + 41*z**2) - 18*x**6*(2*y**4 - 40*y**2*z**2 - 7*z**4) - 14*x**4*(y**6 - 69*y**4*z**2 + 249*y**2*z**4 - 4*z**6) + x**2*(19*y**8 - 392*y**6*z**2 - 630*y**4*z**4 + 2576*y**2*z**6 - 112*z**8) + 3*y**2*(4*y**8 - 199*y**6*z**2 + 994*y**4*z**4 - 952*y**2*z**6 + 112*z**8)) )
            Bz = Bz * -21*rt39*( 1/4*(x**4 - 6*x**2*y**2 + y**4)*z*(5*x**6 + 5*y**6 - 70*y**4*z**2 + 168*y**2*z**4 - 80*z**6 + 5*x**4*(3*y**2 - 14*z**2) + x**2*(15*y**4 - 140*y**2*z**2 + 168*z**4)) + 1j*x*y*(x**2 - y**2)*z*(5*x**6 + 5*y**6 - 70*y**4*z**2 + 168*y**2*z**4 - 80*z**6 + 5*x**4*(3*y**2 - 14*z**2) + x**2*(15*y**4 - 140*y**2*z**2 + 168*z**4)) )
        elif m==5:
            rt39o10 = sqrt(39/10)
            Bx = Bx * 1/2*rt39o10*( z*(-180*x**10 + 5*x**8*(363*y**2 + 419*z**2) + 28*x**6*(90*y**4 - 785*y**2*z**2 - 131*z**4) - 210*x**4*(5*y**6 + 35*y**4*z**2 - 178*y**2*z**4 - 4*z**6) + 5*y**4*(15*y**6 - 125*y**4*z**2 + 28*y**2*z**4 + 168*z**6) - 20*x**2*(75*y**8 - 805*y**6*z**2 + 987*y**4*z**4 + 252*y**2*z**6)) + 1j*-x*y*z*(975*x**8 + 555*y**8 - 5220*y**6*z**2 + 4368*y**4*z**4 + 3360*y**2*z**6 - 300*x**6*(2*y**2 + 37*z**2) - 70*x**4*(51*y**4 - 230*y**2*z**2 - 264*z**4) - 20*x**2*(72*y**6 - 1099*y**4*z**2 + 1904*y**2*z**4 + 168*z**6)) )
            By = By * 1/2*rt39o10*( -x*y*z*(555*x**8 - 180*x**6*(8*y**2 + 29*z**2) + x**4*(-3570*y**4 + 21980*y**2*z**2 + 4368*z**4) - 20*x**2*(30*y**6 - 805*y**4*z**2 + 1904*y**2*z**4 - 168*z**6) + 15*(65*y**8 - 740*y**6*z**2 + 1232*y**4*z**4 - 224*y**2*z**6)) + 1j*z*(75*x**10 - 180*y**10 + 2095*y**8*z**2 - 3668*y**6*z**4 + 840*y**4*z**6 - 125*x**8*(12*y**2 + 5*z**2) - 70*x**6*(15*y**4 - 230*y**2*z**2 - 2*z**4) + 210*x**4*(12*y**6 - 35*y**4*z**2 - 94*y**2*z**4 + 4*z**6) + 5*x**2*(363*y**8 - 4396*y**6*z**2 + 7476*y**4*z**4 - 1008*y**2*z**6)) )
            Bz = Bz * 3/2*rt39o10*( x*(x**4 - 10*x**2*y**2 + 5*y**4)*(5*x**6 + 5*y**6 - 240*y**4*z**2 + 1120*y**2*z**4 - 896*z**6 + 15*x**4*(y**2 - 16*z**2) + 5*x**2*(3*y**4 - 96*y**2*z**2 + 224*z**4)) + 1j*y*(5*x**4 - 10*x**2*y**2 + y**4)*(5*x**6 + 5*y**6 - 240*y**4*z**2 + 1120*y**2*z**4 - 896*z**6 + 15*x**4*(y**2 - 16*z**2) + 5*x**2*(3*y**4 - 96*y**2*z**2 + 224*z**4)) )
        elif m==6:
            rt39o2 = sqrt(39/2)
            Bx = Bx * 3/4*rt39o2*( 1/2*x*(11*x**10 - x**8*(179*y**2 + 426*z**2) + x**6*(-106*y**4 + 7080*y**2*z**2 + 1376*z**4) + 14*x**4*(23*y**6 - 186*y**4*z**2 - 1584*y**2*z**4 - 32*z**6) - y**4*(47*y**6 - 1542*y**4*z**2 + 2912*y**2*z**4 + 2240*z**6) + x**2*(191*y**8 - 8568*y**6*z**2 + 25760*y**4*z**4 + 4480*y**2*z**6)) + 1j*-y*(-36*x**10 + x**8*(83*y**2 + 1371*z**2) + 28*x**6*(7*y**4 - 153*y**2*z**2 - 152*z**4) + y**4*(3*y**6 - 93*y**4*z**2 + 128*y**2*z**4 + 224*z**6) + 2*x**4*(y**6 - 1491*y**4*z**2 + 7840*y**2*z**4 + 560*z**6) - 4*x**2*(18*y**8 - 645*y**6*z**2 + 1512*y**4*z**4 + 560*y**2*z**6)) )
            By = By * 3/4*rt39o2*( 1/2*y*(47*x**10 - 11*y**10 + 426*y**8*z**2 - 1376*y**6*z**4 + 448*y**4*z**6 - x**8*(191*y**2 + 1542*z**2) + x**6*(-322*y**4 + 8568*y**2*z**2 + 2912*z**4) + 2*x**4*(53*y**6 + 1302*y**4*z**2 - 12880*y**2*z**4 + 1120*z**6) + x**2*(179*y**8 - 7080*y**6*z**2 + 22176*y**4*z**4 - 4480*y**2*z**6)) + 1j*-x*(3*x**10 - 36*y**10 + 1371*y**8*z**2 - 4256*y**6*z**4 + 1120*y**4*z**6 - 3*x**8*(24*y**2 + 31*z**2) + 2*x**6*(y**4 + 1290*y**2*z**2 + 64*z**4) + 14*x**4*(14*y**6 - 213*y**4*z**2 - 432*y**2*z**4 + 16*z**6) + x**2*(83*y**8 - 4284*y**6*z**2 + 15680*y**4*z**4 - 2240*y**2*z**6)) )
            Bz = Bz * 17/4*rt39o2*( 1/2*(x**6 - 15*x**4*y**2 + 15*x**2*y**4 - y**6)*z*(15*x**4 + 15*y**4 - 160*y**2*z**2 + 224*z**4 + 10*x**2*(3*y**2 - 16*z**2)) + 1j*x*y*(3*x**4 - 10*x**2*y**2 + 3*y**4)*z*(15*x**4 + 15*y**4 - 160*y**2*z**2 + 224*z**4 + 10*x**2*(3*y**2 - 16*z**2)) )
        elif m==7:
            rt663o2 = sqrt(663/2)
            Bx = Bx * 1/4*rt663o2*( z*(36*x**10 - x**8*(867*y**2 + 251*z**2) + 7*y**6*(3*y**4 - 13*y**2*z**2 - 16*z**4) + 28*x**6*(39*y**4 + 211*y**2*z**2 + 4*z**4) + 42*x**4*(31*y**6 - 285*y**4*z**2 - 40*y**2*z**4) - 28*x**2*(24*y**8 - 131*y**6*z**2 - 60*y**4*z**4)) + 1j*x*y*z*(273*x**8 - 183*y**8 + 888*y**6*z**2 + 672*y**4*z**4 - 84*x**6*(17*y**2 + 22*z**2) - 14*x**4*(15*y**4 - 764*y**2*z**2 - 48*z**4) + 4*x**2*(327*y**6 - 2114*y**4*z**2 - 560*y**2*z**4)) )
            By = By * 1/4*rt663o2*( x*y*z*(183*x**8 - 12*x**6*(109*y**2 + 74*z**2) + 14*x**4*(15*y**4 + 604*y**2*z**2 - 48*z**4) + 28*x**2*(51*y**6 - 382*y**4*z**2 + 80*y**2*z**4) - 21*(13*y**8 - 88*y**6*z**2 + 32*y**4*z**4)) + 1j*z*(-21*x**10 - 36*y**10 + 251*y**8*z**2 - 112*y**6*z**4 + 7*x**8*(96*y**2 + 13*z**2) - 14*x**6*(93*y**4 + 262*y**2*z**2 - 8*z**4) - 42*x**4*(26*y**6 - 285*y**4*z**2 + 40*y**2*z**4) + x**2*(867*y**8 - 5908*y**6*z**2 + 1680*y**4*z**4)) )
            Bz = Bz * 3/4*rt663o2*( -x*(x**6 - 21*x**4*y**2 + 35*x**2*y**4 - 7*y**6)*(x**4 + y**4 - 36*y**2*z**2 + 96*z**4 + 2*x**2*(y**2 - 18*z**2)) + 1j*y*(-7*x**6 + 35*x**4*y**2 - 21*x**2*y**4 + y**6)*(x**4 + y**4 - 36*y**2*z**2 + 96*z**4 + 2*x**2*(y**2 - 18*z**2)) )
        elif m==8:
            rt221 = sqrt(221)
            Bx = Bx * rt221*( 1/8*x*(-11*x**10 + x**8*(361*y**2 + 244*z**2) - 18*x**6*(47*y**4 + 440*y**2*z**2 + 8*z**4) + 3*y**6*(-25*y**4 + 444*y**2*z**2 + 336*z**4) - 42*x**4*(11*y**6 - 588*y**4*z**2 - 72*y**2*z**4) + 3*x**2*(227*y**8 - 4816*y**6*z**2 - 1680*y**4*z**4)) + 1j*-y*(12*x**10 - 3*x**8*(31*y**2 + 87*z**2) + y**6*(y**4 - 17*y**2*z**2 - 18*z**4) + 42*x**6*(y**4 + 52*y**2*z**2 + 3*z**4) + 18*x**4*(6*y**6 - 161*y**4*z**2 - 35*y**2*z**4) + x**2*(-38*y**8 + 720*y**6*z**2 + 378*y**4*z**4)) )
            By = By * rt221*( 1/8*y*(-75*x**10 - 11*y**10 + 244*y**8*z**2 - 144*y**6*z**4 + 3*x**8*(227*y**2 + 444*z**2) - 42*x**6*(11*y**4 + 344*y**2*z**2 - 24*z**4) - 18*x**4*(47*y**6 - 1372*y**4*z**2 + 280*y**2*z**4) + x**2*(361*y**8 - 7920*y**6*z**2 + 3024*y**4*z**4)) + 1j*x*(x**10 - x**8*(38*y**2 + 17*z**2) + 18*x**6*(6*y**4 + 40*y**2*z**2 - z**4) + 3*y**6*(4*y**4 - 87*y**2*z**2 + 42*z**4) + 42*x**4*(y**6 - 69*y**4*z**2 + 9*y**2*z**4) - 3*x**2*(31*y**8 - 728*y**6*z**2 + 210*y**4*z**4)) )
            Bz = Bz * -57*rt221*( 1/8*(x**8 - 28*x**6*y**2 + 70*x**4*y**4 - 28*x**2*y**6 + y**8)*z*(x**2 + y**2 - 6*z**2) + 1j*x*y*(x**6 - 7*x**4*y**2 + 7*x**2*y**4 - y**6)*z*(x**2 + y**2 - 6*z**2) )
        elif m==9:
            rt4199o2 = sqrt(4199/2)
            Bx = Bx * 3/4*rt4199o2*( z*(-4*x**10 + 3*y**8*(y**2 + z**2) - 84*x**6*y**2*(9*y**2 + z**2) + 3*x**8*(57*y**2 + z**2) + 42*x**4*(17*y**6 + 5*y**4*z**2) - 12*x**2*(12*y**8 + 7*y**6*z**2)) + 1j*-x*y*z*(39*x**8 + 31*y**8 + 24*y**6*z**2 - 12*x**6*(37*y**2 + 2*z**2) + 42*x**4*(21*y**4 + 4*y**2*z**2) - 12*x**2*(33*y**6 + 14*y**4*z**2)) )
            By = By * 3/4*rt4199o2*( -x*y*z*(31*x**8 + 39*y**8 - 24*y**6*z**2 + x**6*(-396*y**2 + 24*z**2) + 42*x**4*(21*y**4 - 4*y**2*z**2) + x**2*(-444*y**6 + 168*y**4*z**2)) + 1j*z*(3*x**10 - 4*y**10 + 3*y**8*z**2 + 3*x**8*(-48*y**2 + z**2) + 42*x**6*(17*y**4 - 2*y**2*z**2) - 42*x**4*(18*y**6 - 5*y**4*z**2) + 3*x**2*(57*y**8 - 28*y**6*z**2)) )
            Bz = Bz * 1/4*rt4199o2*( x*(x**8 - 36*x**6*y**2 + 126*x**4*y**4 - 84*x**2*y**6 + 9*y**8)*(x**2 + y**2 - 20*z**2) + 1j*y*(9*x**8 - 84*x**6*y**2 + 126*x**4*y**4 - 36*x**2*y**6 + y**8)*(x**2 + y**2 - 20*z**2) )
        elif m==10:
            rt4199o10 = sqrt(4199/10)
            Bx = Bx * 1/4*rt4199o10*( 1/2*x*(11*x**10 - 5*x**8*(119*y**2 + 2*z**2) - 3*y**8*(37*y**2 + 30*z**2) + 90*x**6*(39*y**4 + 4*y**2*z**2) - 210*x**4*(23*y**6 + 6*y**4*z**2) + 15*x**2*(113*y**8 + 56*y**6*z**2)) + 1j*y*(60*x**10 - 5*y**8*(y**2 + z**2) - 15*x**8*(59*y**2 + 3*z**2) + 84*x**6*(29*y**4 + 5*y**2*z**2) - 90*x**4*(19*y**6 + 7*y**4*z**2) + 20*x**2*(14*y**8 + 9*y**6*z**2)) )
            By = By * 1/4*rt4199o10*( 1/2*y*(111*x**10 - 11*y**10 + 10*y**8*z**2 + x**8*(-1695*y**2 + 90*z**2) + 210*x**6*(23*y**4 - 4*y**2*z**2) - 90*x**4*(39*y**6 - 14*y**4*z**2) + 5*x**2*(119*y**8 - 72*y**6*z**2)) + 1j*-x*(5*x**10 - 60*y**10 + 45*y**8*z**2 + 5*x**8*(-56*y**2 + z**2) + 90*x**6*(19*y**4 - 2*y**2*z**2) - 42*x**4*(58*y**6 - 15*y**4*z**2) + 15*x**2*(59*y**8 - 28*y**6*z**2)) )
            Bz = Bz * 21/4*rt4199o10*( 1/2*(x**10 - 45*x**8*y**2 + 210*x**6*y**4 - 210*x**4*y**6 + 45*x**2*y**8 - y**10)*z + 1j*x*y*(5*x**8 - 60*x**6*y**2 + 126*x**4*y**4 - 60*x**2*y**6 + 5*y**8)*z )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi, n={n} and m is not between -n and n.")

    else:
        print(" n = ", n)
        raise ValueError("In field_xyz.eval_Bi, n>10 but only n=1 to n=10 are supported.")

    Bx = Bx * timeRot
    By = By * timeRot
    Bz = Bz * timeRot

    return Bx, By, Bz
    
#############################################

"""
eval_Bi_Schmidt()
    Same as eval_Bi but using Schmidt semi-normalized harmonics, with particular magnetic moments g_nm, h_nm.
    Usage: `Bx`, `By`, `Bz` = eval_Bi_Schmidt(`n`, `m`, `g_nm`, `h_nm`, `x`, `y`, `z`, `r`, `omega=None`, `t=None`)
    Returns:
        Bx, By, Bz (each): complex, ndarray shape(Nvals). A linear array of field values due to these particular n,m values,
            at each of the specified points. Returns a time sequence if omega and t are passed.
    Parameters:
        n: integer. Degree of magnetic moment to be evaluated.
        m: integer. Order of magnetic moment to be evaluated.
        g_nm, h_nm: complex. Multipole Gauss coefficient of degree and order n,m. Units match the output field.
        x,y,z,r: float, shape(Nvals). Linear arrays of corresponding x,y, and z values. r is redundant but
            saves on computation time to avoid recalculating on every call to this function. If omega and t
            are passed, these quantities are the trajectory locations.
        omega: float (None). Optional oscillation frequency in rads/s for evaluating time series. Requires t to be passed as well.
        t: float, shape(Nvals) (None). Optional time values in TDB seconds since J2000 epoch. Required if omega is passed.
    """
def eval_Bi_Schmidt(n,m,g_nm,h_nm, x,y,z,r, omega=None, t=None):

    if omega is None:
        timeRot = 1.0
    else:
        timeRot = np.exp(-1j*omega*t)

    if n==1:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   Dipole field components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        Bx = 1/r**5 + 0j
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==0:
            Bx = Bx * g_nm * (3*x*z)
            By = By * g_nm * (3*y*z)
            Bz = Bz * g_nm * (-x**2 - y**2 + 2*z**2)
        elif m==1:
            Bx = Bx * ( g_nm*(2*x**2 - y**2 - z**2) + h_nm*(3*x*y) )
            By = By * ( g_nm*(3*x*y) + h_nm*(-x**2 + 2*y**2 - z**2) )
            Bz = Bz * ( g_nm*(3*x*z) + h_nm*(3*y*z) )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi_Schmidt, n={n} and m is not between 0 and n.")

    elif n==2:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   Quadrupole moment components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A2 = sqrt(3/4)
        Bx = A2/r**7 + 0j
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==0:
            rt3 = sqrt(3)
            Bx = Bx * g_nm * -rt3*x*(x**2 + y**2 - 4*z**2)
            By = By * g_nm * -rt3*y*(x**2 + y**2 - 4*z**2)
            Bz = Bz * g_nm * rt3*z*(-3*x**2 - 3*y**2 + 2*z**2)
        elif m==1:
            Bx = Bx * ( g_nm*-2*z*(-4*x**2 + y**2 + z**2) + h_nm*(10*x*y*z) )
            By = By * ( g_nm*(10*x*y*z) + h_nm*-2*z*(x**2 - 4*y**2 + z**2) )
            Bz = Bz * ( g_nm*-2*x*(x**2 + y**2 - 4*z**2) + h_nm*-2*y*(x**2 + y**2 - 4*z**2) )
        elif m==2:
            Bx = Bx * ( g_nm*x*(3*x**2 - 7*y**2 - 2*z**2) + h_nm*-2*y*(-4*x**2 + y**2 + z**2) )
            By = By * ( g_nm*y*(7*x**2 - 3*y**2 + 2*z**2) + h_nm*-2*x*(x**2 - 4*y**2 + z**2) )
            Bz = Bz * ( g_nm*5*(x**2 - y**2)*z + h_nm*(10*x*y*z) )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi_Schmidt, n={n} and m is not between 0 and n.")

    elif n==3:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   Octupole moment components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A3 = sqrt(3/8)
        Bx = A3/r**9 + 0j
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==0:
            rt2o3 = sqrt(2/3)
            Bx = Bx * g_nm * -5*rt2o3*x*z*(3*x**2 + 3*y**2 - 4*z**2)
            By = By * g_nm * -5*rt2o3*y*z*(3*x**2 + 3*y**2 - 4*z**2)
            Bz = Bz * g_nm * rt2o3*(3*x**4 + 3*y**4 - 24*y**2*z**2 + 8*z**4 + 6*x**2*(y**2 - 4*z**2))
        elif m==1:
            Bx = Bx * ( g_nm*(-4*x**4 + y**4 - 3*y**2*z**2 - 4*z**4 - 3*x**2*(y**2 - 9*z**2)) + h_nm*-5*x*y*(x**2 + y**2 - 6*z**2) )
            By = By * ( g_nm*-5*x*y*(x**2 + y**2 - 6*z**2) + h_nm*(x**4 - 4*y**4 + 27*y**2*z**2 - 4*z**4 - 3*x**2*(y**2 + z**2)) )
            Bz = Bz * -5*( g_nm*x*z*(3*x**2 + 3*y**2 - 4*z**2) + h_nm*y*z*(3*x**2 + 3*y**2 - 4*z**2) )
        elif m==2:
            rt10 = sqrt(10)
            Bx = Bx * rt10*( g_nm*x*z*(5*x**2 - 9*y**2 - 2*z**2) + h_nm*-2*y*z*(-6*x**2 + y**2 + z**2) )
            By = By * rt10*( g_nm*y*z*(9*x**2 - 5*y**2 + 2*z**2) + h_nm*-2*x*z*(x**2 - 6*y**2 + z**2) )
            Bz = Bz * -rt10*( g_nm*(x**2 - y**2)*(x**2 + y**2 - 6*z**2) + h_nm*2*x*y*(x**2 + y**2 - 6*z**2) )
        elif m==3:
            rt5o3 = sqrt(5/3)
            Bx = Bx * rt5o3*( g_nm*(4*x**4 + 3*y**2*(y**2 + z**2) - 3*x**2*(7*y**2 + z**2)) + h_nm*x*y*(15*x**2 - 13*y**2 - 6*z**2) )
            By = By * rt5o3*( g_nm*x*y*(13*x**2 - 15*y**2 + 6*z**2) + h_nm*(-3*x**4 - 4*y**4 + 3*y**2*z**2 + 3*x**2*(7*y**2 - z**2)) )
            Bz = Bz * 7*rt5o3*( g_nm*x*(x**2 - 3*y**2)*z + h_nm*-y*(-3*x**2 + y**2)*z )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi_Schmidt, n={n} and m is not between 0 and n.")

    elif n==4:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   Hexadecapole moment components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A4 = sqrt(5/32)
        Bx = A4/r**11 + 0j
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==0:
            rt5o2 = sqrt(5/2)
            Bx = Bx * g_nm * 3*rt5o2*x*(x**4 + y**4 - 12*y**2*z**2 + 8*z**4 + 2*x**2*(y**2 - 6*z**2))
            By = By * g_nm * 3*rt5o2*y*(x**4 + y**4 - 12*y**2*z**2 + 8*z**4 + 2*x**2*(y**2 - 6*z**2))
            Bz = Bz * g_nm * rt5o2*z*(15*x**4 + 15*y**4 - 40*y**2*z**2 + 8*z**4 + 10*x**2*(3*y**2 - 4*z**2))
        elif m==1:
            Bx = Bx * ( g_nm*-2*z*(18*x**4 - 3*y**4 + y**2*z**2 + 4*z**4 + x**2*(15*y**2 - 41*z**2)) + h_nm*-42*x*y*z*(x**2 + y**2 - 2*z**2) )
            By = By * ( g_nm*-42*x*y*z*(x**2 + y**2 - 2*z**2) + h_nm*2*z*(3*x**4 - 18*y**4 + 41*y**2*z**2 - 4*z**4 - x**2*(15*y**2 + z**2)) )
            Bz = Bz * 6*( g_nm*x*(x**4 + y**4 - 12*y**2*z**2 + 8*z**4 + 2*x**2*(y**2 - 6*z**2)) + h_nm*y*(x**4 + y**4 - 12*y**2*z**2 + 8*z**4 + 2*x**2*(y**2 - 6*z**2)) )
        elif m==2:
            rt2 = sqrt(2)
            Bx = Bx * rt2*( g_nm*-x*(5*x**4 - 9*y**4 + 66*y**2*z**2 + 12*z**4 - 2*x**2*(2*y**2 + 23*z**2)) + h_nm*2*y*(-6*x**4 + y**4 - 5*y**2*z**2 - 6*z**4 + x**2*(-5*y**2 + 51*z**2)) )
            By = By * rt2*( g_nm*y*(-9*x**4 + 5*y**4 - 46*y**2*z**2 + 12*z**4 + x**2*(-4*y**2 + 66*z**2)) + h_nm*2*x*(x**4 - 6*y**4 + 51*y**2*z**2 - 6*z**4 - 5*x**2*(y**2 + z**2)) )
            Bz = Bz * -21*rt2*( g_nm*(x**2 - y**2)*z*(x**2 + y**2 - 2*z**2) + h_nm*2*x*y*z*(x**2 + y**2 - 2*z**2) )
        elif m==3:
            rt7 = sqrt(7)
            Bx = Bx * 6*rt7*( g_nm*z*(2*x**4 + y**2*(y**2 + z**2) - x**2*(9*y**2 + z**2)) + h_nm*x*y*z*(7*x**2 - 5*y**2 - 2*z**2) )
            By = By * 6*rt7*( g_nm*x*y*z*(5*x**2 - 7*y**2 + 2*z**2) + h_nm*-z*(x**4 + 2*y**4 - y**2*z**2 + x**2*(-9*y**2 + z**2)) )
            Bz = Bz * 2*rt7*( g_nm*-x*(x**2 - 3*y**2)*(x**2 + y**2 - 8*z**2) + h_nm*y*(-3*x**2 + y**2)*(x**2 + y**2 - 8*z**2) )
        elif m==4:
            rt2 = sqrt(2)
            rt7 = sqrt(7)
            Bx = Bx * rt7*( g_nm/rt2*x*(5*x**4 - 2*x**2*(23*y**2 + 2*z**2) + 3*y**2*(7*y**2 + 4*z**2)) + h_nm*2*rt2*y*(6*x**4 + y**2*(y**2 + z**2) - x**2*(11*y**2 + 3*z**2)) )
            By = By * rt7*( g_nm/rt2*y*(21*x**4 + 5*y**4 - 4*y**2*z**2 + x**2*(-46*y**2 + 12*z**2)) + h_nm*-2*rt2*x*(x**4 + 6*y**4 - 3*y**2*z**2 + x**2*(-11*y**2 + z**2)) )
            Bz = Bz * 9*rt7*( g_nm/rt2*(x**4 - 6*x**2*y**2 + y**4)*z + h_nm*2*rt2*x*y*(x**2 - y**2)*z )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi_Schmidt, n={n} and m is not between 0 and n.")

    elif n==5:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   n=5 moment components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A5 = sqrt(15/16)
        Bx = A5/r**13 + 0j
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==0:
            rt3o5 = sqrt(3/5)
            Bx = Bx * g_nm * 7/2*rt3o5*x*z*(5*x**4 + 5*y**4 - 20*y**2*z**2 + 8*z**4 + 10*x**2*(y**2 - 2*z**2))
            By = By * g_nm * 7/2*rt3o5*y*z*(5*x**4 + 5*y**4 - 20*y**2*z**2 + 8*z**4 + 10*x**2*(y**2 - 2*z**2))
            Bz = Bz * g_nm * -1/2*rt3o5*(5*x**6 + 5*y**6 - 90*y**4*z**2 + 120*y**2*z**4 - 16*z**6 + 15*x**4*(y**2 - 6*z**2) + 15*x**2*(y**4 - 12*y**2*z**2 + 8*z**4))
        elif m==1:
            Bx = Bx * 1/2*( g_nm*(6*x**6 - y**6 + 11*y**4*z**2 + 4*y**2*z**4 - 8*z**6 + x**4*(11*y**2 - 101*z**2) + 2*x**2*(2*y**4 - 45*y**2*z**2 + 58*z**4)) + h_nm*7*x*y*(x**4 + y**4 - 16*y**2*z**2 + 16*z**4 + 2*x**2*(y**2 - 8*z**2)) )
            By = By * 1/2*( g_nm*7*x*y*(x**4 + y**4 - 16*y**2*z**2 + 16*z**4 + 2*x**2*(y**2 - 8*z**2)) + h_nm*(-x**6 + 6*y**6 - 101*y**4*z**2 + 116*y**2*z**4 - 8*z**6 + x**4*(4*y**2 + 11*z**2) + x**2*(11*y**4 - 90*y**2*z**2 + 4*z**4)) )
            Bz = Bz * 7/2*( g_nm*x*z*(5*x**4 + 5*y**4 - 20*y**2*z**2 + 8*z**4 + 10*x**2*(y**2 - 2*z**2)) + h_nm*y*z*(5*x**4 + 5*y**4 - 20*y**2*z**2 + 8*z**4 + 10*x**2*(y**2 - 2*z**2)) )
        elif m==2:
            rt7 = sqrt(7)
            Bx = Bx * rt7*( g_nm*x*z*(-7*x**4 + 11*y**4 - 26*y**2*z**2 - 4*z**4 + x**2*(4*y**2 + 22*z**2)) + h_nm*-2*y*z*(8*x**4 - y**4 + y**2*z**2 + 2*z**4 + x**2*(7*y**2 - 23*z**2)) )
            By = By * rt7*( g_nm*y*z*(-11*x**4 + 7*y**4 - 22*y**2*z**2 + 4*z**4 + x**2*(-4*y**2 + 26*z**2)) + h_nm*2*x*z*(x**4 - 8*y**4 + 23*y**2*z**2 - 2*z**4 - x**2*(7*y**2 + z**2)) )
            Bz = Bz * rt7*( g_nm*(x**2 - y**2)*(x**4 + y**4 - 16*y**2*z**2 + 16*z**4 + 2*x**2*(y**2 - 8*z**2)) + h_nm*2*x*y*(x**4 + y**4 - 16*y**2*z**2 + 16*z**4 + 2*x**2*(y**2 - 8*z**2)) )
        elif m==3:
            rt21o2 = sqrt(21/2)
            Bx = Bx * -1/2*rt21o2*( g_nm*(2*x**6 + y**6 - 7*y**4*z**2 - 8*y**2*z**4 - x**4*(7*y**2 + 23*z**2) + x**2*(-8*y**4 + 90*y**2*z**2 + 8*z**4)) + h_nm*x*y*(7*x**4 - 5*y**4 + 44*y**2*z**2 + 16*z**4 + 2*x**2*(y**2 - 38*z**2)) )
            By = By * 1/2*rt21o2*( g_nm*x*y*(-5*x**4 + 7*y**4 - 76*y**2*z**2 + 16*z**4 + 2*x**2*(y**2 + 22*z**2)) + h_nm*(x**6 + 2*y**6 - 23*y**4*z**2 + 8*y**2*z**4 - x**4*(8*y**2 + 7*z**2) + x**2*(-7*y**4 + 90*y**2*z**2 - 8*z**4)) )
            Bz = Bz * 3/2*rt21o2*( g_nm*-x*(x**2 - 3*y**2)*z*(3*x**2 + 3*y**2 - 8*z**2) + h_nm*y*(-3*x**2 + y**2)*z*(3*x**2 + 3*y**2 - 8*z**2) )
        elif m==4:
            rt21 = sqrt(21)
            Bx = Bx * rt21*( g_nm*1/2*x*z*(7*x**4 + 23*y**4 + 12*y**2*z**2 - 2*x**2*(29*y**2 + 2*z**2)) + h_nm*2*y*z*(8*x**4 + y**2*(y**2 + z**2) - x**2*(13*y**2 + 3*z**2)) )
            By = By * rt21*( g_nm*1/2*y*z*(23*x**4 + 7*y**4 - 4*y**2*z**2 + x**2*(-58*y**2 + 12*z**2)) + h_nm*-2*x*z*(x**4 + 8*y**4 - 3*y**2*z**2 + x**2*(-13*y**2 + z**2)) )
            Bz = Bz * rt21*( g_nm*-1/2*(x**4 - 6*x**2*y**2 + y**4)*(x**2 + y**2 - 10*z**2) + h_nm*-2*x*y*(x**2 - y**2)*(x**2 + y**2 - 10*z**2) )
        elif m==5:
            rt21o10 = sqrt(21/10)
            Bx = Bx * 1/2*rt21o10*( g_nm*(6*x**6 - 5*y**4*(y**2 + z**2) - 5*x**4*(17*y**2 + z**2) + 10*x**2*(8*y**4 + 3*y**2*z**2)) + h_nm*x*y*(35*x**4 + 31*y**4 + 20*y**2*z**2 - 10*x**2*(11*y**2 + 2*z**2)) )
            By = By * 1/2*rt21o10*( g_nm*x*y*(31*x**4 + 5*y**2*(7*y**2 - 4*z**2) + x**2*(-110*y**2 + 20*z**2)) + h_nm*-1*(5*x**6 - 6*y**6 + 5*y**4*z**2 + x**4*(-80*y**2 + 5*z**2) + 5*x**2*(17*y**4 - 6*y**2*z**2)) )
            Bz = Bz * 11/2*rt21o10*( g_nm*x*(x**4 - 10*x**2*y**2 + 5*y**4)*z + h_nm*y*(5*x**4 - 10*x**2*y**2 + y**4)*z )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi_Schmidt, n={n} and m is not between 0 and n.")

    elif n==6:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   n=6 moment components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A6 = sqrt(105/128)
        Bx = A6/r**15 + 0j
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==0:
            rt7o30 = sqrt(7/30)
            Bx = Bx * g_nm * -rt7o30*x*(5*x**6 + 5*y**6 - 120*y**4*z**2 + 240*y**2*z**4 - 64*z**6 + 15*x**4*(y**2 - 8*z**2) + 15*x**2*(y**4 - 16*y**2*z**2 + 16*z**4))
            By = By * g_nm * -rt7o30*y*(5*x**6 + 5*y**6 - 120*y**4*z**2 + 240*y**2*z**4 - 64*z**6 + 15*x**4*(y**2 - 8*z**2) + 15*x**2*(y**4 - 16*y**2*z**2 + 16*z**4))
            Bz = Bz * g_nm * -rt7o30*z*(35*x**6 + 35*y**6 - 210*y**4*z**2 + 168*y**2*z**4 - 16*z**6 + 105*x**4*(y**2 - 2*z**2) + 21*x**2*(5*y**4 - 20*y**2*z**2 + 8*z**4))
        elif m==1:
            rt2o5 = sqrt(2/5)
            Bx = Bx * rt2o5*( g_nm*z*(40*x**6 - 5*y**6 + 15*y**4*z**2 + 12*y**2*z**4 - 8*z**6 + 75*x**4*(y**2 - 3*z**2) + 6*x**2*(5*y**4 - 35*y**2*z**2 + 26*z**4)) + h_nm*3*x*y*z*(15*x**4 + 15*y**4 - 80*y**2*z**2 + 48*z**4 + 10*x**2*(3*y**2 - 8*z**2)) )
            By = By * rt2o5*( g_nm*3*x*y*z*(15*x**4 + 15*y**4 - 80*y**2*z**2 + 48*z**4 + 10*x**2*(3*y**2 - 8*z**2)) + h_nm*-z*(5*x**6 - 40*y**6 + 225*y**4*z**2 - 156*y**2*z**4 + 8*z**6 - 15*x**4*(2*y**2 + z**2) - 3*x**2*(25*y**4 - 70*y**2*z**2 + 4*z**4)) )
            Bz = Bz * -rt2o5*( g_nm*x*(5*x**6 + 5*y**6 - 120*y**4*z**2 + 240*y**2*z**4 - 64*z**6 + 15*x**4*(y**2 - 8*z**2) + 15*x**2*(y**4 - 16*y**2*z**2 + 16*z**4)) + h_nm*y*(5*x**6 + 5*y**6 - 120*y**4*z**2 + 240*y**2*z**4 - 64*z**6 + 15*x**4*(y**2 - 8*z**2) + 15*x**2*(y**4 - 16*y**2*z**2 + 16*z**4)) )
        elif m==2:
            Bx = Bx * ( g_nm*1/2*x*(7*x**6 - 11*y**6 + 210*y**4*z**2 - 240*y**2*z**4 - 32*z**6 + 3*x**4*(y**2 - 50*z**2) - 15*x**2*(y**4 - 4*y**2*z**2 - 16*z**4)) + h_nm*y*(8*x**6 - y**6 + 15*y**4*z**2 - 16*z**6 + 15*x**4*(y**2 - 11*z**2) + 6*x**2*(y**4 - 25*y**2*z**2 + 40*z**4)) )
            By = By * ( g_nm*1/2*y*(11*x**6 - 7*y**6 + 150*y**4*z**2 - 240*y**2*z**4 + 32*z**6 + 15*x**4*(y**2 - 14*z**2) - 3*x**2*(y**4 + 20*y**2*z**2 - 80*z**4)) + h_nm*-x*(x**6 - 8*y**6 + 165*y**4*z**2 - 240*y**2*z**4 + 16*z**6 - 3*x**4*(2*y**2 + 5*z**2) - 15*x**2*(y**4 - 10*y**2*z**2)) )
            Bz = Bz * 3*( g_nm*1/2*(x**2 - y**2)*z*(15*x**4 + 15*y**4 - 80*y**2*z**2 + 48*z**4 + 10*x**2*(3*y**2 - 8*z**2)) + h_nm*x*y*z*(15*x**4 + 15*y**4 - 80*y**2*z**2 + 48*z**4 + 10*x**2*(3*y**2 - 8*z**2)) )
        elif m==3:
            Bx = Bx * ( g_nm*z*(-24*x**6 - 9*y**6 + 15*y**4*z**2 + 24*y**2*z**4 + x**4*(75*y**2 + 95*z**2) + 6*x**2*(15*y**4 - 55*y**2*z**2 - 4*z**4)) + h_nm*-x*y*z*(81*x**4 - 51*y**4 + 140*y**2*z**2 + 48*z**4 + 30*x**2*(y**2 - 10*z**2)) )
            By = By * ( g_nm*x*y*z*(-51*x**4 + 81*y**4 - 300*y**2*z**2 + 48*z**4 + 10*x**2*(3*y**2 + 14*z**2)) + h_nm*z*(9*x**6 + 24*y**6 - 95*y**4*z**2 + 24*y**2*z**4 - 15*x**4*(6*y**2 + z**2) - 3*x**2*(25*y**4 - 110*y**2*z**2 + 8*z**4)) )
            Bz = Bz * ( g_nm*x*(x**2 - 3*y**2)*(3*x**4 + 3*y**4 - 60*y**2*z**2 + 80*z**4 + 6*x**2*(y**2 - 10*z**2)) + h_nm*y*(3*x**2 - y**2)*(3*x**4 + 3*y**4 - 60*y**2*z**2 + 80*z**4 + 6*x**2*(y**2 - 10*z**2)) )
        elif m==4:
            rt3o10 = sqrt(3/10)
            Bx = Bx * -rt3o10*( g_nm*x*(7*x**6 + 23*y**6 - 240*y**4*z**2 - 120*y**2*z**4 - 3*x**4*(17*y**2 + 32*z**2) + x**2*(-35*y**4 + 720*y**2*z**2 + 40*z**4)) + h_nm*4*y*(8*x**6 + y**6 - 9*y**4*z**2 - 10*y**2*z**4 - 5*x**4*(y**2 + 21*z**2) - 6*x**2*(2*y**4 - 25*y**2*z**2 - 5*z**4)) )
            By = By * rt3o10*( g_nm*-y*(23*x**6 + 7*y**6 - 96*y**4*z**2 + 40*y**2*z**4 - 5*x**4*(7*y**2 + 48*z**2) - 3*x**2*(17*y**4 - 240*y**2*z**2 + 40*z**4)) + h_nm*4*x*(x**6 + 8*y**6 - 105*y**4*z**2 + 30*y**2*z**4 - 3*x**4*(4*y**2 + 3*z**2) - 5*x**2*(y**4 - 30*y**2*z**2 + 2*z**4)) )
            Bz = Bz * -11*rt3o10*( g_nm*(x**4 - 6*x**2*y**2 + y**4)*z*(3*x**2 + 3*y**2 - 10*z**2) + h_nm*4*x*y*(x**2 - y**2)*z*(3*x**2 + 3*y**2 - 10*z**2) )
        elif m==5:
            rt33o5 = sqrt(33/5)
            Bx = Bx * rt33o5*( g_nm*z*(8*x**6 - 5*y**4*(y**2 + z**2) + 30*x**2*y**2*(3*y**2 + z**2) - 5*x**4*(21*y**2 + z**2)) + h_nm*x*y*z*(45*x**4 + 33*y**4 + 20*y**2*z**2 - 10*x**2*(13*y**2 + 2*z**2)) )
            By = By * rt33o5*( g_nm*x*y*z*(33*x**4 + 5*y**2*(9*y**2 - 4*z**2) + x**2*(-130*y**2 + 20*z**2)) + h_nm*z*(-5*x**6 + 8*y**6 - 5*y**4*z**2 + x**4*(90*y**2 - 5*z**2) - 15*x**2*(7*y**4 - 2*y**2*z**2)) )
            Bz = Bz * -rt33o5*( g_nm*x*(x**4 - 10*x**2*y**2 + 5*y**4)*(x**2 + y**2 - 12*z**2) + h_nm*y*(5*x**4 - 10*x**2*y**2 + y**4)*(x**2 + y**2 - 12*z**2) )
        elif m==6:
            rt11o5 = sqrt(11/5)
            Bx = Bx * rt11o5*( g_nm*1/2*x*(7*x**6 - 43*y**6 - 30*y**4*z**2 - 3*x**4*(47*y**2 + 2*z**2) + 15*x**2*(15*y**4 + 4*y**2*z**2)) + h_nm*y*(24*x**6 - 3*y**4*(y**2 + z**2) - 5*x**4*(23*y**2 + 3*z**2) + 6*x**2*(11*y**4 + 5*y**2*z**2)) )
            By = By * rt11o5*( g_nm*1/2*y*(43*x**6 - 7*y**6 + 6*y**4*z**2 + x**4*(-225*y**2 + 30*z**2) + 3*x**2*(47*y**4 - 20*y**2*z**2)) + h_nm*x*(-3*x**6 + 3*y**4*(8*y**2 - 5*z**2) + x**4*(66*y**2 - 3*z**2) - 5*x**2*(23*y**4 - 6*y**2*z**2)) )
            Bz = Bz * 13*rt11o5*( g_nm*1/2*(x**6 - 15*x**4*y**2 + 15*x**2*y**4 - y**6)*z + h_nm*x*y*(3*x**4 - 10*x**2*y**2 + 3*y**4)*z )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi_Schmidt, n={n} and m is not between 0 and n.")

    elif n==7:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   n=7 moment components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A7 = sqrt(1/256)
        Bx = A7/r**17 + 0j
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==0:
            Bx = Bx * g_nm * -9*x*z*(35*x**6 + 35*y**6 - 280*y**4*z**2 + 336*y**2*z**4 - 64*z**6 + 35*x**4*(3*y**2 - 8*z**2) + 7*x**2*(15*y**4 - 80*y**2*z**2 + 48*z**4))
            By = By * g_nm * -9*y*z*(35*x**6 + 35*y**6 - 280*y**4*z**2 + 336*y**2*z**4 - 64*z**6 + 35*x**4*(3*y**2 - 8*z**2) + 7*x**2*(15*y**4 - 80*y**2*z**2 + 48*z**4))
            Bz = Bz * g_nm * (6435*z**8 - 12012*z**6*r**2 + 6930*z**4*r**4 - 1260*z**2*r**6 + 35*r**8)
        elif m==1:
            rt7 = sqrt(7)
            Bx = Bx * -1/2*rt7*( g_nm*(40*x**8 - 5*y**8 + 115*y**6*z**2 - 120*y**4*z**4 - 176*y**2*z**6 + 64*z**8 + 5*x**6*(23*y**2 - 247*z**2) + 15*x**4*(7*y**4 - 157*y**2*z**2 + 232*z**4) + x**2*(25*y**6 - 1005*y**4*z**2 + 3360*y**2*z**4 - 1616*z**6)) + h_nm*45*x*y*(x**6 + y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6 + 3*x**4*(y**2 - 10*z**2) + x**2*(3*y**4 - 60*y**2*z**2 + 80*z**4)) )
            By = By * 1/2*rt7*( g_nm*-45*x*y*(x**6 + y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6 + 3*x**4*(y**2 - 10*z**2) + x**2*(3*y**4 - 60*y**2*z**2 + 80*z**4)) + h_nm*(5*x**8 - 40*y**8 + 1235*y**6*z**2 - 3480*y**4*z**4 + 1616*y**2*z**6 - 64*z**8 - 5*x**6*(5*y**2 + 23*z**2) - 15*x**4*(7*y**4 - 67*y**2*z**2 - 8*z**4) + x**2*(-115*y**6 + 2355*y**4*z**2 - 3360*y**2*z**4 + 176*z**6)) )
            Bz = Bz * -9/2*rt7*( g_nm*x*z*(35*x**6 + 35*y**6 - 280*y**4*z**2 + 336*y**2*z**4 - 64*z**6 + 35*x**4*(3*y**2 - 8*z**2) + 7*x**2*(15*y**4 - 80*y**2*z**2 + 48*z**4)) + h_nm*y*z*(35*x**6 + 35*y**6 - 280*y**4*z**2 + 336*y**2*z**4 - 64*z**6 + 35*x**4*(3*y**2 - 8*z**2) + 7*x**2*(15*y**4 - 80*y**2*z**2 + 48*z**4)) )
        elif m==2:
            rt21o2 = sqrt(21/2)
            Bx = Bx * rt21o2*( g_nm*x*z*(135*x**6 + x**4*(75*y**2 - 970*z**2) + x**2*(-255*y**4 + 260*y**2*z**2 + 944*z**4) - 3*(65*y**6 - 410*y**4*z**2 + 272*y**2*z**4 + 32*z**6)) + h_nm*2*y*z*(150*x**6 - 15*y**6 + 65*y**4*z**2 + 32*y**2*z**4 - 48*z**6 + 15*x**4*(19*y**2 - 69*z**2) + 2*x**2*(60*y**4 - 485*y**2*z**2 + 456*z**4)) )
            By = By * rt21o2*( g_nm*y*z*(195*x**6 - 135*y**6 + 970*y**4*z**2 - 944*y**2*z**4 + 96*z**6 + 15*x**4*(17*y**2 - 82*z**2) + x**2*(-75*y**4 - 260*y**2*z**2 + 816*z**4)) + h_nm*-2*x*z*(15*x**6 - 150*y**6 + 1035*y**4*z**2 - 912*y**2*z**4 + 48*z**6 - 5*x**4*(24*y**2 + 13*z**2) + x**2*(-285*y**4 + 970*y**2*z**2 - 32*z**4)) )
            Bz = Bz * -15*rt21o2*( g_nm*(x**2 - y**2)*(x**6 + y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6 + 3*x**4*(y**2 - 10*z**2) + x**2*(3*y**4 - 60*y**2*z**2 + 80*z**4)) + h_nm*2*x*y*(x**6 + y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6 + 3*x**4*(y**2 - 10*z**2) + x**2*(3*y**4 - 60*y**2*z**2 + 80*z**4)) )
        elif m==3:
            rt21 = sqrt(21)
            Bx = Bx * 3/2*rt21*( g_nm*(8*x**8 + 3*y**8 - 57*y**6*z**2 + 20*y**4*z**4 + 80*y**2*z**6 - x**6*(17*y**2 + 207*z**2) + x**4*(-55*y**4 + 585*y**2*z**2 + 420*z**4) - x**2*(27*y**6 - 735*y**4*z**2 + 1320*y**2*z**4 + 80*z**6)) + h_nm*x*y*(27*x**6 - 17*y**6 + 378*y**4*z**2 - 480*y**2*z**4 - 160*z**6 + x**4*(37*y**2 - 678*z**2) + x**2*(-7*y**4 - 300*y**2*z**2 + 1280*z**4)) )
            By = By * 3/2*rt21*( g_nm*x*y*(17*x**6 - 27*y**6 + 678*y**4*z**2 - 1280*y**2*z**4 + 160*z**6 + 7*x**4*(y**2 - 54*z**2) + x**2*(-37*y**4 + 300*y**2*z**2 + 480*z**4)) + h_nm*-1*(3*x**8 + 8*y**8 - 207*y**6*z**2 + 420*y**4*z**4 - 80*y**2*z**6 - 3*x**6*(9*y**2 + 19*z**2) + x**4*(-55*y**4 + 735*y**2*z**2 + 20*z**4) + x**2*(-17*y**6 + 585*y**4*z**2 - 1320*y**2*z**4 + 80*z**6)) )
            Bz = Bz * 55/2*rt21*( g_nm*x*(x**2 - 3*y**2)*z*(3*x**4 + 3*y**4 - 20*y**2*z**2 + 16*z**4 + x**2*(6*y**2 - 20*z**2)) + h_nm*-y*(-3*x**2 + y**2)*z*(3*x**4 + 3*y**4 - 20*y**2*z**2 + 16*z**4 + x**2*(6*y**2 - 20*z**2)) )
        elif m==4:
            rt231 = sqrt(231)
            Bx = Bx * rt231*( g_nm*x*z*(-27*x**6 + x**4*(183*y**2 + 128*z**2) + 5*x**2*(27*y**4 - 176*y**2*z**2 - 8*z**4) + 15*y**2*(-5*y**4 + 16*y**2*z**2 + 8*z**4)) + h_nm*-4*y*z*(30*x**6 + 3*y**6 - 7*y**4*z**2 - 10*y**2*z**4 - 15*x**4*(y**2 + 9*z**2) + x**2*(-42*y**4 + 170*y**2*z**2 + 30*z**4)) )
            By = By * rt231*( g_nm*y*z*(-75*x**6 - 27*y**6 + 128*y**4*z**2 - 40*y**2*z**4 + 15*x**4*(9*y**2 + 16*z**2) + x**2*(183*y**4 - 880*y**2*z**2 + 120*z**4)) + h_nm*4*x*z*(3*x**6 - 7*x**4*(6*y**2 + z**2) - 5*x**2*(3*y**4 - 34*y**2*z**2 + 2*z**4) + 15*(2*y**6 - 9*y**4*z**2 + 2*y**2*z**4)) )
            Bz = Bz * 3*rt231*( g_nm*(x**4 - 6*x**2*y**2 + y**4)*(x**4 + y**4 - 24*y**2*z**2 + 40*z**4 + 2*x**2*(y**2 - 12*z**2)) + h_nm*4*x*y*(x**2 - y**2)*(x**4 + y**4 - 24*y**2*z**2 + 40*z**4 + 2*x**2*(y**2 - 12*z**2)) )
        elif m==5:
            rt231 = sqrt(231)
            Bx = Bx * 1/2*rt231*( g_nm*(-8*x**8 + x**6*(97*y**2 + 127*z**2) + 5*y**4*(y**4 - 11*y**2*z**2 - 12*z**4) + 15*x**4*(y**4 - 103*y**2*z**2 - 4*z**4) + x**2*(-85*y**6 + 1185*y**4*z**2 + 360*y**2*z**4)) + h_nm*x*y*(-45*x**6 - 33*y**6 + 402*y**4*z**2 + 240*y**2*z**4 + x**4*(85*y**2 + 690*z**2) + x**2*(97*y**4 - 1820*y**2*z**2 - 240*z**4)) )
            By = By * 1/2*rt231*( g_nm*x*y*(-33*x**6 + x**4*(97*y**2 + 402*z**2) - 15*y**2*(3*y**4 - 46*y**2*z**2 + 16*z**4) + 5*x**2*(17*y**4 - 364*y**2*z**2 + 48*z**4)) + h_nm*(5*x**8 - 8*y**8 + 127*y**6*z**2 - 60*y**4*z**4 - 5*x**6*(17*y**2 + 11*z**2) + 15*x**4*(y**4 + 79*y**2*z**2 - 4*z**4) + x**2*(97*y**6 - 1545*y**4*z**2 + 360*y**2*z**4)) )
            Bz = Bz * -39/2*rt231*( g_nm*x*(x**4 - 10*x**2*y**2 + 5*y**4)*z*(x**2 + y**2 - 4*z**2) + h_nm*y*(5*x**4 - 10*x**2*y**2 + y**4)*z*(x**2 + y**2 - 4*z**2) )
        elif m==6:
            rt3003o2 = sqrt(3003/2)
            Bx = Bx * 3*rt3003o2*( g_nm*x*z*(3*x**6 - x**4*(57*y**2 + 2*z**2) + 5*x**2*(17*y**4 + 4*y**2*z**2) - 5*(3*y**6 + 2*y**4*z**2)) + h_nm*-2*y*z*(-10*x**6 + y**4*(y**2 + z**2) + 5*x**4*(9*y**2 + z**2) - 2*x**2*(12*y**4 + 5*y**2*z**2)) )
            By = By * 3*rt3003o2*( g_nm*y*z*(15*x**6 - 3*y**6 + 2*y**4*z**2 + x**4*(-85*y**2 + 10*z**2) + x**2*(57*y**4 - 20*y**2*z**2)) + h_nm*-2*x*z*(x**6 + x**4*(-24*y**2 + z**2) + 5*y**4*(-2*y**2 + z**2) + 5*x**2*(9*y**4 - 2*y**2*z**2)) )
            Bz = Bz * -rt3003o2*( g_nm*(x**6 - 15*x**4*y**2 + 15*x**2*y**4 - y**6)*(x**2 + y**2 - 14*z**2) + h_nm*2*x*y*(3*x**4 - 10*x**2*y**2 + 3*y**4)*(x**2 + y**2 - 14*z**2) )
        elif m==7:
            rt429 = sqrt(429)
            Bx = Bx * 1/2*rt429*( g_nm*(8*x**8 + 7*y**6*(y**2 + z**2) + 105*x**4*y**2*(5*y**2 + z**2) - 7*x**6*(31*y**2 + z**2) - 7*x**2*(29*y**6 + 15*y**4*z**2)) + h_nm*x*y*(63*x**6 - 57*y**6 - 42*y**4*z**2 - 7*x**4*(61*y**2 + 6*z**2) + 7*x**2*(59*y**4 + 20*y**2*z**2)) )
            By = By * 1/2*rt429*( g_nm*x*y*(57*x**6 - 63*y**6 + 42*y**4*z**2 - 7*x**4*(59*y**2 - 6*z**2) + 7*x**2*(61*y**4 - 20*y**2*z**2)) + h_nm*-1*(7*x**8 + 8*y**8 - 7*y**6*z**2 + 7*x**6*(-29*y**2 + z**2) + 105*x**4*(5*y**4 - y**2*z**2) - 7*x**2*(31*y**6 - 15*y**4*z**2)) )
            Bz = Bz * 15/2*rt429*( g_nm*x*(x**6 - 21*x**4*y**2 + 35*x**2*y**4 - 7*y**6)*z + h_nm*-y*(-7*x**6 + 35*x**4*y**2 - 21*x**2*y**4 + y**6)*z )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi_Schmidt, n={n} and m is not between 0 and n.")

    elif n==8:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   n=8 moment components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A8 = 3/32
        Bx = A8/r**19 + 0j
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==0:
            Bx = Bx * g_nm * 15/4*x*(7*x**8 + 7*y**8 - 280*y**6*z**2 + 1120*y**4*z**4 - 896*y**2*z**6 + 128*z**8 + 28*x**6*(y**2 - 10*z**2) + 14*x**4*(3*y**4 - 60*y**2*z**2 + 80*z**4) + 28*x**2*(y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6))
            By = By * g_nm * 15/4*y*(7*x**8 + 7*y**8 - 280*y**6*z**2 + 1120*y**4*z**4 - 896*y**2*z**6 + 128*z**8 + 28*x**6*(y**2 - 10*z**2) + 14*x**4*(3*y**4 - 60*y**2*z**2 + 80*z**4) + 28*x**2*(y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6))
            Bz = Bz * g_nm * 3/4*z*(315*x**8 + 315*y**8 - 3360*y**6*z**2 + 6048*y**4*z**4 - 2304*y**2*z**6 + 128*z**8 + 420*x**6*(3*y**2 - 8*z**2) + 126*x**4*(15*y**4 - 80*y**2*z**2 + 48*z**4) + 36*x**2*(35*y**6 - 280*y**4*z**2 + 336*y**2*z**4 - 64*z**6))
        elif m==1:
            Bx = Bx * ( g_nm*z*(-350*x**8 + 35*y**8 - 245*y**6*z**2 + 56*y**4*z**4 + 272*y**2*z**6 - 64*z**8 - 35*x**6*(29*y**2 - 103*z**2) - 7*x**4*(135*y**4 - 995*y**2*z**2 + 872*z**4) + x**2*(-245*y**6 + 3115*y**4*z**2 - 6048*y**2*z**4 + 2032*z**6)) + h_nm*-55*x*y*z*(7*x**6 + 7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6 + 7*x**4*(3*y**2 - 10*z**2) + 7*x**2*(3*y**4 - 20*y**2*z**2 + 16*z**4)) )
            By = By * ( g_nm*-55*x*y*z*(7*x**6 + 7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6 + 7*x**4*(3*y**2 - 10*z**2) + 7*x**2*(3*y**4 - 20*y**2*z**2 + 16*z**4)) + h_nm*z*(35*x**8 - 350*y**8 + 3605*y**6*z**2 - 6104*y**4*z**4 + 2032*y**2*z**6 - 64*z**8 - 245*x**6*(y**2 + z**2) + x**4*(-945*y**4 + 3115*y**2*z**2 + 56*z**4) + x**2*(-1015*y**6 + 6965*y**4*z**2 - 6048*y**2*z**4 + 272*z**6)) )
            Bz = Bz * 5*( g_nm*x*(7*x**8 + 7*y**8 - 280*y**6*z**2 + 1120*y**4*z**4 - 896*y**2*z**6 + 128*z**8 + 28*x**6*(y**2 - 10*z**2) + 14*x**4*(3*y**4 - 60*y**2*z**2 + 80*z**4) + 28*x**2*(y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6)) + h_nm*y*(7*x**8 + 7*y**8 - 280*y**6*z**2 + 1120*y**4*z**4 - 896*y**2*z**6 + 128*z**8 + 28*x**6*(y**2 - 10*z**2) + 14*x**4*(3*y**4 - 60*y**2*z**2 + 80*z**4) + 28*x**2*(y**6 - 30*y**4*z**2 + 80*y**2*z**4 - 32*z**6)) )
        elif m==2:
            rt35o2 = sqrt(35/2)
            Bx = Bx * rt35o2*( g_nm*x*(-9*x**8 + 13*y**8 - 454*y**6*z**2 + 1420*y**4*z**4 - 608*y**2*z**6 - 64*z**8 + x**6*(-14*y**2 + 338*z**2) + 2*x**4*(6*y**4 + 111*y**2*z**2 - 610*z**4) + 10*x**2*(3*y**6 - 57*y**4*z**2 + 20*y**2*z**4 + 80*z**6)) + h_nm*2*y*(-10*x**8 + y**8 - 29*y**6*z**2 + 50*y**4*z**4 + 48*y**2*z**6 - 32*z**8 + x**6*(-29*y**2 + 367*z**2) + x**4*(-27*y**4 + 705*y**2*z**2 - 1270*z**4) + x**2*(-7*y**6 + 309*y**4*z**2 - 1220*y**2*z**4 + 752*z**6)) )
            By = By * rt35o2*( g_nm*y*(-13*x**8 + 9*y**8 - 338*y**6*z**2 + 1220*y**4*z**4 - 800*y**2*z**6 + 64*z**8 + x**6*(-30*y**2 + 454*z**2) - 2*x**4*(6*y**4 - 285*y**2*z**2 + 710*z**4) + 2*x**2*(7*y**6 - 111*y**4*z**2 - 100*y**2*z**4 + 304*z**6)) + h_nm*2*x*(x**8 - 10*y**8 + 367*y**6*z**2 - 1270*y**4*z**4 + 752*y**2*z**6 - 32*z**8 - x**6*(7*y**2 + 29*z**2) + x**4*(-27*y**4 + 309*y**2*z**2 + 50*z**4) + x**2*(-29*y**6 + 705*y**4*z**2 - 1220*y**2*z**4 + 48*z**6)) )
            Bz = Bz * -11*rt35o2*( g_nm*(x**2 - y**2)*z*(7*x**6 + 7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6 + 7*x**4*(3*y**2 - 10*z**2) + 7*x**2*(3*y**4 - 20*y**2*z**2 + 16*z**4)) + h_nm*2*x*y*z*(7*x**6 + 7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6 + 7*x**4*(3*y**2 - 10*z**2) + 7*x**2*(3*y**4 - 20*y**2*z**2 + 16*z**4)) )
        elif m==3:
            rt1155 = sqrt(1155)
            Bx = Bx * rt1155*( g_nm*-z*(-10*x**8 - 3*y**8 + 17*y**6*z**2 + 4*y**4*z**4 - 16*y**2*z**6 + x**6*(19*y**2 + 87*z**2) + x**4*(65*y**4 - 225*y**2*z**2 - 108*z**4) + x**2*(33*y**6 - 295*y**4*z**2 + 312*y**2*z**4 + 16*z**6)) + h_nm*x*y*z*(33*x**6 - 19*y**6 + 138*y**4*z**2 - 96*y**2*z**4 - 32*z**6 + x**4*(47*y**2 - 278*z**2) - 5*x**2*(y**4 + 28*y**2*z**2 - 64*z**4)) )
            By = By * rt1155*( g_nm*x*y*z*(19*x**6 - 33*y**6 + 278*y**4*z**2 - 320*y**2*z**4 + 32*z**6 + x**4*(5*y**2 - 138*z**2) + x**2*(-47*y**4 + 140*y**2*z**2 + 96*z**4)) + h_nm*z*(-3*x**8 - 10*y**8 + 87*y**6*z**2 - 108*y**4*z**4 + 16*y**2*z**6 + x**6*(33*y**2 + 17*z**2) + x**4*(65*y**4 - 295*y**2*z**2 + 4*z**4) + x**2*(19*y**6 - 225*y**4*z**2 + 312*y**2*z**4 - 16*z**6)) )
            Bz = Bz * -rt1155*( g_nm*x*(x**2 - 3*y**2)*(x**6 + y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6 + 3*x**4*(y**2 - 12*z**2) + 3*x**2*(y**4 - 24*y**2*z**2 + 40*z**4)) + h_nm*y*(3*x**2 - y**2)*(x**6 + y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6 + 3*x**4*(y**2 - 12*z**2) + 3*x**2*(y**4 - 24*y**2*z**2 + 40*z**4)) )
        elif m==4:
            rt77 = sqrt(77)
            Bx = Bx * rt77*( g_nm*1/2*x*(9*x**8 - 4*x**6*(13*y**2 + 68*z**2) + x**4*(-106*y**4 + 1728*y**2*z**2 + 664*z**4) - 20*x**2*(y**6 - 68*y**4*z**2 + 212*y**2*z**4 + 8*z**6) + 5*y**2*(5*y**6 - 128*y**4*z**2 + 184*y**2*z**4 + 96*z**6)) + h_nm*2*y*(10*x**8 + y**8 - 23*y**6*z**2 + 16*y**4*z**4 + 40*y**2*z**6 + 5*x**6*(y**2 - 59*z**2) + x**4*(-19*y**4 + 115*y**2*z**2 + 680*z**4) - x**2*(13*y**6 - 387*y**4*z**2 + 760*y**2*z**4 + 120*z**6)) )
            By = By * rt77*( g_nm*1/2*y*(25*x**8 + 9*y**8 - 272*y**6*z**2 + 664*y**4*z**4 - 160*y**2*z**6 - 20*x**6*(y**2 + 32*z**2) + x**4*(-106*y**4 + 1360*y**2*z**2 + 920*z**4) + x**2*(-52*y**6 + 1728*y**4*z**2 - 4240*y**2*z**4 + 480*z**6)) + h_nm*-2*x*(x**8 - x**6*(13*y**2 + 23*z**2) + x**4*(-19*y**4 + 387*y**2*z**2 + 16*z**4) + 5*x**2*(y**6 + 23*y**4*z**2 - 152*y**2*z**4 + 8*z**6) + 5*(2*y**8 - 59*y**6*z**2 + 136*y**4*z**4 - 24*y**2*z**6)) )
            Bz = Bz * 65*rt77*( g_nm*1/2*(x**4 - 6*x**2*y**2 + y**4)*z*(x**4 + y**4 - 8*y**2*z**2 + 8*z**4 + 2*x**2*(y**2 - 4*z**2)) + h_nm*2*x*y*(x**2 - y**2)*z*(x**4 + y**4 - 8*y**2*z**2 + 8*z**4 + 2*x**2*(y**2 - 4*z**2)) )
        elif m==5:
            rt1001 = sqrt(1001)
            Bx = Bx * -5*rt1001*( g_nm*z*(2*x**8 - y**8 + 3*y**6*z**2 + 4*y**4*z**4 - x**6*(23*y**2 + 11*z**2) + x**4*(-5*y**4 + 125*y**2*z**2 + 4*z**4) + x**2*(19*y**6 - 85*y**4*z**2 - 24*y**2*z**4)) + h_nm*x*y*z*(11*x**6 + 7*y**6 - 26*y**4*z**2 - 16*y**2*z**4 - x**4*(19*y**2 + 58*z**2) + x**2*(-23*y**4 + 140*y**2*z**2 + 16*z**4)) )
            By = By * 5*rt1001*( g_nm*-x*y*z*(7*x**6 + 11*y**6 - 58*y**4*z**2 + 16*y**2*z**4 - x**4*(23*y**2 + 26*z**2) + x**2*(-19*y**4 + 140*y**2*z**2 - 16*z**4)) + h_nm*z*(x**8 - 2*y**8 + 11*y**6*z**2 - 4*y**4*z**4 - x**6*(19*y**2 + 3*z**2) + x**4*(5*y**4 + 85*y**2*z**2 - 4*z**4) + x**2*(23*y**6 - 125*y**4*z**2 + 24*y**2*z**4)) )
            Bz = Bz * rt1001*( g_nm*x*(x**4 - 10*x**2*y**2 + 5*y**4)*(x**4 + y**4 - 28*y**2*z**2 + 56*z**4 + 2*x**2*(y**2 - 14*z**2)) + h_nm*y*(5*x**4 - 10*x**2*y**2 + y**4)*(x**4 + y**4 - 28*y**2*z**2 + 56*z**4 + 2*x**2*(y**2 - 14*z**2)) )
        elif m==6:
            rt429o2 = sqrt(429/2)
            Bx = Bx * rt429o2*( g_nm*-x*(3*x**8 - 54*x**6*(y**2 + z**2) + 14*x**4*(2*y**4 + 69*y**2*z**2 + 2*z**4) + 5*y**4*(-3*y**4 + 42*y**2*z**2 + 28*z**4) + 70*x**2*(y**6 - 19*y**4*z**2 - 4*y**2*z**4)) + h_nm*2*y*(-10*x**8 + y**8 - 13*y**6*z**2 - 14*y**4*z**4 + 35*x**6*(y**2 + 5*z**2) + 7*x**4*(3*y**4 - 105*y**2*z**2 - 10*z**4) + x**2*(-23*y**6 + 357*y**4*z**2 + 140*y**2*z**4)) )
            By = By * rt429o2*( g_nm*y*(-15*x**8 + 3*y**8 - 54*y**6*z**2 + 28*y**4*z**4 + 70*x**6*(y**2 + 3*z**2) + 14*x**4*(2*y**4 - 95*y**2*z**2 + 10*z**4) + x**2*(-54*y**6 + 966*y**4*z**2 - 280*y**2*z**4)) + h_nm*2*x*(x**8 - x**6*(23*y**2 + 13*z**2) + 7*x**4*(3*y**4 + 51*y**2*z**2 - 2*z**4) - 5*y**4*(2*y**4 - 35*y**2*z**2 + 14*z**4) + 35*x**2*(y**6 - 21*y**4*z**2 + 4*y**2*z**4)) )
            Bz = Bz * -5*rt429o2*( g_nm*(x**6 - 15*x**4*y**2 + 15*x**2*y**4 - y**6)*z*(3*x**2 + 3*y**2 - 14*z**2) + h_nm*2*x*y*(3*x**4 - 10*x**2*y**2 + 3*y**4)*z*(3*x**2 + 3*y**2 - 14*z**2) )
        elif m==7:
            rt715 = sqrt(715)
            Bx = Bx * rt715*( g_nm*z*(10*x**8 + 7*y**6*(y**2 + z**2) - 7*x**6*(37*y**2 + z**2) + 35*x**4*(17*y**4 + 3*y**2*z**2) - 7*x**2*(31*y**6 + 15*y**4*z**2)) + h_nm*x*y*z*(77*x**6 - 59*y**6 - 42*y**4*z**2 - 7*x**4*(71*y**2 + 6*z**2) + 35*x**2*(13*y**4 + 4*y**2*z**2)) )
            By = By * rt715*( g_nm*x*y*z*(59*x**6 - 77*y**6 + 42*y**4*z**2 - 7*x**4*(65*y**2 - 6*z**2) + 7*x**2*(71*y**4 - 20*y**2*z**2)) + h_nm*z*(-7*x**8 - 10*y**8 + 7*y**6*z**2 + 7*x**6*(31*y**2 - z**2) - 35*x**4*(17*y**4 - 3*y**2*z**2) + 7*x**2*(37*y**6 - 15*y**4*z**2)) )
            Bz = Bz * rt715*( g_nm*-x*(x**6 - 21*x**4*y**2 + 35*x**2*y**4 - 7*y**6)*(x**2 + y**2 - 16*z**2) + h_nm*y*(-7*x**6 + 35*x**4*y**2 - 21*x**2*y**4 + y**6)*(x**2 + y**2 - 16*z**2) )
        elif m==8:
            rt715 = sqrt(715)
            Bx = Bx * rt715*( g_nm*1/4*x*(9*x**8 + 73*y**8 + 56*y**6*z**2 - 4*x**6*(79*y**2 + 2*z**2) + 14*x**4*(77*y**4 + 12*y**2*z**2) - 140*x**2*(5*y**6 + 2*y**4*z**2)) + h_nm*2*y*(10*x**8 + y**6*(y**2 + z**2) - 7*x**6*(13*y**2 + z**2) + 7*x**4*(19*y**4 + 5*y**2*z**2) - x**2*(37*y**6 + 21*y**4*z**2)) )
            By = By * rt715*( g_nm*1/4*y*(73*x**8 + 9*y**8 - 8*y**6*z**2 + x**6*(-700*y**2 + 56*z**2) + 14*x**4*(77*y**4 - 20*y**2*z**2) - 4*x**2*(79*y**6 - 42*y**4*z**2)) + h_nm*-2*x*(x**8 + 10*y**8 - 7*y**6*z**2 + x**6*(-37*y**2 + z**2) + 7*x**4*(19*y**4 - 3*y**2*z**2) + x**2*(-91*y**6 + 35*y**4*z**2)) )
            Bz = Bz * 17*rt715*( g_nm*1/4*(x**8 - 28*x**6*y**2 + 70*x**4*y**4 - 28*x**2*y**6 + y**8)*z + h_nm*2*x*y*(x**6 - 7*x**4*y**2 + 7*x**2*y**4 - y**6)*z )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi_Schmidt, n={n} and m is not between 0 and n.")

    elif n==9:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   n=9 moment components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A9 = sqrt(5)/16
        Bx = A9/r**21 + 0j
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==0:
            rt5 = sqrt(5)
            Bx = Bx * g_nm * 11/8*rt5*x*z*(63*x**8 + 63*y**8 - 840*y**6*z**2 + 2016*y**4*z**4 - 1152*y**2*z**6 + 128*z**8 + 84*x**6*(3*y**2 - 10*z**2) + 126*x**4*(3*y**4 - 20*y**2*z**2 + 16*z**4) + 36*x**2*(7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6))
            By = By * g_nm * 11/8*rt5*y*z*(63*x**8 + 63*y**8 - 840*y**6*z**2 + 2016*y**4*z**4 - 1152*y**2*z**6 + 128*z**8 + 84*x**6*(3*y**2 - 10*z**2) + 126*x**4*(3*y**4 - 20*y**2*z**2 + 16*z**4) + 36*x**2*(7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6))
            Bz = Bz * g_nm * 1/8*rt5*(46189*z**10 - 109395*z**8*r**2 + 90090*z**6*r**4 - 30030*z**4*r**6 + 3465*z**2*r**8 - 63*r**10)
        elif m==1:
            Bx = Bx * 3/8*( g_nm*(70*x**10 - 7*y**10 + 273*y**8*z**2 - 840*y**6*z**4 - 224*y**4*z**6 + 768*y**2*z**8 - 128*z**10 + 21*x**8*(13*y**2 - 163*z**2) + 196*x**6*(2*y**4 - 51*y**2*z**2 + 90*z**4) + 14*x**4*(17*y**6 - 675*y**4*z**2 + 2460*y**2*z**4 - 1424*z**6) + 6*x**2*(7*y**8 - 434*y**6*z**2 + 2660*y**4*z**4 - 3360*y**2*z**6 + 832*z**8)) + h_nm*11*x*y*(7*x**8 + 7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8 + 28*x**6*(y**2 - 12*z**2) + 42*x**4*(y**4 - 24*y**2*z**2 + 40*z**4) + 28*x**2*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6)) )
            By = By * 3/8*( g_nm*11*x*y*(7*x**8 + 7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8 + 28*x**6*(y**2 - 12*z**2) + 42*x**4*(y**4 - 24*y**2*z**2 + 40*z**4) + 28*x**2*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6)) + h_nm*-1*(7*x**10 - 70*y**10 + 3423*y**8*z**2 - 17640*y**6*z**4 + 19936*y**4*z**6 - 4992*y**2*z**8 + 128*z**10 - 21*x**8*(2*y**2 + 13*z**2) + x**6*(-238*y**4 + 2604*y**2*z**2 + 840*z**4) - 14*x**4*(28*y**6 - 675*y**4*z**2 + 1140*y**2*z**4 - 16*z**6) - 3*x**2*(91*y**8 - 3332*y**6*z**2 + 11480*y**4*z**4 - 6720*y**2*z**6 + 256*z**8)) )
            Bz = Bz * 33/8*( g_nm*x*z*(63*x**8 + 63*y**8 - 840*y**6*z**2 + 2016*y**4*z**4 - 1152*y**2*z**6 + 128*z**8 + 84*x**6*(3*y**2 - 10*z**2) + 126*x**4*(3*y**4 - 20*y**2*z**2 + 16*z**4) + 36*x**2*(7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6)) + h_nm*y*z*(63*x**8 + 63*y**8 - 840*y**6*z**2 + 2016*y**4*z**4 - 1152*y**2*z**6 + 128*z**8 + 84*x**6*(3*y**2 - 10*z**2) + 126*x**4*(3*y**4 - 20*y**2*z**2 + 16*z**4) + 36*x**2*(7*y**6 - 70*y**4*z**2 + 112*y**2*z**4 - 32*z**6)) )
        elif m==2:
            rt11o2 = sqrt(11/2)
            Bx = Bx * -3*rt11o2*( g_nm*1/2*x*z*(77*x**8 - 105*y**8 + 1218*y**6*z**2 - 2268*y**4*z**4 + 672*y**2*z**6 + 64*z**8 + 42*x**6*(3*y**2 - 23*z**2) - 42*x**4*(2*y**4 + 17*y**2*z**2 - 50*z**4) - 2*x**2*(119*y**6 - 735*y**4*z**2 + 84*y**2*z**4 + 496*z**6)) + h_nm*y*z*(84*x**8 - 7*y**8 + 63*y**6*z**2 - 42*y**4*z**4 - 80*y**2*z**6 + 32*z**8 + 49*x**6*(5*y**2 - 21*z**2) + 21*x**4*(11*y**4 - 95*y**2*z**2 + 102*z**4) + 3*x**2*(21*y**6 - 301*y**4*z**2 + 700*y**2*z**4 - 304*z**6)) )
            By = By * 3*rt11o2*( g_nm*-1/2*y*z*(105*x**8 - 77*y**8 + 966*y**6*z**2 - 2100*y**4*z**4 + 992*y**2*z**6 - 64*z**8 + 14*x**6*(17*y**2 - 87*z**2) + 42*x**4*(2*y**4 - 35*y**2*z**2 + 54*z**4) - 42*x**2*(3*y**6 - 17*y**4*z**2 - 4*y**2*z**4 + 16*z**6)) + h_nm*x*z*(7*x**8 - 84*y**8 + 1029*y**6*z**2 - 2142*y**4*z**4 + 912*y**2*z**6 - 32*z**8 - 63*x**6*(y**2 + z**2) + x**4*(-231*y**4 + 903*y**2*z**2 + 42*z**4) - 5*x**2*(49*y**6 - 399*y**4*z**2 + 420*y**2*z**4 - 16*z**6)) )
            Bz = Bz * 3*rt11o2*( g_nm*1/2*(x**2 - y**2)*(7*x**8 + 7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8 + 28*x**6*(y**2 - 12*z**2) + 42*x**4*(y**4 - 24*y**2*z**2 + 40*z**4) + 28*x**2*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6)) + h_nm*x*y*(7*x**8 + 7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8 + 28*x**6*(y**2 - 12*z**2) + 42*x**4*(y**4 - 24*y**2*z**2 + 40*z**4) + 28*x**2*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6)) )
        elif m==3:
            rt231o2 = sqrt(231/2)
            Bx = Bx * 1/4*rt231o2*( g_nm*(-10*x**10 + 9*x**8*(y**2 + 49*z**2) + 84*x**6*(y**4 - 9*y**2*z**2 - 23*z**4) + 14*x**4*(7*y**6 - 195*y**4*z**2 + 330*y**2*z**4 + 116*z**6) - 3*y**2*(y**8 - 35*y**6*z**2 + 84*y**4*z**4 + 56*y**2*z**6 - 64*z**8) + 6*x**2*(5*y**8 - 238*y**6*z**2 + 1050*y**4*z**4 - 728*y**2*z**6 - 32*z**8)) + h_nm*-x*y*(33*x**8 - 19*y**8 + 756*y**6*z**2 - 2688*y**4*z**4 + 1120*y**2*z**6 + 384*z**8 + 4*x**6*(20*y**2 - 357*z**2) + 42*x**4*(y**4 - 50*y**2*z**2 + 144*z**4) - 12*x**2*(2*y**6 - 7*y**4*z**2 - 280*y**2*z**4 + 392*z**6)) )
            By = By * 1/4*rt231o2*( g_nm*x*y*(-19*x**8 + 33*y**8 - 1428*y**6*z**2 + 6048*y**4*z**4 - 4704*y**2*z**6 + 384*z**8 + x**6*(-24*y**2 + 756*z**2) + 42*x**4*(y**4 + 2*y**2*z**2 - 64*z**4) + 20*x**2*(4*y**6 - 105*y**4*z**2 + 168*y**2*z**4 + 56*z**6)) + h_nm*(3*x**10 - 15*x**8*(2*y**2 + 7*z**2) + x**6*(-98*y**4 + 1428*y**2*z**2 + 252*z**4) - 42*x**4*(2*y**6 - 65*y**4*z**2 + 150*y**2*z**4 - 4*z**6) - 3*x**2*(3*y**8 - 252*y**6*z**2 + 1540*y**4*z**4 - 1456*y**2*z**6 + 64*z**8) + y**2*(10*y**8 - 441*y**6*z**2 + 1932*y**4*z**4 - 1624*y**2*z**6 + 192*z**8)) )
            Bz = Bz * 13/4*rt231o2*( g_nm*-x*(x**2 - 3*y**2)*z*(7*x**6 + 7*y**6 - 84*y**4*z**2 + 168*y**2*z**4 - 64*z**6 + 21*x**4*(y**2 - 4*z**2) + 21*x**2*(y**4 - 8*y**2*z**2 + 8*z**4)) + h_nm*y*(-3*x**2 + y**2)*z*(7*x**6 + 7*y**6 - 84*y**4*z**2 + 168*y**2*z**4 - 64*z**6 + 21*x**4*(y**2 - 4*z**2) + 21*x**2*(y**4 - 8*y**2*z**2 + 8*z**4)) )
        elif m==4:
            rt1001 = sqrt(1001)
            Bx = Bx * 3*rt1001*( g_nm*1/4*x*z*(11*x**8 + 27*y**8 - 224*y**6*z**2 + 168*y**4*z**4 + 96*y**2*z**6 - 4*x**6*(15*y**2 + 28*z**2) - 42*x**4*(3*y**4 - 16*y**2*z**2 - 4*z**4) - 4*x**2*(7*y**6 - 140*y**4*z**2 + 252*y**2*z**4 + 8*z**6)) + h_nm*y*z*(12*x**8 + y**8 - 7*y**6*z**2 + 8*y**2*z**6 + 7*x**6*(y**2 - 17*z**2) - 7*x**4*(3*y**4 - 5*y**2*z**2 - 24*z**4) - 3*x**2*(5*y**6 - 49*y**4*z**2 + 56*y**2*z**4 + 8*z**6)) )
            By = By * 3*rt1001*( g_nm*1/4*y*z*(27*x**8 + 11*y**8 - 112*y**6*z**2 + 168*y**4*z**4 - 32*y**2*z**6 - 28*x**6*(y**2 + 8*z**2) - 14*x**4*(9*y**4 - 40*y**2*z**2 - 12*z**4) - 12*x**2*(5*y**6 - 56*y**4*z**2 + 84*y**2*z**4 - 8*z**6)) + h_nm*-x*z*(x**8 + 12*y**8 - 119*y**6*z**2 + 168*y**4*z**4 - 24*y**2*z**6 - x**6*(15*y**2 + 7*z**2) - 21*x**4*(y**4 - 7*y**2*z**2) + x**2*(7*y**6 + 35*y**4*z**2 - 168*y**2*z**4 + 8*z**6)) )
            Bz = Bz * -3*rt1001*( g_nm*1/4*(x**4 - 6*x**2*y**2 + y**4)*(x**6 + y**6 - 42*y**4*z**2 + 168*y**2*z**4 - 112*z**6 + 3*x**4*(y**2 - 14*z**2) + 3*x**2*(y**4 - 28*y**2*z**2 + 56*z**4)) + h_nm*x*y*(x**2 - y**2)*(x**6 + y**6 - 42*y**4*z**2 + 168*y**2*z**4 - 112*z**6 + 3*x**4*(y**2 - 14*z**2) + 3*x**2*(y**4 - 28*y**2*z**2 + 56*z**4)) )
        elif m==5:
            rt715o2 = sqrt(715/2)
            Bx = Bx * 3/4*rt715o2*( g_nm*(2*x**10 - 3*x**8*(7*y**2 + 23*z**2) - 28*x**6*(y**4 - 27*y**2*z**2 - 7*z**4) + 14*x**4*(y**6 + 15*y**4*z**2 - 150*y**2*z**4 - 4*z**6) - y**4*(y**6 - 27*y**4*z**2 + 28*y**2*z**4 + 56*z**6) + 6*x**2*(3*y**8 - 98*y**6*z**2 + 210*y**4*z**4 + 56*y**2*z**6)) + h_nm*x*y*(11*x**8 + 7*y**8 - 204*y**6*z**2 + 336*y**4*z**4 + 224*y**2*z**6 - 4*x**6*(2*y**2 + 93*z**2) - 42*x**4*(y**4 - 14*y**2*z**2 - 24*z**4) - 4*x**2*(4*y**6 - 189*y**4*z**2 + 560*y**2*z**4 + 56*z**6)) )
            By = By * 3/4*rt715o2*( g_nm*x*y*(7*x**8 + 11*y**8 - 372*y**6*z**2 + 1008*y**4*z**4 - 224*y**2*z**6 - 4*x**6*(4*y**2 + 51*z**2) - 42*x**4*(y**4 - 18*y**2*z**2 - 8*z**4) - 4*x**2*(2*y**6 - 147*y**4*z**2 + 560*y**2*z**4 - 56*z**6)) + h_nm*-1*(x**10 - 2*y**10 + 69*y**8*z**2 - 196*y**6*z**4 + 56*y**4*z**6 - 9*x**8*(2*y**2 + 3*z**2) - 14*x**6*(y**4 - 42*y**2*z**2 - 2*z**4) + 14*x**4*(2*y**6 - 15*y**4*z**2 - 90*y**2*z**4 + 4*z**6) + 21*x**2*(y**8 - 36*y**6*z**2 + 100*y**4*z**4 - 16*y**2*z**6)) )
            Bz = Bz * 3/4*rt715o2*( g_nm*x*(x**4 - 10*x**2*y**2 + 5*y**4)*z*(15*x**4 + 15*y**4 - 140*y**2*z**2 + 168*z**4 + 10*x**2*(3*y**2 - 14*z**2)) + h_nm*y*(5*x**4 - 10*x**2*y**2 + y**4)*z*(15*x**4 + 15*y**4 - 140*y**2*z**2 + 168*z**4 + 10*x**2*(3*y**2 - 14*z**2)) )
        elif m==6:
            rt429o2 = sqrt(429/2)
            Bx = Bx * rt429o2*( g_nm*-1/2*x*z*(33*x**8 - 141*y**8 + 602*y**6*z**2 + 420*y**4*z**4 - 2*x**6*(285*y**2 + 103*z**2) + 42*x**4*(6*y**4 + 83*y**2*z**2 + 2*z**4) + 42*x**2*(17*y**6 - 105*y**4*z**2 - 20*y**2*z**4)) + h_nm*y*z*(-108*x**8 + 9*y**8 - 33*y**6*z**2 - 42*y**4*z**4 + 21*x**6*(17*y**2 + 31*z**2) + 7*x**4*(33*y**4 - 365*y**2*z**2 - 30*z**4) + x**2*(-225*y**6 + 1113*y**4*z**2 + 420*y**2*z**4)) )
            By = By * rt429o2*( g_nm*1/2*y*z*(-141*x**8 + 33*y**8 - 206*y**6*z**2 + 84*y**4*z**4 + 14*x**6*(51*y**2 + 43*z**2) + 42*x**4*(6*y**4 - 105*y**2*z**2 + 10*z**4) - 6*x**2*(95*y**6 - 581*y**4*z**2 + 140*y**2*z**4)) + h_nm*x*z*(9*x**8 - 3*x**6*(75*y**2 + 11*z**2) + 21*x**4*(11*y**4 + 53*y**2*z**2 - 2*z**4) - 3*y**4*(36*y**4 - 217*y**2*z**2 + 70*z**4) + 7*x**2*(51*y**6 - 365*y**4*z**2 + 60*y**2*z**4)) )
            Bz = Bz * rt429o2*( g_nm*1/2*(x**6 - 15*x**4*y**2 + 15*x**2*y**4 - y**6)*(3*x**4 + 3*y**4 - 96*y**2*z**2 + 224*z**4 + 6*x**2*(y**2 - 16*z**2)) + h_nm*x*y*(3*x**4 - 10*x**2*y**2 + 3*y**4)*(3*x**4 + 3*y**4 - 96*y**2*z**2 + 224*z**4 + 6*x**2*(y**2 - 16*z**2)) )
        elif m==7:
            rt143o2 = sqrt(143/2)
            Bx = Bx * -3/8*rt143o2*( g_nm*(10*x**10 - 3*x**8*(83*y**2 + 67*z**2) + 7*y**6*(y**4 - 15*y**2*z**2 - 16*z**4) + 28*x**6*(12*y**4 + 177*y**2*z**2 + 4*z**4) + 42*x**4*(9*y**6 - 255*y**4*z**2 - 40*y**2*z**4) - 42*x**2*(5*y**8 - 86*y**6*z**2 - 40*y**4*z**4)) + h_nm*x*y*(77*x**8 - 59*y**8 + 936*y**6*z**2 + 672*y**4*z**4 - 84*x**6*(5*y**2 + 18*z**2) - 42*x**4*(y**4 - 220*y**2*z**2 - 16*z**4) + 4*x**2*(99*y**6 - 1974*y**4*z**2 - 560*y**2*z**4)) )
            By = By * 3/8*rt143o2*( g_nm*-x*y*(59*x**8 - 36*x**6*(11*y**2 + 26*z**2) + 42*x**4*(y**4 + 188*y**2*z**2 - 16*z**4) - 7*y**4*(11*y**4 - 216*y**2*z**2 + 96*z**4) + 140*x**2*(3*y**6 - 66*y**4*z**2 + 16*y**2*z**4)) + h_nm*(7*x**10 - 105*x**8*(2*y**2 + z**2) + 14*x**6*(27*y**4 + 258*y**2*z**2 - 8*z**4) + y**6*(10*y**4 - 201*y**2*z**2 + 112*z**4) + 42*x**4*(8*y**6 - 255*y**4*z**2 + 40*y**2*z**4) - 3*x**2*(83*y**8 - 1652*y**6*z**2 + 560*y**4*z**4)) )
            Bz = Bz * 51/8*rt143o2*( g_nm*-x*(x**6 - 21*x**4*y**2 + 35*x**2*y**4 - 7*y**6)*z*(3*x**2 + 3*y**2 - 16*z**2) + h_nm*y*(-7*x**6 + 35*x**4*y**2 - 21*x**2*y**4 + y**6)*z*(3*x**2 + 3*y**2 - 16*z**2) )
        elif m==8:
            rt2431 = sqrt(2431)
            Bx = Bx * 3*rt2431*( g_nm*1/8*x*z*(11*x**8 + 75*y**8 + 56*y**6*z**2 - 4*x**6*(93*y**2 + 2*z**2) + 42*x**4*(29*y**4 + 4*y**2*z**2) - 28*x**2*(27*y**6 + 10*y**4*z**2)) + h_nm*y*z*(12*x**8 + y**6*(y**2 + z**2) - 7*x**6*(15*y**2 + z**2) + 7*x**4*(21*y**4 + 5*y**2*z**2) - 3*x**2*(13*y**6 + 7*y**4*z**2)) )
            By = By * 3*rt2431*( g_nm*1/8*y*z*(75*x**8 + 11*y**8 - 8*y**6*z**2 + x**6*(-756*y**2 + 56*z**2) + 14*x**4*(87*y**4 - 20*y**2*z**2) + x**2*(-372*y**6 + 168*y**4*z**2)) + h_nm*-x*z*(x**8 + 12*y**8 - 7*y**6*z**2 + x**6*(-39*y**2 + z**2) + 21*x**4*(7*y**4 - y**2*z**2) - 35*x**2*(3*y**6 - y**4*z**2)) )
            Bz = Bz * -3*rt2431*( g_nm*1/8*(x**8 - 28*x**6*y**2 + 70*x**4*y**4 - 28*x**2*y**6 + y**8)*(x**2 + y**2 - 18*z**2) + h_nm*x*y*(x**6 - 7*x**4*y**2 + 7*x**2*y**4 - y**6)*(x**2 + y**2 - 18*z**2) )
        elif m==9:
            rt2431o2 = sqrt(2431/2)
            Bx = Bx * 1/8*rt2431o2*( g_nm*(10*x**10 - 9*y**8*(y**2 + z**2) + 252*x**6*y**2*(8*y**2 + z**2) - 9*x**8*(49*y**2 + z**2) - 42*x**4*(47*y**6 + 15*y**4*z**2) + 18*x**2*(23*y**8 + 14*y**6*z**2)) + h_nm*x*y*(99*x**8 + 91*y**8 + 72*y**6*z**2 - 12*x**6*(97*y**2 + 6*z**2) + 126*x**4*(19*y**4 + 4*y**2*z**2) - 36*x**2*(31*y**6 + 14*y**4*z**2)) )
            By = By * 1/8*rt2431o2*( g_nm*x*y*(91*x**8 + 99*y**8 - 72*y**6*z**2 - 36*x**6*(31*y**2 - 2*z**2) + 126*x**4*(19*y**4 - 4*y**2*z**2) - 12*x**2*(97*y**6 - 42*y**4*z**2)) + h_nm*-1*(9*x**10 - 10*y**10 + 9*y**8*z**2 + 9*x**8*(-46*y**2 + z**2) + 42*x**6*(47*y**4 - 6*y**2*z**2) - 126*x**4*(16*y**6 - 5*y**4*z**2) + 63*x**2*(7*y**8 - 4*y**6*z**2 )) )
            Bz = Bz * 19/8*rt2431o2*( g_nm*x*(x**8 - 36*x**6*y**2 + 126*x**4*y**4 - 84*x**2*y**6 + 9*y**8)*z + h_nm*y*(9*x**8 - 84*x**6*y**2 + 126*x**4*y**4 - 36*x**2*y**6 + y**8)*z )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi_Schmidt, n={n} and m is not between 0 and n.")

    elif n==10:
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        #   n=10 moment components
        #@@@@@@@@@@@@@@@@@@@@@@@@@@@@

        A10 = sqrt(11)/64
        Bx = A10/r**23 + 0j
        By = Bx + 0.0
        Bz = Bx + 0.0

        if m==0:
            rt11 = sqrt(11)
            Bx = Bx * g_nm * 3/4*rt11*x*(29393*z**10 - 62985*z**8*r**2 + 46410*z**6*r**4 - 13650*z**4*r**6 + 1365*z**2*r**8 - 21*r**10)
            By = By * g_nm * 3/4*rt11*y*(29393*z**10 - 62985*z**8*r**2 + 46410*z**6*r**4 - 13650*z**4*r**6 + 1365*z**2*r**8 - 21*r**10)
            Bz = Bz * g_nm * 1/4*rt11*z*(88179*z**10 - 230945*z**8*r**2 + 218790*z**6*r**4 - 90090*z**4*r**6 + 15015*z**2*r**8 - 693*r**10)
        elif m==1:
            rt5 = sqrt(5)
            Bx = Bx * 1/2*rt5*( g_nm*z*(756*x**10 - 63*y**10 + 777*y**8*z**2 - 1176*y**6*z**4 - 864*y**4*z**6 + 1024*y**2*z**8 - 128*z**10 + 21*x**8*(141*y**2 - 587*z**2) + 84*x**6*(51*y**4 - 431*y**2*z**2 + 454*z**4) + 18*x**4*(147*y**6 - 1925*y**4*z**2 + 4172*y**2*z**4 - 1712*z**6) + 4*x**2*(126*y**8 - 2499*y**6*z**2 + 8946*y**4*z**4 - 7920*y**2*z**6 + 1504*z**8)) + h_nm*39*x*y*z*(21*x**8 + 21*y**8 - 336*y**6*z**2 + 1008*y**4*z**4 - 768*y**2*z**6 + 128*z**8 + 84*x**6*(y**2 - 4*z**2) + 126*x**4*(y**4 - 8*y**2*z**2 + 8*z**4) + 12*x**2*(7*y**6 - 84*y**4*z**2 + 168*y**2*z**4 - 64*z**6)) )
            By = By * 1/2*rt5*( g_nm*39*x*y*z*(21*x**8 + 21*y**8 - 336*y**6*z**2 + 1008*y**4*z**4 - 768*y**2*z**6 + 128*z**8 + 84*x**6*(y**2 - 4*z**2) + 126*x**4*(y**4 - 8*y**2*z**2 + 8*z**4) + 12*x**2*(7*y**6 - 84*y**4*z**2 + 168*y**2*z**4 - 64*z**6)) + h_nm*-z*(63*x**10 - 756*y**10 + 12327*y**8*z**2 - 38136*y**6*z**4 + 30816*y**4*z**6 - 6016*y**2*z**8 + 128*z**10 - 21*x**8*(24*y**2 + 37*z**2) - 294*x**6*(9*y**4 - 34*y**2*z**2 - 4*z**4) - 18*x**4*(238*y**6 - 1925*y**4*z**2 + 1988*y**2*z**4 - 48*z**6) + x**2*(-2961*y**8 + 36204*y**6*z**2 - 75096*y**4*z**4 + 31680*y**2*z**6 - 1024*z**8)) )
            Bz = Bz * -3/2*rt5*( g_nm*x*(21*x**10 + 21*y**10 - 1260*y**8*z**2 + 8400*y**6*z**4 - 13440*y**4*z**6 + 5760*y**2*z**8 - 512*z**10 + 105*x**8*(y**2 - 12*z**2) + 210*x**6*(y**4 - 24*y**2*z**2 + 40*z**4) + 210*x**4*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6) + 15*x**2*(7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8)) + h_nm*y*(21*x**10 + 21*y**10 - 1260*y**8*z**2 + 8400*y**6*z**4 - 13440*y**4*z**6 + 5760*y**2*z**8 - 512*z**10 + 105*x**8*(y**2 - 12*z**2) + 210*x**6*(y**4 - 24*y**2*z**2 + 40*z**4) + 210*x**4*(y**6 - 36*y**4*z**2 + 120*y**2*z**4 - 64*z**6) + 15*x**2*(7*y**8 - 336*y**6*z**2 + 1680*y**4*z**4 - 1792*y**2*z**6 + 384*z**8)) )
        elif m==2:
            rt15 = sqrt(15)
            Bx = Bx * 1/2*rt15*( g_nm*1/2*x*(77*x**10 + 7*x**8*(29*y**2 - 634*z**2) + 42*x**6*(y**4 - 180*y**2*z**2 + 664*z**4) - 14*x**4*(23*y**6 - 282*y**4*z**2 - 1608*y**2*z**4 + 2896*z**6) + x**2*(-343*y**8 + 12824*y**6*z**2 - 38640*y**4*z**4 + 448*y**2*z**6 + 14464*z**8) - 3*(35*y**10 - 1918*y**8*z**2 + 11088*y**6*z**4 - 13664*y**4*z**6 + 2944*y**2*z**8 + 256*z**10)) + h_nm*-y*(-84*x**10 + 7*y**10 - 329*y**8*z**2 + 1344*y**6*z**4 - 112*y**4*z**6 - 1408*y**2*z**8 + 384*z**10 - 7*x**8*(47*y**2 - 681*z**2) - 28*x**6*(17*y**4 - 499*y**2*z**2 + 1044*z**4) - 42*x**4*(7*y**6 - 317*y**4*z**2 + 1360*y**2*z**4 - 968*z**6) - 4*x**2*(14*y**8 - 945*y**6*z**2 + 6636*y**4*z**4 - 10136*y**2*z**6 + 3264*z**8)) )
            By = By * 1/2*rt15*( g_nm*1/2*y*(105*x**10 - 77*y**10 + 4438*y**8*z**2 - 27888*y**6*z**4 + 40544*y**4*z**6 - 14464*y**2*z**8 + 768*z**10 + 7*x**8*(49*y**2 - 822*z**2) + 14*x**6*(23*y**4 - 916*y**2*z**2 + 2376*z**4) - 42*x**4*(y**6 + 94*y**4*z**2 - 920*y**2*z**4 + 976*z**6) + x**2*(-203*y**8 + 7560*y**6*z**2 - 22512*y**4*z**4 - 448*y**2*z**6 + 8832*z**8)) + h_nm*-x*(7*x**10 - 84*y**10 + 4767*y**8*z**2 - 29232*y**6*z**4 + 40656*y**4*z**6 - 13056*y**2*z**8 + 384*z**10 - 7*x**8*(8*y**2 + 47*z**2) - 42*x**6*(7*y**4 - 90*y**2*z**2 - 32*z**4) - 14*x**4*(34*y**6 - 951*y**4*z**2 + 1896*y**2*z**4 + 8*z**6) + x**2*(-329*y**8 + 13972*y**6*z**2 - 57120*y**4*z**4 + 40544*y**2*z**6 - 1408*z**8)) )
            Bz = Bz * 39/2*rt15*( g_nm*1/2*(x**2 - y**2)*z*(21*x**8 + 21*y**8 - 336*y**6*z**2 + 1008*y**4*z**4 - 768*y**2*z**6 + 128*z**8 + 84*x**6*(y**2 - 4*z**2) + 126*x**4*(y**4 - 8*y**2*z**2 + 8*z**4) + 12*x**2*(7*y**6 - 84*y**4*z**2 + 168*y**2*z**4 - 64*z**6)) + h_nm*x*y*z*(21*x**8 + 21*y**8 - 336*y**6*z**2 + 1008*y**4*z**4 - 768*y**2*z**6 + 128*z**8 + 84*x**6*(y**2 - 4*z**2) + 126*x**4*(y**4 - 8*y**2*z**2 + 8*z**4) + 12*x**2*(7*y**6 - 84*y**4*z**2 + 168*y**2*z**4 - 64*z**6)) )
        elif m==3:
            rt195o2 = sqrt(195/2)
            Bx = Bx * -3*rt195o2*( g_nm*z*(28*x**10 - 7*x**8*(3*y**2 + 59*z**2) - 28*x**6*(8*y**4 - 23*y**2*z**2 - 39*z**4) - 2*x**4*(133*y**6 - 1225*y**4*z**2 + 1218*y**2*z**4 + 332*z**6) + y**2*(7*y**8 - 77*y**6*z**2 + 84*y**4*z**4 + 104*y**2*z**6 - 64*z**8) - 4*x**2*(21*y**8 - 329*y**6*z**2 + 861*y**4*z**4 - 420*y**2*z**6 - 16*z**8)) + h_nm*x*y*z*(91*x**8 - 49*y**8 + 644*y**6*z**2 - 1344*y**4*z**4 + 352*y**2*z**6 + 128*z**8 + 28*x**6*(8*y**2 - 47*z**2) + 14*x**4*(9*y**4 - 142*y**2*z**2 + 240*z**4) - 4*x**2*(14*y**6 + 7*y**4*z**2 - 504*y**2*z**4 + 472*z**6)) )
            By = By * 3*rt195o2*( g_nm*-x*y*z*(49*x**8 - 91*y**8 + 1316*y**6*z**2 - 3360*y**4*z**4 + 1888*y**2*z**6 - 128*z**8 + 28*x**6*(2*y**2 - 23*z**2) - 14*x**4*(9*y**4 - 2*y**2*z**2 - 96*z**4) - 4*x**2*(56*y**6 - 497*y**4*z**2 + 504*y**2*z**4 + 88*z**6)) + h_nm*z*(7*x**10 - 7*x**8*(12*y**2 + 11*z**2) - 14*x**6*(19*y**4 - 94*y**2*z**2 - 6*z**4) + x**4*(-224*y**6 + 2450*y**4*z**2 - 3444*y**2*z**4 + 104*z**6) + x**2*(-21*y**8 + 644*y**6*z**2 - 2436*y**4*z**4 + 1680*y**2*z**6 - 64*z**8) + y**2*(28*y**8 - 413*y**6*z**2 + 1092*y**4*z**4 - 664*y**2*z**6 + 64*z**8)) )
            Bz = Bz * 7*rt195o2*( g_nm*x*(x**2 - 3*y**2)*(x**8 + y**8 - 56*y**6*z**2 + 336*y**4*z**4 - 448*y**2*z**6 + 128*z**8 + 4*x**6*(y**2 - 14*z**2) + 6*x**4*(y**4 - 28*y**2*z**2 + 56*z**4) + 4*x**2*(y**6 - 42*y**4*z**2 + 168*y**2*z**4 - 112*z**6)) + h_nm*y*(3*x**2 - y**2)*(x**8 + y**8 - 56*y**6*z**2 + 336*y**4*z**4 - 448*y**2*z**6 + 128*z**8 + 4*x**6*(y**2 - 14*z**2) + 6*x**4*(y**4 - 28*y**2*z**2 + 56*z**4) + 4*x**2*(y**6 - 42*y**4*z**2 + 168*y**2*z**4 - 112*z**6)) )
        elif m==4:
            rt195 = sqrt(195)
            Bx = Bx * -rt195*( g_nm*1/2*x*(11*x**10 - x**8*(49*y**2 + 556*z**2) - 6*x**6*(31*y**4 - 480*y**2*z**2 - 476*z**4) - 14*x**4*(11*y**6 - 444*y**4*z**2 + 1164*y**2*z**4 + 208*z**6) + 3*y**2*(9*y**8 - 404*y**6*z**2 + 1624*y**4*z**4 - 672*y**2*z**6 - 448*z**8) + x**2*(-y**8 + 1568*y**6*z**2 - 14280*y**4*z**4 + 16576*y**2*z**6 + 448*z**8)) + h_nm*2*y*(12*x**10 + x**8*(19*y**2 - 597*z**2) - 14*x**6*(y**4 + 28*y**2*z**2 - 213*z**4) - 6*x**4*(6*y**6 - 161*y**4*z**2 + 105*y**2*z**4 + 476*z**6) + y**2*(y**8 - 41*y**6*z**2 + 126*y**4*z**4 + 56*y**2*z**6 - 112*z**8) + x**2*(-14*y**8 + 720*y**6*z**2 - 3486*y**4*z**4 + 2576*y**2*z**6 + 336*z**8)) )
            By = By * rt195*( g_nm*-1/2*y*(27*x**10 - x**8*(y**2 + 1212*z**2) + x**6*(-154*y**4 + 1568*y**2*z**2 + 4872*z**4) - 6*x**4*(31*y**6 - 1036*y**4*z**2 + 2380*y**2*z**4 + 336*z**6) + x**2*(-49*y**8 + 2880*y**6*z**2 - 16296*y**4*z**4 + 16576*y**2*z**6 - 1344*z**8) + y**2*(11*y**8 - 556*y**6*z**2 + 2856*y**4*z**4 - 2912*y**2*z**6 + 448*z**8)) + h_nm*2*x*(x**10 - x**8*(14*y**2 + 41*z**2) - 18*x**6*(2*y**4 - 40*y**2*z**2 - 7*z**4) - 14*x**4*(y**6 - 69*y**4*z**2 + 249*y**2*z**4 - 4*z**6) + x**2*(19*y**8 - 392*y**6*z**2 - 630*y**4*z**4 + 2576*y**2*z**6 - 112*z**8) + 3*y**2*(4*y**8 - 199*y**6*z**2 + 994*y**4*z**4 - 952*y**2*z**6 + 112*z**8)) )
            Bz = Bz * -21*rt195*( g_nm*1/2*(x**4 - 6*x**2*y**2 + y**4)*z*(5*x**6 + 5*y**6 - 70*y**4*z**2 + 168*y**2*z**4 - 80*z**6 + 5*x**4*(3*y**2 - 14*z**2) + x**2*(15*y**4 - 140*y**2*z**2 + 168*z**4)) + h_nm*2*x*y*(x**2 - y**2)*z*(5*x**6 + 5*y**6 - 70*y**4*z**2 + 168*y**2*z**4 - 80*z**6 + 5*x**4*(3*y**2 - 14*z**2) + x**2*(15*y**4 - 140*y**2*z**2 + 168*z**4)) )
        elif m==5:
            rt39o2 = sqrt(39/2)
            Bx = Bx * rt39o2*( g_nm*z*(180*x**10 - 5*x**8*(363*y**2 + 419*z**2) - 28*x**6*(90*y**4 - 785*y**2*z**2 - 131*z**4) + 210*x**4*(5*y**6 + 35*y**4*z**2 - 178*y**2*z**4 - 4*z**6) - 5*y**4*(15*y**6 - 125*y**4*z**2 + 28*y**2*z**4 + 168*z**6) + 20*x**2*(75*y**8 - 805*y**6*z**2 + 987*y**4*z**4 + 252*y**2*z**6)) + h_nm*x*y*z*(975*x**8 + 555*y**8 - 5220*y**6*z**2 + 4368*y**4*z**4 + 3360*y**2*z**6 - 300*x**6*(2*y**2 + 37*z**2) - 70*x**4*(51*y**4 - 230*y**2*z**2 - 264*z**4) - 20*x**2*(72*y**6 - 1099*y**4*z**2 + 1904*y**2*z**4 + 168*z**6)) )
            By = By * rt39o2*( g_nm*x*y*z*(555*x**8 - 180*x**6*(8*y**2 + 29*z**2) + x**4*(-3570*y**4 + 21980*y**2*z**2 + 4368*z**4) - 20*x**2*(30*y**6 - 805*y**4*z**2 + 1904*y**2*z**4 - 168*z**6) + 15*(65*y**8 - 740*y**6*z**2 + 1232*y**4*z**4 - 224*y**2*z**6)) + h_nm*z*(-75*x**10 + 125*x**8*(12*y**2 + 5*z**2) + 70*x**6*(15*y**4 - 230*y**2*z**2 - 2*z**4) + y**4*(180*y**6 - 2095*y**4*z**2 + 3668*y**2*z**4 - 840*z**6) - 210*x**4*(12*y**6 - 35*y**4*z**2 - 94*y**2*z**4 + 4*z**6) - 5*x**2*(363*y**8 - 4396*y**6*z**2 + 7476*y**4*z**4 - 1008*y**2*z**6)) )
            Bz = Bz * -3*rt39o2*( g_nm*x*(x**4 - 10*x**2*y**2 + 5*y**4)*(5*x**6 + 5*y**6 - 240*y**4*z**2 + 1120*y**2*z**4 - 896*z**6 + 15*x**4*(y**2 - 16*z**2) + 5*x**2*(3*y**4 - 96*y**2*z**2 + 224*z**4)) + h_nm*y*(5*x**4 - 10*x**2*y**2 + y**4)*(5*x**6 + 5*y**6 - 240*y**4*z**2 + 1120*y**2*z**4 - 896*z**6 + 15*x**4*(y**2 - 16*z**2) + 5*x**2*(3*y**4 - 96*y**2*z**2 + 224*z**4)) )
        elif m==6:
            rt195o2 = sqrt(195/2)
            Bx = Bx * 3/2*rt195o2*( g_nm*1/2*x*(11*x**10 - x**8*(179*y**2 + 426*z**2) + x**6*(-106*y**4 + 7080*y**2*z**2 + 1376*z**4) + 14*x**4*(23*y**6 - 186*y**4*z**2 - 1584*y**2*z**4 - 32*z**6) - y**4*(47*y**6 - 1542*y**4*z**2 + 2912*y**2*z**4 + 2240*z**6) + x**2*(191*y**8 - 8568*y**6*z**2 + 25760*y**4*z**4 + 4480*y**2*z**6)) + h_nm*-y*(-36*x**10 + x**8*(83*y**2 + 1371*z**2) + 28*x**6*(7*y**4 - 153*y**2*z**2 - 152*z**4) + y**4*(3*y**6 - 93*y**4*z**2 + 128*y**2*z**4 + 224*z**6) + 2*x**4*(y**6 - 1491*y**4*z**2 + 7840*y**2*z**4 + 560*z**6) - 4*x**2*(18*y**8 - 645*y**6*z**2 + 1512*y**4*z**4 + 560*y**2*z**6)) )
            By = By * 3/2*rt195o2*( g_nm*1/2*y*(47*x**10 - 11*y**10 + 426*y**8*z**2 - 1376*y**6*z**4 + 448*y**4*z**6 - x**8*(191*y**2 + 1542*z**2) + x**6*(-322*y**4 + 8568*y**2*z**2 + 2912*z**4) + 2*x**4*(53*y**6 + 1302*y**4*z**2 - 12880*y**2*z**4 + 1120*z**6) + x**2*(179*y**8 - 7080*y**6*z**2 + 22176*y**4*z**4 - 4480*y**2*z**6)) + h_nm*-x*(3*x**10 - 36*y**10 + 1371*y**8*z**2 - 4256*y**6*z**4 + 1120*y**4*z**6 - 3*x**8*(24*y**2 + 31*z**2) + 2*x**6*(y**4 + 1290*y**2*z**2 + 64*z**4) + 14*x**4*(14*y**6 - 213*y**4*z**2 - 432*y**2*z**4 + 16*z**6) + x**2*(83*y**8 - 4284*y**6*z**2 + 15680*y**4*z**4 - 2240*y**2*z**6)) )
            Bz = Bz * 17/2*rt195o2*( g_nm*1/2*(x**6 - 15*x**4*y**2 + 15*x**2*y**4 - y**6)*z*(15*x**4 + 15*y**4 - 160*y**2*z**2 + 224*z**4 + 10*x**2*(3*y**2 - 16*z**2)) + h_nm*x*y*(3*x**4 - 10*x**2*y**2 + 3*y**4)*z*(15*x**4 + 15*y**4 - 160*y**2*z**2 + 224*z**4 + 10*x**2*(3*y**2 - 16*z**2)) )
        elif m==7:
            rt3315o2 = sqrt(3315/2)
            Bx = Bx * 1/2*rt3315o2*( g_nm*z*(-36*x**10 + x**8*(867*y**2 + 251*z**2) - 28*x**6*(39*y**4 + 211*y**2*z**2 + 4*z**4) + 7*y**6*(-3*y**4 + 13*y**2*z**2 + 16*z**4) - 42*x**4*(31*y**6 - 285*y**4*z**2 - 40*y**2*z**4) + 28*x**2*(24*y**8 - 131*y**6*z**2 - 60*y**4*z**4)) + h_nm*x*y*z*(-273*x**8 + 183*y**8 - 888*y**6*z**2 - 672*y**4*z**4 + 84*x**6*(17*y**2 + 22*z**2) + 14*x**4*(15*y**4 - 764*y**2*z**2 - 48*z**4) + x**2*(-1308*y**6 + 8456*y**4*z**2 + 2240*y**2*z**4)) )
            By = By * 1/2*rt3315o2*( g_nm*x*y*z*(-183*x**8 + 12*x**6*(109*y**2 + 74*z**2) - 14*x**4*(15*y**4 + 604*y**2*z**2 - 48*z**4) - 28*x**2*(51*y**6 - 382*y**4*z**2 + 80*y**2*z**4) + 21*(13*y**8 - 88*y**6*z**2 + 32*y**4*z**4)) + h_nm*z*(21*x**10 - 7*x**8*(96*y**2 + 13*z**2) + 14*x**6*(93*y**4 + 262*y**2*z**2 - 8*z**4) + y**6*(36*y**4 - 251*y**2*z**2 + 112*z**4) + 42*x**4*(26*y**6 - 285*y**4*z**2 + 40*y**2*z**4) + x**2*(-867*y**8 + 5908*y**6*z**2 - 1680*y**4*z**4)) )
            Bz = Bz * 3/2*rt3315o2*( g_nm*x*(x**6 - 21*x**4*y**2 + 35*x**2*y**4 - 7*y**6)*(x**4 + y**4 - 36*y**2*z**2 + 96*z**4 + 2*x**2*(y**2 - 18*z**2)) + h_nm*-y*(-7*x**6 + 35*x**4*y**2 - 21*x**2*y**4 + y**6)*(x**4 + y**4 - 36*y**2*z**2 + 96*z**4 + 2*x**2*(y**2 - 18*z**2)) )
        elif m==8:
            rt1105 = sqrt(1105)
            Bx = Bx * rt1105*( g_nm*1/4*x*(-11*x**10 + x**8*(361*y**2 + 244*z**2) - 18*x**6*(47*y**4 + 440*y**2*z**2 + 8*z**4) + 3*y**6*(-25*y**4 + 444*y**2*z**2 + 336*z**4) - 42*x**4*(11*y**6 - 588*y**4*z**2 - 72*y**2*z**4) + 3*x**2*(227*y**8 - 4816*y**6*z**2 - 1680*y**4*z**4)) + h_nm*-2*y*(12*x**10 - 3*x**8*(31*y**2 + 87*z**2) + y**6*(y**4 - 17*y**2*z**2 - 18*z**4) + 42*x**6*(y**4 + 52*y**2*z**2 + 3*z**4) + 18*x**4*(6*y**6 - 161*y**4*z**2 - 35*y**2*z**4) + x**2*(-38*y**8 + 720*y**6*z**2 + 378*y**4*z**4)) )
            By = By * rt1105*( g_nm*1/4*y*(-75*x**10 - 11*y**10 + 244*y**8*z**2 - 144*y**6*z**4 + 3*x**8*(227*y**2 + 444*z**2) - 42*x**6*(11*y**4 + 344*y**2*z**2 - 24*z**4) - 18*x**4*(47*y**6 - 1372*y**4*z**2 + 280*y**2*z**4) + x**2*(361*y**8 - 7920*y**6*z**2 + 3024*y**4*z**4)) + h_nm*2*x*(x**10 - x**8*(38*y**2 + 17*z**2) + 18*x**6*(6*y**4 + 40*y**2*z**2 - z**4) + 3*y**6*(4*y**4 - 87*y**2*z**2 + 42*z**4) + 42*x**4*(y**6 - 69*y**4*z**2 + 9*y**2*z**4) - 3*x**2*(31*y**8 - 728*y**6*z**2 + 210*y**4*z**4)) )
            Bz = Bz * -57*rt1105*( g_nm*1/4*(x**8 - 28*x**6*y**2 + 70*x**4*y**4 - 28*x**2*y**6 + y**8)*z*(x**2 + y**2 - 6*z**2) + h_nm*2*x*y*(x**6 - 7*x**4*y**2 + 7*x**2*y**4 - y**6)*z*(x**2 + y**2 - 6*z**2) )
        elif m==9:
            rt20995o2 = sqrt(20995/2)
            Bx = Bx * 3/2*rt20995o2*( g_nm*z*(4*x**10 - 3*y**8*(y**2 + z**2) + 84*x**6*y**2*(9*y**2 + z**2) - 3*x**8*(57*y**2 + z**2) - 42*x**4*(17*y**6 + 5*y**4*z**2) + 12*x**2*(12*y**8 + 7*y**6*z**2)) + h_nm*x*y*z*(39*x**8 + 31*y**8 + 24*y**6*z**2 - 12*x**6*(37*y**2 + 2*z**2) + 42*x**4*(21*y**4 + 4*y**2*z**2) - 12*x**2*(33*y**6 + 14*y**4*z**2)) )
            By = By * 3/2*rt20995o2*( g_nm*x*y*z*(31*x**8 + 39*y**8 - 24*y**6*z**2 + x**6*(-396*y**2 + 24*z**2) + 42*x**4*(21*y**4 - 4*y**2*z**2) + x**2*(-444*y**6 + 168*y**4*z**2)) + h_nm*-z*(3*x**10 - 4*y**10 + 3*y**8*z**2 + 3*x**8*(-48*y**2 + z**2) + 42*x**6*(17*y**4 - 2*y**2*z**2) - 42*x**4*(18*y**6 - 5*y**4*z**2) + 3*x**2*(57*y**8 - 28*y**6*z**2)) )
            Bz = Bz * -1/2*rt20995o2*( g_nm*x*(x**8 - 36*x**6*y**2 + 126*x**4*y**4 - 84*x**2*y**6 + 9*y**8)*(x**2 + y**2 - 20*z**2) + h_nm*y*(9*x**8 - 84*x**6*y**2 + 126*x**4*y**4 - 36*x**2*y**6 + y**8)*(x**2 + y**2 - 20*z**2) )
        elif m==10:
            rt4199o2 = sqrt(4199/2)
            Bx = Bx * 1/2*rt4199o2*( g_nm*1/2*x*(11*x**10 - 5*x**8*(119*y**2 + 2*z**2) - 3*y**8*(37*y**2 + 30*z**2) + 90*x**6*(39*y**4 + 4*y**2*z**2) - 210*x**4*(23*y**6 + 6*y**4*z**2) + 15*x**2*(113*y**8 + 56*y**6*z**2)) + h_nm*y*(60*x**10 - 5*y**8*(y**2 + z**2) - 15*x**8*(59*y**2 + 3*z**2) + 84*x**6*(29*y**4 + 5*y**2*z**2) - 90*x**4*(19*y**6 + 7*y**4*z**2) + 20*x**2*(14*y**8 + 9*y**6*z**2)) )
            By = By * 1/2*rt4199o2*( g_nm*1/2*y*(111*x**10 - 11*y**10 + 10*y**8*z**2 + x**8*(-1695*y**2 + 90*z**2) + 210*x**6*(23*y**4 - 4*y**2*z**2) - 90*x**4*(39*y**6 - 14*y**4*z**2) + 5*x**2*(119*y**8 - 72*y**6*z**2)) + h_nm*-x*(5*x**10 - 60*y**10 + 45*y**8*z**2 + 5*x**8*(-56*y**2 + z**2) + 90*x**6*(19*y**4 - 2*y**2*z**2) - 42*x**4*(58*y**6 - 15*y**4*z**2) + 15*x**2*(59*y**8 - 28*y**6*z**2)) )
            Bz = Bz * 21/2*rt4199o2*( g_nm*1/2*(x**10 - 45*x**8*y**2 + 210*x**6*y**4 - 210*x**4*y**6 + 45*x**2*y**8 - y**10)*z + h_nm*x*y*(5*x**8 - 60*x**6*y**2 + 126*x**4*y**4 - 60*x**2*y**6 + 5*y**8)*z )
        else:
            print(" m = ", m)
            raise ValueError(f"In field_xyz.eval_Bi_Schmidt, n={n} and m is not between 0 and n.")

    else:
        print(" n = ", n)
        raise ValueError("In field_xyz.eval_Bi_Schmidt, n>10 but only n=1 to n=10 are supported.")

    Bx = Bx * timeRot
    By = By * timeRot
    Bz = Bz * timeRot

    return Bx, By, Bz
