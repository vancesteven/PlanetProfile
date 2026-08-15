"""Double-precision ALMA backend vendored for PlanetProfile gravity calculations."""


import numpy as np
import math
import time
from scipy.linalg import lu_solve, lu_factor
from scipy.special import gamma, binom

def initialize(ndigits):
    # dummy function for compatibility with the MP version
    return

def build_model(r_in, rho_in, mu_in, eta_in, rheology, params, norm=True):
    """
    build_model(r_in, rho_in, mu_in, eta_in, rheology, params)
    """

    # Define the Newton constant

    Gnwt   = 6.674e-11

    # Check argument types

    for arg in [r_in, rho_in, mu_in, eta_in, params]:
        if type(arg) not in [list, tuple, np.ndarray]:
            raise TypeError("Invalid argument type")

    if type(rheology) not in [list, tuple]:
        raise TypeError("Invalid argument type")

    # Check argument sizes

    nla = len( rho_in )
    for arg in [r_in, rho_in, mu_in, eta_in, rheology, params]:
        if( len(arg) != nla ):
            raise ValueError("Argument size mismatch")

    # Set model parameters

    r     = np.zeros( nla+1 )
    r[1:] = np.array( r_in )

    rho    = np.array( rho_in )
    mu     = np.array( mu_in )
    eta    = np.array( eta_in )
    rpar   = np.array( params, dtype=float, copy=True )

    # Parse rheology

    rheol = []
    for rcode in rheology:
        if type(rcode)==str:
            if rcode.lower()=='fluid':
                rheol.append(0)
            elif rcode.lower()=='elastic':
                rheol.append(1)
            elif rcode.lower()=='maxwell':
                rheol.append(2)
            elif rcode.lower()=='newton':
                rheol.append(3)
            elif rcode.lower()=='kelvin':
                rheol.append(4)
            elif rcode.lower()=='burgers':
                rheol.append(5)
            elif rcode.lower()=='andrade':
                rheol.append(6)
            elif rcode.lower() == 'andrade_ext':
                rheol.append(7)
            elif rcode.lower() == 'sundcoop':
                rheol.append(8)
            else:
                raise ValueError('Unknown rheology name "'+rcode+'"')
        elif type(rcode)==int:
            if( (rcode>=0) & (rcode<=8) ):
                rheol.append(rcode)
            else:
                raise ValueError('Invalid rheology code '+str(rcode))
        else:
            raise TypeError('Invalid rheology argument')

    # Compute gamma(alpha+1) for Andrade layers unless the caller supplied
    # the second Andrade parameter explicitly.

    for i in range(nla):
        if rheol[i]==6:
            if not np.isfinite(rpar[i,0]):
                raise ValueError("Andrade alpha must be finite.")
            if not np.isfinite(rpar[i,1]) or rpar[i,1] == 0.0:
                rpar[i,1] = gamma(rpar[i,0] + 1)
            elif rpar[i,1] < 0.0:
                raise ValueError("Andrade gamma must be positive when supplied.")

    if any(rcode == 0 for rcode in rheol[1:]):
        raise ValueError(
            "DP ALMA supports 'fluid' rheology only for the innermost layer; "
            "use 'newton' for ocean or fluid shell layers."
        )

    # Build model

    gra    = np.zeros(nla+1)
    mlayer = np.zeros(nla)

    # ----- Compute mass of the planet

    for i in range(nla):
        mlayer[i] = rho[i] * ( r[i+1]**3 - r[i]**3 )
    mlayer = mlayer * 4/3 * np.pi

    mass = sum( mlayer )

    # ----- Compute gravity at the interface boundaries

    for i in range(1,nla+1):
        gra[i] = Gnwt * sum( mlayer[0:i] ) / r[i]**2

    # Normalization

    # ----- Define reference scales

    r0    = r[nla]
    rho0  = max(rho)
    mu0   = max(mu)
    t0    = 1000.0 * 365.25 * 24.0 * 3600.0
    eta0  = mu0 * t0
    mass0 = rho0 * r0**3
    G = Gnwt

    if norm:
        # ------- Normalize model parameters

        r      = r   / r0
        rho    = rho / rho0
        mu     = mu  / mu0
        eta    = eta / eta0

        # ------- Normalize Newton constant

        G      = Gnwt * rho0**2 * r0**2 / mu0

        # ------- Normalized mass

        mass   = mass / mass0
        mlayer = mlayer / mass0

        # ------- Normalized gravity at internal interfaces

        gra    = gra * r0**2 * (G/Gnwt) / mass0

    # ------- Build the model vector

    model = (r, rho, mu, eta, rheol, rpar, r0, rho0, mu0, eta0, t0, mass0, gra, mass, mlayer, G)

    return model

def build_profile(r_profile,model):

    # Normalizes the internal sampling radii, identifies the
    # corresponding layers and computes gravity at the sampling points

    (r, rho, mu, eta, rheol, rpar, r0, rho0, mu0, eta0, t0, mass0, gra, mass, mlayer, G) = model

    r_profile   = np.array(r_profile) / float(r0)
    idx_profile = np.zeros( len(r_profile), dtype=int )

    nla = len(r)-1

    for i in range(1,nla):
        mask = ( r_profile <= r[i+1] ) & ( r_profile >= r[i] )
        idx_profile[mask] = i

    idx_profile = idx_profile.tolist()

    g_profile = np.zeros( len(r_profile) )

    for i in range(len(r_profile)):
        i_layer = idx_profile[i]
        mass =  sum( mlayer[0:i_layer] )
        mass += 4/3 * np.pi * rho[i_layer] * (r_profile[i]**3 - r[i_layer]**3 )
        g_profile[i] = G * mass / r_profile[i]**2

    return r_profile,g_profile,idx_profile

def love_numbers(model,degrees,timesteps,loadtype,mode,order,verbose=False,internal=False,r_profile=[]):
    """
    love_numbers(degrees,timesteps,loadtype,loadfcn,tau,mode,order)

    timesteps are ALMA time values in kyr. Callers with periods in seconds
    should divide by 1000 Julian years before calling this function.
    """

    if loadtype.lower()=='tidal':
        iload=0
    elif loadtype.lower()=='loading':
        iload=1
    else:
        raise ValueError( 'Unknown load type "'+loadtype+'"' )

    if mode.lower()   in [ 'heaviside', 'timedomain' ]:
        itype = 1
    elif mode.lower() in [ 'periodic', 'frequency' ]:
        itype = 2
    else:
        raise ValueError( 'Unknown mode option "'+mode+'"' )

    nt   = len(timesteps)
    ndeg = len(degrees)

    # Check that the sampling radius is of the correct type

    if type(r_profile) not in [list, tuple, np.ndarray]:
        raise TypeError("Invalid argument type for r_profile")

    # Tidal LNs of degree 1 are not allowed

    if (1 in degrees) & (iload==0):
        raise ValueError("Tidal LNs of degree 1 are not allowed")

    # ------- Allocate output arrays

    if internal:
        if len(r_profile)==0:
            nla = len(model[0]) - 1
            arr_size = (ndeg,nt,nla)
            idx_profile=[]
            g_profile=[]
        else:
            arr_size = (ndeg,nt,len(r_profile))
            r_profile, g_profile, idx_profile = build_profile(r_profile,model)
    else:
        arr_size = (ndeg,nt)
        idx_profile=[]
        g_profile=[]

    if itype==1:
        h_love = np.zeros( arr_size, dtype=float )
        l_love = np.zeros( arr_size, dtype=float )
        k_love = np.zeros( arr_size, dtype=float )
    elif itype==2:
        h_love = np.zeros( arr_size, dtype=complex )
        l_love = np.zeros( arr_size, dtype=complex )
        k_love = np.zeros( arr_size, dtype=complex )

    # ------- Compute LNs

    t1 = time.perf_counter()

    idx_n = 0
    idx_t = 0

    if ( itype==1 ):
        zeta = salzer_weights(order)
        for idx_n in range(ndeg):
            n = degrees[idx_n]
            for idx_t in range(nt):
                t = timesteps[idx_t]
                f = np.log(2) / t
                if internal:
                    if len(r_profile)==0:
                        nla = len(model[1])
                        xh = np.zeros(nla)
                        xl = np.zeros(nla)
                        xk = np.zeros(nla)
                    else:
                        xh = np.zeros(len(r_profile))
                        xl = np.zeros(len(r_profile))
                        xk = np.zeros(len(r_profile))
                else:
                    xh = 0.
                    xl = 0.
                    xk = 0.
                for ik in range(1,2*order+1):
                    s = f * ik
                    hh,ll,kk = love_numbers_sampler(model,n,s,iload,internal,r_profile,g_profile,idx_profile)
                    #
                    xh += np.real(hh) * zeta[ik-1] * f / s
                    xl += np.real(ll) * zeta[ik-1] * f / s
                    xk += np.real(kk) * zeta[ik-1] * f / s
                    #
                if internal:
                    h_love[ idx_n, idx_t, : ] = xh
                    l_love[ idx_n, idx_t, : ] = xl
                    k_love[ idx_n, idx_t, : ] = xk
                else:
                    h_love[ idx_n, idx_t ] = xh
                    l_love[ idx_n, idx_t ] = xl
                    k_love[ idx_n, idx_t ] = xk
            if verbose:
                t2 = time.perf_counter()
                print( "Harmonic degree n = " + str(n) + " ( " + str(t2-t1) + " s )" )
                t1 = t2
    elif itype==2:
        for idx_n in range(ndeg):
            for idx_t in range(nt):
                n = degrees[idx_n]
                t = timesteps[idx_t]
                omega = 2.0 * np.pi / t
                s     = 1j * omega
                if internal:
                    hh, ll, kk, y = love_numbers_sampler(model,n,s,iload,internal,
                                                       r_profile,g_profile,
                                                       idx_profile)
                else:
                    hh, ll, kk = love_numbers_sampler(model,n,s,iload,internal,
                                        r_profile,g_profile,
                                        idx_profile)
                if internal:
                    h_love[ idx_n, idx_t, : ] = hh
                    l_love[ idx_n, idx_t, : ] = ll
                    k_love[ idx_n, idx_t, : ] = kk
                else:
                    h_love[ idx_n, idx_t ] = hh
                    l_love[ idx_n, idx_t ] = ll
                    k_love[ idx_n, idx_t ] = kk
            if verbose:
                t2 = time.perf_counter()
                print( "Harmonic degree n = " + str(n) + " ( " + str(t2-t1) + " s )" )
                t1 = t2

    if internal:
        return h_love, l_love, k_love, y
    else:
        return h_love, l_love, k_love


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#
# subroutine 'fluid_core_bc'
# computes the boundary conditions at the CMB for a fluid inviscid core
#
# Initial version DM February 24, 2020
# Modified by DM June 16, 2020 - Converted to type(zm) for complex LNs
#
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def fluid_core_bc(n,r,rho,gra,G):
    """
    fluid_core_bc(n,r,rho,gra)
    """

    assert isinstance(n, (int,np.integer)), "harmonic degree shall be integer"

    b = np.zeros( (6,3) )

    b[0,0] = - 1 / gra * (r**n)
    b[0,2] =   1

    b[1,1] =   1

    b[2,2] =   rho * gra

    b[4,0] =   r**n

    b[5,0] =   2 * ( n - 1 ) * r**(n-1)
    b[5,2] =   4 * np.pi * G * rho

    return b



# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#
# subroutine 'surface_bc'
# computes the boundary conditions at the surface
#
# Initial version DM February 24, 2020
# Implementation of the degree 1 - DM November 29, 2022
#
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def surface_bc(n,r,gra,iload,G):
    """
    surface_bc(n,r,gra,iload)
    """

    assert isinstance(n, (int,np.integer)), "harmonic degree shall be integer"

    bs = np.zeros(3)
    kappa = ( 2 * n + 1) / ( 4 * np.pi * r**2 )

    if n!=1:
        if iload==1:
            bs[0] = - gra * kappa
        bs[2] = -4 * np.pi * G * kappa
    else:
        if iload==1:
            bs[0] = -gra * kappa
        else:
            raise ValueError("Tidal forcing is not defined for n=1")

    return bs



# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#
# subroutine 'salzer_weights'
# computes the weights zeta(i) to be used in the Salzer accelerated
# Post-Widder Laplace inversions scheme
#
# Initial version DM February 24, 2020
#
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+


def salzer_weights(order):
    """
    salzer_weights(order)
    """

    if( type(order) != int ):
        raise TypeError("salzer_weights(): order must be integer")

    m = order
    zeta = np.zeros( 2*m )

    for k in range(1,2*m+1):
        j1 = math.floor( (k+1)/2 )
        j2 = min( k, m )
        for j in range(j1,j2+1):
            fattm = gamma(m+1)
            q1 = binom(m, j)
            q2 = binom(2*j, j)
            q3 = binom(j, k-j)
            zeta[k-1] = zeta[k-1] + j**(m+1) / fattm * q1 * q2 * q3
        if (m+k)%2 != 0:
            zeta[k-1] = -zeta[k-1]

    return zeta

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#
# subroutine 'complex_rigidity'
# computes mu(s) for various rheologies
#
# Initial version DM February 24, 2020
# Modified by DM June 11, 2020 - Burgers and Andrade rheologies
# Modified by DM June 16, 2020 - Complex LNs  - converted to type(zm)
#
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#
#
#
def complex_rigidity(s,mu,eta,code,par):
    """
    complex_rigidity(s,mu,eta,code,par)
    """
    if code==1:                     # Elastic
        mu_s = mu
    elif code==2:                   # Maxwell
        mu_s = mu * s / ( s + mu/eta )
    elif code==3:                   # Newton
        mu_s = eta * s
    elif code==4:                   # Kelvin
        mu_s = mu + eta * s
    elif code==5:                   # Burgers
        mu2  = par[0] * mu
        eta2 = par[1] * eta
        mu_s = mu * s * ( s + mu2/eta2 ) / \
	       ( s**2 + s*(mu/eta + (mu+mu2)/eta2) + (mu*mu2)/(eta*eta2) )
    elif code==6:                   # Andrade
        alpha = par[0]
        gam   = par[1]
        mu_s  = 1/mu + 1/(eta*s) + gam * (1/mu) * (s*eta/mu)**(-alpha)
        mu_s  = 1/mu_s
    elif code == 7:                   # Extended Andrade model (zeta != 1)
        alpha = par[0]
        zeta = par[1]
        mu_s = mu * s / (s + mu / eta + s * mu ** alpha * gamma(alpha + 1) * (s * eta * zeta) ** (-alpha))
    elif code == 8:
        alpha = par[0]
        zeta = par[1]
        mu2  = 5 * mu      # shear modulus defect
        eta2 = 0.02 * eta  # Voigt-Kelvin viscosity

        and_fact = mu ** alpha * gamma(alpha + 1) * (s * eta * zeta) ** (-alpha)

        num = mu * s * (mu2 / eta2 + s)
        den1 = s ** 2 * (1 + and_fact)
        den2 = s * (mu / eta + (mu + mu2) / eta2 + mu2 / eta2 * and_fact) + (mu * mu2) / (eta * eta2)

        mu_s = num / (den1 + den2)
    else:
        raise ValueError("Invalid rheology code: "+str(code))

    return mu_s

# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#
# subroutine 'direct_matrix'
# computes the fundamental matrix in a layer
#
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def direct_matrix(n,r,rho,mu,gra,G):
    """
    direct_matrix(n,r,rho,mu,gra,G)
    computes the fundamental matrix in a given layer
    """

    assert isinstance(n, (int,np.integer)), "harmonic degree shall be integer"

    a1 = 2*n + 3             # 2n+3
    a2 = n + 1               # n+1
    a3 = n + 3               # n+3
    a4 = n**2 - n - 3        # n^2-n-3
    a5 = n + 2               # n+2
    a6 = n - 1               # n-1
    a7 = 2 * n + 1           # 2n+1
    a8 = 2 * n - 1           # 2n-1
    a9 = 2 - n               # 2-n
    a10= n**2 + 3 * n - 1    # n^2+3n-1
    a11= n**2 - 1            # n^2-1

    #Y = mp.matrix(6)
    #Y = np.full( (6, 6), 0.0 )
    Y = np.zeros( (6, 6), dtype=complex )

    Y[0,0] = n/(2 * a1) * r**(n+1)
    Y[1,0] = a3 / ( 2 * a1 * a2 ) * r**(n+1)
    Y[2,0] = ( n * rho * gra * r + 2 * a4 * mu ) / ( 2 * a1 ) * r**n
    Y[3,0] = n * a5 / ( a1 * a2 ) * mu * r**n
    Y[5,0] = 2 * np.pi * G * rho * n / a1 * r**(n+1)

    Y[0,1] = r**(n-1)
    Y[1,1] = 1 / n * r**(n-1)
    Y[2,1] = ( rho * gra * r + 2 * a6 * mu ) * r**(n-2)
    Y[3,1] = 2 * a6 / n * mu * r**(n-2)
    Y[5,1] = 4 * np.pi * G * rho * r**(n-1)

    Y[2,2] = rho * r**n
    Y[4,2] = r**n
    Y[5,2] = a7 * r**(n-1)

    Y[0,3] = a2 / ( 2 * a8 ) * r**(-n)
    Y[1,3] = a9 / ( 2 * n * a8 ) * r**(-n)
    Y[2,3] = ( a2 * rho * gra * r - 2 * a10 * mu ) / ( 2 * a8 ) * r**(-n-1)
    Y[3,3] = a11 / ( n * a8 ) * mu * r**(-n-1)
    Y[5,3] = 2 * np.pi * G * rho * a2 / a8 * r**(-n)

    Y[0,4] = r**(-n-2)
    Y[1,4] = - 1.0 / a2 * r**(-n-2)
    Y[2,4] = ( rho * gra * r - 2 * a5 * mu ) * r**(-n-3)
    Y[3,4] = 2 * a5 / a2 * mu * r**(-n-3)
    Y[5,4] = 4 * np.pi * G * rho * r**(-n-2)

    Y[2,5] = rho * r**(-n-1)
    Y[4,5] = r**(-n-1)

    return Y


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#
# subroutine 'inverse_matrix'
# computes the inverse of the fundamental matrix in a layer
#
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def inverse_matrix(n,r,rho,mu,gra,G):
    """
    inverse_matrix(n,r,rho,mu,gra,G)
    computes the fundamental matrix in a given layer
    """

    assert isinstance(n, (int,np.integer)), "harmonic degree shall be integer"

    a1 = 2 * n + 1           # 2n+1
    a2 = 2 * n - 1           # 2n-1
    a3 = n + 1               # n+1
    a4 = n - 1               # n-1
    a5 = n + 2               # n+2
    a6 = n + 3               # n+3
    a7 = n**2 + 3 * n - 1    # n^2+3n-1
    a8 = n**2 - n - 3        # n^2-n-3
    a9 = 2 - n               # 2-n
    a10= n**2 - 1            # n^2-1
    a11= 2 * n + 3           # 2n+3

    Yinv = np.zeros((6,6), dtype=complex)
    D    = np.zeros((6,6))
    #Yinv = np.full( (6,6), 0.0 )
    #D    = np.full( (6,6), 0.0 )

    D[0,0] = a3 * r**(-(n+1))
    D[1,1] = a3 * n / ( 2 * a2 ) * r**(-(n-1))
    D[2,2] = -r**(-(n-1))
    D[3,3] = n * r**n
    D[4,4] = n * a3 / ( 2 * a11 ) * r**(n+2)
    D[5,5] = r**(n+1)

    D = D / a1

    Yinv[0,0] =   rho * gra * r / mu - 2 * a5
    Yinv[1,0] = - rho * gra * r / mu + 2 * a7 / a3
    Yinv[2,0] =   4 * np.pi * G * rho
    Yinv[3,0] =   rho * gra * r / mu + 2 * a4
    Yinv[4,0] = - rho * gra * r / mu - 2 * a8 / n
    Yinv[5,0] =   4 * np.pi * G * rho * r

    Yinv[0,1] =   2 * n * a5
    Yinv[1,1] = - 2 * a10
    Yinv[3,1] =   2 * a10
    Yinv[4,1] = - 2 * n * a5

    Yinv[0,2] = - r / mu
    Yinv[1,2] =   r / mu
    Yinv[3,2] = - r / mu
    Yinv[4,2] =   r / mu

    Yinv[0,3] =   n * r / mu
    Yinv[1,3] =   a9 * r / mu
    Yinv[3,3] = - a3 * r / mu
    Yinv[4,3] =   a6 * r / mu

    Yinv[0,4] =   rho * r / mu
    Yinv[1,4] = - rho * r / mu
    Yinv[3,4] =   rho * r / mu
    Yinv[4,4] = - rho * r / mu
    Yinv[5,4] =   a1

    Yinv[2,5] = - 1
    Yinv[5,5] = - r

    Yinv = np.matmul(D,Yinv)

    return Yinv



# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#
# subroutine 'love_numbers'
# computes the Love numbers h,l,k in the Laplace domain for
# a given value of s
#
# Initial version DM February 24, 2020
# Modified DM June 11, 2020 - Burgers and Andrade rheologies
# Modified DM June 16, 2020 - Complex LNs
#                             (the 1/s factor is now outside this module)
# Fixed DM October 27, 2020 - Wrong call to complex_rigidty for core BC
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

def love_numbers_sampler(model,n,s,iload,internal,r_ln,g_ln,idx_r_ln):
    """
    love_numbers_sampler
    """

    r, rho, mu, eta, rheol, rpar, r0, rho0, mu0, eta0, t0, mass0, gra, mass, mlayer, G = model

    lam = np.zeros( (6,6) )
    for i in range(6):
        lam[i,i] = 1.0

    #
    # ---- Build the propagator product
    #
    nla = len(rho)

    #for i in range(nla-1,0,-1):
    #    mu_s = complex_rigidity(s,mu[i],eta[i],rheol[i],rpar[i,:])
    #    Ydir = direct_matrix (n,r[i+1],rho[i],mu_s,gra[i+1],G)
    #    Yinv = inverse_matrix(n,r[i  ],rho[i],mu_s,gra[i  ],G)
    #    lam = lam * (Ydir * Yinv)

    prlist=[]

    for i in range(1,nla):
        mu_s = complex_rigidity(s,mu[i],eta[i],rheol[i],rpar[i,:])
        Ydir = direct_matrix (n,r[i+1],rho[i],mu_s,gra[i+1],G)
        Yinv = inverse_matrix(n,r[i  ],rho[i],mu_s,gra[i  ],G)
        prlist.append( np.matmul(Ydir, Yinv) )
        lam = np.matmul(prlist[-1], lam)

    if rheol[0]==0:
        bc = fluid_core_bc (n,r[1],rho[0],gra[1],G)
    else:
        mu_s = complex_rigidity(s,mu[0],eta[0],rheol[0],rpar[0,:])
        Ydir = direct_matrix (n,r[1],rho[0],mu_s,gra[1],G)
        bc   = Ydir[:,0:3]

    bs = surface_bc (n,r[nla],gra[nla],iload,G)

    # ---- Compute the 'R' and 'Q' arrays

    prod = np.matmul(lam, bc)

    rr = np.zeros((3,3),dtype=complex)
    qq = np.zeros((3,3),dtype=complex)

    rr[0,:] = prod[2,:]
    rr[1,:] = prod[3,:]
    if n==1:
        rr[2,:] = prod[4,:]
    else:
        rr[2,:] = prod[5,:]

    #qq[0,:] = prod[0,:]
    #qq[1,:] = prod[1,:]
    #qq[2,:] = prod[4,:]

    # bb = np.linalg.solve( rr, bs )

    lu, piv = lu_factor(rr)
    bb = lu_solve( (lu, piv), bs )
    #x  = qq * bb

    #hh = x[0] * mass / r[nla]
    #ll = x[1] * mass / r[nla]
    #kk = ( - 1 - x[2] * mass / ( r[nla] * gra[nla] ) )

    # ---- Evaluate internal LNs at the layer interfaces, if requested

    if internal:

        if len(r_ln)==0:

            # Compute the solution at the layer interfaces

            hh = np.zeros(nla, dtype=complex)
            ll = np.zeros(nla, dtype=complex)
            kk = np.zeros(nla, dtype=complex)

            y = np.zeros((nla, 6), dtype=complex)

            y[0] = np.matmul(bc, bb)

            hh[0] = y[0][0] * mass / r[nla]
            ll[0] = y[0][1] * mass / r[nla]
            kk[0] = ( - 1 - y[0][4] * mass / ( r[nla] * gra[nla] ) )

            for i in range(1,nla):
                y[i] = np.matmul(prlist[i-1], y[i-1])
                hh[i] = y[i][0] * mass / r[nla]
                ll[i] = y[i][1] * mass / r[nla]
                kk[i] = ( - 1 - y[i][4] * mass / ( r[nla] * gra[nla] ) )

            return hh, ll, kk, y

        else:

            # Compute the solution at the user-supplied radii

            hh = np.zeros(len(r_ln), dtype=complex)
            ll = np.zeros(len(r_ln), dtype=complex)
            kk = np.zeros(len(r_ln), dtype=complex)

            y = np.matmul(bc, bb)     # Initialize solution at the CMB

            # If the inner layer is solid, evaluate the solution in the core

            if rheol[0]!=0:
                if 0 in idx_r_ln:
                    mu_s = complex_rigidity(s,mu[0],eta[0],rheol[0],rpar[0,:])
                    for j in range(len(r_ln)):
                        if idx_r_ln[j]==0:
                            Ydir = direct_matrix (n,r_ln[j],rho[0],mu_s,g_ln[j],G)
                            yr = np.matmul(Ydir[:,0:3] , bb)
                            hh[j] = yr[0] * mass / r[nla]
                            ll[j] = yr[1] * mass / r[nla]
                            kk[j] = ( - 1 - yr[4] * mass / ( r[nla] * gra[nla] ) )

            for i in range(1,nla):

                if i in idx_r_ln:

                    mu_s = complex_rigidity(s,mu[i],eta[i],rheol[i],rpar[i,:])
                    Yinv = inverse_matrix(n,r[i  ],rho[i],mu_s,gra[i  ],G)

                    for j in range(len(r_ln)):

                        if idx_r_ln[j]==i:

                            Ydir = direct_matrix(n,r_ln[j],rho[i],mu_s,g_ln[j],G)
                            yr   = np.matmul( np.matmul(Ydir, Yinv), y)
                            hh[j] = yr[0] * mass / r[nla]
                            ll[j] = yr[1] * mass / r[nla]
                            kk[j] = ( - 1 - yr[4] * mass / ( r[nla] * gra[nla] ) )

                y = np.matmul( prlist[i-1], y )

    else:

        # Compute the solution at the surface

        qq[0,:] = prod[0,:]
        qq[1,:] = prod[1,:]
        qq[2,:] = prod[4,:]

        x  = np.matmul(qq, bb)

        hh = x[0] * mass / r[nla]
        ll = x[1] * mass / r[nla]
        kk = ( - 1 - x[2] * mass / ( r[nla] * gra[nla] ) )

    return hh, ll, kk
