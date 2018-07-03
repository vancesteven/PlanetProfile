#-------------------------------------------------------------------------------
#Name:              tp_rho_liquid
#Purpose:           get ammonia properties for generating the b-spline equation of state
#
#Author:            Steven Fenton Vance
#                   svance@jpl.nasa.gov
#-------------------------------------------------------------------------------

import refprop as rp

def rho_kgm3 = _main(P_MPa,T_K,species='ammonia'):
	rp.setup(u'def',species)
	rho_molar = rp.tprho(T_K, P_MPa * 1000, [1], 1)[u'D']
	mw = rp.wmol([1])
	rho_kgm3 = rho_molar/mw
	
#	inputs:
#    t--temperature [K]
#    p--pressure [kPa]
#    x--composition [array of mol frac]
#    kph--phase flag:
#        1 = liquid
#        2 = vapor
#        0 = stable phase--NOT ALLOWED (use TPFLSH)
#            (unless an initial guess is supplied for rho)
#        -1 = force the search in the liquid phase
#        -2 = force the search in the vapor phase
#    kguess--input flag:
#        1 = first guess for D provided
#        0 = no first guess provided
#    D--first guess for molar density [mol/L], only if kguess = 1
#outputs:
#    D--molar density [mol/L]
#	
