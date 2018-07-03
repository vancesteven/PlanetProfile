#-------------------------------------------------------------------------------
#Name:              rpammonia
#Purpose:           get ammonia properties for generating the b-spline equation of state
#
#Author:            Steven Fenton Vance
#                   svance@jpl.nasa.gov
#-------------------------------------------------------------------------------

import refprop as rp

def _main(P_MPa,T_K):
	rp.setup(u'def',u'ammonia')
	
