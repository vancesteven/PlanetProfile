"""
PPTest48
Andrade hybrid hydrosphere Titan with Yao et al. (2014) spherical convection.

Identical to PPTest46 except Ice Ih convection uses the Yao 2014 3D spherical
shell scaling laws instead of Deschamps & Sotin (2001). This produces thicker
stagnant lids and lower basal heat flux for Titan's geometry (f ~ 0.96).

Used by Test48_mcmc_andrade_yao2014.py for MCMC inference with a heat-flux
consistency constraint.
"""
from copy import deepcopy

from PlanetProfile.Test.PPTest46 import Planet as _PPTest46Planet

Planet = deepcopy(_PPTest46Planet)
Planet.name = 'Test48'

Planet.Do.SPHERICAL_CONVECTION = True
