"""
PPTest49
Andrade hybrid hydrosphere Titan with Yao et al. (2014) spherical convection
and a 4 km clathrate cap.

Identical to PPTest48 except the clathrate cap is 4 km instead of 5 km.
Used to check the sensitivity of the Path B result to clathrate thickness.
"""
from copy import deepcopy

from PlanetProfile.Test.PPTest48 import Planet as _PPTest48Planet

Planet = deepcopy(_PPTest48Planet)
Planet.name = 'Test49'

# Perturbation: 4 km clathrate cap (Path B / Test 48 had 5 km via PPTest42 inheritance)
Planet.Bulk.clathMaxThick_m = 4e3
Planet.Steps.nClath = 30
