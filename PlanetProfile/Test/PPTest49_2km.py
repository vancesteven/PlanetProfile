"""
PPTest49_2km
Andrade hybrid hydrosphere Titan with Yao et al. (2014) spherical convection
and a 2 km clathrate cap.

More consistent with high heat flux (40 mW/m2) per Petricca et al. (2025).
Part of the Test 49 sensitivity series.
"""
from copy import deepcopy

from PlanetProfile.Test.PPTest48 import Planet as _PPTest48Planet

Planet = deepcopy(_PPTest48Planet)
Planet.name = 'Test49_2km'

# Perturbation: 2 km clathrate cap
Planet.Bulk.clathMaxThick_m = 2e3
Planet.Steps.nClath = 30
