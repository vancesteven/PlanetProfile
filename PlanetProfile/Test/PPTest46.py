"""
PPTest46
Andrade hybrid hydrosphere Titan inference prototype.

Thermal/pressure structure identical to PPTest42 (Maxwell ocean Titan); Andrade
rheology parameters added for SetupGravity and for the MCMC harness
(Test46_mcmc_andrade_hybrid_hydro.py).

The Maxwell harness (Test45_mcmc_maxwell_hybrid_hydro.py) uses PPTest45
(identical thermal structure) for its grid template.  The Andrade harness uses
this file (PPTest46) as its template and shares the same grid cache, because
andradExponent / andrade_zeta do NOT affect PP's thermal propagation —
they only enter SetupGravity (tidal Love number calculation).
"""
from copy import deepcopy

from PlanetProfile.Test.PPTest42 import Planet as _PPTest42Planet

Planet = deepcopy(_PPTest42Planet)
Planet.name = 'Test46'

# Andrade rheology parameters (used by SetupGravity; MCMC samples alpha and zeta)
Planet.Gravity.andradExponent = 0.2
Planet.Gravity.andrade_zeta = {
    'Ih':    0.005,
    'III':   0.05,
    'V':     0.05,
    'VI':    0.05,
    'Sil':   0.005,
    'Fe':    1.0,
    'Clath': 1.0,
    'default': 1.0,
}

Planet.hybrid_inference_note = (
    'Andrade hybrid hydrosphere: Tb_K self-consistent in PP, '
    'D_hydro_km and rho_sil controlled at inference time. '
    'alpha and log10_zeta sampled by MCMC (Test46_mcmc_andrade_hybrid_hydro.py).'
)
