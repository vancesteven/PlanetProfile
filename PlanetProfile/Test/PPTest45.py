"""
PPTest45
Hybrid hydrosphere-thickness Titan inference prototype.

This test intentionally reuses the PPTest42 Titan Maxwell ocean setup without
modifying PPTest42.  The corresponding inference cache path keeps the
self-consistent hydrosphere PT calculation from PlanetProfile, but treats total
hydrosphere thickness as an external control parameter.  Global mass and C/MR^2
are then computed as outputs and handled by the inference likelihood rather
than being used by PlanetProfile to choose the hydrosphere/silicate boundary.
"""
from copy import deepcopy

from PlanetProfile.Test.PPTest42 import Planet as _PPTest42Planet


Planet = deepcopy(_PPTest42Planet)
Planet.name = 'Test45'

# Documentation marker consumed by humans; the implementation lives in
# PlanetProfile.Inference.hybrid_structure_cache.
Planet.hybrid_inference_note = (
    'Hydrosphere PT is self-consistent in Tb_K; total hydrosphere thickness is '
    'controlled by D_hydro_km in the inference cache, so Mtot_kg and CMR2mean '
    'are model outputs rather than PlanetProfile closure targets.'
)
