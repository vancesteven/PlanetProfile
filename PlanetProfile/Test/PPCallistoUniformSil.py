"""Metrology template: Callisto with a uniform (constant-density) interior.

NOT a scientific model of Callisto — the production Callisto forward model
(PPCallisto: porous, structured silicate) is unchanged and remains the basis
of all committed Callisto results. This template exists solely so the CMR2
discretization-offset sidecar can be MEASURED on structures whose silicate
shell is uniform, which is the regime where

    offset = CMR2_pp_native - compute_cmr2(no-core reconstruction)

isolates PlanetProfile's radial-integration discretization error rather than
the physical structured-vs-uniform silicate difference (see
PlanetProfile/Inference/build_cmr2_offset_sidecar.py's validity guard and
plans/HANDOFF-2026-07-09-test50-sbi-validation.md, opus design review
2026-07-10; approach ratified by user 2026-07-10).

Usage (measurement grid; same 11-point Tb grid as the production cache):

    python -m PlanetProfile.Inference.build_phase_c1_cache \\
        --config <measurement config with this template + measurement pkl path> \\
        --template PlanetProfile.Test.PPCallistoUniformSil --n-grid 11
"""
from PlanetProfile.Default.Callisto.PPCallisto import Planet

# Uniform interior: PP's CalcMoIConstantRho solves a constant silicate
# density against Bulk.Cmeasured (0.3549, Anderson et al. 2001) given the
# hydrosphere — every silicate cell then shares one density, which is what
# the sidecar producer's uniformity guard requires. Porosity must be off:
# POROUS_ROCK re-enables the structured-density path (see
# Main.AssignPlanetVal, which force-clears CONSTANT_INNER_DENSITY whenever
# porosity parameters are set).
Planet.Do.CONSTANT_INNER_DENSITY = True
Planet.Do.POROUS_ROCK = False
