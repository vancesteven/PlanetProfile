"""Lookup (Nearest2DInterpolator on a P,T grid) vs direct SeaFreeze whichphase:
disagreement rate on random points and the induced error in the Ih/liquid and
liquid/HP boundaries."""
import sys, time, numpy as np
sys.path.insert(0, '.')
import seafreeze.seafreeze as sf
from PlanetProfile.Utilities.DataManip import Nearest2DInterpolator
from PlanetProfile.Thermodynamics.HydroEOS import sfPTgrid, sfPTpairs, SFphase
rng = np.random.default_rng(3)
for dP, dT in ((1.0, 0.5), (0.5, 0.25), (0.25, 0.1)):
    Pg = np.arange(0.1, 1500, dP); Tg = np.arange(240, 300, dT)
    t0 = time.perf_counter(); tab = sf.whichphase(sfPTgrid(Pg, Tg)).astype(np.int8); tb = time.perf_counter() - t0
    lut = Nearest2DInterpolator(Pg, Tg, tab)
    P = rng.uniform(0.1, 1499, 20000); T = rng.uniform(240, 299.9, 20000)
    a = lut(P, T); b = sf.whichphase(sfPTpairs(P, T)).astype(np.int8)
    mism = a != b
    print(f'grid dP={dP} dT={dT}: build {tb:.1f} s ({tab.size} pts); mismatch {mism.mean():.2%} of 20k random points; '
          f'mismatched |dP| to boundary <= {dP/2} MPa by construction')
f = SFphase(0.0, 'PureH2O')
for T in (250.0, 260.0, 268.0):
    # Ih -> liquid along an isotherm, located by bisection with each finder
    def edge(fn):
        lo, hi = 0.1, 400.0
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if int(fn(mid)) == 1: lo = mid
            else: hi = mid
        return 0.5 * (lo + hi)
    Pg = np.arange(0.1, 400, 1.0); Tg = np.arange(240, 300, 0.5)
    lut = Nearest2DInterpolator(Pg, Tg, sf.whichphase(sfPTgrid(Pg, Tg)).astype(np.int8))
    e_direct = edge(lambda p: f(p, T)); e_lut = edge(lambda p: lut(np.array([p]), np.array([T]))[0])
    print(f'T={T} K: Ih upper boundary direct {e_direct:.3f} MPa, lookup(dP=1) {e_lut:.3f} MPa, diff {e_lut-e_direct:+.3f} MPa')
