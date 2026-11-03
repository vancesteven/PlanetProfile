"""SeaFreeze per-property (tdvSpec) cost vs all-property calls vs PP lookup splines. 2026-11-02."""
import time, numpy as np, seafreeze.seafreeze as sf
from scipy.interpolate import RectBivariateSpline
rng = np.random.default_rng(1)
def pts(n, P=(0.1, 200), T=(250, 280)):
    return np.array([(p, t) for p, t in zip(rng.uniform(*P, n), rng.uniform(*T, n))],
                    dtype='f,f').astype(object)   # PP's sfPTpairs format
def t(f, rep=3):
    best = 1e9
    for _ in range(rep):
        t0 = time.perf_counter(); out = f(); best = min(best, time.perf_counter() - t0)
    return best, out
for phase in ('Ih', 'water1'):
    T = (240, 272) if phase == 'Ih' else (273, 300)
    for n in (100, 1000, 10000):
        PT = pts(n, T=T)
        ta, full = t(lambda: sf.seafreeze(PT, phase))
        tr, sub = t(lambda: sf.getProp(PT, phase, sf.defpath, 'rho'))
        trc, sub2 = t(lambda: sf.getProp(PT, phase, sf.defpath, 'rho', 'Cp', 'alpha'))
        ok = np.allclose(sub.rho, full.rho, rtol=0, atol=1e-9)
        print(f'{phase:6s} n={n:5d}: all props {ta*1e3:7.1f} ms | rho only {tr*1e3:6.1f} ms | rho,Cp,alpha {trc*1e3:6.1f} ms | rho identical {ok}')
# whichphase cost (the phase-lookup replacement)
for n in (100, 1000, 10000):
    PT = pts(n, P=(0.1, 400), T=(240, 300))
    tw, ph = t(lambda: sf.whichphase(PT))
    print(f'whichphase n={n:5d}: {tw*1e3:7.1f} ms ({tw/n*1e6:.1f} us/pt)')
# spline lookup vs direct, accuracy and cost (PP's pattern: grid -> RectBivariateSpline)
Pg = np.linspace(0.1, 200, 201); Tg = np.linspace(273.2, 300, 136)
G = np.array([Pg, Tg], dtype=object)
t0 = time.perf_counter(); out = sf.seafreeze(G, 'water1'); tbuild = time.perf_counter() - t0
spl = RectBivariateSpline(Pg, Tg, out.rho)
PT = pts(10000, T=(273.2, 300))
P = np.array([x[0] for x in PT], float); Tt = np.array([x[1] for x in PT], float)
ts, rs = t(lambda: spl(P, Tt, grid=False))
td, rd = t(lambda: sf.getProp(PT, 'water1', sf.defpath, 'rho').rho)
print(f'water1 rho: table build {tbuild*1e3:.0f} ms once; spline eval 10k {ts*1e3:.2f} ms vs direct 10k {td*1e3:.1f} ms; '
      f'max |spline - direct| {np.max(np.abs(rs-rd)):.2e} kg/m3 ({np.max(np.abs(rs-rd)/rd):.1e} rel)')
