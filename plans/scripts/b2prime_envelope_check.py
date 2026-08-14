"""Independent envelope check on the B2' discriminant (Opus, read-only).

Purpose: give me a number to VERIFY the subagent's sweep against, computed
independently. Hydrostatic treatment only (treatment 1), rho_ocean=1020 to
match H&M Fig. 4. Solves zb at libration = 0.120 deg (H&M/Thomas) and
0.092 deg (Park), over a few ocean thicknesses.
"""
import numpy as np
from scipy.optimize import brentq
from PlanetProfile.Gravity.Librations import librations

R_T = 252.22e3
M = 1.08022e20
RHO_ICE = 925.0
RHO_OCEAN = 1020.0   # H&M Fig. 4 nominal
OMEGA = 5.307e-5
ECC = 0.0047
V = 4.0*np.pi/3.0

def stack(zb_km, D_ocean_km):
    R_b = R_T - zb_km*1e3
    R_c = R_b - D_ocean_km*1e3
    if R_c <= 0 or R_b <= 0:
        return None
    m_shell = RHO_ICE*V*(R_T**3 - R_b**3)
    m_ocean = RHO_OCEAN*V*(R_b**3 - R_c**3)
    rho_core = (M - m_shell - m_ocean)/(V*R_c**3)
    if rho_core <= RHO_OCEAN:
        return None
    return np.array([R_c, R_b, R_T]), np.array([rho_core, RHO_OCEAN, RHO_ICE]), rho_core

def lib_deg(zb_km, D_ocean_km, H22_obs_m=None):
    s = stack(zb_km, D_ocean_km)
    if s is None:
        return np.nan
    radial, rho, _ = s
    m = librations(radial, rho, OMEGA, ECC, rigid=True, ocean=True,
                   ocean_idx=1, H22_obs_m=H22_obs_m)
    return float(np.degrees(m/R_T))

print(f"{'D_oc':>6} {'zb@0.120':>9} {'zb@0.092':>9} {'dzb':>7} {'rho_core@0.120':>14}")
for D in [20.0, 30.0, 36.0, 42.0, 50.0, 60.0, 70.0]:
    row = []
    for target in (0.120, 0.092):
        f = lambda z: lib_deg(z, D, None) - target
        lo, hi = 2.0, 120.0
        try:
            flo, fhi = f(lo), f(hi)
            if not (np.isfinite(flo) and np.isfinite(fhi)) or flo*fhi > 0:
                # scan for a bracket
                zs = np.arange(2.0, 120.0, 0.5)
                vals = np.array([f(z) for z in zs])
                ok = np.where(np.isfinite(vals))[0]
                sign = np.where(np.diff(np.sign(vals[ok])) != 0)[0]
                if len(sign) == 0:
                    row.append(np.nan); continue
                i = ok[sign[0]]; lo, hi = zs[i], zs[i+1]
            row.append(brentq(f, lo, hi, xtol=1e-4))
        except Exception:
            row.append(np.nan)
    z120, z092 = row
    rc = stack(z120, D)[2] if np.isfinite(z120) else np.nan
    d = z092 - z120 if np.isfinite(z120) and np.isfinite(z092) else np.nan
    print(f"{D:6.1f} {z120:9.2f} {z092:9.2f} {d:7.2f} {rc:14.1f}")

# local sensitivity at the fiducial for cross-checking the extrapolation
z0 = 25.0
h = 0.25
d1 = (lib_deg(z0+h,36.0)-lib_deg(z0-h,36.0))/(2*h)
print(f"\nd(lib)/d(zb) at zb=25,D=36: {d1:.6f} deg/km  ({d1/0.003:.3f} sigma/km at sigma=0.003)")
print(f"naive extrapolation of 0.028 deg excursion: {0.028/abs(d1):.2f} km thicker")

print("\n=== treatment comparison at BOTH librations (D_ocean sweep) ===")
print(f"{'D_oc':>6} | {'hyd@.120':>9} {'obsfig@.120':>11} | {'hyd@.092':>9} {'obsfig@.092':>11}")
def solve(D, target, H22):
    f = lambda z: lib_deg(z, D, H22) - target
    zs = np.arange(2.0, 120.0, 0.5)
    vals = np.array([f(z) for z in zs])
    ok = np.where(np.isfinite(vals))[0]
    if len(ok) < 2: return np.nan
    sg = np.where(np.diff(np.sign(vals[ok])) != 0)[0]
    if len(sg) == 0: return np.nan
    i = ok[sg[0]]
    try: return brentq(f, zs[i], zs[i+1], xtol=1e-4)
    except Exception: return np.nan

rows=[]
for D in [20.0, 30.0, 36.0, 42.0, 50.0, 60.0, 70.0]:
    a=solve(D,0.120,None); b=solve(D,0.120,857.0)
    c=solve(D,0.092,None); d=solve(D,0.092,857.0)
    rows.append((D,a,b,c,d))
    print(f"{D:6.1f} | {a:9.2f} {b:11.2f} | {c:9.2f} {d:11.2f}")
arr=np.array(rows)
for j,name in [(1,'hydrostatic @0.120'),(2,'observed-figure @0.120'),(3,'hydrostatic @0.092'),(4,'observed-figure @0.092')]:
    v=arr[:,j]; v=v[np.isfinite(v)]
    print(f"{name:26s} range {v.min():6.2f} - {v.max():6.2f} km   overlaps HM 16-22: {v.min()<=22 and v.max()>=16}")
