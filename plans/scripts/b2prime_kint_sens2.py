"""K_int robustness, with pole-aware root finding and |sol| convention."""
import numpy as np
from scipy.optimize import brentq
from PlanetProfile.Gravity.Librations import (librations, compute_eccentricities,
                                             radial2abc, moi_ellps)
from PlanetProfile.Utilities.defineStructs import Constants
G=Constants.G
R_T=252.22e3; M=1.08022e20; RHO_ICE=925.0; RHO_OCEAN=1020.0
OMEGA=5.307e-5; ECC=0.0047; V=4*np.pi/3

def stack(zb,D):
    R_b=R_T-zb*1e3; R_c=R_b-D*1e3
    if R_c<=1e3: return None
    m_s=RHO_ICE*V*(R_T**3-R_b**3); m_o=RHO_OCEAN*V*(R_b**3-R_c**3)
    rc=(M-m_s-m_o)/(V*R_c**3)
    if rc<=RHO_OCEAN: return None
    return np.array([R_c,R_b,R_T]), np.array([rc,RHO_OCEAN,RHO_ICE])

def lib_local(zb,D,H22=None,kint_scale=1.0,oi=1):
    s=stack(zb,D)
    if s is None: return np.nan
    radial,rho=s
    ep,eq=compute_eccentricities(radial,rho,OMEGA)
    a,b,c=radial2abc(radial,ep,eq)
    Ai,Bi,Ci=moi_ellps(a,b,c,rho)
    As,Bs,Cs=Ai[oi+1:].sum(),Bi[oi+1:].sum(),Ci[oi+1:].sum()
    Ac,Bc,Cc=Ai[:oi].sum(),Bi[:oi].sum(),Ci[:oi].sum()
    beta=(a-b)/a
    Bip_Aip=8/15*np.pi*rho[oi]*beta[0]*radial[0]**5
    Bsp_Asp=8/15*np.pi*rho[oi]*beta[oi]*radial[oi]**5
    sc = 1.0 if H22 is None else H22/((a[-1]-b[-1])/6)
    Ks=3*OMEGA**2*((Bs-As)*sc+Bsp_Asp); Kc=3*OMEGA**2*((Bc-Ac)-Bip_Aip)
    K_int=kint_scale*(4*np.pi*G/5*(rho[-1]*beta[-1]+(rho[1]-rho[-1])*beta[1])
                      *((rho[0]-rho[1])*beta[0]*radial[0]**5))
    if Ks<=0 or Kc<=0: return np.nan
    ss=np.sqrt(Ks/Cs); sk=np.sqrt(Kc/Cc)
    A=np.array([[-OMEGA**2*Cs+ss**2*Cs+2*K_int,-2*K_int],
                [-2*K_int,-OMEGA**2*Cc+sk**2*Cc+2*K_int]])
    bb=np.array([2*ECC*ss**2*Cs,2*ECC*sk**2*Cc])
    return float(abs(np.degrees(np.linalg.solve(A,bb)[0]*radial[-1]/R_T)))

def lib_ship(zb,D,H22=None):
    s=stack(zb,D)
    if s is None: return np.nan
    radial,rho=s
    return float(np.degrees(librations(radial,rho,OMEGA,ECC,rigid=True,ocean=True,ocean_idx=1,H22_obs_m=H22)/R_T))

print("=== verification (|local| vs ship), must be ~1e-12 ===")
mx=0.0
for zb in (10.,15.,25.,35.,45.):
    for H22 in (None,857.0):
        s=lib_ship(zb,36.,H22); l=lib_local(zb,36.,H22)
        if np.isfinite(s) and s>0: mx=max(mx,abs(l/s-1))
print(f"  max relative deviation over 10 cases: {mx:.3e}")

def solve_monotone(D,target,H22,ks):
    """Find root on the PHYSICAL monotone-decreasing branch: scan from small zb
    upward, take the FIRST crossing, and reject if the curve is not locally
    decreasing there (that indicates a pole, not a physical solution)."""
    zs=np.arange(3.0,90.0,0.25)
    v=np.array([lib_local(z,D,H22,ks) for z in zs])
    ok=np.isfinite(v)
    for i in range(len(zs)-1):
        if not(ok[i] and ok[i+1]): continue
        if (v[i]-target)*(v[i+1]-target)<0:
            if v[i+1]>=v[i]:   # increasing -> pole side, skip
                continue
            try: return brentq(lambda z: lib_local(z,D,H22,ks)-target, zs[i],zs[i+1],xtol=1e-5)
            except Exception: return np.nan
    return np.nan

KS=8*np.pi/15
print(f"\n=== zb solutions, K_int scale 1.0 vs 8pi/15 ({KS:.4f}) ===")
for target in (0.120,0.092):
    print(f"\n libration = {target} deg")
    print(f"{'D_oc':>6} {'hyd':>7} {'hyd_K':>7} {'obsfig':>7} {'obsfig_K':>9}")
    for D in [20.,30.,36.,42.,50.,60.,70.]:
        r=[solve_monotone(D,target,H,k) for H in (None,857.0) for k in (1.0,KS)]
        print(f"{D:6.1f} {r[0]:7.2f} {r[1]:7.2f} {r[2]:7.2f} {r[3]:9.2f}")
