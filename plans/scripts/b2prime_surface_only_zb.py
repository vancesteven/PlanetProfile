"""Treatments 2 and 3 ARE distinct after all. Re-solve zb for the true
surface-only treatment (scale only the outer-ellipsoid piece of (Bs-As))."""
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
    return np.array([R_c,R_b,R_T]), np.array([rc,RHO_OCEAN,RHO_ICE]), rc

def lib(zb,D,mode='hyd',H22=857.0,oi=1):
    s=stack(zb,D)
    if s is None: return np.nan
    radial,rho,rc=s
    ep,eq=compute_eccentricities(radial,rho,OMEGA)
    a,b,c=radial2abc(radial,ep,eq)
    Ai,Bi,Ci=moi_ellps(a,b,c,rho)
    As,Bs,Cs=Ai[oi+1:].sum(),Bi[oi+1:].sum(),Ci[oi+1:].sum()
    Ac,Bc,Cc=Ai[:oi].sum(),Bi[:oi].sum(),Ci[:oi].sum()
    beta=(a-b)/a
    f=1/5*4/3*np.pi; abc=a*b*c
    bcsq=b**2+c**2; acsq=a**2+c**2
    i=2
    BmA_surf = f*rho[i]*abc[i]*acsq[i]  - f*rho[i]*abc[i]*bcsq[i]
    BmA_base = -(f*rho[i]*abc[i-1]*acsq[i-1] - f*rho[i]*abc[i-1]*bcsq[i-1])
    H22_hyd=(a[-1]-b[-1])/6; sc=H22/H22_hyd
    if mode=='hyd':      BmA = BmA_surf + BmA_base
    elif mode=='whole':  BmA = (BmA_surf + BmA_base)*sc
    elif mode=='surf':   BmA = BmA_surf*sc + BmA_base
    Bip_Aip=8/15*np.pi*rho[oi]*beta[0]*radial[0]**5
    Bsp_Asp=8/15*np.pi*rho[oi]*beta[oi]*radial[oi]**5
    Ks=3*OMEGA**2*(BmA+Bsp_Asp); Kc=3*OMEGA**2*((Bc-Ac)-Bip_Aip)
    if Ks<=0 or Kc<=0: return np.nan
    K_int=4*np.pi*G/5*(rho[-1]*beta[-1]+(rho[1]-rho[-1])*beta[1])*((rho[0]-rho[1])*beta[0]*radial[0]**5)
    ss=np.sqrt(Ks/Cs); sk=np.sqrt(Kc/Cc)
    A=np.array([[-OMEGA**2*Cs+ss**2*Cs+2*K_int,-2*K_int],[-2*K_int,-OMEGA**2*Cc+sk**2*Cc+2*K_int]])
    bb=np.array([2*ECC*ss**2*Cs,2*ECC*sk**2*Cc])
    return float(abs(np.degrees(np.linalg.solve(A,bb)[0]*radial[-1]/R_T)))

# verify hyd + whole against shipped
print("verification vs shipped librations():")
for zb in (15.,25.,35.):
    s=stack(zb,36.); radial,rho,_=s
    sh_h=float(np.degrees(librations(radial,rho,OMEGA,ECC,rigid=True,ocean=True,ocean_idx=1)/R_T))
    sh_w=float(np.degrees(librations(radial,rho,OMEGA,ECC,rigid=True,ocean=True,ocean_idx=1,H22_obs_m=857.0)/R_T))
    print(f"  zb={zb:5.1f} hyd {lib(zb,36.,'hyd'):.10f} vs {sh_h:.10f} | whole {lib(zb,36.,'whole'):.10f} vs {sh_w:.10f}")

def solve(D,target,mode):
    zs=np.arange(3.,90.,0.25); v=np.array([lib(z,D,mode) for z in zs])
    ok=np.isfinite(v)
    for i in range(len(zs)-1):
        if not(ok[i] and ok[i+1]): continue
        if (v[i]-target)*(v[i+1]-target)<0 and v[i+1]<v[i]:
            try: return brentq(lambda z: lib(z,D,mode)-target,zs[i],zs[i+1],xtol=1e-5)
            except Exception: return np.nan
    return np.nan

for target in (0.120,0.092):
    print(f"\n=== libration {target} deg ===")
    print(f"{'D_oc':>6} {'hyd':>7} {'whole':>7} {'surf_only':>10} {'rho_core':>9}")
    rows=[]
    for D in [20.,30.,36.,42.,50.,60.,70.]:
        r=[solve(D,target,m) for m in ('hyd','whole','surf')]
        rc=stack(r[0],D)[2] if np.isfinite(r[0]) else np.nan
        rows.append(r); print(f"{D:6.1f} {r[0]:7.2f} {r[1]:7.2f} {r[2]:10.2f} {rc:9.0f}")
    arr=np.array(rows)
    for j,n in [(0,'hydrostatic'),(1,'whole-difference'),(2,'surface-only')]:
        v=arr[:,j][np.isfinite(arr[:,j])]
        print(f"  {n:18s} {v.min():6.2f} - {v.max():6.2f} km")
    # physically admissible subrange D=30-42
    print("  -- restricted to D_ocean 30-42 km (rho_core in/near HM 2340-2410):")
    for j,n in [(0,'hydrostatic'),(1,'whole-difference'),(2,'surface-only')]:
        v=arr[1:4,j]; v=v[np.isfinite(v)]
        print(f"  {n:18s} {v.min():6.2f} - {v.max():6.2f} km")
