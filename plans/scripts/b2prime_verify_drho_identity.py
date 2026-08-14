"""Verify the reviewer's CRITICAL-1: does the base interface enter Ks at
net weight Delta_rho = rho_ocean - rho_ice, rather than -rho_ice?

Reviewer claim:
  (Bs - As)  == rho_ice * (f_top - f_base)          [exact]
  Bsp_Asp    == rho_ocean * f_base                  [~1.8% linearization]
  => Ks bracket == rho_ice*f_top + (rho_ocean - rho_ice)*f_base
with f_i = 4pi/15 * a_i b_i c_i (a_i^2 - b_i^2).

If true, scaling the WHOLE (Bs-As) perturbs the base term by
-rho_ice*f_base*(r-1) against a net weight of only ~Delta_rho -> order unity.
"""
import numpy as np
from PlanetProfile.Gravity.Librations import compute_eccentricities, radial2abc, moi_ellps
R_T=252.22e3; M=1.08022e20; RHO_ICE=925.0
OMEGA=5.307e-5; V=4*np.pi/3
def stack(zb,D,rho_ocean):
    R_b=R_T-zb*1e3; R_c=R_b-D*1e3
    m_s=RHO_ICE*V*(R_T**3-R_b**3); m_o=rho_ocean*V*(R_b**3-R_c**3)
    rc=(M-m_s-m_o)/(V*R_c**3)
    return np.array([R_c,R_b,R_T]), np.array([rc,rho_ocean,RHO_ICE])

for (zb,D,ro) in [(25.,36.,1005.),(25.,36.,1020.),(20.,36.,1020.),(30.,50.,1020.)]:
    radial,rho=stack(zb,D,ro)
    ep,eq=compute_eccentricities(radial,rho,OMEGA)
    a,b,c=radial2abc(radial,ep,eq)
    Ai,Bi,Ci=moi_ellps(a,b,c,rho)
    oi=1
    BmA = Bi[oi+1:].sum()-Ai[oi+1:].sum()
    beta=(a-b)/a
    Bsp_Asp = 8/15*np.pi*rho[oi]*beta[oi]*radial[oi]**5
    # reviewer's f_i
    f = lambda i: 4*np.pi/15 * a[i]*b[i]*c[i]*(a[i]**2-b[i]**2)
    f_top, f_base = f(2), f(1)
    claim_BmA = RHO_ICE*(f_top - f_base)
    claim_Bsp = ro*f_base
    print(f"\n--- zb={zb} D={D} rho_ocean={ro} ---")
    print(f"  (Bs-As) actual {BmA:.8e}  claim rho_ice*(f_top-f_base) {claim_BmA:.8e}  rel {abs(claim_BmA/BmA-1):.2e}")
    print(f"  Bsp_Asp actual {Bsp_Asp:.8e}  claim rho_oc*f_base       {claim_Bsp:.8e}  rel {abs(claim_Bsp/Bsp_Asp-1):.2e}")
    brack = BmA + Bsp_Asp
    claim_brack = RHO_ICE*f_top + (ro-RHO_ICE)*f_base
    print(f"  Ks bracket actual {brack:.8e}  claim rho_ice*f_top+drho*f_base {claim_brack:.8e}  rel {abs(claim_brack/brack-1):.2e}")
    print(f"  net base weight: drho = {ro-RHO_ICE:.1f} vs the -rho_ice = {-RHO_ICE:.1f} inside (Bs-As)")
    # effective base scale actually applied by shipped whole-difference path
    H22_hyd=(a[-1]-b[-1])/6; r=857.0/H22_hyd
    # shipped: brack_new = BmA*r + Bsp_Asp = rho_ice*f_top*r - rho_ice*f_base*r + rho_oc*f_base
    # write as rho_ice*f_top*r + drho*f_base*s_eff  ->  solve s_eff
    s_eff = (RHO_ICE*f_top*r - RHO_ICE*f_base*r + ro*f_base - RHO_ICE*f_top*r)/((ro-RHO_ICE)*f_base)
    print(f"  r_top = {r:.4f}   effective base scale applied by shipped path = {s_eff:+.3f}  (docstring says 1.000)")
