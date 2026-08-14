"""Does the shipped H22_obs_m path leave the BASE interface untouched?

The reviewer was starting this check when it died. The concern: H22_obs_m
rescales (Bs - As) = the shell's moment DIFFERENCE, which is built from BOTH
the surface AND the shell-base figures (moi_ellps differences consecutive
ellipsoids). So scaling the difference is NOT the same as scaling only the
surface figure -- it implicitly rescales the base contribution too.

Quantify: decompose (Bs-As) into its surface and base pieces and see what
the multiplicative scale does to each.
"""
import numpy as np
from PlanetProfile.Gravity.Librations import (librations, compute_eccentricities,
                                             radial2abc, moi_ellps)
R_T=252.22e3; M=1.08022e20; RHO_ICE=925.0; RHO_OCEAN=1020.0
OMEGA=5.307e-5; ECC=0.0047; V=4*np.pi/3
def stack(zb,D):
    R_b=R_T-zb*1e3; R_c=R_b-D*1e3
    m_s=RHO_ICE*V*(R_T**3-R_b**3); m_o=RHO_OCEAN*V*(R_b**3-R_c**3)
    rc=(M-m_s-m_o)/(V*R_c**3)
    return np.array([R_c,R_b,R_T]), np.array([rc,RHO_OCEAN,RHO_ICE])

zb,D = 25.0, 36.0
radial,rho = stack(zb,D)
ep,eq = compute_eccentricities(radial,rho,OMEGA)
a,b,c = radial2abc(radial,ep,eq)
Ai,Bi,Ci = moi_ellps(a,b,c,rho)

# moi_ellps for the shell layer (index 2) is rho[2]*(abc_2*acsq_2 - abc_1*acsq_1)
# i.e. an OUTER (surface) piece MINUS an INNER (shell base) piece.
f = 1/5*4/3*np.pi
abc = a*b*c
absq=a**2+b**2; bcsq=b**2+c**2; acsq=a**2+c**2
i=2  # shell layer
B_out = f*rho[i]*abc[i]*acsq[i];  B_in = f*rho[i]*abc[i-1]*acsq[i-1]
A_out = f*rho[i]*abc[i]*bcsq[i];  A_in = f*rho[i]*abc[i-1]*bcsq[i-1]
BmA_surface_piece = B_out - A_out
BmA_base_piece    = -(B_in - A_in)
BmA_total = BmA_surface_piece + BmA_base_piece
print(f"(Bs-As) decomposition at zb={zb}, D_ocean={D}:")
print(f"  surface (outer ellipsoid) piece : {BmA_surface_piece:.6e}")
print(f"  base   (inner ellipsoid) piece : {BmA_base_piece:.6e}")
print(f"  total                          : {BmA_total:.6e}   [check vs Bi-Ai: {Bi[2]-Ai[2]:.6e}]")
print(f"  base piece as fraction of total: {BmA_base_piece/BmA_total*100:.2f}%")

H22_hyd = (a[-1]-b[-1])/6
scale = 857.0/H22_hyd
print(f"\nH22_hyd = {H22_hyd:.3f} m ; H22_obs = 857.0 m ; scale = {scale:.6f}")
print("\nWhat the shipped code does:  (Bs-As)*scale")
print(f"   -> surface piece scaled by {scale:.6f}  AND base piece scaled by {scale:.6f}")
print("\nWhat 'surface figure only is observed' would do: scale ONLY the surface piece")
alt = BmA_surface_piece*scale + BmA_base_piece
print(f"   shipped total : {BmA_total*scale:.6e}")
print(f"   surface-only  : {alt:.6e}")
print(f"   ratio shipped/surface-only : {BmA_total*scale/alt:.6f}")
print(f"   i.e. they differ by {(BmA_total*scale/alt-1)*100:.3f}% in (Bs-As)")

# effective scale that a surface-only treatment corresponds to
eff = alt/BmA_total
print(f"\n   surface-only is equivalent to an effective whole-difference scale of {eff:.6f}")
print(f"   (vs shipped {scale:.6f})")
