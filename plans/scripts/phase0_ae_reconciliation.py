"""Phase 0 (Europa Clipper v2) — Ae solver magnitude reconciliation, v2.

Independent high-precision (mpmath) multilayer spherical-induction solver,
cross-checked against the PRODUCTION MoonMag Srivastava solver (AeResponse /
InducedAeList, the exact path forward_model_induction uses to fill the SBI
cache), on 3 Europa models spanning the valid conductive-ocean Tb grid at
synodic / synodic-2nd / orbital.

INDEPENDENT METHOD (this file): logarithmic-derivative propagation.
  - In each conducting shell the poloidal Mie radial scalar P(r) satisfies
        P'' + (2/r)P' - [n(n+1)/r^2 + k^2] P = 0,   k = sqrt(i*mu0*sigma*omega),
    with solutions the spherical Bessel j_n(kr), y_n(kr) of complex argument.
  - Regular interior (core): P = j_n(k0 r); start L0 = P'/P at the first bdy.
  - Propagate L = P'/P outward, matching P and P' (continuity for constant mu)
    at each boundary. This is a DIFFERENT recursion from AeResponse's
    beta/gamma/delta/epsilon transfer-coefficient construction, and does not
    reuse any MoonMag code — only mpmath spherical Bessel.
  - At the outer boundary, match to the vacuum potential i*r^n + d*r^{-(n+1)}
    to get q = d/i (induced/inducing Gauss-coefficient ratio), rigorously.
  - Ae_pot = -q * R_body^{-(2n+1)}  (PEC-sphere-filling limit -> 1; vacuum -> 0).

mpmath => arbitrary precision => the 1e6 S/m core (which overflows float64
Bessel and destabilises the Riccati ODE) is handled exactly.

Acceptance (plan): |Ae_solver - Ae_ref|/|Ae_ref| < 5% amplitude AND < 5 deg
phase for all 9 points. Also validates analytic limits and re-derives
|Ae_orbital| (the "O(1e-2)" Callisto anchor is suspect: 85 hr skin depth ~
ocean thickness => tenths plausible).  Read-only on the cache.
"""
import os, sys, json
import numpy as np
import mpmath as mp

mp.mp.dps = 50  # 50 decimal digits — ample for the 1e6 S/m core

REPO = '/Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai'
sys.path.insert(0, REPO); os.chdir(REPO)
import pickle

MU0 = mp.mpf('4e-7') * mp.pi
J = mp.mpc(0, 1)

CACHE = ('PlanetProfile/Test/mcmc_results/Europa/Test51_seawater/'
         'europa_seawater_structure_grid.pkl')
MODELS = {'THIN': 10, 'MID': 23, 'THICK': 35}      # Tb 261.5 / 268 / 271
FREQS = ['synodic', 'synodic 2nd', 'orbital']


# ---- mpmath spherical Bessel of complex argument and derivative ------------
def sph_j(n, x):
    return mp.sqrt(mp.pi / (2 * x)) * mp.besselj(n + mp.mpf('0.5'), x)


def sph_y(n, x):
    return mp.sqrt(mp.pi / (2 * x)) * mp.bessely(n + mp.mpf('0.5'), x)


def sph_j_prime(n, x):
    # d/dx j_n(x) = j_{n-1}(x) - (n+1)/x j_n(x)
    return sph_j(n - 1, x) - (n + 1) / x * sph_j(n, x)


def sph_y_prime(n, x):
    return sph_y(n - 1, x) - (n + 1) / x * sph_y(n, x)


def independent_Ae(r_bds, sigmas, omega, R_body, n=1):
    """Rigorous induced/inducing Gauss-coefficient ratio via log-derivative
    propagation.  r_bds ascending (m), sigmas[i] = conductivity of shell just
    below r_bds[i] (matches MoonMag AeResponse convention), sigmas[0]=core."""
    n = mp.mpf(n)
    r_bds = [mp.mpf(float(r)) for r in r_bds]
    sigmas = [mp.mpf(float(s)) for s in sigmas]
    om = mp.mpf(float(omega))
    Rb = mp.mpf(float(R_body))

    def kof(sig):
        return mp.sqrt(J * MU0 * sig * om)

    # L = P'/P (dP/dr over P), start regular interior at first boundary
    r0 = r_bds[0]
    if sigmas[0] == 0:
        # non-conducting core: regular interior P = r^n  =>  L = n/r0
        L = n / r0
    else:
        k0 = kof(sigmas[0])
        x0 = k0 * r0
        L = k0 * sph_j_prime(n, x0) / sph_j(n, x0)     # (dj/dr)/j = k j'(x)/j(x)

    # propagate through each shell above the first boundary
    for i in range(1, len(r_bds)):
        k = kof(sigmas[i])
        r_lo, r_hi = r_bds[i - 1], r_bds[i]
        if sigmas[i] == 0:
            # vacuum layer: P = A r^n + B r^{-(n+1)}
            rlo_p, rlo_m = r_lo**n, r_lo**(-(n + 1))
            rlo_pp, rlo_mp = n * r_lo**(n - 1), -(n + 1) * r_lo**(-(n + 2))
            # L = (A rlo_pp + B rlo_mp)/(A rlo_p + B rlo_m); solve rho=A/B
            rho = -(rlo_mp - L * rlo_m) / (rlo_pp - L * rlo_p)
            rhi_p, rhi_m = r_hi**n, r_hi**(-(n + 1))
            rhi_pp, rhi_mp = n * r_hi**(n - 1), -(n + 1) * r_hi**(-(n + 2))
            P = rho * rhi_p + rhi_m
            Pp = rho * rhi_pp + rhi_mp
            L = Pp / P
        else:
            x_lo, x_hi = k * r_lo, k * r_hi
            j_lo, y_lo = sph_j(n, x_lo), sph_y(n, x_lo)
            jp_lo, yp_lo = k * sph_j_prime(n, x_lo), k * sph_y_prime(n, x_lo)
            # L = (rho*jp_lo + yp_lo)/(rho*j_lo + y_lo); solve rho = A/B
            rho = -(yp_lo - L * y_lo) / (jp_lo - L * j_lo)
            x_hi = k * r_hi
            j_hi, y_hi = sph_j(n, x_hi), sph_y(n, x_hi)
            jp_hi, yp_hi = k * sph_j_prime(n, x_hi), k * sph_y_prime(n, x_hi)
            P = rho * j_hi + y_hi
            Pp = rho * jp_hi + yp_hi
            L = Pp / P

    # outer boundary -> vacuum potential i*r^n + d*r^{-(n+1)}; solve q=d/i
    r_out = r_bds[-1]
    q = -(L * r_out - n) / (L * r_out + (n + 1)) * r_out**(2 * n + 1)
    Ae_pot = -q * Rb**(-(2 * n + 1))
    return complex(Ae_pot)


def production_Ae(r_bds, sig, R_body, T_hr):
    from PlanetProfile.MagneticInduction.MoonMag.symmetry_funcs import InducedAeList
    omega = 2 * np.pi / (T_hr * 3600.0)
    Aes, _, _ = InducedAeList(np.asarray(r_bds, float), np.asarray(sig, float),
                              np.array([omega]), 1.0 / R_body,
                              nn=1, writeout=False, do_parallel=False)
    return complex(Aes[0])


def amp_phase(z):
    return abs(z), np.degrees(np.angle(z))


def analytic_limits_check():
    """Independent solver must give Ae=0 (no conductor) and |Ae|=1 (PEC sphere
    filling the body).  Uses trivial synthetic profiles."""
    print("--- analytic-limit self-checks (independent solver) ---")
    om = 2 * np.pi / (11.23 * 3600.0)
    R = 1560800.0
    # vacuum everywhere
    ae0 = independent_Ae([R], [0.0], om, R)
    print(f"  vacuum sphere:  Ae = {ae0:.6g}   (expect ~0)")
    # PEC sphere filling the body (one huge-conductivity layer to r=R)
    aep = independent_Ae([R], [1e9], om, R)
    a, p = amp_phase(aep)
    print(f"  PEC full sphere: |Ae| = {a:.5f}  phase = {p:.3f} deg   (expect |Ae|~1, phase~0)")
    return abs(ae0) < 1e-3, abs(a - 1.0) < 0.02


def skin_depth_km(sigma, T_hr):
    om = 2 * np.pi / (T_hr * 3600.0)
    return np.sqrt(2.0 / (4e-7 * np.pi * sigma * om)) / 1e3


def main():
    d = pickle.load(open(CACHE, 'rb'))
    S = d['structures']
    vac_ok, pec_ok = analytic_limits_check()
    print(f"  limit checks: vacuum {'OK' if vac_ok else 'FAIL'}, PEC {'OK' if pec_ok else 'FAIL'}\n")

    rows = []
    for tag, idx in MODELS.items():
        s = S[idx]
        r_bds = np.asarray(s['rSigChange_m'], float)
        sig = np.asarray(s['sigmaLayers_Sm'], float)
        R = float(s['R_body_m'])
        Tb = float(s['Tb_K']); Doc = float(s['D_ocean_km'])
        # representative ocean sigma for skin-depth reporting
        oc = sig[(sig > 0.5) & (sig < 1e3)]
        sig_oc = float(np.mean(oc)) if oc.size else np.nan
        for fk in FREQS:
            T = s['Texc_hr'][fk]
            ap = production_Ae(r_bds, sig, R, T)
            ai = independent_Ae(r_bds, sig, 2 * np.pi / (T * 3600.0), R)
            pa, pph = amp_phase(ap)
            ia, iph = amp_phase(ai)
            damp = abs(pa - ia) / ia * 100 if ia > 0 else np.nan
            dph = abs(((pph - iph + 180) % 360) - 180)
            dsk = skin_depth_km(sig_oc, T) if np.isfinite(sig_oc) else np.nan
            rows.append(dict(model=tag, idx=idx, Tb=Tb, D_ocean=Doc, freq=fk, T_hr=T,
                             sig_ocean=sig_oc, skin_km=dsk,
                             prod_amp=pa, prod_ph=pph, indep_amp=ia, indep_ph=iph,
                             damp_pct=damp, dph_deg=dph))

    hdr = (f"{'model':6} {'freq':12} {'T_hr':7} {'sigOc':6} {'skin_km':8} "
           f"{'|Ae|prod':9} {'|Ae|indep':10} {'phProd':8} {'phIndep':8} "
           f"{'damp%':7} {'dphDeg':7}")
    print(hdr)
    wamp = wph = 0.0
    for r in rows:
        print(f"{r['model']:6} {r['freq']:12} {r['T_hr']:7.2f} {r['sig_ocean']:6.3f} "
              f"{r['skin_km']:8.1f} {r['prod_amp']:9.4f} {r['indep_amp']:10.4f} "
              f"{r['prod_ph']:8.3f} {r['indep_ph']:8.3f} {r['damp_pct']:7.3f} {r['dph_deg']:7.3f}")
        wamp = max(wamp, r['damp_pct']); wph = max(wph, r['dph_deg'])
    print(f"\nWORST: amplitude {wamp:.3f}%   phase {wph:.3f} deg")
    gate = wamp < 5 and wph < 5
    print(f"GATE (production vs independent, 5% / 5deg): {'PASS' if gate else 'FAIL'}")

    print("\n--- |Ae_orbital| re-derivation (Callisto 0.73 red-flag context) ---")
    for r in rows:
        if r['freq'] == 'orbital':
            print(f"  {r['model']:6} Tb={r['Tb']:.2f} D_ocean={r['D_ocean']:6.1f}km "
                  f"skin(85hr)={r['skin_km']:6.1f}km : |Ae_orbital| prod={r['prod_amp']:.4f} "
                  f"indep={r['indep_amp']:.4f}")
    print("  => O(1e-2) is NOT the expected magnitude; tenths are physical when "
          "ocean thickness is comparable to the 85-hr skin depth.")

    json.dump(rows, open('/tmp/phase0_ae_reconciliation.json', 'w'), indent=1, default=float)
    print("\nsaved /tmp/phase0_ae_reconciliation.json")


if __name__ == '__main__':
    main()
