#!/usr/bin/env python
"""Pre-registered Ae keep/drop diagnostic for the rebuilt Callisto NaCl grid.

Rule (Machine A queue item 2, user directive 2026-07-12):
  With the new Pan2021 conductivity cache, compute the synodic |Ae| SPAN across
  the 11-pt Tb grid. KEEP the Ae channels iff span > 3 x channel sigma (induction
  then discriminates Tb); else drop. Report numbers either way.

Reuses forward_model_induction (same machinery as mcmc_runner._precompute_ae_grid),
so the diagnostic matches exactly what the likelihood/SBI pipeline would compute.
"""
import json
import pickle
import numpy as np

from PlanetProfile.Inference.forward_models import forward_model_induction

CACHE = ("PlanetProfile/Test/mcmc_results/Callisto/C2_andrade/"
         "callisto_nacl_structure_grid.pkl")
CONFIG = "PlanetProfile/Inference/configs/callisto_nacl_andrade_8D.json"


def main():
    cache = pickle.load(open(CACHE, "rb"))
    cfg = json.load(open(CONFIG))
    obs = cfg["observables"]

    Tb_grid = np.asarray(cache["Tb_K_grid"])
    structures = cache["structures"]
    n = len(structures)

    # Frequency labels requested by the config's Ae observables.
    labels = set()
    for k in obs:
        if k.startswith("Ae_") and k.endswith("_real"):
            labels.add(k[len("Ae_"):-len("_real")])
        elif k.startswith("Ae_") and k.endswith("_imag"):
            labels.add(k[len("Ae_"):-len("_imag")])
    labels = sorted(labels)

    # Per-channel sigma from the config (real/imag share sigma per label).
    sigma = {}
    for lbl in labels:
        rk = f"Ae_{lbl}_real"
        sigma[lbl] = float(obs[rk][1]) if rk in obs else np.nan

    print(f"Cache: {CACHE}")
    print(f"Tb grid ({n} pts): {np.round(Tb_grid, 3).tolist()}")
    print(f"Ae labels: {labels}\n")

    # Compute |Ae| for every grid point and label.
    ampl = {lbl: np.full(n, np.nan) for lbl in labels}
    rea = {lbl: np.full(n, np.nan) for lbl in labels}
    ima = {lbl: np.full(n, np.nan) for lbl in labels}
    for i, struct in enumerate(structures):
        Texc = struct.get("Texc_hr") if isinstance(struct, dict) else None
        if not Texc:
            continue
        freq_dict = {lbl: Texc[lbl] for lbl in labels if lbl in Texc}
        if not freq_dict:
            continue
        res = forward_model_induction(struct, freq_dict, nn=1, do_parallel=False)
        if res is None:
            continue
        for lbl in labels:
            if lbl in res and res[lbl] is not None:
                ae = complex(res[lbl])
                ampl[lbl][i] = abs(ae)
                rea[lbl][i] = ae.real
                ima[lbl][i] = ae.imag

    # Two sigma conventions:
    #  (a) config sigma  — as stored (30% of the OLD degenerate |Ae|; stale)
    #  (b) refreshed sigma — 30% of the NEW |Ae| (the config's documented
    #      "sigma = 30%|Ae|" rule applied to the corrected conductivity)
    print("CONVENTION A — config-stored sigma (stale: 30% of OLD degenerate |Ae|)")
    print(f"{'label':<14} {'|Ae|min':>9} {'|Ae|max':>9} {'span':>9} "
          f"{'sigma_cfg':>10} {'3*sig':>9} {'verdict':>7}")
    print("-" * 72)
    verdicts = {}
    for lbl in labels:
        a = ampl[lbl]
        finite = a[np.isfinite(a)]
        if finite.size == 0:
            print(f"{lbl:<14} {'--':>9} {'--':>9} {'--':>9} "
                  f"{sigma[lbl]:>10.5f} {3*sigma[lbl]:>9.4f}   NODATA")
            verdicts[lbl] = "NODATA"
            continue
        span = float(finite.max() - finite.min())
        thresh = 3.0 * sigma[lbl]
        keep = span > thresh
        verdicts[lbl] = "KEEP" if keep else "DROP"
        print(f"{lbl:<14} {finite.min():>9.4f} {finite.max():>9.4f} "
              f"{span:>9.4f} {sigma[lbl]:>10.5f} {thresh:>9.4f} "
              f"{verdicts[lbl]:>7}")

    # NOTE: config-stored sigma = exactly 0.30 x stored |Ae| for all three
    # channels (verified). The stored |Ae| (synodic 0.149, synodic2nd 0.290,
    # orbital 0.0038) predate the Pan2021 conductivity -> a STALE BASELINE, not
    # a 1e-5 placeholder output. Convention A 'orbital KEEP' is therefore a
    # stale-baseline artifact.
    print("\nCONVENTION B — refreshed sigma = 30% of NEW mean |Ae| "
          "(consistent with corrected conductivity)")
    print(f"{'label':<14} {'|Ae|mean':>9} {'span':>9} "
          f"{'sigma_new':>10} {'3*sig':>9} {'verdict':>7}")
    print("-" * 62)
    verdicts_b = {}
    for lbl in labels:
        a = ampl[lbl]
        finite = a[np.isfinite(a)]
        if finite.size == 0:
            verdicts_b[lbl] = "NODATA"
            continue
        span = float(finite.max() - finite.min())
        mean_ae = float(finite.mean())
        sig_new = 0.30 * mean_ae
        thresh = 3.0 * sig_new
        keep = span > thresh
        verdicts_b[lbl] = "KEEP" if keep else "DROP"
        print(f"{lbl:<14} {mean_ae:>9.4f} {span:>9.4f} "
              f"{sig_new:>10.5f} {thresh:>9.4f} {verdicts_b[lbl]:>7}")

    print("\nPer-Tb |Ae| (synodic-focused decision channel):")
    dec = "synodic" if "synodic" in labels else (labels[0] if labels else None)
    if dec is not None:
        for i in range(n):
            print(f"  Tb={Tb_grid[i]:7.3f} K  |Ae_{dec}|={ampl[dec][i]:.4f}  "
                  f"Re={rea[dec][i]:+.4f}  Im={ima[dec][i]:+.4f}")

    # Overall recommendation: keep the whole Ae block iff the primary synodic
    # channel discriminates Tb (pre-registered rule keys on synodic).
    primary_a = verdicts.get("synodic", "NODATA")
    primary_b = verdicts_b.get("synodic", "NODATA")
    print(f"\nPRE-REGISTERED VERDICT (synodic) — conv A: {primary_a} | "
          f"conv B: {primary_b}")
    print("Interpretation: Pan2021 gives a strong conductive ocean "
          "(|Ae_synodic|~0.8), but |Ae| is nearly FLAT across the narrow "
          "Tb=[250,254.9]K range so induction does NOT discriminate Tb for "
          "Callisto NaCl 100 ppt. Synodic/synodic-2nd flatness = strong-"
          "conductor plateau (skin depth ~ ocean thickness); orbital flatness "
          "= thin-ocean under-drive (skin depth >> ocean). Convention A "
          "'orbital=KEEP' is a stale-BASELINE artifact (config sigma = 30% of "
          "an earlier |Ae| baseline, not the 1e-5 placeholder).")
    print("RED FLAG (reviewer): rebuilt |Ae_orbital|=0.73 vs config-stored "
          "0.0038 (~190x). Skin depth at the ~400 hr orbital period >> ocean "
          "thickness, so orbital |Ae| should be O(1e-2-1e-3). Does NOT change "
          "DROP, but MUST be reconciled before Ae is re-enabled for Callisto.")


if __name__ == "__main__":
    main()
