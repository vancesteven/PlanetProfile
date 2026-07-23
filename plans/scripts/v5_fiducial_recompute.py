"""v5 fiducial observable recompute from the v5 cache, at the IDENTICAL v3/v4
fiducial interior (Tb 264.5 K, w 35.165 ppt, core as in v4). Uses the runner's
own _derive_gravity_pair and _induction_channel_value so it is byte-identical to
runtime. Compares to the current v5 config (v4-copied) values and reports deltas.

CMR2/k2/h2 are externally specified conditioning values (Galileo MoI prior /
Mazarico projection), NOT cache-derived, so they are NOT recomputed here.
Only cache-derived channels (C20, C22, 14x Bind_*) can shift with the v5 cache.
"""
import json, numpy as np
from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference import mcmc_runner as mr

CFG = "PlanetProfile/Inference/configs/europa_clipper_v5_geodesy_11D.json"
cfg = InferenceConfig.from_json(CFG)
# v5 samples D_iceIh_km; force the derived-Tb path OFF for this fixed-theta probe
# by feeding Tb_K directly (the fiducial interior is specified in Tb, as in v4).
runner = mr.MCMCRunner(cfg)   # loads v5 cache + Ae sidecar (already built)

# Fiducial interior — identical to v4 (gravity_fiducial_2026_07_19 metadata).
theta = {
    'alpha': 0.31268,
    'log10_zeta_Ih': -0.19696, 'log10_zeta_sil': -0.19696,
    'log10_eta_Ih': 12.938, 'log10_eta_sil': 20.096,
    'Tb_K': 264.5, 'log10_wOcean_ppt': float(np.log10(35.165)),
    'R_core_km': 534.67, 'rho_core_kgm3': 6254.1,
    'dC20_nh': 0.0, 'dC22_nh': 0.0,
}

# --- gravity pair (C20, C22) ---
pair = runner._derive_gravity_pair(theta)
C20, C22 = pair
print(f"[grav] C20 = {C20:.7e}   C22 = {C22:.7e}   -C20/C22 = {-C20/C22:.5f} (expect 3.324)")

# --- induction Bind_* channels via the runner's shared interpolant ---
# _induction_channel_value is a closure inside _generate_sbi_dataset; reproduce
# its two public dependencies directly (same code path): _blended_ae_dict + Be.
from PlanetProfile.Inference.mcmc_runner import _parse_bind_channel
Ae_dict = runner._blended_ae_dict(theta)
assert Ae_dict is not None, "Ae dict None at fiducial — sidecar/label problem"
be = runner._be_excitation
assert be, "Be excitation not loaded"

def bind_val(name):
    label, comp, part = _parse_bind_channel(name)
    Ae = Ae_dict.get(label); Be_comp = be.get(label)
    if Ae is None or Be_comp is None:
        return float('nan')
    Bind = complex(Ae) * Be_comp[comp]
    return Bind.real if part == 'real' else Bind.imag

obs = cfg.observables
recomputed = {'C20': C20, 'C22': C22}
for name in obs:
    if name.startswith('Bind_'):
        recomputed[name] = bind_val(name)

print("\n%-28s %14s %14s %12s" % ("channel", "current(v4-copy)", "v5-recompute", "|delta|"))
maxrel = 0.0; maxrel_ch = None
for name, val in recomputed.items():
    cur = obs[name][0] if isinstance(obs[name], (list, tuple)) else obs[name]['value']
    sig = obs[name][1] if isinstance(obs[name], (list, tuple)) else obs[name]['sigma']
    d = abs(val - cur)
    rel = d / sig if sig else float('inf')
    if rel > maxrel: maxrel = rel; maxrel_ch = name
    flag = "  <-- >0.05sig" if rel > 0.05 else ""
    print("%-28s %14.6e %14.6e %12.3e (%.3f sig)%s" % (name, cur, val, d, rel, flag))

print(f"\nMax deviation: {maxrel:.4f} sigma at {maxrel_ch}")
json.dump({'theta': theta, 'recomputed': recomputed,
           'max_dev_sigma': maxrel, 'max_dev_channel': maxrel_ch},
          open('/tmp/v5_fiducial_recompute.json', 'w'), indent=2)
print("-> /tmp/v5_fiducial_recompute.json")
