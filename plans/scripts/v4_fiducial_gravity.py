"""v4 fiducial (C20,C22) via the runner's own _derive_gravity_pair (same
code as runtime). Ae precompute stubbed (gravity pair needs only the
cached density profile + mass-conservation core, not induction)."""
import json, numpy as np
from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference import mcmc_runner as mr
from PlanetProfile.Inference import gravity_obs as go

mr.MCMCRunner._precompute_ae_grid = lambda self, obs: setattr(self, '_ae_grid_cache', {})

CFG = "PlanetProfile/Inference/configs/europa_seawater_andrade_clipper_v3_8D.json"
cfg = InferenceConfig.from_json(CFG)
cfg.gravity_forward_model = 'clairaut_hydrostatic'
runner = mr.MCMCRunner(cfg)

theta = {
    'alpha': 0.31268, 'log10_zeta': -0.19696,
    'log10_eta_Ih': 12.938, 'log10_eta_sil': 20.096,
    'Tb_K': 264.5, 'log10_wOcean_ppt': float(np.log10(35.165)),
    'R_core_km': 534.67, 'rho_core_kgm3': 6254.1,
    'dC20_nh': 0.0, 'dC22_nh': 0.0,
}
pair = runner._derive_gravity_pair(theta)
print("fiducial w =", 10**theta['log10_wOcean_ppt'], "ppt; pair =", pair)
C20, C22 = pair
print(f"C20_fiducial = {C20:.6e}")
print(f"C22_fiducial = {C22:.6e}")
print(f"-C20/C22 = {-C20/C22:.5f} (expect 3.324)")
q = go.rotation_parameter(2.0479e-05, 1560800.0, 4.8e22)
kf = 4.0*C22/(q*(1560800.0/go.R_REF_GC21_M)**2)
print(f"q_r = {q:.5e}; implied k_f = {kf:.5f}; RD C/MR2 = {go.radau_darwin_cmr2(kf):.5f}")
json.dump({'C20_fiducial': C20, 'C22_fiducial': C22, 'ratio': -C20/C22,
           'q_r': q, 'kf': kf, 'theta': theta},
          open('/tmp/v4_fiducial_gravity.json','w'), indent=2)
print("-> /tmp/v4_fiducial_gravity.json")
