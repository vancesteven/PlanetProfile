"""Joint-anchor limits/W1 gate for Europa Clipper v2 (grid-walk primary + plan-set probes).

Each anchor is a FULL distinct 19-channel x_obs (all frequencies co-vary), so the
stock `validate_sbi limits --anchor-results` single-channel-sweep path does not apply.
This reuses the SAME W1 statistic, tolerance (0.25*sigma_anchor), and prior-box
containment (==1.0) as validate_sbi._run_limits_anchor_check, but conditions the flow
on each anchor config's own observables.

PRIMARY GATE = grid-walk anchors (on-manifold, single physical Tb each): every anchor
must satisfy W1(flow Tb, MCMC-anchor Tb) <= 0.25*sigma_anchor AND 100% prior-box
containment.

PROBES = plan-set independent-frequency anchors: reported with the SAME numbers but
judged by the LOOSER criterion "do flow and MCMC medians rail consistently?" (many are
off-manifold / boundary-railed by design; NOT the calibration gate).
"""
import os, json, glob
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import numpy as np
from PlanetProfile.Inference.validate_sbi import (
    _load_artifact_runner, _normalize_weights, _weighted_quantile,
    _weighted_mean_std, _wasserstein1_weighted,
    LIMITS_ANCHOR_SIGMA_FRAC, LIMITS_CONTAINMENT)
from PlanetProfile.Inference.inference_core import InferenceResult

ART = '/tmp/europa_seawater_andrade_clipper_v2.pt'
PKL_DIR = '/tmp/v2_gates/anchor_pkls'
GW_MANIFEST = 'PlanetProfile/Inference/configs/europa_clipper_v2_anchors_gridwalk/anchor_manifest.json'
PS_MANIFEST = 'PlanetProfile/Inference/configs/europa_clipper_v2_anchors/anchor_manifest.json'
MONOTONE = 'Tb_K'
N_SAMPLES = 2000
SEED = 42

runner = _load_artifact_runner(ART)
pidx = list(runner.param_names).index(MONOTONE)
bounds = np.asarray(runner.artifact_meta['param_bounds'], float)

def eval_anchor(cfg_path, pkl_path, k):
    cfg = json.load(open(cfg_path))
    x_obs = {n: float(cfg['observables'][n][0]) for n in runner.obs_names}
    mcmc = InferenceResult.load(pkl_path)
    aidx = list(mcmc.param_names).index(MONOTONE)
    a = np.asarray(mcmc.samples, float)[:, aidx]
    aw = _normalize_weights(None if mcmc.weights is None else np.asarray(mcmc.weights, float), len(a))
    a_median = float(_weighted_quantile(a, aw, 0.5))
    _, a_std = _weighted_mean_std(a[:, None], aw); a_std = float(a_std[0])
    s = runner.sample_posterior(x_obs, n_samples=N_SAMPLES, seed=SEED + k, reject_outside_prior=False)
    f_median = float(np.median(s[:, pidx]))
    inside = np.all((s >= bounds[:, 0] - 1e-6) & (s <= bounds[:, 1] + 1e-6), axis=1)
    containment = float(np.mean(inside))
    w1 = _wasserstein1_weighted(a, aw, s[:, pidx])
    tol = LIMITS_ANCHOR_SIGMA_FRAC * a_std
    return dict(anchor_median=a_median, anchor_std=a_std, flow_median=f_median,
                w1=float(w1), tol=float(tol), w1_pass=bool(w1 <= tol),
                containment=containment,
                containment_pass=bool(containment >= LIMITS_CONTAINMENT - 1e-12))

# --- PRIMARY: grid-walk ---
gw = json.load(open(GW_MANIFEST))
print("=== GRID-WALK anchors (PRIMARY W1 calibration gate) ===")
print("%-9s %-9s %-9s %-8s %-8s %-6s %-6s" % (
    'Tb_K','|Ae_syn|','anch_med','flow_med','W1','tol','pass'))
gw_rows = []; gw_all_pass = True
for k, a in enumerate(gw['anchors']):
    pkl = f"{PKL_DIR}/gw_anchor_tb{a['tb_K']:.2f}".replace('.', 'p') + '.pkl'
    r = eval_anchor(a['config'], pkl, k)
    ok = r['w1_pass'] and r['containment_pass']
    gw_all_pass = gw_all_pass and ok
    gw_rows.append({**a, **r, 'pass': ok, 'pkl': pkl})
    print("%-9.2f %-9.4f %-9.3f %-8.3f %-8.4f %-6.4f %-6s" % (
        a['tb_K'], a['amp_synodic'], r['anchor_median'], r['flow_median'],
        r['w1'], r['tol'], ok))
    if not r['containment_pass']:
        print(f"    containment {r['containment']:.3f} < 1.0")
print(f"\nGRID-WALK LIMITS GATE: {'PASS' if gw_all_pass else 'FAIL'}")

# --- PROBES: plan-set ---
ps = json.load(open(PS_MANIFEST))
ps_configs = []
for lbl, rows in ps['sweeps'].items():
    for row in rows:
        ps_configs.append(row['config'])
ps_configs.append(ps['corner']['config'])
print("\n=== PLAN-SET anchors (EXTRAPOLATION PROBES; looser criterion) ===")
print("%-26s %-9s %-8s %-8s %-8s %-6s" % ('anchor','anch_med','flow_med','W1','tol','w1<=tol'))
ps_rows = []
for k, cfg in enumerate(sorted(set(ps_configs))):
    base = os.path.basename(cfg)[:-5]
    pkl = f"{PKL_DIR}/ps_{base}.pkl"
    if not os.path.exists(pkl):
        print(f"  {base}: pkl missing, skip"); continue
    r = eval_anchor(cfg, pkl, 100 + k)
    ps_rows.append({'anchor': base, **r})
    print("%-26s %-9.3f %-8.3f %-8.4f %-8.4f %-6s" % (
        base, r['anchor_median'], r['flow_median'], r['w1'], r['tol'], r['w1_pass']))

json.dump({'gridwalk_gate_pass': gw_all_pass, 'gridwalk': gw_rows,
           'planset_probes': ps_rows,
           'gate': f'W1 <= {LIMITS_ANCHOR_SIGMA_FRAC}*sigma_anchor ({MONOTONE}) '
                   f'+ containment=={LIMITS_CONTAINMENT}',
           'note': 'grid-walk = primary calibration gate (on-manifold single-Tb); '
                   'plan-set = extrapolation probes judged by consistency, not the gate'},
          open('/tmp/v2_gates/limits_joint_report.json', 'w'), indent=2)
print("\nreport -> /tmp/v2_gates/limits_joint_report.json")
