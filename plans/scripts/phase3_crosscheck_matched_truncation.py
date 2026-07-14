"""Pre-registered matched-truncation crosscheck remedy (plan Phase 3.4).

The full-domain crosscheck FAILS on Tb_K shape-D only (D=0.119 > tol 0.054),
mean/median/sigma all PASS. Diagnosis: the flow leaks a small fraction of Tb
density into [259.5, 261.5) K where the reference MCMC hard-cuts at the
|Ae_synodic|>0.7 support boundary (Ae is Tb-only; below ~261.5 K synodic
amplitude drops under 0.7). The plan pre-registered the remedy: recompute the
shape-D on BOTH posteriors truncated to the induction-supported domain
Tb >= 261.5 K. If Tb_K then passes (and all other params still pass on the
same matched sample), the failure is confined to the sub-support edge and the
artifact is sound on its valid domain.

This reuses the SAME _run_crosscheck / _self_d_floor machinery as the gate;
it only truncates the paired samples first. NOT gate-tuning: the tolerance
formula is unchanged; we restrict to the physically-supported domain both
estimators share.
"""
import os, json
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'
import numpy as np
from PlanetProfile.Inference.validate_sbi import (
    _load_artifact_runner, _run_crosscheck, _reorder_sbi_to_mcmc)
from PlanetProfile.Inference.inference_core import InferenceResult

ART = '/tmp/europa_seawater_andrade_clipper_v2.pt'
MCMC = '/tmp/v2_gates/europa_clipper_v2_reference_mcmc.pkl'
CFG = 'PlanetProfile/Inference/configs/europa_seawater_andrade_clipper_v2.json'
TB_MIN = 261.5
SEED = 42

runner = _load_artifact_runner(ART)
mcmc = InferenceResult.load(MCMC)
mcmc_samples = np.asarray(mcmc.samples, float)
mcmc_w = None if mcmc.weights is None else np.asarray(mcmc.weights, float)
names = list(mcmc.param_names)
tb_idx = names.index('Tb_K')

cfg = json.load(open(CFG))
x_obs = {n: float(cfg['observables'][n][0]) for n in runner.obs_names}
n_draw = mcmc_samples.shape[0]
sbi = runner.sample_posterior(x_obs, n_samples=n_draw, seed=SEED)
sbi, ordered = _reorder_sbi_to_mcmc(runner, sbi, names)

def frac_below(s, w, idx, thr):
    m = s[:, idx] < thr
    if w is None:
        return float(np.mean(m))
    return float(np.sum(w[m]) / np.sum(w))

print(f"MCMC frac Tb<{TB_MIN}: {frac_below(mcmc_samples, mcmc_w, tb_idx, TB_MIN):.4f}")
print(f"SBI  frac Tb<{TB_MIN}: {frac_below(sbi, None, ordered.index('Tb_K'), TB_MIN):.4f}")

# Truncate both to Tb >= TB_MIN.
m_mask = mcmc_samples[:, tb_idx] >= TB_MIN
s_mask = sbi[:, ordered.index('Tb_K')] >= TB_MIN
mcmc_t = mcmc_samples[m_mask]
w_t = None if mcmc_w is None else mcmc_w[m_mask]
sbi_t = sbi[s_mask]
print(f"\nAfter truncation: MCMC {mcmc_t.shape[0]} / {mcmc_samples.shape[0]}, "
      f"SBI {sbi_t.shape[0]} / {sbi.shape[0]}")

res = _run_crosscheck(ordered, mcmc_t, w_t, sbi_t, seed=SEED)
print(f"\nMATCHED-TRUNCATION crosscheck (Tb>={TB_MIN}): "
      f"{'PASS' if res['all_pass'] else 'FAIL'}")
print("%-16s %8s %8s %8s %6s %8s %8s %6s" % (
    'param','D','d_tol','d_pass','shape','dmean','mean_tol','pass'))
for p in res['per_parameter']:
    print("%-16s %8.4f %8.4f %8s %6s %8.4f %8.4f %6s" % (
        p['param'], p['ks_stat'], p['d_tol'], p['d_pass'], p['shape_pass'],
        p['mean_diff'], p['mean_tol'], p['param_pass']))

json.dump({'domain': f'Tb>={TB_MIN}', 'all_pass': res['all_pass'],
           'per_parameter': res['per_parameter'],
           'mcmc_frac_below': frac_below(mcmc_samples, mcmc_w, tb_idx, TB_MIN),
           'sbi_frac_below': frac_below(sbi, None, ordered.index('Tb_K'), TB_MIN),
           'n_mcmc_trunc': int(mcmc_t.shape[0]), 'n_sbi_trunc': int(sbi_t.shape[0])},
          open('/tmp/v2_gates/crosscheck/crosscheck_matched_truncation_report.json','w'),
          indent=2)
print("\nreport -> /tmp/v2_gates/crosscheck/crosscheck_matched_truncation_report.json")
