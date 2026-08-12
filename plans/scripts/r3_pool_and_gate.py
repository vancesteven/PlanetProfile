"""§0.20 R3 — pool the n_eff=8000 NH3 reference seeds and evaluate the C16
re-ratification gate.

Runs AFTER nh3_diag_matched_reference.py --n-effective 8000 --seeds 72 172 272
has produced the three per-seed pkls in
validation_reports/nh3_diagnosis/matched_reference_neff8000/.

Pooling: RATIFIED treatment A (plans/scripts/pool_nh3_matched_reference.py /
titanG_repool_reference.py) — concatenate samples AND weights, renormalizing
each seed to 1/n_seeds for equal per-seed posterior mass. Ocean fraction =
sum of pooled normalized weights where D_ocean>0.5 km (nan->0), matching
is_correction.py:251.

Re-ratification gate (reviewer required-validation, R2 adjudication 2026-08-11):
  reference side  : pooled ocean fraction at n_eff=8000, 3 seeds (treatment A);
                    between-seed SE = std(per-seed fracs, ddof1)/sqrt(3).
  corrected side  : C13 3-seed corrected fractions (deployed flow) — between-seed
                    SE = std/sqrt(3). NOT the single-seed R2 bootstrap SE.
  combined SE     = sqrt(SE_ref^2 + SE_corr^2).
  GATE            : |corrected_pooled - reference_pooled| <= 2 * combined SE
                    -> C16 re-ratifies. Otherwise the tension is real (gate FAIL).
DIAGNOSTIC assembly; the scientific-reviewer adjudicates the gate. Never tuned
to pass. Output: matched_reference_neff8000/r3_pool_and_gate.json
"""
import json
import pickle
from pathlib import Path

import numpy as np

ROOT = Path('/Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai')
SEEDS = [72, 172, 272]
NEFF = 8000
SRC = ROOT / 'validation_reports/nh3_diagnosis/matched_reference_neff8000'
OUT = SRC / 'r3_pool_and_gate.json'

# Corrected side — deployed-flow 3-seed corrected ocean fractions (C13), EXACT
# committed values (validation_reports/titan_freegrav_nh3_1m/is_correction/):
#   fiducial (seed 72)  = 0.9292539223381953  (R2 raw / R1 committed fiducial)
#   c13_seed1000        = 0.930961106398094
#   c13_seed2000        = 0.9363014238606263
# Committed between-seed std 0.00367 (is_validation_nh3.json adjudication block);
# committed pass_rule: |corrected-reference| <= 2*sqrt(SE_corr^2+SE_ref^2).
C13_CORRECTED = [0.9292539223381953, 0.930961106398094, 0.9363014238606263]

# committed pooled residual (context) and the n_eff=2000 reference for delta.
COMMITTED_POOLED_RESIDUAL = 0.0149
REF_NEFF2000_POOLED_A = 0.9172547835733885


def ocean_frac(doc, w):
    doc = np.nan_to_num(np.asarray(doc, float), nan=0.0)
    w = np.asarray(w, float)
    return float(w[doc > 0.5].sum() / w.sum())


def main():
    per = []
    for sd in SEEDS:
        p = SRC / f'nh3_reference_seed{sd}_neff{NEFF}.pkl'
        if not p.exists():
            raise SystemExit(f'missing seed pkl: {p}')
        o = pickle.load(open(p, 'rb'))
        s = np.asarray(o.samples)
        w = np.asarray(o.weights, float)
        assert w.shape[0] == s.shape[0], \
            f'seed{sd}: weights {w.shape} != samples {s.shape}'
        doc = np.asarray(o.D_ocean_results, float)
        per.append({'seed': sd, 's': s, 'w': w, 'doc': doc,
                    'frac': ocean_frac(doc, w), 'n': int(s.shape[0])})

    # per-seed weighted ocean fractions
    per_frac = np.array([d['frac'] for d in per])

    # treatment-A pooled fraction (concat samples+weights, each seed -> 1/n_seeds)
    n_seeds = len(per)
    pooled_w = np.concatenate([d['w'] / d['w'].sum() / n_seeds for d in per])
    pooled_doc = np.concatenate([d['doc'] for d in per])
    assert abs(pooled_w.sum() - 1.0) < 1e-9, f'wsum={pooled_w.sum()}'
    ref_pooled = ocean_frac(pooled_doc, pooled_w)

    # reference between-seed SE
    ref_between_std = float(per_frac.std(ddof=1))
    ref_se = ref_between_std / np.sqrt(n_seeds)

    # corrected side (C13 between-seed)
    corr = np.array(C13_CORRECTED)
    corr_pooled = float(corr.mean())
    corr_between_std = float(corr.std(ddof=1))
    corr_se = corr_between_std / np.sqrt(len(corr))

    # gate
    residual = corr_pooled - ref_pooled
    combined_se = float(np.sqrt(ref_se ** 2 + corr_se ** 2))
    gate_bound = 2.0 * combined_se
    reratifies = bool(abs(residual) <= gate_bound)

    report = {
        'purpose': 'R3 pooled n_eff=8000 reference + C16 re-ratification gate',
        'n_effective': NEFF,
        'seeds': SEEDS,
        'reference_side': {
            'per_seed_ocean_fraction': {d['seed']: d['frac'] for d in per},
            'per_seed_n': {d['seed']: d['n'] for d in per},
            'pooled_ocean_fraction_treatmentA': ref_pooled,
            'between_seed_std_ddof1': ref_between_std,
            'between_seed_SE': ref_se,
            'delta_vs_neff2000_pooled': ref_pooled - REF_NEFF2000_POOLED_A,
        },
        'corrected_side': {
            'c13_fractions': C13_CORRECTED,
            'pooled_mean': corr_pooled,
            'between_seed_std_ddof1': corr_between_std,
            'between_seed_SE': corr_se,
        },
        'gate': {
            'residual_corrected_minus_reference': residual,
            'combined_SE': combined_se,
            'two_sigma_bound': gate_bound,
            'residual_over_combined_SE': (residual / combined_se
                                          if combined_se > 0 else float('inf')),
            're_ratifies_C16': reratifies,
            'committed_pooled_residual_context': COMMITTED_POOLED_RESIDUAL,
        },
        'reading': (
            'GATE: |corrected_pooled - reference_pooled| <= 2*combined_SE '
            '=> C16 re-ratifies. Both SEs are BETWEEN-SEED pooled SEs (reviewer '
            'R2 required-validation). If the n_eff=8000 reference moves the '
            'reference materially UPWARD (delta_vs_neff2000 positive and large), '
            'the residual shrinks toward the gate; if it stays near the '
            'n_eff=2000 value, the residual persists and the tension is real. '
            'Scientific-reviewer adjudicates; never tuned to pass. '
            'Does NOT re-ratify C16 by itself — the reviewer + manager do.'),
    }
    OUT.write_text(json.dumps(report, indent=2, default=str))
    print(json.dumps(report, indent=2, default=str))
    print(f'\nreport -> {OUT}')


if __name__ == '__main__':
    main()
