"""§0.20 R1 — reference-side ocean-fraction weighting audit (ZERO COMPUTE).

Reads ONLY the committed NH3 matched-resolution reference seed pkls
(72/172/272, n_eff=2000) and reports how the reference ocean fraction
(the C16 comparator) depends on the weighting treatment:

  per seed  : weighted (pocoMC importance weights) vs unweighted (uniform);
              Kish ESS; ocean fraction each way.
  pooled    : (A) CORRECT pooling  = concat samples + weights, each seed
                  renormalised to 1/n_seeds  [what the committed reference uses]
              (B) PRE-REPAIR bug   = last seed's weight array against ALL
                  pooled samples (the titanG_repool_reference.py root-cause;
                  length-mismatch -> here emulated by tiling/truncation the
                  same way the buggy pooled pkl carried only one seed's w)
              (C) UNWEIGHTED pooled = uniform over all pooled samples.

Decision context (manager preregistration): the corrected-vs-reference
residual is +0.0149. If the REFERENCE ocean fraction itself swings by an
amount comparable to that residual across these treatments, the comparator
is weighting-fragile and the reference side is implicated.

No MCMC, no flow, no Test/ writes. Output: /tmp/r1_reference_audit.json
"""
import json
import pickle
from pathlib import Path

import numpy as np

ROOT = Path('/Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai')
SRC = ROOT / 'validation_reports/nh3_diagnosis/matched_reference'
SEEDS = [72, 172, 272]
CORRECTED_FIDUCIAL = 0.9292539223381953  # is_validation_nh3.json branch.ocean.prob_corrected
CORRECTED_RESIDUAL = 0.0149               # committed pooled residual (manager escalation)


def kish_ess(w):
    w = np.asarray(w, float)
    return float((w.sum() ** 2) / np.sum(w ** 2))


def ocean_frac(doc, w):
    doc = np.nan_to_num(np.asarray(doc, float), nan=0.0)
    w = np.asarray(w, float)
    return float(w[doc > 0.5].sum() / w.sum())


def main():
    seeds = {}
    per_s, per_w, per_doc = [], [], []
    for sd in SEEDS:
        o = pickle.load(open(SRC / f'nh3_reference_seed{sd}_neff2000.pkl', 'rb'))
        s = np.asarray(o.samples)
        w = np.asarray(o.weights, float)
        doc = np.nan_to_num(np.asarray(o.D_ocean_results, float), nan=0.0)
        assert w.shape[0] == s.shape[0] == doc.shape[0]
        uni = np.ones_like(w) / len(w)
        seeds[sd] = {
            'n': int(len(w)),
            'kish_ess': kish_ess(w),
            'ess_over_n': kish_ess(w) / len(w),
            'ocean_frac_weighted': ocean_frac(doc, w),
            'ocean_frac_unweighted': ocean_frac(doc, uni),
        }
        per_s.append(s)
        per_w.append(w)
        per_doc.append(doc)

    n_seeds = len(SEEDS)
    pooled_doc = np.concatenate(per_doc)
    pooled_n = int(len(pooled_doc))

    # (A) CORRECT pooling: each seed renormalised to 1/n_seeds
    wA = np.concatenate([w / w.sum() / n_seeds for w in per_w])
    fracA = ocean_frac(pooled_doc, wA)

    # (C) UNWEIGHTED pooled: uniform over all pooled samples
    wC = np.ones(pooled_n) / pooled_n
    fracC = ocean_frac(pooled_doc, wC)

    # (B) PRE-REPAIR bug: the buggy pool carried only the LAST seed's weight
    # array against all samples. Lengths mismatched (last seed 6230 vs pooled
    # 20604), so validate would have raised; but to characterise the STATISTIC
    # the bug implied, emulate "one seed's weights broadcast over the pool":
    # take the last seed's normalised weights, tile/truncate to pooled length,
    # renormalise. This is the closest well-defined reading of the buggy state.
    w_last = per_w[-1] / per_w[-1].sum()
    reps = int(np.ceil(pooled_n / len(w_last)))
    wB = np.tile(w_last, reps)[:pooled_n]
    wB = wB / wB.sum()
    fracB = ocean_frac(pooled_doc, wB)

    # spread of the three per-seed weighted fractions (the committed C13/P1.1
    # reference spread lives on these same seeds)
    seed_fracs = [seeds[sd]['ocean_frac_weighted'] for sd in SEEDS]
    seed_spread = float(np.max(seed_fracs) - np.min(seed_fracs))
    seed_std = float(np.std(seed_fracs, ddof=1))

    treatments = {
        'A_correct_pooled_weighted': fracA,
        'B_prerepair_last_seed_weights': fracB,
        'C_unweighted_pooled': fracC,
    }
    max_treatment = max(treatments.values())
    min_treatment = min(treatments.values())
    treatment_swing = float(max_treatment - min_treatment)

    report = {
        'purpose': 'R1 reference-side ocean-fraction weighting sensitivity '
                   '(zero compute, saved seed pkls only)',
        'seeds': seeds,
        'pooled_n': pooled_n,
        'pooled_ocean_fraction': treatments,
        'committed_reference_value': fracA,
        'corrected_fiducial': CORRECTED_FIDUCIAL,
        'corrected_minus_referenceA': float(CORRECTED_FIDUCIAL - fracA),
        'committed_residual': CORRECTED_RESIDUAL,
        'per_seed_weighted_fraction_spread': seed_spread,
        'per_seed_weighted_fraction_std_ddof1': seed_std,
        'treatment_swing_maxmin': treatment_swing,
        'weighting_fragile': bool(treatment_swing >= 0.5 * CORRECTED_RESIDUAL),
        'reading': (
            'If treatment_swing is comparable to the +0.0149 residual, the '
            'reference comparator is weighting-fragile and the reference side '
            'is implicated. If the correct-pooled (A) and unweighted (C) '
            'fractions nearly coincide, the pocoMC importance weights are '
            'near-uniform on the ocean/no-ocean split and the residual is NOT '
            'a reference-weighting artifact.')
    }
    out = Path('/tmp/r1_reference_audit.json')
    out.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    print(f'\nreport -> {out}')


if __name__ == '__main__':
    main()
