"""§0.20 R4 — record the biased estimator FAIL-ADJUDICATED in the C16 block.

Adds a `r4_resolution_record` sub-block to gates.ocean_fraction. Does NOT flip
`pass` or clear the MANAGER-GATE STOP in `adjudication` — that status change is
the manager's (Machine A) action, per project discipline. This records the R1->R3
outcome + reviewer verdict as provenance so the committed adjudication file is
self-describing when the manager acts.
"""
import json
from pathlib import Path

P = Path('validation_reports/titan_freegrav_nh3_1m/is_correction/is_validation_nh3.json')
d = json.loads(P.read_text())
of = d['gates']['ocean_fraction']

of['r4_resolution_record'] = {
    'sequence': 'section 0.20 R1->R2->R3->R4',
    'biased_estimator': 'reference side (n_eff=2000 pooled MCMC), biased LOW',
    'verdict': 'FAIL-ADJUDICATED',
    'diagnosis': (
        'R3 (n_eff=8000, 3 seeds 72/172/272, all R-hat=1.000, ESS~12.9k) '
        'raised the pooled reference ocean fraction +0.0114 (0.91725->0.92865); '
        'the corrected side (R2, N=100k) is a precisely-resolved fixed point '
        '(~0.933) that moved AWAY from the reference at higher N, falsifying '
        'finite-N corrected-side (SNIS) bias. Therefore the ~0.015 residual was '
        'a REFERENCE-side sampling-resolution artifact: the n_eff=2000 reference '
        'undersampled the higher-dissipation ocean branch (ocean fraction and '
        'weighted |Im_k2| median rose together; |Im_k2| between-seed std '
        'collapsed 0.00066->0.00015).'),
    'gate_at_neff8000': {
        'reference_pooled': 0.9286503,
        'reference_between_seed_SE': 0.0020338,
        'corrected_pooled_C13': 0.9321722,
        'corrected_between_seed_SE': 0.0021226,
        'residual_corrected_minus_reference': 0.0035218,
        'combined_SE': 0.0029397,
        'two_sigma_bound': 0.0058794,
        'residual_over_combined_SE': 1.198,
        're_ratifies_C16': True,
    },
    'reviewer': (
        'scientific-reviewer 2026-08-11 (agent a8d0ea375f5d05048) '
        'PASS-WITH-CONCERNS: independently reproduced the gate to the last digit '
        'from the raw n_eff=8000 seed pkls; confirmed the reference move is '
        'like-for-like (identical config_hash + structure_cache_sha256, '
        'regime-preserving n_active/n_eff ~2); seed-272 inclusion is CONSERVATIVE '
        '(dropping it shrinks the residual to +0.00149). 4 required-validations '
        '(non-blocking) folded into R3_decisive_reference_recompute.md; all point '
        'toward stronger agreement. Verdict: re-ratifies C16, supports releasing '
        'the STOP.'),
    'artifacts': [
        'validation_reports/nh3_diagnosis/R1_ocean_fraction_audit.md',
        'validation_reports/nh3_diagnosis/R2_corrected_100k.md',
        'validation_reports/nh3_diagnosis/R3_decisive_reference_recompute.md',
        'validation_reports/nh3_diagnosis/matched_reference_neff8000/r3_pool_and_gate.json',
    ],
    'status_change_owner': (
        'Flipping `pass` to true and clearing the MANAGER-GATE STOP in '
        '`adjudication` is the manager (Machine A) action, NOT recorded here. '
        'This block records the R1->R3 outcome + reviewer verdict only; the '
        'gate arithmetic passes and the reviewer re-ratifies, but final STOP '
        'release + status flip remain Machine A\'s.'),
}

P.write_text(json.dumps(d, indent=1))
print('R4 record written to', P)
print('pass flag left unchanged:', of['pass'])
