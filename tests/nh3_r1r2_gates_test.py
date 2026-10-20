"""Unit tests for the pure helpers of plans/scripts/nh3_r1r2_gates.py."""
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / 'plans/scripts'))
sys.path.insert(0, str(REPO))
import nh3_r1r2_gates as g  # noqa: E402


def test_c9_stock_pass():
    assert g.c9_decision(True, True, 99.0, 0.0, 0.1) == 'PASS'


def test_c9_pass_by_c9():
    assert g.c9_decision(False, True, 1.05, 1.0, 0.1) == 'PASS_BY_C9'
    assert g.c9_decision(False, True, 1.125, 1.0, 0.125) == 'PASS_BY_C9'  # boundary (exact fp)


def test_c9_fail_modes():
    assert g.c9_decision(False, False, 1.0, 1.0, 0.1) == 'FAIL'   # width wrong
    assert g.c9_decision(False, True, 1.2, 1.0, 0.1) == 'FAIL'    # median off


def _all(v='PASS'):
    return {n: v for n in g.C9_NUISANCE + g.NON_C9}


def test_outcome_lift():
    o = g.outcome_mapping(_all(), True, True)
    assert o['verdict'] == 'LIFT_TAKES_EFFECT' and o['all_13_pass']


def test_outcome_c9_fail_quarantines_only():
    s = _all(); s['log10_zeta'] = 'FAIL'; s['alpha'] = 'PASS_BY_C9'
    o = g.outcome_mapping(s, True, True)
    assert o['verdict'] == 'LIFT_TAKES_EFFECT'
    assert o['quarantined_by_name'] == ['log10_zeta']
    assert o['c9_pass_by_c9'] == ['alpha'] and not o['all_13_pass']


def test_outcome_non_c9_fail_blocks():
    s = _all(); s['Tb_K'] = 'FAIL'
    o = g.outcome_mapping(s, True, True)
    assert o['verdict'] == 'NO_LIFT_REESCALATE' and o['non_c9_fails'] == ['Tb_K']


def test_outcome_r1b_r2_block():
    assert g.outcome_mapping(_all(), False, True)['verdict'] == 'NO_LIFT_REESCALATE'
    assert g.outcome_mapping(_all(), True, False)['verdict'] == 'NO_LIFT_REESCALATE'


def test_pool_and_weights():
    W = g.pool_treatment_A([np.array([1., 3.]), np.array([2., 2., 4.])])
    assert abs(W.sum() - 1) < 1e-15
    assert abs(W[:2].sum() - 0.5) < 1e-12
    w = g.weights_from_log_w(np.array([0., -np.inf, np.log(3.)]))
    assert np.allclose(w, [0.25, 0, 0.75])


def test_max_rel_diff():
    assert g.max_rel_diff([1., 2.], [1., 2.]) == 0.0
    assert g.max_rel_diff([np.nan], [1.]) == float('inf')
