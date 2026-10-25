"""Light tests for plans/scripts/enceladus_fisher_identifiability.py:
linear synthetic forward model -> Fisher inverse == analytic covariance;
a degenerate pair is flagged; a null (insensitive) parameter is reported
unconstrained. No PlanetProfile forward model is involved."""
import importlib.util
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location(
    'fisher_id', REPO / 'plans/scripts/enceladus_fisher_identifiability.py')
fi = importlib.util.module_from_spec(_spec)
sys.modules['fisher_id'] = fi
_spec.loader.exec_module(fi)

NAMES = ['a', 'b', 'c', 'd']


def _cov():
    s = np.array([0.5, 0.2, 1.0, 0.1])
    cov = np.diag(s ** 2)
    cov[0, 1] = cov[1, 0] = 0.47 * s[0] * s[1]
    return cov


def test_fisher_inverse_matches_analytic_covariance():
    rng = np.random.default_rng(3)
    A = rng.normal(size=(4, 4)) + 3 * np.eye(4)          # y = A x, full rank
    cov = _cov()
    F = fi.fisher_from_jacobian(A, cov)
    Ci = np.linalg.inv(A) @ cov @ np.linalg.inv(A).T      # exact cov of x_hat
    inv = fi.invert_fisher(F)
    assert inv['rank'] == 4 and not inv['inf_mask'].any()
    assert np.allclose(inv['C'], Ci, rtol=1e-9, atol=1e-12)
    out = fi.analyze_fisher(F, NAMES, g_cmr2=np.array([1., 0, 0, 0]),
                            i_rho_s=1, i_C2=2)
    assert np.allclose([out['sigma'][n] for n in NAMES], np.sqrt(np.diag(Ci)))
    assert abs(out['sigma_cmr2_derived'] - np.sqrt(Ci[0, 0])) < 1e-12
    # rho_s-fixed == conditional: delete row/col of F, then invert
    keep = [0, 2, 3]
    Cf = np.linalg.inv(F[np.ix_(keep, keep)])
    assert np.allclose([out['sigma_rho_s_fixed'][NAMES[k]] for k in keep],
                       np.sqrt(np.diag(Cf)))


def test_finite_difference_jacobian_exact_for_linear_model():
    A = np.array([[1., 2., 0.], [0., 1., -1.], [3., 0., 1.], [1., 1., 1.]])
    f = lambda x: A @ x
    x0 = np.array([0.3, 0.5, -0.2])
    J, kinds = fi.jacobian(f, x0, [0.1, 0.1, 0.1], [(-1, 1)] * 3)
    assert np.allclose(J, A) and set(kinds) == {'central'}
    J, kinds = fi.jacobian(f, x0, [0.1, 0.1, 0.1], [(0.3, 1)] + [(-1, 1)] * 2)
    assert kinds[0] == 'forward' and np.allclose(J, A)


def test_degenerate_pair_is_flagged_and_covariance_matches_correlated_gaussian():
    # y = [x1 + x2 (+eps tilt), x3, ...]: x1, x2 nearly collinear
    eps = 1e-3
    J = np.array([[1., 1., 0.], [1., 1. + eps, 0.], [0., 0., 1.], [0., 0., 2.]])
    cov = np.diag([0.1, 0.1, 0.2, 0.2]) ** 2
    F = fi.fisher_from_jacobian(J, cov)
    out = fi.analyze_fisher(F, ['rho_s', 'C2', 'x'], i_rho_s=0, i_C2=1)
    assert abs(out['corr_rho_s_C2']) > 0.98
    assert out['flag_C2_not_separately_constrained']
    # a clearly separable pair is not flagged
    J2 = np.eye(3)[:3]
    out2 = fi.analyze_fisher(fi.fisher_from_jacobian(J2, cov[:3, :3]),
                             ['rho_s', 'C2', 'x'], i_rho_s=0, i_C2=1)
    assert not out2['flag_C2_not_separately_constrained']


def test_singular_fisher_reports_inf_not_garbage():
    # 2 observables, 3 parameters -> rank 2, one null direction
    J = np.array([[1., 1., 0.], [0., 1., 1.]])
    F = fi.fisher_from_jacobian(J, np.eye(2) * 0.01)
    out = fi.analyze_fisher(F, ['a', 'b', 'c'], i_rho_s=0, i_C2=2)
    assert out['rank'] == 2 and out['nullity'] == 1
    assert all(out['unconstrained'].values())
    assert out['sigma']['a'] == 'inf'
    # adding a prior on 'a' makes it identifiable
    out2 = fi.analyze_fisher(F + np.diag([1., 0, 0]), ['a', 'b', 'c'])
    assert out2['rank'] == 3 and not any(out2['unconstrained'].values())


def test_zero_sensitivity_column_is_unconstrained_not_poisoning():
    J = np.array([[1., 0., 0.5], [0., 0., 1.], [2., 0., 0.]])
    F = fi.fisher_from_jacobian(J, np.eye(3) * 0.04)
    out = fi.analyze_fisher(F, ['a', 'zb', 'c'])
    assert out['unconstrained']['zb'] and out['sigma']['zb'] == 'inf'
    assert out['sigma']['a'] != 'inf' and out['sigma']['c'] != 'inf'


def test_covariance_builder_matches_likelihood_convention():
    obs = {'C20': (-0.0055, 3.7e-5), 'C22': (0.0015, 1.5e-5),
           'C30': (1.8e-4, 3.3e-5), 'libration_deg': (0.092, 0.003)}
    add = {'C20': 5.3e-6, 'C22': 1.7e-6, 'C30': 4.4e-6, 'libration_deg': 2.5e-4}
    cov, s = fi.build_covariance(obs, add, 0.47)
    assert np.isclose(s[0], np.hypot(3.7e-5, 5.3e-6))
    assert np.isclose(cov[0, 1], 0.47 * s[0] * s[1])
    assert cov[0, 2] == 0 and cov[3, 0] == 0
