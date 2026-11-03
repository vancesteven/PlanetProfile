"""k2 campaign correction C1 (scientific-reviewer 2026-11-02): the result of
apply_parameters must not depend on the order of theta's keys. The ice-Ih
shear-modulus hook used to force rheology_type = 'maxwell', so with Andrade
parameters present {'alpha', 'log10_mu_Ih'} ran Maxwell and
{'log10_mu_Ih', 'alpha'} ran Andrade."""
import itertools
import pickle
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
CACHE = REPO / ('PlanetProfile/Test/mcmc_results/Enceladus/Cassini_isostasy/'
                'enceladus_seawater_zbw_grid_2d.pkl')
pytestmark = pytest.mark.skipif(not CACHE.exists(), reason='cache absent')


@pytest.fixture(scope='module')
def node():
    return pickle.load(open(CACHE, 'rb'))['structures'][1220]


THETA = {'alpha': 0.25, 'log10_zeta': 0.0, 'log10_mu_Ih': np.log10(2.5e9),
         'log10_eta_Ih': 14.0}


def test_rheology_independent_of_hook_order(node):
    from PlanetProfile.Inference.forward_models import apply_parameters
    seen = set()
    for perm in itertools.permutations(THETA):
        out = apply_parameters({k: THETA[k] for k in perm}, node)
        seen.add(out.get('rheology_type'))
        assert out['mu_Ih_override'] == pytest.approx(2.5e9)
    assert seen == {'andrade'}


def test_k2_independent_of_hook_order(node):
    from PlanetProfile.Inference.forward_models import forward_model_k2_flexible
    ks = set()
    for perm in (tuple(THETA), tuple(reversed(tuple(THETA)))):
        Re, Im, *_ = forward_model_k2_flexible({k: THETA[k] for k in perm}, node,
                                               return_heating=False)
        ks.add((round(Re, 12), round(Im, 12)))
    assert len(ks) == 1


def test_mu_only_still_runs_maxwell(node):
    """No Andrade parameter: the solver default (Maxwell) is unchanged, so
    configs that sample log10_mu_Ih without Andrade behave as before."""
    from PlanetProfile.Inference.forward_models import apply_parameters
    out = apply_parameters({'log10_mu_Ih': np.log10(3e9)}, node)
    assert out.get('rheology_type', 'maxwell') == 'maxwell'


def test_andrade_and_maxwell_differ_at_soft_shell(node):
    """Guards that the rheology choice actually reaches the solver."""
    from PlanetProfile.Inference.forward_models import forward_model_k2_flexible
    base = {'log10_mu_Ih': np.log10(3.3e9), 'log10_eta_Ih': 14.0}
    Re_m, Im_m, *_ = forward_model_k2_flexible(dict(base), node, return_heating=False)
    Re_a, Im_a, *_ = forward_model_k2_flexible(
        {**base, 'alpha': 0.2, 'log10_zeta': 0.0}, node, return_heating=False)
    assert abs(Re_a - Re_m) > 1e-4 or abs(Im_a - Im_m) > 1e-4
