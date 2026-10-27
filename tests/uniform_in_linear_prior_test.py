"""prior_type 'uniform-in-linear' (Enceladus run (g1): w uniform in ppt, sampled as log10 w)."""
import numpy as np
import pytest
from scipy import integrate

from PlanetProfile.Inference.mcmc_runner import _UniformInLinearOfLog10
from PlanetProfile.Inference.inference_core import validate_config, InferenceConfig  # noqa: F401


def test_density_normalized_and_proportional_to_10x():
    d = _UniformInLinearOfLog10(-1.0, 2.0)
    tot, _ = integrate.quad(d.pdf, -1.0, 2.0)
    assert tot == pytest.approx(1.0, rel=1e-10)
    x = np.array([-0.5, 0.0, 1.0, 1.9])
    np.testing.assert_allclose(d.pdf(x) / d.pdf(0.0), 10.0 ** x, rtol=1e-12)
    assert np.all(np.isneginf(d.logpdf(np.array([-1.01, 2.01]))))
    assert d.support() == (-1.0, 2.0)


def test_samples_are_uniform_in_linear_value():
    d = _UniformInLinearOfLog10(-1.0, 2.0)
    np.random.seed(0)
    v = 10.0 ** d.rvs(size=200_000)
    assert v.min() >= 0.1 - 1e-12 and v.max() <= 100.0 + 1e-9
    # uniform on [0.1, 100]: mean 50.05, and each decile holds ~10%
    assert v.mean() == pytest.approx(50.05, rel=5e-3)
    counts, _ = np.histogram(v, bins=10, range=(0.1, 100.0))
    np.testing.assert_allclose(counts / v.size, 0.1, atol=3e-3)


def test_rvs_uses_numpy_global_state_like_scipy():
    """pocoMC calls dist.rvs(size=n) with no random_state, so seeding numpy's
    global RandomState must make draws reproducible."""
    d = _UniformInLinearOfLog10(-1.0, 1.3)
    np.random.seed(75)
    a = d.rvs(size=50)
    np.random.seed(75)
    b = d.rvs(size=50)
    np.testing.assert_array_equal(a, b)


def test_cdf_ppf_inverse():
    d = _UniformInLinearOfLog10(-1.0, 1.3)
    q = np.linspace(0, 1, 11)
    np.testing.assert_allclose(d.cdf(d.ppf(q)), q, atol=1e-12)
