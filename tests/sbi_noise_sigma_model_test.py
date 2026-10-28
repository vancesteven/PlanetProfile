"""generate_sbi_dataset(obs_noise=True) must draw noise from the SAME
covariance the likelihood uses: config sigma inflated by
metadata.sigma_model_add (Enceladus B15) and correlated per
metadata.observable_correlations (found 2026-10-27: it used the raw sigma).

Exercises the real generator twice on the same theta (theta_override), once
noise-free and once noisy, so the difference is exactly the injected noise.
Uses the Enceladus smoke cache + campaign observables (same fixture as
tests/enceladus_sbi_dataset_isostatic_test.py)."""
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'tests'))
from PlanetProfile.Inference.mcmc_runner import MCMCRunner  # noqa: E402
import enceladus_sbi_dataset_isostatic_test as E  # noqa: E402

pytestmark = pytest.mark.skipif(not E.SMOKE_CACHE.exists(),
                                reason='Enceladus smoke cache not present')

# Test-only inflation EQUAL to each observable's sigma, so the inflated sd is
# sqrt(2)x the raw one -- large enough that the check discriminates (the real
# campaign inflation, <=1.1%, would not be resolvable at this n).
SMA = {'C20': 36.99e-6, 'C22': 14.70e-6, 'C30': 33.42e-6, 'libration_deg': 0.003}


def _runner(sma):
    cfg = E._config()
    md = dict(cfg.metadata)
    md['observable_correlations'] = {'C20,C22': 0.47}
    if sma:
        md['sigma_model_add'] = dict(sma)
    else:
        md.pop('sigma_model_add', None)
    cfg.metadata = md
    return MCMCRunner(cfg), cfg


def _noise(sma, n=4000):
    r, cfg = _runner(sma)
    th0, _ = r.generate_sbi_dataset(n_samples=60, seed=3, drop_nonfinite=True,
                                    imag_convention='abs')
    th = np.repeat(th0[:1], n, axis=0)
    _, x0 = r.generate_sbi_dataset(theta_override=th, drop_nonfinite=True,
                                   imag_convention='abs')
    _, x1 = r.generate_sbi_dataset(theta_override=th, drop_nonfinite=True,
                                   imag_convention='abs', obs_noise=True,
                                   noise_seed=7)
    sig = np.array([cfg.observables[k][1] for k in E.OBS_NAMES], float)
    return x1 - x0, sig


def test_noise_uses_inflated_sigma_and_correlation():
    e, sig = _noise(SMA)
    target = np.sqrt(sig ** 2 + np.array([SMA[k] for k in E.OBS_NAMES]) ** 2)
    sd = e.std(axis=0, ddof=1)
    np.testing.assert_allclose(sd, target, rtol=0.05)
    rho = np.corrcoef(e[:, 0], e[:, 1])[0, 1]
    assert abs(rho - 0.47) < 0.05
    # raw sigma would give sd/target = 1/sqrt(2) = 0.71: clearly rejected.
    assert np.all(sd / target > 0.9)


def test_no_sigma_model_add_keeps_raw_sigma():
    e, sig = _noise(None)
    np.testing.assert_allclose(e.std(axis=0, ddof=1), sig, rtol=0.05)
