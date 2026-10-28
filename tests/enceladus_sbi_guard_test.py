"""Required refusal-guard test (reviewer adjudication 2026-10-27,
validation_reports/enceladus_isostasy/sbi/RESULT.md): a thin-shard x is
REFUSED; the Cassini datum x_obs is ACCEPTED with 0 % dropped and marked
validated; a non-validated accepted x carries the not-validated label.

Needs the production cache and the trained artifact (Mac scratch); skips
otherwise."""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
SCR = Path(os.environ.get('PP_SCRATCH', '/Users/svance/scratch'))
CACHE = SCR / 'enceladus_zbw_production_cache.pkl'
ART = SCR / 'enc_sbi' / 'enceladus_ocean6D_rhorock_posterior_1m.pt'
THIN = SCR / 'enc_sbi' / 'thin_shard_00.npz'
CFG = REPO / 'PlanetProfile/Inference/configs/enceladus_cassini_isostasy_ocean6D_rhorock.json'

pytestmark = pytest.mark.skipif(
    not (CACHE.exists() and ART.exists() and THIN.exists()),
    reason='production cache / Enceladus SBI artifact not present')


@pytest.fixture(scope='module')
def runners():
    from PlanetProfile.Inference.inference_core import InferenceConfig
    from PlanetProfile.Inference.mcmc_runner import MCMCRunner
    from PlanetProfile.Inference.sbi_runner import SBIRunner
    raw = json.loads(CFG.read_text())
    raw['structure_cache_path'] = str(CACHE)
    mr = MCMCRunner(InferenceConfig(**raw))
    sr = SBIRunner.load_artifact(str(ART))
    obs = {k: float(v[0]) for k, v in raw['observables'].items()}
    return sr, mr, obs, list(raw['observables'])


def test_x_obs_accepted_validated_nothing_dropped(runners):
    from PlanetProfile.Inference.enceladus_sbi_guard import guarded_posterior
    sr, mr, obs, _ = runners
    r = guarded_posterior(sr, mr, obs, validated_x=[obs], n_draws=4000)
    assert r, getattr(r, 'reason', '')
    assert r.validated and r.dropped_fraction == 0.0 and r.frac_thin == 0.0
    assert len(r.samples) == 4000


def test_thin_shard_x_refused(runners):
    from PlanetProfile.Inference.enceladus_sbi_guard import (guarded_posterior,
                                                             Refused)
    sr, mr, obs, names = runners
    d = np.load(THIN)
    i = int(np.argmin(d['theta'][:, 0]))          # thinnest truth in the set
    x = dict(zip(names, map(float, d['x'][i])))
    r = guarded_posterior(sr, mr, x, validated_x=[obs], n_draws=4000)
    assert isinstance(r, Refused) and not r
    assert r.frac_union > 0.01
    assert 'use pocoMC' in r.reason


def test_accepted_non_validated_x_is_labelled(runners):
    from PlanetProfile.Inference.enceladus_sbi_guard import (
        guarded_posterior, NOT_VALIDATED_LABEL)
    sr, mr, obs, _ = runners
    x = dict(obs)
    x['C30'] = obs['C30'] * 1.05                   # off-datum, not a G4 set
    r = guarded_posterior(sr, mr, x, validated_x=[obs], n_draws=4000)
    assert r and not r.validated and r.label == NOT_VALIDATED_LABEL
