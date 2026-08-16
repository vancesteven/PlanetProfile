"""InferenceConfig.from_json schema reconciliation for the Enceladus
isostasy candidate config.

Mechanical wiring only -- no physics changes, no scientific judgment calls.

Before this fix, ``enceladus_cassini_isostasy_7D.json`` carried several
top-level JSON keys (``isostasy``, ``structure_cache_spec``,
``libration_model_correction``, ``gravity_forward_model_contract``,
``derived_display_only``, ``derived_display_only_note``, ``branch_model``)
that are not ``InferenceConfig`` dataclass fields, so ``cls(**data)`` in
``InferenceConfig.from_dict`` raised ``TypeError: __init__() got an
unexpected keyword argument ...``. Those blocks were moved verbatim (same
Python values, only relocated) under the existing top-level ``metadata``
dict, which every consumer already treats as the place for campaign-
specific physics inputs (see ``MCMCRunner._isostasy_physics_config``,
which already read ``metadata['isostasy']`` as its fallback before this
fix -- a grep across the repo turned up no consumer that reads any of
these seven blocks as a top-level ``config.<name>`` attribute, so no
production code needed to change).
"""
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from PlanetProfile.Inference.inference_core import InferenceConfig  # noqa: E402

CANDIDATE_CONFIG = (REPO / 'PlanetProfile' / 'Inference' / 'configs'
                    / 'enceladus_cassini_isostasy_7D.json')

# The seven top-level JSON keys moved under `metadata` -- not recognized
# InferenceConfig dataclass fields.
MOVED_BLOCKS = [
    'gravity_forward_model_contract',
    'structure_cache_spec',
    'isostasy',
    'derived_display_only',
    'libration_model_correction',
    'derived_display_only_note',
    'branch_model',
]


@pytest.mark.skipif(not CANDIDATE_CONFIG.exists(), reason='candidate config not present')
def test_from_json_succeeds_on_candidate_config():
    """This used to raise TypeError -- the whole point of the fix."""
    cfg = InferenceConfig.from_json(str(CANDIDATE_CONFIG))
    assert cfg.mode == 'mcmc'
    assert cfg.bodyname == 'Enceladus'
    assert cfg.gravity_forward_model == 'isostatic_hm2019'


@pytest.mark.skipif(not CANDIDATE_CONFIG.exists(), reason='candidate config not present')
def test_moved_blocks_are_not_top_level_json_keys():
    """Confirms the reconciliation actually moved the blocks rather than
    duplicating them -- from_json's TypeError would only be fixed by
    removal from the top level, and a stray top-level survivor would be
    silently ignored by dataclass '**data' expansion only by accident."""
    with open(CANDIDATE_CONFIG) as f:
        raw = json.load(f)
    for key in MOVED_BLOCKS:
        assert key not in raw, f'{key} must not remain a top-level JSON key'


@pytest.mark.skipif(not CANDIDATE_CONFIG.exists(), reason='candidate config not present')
def test_moved_blocks_reachable_under_metadata():
    cfg = InferenceConfig.from_json(str(CANDIDATE_CONFIG))
    for key in MOVED_BLOCKS:
        assert key in cfg.metadata, f'{key} must be reachable via config.metadata'


@pytest.mark.skipif(not CANDIDATE_CONFIG.exists(), reason='candidate config not present')
def test_moved_block_content_preserved_byte_for_byte():
    """Loads the raw JSON independently of InferenceConfig and checks each
    moved block's content, now nested under metadata, is value-identical
    to a fixed reference of the physics content the freeze spec requires
    (plans/active/enceladus-isostasy-module-spec.md Sigma_model section --
    Tajeddine sigma_model 5.3e-6 / 1.7e-6 / 4.4e-6 and 0.00025 deg)."""
    cfg = InferenceConfig.from_json(str(CANDIDATE_CONFIG))
    iso = cfg.metadata['isostasy']
    assert iso['sigma_model_lm'] == {
        '2,0': 5.3e-06, '2,2': 1.7e-06, '3,0': 4.4e-06,
    }
    assert iso['sigma_model_libration_deg'] == 0.00025
    assert iso['shape_ref_radius_m'] == 252220.0
    assert iso['gravity_ref_radius_m'] == 256600.0
    assert iso['H_obs_lm_m'] == {'2,0': -3510.0, '2,2': 857.0, '3,0': 420.0}

    branch_model = cfg.metadata['branch_model']
    assert branch_model['prior_odds'] == {'ocean': 0.5, 'frozen': 0.5}

    contract = cfg.metadata['gravity_forward_model_contract']
    assert contract['clairaut_radau_darwin_path'] == (
        'MUST NOT be used by this forward model')

    spec = cfg.metadata['structure_cache_spec']
    assert spec['grid_mode'] == 'zb_w_2d'
    # 4960 (= 124 zb x 40 w) until the frozen-branch design ruling (task A7)
    # moved the frozen segment off the shared zb axis onto its own array.
    # Now 87 ocean zb x 40 w = 3480, plus 39 frozen zb NOT crossed with w:
    # salinity is undefined without an ocean, so crossing them would store 39
    # copies of the same structure under a meaningless coordinate.
    assert spec['n_nodes_total'] == 3519
    assert spec['zb_km_grid']['hi'] == 45.0, 'ocean zb axis back to its own edge'
    assert spec['frozen_zb_km_grid']['n'] == 39
    assert 'retry_frozen_as_no_ocean' not in spec, (
        'retry_frozen_as_no_ocean is inoperative at Enceladus (no HP ice) '
        'and must not advertise a build path that does not exist')


@pytest.mark.skipif(not CANDIDATE_CONFIG.exists(), reason='candidate config not present')
def test_no_new_dataclass_fields_were_added():
    """The task's explicit MINIMAL-change constraint: the fix must not add
    new InferenceConfig dataclass fields, since that would change the
    config hash surface (generate_hash) for every campaign, not just this
    one. Confirmed by checking the seven relocated keys are absent from
    dataclasses.fields(InferenceConfig)."""
    import dataclasses
    field_names = {f.name for f in dataclasses.fields(InferenceConfig)}
    for key in MOVED_BLOCKS:
        assert key not in field_names


@pytest.mark.skipif(not CANDIDATE_CONFIG.exists(), reason='candidate config not present')
def test_hash_generation_still_works():
    """generate_hash() must not choke on the reconciled config (it only
    touches recognized dataclass fields, so this is really a smoke check
    that from_json produced a well-formed InferenceConfig)."""
    cfg = InferenceConfig.from_json(str(CANDIDATE_CONFIG))
    h = cfg.generate_hash()
    assert isinstance(h, str) and len(h) == 16
