"""D7 + the r5 freeze-gate double-count correction + D8's budget entry.

r5 (validation_reports/enceladus_isostasy/r5_ADJUDICATION.md):

  D7  ``rho_rock_kgm3`` missing from PARAMETER_REGISTRY -- the FROZEN
      branch's sampling coordinate, and the only one of the campaign's five
      still unregistered, so the frozen branch could not be driven from a
      config at all.
  D8  the first-order H22 term (0.65 sigma one-sided on C22) unbudgeted,
      with Sigma_model's C22 row 5.6x smaller than it.
  FREEZE-GATE CORRECTION (binding): the Tajeddine 0.08 sigma is ALREADY in
      Sigma_model as ``libration_deg`` 0.00025 deg, so it must NOT also
      enter the ``libration_sys_frac`` width -- that would count one shape-
      model uncertainty twice on one channel.

These are record-level repairs (a registry entry and config metadata); the
physics is unchanged, which is exactly what the last test asserts.
"""
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from PlanetProfile.Inference.parameter_registry import (  # noqa: E402
    PARAMETER_REGISTRY, validate_parameter_combination)

CANDIDATE_CONFIG = (REPO / 'PlanetProfile' / 'Inference' / 'configs'
                    / 'enceladus_cassini_isostasy_7D.json')

pytestmark = pytest.mark.skipif(not CANDIDATE_CONFIG.exists(),
                                reason='candidate config not present')


def _meta():
    with open(CANDIDATE_CONFIG) as fh:
        return json.load(fh)['metadata']


# ---------------------------------------------------------------------------
# D7
# ---------------------------------------------------------------------------

def test_registry_rho_rock_kgm3_roundtrip():
    p = PARAMETER_REGISTRY['rho_rock_kgm3']
    assert p.category == 'structure'
    assert p.default_prior == 'uniform'
    assert p.default_bounds == [2200.0, 2600.0]
    assert p.units == 'kg/m^3'
    # zb is DERIVED from it by the analytic bijection, so no per-sample
    # structure rebuild is triggered (the frozen dispatch interpolates).
    assert p.requires_structure_rebuild is False


def test_registry_matches_the_configs_declared_frozen_prior():
    frozen = _meta()['branch_model']['frozen_branch']
    assert frozen['rho_rock_kgm3_prior'] == 'uniform [2200, 2600]'
    assert (PARAMETER_REGISTRY['rho_rock_kgm3'].default_bounds
            == [2200.0, 2600.0])


def test_every_frozen_branch_sampled_parameter_is_registered():
    """The defect in one line: the config declares the frozen branch's
    parameter set, and one of them was not a registered parameter."""
    for name in _meta()['branch_model']['frozen_branch']['sampled']:
        assert name in PARAMETER_REGISTRY, name


def test_frozen_branch_combination_validates():
    result = validate_parameter_combination(
        ['rho_rock_kgm3', 'rho_ice_kgm3', 'libration_sys_frac'])
    assert result['valid'] is True
    # None of the three needs a per-sample structure rebuild: the frozen
    # cache is indexed through the bijection, not rebuilt.
    assert result['requires_rebuild'] is False


# ---------------------------------------------------------------------------
# Freeze-gate double-count correction
# ---------------------------------------------------------------------------

def test_tajeddine_shape_term_is_carried_exactly_once():
    """It lives in Sigma_model (both the isostasy block and the generic
    likelihood-inflation hook), and the freeze gate must NOT ask for it
    again in the nuisance width."""
    meta = _meta()
    assert meta['isostasy']['sigma_model_libration_deg'] == 0.00025
    assert meta['sigma_model_add']['libration_deg'] == 0.00025
    gate = meta['libration_sys_frac_freeze_gate']
    assert 'STRUCK' in gate and 'DOUBLE-COUNT' in gate


def test_the_struck_entry_is_the_same_number_in_both_currencies():
    """0.00025 deg against sigma_obs = 0.003 deg IS the ~0.08 sigma the
    struck freeze-gate entry named -- which is why it is a double count and
    not two different terms."""
    meta = _meta()
    sigma_model = meta['isostasy']['sigma_model_libration_deg']
    sigma_obs = meta['observables']['libration_deg'][1] if 'observables' in meta \
        else 0.003
    assert sigma_model / sigma_obs == pytest.approx(0.0833, abs=5e-4)


def test_freeze_gate_still_names_the_two_live_entries():
    gate = _meta()['libration_sys_frac_freeze_gate']
    assert 'C2 coupling' in gate
    assert 'Bsp_Asp' in gate


# ---------------------------------------------------------------------------
# D8 -- budget entry only, no physics change
# ---------------------------------------------------------------------------

def test_systematics_budget_exists_and_carries_the_H22_term():
    budget = _meta()['systematics_budget']
    entry = budget['H22_first_order_projection']
    assert '0.65 sigma' in entry
    assert '5.6x' in entry
    assert 'NOT applied' in entry


def test_budget_enumerates_every_applied_and_marginalized_term():
    budget = _meta()['systematics_budget']
    for key in ('libration_shape_model', 'libration_formulation_bias_B2',
                'libration_formulation_residual_B2', 'gravity_shape_model',
                'rho_ice_EOS_nuisance', 'rho_ocean_EOS',
                'Bsp_Asp_linearization'):
        assert key in budget, key


def test_D8_changed_no_physics():
    """The H22 term is RECORDED, not applied: Sigma_model's C22 row must be
    untouched at the published H&M Table 1 value."""
    meta = _meta()
    assert meta['isostasy']['sigma_model_lm']['2,2'] == 1.7e-06
    assert meta['sigma_model_add']['C22'] == 1.7e-06


# ---------------------------------------------------------------------------
# Documentation fixes r5 listed
# ---------------------------------------------------------------------------

def test_rescinded_25_80_km_is_never_quoted_as_an_answer():
    """libration_model_correction.effect forbids quoting the rescinded
    option-A numbers, and the grid segment rationale was quoting one as the
    campaign's shifted best fit. The number may still APPEAR -- this repo
    records retractions rather than deleting them -- but only inside an
    explicit retraction or prohibition."""
    meta = _meta()
    seg = [s for s in meta['structure_cache_spec']['zb_km_grid']['segments']
           if s.get('lo') == 22.0][0]['rationale']
    assert 'RESCINDED option-A number' in seg
    assert '25.99' in seg and '27.34' in seg
    # The retracted claim must not stand as a sentence of its own.
    assert "convention) live here" not in seg

    def _strings(node):
        if isinstance(node, str):
            yield node
        elif isinstance(node, dict):
            for v in node.values():
                yield from _strings(v)
        elif isinstance(node, list):
            for v in node:
                yield from _strings(v)

    markers = ('RESCIND', 'MUST NOT', 'CORRECTED', 'DO NOT QUOTE', 'SUSPECT')
    seen = 0
    for s in _strings(meta):
        if '25.80' in s or 'to 25.8 under' in s:
            seen += 1
            assert any(m in s for m in markers), s[:200]
    assert seen >= 1, 'the retraction itself went missing'
    assert 'MUST NOT BE QUOTED' in meta['libration_model_correction']['effect']


def test_status_no_longer_claims_the_code_does_not_exist():
    meta = _meta()
    assert 'STALE CLAIM CORRECTED' in meta['status']
    assert 'code_gaps_status_2026_08_16' in meta['blockers_open']


def test_schema_version_records_that_r5_ran():
    assert 'NOT RATIFIED' in _meta()['schema_version']
