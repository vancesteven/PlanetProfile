"""D1 -- the ocean-branch MoI acceptance window is DECLARED, READ, ENFORCED.

r5 build blocker D1 (validation_reports/enceladus_isostasy/r5_ADJUDICATION.md):
the ocean-branch MoI window was "unspecified/unenforceable -- bulk_overrides is
not a config field; builder defaults to template Cuncertainty 0.001". These
tests pin the repair from both ends: the config declares the window, and
``build_zbw_grid_cache`` actually CONSUMES what the config declares.

No real PlanetProfile build happens here -- ``build_single_structure`` is
monkeypatched so the tests assert on the ``bulk_overrides`` the builder passes
DOWN, which is the thing that was broken. The measured derivation behind the
0.08 value (structural nesting argument + the 3.25e-3 achieved-deviation floor
over 8 real ocean nodes) is recorded in the config's own
``metadata.bulk_overrides_provenance``.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from PlanetProfile.Inference import cache_builder  # noqa: E402
from PlanetProfile.Inference.cache_builder import (  # noqa: E402
    bulk_overrides_from_config, build_zbw_grid_cache)

CANDIDATE_CONFIG = (REPO / 'PlanetProfile' / 'Inference' / 'configs'
                    / 'enceladus_cassini_isostasy_7D.json')

# The template value the builder used to inherit silently. Kept here so this
# test fails loudly if PPEnceladus.py's own window ever changes underneath the
# config's declaration.
TEMPLATE_CUNCERTAINTY = 0.001


def _cfg_dict():
    with open(CANDIDATE_CONFIG) as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# The config end: the window is declared.
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not CANDIDATE_CONFIG.exists(),
                    reason='candidate config not present')
def test_config_declares_the_ocean_branch_moi_window():
    meta = _cfg_dict()['metadata']
    assert 'bulk_overrides' in meta, (
        'D1: the ocean-branch MoI window must be DECLARED in the config, not '
        'inherited from the body template')
    assert meta['bulk_overrides']['Cuncertainty'] == 0.08
    # The number is worthless without its derivation on the record.
    prov = meta['bulk_overrides_provenance']
    assert 'argmin' in prov and '3.2498e-3' in prov


@pytest.mark.skipif(not CANDIDATE_CONFIG.exists(),
                    reason='candidate config not present')
def test_declared_window_is_wider_than_the_template_default():
    """The whole point: the template's 0.001 is measurably conditioning
    (it deletes the zb=5 km corner, achieved deviation 3.2498e-3)."""
    declared = _cfg_dict()['metadata']['bulk_overrides']['Cuncertainty']
    assert declared > TEMPLATE_CUNCERTAINTY
    assert declared > 3.2498e-3, (
        'declared window is below the measured achieved-deviation floor -- '
        'it would silently delete ocean nodes')


# ---------------------------------------------------------------------------
# The resolver.
# ---------------------------------------------------------------------------

def test_resolver_expands_cuncertainty_to_upper_and_lower():
    """PP tests CuncertaintyUpper/Lower, not Cuncertainty (SetupInit.py
    backfills them ONLY when the template left them None), so declaring
    Cuncertainty alone must not be able to become a silent no-op."""
    out = bulk_overrides_from_config({'metadata': {
        'bulk_overrides': {'Cuncertainty': 0.08}}})
    assert out == {'Cuncertainty': 0.08, 'CuncertaintyUpper': 0.08,
                   'CuncertaintyLower': 0.08}


def test_resolver_respects_an_explicit_asymmetric_window():
    out = bulk_overrides_from_config({'metadata': {'bulk_overrides': {
        'Cuncertainty': 0.08, 'CuncertaintyUpper': 0.05}}})
    assert out['CuncertaintyUpper'] == 0.05
    assert out['CuncertaintyLower'] == 0.08


def test_resolver_accepts_a_config_object_and_a_path():
    class _Cfg:
        metadata = {'bulk_overrides': {'Cuncertainty': 0.02}}

    assert bulk_overrides_from_config(_Cfg())['CuncertaintyUpper'] == 0.02
    if CANDIDATE_CONFIG.exists():
        assert bulk_overrides_from_config(
            str(CANDIDATE_CONFIG))['Cuncertainty'] == 0.08


@pytest.mark.parametrize('meta', [
    {},
    {'bulk_overrides': {}},
    {'bulk_overrides': {'Tb_K': 271.0}},
])
def test_resolver_refuses_to_default_silently(meta):
    """Falling through to the template default IS the D1 defect, so it is
    refused rather than defaulted."""
    with pytest.raises(ValueError, match='D1|bulk_overrides'):
        bulk_overrides_from_config({'metadata': meta})


# ---------------------------------------------------------------------------
# The builder end: the declared window actually reaches PlanetProfile.
# ---------------------------------------------------------------------------

def _fake_struct(zb_km, cmr2=0.3345):
    """Minimal ocean-branch structure dict satisfying the builder's checks."""
    r_m = np.array([1.0e5, 2.0e5, 2.521e5])
    return {
        'r_m': r_m,
        'rho': np.array([2350.0, 1005.0, 925.0]),
        'phases': np.array([50, 0, 1]),
        'D_iceIh_km': float(zb_km),
        'D_ocean_km': 30.0,
        'CMR2': float(cmr2),
        'Tb_K': 272.3,
        'Mtot_kg': 1.08022e20,
        'R_body_m': 2.521e5,
    }


@pytest.fixture
def captured_builds(monkeypatch, tmp_path):
    calls = []

    def _fake_build_single_structure(module, Tb_K, **kwargs):
        calls.append(kwargs)
        zb = kwargs['bulk_overrides']['zb_approximate_km']
        return _fake_struct(zb)

    monkeypatch.setattr(cache_builder, 'build_single_structure',
                        _fake_build_single_structure)
    monkeypatch.setattr(cache_builder, '_save_cache_atomic',
                        lambda cache, path: None)
    return calls


def test_builder_consumes_the_config_declared_window(captured_builds, tmp_path):
    """THE D1 PIN: the half-width PlanetProfile is run with comes from the
    config, not from the body template."""
    cache = build_zbw_grid_cache(
        'PlanetProfile.Default.Enceladus.PPEnceladus',
        zb_km_grid=[20.0, 25.0],
        wOcean_ppt_grid=[5.0, 20.0],
        output_path=str(tmp_path / 'x.pkl'),
        config=str(CANDIDATE_CONFIG),
        zb_tol_km=5.0,
        progress=False,
    )
    assert captured_builds, 'builder never ran a node'
    for kw in captured_builds:
        bo = kw['bulk_overrides']
        assert bo['Cuncertainty'] == 0.08
        assert bo['CuncertaintyUpper'] == 0.08
        assert bo['CuncertaintyLower'] == 0.08
        assert bo['Cuncertainty'] != TEMPLATE_CUNCERTAINTY
    win = cache['ocean_moi_window']
    assert win['CuncertaintyUpper'] == 0.08
    assert win['source'] == "config metadata['bulk_overrides']"
    # The floor the build measured for itself, and its margin.
    assert win['max_abs_cmr2_deviation'] == pytest.approx(
        abs(0.3345 - 0.335), abs=1e-12)
    assert win['margin'] == pytest.approx(0.08 - abs(0.3345 - 0.335), abs=1e-12)


def test_explicit_bulk_overrides_win_over_the_config(captured_builds, tmp_path):
    """A probe must still be able to sweep the window deliberately -- but the
    cache records that it did, so a probe cache cannot be mistaken for a
    config-conformant one."""
    cache = build_zbw_grid_cache(
        'PlanetProfile.Default.Enceladus.PPEnceladus',
        zb_km_grid=[20.0, 25.0],
        wOcean_ppt_grid=[5.0, 20.0],
        output_path=str(tmp_path / 'x.pkl'),
        config=str(CANDIDATE_CONFIG),
        bulk_overrides={'Cuncertainty': 0.30, 'CuncertaintyUpper': 0.30,
                        'CuncertaintyLower': 0.30},
        zb_tol_km=5.0,
        progress=False,
    )
    for kw in captured_builds:
        assert kw['bulk_overrides']['Cuncertainty'] == 0.30
    assert cache['ocean_moi_window']['source'].startswith('explicit argument')


def test_builder_without_a_config_is_unchanged(captured_builds, tmp_path):
    """Backwards compatibility: the config= argument is opt-in."""
    build_zbw_grid_cache(
        'PlanetProfile.Default.Enceladus.PPEnceladus',
        zb_km_grid=[20.0, 25.0],
        wOcean_ppt_grid=[5.0, 20.0],
        output_path=str(tmp_path / 'x.pkl'),
        bulk_overrides={'Cuncertainty': 0.08},
        zb_tol_km=5.0,
        progress=False,
    )
    for kw in captured_builds:
        assert kw['bulk_overrides']['Cuncertainty'] == 0.08
        assert 'CuncertaintyUpper' not in kw['bulk_overrides']


def test_builder_refuses_a_config_without_a_window(captured_builds, tmp_path):
    with pytest.raises(ValueError, match='D1|bulk_overrides'):
        build_zbw_grid_cache(
            'PlanetProfile.Default.Enceladus.PPEnceladus',
            zb_km_grid=[20.0, 25.0],
            wOcean_ppt_grid=[5.0, 20.0],
            output_path=str(tmp_path / 'x.pkl'),
            config={'metadata': {}},
            progress=False,
        )
