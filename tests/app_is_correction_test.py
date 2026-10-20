"""Phase 2 GUI wiring of the precomputed NH3 fiducial IS correction
(handoff section 0.32, RATIFICATION 2026-10-20, conditions B5-B9).

Unit tests exercise the weights helpers directly.  AppTest tests drive the
real Inference page (amortized mode, Titan NH3 slot).  The feature flag is
OFF by default; the tests flip ``is_fiducial_bundle.FEATURE_ENABLED`` with
monkeypatch (the GUI reads it at call time).

Status of what this proves: AppTest executes the script headlessly and reads
the element tree; it is NOT visual verification in a browser.
"""
import dataclasses
import json
import sys
import types
from pathlib import Path

import matplotlib
matplotlib.use('Agg')  # AppTest runs the page in a worker thread
import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'PlanetProfileApp'))

from PlanetProfile.Inference import is_fiducial_bundle as fb  # noqa: E402
from Utilities import is_correction_ui as ui  # noqa: E402

PAGE = str(REPO / 'PlanetProfileApp/pages/Inference.py')
MANIFEST = (REPO / 'PlanetProfile/Inference/sbi_artifacts/'
            'titan_freegrav_nh3_is_fiducial.json')

REQUIRED = [
    REPO / 'PlanetProfile/Inference/sbi_artifacts/'
    'titan_freegrav_nh3_posterior_1m.pt',
    REPO / 'PlanetProfile/Test/mcmc_results/Titan/Test54_nh3_ocean/'
    'titan_nh3_joint_structure_grid_2d.pkl',
]
needs_assets = pytest.mark.skipif(
    not all(p.exists() for p in REQUIRED),
    reason='NH3 artifact / structure cache not present')


# ---------------------------------------------------------------------------
# Unit tests: weights helpers
# ---------------------------------------------------------------------------

def _fake_result(n=2000, corrected=True, seed=1):
    from PlanetProfile.Inference.inference_core import InferenceResult
    rng = np.random.default_rng(seed)
    S = rng.normal(size=(n, 3))
    w = np.exp(S[:, 0])           # weights correlated with column 0
    w /= w.sum()
    k2 = np.column_stack([0.6 + 0.01 * S[:, 1], 0.05 + 0.01 * S[:, 2]])
    md = {'sampler': 'sbi', 'heating_indices': np.arange(5),
          'has_ocean': S[:, 1] > 0}
    if corrected:
        md.update({'is_corrected': True, 'is_weights': w,
                   'is_diagnostics': {'ess': 1.0 / np.sum(w ** 2)}})
    return InferenceResult(
        config=types.SimpleNamespace(observables={}), samples=S,
        log_likelihoods=np.zeros(n), param_names=['a', 'b', 'c'],
        param_labels=['a', 'b', 'c'], k2_results=k2,
        D_ocean_results=S[:, 0] * 10.0,
        D_icePhase_results={'III': S[:, 0] + 1.0},
        heating_results=[{'Sil': 1.0}] * 5, metadata=md,
        weights=w if corrected else None), w


def test_result_weights_none_when_uncorrected():
    r, _ = _fake_result(corrected=False)
    assert ui.result_weights(r) is None
    assert ui.display_view(r) is r           # byte-identical path


def test_result_weights_normalized_and_validated():
    r, w = _fake_result()
    got = ui.result_weights(r)
    assert got is not None and np.isclose(got.sum(), 1.0)
    np.testing.assert_allclose(got, w)
    # Wrong-length weights on a result flagged corrected must FAIL LOUDLY: a
    # silent None would render every panel from raw flow draws, unlabelled.
    r.weights = w[:-1]
    r.metadata['is_weights'] = w[:-1]
    with pytest.raises(ValueError, match='IS weights are invalid'):
        ui.result_weights(r)
    r.weights = None
    r.metadata.pop('is_weights')
    with pytest.raises(ValueError, match='carries no IS weights'):
        ui.result_weights(r)


def test_display_view_is_equal_weight_and_not_flagged_corrected():
    r, _ = _fake_result()
    view = ui.display_view(r)
    # The equal-weight view must not re-enter the corrected branch (and must
    # not trip the fail-loudly guard) when helpers are called on it.
    assert ui.result_weights(view) is None
    assert view.metadata['is_display_view_of_corrected'] is True
    assert ui.display_view(view) is view


def test_weighted_summary_matches_numpy_reference():
    rng = np.random.default_rng(0)
    v = rng.normal(size=5000)
    w = np.exp(v)
    w /= w.sum()
    mean, med, std, p16, p84 = ui.weighted_summary(v, w)
    assert mean == pytest.approx(np.sum(w * v))
    assert std == pytest.approx(np.sqrt(np.sum(w * (v - mean) ** 2)))
    # tilted by exp(v): weighted mean ~ +1 vs unweighted ~ 0
    assert mean > 0.8
    assert med > np.median(v) + 0.5
    um, umed, _, _, _ = ui.weighted_summary(v, None)
    assert um == pytest.approx(np.mean(v))


def test_display_view_resamples_consistently_and_reproducibly():
    r, w = _fake_result()
    d1 = ui.display_view(r, seed=7)
    d2 = ui.display_view(r, seed=7)
    d3 = ui.display_view(r, seed=8)
    np.testing.assert_array_equal(d1.samples, d2.samples)
    assert not np.array_equal(d1.samples, d3.samples)
    assert d1.weights is None and ui.result_weights(d1) is None
    # pairing preserved across arrays (same indices everywhere)
    np.testing.assert_allclose(d1.D_ocean_results, d1.samples[:, 0] * 10.0)
    np.testing.assert_allclose(d1.D_icePhase_results['III'],
                               d1.samples[:, 0] + 1.0)
    np.testing.assert_array_equal(d1.metadata['has_ocean'],
                                  d1.samples[:, 1] > 0)
    # heating is dropped (uncorrected subset must not be shown)
    assert d1.heating_results is None
    assert d1.metadata['heating_indices'] is None
    assert d1.metadata['is_display_resample']['seed'] == 7
    # resample reproduces the weighted mean of column 0 (not the raw mean)
    wm = float(np.sum(w * r.samples[:, 0]))
    assert abs(d1.samples[:, 0].mean() - wm) < 0.1
    assert abs(r.samples[:, 0].mean() - wm) > 0.3


def test_weighted_ocean_fraction():
    r, w = _fake_result()
    expect = float(np.sum(w * (r.samples[:, 1] > 0)))
    assert ui.weighted_ocean_fraction(r) == pytest.approx(expect)
    r2, _ = _fake_result(corrected=False)
    assert ui.weighted_ocean_fraction(r2) is None


def test_weighted_csv_has_weight_column_and_manifest_sha():
    r, w = _fake_result(n=50)
    r.metadata['is_manifest_sha256'] = 'abc123'
    r.metadata['is_npz_sha256'] = 'def456'
    csv = ui.weighted_csv(r).splitlines()
    assert 'abc123' in csv[0] and 'def456' in csv[0]
    assert csv[1] == 'a,b,c,weight'
    assert len(csv) == 2 + 50
    assert float(csv[2].split(',')[-1]) == pytest.approx(w[0], rel=1e-6)
    r2, _ = _fake_result(corrected=False)
    assert ui.weighted_csv(r2) is None


def test_warning_text_numbers_come_from_manifest(monkeypatch):
    # B7: changing manifest/bundle numbers changes the text; no literals.
    m = json.loads(MANIFEST.read_text())
    gv = m['gui_gate_values']
    assert ui._ceil2(gv['pushforward_median_gap_sigma_obs']) == 0.15
    assert ui._ceil2(gv['pushforward_ks']) == 0.07
    tmpl = ('{params}|{N}|{ess}|{pareto_k}|{gap}|{ks}|{n_seeds}|{flow_im}|'
            '{corr_im}|{obs_im}|{obs_im_sig}')
    k2 = np.column_stack([np.full(10, 0.6), np.full(10, -0.1)])
    bundle = types.SimpleNamespace(
        manifest={**m, 'n_samples': 12345,
                  'gui_gate_values': {**gv, 'n_seeds': 9}},
        diagnostics={'ess': 777.4, 'pareto_k': 0.31},
        derived={'k2_results': k2}, weights=np.full(10, 0.1))
    txt = ui.corrected_warning(tmpl, bundle, {}, {'Im_k2': 0.135},
                               {'Im_k2': 0.035})
    assert txt == ('Re k₂, |Im k₂||12,345|777|0.31|0.15|0.07|9|0.100|0.10|'
                   '0.135|0.035')
    # g1_passing_params, when present, extends the list (currently absent).
    assert ui._params_text(m, {}) == 'Re k₂, |Im k₂|'
    assert ui._params_text(m, {'g1_passing_params': ['ζ']}) == \
        'Re k₂, |Im k₂|, ζ'


def test_flag_default_off_and_inactive_without_slot_entry(monkeypatch):
    assert fb.FEATURE_ENABLED is False
    assert not ui.correction_active_for({'is_correction': {'manifest': 'x'}})
    monkeypatch.setattr(fb, 'FEATURE_ENABLED', True)
    assert ui.correction_active_for({'is_correction': {'manifest': 'x'}})
    assert not ui.correction_active_for({'sector_warning': 'x'})
    assert ui.attempt_correction({'sector_warning': 'x'}, REPO, None) is None


# ---------------------------------------------------------------------------
# AppTest: end to end through the real page
# ---------------------------------------------------------------------------

def _boot():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(PAGE, default_timeout=900)
    at.session_state['inference_exec_mode'] = 'amortized'
    at.run()
    assert not at.exception, str([e.value for e in at.exception][:2])

    def _sel(pfx):
        return next(s for s in at.selectbox
                    if str(getattr(s, 'key', '')).startswith(pfx))
    an = _sel('amort_analysis')
    an.select(next(o for o in an.options if 'Titan' in o))
    at.run()
    ver = _sel('amort_version_')
    ver.select(next(o for o in ver.options
                    if 'NH₃' in str(o) and 'awaiting' not in str(o)))
    at.run()
    assert not at.exception, str([e.value for e in at.exception][:2])
    return at


def _set_slider(at, key, value):
    next(s for s in at.select_slider if s.key == key).set_value(value)


def _click_run(at):
    btn = next(b for b in at.button if 'Generate Posterior' in b.label)
    btn.click()
    at.run()
    assert not at.exception, str([e.value for e in at.exception][:2])


def _warnings(at):
    return ' '.join(str(w.value) for w in at.warning)


def _all_text(at):
    parts = [str(w.value) for w in at.warning]
    parts += [str(c.value) for c in at.caption]
    parts += [str(c.value) for c in at.info]
    parts += [str(m.value) for m in at.metric]
    parts += [str(m.label) for m in at.metric]
    parts += [str(m.value) for m in at.markdown]
    for t in at.table:
        parts.append(t.value.to_csv())
    return '\n'.join(parts)


@needs_assets
def test_flag_off_warning_unchanged():
    at = _boot()
    _set_slider(at, 'amort_n_post', 1000)
    _set_slider(at, 'amort_n_reeval', 100)
    _set_slider(at, 'amort_n_derived', 500)
    _click_run(at)
    w = _warnings(at)
    assert 'Split ratification' in w and 'under-updates k₂' in w
    assert 'importance-sampling' not in w
    assert 'correction unavailable' not in w
    res = at.session_state['inference_results']
    assert 'is_corrected' not in res.metadata
    assert res.weights is None


@needs_assets
def test_flag_on_default_datum_shows_corrected_weighted_output(monkeypatch):
    monkeypatch.setattr(fb, 'FEATURE_ENABLED', True)
    at = _boot()
    _set_slider(at, 'amort_n_post', 20000)
    _set_slider(at, 'amort_n_derived', 'All')
    _click_run(at)
    w = _warnings(at)
    assert 'Tidal sector at the DEFAULT Cassini datum only' in w, w[:400]
    # G1 passed for all seven rheology params (r1_r2_RESULT.md, 2026-10-20).
    assert ('Re k₂, |Im k₂|, ζ, η (ice Ih, III, V, VI, silicate), α '
            'are shown') in w
    assert 'effective sample size ≈ 1,271' in w
    assert 'Pareto-k 0.19' in w
    assert 'within 0.15 σ_obs' in w and 'KS ≤ 0.07' in w
    assert 'Split ratification' not in w
    res = at.session_state['inference_results']
    assert res.metadata['is_corrected'] is True
    assert res.metadata['is_manifest_sha256']
    assert np.isclose(res.weights.sum(), 1.0)
    # B6: ESS shown, not N
    labels = {m.label: m.value for m in at.metric}
    assert labels['Effective sample size'] == '1,271'
    assert 'Posterior draws' not in labels
    # weighted ocean fraction ~0.929 (not the flow's 0.933)
    assert labels['Ocean fraction (weighted)'] == '0.929'
    text = _all_text(at)
    # weighted |Im k2| median ~0.109 in the pushforward table; the
    # uncorrected flow median (0.045) appears only in the flow column.
    rows = next(t.value for t in at.table
                if 'Corrected median' in t.value.columns)
    im = rows[rows['Quantity'] == '|Im k₂|'].iloc[0]
    assert im['Corrected median'] == '0.109'
    assert im['Uncorrected flow median'] == '0.045'
    # heating explicitly disabled, never silently unweighted
    assert 'deliberately not' in text
    # summary table is flagged as weighted
    assert 'Importance-weighted statistics' in text


@needs_assets
def test_flag_on_modified_datum_falls_back(monkeypatch):
    monkeypatch.setattr(fb, 'FEATURE_ENABLED', True)
    # The loader sees the full-N request (20000 / All) so the ONLY reason to
    # fall back is the modified datum; the uncorrected fallback run itself is
    # shrunk (a full 20000-draw recompute takes minutes) -- the test targets
    # the warning/fallback wiring, not the flow recompute.
    from PlanetProfile.Inference.sbi_runner import SBIRunner
    orig = SBIRunner.infer_from_artifact

    def _small(self, artifact_path, **kw):
        kw.update(n_posterior_samples=500, n_derived=200, n_reeval=50)
        return orig(self, artifact_path, **kw)
    monkeypatch.setattr(SBIRunner, 'infer_from_artifact', _small)
    at = _boot()
    _set_slider(at, 'amort_n_post', 20000)
    _set_slider(at, 'amort_n_derived', 'All')
    reob = next(n for n in at.number_input
                if str(n.key).startswith('amort_obs_')
                and 'Re_k2' in str(n.key))
    reob.set_value(float(reob.value) + 0.01)
    at.run()
    _click_run(at)
    w = _warnings(at)
    assert 'Tidal-sector correction unavailable for this conditioning' in w
    assert 'modified datum' in w
    assert 'about 15%' in w
    assert 'Tidal sector at the DEFAULT' not in w
    res = at.session_state['inference_results']
    assert res.metadata['is_corrected'] is False
    assert 'modified datum' in res.metadata['is_fallback_reason']
    assert res.weights is None


@needs_assets
def test_flag_on_truncated_derived_recompute_falls_back(monkeypatch):
    monkeypatch.setattr(fb, 'FEATURE_ENABLED', True)
    at = _boot()
    _set_slider(at, 'amort_n_post', 20000)
    _set_slider(at, 'amort_n_reeval', 100)
    _set_slider(at, 'amort_n_derived', 1000)
    _click_run(at)
    w = _warnings(at)
    assert 'Tidal-sector correction unavailable' in w
    assert 'full-N likelihoods absent' in w
    res = at.session_state['inference_results']
    assert res.metadata['is_corrected'] is False
    assert res.weights is None


@needs_assets
def test_scope_note_sentence_replaced_only_on_corrected_path(monkeypatch):
    def caps(at):
        return ' '.join(str(c.value) for c in at.caption)
    at = _boot()
    _set_slider(at, 'amort_n_post', 20000)
    _set_slider(at, 'amort_n_derived', 'All')
    at.run()
    assert 'DO NOT quote SBI Re k₂' in caps(at)          # flag off
    monkeypatch.setattr(fb, 'FEATURE_ENABLED', True)
    at.run()
    c = caps(at)
    assert 'DO NOT quote SBI Re k₂' not in c
    assert 'precomputed importance-sampling correction' in c
    _set_slider(at, 'amort_n_derived', 1000)             # fallback config
    at.run()
    assert 'DO NOT quote SBI Re k₂' in caps(at)


@needs_assets
def test_corrected_path_panels_render_without_errors(monkeypatch):
    monkeypatch.setattr(fb, 'FEATURE_ENABLED', True)
    at = _boot()
    _set_slider(at, 'amort_n_post', 20000)
    _set_slider(at, 'amort_n_derived', 'All')
    _click_run(at)
    bad = [str(w.value) for w in at.warning if 'unavailable' in str(w.value)]
    assert not bad, bad
    assert not at.error, [str(e.value) for e in at.error]


def test_fallback_warning_text_hint_and_punctuation():
    tmpl = ('unavailable ({reason}). {scope}; about {fail_pct}% fail.')
    size = ui.fallback_warning(
        tmpl, 'full-N likelihoods absent: n_post=10000 differs from the '
        'validated draw count N=20000.', None)
    assert 'set Posterior samples = 20,000 and Derived subset = All to ' \
           'use the precomputed correction' in size
    assert '.).' not in size and 'modified datum' not in size
    datum = ui.fallback_warning(
        tmpl, 'modified datum: x_obs differs from the validated fiducial '
        '(values or observable set).', None)
    assert 'Posterior samples' not in datum
    assert 'validated only at the default Cassini datum' in datum
    assert '.).' not in datum
