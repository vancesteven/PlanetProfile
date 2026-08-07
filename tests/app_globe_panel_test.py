"""Regression tests for the Inference-page globe panel (promoted from the
2026-07-24 session verify scripts).

Covers:
- radial-profile figure: 6 panels incl. viscosity (log, effective) and
  shear modulus;
- PP-format table: effective eta column, numeric coercion, empty-column
  and wide-range (e-notation) detection logic;
- per-run assumptions markdown (Andrade/Arrhenius and Maxwell variants);
- AppTest against the shipped Europa v4 result: exception-clean render,
  x/y axis selectors, no-ocean axis gating, derived-subset caption,
  PlotWedge caption, no LaTeX (usetex) warnings;
- pp_wedge_exports core suppression note for coreless bodies.

Run from the repo root: pytest tests/app_globe_panel_test.py
"""
import pickle
import sys
import types
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / 'PlanetProfileApp'))

matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')
pytest.importorskip('streamlit')

EUROPA_V4_PKL = (REPO / 'PlanetProfile/Test/mcmc_results/Europa/'
                 'Test53_geodesy_v4/europa_clipper_v4_reference_result.pkl')


def _synth_prof(n=40):
    r = np.linspace(100, 1500, n)
    return {
        'r_km': r,
        'rho': np.linspace(8000, 900, n),
        'T_K': np.linspace(1200, 100, n),
        'P_MPa': np.linspace(5000, 0.1, n),
        'sigma_Sm': np.full(n, 1e-16),
        'eta_Pa_base': np.full(n, 1e14),
        'eta_Pa': np.full(n, 1e15),
        'mu_Pa': np.full(n, 3.5e9),
        'K_Pa': np.full(n, 9e9),
        'phases': np.concatenate([np.full(10, 100), np.full(10, 50),
                                  np.full(10, 0), np.full(10, 1)]),
        'VP_kms': np.full(n, 2.0), 'VS_kms': np.full(n, 1.0),
        'MLayer_kg': np.full(n, 1.2e18), 'VLayer_m3': np.full(n, 4.4e14),
    }


def test_profile_figure_panels():
    from Utilities.radial_profiles import build_profile_figure
    fig = build_profile_figure(_synth_prof(), 'test')
    axes = fig.get_axes()
    assert len(axes) == 6
    labs = [ax.get_xlabel() for ax in axes]
    assert any('eta' in l for l in labs)
    assert any('mu' in l for l in labs)
    eta_ax = axes[3]
    assert eta_ax.get_xscale() == 'log'
    # effective eta (1e15), not base (1e14)
    assert np.allclose(eta_ax.get_lines()[0].get_xdata(), 1e15)


def test_profile_table_effective_eta_and_formats():
    import pandas as pd
    from Utilities.radial_profiles import profile_table
    tbl = profile_table(_synth_prof())
    row0 = tbl[0]
    assert abs(row0['eta (Pa s)'] - 1e15) < 1
    assert abs(row0['GS (GPa)'] - 3.5) < 1e-6

    df = pd.DataFrame(tbl)
    for c in df.columns:
        if c != 'phase ID':
            df[c] = pd.to_numeric(df[c], errors='coerce')
    assert all(df[c].dtype.kind == 'f' for c in df.columns
               if c != 'phase ID')
    empty = [c for c in df.columns if c != 'phase ID'
             and not np.isfinite(df[c].astype(float)).any()]
    assert 'phi (void/solid frac)' in empty and 'QS' in empty
    shown = df.drop(columns=empty)

    def wide(col):
        v = np.abs(shown[col].to_numpy(float))
        v = v[np.isfinite(v) & (v > 0)]
        return v.size and (v.max() >= 1e5 or v.min() < 1e-3)

    wides = [c for c in shown.columns if c != 'phase ID' and wide(c)]
    assert {'eta (Pa s)', 'MLayer (kg)', 'r (m)'} <= set(wides)
    assert 'T (K)' not in wides and 'rho (kg/m3)' not in wides


def _fake_cfg(**over):
    base = dict(
        bodyname='Titan',
        param_space={'log10_eta_Ih': {'prior_type': 'uniform',
                                      'bounds': [12, 17]},
                     'alpha': {'prior_type': 'uniform',
                               'bounds': [0.1, 0.5]}},
        fixed_params={'Tb_K': 250.965},
        param_groups={'log10_eta_HP': ['log10_eta_III', 'log10_eta_V',
                                       'log10_eta_VI']},
        derived_params={},
        observables={'Re_k2': (0.589, 0.075)},
        sampler_settings={},
        arrhenius_params={'activation_energy_kJ_mol': {'Ih': 59.4},
                          'reference_temp_K': 250.965},
        structure_cache_path='foo/titan_cache.pkl',
    )
    base.update(over)
    return types.SimpleNamespace(**base)


def test_assumptions_markdown():
    from Utilities.run_assumptions import describe_assumptions
    res = types.SimpleNamespace(config=_fake_cfg(),
                                metadata={'sampler': 'sbi'})
    md = describe_assumptions(res)
    for frag in ('Andrade', 'Arrhenius', 'POROUS_ROCK = False',
                 'Tb_K` = 250.965', '59.4 kJ/mol', 'titan_cache.pkl',
                 'amortized'):
        assert frag in md, f'missing: {frag}'
    md2 = describe_assumptions(types.SimpleNamespace(
        config=_fake_cfg(arrhenius_params={},
                         param_space={'log10_eta_Ih':
                                      {'prior_type': 'uniform',
                                       'bounds': [12, 17]}}),
        metadata={}))
    assert 'No Arrhenius' in md2 and 'Maxwell' in md2
    assert 'MCMC (pocoMC preconditioned Monte Carlo)' in md2
    assert 'no direct C₂₀/C₂₂ forward model' in md2

    md3 = describe_assumptions(types.SimpleNamespace(
        config=_fake_cfg(gravity_forward_model='clairaut_hydrostatic'),
        metadata={}))
    assert 'direct Clairaut integration' in md3
    assert 'Radau–Darwin' not in md3


@pytest.mark.skipif(not EUROPA_V4_PKL.exists(),
                    reason='Europa v4 reference result not present')
def test_apptest_europa_v4_globe_panel():
    from streamlit.testing.v1 import AppTest
    with open(EUROPA_V4_PKL, 'rb') as f:
        res = pickle.load(f)
    at = AppTest.from_file(
        str(REPO / 'PlanetProfileApp/pages/Inference.py'),
        default_timeout=900)
    at.session_state['inference_results'] = res
    at.session_state['Planet'] = 'Europa'
    at.run()
    assert not at.exception, str([e.value for e in at.exception][:2])
    # No LaTeX/usetex fallout (u panel etc.)
    assert not [w for w in at.warning if 'latex' in str(w.value).lower()]

    bgs = [el for el in at.main if type(el).__name__ == 'ButtonGroup'
           and 'globe_' in str(getattr(el, 'key', ''))]
    keys = [str(w.key) for w in bgs]
    assert any('globe_xaxis' in k for k in keys)
    assert any('globe_yaxis' in k for k in keys)
    xw = next(w for w in bgs if 'globe_xaxis' in str(w.key))
    xopts = [getattr(o, 'content', o) for o in xw.options]
    assert 'Ocean thickness' in xopts  # Europa has an ocean

    caps = [str(c.value) for c in at.main if type(c).__name__ == 'Caption']
    assert any('PlanetProfile PlotWedge for' in c for c in caps)
    md_all = ' '.join(str(getattr(m, 'value', '')) for m in at.main
                      if type(m).__name__ in ('Markdown', 'Caption'))
    assert 'porosity is OFF' in md_all
    assert 'between-node blending' in md_all

    btn_keys = [str(getattr(b, 'key', ''))
                for b in at.get('download_button')]
    assert any('globe_pp_txt' in k for k in btn_keys)
    assert any('globe_pp_csv' in k for k in btn_keys)

    # Heating + geotherm tabs (2026-07-25): totals rendered, provenance
    # captions state the Qrad assumption and Perple_X cache-build-only use.
    mets = [str(getattr(m, 'label', '')) for m in at.main
            if type(m).__name__ == 'Metric']
    assert 'Total tidal power' in mets and 'Total radiogenic power' in mets
    caps_all = ' '.join(str(c.value) for c in at.main
                        if type(c).__name__ == 'Caption')
    assert 'TidalPy radial' in caps_all and 'Q_rad' in caps_all
    assert 'Perple_X at cache-build time only' in caps_all
    assert caps_all.count(
        'between-node structure is bilinearly blended') >= 4
    # wide-range table cells ship as pre-formatted e-notation strings
    dfs = [el for el in at.main if type(el).__name__ == 'Dataframe']
    tbl = next(d.value for d in dfs if 'eta (Pa s)' in d.value.columns)
    assert isinstance(tbl['eta (Pa s)'].iloc[0], str) \
        and 'e+' in tbl['eta (Pa s)'].iloc[0]
    assert isinstance(tbl['MLayer (kg)'].iloc[0], str)


@pytest.mark.skipif(not EUROPA_V4_PKL.exists(),
                    reason='Europa v4 reference result not present')
def test_apptest_heating_inventory_selector():
    # Radiogenic model selector (2026-07-28): body-default vs McDonough+
    # 2020 inventory Q_rad; switching updates the caption source label
    # and the time-rewind slider appears only for inventory modes.
    from streamlit.testing.v1 import AppTest
    with open(EUROPA_V4_PKL, 'rb') as f:
        res = pickle.load(f)
    at = AppTest.from_file(
        str(REPO / 'PlanetProfileApp/pages/Inference.py'),
        default_timeout=900)
    at.session_state['inference_results'] = res
    at.session_state['Planet'] = 'Europa'
    at.run()
    assert not at.exception, str([e.value for e in at.exception][:2])
    sel = next(s for s in at.selectbox
               if str(getattr(s, 'key', '')) == 'heat_qrad_model')
    assert sel.value.startswith('Body default')
    assert not [s for s in at.get('slider')
                if str(getattr(s, 'key', '')) == 'heat_qrad_tGa']
    caps = ' '.join(str(c.value) for c in at.main
                    if type(c).__name__ == 'Caption')
    assert 'body config constant Sil.Qrad_Wkg' in caps

    sel.select('BSE inventory (McDonough+ 2020)')
    at.run()
    assert not at.exception, str([e.value for e in at.exception][:2])
    sliders = [s for s in at.get('slider')
               if str(getattr(s, 'key', '')) == 'heat_qrad_tGa']
    assert sliders, 'time-rewind slider missing in inventory mode'
    caps = ' '.join(str(c.value) for c in at.main
                    if type(c).__name__ == 'Caption')
    assert 'BSE inventory (McDonough+ 2020)' in caps
    # present-day BSE specific heating ~4.9e-12 W/kg must be quoted
    assert '4.91e-12' in caps

    sliders[0].set_value(4.5)
    at.run()
    assert not at.exception
    caps = ' '.join(str(c.value) for c in at.main
                    if type(c).__name__ == 'Caption')
    assert 't = 4.5 Ga' in caps
    assert '2.43e-11' in caps  # rewound heating ~5x present BSE


@pytest.mark.skipif(not EUROPA_V4_PKL.exists(),
                    reason='Europa v4 reference result not present')
def test_apptest_mineralogy_tab():
    # Post-hoc Perple_X mineralogy check (2026-07-28): table of shipped
    # databases with grain rho, mismatch, porosity leverage, verdicts;
    # caption states the post-hoc/Han-2014 assumptions.
    from streamlit.testing.v1 import AppTest
    with open(EUROPA_V4_PKL, 'rb') as f:
        res = pickle.load(f)
    at = AppTest.from_file(
        str(REPO / 'PlanetProfileApp/pages/Inference.py'),
        default_timeout=900)
    at.session_state['inference_results'] = res
    at.session_state['Planet'] = 'Europa'
    at.run()
    assert not at.exception, str([e.value for e in at.exception][:2])
    dfs = [el.value for el in at.main if type(el).__name__ == 'Dataframe']
    mtbl = next((d for d in dfs if 'Verdict' in getattr(d, 'columns', [])),
                None)
    assert mtbl is not None, 'mineralogy dataframe missing'
    assert 'Perple_X table' in mtbl.columns
    assert 'Porosity headroom (%)' in mtbl.columns
    assert len(mtbl) >= 5  # all shipped silicate tables evaluated
    verdicts = ' '.join(mtbl['Verdict'])
    # Europa v4 posterior median under the ±2% band (review 2026-07-28):
    # CI/CM consistent; CV3/CV too dense; 67P too light
    assert 'consistent within' in verdicts
    assert 'too dense by' in verdicts and 'too light by' in verdicts
    # deep mantle: phi_top ill-conditioned (leverage ~0.5% << 8%) —
    # must NOT be quoted in any verdict
    assert 'phi_top' not in verdicts
    # every shipped table must resolve to a physical grain density
    assert all(float(v) > 2000 for v in mtbl['Grain ρ (kg/m³)'])
    caps = ' '.join(str(c.value) for c in at.main
                    if type(c).__name__ == 'Caption')
    assert 'POST-HOC consistency check' in caps
    assert 'Han (2014)' in caps
    assert 'VACUUM pores' in caps
    mets = [str(getattr(m, 'label', '')) for m in at.main
            if type(m).__name__ == 'Metric']
    assert 'Inferred bulk silicate density' in mets
    expanders = [str(getattr(e, 'label', '')) for e in at.get('expander')]
    assert 'Geotherm vs Perple_X table' in expanders
    composition = next(
        s for s in at.selectbox
        if str(getattr(s, 'key', '')) == 'mineralogy_perplex_composition')
    assert composition.value.endswith('.tab')
    download_keys = [str(getattr(b, 'key', ''))
                     for b in at.get('download_button')]
    assert 'mineralogy_perplex_geotherm_pdf' in download_keys
    assert 'mineralogy_perplex_geotherm_png' in download_keys

    # The overlay helper must preserve the exact selected-draw silicate
    # P,T ranges used by the mineralogy summary (posterior median here).
    from Utilities.radial_profiles import (
        build_perplex_geotherm_figure, mineralogy_for_sample,
        perplex_geotherm_overlay_data)
    summary, note = mineralogy_for_sample(res, None, REPO)
    assert summary is not None, note
    overlay, note = perplex_geotherm_overlay_data(
        res, None, REPO, composition.value)
    assert overlay is not None, note
    assert np.allclose(overlay['P_path_range_MPa'],
                       summary['P_sil_range_MPa'])
    assert np.allclose(overlay['T_path_range_K'],
                       summary['T_sil_range_K'])
    assert overlay['rho_grid_kgm3'].shape == (
        len(overlay['P_axis_MPa']), len(overlay['T_axis_K']))
    assert overlay['P_native_range_MPa'] == (
        float(overlay['P_axis_MPa'][0]),
        float(overlay['P_axis_MPa'][-1]))
    assert overlay['T_native_range_K'] == (
        float(overlay['T_axis_K'][0]),
        float(overlay['T_axis_K'][-1]))
    fig = build_perplex_geotherm_figure(overlay, 'test overlay')
    assert len(fig.axes) == 2  # heatmap axis + density colorbar
    assert fig.axes[0].get_lines(), 'selected-draw P,T path missing'
    band = fig.axes[0].patches[0]
    assert band.get_hatch(), 'cold-edge band hatch missing'
    assert band.get_linewidth() > 0, 'cold-edge band outline missing'
    assert band.get_facecolor()[:3] != (1.0, 1.0, 1.0)
    matplotlib.pyplot.close(fig)


@pytest.mark.skipif(not EUROPA_V4_PKL.exists(),
                    reason='Europa v4 reference result not present')
def test_apptest_noocean_axis_gating_and_subset_note():
    from streamlit.testing.v1 import AppTest
    with open(EUROPA_V4_PKL, 'rb') as f:
        res = pickle.load(f)
    n = len(res.samples)
    res.D_ocean_results = np.zeros(n)      # no ocean anywhere
    d_cl = np.full(n, np.nan)              # clathrate: subset, 3 grid nodes
    idx = np.arange(0, n, max(1, n // 200))[:200]
    d_cl[idx] = np.random.default_rng(1).choice([1.5, 2.5, 4.0], idx.size)
    res.D_clath_results = d_cl

    at = AppTest.from_file(
        str(REPO / 'PlanetProfileApp/pages/Inference.py'),
        default_timeout=900)
    at.session_state['inference_results'] = res
    at.session_state['Planet'] = 'Europa'
    at.run()
    assert not at.exception, str([e.value for e in at.exception][:2])
    bgs = [el for el in at.main if type(el).__name__ == 'ButtonGroup'
           and 'globe_' in str(getattr(el, 'key', ''))]
    xw = next(w for w in bgs if 'globe_xaxis' in str(w.key))
    yw = next(w for w in bgs if 'globe_yaxis' in str(w.key))
    xopts = [getattr(o, 'content', o) for o in xw.options]
    yopts = [getattr(o, 'content', o) for o in yw.options]
    assert 'Ocean thickness' not in xopts + yopts
    assert 'Clathrate thickness' in xopts

    xw.set_value('Clathrate thickness').run()
    assert not at.exception
    caps = [str(c.value) for c in at.main if type(c).__name__ == 'Caption']
    notes = [c for c in caps if c.startswith('Note:')]
    assert notes and 'distinct nodes' in notes[0] \
        and 'posterior draws carry' in notes[0]


def test_apptest_titan_slot_readiness_fallback():
    # A just-delivered slot whose structure cache is not present on this
    # machine must never be the default model version (2026-08-03: the NH3
    # joint slot landed cache-less and became the erroring Titan default).
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(
        str(REPO / 'PlanetProfileApp/pages/Inference.py'),
        default_timeout=900)
    at.run()
    radio = next(r for r in at.radio
                 if 'method' in str(getattr(r, 'label', '')).lower())
    radio.set_value(next(o for o in radio.options
                         if 'mort' in str(o) or 'pretrained' in str(o).lower()))
    at.run()

    def _fresh(pfx):
        return next(s for s in at.selectbox
                    if str(getattr(s, 'key', '')).startswith(pfx))
    _fresh('amort_analysis').select(
        next(o for o in _fresh('amort_analysis').options if 'Titan' in o))
    at.run()
    ver = _fresh('amort_version_')
    nh3_cache = (REPO / 'PlanetProfile/Test/mcmc_results/Titan/'
                 'Test54_nh3_ocean/titan_nh3_joint_structure_grid_2d.pkl')
    nh3_opt = next(o for o in ver.options if 'NH' in str(o))
    if nh3_cache.exists():
        # Cache present: newest ready slot (NH3) is the legitimate default.
        assert str(ver.value) == str(nh3_opt)
    else:
        # Cache absent: default falls back to the newest READY slot; the
        # NH3 version stays selectable and shows the awaiting-cache
        # guidance with Generate disabled.
        assert 'Phase A' in str(ver.value)
        _fresh('amort_version_').select(nh3_opt)
        at.run()
        errs = ' '.join(str(e.value) for e in at.error)
        infos = ' '.join(str(i.value) for i in at.info)
        assert 'Structure cache not found' in errs
        assert 'may not have been pushed to this machine yet' in infos
        gen = next(b for b in at.button
                   if 'Generate' in str(getattr(b, 'label', '')))
        assert gen.disabled
    # The NH3 model must be identifiable in the version dropdown (the
    # grouping split eats the label head, so the tail must carry the tag).
    assert 'NH' in str(nh3_opt)
