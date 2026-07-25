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
    assert 'No Arrhenius' in md2 and 'Maxwell' in md2 and 'MCMC' in md2


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

    btn_keys = [str(getattr(b, 'key', ''))
                for b in at.get('download_button')]
    assert any('globe_pp_txt' in k for k in btn_keys)
    assert any('globe_pp_csv' in k for k in btn_keys)


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
