"""Tests for the traceable Fe-S core EOS builder (Thermodynamics/Perple_X/stackFeSTraceable.py).

The builder stacks the five werami tables in Thermodynamics/Perple_X/werami_tables with exact
unit conversions. These tests check that it reproduces the reference product used by the lunar
interior study (same semantic SHA-256), that its conversions are exact, and that its Fe95-S5
slice is physical, unlike the slice of the distributed Fe-S_3D_EOS.mat.
"""
import importlib.util
import os

import numpy as np
import pytest
from scipy.io import loadmat

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
spec = importlib.util.spec_from_file_location(
    'stackFeSTraceable', os.path.join(ROOT, 'Thermodynamics', 'Perple_X', 'stackFeSTraceable.py'))
FeS = importlib.util.module_from_spec(spec); spec.loader.exec_module(FeS)
DIST_MAT = os.path.join(ROOT, 'PlanetProfile', 'Thermodynamics', 'EOStables', 'Perple_X', 'Fe-S_3D_EOS.mat')


@pytest.fixture(scope='module')
def stacked():
    return FeS.StackTraceableFeS()


def test_reproduces_reference_product(stacked):
    assert FeS.SemanticSHA256(stacked) == FeS.SEMANTIC_SHA256


def test_shapes_and_composition_axis(stacked):
    assert stacked['rho_kgm3'].shape == (5, 317, 317)
    assert np.array_equal(stacked['wFe_ppt'], [800, 850, 900, 950, 1000])
    assert np.array_equal(stacked['wFe_ppt'] + stacked['wS_ppt'], np.full(5, 1000.0))


def test_pressure_axis_uses_exact_bar_conversion(stacked):
    src = FeS.ParseWeramiTable(os.path.join(FeS.WERAMI_DIR, '95Fe5S_core.tab'))
    assert np.array_equal(stacked['P_MPa'], src['P_bar'] * 0.1)
    assert stacked['P_MPa'][0] == pytest.approx(4000.0)


def test_fe95s5_slice_is_physical(stacked):
    i = 3  # wFe_ppt = 950
    assert np.median(stacked['KS_GPa'][i]) > 50
    assert np.median(stacked['GS_GPa'][i]) > 10


@pytest.mark.skipif(not os.path.exists(DIST_MAT), reason='distributed Fe-S_3D_EOS.mat not present')
def test_distributed_table_fe95s5_slice_is_broken():
    """Documents why the traceable table exists: the distributed slice has near-zero moduli."""
    d = loadmat(DIST_MAT)
    i = int(np.argmin(np.abs(d['wFe_ppt'][0] - 950)))
    assert np.median(d['KS_GPa'][i]) < 1e-2 * np.median(d['KS_GPa'][-1])


def test_written_file_round_trips(stacked, tmp_path):
    fpath = tmp_path / 'fes.mat'
    FeS.WriteTraceableFeS(stacked, str(fpath))
    d = loadmat(str(fpath))
    for k in FeS.FIELDS:
        assert np.array_equal(np.squeeze(d[k]), np.squeeze(stacked[k]))
