"""Enceladus induction display reference (PlanetProfile/Inference/
enceladus_induction_display.py): the per-sample blend must use the ocean-
branch (zb, log10 w) bilinear weights, reproduce node values on the grid, be
linear between nodes, and give NaN (never a nearest-node fallback) when a
bracketing node is unbuilt. Synthetic Ae values on the real cache geometry."""
import pickle
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
CACHE = REPO / ('PlanetProfile/Test/mcmc_results/Enceladus/Cassini_isostasy/'
                'enceladus_seawater_zbw_grid_2d.pkl')
pytestmark = pytest.mark.skipif(not CACHE.exists(), reason='cache absent')


@pytest.fixture(scope='module')
def setup(tmp_path_factory):
    sd = pickle.load(open(CACHE, 'rb'))
    zb, w = np.asarray(sd['zb_km_grid']), np.asarray(sd['wOcean_ppt_grid'])
    n = len(sd['structures'])
    nw = len(w)

    def val(i):     # linear in (zb, log10 w): bilinear reproduces it exactly
        return complex(0.01 * zb[i // nw], 0.1 * np.log10(w[i % nw]))
    grid = {i: ({'true anomaly': val(i), 'synodic': 2 * val(i)}
                if sd['structures'][i] else None) for i in range(n)}
    p = tmp_path_factory.mktemp('ae') / 'side.pkl'
    pickle.dump({'key': {'n': n}, 'ae_grid': grid, 'ae_grid_core25': grid,
                 'core25_Sm': 25.0}, open(p, 'wb'))
    return sd, zb, w, p


def _run(sd, p, zbs, ws):
    from PlanetProfile.Inference.enceladus_induction_display import (
        induction_display)
    th = np.column_stack([zbs, np.log10(ws)])
    return induction_display(th, ['zb_km', 'log10_wOcean_ppt'], sd, p)


def test_nodes_and_midpoints_reproduce_linear_field(setup):
    sd, zb, w, p = setup
    nw = len(w)
    S = sd['structures']
    # pick an interior cell whose four corners are all built
    for iz in range(40, len(zb) - 1):
        for iw in range(10, nw - 1):
            c = [iz * nw + iw, iz * nw + iw + 1, (iz + 1) * nw + iw,
                 (iz + 1) * nw + iw + 1]
            if all(S[k] for k in c):
                break
        else:
            continue
        break
    zm = 0.5 * (zb[iz] + zb[iz + 1])
    wm = 10 ** (0.5 * (np.log10(w[iw]) + np.log10(w[iw + 1])))
    d = _run(sd, p, [zb[iz], zm], [w[iw], wm])
    o = d['curves']['pp_core']['orbital']
    got = np.asarray(o['Ae_re']) + 1j * np.asarray(o['Ae_im'])
    exp = [complex(0.01 * z, 0.1 * np.log10(x)) for z, x in
           ((zb[iz], w[iw]), (zm, wm))]
    np.testing.assert_allclose(got, exp, rtol=0, atol=1e-12)
    sy = d['curves']['pp_core']['synodic']
    np.testing.assert_allclose(np.asarray(sy['Ae_re']), 2 * got.real, atol=1e-12)


def test_unbuilt_corner_gives_nan_not_fallback(setup):
    sd, zb, w, p = setup
    nw = len(w)
    S = sd['structures']
    k = next(i for i, s in enumerate(S) if s is None
             and 0 < i // nw < len(zb) - 1 and 0 < i % nw < nw - 1)
    iz, iw = k // nw, k % nw
    zq = 0.5 * (zb[iz] + zb[iz + 1])
    wq = 10 ** (0.5 * (np.log10(w[iw]) + np.log10(w[iw + 1])))
    d = _run(sd, p, [zq], [wq])
    assert np.isnan(d['curves']['pp_core']['orbital']['Ae_re'][0])
    assert np.isnan(d['sigma_ocean_Sm'][0])


def test_out_of_box_nan_and_summary(setup):
    from PlanetProfile.Inference.enceladus_induction_display import summarize
    sd, zb, w, p = setup
    d = _run(sd, p, [zb[-1] + 1.0, 27.0], [3.0, 3.0])
    assert np.isnan(d['curves']['pp_core']['orbital']['Ae_re'][0])
    s = summarize(d)
    assert s['pp_core']['orbital']['A']['n'] == 1
    assert s['pp_core']['orbital']['B_nT']['p50'] == pytest.approx(
        5.0 * s['pp_core']['orbital']['A']['p50'])


def test_missing_sidecar_returns_none(setup, tmp_path):
    sd, zb, w, _ = setup
    assert _run(sd, tmp_path / 'absent.pkl', [27.0], [3.0]) is None
