"""Tests for BuildLayeredPerplexTable (PlanetProfile/Thermodynamics/LayeredMantle.py).

Synthetic tables check the row selection and the grid checks. A real-table test checks that
the builder reproduces, byte for byte, a layered lunar table made by the lunar interior study
(moon_w05_preferred_L500A_stx24_unified.tab): rebuilding it from itself above the interface and
from the pristine table below must give it back.
"""
import filecmp
import os

import pytest

from PlanetProfile import _PERPLEXCACHE
from PlanetProfile.Thermodynamics.LayeredMantle import BuildLayeredPerplexTable

HDR = ['|6.6.6', 'NAME_1.tab', '           2', 'T(K)    ', '   300', '   100', '   2', 'P(bar)  ',
       '   1000', '   1000', '   3', '           3', 'T(K)           P(bar)         rho,kg/m3']


def _write(dirpath, name, rho, nodes=None, cols=None):
    nodes = nodes or [(300, 1000), (400, 1000), (300, 2000), (400, 2000), (300, 3000), (400, 3000)]
    hdr = list(HDR); hdr[1] = f'{name[:-4]}_1.tab'
    if cols: hdr[-1] = cols
    rows = [f'{T}  {P}  {r}' for (T, P), r in zip(nodes, rho)]
    with open(os.path.join(dirpath, name), 'w') as f:
        f.write('\n'.join(hdr + rows) + '\n')


def test_rows_split_at_interface(tmp_path):
    _write(tmp_path, 'up.tab', [1, 1, 1, 1, 1, 1]); _write(tmp_path, 'lo.tab', [2, 2, 2, 2, 2, 2])
    out = BuildLayeredPerplexTable('up.tab', 'lo.tab', 250.0, 'lay.tab', EOSdir=str(tmp_path))
    assert (out['nUpperRows'], out['nLowerRows']) == (4, 2)
    rows = open(out['path']).read().splitlines()[13:]
    assert [r.split()[2] for r in rows] == ['1', '1', '1', '1', '2', '2']
    assert open(out['path']).read().splitlines()[1].startswith('lay_1.tab')


def test_interface_at_node_goes_to_lower_table(tmp_path):
    _write(tmp_path, 'up.tab', [1] * 6); _write(tmp_path, 'lo.tab', [2] * 6)
    out = BuildLayeredPerplexTable('up.tab', 'lo.tab', 200.0, 'lay.tab', EOSdir=str(tmp_path))
    assert (out['nUpperRows'], out['nLowerRows']) == (2, 4)


def test_different_grids_raise(tmp_path):
    _write(tmp_path, 'up.tab', [1] * 6)
    _write(tmp_path, 'lo.tab', [2] * 6, nodes=[(300, 1000), (400, 1000), (300, 2000), (400, 2000), (300, 3500), (400, 3500)])
    with pytest.raises(ValueError):
        BuildLayeredPerplexTable('up.tab', 'lo.tab', 250.0, 'lay.tab', EOSdir=str(tmp_path))


def test_columns_must_be_T_and_P_in_bar(tmp_path):
    _write(tmp_path, 'up.tab', [1] * 6, cols='T(K)           P(kbar)        rho')
    _write(tmp_path, 'lo.tab', [2] * 6, cols='T(K)           P(kbar)        rho')
    with pytest.raises(ValueError):
        BuildLayeredPerplexTable('up.tab', 'lo.tab', 250.0, 'lay.tab', EOSdir=str(tmp_path))


def test_existing_output_is_not_overwritten(tmp_path):
    _write(tmp_path, 'up.tab', [1] * 6); _write(tmp_path, 'lo.tab', [2] * 6)
    BuildLayeredPerplexTable('up.tab', 'lo.tab', 250.0, 'lay.tab', EOSdir=str(tmp_path))
    with pytest.raises(FileExistsError):
        BuildLayeredPerplexTable('up.tab', 'lo.tab', 250.0, 'lay.tab', EOSdir=str(tmp_path))


LAY, M0 = 'moon_w05_preferred_L500A_stx24_unified.tab', 'moon_w05_preferred_stx24_unified.tab'
PINT_MPA = 2333.8058874992266  # recorded by the lunar study when it built LAY


@pytest.mark.skipif(not all(os.path.exists(os.path.join(_PERPLEXCACHE, f)) for f in (LAY, M0)),
                    reason='lunar tables not in the Perple_X cache')
def test_reproduces_lunar_layered_table(tmp_path):
    import shutil
    for f in (LAY, M0):
        shutil.copy2(os.path.join(_PERPLEXCACHE, f), tmp_path / ('src_' + f))
    out = BuildLayeredPerplexTable('src_' + LAY, 'src_' + M0, PINT_MPA, LAY, EOSdir=str(tmp_path))
    assert (out['nUpperRows'], out['nLowerRows']) == (41100, 48900)
    assert filecmp.cmp(out['path'], os.path.join(_PERPLEXCACHE, LAY), shallow=False)
