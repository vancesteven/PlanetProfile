"""Deploy snapshot manifest: every file the SBI slots need at serve time must
exist and be git-tracked in the working tree (no network, no git branch ops).
Guards the defect where the Enceladus libration factor sidecar was not staged.
"""
import importlib.util
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
PAGE = REPO / 'PlanetProfileApp/pages/Inference.py'
_spec = importlib.util.spec_from_file_location(
    'extract_sbi_cache_paths', REPO / 'plans/scripts/extract_sbi_cache_paths.py')
ex = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ex)

ENCELADUS_LIB = ('PlanetProfile/Test/mcmc_results/Enceladus/Cassini_isostasy/'
                 'enceladus_seawater_zbw_grid_2d.pkl.lib_factor_sidecar.pkl')


@pytest.fixture(scope='module')
def manifest():
    return ex.build_manifest(PAGE, REPO)


def test_enceladus_libration_sidecar_in_manifest(manifest):
    files = {f['path']: f for fs in manifest.values() for f in fs}
    assert ENCELADUS_LIB in files
    assert files[ENCELADUS_LIB]['kind'] == 'lib_factor_sidecar'
    assert files[ENCELADUS_LIB]['required']
    assert 'enceladus_ocean6D_rhorock_posterior_1m.pt' in manifest


def test_every_registry_cache_is_covered(manifest):
    caches = {f['path'] for fs in manifest.values() for f in fs
              if f['kind'] == 'cache'}
    assert caches == set(ex.extract_cache_paths(PAGE))
    assert len(caches) == 9


def test_required_files_exist_and_are_tracked(manifest):
    flat = ex.flatten_manifest(manifest)
    required = [f for f in flat if f['required']]
    assert required
    missing = [f['path'] for f in required if not (REPO / f['path']).is_file()]
    assert not missing, f'required deploy files missing: {missing}'
    present = [f['path'] for f in flat if (REPO / f['path']).is_file()]
    untracked = []
    for p in present:   # optional files that exist must also ship from git
        r = subprocess.run(['git', 'ls-files', '--error-unmatch', p],
                           cwd=REPO, capture_output=True)
        if r.returncode != 0:
            untracked.append(p)
    assert not untracked, f'deploy files not git-tracked: {untracked}'


def test_is_bundle_companions_listed(manifest):
    kinds = {f['kind'] for f in manifest['titan_freegrav_nh3_posterior_1m.pt']}
    assert {'is_fiducial_manifest', 'is_fiducial_npz'} <= kinds


def test_induction_slots_require_ae_sidecar(manifest):
    req = {f['path'] for fs in manifest.values() for f in fs
           if f['kind'] == 'ae_sidecar' and f['required']}
    assert any('Test51_seawater' in p for p in req)
    assert any('Enceladus' in p for p in req)
    assert not any('Titan' in p for p in req)
