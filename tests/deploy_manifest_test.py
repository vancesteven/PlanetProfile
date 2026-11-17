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


def test_enceladus_slot_needs_no_libration_sidecar(manifest):
    # The Enceladus SBI slot points at the frozen-convention training snapshot
    # (legacy 1.008 libration, no elastic_factor_sidecar), so it must NOT
    # require the repaired-model F sidecar; the artifact itself is deployed.
    files = {f['path']: f for fs in manifest.values() for f in fs}
    assert ENCELADUS_LIB not in files
    assert 'enceladus_ocean6D_rhorock_posterior_1m.pt' in manifest


def test_declared_libration_sidecar_is_required(tmp_path):
    # Machinery check: a slot whose config declares elastic_factor_sidecar
    # must list it as a required lib_factor_sidecar companion.
    cfg = tmp_path / 'cfg.json'
    cfg.write_text('{"metadata": {"libration_model_correction": '
                   '{"elastic_factor_sidecar": {"suffix": ".lib.pkl"}}}}')
    page = tmp_path / 'Inference.py'
    page.write_text("_SBI_ARTIFACT_SLOTS = {'a.pt': {'cache_path': "
                    "'x/c.pkl', 'config_path': 'cfg.json'}}\n")
    m = ex.build_manifest(page, tmp_path)
    lib = [f for f in m['a.pt'] if f['kind'] == 'lib_factor_sidecar']
    assert lib == [{'path': 'x/c.pkl.lib.pkl',
                    'kind': 'lib_factor_sidecar', 'required': True}]


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


# ---------------------------------------------------------------------------
# Review fixes (2026-11-17): Ae sidecar key validation, stage-only end to end
# ---------------------------------------------------------------------------
import os
import pickle
import shutil
import sys

import numpy  # noqa: F401  (the validator unpickles caches)


def _tmp_slot(tmp_path, stored_key_mutator=None):
    """A one-slot registry whose config needs an Ae sidecar, a tiny cache,
    and a legacy sidecar carrying the key the runner would compute."""
    repo = tmp_path / 'repo'
    (repo / 'PlanetProfile/Inference').mkdir(parents=True)
    shutil.copy(REPO / 'PlanetProfile/Inference/ae_sidecar_key.py',
                repo / 'PlanetProfile/Inference/ae_sidecar_key.py')
    (repo / 'x').mkdir()
    (repo / 'cfg.json').write_text(json.dumps(
        {'observables': {'Ae_synodic_real': [1, 1], 'Re_k2': [1, 1]},
         'induction_bounds': {'orbital': [0, 1]}}))
    page = repo / 'Inference.py'
    page.write_text("_SBI_ARTIFACT_SLOTS = {'a.pt': {'cache_path': "
                    "'x/c.pkl', 'config_path': 'cfg.json'}}\n")
    cache = repo / 'x/c.pkl'
    cache.write_bytes(pickle.dumps(
        {'Tb_K_grid': [1.0, 2.0], 'structures': [{}, {}, {}]}))
    key = {'labels': ('orbital', 'synodic'), 'n': 3,
           'cache_size': cache.stat().st_size}
    if stored_key_mutator:
        key = stored_key_mutator(key)
    (repo / 'x/c.pkl.ae_sidecar.pkl').write_bytes(
        pickle.dumps({'key': key, 'ae_grid': {}}))
    return page, repo


import json  # noqa: E402


def test_ae_sidecar_validation_accepts_matching_key(tmp_path):
    page, repo = _tmp_slot(tmp_path)
    assert ex.validate_ae_sidecars(page, repo) == []


@pytest.mark.parametrize('mutate', [
    lambda k: {**k, 'cache_size': k['cache_size'] + 1},      # rebuilt cache
    lambda k: {**k, 'n': 4},
    lambda k: {**k, 'labels': ('synodic',)},                 # wrong labels
])
def test_ae_sidecar_validation_rejects_stale_or_wrong_key(tmp_path, mutate):
    page, repo = _tmp_slot(tmp_path, mutate)
    errs = ex.validate_ae_sidecars(page, repo)
    assert len(errs) == 1
    assert "slot 'a.pt'" in errs[0] and 'recompute' in errs[0]


def test_shipped_ae_sidecars_match_the_keys_their_runners_compute():
    assert ex.validate_ae_sidecars(PAGE, REPO) == []


def test_validate_ae_cli_exits_nonzero_with_message(tmp_path):
    page, repo = _tmp_slot(
        tmp_path, lambda k: {**k, 'cache_size': k['cache_size'] + 1})
    r = subprocess.run(
        [sys.executable, str(REPO / 'plans/scripts/extract_sbi_cache_paths.py'),
         str(page), '--validate-ae', '--repo-root', str(repo)],
        capture_output=True, text=True)
    assert r.returncode == 1 and 'ERROR' in r.stderr and 'recompute' in r.stderr


@pytest.mark.skipif(shutil.which('rsync') is None,
                    reason='rsync is not installed; the deploy script needs it')
def test_build_deploy_branch_stage_only_end_to_end(tmp_path):
    """Runs the real script with a RELATIVE --stage-only dir from a cwd that
    is not the repo root (the two-working-directories defect)."""
    # Stage inside tests/ (cleaned up below) so the relative path means
    # different directories from cwd=tests/ vs the repo root the script
    # cd's to while staging: the old script staged into <repo>/<rel> but
    # copied the SBI files into <cwd>/<rel>.
    cwd = REPO / 'tests'
    rel = f'_stage_only_{os.getpid()}'
    stage = cwd / rel
    assert not stage.exists() and not (REPO / rel).exists()
    try:
        _run_stage_only(cwd, rel, stage)
    finally:
        shutil.rmtree(stage, ignore_errors=True)
        shutil.rmtree(REPO / rel, ignore_errors=True)


def _run_stage_only(cwd, rel, stage):
    env = dict(os.environ, PP_PYTHON=sys.executable)
    r = subprocess.run(
        ['bash', str(REPO / 'plans/scripts/build_deploy_branch.sh'),
         '--stage-only', rel], cwd=cwd, env=env, capture_output=True,
        text=True, timeout=1800)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-2000:]
    assert (stage / 'PlanetProfileApp/pages/Inference.py').is_file()
    assert (stage / 'PlanetProfile/Inference/ae_sidecar_key.py').is_file()
    assert not (stage / '.sbi-file-manifest').exists()
    man = ex.flatten_manifest(ex.build_manifest(PAGE, REPO))
    for f in man:
        if f['required']:
            assert (stage / f['path']).is_file(), f['path']
    staged_test = {p.relative_to(stage).as_posix()
                   for p in (stage / 'PlanetProfile/Test').rglob('*')
                   if p.is_file()}
    expected = {f['path'] for f in man if f['path'].startswith(
        'PlanetProfile/Test/') and (REPO / f['path']).is_file()}
    assert staged_test == expected
    # nothing was staged relative to the repo root instead
    assert not (REPO / rel).exists()
