"""Shared pytest configuration (testing-architecture plan, Phase 0).

plans/active/testing-architecture-plan.md:

* Thread pinning before numpy/torch import, so floating-point results are
  reproducible across runs and machines.
* A "skipped for missing data" summary and a ``--require-data`` mode. Locally,
  a test whose input data is absent is skipped with a reason (as before). In
  CI ``--require-data`` turns every such skip into a failure, so a missing
  cache can never make a tier pass silently. Data skips are recognised by
  their reason ("... absent" / "... not present"), the convention every
  existing data-dependent skip in tests/ already follows; new data-dependent
  tests should use the ``data_asset`` fixture below.
* ``data_asset``: resolves a data path (repo-relative, or under PP_DATA_ROOT),
  skips or fails per the mode, and optionally checks a sha256.
"""
import hashlib
import os
import re
from pathlib import Path

# Pin BLAS/OpenMP threads before numpy is imported by any test module.
for _var in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
             'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ.setdefault(_var, '1')

import pytest  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
_DATA_SKIP = re.compile(r'\b(absent|not present|not found|missing data)\b', re.I)
_data_skips = []


# Tests measured at > ~45 s on the reference Mac (2026-11-02 baseline) that
# lack an explicit @pytest.mark.slow. Kept here, not in the files, so the
# tier boundary lives in one place.
_SLOW_NODEIDS = (
    'tests/enceladus_isostasy_code_gaps_test.py::test_build_zbw_grid_cache_smoke_3x2',
    'tests/enceladus_isostasy_code_gaps_test.py::test_build_zbw_grid_cache_frozen_if2_audit_smoke',
    'tests/coolprop_nh3_test.py::test_single_liquid_band_machine_b',
    'tests/coolprop_nh3_test.py::test_hp_ice_straddle_end_to_end',
)


def pytest_collection_modifyitems(config, items):
    """Tier markers: 'gui' for any test module that drives Streamlit AppTest;
    'slow' for the measured slow tests above."""
    gui_modules = {}
    for item in items:
        mod = item.module
        if mod not in gui_modules:
            try:
                src = Path(mod.__file__).read_text(encoding='utf-8')
            except (OSError, TypeError):
                src = ''
            gui_modules[mod] = 'AppTest' in src
        if gui_modules[mod]:
            item.add_marker(pytest.mark.gui)
        if any(item.nodeid.endswith(n.split('tests/', 1)[-1]) or item.nodeid == n
               for n in _SLOW_NODEIDS):
            item.add_marker(pytest.mark.slow)


def pytest_addoption(parser):
    parser.addoption(
        '--require-data', action='store_true', default=False,
        help='Fail (instead of skip) any test whose input data is absent. '
             'Always set in CI.')


def _is_data_skip(report):
    if not report.skipped:
        return False
    longrepr = report.longrepr
    reason = longrepr[2] if isinstance(longrepr, tuple) else str(longrepr)
    return bool(_DATA_SKIP.search(reason)), reason


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    outcome = yield
    report = outcome.get_result()
    if report.when not in ('setup', 'call') or not report.skipped:
        return
    hit = _is_data_skip(report)
    if not hit or not hit[0]:
        return
    reason = hit[1]
    _data_skips.append((item.nodeid, reason))
    if item.config.getoption('--require-data'):
        report.outcome = 'failed'
        report.longrepr = (f'--require-data: input data missing, test would '
                           f'have been skipped: {reason}')


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    if not _data_skips:
        return
    mode = ('FAILED (--require-data)' if config.getoption('--require-data')
            else 'skipped')
    terminalreporter.section(f'tests {mode} for missing input data')
    for nodeid, reason in _data_skips:
        terminalreporter.line(f'{nodeid}: {reason}')


def _sha256(path, chunk=1 << 20):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(chunk), b''):
            h.update(block)
    return h.hexdigest()


@pytest.fixture
def data_asset(request):
    """Return a callable ``get(relpath, sha256=None) -> Path``.

    Looks for ``relpath`` under the repo, then under ``$PP_DATA_ROOT``. Absent
    -> skip (reason contains "absent", so --require-data fails it). A given
    sha256 that does not match is always a failure: a wrong file is worse
    than a missing one.
    """
    def get(relpath, sha256=None):
        candidates = [REPO / relpath]
        if os.environ.get('PP_DATA_ROOT'):
            candidates.append(Path(os.environ['PP_DATA_ROOT']) / relpath)
        path = next((p for p in candidates if p.exists()), None)
        if path is None:
            pytest.skip(f'data asset absent: {relpath}')
        if sha256 is not None and _sha256(path) != sha256:
            pytest.fail(f'data asset {relpath} sha256 mismatch '
                        f'(expected {sha256[:12]}...)')
        return path
    return get
