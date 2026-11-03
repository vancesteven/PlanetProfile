"""PlanetProfile.BuildTest must run every PPTest body module that exists and
nothing that does not (testing plan T2, 2026-11-03). The old loop imported
PPTest2..PPTest{file count} (54), requesting nonexistent PPTest51-54 and
skipping the suffixed variants (PPTest3NoClath, ...)."""
import importlib.util
import os
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
TEST_DIR = REPO / 'PlanetProfile' / 'Test'


def _bodies():
    from PlanetProfile.BuildTest import ListTestBodies
    return ListTestBodies()


def test_every_listed_body_exists():
    for name in _bodies():
        assert (TEST_DIR / f'{name}.py').exists(), name
        assert importlib.util.find_spec(f'PlanetProfile.Test.{name}') is not None, name


def test_every_body_file_is_listed():
    on_disk = {f[:-3] for f in os.listdir(TEST_DIR)
               if re.fullmatch(r'PPTest.*\.py', f)
               and not any(x in f for x in ('Induct', 'Bayes', 'Explore'))}
    assert set(_bodies()) == on_disk


def test_suffixed_variants_included_and_order():
    names = _bodies()
    for v in ('PPTest3NoClath', 'PPTest41NoClath', 'PPTest46_allice', 'PPTest49_2km'):
        assert v in names
    nums = [int(n[6:]) for n in names if n[6:].isdigit()]
    assert nums == sorted(nums) and nums[0] == 1
    assert all(n[6:].isdigit() for n in names[:len(nums)])
