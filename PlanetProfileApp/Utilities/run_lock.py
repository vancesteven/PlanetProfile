"""Global run lock for compute on shared (public) deployments.

One PlanetProfile forward-model run at a time per container: the Space has
2 vCPUs and a single Streamlit process, so concurrent runs degrade every
visitor's session. Lock = O_CREAT|O_EXCL file in the system temp dir with
a timestamp payload; stale locks (crashed runs) expire after TTL.
"""
import json
import logging
import os
import tempfile
import time
from contextlib import contextmanager

log = logging.getLogger('PlanetProfile')

_LOCK_TTL_S = 15 * 60  # generous: a single PP run is minutes


def _lock_path(name):
    return os.path.join(tempfile.gettempdir(), f'pp_app_lock_{name}.lock')


def acquire(name='pp_run', ttl_s=_LOCK_TTL_S):
    """Try to acquire the named lock. Returns True on success."""
    path = _lock_path(name)
    # Expire stale locks from crashed runs.
    try:
        with open(path) as f:
            stamped = json.load(f).get('t', 0)
        if time.time() - stamped > ttl_s:
            os.remove(path)
            log.warning(f"Removed stale run lock {path}")
    except FileNotFoundError:
        pass
    except Exception:
        # Unreadable lock: treat as stale only if old by mtime.
        try:
            if time.time() - os.path.getmtime(path) > ttl_s:
                os.remove(path)
        except OSError:
            pass
    try:
        fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return False
    with os.fdopen(fd, 'w') as f:
        json.dump({'t': time.time(), 'pid': os.getpid()}, f)
    return True


def release(name='pp_run'):
    try:
        os.remove(_lock_path(name))
    except FileNotFoundError:
        pass


@contextmanager
def run_lock(name='pp_run', ttl_s=_LOCK_TTL_S):
    """Context manager: yields True if the lock was acquired (and releases
    it on exit), False if another run holds it (nothing to release)."""
    got = acquire(name, ttl_s)
    try:
        yield got
    finally:
        if got:
            release(name)
