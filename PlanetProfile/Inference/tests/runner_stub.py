"""Shared MCMCRunner stand-in for the Inference unit tests (import as
``from runner_stub import runner_stub``; not a conftest, so it cannot collide
with tests/conftest.py under the shared module name ``conftest``).

The likelihood tests build a minimal stand-in for MCMCRunner instead of a full
runner (no structure cache, no config file). Each test file used to hand-roll
its own stub; every new predicate the likelihood consults (e.g.
_tidal_solve_needed, 6c929d54) then broke them with AttributeError, and since
these files were never collected nobody noticed. ``runner_stub`` binds the
REAL predicate methods from MCMCRunner, so a stub stays in step with the
runner automatically.
"""
import inspect
import re

from PlanetProfile.Inference.mcmc_runner import MCMCRunner


def _methods_called_by(func):
    """Private ``self._name(`` methods called in ``func``'s source."""
    names = set(re.findall(r'self\.(_[A-Za-z0-9_]+)\(', inspect.getsource(func)))
    return sorted(n for n in names if callable(getattr(MCMCRunner, n, None)))


def _methods_reachable_from(func):
    """All private methods reachable from ``func`` through ``self._x(`` calls."""
    seen, todo = set(), list(_methods_called_by(func))
    while todo:
        name = todo.pop()
        if name in seen:
            continue
        seen.add(name)
        todo.extend(_methods_called_by(getattr(MCMCRunner, name)))
    return tuple(sorted(seen))


# Every method the likelihood closure can reach on ``self``, discovered from
# the real source, so the stub cannot fall behind the runner again.
_BOUND = _methods_reachable_from(MCMCRunner._make_flexible_log_likelihood)


def _init_state_defaults():
    """Private instance attributes that MCMCRunner.__init__ sets to a literal
    (None/False/True/0/{}/[]) -- lazy-cache flags such as
    _cmr2_offset_sidecar_checked that the likelihood path reads."""
    src = inspect.getsource(MCMCRunner.__init__)
    out = {}
    for name, lit in re.findall(
            r'^\s*self\.(_[A-Za-z0-9_]+)\s*(?::[^=\n]*)?=\s*(None|False|True|0|\{\}|\[\])\s*$',
            src, flags=re.M):
        out[name] = {'None': None, 'False': False, 'True': True, '0': 0}.get(lit)
        if lit == '{}':
            out[name] = {}
        elif lit == '[]':
            out[name] = []
    return out


class _RunnerStub:
    pass


def runner_stub(config, **attrs):
    """A stand-in for MCMCRunner with ``config`` and extra attributes, plus the
    real predicate methods bound to it."""
    stub = _RunnerStub()
    for k, v in _init_state_defaults().items():
        setattr(stub, k, v.copy() if isinstance(v, (dict, list)) else v)
    stub.config = config
    stub.param_names = []
    stub.param_groups = {}
    stub.fixed_params = {}
    stub._samples_D_iceIh = False
    stub._TIDAL_OBSERVABLES = MCMCRunner._TIDAL_OBSERVABLES
    for k, v in attrs.items():
        setattr(stub, k, v)
    for name in _BOUND:
        setattr(stub, name, getattr(MCMCRunner, name).__get__(stub, MCMCRunner))
    return stub
