"""Loader/validator for the shipped fiducial importance-sampling bundle.

The Titan NH3 SBI slot ships a precomputed IS correction at the config-
central observation vector (built by plans/scripts/
build_nh3_is_fiducial_bundle.py, ratified in validation_reports/
titan_freegrav_nh3_1m/RATIFICATION.md, 2026-10-20; handoff §0.32,
conditions B1-B4).  ``load_fiducial_correction`` returns a
``CorrectedBundle`` only when ALL of the following hold, else a
``Fallback`` carrying a human-readable reason (it never raises for the
expected fallback cases):

  B2  exact match of every x_obs value (np.isclose rtol 1e-12, identical
      key sets), the observation sigmas, the observable set, and the
      artifact / config / structure-cache hashes;
  B3  diagnostics RECOMPUTED from the shipped log_w satisfy ESS >= 1000,
      Pareto-k <= 0.7, w_max <= 1%, rejected fraction < 50%; the C5.3
      reverse-coverage record is present and passing; npz sha256 matches
      the manifest;
  B4  the request is a full-N conditioning (``n_derived`` None); a
      truncated derived recompute can never receive corrected weights.

``FEATURE_ENABLED`` is the module-level switch the GUI reads (Phase 2);
default False.  The loader itself does not consult it.
"""
from __future__ import annotations

import hashlib
import json
import logging
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Dict, Optional, Union

import numpy as np

from .is_correction import (
    ESS_ABS_FLOOR, PARETO_K_MAX, W_MAX_FRAC_CAP, REJECT_FRAC_HARD_FAIL,
    compute_is_correction)

log = logging.getLogger('PlanetProfile')

# Read by the GUI (Phase 2). Default OFF.
FEATURE_ENABLED = False

# Approved fallback text categories.
R_MODIFIED = 'modified datum'
R_ESS = 'ESS < 1000'
R_PARETO = 'Pareto-k > 0.7'
R_WMAX = 'max weight > 1%'
R_FULLN = 'full-N likelihoods absent'
R_INTEGRITY = 'bundle integrity'

_REQUIRED_ARRAYS = ('samples', 'log_w', 'log_likelihoods', 'flow_log_prob')
_DIAG_RTOL = 1e-9


@dataclass
class Fallback:
    reason: str

    def __bool__(self):  # falsy so `if not result:` reads naturally
        return False


@dataclass
class CorrectedBundle:
    samples: np.ndarray            # (N, n_params), the validated draws
    weights: np.ndarray            # normalized IS weights (PSIS-smoothed if
    #                                the recomputation smoothed), sum to 1
    log_w: np.ndarray
    derived: Dict[str, np.ndarray]  # per-sample derived arrays (length N)
    param_names: list
    diagnostics: Dict[str, Any]    # recomputed from shipped log_w
    manifest: Dict[str, Any] = field(repr=False, default_factory=dict)

    def __bool__(self):
        return True


_sha_cache: Dict[tuple, str] = {}


def _sha256(path) -> str:
    p = Path(path)
    st = p.stat()
    key = (str(p.resolve()), st.st_mtime_ns, st.st_size)
    if key not in _sha_cache:
        h = hashlib.sha256()
        with open(p, 'rb') as f:
            for chunk in iter(lambda: f.read(1 << 20), b''):
                h.update(chunk)
        _sha_cache[key] = h.hexdigest()
    return _sha_cache[key]


def _config_hash(config) -> str:
    if isinstance(config, str):
        return config
    return str(config.generate_hash())


def _same_floats(a: Dict[str, Any], b: Dict[str, Any], rtol=1e-12) -> bool:
    if set(a) != set(b):
        return False
    return all(np.isclose(float(a[k]), float(b[k]), rtol=rtol, atol=0.0)
               for k in a)


def recompute_diagnostics(log_w, derived=None, k2_support_bounds=None):
    """Recompute ESS / Pareto-k / w_max / rejected fraction (and branch
    split) from shipped raw log-weights using the validated
    ``compute_is_correction`` code path.  rejected draws carry -inf, which
    the sentinel mask treats as rejected."""
    derived = derived or {}
    lw = np.asarray(log_w, float)
    # -inf (rejected) -> -1e30 sentinel so compute_is_correction's
    # finite-likelihood guard passes and the rejection mask is identical.
    fake = SimpleNamespace(
        log_likelihoods=np.where(np.isfinite(lw), lw, -1e30),
        metadata={'flow_log_prob': np.zeros(len(lw))},
        D_ocean_results=derived.get('D_ocean_results', np.array([])),
        k2_results=derived.get('k2_results', np.array([])))
    return compute_is_correction(
        fake, k2_support_bounds=k2_support_bounds or None)


def load_fiducial_correction(
        manifest_path, *, artifact_path, config, cache_path, x_obs,
        obs_sigma, n_post=None, n_derived=None
) -> Union[CorrectedBundle, Fallback]:
    """Validate and load the shipped fiducial IS bundle (see module doc).

    Args:
        manifest_path: the bundle's JSON sidecar (npz sits next to it).
        artifact_path / cache_path: files whose sha256 must match.
        config: InferenceConfig (``generate_hash()``) or a hash string.
        x_obs / obs_sigma: {observable: value} dicts of the live request.
        n_post: requested posterior size; the bundle is a fixed-N draw
            set, so anything other than None or the bundle's N falls back.
        n_derived: conditioning derived-recompute size. Must be None (full
            N); a truncated recompute has NaN-padded likelihoods and can
            never be paired with corrected weights (B4).
    """
    try:
        return _load(manifest_path, artifact_path, config, cache_path,
                     x_obs, obs_sigma, n_post, n_derived)
    except (OSError, ValueError, KeyError, json.JSONDecodeError,
            TypeError, AssertionError, RuntimeError) as e:
        return Fallback(f'{R_INTEGRITY}: {type(e).__name__}: {e}')


def _load(manifest_path, artifact_path, config, cache_path, x_obs,
          obs_sigma, n_post, n_derived):
    mp = Path(manifest_path)
    with open(mp) as f:
        m = json.load(f)
    N = int(m['n_samples'])

    # B4: full-N only (checked first: it is a property of the request).
    if n_derived is not None and n_derived < N:
        return Fallback(f'{R_FULLN}: n_derived={n_derived} < N={N}; '
                        'corrected weights require full-N likelihoods '
                        '(n_derived=None).')
    if n_derived is not None and n_derived > N:
        return Fallback(f'{R_FULLN}: n_derived={n_derived} does not '
                        f'correspond to the bundle draw set (N={N}).')
    if n_post is not None and int(n_post) != N:
        return Fallback(f'{R_MODIFIED}: n_post={n_post} differs from the '
                        f'validated draw count N={N}.')

    # B2: exact datum + provenance match.
    if not _same_floats(dict(x_obs), m['x_obs']):
        return Fallback(f'{R_MODIFIED}: x_obs differs from the validated '
                        'fiducial (values or observable set).')
    if not _same_floats(dict(obs_sigma), m['obs_sigma']):
        return Fallback(f'{R_MODIFIED}: observation sigmas differ from the '
                        'validated fiducial.')
    if sorted(x_obs) != sorted(m['observables']):
        return Fallback(f'{R_MODIFIED}: observable set differs.')
    if _config_hash(config) != m['config_hash']:
        return Fallback(f'{R_MODIFIED}: config hash differs from the '
                        'validated fiducial.')
    if _sha256(artifact_path) != m['artifact_sha256']:
        return Fallback(f'{R_MODIFIED}: SBI artifact hash differs from the '
                        'validated fiducial.')
    if _sha256(cache_path) != m['structure_cache_sha256']:
        return Fallback(f'{R_MODIFIED}: structure-cache hash differs from '
                        'the validated fiducial.')

    # Bundle integrity.
    npz_path = mp.parent / m['npz_file']
    if not npz_path.exists():
        return Fallback(f'{R_INTEGRITY}: {npz_path.name} missing.')
    if _sha256(npz_path) != m['npz_sha256']:
        return Fallback(f'{R_INTEGRITY}: npz sha256 does not match the '
                        'manifest.')
    with np.load(npz_path, allow_pickle=False) as z:
        arrs = {k: z[k] for k in z.files}
    for k in _REQUIRED_ARRAYS:
        if k not in arrs or len(arrs[k]) != N:
            return Fallback(f'{R_INTEGRITY}: array {k!r} missing or wrong '
                            'length.')
    log_w = np.asarray(arrs['log_w'], float)
    ll = np.asarray(arrs['log_likelihoods'], float)
    flp = np.asarray(arrs['flow_log_prob'], float)
    if not np.isfinite(ll).all():
        return Fallback(f'{R_FULLN}: shipped log-likelihoods contain '
                        'non-finite (NaN-padded) entries.')
    fin = np.isfinite(log_w)
    if not (np.allclose(log_w[fin], (ll - flp)[fin], rtol=1e-10, atol=0.0)
            and np.array_equal(~fin, ll < -1e29)):
        return Fallback(f'{R_INTEGRITY}: log_w inconsistent with shipped '
                        'log_likelihoods - flow_log_prob.')
    derived = {k: v for k, v in arrs.items()
               if k not in _REQUIRED_ARRAYS and k not in
               ('log_w', 'param_names') and len(v) == N}
    param_names = [str(s) for s in arrs.get('param_names', [])]

    # B3: recompute diagnostics from shipped log_w.
    c53 = m.get('c5_3')
    if not isinstance(c53, dict) or c53.get('pass') is not True:
        return Fallback(f'{R_INTEGRITY}: C5.3 reverse-coverage record '
                        'absent or not passing.')
    k2b = None
    try:
        cfgd = json.load(open(Path(__file__).resolve().parents[2]
                              / m['config_file']))
        k2b = cfgd.get('sampler_settings', {}).get('k2_support_bounds')
    except Exception:
        k2b = None
    try:
        corr = recompute_diagnostics(log_w, derived, k2b)
    except (ValueError, RuntimeError) as e:
        return Fallback(f'{R_INTEGRITY}: diagnostics recomputation failed '
                        f'({e}).')
    if corr.frac_rejected >= REJECT_FRAC_HARD_FAIL:
        return Fallback(f'{R_INTEGRITY}: rejected fraction '
                        f'{corr.frac_rejected:.1%} >= 50%.')
    if not np.isfinite(corr.pareto_k) or corr.pareto_k > PARETO_K_MAX:
        return Fallback(f'{R_PARETO}: recomputed k={corr.pareto_k:.3f}.')
    if corr.ess < ESS_ABS_FLOOR:
        return Fallback(f'{R_ESS}: recomputed ESS={corr.ess:.0f}.')
    if corr.w_max_frac > W_MAX_FRAC_CAP:
        return Fallback(f'{R_WMAX}: recomputed w_max='
                        f'{corr.w_max_frac:.4f}.')
    md = m.get('diagnostics', {})
    for name, val in (('ess', corr.ess), ('pareto_k', corr.pareto_k),
                      ('w_max_frac', corr.w_max_frac),
                      ('frac_rejected', corr.frac_rejected)):
        if name not in md or not np.isclose(md[name], val,
                                            rtol=_DIAG_RTOL, atol=1e-12):
            return Fallback(f'{R_INTEGRITY}: recomputed {name} ({val}) '
                            f'disagrees with manifest ({md.get(name)}).')
    if corr.verdict == 'fail':
        return Fallback(f'{R_INTEGRITY}: recomputed verdict fail '
                        f'({"; ".join(corr.fail_reasons)}).')

    diag = corr.summary()
    return CorrectedBundle(
        samples=np.asarray(arrs['samples'], float), weights=corr.weights,
        log_w=log_w, derived=derived, param_names=param_names,
        diagnostics=diag, manifest=m)
