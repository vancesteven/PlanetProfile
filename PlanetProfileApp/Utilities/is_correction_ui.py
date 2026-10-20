"""GUI-side glue for the precomputed NH3 fiducial IS correction (Phase 2).

Binding spec: plans/MACHINE-B-HANDOFF.md section 0.32 and
validation_reports/titan_freegrav_nh3_1m/RATIFICATION.md (2026-10-20),
conditions B5-B9.  Everything here is inert unless
``PlanetProfile.Inference.is_fiducial_bundle.FEATURE_ENABLED`` is True
(read at CALL time) AND the model slot carries an ``is_correction`` entry.

Design (see Inference.py):
  * On Run the slot's bundle is validated by ``load_fiducial_correction``.
    A ``CorrectedBundle`` is turned into an ``InferenceResult`` DIRECTLY
    from the shipped draws + weights (never re-drawn from the flow: the
    torch RNG is not portable across builds, so fresh draws could not be
    paired with the shipped weights).
  * The raw result keeps the N shipped draws + normalized weights
    (``result.weights``, ``metadata['is_weights']``).  Display panels see
    ``display_view(result)``: a seeded systematic resample (equal weight),
    so every panel that reads ``result.samples`` / derived arrays shows
    the weighted posterior.  Numerical summaries use the raw weights via
    ``result_weights``.
  * B7: every number in warning text comes from the manifest / recomputed
    diagnostics / weighted bundle.  B9: there is no x-keyed path.
"""
from __future__ import annotations

import dataclasses
import hashlib
import json
import math
import time
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np

# Seed recorded in the result for the display resample (sample picker,
# k2 cloud, corner...).  Fixed so a re-render is reproducible.
DISPLAY_SEED = 0

R_TRUNCATED = ('modified datum: truncation bounds are set (the correction '
               'is validated only for the untruncated flow posterior).')


# ---------------------------------------------------------------------------
# Flag / slot gating
# ---------------------------------------------------------------------------

def feature_enabled() -> bool:
    """Read the module flag at call time so tests can monkeypatch it."""
    from PlanetProfile.Inference import is_fiducial_bundle as fb
    return bool(fb.FEATURE_ENABLED)


def correction_active_for(slot_or_spec: Optional[dict]) -> bool:
    return bool(feature_enabled() and slot_or_spec
                and slot_or_spec.get('is_correction'))


# ---------------------------------------------------------------------------
# Weights helpers (B5)
# ---------------------------------------------------------------------------

def result_weights(result) -> Optional[np.ndarray]:
    """Normalized IS weights of a result, or None when the result is not
    IS-corrected (or is an already-resampled display view)."""
    if result is None:
        return None
    md = getattr(result, 'metadata', None) or {}
    if not md.get('is_corrected'):
        return None
    w = getattr(result, 'weights', None)
    if w is None:
        w = md.get('is_weights')
    if w is None:
        return None
    w = np.asarray(w, float)
    s = w.sum()
    if w.ndim != 1 or len(w) != len(result.samples) or not np.isfinite(s) \
            or s <= 0:
        return None
    return w / s


def weighted_mean_std(v, w):
    v = np.asarray(v, float)
    w = np.asarray(w, float)
    ok = np.isfinite(v) & (w > 0)
    v, w = v[ok], w[ok]
    w = w / w.sum()
    mu = float(np.sum(w * v))
    return mu, float(np.sqrt(np.sum(w * (v - mu) ** 2)))


def weighted_summary(v, w):
    """(mean, median, std, p16, p84) weighted (w=None -> unweighted)."""
    from PlanetProfile.Inference.is_correction import weighted_quantile
    v = np.asarray(v, float)
    if w is None:
        return (float(np.mean(v)), float(np.median(v)), float(np.std(v)),
                float(np.percentile(v, 16)), float(np.percentile(v, 84)))
    mu, sd = weighted_mean_std(v, w)
    p16, p50, p84 = weighted_quantile(v, w, [0.16, 0.5, 0.84])
    return mu, float(p50), sd, float(p16), float(p84)


def weighted_ocean_fraction(result) -> Optional[float]:
    """Weighted has_ocean mass (corrected path only)."""
    w = result_weights(result)
    ho = (getattr(result, 'metadata', None) or {}).get('has_ocean')
    if w is None or ho is None or len(ho) != len(w):
        return None
    return float(np.sum(w * np.asarray(ho, float)))


_PER_SAMPLE_META = ('flow_log_prob', 'has_ocean')


def display_view(result, seed: int = DISPLAY_SEED):
    """Equal-weight view of an IS-corrected result for the display panels.

    Systematic resample (seeded, recorded) of ALL per-sample arrays with
    the same indices, so cross-array pairing is preserved.  Heating
    (computed on an uncorrected flow subset) is dropped.  Returns the input
    unchanged for non-corrected results (flag-off path is byte-identical).
    """
    w = result_weights(result)
    if w is None:
        return result
    from PlanetProfile.Inference.is_correction import systematic_resample
    n = len(w)
    idx = systematic_resample(w, n, seed=seed)

    def _take(a):
        if a is None:
            return None
        arr = np.asarray(a)
        return arr[idx] if arr.ndim >= 1 and len(arr) == n else a

    changes = {}
    for f in dataclasses.fields(result):
        if f.name in ('config', 'param_names', 'param_labels', 'metadata',
                      'convergence_metrics', 'heating_results', 'weights'):
            continue
        v = getattr(result, f.name)
        if isinstance(v, dict):
            changes[f.name] = {k: _take(x) for k, x in v.items()}
        else:
            changes[f.name] = _take(v)
    md = dict(result.metadata)
    for k in _PER_SAMPLE_META:
        if k in md:
            md[k] = _take(md[k])
    ind = md.get('induction_Ae')
    if isinstance(ind, dict):
        md['induction_Ae'] = {
            lab: {part: (list(np.asarray(vals)[idx])
                         if len(vals) == n else vals)
                  for part, vals in d.items()} if isinstance(d, dict) else d
            for lab, d in ind.items()}
    md.pop('is_weights', None)   # the view is equal-weight by construction
    diag = md.get('is_diagnostics') or {}
    md.update({
        'heating_indices': None,
        'derived_indices': np.arange(n),
        'is_display_resample': {'seed': int(seed), 'n_out': int(n),
                                'ess': diag.get('ess'),
                                'method': 'systematic'},
    })
    return dataclasses.replace(
        result, heating_results=None, weights=None, metadata=md, **changes)


# ---------------------------------------------------------------------------
# Loader glue
# ---------------------------------------------------------------------------

def _repo_path(repo_root, p) -> Path:
    p = Path(p)
    return p if p.is_absolute() else Path(repo_root) / p


def read_manifest(spec_or_slot, repo_root) -> dict:
    isc = spec_or_slot['is_correction']
    with open(_repo_path(repo_root, isc['manifest'])) as f:
        return json.load(f)


def attempt_correction(spec, repo_root, InferenceConfig):
    """Return None (feature inactive for this spec), a CorrectedBundle, or
    a Fallback.  ``spec`` needs: is_correction, config_path, artifact_path,
    cache_path, x_obs, obs_sigma, n_posterior_samples, n_derived,
    truncate_bounds."""
    if not correction_active_for(spec):
        return None
    from PlanetProfile.Inference import is_fiducial_bundle as fb
    if spec.get('truncate_bounds'):
        return fb.Fallback(R_TRUNCATED)
    try:
        # Hash of the TRAINING config carrying the LIVE observables: GUI-only
        # knobs (seed, n_reeval, cache path) do not change the posterior the
        # bundle encodes, but every observable value / sigma does.
        cfg_dict = json.loads(
            _repo_path(repo_root, spec['config_path']).read_text())
        cfg_dict['mode'] = 'sbi'
        cfg_dict['observables'] = {
            n: [spec['x_obs'][n], spec['obs_sigma'][n]]
            for n in spec['x_obs']}
        cfg_hash = InferenceConfig.from_dict(cfg_dict).generate_hash()
        n_post = int(spec['n_posterior_samples'])
        n_der = spec.get('n_derived')
        if n_der is not None and int(n_der) >= n_post:
            n_der = None
        return fb.load_fiducial_correction(
            _repo_path(repo_root, spec['is_correction']['manifest']),
            artifact_path=_repo_path(repo_root, spec['artifact_path']),
            config=cfg_hash,
            cache_path=_repo_path(repo_root, spec['cache_path']),
            x_obs=spec['x_obs'], obs_sigma=spec['obs_sigma'],
            n_post=n_post, n_derived=n_der)
    except Exception as e:  # never let the correction break the run path
        return fb.Fallback(f'bundle integrity: {type(e).__name__}: {e}')


# ---------------------------------------------------------------------------
# Warning text (B7: numbers from manifest / diagnostics / weighted bundle)
# ---------------------------------------------------------------------------

def _ceil2(x: float) -> float:
    return math.ceil(float(x) * 100.0 - 1e-9) / 100.0


def _params_text(manifest, isc) -> str:
    g1 = (isc or {}).get('g1_passing_params') or manifest.get(
        'g1_passing_params')
    base = ['Re k₂', '|Im k₂|']
    if g1:
        base += list(g1)
    return ', '.join(base)


def corrected_pushforward(bundle):
    """Weighted/unweighted k2 pushforward summaries from the bundle."""
    from PlanetProfile.Inference.is_correction import weighted_quantile
    k2 = np.asarray(bundle.derived['k2_results'], float)
    w = np.asarray(bundle.weights, float)
    out = {}
    for name, v in (('Re_k2', k2[:, 0]), ('Im_k2', np.abs(k2[:, 1]))):
        q = weighted_quantile(v, w, [0.16, 0.5, 0.84])
        fin = v[np.isfinite(v)]
        out[name] = {'corr_p16': float(q[0]), 'corr_median': float(q[1]),
                     'corr_p84': float(q[2]),
                     'flow_median': float(np.median(fin))}
    return out


def corrected_warning(template: str, bundle, isc, x_obs, obs_sigma) -> str:
    m = bundle.manifest
    gv = m.get('gui_gate_values') or {}
    pf = corrected_pushforward(bundle)
    d = bundle.diagnostics
    return template.format(
        params=_params_text(m, isc),
        N=f"{int(m['n_samples']):,}",
        ess=f"{d['ess']:,.0f}",
        pareto_k=f"{d['pareto_k']:.2f}",
        gap=f"{_ceil2(gv['pushforward_median_gap_sigma_obs']):.2f}",
        ks=f"{_ceil2(gv['pushforward_ks']):.2f}",
        n_seeds=int(gv['n_seeds']),
        flow_im=f"{pf['Im_k2']['flow_median']:.3f}",
        corr_im=f"{pf['Im_k2']['corr_median']:.2f}",
        obs_im=f"{float(x_obs['Im_k2']):g}",
        obs_im_sig=f"{float(obs_sigma['Im_k2']):g}")


def fallback_warning(template: str, reason: str, manifest_or_none) -> str:
    gv = (manifest_or_none or {}).get('gui_gate_values') or {}
    frac = gv.get('inrange_prior_predictive_fail_frac')
    pct = f"{100 * float(frac):.0f}" if frac is not None else 'an appreciable'
    return template.format(reason=reason, fail_pct=pct)


# ---------------------------------------------------------------------------
# Result construction from the bundle
# ---------------------------------------------------------------------------

def _sha(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def build_corrected_result(bundle, config, runner, spec, repo_root,
                           t0: float):
    """InferenceResult built FROM the shipped bundle (never re-drawn)."""
    from PlanetProfile.Inference.inference_core import InferenceResult
    m, dv, d = bundle.manifest, bundle.derived, bundle.diagnostics
    n = len(bundle.weights)
    if list(bundle.param_names) != list(runner.param_names):
        raise ValueError('bundle parameter order differs from runner: '
                         f'{bundle.param_names} vs {runner.param_names}')
    manifest_path = _repo_path(repo_root, spec['is_correction']['manifest'])

    # Optional runner-derived provenance (structure radius etc.); best effort
    # because the correction itself never depends on it.
    R_body = Be = ind = None
    try:
        mr = runner._get_mcmc_runner()
        from PlanetProfile.Inference.mcmc_runner import _structure_R_body_km
        R_body = _structure_R_body_km(mr.structure_data)
        Be = mr._be_excitation_metadata()
        ind = mr._collect_posterior_Ae(bundle.samples)
    except Exception:
        pass
    train = {}
    try:
        art = runner._load_artifact_dict(
            spec['artifact_path'],
            validated_version_pairs=spec.get('validated_version_pairs', ()))
        train = {'density_estimator': art.get('density_estimator'),
                 'n_train_effective': art.get('n_train_effective'),
                 'train_seed': art.get('seed'),
                 'artifact_created_utc': art.get('created_utc')}
    except Exception:
        pass

    def g(k):
        return dv.get(k)

    # The loader strips the raw likelihood / flow-density arrays from
    # ``derived``; re-read them from the (sha-verified) npz.  Shipped
    # log_likelihoods are the true Gaussian log-likelihoods (-1e30 sentinel
    # on rejected draws, which carry weight 0).
    with np.load(manifest_path.parent / m['npz_file'],
                 allow_pickle=False) as z:
        ll = np.asarray(z['log_likelihoods'], float)
        flp = np.asarray(z['flow_log_prob'], float)

    has_ocean = g('has_ocean')
    md = {
        'sampler': 'sbi', 'mode': 'amortized',
        'density_estimator': train.get('density_estimator'),
        'elapsed_time_s': time.time() - t0,
        'heating_indices': None, 'derived_indices': np.arange(n),
        'artifact_path': str(spec['artifact_path']),
        'artifact_reused': True,
        'artifact_created_utc': train.get('artifact_created_utc'),
        'n_train_effective': train.get('n_train_effective'),
        'train_seed': train.get('train_seed'),
        'config_hash': m['config_hash'],
        'seed': m.get('effective_seed'),
        'x_obs': {k: float(v) for k, v in spec['x_obs'].items()},
        'imag_convention': m.get('imag_convention'),
        'R_body_km': R_body, 'induction_Ae': ind, 'Be_nT': Be,
        'flow_log_prob': flp,
        'truncate_bounds': None, 'kept_fraction': None,
        # ---- IS correction provenance (B1/B5/B6) ----
        'is_corrected': True,
        'is_weights': np.asarray(bundle.weights, float),
        'is_diagnostics': d,
        'is_manifest': str(manifest_path),
        'is_manifest_sha256': _sha(manifest_path),
        'is_npz_sha256': m.get('npz_sha256'),
        'is_bundle': m.get('bundle'),
        'has_ocean': (np.asarray(has_ocean, bool)
                      if has_ocean is not None else None),
    }
    ph = {k.split('_')[-1]: g(k) for k in
          ('D_icePhase_III', 'D_icePhase_V', 'D_icePhase_VI')}
    conv = {'n_samples': int(n), 'sampler': 'sbi_npe_is_corrected',
            'ess': float(d['ess']), 'acceptance_rate': None, 'r_hat': None,
            'pareto_k': float(d['pareto_k']),
            'w_max_frac': float(d['w_max_frac'])}
    return InferenceResult(
        config=config, samples=np.asarray(bundle.samples, float),
        log_likelihoods=ll,
        param_names=list(runner.param_names),
        param_labels=list(runner.param_labels),
        k2_results=g('k2_results'), h2_results=g('h2_results'),
        cmr2_results=g('cmr2_results'),
        D_ocean_results=g('D_ocean_results'),
        D_clath_results=g('D_clath_results'),
        D_iceIh_results=g('D_iceIh_results'),
        D_hsphere_results=g('D_hsphere_results'),
        D_iceHP_results=g('D_iceHP_results'),
        D_icePhase_results=ph,
        c20_results=g('c20_results'), c22_results=g('c22_results'),
        cmr2_hydro_results=g('cmr2_hydro_results'),
        heating_results=None, convergence_metrics=conv, metadata=md,
        weights=np.asarray(bundle.weights, float))


def weighted_csv(result) -> Optional[str]:
    """CSV of the shipped draws with a weights column + manifest sha."""
    md = result.metadata or {}
    w = result_weights(result)
    if w is None:
        return None
    lines = [f"# IS-corrected posterior draws; weights sum to 1; "
             f"manifest_sha256={md.get('is_manifest_sha256')}; "
             f"npz_sha256={md.get('is_npz_sha256')}",
             ','.join(list(result.param_names) + ['weight'])]
    S = np.asarray(result.samples, float)
    for row, wi in zip(S, w):
        lines.append(','.join(f"{v:.10g}" for v in row) + f",{wi:.10e}")
    return '\n'.join(lines) + '\n'
