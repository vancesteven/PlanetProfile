"""Pure (Streamlit-free) helpers behind the Inference page's MCMC config loader.

The page keeps the FULL loaded config dict as the single source of truth and
edits only a working copy of the user-editable blocks; ``build_run_config``
folds those edits back into the loaded dict so every non-edited block
(fixed_params, derived_params, param_groups, metadata, induction_bounds,
gravity_forward_model, arrhenius_params, ocean_overrides, bodyname, ...)
reaches the runner unchanged.
"""

from __future__ import annotations

import copy
import json
import math
import pickle
from collections import OrderedDict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from .inference_core import InferenceConfig
from .structure_cache import resolve_cache_path

CONFIGS_DIR = Path(__file__).resolve().parent / "configs"
REPO_ROOT = Path(__file__).resolve().parents[2]

GUI_STATUSES = ("production", "sensitivity", "smoke", "candidate", "retired")
PRODUCTION_STATUS = "production"

SAMPLER_EDIT_FIELDS = ("n_effective", "n_reeval", "random_state",
                       "checkpoint_interval", "max_iterations")
SAMPLER_FORM_DEFAULTS = {
    "n_effective": 500, "random_state": 42, "n_reeval": 500,
    "checkpoint_interval": 100, "max_iterations": 10000,
}

# prior_type -> the numeric fields the form edits for it. Semantics are the
# runner's (mcmc_runner._build_prior): log-uniform and uniform-in-linear
# bounds are log10 EXPONENTS; truncated_gaussian bounds are in linear units.
PRIOR_FIELDS = {
    "uniform": ("bounds",),
    "log-uniform": ("bounds",),
    "uniform-in-linear": ("bounds",),
    "normal": ("mean", "std"),
    "truncated_gaussian": ("mean", "std", "bounds"),
}

# Cache axis -> (axis key in the pickle, sampled-parameter name, transform
# from axis units to the parameter's units).
_GRID_AXES = (
    ("Tb_K_grid", "Tb_K", lambda v: v),
    ("wOcean_ppt_grid", "log10_wOcean_ppt", math.log10),
    ("zb_km_grid", "zb_km", lambda v: v),
)

BUILD_COMMAND = ("python -m PlanetProfile.Inference.build_phase_c1_cache "
                 "--config {config} --n-grid <N_TB>{wflag}")


# --------------------------------------------------------------------------
# Config catalogue / status
# --------------------------------------------------------------------------

def gui_info(cfg: Dict[str, Any]) -> Dict[str, Any]:
    """``metadata.gui`` as {status, note, runnable}; unlabeled -> candidate."""
    gui = (cfg.get("metadata") or {}).get("gui") or {}
    status = gui.get("status")
    if status not in GUI_STATUSES:
        status = "candidate"
        note = gui.get("note") or "No metadata.gui status recorded."
    else:
        note = gui.get("note", "")
    return {"status": status, "note": note,
            "runnable": gui.get("runnable", True) is not False}


def config_catalogue(configs_dir: Optional[Path] = None) -> List[Dict[str, Any]]:
    """One entry per config JSON: path, stem, body, status, note, runnable."""
    entries = []
    for fp in sorted(Path(configs_dir or CONFIGS_DIR).glob("*.json")):
        try:
            data = json.loads(fp.read_text(encoding="utf-8"))
        except Exception:
            continue
        body = (data.get("bodyname") or data.get("body")
                or data.get("planet_name") or data.get("Body") or "Unknown")
        entries.append({"path": str(fp), "stem": fp.stem, "body": body,
                        **gui_info(data)})
    return entries


# --------------------------------------------------------------------------
# Run-config assembly
# --------------------------------------------------------------------------

def _pair(v) -> tuple:
    return (float(v[0]), float(v[1]))


def observable_pair(v) -> Optional[tuple]:
    """``(value, sigma)`` floats if ``v`` is a 2-element numeric sequence
    (the only shape the form renders), else None."""
    if isinstance(v, (list, tuple)) and len(v) == 2:
        try:
            return _pair(v)
        except (TypeError, ValueError):
            return None
    return None


def merge_observables(loaded_obs: Dict[str, Any],
                      form_obs: Dict[str, Any]) -> Dict[str, Any]:
    """Loaded observables with the form's edits applied.

    Entries whose (value, sigma) match the loaded ones keep the loaded
    object verbatim (so the loaded key spelling, e.g. ``Im_k2``, and list
    type are preserved); edited ones keep the loaded key; keys the form
    added are appended as [value, sigma]. A loaded observable the form
    could not render (not a 2-element [value, sigma]) is passed through
    unchanged; only a rendered observable the user removed is dropped.
    """
    out: Dict[str, Any] = {}
    for key, val in form_obs.items():
        orig = key
        if key not in loaded_obs and key == "abs_Im_k2" and "Im_k2" in loaded_obs:
            orig = "Im_k2"
        if orig in loaded_obs:
            same = observable_pair(loaded_obs[orig]) == observable_pair(val) \
                and observable_pair(val) is not None
            out[orig] = loaded_obs[orig] if same else [val[0], val[1]]
        else:
            out[key] = [val[0], val[1]]
    for key, val in loaded_obs.items():
        if observable_pair(val) is None and key not in out:
            out[key] = val              # not rendered by the form: untouched
    # preserve the loaded ordering
    ordered = {k: out[k] for k in loaded_obs if k in out}
    ordered.update({k: v for k, v in out.items() if k not in ordered})
    return ordered


def build_run_config(loaded: Optional[Dict[str, Any]],
                     *,
                     param_space: Dict[str, Any],
                     observables: Dict[str, Any],
                     sampler_settings: Dict[str, Any],
                     structure_cache_path: Optional[str],
                     selected_params: Optional[Sequence[str]] = None,
                     fallback_bodyname: str = "Custom") -> InferenceConfig:
    """InferenceConfig for a run: the loaded dict plus the form's edits.

    With no edits the result equals ``InferenceConfig.from_dict(loaded)``.
    With ``loaded=None`` (no config file loaded) a bare config is built from
    the form fields alone (legacy preset/custom path).
    """
    if selected_params is not None:
        ps = {k: copy.deepcopy(param_space[k]) for k in selected_params
              if k in param_space}
    else:
        ps = copy.deepcopy(param_space)

    if loaded is None:
        return InferenceConfig(
            mode="mcmc", bodyname=fallback_bodyname, param_space=ps,
            observables=observables, sampler_settings=sampler_settings,
            structure_cache_path=structure_cache_path,
            random_state=sampler_settings["random_state"])

    cfg = copy.deepcopy(loaded)
    cfg["param_space"] = ps
    cfg["observables"] = merge_observables(loaded.get("observables", {}),
                                           observables)

    ss = dict(cfg.get("sampler_settings") or {})
    loaded_ss = loaded.get("sampler_settings") or {}
    for field in SAMPLER_EDIT_FIELDS:
        if field not in sampler_settings:
            continue
        if field == "random_state":
            base = loaded_ss.get("random_state", loaded.get("random_state",
                                 SAMPLER_FORM_DEFAULTS["random_state"]))
        else:
            base = loaded_ss.get(field, SAMPLER_FORM_DEFAULTS[field])
        if sampler_settings[field] != base:
            ss[field] = sampler_settings[field]
            if field == "random_state":
                cfg["random_state"] = sampler_settings[field]
    cfg["sampler_settings"] = ss

    if structure_cache_path and structure_cache_path != loaded.get(
            "structure_cache_path"):
        cfg["structure_cache_path"] = structure_cache_path
    return InferenceConfig.from_dict(cfg)


def run_config_from_state(state, fallback_bodyname: str = "Custom"
                          ) -> InferenceConfig:
    """``build_run_config`` fed from the page's session-state keys."""
    return build_run_config(
        state.get("inference_loaded_config"),
        param_space=state["inference_param_space"],
        observables=state["inference_observables"],
        sampler_settings=state["inference_sampler_settings"],
        structure_cache_path=state["inference_structure_cache_path"],
        selected_params=state["inference_selected_params"],
        fallback_bodyname=fallback_bodyname)


# --------------------------------------------------------------------------
# Prior editing helpers
# --------------------------------------------------------------------------

# Prior types whose bounds are log10 EXPONENTS (see mcmc_runner._build_prior).
EXPONENT_TYPES = ("log-uniform", "uniform-in-linear")
_RT_REL, _RT_ABS = 1e-9, 1e-12


def _close(a, b) -> bool:
    return abs(a - b) <= _RT_ABS + _RT_REL * max(abs(a), abs(b))


def retype_prior_entry(param_def, prior_type: str,
                       from_entry: Optional[Dict[str, Any]] = None,
                       origin: Optional[Dict[str, Any]] = None):
    """Prior dict of ``prior_type`` for a newly added/retyped parameter.

    Returns ``(entry, notice)``; ``notice`` is None or a sentence for the UI.

    Conversion preserves the SAMPLED RANGE of the parameter exactly as
    ``mcmc_runner._build_prior`` builds it:

    - ``uniform`` / ``truncated_gaussian``: bounds are linear values.
    - ``log-uniform``: ``loguniform(10**lo, 10**hi)``, so the sampled range
      is ``[10**lo, 10**hi]``; uniform<->log-uniform converts with
      ``log10`` / ``10**``.
    - ``uniform-in-linear``: the sampled variable is ``x = log10(v)`` and its
      bounds ``[lo, hi]`` are the range of x itself, so to/from the linear
      types the bounds are unchanged; to/from ``log-uniform`` (both are
      exponent-valued) they are also unchanged.

    A linear -> log-uniform conversion needs both bounds > 0; otherwise the
    registry default for the new type is used and the notice says so.
    ``origin`` (the entry before the first retype, if any) is restored
    exactly when converting back to its type reproduces it within float
    tolerance, so a round trip is not an edit.
    """
    entry: Dict[str, Any] = {"prior_type": prior_type}
    fields = PRIOR_FIELDS[prior_type]
    notice = None
    src_type = (from_entry or {}).get("prior_type")
    linear_types = ("uniform", "truncated_gaussian")
    pname = getattr(param_def, "id", None) or getattr(param_def, "label", "")
    bounds = None
    if from_entry and "bounds" in from_entry:
        lo, hi = float(from_entry["bounds"][0]), float(from_entry["bounds"][1])
        if src_type == "log-uniform" and prior_type in linear_types:
            try:
                bounds = [10.0 ** lo, 10.0 ** hi]
            except OverflowError:
                notice = (f"`{pname}`: the exponent bounds are too large to "
                          "convert to linear values; the registry default "
                          f"for {prior_type} was used instead.")
        elif src_type in linear_types and prior_type == "log-uniform":
            if lo > 0 and hi > 0:
                bounds = [math.log10(lo), math.log10(hi)]
            else:
                notice = (f"`{pname}` has a non-positive bound, which a "
                          "log-uniform prior cannot represent; the registry "
                          "default for log-uniform was used instead.")
        else:
            bounds = [lo, hi]   # same meaning (linear<->linear, exponent<->exponent, uniform-in-linear<->linear)
    if bounds is None:
        bounds = list(param_def.default_bounds)
        if prior_type == "log-uniform":
            if bounds[0] > 0 and bounds[1] > 0:
                bounds = [math.log10(bounds[0]), math.log10(bounds[1])]
            else:
                bounds = [-1.0, 1.0]
                notice = ((notice or "") + " The registry default bounds are "
                          "non-positive too, so placeholder exponents "
                          "[-1, 1] (i.e. 0.1 to 10) are shown; set the "
                          "bounds manually.").strip()
    if "bounds" in fields:
        entry["bounds"] = [float(bounds[0]), float(bounds[1])]
    if "mean" in fields:
        mean = param_def.default_mean
        if mean is None:
            mean = 0.5 * (bounds[0] + bounds[1])
        std = param_def.default_std
        if not std:
            std = 0.25 * abs(bounds[1] - bounds[0]) or 1.0
        entry["mean"] = float(mean)
        entry["std"] = float(std)
    if origin and origin.get("prior_type") == prior_type and notice is None:
        keys = [k for k in ("bounds", "mean", "std") if k in entry]
        try:
            same = all(
                _close(float(x), float(y))
                for k in keys for x, y in zip(
                    entry[k] if k == "bounds" else [entry[k]],
                    origin[k] if k == "bounds" else [origin[k]]))
        except (KeyError, TypeError, ValueError):
            same = False
        if same:
            return copy.deepcopy(origin), None
    return entry, notice


def default_prior_entry(param_def, prior_type: str,
                        from_entry: Optional[Dict[str, Any]] = None
                        ) -> Dict[str, Any]:
    """``retype_prior_entry`` without the notice."""
    return retype_prior_entry(param_def, prior_type, from_entry)[0]


def prior_number_format(*vals) -> str:
    """printf format (for st.number_input) that shows every given value
    without losing digits: fixed notation with at least 4 decimals and as
    many more as the values need (up to 12), or scientific with enough
    mantissa digits for tiny/huge magnitudes."""
    mags = []
    for v in vals:
        try:
            f = float(v)
        except (TypeError, ValueError):
            continue
        if math.isfinite(f) and f != 0.0:
            mags.append(abs(f))
    sci = bool(mags) and (max(mags) < 1e-2 or max(mags) >= 1e5)
    n = 4 if not sci else 3
    for v in vals:
        try:
            f = float(v)
        except (TypeError, ValueError):
            continue
        if not math.isfinite(f) or f == 0.0:
            continue
        if sci:
            mant = format(f, ".12e").split("e")[0].rstrip("0").rstrip(".")
            digits = len(mant.split(".")[1]) if "." in mant else 0
            n = max(n, min(digits, 12))
        else:
            dec = format(f, ".12f").rstrip("0")
            digits = len(dec.split(".")[1]) if "." in dec else 0
            n = max(n, min(digits, 12))
    return f"%.{n}{'e' if sci else 'f'}"


# --------------------------------------------------------------------------
# Structure-grid check
# --------------------------------------------------------------------------

# Per-node 1-D profile arrays that must share one length in a valid node.
_NODE_PROFILE_KEYS = ("r_m", "rho", "K_Pa", "mu_Pa", "T_K", "P_MPa")

_AXES_CACHE: "OrderedDict[tuple, Dict[str, Any]]" = OrderedDict()
_AXES_CACHE_MAX = 8


def _axis_info(grid, fn, axis_key):
    vals = [float(v) for v in grid]
    zero = any(v <= 0.0 for v in vals) and fn is math.log10
    if fn is math.log10:
        vals = [v for v in vals if v > 0.0]
    if not vals:
        return None
    info = {"axis": axis_key, "range": (fn(min(vals)), fn(max(vals))),
            "n": len(grid)}
    if zero:
        info["zero_node"] = True     # w = 0 has no log10; excluded from range
    return info


def _load_axes(path: str) -> Dict[str, Any]:
    with open(path, "rb") as fh:
        data = pickle.load(fh)
    axes: Dict[str, Any] = {}
    quality = {"n_nodes": 0, "n_none": 0, "n_inconsistent": 0,
               "inconsistent_examples": []}
    if isinstance(data, dict):
        for axis_key, param, fn in _GRID_AXES:
            grid = data.get(axis_key)
            if grid is None or len(grid) == 0:
                continue
            info = _axis_info(grid, fn, axis_key)
            if info is not None:
                axes[param] = info
        # The frozen branch's zb axis (zb > ocean-branch maximum) is its own
        # axis: a sampled zb_km prior is range-checked against the ocean
        # axis only.
        fz = data.get("frozen_zb_km_grid")
        if fz is not None and len(fz):
            info = _axis_info(fz, lambda v: v, "frozen_zb_km_grid")
            if info is not None:
                axes["frozen_zb_km"] = info
        nodes = data.get("structures")
        if nodes is not None:
            try:
                quality["n_nodes"] = len(nodes)
                for i, node in enumerate(nodes):
                    if node is None:
                        quality["n_none"] += 1
                    elif isinstance(node, dict):
                        lens = {len(node[k]) for k in _NODE_PROFILE_KEYS
                                if k in node and hasattr(node[k], "__len__")}
                        if len(lens) > 1:
                            quality["n_inconsistent"] += 1
                            if len(quality["inconsistent_examples"]) < 3:
                                quality["inconsistent_examples"].append(i)
            except TypeError:
                pass
    return {"axes": axes, "quality": quality}


def _cache_key(path: Path) -> tuple:
    st = path.stat()
    return (str(path), st.st_mtime, st.st_size)


def axes_cached(path: Path) -> bool:
    """True when ``cache_axes(path)`` will not re-read the pickle."""
    try:
        return _cache_key(path) in _AXES_CACHE
    except OSError:
        return False


def cache_report(path: Path) -> Dict[str, Any]:
    """``{'axes', 'quality'}`` for the cache pickle, memoized on (resolved
    path, mtime, size): Streamlit reruns never re-unpickle a large cache."""
    key = _cache_key(path)
    hit = _AXES_CACHE.get(key)
    if hit is None:
        hit = _load_axes(str(path))
        _AXES_CACHE[key] = hit
        while len(_AXES_CACHE) > _AXES_CACHE_MAX:
            _AXES_CACHE.popitem(last=False)
    else:
        _AXES_CACHE.move_to_end(key)
    return hit


def cache_axes(path: Path) -> Dict[str, Any]:
    return cache_report(path)["axes"]


def _prior_extent(spec: Dict[str, Any]) -> Optional[tuple]:
    ptype = spec.get("prior_type")
    if ptype in ("uniform", "uniform-in-linear", "truncated_gaussian"):
        b = spec.get("bounds")
        return (float(b[0]), float(b[1])) if b else None
    if ptype == "normal":
        m, s = spec.get("mean"), spec.get("std")
        return (m - 3 * s, m + 3 * s) if m is not None and s else None
    return None


def resolve_check_path(cache_path: str) -> Path:
    """The file ``check_structure_cache`` will read for ``cache_path``."""
    resolved = resolve_cache_path(cache_path)
    if not resolved.is_absolute() and not resolved.exists():
        resolved = REPO_ROOT / resolved
    return resolved


def check_structure_cache(param_space: Dict[str, Any],
                          fixed_params: Optional[Dict[str, Any]],
                          cache_path: Optional[str],
                          config_path: Optional[str] = None
                          ) -> Dict[str, Any]:
    """Real structure-grid validation.

    Returns ``{'status': 'unset'|'missing'|'unreadable'|'ok', 'resolved',
    'axes', 'quality', 'issues', 'command'}``. ``quality`` counts None and
    internally inconsistent nodes (profile arrays of unequal length). ``issues`` lists each sampled/fixed
    parameter whose prior (or fixed value) extends outside the cached axis,
    as ``{'param', 'prior_range', 'cache_range'}``.
    """
    if not cache_path:
        return {"status": "unset", "resolved": None, "axes": {}, "issues": [],
                "command": None}
    resolved = resolve_check_path(cache_path)
    needs_w = "log10_wOcean_ppt" in (param_space or {})
    cfg_ref = config_path or "<config.json>"
    command = BUILD_COMMAND.format(
        config=cfg_ref, wflag=" --n-wgrid <N_W>" if needs_w else "")
    out = {"status": "missing", "resolved": str(resolved), "axes": {},
           "issues": [], "command": command}
    if not resolved.exists():
        return out
    try:
        rep = cache_report(resolved)
        axes = rep["axes"]
    except Exception as exc:  # corrupt / unreadable pickle
        out.update(status="unreadable", error=str(exc))
        return out
    out["status"] = "ok"
    out["axes"] = axes
    out["quality"] = rep["quality"]
    extents = {k: _prior_extent(v) for k, v in (param_space or {}).items()}
    for k, v in (fixed_params or {}).items():
        try:
            extents[k] = (float(v), float(v))
        except (TypeError, ValueError):
            pass
    for param, info in axes.items():
        if param == "frozen_zb_km":     # reported separately, never checked
            continue
        ext = extents.get(param)
        if ext is None:
            continue
        lo, hi = info["range"]
        tol = 1e-6 * max(abs(hi - lo), 1e-12)
        if ext[0] < lo - tol or ext[1] > hi + tol:
            out["issues"].append({"param": param, "prior_range": ext,
                                  "cache_range": (lo, hi)})
    return out
