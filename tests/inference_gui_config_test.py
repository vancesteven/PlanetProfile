"""Inference-page MCMC config loading: priors survive, run config is the
loaded config, structure-grid messaging is real, configs carry a GUI status.

Defects covered (audit 2026-11): loaded priors overwritten by registry
defaults on the first rerun; Run built a bare InferenceConfig (dropping
fixed_params, derived_params, metadata, ...); a false "Structure Grid
Required" banner pointing at nonexistent scripts; unlabeled smoke/retired/
non-runnable configs.
"""

import copy
import dataclasses
import json
import math
import pickle
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "PlanetProfileApp"))
PAGE = REPO / "PlanetProfileApp/pages/Inference.py"
CONFIGS = REPO / "PlanetProfile/Inference/configs"
CONFIG_FILES = sorted(CONFIGS.glob("*.json"))
STEMS = [p.stem for p in CONFIG_FILES]

from PlanetProfile.Inference import gui_config  # noqa: E402
from PlanetProfile.Inference.inference_core import InferenceConfig  # noqa: E402

pytest.importorskip("streamlit")


def _file(stem):
    return json.loads((CONFIGS / f"{stem}.json").read_text(encoding="utf-8"))


def _diff_blocks(a, b):
    """Names of InferenceConfig fields that differ."""
    da, db = dataclasses.asdict(a), dataclasses.asdict(b)
    return sorted(k for k in da if da[k] != db[k])


# ---------------------------------------------------------------------------
# metadata.gui and hashing (no AppTest)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("stem", STEMS)
def test_every_config_has_gui_status(stem):
    gui = _file(stem)["metadata"]["gui"]
    assert gui["status"] in gui_config.GUI_STATUSES
    assert gui["note"].strip()


def test_non_runnable_configs_are_marked_with_reason():
    marked = {s for s in STEMS if _file(s)["metadata"]["gui"].get("runnable") is False}
    assert marked == {"callisto_mgso4_andrade_8D", "enceladus_cassini_isostasy_7D"}
    for s in marked:
        assert len(_file(s)["metadata"]["gui"]["note"]) > 30


@pytest.mark.parametrize("stem", STEMS)
def test_metadata_gui_does_not_change_hash(stem):
    """metadata.gui is display-only: the config hash must not depend on it
    (present == absent == a different gui block). No git, no skip."""
    data = _file(stem)
    assert "gui" in data["metadata"]
    h_with = InferenceConfig.from_dict(data).generate_hash()
    stripped = copy.deepcopy(data)
    stripped["metadata"].pop("gui")
    assert h_with == InferenceConfig.from_dict(stripped).generate_hash()
    other = copy.deepcopy(data)
    other["metadata"]["gui"] = {"status": "retired", "runnable": False,
                                "note": "a different gui block"}
    assert h_with == InferenceConfig.from_dict(other).generate_hash()
    # and the hash is not trivially constant: a real field change moves it
    changed = copy.deepcopy(data)
    changed["random_state"] = int(changed.get("random_state", 42)) + 1
    assert h_with != InferenceConfig.from_dict(changed).generate_hash()


def test_catalogue_matches_files():
    cat = gui_config.config_catalogue()
    assert [e["stem"] for e in cat] == STEMS
    assert {e["status"] for e in cat} <= set(gui_config.GUI_STATUSES)
    assert any(e["status"] == "production" for e in cat)


# ---------------------------------------------------------------------------
# build_run_config (pure)
# ---------------------------------------------------------------------------

def _state_from_file(data):
    obs = {}
    for k, v in data["observables"].items():
        obs["abs_Im_k2" if k == "Im_k2" else k] = (float(v[0]), float(v[1]))
    ss = dict(gui_config.SAMPLER_FORM_DEFAULTS)
    for f in gui_config.SAMPLER_EDIT_FIELDS:
        if f in data.get("sampler_settings", {}):
            ss[f] = data["sampler_settings"][f]
    if "random_state" not in data.get("sampler_settings", {}):
        ss["random_state"] = data.get("random_state", 42)
    return dict(
        inference_loaded_config=data,
        inference_param_space=copy.deepcopy(data["param_space"]),
        inference_observables=obs,
        inference_sampler_settings=ss,
        inference_structure_cache_path=data["structure_cache_path"],
        inference_selected_params=list(data["param_space"]),
    )


@pytest.mark.parametrize("stem", STEMS)
def test_unedited_form_state_reproduces_file_config(stem):
    data = _file(stem)
    built = gui_config.run_config_from_state(_state_from_file(data))
    ref = InferenceConfig.from_dict(data)
    assert _diff_blocks(built, ref) == []
    assert built.generate_hash() == ref.generate_hash()
    assert built.metadata == ref.metadata


def test_run_config_keeps_non_form_blocks():
    data = _file("test50_titan_noocean_andrade_5D")
    built = gui_config.run_config_from_state(_state_from_file(data))
    assert built.fixed_params == {"Tb_K": 250.965}
    assert built.bodyname == "Titan"
    data = _file("europa_clipper_v5_geodesy_11D")
    built = gui_config.run_config_from_state(_state_from_file(data))
    assert built.derived_params and built.induction_bounds
    assert built.gravity_forward_model == data["gravity_forward_model"]
    assert built.param_space["D_iceIh_km"] == {"prior_type": "uniform",
                                               "bounds": [5.0, 80.0]}


def test_no_loaded_config_keeps_bare_path():
    data = _file("test50_titan_noocean_andrade_8D")
    st = _state_from_file(data)
    st["inference_loaded_config"] = None
    cfg = gui_config.run_config_from_state(st, fallback_bodyname="Titan")
    assert cfg.bodyname == "Titan" and cfg.metadata == {}
    assert cfg.param_space == data["param_space"]


# ---------------------------------------------------------------------------
# Structure-grid check (pure)
# ---------------------------------------------------------------------------

@pytest.fixture
def grid_cache(tmp_path):
    p = tmp_path / "grid.pkl"
    p.write_bytes(pickle.dumps({
        "Tb_K_grid": [250.0, 255.0, 260.0],
        "wOcean_ppt_grid": [0.1, 1.0, 10.0, 100.0],
        "structures": [],
    }))
    return p


def test_grid_check_in_range_no_issue(grid_cache):
    rep = gui_config.check_structure_cache(
        {"Tb_K": {"prior_type": "uniform", "bounds": [250.0, 260.0]},
         "log10_wOcean_ppt": {"prior_type": "uniform", "bounds": [-1.0, 2.0]}},
        {}, str(grid_cache))
    assert rep["status"] == "ok" and rep["issues"] == []
    assert rep["axes"]["Tb_K"]["range"] == (250.0, 260.0)
    assert rep["axes"]["log10_wOcean_ppt"]["range"] == pytest.approx((-1.0, 2.0))


def test_grid_check_outside_range_names_both_ranges(grid_cache):
    rep = gui_config.check_structure_cache(
        {"Tb_K": {"prior_type": "uniform", "bounds": [248.0, 262.0]}},
        {}, str(grid_cache))
    (issue,) = rep["issues"]
    assert issue["param"] == "Tb_K"
    assert issue["prior_range"] == (248.0, 262.0)
    assert issue["cache_range"] == (250.0, 260.0)


def test_grid_check_fixed_param_outside(grid_cache):
    rep = gui_config.check_structure_cache(
        {"alpha": {"prior_type": "uniform", "bounds": [0.1, 0.5]}},
        {"Tb_K": 270.0}, str(grid_cache))
    assert [i["param"] for i in rep["issues"]] == ["Tb_K"]


def test_grid_check_missing_names_path_and_real_command(tmp_path):
    missing = tmp_path / "nope.pkl"
    rep = gui_config.check_structure_cache(
        {"Tb_K": {"prior_type": "uniform", "bounds": [1, 2]}}, {},
        str(missing), "PlanetProfile/Inference/configs/x.json")
    assert rep["status"] == "missing" and rep["resolved"] == str(missing)
    assert "build_phase_c1_cache --config PlanetProfile/Inference/configs/x.json" \
        in rep["command"]
    assert "scripts/prepare" not in rep["command"]


def test_registry_no_longer_emits_false_grid_warning():
    from PlanetProfile.Inference.parameter_registry import (
        validate_parameter_combination)
    r = validate_parameter_combination(["Tb_K", "alpha", "log10_zeta"])
    assert r["requires_rebuild"] is True
    assert not any("prepare_structure" in w or "structure grid" in w
                   for w in r["warnings"])


def test_page_has_no_nonexistent_script_paths():
    src = PAGE.read_text(encoding="utf-8")
    assert "scripts/prepare_structure" not in src
    assert "Structure Grid Required" not in src


# ---------------------------------------------------------------------------
# AppTest: load via the selector, rerun, compare
# ---------------------------------------------------------------------------

def _app():
    from streamlit.testing.v1 import AppTest
    return AppTest.from_file(str(PAGE), default_timeout=300)


def _w(at, kind, name):
    """Form widget by base name: loader-reset widgets carry a load-generation
    suffix (``name__g<N>``) so each config load creates new widget identities."""
    gen = at.session_state["inference_load_gen"] \
        if "inference_load_gen" in at.session_state else 0
    return getattr(at, kind)(key=name if not gen else f"{name}__g{gen}")


def _ok(at):
    assert not at.exception, [e.value for e in at.exception]


def _load_via_selector(at, stem, reruns=2, data=None):
    at.run()
    _ok(at)
    at.checkbox(key="cfg_selector_show_all").check()
    at.run()
    body = (data or _file(stem))["bodyname"]
    at.selectbox(key="cfg_selector_body").select(body)
    at.run()
    sel = at.selectbox(key="cfg_selector_variant")
    label = next(o for o in sel.options if str(o).startswith(f"{stem}  ["))
    sel.select(label)
    at.run()
    at.button(key="cfg_selector_apply").click()
    at.run()
    _ok(at)
    for _ in range(reruns):
        at.run()
        _ok(at)


def _state(at):
    keys = ("inference_loaded_config", "inference_param_space",
            "inference_observables", "inference_sampler_settings",
            "inference_structure_cache_path", "inference_selected_params")
    return {k: at.session_state[k] for k in keys}


def _norm(x):
    return json.loads(json.dumps(x))


def _texts(at):
    out = []
    for kind in ("error", "warning", "success", "info", "caption", "markdown"):
        out += [str(e.value) for e in getattr(at, kind)]
    return "\n".join(out)


@pytest.mark.gui
@pytest.mark.parametrize("stem", STEMS)
def test_loaded_priors_survive_reruns_and_run_config_matches_file(stem):
    data = _file(stem)
    at = _app()
    _load_via_selector(at, stem)
    st = _state(at)
    assert _norm(st["inference_param_space"]) == _norm(data["param_space"])
    built = gui_config.run_config_from_state(st)
    ref = InferenceConfig.from_dict(data)
    assert _diff_blocks(built, ref) == []
    assert built.generate_hash() == ref.generate_hash()
    assert built.metadata == ref.metadata
    assert built.bodyname == data["bodyname"]


@pytest.mark.gui
def test_v5_D_iceIh_stays_uniform_5_80():
    at = _app()
    _load_via_selector(at, "europa_clipper_v5_geodesy_11D", reruns=3)
    assert at.session_state["inference_param_space"]["D_iceIh_km"] == {
        "prior_type": "uniform", "bounds": [5, 80]}
    assert _w(at, "selectbox", "prior_type_D_iceIh_km").value == "uniform"


@pytest.mark.gui
def test_truncated_gaussian_prior_survives_and_is_editable():
    stem = "enceladus_cassini_isostasy_ocean6D_rhorock"
    data = _file(stem)
    at = _app()
    _load_via_selector(at, stem)
    tg = {k: v for k, v in data["param_space"].items()
          if v["prior_type"] == "truncated_gaussian"}
    assert tg
    for k, v in tg.items():
        assert _w(at, "selectbox", f"prior_type_{k}").value == "truncated_gaussian"
        assert at.session_state["inference_param_space"][k] == v


@pytest.mark.gui
def test_single_bound_edit_changes_only_that_field():
    stem = "europa_clipper_v5_geodesy_11D"
    data = _file(stem)
    at = _app()
    _load_via_selector(at, stem)
    _w(at, "number_input", "lower_D_iceIh_km_uniform").set_value(10.0)
    at.run()
    _ok(at)
    built = gui_config.run_config_from_state(_state(at))
    ref = InferenceConfig.from_dict(data)
    assert _diff_blocks(built, ref) == ["param_space"]
    for k in ref.param_space:
        if k == "D_iceIh_km":
            assert built.param_space[k] == {"prior_type": "uniform",
                                            "bounds": [10.0, 80]}
        else:
            assert built.param_space[k] == ref.param_space[k]
    assert built.generate_hash() != ref.generate_hash()


@pytest.mark.gui
def test_observable_and_sampler_edits_are_isolated():
    stem = "europa_clipper_v4_geodesy_11D"
    data = _file(stem)
    at = _app()
    _load_via_selector(at, stem)
    _w(at, "number_input", "Re_k2_value").set_value(0.31)
    at.run()
    _ok(at)
    built = gui_config.run_config_from_state(_state(at))
    ref = InferenceConfig.from_dict(data)
    assert _diff_blocks(built, ref) == ["observables"]
    for k in ref.observables:
        if k == "Re_k2":
            assert built.observables[k][0] == 0.31
            assert built.observables[k][1] == ref.observables[k][1]
        else:
            assert built.observables[k] == ref.observables[k]
    _w(at, "slider", "n_eff_slider").set_value(900)
    at.run()
    _ok(at)
    built = gui_config.run_config_from_state(_state(at))
    assert _diff_blocks(built, ref) == ["observables", "sampler_settings"]
    assert built.sampler_settings["n_effective"] == 900
    assert {k: v for k, v in built.sampler_settings.items()
            if k != "n_effective"} == {
        k: v for k, v in ref.sampler_settings.items() if k != "n_effective"}


@pytest.mark.gui
def test_grid_messaging_in_range_no_warning():
    at = _app()
    _load_via_selector(at, "europa_clipper_v4_geodesy_11D")
    txt = _texts(at)
    assert "Cache file found" in txt and "Tb_K 259.5 to 271" in txt
    assert "extends outside the cached grid" not in txt
    assert "Structure Grid Required" not in txt
    assert "scripts/prepare" not in txt


@pytest.mark.gui
def test_grid_messaging_outside_range_warns_with_both_ranges():
    at = _app()
    _load_via_selector(at, "europa_clipper_v4_geodesy_11D")
    _w(at, "number_input", "upper_Tb_K_uniform").set_value(275.0)
    at.run()
    _ok(at)
    warns = [str(w.value) for w in at.warning
             if "extends outside the cached grid" in str(w.value)]
    assert len(warns) == 1 and "Tb_K" in warns[0]
    assert "259.5" in warns[0] and "271" in warns[0] and "275" in warns[0]


@pytest.mark.gui
def test_grid_messaging_missing_cache_single_error_with_command():
    stem = "europa_clipper_v4_geodesy_11D"
    at = _app()
    _load_via_selector(at, stem)
    bogus = "PlanetProfile/Test/mcmc_results/NoSuch/missing_grid.pkl"
    _w(at, "text_input", "cache_path_input").set_value(bogus)
    at.run()
    _ok(at)
    errs = [str(e.value) for e in at.error if "cache not found" in str(e.value)]
    assert len(errs) == 1
    assert bogus in errs[0]
    assert f"build_phase_c1_cache --config PlanetProfile/Inference/configs/{stem}.json" \
        in errs[0]
    assert "--n-wgrid" in errs[0]
    assert "scripts/prepare" not in _texts(at)
    assert "Cache file found" not in _texts(at)


@pytest.mark.gui
def test_selector_defaults_to_production_and_toggle_reveals_rest():
    at = _app()
    at.run()
    _ok(at)
    shown = []
    for body in at.selectbox(key="cfg_selector_body").options:
        at.selectbox(key="cfg_selector_body").select(body)
        at.run()
        shown += list(at.selectbox(key="cfg_selector_variant").options)
    assert shown and all("[production" in s for s in shown)
    assert not any("smoke" in s or "retired" in s or "candidate" in s for s in shown)
    at.checkbox(key="cfg_selector_show_all").check()
    at.run()
    all_labels = []
    for body in at.selectbox(key="cfg_selector_body").options:
        at.selectbox(key="cfg_selector_body").select(body)
        at.run()
        all_labels += list(at.selectbox(key="cfg_selector_variant").options)
    assert len(all_labels) == len(STEMS)
    assert any("[retired · cannot run]" in s for s in all_labels)


@pytest.mark.gui
@pytest.mark.parametrize("stem", ["enceladus_cassini_isostasy_7D",
                                  "callisto_mgso4_andrade_8D"])
def test_non_runnable_config_run_disabled_with_reason(stem):
    at = _app()
    _load_via_selector(at, stem)
    note = _file(stem)["metadata"]["gui"]["note"]
    assert any(note in str(e.value) for e in at.error)
    run_btn = next(b for b in at.button if "Run MCMC" in str(b.label))
    assert run_btn.disabled


# ---------------------------------------------------------------------------
# Headless: every config constructs a runner, or is marked non-runnable
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def cache_copies(tmp_path_factory):
    """Copy of each cache (+ sidecars) so runner construction can never
    rewrite a tracked file under PlanetProfile/Test/."""
    root = tmp_path_factory.mktemp("cache_copies")
    done = {}

    def get(rel):
        if rel not in done:
            src = REPO / rel
            dst_dir = root / f"c{len(done)}"
            dst_dir.mkdir()
            stem = src.name[:-len(".pkl")]
            for f in src.parent.glob(f"{stem}.*"):
                if f.is_file():
                    shutil.copy2(f, dst_dir / f.name)
            done[rel] = dst_dir / src.name
        return str(done[rel])
    return get


@pytest.mark.parametrize("stem", STEMS)
def test_config_constructs_runner_or_is_marked_non_runnable(stem, cache_copies):
    data = _file(stem)
    gui = data["metadata"]["gui"]
    if gui.get("runnable") is False:
        assert gui["note"].strip()
        return
    pytest.importorskip("pocomc")
    from PlanetProfile.Inference.mcmc_runner import MCMCRunner
    data["structure_cache_path"] = cache_copies(data["structure_cache_path"])
    runner = MCMCRunner(InferenceConfig.from_dict(data))
    assert runner.param_names == list(data["param_space"])


# ---------------------------------------------------------------------------
# Review fixes (2026-11-17): prior retype, k2 channels, non-pair observables,
# carry-over between configs, cache-check axes/zero node/quality/memoization
# ---------------------------------------------------------------------------

from PlanetProfile.Inference.parameter_registry import (  # noqa: E402
    PARAMETER_REGISTRY as _REG)

_BASE = "test50_titan_noocean_andrade_5D"


@pytest.fixture
def fx(tmp_path, monkeypatch):
    """Fixture configs in a tmp configs dir that the selector lists."""
    base = _file(_BASE)
    out = {}

    def add(stem, mut):
        d = copy.deepcopy(base)
        d["metadata"]["gui"] = {"status": "candidate", "note": "fixture"}
        mut(d)
        (tmp_path / f"{stem}.json").write_text(json.dumps(d))
        out[stem] = d

    def full(d):
        d["observables"] = {"Re_k2": [0.608, 0.048], "Im_k2": [0.135, 0.035],
                            "CMR2": [0.3547, 0.0024]}
        d["sampler_settings"] = {"n_effective": 77, "n_reeval": 33,
                                 "random_state": 7}
        d["random_state"] = 7

    def re_only(d):
        d["observables"] = {"Re_k2": [0.608, 0.048]}

    def im_only(d):
        d["observables"] = {"Im_k2": [0.135, 0.035]}

    def obs3(d):
        d["observables"] = {"Re_k2": [0.608, 0.048],
                            "Re_h2": [1.2, 0.1, 0.2]}

    def bare(d):
        d["observables"] = {"Re_k2": [0.5, 0.05]}
        d["sampler_settings"] = {}
        d.pop("structure_cache_path", None)

    def prec(d):
        d["param_space"]["log10_zeta"]["bounds"] = [-3.0, 1.845098]

    def posneg(d):
        d["param_space"]["alpha"]["bounds"] = [-0.1, 0.45]

    def noobs(d):
        d["observables"] = {}

    for stem, mut in (("fx_noobs", noobs), ("fx_full", full), ("fx_re_only", re_only),
                      ("fx_im_only", im_only), ("fx_obs3", obs3),
                      ("fx_bare", bare), ("fx_prec", prec),
                      ("fx_posneg", posneg)):
        add(stem, mut)
    monkeypatch.setattr(gui_config, "CONFIGS_DIR", tmp_path)
    return out


def _built_diff(at, data):
    built = gui_config.run_config_from_state(_state(at))
    return _diff_blocks(built, InferenceConfig.from_dict(data)), built


# -- finding 2: prior retype -------------------------------------------------

def _sampled_range(entry):
    """Range of the SAMPLED variable exactly as mcmc_runner._build_prior
    builds it (scipy support), for the non-normal types."""
    from scipy.stats import loguniform, truncnorm, uniform
    t, b = entry["prior_type"], entry.get("bounds")
    if t == "uniform":
        d = uniform(loc=b[0], scale=b[1] - b[0])
    elif t == "log-uniform":
        d = loguniform(a=10 ** b[0], b=10 ** b[1])
    elif t == "uniform-in-linear":
        return (b[0], b[1])               # sampled x = log10(v) in [lo, hi]
    elif t == "truncated_gaussian":
        return (b[0], b[1])
    return tuple(d.support())


def test_retype_log_uniform_to_linear_types_preserves_sampled_range():
    pd = _REG["Tb_K"]
    lu = {"prior_type": "log-uniform", "bounds": [2.414, 2.433]}
    want = _sampled_range(lu)
    for t in ("uniform", "truncated_gaussian"):
        got = gui_config.default_prior_entry(pd, t, from_entry=lu)
        assert _sampled_range(got) == pytest.approx(want, rel=1e-12)
        assert got["bounds"][0] > 200         # linear K, not exponents
    # exponent-valued <-> exponent-valued and uniform-in-linear (sampled x
    # is itself log10 v): bounds are unchanged
    uil = {"prior_type": "uniform-in-linear", "bounds": [-1.0, 2.0]}
    assert gui_config.default_prior_entry(
        _REG["log10_wOcean_ppt"], "uniform", uil)["bounds"] == [-1.0, 2.0]
    assert gui_config.default_prior_entry(
        _REG["log10_wOcean_ppt"], "uniform-in-linear",
        {"prior_type": "uniform", "bounds": [-1.0, 2.0]})["bounds"] == [-1.0, 2.0]


def test_retype_uniform_to_log_uniform_preserves_sampled_range():
    pd = _REG["Tb_K"]
    u = {"prior_type": "uniform", "bounds": [259.5, 271.0]}
    entry, notice = gui_config.retype_prior_entry(pd, "log-uniform", u)
    assert notice is None
    assert _sampled_range(entry) == pytest.approx((259.5, 271.0), rel=1e-12)


def test_retype_non_positive_falls_back_to_registry_default_with_notice():
    # registry default for alpha is positive: used (as exponents)
    pd = _REG["alpha"]
    entry, notice = gui_config.retype_prior_entry(
        pd, "log-uniform", {"prior_type": "uniform", "bounds": [-0.1, 0.45]})
    assert notice and "non-positive" in notice
    lo, hi = pd.default_bounds
    assert entry["bounds"] == pytest.approx([math.log10(lo), math.log10(hi)])
    # registry default also non-positive: explicit placeholder + notice
    entry, notice = gui_config.retype_prior_entry(
        _REG["dC20_nh"], "log-uniform",
        {"prior_type": "uniform", "bounds": [-2e-5, 2e-5]})
    assert entry["bounds"] == [-1.0, 1.0]
    assert "placeholder" in notice and "set the bounds manually" in notice


def test_retype_round_trip_returns_original_and_restores_origin_exactly():
    pd = _REG["Tb_K"]
    u = {"prior_type": "uniform", "bounds": [259.5, 271.0]}
    lu, _ = gui_config.retype_prior_entry(pd, "log-uniform", u)
    back_no_origin, _ = gui_config.retype_prior_entry(pd, "uniform", lu)
    assert back_no_origin["bounds"] == pytest.approx(u["bounds"], rel=1e-12)
    back, notice = gui_config.retype_prior_entry(pd, "uniform", lu, origin=u)
    assert back == u and notice is None
    # a user edit in the log-uniform stage is NOT overridden by the origin
    lu_edit = {"prior_type": "log-uniform", "bounds": [2.42, 2.433]}
    back, _ = gui_config.retype_prior_entry(pd, "uniform", lu_edit, origin=u)
    assert back["bounds"] == pytest.approx([10 ** 2.42, 10 ** 2.433])


@pytest.mark.gui
def test_retype_round_trip_in_page_is_not_an_edit_and_notice_is_visible(fx):
    data = fx["fx_full"]
    at = _app()
    _load_via_selector(at, "fx_full", data=data)
    _w(at, "selectbox", "prior_type_alpha").select("log-uniform")
    at.run()
    _ok(at)
    lu = at.session_state["inference_param_space"]["alpha"]
    assert lu["prior_type"] == "log-uniform"
    assert lu["bounds"] == pytest.approx([math.log10(0.15), math.log10(0.45)])
    _w(at, "selectbox", "prior_type_alpha").select("uniform")
    at.run()
    at.run()
    _ok(at)
    assert at.session_state["inference_param_space"]["alpha"] == \
        data["param_space"]["alpha"]
    diff, _ = _built_diff(at, data)
    assert diff == []
    # impossible conversion: log10_zeta [-3, 2] has non-positive bounds
    _w(at, "selectbox", "prior_type_log10_zeta").select("log-uniform")
    at.run()
    _ok(at)
    assert any("non-positive" in str(w.value) for w in at.warning)


# -- finding 3: k2 channel leakage -------------------------------------------

@pytest.mark.gui
@pytest.mark.parametrize("stem,have", [("fx_re_only", "Re_k2"),
                                       ("fx_im_only", "Im_k2")])
def test_single_k2_channel_config_round_trips_without_gaining_the_other(
        fx, stem, have):
    data = fx[stem]
    at = _app()
    _load_via_selector(at, stem, data=data, reruns=3)
    diff, built = _built_diff(at, data)
    assert diff == []
    assert list(built.observables) == [have]
    # explicitly adding the missing channel is honoured
    other_box = "use_im_k2_obs" if have == "Re_k2" else "use_re_k2_obs"
    assert _w(at, "checkbox", other_box).value is False
    _w(at, "checkbox", other_box).check()
    at.run()
    _ok(at)
    _, built = _built_diff(at, data)
    assert set(built.observables) == {"Re_k2", "abs_Im_k2" if have == "Re_k2"
                                      else "Im_k2"} or \
        set(built.observables) == {"Re_k2", "abs_Im_k2"}


@pytest.mark.gui
def test_k2_master_toggle_off_on_restores_loaded_values(fx):
    data = fx["fx_full"]
    at = _app()
    _load_via_selector(at, "fx_full", data=data)
    _w(at, "checkbox", "use_k2_obs").uncheck()
    at.run()
    _, built = _built_diff(at, data)
    assert "Re_k2" not in built.observables and "Im_k2" not in built.observables
    _w(at, "checkbox", "use_k2_obs").check()
    at.run()
    at.run()
    diff, _ = _built_diff(at, data)
    assert diff == []                    # loaded values, not body defaults


# -- finding 4: non-pair observables -----------------------------------------

def test_merge_observables_passes_through_unrendered_and_drops_removed():
    loaded = {"Re_k2": [0.6, 0.05], "Re_h2": [1.2, 0.1, 0.2],
              "Im_k2": [0.13, 0.03]}
    form = {"Re_k2": (0.6, 0.05)}               # user removed Im_k2
    out = gui_config.merge_observables(loaded, form)
    assert out == {"Re_k2": [0.6, 0.05], "Re_h2": [1.2, 0.1, 0.2]}


@pytest.mark.gui
def test_three_element_observable_survives_load_and_run(fx):
    data = fx["fx_obs3"]
    at = _app()
    _load_via_selector(at, "fx_obs3", data=data)
    diff, built = _built_diff(at, data)
    assert diff == []
    assert built.observables["Re_h2"] == [1.2, 0.1, 0.2]
    _w(at, "checkbox", "use_re_k2_obs").uncheck()      # explicit removal
    at.run()
    _ok(at)
    _, built = _built_diff(at, data)
    assert "Re_k2" not in built.observables
    assert built.observables["Re_h2"] == [1.2, 0.1, 0.2]


# -- finding 5: carry-over between configs -----------------------------------

@pytest.mark.gui
def test_loading_empty_blocks_after_another_config_does_not_carry_over(fx):
    at = _app()
    _load_via_selector(at, "fx_full", data=fx["fx_full"])
    full_diff, _ = _built_diff(at, fx["fx_full"])
    assert full_diff == []
    data = fx["fx_bare"]
    _load_via_selector(at, "fx_bare", data=data)
    diff, built = _built_diff(at, data)
    assert diff == []
    assert built.observables == {"Re_k2": [0.5, 0.05]}
    assert not built.structure_cache_path
    assert built.sampler_settings == InferenceConfig.from_dict(
        data).sampler_settings
    assert "n_effective" not in built.sampler_settings or \
        built.sampler_settings["n_effective"] != 77
    # an EMPTY observables block also clears the form (no carry-over)
    _load_via_selector(at, "fx_full", data=fx["fx_full"])
    _load_via_selector(at, "fx_noobs", data=fx["fx_noobs"])
    assert at.session_state["inference_observables"] == {}
    assert _w(at, "checkbox", "use_k2_obs").value is False
    assert _w(at, "checkbox", "use_cmr2_obs").value is False


# -- finding 2/10: precision of loaded bounds --------------------------------

def test_prior_number_format_is_lossless_for_shown_values():
    assert gui_config.prior_number_format(-3.0, 1.845098) == "%.6f"
    assert gui_config.prior_number_format(0.15, 0.45) == "%.4f"
    assert gui_config.prior_number_format(-2e-5, 2e-5) == "%.3e"
    assert gui_config.prior_number_format(1.2345678e-5) == "%.7e"
    # the format reproduces the value exactly for display
    for v in (1.845098, 259.5, 0.123456789):
        fmt = gui_config.prior_number_format(v)
        assert float(fmt % v) == pytest.approx(v, rel=1e-9)


@pytest.mark.gui
def test_six_decimal_bound_survives_load_rerun_build_and_is_not_an_edit(fx):
    data = fx["fx_prec"]
    at = _app()
    _load_via_selector(at, "fx_prec", data=data, reruns=4)
    assert at.session_state["inference_param_space"]["log10_zeta"][
        "bounds"] == [-3.0, 1.845098]
    diff, built = _built_diff(at, data)
    assert diff == []
    assert built.param_space["log10_zeta"]["bounds"] == [-3.0, 1.845098]
    ni = _w(at, "number_input", "upper_log10_zeta_uniform")
    assert ni.value == 1.845098
    # AppTest applies no browser formatting, so inspect the widget protocol:
    # the browser shows 6 decimals ("%.6f") -> "1.845098"; the stored value
    # is submitted untouched unless the user edits the field.
    assert ni.proto.format == "%.6f"
    assert ni.proto.step == pytest.approx(0.01)
    assert _w(at, "number_input", "lower_log10_zeta_uniform").proto.format == "%.6f"


# -- findings 6, 7, 10: cache check -------------------------------------------

def _write_cache(tmp_path, **kw):
    p = tmp_path / "c.pkl"
    p.write_bytes(pickle.dumps(kw))
    return p


def test_cache_check_reports_frozen_zb_axis_separately(tmp_path):
    p = _write_cache(
        tmp_path, Tb_K_grid=[250.0, 260.0], zb_km_grid=[5.0, 25.0, 45.0],
        frozen_zb_km_grid=[46.5, 55.5, 65.5], structures=[{}] * 3)
    ps = {"zb_km": {"prior_type": "uniform", "bounds": [5.0, 60.0]}}
    rep = gui_config.check_structure_cache(ps, {}, str(p))
    assert rep["axes"]["zb_km"]["range"] == (5.0, 45.0)
    assert rep["axes"]["zb_km"]["n"] == 3
    assert rep["axes"]["frozen_zb_km"]["range"] == (46.5, 65.5)
    assert [i["param"] for i in rep["issues"]] == ["zb_km"]   # 60 > ocean max
    ok = gui_config.check_structure_cache(
        {"zb_km": {"prior_type": "uniform", "bounds": [5.0, 45.0]}}, {},
        str(p))
    assert ok["issues"] == []


def test_enceladus_cache_ocean_axis_is_not_merged_with_frozen_axis():
    cfg = _file("enceladus_cassini_isostasy_ocean6D_rhorock")
    cache = REPO / cfg["structure_cache_path"]
    if not cache.exists():
        pytest.skip("Enceladus cache not present")
    ps = dict(cfg["param_space"])
    ps["zb_km"] = {"prior_type": "uniform", "bounds": [5.0, 60.0]}
    rep = gui_config.check_structure_cache(ps, cfg.get("fixed_params"),
                                           cfg["structure_cache_path"])
    assert rep["axes"]["zb_km"]["range"] == (5.0, 45.0)
    assert rep["axes"]["frozen_zb_km"]["range"][0] > 45.0
    assert any(i["param"] == "zb_km" for i in rep["issues"])


def test_salinity_axis_with_zero_node_is_readable(tmp_path):
    p = _write_cache(tmp_path, Tb_K_grid=[250.0, 260.0],
                     wOcean_ppt_grid=[0.0, 10.0, 100.0], structures=[None] * 6)
    rep = gui_config.check_structure_cache(
        {"log10_wOcean_ppt": {"prior_type": "uniform",
                              "bounds": [1.2, 1.8]}}, {}, str(p))
    assert rep["status"] == "ok"
    w = rep["axes"]["log10_wOcean_ppt"]
    assert w["range"] == pytest.approx((1.0, 2.0)) and w["zero_node"] is True
    assert w["n"] == 3
    assert rep["issues"] == [] and rep["quality"]["n_none"] == 6


def test_cache_quality_counts_none_and_inconsistent_nodes(tmp_path):
    good = {"r_m": [1, 2, 3], "P_MPa": [1, 2, 3]}
    bad = {"r_m": [1, 2, 3, 4], "P_MPa": [1, 2, 3]}
    p = _write_cache(tmp_path, Tb_K_grid=[250.0], structures=[good, None, bad, bad])
    q = gui_config.check_structure_cache({}, {}, str(p))["quality"]
    assert (q["n_nodes"], q["n_none"], q["n_inconsistent"]) == (4, 1, 2)
    assert q["inconsistent_examples"] == [2, 3]


def test_callisto_mgso4_cache_reports_inconsistent_nodes_read_only():
    rel = "PlanetProfile/Test/mcmc_results/Callisto/C2_andrade/callisto_mgso4_structure_grid.pkl"
    if not (REPO / rel).exists():
        pytest.skip("callisto_mgso4 cache not present")
    q = gui_config.check_structure_cache({}, {}, rel)["quality"]
    assert q["n_inconsistent"] == 5 and q["n_none"] == 0


def test_cache_read_is_memoized_on_path_mtime_size(tmp_path, monkeypatch):
    p = _write_cache(tmp_path, Tb_K_grid=[250.0, 260.0], structures=[])
    calls = []
    real = gui_config._load_axes
    monkeypatch.setattr(gui_config, "_load_axes",
                        lambda path: (calls.append(path), real(path))[1])
    assert not gui_config.axes_cached(p)
    for _ in range(3):
        gui_config.check_structure_cache({}, {}, str(p))
    assert len(calls) == 1 and gui_config.axes_cached(p)
    p.write_bytes(pickle.dumps({"Tb_K_grid": [1.0, 2.0, 3.0], "structures": []}))
    gui_config.check_structure_cache({}, {}, str(p))     # size changed
    assert len(calls) == 2


@pytest.mark.gui
def test_page_does_not_reread_cache_on_rerun(fx, monkeypatch):
    data = fx["fx_full"]
    calls = []
    real = gui_config._load_axes
    monkeypatch.setattr(gui_config, "_load_axes",
                        lambda path: (calls.append(path), real(path))[1])
    at = _app()
    _load_via_selector(at, "fx_full", data=data, reruns=1)
    n = len(calls)
    for _ in range(3):
        at.run()
        _ok(at)
    assert len(calls) == n


# -- finding 11: widgets must not keep the previous config's frontend value --

def _assert_loaded_params_hold(at, stem):
    data = _file(stem)
    ms = _w(at, "multiselect", "param_multiselect")
    assert list(ms.value) == list(data["param_space"])
    assert at.session_state["inference_selected_params"] == list(
        data["param_space"])
    built = gui_config.run_config_from_state(_state(at))
    assert _norm(built.param_space) == _norm(
        InferenceConfig.from_dict(data).param_space)


@pytest.mark.gui
def test_widgets_get_new_identity_on_every_config_load():
    """Streamlit 1.59 keys these widgets' identity on the key alone, so in a
    real browser a popped key keeps the OLD frontend value and re-sends it on
    the next interaction (AppTest has no frontend and cannot show that).
    The observable, harness-independent property: every Apply moves the
    form widgets to fresh keys (new identity), so no stale value can win."""
    at = _app()
    _load_via_selector(at, "test50_titan_noocean_andrade_5D")
    g1 = at.session_state["inference_load_gen"]
    keys1 = {w.key for w in at.multiselect} | {w.key for w in at.slider}
    _load_via_selector(at, "europa_clipper_v5_geodesy_11D")
    g2 = at.session_state["inference_load_gen"]
    keys2 = {w.key for w in at.multiselect} | {w.key for w in at.slider}
    assert g2 > g1
    assert keys1.isdisjoint(keys2)
    assert not any(k in keys2 for k in ("param_multiselect", "n_eff_slider"))


@pytest.mark.gui
@pytest.mark.parametrize("stems", [
    ["europa_clipper_v5_geodesy_11D"],
    ["test50_titan_noocean_andrade_5D", "europa_clipper_v5_geodesy_11D",
     "europa_clipper_v6_freegrav_11D"]])
def test_params_hold_after_widget_interaction_following_apply(stems):
    at = _app()
    for stem in stems:
        _load_via_selector(at, stem)
        _assert_loaded_params_hold(at, stem)
        # a real widget interaction after Apply (the browser repro: moving
        # the ESS slider made the stale multiselect value win)
        _w(at, "slider", "n_eff_slider").set_value(900)
        at.run()
        _ok(at)
        _assert_loaded_params_hold(at, stem)
        space = at.session_state["inference_param_space"]
        data = _file(stem)
        assert _norm(space) == _norm(data["param_space"])
        _w(at, "number_input", "random_seed").set_value(11)
        at.run()
        _ok(at)
        _assert_loaded_params_hold(at, stem)
