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
    data = _file(stem)
    stripped = copy.deepcopy(data)
    stripped["metadata"].pop("gui")
    h_with = InferenceConfig.from_dict(data).generate_hash()
    assert h_with == InferenceConfig.from_dict(stripped).generate_hash()
    try:
        head = subprocess.run(
            ["git", "show", f"HEAD:PlanetProfile/Inference/configs/{stem}.json"],
            cwd=REPO, capture_output=True, text=True, check=True).stdout
    except Exception:
        pytest.skip("git HEAD copy unavailable")
    assert h_with == InferenceConfig.from_dict(json.loads(head)).generate_hash()


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


def _ok(at):
    assert not at.exception, [e.value for e in at.exception]


def _load_via_selector(at, stem, reruns=2):
    at.run()
    _ok(at)
    at.checkbox(key="cfg_selector_show_all").check()
    at.run()
    body = _file(stem)["bodyname"]
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
    assert at.selectbox(key="prior_type_D_iceIh_km").value == "uniform"


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
        assert at.selectbox(key=f"prior_type_{k}").value == "truncated_gaussian"
        assert at.session_state["inference_param_space"][k] == v


@pytest.mark.gui
def test_single_bound_edit_changes_only_that_field():
    stem = "europa_clipper_v5_geodesy_11D"
    data = _file(stem)
    at = _app()
    _load_via_selector(at, stem)
    at.number_input(key="lower_D_iceIh_km_uniform").set_value(10.0)
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
    at.number_input(key="Re_k2_value").set_value(0.31)
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
    at.slider(key="n_eff_slider").set_value(900)
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
    at.number_input(key="upper_Tb_K_uniform").set_value(275.0)
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
    at.text_input(key="cache_path_input").set_value(bogus)
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
