"""AppTest coverage for the scope-bound Cassini-Enceladus SBI slot.

Reviewer ruling 2026-10-27 (validation_reports/enceladus_isostasy/sbi/
RESULT.md): the slot shows the verbatim deploy-scope text, routes every run
through enceladus_sbi_guard.guarded_posterior, labels the Cassini datum
"validated", and REFUSES a thin-shell query with no posterior. Replaces the
awaiting-artifact placeholder test (tests/enceladus_placeholder_test.py).
"""

import ast
import sys
from pathlib import Path

import numpy as np
import pytest


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "PlanetProfileApp"))
PAGE = REPO / "PlanetProfileApp/pages/Inference.py"
ART = "enceladus_ocean6D_rhorock_posterior_1m.pt"
SLOT_KEY = "enceladus_cassini"

pytest.importorskip("streamlit")


def _registry():
    tree = ast.parse(PAGE.read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == "_SBI_ARTIFACT_SLOTS"
                for t in node.targets):
            return ast.literal_eval(node.value)
    raise AssertionError("_SBI_ARTIFACT_SLOTS not found")


SLOT = _registry()[ART]
_READY = ((REPO / "PlanetProfile/Inference/sbi_artifacts" / ART).exists()
          and (REPO / SLOT["cache_path"]).exists())


def test_scope_note_is_verbatim_guard_text():
    from PlanetProfile.Inference.enceladus_sbi_guard import SCOPE_TEXT
    assert SLOT["scope_note"] == SCOPE_TEXT
    assert SLOT["scope_guard"] == "enceladus_ocean_branch"
    assert len(SLOT["validated_x"]) == 7          # x_obs + six G4 sets
    assert SLOT["validated_x"][0] == SLOT["default_obs"]


def _open_enceladus(app):
    app.session_state["inference_exec_mode"] = "amortized"
    app.run()
    assert not app.exception, [e.value for e in app.exception]
    analysis = next(s for s in app.selectbox
                    if str(getattr(s, "key", "")) == "amort_analysis")
    opt = next(o for o in analysis.options if "Enceladus" in str(o))
    analysis.select(opt)
    app.run()
    assert not app.exception, [e.value for e in app.exception]


def _generate(app):
    btn = next(b for b in app.button
               if "Generate Posterior" in str(getattr(b, "label", "")))
    assert not btn.disabled
    btn.click()
    app.run()
    assert not app.exception, [e.value for e in app.exception]


def _text(app):
    out = []
    for kind in ("error", "warning", "success", "info", "caption"):
        out += [str(e.value) for e in getattr(app, kind)]
    return " ".join(out)


@pytest.mark.skipif(not _READY, reason="Enceladus artifact / cache absent")
def test_cassini_datum_accepted_validated_with_scope_text():
    from streamlit.testing.v1 import AppTest
    from PlanetProfile.Inference.enceladus_sbi_guard import SCOPE_TEXT
    app = AppTest.from_file(str(PAGE), default_timeout=900)
    _open_enceladus(app)
    assert SCOPE_TEXT in _text(app)          # pre-run, on the config panel
    # Persistent superseded-simulator warning (frozen libration convention).
    assert SLOT["artifact_warning"] in _text(app)
    assert "SUPERSEDED SIMULATOR" in _text(app)
    assert "disabled for this model" in _text(app)
    app.select_slider(key="amort_n_post").set_value(1000)
    app.select_slider(key="amort_n_derived").set_value(500)
    _generate(app)
    res = app.session_state["inference_results"]
    sg = res.metadata["scope_guard"]
    assert sg["validated"] and sg["dropped_fraction"] == 0.0
    assert sg["scope_text"] == SCOPE_TEXT
    assert len(res.samples) == 1000
    assert "validated against pocoMC at this datum" in _text(app)
    assert "frozen-libration-convention simulator only" in _text(app)
    assert SLOT["artifact_warning"] in _text(app)   # also on the result
    zb = float(np.median(res.samples[:, list(res.param_names).index("zb_km")]))
    print("Cassini-datum zb median (km):", zb)
    # Frozen-convention flow + frozen-convention recompute: ~26.75-26.8 km
    # (NOT the repaired headline 26.5).
    assert 26.5 < zb < 27.1
    # No tidal channel: the k2 and heating panels explain the absence
    # instead of comparing NaN placeholders (user report 2026-10-30).
    assert res.metadata["tidal_solve_skipped"] is True
    assert _text(app).count("This model has no tidal channel") >= 2
    assert "k₂ scatter unavailable" not in _text(app)
    # Induction display reference (display only): both core curves present
    # when the Ae sidecar ships next to the cache.
    from PlanetProfile.Inference.enceladus_induction_display import (
        CAPTION, sidecar_path_for)
    ind = res.metadata["enceladus_induction"]
    if sidecar_path_for(REPO / SLOT["cache_path"]).exists():
        assert set(ind["curves"]) == {"pp_core", "core_25Sm"}
        assert len(ind["curves"]["pp_core"]["orbital"]["Ae_re"]) == 1000
        assert CAPTION in _text(app)
    else:
        assert "missing_sidecar" in ind


@pytest.mark.skipif(not _READY, reason="Enceladus artifact / cache absent")
def test_thin_shell_query_refused_without_posterior():
    from streamlit.testing.v1 import AppTest
    app = AppTest.from_file(str(PAGE), default_timeout=900)
    _open_enceladus(app)
    # The thinnest truth of the held-out thin shard (zb 5.1 km): its large
    # libration (1.27 deg) puts the flow in the thin-shell regime.
    app.number_input(key=f"amort_obs_{SLOT_KEY}_libration_deg").set_value(
        1.2703469626030577)
    app.select_slider(key="amort_n_post").set_value(1000)
    _generate(app)
    assert app.session_state["inference_results"] is None
    txt = _text(app)
    assert "Refused by the scope guard" in txt and "use pocoMC" in txt


# ---------------------------------------------------------------------------
# Pairing-sensitive check (reviewer binding condition): zb comes from the flow
# alone, so the zb window above cannot detect a flow/simulator mismatch.  The
# pushforward libration (flow draws -> forward model's libration) can: the
# frozen-convention flow pushed through the frozen snapshot config lands about
# -0.05 sigma_infl from the observed 0.092 deg (gate record <= 0.056); through
# the live repaired config it sits ~ -0.37 sigma.  The median's Monte-Carlo
# error at 2000 draws is ~0.03 sigma and the draws are not bit-reproducible
# across torch builds, so the two bounds sit between the cases with margin on
# both sides instead of at the gate's 0.06.
# ---------------------------------------------------------------------------
LIB_TOL_SIGMA = 0.15
LIB_MISPAIRED_MIN_SIGMA = 0.25
_LIVE_CFG = ("PlanetProfile/Inference/configs/"
             "enceladus_cassini_isostasy_ocean6D_rhorock.json")


def _libration_pushforward_offset(config_rel, n_draws=2000, seed=75):
    """(median pushforward libration - observed) / sigma_infl at the Cassini
    datum.  Same flow draws (guard seed/sampler), same forward-model path the
    likelihood uses (MCMCRunner._derive_libration_deg), hole draws dropped as
    the guard drops them; sigma_infl = hypot(sigma_obs, sigma_model_add) read
    from the config as the runner's _inflate_sigma does."""
    import json
    from PlanetProfile.Inference.inference_core import InferenceConfig
    from PlanetProfile.Inference.sbi_runner import SBIRunner
    from PlanetProfile.Inference.enceladus_sbi_guard import is_hole

    d = json.loads((REPO / config_rel).read_text(encoding="utf-8"))
    d["mode"] = "sbi"
    # Observed values = the slot's Cassini datum; sigmas stay the config's.
    d["observables"] = {n: [SLOT["default_obs"][n], d["observables"][n][1]]
                        for n in SLOT["default_obs"]}
    d["structure_cache_path"] = SLOT["cache_path"]
    cfg = InferenceConfig.from_dict(d)
    mcmc = SBIRunner(cfg)._get_mcmc_runner()
    flow = SBIRunner.load_artifact(REPO / "PlanetProfile/Inference/sbi_artifacts" / ART)
    x_obs = dict(SLOT["default_obs"])
    th = np.asarray(flow.sample_posterior(
        x_obs, n_samples=n_draws, seed=seed, reject_outside_prior=True), float)
    lib = []
    for t in th:
        if is_hole(mcmc, t):
            continue
        v = mcmc._derive_libration_deg(mcmc._expand_theta(t))
        if v is not None and np.isfinite(v):
            lib.append(v)
    assert len(lib) > 0.9 * n_draws
    sig_obs = float(cfg.observables["libration_deg"][1])
    sig_add = float(((cfg.metadata or {}).get("sigma_model_add") or {})
                    .get("libration_deg", 0.0))
    sigma_infl = float(np.hypot(sig_obs, sig_add))
    return (float(np.median(lib)) - x_obs["libration_deg"]) / sigma_infl


@pytest.mark.skipif(not _READY, reason="Enceladus artifact / cache absent")
def test_pushforward_libration_matches_snapshot_pairing():
    off = _libration_pushforward_offset(SLOT["config_path"])
    print("pushforward libration offset, snapshot config (sigma_infl):", off)
    assert abs(off) <= LIB_TOL_SIGMA


@pytest.mark.skipif(not _READY, reason="Enceladus artifact / cache absent")
def test_pushforward_libration_detects_wrong_pairing():
    """Negative control: the same flow through the live REPAIRED config must
    fall outside the tolerance, so the check above can see a mispairing."""
    assert _LIVE_CFG != SLOT["config_path"]
    off = _libration_pushforward_offset(_LIVE_CFG)
    print("pushforward libration offset, live repaired config (sigma_infl):", off)
    assert abs(off) > LIB_MISPAIRED_MIN_SIGMA
