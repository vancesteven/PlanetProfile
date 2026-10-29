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
