"""AppTest coverage for the Cassini-Enceladus awaiting-artifact slot."""

import importlib
import sys
from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "PlanetProfileApp"))

pytest.importorskip("streamlit")


def test_apptest_enceladus_placeholder_is_listed_nondefault_and_rendered():
    """The placeholder is visible, nondefault, and stops before loading."""
    from streamlit.testing.v1 import AppTest

    app = AppTest.from_file(
        str(REPO / "PlanetProfileApp/pages/Inference.py"),
        default_timeout=900,
    )
    app.session_state["inference_exec_mode"] = "amortized"
    app.run()
    assert not app.exception, str([error.value for error in app.exception])

    analysis = next(
        selector
        for selector in app.selectbox
        if str(getattr(selector, "key", "")) == "amort_analysis"
    )
    enceladus_option = next(
        option for option in analysis.options if "Enceladus" in str(option)
    )
    assert "Enceladus" not in str(analysis.value)

    analysis.select(enceladus_option)
    app.run()
    assert not app.exception, str([error.value for error in app.exception])

    version = next(
        selector
        for selector in app.selectbox
        if str(getattr(selector, "key", "")).startswith("amort_version_")
    )
    assert len(version.options) == 1
    assert "awaiting artifact" in str(version.value).lower()

    info_text = " ".join(str(info.value) for info in app.info)
    assert "Awaiting artifact" in info_text
    assert "this model slot is scaffolded" in info_text

    captions = " ".join(
        str(element.value)
        for element in app.main
        if type(element).__name__ == "Caption"
    )
    for expected in (
        "C₂₀/C₂₂/C₃₀",
        "H&M shape-input isostasy",
        "zb × Seawater salinity",
        "induction-response bands",
        "plans/active/enceladus-config-freeze.md",
    ):
        assert expected in captions

    generate = next(
        button
        for button in app.button
        if "Generate Posterior" in str(getattr(button, "label", ""))
    )
    assert generate.disabled

    inference_page = importlib.import_module("pages.Inference")
    slot = inference_page._SBI_ARTIFACT_SLOTS[
        "enceladus_cassini_awaiting_artifact"
    ]
    assert slot["slot_id"] == "enceladus_cassini"
    assert slot["bodyname"] == "Enceladus"
    assert slot["artifact_status"] == "awaiting_artifact"
    assert slot["artifact_filename"] is None
    assert slot["config_path"] is None
    assert slot["cache_path"] is None
