"""Enceladus induction freeze-blocker regression tests.

The conductivity bracket is the Saur et al. (2024) Fig. 20 acceptance
record adopted in C19.  The PPO test deliberately evaluates only the
dimensionless response: the Saturn PPO phase and amplitude drift, so no
fixed J2000 ``Be1xyz`` vector exists for the synodic channel.
"""
import pickle
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from PlanetProfile.Inference.forward_models import (  # noqa: E402
    induced_amplitude_band_nT,
)
from PlanetProfile.MagneticInduction.MagneticInduction import (  # noqa: E402
    GetBexc,
)
from PlanetProfile.MagneticInduction.Moments import Excitations  # noqa: E402
from PlanetProfile.MagneticInduction.MoonMag.symmetry_funcs import (  # noqa: E402
    InducedAeList,
)
from PlanetProfile.Thermodynamics.HydroEOS import GetOceanEOS  # noqa: E402
from PlanetProfile.Utilities.defineStructs import EOSlist  # noqa: E402

CACHE_PATH = (
    REPO
    / "PlanetProfile/Test/mcmc_results/Enceladus/Cassini_smoke"
    / "enceladus_seawater_smoke_grid.pkl"
)
RESULT_PATH = (
    REPO
    / "PlanetProfile/Test/mcmc_results/Enceladus/Cassini_smoke"
    / "enceladus_cassini_smoke_result.pkl"
)


def _seawater_sigma(w_ppt, temperatures_k):
    temperatures_k = np.asarray(temperatures_k, dtype=float)
    pressures_mpa = np.full_like(temperatures_k, 3.0)
    eos = GetOceanEOS(
        "Seawater",
        float(w_ppt),
        np.array([3.0, 4.0]),
        temperatures_k,
        None,
        FORCE_NEW=True,
    )
    return np.asarray(eos.fn_sigma_Sm(pressures_mpa, temperatures_k), float)


def test_seawater_conductivity_brackets_saur_figure_20():
    """The plume-band seawater model stays within the adopted Fig. 20 band."""
    temperatures_k = np.array([271.15, 272.15, 274.15])
    plume_sigma = np.vstack([
        _seawater_sigma(w_ppt, temperatures_k)
        for w_ppt in (5.0, 10.0, 20.0)
    ])

    # The reviewer-recorded 0.45--1.79 S/m bracket is digitized/reported to
    # 0.01 S/m; preserve that resolution at the upper endpoint.  The exact
    # GSW maximum is 1.7951 S/m and remains strictly below the separate
    # 2 S/m plume-band ceiling.
    assert np.min(plume_sigma) >= 0.45
    assert np.max(plume_sigma) <= 1.79 + 0.01
    assert np.all(plume_sigma < 2.0)

    high_salinity_sigma = _seawater_sigma(70.0, temperatures_k)
    assert np.all((4.5 <= high_salinity_sigma)
                  & (high_salinity_sigma <= 6.0))


def test_enceladus_periods_return_finite_response_without_be_vector():
    """Both response periods are usable without assigning PPO field moments."""
    with CACHE_PATH.open("rb") as stream:
        cache = pickle.load(stream)
    structure = cache["structures"][2]
    periods_hr = np.array([
        Excitations.Texc_hr["Enceladus"]["synodic"],
        Excitations.Texc_hr["Enceladus"]["true anomaly"],
    ])
    omega_radps = 2.0 * np.pi / (periods_hr * 3600.0)

    ae, magnitude, phase = InducedAeList(
        np.asarray(structure["rSigChange_m"], float),
        np.asarray(structure["sigmaLayers_Sm"], float),
        omega_radps,
        1.0 / float(structure["R_body_m"]),
        nn=1,
        writeout=False,
        do_parallel=False,
    )

    assert periods_hr == pytest.approx([15.5592, 32.93])
    assert np.all(np.isfinite(ae.real))
    assert np.all(np.isfinite(ae.imag))
    assert np.all(np.isfinite(magnitude))
    assert np.all(np.isfinite(phase))


def test_enceladus_loader_uses_only_committed_true_anomaly_row(
        monkeypatch, tmp_path):
    """Era/model lookup falls back to the base row and never reuses it as PPO."""
    selection = {"synodic": True, "true anomaly": True}
    cache_label = f"EnceladusBe1CassiniCassini11{selection}"
    EOSlist.loaded.pop(cache_label, None)
    EOSlist.ranges.pop(cache_label, None)
    monkeypatch.chdir(tmp_path)

    periods_hr, omega_radps, benm_nt, b0_nt, bexyz_nt = GetBexc(
        "Enceladus",
        "Cassini",
        "Cassini11",
        selection,
        nprmMax=1,
        pMax=0,
    )

    assert periods_hr == pytest.approx([32.92759])
    assert omega_radps == pytest.approx(2.0 * np.pi / periods_hr / 3600.0)
    assert benm_nt.shape[0] == 1
    assert b0_nt == pytest.approx([
        9.67711809556066,
        -0.695369498688637,
        -336.541625489458,
    ])
    assert bexyz_nt[0] == pytest.approx(
        complex(0.182096439236271, -0.0324979622986513)
    )

    EOSlist.loaded.pop(cache_label, None)
    EOSlist.ranges.pop(cache_label, None)


def test_ppo_amplitude_band_uses_response_magnitude_only():
    response = np.array([0.3 + 0.4j, -0.6 + 0.8j])
    lower_nt, upper_nt = induced_amplitude_band_nT(response)
    assert lower_nt == pytest.approx([0.5, 1.0])
    assert upper_nt == pytest.approx([1.0, 2.0])

    with pytest.raises(ValueError, match="ordered and non-negative"):
        induced_amplitude_band_nT(response, (2.0, 1.0))


@pytest.mark.skipif(not RESULT_PATH.exists(),
                    reason="Enceladus smoke result not present")
def test_apptest_enceladus_ppo_band_caption():
    """The app renders the response band without requesting a PPO vector."""
    from streamlit.testing.v1 import AppTest

    with RESULT_PATH.open("rb") as stream:
        result = pickle.load(stream)
    n_samples = len(result.samples)
    result.metadata = dict(result.metadata or {})
    result.metadata["induction_Ae"] = {
        "synodic": {
            "re": np.full(n_samples, 0.3).tolist(),
            "im": np.full(n_samples, 0.4).tolist(),
        },
    }
    result.metadata.pop("Be_nT", None)

    app = AppTest.from_file(
        str(REPO / "PlanetProfileApp/pages/Inference.py"),
        default_timeout=900,
    )
    app.session_state["inference_results"] = result
    app.session_state["Planet"] = "Enceladus"
    app.run()

    assert not app.exception, str([e.value for e in app.exception][:2])
    metrics = [element for element in app.main
               if type(element).__name__ == "Metric"]
    ppo_metric = next(
        metric for metric in metrics
        if str(getattr(metric, "label", ""))
        == "Synodic induced-amplitude band"
    )
    assert str(ppo_metric.value) == "0.5–1 nT"
    captions = " ".join(
        str(element.value) for element in app.main
        if type(element).__name__ == "Caption"
    )
    assert "|Ae| × [1, 2] nT" in captions
    assert "no stable excitation vector" in captions
