"""Focused provenance-schema regression for the Titan salt gate runner."""

import hashlib

from plans.scripts import titanG_ocean_run_gates as gates


def test_gate_manifest_records_complete_provenance_and_merges_partial_runs():
    spec = gates.COMPS["MgSO4"]
    artifact = gates.ROOT / (
        "PlanetProfile/Inference/sbi_artifacts/"
        "titan_freegrav_mgso4_posterior_1m.pt"
    )
    config = gates.ROOT / spec["cfg"]
    reports = gates.ROOT / "validation_reports/titan_freegrav_mgso4_1m"
    reference = reports / "reference/titan_freegrav_mgso4_reference_pooled.pkl"
    previous = {"gate_exit_codes": {"sbc": 0}}

    manifest = gates._build_manifest(
        comp="MgSO4", spec=spec, art=artifact, cfg=config,
        reports=reports, ref_mcmc=reference, rc={"limits": 2},
        previous=previous,
    )

    assert manifest["schema_version"] == 2
    assert manifest["gate_exit_codes"] == {
        "sbc": 0, "limits": 2, "crosscheck": "NOT_RUN",
    }
    assert manifest["artifact"] == artifact.relative_to(gates.ROOT).as_posix()
    assert manifest["config"] == config.relative_to(gates.ROOT).as_posix()
    assert manifest["reference_mcmc"] == reference.relative_to(gates.ROOT).as_posix()
    assert set(manifest["gate_report_paths"]) == set(gates._GATE_NAMES)
    assert set(manifest["gate_log_paths"]) == set(gates._GATE_NAMES)
    assert all(not value.startswith("/") for key, value in manifest.items()
               if key in {"artifact", "config", "reference_mcmc", "reports"})
    assert len(manifest["git_sha"]) == 40
    assert manifest["versions"]["torch"]
    assert manifest["versions"]["sbi"]
    assert manifest["artifact_sha256"] == hashlib.sha256(
        artifact.read_bytes()
    ).hexdigest()
