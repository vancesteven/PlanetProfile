"""End-to-end provenance and sampling regression for deployed SBI slots.

The GUI registry is intentionally read from its source assignment rather than
imported: importing the Streamlit page executes its UI body.  The registry is
a literal dictionary by design, so ``ast.literal_eval`` gives this test the
same slot records without starting the app.
"""

import ast
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
import torch

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.sbi_runner import SBIRunner


REPO_ROOT = Path(__file__).resolve().parents[1]
REGISTRY_PATH = REPO_ROOT / "PlanetProfileApp/pages/Inference.py"
ARTIFACT_DIR = REPO_ROOT / "PlanetProfile/Inference/sbi_artifacts"
REFERENCE_PATH = REPO_ROOT / "tests/data/slot_reproduction_test50.json"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_registry() -> dict:
    tree = ast.parse(REGISTRY_PATH.read_text(encoding="utf-8"))
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if any(
            isinstance(target, ast.Name) and target.id == "_SBI_ARTIFACT_SLOTS"
            for target in node.targets
        ):
            return ast.literal_eval(node.value)
    raise AssertionError(f"_SBI_ARTIFACT_SLOTS not found in {REGISTRY_PATH}")


def _slot_paths(slot_name: str, slot: dict) -> tuple[Path, Path, Path]:
    artifact_name = slot.get("artifact_filename", slot_name)
    return (
        ARTIFACT_DIR / artifact_name,
        REPO_ROOT / slot["config_path"],
        REPO_ROOT / slot["cache_path"],
    )


SLOTS = _load_registry()
READY_SLOTS = {
    name: slot
    for name, slot in SLOTS.items()
    if slot.get("config_path")
    and slot.get("cache_path")
    and all(path.exists() for path in _slot_paths(name, slot))
}


@pytest.mark.parametrize("slot_name", READY_SLOTS, ids=READY_SLOTS)
def test_ready_slot_config_hash_matches_artifact(slot_name):
    """Every locally complete slot resolves to its artifact's SBI config."""
    slot = READY_SLOTS[slot_name]
    artifact_path, config_path, cache_path = _slot_paths(slot_name, slot)

    config = InferenceConfig.from_json(config_path)
    assert (REPO_ROOT / config.structure_cache_path).resolve() == cache_path.resolve()

    # Committed source configs default to the MCMC execution path.  Training
    # changes only this execution-mode tag before hashing/saving the SBI flow.
    config.mode = "sbi"
    artifact = torch.load(artifact_path, map_location="cpu", weights_only=False)
    documented_hash = slot.get("artifact_config_hash")
    accepted_hashes = {config.generate_hash()}
    if documented_hash:
        accepted_hashes.add(documented_hash)

    assert artifact["config_hash"] in accepted_hashes


def _committed_cache_manifests() -> list[tuple[Path, dict]]:
    manifests = []
    roots = (
        REPO_ROOT / "validation_reports",
        ARTIFACT_DIR / "validation_reports",
    )
    for root in roots:
        for path in sorted(root.rglob("*manifest*.json")):
            payload = json.loads(path.read_text(encoding="utf-8"))
            if payload.get("structure_cache_path") and payload.get(
                "structure_cache_sha256"
            ):
                manifests.append((path, payload))
    return manifests


def test_ready_slot_cache_hashes_match_committed_manifests():
    """Validate each applicable manifest against the cache used by a slot."""
    ready_cache_paths = {
        slot["cache_path"]
        for slot in READY_SLOTS.values()
    }
    checked = []
    computed = {}
    for manifest_path, manifest in _committed_cache_manifests():
        cache_rel = Path(manifest["structure_cache_path"]).as_posix()
        if cache_rel not in ready_cache_paths:
            continue
        cache_path = REPO_ROOT / cache_rel
        actual = computed.setdefault(cache_rel, _sha256(cache_path))
        assert actual == manifest["structure_cache_sha256"], manifest_path
        checked.append(manifest_path)

    # Guards against a broken search silently turning this into a no-op while
    # allowing future manifests to join the same automatic scan.
    assert checked


@pytest.mark.slow
def test_test50_fixed_seed_conditioning_reproduces_reference():
    """The designated fast slot reproduces its committed 1k-draw summary."""
    reference = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))
    slot_name = reference["slot"]
    slot = READY_SLOTS[slot_name]
    artifact_path, _, cache_path = _slot_paths(slot_name, slot)

    assert artifact_path == REPO_ROOT / reference["artifact_path"]
    assert cache_path == REPO_ROOT / reference["cache_path"]
    assert _sha256(artifact_path) == reference["artifact_sha256"]
    assert _sha256(cache_path) == reference["cache_sha256"]
    assert slot["default_obs"] == reference["conditioning"]

    runner = SBIRunner.load_artifact(
        artifact_path,
        validated_version_pairs=slot.get("validated_version_pairs"),
    )
    assert runner.artifact_meta["config_hash"] == reference["artifact_config_hash"]
    assert runner.param_names == reference["param_names"]

    samples = runner.sample_posterior(
        reference["conditioning"],
        n_samples=reference["n_draws"],
        seed=reference["seed"],
    )
    assert samples.shape == (reference["n_draws"], len(runner.param_names))

    actual = {}
    for index, name in enumerate(runner.param_names):
        values = samples[:, index]
        actual[name] = {
            "mean": float(np.mean(values)),
            "std": float(np.std(values, ddof=1)),
            "median": float(np.median(values)),
            "q05": float(np.quantile(values, 0.05)),
            "q95": float(np.quantile(values, 0.95)),
        }

    tolerance = reference["tolerance"]
    for name, expected_stats in reference["statistics"].items():
        for statistic, expected in expected_stats.items():
            np.testing.assert_allclose(
                actual[name][statistic],
                expected,
                rtol=tolerance["rtol"],
                atol=tolerance["atol"],
                err_msg=f"{name}.{statistic}",
            )
