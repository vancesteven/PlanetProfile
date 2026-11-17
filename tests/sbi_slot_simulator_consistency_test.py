"""Slot config vs. training-time simulator consistency.

``InferenceConfig.generate_hash()`` ignores ``metadata``, but the runner reads
simulator-defining blocks from it (libration formulation/correction, isostasy
closure, gravity reference, noise model...).  A slot whose ``config_path`` was
edited after training therefore keeps a matching hash while silently pairing
the flow with a different forward model (Enceladus, 2a75a910).

What artifacts store: schema-1 ``.pt`` files embed ``config_hash`` and
``git_sha`` but NOT the training config or its metadata.  So the baseline here
is ``git show <artifact git_sha>:<slot config_path>`` (no baseline is invented
from the current files).  Frozen snapshots use their declared ``source_file`` at the same git_sha.
Slots whose config does not exist at the artifact's git_sha fall back to the
FIRST-ADD baseline: the commit that first added the slot's config file
(``git log --diff-filter=A``), accepted ONLY if that version's
``InferenceConfig`` (mode='sbi') hashes to the artifact's stored
``config_hash``; results are labelled "first-add baseline (hash-matched)".
Only a slot with no hash-matching commit is skipped (never passed).

Compared: the ``metadata`` keys the runner reads (RUNNER_METADATA_KEYS) PLUS
the top-level fields that ``generate_hash()`` does not cover but that reach
the forward model (UNHASHED_SIMULATOR_FIELDS).

LIMITATION: this test compares config CONTENT only.  It cannot see
forward-model CODE drift since training (a changed ``mcmc_runner`` /
``forward_models`` / ``Librations`` with an unchanged config passes).
"""
import ast
import json
import subprocess
from pathlib import Path

import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("sbi")

REPO = Path(__file__).resolve().parents[1]
PAGE = REPO / "PlanetProfileApp/pages/Inference.py"
ARTIFACT_DIR = REPO / "PlanetProfile/Inference/sbi_artifacts"

# Top-level ``config.metadata`` keys read by the runner / forward model /
# cache path (grep of ``metadata`` in mcmc_runner.py, forward_models.py,
# cache_builder.py, enceladus_induction_display.py):
RUNNER_METADATA_KEYS = (
    "libration_model_correction",   # formulation, factor, sidecar, figure ruling
    "isostasy",                     # stack_mass_closure, moi_target, H_obs...
    "bulk_overrides",               # Cmeasured
    "gravity_ref_radius_m",
    "gravity_j2_over_c22",
    "observable_correlations",      # likelihood / training noise
    "sigma_model_add",              # inflated sigmas
    "branch_model",
    "structure_cache_spec",
    "induction_provenance",         # presence drives the Ae display sidecar
)


# Top-level InferenceConfig fields NOT covered by generate_hash() (which hashes
# only mode, bodyname, param_space, observables, sampler_settings, random_state,
# arrhenius_params, induction_bounds, gravity_forward_model) that reach the
# forward model / structure lookup.
UNHASHED_SIMULATOR_FIELDS = (
    "fixed_params",
    "derived_params",
    "param_groups",
    "planet_template_module",
    "ocean_overrides",
    "structure_cache_path",
)


def _registry():
    tree = ast.parse(PAGE.read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == "_SBI_ARTIFACT_SLOTS"
                for t in node.targets):
            return ast.literal_eval(node.value)
    raise AssertionError("_SBI_ARTIFACT_SLOTS not found")


SLOTS = {n: s for n, s in _registry().items()
         if s.get("config_path") and s.get("artifact_filename", n)
         and (ARTIFACT_DIR / s.get("artifact_filename", n)).exists()}


def _baseline(config_path, sha):
    r = subprocess.run(["git", "show", f"{sha}:{config_path}"], cwd=REPO,
                       capture_output=True, text=True)
    return json.loads(r.stdout) if r.returncode == 0 else None


def _sbi_hash(cfg_dict):
    from PlanetProfile.Inference.inference_core import InferenceConfig
    cfg = InferenceConfig.from_dict(cfg_dict)
    cfg.mode = "sbi"
    return cfg.generate_hash()


def _first_add_baseline(config_path, config_hash):
    """(sha, config dict) of the commit that first added ``config_path``, if
    its content hashes (mode='sbi') to the artifact's ``config_hash``."""
    r = subprocess.run(
        ["git", "log", "--diff-filter=A", "--format=%H", "--", config_path],
        cwd=REPO, capture_output=True, text=True)
    shas = [x for x in r.stdout.split() if x]
    if r.returncode != 0 or not shas:
        return None, None
    sha = shas[-1]                      # earliest addition
    base = _baseline(config_path, sha)
    if base is None:
        return None, None
    try:
        ok = _sbi_hash(base) == config_hash
    except Exception:
        ok = False
    return (sha, base) if ok else (None, None)


def _differences(slot_cfg, base_cfg):
    cur, old = slot_cfg.get("metadata") or {}, base_cfg.get("metadata") or {}
    bad = [f"metadata.{k}" for k in RUNNER_METADATA_KEYS
           if cur.get(k) != old.get(k)]
    bad += [k for k in UNHASHED_SIMULATOR_FIELDS
            if slot_cfg.get(k) != base_cfg.get(k)]
    return bad


@pytest.mark.parametrize("slot_name", SLOTS, ids=SLOTS)
def test_slot_config_runner_metadata_matches_training(slot_name, record_property):
    slot = SLOTS[slot_name]
    art = torch.load(ARTIFACT_DIR / slot.get("artifact_filename", slot_name),
                     map_location="cpu", weights_only=False)
    sha = art.get("git_sha")
    cur = json.loads((REPO / slot["config_path"]).read_text(encoding="utf-8"))
    base = _baseline(slot["config_path"], sha) if sha else None
    kind = f"config at artifact git_sha {sha}"
    if base is None and sha:
        # A frozen snapshot is a new file; its baseline is the ORIGINAL
        # config it copies, read at the artifact's own training git_sha
        # (not at the snapshot's self-declared commit).
        src = (cur.get("metadata") or {}).get("frozen_snapshot", {}).get(
            "source_file")
        base = _baseline(src, sha) if src else None
        kind = f"snapshot source {src} at artifact git_sha {sha}"
    if base is None:
        fsha, base = _first_add_baseline(slot["config_path"],
                                         art.get("config_hash"))
        if base is not None:
            kind = f"first-add baseline (hash-matched) {fsha}"
    if base is None:
        pytest.skip(f"NO TRUSTWORTHY BASELINE: {slot['config_path']} does not "
                    f"exist at artifact git_sha {sha} and no first-add commit "
                    "hash-matches the artifact config_hash")
    record_property("baseline", kind)
    print(f"[baseline] {slot_name}: {kind}")
    bad = _differences(cur, base)
    assert not bad, (
        f"slot {slot_name} ({kind}): config {slot['config_path']} differs from "
        f"the training baseline in simulator-defining fields {bad}; flow and "
        "forward recomputation would use different simulators. Point the "
        "slot at a frozen training snapshot.")


def test_detector_flags_known_enceladus_mismatch():
    """Guard against a vacuous test: the repaired MCMC config must differ from
    the artifact's training commit in libration_model_correction."""
    cfg = "PlanetProfile/Inference/configs/enceladus_cassini_isostasy_ocean6D_rhorock.json"
    base = _baseline(cfg, "53f333fa")
    if base is None:
        pytest.skip("training commit not available in this checkout")
    cur = json.loads((REPO / cfg).read_text(encoding="utf-8"))
    assert _differences(cur, base) == ["metadata.libration_model_correction"]


FROZEN = ("PlanetProfile/Inference/configs/"
          "enceladus_cassini_isostasy_ocean6D_rhorock_sbi_frozenconv.json")
ENC_ART = "enceladus_ocean6D_rhorock_posterior_1m.pt"


def test_enceladus_slot_uses_frozen_snapshot_with_legacy_libration():
    slot = _registry()[ENC_ART]
    assert slot["config_path"] == FROZEN
    assert "SUPERSEDED SIMULATOR" in slot["artifact_warning"]
    cur = json.loads((REPO / FROZEN).read_text(encoding="utf-8"))
    meta = cur["metadata"]
    assert meta["gui"]["status"] == "retired" and meta["gui"]["runnable"] is True
    assert meta["frozen_snapshot"]["source_commit"] == "1af55a95"
    lib = meta["libration_model_correction"]
    assert lib["multiplicative_factor"] == 1.008
    assert "formulation" not in lib and "elastic_factor_sidecar" not in lib


def test_frozen_snapshot_is_verbatim_1af55a95_plus_labels():
    old = _baseline("PlanetProfile/Inference/configs/"
                    "enceladus_cassini_isostasy_ocean6D_rhorock.json",
                    "1af55a95")
    if old is None:
        pytest.skip("commit 1af55a95 not available in this checkout")
    cur = json.loads((REPO / FROZEN).read_text(encoding="utf-8"))
    for k in ("gui", "frozen_snapshot"):
        cur["metadata"].pop(k)
    assert cur == old


def test_frozen_snapshot_hash_equals_artifact_hash():
    from PlanetProfile.Inference.inference_core import InferenceConfig
    art = torch.load(ARTIFACT_DIR / ENC_ART, map_location="cpu",
                     weights_only=False)
    cfg = InferenceConfig.from_json(REPO / FROZEN)
    cfg.mode = "sbi"
    assert cfg.generate_hash() == art["config_hash"]
