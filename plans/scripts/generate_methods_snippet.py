#!/usr/bin/env python
"""Generate source-backed methods snippets for deployed SBI campaigns.

Each output sentence is assembled from committed repository sources and is
immediately followed by an HTML source-pointer comment.  ``--verify`` reloads
those sources, regenerates the expected Markdown, checks cross-source hashes
and names, and renders the Markdown with Python-Markdown.

Usage:
  python plans/scripts/generate_methods_snippet.py nh3
  python plans/scripts/generate_methods_snippet.py all
  python plans/scripts/generate_methods_snippet.py all --verify
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import re
import sys
import warnings
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[2]
CONFIG_DIR = REPO / "PlanetProfile/Inference/configs"
ARTIFACT_DIR = REPO / "PlanetProfile/Inference/sbi_artifacts"
INDEX_PATH = ARTIFACT_DIR / "INDEX.md"


CAMPAIGNS = {
    "nh3": {
        "title": "Titan free-gravity NH3",
        "config": "test54_titan_nh3_freegrav.json",
        "artifact": "titan_freegrav_nh3_posterior_1m.pt",
        "report_dir": "validation_reports/titan_freegrav_nh3_1m",
        "cache_manifest": (
            "validation_reports/titan_freegrav_nh3_1m/reference/"
            "titanG_nh3_reference_manifest.json"
        ),
        "cache_sha_key": "structure_cache_sha256",
        "gate_kind": "ratification",
        "gate_source": "validation_reports/titan_freegrav_nh3_1m/RATIFICATION.md",
    },
    "mgso4": {
        "title": "Titan free-gravity MgSO4",
        "config": "test54_titan_mgso4_freegrav.json",
        "artifact": "titan_freegrav_mgso4_posterior_1m.pt",
        "report_dir": "validation_reports/titan_freegrav_mgso4_1m",
        "cache_manifest": "validation_reports/titan_freegrav_mgso4_1m/gen_manifest.json",
        "cache_sha_key": "structure_cache_sha256",
        "gate_kind": "ratification",
        "gate_source": "validation_reports/titan_freegrav_mgso4_1m/RATIFICATION.md",
    },
    "nacl": {
        "title": "Titan free-gravity NaCl",
        "config": "test54_titan_nacl_freegrav.json",
        "artifact": "titan_freegrav_nacl_posterior_1m.pt",
        "report_dir": "validation_reports/titan_freegrav_nacl_1m",
        "cache_manifest": "validation_reports/titan_freegrav_nacl_1m/gen_manifest.json",
        "cache_sha_key": "structure_cache_sha256",
        "gate_kind": "ratification",
        "gate_source": "validation_reports/titan_freegrav_nacl_1m/RATIFICATION.md",
    },
    "v5": {
        "title": "Europa Clipper v5 geodesy",
        "config": "europa_clipper_v5_geodesy_11D.json",
        "artifact": "europa_clipper_v5_geodesy_11D_posterior_1m.pt",
        "report_dir": "validation_reports/europa_clipper_v5_baseline_1m",
        "cache_manifest": (
            "PlanetProfile/Test/mcmc_results/Europa/Test52_seawater_v5/"
            "v5_reference_manifest.json"
        ),
        "cache_sha_key": "structure_cache_sha256",
        "gate_kind": "reviewer_json",
        "gate_source": "validation_reports/v5_gate_summary.json",
    },
    "v6": {
        "title": "Europa Clipper v6 free-gravity",
        "config": "europa_clipper_v6_freegrav_11D.json",
        "artifact": "europa_clipper_v6_freegrav_11D_posterior_1m.pt",
        "report_dir": "validation_reports/europa_clipper_v6_baseline_1m",
        "cache_manifest": (
            "PlanetProfile/Test/mcmc_results/Europa/Test53_seawater_v6/"
            "v6_reference_manifest.json"
        ),
        "cache_sha_key": "structure_cache_sha256",
        "gate_kind": "reviewer_json",
        "gate_source": "validation_reports/v6_gate_summary.json",
    },
    "v1p1": {
        "title": "Europa Galileo v1.1",
        "config": "europa_galileo_v1p1_8D.json",
        "artifact": "europa_galileo_v1p1_8D_posterior_1m.pt",
        "report_dir": "validation_reports/europa_galileo_v1p1_1m",
        "cache_manifest": (
            "validation_reports/europa_galileo_v1p1_1m/reference/"
            "titanG_reference_manifest.json"
        ),
        "cache_sha_key": "structure_cache_sha256",
        "gate_kind": "overall_heading",
        "gate_source": "validation_reports/europa_galileo_v1p1_1m/PHASE_GATES.md",
    },
    "v4": {
        "title": "Europa Clipper v4 geodesy",
        "config": "europa_clipper_v4_geodesy_11D.json",
        "artifact": "europa_clipper_v4_geodesy_11D_posterior_1m.pt",
        "report_dir": "validation_reports/europa_clipper_v4_1m",
        "cache_manifest": None,
        "cache_sha_key": None,
        "gate_kind": "phase_table",
        "gate_source": "validation_reports/europa_clipper_v4_1m/PHASE_GATES.md",
        "deployment_extra": (
            "Deployment stands (user-ratified) but the Tb/log10_w SBC defect "
            "is part of the persistent D–w degeneracy-direction signal under "
            "investigation (Machine B follow-ups B2/B5); whether v4 needs "
            "re-adjudication is a USER decision."
        ),
    },
    "test50": {
        "title": "Titan Test50 no-ocean Andrade",
        "config": "test50_titan_noocean_andrade_8D.json",
        "artifact": "titan_andrade_noocean_posterior.pt",
        "report_dir": "validation_reports/test50",
        "cache_manifest": "tests/data/slot_reproduction_test50.json",
        "cache_sha_key": "cache_sha256",
        "gate_kind": "index_excerpt",
        "gate_source": "PlanetProfile/Inference/sbi_artifacts/INDEX.md",
        "gate_excerpt": "ALL GREEN within domain (amended rules, 2026-07-11)",
    },
}


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPO).as_posix()


def _source_comment(*pointers: str) -> str:
    return f"<!-- sources: {'; '.join(pointers)} -->"


def _fmt(value: Any) -> str:
    if isinstance(value, float):
        return format(value, ".10g")
    return str(value)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_artifact(path: Path) -> dict[str, Any]:
    try:
        import torch
    except ImportError as exc:  # pragma: no cover - environment diagnostic
        raise RuntimeError("torch is required to read committed SBI metadata") from exc
    return torch.load(path, map_location="cpu", weights_only=False)


def _split_md_row(line: str) -> list[str]:
    sentinel = "\x00PIPE\x00"
    return [cell.strip().replace(sentinel, r"\|") for cell in
            line.strip().strip("|").replace(r"\|", sentinel).split("|")]


def _markdown_table_after(path: Path, heading: str) -> tuple[list[str], list[list[str]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    try:
        start = next(i for i, line in enumerate(lines)
                     if line.strip() == heading)
    except StopIteration as exc:
        raise ValueError(f"heading not found in {_rel(path)}: {heading}") from exc
    table_lines: list[str] = []
    for line in lines[start + 1:]:
        if line.startswith("|"):
            table_lines.append(line)
        elif table_lines:
            break
    if len(table_lines) < 3:
        raise ValueError(f"table not found after {heading} in {_rel(path)}")
    header = _split_md_row(table_lines[0])
    rows = [_split_md_row(line) for line in table_lines[2:]]
    return header, rows


def _index_row(artifact_name: str) -> dict[str, str]:
    lines = INDEX_PATH.read_text(encoding="utf-8").splitlines()
    for line in lines:
        if line.startswith(f"| {artifact_name} |"):
            cells = _split_md_row(line)
            if len(cells) == 8:
                keys = ("artifact", "config", "config_hash", "git_sha",
                        "seed", "n_train", "gates", "deployed")
            elif len(cells) == 7:
                keys = ("artifact", "config", "config_hash", "seed",
                        "n_train", "gates", "deployed")
            else:
                raise ValueError(
                    f"unexpected INDEX row width for {artifact_name}: {len(cells)}")
            return dict(zip(keys, cells))
    raise ValueError(f"artifact row not found in INDEX.md: {artifact_name}")


def _config_sentences(config_path: Path, config: dict[str, Any]) -> list[str]:
    param_specs = []
    for name, spec in config["param_space"].items():
        prior = spec.get("prior_type", "unspecified")
        bounds = spec.get("bounds")
        if bounds is None:
            param_specs.append(f"{name} ({prior})")
        else:
            param_specs.append(
                f"{name} ({prior} [{_fmt(bounds[0])}, {_fmt(bounds[1])}])")
    observables = [
        f"{name}={_fmt(value[0])} ± {_fmt(value[1])}"
        for name, value in config["observables"].items()
    ]
    ptr = _rel(config_path)
    return [
        (f"The {config['bodyname']} campaign sampled {len(param_specs)} "
         f"parameters: {', '.join(param_specs)}.\n"
         f"{_source_comment(ptr + '#/param_space')}"),
        (f"The likelihood used {len(observables)} configured observables: "
         f"{', '.join(observables)}.\n"
         f"{_source_comment(ptr + '#/observables')}"),
    ]


def _cache_sentence(spec: dict[str, Any], config: dict[str, Any]) -> str:
    cache_path = REPO / config["structure_cache_path"]
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"numpy\.core\.numeric is deprecated",
            category=DeprecationWarning,
        )
        with cache_path.open("rb") as stream:
            cache = pickle.load(stream)
    if not isinstance(cache, dict):
        raise ValueError(f"cache payload is not a dict: {_rel(cache_path)}")

    axes = []
    for key, values in cache.items():
        if key.endswith("_grid"):
            vals = list(values)
            axes.append(
                f"{key}={len(vals)} points [{_fmt(float(vals[0]))}, "
                f"{_fmt(float(vals[-1]))}]"
            )
    n_structures = len(cache.get("structures", []))
    composition = cache.get("ocean_comp")
    if composition is None:
        composition = (config.get("metadata") or {}).get("ocean_composition")
    if composition is None:
        composition = (config.get("sampler_settings") or {}).get(
            "phase_stability", {}).get("enforce", "not recorded")

    flags: list[str] = []
    for key, value in cache.items():
        if isinstance(value, bool):
            flags.append(f"{key}={value}")

    manifest_path = (REPO / spec["cache_manifest"]
                     if spec["cache_manifest"] else None)
    manifest = None
    if manifest_path:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for key, value in manifest.items():
            if isinstance(value, bool) and key not in {"validate"}:
                item = f"{key}={value}"
                if item not in flags:
                    flags.append(item)

    cache_sha = _sha256(cache_path)
    if manifest is not None:
        expected = str(manifest[spec["cache_sha_key"]])
        if cache_sha != expected:
            raise ValueError(
                f"cache SHA mismatch for {_rel(cache_path)}: "
                f"{cache_sha} != {expected}")

    schema = cache.get("schema_version")
    schema_text = f", schema {schema}" if schema is not None else ""
    flag_text = ", ".join(flags) if flags else "none recorded"
    pointers = [f"{_rel(cache_path)}#cache-metadata"]
    if manifest_path:
        pointers.append(
            f"{_rel(manifest_path)}#/{spec['cache_sha_key']}")
    return (
        f"The committed {composition} cache {_rel(cache_path)} contains "
        f"{n_structures} structure entries on {', '.join(axes)}{schema_text}, "
        f"with builder flags {flag_text} and SHA-256 `{cache_sha}`.\n"
        f"{_source_comment(*pointers)}"
    )


def _artifact_sentence(spec: dict[str, Any], config: dict[str, Any],
                       artifact: dict[str, Any]) -> str:
    artifact_path = ARTIFACT_DIR / spec["artifact"]
    param_names = list(config["param_space"])
    obs_names = list(config["observables"])
    if artifact.get("param_names") != param_names:
        raise ValueError(f"parameter-name mismatch in {spec['artifact']}")
    if artifact.get("obs_names") != obs_names:
        raise ValueError(f"observable-name mismatch in {spec['artifact']}")
    return (
        f"The committed artifact {_rel(artifact_path)} retained "
        f"{int(artifact['n_train_effective']):,} training simulations and "
        f"used seed {artifact['seed']}, architecture "
        f"{artifact['density_estimator']}, torch {artifact['torch_version']}, "
        f"and sbi {artifact['sbi_version']}.\n"
        f"{_source_comment(_rel(artifact_path) + '#artifact-metadata')}"
    )


def _gate_sentences(spec: dict[str, Any], index: dict[str, str]) -> list[str]:
    source = REPO / spec["gate_source"]
    kind = spec["gate_kind"]
    ptr = _rel(source)
    if kind == "ratification":
        header, rows = _markdown_table_after(source, "## Gate outcomes")
        verdict_col = "Raw verdict" if "Raw verdict" in header else "Verdict"
        gate_i, verdict_i = header.index("Gate"), header.index(verdict_col)
        pairs = [f"{row[gate_i]} = {row[verdict_i]}" for row in rows]
        out = [
            f"The committed raw gate outcomes are {'; '.join(pairs)}.\n"
            f"{_source_comment(ptr + '#gate-outcomes')}"
        ]
        if "FAIL-ADJUDICATED-ACCEPTABLE" in source.read_text(encoding="utf-8"):
            out.append(
                "The committed adjudication vocabulary includes "
                "`FAIL-ADJUDICATED-ACCEPTABLE`.\n"
                f"{_source_comment(ptr + '#adjudications-and-required-caveats')}"
            )
        return out
    if kind == "reviewer_json":
        data = json.loads(source.read_text(encoding="utf-8"))
        adjudication = data["reviewer_adjudication_2026_08_06"]
        return [
            f"The committed gate adjudication states: “{adjudication[key]}”.\n"
            f"{_source_comment(ptr + '#/reviewer_adjudication_2026_08_06/' + key)}"
            for key in ("verdict", "sbc_baseline", "crosscheck_baseline",
                        "limits_baseline")
        ]
    if kind == "overall_heading":
        text = source.read_text(encoding="utf-8")
        match = re.search(r"^## Overall: \*\*(.+?)\*\* \((.+?)\)$",
                          text, flags=re.MULTILINE)
        if not match:
            raise ValueError(f"overall gate heading not found in {ptr}")
        return [
            f"The committed overall gate outcome is {match.group(1)} "
            f"({match.group(2)}).\n"
            f"{_source_comment(ptr + '#overall')}"
        ]
    if kind == "phase_table":
        header, rows = _markdown_table_after(source, "## Gate results")
        gate_i, verdict_i = header.index("Gate"), header.index("Verdict")
        pairs = [f"{row[gate_i]} = {row[verdict_i]}" for row in rows]
        return [
            f"The committed raw gate outcomes are {'; '.join(pairs)}.\n"
            f"{_source_comment(ptr + '#gate-results')}"
        ]
    if kind == "index_excerpt":
        excerpt = spec["gate_excerpt"]
        if excerpt not in index["gates"]:
            raise ValueError(f"gate excerpt is no longer present for {spec['artifact']}")
        return [
            f"The committed gate record states: “{excerpt}”.\n"
            f"{_source_comment(ptr + '#deployed-artifacts')}"
        ]
    raise ValueError(f"unknown gate kind: {kind}")


def _deployment_sentences(spec: dict[str, Any], index: dict[str, str]) -> list[str]:
    deployed = index["deployed"]
    match = re.search(r"\*\*(.+?)\*\*", deployed)
    excerpt = (match.group(1) if match else deployed).rstrip(". ")
    out = [
        f"The artifact index records deployment state as: “{excerpt}”.\n"
        f"{_source_comment(_rel(INDEX_PATH) + '#artifact-' + spec['artifact'])}"
    ]
    extra = spec.get("deployment_extra")
    if extra:
        if extra not in deployed:
            raise ValueError(
                f"deployment adjudication changed for {spec['artifact']}")
        out.append(
            f"“{extra}”\n"
            f"{_source_comment(_rel(INDEX_PATH) + '#artifact-' + spec['artifact'])}"
        )
    return out


def render_campaign(key: str) -> str:
    spec = CAMPAIGNS[key]
    config_path = CONFIG_DIR / spec["config"]
    artifact_path = ARTIFACT_DIR / spec["artifact"]
    config = json.loads(config_path.read_text(encoding="utf-8"))
    artifact = _load_artifact(artifact_path)
    index = _index_row(spec["artifact"])

    if artifact["config_hash"] not in index["config_hash"]:
        raise ValueError(f"artifact config hash absent from INDEX for {key}")
    if index["config"] != spec["config"]:
        raise ValueError(f"INDEX config mismatch for {key}")

    blocks = [
        f"# Methods snippet: {spec['title']}",
        "<!-- generated by plans/scripts/generate_methods_snippet.py; do not edit -->",
        *_config_sentences(config_path, config),
        _cache_sentence(spec, config),
        _artifact_sentence(spec, config, artifact),
        *_gate_sentences(spec, index),
        *_deployment_sentences(spec, index),
    ]
    return "\n\n".join(blocks) + "\n"


def output_path(key: str) -> Path:
    return REPO / CAMPAIGNS[key]["report_dir"] / "methods_snippet.md"


def generate(keys: list[str]) -> None:
    for key in keys:
        path = output_path(key)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(render_campaign(key), encoding="utf-8")
        print(f"generated {key}: {_rel(path)}")


def verify(keys: list[str]) -> None:
    try:
        import markdown
    except ImportError as exc:  # pragma: no cover - environment diagnostic
        raise RuntimeError("python-markdown is required for render verification") from exc

    for key in keys:
        path = output_path(key)
        expected = render_campaign(key)
        if not path.is_file():
            raise ValueError(f"missing generated snippet: {_rel(path)}")
        actual = path.read_text(encoding="utf-8")
        if actual != expected:
            raise ValueError(f"generated snippet is stale: {_rel(path)}")
        html = markdown.markdown(actual, extensions=["tables"])
        if "<h1>Methods snippet:" not in html or "<p>" not in html:
            raise ValueError(f"Markdown render check failed: {_rel(path)}")
        sentence_blocks = [block for block in actual.split("\n\n")
                           if block and not block.startswith("#")
                           and not block.startswith("<!-- generated")]
        if not all("<!-- sources:" in block for block in sentence_blocks):
            raise ValueError(f"sentence lacks source pointer: {_rel(path)}")
        print(f"verified {key}: source values, hashes, names, and Markdown render")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("campaign", nargs="?", default="all",
                        choices=[*CAMPAIGNS, "all"])
    parser.add_argument("--verify", action="store_true",
                        help="compare generated files with committed sources")
    args = parser.parse_args()
    keys = list(CAMPAIGNS) if args.campaign == "all" else [args.campaign]
    if args.verify:
        verify(keys)
    else:
        generate(keys)
    return 0


if __name__ == "__main__":
    sys.exit(main())
