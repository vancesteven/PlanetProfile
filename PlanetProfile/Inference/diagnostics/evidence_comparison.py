"""
evidence_comparison.py — Pairwise log-Bayes factor comparison for MCMC runs.

Computes ln B₁₂ = log_Z₁ - log_Z₂ for matched pairs of MCMC runs and
classifies the strength of evidence per Kass & Raftery (1995).

Kass & Raftery (1995) classification (on |ln B|):
    |ln B| < 1     : not worth more than a bare mention
    1 ≤ |ln B| < 3 : positive evidence
    3 ≤ |ln B| < 5 : strong evidence
    |ln B| ≥ 5     : very strong evidence

MC error on ln B:
    σ(ln B) = sqrt(σ_Z1² + σ_Z2²)

A pair is flagged "MC error too large to resolve" when σ(ln B) > 1.0.

Usage
-----
    python3 PlanetProfile/Inference/diagnostics/evidence_comparison.py
        [--results-dir DIR]
        [--pairs-file FILE]   # optional YAML/JSON with custom pairs
        [--output PATH]

All Titan comparison pairs are built in and run by default.
"""

from __future__ import annotations

import argparse
import datetime
import json
import math
import os
import pickle
import sys
from pathlib import Path
from typing import Any, Optional

import numpy as np

# ---------------------------------------------------------------------------
# Project root on sys.path so InferenceResult can be unpickled
# ---------------------------------------------------------------------------
os.environ.setdefault("NUMBA_DISABLE_JIT", "1")

_this_file = Path(__file__).resolve()
_project_root = _this_file.parent.parent.parent.parent  # …/planetprofile-genai
if str(_project_root) not in sys.path:
    sys.path.insert(0, str(_project_root))

# ---------------------------------------------------------------------------
# Built-in comparison pairs
# Format: (body, run_name_1, run_name_2, label)
# Run names must match subdirectory names under mcmc_results/<body>/
# ---------------------------------------------------------------------------

TITAN_PAIRS: list[tuple[str, str, str, str]] = [
    (
        "Titan",
        "Test43_andrade_arrhenius_no_ocean",
        "Test44_maxwell_arrhenius_ocean",
        "Andrade no-ocean vs Maxwell ocean",
    ),
    (
        "Titan",
        "Test45_maxwell_hybrid_hydro",
        "Test46_allice",
        "Maxwell hybrid-hydro vs all-ice",
    ),
    (
        "Titan",
        "Test43_andrade_arrhenius_no_ocean",
        "Test46_allice",
        "Andrade no-ocean vs all-ice",
    ),
    (
        "Titan",
        "Test44_maxwell_arrhenius_ocean",
        "Test45_maxwell_hybrid_hydro",
        "Maxwell ocean vs Maxwell hybrid-hydro",
    ),
]

# Additional pairs can be added here or via --pairs-file
DEFAULT_PAIRS = TITAN_PAIRS


# ---------------------------------------------------------------------------
# Kass & Raftery (1995) classification
# ---------------------------------------------------------------------------

def kr_grade(abs_ln_b: float) -> str:
    """Return K&R (1995) grade string for |ln B|."""
    if abs_ln_b < 1.0:
        return "not worth mentioning"
    elif abs_ln_b < 3.0:
        return "positive"
    elif abs_ln_b < 5.0:
        return "strong"
    else:
        return "very strong"


# ---------------------------------------------------------------------------
# pkl loading and log_Z extraction
# ---------------------------------------------------------------------------

def _safe_get(obj: Any, *keys, default=None):
    """Traverse nested attrs / dict keys, returning default on any failure."""
    cur = obj
    for k in keys:
        try:
            if isinstance(cur, dict):
                cur = cur[k]
            else:
                cur = getattr(cur, k)
        except (KeyError, AttributeError, TypeError):
            return default
    return cur


def find_result_pkl(run_dir: Path) -> Optional[Path]:
    """
    Find the best candidate result pkl in run_dir.

    Preference order:
      1. *result*.pkl files that are NOT structure_cache or conflicted copies
      2. *mcmc*.pkl files (legacy format)
      3. Any *.pkl that is not a structure cache or conflicted copy
    """
    if not run_dir.exists():
        return None

    def _is_valid(p: Path) -> bool:
        name = p.name.lower()
        return "conflicted" not in name and "structure" not in name

    # Priority 1: *result*.pkl
    candidates = sorted(p for p in run_dir.glob("*result*.pkl") if _is_valid(p))
    if candidates:
        # Prefer the most recent by mtime
        return max(candidates, key=lambda p: p.stat().st_mtime)

    # Priority 2: *mcmc*.pkl
    candidates = sorted(p for p in run_dir.glob("*mcmc*.pkl") if _is_valid(p))
    if candidates:
        return max(candidates, key=lambda p: p.stat().st_mtime)

    # Priority 3: any valid pkl
    candidates = sorted(p for p in run_dir.glob("*.pkl") if _is_valid(p))
    if candidates:
        return max(candidates, key=lambda p: p.stat().st_mtime)

    return None


def load_log_z(pkl_path: Path) -> tuple[Optional[float], Optional[float]]:
    """
    Load a pkl and return (log_Z, log_Z_err).

    Returns (None, None) when the pkl cannot be read or log_Z is not stored.
    """
    try:
        with open(pkl_path, "rb") as fh:
            result = pickle.load(fh)
    except Exception as exc:
        raise RuntimeError(f"Cannot load pkl {pkl_path}: {exc}") from exc

    log_Z: Optional[float] = None
    log_Z_err: Optional[float] = None

    # Branch A: modern InferenceResult (has .metadata dict attribute)
    if hasattr(result, "metadata") and isinstance(result.metadata, dict):
        meta = result.metadata
        raw_Z = meta.get("log_Z") or meta.get("log_evidence")
        raw_Zerr = meta.get("log_Z_err") or meta.get("log_evidence_err")
        if raw_Z is not None:
            try:
                log_Z = float(raw_Z)
            except (TypeError, ValueError):
                log_Z = None
        if raw_Zerr is not None:
            try:
                log_Z_err = float(raw_Zerr)
            except (TypeError, ValueError):
                log_Z_err = None

    # Branch B: legacy plain-dict format
    elif isinstance(result, dict):
        raw_Z = result.get("log_Z") or result.get("log_evidence")
        raw_Zerr = result.get("log_Z_err") or result.get("log_evidence_err")
        # Also check nested metadata key
        if raw_Z is None:
            raw_Z = _safe_get(result, "metadata", "log_Z")
        if raw_Zerr is None:
            raw_Zerr = _safe_get(result, "metadata", "log_Z_err")
        if raw_Z is not None:
            try:
                log_Z = float(raw_Z)
            except (TypeError, ValueError):
                log_Z = None
        if raw_Zerr is not None:
            try:
                log_Z_err = float(raw_Zerr)
            except (TypeError, ValueError):
                log_Z_err = None

    return log_Z, log_Z_err


# ---------------------------------------------------------------------------
# Per-pair computation
# ---------------------------------------------------------------------------

def compute_pair(
    body: str,
    run1: str,
    run2: str,
    label: str,
    results_dir: Path,
) -> dict:
    """
    Compute ln B₁₂ for a single (run1, run2) pair.

    Returns a dict with fields:
        body, run1, run2, label,
        log_Z_1, log_Z_err_1, pkl_1, found_1,
        log_Z_2, log_Z_err_2, pkl_2, found_2,
        ln_B, sigma_lnB, kr_grade,
        skip_reason, flags
    """
    out: dict[str, Any] = dict(
        body=body,
        run1=run1,
        run2=run2,
        label=label,
        log_Z_1=None,
        log_Z_err_1=None,
        pkl_1=None,
        found_1=False,
        log_Z_2=None,
        log_Z_err_2=None,
        pkl_2=None,
        found_2=False,
        ln_B=None,
        sigma_lnB=None,
        kr_grade=None,
        skip_reason=None,
        flags=[],
    )

    body_dir = results_dir / body

    # --- Load run 1 ---
    run1_dir = body_dir / run1
    pkl1 = find_result_pkl(run1_dir)
    if pkl1 is None:
        out["skip_reason"] = f"pkl not found for {run1}"
        return out
    out["pkl_1"] = str(pkl1)
    out["found_1"] = True

    try:
        log_Z_1, log_Z_err_1 = load_log_z(pkl1)
    except RuntimeError as exc:
        out["skip_reason"] = str(exc)
        return out
    out["log_Z_1"] = log_Z_1
    out["log_Z_err_1"] = log_Z_err_1

    # --- Load run 2 ---
    run2_dir = body_dir / run2
    pkl2 = find_result_pkl(run2_dir)
    if pkl2 is None:
        out["skip_reason"] = f"pkl not found for {run2}"
        return out
    out["pkl_2"] = str(pkl2)
    out["found_2"] = True

    try:
        log_Z_2, log_Z_err_2 = load_log_z(pkl2)
    except RuntimeError as exc:
        out["skip_reason"] = str(exc)
        return out
    out["log_Z_2"] = log_Z_2
    out["log_Z_err_2"] = log_Z_err_2

    # --- Check for missing log_Z ---
    missing = []
    if log_Z_1 is None:
        missing.append(run1)
    if log_Z_2 is None:
        missing.append(run2)
    if missing:
        out["skip_reason"] = (
            "log_Z not available (pre-fix run — rerun to populate): "
            + ", ".join(missing)
        )
        return out

    # --- Compute ln B₁₂ ---
    ln_B = log_Z_1 - log_Z_2
    out["ln_B"] = ln_B

    # --- Compute σ(ln B) ---
    if log_Z_err_1 is not None and log_Z_err_2 is not None:
        sigma = math.sqrt(log_Z_err_1 ** 2 + log_Z_err_2 ** 2)
        out["sigma_lnB"] = sigma
        if sigma > 1.0:
            out["flags"].append("MC error too large to resolve")
    else:
        out["sigma_lnB"] = None
        out["flags"].append("σ(ln B) unknown — log_Z_err missing")

    # --- K&R grade ---
    out["kr_grade"] = kr_grade(abs(ln_B))

    return out


# ---------------------------------------------------------------------------
# Formatting
# ---------------------------------------------------------------------------

def _fmt_float(val: Optional[float], fmt: str = ".3f") -> str:
    if val is None:
        return "—"
    try:
        return format(val, fmt)
    except (TypeError, ValueError):
        return str(val)


def build_markdown_table(pairs: list[dict]) -> str:
    """Render comparison results as a Markdown table."""
    header = (
        "| Model 1 | Model 2 | ln B₁₂ | σ(ln B) | K&R grade | Flag |"
    )
    divider = "|---|---|---:|---:|---|---|"
    lines = [header, divider]

    for p in pairs:
        if p["skip_reason"]:
            # One row with skip reason spanning the numeric columns
            ln_b_str = "—"
            sigma_str = "—"
            grade_str = "—"
            flag_str = p["skip_reason"]
        else:
            ln_b_str = _fmt_float(p["ln_B"], "+.3f")
            sigma_str = _fmt_float(p["sigma_lnB"], ".3f")
            grade_str = p["kr_grade"] or "—"
            flag_str = "; ".join(p["flags"]) if p["flags"] else "ok"

        lines.append(
            f"| {p['run1']} "
            f"| {p['run2']} "
            f"| {ln_b_str} "
            f"| {sigma_str} "
            f"| {grade_str} "
            f"| {flag_str} |"
        )

    return "\n".join(lines)


def build_report(pairs: list[dict], results_dir: Path) -> str:
    """Build the full Markdown evidence comparison report."""
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    n_available = sum(
        1 for p in pairs
        if p["ln_B"] is not None
    )
    n_skipped = len(pairs) - n_available

    lines = [
        "# MCMC Evidence Comparison Report",
        "",
        f"Generated: {now}  ",
        f"Results dir: `{results_dir}`  ",
        f"Pairs evaluated: {len(pairs)}  ",
        f"Pairs with log_Z available: {n_available}  ",
        f"Pairs skipped (log_Z unavailable): {n_skipped}  ",
        "",
        "## Kass & Raftery (1995) classification",
        "",
        "| |ln B₁₂| | Grade |",
        "|---|---|",
        "| < 1 | not worth mentioning |",
        "| 1 – 3 | positive evidence |",
        "| 3 – 5 | strong evidence |",
        "| ≥ 5 | very strong evidence |",
        "",
        "ln B₁₂ = log Z₁ − log Z₂ > 0 favours Model 1.  ",
        "σ(ln B) = √(σ_Z1² + σ_Z2²).  ",
        "Pairs with σ(ln B) > 1.0 are flagged as unresolvable.  ",
        "",
        "## Pairwise comparisons",
        "",
        build_markdown_table(pairs),
        "",
    ]

    # Per-pair detail block for skipped runs
    skipped = [p for p in pairs if p["skip_reason"]]
    if skipped:
        lines += [
            "## Skipped pairs — detail",
            "",
        ]
        for p in skipped:
            lines.append(
                f"- **{p['label']}** ({p['run1']} vs {p['run2']}):  "
            )
            lines.append(f"  {p['skip_reason']}")
            if p.get("pkl_1"):
                lines.append(f"  pkl 1: `{p['pkl_1']}`")
            if p.get("pkl_2"):
                lines.append(f"  pkl 2: `{p['pkl_2']}`")
            lines.append("")

    # Footer
    lines += [
        "---",
        "",
        "> **Note:** All comparisons require production-quality runs "
        "(n_eff ≥ 500, ESS ≥ 500).  ",
        "> Smoke runs (n_eff = 100) are listed for pipeline validation only "
        "and **must not** be used for Bayes-factor conclusions.  ",
        "> Runs flagged ⚠ pre-bugfix (before commit ca1b600, 2026-05-23) "
        "have unreliable importance weights.  ",
        "",
    ]

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def load_pairs_from_file(path: Path) -> list[tuple[str, str, str, str]]:
    """
    Load comparison pairs from a JSON or YAML file.

    Expected format (JSON example):
        [
            {"body": "Titan", "run1": "TestXX_...", "run2": "TestYY_...",
             "label": "Human-readable description"},
            ...
        ]
    """
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() in (".yaml", ".yml"):
        try:
            import yaml  # type: ignore
            data = yaml.safe_load(text)
        except ImportError:
            raise ImportError(
                "PyYAML is required for YAML pair files. "
                "Install with: pip install pyyaml"
            )
    else:
        data = json.loads(text)

    pairs = []
    for item in data:
        pairs.append((
            item["body"],
            item["run1"],
            item["run2"],
            item.get("label", f"{item['run1']} vs {item['run2']}"),
        ))
    return pairs


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Compute pairwise log-Bayes factors (ln B₁₂) between MCMC runs "
            "and write a Markdown evidence comparison report."
        )
    )
    parser.add_argument(
        "--results-dir",
        default=None,
        help=(
            "Root directory for MCMC results (expected layout: "
            "<results-dir>/<Body>/<RunName>/*.pkl). "
            "Default: PlanetProfile/Test/mcmc_results relative to project root."
        ),
    )
    parser.add_argument(
        "--pairs-file",
        default=None,
        help=(
            "Path to a JSON or YAML file defining custom comparison pairs. "
            "If not provided, the built-in Titan pairs are used."
        ),
    )
    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Path to write the evidence_comparison_report.md. "
            "Default: PlanetProfile/Inference/diagnostics/evidence_comparison_report.md"
        ),
    )
    args = parser.parse_args()

    # Resolve paths
    if args.results_dir:
        results_dir = Path(args.results_dir).resolve()
    else:
        results_dir = _project_root / "PlanetProfile" / "Test" / "mcmc_results"

    if args.output:
        output_path = Path(args.output).resolve()
    else:
        output_path = _this_file.parent / "evidence_comparison_report.md"

    if not results_dir.exists():
        print(
            f"ERROR: results directory does not exist: {results_dir}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Load pairs
    if args.pairs_file:
        pairs_file = Path(args.pairs_file)
        print(f"Loading comparison pairs from {pairs_file}", flush=True)
        raw_pairs = load_pairs_from_file(pairs_file)
    else:
        raw_pairs = DEFAULT_PAIRS
        print(
            f"Using {len(raw_pairs)} built-in Titan comparison pair(s).",
            flush=True,
        )

    # Process each pair
    print(f"Results directory: {results_dir}\n", flush=True)
    pair_results = []
    for body, run1, run2, label in raw_pairs:
        print(f"  Pair: {label}", flush=True)
        print(f"    Model 1: {run1}", flush=True)
        print(f"    Model 2: {run2}", flush=True)
        result = compute_pair(body, run1, run2, label, results_dir)
        pair_results.append(result)

        if result["skip_reason"]:
            print(f"    SKIP: {result['skip_reason']}", flush=True)
        else:
            flags = ", ".join(result["flags"]) if result["flags"] else "ok"
            print(
                f"    ln B₁₂ = {result['ln_B']:+.4f}  "
                f"σ = {_fmt_float(result['sigma_lnB'], '.4f')}  "
                f"grade = {result['kr_grade']}  "
                f"flags = {flags}",
                flush=True,
            )
        print("", flush=True)

    # Build and write report
    report = build_report(pair_results, results_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(report, encoding="utf-8")
    print(f"Evidence comparison report written to: {output_path}\n")
    print("=" * 70)
    print(report)


if __name__ == "__main__":
    main()
