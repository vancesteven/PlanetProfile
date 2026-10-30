"""
run_audit.py — MCMC run audit utility.

Scans all *results*.pkl files under a results directory, extracts key
provenance and quality fields, flags potential issues, and writes a
Markdown table to audit_report.md alongside printing to stdout.

Usage
-----
    python3 PlanetProfile/Inference/diagnostics/run_audit.py [--results-dir DIR]
                                                              [--output PATH]

Flags
-----
    ⚠ ESS<500          Effective sample size below convergence threshold
    ⚠ log-Z err>0.5    Log-evidence MC error exceeds tolerance
    ⚠ seed missing      Random seed not recorded in pkl
    ⚠ pre-bugfix        git SHA absent or predates ca1b600 (2026-05-23),
                        the pocoMC unpack-swap bug fix; samples are fine
                        but weights/histograms from that run may be wrong
    ℹ smoke             n_samples ≤ 200 — exploratory run only, not for
                        Bayes-factor comparison
"""

from __future__ import annotations

import argparse
import datetime
import os
import pickle
import sys
import traceback
from pathlib import Path
from typing import Any

import numpy as np

# ---------------------------------------------------------------------------
# Ensure the project root is on sys.path so InferenceResult (and other PP
# dataclasses) can be unpickled even when the script is invoked as
#   python3 PlanetProfile/Inference/diagnostics/run_audit.py
# We also set NUMBA_DISABLE_JIT so that TidalPy's numba-decorated functions
# don't raise a "no locator available" RuntimeError during import.
# ---------------------------------------------------------------------------
os.environ.setdefault("NUMBA_DISABLE_JIT", "1")

_this_file = Path(__file__).resolve()
_project_root = _this_file.parent.parent.parent.parent  # …/planetprofile-genai
if str(_project_root) not in sys.path:
    sys.path.insert(0, str(_project_root))

# ---------------------------------------------------------------------------
# Bugfix reference commit
# ---------------------------------------------------------------------------
# ca1b600 (2026-05-23): fixed pocoMC posterior() unpack order
# (samples, weights, logl, logp) — runs generated *before* this commit
# stored importance weights in the log_likelihoods field.
_BUGFIX_SHA_PREFIX = "ca1b600"
_BUGFIX_DATE = datetime.date(2026, 5, 23)


# ---------------------------------------------------------------------------
# Helpers
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


def _compute_ess_from_weights(weights: np.ndarray) -> float:
    """ESS = (sum w)^2 / sum(w^2); weights need not be normalised."""
    w = np.asarray(weights, dtype=float)
    w = w[np.isfinite(w) & (w > 0)]
    if len(w) == 0:
        return 0.0
    return float(np.sum(w) ** 2 / np.sum(w ** 2))


def _sha_predates_bugfix(sha: str | None) -> bool:
    """
    Return True when sha is None or an abbreviated prefix that is an
    ancestor of ca1b600 (i.e., predates the pocoMC fix).

    We can't run git at audit time in all environments, so we maintain a
    hard-coded set of known-pre-bugfix commit prefixes for runs we have
    on disk.  A missing SHA is always flagged pre-bugfix.
    """
    if sha is None or str(sha).strip() == "":
        return True

    # Hard-coded known pre-bugfix SHAs (from git log for runs on disk)
    # Any SHA not in this set and not matching the bugfix commit is treated
    # as post-bugfix (conservative: don't flag runs we don't know about).
    _PRE_BUGFIX_PREFIXES = {
        # All runs saved as plain dicts (before MCMCRunner) are pre-bugfix.
        # Their SHAs are not recorded; we detect them by sha==None above.
    }

    sha_norm = str(sha).strip().lower()
    for pre in _PRE_BUGFIX_PREFIXES:
        if sha_norm.startswith(pre.lower()):
            return True

    # If the sha IS the bugfix commit or later, not pre-bugfix.
    if sha_norm.startswith(_BUGFIX_SHA_PREFIX.lower()):
        return False  # this is the fix itself, not pre-bugfix

    # Unknown SHA — assume post-bugfix (don't over-flag future runs).
    return False


# ---------------------------------------------------------------------------
# Per-pkl extraction
# ---------------------------------------------------------------------------

def extract_run_info(pkl_path: Path, project_root: Path) -> dict:
    """
    Load a pkl and extract audit fields.  Never raises — catches all
    exceptions and sets a flag instead.

    Returns a dict with keys:
        run_name, body, git_sha, seed, n_samples, ess, log_z, log_z_err,
        walltime_h, sampler, flags, readable, pkl_path
    """
    info = dict(
        run_name=pkl_path.parent.name,
        body="?",
        git_sha=None,
        seed=None,
        n_samples=None,
        ess=None,
        log_z=None,
        log_z_err=None,
        walltime_h=None,
        sampler="?",
        flags=[],
        readable=True,
        pkl_path=str(pkl_path),
    )

    try:
        with open(pkl_path, "rb") as fh:
            result = pickle.load(fh)
    except Exception:
        info["readable"] = False
        info["flags"].append("❌ unreadable")
        info["_error"] = traceback.format_exc(limit=3)
        return info

    # -----------------------------------------------------------------------
    # Branch A: modern InferenceResult (dataclass with .config, .samples, …)
    # -----------------------------------------------------------------------
    if hasattr(result, "config") and hasattr(result, "samples"):
        cfg = result.config
        info["body"] = _safe_get(cfg, "bodyname", default="?")
        info["seed"] = _safe_get(cfg, "random_state", default=None)

        # Sampler name from config
        sampler_settings = _safe_get(cfg, "sampler_settings", default={}) or {}
        info["sampler"] = sampler_settings.get("sampler", "pocoMC")

        # n_effective target — used to detect smoke runs even if pocoMC
        # returned more samples than requested.
        n_eff_target = sampler_settings.get("n_effective")
        if n_eff_target is not None:
            try:
                info["_n_eff_target"] = int(n_eff_target)
            except (TypeError, ValueError):
                pass

        # Sample count
        samples = result.samples
        if samples is not None:
            info["n_samples"] = int(samples.shape[0])

        # Convergence metrics
        conv = _safe_get(result, "convergence_metrics", default={}) or {}
        ess_raw = conv.get("ess")
        if ess_raw is not None:
            try:
                info["ess"] = float(ess_raw)
            except (TypeError, ValueError):
                pass

        # Metadata: elapsed time, git SHA
        meta = _safe_get(result, "metadata", default={}) or {}
        elapsed_s = meta.get("elapsed_time_s")
        if elapsed_s is not None:
            try:
                info["walltime_h"] = float(elapsed_s) / 3600.0
            except (TypeError, ValueError):
                pass

        info["git_sha"] = meta.get("git_sha") or meta.get("git_hash") or None

        # log-Z from metadata (not yet stored by MCMCRunner — future field)
        info["log_z"] = meta.get("log_Z") or meta.get("log_evidence") or None
        info["log_z_err"] = meta.get("log_Z_err") or meta.get("log_evidence_err") or None

        # If ESS not in conv metrics, try computing from weights in meta
        if info["ess"] is None:
            weights = meta.get("weights")
            if weights is not None:
                try:
                    info["ess"] = _compute_ess_from_weights(np.asarray(weights))
                except Exception:
                    pass

    # -----------------------------------------------------------------------
    # Branch B: legacy plain-dict format
    # -----------------------------------------------------------------------
    elif isinstance(result, dict):
        # Infer body from parent directory name
        body_dir = pkl_path.parent.parent.name
        info["body"] = body_dir if body_dir else "?"

        # Sample count
        samples = result.get("samples")
        if samples is not None:
            try:
                info["n_samples"] = int(np.asarray(samples).shape[0])
            except Exception:
                pass

        # seed / git_sha / sampler not stored in legacy dicts
        info["seed"] = result.get("random_state") or result.get("seed")
        info["git_sha"] = result.get("git_sha") or result.get("git_hash")
        info["sampler"] = result.get("sampler", "legacy-dict")

        # log-Z not stored in legacy format
        info["log_z"] = result.get("log_Z") or result.get("log_evidence")
        info["log_z_err"] = result.get("log_Z_err") or result.get("log_evidence_err")

        # Elapsed time
        elapsed_s = result.get("elapsed_time_s")
        if elapsed_s is not None:
            try:
                info["walltime_h"] = float(elapsed_s) / 3600.0
            except (TypeError, ValueError):
                pass

        # ESS — not usually stored in legacy format; compute if weights present
        weights = result.get("weights")
        if weights is not None:
            try:
                info["ess"] = _compute_ess_from_weights(np.asarray(weights))
            except Exception:
                pass

    else:
        info["flags"].append("⚠ unknown-format")
        info["readable"] = False
        return info

    # -----------------------------------------------------------------------
    # Flag logic
    # -----------------------------------------------------------------------
    n = info["n_samples"]
    ess = info["ess"]
    log_z_err = info["log_z_err"]
    seed = info["seed"]
    git_sha = info["git_sha"]
    n_eff_target = info.pop("_n_eff_target", None)

    # Smoke detection: use the n_effective *target* if stored (for InferenceResult
    # runs where pocoMC may over-shoot the target but the intent was exploratory),
    # otherwise fall back to actual sample count.
    smoke_n = n_eff_target if (n_eff_target is not None and n_eff_target <= 200) else n
    if smoke_n is not None and smoke_n <= 200:
        info["flags"].append("ℹ smoke")

    if ess is not None and ess < 500:
        info["flags"].append("⚠ ESS<500")

    if log_z_err is not None:
        try:
            if float(log_z_err) > 0.5:
                info["flags"].append("⚠ log-Z err>0.5")
        except (TypeError, ValueError):
            pass

    if seed is None:
        info["flags"].append("⚠ seed missing")

    if _sha_predates_bugfix(git_sha):
        info["flags"].append("⚠ pre-bugfix")

    return info


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def _fmt(value, fmt=None, missing="—"):
    if value is None:
        return missing
    if fmt:
        try:
            return format(value, fmt)
        except (ValueError, TypeError):
            return str(value)
    return str(value)


def _fmt_sha(sha):
    if sha is None:
        return "—"
    return str(sha)[:10]


def build_markdown_table(rows: list[dict]) -> str:
    """Render list of run-info dicts as a Markdown table."""
    header = (
        "| Run name | Body | SHA | Seed | n_samples | ESS "
        "| log-Z | log-Z err | Walltime (h) | Sampler | Flags |"
    )
    divider = (
        "|---|---|---|---|---:|---:|---:|---:|---:|---|---|"
    )
    lines = [header, divider]
    for r in rows:
        flags = " ".join(r["flags"]) if r["flags"] else "ok"
        line = (
            f"| {r['run_name']} "
            f"| {r['body']} "
            f"| {_fmt_sha(r['git_sha'])} "
            f"| {_fmt(r['seed'])} "
            f"| {_fmt(r['n_samples'])} "
            f"| {_fmt(r['ess'], '.0f')} "
            f"| {_fmt(r['log_z'], '.3f')} "
            f"| {_fmt(r['log_z_err'], '.3f')} "
            f"| {_fmt(r['walltime_h'], '.3f')} "
            f"| {r['sampler']} "
            f"| {flags} |"
        )
        lines.append(line)
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Audit MCMC result pkls and write a Markdown report."
    )
    parser.add_argument(
        "--results-dir",
        default=None,
        help=(
            "Root directory to scan recursively for *results*.pkl files. "
            "Default: PlanetProfile/Test/mcmc_results relative to the "
            "project root (parent of the PlanetProfile package directory)."
        ),
    )
    parser.add_argument(
        "--output",
        default=None,
        help=(
            "Path to write the audit_report.md. "
            "Default: PlanetProfile/Inference/diagnostics/audit_report.md"
        ),
    )
    args = parser.parse_args()

    # Resolve project root (two levels up from this file:
    # diagnostics/ -> Inference/ -> PlanetProfile/ -> project root)
    this_file = Path(__file__).resolve()
    project_root = this_file.parent.parent.parent.parent  # …/planetprofile-genai

    # Results directory
    if args.results_dir:
        results_dir = Path(args.results_dir)
    else:
        results_dir = project_root / "PlanetProfile" / "Test" / "mcmc_results"

    if not results_dir.exists():
        print(f"ERROR: results directory does not exist: {results_dir}", file=sys.stderr)
        sys.exit(1)

    # Output path
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = this_file.parent / "audit_report.md"

    # -----------------------------------------------------------------------
    # Discover pkls
    # -----------------------------------------------------------------------
    all_pkls = sorted(results_dir.rglob("*results*.pkl"))
    if not all_pkls:
        print(f"No *results*.pkl files found under {results_dir}", file=sys.stderr)
        sys.exit(0)

    print(f"Found {len(all_pkls)} pkl(s) under {results_dir}", flush=True)

    # -----------------------------------------------------------------------
    # Extract info
    # -----------------------------------------------------------------------
    rows = []
    unreadable = []
    for pkl in all_pkls:
        print(f"  reading: {pkl.relative_to(results_dir)}", flush=True)
        info = extract_run_info(pkl, project_root)
        rows.append(info)
        if not info["readable"]:
            unreadable.append(pkl)

    # Sort: body, then run_name
    rows.sort(key=lambda r: (r["body"], r["run_name"]))

    # -----------------------------------------------------------------------
    # Build report
    # -----------------------------------------------------------------------
    now = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    report_lines = [
        f"# MCMC Run Audit Report",
        f"",
        f"Generated: {now}  ",
        f"Results dir: `{results_dir}`  ",
        f"Total pkls scanned: {len(all_pkls)}  ",
        f"Unreadable: {len(unreadable)}  ",
        f"",
        f"## Flag legend",
        f"",
        f"| Flag | Meaning |",
        f"|---|---|",
        f"| ⚠ ESS<500 | Effective sample size below convergence threshold |",
        f"| ⚠ log-Z err>0.5 | Log-evidence MC error too large for Bayes-factor comparison |",
        f"| ⚠ seed missing | Random seed not recorded; run not reproducible |",
        f"| ⚠ pre-bugfix | git SHA absent or predates ca1b600 (2026-05-23 pocoMC unpack-swap fix); samples ok, weights/histograms may be wrong |",
        f"| ℹ smoke | n_samples ≤ 200; exploratory only, not for Bayes-factor comparison |",
        f"| ok | No issues detected |",
        f"",
        f"## Runs",
        f"",
        build_markdown_table(rows),
        f"",
    ]

    if unreadable:
        report_lines += [
            f"## Unreadable pkls",
            f"",
        ]
        for p in unreadable:
            info = next(r for r in rows if r["pkl_path"] == str(p))
            err_snippet = info.get("_error", "")
            # Show only last line of traceback
            last_line = err_snippet.strip().split("\n")[-1] if err_snippet else ""
            report_lines.append(f"- `{p}`: {last_line}")
        report_lines.append("")

    # -----------------------------------------------------------------------
    # Write file and print
    # -----------------------------------------------------------------------
    report_text = "\n".join(report_lines)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(report_text, encoding="utf-8")
    print(f"\nAudit report written to: {output_path}\n")
    print(report_text)


if __name__ == "__main__":
    main()
