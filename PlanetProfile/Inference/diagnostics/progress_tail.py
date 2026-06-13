#!/usr/bin/env python3
"""Tail an MCMC progress JSONL file and print human-readable progress lines.

Usage
-----
    python progress_tail.py <path/to/progress.jsonl>

Waits up to 60 s for the file to appear, then follows new lines as they are
appended (poll interval: 2 s).  Press Ctrl-C to exit cleanly.

Output format (one line per JSONL record):
    [HH:MM:SS] iter=N  log_Z=-45.2 ± 0.30  ESS=120  n_acc=85/500  elapsed=142s
"""

import argparse
import json
import sys
import time
from pathlib import Path


def _fmt(value, fmt='.2f', fallback='N/A'):
    """Format a value, returning fallback string when value is None."""
    if value is None:
        return fallback
    try:
        return format(value, fmt)
    except (TypeError, ValueError):
        return str(value)


def _format_record(record: dict) -> str:
    """Return a human-readable single-line summary of a JSONL progress record."""
    ts = record.get('timestamp', '')
    # Extract HH:MM:SS from ISO-8601 timestamp if present
    if 'T' in ts:
        hms = ts.split('T')[1].rstrip('Z')[:8]
    else:
        hms = time.strftime('%H:%M:%S')

    iteration = record.get('iteration')
    iter_str = str(iteration) if iteration is not None else 'N/A'

    log_z = record.get('log_Z')
    log_z_err = record.get('log_Z_err')
    if log_z is not None and log_z_err is not None:
        logz_str = f'{log_z:.2f} ± {log_z_err:.2f}'
    elif log_z is not None:
        logz_str = f'{log_z:.2f} ± N/A'
    else:
        logz_str = 'N/A'

    ess = record.get('ESS')
    ess_str = str(ess) if ess is not None else 'N/A'

    n_accepted = record.get('n_accepted')
    n_total = record.get('n_total')
    if n_accepted is not None and n_total is not None:
        nacc_str = f'{n_accepted}/{n_total}'
    elif n_total is not None:
        nacc_str = f'N/A/{n_total}'
    else:
        nacc_str = 'N/A'

    elapsed = record.get('elapsed_s')
    elapsed_str = f'{int(elapsed)}s' if elapsed is not None else 'N/A'

    return (
        f'[{hms}] iter={iter_str:>4}  '
        f'log_Z={logz_str:>14}  '
        f'ESS={ess_str:>6}  '
        f'n_acc={nacc_str:>10}  '
        f'elapsed={elapsed_str}'
    )


def tail_jsonl(path: str, poll_interval: float = 2.0,
               wait_timeout: float = 60.0) -> None:
    """Follow a JSONL file, printing each new line as it appears.

    Parameters
    ----------
    path:
        Path to the JSONL progress file.
    poll_interval:
        Seconds between polls when no new data is available.
    wait_timeout:
        Maximum seconds to wait for the file to appear before giving up.
    """
    jsonl_path = Path(path)

    # Wait for the file to appear
    waited = 0.0
    while not jsonl_path.exists():
        if waited >= wait_timeout:
            print(
                f'[progress_tail] File not found after {wait_timeout:.0f}s: {path}',
                file=sys.stderr,
            )
            sys.exit(1)
        print(
            f'[progress_tail] Waiting for {path} ... ({waited:.0f}s)',
            file=sys.stderr,
        )
        time.sleep(poll_interval)
        waited += poll_interval

    print(f'[progress_tail] Following {path}', file=sys.stderr)

    with open(jsonl_path, 'r') as fh:
        # Replay any lines already present
        for raw in fh:
            raw = raw.strip()
            if not raw:
                continue
            try:
                record = json.loads(raw)
                print(_format_record(record))
            except json.JSONDecodeError:
                print(f'[progress_tail] Malformed line: {raw}', file=sys.stderr)

        # Follow new lines
        while True:
            raw = fh.readline()
            if raw:
                raw = raw.strip()
                if raw:
                    try:
                        record = json.loads(raw)
                        print(_format_record(record), flush=True)
                    except json.JSONDecodeError:
                        print(
                            f'[progress_tail] Malformed line: {raw}',
                            file=sys.stderr,
                        )
            else:
                time.sleep(poll_interval)


def main() -> None:
    parser = argparse.ArgumentParser(
        description='Tail an MCMC progress JSONL file and print human-readable progress lines.'
    )
    parser.add_argument(
        'jsonl_path',
        metavar='PATH',
        help='Path to the JSONL progress file to tail.',
    )
    parser.add_argument(
        '--interval',
        type=float,
        default=2.0,
        metavar='SECONDS',
        help='Poll interval in seconds (default: 2.0).',
    )
    parser.add_argument(
        '--wait',
        type=float,
        default=60.0,
        metavar='SECONDS',
        help='Maximum seconds to wait for the file to appear (default: 60).',
    )
    args = parser.parse_args()

    try:
        tail_jsonl(args.jsonl_path, poll_interval=args.interval,
                   wait_timeout=args.wait)
    except KeyboardInterrupt:
        print('\n[progress_tail] Interrupted.', file=sys.stderr)
        sys.exit(0)


if __name__ == '__main__':
    main()
