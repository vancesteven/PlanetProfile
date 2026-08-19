#!/usr/bin/env python
"""Write a seed-patched copy of an inference config into scratch.

`run_inference_cli` has no `--seed` flag -- the seed comes from the
config's `random_state` (and the sampler knobs from `sampler_settings`).
The multi-seed protocols this project runs (B3 reference wander, r7's
"3 seeds each" two-run branch comparison) therefore need one config per
(config, seed) pair. Rather than committing 6 near-identical configs per
campaign, this emits them into scratch at submit time, so the committed
config stays the single source of truth and the only difference on disk
is the seed it was stamped with.

Records what it changed in `metadata.cluster_seed_patch` so a result
pickle can always be traced back to the committed parent config.

Usage:
  python plans/scripts/cluster_seeded_config.py \
      --config PlanetProfile/Inference/configs/foo.json \
      --seed 1234 --out /scratch/foo_seed1234.json \
      [--set sampler_settings.n_effective=8000] [--set ...]
"""
import argparse
import json
from pathlib import Path


def _coerce(text):
    """JSON-parse a --set value, falling back to the raw string."""
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return text


def _assign(cfg, dotted, value):
    """Set a dotted path, refusing to CREATE a new leaf.

    A typo like `sampler_settings.n_efective=8000` must not silently add a
    dead key that the sampler ignores -- that is exactly the class of
    "sampled but never used" defect this campaign has already had to
    adjudicate (dead-parameter criterion). Parent containers must exist and
    the leaf must already be present.
    """
    parts = dotted.split('.')
    node = cfg
    for p in parts[:-1]:
        if not isinstance(node, dict) or p not in node:
            raise SystemExit(f'--set {dotted}: no such section {p!r} in the '
                             f'config; refusing to create it')
        node = node[p]
    leaf = parts[-1]
    if not isinstance(node, dict) or leaf not in node:
        raise SystemExit(f'--set {dotted}: key {leaf!r} does not exist in '
                         f'the config; refusing to create it (a typo here '
                         f'would be a silently-ignored setting)')
    old = node[leaf]
    node[leaf] = value
    return old


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', required=True)
    ap.add_argument('--seed', type=int, required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--set', action='append', default=[],
                    metavar='DOTTED.KEY=VALUE',
                    help='Override an EXISTING config key (repeatable).')
    args = ap.parse_args()

    src = Path(args.config)
    with open(src) as f:
        cfg = json.load(f)

    changes = {}
    if 'random_state' not in cfg:
        raise SystemExit(f'{src} has no top-level random_state to patch')
    changes['random_state'] = {'from': cfg['random_state'],
                              'to': args.seed}
    cfg['random_state'] = args.seed

    for item in args.set:
        if '=' not in item:
            raise SystemExit(f'--set {item!r} is not DOTTED.KEY=VALUE')
        dotted, raw = item.split('=', 1)
        value = _coerce(raw)
        old = _assign(cfg, dotted, value)
        changes[dotted] = {'from': old, 'to': value}

    cfg.setdefault('metadata', {})['cluster_seed_patch'] = {
        'parent_config': str(src),
        'generated_by': 'plans/scripts/cluster_seeded_config.py',
        'changes': changes,
        'note': ('Scratch-only derivative for one array task. The committed '
                 'parent config is the source of truth; only the listed '
                 'keys differ.'),
    }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, 'w') as f:
        json.dump(cfg, f, indent=1)
    print(f'{src.name} + seed {args.seed} -> {out}')
    for k, v in changes.items():
        print(f'  {k}: {v["from"]} -> {v["to"]}')


if __name__ == '__main__':
    main()
