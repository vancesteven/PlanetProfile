#!/usr/bin/env python3
"""
Train an SBI posterior artifact from a pre-generated (theta, x) dataset.

Reusable entry point for the SBI training campaigns (Test50 Titan and
subsequent bodies). Decouples training from dataset generation so a single
large ``.npz`` (produced by ``SBIRunner.generate_training_set`` and saved with
``np.savez``) can be trained repeatedly under different estimators/seeds
without re-running the forward model.

The dataset ``.npz`` must contain ``theta`` (n, D) and ``x`` (n, K) arrays in
the config's ``param_names`` / ``obs_names`` order. If it also contains a
``stats`` field (JSON string of the generate_training_set rejection stats),
it is carried into the artifact's ``rejection_stats`` provenance so the saved
``.pt`` records how its training data was drawn.

Usage:
    python -m PlanetProfile.Inference.train_sbi_artifact \\
        --dataset /tmp/train_1m.npz \\
        --config PlanetProfile/Inference/configs/test50_titan_noocean_andrade_8D.json \\
        --density-estimator nsf --seed 42 \\
        --output PlanetProfile/Inference/sbi_artifacts/titan_andrade_noocean_posterior_1m.pt

Author: PlanetProfile Team
Date: 2026-07-10
"""
import argparse
import json
import logging
import sys
import time

import numpy as np

from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.sbi_runner import SBIRunner

log = logging.getLogger(__name__)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--dataset', required=True,
        help='Pre-generated (theta, x) .npz from generate_training_set + np.savez.',
    )
    parser.add_argument(
        '--config', required=True,
        help='InferenceConfig JSON (forced to mode="sbi" for training).',
    )
    parser.add_argument(
        '--output', required=True,
        help='Output artifact .pt path (parent dirs created).',
    )
    parser.add_argument(
        '--seed', type=int, default=42,
        help='torch.manual_seed for reproducible training (default: 42).',
    )
    parser.add_argument(
        '--density-estimator', default='nsf', choices=['maf', 'nsf'],
        help='Normalizing-flow family (default: nsf).',
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')

    # Load dataset (theta, x, optional stats provenance).
    with np.load(args.dataset, allow_pickle=True) as npz:
        theta = np.asarray(npz['theta'], dtype=np.float64)
        x = np.asarray(npz['x'], dtype=np.float64)
        stats = None
        if 'stats' in npz.files:
            raw = npz['stats']
            try:
                stats = json.loads(str(raw.item()) if hasattr(raw, 'item') else str(raw))
            except (ValueError, TypeError):
                stats = None
    log.info(f"Loaded dataset {args.dataset}: theta={theta.shape}, x={x.shape}")

    # Build the SBI runner from the config.
    cfg = json.load(open(args.config))
    cfg['mode'] = 'sbi'
    runner = SBIRunner(InferenceConfig.from_dict(cfg))

    # Train.
    t0 = time.time()
    runner.train(
        theta, x, seed=args.seed,
        density_estimator=args.density_estimator,
    )
    # Carry the dataset's rejection stats into the artifact provenance.
    if stats is not None:
        runner._train_info['rejection_stats'] = stats

    runner.save_artifact(args.output)
    log.info(
        f"Trained {args.density_estimator} (seed {args.seed}) on {len(theta)} sims "
        f"in {time.time() - t0:.0f}s -> {args.output}"
    )


if __name__ == '__main__':
    sys.exit(main())
