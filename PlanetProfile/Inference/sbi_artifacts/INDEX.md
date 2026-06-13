# SBI Artifact Index

Last audited: 2026-06-12

## Status

No pre-trained SBI artifacts found in repository.

The `sbi_artifacts/` directory was created (commit history shows it was added 2026-06-12 at 13:39)
but contains no files. No `.pt`, `.pth`, `.pkl`, or `.npz` SBI training datasets exist anywhere
in the repository tree under this directory or elsewhere under a SBI-specific naming convention.

## Expected artifacts (for AmortizedInference page)

The `PlanetProfileApp/pages/AmortizedInference.py` page is scaffolded but currently **disabled**
(the "Generate Instant Posterior" button is `disabled=True` and a status banner reads
"This tab is currently under development").

The page UI offers three model slots that imply three eventual artifact files:

| UI Label | Implied artifact |
|---|---|
| Titan (Andrade, No-Ocean) | Normalizing-flow posterior estimator trained on Titan no-ocean simulations with Andrade rheology |
| Titan (Maxwell, Ocean) | Normalizing-flow posterior estimator trained on Titan ocean-bearing simulations with Maxwell rheology |
| Europa (Andrade, Thin Shell) | Normalizing-flow posterior estimator trained on Europa thin-shell simulations with Andrade rheology |

The page accepts two observables as inputs: `Re(k2)` and `Im(k2)`.

No explicit file-load call exists in `AmortizedInference.py` yet — artifact paths and a
`torch.load` / `sbi` posterior call have not been wired in.

## SBI training pipeline (not yet run)

- `PlanetProfile/Inference/mcmc_runner.py` exposes `MCMCRunner.generate_sbi_dataset(n_samples, output_path)`
  which produces `(theta, x)` NumPy arrays and saves them as `.npz` when `output_path` is supplied.
- `inference_core.py` references `from .sbi_runner import SBIRunner` (mode `'sbi'`) but
  `sbi_runner.py` does not exist in the `Inference/` directory — the SBI training runner is
  not yet implemented.
- The `sbi` Python package (and `torch`) must be installed separately (`pip install sbi torch`)
  before any training can occur.

## Related files

- `plans/titan-mcmc-sbi-plan.md` — roadmap; Phase 2 (SBI pre-training) is marked complete for
  dataset-export scaffolding but no training run has been executed.
- `PlanetProfile/Inference/mcmc_runner.py` — `generate_sbi_dataset` method (line 722)
- `PlanetProfileApp/pages/AmortizedInference.py` — GUI scaffold (no artifact load yet)

## Notes

Training is deferred. The current priority is MCMC validation on Titan and Europa (see
`plans/titan-mcmc-sbi-plan.md` and the MCMC task list in memory). Once MCMC posteriors are
stable, the `generate_sbi_dataset` pipeline can be used to produce training data, and
`sbi_runner.py` can be implemented to train and serialize the normalizing-flow estimators
into this directory.

Suggested artifact naming convention (not yet enforced):
```
sbi_artifacts/
  titan_andrade_noocean_posterior.pt
  titan_maxwell_ocean_posterior.pt
  europa_andrade_thinshell_posterior.pt
```
