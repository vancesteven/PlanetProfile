# MoonMag (maintained in-tree)

Magnetic induction forward modeling for icy moons: spherically symmetric and asymmetric
induced-field responses, excitation-moment data, and spacecraft-trajectory evaluation.

## Provenance and maintenance

MoonMag was originally developed as a standalone package by M. J. Styczinski. It is
**fully vendored here and maintained as first-party PlanetProfile code** — the copy in this
directory is the authoritative one for PlanetProfile; it is not synced from, and does not
depend on, any external repository or package. Fixes and modernizations are applied directly
in-tree (e.g., NumPy 2.0 dtype migration; SciPy ≥1.17 `sph_harm` → `sph_harm_y` migration in
`asymmetry_funcs.py`, verified equivalent to machine precision).

## Citation

When publishing results that use the induction calculations in this directory, cite the
method paper:

> Styczinski, M. J., Vance, S. D., Harnett, E. M., & Cochrane, C. J. (2022). A perturbation
> method for evaluating the magnetic field induced from an arbitrary, asymmetric ocean world
> analytically. *Icarus*, 376(1), 114840. https://doi.org/10.1016/j.icarus.2021.114840

and the PlanetProfile framework paper describing its integration:

> Styczinski, M. J., Vance, S. D., & Melwani Daswani, M. (2023). PlanetProfile: Self-consistent
> interior structure modeling for ocean worlds and rocky dwarf planets in Python.
> *Earth and Space Science* (see `papers/styczinski2023planetprofile.pdf`).

## Layout

- `symmetry_funcs.py` / `asymmetry_funcs.py` — induced-moment calculations (symmetric /
  asymmetric boundaries)
- `eval_induced_field.py`, `field_xyz.py`, `trajec_analysis.py` — field evaluation
- `excitation/` — excitation-moment data tables per body (packaged with PlanetProfile)
- `induced/`, `interior/` — outputs and interior-model inputs
- `plotting_funcs.py`, `run_all.py`, `config.py` — plotting and standalone-run scaffolding

Smoke coverage: `tests/moonmag_smoke_test.py` guards the import chain and core spherical-
harmonic evaluation against dependency API drift (e.g., the SciPy 1.17 `sph_harm` removal).
