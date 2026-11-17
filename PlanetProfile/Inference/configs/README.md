# Inference BodyConfigs

JSON configuration files consumed by `MCMCRunner(InferenceConfig.from_json(path))`.

Each file fully specifies an MCMC inference problem: which body, which parameters
to sample, the priors, any parameters that are derived (not sampled) from
mass-balance, the observables and their uncertainties, sampler settings, and where
to find the precomputed structure cache.

## Config status (`metadata.gui`)

Every config carries `metadata.gui = {"status": ..., "note": "..."}` (plus
`"runnable": false` where it cannot run). `metadata` is excluded from
`InferenceConfig.generate_hash`, so the block does not affect artifact hash
checks. The Streamlit Inference page lists only `production` configs by
default; a "Show non-production configs" checkbox reveals the rest, and a
config marked `runnable: false` is listed but its Run button is disabled with
the note as the reason.

| status | meaning |
|---|---|
| `production` | ratified or deployed result exists (STATUS.md / `sbi_artifacts/INDEX.md`) |
| `sensitivity` | ablation, slice or variant of a production config; delivered, not ratified |
| `smoke` | reduced-size config for fast sanity runs |
| `candidate` | not ratified, or the evidence is unclear; do not quote results as production |
| `retired` | superseded, vetoed or deprecated |

Statuses are judged from each config's own metadata, `plans/STATUS.md` and
`sbi_artifacts/INDEX.md`; when evidence was unclear the config is `candidate`.

| File | Body | Sampled dim | Tb | Status |
|---|---|---|---|---|
| `test50_titan_noocean_andrade_8D` | Titan (no ocean) | 8 | sampled | production |
| `test50_titan_noocean_andrade_5D` | Titan (no ocean) | 5 | fixed 250.965 K | sensitivity |
| `test52_titan_noocean_andrade_10D` | Titan (no ocean, core) | 10 | sampled | candidate |
| `test52_titan_noocean_andrade_10D_smoke` | Titan | 10 | sampled | smoke |
| `titan_freegrav_noocean` | Titan (no ocean, free gravity) | 12 | sampled | candidate |
| `test54_titan_nh3_freegrav` | Titan (NH3 joint) | 13 | sampled | production (k2 sector not verified) |
| `test54_titan_mgso4_freegrav` | Titan (MgSO4 joint) | 13 | sampled | production |
| `test54_titan_nacl_freegrav` | Titan (NaCl joint) | 13 | sampled | production |
| `europa_galileo_v1p1_8D` | Europa | 8 | sampled | production |
| `europa_clipper_v4_geodesy_11D` | Europa | 11 | sampled | production |
| `europa_clipper_v5_geodesy_11D`, `_v5_noinduction_7obs`, `_v5_nok2_17obs` | Europa | 11 | via `D_iceIh_km` | sensitivity |
| `europa_clipper_v6_freegrav_11D`, `_noinduction_6obs`, `_nok2_16obs` | Europa | 11 | via `D_iceIh_km` | sensitivity |
| `europa_clipper_v7_openae_11D` | Europa | 11 | via `D_iceIh_km` | sensitivity |
| `europa_seawater_andrade_7D`, `_clipper_v2`, `_clipper_v3_8D` | Europa | 7-8 | sampled | retired |
| `europa_seawater_smoke` | Europa | 7 | sampled | smoke |
| `ganymede_pureh2o_andrade_8D` | Ganymede | 8 | sampled | candidate |
| `callisto_nacl_andrade_8D` | Callisto | 8 | sampled | candidate |
| `callisto_nacl_andrade_8D_cmr2_minus5pct`, `_minus10pct` | Callisto | 8 | sampled | sensitivity |
| `callisto_mgso4_andrade_8D` | Callisto | 8 | sampled | retired, **cannot run** (defective cache) |
| `enceladus_cassini_isostasy_ocean6D_rhorock` | Enceladus | 6 | via `zb_km` x salinity | production |
| `enceladus_cassini_isostasy_7D` | Enceladus | 5 | via `zb_km` | candidate, **cannot run** |
| `enceladus_cassini_smoke_6D` | Enceladus | 6 | sampled | smoke |

All of the structure caches named by these configs exist in the repository
under `PlanetProfile/Test/mcmc_results/`. The two non-runnable configs are:

- `callisto_mgso4_andrade_8D`: its cache has malformed nodes (about one third
  of prior draws raise a `_blend_b_layered ... P_MPa has length ...` error)
  and the config is deprecated in favour of the NaCl one.
- `enceladus_cassini_isostasy_7D`: a frozen-config candidate that has never
  been run and cannot construct a runner from the CLI; the runnable
  Enceladus config is `enceladus_cassini_isostasy_ocean6D_rhorock`.

### Prior-type semantics (as the runner interprets them)

| `prior_type` | fields | meaning of `bounds` |
|---|---|---|
| `uniform` | `bounds` | linear limits in the parameter's units |
| `log-uniform` | `bounds` | log10 **exponents**: log-uniform on `[10^lo, 10^hi]` |
| `uniform-in-linear` | `bounds` | log10 exponents of the sampled variable, prior uniform in the linear value on `[10^lo, 10^hi]` |
| `normal` | `mean`, `std` | none |
| `truncated_gaussian` | `mean`, `std`, `bounds` | linear limits of the truncation |

Parameters sampled as `log10_*` (e.g. `log10_wOcean_ppt`) are themselves
exponents, with `uniform` bounds on the exponent.

### How the GUI uses a loaded config

Loading a config keeps the whole file as the base. The form edits only the
priors, observable values/sigmas, sampler fields and cache path; the Run
button rebuilds the config with `InferenceConfig.from_dict(<file> + edits)`, so
`bodyname`, `fixed_params`, `derived_params`, `param_groups`, `metadata`,
`induction_bounds`, `gravity_forward_model`, `arrhenius_params` and
`ocean_overrides` pass through unchanged. Untouched, the run config equals
`InferenceConfig.from_json(file)`. Logic lives in
`PlanetProfile/Inference/gui_config.py`.

### Structure caches

Build or extend a Tb (or Tb x salinity) cache from a config with

```bash
python -m PlanetProfile.Inference.build_phase_c1_cache \
    --config PlanetProfile/Inference/configs/<config>.json --n-grid <N> [--n-wgrid <M>]
```

(`--n-wgrid` adds the log-spaced salinity axis for `log10_wOcean_ppt`
configs.) The GUI checks that the selected cache exists and that each sampled
grid-axis prior (`Tb_K`, `log10_wOcean_ppt`, `zb_km`) and each fixed
`Tb_K` lies inside the cached range, and warns otherwise.

## Schema (v2 — adds `derived_params` and core sampling)

```jsonc
{
  "mode": "mcmc",
  "bodyname": "Body name (matches PP body directories)",
  "param_space": {
    "<param_id>": {"prior_type": "uniform" | "log-uniform" | "normal", "bounds": [...]}
    // Rheology (any subset):
    //   alpha, log10_zeta, log10_eta_Ih, log10_eta_HP, log10_eta_sil
    //   per-phase overrides: log10_eta_III, log10_eta_V, log10_eta_VI,
    //                        log10_zeta_Ih, log10_zeta_HP, log10_zeta_sil
    // Structure (v2):
    //   Tb_K              — basal temperature (K); requires Tb-grid cache
    //   R_core_km         — core radius (km); 0.0 admits the no-core limit
    //   rho_core_kgm3     — core density (kg/m^3)
  },
  "param_groups": {
    // Optional: collapse several parameters onto one sampled scalar.
    // Example: lock the three HP ice viscosities together.
    "log10_eta_HP": ["log10_eta_III", "log10_eta_V", "log10_eta_VI"]
  },
  "fixed_params": {
    // Optional: parameters injected into every forward-model call without
    // entering the prior. Use to pin Tb_K or salinity for slice studies.
    "Tb_K": 250.965
  },
  "derived_params": {
    // Optional (v2): parameters whose value is computed from other sampled
    // parameters via a physical constraint. Currently supported:
    //   rho_sil_kgm3: mass_conservation
    "rho_sil_kgm3": {
      "derivation": "mass_conservation",
      "bounds": [2200.0, 3500.0],
      "reject_if_outside_bounds": true
    }
  },
  "observables": {
    // Each entry is [value, sigma]. Supported keys (handled by mcmc_runner):
    //   Re_k2, Im_k2, abs_Im_k2, k2 (= sqrt(Re_k2^2 + Im_k2^2))
    //   CMR2, Mtot_kg
    "CMR2": [0.3115, 0.0028]
  },
  "sampler_settings": {
    "n_effective": 500,
    "n_reeval": 500,
    "rheology": "andrade" | "maxwell",
    "arrhenius_params": {
      "activation_energy_kJ_mol": {"Ih": 60.0},
      "R_J_mol_K": 8.314462
    },
    // Optional: phase-stability gate. Currently only "no_ocean_Ih" is wired.
    "phase_stability": {"enforce": "no_ocean_Ih", "margin_K": 0.1}
  },
  "structure_cache_path": "PlanetProfile/Test/mcmc_results/<body>_<comp>_structure_grid.pkl",
  "random_state": 42,
  "metadata": {
    // Free-form. Use schema_version, description, primary_constraint, etc.
    // Phase C1 v2 configs additionally include Tb_bounds_rationale,
    // core_bounds_rationale, stage2_dependencies, stage3_dependencies.
  }
}
```

### `derived_params` semantics (v2)

`derived_params.rho_sil_kgm3.derivation = "mass_conservation"` instructs the
runner to solve, at each MCMC sample:

```
rho_sil = (M_total - M_ice - M_ocean - M_core)
          / V_silicate

where  M_core      = (4/3) π R_core_km^3 ρ_core_kgm3
       V_silicate  = (4/3) π (R_oceanbot_m^3 - R_core_m^3)
       M_ice, M_ocean, R_oceanbot_m, M_total  come from the structure cache
```

If `rho_sil` falls outside `bounds = [2200.0, 3500.0]` and
`reject_if_outside_bounds = true`, the sample is assigned `log L = -∞` and
rejected. This is a physically meaningful prior cut: silicate densities below
2200 or above 3500 kg/m³ are inconsistent with plausible silicate
mineralogies for icy moons.

CMR² is then recomputed at sample time from the full layered profile:
sampled `(R_core, ρ_core)` + derived `ρ_sil` + cached ice/ocean profile.

## Structure cache formats

The runner accepts two cache shapes, picked automatically based on whether
`Tb_K` is in the sampled `param_space`, in `fixed_params`, or absent:

1. **Single structure** (used when `Tb_K` is fixed or not present): a pickled
   dict with keys like `r_m`, `rho_kgm3`, `T_K`, `P_MPa`, `phase`, `CMR2`,
   `Mtot_kg`, etc.
2. **Tb grid** (used when `Tb_K` is sampled): a pickled dict
   `{'Tb_K_grid': [...], 'structures': [structure_dict, ...]}` with one
   structure per grid point. Linear interpolation between bracketing grid
   points is bit-for-bit identical to Test 50's original `_interp_structure`.

**v2 cache schema additions** (required by `derived_params.mass_conservation`):
each `structure_dict` must additionally expose the scalars `M_ice_kg`,
`M_ocean_kg`, `R_oceanbot_m`, and `M_total_kg`. The Test 50 caches predate
this and won't carry these fields; that's fine — those configs don't use
`derived_params`.

## Adding a new BodyConfig

1. Pick a name: `<body>_<composition>_<rheology>_<dim>D.json`.
2. Copy a sibling config (Phase C1 v2 if your body has an ocean and a core;
   Test50 if it doesn't). Edit `bodyname`, prior bounds, `observables`,
   `structure_cache_path`, and `metadata` (including `metadata.gui`).
3. Decide which body radius caps `R_core_km` (use `0.5 × R_body` as a
   conservative upper bound unless you have a tighter physical constraint).
4. Build the structure cache with `build_phase_c1_cache` (see above).
5. Smoke-test the parse:
   ```bash
   /Users/svance/mamba/envs/PPcl/bin/python -c "from PlanetProfile.Inference.inference_core import InferenceConfig; \
     c = InferenceConfig.from_json('PlanetProfile/Inference/configs/<your_config>.json'); \
     print(c.bodyname, list(c.param_space.keys()))"
   ```
6. Run the MCMC:
   ```bash
   python -m PlanetProfile.Inference.run_inference_cli \
       --config PlanetProfile/Inference/configs/<your_config>.json
   ```

## Roadmap

The v2 toolkit (derived parameters, core sampling, mass conservation, 1D/2D
caches) is implemented. Open work and campaign status are tracked in
`plans/STATUS.md`; deployed SBI artifacts in `sbi_artifacts/INDEX.md`.
