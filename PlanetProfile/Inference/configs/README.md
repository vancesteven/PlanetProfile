# Inference BodyConfigs

JSON configuration files consumed by `MCMCRunner(InferenceConfig.from_json(path))`.

Each file fully specifies an MCMC inference problem: which body, which parameters
to sample, the priors, any parameters that are derived (not sampled) from
mass-balance, the observables and their uncertainties, sampler settings, and where
to find the precomputed structure cache.

## Available configs

| File | Body | Dim | `Tb_K` | Status |
|---|---|---|---|---|
| `test50_titan_noocean_andrade_8D.json` | Titan (no-ocean) | 8 | sampled `[249.0, 250.965]` | Runnable (cache exists) |
| `test50_titan_noocean_andrade_5D.json` | Titan (no-ocean) | 5 | fixed `250.965` | Runnable (cache exists) |
| `europa_seawater_andrade_7D.json` | Europa (seawater) | 7 | sampled `[260.5, 270.0]` | **Schema v2 — needs Stage 2 toolkit + Stage 3 cache** |
| `ganymede_pureh2o_andrade_8D.json` | Ganymede (pure H₂O) | 8 | sampled `[252.2, 261.0]` | **Schema v2 — needs Stage 2 toolkit + Stage 3 cache** |
| `callisto_mgso4_andrade_8D.json` | Callisto (MgSO₄ 100 ppt) | 8 | sampled `[252.2, 263.0]` | **Schema v2 — needs Stage 2 toolkit + Stage 3 cache** |

The Titan configs were the Phase B reference and use the v1 schema (rheology +
optional `Tb_K` only). The three ocean-world configs are Phase C1 v2: they add
sampled core parameters and a derived (mass-conserving) silicate density.
They parse cleanly through `InferenceConfig.from_json()` but are unrunnable
end-to-end until the toolkit's mass-conservation solver lands (Stage 2) and the
per-body structure caches are built (Stage 3).

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
   `structure_cache_path`, and `metadata`.
3. Decide which body radius caps `R_core_km` (use `0.5 × R_body` as a
   conservative upper bound unless you have a tighter physical constraint).
4. Build the structure cache (Phase C2 / C3 work — not yet automated).
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

### Stage 2 — toolkit support for v2 schema

- Parse `derived_params` block in `InferenceConfig`.
- Add parameter hooks for `R_core_km`, `rho_core_kgm3` in
  `apply_*_params` / runner.
- Implement `mass_conservation` derivation for `rho_sil_kgm3` inside
  `MCMCRunner._make_flexible_log_likelihood`, with hard-bound rejection.
- Recompute CMR² at sample time from the layered profile assembled from
  cached ice/ocean + sampled core + derived silicate.
- Extend the structure-cache schema to expose `M_ice_kg`, `M_ocean_kg`,
  `R_oceanbot_m`, `M_total_kg` per Tb point.
- Test-driven: write unit tests for the mass-conservation solver against
  hand-computed cases (e.g. Earth-like proportions) before wiring into the
  runner.

### Stage 3 — per-body cache builds

- Probe Tb range for each body: run PP at the bracket Tb values and confirm
  the resulting Ih shell thickness matches the design target
  (Europa 3-100 km; Ganymede/Callisto ~50 km at the thin end, near no-ocean
  at the thick end).
- Build Tb-grid caches with the v2 schema additions.
- Run a small smoke MCMC (~50 effective samples) per config to verify the
  end-to-end path.

### Beyond Phase C1

- **Per-composition runs** for Ganymede and Callisto: 0.1 molal increments
  of NaCl in the ocean using the SeaFreeze NaCl representation. Builds a
  composition axis on top of the Tb axis. Deferred — not part of Phase C1.
- **Induction observables** (Galileo: Khurana et al. 1997 / Kivelson et al.
  2002 / Zimmer et al. 2000) — requires extending the runner's observable
  switch beyond `Re_k2 / Im_k2 / k2 / abs_Im_k2 / CMR2 / Mtot_kg`.
- **Gravity J₂ / C₂₂ constraints** (Anderson et al. 1996 for Ganymede,
  Anderson et al. 2001 for Callisto, Schubert et al. for Europa) — same
  toolkit extension.
- **No-ocean Callisto variant** as a separate config, analogous to Test 50
  for Titan.
- **Posterior viewer** in PlanetProfileApp that loads multiple `.pkl`
  results and overlays them.
