# PlanetProfile MCMC Inference Toolkit

Bayesian inference of interior-structure parameters for ocean worlds and rocky
dwarf planets, given observables such as the moment-of-inertia factor (CMR²),
tidal Love numbers (Re k₂, Im k₂), or total mass (Mtot_kg).  Parameters that
can be constrained include rheological constants (Andrade creep exponent α,
fluidity ζ, phase-specific viscosities η), basal ice temperature Tb_K, core
radius, core density, and—via mass conservation—silicate mantle density.  The
toolkit is designed for HPC-scale throughput: a single MCMC run typically costs
minutes to tens of minutes on a workstation, not hours.

## Two-tier architecture

The toolkit separates an expensive one-time step from a cheap many-times step.

**Tier 1 — cache build.** PlanetProfile's full layer-propagation engine is
expensive: a single run at one Tb_K takes on the order of one to a few minutes
depending on EOS table coverage.  Running PP at every MCMC sample is
therefore impractical.  Instead, `cache_builder.build_tb_grid_cache` runs PP
at a small grid of Tb_K values (typically 9–15 points spanning the prior
bounds) and serialises each resulting per-cell structure to a pickle file.
Total wall time for a cache build is roughly 30–60 minutes.  The cache needs
to be built only once per body/composition combination; it is reused across
every subsequent MCMC run.

**Tier 2 — MCMC sampling.** `MCMCRunner` loads the pre-built cache and evaluates
the log-likelihood entirely within Python/NumPy by interpolating between cached
structures (the B-layered blend, described in section 4) and then computing
observables—CMR², k₂—from the resulting structure.  Each likelihood evaluation
costs a few milliseconds, giving a throughput of roughly 10 evaluations per
second.  At that rate, 10⁵ total evaluations (the typical pocoMC workload for
500 effective samples in an 8D problem) complete in minutes.

## Key entry points

**Programmatic:**

```python
from PlanetProfile.Inference.inference_core import InferenceConfig
from PlanetProfile.Inference.mcmc_runner import MCMCRunner

config = InferenceConfig.from_json("PlanetProfile/Inference/configs/callisto_nacl_andrade_8D.json")
runner = MCMCRunner(config)
result = runner.run()
result.save("callisto_run1.pkl")
```

**Headless / HPC:**

```bash
# Validate config without running
python -m PlanetProfile.Inference.run_inference_cli \
    --config PlanetProfile/Inference/configs/callisto_nacl_andrade_8D.json \
    --validate-only

# Full run, output to a named file
python -m PlanetProfile.Inference.run_inference_cli \
    --config PlanetProfile/Inference/configs/callisto_nacl_andrade_8D.json \
    --output results/callisto_run1.pkl
```

**Cache build:**

```bash
python -m PlanetProfile.Inference.build_phase_c1_cache \
    --config PlanetProfile/Inference/configs/callisto_nacl_andrade_8D.json \
    --n-grid 9 --force
```

---

## 1. pocoMC: preconditioned Monte Carlo in this framework

### What pocoMC is

pocoMC (Karamanis et al. 2022) is a Sequential Monte Carlo (SMC) sampler that
uses a normalising flow to *precondition* the posterior before each SMC
tempering step.  Unlike vanilla SMC—which draws proposals in the original
parameter space—pocoMC trains a normalising flow on the current population of
weighted particles, then uses the flow to map the parameter space into an
approximately Gaussian representation.  Subsequent MCMC steps happen in that
transformed space, where standard isotropic proposals mix efficiently.  After
each step the flow is retrained on the updated population.

This architecture makes pocoMC particularly effective when the posterior is
multimodal or has strong parameter correlations—both common in ice-shell
inference, where an ice-dominated regime and an ocean-dominated regime can
produce similar CMR² values through different combinations of layer thickness
and density.

### Sequential Monte Carlo tempering

pocoMC starts from the prior (inverse temperature β = 0, meaning all samples
are equally weighted by prior probability) and walks β toward 1 (the full
posterior) in adaptive steps.  At each β level it maintains a population of
weighted particles.  Particles are resampled to have uniform weight and then
dispersed via short MCMC chains in the *preconditioned* (flow-transformed)
space.  The β schedule is adapted so that each step's effective sample size
(ESS) does not drop below a user-controllable floor.

This tempering strategy avoids the "cold start" failure mode of plain MCMC:
a chain initialised from the prior will reach the high-probability region
regardless of where it starts, because the posterior is introduced gradually
rather than all at once.  It also sidesteps mode-switching problems that
afflict fixed-temperature samplers on multimodal posteriors.

### n_effective: target effective sample size

`n_effective` sets the ESS at which pocoMC stops.  Termination is not governed
by a fixed number of iterations or likelihood evaluations; the sampler runs
until the posterior particle population has at least `n_effective` independent
draws.  ESS counts independent samples: a population of 1000 strongly
correlated particles might have ESS = 50, meaning only 50 independent draws
worth of information.

Practical values used in this toolkit:
- `n_effective = 100` — smoke test; confirms the end-to-end path works; not
  suitable for scientific conclusions.
- `n_effective = 500` — production default, set in
  `callisto_nacl_andrade_8D.json` and the other Phase C1 configs.  Yields
  reliable marginal distributions and meaningful Bayes factors for 5–10D
  problems.
- `n_effective = 1000` or higher — warranted when the posterior has sharp
  secondary modes or very high parameter correlations, at roughly 2× wall
  time.

The `validate_config` function in `inference_core.py` enforces a minimum of
100; configs with lower values are rejected at parse time (`inference_core.py`
line 379 region).

### n_reeval and the hardcoded n_total

`n_reeval` in `sampler_settings` is the intended refinement particle count
(pocoMC's `n_total` argument), which controls the accuracy of the evidence
integral and tail coverage.  It should be ≥ `n_effective`.

However, `mcmc_runner.py:396` currently calls:

```python
sampler.run(n_total=4096, progress=True)
```

The `4096` is hardcoded and does not read `n_reeval` from the config.  This
means `n_reeval` in the JSON has no effect on the actual run until this line is
updated to `sampler.run(n_total=self.n_reeval, progress=True)`.  Treat
`n_reeval` as a placeholder that documents intent; the actual refinement
particle count is 4096 in all current runs.

### What the runner extracts from pocoMC

After `sampler.run()` completes, `mcmc_runner.py` calls:

```python
samples, log_likes, log_post, _ = sampler.posterior()
```

This returns the weighted posterior particle population.  `samples` is an array
of shape `(n_samples, n_params)`.  `log_likes` contains the log-likelihood at
each sample.  The runner then recomputes k₂ for a diagnostic subset of
posterior samples using the full rheology forward model (`forward_models.py`),
and packages everything into an `InferenceResult` dataclass defined in
`inference_core.py:171`.  Fields include:

| Field | Shape / type | Description |
|---|---|---|
| `samples` | `(n_samples, n_params)` | Posterior samples, column order matches `param_names` |
| `log_likelihoods` | `(n_samples,)` | Log-likelihood at each sample |
| `param_names` | `List[str]` | Parameter names in sampling order |
| `k2_results` | `(n_eval, 2)` or None | Re k₂ and Im k₂ for a diagnostic subset |
| `cmr2_results` | `(n_samples,)` or None | CMR² for all posterior samples |
| `D_iceIh_results` | `(n_samples,)` or None | Ice Ih shell thickness (km) |
| `convergence_metrics` | dict | ESS, acceptance rate, R-hat proxy |

The `InferenceResult.get_summary_stats()` method (`inference_core.py:263`)
returns median, mean, std, and 1σ/2σ credible intervals per parameter as a
dict of arrays.

### Acceptance rate and convergence proxy

pocoMC does not expose individual chains the way emcee or PyMC do, so
the classical Gelman–Rubin R-hat diagnostic is not directly applicable.  The
runner approximates acceptance rate from the progress bar's `acc` field
(`mcmc_runner.py:524`):

```python
acceptance_rate = sampler.pbar.info.get('acc')
```

This is a sampler health metric—values near 0 or 1 indicate the flow
transformation has failed—not a convergence diagnostic in the posterior sense.
The runner stores it in `convergence_metrics['acceptance_rate']` and logs it
after the run.  Treat it as a red flag for degenerate runs, not as a
sufficiency condition for trustworthy posteriors.  The primary convergence
indicator is whether `ESS >= n_effective` (guaranteed by pocoMC's termination
condition).

### When pocoMC is the right choice

pocoMC is well-suited to this problem because:

- Likelihood throughput is ~10 evaluations/second with the cached structures,
  which matches the SMC regime where total evaluations are O(10⁴)–O(10⁵)
  rather than O(10⁶)–O(10⁷).
- The parameter space is 5–10D with bounded uniform priors—a regime where
  normalising flows are well-conditioned and the tempering schedule is short.
- Multimodal posteriors arise naturally in ice-shell inference (e.g.,
  ice-dominated vs. ocean-dominated solutions); SMC tempering handles these
  without requiring mode-jumping proposals.

**Compared to emcee:** emcee is better understood and easier to debug for
unimodal posteriors, but struggles with multimodality and correlated parameters
without careful ensemble tuning.  pocoMC's preconditioned proposals handle
both at comparable cost.

**Compared to nested sampling (dynesty, nautilus):** nested samplers provide
accurate evidence estimates and naturally handle multimodality, but their
wall-time advantage is strongest when the prior volume greatly exceeds the
posterior volume.  For the bounded uniform priors used here, pocoMC typically
reaches `n_effective = 500` in fewer likelihood evaluations.

**Caveats:**
- The normalising flow is an extra hyperparameter not exposed in our configs;
  the pocoMC defaults are appropriate for our dimensionality.
- For problems above ~20D the flow becomes harder to train; all current
  configs are ≤ 10D and well below this limit.
- For sharply-peaked posteriors the adaptive β schedule may take many
  tempering steps; `n_effective` gates termination, so wall time is
  indeterminate but correctness is not compromised.

---

## 2. Configuration: InferenceConfig and the JSON schema

`InferenceConfig` is a Python dataclass (`inference_core.py:20`) that fully
specifies one inference problem.  It can be instantiated directly or loaded
from a JSON file via `InferenceConfig.from_json(path)`.

### Core fields

```jsonc
{
  "mode": "mcmc",                // "mcmc" or "sbi"
  "bodyname": "Callisto",        // Must match a PP body directory
  "param_space": { ... },        // Parameters to sample and their priors
  "param_groups": { ... },       // Optional: collapse params onto shared scalar
  "fixed_params": { ... },       // Optional: constants injected into every call
  "derived_params": { ... },     // Optional (v2): mass-conservation derivations
  "ocean_overrides": { ... },    // Optional (v2.1): Planet.Ocean.* overrides
  "observables": { ... },        // Observed values and uncertainties
  "sampler_settings": { ... },   // n_effective, n_reeval, rheology, arrhenius_params
  "structure_cache_path": "...", // Path to pre-built cache pickle
  "arrhenius_params": { ... },   // Optional (top-level): Arrhenius viscosity params, e.g. {"E_Ih_J_per_mol": 60000}.
                                 // Preferred over sampler_settings["arrhenius_params"]; stored on the
                                 // saved InferenceResult pickle so post-hoc reanalyze_k2_from_pickle
                                 // calls are self-describing without the caller re-supplying this dict.
  "planet_template_module": "PlanetProfile.Test.PPTest48",  // Optional: importable PPTest module used
                                 // by plot_structure_wedge_pp to re-run PlanetProfile at the posterior
                                 // median and produce a canonical PP wedge plot. Must be set here (or
                                 // passed explicitly) for plot_structure_wedge_pp to work.
  "random_state": 42,
  "metadata": { ... }
}
```

### param_space

Each entry maps a parameter name to a prior specification:

```jsonc
"param_space": {
  "alpha":        {"prior_type": "uniform", "bounds": [0.15, 0.45]},
  "log10_eta_Ih": {"prior_type": "uniform", "bounds": [10.0, 16.0]},
  "Tb_K":         {"prior_type": "uniform", "bounds": [252.2, 263.0]}
}
```

Supported prior types: `"uniform"`, `"log-uniform"`, `"normal"`.  The sampled
parameters become the columns of `InferenceResult.samples` in declaration order.

### param_groups

Collapses several parameters onto a single sampled scalar.  The group key
must appear in `param_space`; the member names must not:

```jsonc
"param_groups": {
  "log10_eta_HP": ["log10_eta_III", "log10_eta_V", "log10_eta_VI"]
}
```

At every likelihood evaluation, the runner expands `log10_eta_HP` into
`log10_eta_III = log10_eta_V = log10_eta_VI = <sampled value>` before calling
the forward model.  This reduces the dimensionality when HP-ice viscosities are
physically expected to be similar, as in `callisto_nacl_andrade_8D.json`.

### derived_params

Declares parameters computed from sampled values via a physical constraint rather
than sampled directly.  Currently only `rho_sil_kgm3` with
`derivation = "mass_conservation"` is wired in the runner:

```jsonc
"derived_params": {
  "rho_sil_kgm3": {
    "derivation": "mass_conservation",
    "bounds": [2200.0, 3500.0],
    "reject_if_outside_bounds": true
  }
}
```

At each sample, the runner solves:

```
rho_sil = (M_total - M_ice - M_ocean - M_core) / V_silicate
```

where M_core and V_silicate are computed from the sampled `R_core_km` and
`rho_core_kgm3`, and M_ice, M_ocean, R_oceanbot_m, M_total come from the
structure cache.  If `rho_sil` falls outside the stated bounds and
`reject_if_outside_bounds` is true, the sample is assigned log L = −∞.  This
implements a physically motivated prior cut: silicate densities outside
[2200, 3500] kg/m³ are inconsistent with plausible icy-moon mineralogies.

### Observables

```jsonc
"observables": {
  "CMR2": [0.3549, 0.0042]
}
```

Each entry is `[value, 1σ_uncertainty]`.  Supported keys: `Re_k2`, `Im_k2`,
`abs_Im_k2`, `k2` (= |k₂|), `CMR2`, `Mtot_kg`.  The log-likelihood is a sum
of Gaussian terms:

```
log L = -0.5 * Σ  [(obs_model - obs_measured) / sigma]²
```

implemented in `MCMCRunner._make_flexible_log_likelihood` (around
`mcmc_runner.py:340–357`).

### Validation

`validate_config(config)` (`inference_core.py:347`) returns `(is_valid, errors)`.
The `--validate-only` CLI flag runs this check and exits without sampling.
`InferenceConfig.__post_init__` (`inference_core.py:87`) enforces stricter
structural constraints at construction time (non-overlapping param_space /
param_groups / fixed_params; positive observable uncertainties).

For the complete field-by-field JSON schema see `configs/README.md`.

---

## 3. v2.1 cache schema and transition refinement

### What the cache pickle contains

A v2.1 structure cache is a Python dict with four top-level keys:

```python
{
    'Tb_K_grid':      np.ndarray,          # shape (N,), sorted ascending, units K
    'structures':     [dict, dict, ...],   # N structure dicts, one per grid point
    'transitions':    [dict, dict, ...],   # metadata for detected layer-set changes
    'schema_version': 'v2.1',
}
```

Each `structures[i]` is a per-cell dict built by `cache_builder.build_single_structure`
(`cache_builder.py:78`).  The per-cell arrays are parallel: index k of `r_m`,
`rho`, `K_Pa`, `mu_Pa`, `eta_Pa_base`, `T_K`, `P_MPa`, `bulk_visc`, and
`phases` all refer to the same radial cell.  Metadata fields on each structure
dict:

| Key | Type | Description |
|---|---|---|
| `r_m` | `np.ndarray` | Cell radii (m), ascending |
| `rho` | `np.ndarray` | Density (kg/m³) |
| `K_Pa` | `np.ndarray` | Bulk modulus (Pa) |
| `mu_Pa` | `np.ndarray` | Shear modulus (Pa) |
| `eta_Pa_base` | `np.ndarray` | Baseline viscosity (Pa·s) |
| `T_K` | `np.ndarray` | Temperature (K) |
| `P_MPa` | `np.ndarray` | Pressure (MPa) |
| `bulk_visc` | `np.ndarray` | Bulk viscosity (Pa·s), zeros by default |
| `phases` | `np.ndarray` | PP integer phase codes per cell |
| `changeIndices` | `np.ndarray` | Layer boundary indices into per-cell arrays |
| `n_layers` | `int` | Number of distinct layers |
| `layer_types` | `tuple` | Layer role labels |
| `region_phases` | `list[str]` | Phase string per layer (e.g. `'Ih'`, `'Ih_conv'`, `'0'`) |
| `layer_upper_radii` | `tuple` | Upper radius of each layer (m) |
| `R_body_m` | `float` | Surface radius (m) |
| `Mtot_kg` | `float` | Total body mass (kg) |
| `omega` | `float` | Mean motion (rad/s) |
| `eccentricity` | `float` | Orbital eccentricity |
| `host_mass` | `float` | Host planet mass (kg) |
| `a_m` | `float` | Semi-major axis (m) |
| `CMR2` | `float` | Moment-of-inertia factor from PP |
| `Tb_K` | `float` | Basal temperature at which this structure was built |
| `rhoSil_kgm3` | `float` | Mean silicate density (kg/m³) |
| `D_iceIh_km` | `float` | Ice Ih shell thickness (km) |
| `D_iceIII_km`, `D_iceV_km`, `D_iceVI_km` | `float` | HP ice thicknesses (km) |

### Why a Tb grid

PP's layer-build is the computational bottleneck: a single run takes O(1 min).
MCMC samples Tb_K at every likelihood evaluation—O(10⁵) times per run.  The
cache converts that cost from O(10⁵ × 1 min) = months to O(9 × 1 min) = minutes
for the build, plus O(10⁵ × 1 ms) = minutes for the MCMC itself.  The
interpolation error introduced by the grid is controlled by the B-layered blend
(section 4) and by transition refinement (this section).

### What "transition" means

Between two adjacent grid points, the set of physical layers present—encoded in
`region_phases`—sometimes changes.  For example, as Tb_K rises toward the
Ice Ih–ocean melting curve, the ocean opens and an ocean layer (`'0'`) appears
in the sequence; at a lower transition, a convective sublayer in Ice Ih may
appear or collapse.  These are genuine phase-topology changes, not numerical
artifacts.

A concrete example from the Callisto NaCl cache: at Tb ≈ 255.535 K a silicate
sublayer disappears from the top of the sequence:

```
regions below: ['Sil', 'Sil', 'Sil', '0', 'Ih', 'Ih_conv', 'Ih']
regions above: ['Sil', 'Sil',        '0', 'Ih', 'Ih_conv', 'Ih']
```

The element-wise blend that worked at interior grid points is not meaningful
across this boundary: an ice layer cannot be linearly blended with a silicate
layer.  The B-layered blend's precondition (`_regions_match`) catches this and
falls back to nearest-neighbour selection instead.

### Transition refinement via bisection

`build_tb_grid_cache` (`cache_builder.py:474`) operates in two phases:

1. **Regular grid.** PP is run at each user-supplied Tb value.  Per-grid-point
   results are stored in order.
2. **Transition detection and bisection.** Adjacent pairs are compared via
   `_regions_match` (`cache_builder.py:389`).  For any pair where the layer
   set changes, `_bisect_transition` (`cache_builder.py:397`) is called.

`_bisect_transition` halves the interval repeatedly, running PP at each
midpoint and comparing `region_phases` to the low and high endpoints.  The loop
continues until `Tb_hi - Tb_lo < eps_T` (default 0.01 K) or `max_iter` (20)
is reached.  All midpoints—whether they fall on the low-side or high-side of
the transition—are inserted into the grid, increasing resolution around the
transition boundary.  If a third distinct layer set is found at a midpoint,
the function recurses into both sub-intervals (`cache_builder.py:455–467`).

The result is that every transition in the cache is flanked by two grid points
separated by less than 0.01 K.  The `transitions` list records each such pair:

```python
{
    'Tb_lo':      255.530,
    'Tb_hi':      255.540,
    'regions_lo': ['Sil', 'Sil', 'Sil', '0', 'Ih', 'Ih_conv', 'Ih'],
    'regions_hi': ['Sil', 'Sil', '0', 'Ih', 'Ih_conv', 'Ih'],
}
```

The runner does not need to read this list explicitly: the B-layered dispatch
in `apply_bottom_temperature` calls `_regions_match` on the two bracketing
structures at runtime, so it naturally falls back to nearest-neighbour whenever
the bracket straddles a transition.  The `transitions` list is available for
diagnostics and for any runner extension that wants to pre-check before
dispatching.

### Schema versioning

v2.0 caches (built before transition detection was added) lack the `transitions`
field.  The runner tolerates them: `apply_bottom_temperature`
(`forward_models.py:374`) reads the `'Tb_K_grid'` / `'structures'` format and
calls `_regions_match` at runtime, which is safe on any cache that has
`region_phases` populated.  A v2.0 cache will simply never hit the transition
fast path because no `transitions` list was pre-computed; the B-layered blend
either succeeds or falls back per the runtime check.  The string
`'v2.1'` in the `schema_version` field is informational; the runner does not
branch on it.

---

## 4. B-layered blend across the Tb grid

When `Tb_K` is sampled, the runner must reconstruct a structure at an arbitrary
Tb between two grid points.  The B-layered blend (`forward_models._blend_b_layered`,
`forward_models.py:246`) does this correctly when the two bracketing structures
share the same layer topology.

### Why naive element-wise blending is wrong

Any ocean-bearing body has layer boundaries (inner surface of Ice Ih, top of
ocean, bottom of ocean) that move with Tb.  At a lower Tb the ice shell is
thicker; at a higher Tb it is thinner.  If adjacent grid points have the same
total cell count, a naive element-wise blend averages, say, cell 1400 of s0
(which is a deep-ocean cell) with cell 1400 of s1 (which is an ice cell at the
same index due to different layer boundary positions).  The result is a
"mush" cell with density ~1010 kg/m³ and intermediate elastic moduli that
correspond to no physical material—and an incorrect CMR² derived from it.

### The B-layered blend step by step

`apply_bottom_temperature` (`forward_models.py:375`) dispatches the blend:

```python
j = np.searchsorted(Tb_grid, Tb_clamped)
s0, s1 = structs[j-1], structs[j]
w = (Tb_clamped - Tb_grid[j-1]) / (Tb_grid[j] - Tb_grid[j-1])

if _regions_match(s0, s1):
    out = _blend_b_layered(s0, s1, w)
else:
    nearer = s0 if (Tb_clamped - t0_v) <= (t1_v - Tb_clamped) else s1
    out = _copy_struct(nearer)
```

`_regions_match` checks `s0['region_phases'] == s1['region_phases']` and
`s0['n_layers'] == s1['n_layers']` (`forward_models.py:238`).  When the
regions match, `_blend_b_layered` is called with blend weight `w`.

Inside `_blend_b_layered` (`forward_models.py:246`), the algorithm is:

**Step 1: Precondition check.** Every continuous per-cell field in s0 and s1
must have length equal to `len(r_m)`.  An assertion loop checks this for each
field in `_BLEND_CONT_FIELDS = ('rho', 'K_Pa', 'mu_Pa', 'eta_Pa_base', 'T_K',
'P_MPa', 'bulk_visc')` before any computation, raising a detailed
`ValueError` naming the offending field if the cache has a padding inconsistency
(`forward_models.py:277–294`).

**Step 2: Per-layer boundary blending.** For each layer index k, extract the
layer's lower and upper boundary radii from both structures:

```python
r_lo = (1 - w) * s0_r[s0_ci[k]]     + w * s1_r[s1_ci[k]]
r_hi = (1 - w) * s0_r[s0_ci[k+1]-1] + w * s1_r[s1_ci[k+1]-1]
```

The blended layer spans `[r_lo, r_hi]` in the output structure.  This is the
physically correct treatment: the ice–ocean boundary is at radius `r_hi` of the
ice layer, and it moves linearly with Tb between the two grid values.

**Step 3: Target cell grid.** The target intra-layer cell positions are derived
from s0's normalised positions:

```python
t0 = (s0_r[s0_lo:s0_hi] - r0_lo_b) / (r0_hi_b - r0_lo_b)
out_r[s0_lo:s0_hi] = r_lo + t0 * (r_hi - r_lo)
```

s0's cell count and relative distribution are preserved; cells are simply
rescaled to the new layer extent.

**Step 4: Resampling s1's values.** For each continuous field, s1's values are
resampled onto the target normalised positions via `np.interp`, then
linearly blended with s0's values cell-wise:

```python
s1_resampled = np.interp(t0, t1_src, s1_val)
out[field][s0_lo:s0_hi] = (1 - w) * s0_val + w * s1_resampled
```

This ensures every blended cell contains a material-consistent property value:
an ice cell receives a blend of two ice-cell properties at two nearby Tb values,
not a blend of ice and ocean.

**Step 5: Discrete and scalar fields.**

- Phase codes, `changeIndices`, `n_layers`, `layer_types`, `region_phases`: copied
  from s0 (identical to s1 by the matching-regions precondition).
- Body-constant scalars (`omega`, `R_body_m`, `Mtot_kg`, etc.): copied from s0
  (Tb-invariant by construction).
- Per-Tb scalars (`CMR2`, `rhoSil_kgm3`, `D_iceIh_km`, etc.): linearly blended
  (`forward_models.py:349–354`).
- `layer_upper_radii`: element-wise linear blend of the tuple
  (`forward_models.py:357–364`).

### Fallback for transition brackets

When `_regions_match` returns False, the blend is skipped and the nearer grid
point is returned unmodified.  Because `_bisect_transition` has already
narrowed each transition interval to < 0.01 K, the maximum interpolation error
from this fallback is bounded by the change in the observable across a 0.01 K
Tb step—typically well below the observational uncertainty in CMR² (a few
parts per thousand).

---

## 5. ocean_overrides and end-to-end workflow

### ocean_overrides

`ocean_overrides` is an optional top-level key in the v2.1 BodyConfig JSON.
During the cache build, `cache_builder.build_single_structure` (`cache_builder.py:78`)
deep-copies the body's PP template module, overrides `Planet.Bulk.Tb_K`, and
then applies each entry in `ocean_overrides` via `setattr(Planet.Ocean, key,
value)` (`cache_builder.py:134–141`).  The source template module is not
mutated.

Typical use cases:
- Switch composition: `{"comp": "NaCl", "wOcean_ppt": 100.0}`
- Sweep concentration without a new template file.
- Pin any `Planet.Ocean.*` attribute that differs from the body's PP default.

### Callisto NaCl example

`callisto_nacl_andrade_8D.json` sets:

```json
"ocean_overrides": {
  "comp": "NaCl",
  "wOcean_ppt": 100.0
}
```

The PP default for Callisto uses MgSO₄.  The switch to NaCl is motivated by
EOS coverage: SeaFreeze's NaCl EOS extends to P = 5000 MPa and T = 501 K,
while the MgSO₄ EOS covers only P < 800 MPa.  Callisto's deep ocean at high
Tb values exceeds this MgSO₄ limit, forcing PP to generate slow on-the-fly
extrapolations at each cache-build grid point.  The NaCl override eliminates
that bottleneck without requiring a separate Callisto template module for each
composition variant.  Future composition sweeps (e.g., wOcean_ppt in 5–10 ppt
increments) can be run by copying this config and editing `wOcean_ppt` alone.

### End-to-end workflow

The numbered steps below constitute the standard workflow for a new body or
composition.

**1. Choose or write a v2.1 BodyConfig JSON.**

Place it under `PlanetProfile/Inference/configs/`, using the naming convention
`<body>_<composition>_<rheology>_<dim>D.json`.  Start from
`callisto_nacl_andrade_8D.json` if the body has a sampled ocean and core;
start from `test50_titan_noocean_andrade_8D.json` if it does not.  Set
`param_space`, `observables`, `structure_cache_path`, and `ocean_overrides`
(if needed).  See `configs/README.md` for the full field-by-field schema.

**2. Build the structure cache.**

```bash
python -m PlanetProfile.Inference.build_phase_c1_cache \
    --config PlanetProfile/Inference/configs/<body>.json \
    --n-grid 9 \
    --force
```

Expected wall time: 30–60 minutes.  The `--force` flag overwrites an existing
cache at `structure_cache_path`; omit it to skip if the cache already exists.
With `--n-grid 9`, PP runs at 9 regular Tb points plus any bisection midpoints
needed to resolve transitions.  A total of 12–20 PP runs is typical.

**3. Validate the config.**

```bash
python -m PlanetProfile.Inference.run_inference_cli \
    --config PlanetProfile/Inference/configs/<body>.json \
    --validate-only
```

This runs `validate_config` (`inference_core.py:347`) and reports any schema
errors, prior inconsistencies, or missing cache file without starting the
sampler.  Fix any errors before proceeding.

**4. Smoke MCMC (n_effective = 100).**

Copy the config to a smoke variant, lower `n_effective` to 100, and run:

```bash
python -m PlanetProfile.Inference.run_inference_cli \
    --config PlanetProfile/Inference/configs/<body>_smoke.json \
    --output /tmp/<body>_smoke.pkl
```

A smoke run should complete in 1–3 minutes on a single core.  Inspect the
output to confirm the sampler terminates, the posterior samples are within
prior bounds, and key observables (CMR², k₂) are in the right ballpark before
committing to a full run.

**5. Production MCMC (n_effective = 500).**

Use the original config:

```bash
python -m PlanetProfile.Inference.run_inference_cli \
    --config PlanetProfile/Inference/configs/<body>.json \
    --output results/<body>_run1.pkl
```

Expected wall time on a workstation: 10–30 minutes.  On a shared HPC node the
job can be submitted with a 2-hour wall time limit.

**6. Inspect the result.**

```python
import pickle
import numpy as np

with open("results/<body>_run1.pkl", "rb") as f:
    result = pickle.load(f)

# Posterior summary statistics
stats = result.get_summary_stats()
for i, name in enumerate(result.param_names):
    print(f"{name}: {stats['median'][i]:.4g} ± {stats['std'][i]:.2g}")

# Best-fit sample
best = result.get_best_fit()

# CMR² posterior
if result.cmr2_results is not None:
    print(f"CMR2: {np.median(result.cmr2_results):.4f} ± "
          f"{np.std(result.cmr2_results):.4f}")

# Convergence diagnostics
print(result.convergence_metrics)
```

`InferenceResult.save` and `InferenceResult.load` (`inference_core.py:231,244`)
wrap `pickle.dump` / `pickle.load` with atomic writes and directory creation.

### Cross-references

- `configs/README.md` — complete JSON schema field-by-field, observable keys,
  derived_params semantics, and roadmap for Stage 2/3 features.
- `MCMC_INFERENCE_GUIDE.md` (repo root) — narrative walkthrough of the
  Titan-specific Test41–50 sequence; context for the no-ocean vs. ocean-bearing
  model comparison and the progression from fixed-Tb to sampled-Tb configs.

## Wedge plot (PP canonical)

```python
plot_structure_wedge_pp(
    result, grid_cache, output_path,
    *,
    planet_template_module=None,
    param_overrides=None,
    use='median',
    sample_index=None,
    strict_validate=False,
    fig_format=None,
)
```

`mcmc_plots.plot_structure_wedge_pp(...)` re-runs PlanetProfile's full forward
model at a chosen posterior point (`use='median'`, `'best_fit'`, or
`'sample'`) and saves the canonical `PlotWedge` figure to `output_path`.
This bypasses the cached-structure interpolation used during sampling and
produces a manuscript-quality wedge that exactly reflects what PP itself
would draw for the inferred parameters.

Three reproducibility-relevant pitfalls are handled internally:

1. **`FigMisc.figFormat` controls the saved byte stream.**  PP's `PlotWedge`
   calls `fig.savefig(path, format=FigMisc.figFormat, ...)`, so the format
   on disk is governed by the PP-wide config (default `'pdf'`), not by the
   extension on the filename.  By default (`fig_format=None`)
   `plot_structure_wedge_pp` honors whatever `FigMisc.figFormat` is currently
   set to and rewrites `output_path`'s extension to match — so the filename
   always agrees with the bytes.  Pass `fig_format='png'` (or `'pdf'`,
   `'svg'`, `'eps'`, `'jpg'`, `'tif'`) to request a specific format
   regardless of the PP-wide setting; the wrapper temporarily sets
   `FigMisc.figFormat` and `FigMisc.xtn` for the duration of the call and
   restores them in a `finally` block.
2. **`Params.SKIP_PLOTS` leaks from the forward run.**  `_run_pp_with_overrides`
   sets `SKIP_PLOTS=True` so the forward call is silent; the wrapper resets it
   to `False` before invoking `PlotWedge`.
3. **`tight_layout` chokes on zero-width Wedge patches.**  Both `tight_layout`
   and `savefig` are locally monkey-patched to inject `bbox_inches='tight'`
   and to swallow the `ValueError` matplotlib raises on degenerate patches.
   Original methods are restored in a `finally` block.

Wedge fill colors come from PP's `Color` config via the `_wedge_color_map()`
helper rather than hardcoded fallbacks, so the inferred-structure wedge
matches the style of any other PP wedge for the same body.

Example:

```python
plot_structure_wedge_pp(
    result, grid_cache,
    "mcmc_results/Titan/Test48_andrade_yao2014/wedge_pp.png",
    use='median',
    fig_format='png',
)
```
