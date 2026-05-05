# Tidal Love Numbers in PlanetProfile

PlanetProfile computes tidal Love numbers (k₂, h₂, l₂) from the 1D interior structure model.
Love numbers quantify how much a body deforms under tidal forcing — they are essential for
interpreting spacecraft gravity measurements and for estimating interior dissipation.

Two solvers are available: **PyALMA3** (default) and **TidalPy**.

---

## 1. Overview

Given a fully computed interior structure, PlanetProfile assembles a radial profile of density,
elastic moduli, and viscosity, then solves the linearized viscoelastic tidal deformation problem
to obtain the complex Love numbers:

- **k₂**: potential Love number (relates tidal potential to gravitational response)
- **h₂**: radial displacement Love number
- **l₂**: tangential displacement Love number
- **delta = 1 + k₂ - h₂**: auxiliary combination

Both the amplitude (`Planet.Gravity.kAmp`, `hAmp`, `lAmp`, `deltaAmp`) and phase delay
(`Planet.Gravity.kPhase`, `hPhase`, `lPhase`, `deltaPhase`, in degrees) are stored.

---

## 2. Enabling Love Number Calculations

Love numbers are computed automatically when all four of these conditions are met
in `configPP.py`:

```python
Params.CALC_NEW_GRAVITY = True   # Recalculate Love numbers (not reload from file)
Params.CALC_SEISMIC     = True   # Seismic properties required as input
Params.CALC_VISCOSITY   = True   # Viscosity profile required as input
Params.SKIP_INNER       = False  # Must not skip interior calculations
```

To reload previously saved results without recalculating:

```python
Params.CALC_NEW_GRAVITY = False  # Reloads from gravityParameters.txt if it exists
```

---

## 3. PyALMA3 (Default Solver)

PyALMA3 uses a propagator matrix method for viscoelastic radial deformation. It treats the
planet as a stack of Maxwell or Andrade layers and uses Laplace-transform techniques to
produce complex Love numbers at a specified tidal period.

**Installation:**
```bash
pip install PyALMA3
```

**No extra configuration is needed** — PyALMA3 runs whenever Love numbers are enabled.
The default rheology assigned per layer type is controlled automatically by `SetupGravity`.

**Andrade exponent and zeta (PyALMA3):**

```python
# In PPBody.py or a configPP-style script
Planet.Gravity.andradExponent = 0.2       # Andrade creep exponent (default: 0.2)

# Uniform zeta for all layers:
Planet.Gravity.andrade_zeta = 1.0         # Default; <1 amplifies transient creep

# Or per-phase zeta (Petricca 2025 range: 0.01 – 100):
Planet.Gravity.andrade_zeta = {
    'Ih':  1.0,
    'III': 0.1,
    'V':   0.1,
    'VI':  0.1,
    'Sil': 1.0,
}
```

**Reference:** Spada et al. (PyALMA / ALMA3 series).

---

## 4. TidalPy (Alternative Solver)

TidalPy uses a radial solver approach and supports multiple rheology models per layer.
It additionally computes per-layer volumetric tidal heating, making it the preferred backend
for detailed dissipation analysis.

**Installation:**
```bash
pip install TidalPy
```

**Enable TidalPy backend** in `configPP.py` or at the top of your run script:

```python
Params.Gravity.backend = 'tidalpy'
```

### Rheology options

TidalPy supports four shear rheologies. Specify them per layer region using
`Params.Gravity.rheology_models`, a dict keyed by layer-type string:

| Key string | Represents |
|------------|------------|
| `'Ih'`, `'Ih_conv'` | Surface ice Ih (conductive / convective) |
| `'III'`, `'III_conv'` | High-pressure ice III |
| `'V'`, `'V_conv'` | High-pressure ice V |
| `'VI'` | High-pressure ice VI |
| `'Clath'`, `'Clath_conv'` | Clathrate layers |
| `'Sil'` | Silicate mantle |
| `'Fe'` | Iron core |
| `'0'` | Ocean / fluid (always Newton) |

Supported rheology values per key:

| Value string | Rheology |
|---|---|
| `'andrade'` | Andrade power-law (uses `andradExponent` and `andrade_zeta`) |
| `'maxwell'` | Maxwell viscoelastic |
| `'elastic'` | Purely elastic (no dissipation) |
| `'newton'` | Purely viscous (fluid) |

**Example — Andrade ice, elastic core:**
```python
Params.Gravity.backend = 'tidalpy'
Params.Gravity.rheology_models = {
    '0':          'newton',
    'Ih':         'andrade',
    'Ih_conv':    'andrade',
    'III':        'andrade',
    'III_conv':   'andrade',
    'V':          'andrade',
    'V_conv':     'andrade',
    'VI':         'andrade',
    'Sil':        'andrade',
    'Fe':         'elastic',
    'Clath':      'elastic',
    'Clath_conv': 'maxwell',
}
Planet.Gravity.andradExponent = 0.2
Planet.Gravity.andrade_zeta   = {'Ih': 1.0, 'III': 0.1, 'V': 0.1, 'VI': 0.1, 'Sil': 1.0}
```

**Example — Maxwell rheology everywhere:**
```python
Params.Gravity.rheology_models = {
    '0': 'newton', 'Ih': 'maxwell', 'Ih_conv': 'maxwell',
    'III': 'maxwell', 'III_conv': 'maxwell',
    'V': 'maxwell', 'V_conv': 'maxwell',
    'VI': 'maxwell', 'Sil': 'maxwell',
    'Fe': 'elastic', 'Clath': 'elastic',
}
```

**Reference:** Renaud & Henning (2018), *JGR: Planets*, doi:10.1029/2017JE005450.

---

## 5. Configuration Example

A minimal `configPP.py` block for a tidal Love number run with TidalPy:

```python
Params.CALC_NEW_GRAVITY = True
Params.CALC_SEISMIC     = True
Params.CALC_VISCOSITY   = True
Params.SKIP_INNER       = False
Params.Gravity.backend  = 'tidalpy'
```

Required orbital parameters in `PPBody.py` (needed for tidal heating diagnostics):

```python
Planet.Bulk.eccentricity      = 0.0094    # Orbital eccentricity (dimensionless)
Planet.Bulk.meanMotion_radps  = 2.048e-5  # Mean orbital motion in rad/s
                                           # = 2*pi / orbital_period_s
```

> `Planet.Bulk.meanMotion_radps` is required by TidalPy; it is also used by the tidal
> heating consistency check in both backends. If it is `None`, tidal heating diagnostics
> are skipped silently.

---

## 6. Output Fields

After a successful run, results are stored in `Planet.Gravity`:

| Attribute | Description |
|---|---|
| `Planet.Gravity.k` | Complex k₂ Love number |
| `Planet.Gravity.h` | Complex h₂ Love number |
| `Planet.Gravity.l` | Complex l₂ Love number |
| `Planet.Gravity.delta` | Complex delta = 1 + k₂ - h₂ |
| `Planet.Gravity.kAmp` | |k₂| (amplitude) |
| `Planet.Gravity.hAmp` | |h₂| (amplitude) |
| `Planet.Gravity.lAmp` | |l₂| (amplitude) |
| `Planet.Gravity.deltaAmp` | |delta| (amplitude) |
| `Planet.Gravity.kPhase` | Phase delay of k₂ (degrees, sign convention: positive = lag) |
| `Planet.Gravity.tidalpy_Htidal_Wm3` | Per-radial-node volumetric heating, shape (N,) W/m³ (TidalPy only) |
| `Planet.Gravity.tidalpy_Htidal_perPhase_W` | Dict: phase label → total power in W (TidalPy only) |
| `Planet.Gravity.tidalpy_Htidal_perPhase_Wm3` | Dict: phase label → volume-averaged W/m³ (TidalPy only) |

Results are also written to `<body>/figures/<body>_gravityParameters.txt`.

### Tidal heating from k₂

The total tidal dissipation power implied by Im(k₂) is:

    E_dot = (21/2) * n^5 * R^5 * e^2 * |Im(k₂)| / G

PlanetProfile evaluates this automatically and stores the implied ice-layer heating rate
in `Planet.Ocean.HtidalIce_Wm3_computed` for comparison with the user-specified value.

To have PlanetProfile override the user-specified `Ocean.HtidalIce_Wm3` with the
self-consistent TidalPy value:

```python
Planet.Do.DO_SELF_CONSISTENT_HTIDAL = True  # Requires Params.Gravity.backend = 'tidalpy'
```

---

## 7. Solver Comparison

| Feature | PyALMA3 | TidalPy |
|---|---|---|
| Default backend | Yes | No (opt-in) |
| Rheology options | Andrade, Maxwell (via layer type) | Andrade, Maxwell, elastic, Newton |
| Per-layer rheology assignment | Via `andrade_zeta` dict | Full `rheology_models` dict |
| Per-layer heating output | No | Yes |
| Self-consistent Htidal update | No | Yes |
| Speed | Fast | Moderate |
| Best for | Quick k₂ / h₂ estimate | Detailed dissipation analysis, MCMC |

---

## 8. Connection to MCMC Inference

When using Bayesian inference to constrain rheology from observed Love numbers
(e.g., Europa Clipper k₂ measurement), TidalPy is required as the forward model backend.
The MCMC sampler varies parameters such as `andradExponent`, `andrade_zeta`, and
`etaRock_Pas`, calling TidalPy for each sample to compute the predicted k₂.

See [MCMC_INFERENCE_GUIDE.md](MCMC_INFERENCE_GUIDE.md) for full details on
setting up inference runs (Tests 41–44).

---

## 9. References

- **PyALMA3**: Spada, G. et al. — ALMA3 propagator matrix code for viscoelastic Earth/planetary models.
- **TidalPy**: Renaud, J. P. & Henning, W. G. (2018). Increased tidal dissipation using advanced rheological models. *JGR: Planets*, 123(9), 2310–2337. https://doi.org/10.1029/2017JE005450
- **Tidal theory**: Tobie, G., Mocquet, A., & Sotin, C. (2005). Tidal dissipation within large icy satellites. *Icarus*, 177(2), 534–549.
- **k₂ from Cassini/Clipper**: Iess, L. et al. (2012, 2014); Mazarico, E. et al. (Europa Clipper mission).
