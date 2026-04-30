# Julia ALMA3 Integration Assessment

## Current PyALMA3 Usage

PlanetProfile calls exactly two functions from `alma` (PyALMA3):

1. **`build_model(r, rho, mu, vis, rheology, params, ndigits, verbose, parallel)`**
   - Inputs: radial, density, shear modulus, viscosity arrays + rheology type per layer + Andrade/Burgers parameters
   - Returns: `model_params` (opaque structure for love_numbers)

2. **`love_numbers(harmonic_degrees, time_log_kyrs, loading_type, time_history_function, tau, model_params, output_type, gorder, verbose, parallel)`**
   - Returns: `(h, l, k)` complex Love number arrays, shape `(n_degrees, n_times)`

Called from `Gravity.py:24-33`. Output is `h, l, k` complex arrays converted to amplitude+phase.

## Julia ALMA.jl

The Julia implementation lives at `github.com/drsaikirant88/PyALMA3/tree/main/JULIA`. 

### Integration approach

Create `PlanetProfile/Gravity/AlmaJulia.py`:
```python
from juliacall import Main as jl
jl.include("path/to/ALMA.jl")

def build_model_julia(...):
    # Convert numpy arrays to Julia arrays
    # Call jl.build_model(...)
    return model_params

def love_numbers_julia(...):
    # Call jl.love_numbers(...)
    # Convert result back to numpy complex arrays
    return h, l, k
```

Dispatch in `Gravity.py`:
```python
if Params.Gravity.USE_JULIA_ALMA:
    from PlanetProfile.Gravity.AlmaJulia import build_model_julia as build_model, ...
else:
    from alma import build_model, love_numbers
```

### Requirements

- Julia >= 1.6 installed
- `juliacall` Python package (`pip install juliacall`)
- ALMA.jl source file (vendored or fetched)

### Risk assessment

| Factor | Impact |
|--------|--------|
| New dependency (Julia runtime) | High — users must install Julia separately |
| API compatibility | Medium — ALMA.jl may have different function signatures |
| Precision matching | Medium — Julia's arbitrary precision via `ArbNumerics.jl` vs Python's `mpmath` |
| Startup overhead | Medium — Julia JIT compilation on first call (~30s) |
| Maintenance burden | Low — single-file wrapper, optional backend |

### Recommendation

**Defer until benchmarking justifies it.** The PyALMA3 call is typically ~1-5 seconds for a single (degree-2, single-time) calculation. Julia's advantage is mainly for:
- Multi-degree calculations (`harmonic_degrees = [2, 3, 4, ...]`)
- Large exploreograms where love_numbers is called thousands of times
- The JIT startup cost (~30s) dominates for single runs

**If pursued:**
1. Add `Params.Gravity.USE_JULIA_ALMA = False` to `defaultConfigGravity.py`
2. Create `AlmaJulia.py` with try/except import of juliacall
3. Vendor or auto-download ALMA.jl source
4. Benchmark against PyALMA3 for Europa degree-2 and multi-degree cases
5. Verify Love numbers match to 6+ significant figures

### Config flag (ready to add)

```python
GravityParams.USE_JULIA_ALMA = False  # Use Julia ALMA.jl backend for Love number calculations (requires Julia + juliacall)
```
