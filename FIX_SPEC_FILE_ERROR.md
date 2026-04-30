# Fix for AttributeError: 'ParamsStruct' object has no attribute 'SPEC_FILE'

## Problem

When running Ganymede (or any body) through PlanetProfileApp, users encountered this error:
```
AttributeError: 'ParamsStruct' object has no attribute 'SPEC_FILE'
```

## Root Cause

The `SPEC_FILE` attribute was being set conditionally in `PlanetProfile/Main.py` (line 506) when specific files were passed to the `run()` function:

```python
if fNames is not None:
    fNames = list(fNames)
    Params.SPEC_FILE = True  # Only set here!
```

However, this attribute was **not initialized** in the default configuration (`PlanetProfile/defaultConfig.py`).

When the GUI runs PlanetProfile, it doesn't always go through the conditional logic that sets `SPEC_FILE`, so when the code later tries to check this attribute (lines 2043 and 2074 in Main.py), it doesn't exist and throws an `AttributeError`.

```python
# Line 2043 in Main.py - tries to check SPEC_FILE
if Params.SPEC_FILE and fNames is not None:
    # ... code ...

# Line 2074 in Main.py - also checks SPEC_FILE
if (Params.RUN_ALL_PROFILES and Params.COMPARE) or Params.SPEC_FILE:
    # ... code ...
```

## Solution

Added initialization of `SPEC_FILE` to `PlanetProfile/defaultConfig.py`:

```python
Params.SPEC_FILE = False  # Whether specific files were passed to run() function
```

This ensures the attribute always exists with a sensible default value, whether or not it gets set to `True` later.

## File Modified

**PlanetProfile/defaultConfig.py** (line 46):
```python
Params.RUN_ALL_PROFILES = False  # Whether to run all PPBody.py files for the named body and plot together
Params.COMPARE =          False  # Whether to plot each new run against other runs from the same body
Params.SPEC_FILE =        False  # Whether specific files were passed to run() function  [NEW]
```

## Testing

Created test script `test_ganymede_fix.py` to verify:
- ✅ `SPEC_FILE` attribute exists in default Params
- ✅ Ganymede profile loads without errors
- ✅ Main module imports successfully
- ✅ All required Params attributes present

Test results:
```
✅ All tests passed! Issue is fixed.
```

## Impact

- **CLI**: No change - already worked because it goes through the `run()` function
- **GUI**: Now works correctly - `SPEC_FILE` always initialized
- **Backwards compatible**: Default value of `False` matches expected behavior

## Related Attributes

While investigating, confirmed these attributes are properly initialized:
- `RUN_ALL_PROFILES` ✓
- `COMPARE` ✓
- `SKIP_INDUCTION` ✓
- `SKIP_PLOTS` ✓
- `TIME_AND_DATE_LABEL` ✓
- `NO_SAVEFILE` ✓
- `DO_PARALLEL` ✓
- `CALC_NEW` ✓
- `CALC_NEW_INDUCT` ✓

## Verification

To verify the fix works:
```bash
python test_ganymede_fix.py
```

To test Ganymede specifically through the app:
```bash
streamlit run PlanetProfileApp/PlanetProfileApp.py
# Then select Ganymede and run with defaults
```

---

**Status**: ✅ Fixed
**Date**: 2026-03-10
**Affects**: PlanetProfile 3.1.0+
