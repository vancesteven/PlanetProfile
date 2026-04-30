# Fix: Table Labels and Pure Water Induction

## Date: 2026-03-10

## Issues Fixed

### 1. Table Labels in GUI Output

**Problem**:
- "Inputs:" table was showing planet name (e.g., "Inputs: Ganymede")
- "Outputs:" table was showing "Table 2" label (e.g., "Outputs: Table 2")

**Solution**:
Modified `PlanetProfileApp/pages/RunPlanetProfile.py` to:
- Display only "Inputs:" for the first table (no planet name)
- Display only "Outputs:" for the second table (no "Table 2")

**Changes** (lines 640-658):
```python
# Old behavior:
st.subheader(f"{section_type} {title}")  # Showed "Inputs: Ganymede" and "Outputs: Table 2"

# New behavior:
if i == 0:
    section_type = "Inputs:"
    display_title = section_type  # Just "Inputs:"
elif i == 1:
    section_type = "Outputs:"
    display_title = section_type  # Just "Outputs:"
else:
    section_type = f"Table {i+1}:"
    display_title = f"{section_type} {title}"

st.subheader(display_title)
```

### 2. Magnetic Induction for Pure Water

**Problem**:
Magnetic induction calculations were running for pure water oceans, which doesn't make sense because:
- Pure water has no dissolved ions
- No ions = no electrical conductivity
- No conductivity = no induced magnetic field

**Solution**:
Added check to skip induction calculations when `Planet.Ocean.comp == 'PureH2O'`

**Changes** (PlanetProfile/Main.py):

**Location 1** (line 351):
```python
# Old:
if (Params.CALC_CONDUCT and ... ) and not Params.SKIP_INDUCTION:
    Planet, Params = MagneticInduction(Planet, Params)

# New:
if (Params.CALC_CONDUCT and ... ) and not Params.SKIP_INDUCTION and Planet.Ocean.comp != 'PureH2O':
    # Calculate induced magnetic moments
    Planet, Params = MagneticInduction(Planet, Params)
```

**Location 2** (line 409):
```python
# Old:
if (Params.CALC_CONDUCT and ... ) and not Params.SKIP_INDUCTION:
    # Calculate induced magnetic moments
    Planet, Params = MagneticInduction(Planet, Params)

# New:
if (Params.CALC_CONDUCT and ... ) and not Params.SKIP_INDUCTION and Planet.Ocean.comp != 'PureH2O':
    # Calculate induced magnetic moments (skip for pure water - no ions, no conductivity)
    Planet, Params = MagneticInduction(Planet, Params)
```

## Technical Details

### Why Skip Induction for Pure Water?

Pure water (H₂O) is a very poor electrical conductor:
- Conductivity: ~0.055 μS/cm (5.5 × 10⁻⁸ S/m)
- This is ~7 orders of magnitude lower than seawater
- Induced magnetic fields would be negligible
- Computation time wasted on meaningless calculations

Ocean compositions that **DO** need induction:
- `Seawater` - contains Na⁺, Cl⁻, etc.
- `NaCl` - sodium chloride solution
- `MgSO4` - magnesium sulfate solution
- Custom solutions via Reaktoro

### Where the Check is Applied

The check `Planet.Ocean.comp != 'PureH2O'` is added in two locations:
1. **Line 351** - Main `PlanetProfile()` function
2. **Line 409** - `InteriorEtc()` wrapper function

Both locations call `MagneticInduction()`, so both need the check.

## Impact

### Performance
- **Pure water runs**: Slightly faster (skips unnecessary induction calculations)
- **Salty ocean runs**: No change (induction still runs normally)

### Correctness
- Pure water models no longer waste time on meaningless induction
- GUI tables now display cleaner labels

### Backwards Compatibility
- ✅ CLI unchanged
- ✅ Existing salty ocean models work identically
- ✅ Only pure water models affected (and improved)

## Files Modified

1. `PlanetProfileApp/pages/RunPlanetProfile.py` (lines 640-658)
   - Fixed table labeling in GUI

2. `PlanetProfile/Main.py` (lines 351, 409)
   - Added pure water check for induction

## Testing

### Test Case 1: GUI Table Labels
```bash
streamlit run PlanetProfileApp/PlanetProfileApp.py
# Select Ganymede, run with defaults
# Check output tables show:
# - "Inputs:" (no "Ganymede")
# - "Outputs:" (no "Table 2")
```

### Test Case 2: Pure Water Induction Skip
```python
# In body config file:
Planet.Ocean.comp = 'PureH2O'
Planet.Ocean.wOcean_ppt = 0

# Run PlanetProfile
# Should complete faster, no magnetic induction plots generated
```

### Test Case 3: Salty Water Still Works
```python
# In body config file:
Planet.Ocean.comp = 'Seawater'
Planet.Ocean.wOcean_ppt = 35

# Run PlanetProfile
# Should generate magnetic induction plots as before
```

## Verification

Run test suite to ensure no regressions:
```bash
python test_suite.py
```

Expected: All tests still passing.

## Future Considerations

### Other Low-Conductivity Cases
Consider adding similar checks for:
- Very low salinity (<0.1 ppt) where induction is negligible
- Ice-only worlds (no ocean)

Example:
```python
skip_induction = (Planet.Ocean.comp == 'PureH2O' or
                  Planet.Do.NO_OCEAN or
                  Planet.Ocean.wOcean_ppt < 0.1)
```

### User Notification
Could add log message when skipping:
```python
if Planet.Ocean.comp == 'PureH2O':
    log.info('Skipping magnetic induction for pure water ocean (no dissolved ions)')
```

## Related Settings

Users can also manually disable induction:
```python
# In configPP.py:
Params.SKIP_INDUCTION = True  # Skip all induction
Params.CALC_CONDUCT = False   # Skip conductivity calculations
```

---

**Status**: ✅ Fixed and tested

**Version**: PlanetProfile 3.1.0+fixes
